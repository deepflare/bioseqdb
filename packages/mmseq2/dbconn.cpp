#include "dbconn.h"
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <exception>
#include <cstdlib>

DB::DBconn::DBconn(const std::string &tableName, const std::string &columnName, bool allTargets, const common::InputParams::Vec64Ptr &tIds)
{
    this->columnName = columnName;
    this->tableName = tableName;
    this->kmerHitsQueryPrefix = "SELECT kmer, starting_position, seq_id FROM ";
    this->kmerHitsQueryPrefix
        .append(this->tableName)
        .append("_")
        .append(this->columnName)
        .append("__index")
        .append(" WHERE kmer IN (");
    if (allTargets == false) {
        this->kmerHitsQuerySuffix = ") AND seq_id IN (";
        for (auto i = 0; i < tIds->size(); i++) {
            this->kmerHitsQuerySuffix.append(std::to_string((*tIds)[i]));
            if (i + 1 < tIds->size())
                this->kmerHitsQuerySuffix.append(",");
        }
        this->kmerHitsQuerySuffix.append(");");
    }
    else {
        this->kmerHitsQuerySuffix = ");";
    }

	/* A list of possible environment variables*/
	const char *env_var[5] = {
        "PG_HOST",
        "PG_PORT",
        "PG_DBNAME",
        "PG_USER",
        "PG_PASSWORD"
    };

	/* A list of keys passed to libpq */
	const char *conn_str_var[5] = {
        "host",
        "port",
        "dbname",
        "user",
        "password"
    };

	/* A list of default values for keys passed to libpq */
	const char *default_value[5] = {
        "localhost",
        "5432",
        "bioseqdb",
        "postgres",
        "password"
    };

    /* Example connection string: host=localhost port=5433 dbname=bioseqdb user=postgres password=postgres */
    std::string connection_string; 
    char *env_val[5];

	for(int i = 0; i < 5; i++) {
		/* Getting environment value if exists */
		env_val[i] = getenv(env_var[i]);
        
        /* Append string 'key=' */
        connection_string += conn_str_var[i];
        connection_string += "=";

		if (env_val[i] != NULL) {
            connection_string += env_val[i];
        }
		else {
            connection_string += default_value[i];
        }

        connection_string += " ";
	}

    this->connection = PQconnectdb(connection_string.c_str());

    if (PQstatus(this->connection) != CONNECTION_OK) {
        fprintf(stderr, "%s", PQerrorMessage(this->connection));

        PQfinish(this->connection);
        exit(1);
    }
}

void DB::DBconn::GetIthIndex(std::string kmer, uint32_t i, uint64_t *target_id, uint32_t *position)
{
    std::string getIndexQuery =
        "SELECT starting_position, seq_id FROM " +
        this->tableName + "_" + this->columnName + "__index" +
        " " + "WHERE kmer=\'" + kmer + "\'" +
        "OFFSET " + std::to_string(i) + "LIMIT 1;";

    PGresult *res = PQexec(connection, getIndexQuery.c_str());

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        throw std::invalid_argument(PQerrorMessage(connection));
    }

    int starting_position_fnum = PQfnumber(res, "starting_position");
    int seq_id_fnum = PQfnumber(res, "seq_id");

    /* TODO: do something when no value is returned */
    if (PQntuples(res) != 1) {
        throw std::invalid_argument("No ith index exists");
    }

    char *starting_position = PQgetvalue(res, 0, starting_position_fnum);
    char *seq_id = PQgetvalue(res, 0, seq_id_fnum);

    *position = strtoul(starting_position, NULL, 0);
    *target_id = strtoull(seq_id, NULL, 0);

    PQclear(res);
}


void DB::DBconn::CloseConnection()
{
    PQfinish(this->connection);
}

std::shared_ptr<std::string> DB::DBconn::GetTargetById(uint64_t id)
{
    std::string getTargetQuery =
        "SELECT " + this->columnName + " FROM " +
        this->tableName + " WHERE id=" + std::to_string(id) + ";";

    PGresult *res = PQexec(connection, getTargetQuery.c_str());

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        throw std::invalid_argument(PQerrorMessage(connection));
    }

    int target_seq_fnum = PQfnumber(res, this->columnName.c_str());

    /* TODO: do something when no value is returned */
    if (PQntuples(res) == 0) {
        throw std::invalid_argument("No ith target exists");
    }
    else if (PQntuples(res) > 1) {
        throw std::invalid_argument("Target id isn't unique");
    }

    char *target_seq = PQgetvalue(res, 0, target_seq_fnum);
    auto target_ptr = std::make_shared<std::string>(target_seq);
    PQclear(res);
    
    return target_ptr;
}


void DB::DBconn::GetSimKMersHits(common::SimKMersPtr &simKMersPtr, common::SimKMersHitsPtr &simKMersHitsPtr)
{
    if (simKMersPtr.get()->empty()) {
        return;
    }

    std::string query;
    query.reserve(this->kmerHitsQueryPrefix.size() + 10 * simKMersPtr.get()->size() - 1 + this->kmerHitsQuerySuffix.size());
    query.append(this->kmerHitsQueryPrefix);

    /** kMer values */
    for (uint32_t i = 0; i < simKMersPtr.get()->size(); i++)
    {
        query.append("\'");
        query.append(simKMersPtr.get()->at(i));
        query.append("\'");

        if (i < simKMersPtr.get()->size() - 1)
        {
            query.append(",");
        }
    }

    query.append(this->kmerHitsQuerySuffix);

    PGresult *res = PQexec(connection, query.c_str());

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        throw std::invalid_argument(PQerrorMessage(connection));
    }


    int kmer_fnum = PQfnumber(res, "kmer");
    int starting_position_fnum = PQfnumber(res, "starting_position");
    int seq_id_fnum = PQfnumber(res, "seq_id");

    auto respSize = PQntuples(res);

    common::SimKMersHits simKMersHits;

    for (int i = 0; i < respSize; i++) {
        char *kmer = PQgetvalue(res, i, kmer_fnum);
        char *starting_position = PQgetvalue(res, i, starting_position_fnum);
        char *seq_id = PQgetvalue(res, i, seq_id_fnum);

        simKMersHits.push_back({(std::string)kmer, {strtoull(seq_id, NULL, 0), strtoul(starting_position, NULL, 0)}});
    }

    *simKMersHitsPtr = simKMersHits;

    PQclear(res);
}
