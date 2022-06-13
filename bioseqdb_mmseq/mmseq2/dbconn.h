#ifndef BIOSEQDB_DBCONN_H
#define BIOSEQDB_DBCONN_H

extern "C" {
#include <libpq-fe.h>
}
#include <string>
#include <memory>
#include "../common/mmseq2lib.h"

namespace DB
{
    class DBconn
    {
    public:
        DBconn(const std::string&, const std::string&, bool, const common::InputParams::Vec64Ptr&);
        void GetIthIndex(std::string, uint32_t, uint64_t *, uint32_t *);
        void GetSimKMersHits(common::SimKMersPtr&, common::SimKMersHitsPtr&);
        void CloseConnection();
        std::shared_ptr<std::string> GetTargetById(uint64_t);


    private:
        PGconn *connection;
        std::string tableName;
        std::string columnName;
        std::string kmerHitsQueryPrefix;
        std::string kmerHitsQuerySuffix;
    };
}

#endif // BIOSEQDB_DBCONN_H
