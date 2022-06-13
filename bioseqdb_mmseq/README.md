# Building the PostgreSQL extension
To build the extension run:
```
mkdir build
cd build
cmake -D PostgreSQL_TYPE_INCLUDE_DIR=<postgres include directory> ..
make
sudo make install
```
in the pg directory. You can also add optional variables in cmake:
- MMSEQ\_HOSTNAME - "localhost" by default
- MMSEQ\_PORT - 8080 by default
- MMSEQ\_TIMEOUT - timeout in milliseconds, 30000 by default (30 seconds)

For example:
```
cmake -D PostgreSQL_TYPE_INCLUDE_DIR=/usr/include/postgresql/14/server/ -D MMSEQ_TIMEOUT=60000 ..
```
Then run ``` CREATE EXTENSION mmseq2; ``` in a PostgreSQL database.

# Building the MMSeqs2 computational service
To build the MMSeqs2 computational service run:
```
mkdir build
cd build
cmake -D PostgreSQL_TYPE_INCLUDE_DIR=<postgres include directory> ..
make
```
in the mmseq2 directory. You can also optionally set the following environment variables:
- PG\_HOST - "localhost" by default
- PG\_PORT - 5432 by default
- PG\_DBNAME - "bioseqdb" by default
- PG\_USER - "postgres" by default
- PG\_PASSWORD - "password" by default

to specify for MMSeqs2 from what database to fetch the data such as the target sequences and the indexes.

Then you can run MMSeqs2 by executing the command:
```
echo <MMSeqs2 port> | ./mmseq2
```

# Running the project as a containerized application
Another way to build the project is to use the ```ops/bin/up.sh``` shell script
which creates two docker containers - one for the PostgreSQL database
and one for the MMSeqs2 computational service.
