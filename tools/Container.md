### Build
```sh
docker build . -t bioseqdb
```
### Run
Assuming you want database to be stored in `data` folder run this command:
```sh
docker run -d \
	--name bioseqdb-instance \
	-e POSTGRES_PASSWORD=mypassword \
	-e PGDATA=/var/lib/postgresql/data/pgdata \
	-v data:/var/lib/postgresql/data \
	-p 5432:5432 \
	bioseqdb
```
Connect to it using your favourite client with user `postgres` and password `mypassword` on port `5432`.
