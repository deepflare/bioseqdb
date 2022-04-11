# BioSeqDB

**Note: This project is work in progress (very early stage), not intended for production use**

## Development

The relevant files and directories are `bioseqdb_pg/` and `CMakeLists.txt`. To work on the extension, install Postgres (including additional libraries, `apt install postgresql postgresql-server-dev-all` on Ubuntu, `pacman -S postgresql` on Arch) and CMake, and set up a Postgres instance (see [Ubuntu docs](https://ubuntu.com/server/docs/databases-postgresql) or [Arch docs](https://wiki.archlinux.org/title/PostgreSQL)).

To build and install the extension, create a `build/` directory and run `cmake ..` from it. You can now build the extension by running the `make` command in the build directory, and install it with `sudo make install`. A typical development flow is running `make && sudo make install && sudo systemctl restart postgresql`. The entire process should take about a second.

After first installing the extension, you need to run `CREATE EXTENSION bioseqdb;` to load the additional types. If you modify the definitions of any SQL functions or types, remember to drop any affected tables, `DROP EXTENSION bioseqdb CASCADE;` and repeate the `CREATE EXTENSION` command.

## Building

### Requiements

First please install the following tools:

- [Python (version 3.8)](https://www.python.org/downloads/release/python-382/)
- [Poetry](https://python-poetry.org/docs/)
- [CMake](https://cmake.org/download/)
- [Docker](https://docs.docker.com/get-docker/)

### Building and running locally (Docker)

Please firstly login into the [Github Docker registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry):

```bash
    export GIT_TOKEN="your_pat_token" && \
    export GIT_USER="your_git_username" && \
    make docker-login
```

You need to specify PAT (Personal Access Token) for github and github username.
PAT [can be generated here](https://github.com/settings/tokens/new)
Please write down your PAT somewhere.

Now you can build the extension:

```bash
    make docker-build-dev
```

This command performs build inside the Docker container and runs smoke test specified in `./tests/smoke_test.sql`

To run the entire local enviornment please type:
```bash
    make docker-run-dev
```

This command spawns local build using [docker-compose](https://docs.docker.com/compose/gettingstarted/).
You can access [localhost:5050](http://localhost:5050/)
Now type any password for pgadmin.

Please select servers tab and add new server:
```
    Host:      postgres_container
    User:      postgres
    Password:  changeme
    Port:      5432 (default)
```

Now you can access the database using pgadmin interface.

### Building locally (Native)

To build the code without using containerization features please do:
```
    make build
```

### Few words about Docker building

For build Docker uses custom image that is defined by `./docker/base.Dockerfile` 
To build and release this image please do:
```bash
    make docker-login && \
    make docker-build-base && \
    make docker-push-base && \
    echo "DONE"
```

## Project layout

* `bioseq_lib` - Subproject that defines static C++ library with the Bioseq algoithms implementation (that is data source agnostic and does not rely on aly PostgreSQL internals)
* `bioseq_pg` - Postgre extension written in C++ that uses interface of `bioseq_lib` and is linked against it
* `docker` - Directory with definitions of all useful images used for releases and building
* `scripts` - Directory with handy build scripts
* `sql` - Contains code that initializes the module on installation
* `tests` - Tests cases

The project uses CMake for building and python with poetry for running tests, additional scripts etc.

