# BioSeqDB

**Note: This project is work in progress (very early stage), not intended for production use**

## Development

The relevant files and directories are `bioseqdb_pg/` and `CMakeLists.txt`. To work on the extension, install Postgres (including additional libraries, `apt install postgresql postgresql-server-dev-all` on Ubuntu, `pacman -S postgresql` on Arch) and CMake, and set up a Postgres instance (see [Ubuntu docs](https://ubuntu.com/server/docs/databases-postgresql) or [Arch docs](https://wiki.archlinux.org/title/PostgreSQL)).

Also, you need to compile and install SeqLib, which is used as an intermediate layer for BWA implementation. Follow the instructions in the [.github/workflows/ci.yml](.github/workflows/ci.yml) file, which should work out exactly on Ubuntu. On Arch, you may need to delete `const uint8_t rle_auxtab[8];` lines in `bwa/` and `fermi-lite/` dependencies inside `SeqLib/` to fix mutiple symbol definitions errors.

To build and install the extension, create a `build/` directory and run `cmake ..` from it. You can now build the extension by running the `make` command in the build directory, and install it with `sudo make install`. A typical development flow is running `make && sudo make install && sudo systemctl restart postgresql`. The entire process should take about a second.

After first installing the extension, you need to run `CREATE EXTENSION bioseqdb;` to load the additional types. If you modify the definitions of any SQL functions or types, remember to drop any affected tables, `DROP EXTENSION bioseqdb CASCADE;` and repeate the `CREATE EXTENSION` command.
