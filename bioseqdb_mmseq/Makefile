.PHONY : all clean

CXX=g++ -std=c++17 -O2
DIR=mmseq2/

all: mock_structures.o mmseq2.o main.o
	$(CXX) -pthread mock_structures.o mmseq2.o main.o -o main

mock_structures.o: $(DIR)mock_structures.h $(DIR)mock_structures.cpp
	$(CXX) $(DIR)mock_structures.cpp -c -o mock_structures.o

mmseq2.o: $(DIR)mock_structures.h $(DIR)mmseq2.h $(DIR)mmseq2.cpp
	$(CXX) -pthread $(DIR)mmseq2.cpp -c -o mmseq2.o

main.o: $(DIR)mock_structures.h $(DIR)mmseq2.h $(DIR)main.cpp
	$(CXX) $(DIR)main.cpp -c -o main.o

clean:
	rm -rf mock_structures.o mmseq2.o main.o main