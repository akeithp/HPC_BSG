CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
INCLUDES=-I./include/
SOURCES=./include/BasicCDS.cpp binarySearchP.cpp

unary:
	@echo " [BLD] Binary Searches with openMP improvements"
	@$(CPP) -fopenmp $(CPPFLAGS) $(INCLUDES) $(SOURCES) -o bSearch
