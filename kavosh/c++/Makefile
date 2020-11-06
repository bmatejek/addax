CC=g++
CFLAGS=-I.

all: cpp-enumerate

cpp-graph.o: cpp-graph.cpp
	$(CC) -c $(CFLAGS) cpp-graph.cpp

cpp-enumerate: cpp-enumerate.o cpp-graph.o
	$(CC) -o cpp-enumerate $(CFLAGS) cpp-enumerate.o cpp-graph.o

clean:
	rm -f *.o \# *~
