CC=g++
CFLAGS= -fPIC
DEPS = cpp-graph.h

INC_DIR = -I. -I/home/brian/software/nauty
LIBS = /home/brian/software/nauty/nauty.a /usr/lib/x86_64-linux-gnu/libbz2.a


# Makefile syntax
# $@ = object to the left of the colon
# $^ = object to the right of the colon

cpp-enumerate: cpp-graph.o cpp-nauty.o cpp-enumerate.o
	$(CC) $(INC_DIR) -o $@ $^ $(LIBS) $(CFLAGS)

cpp-graph.o: cpp-graph.cpp cpp-graph.h
	$(CC) $(INC_DIR) -c -o $@ $< $(CFLAGS)

cpp-nauty.o: cpp-nauty.cpp cpp-graph.cpp cpp-nauty.h cpp-graph.h
	$(CC) $(INC_DIR) -c -o $@ $< $(CFLAGS)

cpp-enumerate.o: cpp-enumerate.cpp cpp-graph.h cpp-nauty.h
	$(CC) $(INC_DIR) -c -o $@ $< $(CFLAGS)




clean:
	rm -f *.o \# *~
