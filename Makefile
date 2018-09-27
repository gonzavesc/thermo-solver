CC = g++
FLAG = -std=c++11 -c -Wall -O3

all: init iteration export heat
	$(CC) *.o -o heat_1.2
init:
	$(CC) $(FLAG) init.cpp

iteration:
	$(CC) $(FLAG) iteration.cpp
export:
	$(CC) $(FLAG) export.cpp
heat:
	$(CC) $(FLAG) heat_1.2.cpp
clean:
	rm *.o