CC = g++
FLAG = -std=c++11 -c -Wall

all: test iteration export heat
	$(CC) *.o -o heat_1.2
test:
	$(CC) $(FLAG) test.cpp

iteration:
	$(CC) $(FLAG) test.cpp
export:
	$(CC) $(FLAG) export.cpp
heat:
	$(CC) $(FLAG) heat_1.2.cpp
clean:
	rm *.o