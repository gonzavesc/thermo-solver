CC = g++
FLAG = -std=c++11 -c  

all: test iteration export heat
	$(CC) *.o -o heat_1.2
test:
	$(CC) $(FLAG) test.cpp

iteration:
	$(CC) $(FLAG) test.cpp
export:
	$(CC) $(FLAG) test.cpp
heat:
	$(CC) $(FLAG) heat_1.2.cpp
clean:
	rm *.o