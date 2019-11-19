all: main

main: handle_database.o
	g++ -std=c++11 -Wall -o main handle_database.o main.cpp

handle_database.o: handle_database.h 
	g++ -std=c++11 -Wall -c handle_database.h handle_database.cpp

clean:  
	rm -f main *.o
