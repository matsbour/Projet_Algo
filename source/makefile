all: main

main: handle_database.o myProtein.o main.cpp
	g++ -std=c++11 -Wall -o main handle_database.o myProtein.o main.cpp

handle_database.o: handle_database.h
	g++ -std=c++11 -Wall -c handle_database.h handle_database.cpp

myProtein.o: myProtein.h 
	g++ -std=c++11 -Wall -c myProtein.h myProtein.cpp

clean:  
	rm -f main *.o
