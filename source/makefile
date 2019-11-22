all: main

main: handle_database myProtein main.cpp
	g++ -std=c++11 -Wall -o main handle_database.o myProtein.o main.cpp

handle_database: handle_database.h 
	g++ -std=c++11 -Wall -c handle_database.h handle_database.cpp

myProtein: myProtein.h 
	g++ -std=c++11 -Wall -c myProtein.h myProtein.cpp

clean:  
	rm -f main *.o
