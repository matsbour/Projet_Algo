all: main

main: handle_database.o smith_waterman.o myProtein.o main.cpp
	g++ -std=c++11 -Wall -o main handle_database.o myProtein.o smith_waterman.o main.cpp

handle_database.o: handle_database.h
	g++ -std=c++11 -Wall -c handle_database.h handle_database.cpp

myProtein.o: myProtein.h 
	g++ -std=c++11 -Wall -c myProtein.h myProtein.cpp

smith_waterman.o: smith_waterman.h 
	g++ -std=c++11 -Wall -c smith_waterman.h smith_waterman.cpp

clean:  
	rm -f main *.o
