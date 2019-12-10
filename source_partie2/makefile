all: main


main: handle_database.o myProtein.o smith_waterman.o main.cpp
	g++ -std=c++11 -O2 -Wall -o main handle_database.o myProtein.o smith_waterman.o main.cpp #-O2 est une optimisation

handle_database.o: handle_database.h handle_database.cpp
	g++ -std=c++11 -O2 -Wall -c handle_database.h handle_database.cpp

myProtein.o: myProtein.h myProtein.cpp
	g++ -std=c++11 -O2 -Wall -c myProtein.h myProtein.cpp

smith_waterman.o: smith_waterman.h smith_waterman.cpp
	g++ -std=c++11 -O2 -Wall -c smith_waterman.h smith_waterman.cpp

clean:  
	rm -f main *.o *gch
