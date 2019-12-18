Flag = -std=c++11 -lpthread -O2 -Wall -Ofast #-O2 est une optimisation


all: main

main: handle_database.o myProtein.o smith_waterman.o main.o
	g++ ${Flag} -o projet handle_database.o myProtein.o smith_waterman.o main.o
	
main.o:
	g++ ${Flag} -c main.cpp

handle_database.o: handle_database.h handle_database.cpp
	g++ ${Flag} -c handle_database.h handle_database.cpp

myProtein.o: myProtein.h myProtein.cpp
	g++ ${Flag} -c myProtein.h myProtein.cpp

smith_waterman.o: smith_waterman.h smith_waterman.cpp
	g++ ${Flag} -c smith_waterman.h smith_waterman.cpp


thread: handle_database.o myProtein.o smith_waterman_thread.o main_thread.o
	g++ ${Flag} -o projet_thread handle_database.o myProtein.o smith_waterman_thread.o main_thread.o

main_thread.o:
	g++ ${Flag} -c main_thread.cpp

smith_waterman_thread.o: smith_waterman_thread.h smith_waterman_thread.cpp
	g++ ${Flag} -c smith_waterman_thread.h smith_waterman_thread.cpp


clean:  
	rm -f main thread *.o *gch
