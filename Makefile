default: main

#mettre tous les .o
main: 
	g++ $^ -o $@

%.o: %.c
	g++ -c $< -o $@
  
clean:
	rm -f main
	rm -f *.o
