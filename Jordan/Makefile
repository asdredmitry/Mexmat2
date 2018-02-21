all : main

help.o : help.cpp 
	mpic++ -c help.cpp -o help.o
sol.o : sol.cpp 
	mpic++ -c sol.cpp -o sol.o
main.o : main.cpp 
	mpic++ -c main.cpp -o main.o
main : main.o sol.o help.o
	mpic++ main.o sol.o help.o -o main
	
clean : 
	rm -rf main.o sol.o help.o main