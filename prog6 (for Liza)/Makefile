all: main.o rate.o
	clear
	mpicc main.o rate.o -o main.out -lm
	mpirun -n 1 ./main.out 100
	mpirun -n 2 ./main.out 100
	mpirun -n 3 ./main.out 100
	mpirun -n 4 ./main.out 100
	rm *.o *.out

main.o: main.c matrix.h
	mpicc -c main.c -o main.o

rate.o: rate.c matrix.h
	gcc -c rate.c -o rate.o
