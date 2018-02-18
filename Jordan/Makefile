NAME	= MPI

DEBUG	= -O3
CC	= mpicc -c
LD	= mpicc
CFLAGS	= $(DEBUG) -W -Wall
LIBS	= -lm
LDFLAGS	= $(DEBUG)

OBJS = main.o help.o sol.o

all: $(NAME)

$(NAME) : $(OBJS)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJS) $(NAME)

main.o		: main.c help.h sol.h
help.o		: help.c help.h
sol.o		: sol.c sol.h
