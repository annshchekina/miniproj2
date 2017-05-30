CFLAGS=-Wall -O3 -g -pg -std=gnu99 
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 
LIBS = -lm -lpthread

OBJS = fd.o hamiltonian.o linalg.o plotfunc.o powermethod.o time_mes.o

fd: $(OBJS)
	gcc -o fd $(OBJS) $(LDFLAGS) $(LIBS)

fd.o: fd.c hamiltonian.h linalg.h plotfunc.h powermethod.h time_mes.h
	gcc $(CFLAGS) $(INCLUDES) -c fd.c

hamiltonian.o: hamiltonian.c hamiltonian.h
	gcc $(CFLAGS) $(INCLUDES) -c hamiltonian.c

linalg.o: linalg.c linalg.h
	gcc $(CFLAGS) $(INCLUDES) -c linalg.c

plotfunc.o: plotfunc.c plotfunc.h
	gcc $(CFLAGS) $(INCLUDES) -c plotfunc.c

powermethod.o: powermethod.c powermethod.h linalg.h
	gcc $(CFLAGS) $(INCLUDES) -c powermethod.c

time_mes.o: time_mes.c time_mes.h
	gcc $(CFLAGS) $(INCLUDES) -c time_mes.c

clean:
	rm -f *.o fd *~
