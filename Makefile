CFLAGS=-Wall -g -pg -std=gnu99 
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 
LIBS = -lm -lpthread

OBJS = fd.o linalg.o plotfunc.o powermethod.o 

fd: $(OBJS)
	gcc -o fd $(OBJS) $(LDFLAGS) $(LIBS)

fd.o: fd.c hamiltonian.h linalg.h plotfunc.h powermethod.h time_mes.h
	gcc $(CFLAGS) $(INCLUDES) -c fd.c

linalg.o: linalg.c linalg.h
	gcc $(CFLAGS) $(INCLUDES) -c linalg.c

plotfunc.o: plotfunc.c plotfunc.h
	gcc $(CFLAGS) $(INCLUDES) -c plotfunc.c

powermethod.o: powermethod.c powermethod.h linalg.h
	gcc $(CFLAGS) $(INCLUDES) -c powermethod.c

clean:
	rm -f *.o fd *~
