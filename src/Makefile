CC=gcc
CFLAGS= -std=c99  -fPIC -Wall
CFLAGS+= -L/usr/local/lib -lpcre2-8

builder: read_sdf.o read_bng.o  eprintf.o
	$(CC) ${CFLAGS} *.o -o builder.x

read_bng.o: read_bng.c
	${CC} ${CFLAGS} -c read_bng.c

read_sdf.o: read_sdf.c
	${CC} ${CFLAGS} -c read_sdf.c

eprintc.o: eprintf.c
	${CC} ${CFLAGS} -c eprintf.c


clean:
	rm -rf *.o *.x
