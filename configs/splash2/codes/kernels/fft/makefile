
PARMACS_DIR = ../../../../parmacs.upc.3/pthreads

TARGET = FFT
OBJS = fft.o
#CFLAGS = -O2 -Olimit 2000 -w 
#LDFLAGS = -lmpc -lm
#MACROS = ../../null_macros/c.m4.null
CFLAGS = $(OPT) -I$(PARMACS_DIR)
LDFLAGS = -non_shared -lm -L$(PARMACS_DIR) -lparmacs_pthreads -lpthread -lexc
MACROS = $(PARMACS_DIR)/c.m4.pthreads

x = *

$(TARGET): $(OBJS)
	cc $(OBJS) $(CFLAGS) -o $(TARGET) $(LDFLAGS)

.SUFFIXES:
.SUFFIXES:	.o .c .C .h .H

.H.h:
	m4 ${MACROS} $*.H > $*.h

.C.c:
	m4 $(MACROS) $*.C > $*.c

.c.o:
	cc -c $(CFLAGS) $*.c

.C.o:
	m4 $(MACROS) $*.C > $*.c
	cc -c $(CFLAGS) $*.c

clean:
	rm $(OBJS) fft.c
