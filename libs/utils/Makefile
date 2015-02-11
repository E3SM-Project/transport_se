.SUFFIXES: .c .o

ARCH=$(shell uname -s)

BUILDDIR = .
TARGET1 = $(BUILDDIR)/libutils.a 
TARGET2 = $(BUILDDIR)/nbits

include ../../Params.inc
include ../../bld/Makefile.$(ARCH)

OBJS = timer.o nbits.o precision.o derfc.o

all: $(TARGET1) $(TARGET2)

thread_bind.o: thread_bind.c
	$(MPCC) $(CFLAGS) $(cdefs) thread_bind.c -c
#
# Build unit testers and the library
#

$(TARGET2): nbits.o testnbits.o $(TARGET1)
	$(CC) $(CFLAGS) $^ -o $@

$(TARGET1): $(OBJS) $(aobjs)
	$(AR) $(ARFLAGS)  $@ $?
	$(RANLIB) $@

clean:
	rm -f *.o libutils.a nbits

realclean:
	rm -f *.o libutils.a 



