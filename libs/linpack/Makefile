#
# linpack library Makefile...
#
ARCH=$(shell uname -s)

BUILDDIR = .

TARGET = $(BUILDDIR)/liblinpack.a

SRC = dgedi.f dgesl.f dgefa.f  dgeco.f
OBJS = $(SRC:.f=.o)


.SUFFIXES: .f .o

include ../../Params.inc
include ../../bld/Makefile.$(ARCH)


%.o:%.f
	$(F90) $(FFLAGS) $(FIXED) $(OMP_FLAGS) -c $<

all: $(TARGET)

#
# This method of specifying the library target allows a parallel make 
# but forces a single execution of the AR command (which will fail in parallel)
#

$(TARGET) : $(SRC)
	$(MAKE) $(?:.f=.o)
	$(AR) $(ARFLAGS) $@ $(?:.f=.o)
	$(RANLIB) $@
	$(RM) $(?:.f=.o)


clean:
	$(RM) $(OBJS) $(TARGET)

realclean:
	$(RM) $(OBJS) ; $(RM) $(TARGET)

