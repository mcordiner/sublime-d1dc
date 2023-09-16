# Makefile
# 
# This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
#
# Copyright (C) 2006-2014 Christian Brinch
# Copyright (C) 2015-2017 The LIME development team
# Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)

# Platform-dependent stuff:
include Makefile.defs

##
## Make sure to put the correct paths.
##
PREFIX  =  ${PATHTOLIME}

# Paths:
srcdir		= ${CURDIR}/src
docdir		= ${CURDIR}/doc
exampledir	= ${CURDIR}/example

ifneq (,$(wildcard ${PREFIX}/lib/.))
    LIBS += -L${PREFIX}/lib
endif
ifneq (,$(wildcard ${HOME}/lib/.))
    LIBS += -L${HOME}/lib
endif
ifneq (,$(wildcard /opt/local/lib/.))
    LIBS += -L/opt/local/lib
endif
ifneq (,$(wildcard /sw/lib/.))
    LIBS += -L/sw/lib
endif
ifneq (,$(wildcard /usr/local/lib/.))
    LIBS += -L/usr/local/lib
endif

CPPFLAGS += -I${PREFIX}/include \
            -I${PREFIX}/src \
            ${EXTRACPPFLAGS}


# Names of source files included:
include Makefile.srcs

CCFLAGS += -O3 -falign-loops=16 -fno-strict-aliasing
LDFLAGS += -lgsl -lgslcblas -lcfitsio -lncurses -lsundials_cvode -lsundials_nvecserial -lsundials_nvecmanyvector -lm

ifeq (${DOTEST},yes)
  CCFLAGS += -DTEST
#  CC += -g -Wunused -Wno-unused-value -Wformat -Wformat-security
  CC += -g -Wall
endif

ifeq (${VERBOSE},no)
  CCFLAGS += -DNOVERBOSE
endif

ifeq (${USEHDF5},yes)
  CPPFLAGS += -DUSEHDF5
  CCFLAGS += -DH5_NO_DEPRECATED_SYMBOLS
  LDFLAGS += -lhdf5_hl -lhdf5 -lz
  CORESOURCES += ${HDF5SOURCES}
  CONVSOURCES += ${HDF5SOURCES}
  COREINCLUDES += ${HDF5INCLUDES}
endif

SRCS = ${CORESOURCES} ${STDSOURCES}
INCS = ${COREINCLUDES}
OBJS = $(SRCS:.c=.o)

CONV_OBJS = $(CONVSOURCES:.c=.o)

.PHONY: all doc docclean objclean limeclean clean distclean

all:: ${TARGET}

# Implicit rules:
%.o : %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

${TARGET}: ${OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}

${OBJS} : ${INCS}
${CONV_OBJS} : ${CONVINCLUDES}

objclean::
	rm -f ${srcdir}/*.o

clean:: objclean
	rm -f *~ ${srcdir}/*~
