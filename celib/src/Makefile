# Makefile for CELib.
.SUFFIX: .c .o

CC = clang

CFLAGS = -std=c99 -O3 -g -Wall -ffast-math  #-fPIC

GIT_VERSION = $(shell git log -1 --format="%H")
CFLAGS += -DGITVERSION=\"$(GIT_VERSION)\"

LDFLAGS =  -lm #-lefence

OBJS = IMFTypes.o Info.o WriteIMF.o Init.o YieldsInfo.o \
	   LifeTime.o IMFFunctions.o Integral.o SNIIYields.o SNIaYields.o \
	   SNIIRate.o SNIaRate.o NSMYields.o NSMRate.o \
	   InitialMetallicity.o \
	   AGBMassLoss.o UnifiedAPIs.o \
	   RunParameters.o CheckRunParameters.o \
	   Utilities.o

.PHONY: all
all:	depend libCELib.a  test checkyields

.PHONY: clean
clean:
	rm *.o *.s *.a *.swp depend.inc libCELib.a -f 

libCELib.a: $(OBJS) 
	$(AR) -r $@ $^
	ranlib $@

test:
	$(CC) $(CFLAGS) ./main.c WriteIMF.c Utilities.c -o chtest.out -L./ -lCELib $(LDFLAGS)

checkyields:
	$(CC) $(CFLAGS) ./checkyields.c Utilities.c -o yields.out -L./ -lCELib $(LDFLAGS)

tag:
	ctags *.[c,h] 

DATE = $(shell date '+%Y%m%d%H%M')
pkg:
	tar cv *.[c,h] Makefile \
	../data \
	../test | gzip > libCELib.$(DATE).tar.gz

doc:
	platex doc/readme.tex

.PHONY: depend
depend: $(OBJS:.o=.c)
	-@ $(RM) depend.inc
	-@ for i in $^; do cpp -MM $$i | sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.c//g" >> depend.inc; done

-include depend.inc
