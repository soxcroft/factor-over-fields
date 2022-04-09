#
# Makefile for polynomial factorizer
#

# compiler flags
DEBUG    = -ggdb
OPTIMISE = -O0
WARNINGS = -Wall -Wextra -Wno-variadic-macros -Wno-overlength-strings -pedantic
CFLAGS   = $(DEBUG) $(OPTIMISE) $(WARNINGS)
#DFLAGS = -TODO

CC       = clang
RM       = rm -f
COMPILE  = $(CC) $(CFLAGS) $(DFLAGS)
INSTALL  = install

# files
EXES = factor testberlekamp testeuclid testlift

BINDIR = ../bin
LOCALBIN = ~/.local/bin

# RULES

# executables

factor: factor.c euclid.o berlekamp.o lift.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

testlift: testlift.c lift.o euclid.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

testberlekamp: testberlekamp.c euclid.o berlekamp.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

testeuclid: testeuclid.c euclid.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

# units

lift.o: lift.c euclid.h lift.h
	$(COMPILE) -c $<

berlekamp.o: berlekamp.c berlekamp.h euclid.h
	$(COMPILE) -c $<

euclid.o: euclid.c euclid.h
	$(COMPILE) -c $<

# PHONY TARGETS

all: factor

clean:
	$(RM) $(foreach EXEFILE, $(EXES), $(BINDIR)/$(EXEFILE))
	$(RM) *.o
	$(RM) -rf $(BINDIR)/*.dSYM
