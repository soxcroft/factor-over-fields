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

# files TODO add euclid
EXES = testberlekamp testeuclid

BINDIR = ../bin
LOCALBIN = ~/.local/bin

# RULES

# executables

testberlekamp: testberlekamp.c euclid.o berlekamp.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

testeuclid: testeuclid.c euclid.o | $(BINDIR)
	$(COMPILE) -o $(BINDIR)/$@ $^

# units

berlekamp.o: berlekamp.c euclid.h berlekamp.h
	$(COMPILE) -c $<

euclid.o: euclid.c euclid.h
	$(COMPILE) -c $<

# PHONY TARGETS

# all: TODO

clean:
	$(RM) $(foreach EXEFILE, $(EXES), $(BINDIR)/$(EXEFILE))
	$(RM) *.o
	$(RM) -rf $(BINDIR)/*.dSYM
