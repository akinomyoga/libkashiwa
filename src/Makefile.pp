# -*- mode:makefile-gmake -*-

all:
.PHONY: all

Makefile: Makefile.pp

#------------------------------------------------------------------------------
# Host Configurations

# default configuration
opt_isystem  := -isystem
OTFFT_PREFIX := $(HOME)/opt/otfft-6.4
EIGEN_PREFIX := $(HOME)/opt/Eigen-3.2.8

host:=unknown

ifeq ($(HOSTNAME),padparadscha)
  host:=found

  OTFFT_PREFIX   := $(HOME)/work/idt/otfft/out/bin
  CXX:=cxx
  CXXFLAGS:= -Wall -std=gnu++14 -O3 -march=native
  LDFLAGS:= -Wall -std=gnu++14 -O3 -march=native
endif

ifeq ($(USER)@$(HOSTNAME),murase@tkyntn.phys.s.u-tokyo.ac.jp)
  host:=found

  OTFFT_PREFIX :=$(HOME)/opt/otfft-6.4
  CXX:=g++-5.3.0
  CXXFLAGS:= -Wall -std=gnu++14 -O3 -march=native
  LDFLAGS:= -Wall -std=gnu++14 -O3 -march=native

  LIBMWG_PREFIX:=$(HOME)/opt/libmwg-20160527
  CXXFLAGS+= \
    -isystem $(LIBMWG_PREFIX)/include \
    -isystem $(LIBMWG_PREFIX)/include/x86_64-unknown-linux-gnu-gcc-5.3.0+default
  LDFLAGS+= \
    -L $(LIBMWG_PREFIX)/lib/x86_64-unknown-linux-gnu-gcc-5.3.0+default
endif

ifeq ($(host),unknown)
  CXX:=g++
  CXXFLAGS:= -Wall -std=gnu++0x -O3 -march=native
endif

#------------------------------------------------------------------------------

OTFFT_LDFLAGS  := -L $(OTFFT_PREFIX)/lib -Wl,-rpath,$(OTFFT_PREFIX)/lib
OTFFT_CXXFLAGS := $(opt_isystem) $(OTFFT_PREFIX)/include
OTFFT_LIBS     := -lotfft

OUTDIR := ../out
OBJDIR := $(OUTDIR)/bin
directories := $(OUTDIR) $(OBJDIR)

#
# Libraries
#

# ksh

directories += $(OBJDIR)/ksh

objectfiles += $(OBJDIR)/ksh/embedded_runge_kutta.o
-include $(OBJDIR)/ksh/embedded_runge_kutta.d
$(OBJDIR)/ksh/embedded_runge_kutta.o: ksh/embedded_runge_kutta.cpp | $(OBJDIR)/ksh
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ $<

objectfiles += $(OBJDIR)/ksh/integrator.o
-include $(OBJDIR)/ksh/integrator.d
$(OBJDIR)/ksh/integrator.o: ksh/integrator.cpp | $(OBJDIR)/ksh
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ -I ksh $<

objectfiles += $(OBJDIR)/ksh/linear_lu.o
-include $(OBJDIR)/ksh/linear_lu.d
$(OBJDIR)/ksh/linear_lu.o: ksh/linear_lu.cpp | $(OBJDIR)/ksh
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ -I ksh $<

all: $(OUTDIR)/libksh.a
$(OUTDIR)/libksh.a: $(objectfiles) | $(OUTDIR)
	ar crs $@ $^

#
# Binaries
#

all: rktest.exe
-include $(OBJDIR)/rktest.d
-include $(OBJDIR)/rktest_erk.d
$(OBJDIR)/rktest.o: rktest.cpp ksh/embedded_runge_kutta.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ $<
$(OBJDIR)/rktest_erk.o: rktest_erk.cpp ksh/erk.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ $<
rktest.exe: $(OBJDIR)/rktest.o $(OBJDIR)/rktest_erk.o $(OBJDIR)/ksh/embedded_runge_kutta.o
	$(CXX) $(LDFLAGS) -o $@ $^

$(OUTDIR)/rk/rkeuler.txt: | $(OUTDIR)/rk
	./rktest.exe

directories += $(OUTDIR)/rk
all: $(OUTDIR)/rk/rktest.pdf
$(OUTDIR)/rk/rktest.pdf: rktest.gp $(OUTDIR)/rk/rkeuler.txt | $(OUTDIR)/rk
	gnuplot rktest.gp


#------------------------------------------------------------------------------
#
# Check
#

check_LDFLAGS := $(LDFLAGS) -L $(OUTDIR)
check_LIBS    := $(LIBS) -lksh
check:
.PHONY: check

directories += $(OBJDIR)/test

check: linear
.PHONY: linear
linear: test/linear.exe
	./$<
-include $(OBJDIR)/test/linear.d
$(OBJDIR)/test/linear.o: test/linear.cpp | $(OBJDIR)/test
	$(CXX) $(CXXFLAGS) -I . -MD -MF $(@:.o=.d) -c -o $@ $<
test/linear.exe: $(OBJDIR)/test/linear.o $(OUTDIR)/libksh.a
	$(CXX) $(check_LDFLAGS) -o $@ $^ $(check_LIBS)

clean:
	-rm -rf *.o *.d

$(directories):
	mkdir -p $@
