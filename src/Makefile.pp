# -*- mode:makefile-gmake -*-

all:
.PHONY: all

Makefile: Makefile.pp
	mwg_pp.awk $< > $@

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

  # libmwg
  LIBMWG_PREFIX := $(HOME)/opt/libmwg-201509
  LIBMWG_BUILD  := i686-pc-linux-gnu-gcc-6.3.1+cxx11-release
  CXXFLAGS += \
    -isystem $(LIBMWG_PREFIX)/include \
    -isystem $(LIBMWG_PREFIX)/include/$(LIBMWG_BUILD)
  LDFLAGS += \
    -L $(LIBMWG_PREFIX)/lib/$(LIBMWG_BUILD)
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

ifeq ($(HOSTNAME),chatoyancy)
  host:=found

  OTFFT_PREFIX   := $(HOME)/opt/otfft-6.4
  CXX := g++
  CXXFLAGS := -Wall -std=gnu++14 -O3 -march=native
  LDFLAGS := -Wall -std=gnu++14 -O3 -march=native

  # libmwg
  LIBMWG_PREFIX := $(HOME)/opt/libmwg-20170824
  LIBMWG_BUILD  := x86_64-pc-linux-gnu-gcc-7.3.1+cxx98-debug
  CXXFLAGS += \
    -isystem $(LIBMWG_PREFIX)/include \
    -isystem $(LIBMWG_PREFIX)/include/$(LIBMWG_BUILD)
  LDFLAGS += \
    -L $(LIBMWG_PREFIX)/lib/$(LIBMWG_BUILD)
endif

ifeq ($(host),unknown)
  CXX:=g++
  CXXFLAGS:= -Wall -std=gnu++14 -O3 -march=native
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

#%m register_object
#%%x
$"target"_objects += $(OBJDIR)/%name%.o
-include $(OBJDIR)/%name%.d
$(OBJDIR)/%name%.o: %name%.cpp | $(OBJDIR)/ksh
	$(CXX) $(CXXFLAGS) -MD -MF $(@:.o=.d) -c -o $@ $<
#%%end.i
#%end
#%m register_binary
##%x
%target%.exe: $($"target"_objects)
	$(CXX) $(LDFLAGS) -o $@ $^
##%end.i
#%end


all: $(OUTDIR)/libksh.a
directories += $(OBJDIR)/ksh
#%[target="libksh"]
#%x register_object.r|%name%|ksh/embedded_runge_kutta|
#%x register_object.r|%name%|ksh/integrator|
#%x register_object.r|%name%|ksh/linear_lu|
$(OUTDIR)/libksh.a: $(libksh_objects) | $(OUTDIR)
	ar crs $@ $^

#
# rktest
#

all: rktest.exe
#%[target="rktest"]
#%x register_object.r|%name%|rktest|
#%x register_object.r|%name%|rktest_erk|
rktest.exe: $(rktest_objects) $(OBJDIR)/ksh/embedded_runge_kutta.o
	$(CXX) $(LDFLAGS) -o $@ $^

directories += $(OUTDIR)/rk
$(OUTDIR)/rk/rkeuler.txt: rktest.exe | $(OUTDIR)/rk
	./rktest.exe
all: $(OUTDIR)/rk/rktest.pdf
$(OUTDIR)/rk/rktest.pdf: rktest.gp $(OUTDIR)/rk/rkeuler.txt | $(OUTDIR)/rk
	gnuplot rktest.gp

.PHONY: rktest
rktest: rktest.gp rktest.exe | $(OUTDIR)/rk
	./rktest.exe
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

#%m register_test
check: %name%
.PHONY: %name%
%name%: test/test_%name%.exe
	./$<
-include $(OBJDIR)/test/test_%name%.d
$(OBJDIR)/test/test_%name%.o: test/test_%name%.cpp | $(OBJDIR)/test
	$(CXX) $(CXXFLAGS) -I . -MD -MF $(@:.o=.d) -c -o $@ $<
test/test_%name%.exe: $(OBJDIR)/test/test_%name%.o $(OUTDIR)/libksh.a
	$(CXX) $(check_LDFLAGS) -o $@ $^ $(check_LIBS)
#%end
#%x register_test.r/%name%/integrator/
#%x register_test.r/%name%/polynomial/
#%x register_test.r/%name%/rational/
#%x register_test.r/%name%/big_integer/

experiment: multi_precision
.PHONY: experiment

#%[target="test_multi_precision"]
#%x register_object.r|%name%|test/multi_precision|
#%x register_binary.r|%out%|test/multi_precision|

#------------------------------------------------------------------------------
#
# Samples
#

sample_LDFLAGS := $(LDFLAGS) -L $(OUTDIR)
sample_LIBS    := $(LIBS) -lksh
sample:
.PHONY: sample
directories += $(OBJDIR)/sample

%.sample: $(OBJDIR)/sample/%.exe
	@echo SAMPLE $@
	@./$<
$(OBJDIR)/sample/%.exe: $(OBJDIR)/sample/%.o $(OUTDIR)/libksh.a
	$(CXX) $(sample_LDFLAGS) -o $@ $^ $(sample_LIBS)
$(OBJDIR)/sample/%.o: sample/%.cpp | $(OBJDIR)/sample
	$(CXX) $(CXXFLAGS) -I . -MD -MF $(@:.o=.d) -c -o $@ $<

sample-names += legendre_polynomial
sample-names += ode_dop853

-include $(sample-names:%=$(OBJDIR)/sample/%.d)
sample: $(sample-names:%=%.sample)

.SECONDARY:

#------------------------------------------------------------------------------

clean:
	-find $(OBJDIR) -name \*.d -o -name \*.o | xargs rm -f

$(directories):
	mkdir -p $@
