# -*- mode:makefile-gmake -*-

all:
.PHONY: all

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
  CXXFLAGS:= -Wall -std=gnu++17 -O3 -march=native
  LDFLAGS:= -Wall -std=gnu++17 -O3 -march=native

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
  CXXFLAGS:= -Wall -std=gnu++17 -O3 -march=native
  LDFLAGS:= -Wall -std=gnu++17 -O3 -march=native

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
  CXXFLAGS := -Wall -std=gnu++17 -O3 -march=native
  LDFLAGS := -Wall -std=gnu++17 -O3 -march=native

  # libmwg
  LIBMWG_PREFIX := $(HOME)/opt/libmwg-20170824
  LIBMWG_BUILD  := x86_64-pc-linux-gnu-gcc-7.3.1+cxx98-debug
  CXXFLAGS += \
    -isystem $(LIBMWG_PREFIX)/include \
    -isystem $(LIBMWG_PREFIX)/include/$(LIBMWG_BUILD)
  LDFLAGS += \
    -L $(LIBMWG_PREFIX)/lib/$(LIBMWG_BUILD)
endif

ifeq ($(HOSTNAME),letsnote2025)
  host:=found

  OTFFT_PREFIX   := $(HOME)/opt/otfft-6.4
  CXX := g++
  CXXFLAGS := -Wall -std=gnu++17 -O3 -march=native
  LDFLAGS := -Wall -std=gnu++17 -O3 -march=native

  # libmwg
  LIBMWG_PREFIX := $(HOME)/.opt/libmwg/master
  LIBMWG_BUILD  := x86_64-cygwin-gcc-13.4.0+cxx98-debug
  CXXFLAGS += \
    -isystem $(LIBMWG_PREFIX)/include \
    -isystem $(LIBMWG_PREFIX)/include/$(LIBMWG_BUILD)
  LDFLAGS += \
    -L $(LIBMWG_PREFIX)/lib/$(LIBMWG_BUILD)
endif

ifeq ($(host),unknown)
  CXX:=g++
  CXXFLAGS:= -Wall -std=gnu++17 -O3 -march=native
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

define register_object
$(register_target)_objects += $$(OBJDIR)/$1.o
-include $$(OBJDIR)/$1.d
$$(OBJDIR)/$1.o: $1.cpp | $$(OBJDIR)/ksh
	$$(CXX) $$(CXXFLAGS) -MP -MD -MF $$(@:.o=.d) -c -o $$@ $$<
endef

define register_binary
$(register_target).exe: $$($(register_target)_objects)
	$$(CXX) $$(LDFLAGS) -o $$@ $$^
endef

all: $(OUTDIR)/libksh.a
directories += $(OBJDIR)/ksh

register_target = libksh
$(eval $(call register_object,ksh/embedded_runge_kutta))
$(eval $(call register_object,ksh/integrator))
$(eval $(call register_object,ksh/linear_lu))
$(OUTDIR)/libksh.a: $(libksh_objects) | $(OUTDIR)
	ar crs $@ $^

#
# rktest
#

all: rktest.exe
register_target = rktest
$(eval $(call register_object,rktest))
$(eval $(call register_object,rktest_erk))
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
	$(CXX) $(CXXFLAGS) -I . -MP -MD -MF $(@:.o=.d) -c -o $@ $<
test/linear.exe: $(OBJDIR)/test/linear.o $(OUTDIR)/libksh.a
	$(CXX) $(check_LDFLAGS) -o $@ $^ $(check_LIBS)

define register_test
check: $1
.PHONY: $1
$1: test/test_$1.exe
	./$$<
-include $$(OBJDIR)/test/test_$1.d
$$(OBJDIR)/test/test_$1.o: test/test_$1.cpp | $$(OBJDIR)/test
	$$(CXX) $$(CXXFLAGS) -I . -MP -MD -MF $$(@:.o=.d) -c -o $$@ $$<
test/test_$1.exe: $$(OBJDIR)/test/test_$1.o $$(OUTDIR)/libksh.a
	$$(CXX) $$(check_LDFLAGS) -o $$@ $$^ $$(check_LIBS)
endef
$(eval $(call register_test,integrator))
$(eval $(call register_test,polynomial))
$(eval $(call register_test,rational))
$(eval $(call register_test,big_integer))

experiment: multi_precision
.PHONY: experiment

register_target = test_multi_precision
$(eval $(call register_object,test/multi_precision))
$(eval $(call register_binary,test/multi_precision))

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
