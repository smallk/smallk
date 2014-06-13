
# locations of source files separated by ':'
VPATH=src

SRC = \
	smallk_test.cpp

TARGET=smallk_tester

# paths to all include folders
INCLUDES = \
	-Iinclude

# object files replace the .cpp extension with .o
OBJS = $(patsubst %.cpp, objs/%.o, $(SRC))

all: smallk_bin

smallk_bin: bin/$(TARGET)

# link the test program with libsmallk

bin/$(TARGET) : $(OBJS) ../build/bin/$(LIBNAME).a
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ -L../build/bin $(LDFLAGS) -lsmallk $(LIBS)


# The -MMD switch causes several things to happen.  First, a makefile rule is
# generated that contains the the compiler to generate dependency files for 
# each file being compiled.  The dependency files are given a '.d' extension 
# and are written to the default location for preprocessed output.  The 
# 'include' statement below includes them in the makefile.  The -MP option 
# causes a phony target to be generated for each dependency other than the 
# main CPP file being compiled (such as for header files).  This feature 
# prevents the problem of having to include headers explicitly in the makefile
# and keeping them updated.
objs/%.o: %.cpp
	@echo "PWD = $(shell pwd)"
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -MMD -MP -o $@ $<

check: smallk_bin
	sh ../tests/scripts/test_smallk.sh

clean:
	@rm -f *.csv
	@rm -f *.xml
	@rm -rf bin* objs* 

# include all dependency files as individual targets
ifneq ($(MAKECMDGOALS),clean)
-include $(wildcard objs/*.d)
endif

.PHONY: all clean check
