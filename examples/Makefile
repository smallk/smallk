###############################################################################
#
# This is a sample makefile for use with the smallk library.
#
# IMPORTANT: you must first build and install the smallk library!
#
# If you do not install smallk into the default location (~/smallk), then 
# either create an environment variable called SMALLK_INSTALL_DIR that 
# specifies the path to the top-level folder of your smallk installation, or
# define the same variable on the make command line.
#
###############################################################################

# find the smallk config file and include it; error if not found
SMALLK_INSTALL_DIR ?= /usr/local/smallk
SMALLK_VAR_FILE := $(SMALLK_INSTALL_DIR)/conf/smallkvariables
ifneq (,$(wildcard $(SMALLK_VAR_FILE)))
include $(SMALLK_VAR_FILE)
else
$(error Cannot find the smallkvariables file)
endif

# locations of source files separated by ':'
VPATH=src

SRC = \
	smallk_example.cpp

TARGET=example

# paths to all include folders
INCLUDES = \
	$(SMALLK_INC)

# object files replace the .cpp extension with .o
OBJS = $(patsubst %.cpp, objs/%.o, $(SRC))

all: example_bin

example_bin: bin/$(TARGET)

# link the test program with libsmallk
bin/$(TARGET) : $(OBJS)
	@echo "SMALLK_LIBS: $(SMALLK_LIBS)"
	@mkdir -p $(dir $@)
	$(CXX) $(SMALLK_CXXFLAGS) $(INCLUDES) -o $@ $^ $(SMALLK_LDFLAGS) $(SMALLK_LIBS)

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
	@mkdir -p $(dir $@)
	$(CXX) -c $(SMALLK_CXXFLAGS) $(INCLUDES) -MMD -MP -o $@ $<

clean:
	@rm -rf bin* objs* 

# include all dependency files as individual targets
ifneq ($(MAKECMDGOALS),clean)
-include $(wildcard objs/*.d)
endif

.PHONY: all clean example_bin
