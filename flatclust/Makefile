include $(ELEMVARS)

# locations of source files separated by ':'
VPATH=src:../common/src

SRC = \
	main.cpp                  \
	command_line.cpp          \
	matrix_market_file.cpp    \
	utils.cpp                 \
	file_loader.cpp           \
	flat_clust_options.cpp    \
	nmf.cpp                   \
	nmf_options.cpp           \
	$(SYSTEM_FILE)            \
	delimited_file.cpp        \
	assignments.cpp           \
	flat_clust.cpp            \
	bit_matrix.cpp            \
	bit_matrix_ops.cpp        \
	spooky_v2.cpp             \
	nnls.cpp                  \
	constants.cpp             \
	flat_clust_output.cpp     \
	flatclust_xml_writer.cpp  \
	flatclust_json_writer.cpp

# object files replace the .cpp extension with .o
OBJS = $(patsubst %.cpp, objs/%.o, $(SRC))

TARGET=flatclust

# # paths to all include folders
INCLUDES = \
	-Iinclude \
	-I../common/include \
	-I$(ELEMENTAL_INCLUDES)

all: flatclust_bin

flatclust_bin: bin/$(TARGET)

# link the target binary
bin/$(TARGET) : $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

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
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -MMD -MP -o $@ $<

clean:
	@rm -rf bin* objs*

install:
	install -m 0755 bin/$(TARGET) $(DESTDIR)$(bindir)

uninstall:
	-rm $(DESTDIR)$(bindir)/$(TARGET) >/dev/null 2>&1

# include all dependency files as individual targets
ifneq ($(MAKECMDGOALS),clean)
-include $(wildcard objs/*.d)
endif

.PHONY: clean all flatclust_bin install uninstall

