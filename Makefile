###############################################################################
#
# This is the master makefile for the smallk project.  Much of this is based
# on information in the books 'Autotools: A Practitioner's Guide to GNU
# Autoconf, Automake, and Libtool', by John Calcote, and 'Managing Projects
# with GNU Make', by Robert Mecklinburg.
#
# Eventually we will have a build system based either on the autotools or
# CMake.
#
###############################################################################

###############################################################################
#
# The default install location for smallk is /usr/local/smallk, unless 
# SMALLK_INSTALL_DIR exists, in which case smallk will be installed there.
#
###############################################################################

SMALLK_INSTALL_DIR ?= /usr/local/smallk
export SMALLK_INSTALL_DIR

###############################################################################
#
# Elemental
#
# Sets the version of Elemental to link against, as well as the appropriate
# compiler and linker flags from the Elemental build.
#
# If multiple versions of Elemental have been installed, select the version
# to link against with the make command-line argument 'ELEMVER'.
#
# If ELEMENTAL_INSTALL_DIR exists, use that to locate the appropriate 
# Elemental config files.  If it does not exist, assume that Elemental has 
# been installed to /usr/local/elemental.
#
###############################################################################

# default to Elemental version 0.84
ELEMVER ?= 0.84

# name of the elemental config file
ELEMVARIABLES_FILE := ElemVars

# override if command-line arg 'ELEMVER' has been specified
ifeq ($(ELEMVER), 0.81)
ELEMENTAL_VERSION := 0.81
ELEM_INT_VER := 81
ELEMVARIABLES_FILE := elemvariables
endif

ifeq ($(ELEMVER), 0.82)
$(error Elemental release 0.82 is not supported)
endif

ifeq ($(ELEMVER), 0.82-p1)
$(error Elemental release 0.82-p1 is not supported)
endif

ifeq ($(ELEMVER), 0.83)
ELEMENTAL_VERSION := 0.83
ELEM_INT_VER := 83
ELEMVARIABLES_FILE := elemvariables
endif

ifeq ($(ELEMVER), 0.84)
ELEMENTAL_VERSION := 0.84
ELEM_INT_VER := 84
endif

# find the top-level Elemental install directory
ELEMENTAL_INSTALL_DIR ?= /usr/local/elemental
ELEMENTAL_BASE = $(ELEMENTAL_INSTALL_DIR)/$(ELEMENTAL_VERSION)

# locate the appropriate Elemental config file based on the build type
ifeq ($(CFG), purerelease)
ELEMVARS := $(ELEMENTAL_BASE)/PureRelease/conf/$(ELEMVARIABLES_FILE)
else
ELEMVARS := $(ELEMENTAL_BASE)/HybridRelease/conf/$(ELEMVARIABLES_FILE)
endif

# find the file and include it; error if not found
ifneq (,$(wildcard $(ELEMVARS)))
include $(ELEMVARS)
else
$(error Cannot find the elemvariables file)
endif
export ELEMVARS

###############################################################################
#
# Other
#
###############################################################################
package     = libsmallk
version     = 1.1.0
tarname     = $(package)
distdir     = $(tarname)-$(version)
prefix      = $(SMALLK_INSTALL_DIR)
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
includedir  = $(prefix)/include
configdir   = $(prefix)/conf
libdir      = $(exec_prefix)/lib
export prefix
export exec_prefix
export includedir
export configdir
export libdir
export bindir

# uncomment to use AddressSanitizer (must be in debug mode)
#ASAN_CXXFLAGS = -fsanitize=address -fno-omit-frame-pointer
#ASAN_LDFLAGS = -fsanitize=address

# detect the operating system
# returns either 'Darwin' (MacOSX) or 'Linux'
UNAME := $(shell uname)

# Linux needs to explicitly link with librt
SYSTEM_FILE := system_mac.cpp
ifeq ($(UNAME), Linux)
LIBRT := -lrt
SYSTEM_FILE := system_posix.cpp
endif

# macros passed directly to the compiler on the command line
DEFINES := -D ELEM_VER=$(ELEM_INT_VER)

# compiler flags, linker flags, and libs for all configurations
CXXFLAGS += $(DEFINES) -Wall $(ELEM_COMPILE_FLAGS)
LDFLAGS += -v $(ELEM_LINK_FLAGS)
LIBS = $(ELEM_LIBS)

# add support for OpenMP if not a purerelease build
ifneq ($(CFG), purerelease)
CFG=release
CXXFLAGS += -fopenmp
LIBS += -lgomp
endif

# name of the smallk config file, similar to Elemental's elemvariables file
SMALLK_CONFIGFILE=smallkvariables

# locations of source files separated by ':'
#VPATH=

SMALLK_SRC = \
	smallk/src/smallk.cpp \
	common/src/utils.cpp  \
	common/src/matrix_market_file.cpp \
	common/src/delimited_file.cpp \
	common/src/constants.cpp \
	common/src/nmf.cpp \
	common/src/nmf_options.cpp \
	common/src/bit_matrix.cpp \
	common/src/bit_matrix_ops.cpp \
	common/src/spooky_v2.cpp \
	common/src/nnls.cpp \
	common/src/$(SYSTEM_FILE) \
	common/src/file_loader.cpp \
	common/src/flat_clust_output.cpp \
	common/src/flatclust_json_writer.cpp \
	common/src/flatclust_xml_writer.cpp \
	common/src/assignments.cpp \
	hierclust/src/clust.cpp \
	hierclust/src/tree.cpp \
	hierclust/src/node.cpp \
	hierclust/src/clust_options.cpp \
	hierclust/src/hierclust_writer_factory.cpp \
	hierclust/src/hierclust_json_writer.cpp \
	hierclust/src/hierclust_xml_writer.cpp

# object files replace the .cpp extension with .o
LIB_OBJS  = $(patsubst %.cpp, build/objs/%.o, $(SMALLK_SRC))

# paths to all include folders
INCLUDES = \
	$(MPI_CXX_INCLUDE_STRING) \
	-Ismallk/include \
	-Icommon/include \
	-Ihierclust/include \
	-I$(ELEM_INC)

# libraries - add install location of libsmallk.a to this - TBD
LIBS += $(LIBRT) -lm

# includes for the smallkvariables file
SMALLK_INC = \
	$(MPI_CXX_INCLUDE_STRING) \
	-I$(includedir) \
	-I$(ELEM_INC)

export CXX
export CFG
export LIBS
export LDFLAGS
export INCLUDES
export CXXFLAGS
export SYSTEM_FILE

LIB_TARGET=libsmallk
LIBNAME=$(LIB_TARGET)
export LIBNAME

all: inform libsmallk smallk_tester matrixgen preprocessor nmf hierclust flatclust tests

smallk_tester: build/bin/$(LIBNAME).a
	cd smallk && $(MAKE) smallk_bin

matrixgen:
	cd matrixgen && $(MAKE) all

preprocessor: inform
	cd preprocessor && $(MAKE) all

nmf: inform
	cd nmf && $(MAKE) all

hierclust: inform
	cd hierclust && $(MAKE) all

flatclust: inform
	cd flatclust && $(MAKE) all

tests: inform
	cd tests && $(MAKE) all

libsmallk: build/bin/$(LIBNAME).a

build/bin/$(LIBNAME).a : $(LIB_OBJS)
	@mkdir -p $(dir $@)
	ar rs build/bin/$(LIBNAME).a $(LIB_OBJS)
	@echo "ELEMVARS = $(ELEMVARS)" > build/$(SMALLK_CONFIGFILE)
	@echo "CFG = $(CFG)" >> build/$(SMALLK_CONFIGFILE)
	@echo "CC = $(CC)" >> build/$(SMALLK_CONFIGFILE)
	@echo "CXX = $(CXX)" >> build/$(SMALLK_CONFIGFILE)

install: all
	@echo "SMALLK_INC = $(SMALLK_INC)" >> build/$(SMALLK_CONFIGFILE)
	@echo "SMALLK_CXXFLAGS = $(CXXFLAGS)" >> build/$(SMALLK_CONFIGFILE)
	@echo "SMALLK_LDFLAGS = -L$(libdir) $(LDFLAGS)" >> build/$(SMALLK_CONFIGFILE)
	@echo "SMALLK_LIBS = -lsmallk $(LIBS)" >> build/$(SMALLK_CONFIGFILE)
	install -d $(DESTDIR)$(includedir)
	install -d $(DESTDIR)$(configdir)
	install -d $(DESTDIR)$(libdir)
	install -d $(DESTDIR)$(bindir)
	install smallk/include/smallk.hpp $(DESTDIR)$(includedir)
	install build/$(SMALLK_CONFIGFILE) $(DESTDIR)$(configdir)
	install build/bin/$(LIBNAME).a $(DESTDIR)$(libdir)
	cd matrixgen && $(MAKE) install
	cd preprocessor && $(MAKE) install
	cd nmf && $(MAKE) install
	cd hierclust && $(MAKE) install
	cd flatclust && $(MAKE) install

uninstall:
	@rm -f $(DESTDIR)$(includedir)/smallk.hpp
	@rm -f $(DESTDIR)$(configdir)/$(SMALLK_CONFIGFILE)
	@rm -f $(DESTDIR)$(libdir)/$(LIBNAME).a
	cd matrixgen && $(MAKE) uninstall
	cd preprocessor && $(MAKE) uninstall
	cd nmf && $(MAKE) uninstall
	cd hierclust && $(MAKE) uninstall
	cd flatclust && $(MAKE) uninstall

# The -MMD switch causes several things to happen.  First, a makefile rule is
# generated that contains the the compiler to generate dependency files for 
# each file being compiled.  The dependency files are given a '.d' extension 
# and are written to the default location for preprocessed output.  The 
# 'include' statement below includes them in the makefile.  The -MP option 
# causes a phony target to be generated for each dependency other than the 
# main CPP file being compiled (such as for header files).  This feature 
# prevents the problem of having to include headers explicitly in the makefile
# and keeping them updated.
build/objs/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -MMD -MP -o $@ $<

clean:
	@rm -rf build
	cd smallk && $(MAKE) clean
	cd matrixgen && $(MAKE) clean
	cd preprocessor && $(MAKE) clean
	cd nmf && $(MAKE) clean
	cd hierclust && $(MAKE) clean
	cd flatclust && $(MAKE) clean
	cd tests && $(MAKE) clean

distclean: clean
	rm -rf $(distdir)
	rm -f $(distdir).tar.gz

dist: $(distdir).tar.gz

$(distdir).tar.gz: $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $@
	rm -rf $(distdir)

# generate the folders for the distribution
projects := common flatclust hierclust matrixgen nmf preprocessor smallk tests
proj_inc := $(addsuffix /include,$(projects))
proj_src := $(addsuffix /src,$(projects))

$(distdir): FORCE
	mkdir -p $(addprefix $(distdir)/,$(proj_inc))
	mkdir -p $(addprefix $(distdir)/,$(proj_src))
	mkdir -p $(distdir)/examples
	mkdir -p $(distdir)/data/test
	mkdir -p $(distdir)/tests/scripts
	mkdir -p $(distdir)/doc
	for d in $(proj_inc); \
	do \
	cp -r $$d/ $(distdir)/$$d; \
	done
	for d in $(proj_src); \
	do \
	cp -r $$d/ $(distdir)/$$d; \
	done
	cp Makefile $(distdir)
	cp README.txt $(distdir)
	cp LICENSE-2.0.txt $(distdir)
	cp smallk/Makefile $(distdir)/smallk
	cp examples/Makefile $(distdir)/examples
	cp examples/smallk_example.cpp $(distdir)/examples
	cp preprocessor/Makefile $(distdir)/preprocessor
	cp nmf/Makefile $(distdir)/nmf
	cp hierclust/Makefile $(distdir)/hierclust
	cp flatclust/Makefile $(distdir)/flatclust
	cp matrixgen/Makefile $(distdir)/matrixgen
	cp tests/Makefile $(distdir)/tests
	cp data/matrix.mtx $(distdir)/data
	cp data/dictionary.txt $(distdir)/data
	cp data/documents.txt $(distdir)/data
	cp data/test/reduced_matrix_20news.mtx $(distdir)/data/test
	cp data/test/reduced_dictionary_20news.txt $(distdir)/data/test
	cp data/test/reduced_documents_20news.txt $(distdir)/data/test
	cp data/reuters.mtx $(distdir)/data
	cp data/reuters_dictionary.txt $(distdir)/data
	cp data/nmf_init_w.csv $(distdir)/data
	cp data/nmf_init_h.csv $(distdir)/data
	cp data/nmf_rank2_init_w.csv $(distdir)/data
	cp data/nmf_rank2_init_h.csv $(distdir)/data
	cp data/test/nmf_result_w.csv $(distdir)/data/test
	cp data/test/nmf_result_h.csv $(distdir)/data/test
	cp data/hierclust_init_w.csv $(distdir)/data
	cp data/hierclust_init_h.csv $(distdir)/data
	cp data/test/reuters_tree_5.xml $(distdir)/data/test
	cp data/test/reuters_assignments_5.csv $(distdir)/data/test
	cp data/rnd_256_256.csv $(distdir)/data
	cp data/flatclust_init_w.csv $(distdir)/data
	cp data/flatclust_init_h.csv $(distdir)/data
	cp data/test/flatclust_rnd_clusters_16.xml $(distdir)/data/test
	cp data/test/flatclust_rnd_assignments_16.csv $(distdir)/data/test
	cp tests/scripts/* $(distdir)/tests/scripts
	cp doc/smallk_readme.pdf $(distdir)/doc

distcheck: $(distdir).tar.gz
	gzip -cd $(distdir).tar.gz | tar xvf -
	cd $(distdir) && $(MAKE) all
	cd $(distdir) && $(MAKE) check
	cd $(distdir) && $(MAKE) DESTDIR=$${PWD}/_inst install
	cd $(distdir) && $(MAKE) DESTDIR=$${PWD}/_inst uninstall
	@remaining="`find $${PWD}/$(distdir)/_inst -type f | wc -l`"; \
	if test "$${remaining}" -ne 0; then \
	  echo "distcheck: uninstall left $${remaining} file(s) in the temp directory!"; \
	  exit 1; \
	fi
	cd $(distdir) && $(MAKE) clean
	rm -rf $(distdir)
	@echo "The tarball $(distdir).tar.gz is ready."

check:
	sh tests/scripts/test_smallk.sh | tee smallk_test_results.txt
	@count="`grep -i 'failed' smallk_test_results.txt | wc -l`"; \
	if test "$${count}" -ne 0; then \
	echo "***** $${count} tests failed, exiting... *****"; \
	exit 1; \
	else \
	echo "***** All tests passed. *****"; \
	fi
	@rm smallk_test_results.txt

FORCE:
	-rm $(distdir).tar.gz > /dev/null 2>&1
	-rm -rf $(distdir) > /dev/null 2>&1

# include all dependency files as individual targets
ifneq ($(MAKECMDGOALS),clean)
-include $(wildcard build/objs/*.d)
endif

# tell the user the valid config options
inform:
ifneq ($(CFG),debug)
ifneq ($(CFG),purerelease)
ifneq ($(CFG),release)
	@echo "Invalid configuration: "$(CFG)
	@echo "Allowed values of CFG are 'debug', 'purerelease', and 'release'"
	@exit 1
endif
endif
endif
	@echo "Build configuration: "$(CFG)

.PHONY: FORCE inform all clean distclean check dist distcheck install uninstall \
	smallk_tester matrixgen preprocessor nmf hierclust flatclust tests
