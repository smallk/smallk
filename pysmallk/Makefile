
PY = python

TARGET = pysmallk.so

SRC = setup.py

PFLAGS = build_ext --inplace

TEST_DIR = tests/

all: pysmallk

pysmallk: 
	$(PY) $(SRC) $(PFLAGS)

install: 
	cp $(TARGET) $(INSTALLATION_DIR)

check:
	sh ../tests/scripts/test_pysmallk.sh

clean:
	@rm -r -f *.so
	@rm -f *.csv
	@rm -f *.xml
	@rm -rf bin* objs* build*
	@rm -rf interface/*.cpp


uninstall:
	rm -f $(INSTALLATION_DIR)$(TARGET) >/dev/null 2>&1

.PHONY: all clean check uninstall
