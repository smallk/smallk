Notes for the RC3 release:

      0. When running 'make check', the smallk_data repository must be cloned first.
            See notes at:
            https://github.com/smallk/smallk_data
            
            Use 'make check DATA_DIR=../smallk_data
            
      1. The cython interface (pysmallk) is currently broken.
      2. Binaries of all compiled C++ smallk programs can be called from
            Python with an 'exec' call.

      These will be updated in a future release.

