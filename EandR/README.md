This archive contains code to perform the likelihood calculations
described in Terhorst, Schlötterer & Song (2015).

### Compilation Instructions
There are two dependencies required to build this package which are not
standard on all UNIX-like systems:

  * A C++11-compliant compiler (such as GCC 4.8 or higher). The clang
    compiler shipped with recent versions of Mac OS X should also work.

  * The google-sparsehash library (https://code.google.com/p/sparsehash/).

To compile the package, type `make` in the directory containing this
README file.

### Usage Instructions
The file `example.py` is example of how to use this library to perform
inference. It analyzes a test data set (also included) in order to
test for and estimate selection. Refer to the comments in that file
for further details. Execute `example.py` using the wrapper script
`example.sh` in order to ensure that the threelocus library can be
located on your system's library path.
