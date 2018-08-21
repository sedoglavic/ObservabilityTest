# ObservabilityTest

A maple package that test observability/identifiability of ordinary differential systems in polynomial time.

## Getting Started

### Prerequisites

This package works with maple computer system (with Maple 2017 and previous versions). 

### Installing 

Clone the repository.

This creates: 
* a directory `ObservabilityTest` (source code) and in it:
* a directory `examples` 
* a directory `release` with three files: `maple.lib`, `maple.ind`, `maple.hdb`.

### Deployment
To make the observabilityTest procedure available to your Maple system, the absolute path to the newly created `ObservabilityTest/release` directory has to be put in the Maple global variable `libname`.
This is either done once forever in your `.mapleinit` file or has to be done at the beginning of each Maple session by the command
```
 libname:="<ObservabilityTest release path>",libname:
```
where `<ObservabilityTest release path>` depends on your system type.

#### Under Unix (like):
The `ObservabilityTest/release` path is something like 
```
      /foo/bar/.../ObservabilityTest/release
```
 so that your command reads:
```
libname:="/foo/bar/.../ObservabilityTest/release",libname:
```
#### Under Windows:
The `ObservabilityTest/release` path is something like 
```
C:\ObservabilityTest\release
``` 
so that your command reads:
```
libname:="C:\ObservabilityTest\release",libname:
```

## Running the tests

Some examples are available using a command like
```
read "/foo/bar/../ObservabilityTest/examples/CGV1990.mpl";
```
in maple or
```
maple -i examples/CGV1990.mpl
```
in your bash (see the files in directory examples).

## Remarks

The given Makefile is adapted to Unix-like environment

## Warnings

This package comes with ! NO WARRANTY !

## Reference (to cite this package):

Sedoglavic, A. (2002), A probabilistic algorithm to test local algebraic observability in polynomial time. Journal of Symbolic Computation, 33(5), 735-755.

## Author
Any support request, comments, bug reports, critics are welcome!
Please send to: email: Alexandre.Sedoglav@univ-lille.fr

## Licence
GNU GPL v2
