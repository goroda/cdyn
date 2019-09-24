# Dynamical Systems in C

A very simple library in C that has some utilities for monitoring and time stepping for ordinary differential equations. It is mostly used as a standalone library in my other packages that work with dynamical systems

## Prerequisites
  * BLAS
  * LAPACK
  
## Installation instructions

```
cd <parent-directory-for-library>  
git clone https://github.com/goroda/cdyn.git  
cd cdyn  
mkdir build && cd build  
cmake -DCMAKE_INSTALL_PREFIX:PATH=<location> ..  
make  
make install  
```

Notes:  
1)  If no install prefix is specified then some default installation will be used, typically  <kbd> </usr/local> </kbd> in unix systems  
2)  you make need to use do <kbd>sudo make install</kbd> for proper permissions


## Known problems

### error while loading shared libraries: libcdyn.so: cannot open shared object file: No such file or directory

This error occurs because the system did not properly configure and link the shared library when it was created. It can be resolved (found solution [here](https://stackoverflow.com/questions/480764/linux-error-while-loading-shared-libraries-cannot-open-shared-object-file-no-s)) by running the following command

```
sudo ldconfig -v
```


## More information

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2016 Massachusetts Institute of Technology  
Copyright (c) 2018 University of Michigan  
License: MIT
