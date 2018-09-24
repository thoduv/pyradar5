# pyradar5

The radar5 module is a numerical solver for delay differential equations. It is a wrapper around an original FORTRAN code, that uses a colocation method on Radau nodes for an efficient integration of stiff problems. Please note that all the capabilities of the FORTRAN code are not wrapped yet. Here is what is currently available:

* Integration of DDE specified by Python function or runtime-compiled C code.
* Constant or time-dependant initial conditions, with interpolation if needed.
* Constant, variable-dependant, or time-dependant delay.

Here is what is not implemented:

* User-specified Jacobian, and delay-components Jacobian.
* Implicit systems and mass matrix
* Advanced breakpoints detection

The original FORTRAN code is the work of Nicola Guglielmi and Ernst Hairer, and can be found on: http://www.unige.ch/~hairer/software.html

## Installation

The module has been successfully tested on Windows and Linux x86_64. It has not been tested on OSX yet. 

The simplest installation process should be to use the command:

```
pip install radar5
```

On _Windows_, this will install binary packages. Please note that they are build against a _recent_ version of numpy. Updating your version of numpy can be required if you run into the error `RuntimeError: module compiled against API version 0xc but this version of numpy is 0x9`. If you have installed Python using Αnaconda, you can do this by running the `conda upgrade numpy` command. If you have only used pip, run `pip install --upgrade numpy`. Finally, if you don't know much about how your Python has been installed, a foolproof solution can be to copy the binary files into your working directory. Such files can be found here: https://github.com/thoduv/pyradar5/tree/master/windows_binaries

On _Linux systems_, no binary package are provided, but building is easy as long as you have a FORTRAN compiler. Grab any version of `gfortran` on your package manager. Be careful if you have installed numpy by your own means (probably through the system package manager). In that case I strongly suggest that you do not install a concurrent version through pip. This can be avoided using the following command.

```
pip install --no-deps radar5
```

This will build the package against the systemwide version of numpy version, which is great.

## Usage

A very simple example can be run straight away to check if the installation went well.
```python
import radar5
radar5.test()
```

If you have `matplotlib`, a window with two oscillating curves should pop up.

Then, more useful examples can be found in the example folder of this repository.

## Changelog

0.1 - 21 September 2018
-- Initial release.


