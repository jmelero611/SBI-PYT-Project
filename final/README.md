# Reconstruction of a macro-complex using interacting subunits

This project aims to construct a whole protein complex using pair of interacting chains.

## Getting Started

These instructions will get you how to install the package and use it.

### Prerequisites

There are package required to run the program. These packages are Biopython (Bio), re, sys, os, subprocess, argparse, string, copy and numpy. If you do not have them installed, you can install them using the following command:

```
pip install biopython
```
You must find the correct name to install the package. All of them offer alternatives of installations.

### Installing

To install the package first unzip the folder:

```
tar -xzfcv complex_reconstr.tar.gz
```

Change the directory to the folder and execute setup script:

```
sudo python setup.py install
```

At the end, you will be able to use it as a package and also as a script.

## Running the program

To run the program use the following command:

```
python3 complex_reconstruction.py -i [inputs] -o [output] -v -vz
```
Make sure that if you execute the program as a script, have the modules homodimers, heterodimers and common_functions in the same directory as the script.

All options are optionals, but it is important to write inputs and outputs. 

You will find more information about the program, how to run it, tests and options explained in the written project.


## Authors

* **Lydia Fortea** - *Mainly in charge of the script*
* **Juan Luis Melero** - *Mainly in charge of the project*



