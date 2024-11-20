
README file for CADS/UVS Stellar Spectrum Calculator
----------------------------------------------------

1. BUILD/INSTALL
----------------

  Follow these 4 steps:

  1. `cd' to the directory containing the package's source code and type
     `./configure' to configure the package for your system.

     Running `configure' takes a few seconds.  While running, it prints some
     messages telling which features it is checking for.

  2. Type `make' to compile the package.

  3. Type `make install' to install the programs and any data files and
     documentation. By default, the binary is copied to /usr/local/bin.
     You may change this destination by passing on the following
     argument to configure:

             --prefix=/your/chosen/destination

     Which would result in the binary in /your/chosen/destination/bin
     You may need to ensure that in $PATH to run the program.

  4. You can remove the program binaries and object files from the
     source code directory by typing `make clean'.

You may read more about configure script and others in the accompanying
file named `INSTALL'.


2. USAGE
--------
All input parameters are read from an input parameter file named
stellarspec_initparams.txt. The general usage pattern is:

  uvs_stellar_spectrum [options]

  Options are:
    -h help
    -g Generates a well documented default parameter file

Once you have the default parameter file, you can edit it as per your needs
and run the program without any arguments to generate the stellar spectrum.

NOTE: This program requires the kurucz model datafiles to simulate stellar 
spectrum. If you have installed it right, these files would be installed in 
<prefix>/share/cads/data/ and the default parameter file is aware of the
location. These datafiles are also available in the subdirectory named
"spectral_data" inside the source distribution as well. 


3. DOCUMENTATION
----------------
Documentation and test cases are available at
http://cads.iiap.res.in/software


4. LICENSE
----------
GPL [See the file COPYING.txt for details]


5. DISCLAIMER
-------------
You may encounter bugs in this software. If you do, please
report them. Your bug reports are valuable contributions,
since they allow us to notice and fix problems on machines
we don't have and/or remained un-noticed.


6. REPORTING BUGS
-----------------
You can register with the cads bug reporting tool:
http://cads.iiap.res.in/bugzilla/
and file a bug report.

If you are too lazy, drop in an email to: cads_AT_iiap.res.in

Either way, please include as many details as possible.

-----------------------------------------------------------
Reks, 01 July 2012. <reks_at_iiap.res.in>

