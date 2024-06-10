

README file for CADS/UVS Zodiacal light model software
------------------------------------------------------

0. Dependencies
----------------

UVS Zodiacal light model requires libnova to perform coordinate
transformation and stuff like that. If you are using linux, then
a binary package of libnova should be available in official OS
repository since popular sofware like kstars desktop planetarium
used libnova. You need to install libnova-devel package so that
we have all the header files for compiling this one.


1. BUILD/INSTALL
----------------

Compiling this code is not difficult task after all!
On unix/linux machines, you might even get away with
commands like:

cd src && cc *.c -o your_binary_file -lm -L/path/to/libnova -lnova

The standard route, however, is to read the file "INSTALL"
:-) Once you are done with it, follow these 4 steps:


  i. `cd' to the directory containing the package's source code and type
     `./configure' to configure the package for your system.  If you're
     using `csh' on an old version of System V, you might need to type
     `sh ./configure' instead to prevent `csh' from trying to execute
     `configure' itself.

     Running `configure' takes a few seconds.  While running, it prints
     some messages telling which features it is checking for. It will try
     to locate libnova distribution but if it fails, you can manually
     direct configure script to the location of your libnova by:
     
     ./configure --with-libnova=path/to/libnova
     
     Script will search libnova header files in:
         1. path/to/libnova
         2. path/to/libnova/include
         3. path/to/libnova/include/libnova
         
     and library (libnova.so or libnova.a) in
         1. path/to/libnova
         2. path/to/libnova/lib
         3. path/to/libnova/lib64
         
 ii. If configure ended successfully, Type `make' to compile the package.

iii. Type `make install' to install the programs and any data files and
     documentation. By default, the binary is copied to /usr/local/bin.
     You may change this destination by passing on the following
     argument to configure:

             --prefix=/your/chosen/destination

     Which would result in the binary in /your/chosen/destination/bin
     You may need to ensure that in $PATH to run the program.

 iv. You can remove the program binaries and object files from the
     source code directory by typing `make clean'.

You may read more about configure script and others in the accompanying
file named `INSTALL'.


2.USAGE
-------

If you have installed the package, you can invoke the binary: 

<prefix>/bin/uvs_zodiacal_model 

If you run the program without any arguments, and there is no parameter 
file (zodiacal_initparams.txt) in the working directory, the program will 
generate a new parameter file, which is pretty well documented internally. 
You may edit the values in parameter file. 

Program requires two data files -
         1. leinert_dist.txt
         2. zodiacal_spec.txt

These two files are packed along in the src subdirectory. If you have 
followed the standard installation procedure listed above, they will be
installed in <prefix>/share/cads/data/ and the default paramenter file is 
aware of the location. 

By default, output will be written to zodiacal_output.txt.


3. DOCUMENTATION
----------------
Documentation : cads.iiap.res.in/tools/zodiacalCalc/Documentation
Test reports  : Not yet :-)


4. LICENSE: GPL [See the file COPYING.txt for details]
-----------


5. DISCLAIMER
-------------
You may encounter bugs in this software. If you do, please
report them. Your bug reports are valuable contributions,
since they allow us to notice and fix problems on
machines/platforms we don't have, and/or remained un-noticed.


6. REPORTING BUGS
-----------------
You can register with the cads bug reporting tool: 
http://cads.iiap.res.in/bugzilla/ 
and file a bug report.

If you are too lazy, drop in an email to: cads_AT_iiap.res.in

Either way, please include as many details as possible. 

-----------------------------------------------------------
Reks, 22 June 2012 <reks_at_iiap.res.in>
      Last modified: 23 June 2012

