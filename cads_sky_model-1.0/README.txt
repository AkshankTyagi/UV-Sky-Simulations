

README file for CADS sky_model software
-------------------------------------------


0. FEATURES
-----------

Sky model is a tool to predict the appearance of the sky at any wavelength.
It comprises three parts:

   i. Zodiacal light background (using Leinert's model)
  ii. Diffuse galactic light (from Galex data)
 iii. Stars (hipparcos catalog).

Provision is made for different filters.
Output will be in fits format
http://en.wikipedia.org/wiki/FITS


1. BUILD/INSTALL
----------------

Pre-requisites:

 a) CFITSIO
    CFITSIO is a library and headers to read and write FITS files.
    If you do not have cfitsio installed, get it from:

    http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html

 b) LIBNOVA
    libnova library and headers.
    Once you have cfitsio and libnova installed,
    you may proceed with compiling sky_model.

    Most linux distributions have cfitsio and libnova in their deb/yum/rpm
    repositories.

    Follow these 4 steps to compile and install the code:


  1. `cd' to the directory containing the package's source code and type
     `./configure' to configure the package for your system.

     If cfitsio library and header files are not in the standard path,
     then you may need to pass on the following argument to configure:

        --with-cfitsio=path/to/cfitsio

         Script will search fitsio header files in:
             1. path/to/cfitsio
             2. path/to/cfitsio/include
             3. path/to/cfitsio/include/cfitsio

        and library (libcfitsio.so or libcfitsio.a) in
             1. path/to/cfitsio
             2. path/to/cfitsio/lib
             3. path/to/cfitsio/lib64


     Similarly, for libnova (incase configure failed to detect it)

              --with-libnova=/path/to/libnova

     Running `configure' takes a few seconds. While running, it prints some
     messages telling which features it is checking for.

  2. Type `make' to compile the package.

  3. Type `make install' to install the programs and any data files and
     documentation. By default, the binary is copied to /usr/local/bin. You may
     change this destination by passing on the following argument to configure:

             --prefix=/your/chosen/destination

     Which would result in the binary in /your/chosen/destination/bin
     You may need to ensure that in $PATH to run the program.

  4. You can remove the program binaries and object files from the source code
     directory by typing `make clean'.

You may read more about configure script and others in the accompanying file
named `INSTALL.txt'.


2.USAGE
-------
This software is governed by an input parameter file called

diffuse_initparams.txt.

When you invoke the program, sky_model will generate a default parameter file
if the file is not present in the working directory. This file is well
documented internally and it should help you in choosing the input values.


3. DOCUMENTATION
----------------
If you have downloaded the source package, a detailed user manual
(skymodel_manual.pdf) can be found in the doc subdirectory after you've gone
through configure and make (section 1, see above).

If you installed a binary package, you will find the documentation alongside
this file.


4. LICENSE: GPL [See the file COPYING.txt for details]
-----------


5. DISCLAIMER
-------------
You may encounter bugs in this software. If you do, please report them. Your bug
reports are valuable contributions, since they allow us to notice and fix
problems on machines/platforms we don't have, and/or remained un-noticed.


6. REPORTING BUGS
-----------------
You can register with the cads bug reporting tool:
http://cads.iiap.res.in/bugzilla/
and file a bug report.

If you are too lazy, drop in an email to: cads_AT_iiap.res.in

Either way, please include as many details as possible.


-----------------------------------------------------------
Reks, 31 Oct 2012 <reks_at_iiap.res.in>
      Last modified: 08 Nov 2012
