=================
parallelized blat
=================
blat with multi-threads support
-------------------------------
.. image:: https://travis-ci.org/icebert/pblat.svg?branch=master
    :target: https://travis-ci.org/icebert/pblat


When the query file format is fasta, you can specify many threads to process it.
It can reduce run time linearly, and use almost equal memory as the original blat
program. This is useful when you blat a big query file to a huge reference like
human whole genome sequence.

The program is based on the original blat program which was written by Jim Kent.

pblat can run on Linux and Mac OS.

----

Install
---------------
To compile the source code, simply enter the source code directory in terminal
and issue the "make" command. When the compiling finished, the executable pblat
will locate in the same directory. Then it can be moved to where you want.

Run
---------------
To run with multiple threads, add option "-threads=<number of threads>" in the
command line.

Licence
---------------
pblat is modified from blat, the licence is the same as blat. The source code and
executables are freely available for academic, nonprofit and personal use. Commercial
licensing information is available on the Kent Informatics website.

----

Copyright (C) 2012 - 2017 Wang Meng

Contact me: wangm@mail.cbi.pku.edu.cn 
