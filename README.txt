This software requires zlib both regular and development versions to be installed on your machine.

To compile, from the top level directory, type: make
To test (after compiling), from the top level directory, type: make test

Under the main statgen repository, there are: 
* lib - the library code
** The library code is compiled into libStatGen.a which, after compiling is located at: statgen/lib/libStatGen.a
** After compiling, library headers can all be found in: statgen/lib/include
* src - the tools we developed
** After compling, the executables are found in: statgen/src/bin
* scripts - the scripts we developed

Makefiles
---------
statgen/Makefile.include should contain the definitions that you need for creating software using this library.

statgen/lib/Makefile.lib and statgen/lib/Makefile.src can be used as templates for creating Makefiles for new software.  If possible, just link to them.  They look for a file called Makefile.tool that should be all you need to update for your specific software.  (both Makefiles automatically include Makefile.include)

A similar setup should be used for test code, linking your Makefile to statgen/Makefile.test and defining Makefile.testtool for your specific test.

Other Notes
-----------
* Typically the .o files are compiled into their own directory called obj.
