Normally, you can just enter the directory and run 
> make
This will create a binary named Multrec, which you can then execute.


If that does not work for you, here are the specifics.
Qt was used to develop this program, although there are no dependencies to the 
Qt library.  Qt was used for the convenience of its QtCreator IDE and its 
ability to generate a Makefile.  This Makefile was generated with the command

> qmake Multrec.pro CONFIG+=RELEASE

If you have Qt installed.  You should be able to run the above command, 
and then run 
> make

And otherwise, the simplest way is to include all .h and .cpp from this directory and use your favorite compiler, e.g. g++.
As mentioned before, there should be no dependencies to a fancy external library.
However, whichever environment is used, the user must ensure that the equivalent to the following Qt
compiler directives are handled : 

QMAKE_CXXFLAGS += -std=c++0x

The project requires the use of c++0x (for the unordered_set and the unordered_map classes).
