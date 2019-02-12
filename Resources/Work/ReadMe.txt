This directory is meant to be used as temporary working directory for 
converting Fortran .f files to C .c files via the f2c command line 
application. The original Fortran source functions are grabbed from:

lapack-3.8.0\BLAS\SRC   for BLAS sources
lapack-3.8.0\SRC        for LAPACK sources

the f2c Fortran to C converter can be donwloaded here:
http://www.netlib.org/f2c/
http://www.netlib.org/f2c/mswin/index.html Windows version

the fable Fortran to C++ converter can be downloaded here:
http://cci.lbl.gov/fable/
http://cci.lbl.gov/cctbx_build/
..i didn't try this yet - the installer is HUUUGE - WTF?



in the directory:
Lapackpp1.1a\include
there's a "modified" f2c.h file - such a file is missing in
f2c\libf2c
but i think it's supposed to be there as f2c.h is included by all code that
f2c generates. 

I think, i should first set up a fresh console-application project, drag in all 
the libf2c stuff and see, if i can get the converted daxpy function to run.

To convert .f files, copy the files and f2c.exe into Temp and run .\f2c *.f 
from the console while being in that Temp directory. This converts them
all at once (i got an error for xerbla.f, though and it generated an empty
xerbla.c file - but for a few others, it seemed to work - i may have to 
manually translate files where such errors occur - i hope such errors are very 
rare :-O). Or: the error was an unknown intrinsic functions - modifying the
Fortran source in a way that the missing function is treated as an external
function by replacing INTRINSIC by EXTERNAL in the xerbla.f fixed it. Of 
course, it means i must provide the missing function "len_trim__" in this case.
...we'll see how this works out...
btw. it alsoe works to use Start GIT bash here in Windows explorer and then
use ./f2c *.f (using the forward slash instead of the backslash)

After a rotuine has been converted, the conversion result goes into
Source/GeneratedByF2c and the original .f file goes into the Done folder. If 
an original routine has to be modified to be digestible for the translator, the
original version will be kept with the appendix _original appended to its 
filename. 


Maybe after a converted file, like "daxpy.c", works, make a templated version 
"taxpy.cpp" or just "axpy.cpp" ...or combine a lot of these functions into a 
single file. The comments for how the routine is supposed to be used go into 
the header file - make them doxygen compliant. Maybe wrap the C++ version into
a namespace - maybe use namespaces BLAS, LaPack, later maybe also LinPack, 
EisPack

Notes on real vs complex valued functions: i think, the complex case should be
used for conversion - it's more general than the real case, because many matrix 
formulas include a complex conjugation when in the real case it's just 
transposition - that conjugation just boils down to an identity for real 
inputs. -> for the already converted "d"-routines, check carefully, if they can
be used also for the complex case - for new translations, start with the 
complex version from the start.

already converted functions, that should be safe to complexify:
daxpy, dcopy, dger (but double-check!), dlaswp, dscal(check!), dswap, idamax



questionable:

dgbrfsx vs zgbrsfx: 

-the z-version does not translate - f2c produces an error related to the 
intrinsic "TRANSFER" function declaration - changing that to external does not 
fix the problem in this case (well, it fixes the error for the declartion, but 
there are still syntax errors in lines 650 and 658, where
ZLA_GBRFSX_EXTENDED is called with TRANSFER in the argument list

-in line 556, there's this assignment: NOTRAN = LSAME( TRANS, 'N' ) which is 
later used in a conditional to decide wether ZLA_GBRFSX_EXTENDED is called with
"...,COLEQU, C,..." or with "...,ROWEQU, R,..." - so NOTRAN probably has to with
matrix-transposition or not (not with "transfer or not" - whatever that means)

-the difference between dgbrfsx and zgbrsfx is that dgbrfsx calls 
ZLA_GBRFSX_EXTENDED with WORK( 1 ) where zgbrsfx uses 
TRANSFER (RWORK(1:2*N), (/ (ZERO, ZERO) /), N) as argument

-i need to figure out what that "TRANSFER" means - probably some sort of 
temporary data buffer or related function...maybe install a fortran compiler 
and try to compile the .f files
it's probably this one:
http://ergodic.ugr.es/cphysn/LECCIONES/FORTRAN/f77to90/a5.html#section11
"TRANSFER(SOURCE, MOULD, size) specifies that the physical representation of the first argument 
SOURCE shall be treated as if it had type and parameters as the second argument MOULD, but without 
converting it. The purpose is to give a possibility to move a quantity of a certain type via a 
routine that does not have exactly that data type."
...so that baiscally seems to be something like a reinterpret_cast? maybe it's because it's a 
newer (Fortran90) function that f2c does not yet know because it supports only Fortran77
...but why is it then a syntax error when we decalre it as EXTERNAL function?

-i think, in this statement TRANSFER (RWORK(1:2*N), (/ (ZERO, ZERO) /), N),
 RWORK(1:2*N) creates a copy of element 1..2N of the array RWORK ...and that 
 copy shall be reinterpreted as..



-my guess is that is has to do with promoting an array of real numbers to their
corresponding complex counterparts, i.e. just the same array with zeros 
interspersed for the imaginary parts

-in the translated dgbrfsx.c, the corresponding call of dla_gbrfsx_extended__ 
is in line 710

-dla_.. or zla_gbrfsx_extended__ is an iterative refinement procedure




dgbmv, 
dgbrfs (in contains an optiona to pass 'C' for complex conjugation),
dgbsv (the doc explicitly says a "real" system), 
dgbtf2, 
dgbtrf,
dgbtrs (has TRANS='C' option), 
dgemm (..has 'C' option - boils down to 'T'),
dgemv ('C' -> 'T'), 
dtbsv, 
dtrsm

compare with cooversion results of corresponding complex functions

not possible:


