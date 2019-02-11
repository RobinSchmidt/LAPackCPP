# LAPackCPP

With this repo, i try to translate the original LAPACK Fortran routines (from here http://www.netlib.org/lapack/) to C++, by letting them first be automatically translated to C by the f2c.exe commandline tool (from here: https://www.netlib.org/f2c/) and then manually editing the results afterwards. The motivation was that at one point, i needed a solver routine for pentadiagonal linear systems for my RAPT library. After looking into various options:  

- the LinPack port here https://people.sc.fsu.edu/~jburkardt/cpp_src/linpack/linpack.html looked really good for my purpose but has the annoying LGPL license
- i wrote my own simple solver (without pivot-search) which looks pretty neat: https://github.com/RobinSchmidt/RS-MET/issues/225#issuecomment-462045482 but seemed to be numerically not-so-good - and writing my own general band-diagonal solver with pivot search seemed to be moderately complex and i think, i have re-invented far too many wheels already in my life as a programmer
- projects like http://lapackpp.sourceforge.net/ are not really a translation but rather an interface which still requires the original code and a Fortran compiler (i cannot bother clients who use my library to install a Fortran compiler - i can't even bother myself to do that)
- more modern linear algebra libraries like Eigen have a - to my taste - too fancy interface. Curiously recurring template pattern? Great but really? I just wanted some raw and straightforward number crunching routines! :-P Others like Armadillo even have a dependency on the original LAPACK (and therefore on the availability of a Fortran compiler)
- it goes without saying that commercial, closed source solutions like Intel's MKL are not an option as a dependency for an open source project like my RAPT library (which is part of the RS-MET codebase, by the way)

i decided that porting the original LAPACK routines is the way to go. I hope that, over time, this will grow into a useful library. After all the fiddly setup work is done (which was, in fact, much less fiddly than i feared it to be :-)), i expect the rest to be mostly grunt work. A big pile of grunt work, though. But as a chinese saying goes: the beginning is half of the whole - or something.

I just started with this, so at the moment, not very much has been translated yet. But i hope that over time, this repo may evolve to a more and more complete translation of this venerable, classic and awesome library. I don't know, how fast this will go - it will probably depend on what specific functionality from LAPACK i will need for my other projects (RAPT in particular).

What bothers me about LAPACK (and also LINPACK, EISPACK, etc.) is that you have the same routines 4 times - for single and double precision and real and complex numbers. Code duplication is the mother of all evil. In my translation, i don't want to use daxpy, saxpy, caxpy and zaxpy - instead, i want to use a type-independent, i.e. templatized, axpy. Now the design goal of LAPACK was that the user may drop in replacements for the BLAS routines (like daxpy) with implementations that are optimized for the particular target hardware platform...i'm not yet sure how to solve the problem of using templatized BLAS routines while still alowing them to be replaced. To get the ball rolling, i'll first just templatize the reference BLAS implementations that come packaged with LAPACK. Later i may probably try some sort of compile-time dispatch in axpy with some sort of explicit specializations for daxpy etc. - which may then be replaced - we'll see...

I'm very indebted to the authors of the original Fortran sources:
http://www.netlib.org/lapack/contributor-list.html
for making this library available to the public under such liberal licensing conditions. Of course, i apply the same licensing scheme to this translation - anything else would be very unfair.

...tbc...
