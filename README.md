# LAPackCPP

With this repo, i try to translate the original LAPACK Fortran routines (from here http://www.netlib.org/lapack/) to C++, by letting them first be automatically translated to C by the f2c.exe commandline tool (from here: https://www.netlib.org/f2c/) and then manually editing the results afterwards. The motivation was that at one point, i needed a solver routine for pentadiagonal linear systems for my RAPT library. After looking into various options:  

- the LinPack port here https://people.sc.fsu.edu/~jburkardt/cpp_src/linpack/linpack.html looked really good for my purpose but has the annoying LGPL license
- i wrote my own simple solver (without pivot-search) which looks pretty neat: https://github.com/RobinSchmidt/RS-MET/issues/225#issuecomment-462045482 but seemed to be numerically not-so-good
- projects like http://lapackpp.sourceforge.net/ are not really a translation but rather an interface which still requires the original code and a Fortran compiler (i cannot bother clients who use my library to install a Fortran compiler - i can't even bother myself to do that)
- more modern linear algebra libraries like Eigen have a - to my taste - too fancy and bloated interface. Others like Armadillo even have a dependency on the original LAPACK (and therefore on the availability of a Fortran compiler)

i decided that porting the original LAPACK routines is the way to go. I hope that, over time, this will grow into a useful library.

I just started with this, so at the moment, not very much has been translated yet. But i hope that over time, this repo may evolve to a more and more complete translation of this venerable, classic and awesome library. I don't know, how fast this will go - it will probably depend on what specific functionality from LAPACK i will need for my other projects (RAPT in particular).

For the licensing, i apply the same liberal licensing scheme as applies to the original Fortran code - anything else would be very unfair.

...tbc...
