#
#    @@@@@@@  @            @@@    @@@@@@@@             @@@@@@@        @@@@@  
#   @       @ @           @   @   @       @           @       @      @    @@ 
#   @         @          @     @  @       @                   @     @    @  @
#    @@@@@@@  @         @       @ @@@@@@@@              @@@@@@      @   @   @
#           @ @         @@@@@@@@@ @                    @            @  @    @
#   @       @ @         @       @ @                   @         @@@  @@    @ 
#    @@@@@@@  @@@@@@@@@ @       @ @                   @@@@@@@@@ @@@   @@@@@  
#                                                                          
#
#	This file defines the SLAP version 2.0 package and it's test routines.
#	To compile the SLAP version 2.0 quick check do the following:
#		1) Read the headder in the source file "slapqc.f" and
#		   make the appropriate system dependent changes.
#		2) Edit the makefile and change the FFLAGS, LFLAGS,
#		   OWNER, GROUP and LIBDIR macros for your system.
#		3) Compile the SLAP code and the driver with "make slapqc".
#		4) If you also wish to install the SLAP code as a Unix library
#		   then do "make install"
#	To run the quick check simply type "slapqc" at the command line and
#	give an input value for KPRINT of 0, 1, 2, or 3 in I1 format.
#	KPRINT=2 recomended.  Results of the test are printed to the STANDARD
#	OUT.
#
FFLAGS=-g -ffpa
LFLAGS=-L$(HOME)/lib -lmksF -lm
OWNER=seager
GROUP=src
LIBDIR=$(HOME)/lib

NSRCS=sir.f scg.f scgn.f sbcg.f scgs.f somn.f sgmres.f smset.f smvops.f slaputil.f
NOBJS=sir.o scg.o scgn.o sbcg.o scgs.o somn.o sgmres.o smset.o smvops.o slaputil.o
SRCS=slapqc.f tstutil.f
USRCS=xersla.f mach.f blas.f
#UOBJS=xersla.o mach.o blas.o
OBJS=$(NOBJS) $(UOBJS)

.f.o:
	$(FC) ${FFLAGS} -c $<

slapqc : slapqc.o ${OBJS}
	$(FC) ${FFLAGS} slapqc.o ${OBJS} -o slapqc ${LFLAGS}
slapcmp : slapcmp.o tstutil.o ORIG/slap.o ${OBJS}
	$(FC) ${FFLAGS} slapcmp.o tstutil.o ${OBJS} ORIG/slap.o -o slapcmp ${LFLAGS}
ORIG/slap.o : ORIG/slap.f
	(cd ORIG; $(FC) $(FFLAGS) -c slap.f)

slapqc.o : slapqc.f
slapcmp.o : slapcmp.f
tstutil.o : tstutil.f

sir.o : sir.f
scg.o : scg.f
scgn.o : scgn.f
sbcg.o : sbcg.f
scgs.o : scgs.f
somn.o : somn.f
sgmres.o : sgmres.f
smset.o : smset.f
smvops.o : smvops.f
slaputil.o : slaputil.f
xersla.o : xersla.f
mach.o : mach.f
blas.o : blas.f

DOBJS = dir.o dcg.o dcgn.o dbcg.o dcgs.o domn.o dgmres.o\
	dmvops.o dmset.o dlaputil.o ${UOBJS}
DSRCS = dir.f dcg.f dcgn.f dbcg.f dcgs.f domn.f dgmres.f\
	dmvops.f dmset.f dlaputil.f
dlapqc : dlapqc.o ${DOBJS}
	$(FC) ${FFLAGS} -o dlapqc dlapqc.o ${DOBJS} ${LFLAGS}
dlapqc.o : dlapqc.f
dir.o : dir.f
dcg.o : dcg.f
dmeth.o : dmeth.f
dgmres.o : dgmres.f
dmvops.o : dmvops.f
dlaputil.o : dlaputil.f
#
#	This list makes the Unix library "libslap.a"
lib libslap.a : $(NOBJS)
	rm -f libslap.a
	ar cr libslap.a $(NOBJS)
	ranlib libslap.a

#
#	This list installs the library in $(LIBDIR)
install : libslap.a
	install -o $(OWNER) -g $(GROUP) -c -m 644 libslap.a $(LIBDIR); ranlib $(LIBDIR)/libslap.a
	
#
#	This list generates the distribution "shar" files	
slap.shar : READ.ME makefile $(NSRCS)
	shar -a READ.ME makefile $(NSRCS) > slap.shar
slapqc.shar :  $(SRCS) $(USRCS) slapdoc.f
	shar -a slapdoc.f $(SRCS) $(USRCS) > slapqc.shar
dlap.shar : READ.ME makefile $(DSRCS)
	shar -a READ.ME makefile $(DSRCS) > dlap.shar
dlapqc.shar : dlapdoc.f dlapqc.f $(USRCS)
	shar -a dlapdoc.f dlapqc.f $(USRCS) > dlapqc.shar
#
#	This makes the SLATEC slap.f file.
slap.f : slapdoc.f $(NSRCS)
	cat slapdoc.f $(NSRCS) > slap.f
slapqcqc.f : $(SRCS)
	cat $(SRCS) > slapqcqc.f
dlap.f : dlapdoc.f $(DSRCS)
	cat dlapdoc.f $(DSRCS) > dlap.f

slaptomatlab : slaptomatlab.o
	$(FC) $(FFLAGS) slaptomatlab.o -o slaptomatlab -lslap $(LFLAGS)
clean : 
	rm -f *.o slaptst gmon.out core *.a

TAGS : slaptst.f slap.f newslap.f scg.f silu.f sbcg.f
	etags slaptst.f slap.f newslap.f scg.f silu.f sbcg.f
