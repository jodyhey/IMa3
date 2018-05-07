/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* This is the main header file for pretty much everything, nearly all files include this file */


/* sections in this file */
/******** INCLUDES *******************/
/******** COMPILATION MACROS *******************/
/******** SIMPLE DEFINITIONS MACRO ************/
/******** CODING MACROS *******************/
/****** ENUMERATED THINGS ************/
/***** GLOBAL STRUCTURES  ************/
/*******     GLOBAL VARIABLES     *********/
/***** GLOBAL FUNCTION PROTOTYPES *********/


/*************************************/
/******** INCLUDES *******************/
/*************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stack>
#include <sstream>


/***********************************************/
/******** COMPILATION MACROS *******************/
/***********************************************/

/*  compilation macro notes

  MPI_ENABLED must be defined if compiling for MPI
     this definition is usually handled at the command line or in the development environment

  When debugging
      DEBUG should be defined, usually handled at the command line or in the development environment
      if DEBUG is defined,  TURNONCHECKS will also be defined, which leads to lots of things being checked and makes it quite slow
      To run in debug mode without all the checks,  uncomment #define TURNONCHECKS

  When the code is going to be run on the testbed
   - compile with STDTEST defined, this will cause RANDOM_NUMBERS_FROM_FILE to be defined
   -  if STDTEST is defined TURNONCHECKS will be undefined so as not to take too long if running standard tests in debug mode

  Whenever code has reached the point that it is a definite release version:
   - increment  IMA3RELEASEVERSION
   - uncomment  #define IMA3RELEASE
   - put this copy of the code (with IMA3RELEASE defined)  in a safe spot, clearly labeled as release code
*/

//#define MPI_ENABLED  // usually done in compiler settings or the command line
//#undef MPI_ENABLED

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#ifdef XMLOUTPUT
#include "tinyxml.h"
#include "tinystr.h"

#endif

#define IMA3RELEASEVERSION  "0.0"  // update only when a release is made
//#define IMA3RELEASE   // uncomment  if this is release code,  use sparingly with updates of IMA3RELEASEVERSION


#ifdef IMA3RELEASE
#undef INDEVELOPMENT
#else
#undef IMA3RELEASE
#define INDEVELOPMENT
#endif

//#define STDTEST //  in visual studio, uncomment this to compile for testbed


#ifndef STDTEST //  STDTEST can be defined at compile time or in the code
#undef STDTEST
#undef RANDOM_NUMBERS_FROM_FILE
#else
#define RANDOM_NUMBERS_FROM_FILE  // STDTEST and RANDOM_NUMBERS_FROM_FILE should both be defined or both not
#endif

/* #define RANDOM_NUMBERS_FROM_FILE   coded by Janeen, used to be called SANITY_TEST
 * for testing only. Usually used on testbed.  Avoids use of random number generator to make things more repeatable.
 * read a random number from a file "randNums" which contains a list of random numbers (in ascii)
 * then convert the ascii to unsigned long and return.*/

#if defined(_DEBUG) && !defined(DEBUG) // defines DEBUG by default,  _DEBUG is visual studio macro // will be overridden by compiling with -D NDEBUG in release mode
#define DEBUG
#endif /* _DEBUG or DEBUG */
//#define NDEBUG  usage
// invoked NDEBUG will cause program to compile without use of asserts.
// visual studio compilers include NDEBUG in release configurations when compiling
//when compiling under linux NDEBUG must be specifically invoked
//  e.g.  mpicxx *.cpp -D MPI_ENABLED -D NDEBUG -O3 -o IMa3
//on linux using -D NDEBUG makes program about 5% smaller
// and it runs about 8% faster


#ifdef DEBUG
#define TURNONCHECKS
#endif


#undef TURNONCHECKS   // turn off debugging check functions (mostly in treeprint.cpp and chainprint.cpp)


#ifdef STDTEST
#undef TURNONCHECKS   // turn off debugging check when running standard tests
#endif
#ifdef CLTURNONCHECKS // allow for command line turn on of TURNONCHECKS, i.e. -D CLTURNONCHECKS // debugging files must be included if using this
#define TURNONCHECKS
#endif


// vld - visual leak dectector for finding memory leaks,  only runs under visual studio on windows
// output is not appearing in the IDE,  so need to look in  memory_leak_report.txt
#ifdef  _MSC_VER  // a built-in visual studio macro
//#define USEVLD   // turn on visual leak detector vld  // usually not on // works better in 32bit mode
#undef USEVLD  //usually this is undefined,  but sometimes turn it on in debugging mode to check
#else
#undef USEVLD  // can't use vld outside windows,  or with MPI
#endif /*  _MSC_VER */

#ifdef USEVLD
// vld is a leak detection library,  only invoked in MVC++ debug mode and if vld is installed
// libs for both 32 or 64 are provided, but must set this in configuration
#include <vld.h> // for windows, a useful library for checking memory leaks under microsoft visual c++ 10, otherwise do not use this
#undef MPI_ENABLED  // does not work to run mpi with vld
#endif


/********* DIRENT USAGE **************/
/* Microsoft compilers use _findfirst and _findnext to search directories for files.
Most other compilers use a library called dirent.  It is not a formal standard
but is very widespread.  dirent does not have _findfirst or _findnext but instead
uses functions called  opendir() closedir() readdir() and rewinddir().
In order for this code to be portable,  I've included a file called sdirent that
maps the usual dirent functions calls onto the microsoft functions.  If the compiler
is not microsoft, then this file is not included, and the regular dirent is used */

#ifdef _MSC_VER
#include <io.h>
#include "msdirent.hpp"
#include <errno.h>              /* _findfirst and _findnext set errno iff they return -1 so these must be included */
#else
#include <dirent.h>
#endif

/*  forceinline was used for  labelgtree(), which uses it for HKY model likelihood calculations
but it does not seem to work when compiled on linux, changed _forceinline  to inline
 Ended up just commented this out, was causing problems compiling under linux, not sure if it was very very helpful JH 6/8/2016
#ifndef _MSC_VER
#ifndef __forceinline
#define __forceinline __attribute__((__always_inline__)) inline
#endif
#endif
*/

/* microsft visual studio compiler stuff */
/* disable deprecate warning  C4996 */
#ifdef _MSC_VER
#pragma warning( disable : 4996)
#define _CRT_SECURE_NO_WARNINGS
#endif /* _MSC_VER */

/* By turning on conditional compiler definition MORESTABLE,
 * we use LogDiff or LogSum2. Otherwise, we use the original form of equations.
 * Search codes for LogDiff and LogSum2 for example. */
#define MORESTABLE

/* splittime updateing  */
#define DO_NWUPDATE    // only turn off when debugging something
//#undef DO_NWUPDATE   // use undef to turn this update off
#define DO_RY1UPDATE   // only turn off when debugging something
//#undef DO_RY1UPDATE  // use undef to turn this update off


/**********************************************/
/******** SIMPLE DEFINITION MACROS ************/
/**********************************************/

/* CONSTANTS */
#define MAXLOCI  201            //100
#define MAXGENES 1000           // maximum sample size for a locus
#define DEFAULTNUMCHAINS 1
#define SWAPDIST 7 //10         // maximum distance in chain array of two chains with betas being swapped, not sure how much this matters
#define MINNUMCHAINSPERPROCESSOR 2 // changed to 2 after some testing on 10/4/2017  4  // must have at least 4 for swaps within chains
#define MAXPOPS 9 // 10     set this back to 9 , have not tested it with this many        // MAXPOPS cannot exceed 10 because the treestring functions assume populations and nodes are represented by single integers
#define MAXPOPS_PHYLOGENYESTIMATION 9  // actually 8,  but 9 prevents a crash  (with ghost it becomes 9) for larger numbers the list of possble trees just gets too large.
#define MAXPERIODS (MAXPOPS+1)
#define MAXTREEPOPS  (2*MAXPOPS - 1)
#define MAXUPPOPS    2           /* added for hidden genealogy stuff.  maximum number of populations that can come together in the tree at one time */
#define MAXLINKED 15            // largest number of linked loci with the same genealogy - each neads its own mutation rate
#define TIMEMAX 1000000.0       // no branch can have a bottom time greater than this
#define MINPARAMVAL 0.0000001   //0.0001      // smallest parameter value for parameters that are in the MCMC
#define LOGFACTMAX  100001   // large number for logarithm of factorial,  needed in various contexts,  but must be big to handle large sample sizes if they come up
#define MINPRIORFROMHYPERPRIOR 0.0001  // only used for selecting new hyperpriorvalue, as values to close to 0 cause integration problems.
#define MINPOPSFORDISTANCE 6 // minimum number of npops for which to use topology distance value rather than just topology number
/*
  MIGMAX is the maximum number of migration events on an edge.
  MIGMAX has a pretty big effect on how much memory is used.
  a max of 20 per edge seems like a lot but more will be needed when there are more populations
  and we will want to be able to have lots of migration in heated models
  particularly when updating topologies
  settle on 50 6/7/2016
*/
#define MIGMAX 50              // number of migration events that can be in .mig[] and .mighg[] arrays
#define MIGARRAYSIZE  (MIGMAX + 1)  //  +1 allows for last migration element at end of array to be terminator,  with -1 for time
#define MPRIORMIN 0.000001      /* small value for setting upper bound on migration to (near) zero */
#define KAPPAMAX  100           /* maximum value of HKY parameter */
#define RECORDINTERVALDEFAULT 10     // default # of steps between recording values of things - used to call record()
#define SAMPLEGENEALOGYINTERVALDEFAULT 100  // default # of steps between recording information about the genealogies
#define MINSTRLENGTH  3         // minimum allowed number of STR repeats, so users don't use data that doesn't fit // JH changed to 3  11/29/2010
#define GRIDSIZE 1000           // # of bins in histogram
#define TRENDDIM 500            // number of points saved for the trendline plots
#define MARGIN2DGRIDSIZE 50     // # of bins on each axis in 2D histograms
#define PROFILEGRIDSIZE 50      // 100 //30  // # points along single dimension profile curve
#define MAXGENEALOGIESTOSAVE 300000 // the max to save in ram during a run - some bug,  setting this to 500,000 caseus memory management errors
#define PRINTINTDEFAULT 10000   // default # steps between writing to the screen
#define MAXLOADFILES 500        // max # of files with genealogy information that can be loaded
#define UMAX  10000.0           // highest value for u scalars  - can differ by UMAX^2 fold
#define HMAX  20                // highest value for h scalars, if HMAX=20 it gives a range of ratios for h scalars from 1/20  to 20  (i.e. they can differ by up to 400 fold)
#define TIMEPRIORMULTIPLIER 10  // useful for plotting TMRCAs
#define EXPOMIGPLOTSCALE  10 // changed to 10 when adding hyperprior   20 //with exponential prior on migration, this number times the given prior mean value sets peak search interval and plot scale
#define FORCEREJECTIONCONSTANT  -1e100      // a very low value, indicates failure of IS model
#define REJECTMIGRATIONPROPOSAL  -1000000000.0      // a very low value, indicates failure of migration proposal because too many were proposed
#define NUMTARRAYBINS  100      // number of bins for multidimensional peak estimation of splitting times
#define ARBCHAIN 0   // used as a chain index for those cases when any chain will do because they are all the same for the info referred to
#define HEADNODE 0 //the headnode is the cpu that does writing of output files
#define BURNTRENDSTARTDELAYDEFAULT  10000
#define GENERATIONTIMEDEFAULT   1
#define DEFAULTNUMGENEALOGIES 10000
#define CHECKINTERVALSTEPS 500 //1000
#define MPRIORFRACFORHG  0.1  // used when hiddenoptions[UPDATEMRATEFORHGUPDATE]==0 to fix the migration rate used when updating hg
#define NUMPCOUNTARRAYS 10 // # maximum # of phylogeny count arrays in the burntrend file

/* some things that were once command line options but are now just turned on */
#define PRINTDEMOGHIST 1 //- print out distributions on demographic scales - requires mutation rates and generation times


/*constants used for specifying character array lengths*/
#define PARAMSTRLEN  12         /* length of string of parameter name */
#define PARAMSTRLENSHORT  7     /* length of string of parameter name */
#define UPDATELABELLEN 18       /* length of string for an update name */
#define FNSIZE  1000 // increased 8/22/2016 to allow for long paths  500           // max file name length
#define POPTREESTRINGLENGTHMAX  100
#define POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION  60  /* with npops=8 + outgroup,  the maximum treestring length is apparently 56  */
#define COMMANDLINESTRINGLENGTHMAX 1001
#define BANNERMAXLENGTH 102  // just for printing headings in output file
#define NAMELENGTH 151          // max length of population names and locus names
#define GENENAMELENGTH 10          /* a gene name can be up to 10 */
#define POPTOPOLOGYSEQUENCELENGTH 1000 // starting length of array of saved topology and distance values
#define FPSTRIMAXLENGTH 200000 // max length of input file header needs to be really big,  sprintf is used a lot on a string of this length, a bit risky

/* autocorrelation estimation constants */
#define AUTOCTERMS 12           // the number of lag values for which autocorrelations are recorded
#define CHECKAUTOCWAIT 10000   /* 100000 */     // number of steps before recording autocorrelation values
#define AUTOCINT 1000           // interval between measurements for autocorrelations - cannot be changed easily
#define AUTOCNEXTARRAYLENGTH 1000       /*  = (largest value in autoc_checkstep[] divided by AUTOCINT) */
#define AUTOCCUTOFF 100  /*500 */       /* minimum number of measurements to have for autocorrelation before printing results */
#define AUTOCSTEPSCALAR 1//50//5       /* scale over which autocorrelation is measured, 1 means the scale is steps */

// some alternative values when debugging autoc stuff to speed things up
//#define CHECKAUTOCWAIT 10 //10000
//#define AUTOCCUTOFF 2//100
//#define BURNTRENDSTARTDELAYDEFAULT 10// 10000
//#define CHECKINTERVALSTEPS 10 //1000

/******************************************/
/******** CODING MACROS *******************/
/******************************************/
/* nan, inf macros */
  /* some C libraries have things like this built in, but others don't,  so these are included.  for isnan end name in '_' so as not to conflict with same
  names in some math.h files */
#define isnan_(x) ((x) != (x))
#define isnotnan(x) ((x)==(x))
#define isinf_DBL(x) (((x) > DBL_MAX) || ((x) < -DBL_MAX))  // only for doubles
#define isninf_DBL(x) ((x) < -DBL_MAX)  // only for doubles
#define ispinf_DBL(x) ((x) > DBL_MAX)  // only for doubles
#define isnotinf_DBL(x) (((x) <= DBL_MAX) && ((x) >= -DBL_MAX))  // only for doubles
#define isinf_FLT(x) (((x) > FLT_MAX) || ((x) < -FLT_MAX))  // only for float
#define isninf_FLT(x) ((x) < -FLT_MAX)  // only for floats
#define ispinf_FLT(x) ((x) > FLT_MAX)  // only for floats
#define isnotinf_FLT(x) (((x) <= FLT_MAX) && ((x) >= -FLT_MAX))  // only for floats



#ifdef TURNONCHECKS  // these are used for checking if floating point values are 'equal' ,  very crude
//#define checkfloatsclose(x,y) ( ((x)==(y)) || ((fabs((x) - (y)) / DMAX(fabs(x), fabs(y))) < (1e-7) ) ) // evaluates to 1 if floats x and y are near each other
#define checkfloatsclose(x,y) ( ((x)==(y)) || ( (fabs(x) - fabs(y)) < 1e-12 )  || ( ((x)==0.0)  && (fabs(y) < 1e-12) ) || ( ((y)==0.0)  && (fabs(x) < 1e-12) ) || ((fabs((x) - (y)) / DMAX(fabs(x), fabs(y))) < (1e-7) )  )
#define checkfloatsortaclose(x,y) ( ((x)==(y)) || ((fabs((x) - (y)) / DMAX(fabs(x), fabs(y))) < (1e-4) ) ) // evaluates to 1 if floats x and y are near each other
#endif

#define MYDBL_MAX DBL_MAX/1e10 // these avoid some over- under-flow issues, I think.  Not used much. May not be necessary
#define MYDBL_MIN DBL_MIN * 1e10
#define IM_BESSI_MIN (-1e+100) /* very small value that would rejct any update with stepwise muation model */
#define POSROUND(a) (long) ((a)+0.5)    // simple rounding to positive integers
#define INTEGERROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))        // simple rounding to integer
#define ODD(a) ( (((a) & 1) == 1) ? 1 : 0 )     /* for nonnegative integers , if odd returns 1, else 0 */
#define ISEVEN(x) (!ODD(x))

#define FP fprintf(outfile,     // handy way to avoid retyping the same thing again and again
#define SP *fpstri += sprintf(&fpstr[*fpstri],  // widely used to add text to the long string (called 'fpstr') that goes at the beginning of the output files

//#define f_close(a)  fclose(a); (a) = NULL       //regular close() does not make the pointer null,  but it is useful to have this set to null

#define FCLOSE(fp)  ((fp) ? fclose(fp) : 0, (fp) = 0)  // found this online //regular close() does not make the pointer null,  but it is useful to have this set to null

#define XFREE(p) do { free((p)); (p) = NULL; } while(0)  // replaces free(p) to declare the pointer NULL

#define LOG10  2.3025850929940456840
#define LOG2 0.69314718055994530941723212146
#define LOG2HALF 0.34657359027997265470861606073   // half of LOG2, used in update_t_NW
#define LOG_DBL_MAX 7.0978271289338397e+02

/* these are supposedly bulletproof macros suggested by Melissa Hibnuz
 * but still, do not nest calls to these min and max macros
  - later found that these were copied for nrutil.h */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static double dminarg1, dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static float minarg1, minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
static long lmaxarg1, lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))
static long lminarg1, lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))
static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

/* Macro LogDiff approximately computes the value of v of the following
 * equation:
 *  v = Log
 * exp(v) = exp(a) - exp(b)
 * where a > b.
 *  typically use when a and b are logarithms and we want the logarithm of
 *  the difference between what they are logarithms of.
 *  Use if you want the logarithm of the difference, ax-bx
 *  and a is the log of ax  and b is the log of bx
 */
 // added parentheses around (a) - (b)  4/27/2017 just to be sure
#define LogDiff(v,a,b) \
    if ((a) <= (b)) \
      IM_err (IMERR_LOGDIFF, " inside LogDiff() macro, a %lf  b %lf",(a),(b)); \
    if ( ((a) - (b)) < 7.0978271289338397e+02) \
      (v) = (b) + log (exp ((a) - (b)) - 1.0); \
    else \
      (v) = (a)

/* Macro LogSum2 approximately computes the value of v of the following
 * equation:
 * exp(v) = exp(a) + exp(b) */
 // added parentheses around (a) - (b)  4/27/2017 just to be sure
#define LogSum2(v,a,b) \
    if ((a) > (b)) \
      if ( ((a) - (b)) < 7.0978271289338397e+02) \
        (v) = (b) + log (exp ((a) - (b)) + 1.0); \
      else \
        (v) = (a); \
    else \
      if ( ((b) - (a)) < 7.0978271289338397e+02) \
        (v) = (a) + log (exp ((b) - (a)) + 1.0); \
      else \
        (v) = (b)

/* Macros for using IM_errloc */
#define STRINGIFY(x) #x  // stringifcation, whatever x is it is replaced with a string literal of x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)



/* SET MACROS */
/* below are some routines for using sets ala pascal.
	I got this from "C reference manual" by Harbison p 175. useful for sets
	of small numbers from 0(inclusive) to the size of an unsigned long
	integer (32, exclusive). A number is in the set if the bit in the place
	position of that number (plus 1, because we permit the number 0 to be in
	the set) is a 1. If the bit is zero then it is not in the set.  Thus
	the numbers 0-31 can be in the set. */
/* any program using this must also declare typedef unsigned long set; */
/* changed to unsigned short to save memory - won't need more than 16 populations */
typedef unsigned int SET;
#define SET_BITS 32             // make sure this is correct for the compiler
#define EMPTYSET   ((SET) 0)    /* produces a SET of value 0 */
#define ISEMPTY(setj)  ((setj) == 0)    /*true if setj has no elements */
#define SINGLESET(i)  (((SET) 1) << (i))        /* makes a set of the integer i, CAREFUL!! only works for i values from 0 thru 31 !! */
#define SETADD(setj,i)  ((setj) | SINGLESET (i))        /* setj with i added to it */
#define ISELEMENT(i,setj)   (SINGLESET((i))  & (setj) ) /*true if i in the setj */
#define INTERSECT(set1,set2)   ((set1) & (set2))
#define UNION(set1,set2)   ((set1) | (set2))
#define SETDIFF(set1,set2)   ((set1) ^ (set2))
#define SETREMOVE(setj,i)  ((setj) ^ SINGLESET(i))      /* setj with i removed from it , be sure that is is aleady in it */

#define FORALL(i,setj) \
  for ((i) = 0; (i) < SET_BITS; ++(i)) \
  if (ISELEMENT ((i), (setj)))

      /* forall : permits an action to be applied to all members of a set
         for example
         int k;
         forall(k,z) printf("%d ",k);   */

#define GETLOW(setj)  ( (setj) & -(setj) )
      /*returns the lowest member of the set, as a set
         can be used in conjuntion with setdiff to reduce sets an item at a time */

/*  for counting # of elements in a SET,  must copy this somewhere
int cardinality(SET x){
	int count = 0;
	while (x != EMPTYSET) {
		x ^= (x & -x); ++count;
		}
	return(count);
	}  */


/*************************************/
/****** ENUMERATED THINGS ************/
/*************************************/

/* run modes  runmode is set to one of these

this is an internal organizing system to help keep track of some things

these modes are only referenced in the code and are not used in output

they pertain to 4 variables:
sampling phylogenies vs not:    0,1  versus 2,3,4,5,6,
sample priors from hyperprior vs not:  0,2,5  versus  1,3,4,6
use hidden genealogy with fixed topology vs not: 5,6 vs 0,1,2,3,4
load genealogies vs not:   4 vs 0,1,2,3,5,6

*/


enum{
POPTREEHYPERPRIORmode0 = 0,  //sample poptree topology, sample priors given hyperprior   -j03
POPTREEmode1 = 1,  //sample poptree topology         -j0
GHYPERPRIORmode2 = 2,   // topology is given, sample priors given hyperprior,    -j3
Gmode3 = 3,  //topology is given, sample genealogies, aka M mode   no model options (no -j)
LOADGmode4	= 4, //topology is given, load genealogies, aka L mode  -r0   (no -j)
HGHYPERPRIORmode5 = 5, //topology is given, sample  priors given hyperprior, use hidden genealogy -jh2 -j3
HGmode6 = 6, //topology is given, sample genealogies, use hidden genealogies - like mode 3 with with hidden genealogies  -jh2  (no -j)
NUMRUNMODES = 7};

/*  hidden options using -jh flag  (formlerly called jheyoptions[])
  these do not appear on the command line
  developer can use them
  some are invoked indirectly from command line
  e.g. -j0 invokes hiddenoptions[HIDDENGENEALOGY]
	   SWINPUTOPTION   alternate way of inputing SW data when there is just one SW portion for a given genealogy
    WRITEMIGRATIONNAME  -jh1filename  (i.e. the numberid symbol for WRITEMIGRATIONNAME is immediately followed by
      two integers,  the from and to population numbers,  and then immediately followed by the filename.
      e.g. '-jh104migrate0to4.out')
    HIDDENGENEALOGY - use hidden genealogy for updating genealogy and population tree
    GSAMPINFOEXTRA  - save a .ti file in -j0 mode,  include probhgg and population id info
    UPDATEMRATEFORHGUPDATE - use a migration rate for hidden genealogy updating drawn from migration prior, update this in mcmc
    FIXMUTATIONSCALARUPATES  mutation rate scalars are fixed to relative mutation scalar values, as found in input file
    CALCGEWEKEZ  calculate gewekez stat when updating phylogenies
    PRINT2NMASCIICURVES  prints ascii curves for 2Nm,  takes some time so off by default
    SKIPMOSTUSCALAROUTPUT  stop printing of  most mutation scalar stuff, use -jha  or -jhA on command line
    STOPMOSTINTERVALOUTPUT  stop writing most interval output to the screen  use -jhb or -jhB on command line

  some options tried, but proved unnecessary or harmful:
   - SUPERHEATING    in chain swapping raised the mh ratio term to the largest of the two betas,  did not work at all
   - tried having an update in which topology was never updated within a chain,  but only updated by swapping in chains
  and then at the beginning all tree numbers were used evenly across chains.   Tried it to see if bugs were in the tree
  updating code.  Dumb idea. Did not work.  This forces a uniform distribution of trees among chains,  and the different
  chains no longer share the same state space.

    MIGRATEPAIRHYPERPRIOR = 9, // hyperprior for pairs of migration rates  // little point, just use regular hyperprior
    POPTREETOPOLOGYUPDATE=3,  11/1/2017 stop using this and only use modeloptions[POPTREETOPOLOGYUPDATE]


  */
enum
{ SWINPUTOPTION = 0,
  WRITEMIGRATIONNAME=1,   // 1_22_2018  not sure what this is for
  HIDDENGENEALOGY=2,
  GSAMPINFOEXTRA=3,
  UPDATEMRATEFORHGUPDATE=4,
  NOMUTATIONSCALARUPATES = 5, // no mutation rate scalar updating
  FIXMUTATIONSCALARUPATES = 6, // mutation rate scalars are fixed to relative mutation scalar values, as found in input file
  CALCGEWEKEZ = 7, //calculate gewekez stat when updating phylogenies
  PRINT2NMASCIICURVES = 8,
  SKIPMOSTUSCALAROUTPUT = 9, // stop printing of  most mutation scalar stuff, use -jha  or -jhA on command line
  STOPMOSTINTERVALOUTPUT = 10,// stop writing most interval output to the screen  use -jhb or -jhB on command line
  READOLDMCFFILE = 11, // -jhC on command line, read mcf files generated before 1/17/2018  use -jhb or -jhC on command line
  HIDDENOPTIONNUMBER=12};

/* calcoptions -c
    DONTCALCLIKELIHOODMUTATION - don't calculate p(Data|G)  if set to 1 then data likelihood functions return a constant
	   MUTATIONPRIORRANGE - include pior ranges on mutation rate scalrs, as included in input file
    FINDJOINTPOSTERIOR - find joint posterior for full model and/or nested models loaded in nestedmodelfile
    LOADPRIORSFROMFILE - read a file with prior's for each parameter
*/
enum
{
  DONTCALCLIKELIHOODMUTATION = 0,
  MUTATIONPRIORRANGE = 1,
  FINDJOINTPOSTERIOR = 2,
  LOADPRIORSFROMFILE = 3,
  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  CALCMARGINALLIKELIHOOD = 4,  // as of 6/17/2017 no confidence that this is useful.  it does not appear in help menu
  DONTCALCGENEALOGYPRIOR = 5, // JH added 4/27/2017,  causes the treeweight and integrate functions to do nothing , not checked much
  CALOPTIONSNUMBER = 6
  /*, FINDJOINTPEAK, PRINTPEAKFINDDETAILS */
};

/*  modeloptions  using the -j flag
 POPTREETOPOLOGYUPDATE
 ADDGHOSTPOP - add a non-sampled ghost population to the model
 EXPOMIGRATIONPRIOR - exponential distribution on migration prior -m gives mean value
 POPSIZEANDMIGRATEHYPERPRIOR - use hyperpriors for population size and migration parameters
	NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
	SINGLEMIGRATIONBOTHDIRECTIONS
	MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
	PARAMETERSBYPERIOD  - every population size and migratino parameter applies for only 1 period
	NOMIGRATION  - no migration
 ONEMIGRATIONPARAMETER - just one migration parameter for the entire model
 ONEPOPSIZEPARAMETER
    these were removed:
      ANCESTRALPOPSIZESSETTOLARGESTSAMPLED - ancestral population size parameter is set to be the same as its largest sampled descendant population
	    SINGLEPOPULATION - single population, no migration, no population tree, just a single theta      ---> replaced by npops == 1

	*/

enum
{
  POPTREETOPOLOGYUPDATE = 0,
  ADDGHOSTPOP = 1,  //1
  EXPOMIGRATIONPRIOR, //2
  POPSIZEANDMIGRATEHYPERPRIOR,//3
  NOMIGBETWEENNONSISTERS, //4
  SINGLEMIGRATIONBOTHDIRECTIONS,//5
  MIGRATIONBETWEENSAMPLED, //6
  PARAMETERSBYPERIOD,//7
  NOMIGRATION, //8
  ONEMIGRATIONPARAMETER,//9
  ONEPOPSIZEPARAMETER,//10
  MODELOPTIONSNUMBER
  };

/* outputoptions -p
   DONTPRINTASCIITREND - don't print ascii trend lines
   DONTPRINTASCIICURVE - don't print ascii curves
   PRINTTMRCA - print a table of the TMRCA distribution
   THISTDIVIDEBYPRIOR - print histogram table of splittime values, divided by the prior
   NOPOPMIGPARAMHIST - do not print the histograms and plots for the population migration rate parameters
   PARAMGREATERTHAN - print probabilities of pairwise greater than probabilities
   MIGRATEHIST - print out distributions of numbers of migration events
   PRINTJOINTTEST  - print the joint estimate for splitting times,  for more then 1 splitting time,  not in LOADRUN mode
   // options removed just turn these on by default
    reuse these names as constants and replace e.g. outputoptions[PRINTDEMOGHIST] with  PRINTDEMOGHIST
    POPMIGPARAMHIST - print the histograms and plots for the population migration rate parameters
      this slows things down quite a bit at the end of a run as it can take a while to do peakfinding in 2Nm space
    PRINTDEMOGHIST - print out distributions on demographic scales - requires mutation rates and generation times
*/
enum
{
  DONTPRINTASCIITREND=0,
  DONTPRINTASCIICURVE=1,
  PRINTTMRCA=2,
  THISTDIVIDEBYPRIOR=3,
  NOPOPMIGPARAMHIST = 4,
  PARAMGREATERTHAN =5,
  MIGRATEHIST=6,
  PRINTJOINTTEST=7,
  OUTPUTOPTIONSNUMBER=8
};

/* runoptions -r
	LOADRUN  - load genealogy and t information saved in a previous run
	DONTSAVEGENEALOGIES  - save gtree information for analytical integration
 SAVEMCSTATEFILE - save the state of the markvon chain in a file
 LOADMCSTATE - load a previously saved markov chain state
 PRINTMUTATIONUPDATESTOSCREEN - print all of the mutation related updates (scalars, kappas, A values) and values to stdout during the run
    PRINTBURNTREND - print trend blot during burnin
 MCFLOADONLYSTATESPACE - do not load all the data on sampled values and acceptance rates
 SAVELOADSAMEMCFFILE - load and save mcf files based on the outfile name, if files not present at start, begin run without them
    - skip burn period upon subsequent starts,  if mcf files are present
    - use for reusing same command line

*/
enum
{
  LOADRUN = 0,
  DONTSAVEGENEALOGIES, //1
  SAVEMCSTATEFILE,  //2
  LOADMCSTATE,  //3
  PRINTMUTATIONUPDATESTOSCREEN, //4
  PRINTBURNTREND,  //5
  MCFLOADONLYSTATESPACE,  //6
  SAVELOADSAMEMCFFILE,  //7
  UNIQUEBURNRUN, //8
  RUNOPTIONSNUMBER
};

/* mutation models
 * ---------------
 * INFINITESITES
 * HKY
 * STEPWISE    - one or more linked stepwise loci, no IS portion
 * JOINT_IS_SW - one IS portion and one or more STEPWISE portions
 */
enum
{
  INFINITESITES,
  HKY,
  STEPWISE,
  JOINT_IS_SW,
  IM_MODEL_NUMBER
};

/*  different ways that the length of the run can be specified using the '-B' and '-L' flags
	TIMESTEPS - time counted in # of steps
    TIMEINF  - run until first charcater in 'IMrun' is not 'y'
	INDEFINITE - used for burnin */
enum
{
  TIMESTEPS,
  TIMEINF,
  INDEFINITE
};

/* heating modes */
/* 6/14/2017  stopped using all but HGEOMETRIC
  made use of -hfg on command line optional */
enum
{
  HLINEAR,
  HGEOMETRIC,
  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  HEVEN,
  HFULL   //JH added to deal with hidden genealogies and topology updating - use -hfs
};

/* kinds of acceptance of genealogy update */
enum
{
  IM_UPDATE_GENEALOGY_ANY,
  IM_UPDATE_GENEALOGY_TOPOLOGY,
  IM_UPDATE_GENEALOGY_TMRCA,
  IM_UPDATE_GENEALOGY_NUMBER
};

/* kinds of acceptance of poptree updates */
enum
{
  IM_UPDATE_POPTREE_ANY,
  IM_UPDATE_POPTREE_TOPOLOGY,
  IM_UPDATE_POPTREE_TMRCA,
  IM_UPDATE_POPTREE_NUMBER
};


/* update method of tree */
enum
{
  IM_UPDATE_TREE_WCH,
  IM_UPDATE_TREE_NUMBER
};

/* update methods for migration rate priors */

#define NUM_PRIOR_UPDATE_RECORD_TYPES 1 // 1 if just MPRIOR_ALONEE (/to skip recording stuff about the the "WithTopol" update,  which is not very useful) else 2.

enum
{
  PRIOR_UPDATE,
};

/* error code enumeration - these are used by IM_err() and must correspond to the text messages in simerrmsg[]  in utilities.cpp*/
enum
{
  IMERR_SUCCESS = 0,
  IMERR_READFILEOPENFAIL = 1,
  IMERR_MEM = 2,
  IMERR_TIFILE = 3,
  IMERR_CREATEFILEFAIL = 4,
  IMERR_APPENDFILEFAIL = 5,
  IMERR_OUTPUTFILECHECK = 6,
  IMERR_READFILEFAIL = 7,

  IMERR_COMMANDLINEINCOMPAT = 8,
  IMERR_MISSINGCOMMANDINFO = 9,
  IMERR_COMMANDLINEFORMAT = 10,
  IMERR_COMMANDLINEHEATINGTERMS = 11,

  IMERR_MISSINGPOPSTRING = 12,
  IMERR_POPTREESTRINGFAIL = 13,
  IMERR_OPTIONNUMBERFAIL = 14,

  IMERR_INFILEFAIL_NLOCI = 15,
  IMERR_NESTEDMODELLSPECIFYLFAIL = 16,

  IMERR_MUTSCALARPRIORRANGEFAIL = 18,
  IMERR_MUTSCALEPRODUCTFAIL = 19,
  IMERR_MUTSCALEPROBLEM = 20,

  IMERR_TOOMANYMIG = 21,
  IMERR_MIGARRAYTOOBIG = 22,
  IMERR_HPD95 = 23,
  IMERR_HYPERPROB = 24,
  IMERR_DATAREADOVERRUN = 25,
  IMERR_DATAERROR = 26,
  IMERR_LOCUSERROR  = 27,

  IMERR_MCFREADFAIL = 28,
  IMERR_MCFSPLITTIMEPROB = 29,
  IMERR_MCFWRITEFAIL = 30,

  IMERR_ROOTTIMEMAXFAIL = 33,

  IMERR_INPUTFILEINVALID = 35,
  IMERR_INFINITESITESFAIL = 36,

  IMERR_SWCHECK = 38,

  IMERR_LOWERGAMMA = 41,
  IMERR_UPPERGAMMA = 42,
  IMERR_LOGDIFF = 43,

  IMERR_MULTITPRIOR = 45,
  IMERR_PRIORFILEVALS = 46,
  IMERR_MIGRATIONPRIOR0 = 47,

  IMERR_NUMERICALRECIPES = 51,
  IMERR_GSL = 52,
  IMERR_ASSERT = 55,
  IMERR_GENENAME = 56,
  IMERR_ASN = 57,
  IMERR_MIGPATHPROB = 60,
  IMERR_MISCPROBPROBLEM = 62,
  IMERR_CHAINNUM = 70,
  IMERR_MISCELLANEOUS = 72
};



typedef char strn[PARAMSTRLEN];
typedef char strnl[UPDATELABELLEN];

/*************************************/
/***** GLOBAL STRUCTURES  ************/
/************************************/

/*Kernigan Ritchie dictionary item */
struct dictionary_node_kr {
    struct dictionary_node_kr *next; /* next entry in chain */
    char *name; /* defined name */
    double priorval;
};


typedef struct  // used in  sort_and_print_alltreestrings()
{
  char treestr[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  char treestrnoghost[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  int count;
  double freqset1;
  double freqset2;
  double freqall;
  double  ppcp; // product of the posterior clade probabilities
  int origi; // original index
}  foralltreestringsort;

struct topolseq // holds the entire sampled sequences of tree numbers and Robinson Foulds distances (from tree 0).  Kept only on cpu with rank 0.
{
  int  *vals; // the current tree number
  double *disvals; // Robinson Foulds distances between tree 0 and the current tree
  int currentlength;
  int maxlength;
};


/* extendnum is used for calls to eexp() */
struct extendnum
{
  double m;
  int z;
};

/* eevent contains info needed to calculate the mean and variance, correlations and autocorrelations over course of run */
struct eevent
{
  double s;                     /* sum of times */
  double ss;                    /* in case need to sum a second variable */
  double s2;                    /* sum of squares of times */
  int n;                        /* number of events */
};

//typedef struct eevent im_eevent;

/* for calculating autocorrelations */

struct autoc
{
  struct eevent cov;
  struct eevent var[2];
  double vals[AUTOCNEXTARRAYLENGTH];
};

//typedef struct autoc im_autoc;

// used for calculating 90% HPD intervals
struct hlists
{
  double v;
  double p;
};

//typedef struct hlists im_hlists;

// terms needed for loci with HKY mutation model
struct hkyinfo
{
  double **frac;
  double **newfrac;
  double *scalefactor;
  double *oldscalefactor;
};

/* a little structure for sorting population tree node numbers and times together */
struct pnodetime
{
  int ptreepos;
  double time;
};

//typedef struct hkyinfo im_hkyinfo;

/*Main data structures: edge, locus, parameter */
/* an edge is a branch in the genealogy
each edge gets a number.  the n sequences are numbers 0 thru n-1.  The remaining n-1 edges are numbered after that.
For edge i.
up[2] - contains the numbers of the edges to which edge i connects to  (i.e. the numbers of its descendants in the genealogy).
down - contains the number of the edge to which edge i connects down to  (i.e. the immediate ancestor in the genealogy
time - contains the time at the bottom of the edge, that is the time at the top of the ancestral edge (down).
      the root edge has a negative value for time.
*mig - pointer to arrya that contains the times of migration, and identity of pops migrated to, on edge i, these times are on the same scale as time,
   and thus must be less than the population splitting time
cmm - number of migration events in array (position in array of last migration event +1)
mut - used for labeling under INFINITESITES
A - an allele state for stepwise or other allelic model.  This is the state of the node at the top.
dlike - the likelihood of the distance between A, the allele state at the top of the node, and the allele state at the top
of the down edge
pop - contains the population that the edge is in at its top, which may be different than
    the one it is in at the bottom, depending on migration.
frac, newfrac, scalefactor, oldscalefactor - used for the HKY mutation model
*/
struct migstruct
{
  double mt; /* for the time                                                */
  int    mp; /* for the population a migration went to (in the coalescent); */
};

//typedef struct migstruct im_migstruct;

/* struct edgmiginfo holds info about the genealogy edges being updated
 * info is different than what is easily available in the genealogy itself
 * particularly in the arrays mtimeavail and mp which hold the times and number
 * of migration events in each period for the branch.
 * Comment on [[e]]: time period in which the edge ends  e >= b,  DO NOT CONFUSE
 * these values with poptree b and e  which refer to populations.
 *
 * When we move a branch, we could move either only the branch or the branch
 * with its sister. We may know whether a sister branch is involved in a branch
 * using member edgeid. If edgeid for a sister branch is negative, then we may
 * assume that the moving branch is not attached to a root.
 *
 * We may set the following members before using an edgemiginfo variable: li,
 * sisid, edgeid, upt, dnt, pop, temppop, fpop. By calling function
 * [[fillmiginfoperiods]] we determine the values of b, e, mtimeavail, and
 * mtall. We set mp, mpall, and mig by migration simulation.
 * During migration simulations it is necessary to determine which population
 * the edge is currently in.  temppop is used for this */
struct edgemiginfo
{
  /* by manual */
  int li;
  int sisid;
  int edgeid;
  double upt;
  double dnt;
  int pop;
  int temppop;
  int fpop;
    /* by fillmiginfoperiods */
  int b;                        /* time period in which the edge begins */
  int e;
  double *mtimeavail;
  double mtall;
  /* by migration simulation or copying */
  int *mp;  // migration count in period
  int mpall;
  struct migstruct mig[MIGARRAYSIZE];
  int cmm;   // totalnumber of migration events in mig
};


struct edge
{
  /* key edge members */
  int up[2];             /* daughter edge IDs: -1 for leaves       */
  int down;              /* down edge ID: -1 for root edge         */
  double time;           /* time at the bottom! of the edge                */
  /* migration arrays  mig and mighg
    cmm and cmmhg are the respective counts
    the last migration event in mig[] is at position cmm-1
    and the last migration event in mighg[] is at position cmmhg-1
    cmm can equal MIGMAX,  but not exceed it
    same goes for cmmhg
    at mig[cmm] there should be a false event with time of -1 i.e. mig[cmm].mt == -1
    at mighg[cmmhg] there should be a false event with time of -1 i.e. mighg[cmmhg].mt == -1
    these are indicators of being at the end of the array
    these indicators require that the length of mig and mighg should be MIGMAX + 1
  */
  struct migstruct mig[MIGARRAYSIZE]; /* migration events                       */
  struct migstruct mighg[MIGARRAYSIZE]; /* hgstuff migration events                       */
  int cmm;               /*number of migration events in the array */
  int cmmhg;              /* number of migration events in the array mighg */
  /* supplementary members */
  int pop;               /* population the edge is in at its top   */
  int pophg; /*hgstuff */
  int mut;               /* number of mutations on the edge        */
  int *A;
  double *dlikeA;
  struct hkyinfo hkyi;   /* only use if the mutation model for the */
                         /* locus is HKY                           */
  int i;                 /* identifier of individual whose         */
  int ei;                /* index at gtree for saving purposes     */
                         /* changed from gi to ei  1/9/09 is this  */
                         /* used?                                  */
  //char exist;    not in use 5/12/2016  /* 'F' if edge is detached from a  genealogy                              */
  int *seq;              /* sequences at the top of edge           */
  int fpop;
  int fpophg; /*hgstuff */
};

/* struct gtreeevent  is used in treeweight() which sorts all of the events in a genealogy */
struct gtreeevent
{
  double time;
  int pop;     /* population the edge is in when the event happens        */
  int topop;   /* population the edge goes to if it is a migration event  */
  int cmt;     /* coalesce =0 or  migrate=1  or process reaches divT = -1 */
  int periodi; /* period the event is in                                  */
               /* noticed 2/16/09  could probably delete periodi as it    */
               /* seems not to get used                                   */
};

struct upairlist
{
  int l;
  int u;
};

struct priorvalues
{
  double min; // always 0 ?   only for uniform priors
  double max;  // only for uniform priors
  double expomean;  // only for exponential priors
};

struct plotpoint                /* (x axis are the possible values,  y axis is the counts) */
{
  double x;
  double y;
};

struct update_rate_calc
{
  double accp;
  double tries;
};

struct acceptancerate
{
  struct update_rate_calc *au_perturb;
};

/* weightposition holds information on which values in the cc, fc, mc
 * and fm arrays of genealogy_weights to sum for doing the integration */
struct weightposition
{
  int n;                        /* # of terms summed to calculate the weight                 */
  int *p;                       /* period - array of length n                                */
  int *r;                       /* row  - array of length n                                  */
  int *c;                       /* col  - array of length n (only use for migration weights) */
};

/* updatescalarinfo for setting the value used in proposals for the distribution used for sliding distances when updating population trees */
/* careful - can't change this without also changing MPI_updatescalar */
struct updatescalarinfo
{
  double updatescalarval;  // the current scalar
  double updatescalaradjustval;  //  amount to adjust the scalar
  double updatescalarvalmin;
  double updatescalarvalmax;
  double targetupdaterate;  // target update rate
  int recentaccp;
  int recenttries;
  int allaccp;
  int alltries;
  int numattemptscheck;   // the number of attempts to use for calculating the recent acceptance rate
};

struct migdir
{
  int from;
  int to;
};

struct i_param
{
  struct priorvalues pr;
  struct plotpoint *xy;         /* if used points to an array of gridize * elements, each a plotpoint */
  strn str;
  int b;                        /* period when this parameter first appears */
  int e;                        /* period when this ends */
  struct weightposition wp;     /* info on where the weights are for
                                 * calculating the integration */
  /* for migration hyperprior stuff, not used for popsize params */
  struct migdir md;
  strn descstr;  // holds a string about the relevant descendant populations,  used only with hyperpriors
  int dir;  //dir is 0 if 'md.from' pops are on the left side of descstr and 'md.to' pops are on the right, else dir is 1
};

/*
struct value_record  for recording stuff needed for posterior plots, trends and ESS
all purpose
there are many differentn instances of these
they are only used on the head node because that is the node that writes output
so with mpi  there is a lot of MPI_Send and MPI_Recv of the data that gets putinto a value_record

char  str[PARAMSTRLN];  name
struct plotpoint *xy;   recorded values that are can later be plotted,an array of gridize elements, each a plotpoint
double *trend;    record of trend
struct  autoc ac[AUTOCTERMS];  autocorrelation calculations
int beforemin;   number of values  below the prior
int aftermax;    number of values  above the prior
*/
struct value_record
{
  strnl str;
  strn strshort;
  int do_xyplot;
  int do_logplot;
  struct plotpoint *xy;
  struct priorvalues plotrange;
  double plotrescale;           // 1 unless needed for some reason
  int do_trend;
  double *trend;
  int do_autoc;
  struct autoc ac[AUTOCTERMS];
  int beforemin;
  int aftermax;
};

/*
1/9/09
 struct chainstate_record_updates_and_values

 this structure contains a bunch of stuff related to something that can be measured about the state of the Chain 0
 examples include splitting times, mutation rate scalars,  features of genealogies etc etc
 An instance of chainstate_record_updates_and_values can hold:
	info on priors and updating window width
	numbers of update types
	acceptance rates for each type
	pointers to another structure (value_record) that hold information on recorded values

 components of chainstate_record_updates_and_values:
 char str[PARAMSTRLEN]  name of the thing being update
 struct priorvalues pr;
 double win    window width, may be used for updating
 int num_uptypes   number of different update types
 char **upnames   array of length num_uptypes of names of different update types
 struct update_rate_calc *upinf   - array of length num_uptypes of calculations of update rates
 int num_vals    Number of different measurements made -  for some things this will be zero as often we just want to know about update rates, not values
 struct  value_record  *v    pointer to array of length num_vals
*/
struct chainstate_record_updates_and_values
{
  strn str;
  struct priorvalues pr;
  double win;
  int num_uptypes;
  strnl *upnames;
  struct update_rate_calc *upinf;
  int num_vals; // set to 0 if *v is not going to be used
  struct value_record *v;
};

struct popedge
{
  /* int id;    turns out we don't need id# as the position in the list acts as the id # */
  int numup;                    /* number of descendant nodes                 */
  int up[MAXUPPOPS];           /*jh changed to a fixed array length when adding hidden genealogies*/ /* the identities of the descendant nodes */
  int down;                     /* the down node                              */

  /* int ti; the time interval at the bottom of the branch    */
  double time;                  /* the time on the branch to the down node  (not time interval, but time)  */
  int b;                        /* time period in which the population starts */
  int e;                        /* time period that begins where the population
                                 * ends, e > b    e == -1  if the population is
                                 * the root */
};

struct probcalc
{
  double *qintegrate;           /*  integration of theta terms                   */
  double *mintegrate;           /* integration of migration terms                */
  double pdg;                   /* probability of data given all genealogies     */
  double probg;                 /* prior prob of current genealogy  based on
                                 * integration over parameter                    */
  double probhgg;              /* HGSTUFF prior prob of hidden genealogy given genealogy*/
};

/* genealogy_weights holds the information needed to do integrations
 * of p(Genealogy) for the coalescent-related quantities (cc, hcc, fc),
 * there is just an array with positions corresponding to population numbers
 * for the migration-related quantitites there is a 3d array,
 * with each layer being a period in the poptree, and within each period
 * there is a 2D array of migration rates.
 * The indexing of these does not follow population numbers,
 * but rather follows the position of the population numbers that is
 * given in plist */
struct genealogy_weights
{
  int **cc;                     /* coalescent counts for pops in period i            */
  double **hcc;                 /* inheritance weights for population in period i    */
  double **fc;                  /*coalescent weights for pops in period i            */
  int ***mc;                    /* migration counts between populations in period i,
                                 * the positions in layer k of mc (i.e. mc[k] follow
                                 * the same order and listing in plist[k]            */
  double ***fm;                 /* migration weights between populations in period i */
};


struct locus
{
  char name[NAMELENGTH]; // name can be longer than gNames
  int pairs[MAXGENES];
  char gNames[MAXGENES][GENENAMELENGTH+1];
  int numgenesknown;   // 5/17/2017 don't know what this is for
  int numgenesunknown;  // 5/17/2017 don't know what this is for
    int numgenes;   // sample size,  the total number of sampled gene copies at a locus
  int samppop[MAXPOPS];
  int numlines;
  int model;                    /* the overall mutation model for the locus */
  double hval;
  int nlinked;                  /* # of linked portions  = 1 + # linked SW portions */
  int nAlinked;                 /* # of linked SW portions */
  int **A;                      /* points to microsat allele values, if nlinked > 1,  2D array */
  int maxA[MAXLINKED];          /* array of maximum allele size */
  int minA[MAXLINKED];          /* array of minimum allele size */
  int umodel[MAXLINKED];        /* array of locus specific model */
  int numbases;
  int numsites;
  int totsites;       /* added this 5/15/09,  fixed bad HKY bug, not sure why it was missing*/
  int **seq;
  int *mult;
  int *badsite;

  // mutation rate per year values and priors
  double uperyear_vals[MAXLINKED];
  struct priorvalues uperyear_prior[MAXLINKED];
  int uii[MAXLINKED];

  // records for mutation scalars
  struct chainstate_record_updates_and_values *u_rec;
  struct chainstate_record_updates_and_values *kappa_rec;

  // records for microsat allele state updates
  struct chainstate_record_updates_and_values *A_rec;

  // records for genealogy measurements
  struct chainstate_record_updates_and_values *g_rec;

 // struct chainstate_record_updates_and_values *a_rec; does not appear to be used
};

/* renamed form struct locus to struct genealogy  on 1/9/09
 moved a bunch of stuff out to the new struct locus

*gtree  pointer to the genealogy
root - edge number of root
mignum - current # of migration events in gtree
roottime - time of root of genetree
length - total length of gtree
tlength -  length of gtree except for time in last period, used for calculating a migration rate for updating genealogies
pdg - current total P(X|G)
pdg_a array of P(X|G) for each part of locus
*gtree  pointer to array of branch info
genealogy_weights - array of info that is saved when the genealogy is sampled
pi - used for HKY model
*/

struct genealogy
{
  /* genealogy stuff - different for each chain */
  struct edge *gtree;
  int root;             /* Should be updated in each genealogy update.                            */
  int mignum;           /* Is computed in treeweight.                                             */
  double roottime;      /* Should be updated in each genealogy update.                            */
  double length;        /* Is computed in treeweight.                                             */
  double tlength;       /* Is computed in treeweight.                                             */
  double pdg;           /* the total probability of the data given the genealogy (sum of pdg_a)   */
  double *pdg_a;        /* points to an array of of length nlinked  used for multiple linked loci */
  double *uvals;        /* points to an array of length nlinkded - the mutation scalar values     */
  double kappaval;      /* value of HKY parameter if needed */
  double pi[4];         /* for HKY model */
  double hilike;
#ifdef DEBUG
  double hiprob; // debugging 8/18/2016
#endif
  struct genealogy_weights gweight;
  double asn;
  int *mut;
  double hgprob; // probability of hidden genealogy given genealogy
  double mhg;  // the current sampled migration rate for hg updates // fixed if hiddenoptions[UPDATEMRATEFORHGUPDATE]==0
};

struct chain
{
  char chainpoptreestring[POPTREESTRINGLENGTHMAX];
  SET periodset[MAXPERIODS];
  int addpop[MAXPERIODS];
  int droppops[MAXPERIODS][2];
  struct i_param *itheta; // moved here for hg stuff
  struct i_param *imig;  // moved here for hg stuff
  int **plist;
  int **ancplist;  /*hgstuff */
  int poptreenum; /*hgstuff */
  struct popedge *poptree;
  double *tvals;                // array of time values
  char name[10];                /* chain name or id */
  struct updatescalarinfo branchslideinfo;   /* for poptree branch sliding updates */
  struct updatescalarinfo *RYwidthinfo; /* for Rannala Yang updates */
  struct updatescalarinfo *NWwidthinfo; /* for Neilsen Wakeley Yang updates */
  struct genealogy_weights allgweight;  // accumulated stuff used to calculate p(G)
  struct probcalc allpcalc;     // prob(Data|G) and prob(G) stuff
  struct genealogy *G;          // points to an array of genealogy
  int currallbetapos;  // jh added 6/29/2016 current position in allbetas of the beta value associated with this chain

//hyperprior stuff
   SET descendantpops[MAXTREEPOPS];
      // the ancestral pops currently in the tree  for  0 < i < npops descdendantpops[i] is just a SET made of  i
    	 // for npops <= i < numtreepops  descdendantpops[i] is the SET of all the populations descendant from node i
      // descendantpops[numtreepops-1] is always the full SET of sampled populations

  // popsize parameter allocations if modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]
    double *qhpriors; //  length 2^n, index position is the set of subpopulations with that prior
	   // make a SET of subpopulations, subpopset,   then qpriors[ (int) subpopset] gives the prior
  	 // position 0 cast as a set is the null set,  so this will not be used
   	// position 2^npops - 1   cast as a set is the full set, used for the ancestral pop

  //migration parameter allocations if using hyperprior
    struct dictionary_node_kr **mltorhpriors; // used for imig[i].dir == 0  (i.e. md.from is on left side of imig[i].descstr)
    struct dictionary_node_kr **mrtolhpriors; // used for imig[i].dir == 1  (i.e. md.to is on left side of imig[i].descstr)

};

/******************************************/
/*******     GLOBAL VARIABLES     *********/
/******************************************/

/* if GLOBVARS is defined then gextern is ignored. This
causes the corresponding variable declarations to be definitions.
GLOBVARS is only defined in ima_main.cpp.  If GLOBVARS is not defined,
as is the case at the beginning of the other files, then gextern
gets replaced by extern and the variable declarations are simply
declarations */

#ifdef GLOBVARS
#define gextern
#else /*  */
#define gextern extern
#endif /*  */

gextern char **debug_ti_addinfo; //strings to contain stuff to add to the ti file to look at genealogies when sampling trees
gextern int runmode; // takes on one of the values from 0 to NUMRUNMODES
gextern double thetaprior, mprior, tprior;  // by default this are set to -1 to begin with, they may be set using command line or priorfile
gextern double hyperprior_expo_m_mean; //set equal to mprior, for simplicity
gextern double expo_m_mean; //set equal to mprior, for simplicity
gextern double hyperprior_uniform_m_max; // set equal to mprior for simplicity
gextern double hyperprior_uniform_q_max; // set equal to mprior for simplicity
gextern double m_max; // set equal to mprior for simplicity
gextern double q_max; // set equal to thetaprior for simplicity
gextern double tperiodpriors[MAXPOPS-1];
gextern struct chain **C;  //points to an array of pointers to chains
gextern struct chainstate_record_updates_and_values *T;
gextern struct chainstate_record_updates_and_values *poptreeuinfo;   // only used with topology updating
gextern struct chainstate_record_updates_and_values *mh;  // only used with a hyperprior on migration rate parameters
gextern struct chainstate_record_updates_and_values *mhnit; // only used with a hyperprior on migration rate parameters
gextern struct chainstate_record_updates_and_values *qh;  // only used with a hyperprior on popsize rate parameters
gextern struct chainstate_record_updates_and_values *qhnit; // only used with a hyperprior on popsize rate parameters
gextern struct locus *L;

gextern struct value_record *lpgpd_v; // record for likelihood measurements
gextern struct value_record **migration_counts;   // used if outputoptions[MIGRATEHIST]
//gextern struct value_record **migration_counts_times; // no longer in use
gextern int step;
gextern int nurates;
gextern int nkappas;
gextern int numchainspp;  // numchains per process (i.e. per cpu)
gextern int numchainstotal;  // numchainspp * numprocesses
gextern int numprocesses; // # cpus
gextern int nloci;
gextern int npops;
gextern int numtreepops;        /* number of distinct populations in the tree */
gextern int numpopsizeparams;   /* number of distinct population size parameters in the model */
gextern int nummigrateparams;   /* number of distinct migration rate parameters in the model */
gextern int nummigrateparampairs; //nummigrateparams/2
gextern int numpopsets; /* number of subsets of sampled populations, including empty and full  numpopsets-2 is the number that can actually occur  */
gextern int nummigdirs;

gextern int hiddenoptions[HIDDENOPTIONNUMBER];
gextern int modeloptions[MODELOPTIONSNUMBER];
gextern int calcoptions[CALOPTIONSNUMBER];
gextern int outputoptions[OUTPUTOPTIONSNUMBER];
gextern int runoptions[RUNOPTIONSNUMBER];
gextern double *beta;
gextern double *allbetas; // an array to hold all the betas, regardless of cpu, should never change once its initialized
gextern int genealogysamples;
gextern int burninsteps;
gextern int runsteps;
gextern int numpriormcfruns;
gextern time_t totaltime;
gextern int mcf_was_read;
gextern int somestepwise;
gextern int countuprior;
gextern int counturateperyear;
gextern struct upairlist ul[2 * MAXLOCI];       /* listing of mutation rate scalars  not clear how big to make it,  because not clear how many portions each locus will have */
gextern int lastperiodnumber;   // just numsplittimes, but it comes up a lot
gextern int numsplittimes;      // same as lastperiodnumber, but it comes up a lot
gextern float **gsampinf;
/* position markers for gsampinf */
gextern int gsamp_ccp;
gextern int gsamp_fcp;
gextern int gsamp_hccp;
gextern int gsamp_mcp;
gextern int gsamp_fmp;
gextern int gsamp_qip;
gextern int gsamp_mip;
gextern int gsamp_pdgp;
gextern int gsamp_probgp;
gextern int gsamp_tp;
gextern double logfact[LOGFACTMAX];
gextern struct weightposition nomigrationchecklist;
gextern struct extendnum *eexpsum;

gextern int total_numgenes;  /* sum of sample sizes across all loci */
gextern double max_m_from_priorfile;
gextern double *nnminus1; // moved this here, out of update_gtree_common.cpp  when adding hidden genealogy stuff
gextern char **alltreestrings;  /* hgstuff holds all possible population tree strings - including ghosts if ADDGHOSTPOP is invoked */
gextern char **alltreestrings_noghost;  /* all possible population tree strings but without ghost */
gextern int  *poptopologycounts;  /* hgstuff,  holds observed counts of different population trees */
gextern  unsigned short *RFtreedis;  // used in getting distances between trees
gextern int numpoptopologies;   /* hgstuff  total number of possible distinct topologies*/
gextern struct edgemiginfo oldedgemig;
gextern struct edgemiginfo oldsismig;
gextern struct edgemiginfo newedgemig;
gextern struct edgemiginfo newsismig;
gextern int *poptopologyproposedlist; // checklist of all trees proposed
gextern struct topolseq poptopologysequence;
gextern int totaltopolupdates;
gextern int chain0topolupdates;
gextern int chain0topolswaps;
gextern char popnames[MAXPOPS][NAMELENGTH];
gextern char *poptreenewickstring;
gextern double *topologypriors;  // holds log of prior for each topology
gextern int usetopologypriors;
gextern double *popsizeprior_fromfile;
gextern double **mprior_fromfile;
gextern int poptopologiessampled;
gextern char startpoptreestring[POPTREESTRINGLENGTHMAX];

// hyperprior stuff
gextern struct i_param *holdimig;
gextern char **poppairs;
  	 //list of all possible poppairs as strings,  also extern in update_migration_priors.cpp
	   //length numdistinctpopulationpairs[]
	   //never changes
gextern int hashsize;

//These are set in start()
gextern int doRYupdate;
gextern int doNWupdate;
gextern int domutationscalarupdate;
gextern int time_update_type_count;
gextern int update_type_RY;
gextern int update_type_NW;

gextern int scanfval;  // return value by scanf functions, to avoid warnings mostly
gextern char* fgetval;  // return value by fget functions, to avoid warnings mostly

#ifdef TURNONCHECKS
gextern int currentid_debug;  // turn this on when debugging as it makes it easier
#endif


/******************************************/
/***** GLOBAL FUNCTION PROTOTYPES *********/
/******************************************/

/**** GLOBAL FUNCTIONS IN ima_main.cpp ****/
int whichiscoldchain(void);

/**** GLOBAL FUNCTIONS IN autoc.cpp ****/
void init_autoc_pointers (void);
void free_autoc_pointers (void);
void checkautoc (int start_autocorrelations, int burndone, int burninsteps, int currentid);
void callprintautoctable (FILE * outto/*, int step*/);

/**** GLOBAL FUNCTIONS IN update_poptree.cpp ****/

void init_change_poptree(char topologypriorinfostring[]);
void free_change_poptree(void);
int change_poptree (int ci,int *trytopolchange, int *topolchange, int *trytmrcachange, int *tmrcachange, int topologychangeallowed);

/**** GLOBAL FUNCTIONS IN update_hg.cpp ****/

void makegenealogy_from_hiddengenealogy(int ci,int li);
double prob_hg_given_g(int ci,int li);
int  update_hidden_genealogy(int ci,  int li, int *topolchange, int *tmrcachange);

/**** GLOBAL FUNCTIONS IN alltreestrings.cpp ****/
void filldescendantpops(int ci);
int makepairstring(SET d0,SET d1, char s[]);
int pairhash(char *str);
int fillmigratepairs(void);
char **allocalltreestrings(void);
void freepoptreestringarrays(void);
int buildpoptreestringarray(void);
void printnewickstring(FILE * outto, char *ps, double *tvals, int ghostintree);
void getnodecounts(char **na,char *tstr,int count,int *nc, int *nunique);
double calcppcp(char **na,char *tstr, int *nc, int totaltreecount,int nunique);
int foralltreestringsort_comp (const void * a, const void * b);
void printallpoptreesamples (FILE * outfile, int *poptopologycounts,foralltreestringsort *fa, int *poptreeproposed, int uniformprior );
void init_RF_nodeinfo(void);
/**** GLOBAL FUNCTIONS IN build_gtree.cpp ****/
void makeHKY (int ci, int li, int nosimmigration);
void makeIS (int ci, int li, int nosimmigration);
void makeJOINT_IS_SW (int ci, int li, int nosimmigration);
void makeSW (int ci, int li, int nosimmigration);

/**** GLOBAL FUNCTIONS IN build_poptree.cpp ****/
void testtreewrite (char poptreestring[]);
void add_ghost_to_popstring (char poptreestring[]);
void setup_poptree (int ci, char poptreestring[]);
void reset_poptree (int ci, char poptreestring[]);
void set_poptree_update_record(void);
int getpoptreestringnum(char *s);   /* added with hidden genealogy stuff */
void rewrite (char *substr);

/**** GLOBAL FUNCTIONS IN calc_prob_data.cpp ****/
void labelgtree (int ci, int li, int edge);
double likelihoodHKY (int ci, int li, double mutrate, double kappa, int e1,int e2, int e3, int e4);

void calc_sumlogk (int ci, int li, double *sumlogk);
void free_sumlogk (void);
double likelihoodIS (int ci, int li, double mutrate);
double likelihoodSW (int ci, int li, int ai, double u, double tr);
void checklikelihoodSW (int ci, int li, int ai, double u);

/**** GLOBAL FUNCTIONS IN ginfo.cpp ****/
void init_genealogy_weights (struct genealogy_weights *gweight);
void setzero_genealogy_weights (struct genealogy_weights *gweight);
void free_genealogy_weights (struct genealogy_weights *gweight);
void init_probcalc (struct probcalc *pcalc);
void free_probcalc (struct probcalc *pcalc);
void copy_treeinfo (struct genealogy_weights *dest,
                    struct genealogy_weights *srce);
void copy_probcalc (struct probcalc *dest, struct probcalc *srce);
void sum_subtract_treeinfo (struct genealogy_weights *addup,
                            struct genealogy_weights *addtoplus,
                            struct genealogy_weights *addtominus);
int calc_gsampinf_length (void);
//AS: Added the index of cold chain to this as on Mon Mar 31 16:20:45 EDT 2014
void savegsampinf (float *g, int z);
void savegsampinf_debug_ti_addinfo (float *g, int z, char *a); // used when if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
void sum_treeinfo (struct genealogy_weights *addup,
                   struct genealogy_weights *addto);

/**** GLOBAL FUNCTIONS IN histograms.cpp ****/
void printhistograms (FILE * outfile, long int mcmcrecords,
                      double generationtime,int usegenerationtimedefault, double scaleumeaninput,char priorfilename[]);
void printmigrationhistograms (FILE * outfile, long int mcmcrecords);

/**** GLOBAL FUNCTIONS IN initialize.cpp ****/
void set_iparam_poptreeterms (int ci);// also called from update_poptree
void setup (char infilename[], int *fpstri, char fpstr[], char priorfilename[],char topologypriorinfostring[],int currentid);
/**** GLOBAL FUNCTIONS IN output.cpp ****/
void closeopenout (FILE ** p_to_file, char fname[]);
char* outputbanner(const char *bannertext);
void checkoutfileclosed (FILE ** outfile, char outfilename[]);
void printrunbasics (FILE * outfile, int loadrun, char fpstr[],
                      int burninsteps,int burninsteps_old,int runsteps_old, int mcmcrecords_old,int genealogysamples_old,
                     int recordint, int mcmcrecords,
                     int savegenealogyint,double hilike, double hiprob);
void asciitrend (FILE * outfile, struct value_record *v, int trenddoublepoint,int trendspot);
void asciicurve (FILE * outfile, struct plotpoint *a, char *qlabel,
                 int logscale, int mcmcrecords);
void printacceptancerates (FILE * outto, int numrec,
                           struct chainstate_record_updates_and_values *rec[],
                           const char *printstring);
void callprintacceptancerates (FILE * outto, int currentid);
void printcurrentvals (FILE * outto);
void printcurrent_tvals (FILE * outto,int  currentid);
void savegenealogyfile (char genealogyinfosavefilename[], FILE * genealogyinfosavefile,
                   int *lastgenealogysavedvalue, int gsampinflength);
void preparehistogram (FILE * outfile, int mode, long int mcmcrecords,
                       double scaleumeaninput, double generationtime);
void printmigratehist (FILE * outfile, int mcmcrecords);
void printtmrcahist (FILE * outfile, int mcmcrecords);
void print_means_variances_correlations (FILE * outfile);
void sort_and_print_alltreestrings(FILE * outfile, int *poptopologycounts, int *poptopologyproposedlist_rec, char *topologypriorinfostring); /* for printing population tree posteriors,  sorts and prints */

/**** GLOBAL FUNCTIONS IN readata.cpp ****/
void read_datafile_top_lines (char infilename[], int *fpstri, char fpstr[]);
void readdata (char infilename[], int *fpstri,
               char fpstr[], int **numsitesIS, int currentid);
int imaInfileNpops (const char *fn);
void callasciitrend (FILE * outtofile,int trenddoublepoint,int trendspot);
void callasciicurves (FILE *outfile,int mcmcrecords);
void printsteps (FILE * outto, double like, double probg,int burndone, int burninsteps);
/**** GLOBAL FUNCTIONS IN surface_call_functions.cpp ****/
double margincalc (double x, double yadjust, int pi, int logi);
void marginalopt (int firsttree, int lasttree, double *mlval,double *peakloc);
double margin95 (double mlval[], double peakloc[], int pi, int UL);
void findmarginpeaks (FILE * outfile, float *holdpeakloc);

/**** GLOBAL FUNCTIONS IN surface_search_functions.cpp ****/
void bracket (int ndim, int firsttree, int lasttree, double *axval, double *bxval,
           double *cxval, double *funca, double *funcb, double *funcc,
           double (*func) (int, int, int, double, int), int ifneeded);
double goldenmin (int ndim, int firsttree, int lasttree, double axval, double bxval,
           double cxval, double tol_, double *xmin, double (*func) (int, int,
                                                                 int, double, int), int ifneeded);

/**** GLOBAL FUNCTIONS IN popmig.cpp ****/
double calc_popmig (int thetai, int mi, double x, int prob_or_like);
double calc_pop_expomig (int thetai, int mi, double x, int prob_or_like);
void marginalopt_popmig (int firsttree, int lasttree, double *mlval, double *peakloc, int *mpop, int *mterm);
double marginpopmig (int mi, int firsttree, int lasttree, double x, int thetai);

/**** GLOBAL FUNCTIONS IN gtint.cpp ****/
void print_greater_than_tests (FILE * outfile);

/**** GLOBAL FUNCTIONS IN swapchains.cpp ****/
void setheat (double hval1, double hval2, int heatmode, int currentid);
void freeswapstuff(void);

int setswaptries(void);
int swapchains (int nargs, ...);
void swapchains_bwprocesses(int currentid, int swaptries,int *numattemptwithin,int *numattemptbetween,int *numsuccesswithin,int *numsuccessbetween);
void printchaininfo (FILE * outto, int heatmode,
                     double hval1, double hval2, int currentid);

/**** GLOBAL FUNCTIONS IN treeprint.cpp  (for debugging only) ****/
#ifdef TURNONCHECKS
void gtreeprint (int ci, int li/*, int step , int callsource */ );
void poptreeprint (int ci );
void poptreeprint_frompointer (int ci, struct popedge *poptree/*, int step*/);
void checkgenealogy(int ci, int li, int mode );/* 1 check migs, but not length, 2 check mighg but not  mig, 3 check neither mig or mighg, 0 check everything including lengths*/
void checkpoptree(int ci, int mode); // mode: 0 do all checks, 1 skip period check
void checkprobs(int ci, int locusi);
void printgenealogyweights(int ci,int li);
void checkgenealogyweights(int ci);
void check_hgprob_sums(int ci);
void checkdetailedbalance(double newlikelihood, double oldlikelihood, double newprior, double oldprior, double propose_old_given_new, double propose_new_given_old, double beta);
void checkdetailedbalance_chainswap(double likelihood_i, double likelihood_k, double prior_i, double prior_k, double beta_i, double beta_k);
void compare_genealogy_weights (struct genealogy_weights *gw1,struct genealogy_weights *gw2);
void pcheck(int ci, int fromcode);
#endif


/**** GLOBAL FUNCTIONS IN update_gtree_common.cpp ****/

int findperiod (int ci, double t);
int findtindex (int ci, double t);
int nowedgepop (int ci, struct edge *gtree, double ptime);
void init_treeweight (void);
void free_treeweight (void);
void init_holdgtree (struct genealogy *g, int numgenes);
void free_holdgtree (struct genealogy *g, int numgenes);
int treeweight (int ci, int li, int callcode);

#ifdef TURNONCHECKS
/* for debugging genealogy structure
checks some things that checkgenealogy() does not
This is in update_gtree_common.cpp in order to use sgEvent and related stuff
that is also used by treeweight() */
void checktreeweight (int ci, int li); // for debugging genealogy structure
#endif //TURNONCHECKS

void initialize_integrate_tree_prob (int ci,
                                     struct genealogy_weights *gweight,
                                     struct probcalc *pcalc);
void integrate_tree_prob (int ci, struct genealogy_weights *gweight,
                          struct genealogy_weights *holdgweight,
                          struct probcalc *pcalc, struct probcalc *holdpcalc);
void copyfraclike (int ci, int li);
void storescalefactors (int ci, int li);
void restorescalefactors (int ci, int li);
double finishSWupdateA (int ci, int li, int ai, int edge, int downedge,
                        int sisedge, int newsisedge, double u, double *Atermnum, double *Atermdenom);
double oldfinishSWupdateA (int ci, int li, int ai, int edge,
                 int downedge, int sisedge, int newsisedge,
                 double u, double *Aterm);
double finishSWupdateAD (int ci, int li, int ai, int edge, int downedge,
                         int sisedge, int newsisedge, double u,
                         double *Aterm);
#if 0
//  CR 110825.1 updateA() not used and no longer compiled into code
double updateA (int ci, int li, int ai, double u, int *count);
#endif   // if 0



/**** GLOBAL FUNCTIONS IN update_gtree.cpp ****/

void init_updategenealogy (void);
void free_updategenealogy (void);
int updategenealogy (int ci, int li, int *topolchange, int *tmrcachange);


/******* global in update_t_NW **********/
void init_t_NW (void);
void free_t_NW (void);
int changet_NW (int ci, int timeperiod);

/******* global in update_t_RY **********/
void init_t_RY (void);
void free_t_RY (void);

int changet_RY1 (int ci, int p);

/******* global in update_t_RYhg **********/
void init_t_RYhg (void);
void free_t_RYhg (void);
int changet_RYhg (int ci, int p);

/**** GLOBAL FUNCTIONS IN update_mc_params.cpp ****/

int changeu (int ci, int j, int *k);
int changekappa (int ci);

/**** GLOBAL FUNCTIONS IN utilities.cpp ****/
int imaDirBase (const char *a, char **b);
char *shorten_e_num (char *s);
int skip_a_line (FILE *fp);
int read_int (FILE *fp);
/* Print error located at a source file.
  use AT macro to print out file name and line number
 eg  IM_errloc (AT, "Error %d", d);
 */
void IM_errloc (const char *loc, const char *fmt, ...);
/* used with error messages  in   simerrmsg[] */
void IM_err (int i, const char *fmt, ...);
void errr (int ci, int li, int i, double val1, double val2);
void nrerror (const char error_text[]);
double mylogcosh (double x);
double mylogsinh (double x);
void setseeds (int seed);
void resetseeds (int seed);     /* use for debugging, provides control over sequence of random numbers */
void unsetseeds ();
double uniform ();
double uniforminterval(double lower,double upper);
double expo (double c);
double normprob (double mean, double stdev, double val);
double normdev (double mean, double stdev);
int poisson (double param, int condition,int cmm);
int geometric (double p);
void hpsortmig (struct migstruct *lptr, int n);
void shellpnodetimes (struct pnodetime *hptr, int length); /* added for updating population tree using hidden genealogies hgstuff */
void shellhist (struct hlists *hptr, int length);
void indexx (unsigned long n, struct gtreeevent *arr, unsigned long *indx);
void gcf (double *gammcf, double a, double x, double *gln);
void gser (double *gamser, int a, double x, double *gln);
void checktreestring (char *t);
int put_spaces_in_filepaths(char *pathstr);
void convertToUpperCase(char *sPtr);
/* double expint(int n, double x); */
double expint (int n, double x, int *islog);
double uppergamma (int a, double x);
double lowergamma (int a, double x);
int isemptystring(char *s);
void strdelete (char *s, int pos, int len);
void strinsert (char *dest, char *source, int pos);
void strtrunc (char *s, char c);
char *nextwhite (char *c);
int allwhitespace (char *c);
char *nextnonspaceafterspace (char *textline);
char *nextnonspace (char *textline);
void setlogfact (void);
void ieevent (struct eevent *a);
double bessi (int n, double x);
int **alloc2Dint (int rows, int cols);
double **orig2d_alloc2Ddouble (int rows, int cols);
void orig2d_free2D (void **a, int rows);
extern void eexp (double x, double *m, int *z);
extern void copymig (struct migstruct *m1, struct migstruct *m2);
extern int bitran (void);
extern int randposint (int lessthanval);
double **alt2d_alloc2Ddouble(long m, long n);
int **alt2d_alloc2Dint(long m, long n);
double **allocAndInit2Ddouble (int rows, int cols); /* cr 110907.1 */
void alt2d_free2D(double **x/*, long m , long n */);
void alt2d_free2Dint(int **x/*, long m , long n */);
int **alloc2Dint_contiguous(int **arrdata, int rows, int columns);
void free_2Dint_contiguous(int **arr, int *arrdata);
void setupdatescalarinfo(struct updatescalarinfo *a, double sav, double saav, double savmin, double savmax, double target, int numattemptscheck);
void checkupdatescalarer(struct updatescalarinfo *a);
void resetupdatescalarer(struct updatescalarinfo *a);
int  file_exists(const char * filename);

struct dictionary_node_kr *dictionary_install(char *name, double val,struct dictionary_node_kr **hashtab);  // returns an nlist pointer
double getvalue(char *key, struct dictionary_node_kr **hashtab);
int metropolishastingsdecide(double logmhratio,int othercriteria);
char* timestring(time_t seconds);

/* GLOBAL Functions in multi_t_bins */
void setup_multi_t_arrays (int z);
void free_multi_t_arrays (void);
void return_joint_t (double tvals[]);
double joint_t_prob (double *tvals);


/* freemem.cpp */
void freeanymemory (void);

/**** Global Functions in File: mcmcfile */
void writemcf (char mcffilename[],char commandline[],int mcmcrecords,int mcmcrecords_old,int genealogysamples_old,int burninsteps_old,int runsteps_old, double hilike,double hiprob,int currentid);
void readmcf (char mcffilename[],int *mcmcrecords,double *hilike,double *hiprob, int currentid);
#ifdef TURNONCHECKS
void readima2mcf (char ima2mcffilename[]);
#endif

/**** Global Functions in File: jointfind */
void findjointpeaks(FILE **outfile,char *outfilename, char *nestfname,int number_of_parameters);

/**** Global Functions in File: readpriorfile*/
void readpriorfile(char priorfilename[],double *popsizepriorvals, double **mpriorvals);


/* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
/**** Global Functions in File: marglike.cpp */
void initmarginlikecalc();
void freemarginlikecalc();
void move_calcmarglike_vals(void);
void summarginlikecalc();
void thermo_marginlike_calc(int n, double *estimator);
// 8/22/2016  shut down most marginal likelihood stuff because it was not working well with non-constant beta intervals
//void stepstone_get_Lmax();
//double harmonicmarginlikecalc(int k);
//double thermomarginlikecalc(int k);
//void steppingstone_marginlike_calc2(int n,double *estimator2, double *stdev2);


/****** Global Functions in File update_migration_priors.cpp */

int update_migration_prior_intree(int ci, int mi);
void update_migration_prior_not_intree(int ci, int *attempted, int *accepted);
int update_popsize_prior_intree(int ci, int qi);
void update_popsize_prior_not_intree(int ci, int *attempted, int *accepted);
void init_migration_prior_update();
void free_migration_prior_update();
void init_hyperprior_arrays(int ci);
void free_hyperprior_arrays(int ci);

/**** Global Functions in chainprint.cpp (for debugging only) ******/
#ifdef TURNONCHECKS
void chaininfo_print(int currentid,int recordint);
#endif //TURNONCHECKS

/***** Global Functions in burnphylogenycounts.cpp ********/
void init_burn_phylogeny_counts();
void free_burn_phylogeny_counts();
void recordburntopology(void);
void outputburntopologycounts(FILE *burnfile, int currentid);

/***** Global Functions in gewekeIM.cpp ********/
double twosidedp(double z);
double gewekez(double *vars1, double *vars2, int lenvars1, int lenvars2);
void printgewekez(FILE * outto);

/**** Global Functions in File: msdirent.cpp */
// nothing needed here,  prototypes are in msdirent.hpp


//gextern int badmigrationupdate;   used sometimes for debugging
