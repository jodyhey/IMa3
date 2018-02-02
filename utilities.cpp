/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

#undef GLOBVARS
#include "ima.hpp"
//#include "utilities.hpp" moved stuff out of here and deleted it

/*********** LOCAL STUFF **********/
//double loghalffact[50 * ABSMIGMAX];   //leftover from assignment stuff, required xtrabbits
//int BITNUMBERTRUE[256];
unsigned long *iseed;

/* local function prototypes */

char *strdup(char *s); /* make a duplicate of s */
struct dictionary_node_kr *dictionarylookup(char *key, struct dictionary_node_kr **hashtab);
unsigned hashkr(char *s);

static int *ivector (long nl, long nh);
static void free_ivector (int *v, long nl/* , long nh*/);
static double bessi0 (double x);
static double bessi1 (double x);

/* these could be deleted probably, as they only turn up in this file */
int comma_exists (char *s);  // only in utilities.cpp
double logsum (int n, ...);  // only in utilities.cpp
double logsuma (int n, double *d);   // only in utilities.cpp
double read_double (FILE *fp);  // only in utilities.cpp

/* try a new random number generator in rand.c */
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);

#ifdef RANDOM_NUMBERS_FROM_FILE
/* used in testing only */
unsigned long getKnownRandom(void);
#endif /* RANDOM_NUMBERS_FROM_FILE */

// stuff called by indexx() - sorting an index,  from NR 
//#define NRANSI
#define NR_END 1l
#define FREE_ARG char*
int *
ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  v = (int *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (int)));
  if (!v)
  {
    nrerror ("allocation failure in ivector()");
  }
  return v - nl + NR_END;
}

void
free_ivector (int *v, long nl/*, long nh*/)
/* free an int vector allocated with ivector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
} 

double
bessi0 (double x)
{
  double ax, ans;
  double y;
  if ((ax = fabs (x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
      + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2))))); 
  }
  else
  {
    y = 3.75 / ax;
    ans = (exp (ax) / sqrt (ax)) * (0.39894228 + y * (0.1328592e-1
      + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
  }
  return ans;
}

double
bessi1 (double x)
{
  double ax, ans;
  double y;
  if ((ax = fabs (x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
      + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
  }
  else
  {
    y = 3.75 / ax;
    ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
    ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 +
      y * (-0.1031555e-1 + y * ans))));
    ans *= (exp (ax) / sqrt (ax));
  }
  return x < 0.0 ? -ans : ans;
}

/********** GLOBAL FUNCTIONS ***********/

/* Here is how to add an error message.
 * 1. Add an element to enum in utilities.h with prefix of "IMERR_XXX"
 * 2. Add the corresponding error message string to simerrmsg.
 * 3. Call function IM_err (IMERR_XXX, ...) where ... is like the way that
 *    printf function arguments are used.
 */
static const char *simerrmsg[] = {
  /* 0 */ "no error",
  /* 1 */ "cannot open a file for reading",
  /* 2 */ " memory error ",
  /* 3 */ "problem finding .ti file(s)",
  /* 4 */ "cannot create file",
  /* 5 */ "cannot open file for appending",
  /* 6 */ " problem with opening or closing output file",
  /* 7 */ " problem reading file",
  /* 8 */ "incompatibility on command line",
  /* 9 */ "not enough information provided on command line",
  /* 10 */ "problem with command line formatting ",
  /* 11 */ "problem with heating terms in command line",
  /* 12 */ "missing population string in input file",
  /* 13 */ "some problem in the population tree string ",
  /* 14 */ "option number on command line not allowed",
  /* 15 */ "problem w/ number of loci indicated in data file",
  /* 16 */ "problem in specifying nested models for LLR tests",
  /* 17 */ "",
  /* 18 */ "mutation range priors to constraining - not able to set starting values",
  /* 19 */ "product of mutation scalars not equal to 1",
  /* 20 */ "problem with mutation rate scalars",
  /* 21 */ "Too many migration events found when storing edges",
  /* 22 */ "Too much migration called for in checkmig - reduce maximum value of migration parameter",
  /* 23 */ "problem in calculating HPD interval",
  /* 24 */ "problem with hyperpriors",
  /* 25 */ "error reading data, too many lines or line too long",
  /* 26 */ "error in data",
  /* 27 */ "error in locus information",
  /* 28 */ "problem reading mcf file",
  /* 29 */ "cannot load mcf file, split times in file not compatiable with t prior",
  /* 30 */ "problem writing mcf file",
  /* 31 */ "",
  /* 32 */ "",
  /* 33 */ "problem with root ",
  /* 34 */ "",
  /* 35 */ "input file invalid ",
  /* 36 */ "data not compatible with infinite sites model",
  /* 37 */ "",
  /* 38 */ "likelihoods do not add up for stepwise model",
  /* 39 */ "",
  /* 40 */ "",
  /* 41 */ "error in lowergamma function",
  /* 42 */ "error in uppergamma function",
  /* 43 */ "error using LogDiff macro,  difference is non-positive",
  /* 44 */ "",
  /* 45 */ "error calculating prior of t  in multi_t_prior_func",
  /* 46 */ "problem in the values specified in a file with parameter priors",
  /* 47 */ "model calls for 0 migration parameters",
  /* 48 */ "",
  /* 49 */ "",
  /* 50 */ "",
  /* 51 */ "error in NR functions",
  /* 52 */ "GSL error",
  /* 53 */ "",
  /* 54 */ "",
  /* 55 */ "Potential bug",
  /* 56 */ "Invalid input file: gene name",
  /* 57 */ "Assignment Error",
  /* 58 */  "",
  /* 59 */  "",
  /* 60 */  "problem calculating migration path probability when updating genealogy",
  /* 61 */  "",
  /* 62 */  "some problem with prior, or likelihood",
  /* 63 */  "",
  /* 64 */  "",
  /* 65 */  "",
  /* 66 */  "",
  /* 67 */  "",
  /* 68 */  "",
  /* 69 */  "",
  /* 70 */   "problem in the number of chains per CPU",
  /* 71 */  "",
  /* 72 */ "miscellaneous error"
};

void
IM_err (int i, const char *fmt, ...)
{
  va_list args;
  va_start (args, fmt);
  fprintf (stderr, "IMa3: %s - ", simerrmsg[i]);
  vfprintf (stderr, fmt, args);
  fprintf (stderr, "\n");
  va_end (args);
  exit (i);
}

/* Print error located at a source file. 
  use AT macro to print out file name and line number 
 eg  IM_errloc (AT, "Error %d", d);
 */
void
IM_errloc (const char *loc, const char *fmt, ...)   // don't use this much.  mostly use IM_err
{
  va_list args;
  va_start (args, fmt);
  /* fprintf (stderr, "ima at %s: %s - ", loc, simerrmsg[i]); */
  fprintf (stderr, "ima at %s - ", loc);
  vfprintf (stderr, fmt, args);
  fprintf (stderr, "\n");
  va_end (args);
  exit (1);
}

void
nrerror (const char error_text[])
/* Numerical Recipes standard error handler */
{
  IM_err (IMERR_NUMERICALRECIPES, " error text: %s", error_text);
}

/* for a moderately large value x,  cosh[x] = sinh[x] 
also (cosh[x] + sinh[x] = exp[x]  so for a large value of x  cosh[x] = sinh[x] = Exp[x]/2
so rather than return floating error,  for the log of a cosh of a large number  just return x - log[2]  
used to have the value at 100,  but changed to 20, for which they are identical to 17 digits
*/

#define SINHCOSHXBIG 20
double
mylogcosh (double x)
{
  if (x < SINHCOSHXBIG)
    return log (cosh (x));
  else
    return x - LOG2;
}

double
mylogsinh (double x)
{
  if (x < SINHCOSHXBIG)
    return log (sinh (x));
  else
    return x - LOG2;
}
#undef SINHCOSHXBIG

/***********************************/
/* RANDOM NUMBER RELATED FUNCTIONS */
/***********************************/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

#ifdef RANDOM_NUMBERS_FROM_FILE
 
/* for testing only.  read a random number from a file
 * containing a known list of random numbers (in ascii)
 * then convert the ascii to unsigned long and return.
 */ 
#define ranNumFile "randNums"
#define MAXLINE 120

unsigned long getKnownRandom(void)
{
  static FILE *ranNum = NULL;
  static int firstTime = 1;

  unsigned long y;
  char line[MAXLINE];
  char comment = '#';

  while ( 1 ) 
  {
  	if ( firstTime )
  	{
	  	firstTime = 0;
    if ((ranNum = fopen(ranNumFile, "r")) == NULL )
    {
      IM_err (IMERR_CREATEFILEFAIL, "Error opening randNums file");
    }
	  	//if (( ranNum = fopen(ranNumFile, "r")) == NULL )
			//assert (0);
	  	if (( fgets(line, MAXLINE,  ranNum)) == NULL)
		  	assert (0);
	  	/* skip the beginning comment lines started with # */
	  	while( line[0] == comment )
	  	{
		  	fgetval = fgets(line, MAXLINE,  ranNum);
	  	}
		break;
  	}
  	else 
	  	if ( fgets(line, MAXLINE,  ranNum) == NULL )
	  	{
		  	firstTime = 1;
		  	fclose(ranNum);
	  	}
		else
			break;
  	}  /* end while */
  y = strtoul(line, NULL, 0);
  return y;
}

#endif /* RANDOM_NUMBERS_FROM_FILE */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

// hook for sanity testing with a known list of random numbers - janeen
#ifdef RANDOM_NUMBERS_FROM_FILE
    y = getKnownRandom();
#endif /* RANDOM_NUMBERS_FROM_FILE */
    return y;
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    /*  CR 110715.2
     *  converted 1/2^32 to a const instead of performing a calculation 
     *  each time function is called. 
     *       n * 1/2^32 = n * 2.32830643653869628906e-10
     */
    return (((double)genrand_int32()) + 0.5) * 2.32830643653869628906e-10;
    /* divided by 2^32 */
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

void
setseeds (int seed)
{
  init_genrand((unsigned long) seed);
  /* try a new random number generator 9_19_10  don;t use z_rndu  */
 // z_rndu = 170 * (seed % 178) + 137;
  iseed = static_cast<unsigned long *> (malloc (sizeof (unsigned long)));
  *iseed = ULONG_MAX / 2 + (unsigned long) seed;        // just set it so that iseed points to a large number - probably not necessary
  
}

void
unsetseeds ()
{
  XFREE (iseed);
  iseed = NULL;
  return;
}

double
uniform ()
{
    double r;
 /* started using a new random number generator 9_19_10 */
 /* see rand.c */
  r = genrand_real3();
  return r;
}

double uniforminterval(double lower,double upper)
{
  return lower + uniform() * (upper - lower);
}

/* c is the parameter,  i.e. 1/c  is the mean */ 
double
expo (double c)
{
  return -log (uniform ()) / c;
}
int
randposint (int lessthanval) // int, range is 0 thru (lessthanval-1) inclusive (i.e. can't get lessthanval). actually should be rand_nonneg_int because 0 is allowe
{
  return (int) floor (uniform () * lessthanval);
}

/* generate a sorted list of random uniform values between a and b
- generate n+1 independent exponentially-distributed variates (parameter does not matter so use 1). 
- Normalize their cumulative sum to the range (0,1) by dividing by the sum. 
- Drop the largest value (which will always be 1). Rescale to the range (a,b).
*/
void randomsorted(double a, double b, int n, double r[])
{
  int i;
  double sum = 0.0;
  for (i=0;i<n;i++)
  {
    r[i] = (-log(uniform ())) + sum;
    sum += r[i];
  }
  sum +=  expo(1.0);
  for (i=0;i<n ;i++)
  {
	r[i] = (r[i] * ((b-a)/sum)) + a;
    //r[i]/= sum;
    //r[i] *= (b-a);
    //r[i] += a;

  }
}


/* for binary random numbers  -	quick - based on Press et al irbit2() */
#define IB1 1
#define IB18 131072
#define MASK 19
int
bitran (void)
{
#ifdef RANDOM_NUMBERS_FROM_FILE
  if (uniform() < 0.5)
    return 0;
  else
    return 1;
#endif
  if (*iseed & IB18)
  {
    *iseed = ((*iseed ^ MASK) << 1) | IB1;
    return 1;
  }
  else
  {
    *iseed <<= 1;
    return 0;
  }
}

#undef MASK
#undef IB18
#undef IB1

#define ONEOVERSQRT2PI  0.3989422803
double
normprob (double mean, double stdev, double val)
{
  return ONEOVERSQRT2PI * exp (-SQR ((val - mean) / stdev) / 2) / stdev;
}

#undef ONEOVERSQRT2PI
double
normdev (double mean, double stdev)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  double rescale;
  if (iset == 0)
  {
    do
    {
      v1 = 2.0 * uniform () - 1.0;
      v2 = 2.0 * uniform () - 1.0;
      rsq = v1 * v1 + v2 * v2;
    }
    while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt (-2.0 * log (rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    rescale = ((v2 * fac * stdev) + mean);
    return rescale;
  }
  else
  {
    iset = 0;
    rescale = ((gset * stdev) + mean);
    return rescale;
  }
}                               /* normdev */

/* gets a poisson random variable.  checked it out with simulations 
if condition is 0, pick an even number
if condition is 1, pick an odd number
if condition is 2,  pick any number other than 0 (because they must be different)
if condition is 3, pick any number other than 1
if condition is -1 pick any number 
If param > minpp the normal distribution is used as an approximation
Also stuck in some stuff to deal with the case when condition is 1 and param << 1

quickly returns -1 if the current migration count cmm  is too big 

If the new poisson variable causes the total  (cmm + i)  to be too big 
If cmm + i >= MIGMAX  return the largest value of i consistent with  cmm + i <= MIGMAX and condition

*/

#define USENORMAL  100.0  // only for very large parameter values, because it is not a perfect approximation
#define MINPP 0.25
int
poisson (double param, int condition,int cmm)
{
  long double randnum, raised, rcheck;
  int i;
  int stop;

  switch (condition) // check incoming cmm values to see if they allow room
  {
  case -1:
    if (cmm > MIGMAX)
      return -1; // update must be rejected
    break;
  case 0:
    if (cmm > MIGMAX)
      return -1; // update must be rejected
    break;
  case 1:
    if (cmm > MIGMAX-1)
      return -1; // update must be rejected
    break;
  case 2:
    if (cmm > MIGMAX-1)
      return -1; // update must be rejected
    break;                    //*
  case 3:
    if (cmm > MIGMAX)
      return -1; // update must be rejected
    break;
  }
  if (param < MINPP && condition == 1)  /* i.e. need an odd number but the parameter is a small value */
  {
    randnum = uniform ();
    rcheck = param / sinh (param);
    if (randnum < rcheck)       /* param / sinh(param) is prob of # being 1, given it must be odd */
      i = 1;
    else
    {
      rcheck = param * (6 + param * param) / (6 * sinh (param));
      if (randnum < rcheck)
        i = 3;
      else
        i = 5;
    }
    return (i);
  }
  /*need any number not zero,  but the parameter is a small value  */
  /* seems to work - check  file Poisson_low_value_parameter_simulation_check.nb */
  if (param < MINPP && condition == 2)
  {
    raised = exp (-param);
    rcheck = raised = param * raised / (1 - raised);
    randnum = uniform ();
    i = 1;
    while (randnum > rcheck)
    {
      raised *= param / (i + 1);
      rcheck += raised;
      i++;
      //trap cases where too large a value is being simulated
      if (cmm + i >= MIGMAX)
        break;
    }
    return (i);
  }
  /* can only get here if either param > MINPP or  param < MINPP but condition not in [1,2] */
  do
  {
    if (param >= USENORMAL)     // used normal approximation
    {
      i = IMAX (0, POSROUND (normdev (param, param)));
      /* I checked this out extensively, and although the mean and variance of the normal and the poisson come to be 
         very close for parameter values of 10 or or,  non-trivial differences persist for the tails of the distribution
         even for very large parameter values*/
    }
    else
    {
      raised = exp (-param);
      randnum = uniform ();
      for (i = 0; randnum > raised; i++)
      {
        randnum *= uniform ();
        //trap cases where too large a value is being simulated
        if (cmm + i >= MIGMAX)
        {
          // MIGMAX > 1  so conditions -1,2 and 3 are satisified 
          // have to make sure conditions 0 and 1 are satisfied
          if ( (condition == 1 && !ODD(i)) ||(condition == 0 && ODD(i))) // i is at the max,  but must make it fit with condition
              i -= 1;
          return (i);  // go ahead and return at this point
          break;
        }
      }
    }
    switch (condition)
    {
    case -1:
      stop = 1;
      break;
    case 0:
      stop = !ODD (i);
      break;
    case 1:
      stop = ODD (i);
      break;
    case 2:
      stop = (i != 0);
      break;                    //*
    case 3:
      stop = (i != 1);
      break;
    }
    stop = stop && (cmm + i) <= MIGMAX;
  } while (!stop); 
  return (i);
}                               /* pickpoisson */

#undef USENORMAL
#undef MINPP
int
geometric (double p)
/* returns a geometrically distributed random integer variable >= 0 */
/* this distribution is given by prob(k) = p*(1-p)^(k-1)  where k = 1,2,...  which has a mode of 1 and an expectation of 1/p */
/* it is a bit different from prob(k) = p*(1-p)^k   where k  = 0,1, 2... , which has a mode of 0 and an expectation of (1-p)/p */
/* the variance of these two different versions is the same, i.e. (1-p)/p^2 */
/* checked this in various ways - seems ok */
{
  return (int) ceil (log (uniform ()) / log (1.0 - p));
}

/***********************************/
/* SORTING FUNCTIONS */
/***********************************/

void
hpsortmig (struct migstruct *lptr, int n)
{
  unsigned long i, ir, j, l;
  double t;
  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
    {
      t = (lptr + --l)->mt;
    }
    else
    {
      t = (lptr + ir)->mt;
      (lptr + ir)->mt = (lptr + 1)->mt;
      if (--ir == 1)

      {
        (lptr + 1)->mt = t;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir)
    {
      if (j < ir && (lptr + j)->mt < (lptr + (j + 1))->mt)
        j++;
      if (t < (lptr + j)->mt)
      {
        (lptr + i)->mt = (lptr + j)->mt;
        i = j;
        j <<= 1;
      }
      else
      {
        j = ir + 1;
      }
    }
    (lptr + i)->mt = t;
  }
}                               /*hpsortmig */

/* added for updating population tree using hidden genealogies hgstuff */
void
shellpnodetimes (struct pnodetime *hptr, int length)
{
  double aln = 1.442695022, tiny = 1.0e-5;
  struct pnodetime h;
  static int nn, m, lognb2, i, j, k, l;
  lognb2 = (int) floor (log ((double) length) * aln + tiny);
  m = length;
  for (nn = 1; nn <= lognb2; nn++)
  {
    m = m / 2;
    k = length - m;
    for (j = 0; j <= k - 1; j++)
    {
      i = j;
    reloop:l = i + m;
      if ((hptr + l)->time < (hptr + i)->time)
      {
        h = *(hptr + i);
        *(hptr + i) = *(hptr + l);
        *(hptr + l) = h;
        i = i - m;
        if (i >= 0)
          goto reloop;
      }
    }
  }
}                               /* shellpnodetimes */

void
shellhist (struct hlists *hptr, int length)
{
  double aln = 1.442695022, tiny = 1.0e-5;
  struct hlists h;
  static int nn, m, lognb2, i, j, k, l;
  lognb2 = (int) floor (log ((double) length) * aln + tiny);
  m = length;
  for (nn = 1; nn <= lognb2; nn++)
  {
    m = m / 2;
    k = length - m;
    for (j = 0; j <= k - 1; j++)
    {
      i = j;
    reloop:l = i + m;
      if (((hptr + l)->p < (hptr + i)->p)
          || (((hptr + l)->p == (hptr + i)->p)
              && (((hptr + l)->v < (hptr + i)->v))))
      {
        h = *(hptr + i);
        *(hptr + i) = *(hptr + l);
        *(hptr + l) = h;
        i = i - m;
        if (i >= 0)
          goto reloop;
      }
    }
  }
}                               /* shellhist */

/* quicksort of an index of locations */
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp
#define M 7
#define NSTACK 50
void
indexx (unsigned long n, struct gtreeevent *arr, unsigned long *indx)
{
  unsigned long i, indxt, ir = n, itemp, j, k, l = 1;
  int jstack = 0, *istack;
  double a;
  istack = ivector (1, NSTACK);
  for (j = 1; j <= n; j++)
    indx[j] = j;
  for (;;)
  {
    if (ir - l < M)
    {
      for (j = l + 1; j <= ir; j++)
      {
        indxt = indx[j];
        a = arr[indxt].time;
        for (i = j - 1; i >= l; i--)
        {
          if (arr[indx[i]].time <= a)
            break;
          indx[i + 1] = indx[i];
        }
        indx[i + 1] = indxt;
      }
      if (jstack == 0)
        break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else
    {
      k = (l + ir) >> 1;
      SWAP (indx[k], indx[l + 1]);
      if (arr[indx[l]].time > arr[indx[ir]].time)
      {
        SWAP (indx[l], indx[ir]);
      }
      if (arr[indx[l + 1]].time > arr[indx[ir]].time)
      {
        SWAP (indx[l + 1], indx[ir]);
      }
      if (arr[indx[l]].time > arr[indx[l + 1]].time)
      {
        SWAP (indx[l], indx[l + 1]);
      }
      i = l + 1;
      j = ir;
      indxt = indx[l + 1];
      a = arr[indxt].time;
      for (;;)
      {
        do
          i++;
        while (arr[indx[i]].time < a);

        do
          j--;
        while (arr[indx[j]].time > a);

        if (j < i)
          break;
        SWAP (indx[i], indx[j]);
      }
      indx[l + 1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK)
        nrerror ("NSTACK too small in indexx.");
      if (ir - i + 1 >= j - l)
      {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else
      {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
  free_ivector (istack, 1/*, NSTACK*/);
}

#undef SWAP

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

#undef SWAP
#undef M
#undef NSTACK
#undef SWAP

/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define ITMAX 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void
gcf (double *gammcf, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = logfact[(int) a - 1];
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++)
  {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs (d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs (c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    IM_err (IMERR_LOWERGAMMA, " too many iterations within gcf()");
  *gammcf = exp (-x + a * log (x) - (*gln)) * h;
}

// gcflog is the same as gcf but returns the logarithm of gammcf
// saves a little time and stops some underflows
void
gcflog (double *gammcflog, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = logfact[(int) a - 1];
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++)
  {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs (d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs (c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs (del - 1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    IM_err (IMERR_UPPERGAMMA, " too many iterations within gcflog()");

  //*gammcf=exp(-x+a*log(x)-(*gln))*h;
  *gammcflog = (-x + a * log (x) - (*gln)) + log (h);
}                               //gcflog

#undef FPMIN

/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
void
gser (double *gamser, int a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln = logfact[a - 1];
  if (x <= 0.0)
  {
    if (x < 0.0)
      IM_err (IMERR_UPPERGAMMA,
              " gser() called with negative x value, x: %lf ", x);
    *gamser = 0.0;
    return;
  }
  else
  {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs (del) < fabs (sum) * EPS)
      {
        *gamser = sum * exp (-x + a * log (x) - (*gln));
        return;
      }
    }
    IM_err (IMERR_UPPERGAMMA, " too many iterations within gser() ");
    return;
  }
}

// gserlog is the same as gser but returns the logarith of gamser
// saves a little time and stops some underflows
void
gserlog (double *gamserlog, int a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln = logfact[a - 1];
  if (x <= 0.0)
  {
    if (x < 0.0)
      IM_err (IMERR_LOWERGAMMA,
              " function gserlog() called from lowergamma() with negative x value, x: %lf",
              x);
    *gamserlog = 0.0;
    return;
  }
  else
  {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++)
    {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs (del) < fabs (sum) * EPS)
      {
        //*gamser=sum*exp(-x+a*log(x)-(*gln));
        *gamserlog = log (sum) + (-x + a * log (x) - (*gln));
        return;
      }
    }
    IM_err (IMERR_LOWERGAMMA, " too many iterations within gserlog() ");
    return;
  }
}

#undef ITMAX
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
double
expint (int n, double x, int *islog)
{
  int i, ii, nm1;
  double a, b, c, d, del, fact, h, psi, ans;
  *islog = 0;
  nm1 = n - 1;
  if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
  {
    IM_err (IMERR_UPPERGAMMA,
            " expint() called from uppergamma() with bad value(s). n %d  x %lf",
            n, x);
  }
  else
  {
    if (n == 0)
    {
      ans = exp (-x) / x;
    }
    else
    {
      if (x == 0.0)
      {
        ans = 1.0 / nm1;
      }
      else
      {
        if (x > 1.0)
        {
          b = x + n;
          c = 1.0 / FPMIN;
          d = 1.0 / b;
          h = d;
          for (i = 1; i <= MAXIT; i++)
          {
            a = -i * (nm1 + i);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            del = c * d;
            h *= del;
            if (fabs (del - 1.0) < EPS)
            {
              *islog = 1;

              //ans = h*exp(-x);
              ans = log (h) - x;
              return ans;
            }
          }

          IM_err (IMERR_UPPERGAMMA,
                  " too many iterations in expint() called from uppergamma()");
        }
        else
        {
          ans = (nm1 != 0 ? 1.0 / nm1 : -log (x) - EULER);
          fact = 1.0;
          for (i = 1; i <= MAXIT; i++)
          {
            fact *= -x / i;
            if (i != nm1)
            {
              del = -fact / (i - nm1);
            }
            else
            {
              psi = -EULER;
              for (ii = 1; ii <= nm1; ii++)
                psi += 1.0 / ii;
              del = fact * (-log (x) + psi);
            }
            ans += del;
            if (fabs (del) < fabs (ans) * EPS)
              return ans;
          }
          IM_err (IMERR_UPPERGAMMA,
                  " too many iterations in expint() called from uppergamma()");
        }
      }
    }
  }
  return ans;
}                               //expint
/*
incomplete gamma function confusion

According to Boost:
There are four incomplete gamma functions: two are normalised versions (also known as regularized incomplete gamma functions) that return values in the range [0, 1], and two are non-normalised and return values in the range [0, Γ(a)]. Users interested in statistical applications should use the normalised versions (gamma_p and gamma_q).
All of these functions require a > 0 and x >= 0, 
*/


#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */
double
uppergamma (int a, double x)
/*  Returns the log of what Mathematica calls the incomplete gamma function Gamma[a, x]
This is a unnormalized function.   Modified from numerical recipes gammaq() in two ways:
1.return the unnormalized function (gammaq() returns normalized version)
2.return the logarithm */
{
  int logindicator;
  double gamser, gammcf, gln, p;
  double temp;
  if (x < 0.0 || a < 0.0)
    IM_err (IMERR_UPPERGAMMA, "  step %d, uppergamma arguments: a  %d, x %lf",
            step, a, x);
  if (a == 0)
  {
    /* JH 1/9/2018
      rarely get a bug here  when x = 0.0 
      would like to trap this without triggering error expint()
      but what should temp be in this case?? */ 
    temp = expint (1, x, &logindicator);
    if (!logindicator)
      p = log (temp);

    else
      p = temp;
  }
  else
  {
    if (x < (a + 1.0))
    {
      gser (&gamser, a, x, &gln);
      p = gln + log (1.0 - gamser);
    }
    else
    {
      //gcf(&gammcf,a,x,&gln);
      //p = gln + log(gammcf);
      gcflog (&gammcf, (double) a, x, &gln);
      p = gln + gammcf;
    }
  }
  if (p < -1e200)
    p = -1e200;
  return p;
}                               //uppergamma

double
lowergamma (int a, double x)
/*  Returns the log of  what mathematica would call the complement of 
 the incomplete gamma function.    It is the log of what mathematica implements as Gamma[a, 0, z] 
This is a unnormalized function.   Modified from numerical recipes gammap() in two ways:
1.return the unnormalized function (gammap() returns normalized version)
2.return the logarithm */
{
  double gamser, gammcf, gln, p;
  if (x < 0.0 || a < 0.0)
    IM_err (IMERR_LOWERGAMMA, "  step %d, lowergamma arguments: a  %d, x %lf",
            step, a, x);
  if (a == 0)
  {
    IM_err (IMERR_LOWERGAMMA, "  step %d, lowergamma arguments: a  %d, x %lf",
            step, a, x);
  }
  if (x < (a + 1.0))
  {

    //gser(&gamser,a,x,&gln);
    //p = gln + log(gamser);
    gserlog (&gamser, a, x, &gln);
    p = gln + gamser;
  }
  else
  {
    gcf (&gammcf, (double) a, x, &gln);
    p = gln + log (1 - gammcf);
  }
  if (p < -1e200)
    p = -1e200;
  return p;
}

/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

/***********************************/
/* TEXT, CHARACTER RELATED FUNCTIONS */
/***********************************/

int isemptystring(char *s)
{
  return s[0] == '\0';
}

/* Delete the substring of length "len" at index "pos" from "s". Delete less if out-of-range. */
/* note the counting for pos starts at 1,  not 0,  so might have to call with pos being 1 higher than you think */ 
void
strdelete (char *s, int pos, int len)
{
  int slen;
  if (--pos < 0)
    return;
  slen = (int) (strlen (s) - pos);
  if (slen <= 0)
    return;
  s += pos;
  if (slen <= len)
  {
    *s = 0;
    return;
  }
  while ((*s = s[len]))
    s++;
}

/* insert source into dest at pos  */
void
strinsert (char *dest, char *source, int pos)
{
  char temp[POPTREESTRINGLENGTHMAX];
  temp[0] = '\0';
  if (pos > 0)
    strncpy (temp, dest, (size_t) pos);
  temp[pos] = '\0';
  strcat (temp, source);
  if ((int) strlen (dest) > pos)
    strcat (temp, &dest[pos - 1]);
  strcpy (dest, temp);
}

/* truncate a string at the first instance of c */
void
strtrunc (char *s, char c)
{
  char *cspot;
  if (strchr (s, c))
  {
    cspot = strchr (s, c);
    *cspot = '\0';
  }
}

/* return -1 if empty string,  0 if not all whitespace,  1 if all whitespace */
int allwhitespace (char *c)
{
  size_t i;
  int iswhite = 1;
  char *ci;
  i = strlen(c);
  if (i==0)
    return -1;
  ci = c;
  while (*ci != 0 && iswhite )
  {
    iswhite &= (isspace(*ci) != 0);
    ci++;
  }
  return iswhite;
}

/* find next whitespace, after next non-whitespace */
char *
nextwhite (char *c)
{
  int nonw;
  nonw = !isspace (*c);
  while (*c != '\0')
  {
    while ((nonw == 1 && !isspace (*c)) || (nonw == 0 && isspace (*c)))
    {
      c++;
      if (nonw == 0 && !isspace (*c))
        nonw = 1;
    }
    return c;
  }
  return c;
}                               /* nextwhite */

/* finds the next non-space character after the next space */
char *
nextnonspaceafterspace (char *textline)
{
  char *cc;
  if (textline == NULL)
    return NULL;
  cc = textline;
  while (*cc != ' ' && *cc != '\0')
    cc++;
  while (*cc == ' ')
    cc++;
  if (*cc == '\0')
    return NULL;
  else
    return cc;
}

/* finds the next non-space character */
char *
nextnonspace (char *textline)
{
  char *cc;
  if (textline == NULL)
    return NULL;
  cc = textline;
  while (*cc == ' ')
    cc++;
  if (*cc == '\0')
    return NULL;
  else
    return cc;
}


int
comma_exists (char *s)
{
  size_t l;
  size_t i;
  int v;

  v = 0;
  l = strlen (s);
  for (i = 0; i < l; i++)
  {
    if (s[i] == ',' || s[i] == ';')
    {
      v = 1;
    }
  }
  return v;
}

/* remove spaces and count parentheses */
void
checktreestring (char *t)
{
  int closepcheck, openpcheck, i;

  for (closepcheck = 0, openpcheck = 0, i = 0; (unsigned) i < strlen (t); i++)  // check if the number of open parenthese match the close parentheses
  {
    closepcheck += t[i] == ')';
    openpcheck += t[i] == '(';
    if ((t[i] == ' ') || (t[i] == '\t') || (t[i] == '\n') || (t[i] == '\r'))    // remove spaces
    {
      t[i] = '\0';
      strcat (t, &(t[i + 1]));
      i--;
    }
  }
  if (closepcheck != openpcheck)
  {
    //  errr (-1, -1, 17, closepcheck, openpcheck);
    IM_err (IMERR_POPTREESTRINGFAIL,
            " wrong number of parentheses, string %s,  open '(' %d  close ')' , step %d",
            t, openpcheck, closepcheck, step);
  }
}                               /* checktreestring */

void
copymig (struct migstruct *m1, struct migstruct *m2)
{
  m1->mp = m2->mp;
  m1->mt = m2->mt;
}

/************************/
/* Other Misc Functions */
/************************/

/** 
 * logfact is absurdly large, but this is because it must be usable 
 * for all of the migration events in a * data set, which might be very large
 * for a large data set with large upper bounds on migration rates.
 * We need log factorial of n + 1/2 where n a is nonnegative integer.
 * loghalffact (n) := log (n + 1/2)! where n starts from 0 to 50 * ABSMIGMAX.
 */
/* copy to utilities.h
#ifndef M_LNPI
#  define M_LNPI 1.14472988584940017414342735135
#endif
#ifndef M_LN2
#  define M_LN2 0.69314718055994530941723212146
#endif */

void
setlogfact (void)
{
  int i;
  //UByteP a;   leftover from assignment stuff, required xtrabbits
  // int bitj;   leftover from assignment stuff, required xtrabbits
  //int ic;
  logfact[0] = 0;
  for (i = 1; i <= LOGFACTMAX-1; i++)
    logfact[i] = ((double) logfact[i - 1]) + log ((double) i);
  /*
  for (i = 0; i < 50 * ABSMIGMAX; i++)
  {
    loghalffact[i] = M_LNPI / 2
      + logfact[2 * i + 2] - logfact[i + 1] - (2 * i + 2) * M_LN2;
  }
  a = (UByteP) malloc (sizeof (unsigned char));    //leftover from assignment stuff, required xtrabbits
  for (ic = 0; ic < 256; ic++)
  {
    BitZero (a, 1);
    a[0] = (unsigned char) ic;
    BITNUMBERTRUE[ic] = 0;
    for (i = 0; i < 8; i++)
    {
      if (BitIsTrue (a, i))
      {
        BITNUMBERTRUE[ic] += 1;
      }
    }
  }
  XFREE (a);
  a = NULL; */
}

void
ieevent (struct eevent *a)      // initialize structure
{
  a->n = 0;
  a->s = 0;
  a->ss = 0;
  a->s2 = 0;
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
double
bessi (int n, double x)
{
  int j;
  double bi, bim, bip, tox, ans;
  assert (x >= 0);
  n = abs (n);
  if (x > 700)
    return (MYDBL_MAX);
  if (n == 0)
    return bessi0 (x);
  if (n == 1)
    return bessi1 (x);
  if (x == 0.0)
  {
    return 0.0;
  }
  else
  {
    tox = 2.0 / fabs (x);
    bip = ans = 0.0;
    bi = 1.0;
    for (j = 2 * (n + (int) sqrt (ACC * n)); j > 0; j--)
    {
      bim = bip + j * tox * bi;
      bip = bi;
      bi = bim;
      if (fabs (bi) > BIGNO)
      {
        ans *= BIGNI;
        bi *= BIGNI;
        bip *= BIGNI;
      }
      if (j == n)
        ans = bip;
    }
    ans *= bessi0 (x) / bi;
    return x < 0.0 && (n & 1) ? -ans : ans;
  }
}

#undef ACC
#undef BIGNO
#undef BIGNI

/* eexp is used, generally together with struct extendnum to calculated exponentials of sums of numbers 
that may have very large exponents */ 
#define LOG_10_2  0.30102999566398119521
static   double us[10] = { 1.0, 0.5, 0.16666666666666666666666666667,
  0.04166666666666666666666666667, 0.00833333333333333333333333333,
  0.001388888888888888888888888889, 0.000198412698412698412698412698,
  0.000024801587301587301587301587301, 2.75573192239858906525573192239859e-6,
  2.75573192239858906525573192239e-7};
void
eexp (double x, double *m, int *z)
{
  double u, zr, temp;
  int n;
  n = (int) floor (x / LOG2);
  zr = LOG_10_2 * (double) n;
  *z = (int) zr;
  zr -= (double) *z;
  u = x - (((double) n) * LOG2);
  temp =
    1 + u * (us[0] +
             u * (us[1] +
                  u * (us[2] +
                       u * (us[3] +
                            u * (us[4] +
                                 u * (us[5] +
                                      u * (us[6] +
                                           u * (us[7] +
                                                u * (us[8] +
                                                     u * (us[9]))))))))));
  *m = temp * pow (10.0, zr);
  if (fabs (*m) > 10)
  {
    *m = *m / 10.0;
    *z = *z + 1;
  }

  if (fabs (*m) < 1)
  {
    *m = *m * 10.0;
    *z = *z - 1;
  }
}                               //eexp

int **
alloc2Dint (int rows, int cols)
{
  int **a;
  int i;
  if ((a = static_cast<int **> (malloc (rows * sizeof (*a)))) == NULL)
    IM_err (IMERR_MEM, "  alloc2Dint main malloc did not work.  step %d",
            step);
  for (i = 0; i < rows; i++)
    if ((a[i] = static_cast<int *> (malloc (cols * sizeof (*a[i])))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc2Dint loop malloc did not work. rows %d,  cols %d, step %d",
              rows, cols, step);
  return a;
}

double **
orig2d_alloc2Ddouble (int rows, int cols)
{
  double **a;
  int i;
  if ((a = static_cast<double **> (malloc (rows * sizeof (*a)))) == NULL)
    IM_err (IMERR_MEM, "  alloc2Ddouble main malloc did not work.  step %d",
            step);
  for (i = 0; i < rows; i++)
    if ((a[i] = static_cast<double *> (malloc (cols * sizeof (*a[i])))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc2Ddouble loop malloc did not work. rows %d,  cols %d, step %d",
              rows, cols, step);
  return a;
} 

double **alt2d_alloc2Ddouble(long m, long n)
// got this from http://www.materialography.de/Book/sourcecodes/malloc2d.c 
/* a double matrix x[0..m-1][0..n-1] is allocated */ 
{  long i;
   double **x;
   
   x=(double **)malloc((size_t)(m*sizeof(double *)));
   x[0]=(double *)malloc((size_t)(m*n*sizeof(double)));
   for(i=1l;i<m;i++) x[i]=x[i-1]+n;
   return x;
}

int **alt2d_alloc2Dint(long m, long n)
// got this from http://www.materialography.de/Book/sourcecodes/malloc2d.c 
/* a int matrix x[0..m-1][0..n-1] is allocated */ 
{  long i;
   int **x;
   
   x=(int **)malloc((size_t)(m*sizeof(int *)));
   x[0]=(int *)malloc((size_t)(m*n*sizeof(int)));
   for(i=1l;i<m;i++) x[i]=x[i-1]+n;
   return x;
}


/* cr 110907.1    added function
 *  allocate the space for a 2 dimensional array (rows x cols) of 
 *  doubles, and init each element to 0.0
 */
double **
allocAndInit2Ddouble (int rows, int cols)
{

  double **a;
  int i,j;
  if ((a = static_cast<double **> (malloc (rows * sizeof (double *)))) == NULL)
    IM_err (IMERR_MEM, 
        "  allocAndInit2Ddouble main malloc failed!  step %d", step);
  for (i = 0; i < rows; i++)
  {
    if ((a[i] = static_cast<double *> 
                (malloc (cols * sizeof (double)))) == NULL)
      IM_err (IMERR_MEM,
        "  allocAndInit2Ddouble loop malloc failed! rows %d,  cols %d, step %d",
              rows, cols, step);
    for (j = 0; j < cols; j++)
    {
      a[i][j] = 0.0;
    }
  }
  return a;
} 

void
orig2d_free2D (void **a, int rows)
{
  int i;
  for (i = 0; i < rows; i++)
    XFREE (a[i]);
  XFREE (a);
} 

void alt2d_free2D(double **x/*, long m , long n */)
// modified this from http://www.materialography.de/Book/sourcecodes/malloc2d.c 
/* free the double matrix x[0..m-1][0..n-1] */ 
{  XFREE(x[0]); XFREE(x); } 

void alt2d_free2Dint(int **x/*, long m , long n */)
// modified this from http://www.materialography.de/Book/sourcecodes/malloc2d.c 
/* free the double matrix x[0..m-1][0..n-1] */ 
{  XFREE(x[0]); XFREE(x); } 

/* make a 2D matrix that is truly contiguous
  requires two variables both of which must be malloced and later freed.
  the main variable is the address of the 2D matrix
  the aardata variable is used so that all the memory for the ints is contiguous
  aardata could be an *int  (not that *arrdata is malloced) but this address must be passed back to the calling function and later used for freeing,  so it is **int 
  use free_2Dint_contiguous() */

int **alloc2Dint_contiguous(int **arrdata, int rows, int columns)
{
  int **arr;
  if ((arr = static_cast<int **> (malloc (rows * sizeof (*arr)))) == NULL)
    IM_err (IMERR_MEM, "  alloc2Dint_contiguous main malloc did not work.  step %d",
            step);
  if ((*arrdata = static_cast<int *> (malloc (columns * rows * sizeof (**arrdata)))) == NULL)
      IM_err (IMERR_MEM,
              "  alloc2Dint_contiguous did not work. rows %d,  cols %d, step %d",
              rows, columns, step);
  for(int i = 0; i < rows; i++)
     arr[i]  = *arrdata + i * columns ;
  return arr;
}

void free_2Dint_contiguous(int **arr, int *arrdata)
{
  XFREE(arr);
  XFREE(arrdata);
}

/* not used as of 4/27/2017,  not sure why it is in here
 * Function logsum approximately computes the value of v of the following
 * equation:
 * exp(v) = exp(a_1) + exp(a_2) + ...
 * The first arguement is the number of elements of the sum, and the rest of the
 * argements are the n elements. For example, we can compute 
 * exp(3.0) + exp(1.5) + exp(2.5), which is exp(3.604131) by using function call
 * logsum (3, 3.0, 1.5, 2.5). The returned value should be 3.604131.
 */
double logsum (int n, ...)
{
  register int i;
  double a;
  double b;
  va_list ap;

  assert (n > 1);
  va_start (ap, n);
  a = va_arg (ap, double);
  for (i = 1; i < n; i++)
  {
    b = va_arg (ap, double);

    if (a > b)
    {
      if (a - b < LOG_DBL_MAX)
        a = b + log (exp (a - b) + 1.0);
    }
    else
    {
      if (b - a < LOG_DBL_MAX)
        a = a + log (exp (b - a) + 1.0);
      else
        a = b;
    }
  }
  va_end (ap);
  return a;
}

/*  version of logsum that uses an array 
 not used as of 4/27/2017,  not sure why it is in here
*/
double logsuma (int n, double *d) 
{
  register int i;
  double a;
  double b;

  assert (n > 1);
  a = d[0];
  for (i = 1; i < n; i++)
  {
    b = d[i];

    if (a > b)
    {
      if (a - b < LOG_DBL_MAX)
        a = b + log (exp (a - b) + 1.0);
    }
    else
    {
      if (b - a < LOG_DBL_MAX)
        a = a + log (exp (b - a) + 1.0);
      else
        a = b;
    }
  }
  return a;
}

/* Find the default directory based on loadfilebase.
  * /this/directory/a.out -> defaultdir is /this/directory/
  * /a.out                -> /
  * a.out                 -> ./
  */
int
imaDirBase (const char *a, char **b)
{
  size_t len;
  int i;  // 8/26/2012  was size_t   but needs to be negative or zero 
  
  len = strlen (a);
  for (i = (int) (len - 1); i >= 0; i--)
    {
      if (a[i] == '/' || a[i] == '\\' )
        {
          break;
        }
    }
  if (i < 0)
    {
      *b = static_cast<char *> (malloc (3 * sizeof (char)));
      strcpy (*b, "./");
    }
  else
    {
      *b = static_cast<char *> (malloc ((i+2) * sizeof (char)));
      strncpy (*b, a, i + 1);      
      (*b)[i+1] = '\0';
    }
  
  return 0;
} 

/* paths that have a space (' ') in the file system cannot be read on the command line
instead the user must change all spaces to "%20" when preparing the command line.  
This function puts the space back */ 
int put_spaces_in_filepaths(char *pathstr)
{
  char *c;
  while (c = strstr(pathstr,"%20"))
  {
    *c = ' ';
    strdelete(pathstr,c-pathstr + 2,2);
  }
  return 0;
}

/* seems to shorten a string of a number in exponential notation */ 
char *
shorten_e_num (char *s)
{
  char *pos;
  char e = 'e';
  int i;
  pos = strchr (s, e);
  i = (int) (pos-s);
  while (s[i+1] == '0' ||s[i+1] == '+')
  {
    strdelete (s, i+2,1);
  }
  return s;
}

/* some functions sang chul wrote for reading from a text file */

int
skip_a_line (FILE *fp)
{
  char c;
  c = 'x';
  while (c != '\n')
    {
      c = fgetc(fp);
    }
  return 0;
}

int
read_int (FILE *fp)
{
  int v;
  char c;
  char *buf;
  int len;
  int i;
  
  len = 5;
  buf = static_cast<char *> (malloc (len * sizeof (char)));
  
  c = fgetc (fp);
  i = 0;
  while (c != '\t' && c != '\n')
    {
      buf[i++] = c;
      if (i == len)
        {
          len += 5;
          buf = static_cast<char *> (realloc (buf, len * sizeof (char)));
        }
      c = fgetc (fp);
    }
  buf[i] = '\0';
  v = atoi (buf);
  
  free (buf);
  buf = NULL;
  return v;
}

double 
read_double (FILE *fp)
{
  double v;
  char c;
  char *buf;
  int len;
  int i;
  
  len = 5;
  buf = static_cast<char *> (malloc (len * sizeof (char)));
  
  c = fgetc (fp);
  i = 0;
  while (c != '\t' && c != '\n')
    {
      buf[i++] = c;
      if (i == len)
        {
          len += 5;
          buf = static_cast<char *> (realloc (buf, len * sizeof (char)));
        }
      c = fgetc (fp);
    }
  buf[i] = '\0';
  v = atof (buf);
  
  free (buf);
  buf = NULL;
  
  return v;
}

void convertToUpperCase(char *sPtr)
{
  while(*sPtr != '\0')
  {
    *sPtr = toupper((unsigned char)*sPtr);
    ++sPtr;
  }
}

void setupdatescalarinfo(struct updatescalarinfo *a, double sav, double saav, double savmin, double savmax, double target, int numattemptscheck)
{
  a->updatescalarval = sav;
  a->updatescalaradjustval = saav;
  a->updatescalarvalmin = savmin;
  a->updatescalarvalmax = savmax;
  a->targetupdaterate = target;
  a->numattemptscheck = numattemptscheck;
  a->recentaccp = 0;
  a->recenttries = 0;
  a->allaccp = 0;
  a->alltries = 0;
}

/* checks to see if the udpate rate is close to the target and if it is not change it */
void checkupdatescalarer(struct updatescalarinfo *a)
{
  static double minfactor= 0.9;
  static double maxfactor= 1.1;
  
  double rate; 
  if (a->recenttries >= a->numattemptscheck) // && a->alltries < 10000)
  {
    rate =   static_cast<double>(a->recentaccp)/static_cast<double>(a->recenttries);
    if (rate < minfactor * a->targetupdaterate)// rate too low so reduce the scalar 
    {
      a->updatescalarval = DMAX(a->updatescalarvalmin, a->updatescalarval/a->updatescalaradjustval);
    }
    else if (rate > maxfactor * a->targetupdaterate)// rate too high so increase the scaler 
    {
       a->updatescalarval = DMIN(a->updatescalarvalmax, a->updatescalarval*a->updatescalaradjustval);
    }
    a->recenttries = 0;
    a->recentaccp = 0;
  }
  a->recenttries += 1;
  a->alltries += 1; 
  
} // checkupdatescalarer

/*void resetupdatescalarer(struct updatescalarinfo *a)
{
  if (a->c == 1)
  {
    a->updatescalarval   =  a->holdupdatescalarval;
    a->recenttries = a->holdrecenttries;
    a->recentaccp = a->holdrecentaccp;
    a->c = 0;
  }
} //resetupdatescalarer */

int  file_exists(const char * filename) {
    if (FILE * file = fopen(filename, "r")) {
        fclose(file);
        return 1;
    }
    return 0;
}

/*********** dictionary code modified from Kernigan & Ritchie 2nd ed sec 6.6 **********/
// use a list of pointers to nodes,  
// for collisions each node becomes the head of a linked list of nodes that share that hash value 
// see struct dictionary_node_kr  in ima.hpp
// see also  void freekrlinkedlist(struct dictionary_node_kr* head) in update_priors.cpp for freeing things

/* hash: form hash value for string s  - given a string s returns a number from 0 to hashsize
simple function out of K&R 2nd ed
*/

unsigned hashkr(char *s)
{
    unsigned hashval;
    for (hashval = 0; *s != '\0'; s++)
      hashval = *s + 31 * hashval;
    return hashval % hashsize ;
}

/* change this so hashtab is passed to the function  */
/* the same as dictionarylookup()  but returns the value */ 
double getvalue(char *key, struct dictionary_node_kr **hashtab)
{
    struct dictionary_node_kr *np;
    for (np = hashtab[hashkr(key)]; np != NULL; np = np->next)
        if (strcmp(key, np->name) == 0)
        {
          return np->priorval;
        }
    return -1;
}

/* dictionarylookup: look for s in hashtab  - change this so hashtab is passed to the function
  this is used when building the dictionary. use this for simple lookup 
*/
struct dictionary_node_kr *dictionarylookup(char *key, struct dictionary_node_kr **hashtab)
{
    struct dictionary_node_kr *np;
    for (np = hashtab[hashkr(key)]; np != NULL; np = np->next)
        if (strcmp(key, np->name) == 0)
          return np; /* found */
    return NULL; /* not found */
}

char *strdup(char *s) /* make a duplicate of s */
{
    char *p;
    assert(strlen(s) < 100);
    p = (char *) malloc(strlen(s)+1); /* +1 for ’\0’ */
    if (p != NULL)
       strcpy(p, s);
    return p;
}

/* dictionary_install: put (name, i) in hashtab */
/* change this so hashtab is passed to the function  */
/*  collisions are handled by generating a linked list at position hashval  in hashtab_i of all things that share that hashval */
struct dictionary_node_kr *dictionary_install(char *name, double val,struct dictionary_node_kr **hashtab)  // returns an nlist pointer
{
    struct dictionary_node_kr *np;
    unsigned hashval;

    if ((np = dictionarylookup(name,hashtab)) == NULL) { /* not found */
        np = (struct dictionary_node_kr *) malloc(sizeof(*np));
        if (np == NULL || (np->name = strdup(name)) == NULL)
          return NULL;
        hashval = hashkr(name);
        np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
    np->priorval = val;
    return np;
}

/*********** end of dictionary code modified from Kernigan & Ritchie 2nd ed sec 6.6 **********/

/* 
  metropolishastingsdecide() returns 0 (reject) or 1 (accept) the update
  uses the logarithm of the metropolishastings ratio
   if logmhratio is out of bounds, reject
   othercriteria allows for various other things to also be true
   othercriteria is 0 or 1,  if 0 reject
 */
int metropolishastingsdecide(double logmhratio,int othercriteria)
{
  double logU;
  if (isinf_DBL(logmhratio) || othercriteria == 0)
    return 0;
  if (logmhratio >= 0.0) // mh ratio > 1
    return 1;
   // could either take log of random value or exp of logmhratio 
  logU = log(uniform());
  if (logU < logmhratio)
    return 1;
  else
    return 0;
}

/* return point to a string with time info
  should be safe since ts is static and should not dissappear after exiting */ 
char* timestring(time_t seconds)
{
  static char ts[50];
  sprintf(ts,"%d hours, %d minutes, %d seconds",
        (int) seconds / (int) 3600,
          ((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),
              (int) seconds - (int) 60 *((int) seconds / (int) 60));
   ts[sizeof(ts)-1] = '\0';  // just in case it has written past the end 
   return ts;
}