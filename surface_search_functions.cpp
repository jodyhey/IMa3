/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* this file is numerical recipes stuff with some small modifications */
#undef GLOBVARS
#include "ima.hpp"              // may not really be needed in this particular file

/*********** LOCAL STUFF **********/
// these 3 defines don't seem to get used 
//#define NRANSI
//#define NR_END 1
//#define FREE_ARG char*

/* local function prototypes */

static double sign (double a, double b);

/* local functions */

static double
sign (double a, double b)
{
  if (b > 0.0)
    return fabs (a);

  else
    return (-fabs (a));
}

/************ GLOBAL FUNCTIONS ******************/

/* NR stuff for finding minimum using golden section */

#define goldratio  1.61803399
#define x100       100.0
#define smallval   1.0e-20
#define grless1    0.61803399

void
bracket (int ndim, int firsttree, int lasttree, double *axval, double *bxval,
           double *cxval, double *funca, double *funcb, double *funcc,
           double (*func) (int, int, int, double, int), int ifneeded)
{

  /* Programs using routine MNBRAK must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL FOR which a minimum is TO be found */
  /* I converted this to a function and included code to check for
     a failure to find a braket */
  double ulim, uval, rr, q, funcu, dum;
  *funca = func (ndim, firsttree, lasttree, *axval, ifneeded);
  *funcb = func (ndim, firsttree, lasttree, *bxval, ifneeded);
  if (*funcb > *funca)
  {
    dum = *axval;
    *axval = *bxval;
    *bxval = dum;
    dum = *funcb;
    *funcb = *funca;
    *funca = dum;
  }
  *cxval = fabs (*bxval + goldratio * (*bxval - *axval));
  *funcc = func (ndim, firsttree, lasttree, *cxval, ifneeded);
_L1:
  /* need to add some excapes for floating point problems */
  if (*funcb >= *funcc && *funcb > -DBL_MAX && !(*funcb == 0 && *funcc == 0))
  {
    rr = (*bxval - *axval) * (*funcb - *funcc);
    q = (*bxval - *cxval) * (*funcb - *funca);
    uval = (double) *bxval - ((*bxval - *cxval) * q - (*bxval - *axval) * rr) /
      (2.0 * sign (FMAX (fabs (q - rr), (float) smallval), q - rr));
    ulim = *bxval + x100 * (*cxval - *bxval);
    if ((*bxval - uval) * (uval - *cxval) > 0.0)
    {
      funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
      if (funcu < *funcc)
      {
        *axval = *bxval;
        *funca = *funcb;
        *bxval = uval;
        *funcb = funcu;
        goto _L1;
      }
      if (funcu > *funcb)
      {
        *cxval = uval;
        *funcc = funcu;
        goto _L1;
      }
      uval = *cxval + goldratio * (*cxval - *bxval);
      funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
    }
    else if ((*cxval - uval) * (uval - ulim) > 0.0)
    {
      funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
      if (funcu < *funcc)
      {
        *bxval = *cxval;
        *cxval = uval;
        uval = *cxval + goldratio * (*cxval - *bxval);
        *funcb = *funcc;
        *funcc = funcu;
        funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
      }
    }
    else if ((uval - ulim) * (ulim - *cxval) >= 0.0)
    {
      uval = ulim;
      funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
    }
    else
    {
      uval = *cxval + goldratio * (*cxval - *bxval);
      funcu = func (ndim, firsttree, lasttree, uval, ifneeded);
    }
    *axval = *bxval;
    *bxval = *cxval;
    *cxval = uval;
    *funca = *funcb;
    *funcb = *funcc;
    *funcc = funcu;
    goto _L1;
  }
  //if (*funcb < -DBL_MAX)
    //printf (" minus infinity in bracket \n");
  assert (!isninf_DBL(*funcb));
      
}

#undef goldratio
#undef x100
#undef smallval

double
goldenmin (int ndim, int firsttree, int lasttree, double axval, double bxval,
           double cxval, double tol_, double *xmin, double (*func) (int, int,
                                                                 int, double, int), int ifneeded)
{

  /* Programs using routine GOLDEN must supply an EXTERNAL
     FUNCTION func(x:REAL):REAL whose minimum is TO be found. */
  double Result, f1, f2, c, x0, x1, x2, x3;
  c = 1.0 - grless1;
  x0 = axval;
  x3 = cxval;
  if (fabs (cxval - bxval) > fabs (bxval - axval))
  {
    x1 = bxval;
    x2 = bxval + c * (cxval - bxval);
  }
  else
  {
    x2 = bxval;
    x1 = bxval - c * (bxval - axval);
  }
  f1 = func (ndim, firsttree, lasttree, x1, ifneeded);
  f2 = func (ndim, firsttree, lasttree, x2, ifneeded);
  while (fabs (x3 - x0) > tol_ * (fabs (x1) + fabs (x2)) && (x3 > -DBL_MAX))
  {
    if (f2 < f1)
    {
      x0 = x1;
      x1 = x2;
      x2 = grless1 * x1 + c * x3;
      f1 = f2;
      f2 = func (ndim, firsttree, lasttree, x2, ifneeded);
      continue;
    }
    x3 = x2;
    x2 = x1;
    x1 = grless1 * x2 + c * x0;
    f2 = f1;
    f1 = func (ndim, firsttree, lasttree, x1, ifneeded);
  }
  if (f1 < f2)
  {
    Result = f1;
    *xmin = x1;
  }
  else
  {
    Result = f2;
    *xmin = x2;
  }
  if (x3 < -DBL_MAX)
    printf (" minus infinity in goldenmin \n");
  return Result;
}
#undef grless1
