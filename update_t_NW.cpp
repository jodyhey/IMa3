/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

#undef GLOBVARS

#include "ima.hpp"
#include "update_gtree_common.hpp"
#include "update_gtree.hpp"

/* this file holds stuff for changet_NW()
This is a t update modeled on the way Nielsen and Wakeley  (2001)did it as implemented in the midv program,  by moving a split time up or down
and adding or erasing migration events. 
check out mathematica file: newt_t_updating_8_21_08.nb 
*/

/*********** local to this file  ***************/

/* struct migrationinfo_tNW  
holds info on each of the simulations that are done in the update
each element can hold info on one or two edges 

a static array, large enough to hold the largest of genealogies in the data set, 
is set up the first time through 
*/

struct migrationinfo_tNW
{
  int n;                            // 1 or 2 for # of edges
  int edgeid[2];                    // edge numbers
  int upa[2];                      // the population state at the top of (but just below) the split interval, after the update
  int upb[2];                      // the population state at the top of (but just below) the split interval, before the update
  int da;                       // the population state just before the bottom of the split interval, after the update
  int db;                             // the population state just before the bottom of the split interval, before the update
  double uptime[2];               // time at top of edge
  double mtime[2];                 // time of edge that is in the interval between the new and old split times
  int mcount[2];                   // # of migrations in that time interval before the update
  int mnew[2];                     // # of migrations in that interval after the uptdate
  int npopsb;                   // # of populations in period before
  int npopsa;                   // # of populations in period after
  int cm2_b[2];                    // population state before the second to last migration before the update
  int cm2_a[2];                    // population state before the second to last migration after  the update
  double mrate[2];                 // migration rate over that interval before udpate 
  double mrate_r[2];               // migration rate over that interval after udpate 
  double logpfpop;                 // probability of simulating da for the update
  double logpfpop_r;               // probability of simulating da for the reverse update 
  int cmm[2];                     // # of migrations on edges, not counting mcount
} *minfo;

static struct genealogy *holdallgtree_t_NW;
static int largestsamp;
static int ci;
static struct genealogy *G;
static double oldt, newt, tu, td;

static struct genealogy_weights holdallgweight_t_NW;
static struct genealogy_weights holdgweight_t_NW[MAXLOCI];
static struct probcalc holdallpcalc_t_NW;

static int mforward, mreverse;
/******* local prototypes ********/
static void copy_all_gtree (int mode);
static int simmpath_tNW (int edge, struct edge *gtree, int numm, int lastm,
                         double timein, double upt, int period, int pop,
                         int constrainpop);
void initminfo(int j);
static void addmigration_NW(struct edge *gtree, int period, int numupdates);
static double getmprob_NW(int period, double mrate, double mtime, int mcount, int uppop, int dpop, int cm2pop, int numpops);
static double update_mig_tNW (int li, struct edge *gtree, int period);

// copy_all_gtree  makes a copy of all genealogies in a given chain 
// just like copy_gtree() but for all loci in a chain
// could set it up so it calls copy_gtree() but that function 
// uses a local holdgtree variable - could change this to cut down on code
// if mode==1  copy *G into holdgtree  if mode==0 copy holdgtree into *G
// used only for NW update
void
copy_all_gtree (int mode)
{
  struct edge *togtree, *fromgtree;
  struct genealogy *fromG, *toG;
  int li, i, j, ai;

  for (li = 0; li < nloci; li++)
  {
    if (mode)
    {
      toG = &holdallgtree_t_NW[li];
      fromG = &C[ci]->G[li];
    }
    else
    {
      toG = &C[ci]->G[li];
      fromG = &holdallgtree_t_NW[li];
    }
    togtree = toG->gtree;
    fromgtree = fromG->gtree;
    toG->length = fromG->length;
    toG->mignum = fromG->mignum;
    toG->root = fromG->root;
    toG->roottime = fromG->roottime;
    toG->tlength = fromG->tlength;
    for (i = 0; i < L[li].numlines; i++)
    {
      togtree[i].down = fromgtree[i].down;
      j = -1;
      do
      {
        j++;
        togtree[i].mig[j] = fromgtree[i].mig[j];
      }
      while (fromgtree[i].mig[j].mt > -0.5);
      togtree[i].cmm = j;
      togtree[i].time = fromgtree[i].time;
      togtree[i].pop = fromgtree[i].pop;
      togtree[i].fpop = fromgtree[i].fpop; // fix a bug ? 
      togtree[i].up[0] = fromgtree[i].up[0];
      togtree[i].up[1] = fromgtree[i].up[1];
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
        {
          togtree[i].A[ai] = fromgtree[i].A[ai];
          togtree[i].dlikeA[ai] = fromgtree[i].dlikeA[ai];
        }
    }
  }
}                               /* copy_all_gtree */

/* simulate the migration path  -  very similar to simmpath() */
int
simmpath_tNW (int edge, struct edge *gtree, int numm, int lastm,
              double timein, double upt, int period, int pop,
              int constrainpop)
{
  int i, lastpop, startm, lastpop_2;
                    /* CR 110715.1 */
  int dupCheck;     /* flag turns off dup migration time check of mig events */
  int migIndex;     /* index used to look for duplicate migration times */

  assert (numm > 0);
  startm = lastm + 1;
  lastm = lastm + numm;

  do
  {
    for (i = startm; i <= lastm; i++)
    {
      gtree[edge].mig[i].mt = upt + uniform () * timein;

    }
    gtree[edge].mig[i].mt = -1;
    gtree[edge].cmm = i;
    dupCheck=0;
    if (numm > 1)
    {
      hpsortmig (&gtree[edge].mig[startm] - 1, numm);

      /* look for duplicate migration times in sorted list.
       * This solves a charateristic of the Mersennes Twister 
       * random number generator in which identical random numbers may be
       * returned from the random number sequence in a very small number
       * of calls.  With some seeds it was noted as small as within 4 calls.
       */
      for (migIndex = startm; migIndex < lastm; ++migIndex)
      {
        if  ( gtree[edge].mig[migIndex].mt != gtree[edge].mig[migIndex + 1].mt )
        {
           continue;
        }
        else
        {  /* if duplicate time found, a new migration path must be simulated */
           dupCheck=1;
           break;
        }
      }
    }
    else
    {
      dupCheck=0;
    }
  } while (dupCheck == 1);    /* when no duplicate times found, exit loop  */

  lastpop = pop;
  lastpop_2 = -1;
  if (constrainpop < 0)
  {
    for (i = startm; i <= lastm; i++)
    {
      gtree[edge].mig[i].mp =
        picktopop (lastpop, C[ci]->plist[period], npops - period);
      lastpop = gtree[edge].mig[i].mp;
if (step==1542608) printf("%d %d %.3lf %d\t",edge,i,gtree[edge].mig[i].mt,lastpop);
    }
  }
  else
  {
    if (numm >= 2)
    {
      i = startm;
      while (i < lastm - 1)
      {
        gtree[edge].mig[i].mp =
          picktopop (lastpop, C[ci]->plist[period], npops - period);
        lastpop = gtree[edge].mig[i].mp;
if (step==1542608) printf("%d %d %.3lf %d\t",edge,i,gtree[edge].mig[i].mt,lastpop);
        i++;
      }
      lastpop_2 = lastpop;
      
      gtree[edge].mig[lastm - 1].mp =
        picktopop2 (lastpop, C[ci]->plist[period], npops - period,
                    constrainpop);
      gtree[edge].mig[lastm].mp = constrainpop;
if (step==1542608) printf("%d %d %.3lf %d\t",edge,lastm,gtree[edge].mig[lastm].mt,constrainpop);
    }
    else
    {
      if (numm == 1)
      {
        gtree[edge].mig[lastm].mp = constrainpop;
        lastpop = constrainpop;
if (step==1542608) printf("%d %d %.3lf %d\t",edge,lastm,gtree[edge].mig[lastm].mt,constrainpop);
      }
    }
  }
  return lastpop_2;
}                               /* simmpath_tNW */

 
void addmigration_NW(struct edge *gtree, int period, int numupdates)
{
  int i, j, k, ei, mi, mii;
  int lastabovem, numskip, numheld;
  struct migstruct holdmig[MIGARRAYSIZE]; 

  for (j=0;j<numupdates;j++)
  {
    for (i=0;i<minfo[j].n;i++)
    {
      ei = minfo[j].edgeid[i];
      if (ei != G->root)
      {
       /* determine position in mig array of the first migration event in the relevant period */
        mi = 0;
        if (minfo[j].uptime[i] < tu)
          while (gtree[ei].mig[mi].mt > -0.5 && gtree[ei].mig[mi].mt < tu)
          {
            mi++;
          }
        lastabovem = mi - 1;           // position to start simulating

        /* count how many migrations must be erased for the update */
        mii = mi;
        while (gtree[ei].mig[mii].mt > -0.5 && gtree[ei].mig[mii].mt < td)
        {
          mii++;
        }
        numskip = mii - mi;

        /* count and save the migrations after the relevant period  */
        if (gtree[ei].mig[mii].mt < 0)
        {
          numheld = 0;
        }
        else
        {
          k = 0;
          while (gtree[ei].mig[mii].mt > -0.5)
          {
            holdmig[k] = gtree[ei].mig[mii];
            mii++;
            k++;
          }
          holdmig[k].mt = -1;
          numheld = k;
        }
        if (newt < td && period == lastperiodnumber)  // do not add migration
        {
          assert(minfo[j].mnew[i] == 0);
          gtree[ei].mig[lastabovem + 1].mt = -1;
          gtree[ei].cmm = lastabovem + 1; // fix bug, this was missing 
        }
        else
        {
          if (minfo[j].mnew[i] > 0)
          {
            minfo[j].cm2_a[i] =
              simmpath_tNW (ei, gtree, minfo[j].mnew[i], lastabovem, minfo[j].mtime[i],
                            DMAX (tu, minfo[j].uptime[i]), period, minfo[j].upa[i],
                            minfo[j].da);
            /* put back stored migration events */
            if (numheld > 0)
            {
              mi = lastabovem + minfo[j].mnew[i] + 1;
              for (k = 0; k < numheld; k++, mi++)
              {
                gtree[ei].mig[mi] = holdmig[k];
              }
              gtree[ei].mig[mi].mt = -1;
              gtree[ei].cmm = mi;
            }
          }
          else
          {
            if (numskip > 0)
            {
              /* put back stored migration events */
              if (numheld > 0)
              {
                for (mi = lastabovem + 1, k = 0; k < numheld; k++, mi++)
                  gtree[ei].mig[mi] = holdmig[k];
                gtree[ei].mig[mi].mt = -1;
                gtree[ei].cmm = mi;
              }
              else
              {
                gtree[ei].mig[lastabovem + 1].mt = -1;
                gtree[ei].cmm = lastabovem + 1;
              }
            }
          }
        }

      }
    }
  }
}  // addmigration_NW

double 
getmprob_NW(int period, double mrate, double mtime, int mcount, int uppop, int dpop, int cm2pop, int numpops)
{
  double logb, logs, lognp; 
  if (period ==lastperiodnumber)
    return 0;
  if (period == lastperiodnumber - 1)
  {
    if (ODD(mcount))
      return  mcount * log(mrate/mtime) - mylogsinh (mrate);
    else
      return  mcount * log(mrate/mtime) - mylogcosh (mrate);
  }
  // else period > lastperiodnumber - 1
  assert(period < lastperiodnumber - 1);
  if (uppop == dpop)
    logs = log (1 - mrate * exp (-mrate));
  else
    logs = log (1 - exp (-mrate));
  switch (mcount)
  {
  case  0 : 
    return -mrate - logs;
  case  1 :
    return log(mrate/mtime) - mrate - logs;
  default :
    {
      assert(cm2pop >= 0);
      assert(mcount > 0);
      lognp = log(static_cast<double>(numpops)-1);
      if (cm2pop == dpop)
        logb = - lognp;
      else
        logb = - log(static_cast<double>(numpops)- 2);
      return mcount * log(mrate/mtime) + (2- mcount)* lognp + logb - mrate - logs;

    }
  }
} //getmprob_NW

void initminfo(int j)
{
  minfo[j].n = minfo[j].edgeid[0] = minfo[j].edgeid[1] = -1;
  minfo[j].upa[0] = minfo[j].upa[1] = minfo[j].upb[0] = minfo[j].upb[1] = -1;
  minfo[j].db = minfo[j].da = -1;
  minfo[j].uptime[0] =minfo[j].uptime[1] = -1;
  minfo[j].mtime[0] = minfo[j].mtime[1] = -1;
  minfo[j].mcount[0] = minfo[j].mcount[1] = minfo[j].mnew[0] =minfo[j].mnew[1] = -1;
  minfo[j].npopsb = minfo[j].npopsa  = -1;
  minfo[j].cm2_b[0] = minfo[j].cm2_b[1] = minfo[j].cm2_a[0] =minfo[j].cm2_a[1] = -1;
  minfo[j].mrate[0] = minfo[j].mrate[1] = minfo[j].mrate_r[0] = minfo[j].mrate_r[1] = -1;
  minfo[j].logpfpop = minfo[j].logpfpop_r = -1;
}

#define  MIGSIMFRAC  0.999  /* fraction of updated edges with a final population the same as a starting population,  
  in cases when there are two possible final populations */

double
update_mig_tNW (int li, struct edge *gtree, int period)
{
  double uptime;                //, tu, td;
  int i, j, k,  sis;
  double mproposedenom, mproposenum;
  int period_a, period_b, periodp1;   
  int numupdates, ei,mi;
  int mstart, checkup0, checkup1;
  int setfpop[2*MAXGENES-1]; // contains the position in the minfo[] array in which the update information for an edge is contained
  if (newt < oldt)
  {
    tu = newt;
    td = oldt;
    period_a = period + 1;
    period_b = period;
  }
  else
  {
    tu = oldt;
    td = newt;
    period_a = period;
    period_b = period + 1;
  }
  periodp1 = period+1;
  for (i = 0; i < L[li].numlines; i++)
    setfpop[i] = -1;

    /*********** Setup minfo - indicate which branches need attention  ***********/
  // read about struct migrationinfo_tNW at the top of this file
  for (i = 0, j=0; i < L[li].numlines; i++)
  {
    uptime = (i < L[li].numgenes) ? 0 : gtree[gtree[i].up[0]].time;
    // identify edges that need updating
    //  identify the population state of the bottoms of edges before and after the update
    if ((gtree[i].time > tu && uptime <= td) /*&& i != G->root */) // bottom falls in interval or edge passes through interval
    {
      if (setfpop[i] == -1)
      {
        initminfo(j);
        minfo[j].uptime[0] = uptime;
        /* set the population states at the bottom of the time interval,  db and da */
        /* set the pfpop (logpfpop and logpfpop_r)  terms */
        if (gtree[i].time > td)  // bottom of edge crosses lower time,  so only a single edge is done, fpop is unchanged
        {
          setfpop[i] = j;
          minfo[j].n = 1;
          minfo[j].edgeid[0] = i;
          minfo[j].db = nowedgepop (ci, &gtree[i], td);
          if (oldt < newt)
          {
            if (minfo[j].db == C[ci]->addpop[periodp1])
            {
 /*             minfo[j].logpfpop = -LOG2;
              minfo[j].logpfpop_r = 0.0;
              minfo[j].da = C[ci]->droppops[periodp1][bitran()];
              // see what upa is,   */

              minfo[j].logpfpop_r = 0.0;
              checkup0 = nowedgepop (ci, &gtree[i], tu);
              if (uptime < tu && (checkup0 == C[ci]->droppops[periodp1][0] || checkup0 == C[ci]->droppops[periodp1][1]))  // then we know what upa will be and can use it to simulate db;
              {
                if (uniform() < MIGSIMFRAC)
                {
                  minfo[j].da = checkup0;
                  minfo[j].logpfpop = log(MIGSIMFRAC);
                }
                else
                {
                   if (checkup0 == C[ci]->droppops[periodp1][0])
                    minfo[j].da = C[ci]->droppops[periodp1][1];
                  else
                    minfo[j].da = C[ci]->droppops[periodp1][0];
                  minfo[j].logpfpop = log(1.0 - MIGSIMFRAC);
                }
              }
              else
              {
                minfo[j].logpfpop = -LOG2;
                minfo[j].da = C[ci]->droppops[periodp1][bitran()];
              }
            }
            else
            {
              minfo[j].da = minfo[j].db;
              minfo[j].logpfpop = minfo[j].logpfpop_r =  0.0;
            }
          }
          else //oldt > newt
          {
            if (minfo[j].db == C[ci]->droppops[periodp1][0] ||minfo[j].db == C[ci]->droppops[periodp1][1])
            {
              minfo[j].logpfpop = 0.0;
              minfo[j].da = C[ci]->addpop[periodp1];
              checkup0 = nowedgepop (ci, &gtree[i], tu);
              if (uptime < tu && (checkup0 == C[ci]->droppops[periodp1][0] || checkup0 == C[ci]->droppops[periodp1][1]))  // then we know what upb will be and can use when imagining simulating db;
              {
                if (checkup0 == minfo[j].db)
                  minfo[j].logpfpop_r = log(MIGSIMFRAC);
                else
                  minfo[j].logpfpop_r = log(1.0 - MIGSIMFRAC);
              }
              else
              {
                minfo[j].logpfpop_r = -LOG2;
              }
            }
            else
            {
              minfo[j].da = minfo[j].db;
              minfo[j].logpfpop = minfo[j].logpfpop_r = 0.0;
            }
          }
        }
        else // gtree[i].time <=  td   two edges need to be done
        {
          setfpop[i] = j;
          if (i== gtree[gtree[i].down].up[0])
            sis = gtree[gtree[i].down].up[1];
          else
            sis = gtree[gtree[i].down].up[0];
          minfo[j].uptime[1] = (sis < L[li].numgenes) ? 0 : gtree[gtree[sis].up[0]].time;
          setfpop[sis] = j;
          minfo[j].n = 2;
          minfo[j].edgeid[0] = i;
          minfo[j].edgeid[1] = sis;
          minfo[j].db = nowedgepop (ci, &gtree[i], gtree[i].time);
          if (oldt < newt)
          {

            minfo[j].logpfpop_r = 0.0;
            if (minfo[j].db == C[ci]->addpop[periodp1])
            {
              checkup0 = nowedgepop (ci, &gtree[i], tu);
              checkup1 = nowedgepop (ci, &gtree[sis], tu);
              if (minfo[j].uptime[0] < tu && minfo[j].uptime[1] < tu && checkup0 == checkup1 &&
                (checkup0 == C[ci]->droppops[periodp1][0] || checkup0 == C[ci]->droppops[periodp1][1]))
              {
                if (uniform() < MIGSIMFRAC)
                {
                  minfo[j].da = checkup0;
                  minfo[j].logpfpop = log(MIGSIMFRAC)/2.0;
                }
                else
                {
                   if (checkup0 == C[ci]->droppops[periodp1][0])
                    minfo[j].da = C[ci]->droppops[periodp1][1];
                  else
                    minfo[j].da = C[ci]->droppops[periodp1][0];
                  minfo[j].logpfpop = log(1.0 - MIGSIMFRAC)/2.0;
                }

              }
              else
              {
                minfo[j].logpfpop = - LOG2HALF;
                minfo[j].da = C[ci]->droppops[periodp1][bitran()];
              }
            }
            else
            {
              minfo[j].da = minfo[j].db;
              minfo[j].logpfpop = 0.0;
            }

          }
          else  //newt < oldt 
          {
            minfo[j].logpfpop = 0.0;
            if (minfo[j].db == C[ci]->droppops[periodp1][0] ||minfo[j].db == C[ci]->droppops[periodp1][1])
            {
              checkup0 = nowedgepop (ci, &gtree[i], tu);
              checkup1 = nowedgepop (ci, &gtree[sis], tu);
              minfo[j].da = C[ci]->addpop[periodp1];
              if (minfo[j].uptime[0] < tu && minfo[j].uptime[1] < tu && checkup0 == checkup1 &&
                (checkup0 == C[ci]->droppops[periodp1][0] || checkup0 == C[ci]->droppops[periodp1][1]))
              {
                if (checkup0 == minfo[j].db)
                  minfo[j].logpfpop_r = log(MIGSIMFRAC)/2.0;
                else
                  minfo[j].logpfpop_r = log(1.0 - MIGSIMFRAC)/2.0;
              }
              else
              {
                minfo[j].logpfpop_r = -LOG2HALF;
              }
            }
            else
            {
              minfo[j].da = minfo[j].db;
              minfo[j].logpfpop_r = 0.0;
            }
          }
          gtree[i].fpop = gtree[sis].fpop = minfo[j].da;
        }
        j++;
      }
    }
  }
  numupdates = j;
  for (j=0;j<numupdates;j++)
  {
    // identify population states at the top of the edge within the interval, before and after the update
    for (i=0;i<minfo[j].n;i++)
    {
      ei = minfo[j].edgeid[i];

      if (newt < oldt)
      {
        if (minfo[j].uptime[i] < tu)
        {
          minfo[j].upb[i] = nowedgepop (ci, &gtree[ei], tu);
          if (minfo[j].upb[i] == C[ci]->droppops[periodp1][0] ||minfo[j].upb[i] == C[ci]->droppops[periodp1][1])
            minfo[j].upa[i] = C[ci]->addpop[periodp1];
          else
            minfo[j].upa[i] =  minfo[j].upb[i];
          
        }
        else  // uptime falls in the interval
        {
          minfo[j].upb[i] = minfo[setfpop[gtree[ei].up[0]]].db;
          assert(minfo[j].upb[i] == nowedgepop (ci, &gtree[ei], minfo[j].uptime[i]));
          minfo[j].upa[i] = minfo[setfpop[gtree[ei].up[0]]].da;
          gtree[ei].pop = minfo[j].upa[i];
        }
      }
      else // (newt > oldt)
      {
        if (minfo[j].uptime[i] < tu)
        {
          minfo[j].upb[i] = nowedgepop (ci, &gtree[ei], tu * (1 + DBL_EPSILON)); // just below tu 
          if (minfo[j].upb[i] == C[ci]->addpop[periodp1])
            minfo[j].upa[i] = nowedgepop (ci, &gtree[ei], tu); // just at tu
          else 
            minfo[j].upa[i] = minfo[j].upb[i];
        }
        else  // uptime falls in the interval
        {
          minfo[j].upb[i] = minfo[setfpop[gtree[ei].up[0]]].db;
          assert(minfo[j].upb[i] == nowedgepop (ci, &gtree[ei], minfo[j].uptime[i]));
          minfo[j].upa[i] = minfo[setfpop[gtree[ei].up[0]]].da;
          gtree[ei].pop = minfo[j].upa[i];
        }
      }
    }
  }

  /* determine migration counts and rates */
  for (j=0;j<numupdates;j++)
  {
    minfo[j].npopsb = npops - period_b;
    minfo[j].npopsa = npops - period_a;
    for (i=0;i<minfo[j].n;i++)
    {
      ei = minfo[j].edgeid[i];
      if (ei != G->root)
      {
        minfo[j].mtime[i] = DMIN (td, gtree[ei].time) - DMAX (tu, minfo[j].uptime[i]);
        assert(minfo[j].mtime[i] > 0);
        /* determine mcount */
        mi = 0;
        k = 0;
        mstart  = -1;
        while (gtree[ei].mig[k].mt > -0.5 && gtree[ei].mig[k].mt < td)
        {
          if (gtree[ei].mig[k].mt > tu)
          {
            if (mi==0)
              mstart = k;
            mi++;
          }
          k++;
        }
        if (k >= 2 && mstart >= 0 && (k- mstart >= 2))
        {
          if (k==2)
            minfo[j].cm2_b[i] = minfo[j].upb[i];
          else
            //minfo[j].cm2_b[i] = gtree[ei].mig[k-3].mp;
            minfo[j].cm2_b[i] = gtree[ei].mig[k-2].mp;
        }
        else
          minfo[j].cm2_b[i] = -1;
        assert(mi >= 0);
        minfo[j].mcount[i] = mi;
        minfo[j].cmm[i] = gtree[ei].cmm - mi;
        assert(period_b < lastperiodnumber || mi==0);
        /* set mrate */ 
        if (period_a < lastperiodnumber)
        {
          minfo[j].mrate[i] = calcmrate(minfo[j].mcount[i],minfo[j].mtime[i])*minfo[j].mtime[i];
        }
        else
        {
          minfo[j].mrate[i] = 0;
        }
        /* determine mnew */
        if (minfo[j].npopsa == 1)
          minfo[j].mnew[i] = 0;
        else
        {
          if (minfo[j].npopsa == 2)
          {
            if (minfo[j].upa[i] == minfo[j].da)
              minfo[j].mnew[i] = poisson (minfo[j].mrate[i],0,minfo[j].cmm[i]);
            else
              minfo[j].mnew[i] = poisson (minfo[j].mrate[i],1,minfo[j].cmm[i]);
          }
          else
          {
            if (minfo[j].upa[i] == minfo[j].da)
              minfo[j].mnew[i] = poisson (minfo[j].mrate[i],3,minfo[j].cmm[i]);
            else
              minfo[j].mnew[i] = poisson (minfo[j].mrate[i],2,minfo[j].cmm[i]);
          }
          if (minfo[j].mnew[i] < 0)
            return FORCEREJECTIONCONSTANT;
        }
        /* set the reverse update rate */
        if (period_b < lastperiodnumber)
        {
          minfo[j].mrate_r[i] = calcmrate(minfo[j].mnew[i],minfo[j].mtime[i])*minfo[j].mtime[i];
        }
        else
        {
          minfo[j].mrate_r[i] = 0.0;
        }
        mforward += minfo[j].mnew[i];
        mreverse += minfo[j].mcount[i];
      }
    }
  }
  /* simulate migration events */
  addmigration_NW(gtree, period_a, numupdates);
  /* calculate the Hastings term associated with the migration events before and after the update */
  mproposedenom =  mproposenum = 0;
  for (j=0;j<numupdates;j++)
  {
    for (i=0;i<minfo[j].n;i++)
    {
      assert(minfo[j].mrate[i]!=0 || minfo[j].logpfpop == 0);
      if (minfo[j].mrate[i] > 0)
      {
        mproposedenom += (minfo[j].logpfpop + 
          getmprob_NW(period_a,minfo[j].mrate[i],minfo[j].mtime[i],minfo[j].mnew[i],minfo[j].upa[i],minfo[j].da, minfo[j].cm2_a[i], minfo[j].npopsa) );
      }
      assert(minfo[j].mrate_r[i]!=0 || minfo[j].logpfpop_r == 0);
      if (minfo[j].mrate_r[i] > 0)
      {
        mproposenum += (minfo[j].logpfpop_r + 
          getmprob_NW(period_b,minfo[j].mrate_r[i],minfo[j].mtime[i],minfo[j].mcount[i],minfo[j].upb[i],minfo[j].db, minfo[j].cm2_b[i], minfo[j].npopsb) );
      }
    }
  }
   return mproposenum - mproposedenom; 
}                               /* update_mig_tNW */

/**************************************/
/*  GLOBAL FUNCTIONS in update_t_NW.c */
/**************************************/

void
init_t_NW (void)
{
  int li, j;
  init_genealogy_weights (&holdallgweight_t_NW);

  /* CR:110331.1 
   * Incorrect structure used when allocating memory to a genealogy pointer.
   */
  holdallgtree_t_NW = static_cast<struct genealogy *> 
                      (calloc ((size_t) nloci, sizeof (struct genealogy)));
  for (li = 0; li < nloci; li++)
  {
    init_genealogy_weights (&holdgweight_t_NW[li]);
    init_holdgtree (&holdallgtree_t_NW[li], L[li].numgenes);
  }
  init_probcalc (&holdallpcalc_t_NW);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  minfo = static_cast<struct migrationinfo_tNW *> 
          (malloc ((2 * largestsamp - 1) * sizeof (struct migrationinfo_tNW)));
}                               /* init_t_NW */

void
free_t_NW (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_NW);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_NW[li]);
    free_holdgtree (&holdallgtree_t_NW[li], L[li].numgenes);
  }
  XFREE (holdallgtree_t_NW);
  free_probcalc (&holdallpcalc_t_NW);
  XFREE (minfo);
}                               /* init_t_NW */

/* changet_NW
   this is modelled on the update of t in Nielsen and Wakeley (2001) 
   does nothing whatsover to branch lengths  so no change in pdg
   
   This update applies to all genealogies in the chain 
   */

/*  The update is done to the actual genealogies
    Could improve speed by doing the update to a copy of the genealogies, 
    This is because most updates are rejected, so if the copy is being changed and it is rejected 
    then there would be no need to restore it 

    to make this change would have to change a few lines of code (below)
    also have to pass treeweight a genealogy pointer and not just the genealogy number
    there will no doubt be other stuff
    spent a couple hours trying this on 8/28/08,  but wasn't working so gave up.  maybe try later

        
        Current Setup -  work on the original:
    --------------------------------------

        copy_all_gtree(1);  -  copy each C[ci]->G  to holdallgtree_t_NW
        copy_treeinfo (&holdallgweight_t_NW, &C[ci]->allgweight);  - copy allgweight to holdallgweight_t_NW - then empty allgweight
        copy_probcalc (&holdallpcalc_t_NW, &C[ci]->allpcalc);  - copy allpcalc to holdallpcalc_t_NW  (stuff for integrating)
        setzero_genealogy_weights (&C[ci]->allgweight); - empty allgweight

        
        going thru the loop:
            copy_treeinfo (&holdgweight_t_NW[li], &G->gweight);  - copy gweight to holdgweight_t_NW for each locus
            update C[ci]->G
            setzero_genealogy_weights (&G->gweight);   -  set gweight to zero
            treeweight (ci, li);   -   calculate weights
            sum_treeinfo (&C[ci]->allgweight, &G->gweight);  -     sum allgweight
        
        integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc); - integrate and reset allpcalc 

        if accept, 
            everything is already in gweight, allgweight and C[ci]->G 
        if reject:
            copy_all_gtree(0);  -  copy holdallgtree_t_NW back to each C[ci]->G
            copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);  -  copy holdallgweight_t_NW back to allgweight
            copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);   -  copy holdallpcalc_t_NW back to allpcalc
            for (li = 0; li < nloci; li++)   - loop thru loci
            {
                copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);  - copy each locus's from holdgweight_t_NW  back to gweight
            }

       If we reverse it :
       ------------------------------------
        copy_all_gtree(1);  -  copy each C[ci]->G  to holdallgtree_t_NW
        setzero_genealogy_weights (&holdallgweight_t_NW); - empty holdallgweight_t_NW

        going thru the loop:
            update holdgallgtree[li]
            setzero_genealogy_weights (&holdgweight_t_NW[li]);   -  set holdgweight_t_NW to zero
            local_treeweight (ci, li, holdallgtree_t_NW[li]);   -   calculate weights
            sum_treeinfo (&holdallgweight_t_NW, &holdgweight_t_NW[li]);  -     sum allgweight

        integrate_tree_prob (ci, &holdallgweight_t_NW, &holdallpcalc_t_NW); - integrate and reset allpcalc 

        if accept, 
            copy_all_gtree(0);  -  copy holdallgtree_t_NW to each C[ci]->G
            copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);  -  copy holdallgweight_t_NW into allgweight
            copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);   -  copy holdallpcalc_t_NW to allpcalc
            for (li = 0; li < nloci; li++)   - loop thru loci
            {
                copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);  - copy each locus's from holdgweight_t_NW  back to gweight
            }

        if reject:
            don't do anything because nothing in C[ci]->G, allgweight or allpcalc was changed.
 */

/* let u refer to the more recent time  and d to the older time  */

/* 9/25/08  updated this
revised genealogy updating so that the current migration rate is based on the current number of migration events and the 
current length of the branch that is being updated 
reasoned that this might work better than using the rate that occurs for the entitre genealogy - e.g. help to avoid promoting
correlations and improve mixing  */
/* only works for nonzero migration priors */

int
changet_NW (int chain, int timeperiod)
{

  double metropolishastingsratio;  
  double likelihoodratio = 0.0,priorratio,proposalratio;
  int li;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double migweight;

  checkupdatescalarer(&C[ci]->NWwidthinfo[timeperiod]);
  mforward = mreverse = 0;
  ci = chain;
  /* select a new time */
  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d = (timeperiod == (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];
  //newt = getnewt (t_u_prior, t_d_prior, oldt, 0);
  newt = getnewt (t_u_prior, t_d_prior, oldt,C[ci]->NWwidthinfo[timeperiod].updatescalarval);
  assert (newt < T[timeperiod].pr.max);
  /* store stuff and prepare for adding to storing allgweight */

  copy_all_gtree (1);
  copy_treeinfo (&holdallgweight_t_NW, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc_t_NW, &C[ci]->allpcalc);
  setzero_genealogy_weights (&C[ci]->allgweight);

  /* initialize migration hastings term and loop thru the loci */
  migweight = 0;
  int tw = 0, ttw;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_NW[li], &G->gweight);
    /* if the root of the genealogy is younger than oldt and newt no updating of this genealogy is needed */
    if ((newt > oldt && G->roottime > oldt)
        || (newt < oldt && G->roottime > newt))
    {
      gtree = G->gtree;
      migweight += update_mig_tNW (li, gtree, timeperiod);
      C[ci]->tvals[timeperiod] = newt;  // reset for treeweight() calculations
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = \
        C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;
      setzero_genealogy_weights (&G->gweight);
#ifdef TURNONCHECKS
      checkgenealogy(0,0,4);
#endif //TURNONCHECKS

      int treeweightcallcode = 4;  // a debugging code,  if treeweight has an error results are written to an output file with this code 
      ttw = treeweight (ci, li,treeweightcallcode);
      if (ttw < 0)
        tw = -1;
     
      C[ci]->tvals[timeperiod] = oldt;  // put old value back for now because update_mig_tNW() depends on old value 
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = \
        C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = oldt;
    }
    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
  }
  if (tw < 0)
     goto rejectjump;  // don't like to do this kind of jump,  but need to reject if treeweight() fails 


  /* calculate the prior for the new genealogy - needed even if genealogy not changed */
  C[ci]->tvals[timeperiod] = newt;      // reset for integrate calculations
  C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time = \
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;
//initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc); 
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_NW,&C[ci]->allpcalc, &holdallpcalc_t_NW);

  /* calculate the MH term - depends ratio of prior probabilities and hastings term for simulated migration events */
  /* does not depend on P(D|G)  because branch lengths are not changed 
  so likelihoodratio = 0.0 */
  priorratio = C[ci]->allpcalc.probg - holdallpcalc_t_NW.probg;
  proposalratio = migweight;
  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  /* the likelihood does not appear in this update so beta is not used when calculating the marginallikelihood */

  if (calcoptions[CALCMARGINALLIKELIHOOD]) 
  {
    metropolishastingsratio = priorratio + proposalratio;
  }
  else
  {
      metropolishastingsratio = beta[ci] * priorratio + proposalratio;
  }
  if (metropolishastingsdecide(metropolishastingsratio,1))
  {
    C[ci]->NWwidthinfo[timeperiod].recentaccp += 1;
    C[ci]->NWwidthinfo[timeperiod].allaccp += 1;
#ifdef TURNONCHECKS
    for (li = 0; li < nloci; li++)
    //checkprobs(ci,li);
#endif //TURNONCHECKS
    return 1;
  }
  else
  {
    /* copy stuff back from where it was held */
rejectjump:
    C[ci]->tvals[timeperiod] = oldt;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = oldt;
    copy_all_gtree (0);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_NW);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_NW);
  
//initialize_integrate_tree_prob (ci,&C[ci]->allgweight,&C[ci]->allpcalc);

    for (li = 0; li < nloci; li++)
    {
      copy_treeinfo (&C[ci]->G[li].gweight, &holdgweight_t_NW[li]);
      //      assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);   
#ifdef TURNONCHECKS
  //checkprobs(ci,li);
#endif //TURNONCHECKS
    }

    return 0;
  }
}                               /* changet_NW  - after Nielsen and Wakeley (2001) */
