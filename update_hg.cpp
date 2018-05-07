/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

// 3/31/2016 version of this file.   Had a series of versions with various bugs. this one is best so far
// rnamed this to just update_hg.cpp  (i.e. without the date in the name

#undef GLOBVARS
#include "ima.hpp"

#include "update_gtree_common.hpp"
#include "update_gtree.hpp"

extern struct genealogy_weights holdgweight_updategenealogy;
extern struct genealogy_weights holdallgweight_updategenealogy;
extern struct probcalc holdallpcalc_updategenealogy;

struct edge *copyedgehg; // holds copies of genealogy edges during update in case update is rejected
static double holdsisdlikeA[MAXLINKED];

/*********** local to this file  ***************/

extern int rootmove;  //declared in  update_gtree.cpp

static int holddownA[MAXLINKED];
static double forcereject = -1e100;
#define SLIDESTDVMAX 20   // not sure what's best, copied from update_ptree.cpp
#define  MIGCLOSEFRAC  0.9  // when both edges updating, chance of zero or 1 migration event in last migration period, copied from update_gtree_common.cpp

struct hgcalcstruct{
  int nummighg;  // number of hidden migrations in the segment.
  int ndpops; // number of descendant pops to the current population the edge is in;
  double l;  // time span of interval
  double mrate;  // instantaneous migration rate
  double logprobendpop;  // probability of a particular target population if interval ends in a regular migration event
};

static struct hgcalcstruct hgczero = {0,0,0.0,0.0,0.0}; // for initialization

/******* local prototypes ********/

static void edgemask(int ci, struct edge *e,double upt);
#ifdef TURNONCHECKS
static void edgemask_debugprint(int ci, struct edge *e, double upt);
#endif
static void IMA_initmemory_edgemiginfohg (struct edgemiginfo *em);
static void storegenealogystatshg (int ci, int li, int mode);
void  init_gtreecommonhg (void);         // initialize copyedgehg
void free_gtreecommonhg (void);
static void fillmiginfohg (int ci, int li, struct edge *gtree, int edge, int sisedge);
static void storeoldedgeshg (int ci, int li, int edge, int sisedge, int downedge);
static void restoreedgeshg (int ci, int li, int edge, int sisedge, int downedge,  int newsisedge);
static double addmigrationhg (int ci,int  li);

static int mwork_single_edgehg (int ci, struct edgemiginfo *edgem,struct edgemiginfo *oldedgem, double mhg);
static int mwork_two_edgeshg(int ci, struct edgemiginfo *edgem, struct edgemiginfo *sisem,  struct edgemiginfo *oldedgem,  struct edgemiginfo *oldsisem,double mhg);
static int getmhg (int ci, struct edgemiginfo *edgem,struct edgemiginfo *sisem, struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem, double mhg);
double getmprobedge(struct edgemiginfo *e, double mhg);
double getmprobhg(int ci, struct edgemiginfo *edgem,struct edgemiginfo *sisem,struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem,double mhg);
int simmpathhg (int ci, struct edgemiginfo *edgem, int numm, double timein, double upt, int *pop, int constrainpop);
double hgedgecalcprob(struct hgcalcstruct hgc);
double hgcalccoalprob(int nm[], int cpop[], int thirdlastpop[], double l[], int fpop, double mrate,int ndpops);
int countnummighg(int ci,int li,int ei,double t, double endt,int *lastpop, int *thirdlastpop);
/******* local functions ********/


/* store a few basic tree statistics, just a copy of storegenealogystats() */
void
storegenealogystatshg (int ci, int li, int mode)
{
  static double holdlength, holdtlength;
  static double holdroottime;
  static int holdroot;
  static int holdmig;
  static double holdhgprob;
  if (mode == 0)
  {
    holdlength = C[ci]->G[li].length;
    holdtlength = C[ci]->G[li].tlength;
    holdroottime = C[ci]->G[li].roottime;
    holdroot = C[ci]->G[li].root;
    holdmig = C[ci]->G[li].mignum;
    holdhgprob =  C[ci]->G[li].hgprob;
  }
  else
  {
    C[ci]->G[li].length = holdlength;
    C[ci]->G[li].tlength = holdtlength;
    C[ci]->G[li].mignum = holdmig;
    C[ci]->G[li].roottime = holdroottime;
    C[ci]->G[li].root = holdroot;
    C[ci]->G[li].hgprob = holdhgprob;
  }
}   /* storegenealogystatshg */

void
storeoldedgeshg (int ci, int li, int edge, int sisedge, int downedge)
{
  int i,j,ai;
  int a[3];
  struct edge *gtree = C[ci]->G[li].gtree;
  a[0] = edge;
  a[1] = sisedge;
  a[2] = downedge;

  for (j = 0;j<3;j++)
  {
    copyedgehg[j].down = gtree[a[j]].down;
    copyedgehg[j].cmm = gtree[a[j]].cmm;
    copyedgehg[j].time = gtree[a[j]].time;
    copyedgehg[j].pop = gtree[a[j]].pop;
    copyedgehg[j].fpop = gtree[a[j]].fpop;
    copyedgehg[j].up[0] = gtree[a[j]].up[0];
    copyedgehg[j].up[1] = gtree[a[j]].up[1];
    copyedgehg[j].cmmhg = gtree[a[j]].cmmhg;
    copyedgehg[j].pophg = gtree[a[j]].pophg;
    copyedgehg[j].fpophg = gtree[a[j]].fpophg;
    i = -1;
    do
    {
      i++;
      if (i > MIGMAX)
      {
        IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] a[j] [%d] mig %d > %d",
                step, li, a[j], MIGMAX, i);
      }
      copyedgehg[j].mig[i] = gtree[a[j]].mig[i];
    } while (copyedgehg[j].mig[i].mt > -0.5);
    copyedgehg[j].cmm = i;
    i = -1;
    do
    {
      i++;
      if (i > MIGMAX)
      {
        IM_err (IMERR_TOOMANYMIG, "step %d: locus [%d] a[j] [%d] mig %d > %d",
                step, li, a[j], MIGMAX, i);
      }
      copyedgehg[j].mighg[i] = gtree[a[j]].mighg[i];
    } while (copyedgehg[j].mighg[i].mt > -0.5);
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
  {
    for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
    {
      copyedgehg[0].A[ai] = gtree[edge].A[ai];
      copyedgehg[0].dlikeA[ai] = gtree[edge].dlikeA[ai];
      copyedgehg[1].A[ai] = gtree[sisedge].A[ai];
      copyedgehg[1].dlikeA[ai] = gtree[sisedge].dlikeA[ai];
      copyedgehg[2].A[ai] = gtree[downedge].A[ai];
      if (gtree[downedge].down != -1)   // added this 5/27/08
      {
        copyedgehg[2].dlikeA[ai] = gtree[downedge].dlikeA[ai];
        holddownA[ai] = gtree[gtree[downedge].down].A[ai];      // changed this 5/27/08
      }
      else
      {
        copyedgehg[2].down = -1;  // inserted this 5/27/08
        holddownA[ai] = -1;
        copyedgehg[2].dlikeA[ai] = 0;
      }
    }
  }
}                               /* storoldegeshg */

/* set the gtree back to the way it was */
void
restoreedgeshg (int ci, int li, int edge, int sisedge, int downedge,  int newsisedge)
/*all this can be optimized some*/
{
  int i, j, ai, down;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (newsisedge != sisedge)
  {
    down = gtree[downedge].down;
    if (down != -1)
    {
      if (gtree[down].up[0] == downedge)
        gtree[down].up[0] = newsisedge;
      else
        gtree[down].up[1] = newsisedge;
    }
    else
    {
      C[ci]->G[li].root = newsisedge;
      C[ci]->G[li].roottime = gtree[gtree[newsisedge].up[0]].time;
      assert (C[ci]->G[li].roottime <= TIMEMAX);
    }
    gtree[newsisedge].down = down;
    if (down != -1)
    {
      i = 0;
      while (gtree[newsisedge].mig[i].mt > -0.5)
        i++;
      j = -1;
      do
      {
        j++;
        gtree[newsisedge].mig[i + j] = gtree[downedge].mig[j];
      } while (gtree[downedge].mig[j].mt > -0.5);
      gtree[newsisedge].cmm = i+j;
      i = 0;
      while (gtree[newsisedge].mighg[i].mt > -0.5)
        i++;
      j = -1;
      do
      {
        j++;
        gtree[newsisedge].mighg[i + j] = gtree[downedge].mighg[j];
      } while (gtree[downedge].mighg[j].mt > -0.5);
      gtree[newsisedge].cmmhg = i + j;
      gtree[newsisedge].fpop = C[ci]->ancplist[gtree[newsisedge].fpop][findperiod(ci,gtree[downedge].time)]; // do we need this ?
      gtree[newsisedge].fpophg = gtree[downedge].fpophg;
    }
    else
    {
      gtree[newsisedge].mig[0].mt = -1;
      gtree[newsisedge].mighg[0].mt = -1;
      gtree[newsisedge].cmm = gtree[newsisedge].cmmhg = 0;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
          gtree[newsisedge].dlikeA[ai] = 0;
    }
    gtree[newsisedge].time = gtree[downedge].time;
  }
  gtree[edge].down = copyedgehg[0].down;
  i = -1;
  do
  {
    i++;
    gtree[edge].mig[i] = copyedgehg[0].mig[i];
  } while (gtree[edge].mig[i].mt > -0.5);
  gtree[edge].cmm = i;
  i = -1;
  do
  {
    i++;
    gtree[edge].mighg[i] = copyedgehg[0].mighg[i];
  } while (gtree[edge].mighg[i].mt > -0.5);
  gtree[edge].cmmhg = i;
  gtree[edge].pophg = copyedgehg[0].pophg;
  gtree[edge].fpophg = copyedgehg[0].fpophg;
  gtree[edge].time = copyedgehg[0].time;
  gtree[edge].pop = copyedgehg[0].pop;
  gtree[edge].fpop = copyedgehg[0].fpop;
  down = gtree[sisedge].down;
  gtree[sisedge].down = copyedgehg[1].down;
  if (down != -1)
  {
    if (gtree[down].up[0] == sisedge)
      gtree[down].up[0] = downedge;

    else
      gtree[down].up[1] = downedge;
  }
  i = -1;
  do
  {
    i++;
    gtree[sisedge].mig[i] = copyedgehg[1].mig[i];
  } while (gtree[sisedge].mig[i].mt > -0.5);
  gtree[sisedge].cmm = i;
  i = -1;
  do
  {
    i++;
    gtree[sisedge].mighg[i] = copyedgehg[1].mighg[i];
  } while (gtree[sisedge].mighg[i].mt > -0.5);
  gtree[sisedge].cmmhg = i;
  gtree[sisedge].pophg = copyedgehg[1].pophg;
  gtree[sisedge].fpophg = copyedgehg[1].fpophg;
  gtree[sisedge].time = copyedgehg[1].time;
  gtree[sisedge].pop = copyedgehg[1].pop;
  gtree[sisedge].fpop = copyedgehg[1].fpop;
  gtree[downedge].down = copyedgehg[2].down;
  i = -1;
  do
  {
    i++;
    gtree[downedge].mig[i] = copyedgehg[2].mig[i];
  } while (gtree[downedge].mig[i].mt > -0.5);
  gtree[downedge].cmm = i;
  i = -1;
  do
  {
    i++;
    gtree[downedge].mighg[i] = copyedgehg[2].mighg[i];
  } while (gtree[downedge].mighg[i].mt > -0.5);
  gtree[downedge].cmmhg = i;
  gtree[downedge].pophg = copyedgehg[2].pophg;
  gtree[downedge].fpophg = copyedgehg[2].fpophg;
  gtree[downedge].time = copyedgehg[2].time;
  gtree[downedge].pop = copyedgehg[2].pop;
  gtree[downedge].fpop = copyedgehg[2].fpop;
  gtree[downedge].up[0] = copyedgehg[2].up[0];
  gtree[downedge].up[1] = copyedgehg[2].up[1];
  if (gtree[downedge].down == -1)
  {
    C[ci]->G[li].roottime = gtree[gtree[downedge].up[0]].time;
    assert (C[ci]->G[li].roottime <= TIMEMAX);
    C[ci]->G[li].root = downedge;
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    for (ai = (L[li].model == JOINT_IS_SW); ai < L[li].nlinked; ai++)
    {
      gtree[edge].A[ai] = copyedgehg[0].A[ai];
      gtree[sisedge].A[ai] = copyedgehg[1].A[ai];
      gtree[downedge].A[ai] = copyedgehg[2].A[ai];
      gtree[edge].dlikeA[ai] = copyedgehg[0].dlikeA[ai];
      gtree[sisedge].dlikeA[ai] = copyedgehg[1].dlikeA[ai];
      gtree[downedge].dlikeA[ai] = copyedgehg[2].dlikeA[ai];
      if (holdsisdlikeA[ai] != 0)
        gtree[newsisedge].dlikeA[ai] = holdsisdlikeA[ai];
    }
}                               /* restoreedgeshg  */

void
init_gtreecommonhg (void)         // initialize copyedgehg
{
  int i;
  copyedgehg = static_cast<struct edge *> (malloc (3 * (sizeof (struct edge))));
  if (somestepwise)
    for (i = 0; i < 3; i++)
    {
      copyedgehg[i].A = static_cast<int *> (calloc (MAXLINKED, sizeof (int)));
      copyedgehg[i].dlikeA = static_cast<double *>
                            (calloc (MAXLINKED, sizeof (double)));
    }
}                               /* init_gtreecommonhg */


void free_gtreecommonhg (void)
{
  int i;
  for (i = 0; i < 3; i++)
  {
    if (somestepwise)
    {
      XFREE (copyedgehg[i].A);
      XFREE (copyedgehg[i].dlikeA);
    }
  }
  XFREE (copyedgehg);
  return;
}                               //free_gtreecommonhg


/* copy info into oldedgemig and oldsismig,  used for calculating probability of updating the genealogy
these are instances of edgemiginfo.  They apply to the hidden genealogy, but the code was copied from work on the genealogy, and so there is
nothing here that refers explicitly to the hidden genealogy.
oldsismig is only used when sisedge != -1,  which happens when oldedgemig connects to the root.
*/
void
fillmiginfohg (int ci, int li, struct edge *gtree, int edge, int sisedge)
{
  int i, k;
  double uptime;
  struct edgemiginfo *ep[2];  // short array that holds two pointers
  int a[2];

  ep[0] = &oldedgemig;
  ep[1] = &oldsismig;
  a[0] = edge;
  a[1] = sisedge;

  for (k=0;k< 2;k++)
  {
    ep[k]->edgeid = a[k];
    if (a[k] >= 0) // sisedge (aj[1]) can be -1
    {
      IMA_reset_edgemiginfo (ep[k]);   // no hg version of this
      ep[k]->edgeid = a[k];
      ep[k]->li = li;
      if (a[k] < L[li].numgenes)
        uptime = 0;
      else
        uptime = gtree[gtree[a[k]].up[0]].time;
      ep[k]->dnt = gtree[a[k]].time;
      ep[k]->upt = uptime;
      ep[k]->b = 0;
      while (ep[k]->upt > C[ci]->tvals[ep[k]->b])
        ep[k]->b++;
      ep[k]->e = ep[k]->b;
      while (ep[k]->dnt > C[ci]->tvals[ep[k]->e])
        ep[k]->e++;
      ep[k]->mtall = ep[k]->dnt - ep[k]->upt;
      assert(ep[k]->mtall > 0);
      ep[k]->pop = gtree[a[k]].pophg;
      ep[k]->fpop = gtree[gtree[a[k]].down].pophg;
      assert (ep[k]->mig != NULL);
      i = 0;
      while (gtree[a[k]].mighg[i].mt > -0.5)
      {
        ep[k]->mig[i] = gtree[a[k]].mighg[i];
        i++;
      }
      ep[k]->mig[i].mt = -1;
      ep[k]->mpall = i;
      ep[k]->cmm = i;
    }
    else
    {
      ep[k]->mtall = 0.0;
      ep[k]->mpall = 0;
      ep[k]->mig[0].mt = -1;
      ep[k]->cmm = 0;
    }
  }
}                               /* fillmiginfohg */

/* simulate the migration path along the moved edge
unlike simmpath(),  no use of time periods in this code
numm is the number of migration events to add
if it must end up in a particular population,  then constrainpop is that population */
int
simmpathhg (int ci, struct edgemiginfo *edgem, int numm, double timein, double upt, int *pop, int constrainpop)
{
  int i, lastpop,lastm;
                    /* CR 110715.1 */
  int dupCheck;     /* flag turns off dup migration time check of mig events */
  int migIndex;     /* index used to look for duplicate migration times */
  assert (numm > 0);
  lastm = numm-1; // position in mig array of the last migration event to be added
  do
  {
    for (i = 0; i <= lastm; i++)
        edgem->mig[i].mt = upt + uniform () * timein;
    edgem->mig[i].mt = -1;
    dupCheck=0;
    if (numm > 1)
    {
      hpsortmig (&edgem->mig[0] - 1, numm);

      /* CR 110715.1
       * look for duplicate migration times in sorted event list.
       * This solves a charateristic of the Mersennes Twister
       * random number generator in which identical random numbers may be
       * returned from the random number sequence in a very small number
       * of calls.  With some seeds it was noted as small as within 4 calls.
       */
      for (migIndex = 0; migIndex < lastm; ++migIndex)
      {
        if  (edgem->mig[migIndex].mt != edgem->mig[migIndex + 1].mt)
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
  } while (dupCheck == 1); /* when no dup times found, exit loop   */

  lastpop = *pop;
  for (i = 0; i <= lastm; i++)
  {
    if (constrainpop >= 0 && i >= lastm-1)
    {
      if (i==lastm-1)
        edgem->mig[i].mp =  picktopop2 (lastpop, C[ci]->plist[0], npops,constrainpop); // pick a population other than constrainpop
      else //i==lastm
        edgem->mig[i].mp = constrainpop;
    }
    else
      edgem->mig[i].mp = picktopop (lastpop, C[ci]->plist[0], npops);
    lastpop = edgem->mig[i].mp;
  }
  *pop = lastpop;
  return lastm;
} // simmpathhg

/*
  simulates migration times for single edges,  updates temppop (population id of the edge at the time under consideration) as needed
*/
int
mwork_single_edgehg (int ci, struct edgemiginfo *edgem,struct edgemiginfo *oldedgem,double mhg)
{
  int mpall,lastm;
  double r;

  edgem->temppop = edgem->pop;
  r = mhg *edgem->mtall;
  assert(edgem->cmm == 0);
  if (npops == 2)
  {
    if (edgem->temppop == edgem->fpop)      //even
    {
      mpall = poisson (r, 0,edgem->cmm);
    }
    else                    // odd
    {
      mpall = poisson (r, 1,edgem->cmm);
    }
    edgem->cmm = mpall;
    if (mpall > 0)
    {
      lastm = simmpathhg (ci, edgem, mpall,edgem->mtall,edgem->upt, &edgem->temppop, -1);
    }
    else
      lastm = -1;
  }
  else
  {
    if (edgem->temppop == edgem->fpop)      //cannot be just one migration event
    {
      mpall = poisson (r, 3,edgem->cmm);
    }
    else                    // cannot be zero migration events
    {
      mpall = poisson (r, 2,edgem->cmm);
    }
    edgem->cmm = mpall;
    if (mpall < 0)
      return -1;
    if (mpall > 0)
    {
      lastm = simmpathhg (ci, edgem, mpall,edgem->mtall, edgem->upt,&edgem->temppop, edgem->fpop);
    }
    else
      lastm = -1;
  }
  assert(lastm < 0 || edgem->temppop == edgem->mig[lastm].mp);
  return mpall;
}                               /* mwork_single_edgehg */

/* return integer is just for reporting two many migrations */

int mwork_two_edgeshg(int ci, struct edgemiginfo *edgem, struct edgemiginfo *sisem,  struct edgemiginfo *oldedgem,  struct edgemiginfo *oldsisem,double mhg)
{
  int lastm[2];
  int ii;
  double r;
  struct edgemiginfo *mm;//, *oldmm; not used

  assert (edgem->e == sisem->e);
  assert(edgem->e >= edgem->b);
  assert(sisem->e >= sisem->b);
  edgem->temppop = edgem->pop; // should already be set to this
  sisem->temppop = sisem->pop; // should already be set to this
/* do constrained migration to each of two portions of both edges that are in lastmigperiod */
  if (edgem->temppop == sisem->temppop)
  {
    if (uniform() < MIGCLOSEFRAC)
    {
      edgem->fpop = sisem->fpop = edgem->temppop;
    }
    else
    {
      edgem->fpop = sisem->fpop = picktopop(edgem->temppop,C[ci]->plist[0],npops);
    }
  }
  else
  {
    if (npops == 2)
    {
      if (uniform() < 0.5)
        edgem->fpop = sisem->fpop = edgem->temppop;
      else
        edgem->fpop = sisem->fpop = sisem->temppop;
    }
    else
    {
      if (uniform() < MIGCLOSEFRAC)
      {
        if (uniform() < 0.5)
          edgem->fpop = sisem->fpop = edgem->temppop;
        else
          edgem->fpop = sisem->fpop = sisem->temppop;
      }
      else // pick something besides edgem->temppop and sisem->temppop
      {
        edgem->fpop = sisem->fpop = picktopop2(edgem->temppop,C[ci]->plist[0],npops,sisem->temppop);
      }
    }
  }
  for (ii=0;ii<2;ii++)
  {
    mm = (ii==0) ? edgem : sisem;
    //oldmm = (ii==0) ? oldedgem : oldsisem; not used
    r = mhg*mm->mtall;
    if (npops == 2)
    {
      if (mm->temppop == mm->fpop)      //even
      {
        mm->mpall = poisson (r, 0,mm->cmm);
        if (mm->mpall < 0)
          return -1;
        mm->cmm = mm->mpall;
      }
      else                    // odd
      {
        mm->mpall = poisson (r, 1,mm->cmm);
        if (mm->mpall < 0)
          return -1;
        mm->cmm = mm->mpall;
      }
      if (mm->mpall > 0)
      {
        lastm[ii] =
          simmpathhg (ci, mm, mm->mpall,mm->mtall, mm->upt, &mm->temppop, -1);  //lastm is not used for anything ?
      }
    }
    else
    {
      if (mm->temppop == mm->fpop)      //cannot be just one migration event
      {
        mm->mpall = poisson (r, 3,mm->cmm);
        if (mm->mpall < 0)
          return -1;
        mm->cmm = mm->mpall;
      }
      else                    // cannot be zero migration events
      {
        mm->mpall = poisson (r, 2,mm->cmm);
        if (mm->mpall < 0)
          return -1;
        mm->cmm += mm->mpall;
      }
      if (mm->mpall > 0)
      {
        lastm[ii] =
          simmpathhg (ci, mm, mm->mpall,mm->mtall,mm->upt,&mm->temppop, mm->fpop); // lastm is not used for anything
      }
    }
  }
  return 0;
}                               /* mwork_two_edgeshg*/


/* return integer is in case update should be rejected because of two much migration */
int  getmhg (int ci,struct edgemiginfo *edgem,struct edgemiginfo *sisem, struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem, double mhg)
{
  if ( sisem->edgeid == -1  /* sisem->mtall <= 0*/)  // no sister edge, or sister edge not in a period where migration can occur
  {
    assert(sisem->mtall == 0);
    edgem->mpall = mwork_single_edgehg (ci, edgem,oldsisem,mhg);
    if (edgem->mpall < 0)
      return -1;
  }
  else
  {
    assert (edgem->e == sisem->e);
  // both edge and sis have length in periods with migration
    if (mwork_two_edgeshg(ci, edgem, sisem,oldedgem,oldsisem,mhg) <0)
      return -1;
   }
  if (edgem->mpall > MIGMAX)
    return -1;
  else
    return 0;
}  //getmhg

/* copy gtree info from genealogy into newedgemig and newsismig.
and then simulate migration events on those.
have to fill up newdgemig and newsismig.  These are instances of edgemiginfo,  and as such are not unique to hidden genealogy work
use pop,fpop and mig[] parts of these to hold info about the hidden genealogy
*/
static double addmigrationhg (int ci,int  li)
{
  int newsis, edge;
  double weight;
  double mproposenum, mproposedenom;
  struct edge *gtree = C[ci]->G[li].gtree;
  int nm;
  double mhgproposeratio;
  double U, windowsize, oldmhg, newmhg,maxval,troldmhg,trnewmhg;

  /* code for updating the mhg value  */
  if (hiddenoptions[UPDATEMRATEFORHGUPDATE] && uniform() < 0.2) // do a migration prior update about 20% of the time
  {
    if (modeloptions[EXPOMIGRATIONPRIOR])
    {

     // C[ci]->G[li].mhg = expo(mprior);
      oldmhg =  C[ci]->G[li].mhg;
      troldmhg =  1.0 - exp(-C[ci]->G[li].mhg/mprior);  // transformed old mhg
      assert (0.0 < troldmhg  && troldmhg < 1.0);

      // make a value transformed to uniform[0,1] tra = 1- Exp[-a /m]
      // tra now goes from 0 to 1
      // get a new tra in a window on an interval of [0,1] with reflecting boundaries
      // back transform newa =  -Log[1- newtra] * m

      U = uniform();
      //maxval =  1.0;
      windowsize = DMIN(0.1,0.1/mprior);  // the bigger the mean of the exponential distribution the smaller the window size


      if (U > 0.5)
          {
          trnewmhg = troldmhg + (2.0*U)*windowsize;
          if (trnewmhg > 1.0)
              trnewmhg = 2.0 - trnewmhg;
          }
      else
          {
          trnewmhg = troldmhg - windowsize*U*2.0;
          if (trnewmhg < 0.0)
              trnewmhg = - trnewmhg;
          }
      C[ci]->G[li].mhg = - log(1 - trnewmhg) * mprior;  // back transform
      assert(C[ci]->G[li].mhg  > 0.0);

     // holdmhg = C[ci]->G[li].mhg;
     // C[ci]->G[li].mhg = expo(mprior);
      mhgproposeratio = (oldmhg - C[ci]->G[li].mhg)/ mprior;
    }
    else
    {
      //C[ci]->G[li].mhg = uniform() * mprior;
      oldmhg = C[ci]->G[li].mhg;
      U = uniform();
      maxval =  mprior;
      windowsize = 0.1; // fixed small window size

      if (U > 0.5)
        {
        newmhg = oldmhg + (2.0*U-1.0)*windowsize;
        if (newmhg > maxval)
            newmhg = 2.0*maxval - newmhg;
        }
      else
        {
        newmhg = oldmhg - windowsize*U*2.0;
        if (newmhg < 0)
            newmhg = - newmhg;
        }
      C[ci]->G[li].mhg = newmhg;
      mhgproposeratio  = 0.0;
    }
  }
  else
  {
    mhgproposeratio  = 0.0;
    oldmhg = C[ci]->G[li].mhg;
  }

  IMA_reset_edgemiginfo (&newedgemig); // non-hg version of this
  IMA_reset_edgemiginfo (&newsismig);  // non-hg version of this
  newedgemig.edgeid = edge = oldedgemig.edgeid;
  newedgemig.li = li;
  if (edge < L[li].numgenes)
    newedgemig.upt = 0;
  else
    newedgemig.upt = gtree[gtree[edge].up[0]].time;
  newedgemig.dnt = gtree[edge].time;
  newedgemig.mtall =  newedgemig.dnt - newedgemig.upt;
  newedgemig.mig[0].mt = -1;
  newedgemig.mpall = 0;
  newedgemig.cmm = 0;
  newedgemig.pop = newedgemig.temppop =  gtree[edge].pophg;
  newedgemig.fpop =gtree[gtree[edge].down].pophg;
  if (gtree[edge].down == C[ci]->G[li].root)    /* simulate migration on the sister branch as well */
  {
    newedgemig.fpop = -1;       //pop unknown, as edge must be determined by migration
    if (gtree[gtree[edge].down].up[0] == edge)
      newsis = gtree[gtree[edge].down].up[1];
    else
      newsis = gtree[gtree[edge].down].up[0];
    if (newsis < L[li].numgenes)
      newsismig.upt = 0;
    else
      newsismig.upt = gtree[gtree[newsis].up[0]].time;
    newsismig.edgeid = newsis;
    newsismig.li = li;
    newsismig.dnt = gtree[newsis].time;
    newsismig.mig[0].mt = -1;
    newsismig.cmm = 0;
    newsismig.mpall = 0;
    newsismig.pop = newsismig.temppop = gtree[newsis].pophg;
    newsismig.fpop = -1;    // this is not known until migration events are simulated
    newsismig.mtall = newsismig.dnt - newsismig.upt;
  }
  else
  {
    newsismig.edgeid = -1;
    newsismig.mtall = 0;
    newsismig.b = newsismig.e = -1;// no second edge to deal with
  }

  assert((newsismig.mtall > 0 && newsismig.edgeid >= 0) || newsismig.mtall == 0);
  nm = getmhg (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig,C[ci]->G[li].mhg);  // simulate the migration events
  if (nm < 0)
  {
    //badmigrationupdate += 1;
    return forcereject;
  }

  assert (newedgemig.fpop >= 0);

  mproposedenom = getmprobhg (ci, &newedgemig, &newsismig, &oldedgemig, &oldsismig,oldmhg);
  if (isnan_(mproposedenom) || isinf_DBL(mproposedenom))
  {
    //badmigrationupdate += 1;
    return forcereject;
  }
  /* calculate probability of reverse update    */
  mproposenum = getmprobhg (ci, &oldedgemig, &oldsismig, &newedgemig, &newsismig,C[ci]->G[li].mhg);
  if (isnan_(mproposenum) || isinf_DBL(mproposenum))
  {
    //badmigrationupdate += 1;
    return forcereject;
  }
  weight = mproposenum - mproposedenom + mhgproposeratio;
  return weight;
} /* addmigrationhg */

double getmprobedge(struct edgemiginfo *e, double mhg)
{
  double lr = mhg * e->mtall;
  double tempp, pathlog,popc;
  int lastm_2_pop;

  if (npops == 2)
  {
    if (ODD (e->mpall))
    {
      tempp = e->mpall * log(mhg) - mylogsinh(lr);  // why was this tempp +=  ??
    }
    else
    {
      tempp = e->mpall * log(mhg) - mylogcosh(lr);  // why was this tempp +=  ??
    }
  }
  else                      // 3 or more pops in the last period
  {
    popc = (double) (npops - 1);
    switch (e->mpall)
    {
    case 0:
            tempp = -log(exp(lr)-lr);
            break;
    case 1:
            tempp = log(mhg) - log(exp(lr)-1.0);
            break;
    default:
          {
            if (e->mpall==2)
              lastm_2_pop = e->pop;
            else
              lastm_2_pop = e->mig[e->mpall-3].mp;
            if (lastm_2_pop == e->fpop)
              pathlog = (e->mpall-1) * (-log(popc)) ;
            else
              pathlog = (e->mpall-2) * (-log(popc))  - log(popc - 1.0);
            if (e->pop == e->fpop)
              tempp = pathlog + e->mpall * log(mhg) -log(exp(lr)-lr);
            else
              tempp = pathlog + e->mpall * log(mhg) -log(exp(lr)-1);
          }
    }
  }
  return tempp;
}  //getmprobedge


/* calculate propobility of proposing edgem and sisem
   based simply on simulation probabilities
   note that struct edgemiginfo does not have explicit hg terms
   (the same structure is also used in update_gtree.cpp)
   so pop, fpop, mig in these structures all refer to hidden genealogies in this file
*/
double
getmprobhg(int ci, struct edgemiginfo *edgem,
          struct edgemiginfo *sisem,struct edgemiginfo *oldedgem,
          struct edgemiginfo *oldsisem, double mhg)
{
  double tempp;
  //int lastmigrationperiod; not used
  double r;

  tempp = 0;
  //lastmigrationperiod = IMIN(edgem->e, lastperiodnumber); not uised
  if (sisem->mtall <= 0)       // only deal with edgem
  {
    r = mhg;
    tempp = getmprobedge(edgem,r);
  }
  else  // both edges
  {
    assert(edgem->mtall > 0 && sisem->mtall > 0);
    assert(edgem->e ==sisem->e);
    assert(edgem->fpop == sisem->fpop);
    // probability of simulating fpop
    if (edgem->pop == sisem->pop)
    {
      if (edgem->pop == edgem->fpop)
        tempp = log(MIGCLOSEFRAC);
      else
        tempp = log((1.0-MIGCLOSEFRAC) / (double) (npops-1));
    }
    else
    {
      if (npops == 2)
        tempp = log(0.5);// 50:50 chance of fpop being either
      else
      {
        if (edgem->fpop == edgem->pop || edgem->fpop == sisem->pop)
          tempp = log(0.5* MIGCLOSEFRAC);
        else
          tempp = log((1.0 - MIGCLOSEFRAC)/(double) (npops-2));
      }
    }
    r = mhg;
    tempp += getmprobedge(edgem,r);
    tempp += getmprobedge(sisem,r);
  }
  assert (isnotnan(tempp));
  return tempp;
}                               /* getmprobhg*/


/* set the genealogy edge,  given the hidden genealogy edge and the phylogeny*/
void edgemask(int ci, struct edge *e, double upt)
{
  int mi,ami,cpop,cpophg,topop,topophg;
  int fp;
  cpophg = e->pophg;
  assert (0<=cpophg && cpophg < numtreepops);
  fp = findperiod(ci,upt);
  assert (0<=fp && fp < npops);
  e->pop = C[ci]->ancplist[cpophg][fp]; // the ancestral pop of cpophg in period fp, i.e. pop that includes the pophg in the fp
  cpop = e->pop;
  mi = 0;
  ami = 0;
  while (e->mighg[mi].mt > -0.5)
  {
    topophg = e->mighg[mi].mp;
    assert (0<=topophg && topophg < npops);
    assert(topophg != cpophg);
    fp =  findperiod(ci,e->mighg[mi].mt) ;
    assert (0<=  fp  && fp < npops);
    assert(0<= cpop && cpop< numtreepops);
    cpop = C[ci]->ancplist[cpop][fp];
    assert(0<=topophg && topophg < numtreepops);
    topop = C[ci]->ancplist[topophg][fp];
    if (topop != cpop) // not a hidden migration
    {
      assert (0<=  fp  && fp < npops-1); // can't be an oldest ancestral population so can't be in last period
      e->mig[ami].mp = topop;
      assert (e->mighg[mi].mt >= 0.0);
      e->mig[ami].mt = e->mighg[mi].mt;
      cpop = topop;
      assert (e->mig[ami].mt < C[ci]->tvals[fp]);
      ami += 1;

    }
    mi += 1;
    cpophg = topophg;
  }
  assert (e->mighg[mi].mt < -0.5);
  e->mig[ami].mt = -1;
  e->cmm = ami;
  fp =  findperiod(ci,e->time) ;
  assert (0<= fp && fp < npops);
  //e->fpop = C[ci]->ancplist[e->fpophg][fp]; // this seems likely to be a bug because e->fpophg has not been set yet,  probably minor
  e->fpop = C[ci]->ancplist[cpophg][fp];
  e->fpophg = cpophg;
}  /* edgemask */

/* debug version*/
#ifdef TURNONCHECKS
void edgemask_debugprint(int ci, struct edge *e, double upt)
{
  FILE *ef;
  char debugfilename[] = "edgemaskprint.out";
#define EFP fprintf(ef,
  int mi,ami,cpop,cpophg,topop,topophg;
  int fp;
  if (!(ef = fopen (debugfilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening edgemask debug file for appending");
  }
 EFP "step %d chain %d  edgenum %d  upt %.4lf\n",step,ci,e->ei,upt);
  cpophg = e->pophg;
  assert (0<=cpophg && cpophg < numtreepops);
  fp = findperiod(ci,upt);
  assert (0<=fp && fp < npops);
  e->pop = C[ci]->ancplist[cpophg][fp]; // the ancestral pop of cpophg in period fp, i.e. pop that includes the pophg in the fp
  cpop = e->pop;
  mi = 0;
  ami = 0;
EFP "cpophg %d fp %d cpop %d \n",cpophg,fp,cpop);
  while (e->mighg[mi].mt > -0.5)
  {
    topophg = e->mighg[mi].mp;
    assert (0<=topophg && topophg < npops);
    assert(topophg != cpophg);
    fp =  findperiod(ci,e->mighg[mi].mt) ;
    assert (0<=  fp  && fp < npops);
    assert(0<= cpop && cpop< numtreepops);
    cpop = C[ci]->ancplist[cpop][fp];
    assert(0<=topophg && topophg < numtreepops);
    topop = C[ci]->ancplist[topophg][fp];
  EFP " mi  %d  mighg[mi].mt %.4lf period  %d topophg %d cpop %d topop %d\n",mi,e->mighg[mi].mt,fp,topophg,cpop,topop);
    if (topop != cpop) // not a hidden migration
    {
      assert (0<=  fp  && fp < npops-1); // can't be an oldest ancestral population so can't be in last period
      e->mig[ami].mp = topop;
      assert (e->mighg[mi].mt >= 0.0);
      e->mig[ami].mt = e->mighg[mi].mt;
      cpop = topop;
      assert (e->mig[ami].mt < C[ci]->tvals[fp]);
  EFP "not hidden  ami %d  e->mig[ami].mt %.4lf  e->mig[ami].mp %d \n",ami, e->mig[ami].mt,e->mig[ami].mp);
      ami += 1;

    }
    mi += 1;
    cpophg = topophg;
  }
  assert (e->mighg[mi].mt < -0.5);
  e->mig[ami].mt = -1;
  e->cmm = ami;
  fp =  findperiod(ci,e->time) ;
  assert (0<= fp && fp < npops);
  //e->fpop = C[ci]->ancplist[e->fpophg][fp]; // this seems likely to be a bug because e->fpophg has not been set yet,  probably minor
  e->fpop = C[ci]->ancplist[cpophg][fp];
  e->fpophg = cpophg;
EFP" period %d fpop %d fpophg %d\n",fp,e->fpop,e->fpophg);
FCLOSE (ef);
}  /* edgemask_debugprint */
#endif //TURNONCHECKS

/*********** global functions *********/
void makegenealogy_from_hiddengenealogy(int ci,int li)
{
  int ei;
  double upt;
  struct genealogy *G;
  struct edge *gtree;

  G = &C[ci]->G[li];
  gtree = G->gtree;

  for (ei =0; ei < L[li].numlines; ei++)
  {
    if (ei >= L[li].numgenes)
      upt = gtree[gtree[ei].up[0]].time;
    else
      upt = 0.0;
    assert (upt >= 0.0);
    assert (&gtree[ei]);
    edgemask(ci,&gtree[ei],upt);
    //edgemask_debugprint(ci,&gtree[ei],upt);
  }
}

/* for each population,  count the number of sampled populations that are descendants */
void calcnumancestralpops(int ci, int counts[])
{
  int i,j;
  for (i=0;i<numtreepops;i++)
    counts[i] = 0;
  for (i=0;i<npops;i++)
  {
    j = i;
    while  (C[ci]->poptree[j].down != -1)
    {
      j = C[ci]->poptree[j].down;
      counts[j] += 1;
    }
  }
}

/* go through the mighg array, count hidden migrations on an edge that are between two time points (not inclusive) */
int countnummighg(int ci,int li,int ei,double t, double endt,int *lastpop, int *thirdlastpop)
{
  int c,mi;

  // first go through events up to and including t (t might be at a migration time and we don't want it included)
  mi = 0;
  while (C[ci]->G[li].gtree[ei].mighg[mi].mt <= t && C[ci]->G[li].gtree[ei].mighg[mi].mt > -0.5)
    mi += 1;
  // now go through events up to, but not including endt
  c = 0;
  while (C[ci]->G[li].gtree[ei].mighg[mi].mt < endt && C[ci]->G[li].gtree[ei].mighg[mi].mt > -0.5)
  {
    c += 1; // (C[ci]->G[li].gtree[ei].mighg[mi].mt > t && C[ci]->G[li].gtree[ei].mighg[mi].mt < endt);
    mi += 1;
  }
  if (c>=1)
    *lastpop = C[ci]->G[li].gtree[ei].mighg[mi-1].mp;
  else
      *lastpop = -1;
  if (c>=3)
    *thirdlastpop =C[ci]->G[li].gtree[ei].mighg[mi-3].mp;
  else
    *thirdlastpop = -1;
  return c;
}

/* probability of hidden migration events in an interval of an edge
  If interval ends in a regular migration,  then this calculation includes the
  probability of the hidden population state */
double hgedgecalcprob(struct hgcalcstruct hgc)
{
  double plog;
  plog = hgc.logprobendpop;
  if (hgc.ndpops >=2)
  {
    plog += hgc.nummighg * log(hgc.mrate/(hgc.ndpops-1 ) ) - (hgc.l * hgc.mrate);
  }
  assert(isnotnan(plog));
  assert(isnotinf_DBL(plog));
  return plog;
}  /* hgedgecalcprob */

double hgcalccoalprob(int nm[], int cpop[], int thirdlastpop[], double l[], int fpop, double mrate,int ndpops)
{
  int ii;
  int lastm_2_pop; // population state just prior to the second to last migration
  double p, plog,pathc;

  if (cpop[0]== cpop[1])
  {
    if (cpop[0] == fpop)
      p = MIGCLOSEFRAC;     // prob of final pop being the same as two current pops
    else
      p = (1.0 - MIGCLOSEFRAC) * (1.0/(ndpops-1)); // prob of particular final pop being different from two sampled pops
  }
  else
  {
    if (ndpops == 2)
    {
      p = 0.5;
    }
    else
    {
      if (cpop[0] == fpop || cpop[1] == fpop)
        p = 0.5*MIGCLOSEFRAC;
      else
        p = (1.0 - MIGCLOSEFRAC)/(ndpops-2); // prob of particular final pop being different from both sampled pops
    }
  }
  plog = log(p);
  for (ii=0;ii<2;ii++)
  {
    if (ndpops == 2)
    {
      if (cpop[ii] == fpop)      //even
      {
        assert (ISEVEN(nm[ii]));
        plog += nm[ii] * log(mrate) -  mylogcosh(mrate * l[ii]);
      }
      else                    // odd
      {
        assert (ODD(nm[ii]));
        plog += nm[ii] * log(mrate) -  mylogsinh(mrate * l[ii]);
      }
    }
    else
    {
      if (cpop[ii] == fpop)
      {
        assert(nm[ii] != 1);
        plog += nm[ii] * log(mrate) - log( exp (mrate * l[ii]) -  mrate * l[ii]);
      }
      else
      {
        assert(nm[ii] > 0);
        plog += nm[ii] * log(mrate) - log( exp (mrate * l[ii]) -  1.0);
      }
      if (nm[ii] >= 2 )
      {
        if (nm[ii] == 2)
        {
          lastm_2_pop = cpop[ii];
          assert(lastm_2_pop == thirdlastpop[ii]);
        }
        else
        {
          lastm_2_pop = thirdlastpop[ii];// second to last
        }
        if (lastm_2_pop == fpop)// can go to any pop except the one its in (which is fpop)
          pathc = -log(ndpops-1.0);
        else
          pathc = -log(ndpops-2.0); // can go to any pop except the one its in and fpop
        pathc -= (nm[ii]-2) * log(ndpops-1.0);
        plog += pathc;
      }
    }
  }
 assert(isnotnan(plog));
 assert(isnotinf_DBL(plog));
 return plog;
}  /* hgcalccoalprob */

/*
prob_hg_given_g calculate the probability of the hidden genealogy given the genealogy
Nothing is being set, just calculating the probability of what we find when observing
what the hidden genealogy has,  condition on on what the genealogy has.
*/
double prob_hg_given_g(int ci,int li)
{
  int i,ui,ei,gmi,hgmi,currhgpop,currgpop;
  double uptg,tmig,stime,ctime;
  struct genealogy *G;
  struct edge *gtree;
  double logpsum,mrate;
  struct hgcalcstruct hgc;
  int descounts[MAXTREEPOPS];
  int lastpophg,popbeforesecondtolasthg, reachedcoalescent[2];
  double temp;
  int cnummighg[2];
  int ccurrhgpop[2], ccurrgpop[2];
  int thirdlastpop[2], ndpops;
  double l[2];
  int fpophg;
  //FILE *debugprintfile;
  int debug_probhg_fprint = 0;
  double log_mighgprior;
  int cmmhgcheck[2];
  int nexteventtype,eii[2];

  if (calcoptions[DONTCALCGENEALOGYPRIOR]==1)
    return 0.0;

  /*
  to print to a file a table of numbers used to calculate probhg
  usually use this when checkprobs()  finds a problem with prob_hg_given_g()
  Must set these vars to the step, chain and locus that showed the problem
  int steptoprint = 887;//-1;
  int chaintoprint =4;// -1;
  int locustoprint =0;// -1;
  debug_probhg_fprint = (steptoprint == step && chaintoprint ==ci  && locustoprint == li);
  */
#ifdef TURNONCHECKS
  //gtreeprint(ci,li);
  checkgenealogy(ci,li,0);
#endif //TURNONCHECKS
  /*if (debug_probhg_fprint)  //for tracking down a specific bug,  not sure where this is at  as of 3/31/2017
  {
    char tempdebugfilename[] = "hiddengenealogydebug.out";
    debugprintfile= fopen (tempdebugfilename, "a");
    fprintf(debugprintfile,"\nstep=%d chain=%d locus= %d\n",step,ci,li);
  } */
  G = &C[ci]->G[li];
  gtree = G->gtree;
  logpsum = 0.0;
  mrate = G->mhg;
  if (modeloptions[EXPOMIGRATIONPRIOR])
    log_mighgprior = -(mrate/mprior) - log(mprior);
  else
    log_mighgprior = 0.0;  // technically log (1.0/(C[ci]->imig[li].pr.mean )) but this always cancels
  calcnumancestralpops(ci,descounts); // fill descounts  descounts[i] is the number of descendant populations that population i has

  for (i=L[li].numgenes;i< L[li].numlines;i++) // over internal nodes
  {
    for (ui = 0;ui < 2; ui++)  // over each up edge
    {
      ei = gtree[i].up[ui];
      eii[ui] = ei;
      if (ei < L[li].numgenes) // edge ei is an external edge
        uptg = 0.0;
      else
        uptg = gtree[gtree[ei].up[0]].time;
      ctime = gtree[ei].time; // time when edge coalesces
      currhgpop = gtree[ei].pophg; // hidden pop state at top of edge
      currgpop = gtree[ei].pop;  // pop state at top of edge
      gmi = 0; // indexes regular migration events
      reachedcoalescent[ui] = 0; // array for each edge.
      cmmhgcheck[ui] = 0;  // used to check on counts of hidden migrations
      do // move down the edge, identify all intervals between events in G that could have events in GH
      {
        hgc = hgczero;  // initialize to an empty structure  // this sets logprobendpop to 0.0
        tmig = (gtree[ei].mig[gmi].mt > -0.5) ? gtree[ei].mig[gmi].mt : TIMEMAX;
        stime = C[ci]->poptree[currgpop].time;  // next split time
        hgc.ndpops = descounts[currgpop]; // # descendant populations of current population
        hgc.mrate = mrate;
        nexteventtype = 0;
        if (stime <= DMIN(tmig,ctime)) // current interval ends in a split
        {
          hgc.nummighg = countnummighg(ci,li,ei,uptg,stime,&lastpophg,&popbeforesecondtolasthg);
          hgc.l = stime - uptg;
          assert (hgc.l > 0.0);
          uptg = stime;
          if (lastpophg == -1)  // no hidden migrations,  currhgpop does not change
          {
            assert (hgc.nummighg == 0);
          }
          else// currhgpop changes
          {
            currhgpop = lastpophg;
            assert (hgc.nummighg > 0);
          }
          // when we get to the split the current population changes
          currgpop = C[ci]->poptree[currgpop].down;
          nexteventtype = 1;
        }
        if (tmig <= DMIN(stime,ctime))  // current interval ends in a migration
        {
          assert(tmig==gtree[ei].mig[gmi].mt);
          hgc.logprobendpop = (gtree[ei].mig[gmi].mp < npops) ? 0.0 : -log((double) descounts[gtree[ei].mig[gmi].mp]); // if not a sampled pop,  log of 1 over number ancestral pops
          hgc.nummighg = countnummighg(ci,li,ei,uptg,tmig,&lastpophg,&popbeforesecondtolasthg);
          hgc.l = tmig - uptg;
          assert (hgc.l > 0.0);
          uptg = tmig;
          hgmi = 0; // find the migration in mighg that is the same event as tmig in order to set currhgpop for the next interval
          while (gtree[ei].mighg[hgmi].mt != tmig && gtree[ei].mighg[hgmi].mt > -0.5)
            hgmi++;
          assert(gtree[ei].mighg[hgmi].mt > -0.5); // would only be < 0 if it was not found
          assert(gtree[ei].mighg[hgmi].mt == tmig);
          currhgpop = gtree[ei].mighg[hgmi].mp;
          assert(currhgpop < npops);
          currgpop = gtree[ei].mig[gmi].mp;
          gmi += 1;
          cmmhgcheck[ui] += 1;
          nexteventtype = 2;
        }
        if (ctime <= DMIN(stime,tmig) ||(stime >= TIMEMAX && tmig >= TIMEMAX))  // current interval ends in a coalescent
        {
          l[ui] = ctime - uptg;
          assert (l[ui] > 0.0);
          cnummighg[ui] = countnummighg(ci,li,ei,uptg,ctime,&lastpophg,&popbeforesecondtolasthg);
          ccurrgpop[ui] = currgpop;
          ccurrhgpop[ui] = currhgpop;
          if (cnummighg[ui]>=3)
          {
            assert( popbeforesecondtolasthg >= 0);
            thirdlastpop[ui] = popbeforesecondtolasthg;
          }
          else
          {
            if (cnummighg[ui]==2)
              thirdlastpop[ui] = currhgpop;
            else
              thirdlastpop[ui] = -1;
          }
          fpophg = gtree[ei].fpophg;
          assert(lastpophg==-1 || fpophg == lastpophg);
          reachedcoalescent[ui] = 1;
          nexteventtype = 3;
        }
        assert (nexteventtype != 0);
        if ( reachedcoalescent[ui] == 0)
        {
           // BAD BUG  3/27/2017  found by Yujin.
            //hgc.logprobendpop can be nonzero even when ndpops is < 2 because migration from non-ancestor to ancestor  has multiple possibilities
            //  Even though there can be no hidden migrations in this case,  there is still a nonzero probability for the interval.
            // edited hgedgecalcprob() to handle case when ndpops < 2
          cmmhgcheck[ui] += hgc.nummighg;
          temp = hgedgecalcprob(hgc);
          logpsum += temp;
        }
      }while (reachedcoalescent[ui] == 0);
    }
    assert (ccurrgpop[0] == ccurrgpop[1]);
    assert(l[0] > 0.0 && l[1] > 0.0);
    if (ccurrgpop[0] >=  npops) // else do nothing as the probability of the hidden genealogy associated with an edge that coalesces before the first time split  must be 1
    {
      ndpops = descounts[ccurrgpop[0]];
      if (ndpops >= 0)
      {
        temp = hgcalccoalprob(cnummighg,ccurrhgpop,thirdlastpop,l,fpophg,mrate,ndpops);
        assert (temp != 0.0 && temp==temp &&  !isinf_DBL(temp));
#ifdef DEBUG
        cmmhgcheck[0] += cnummighg[0];
        cmmhgcheck[1] += cnummighg[1];
        assert (cmmhgcheck[0]==C[ci]->G[li].gtree[eii[0]].cmmhg);
        assert (cmmhgcheck[1]==C[ci]->G[li].gtree[eii[1]].cmmhg);
#endif
        logpsum += temp;
       /*if (debug_probhg_fprint)
        {
          fprintf(debugprintfile,"i cnummighg[0] cnummighg[1] ccurrhgpop[0] ccurrhgpop[1] thirdlastpop[0] thirdlastpop[1] l[0] l[1] fpophg mrate ndpops logpsum\n");
          fprintf(debugprintfile,"%d %d %d %d %d %d %d %.6lf %.6lf %d %.4f %d %.6lf\n",i,cnummighg[0],cnummighg[1],ccurrhgpop[0],ccurrhgpop[1],thirdlastpop[0],thirdlastpop[1],l[0],l[1],fpophg,mrate,ndpops,logpsum);
        } */
      }
    }
  }

  /*if (debug_probhg_fprint)
  {
    fprintf(debugprintfile,"final logpsum %.6lf\n",logpsum);
    FCLOSE (debugprintfile);    
  } */
  assert( isnotinf_DBL(logpsum));
  assert(isnotnan(logpsum));
  return (logpsum  + log_mighgprior);
}   /*prob_hg_given_g*/

/*
  save all info about Gh and G
  make a change to Gh
  mask this change onto G
  MH term includes:
    likelihood for locus
    Prog_Gh_given_G
    Prog_G
    Proposal prob ratio for Gh

main steps in updateHG():
  copy_hginfo() [locus li, chain ci]
  copy_hginfo()  [chain ci, overall]
  storehgstats()
  do sliding update
  rebuild hidden genealogy with update
 *maskgenealogy()     [make the new genealogy, given the new hg]
  treeweight(ci,li)
  sum_subtract_treeinfo()    [updates the summary stats for the joint genealogies]
  calculate the likelihood on the genealogy depending on mutation model
  copy_probcalc()
  calculate the MH term, including:
      ratio of locus specific likelihoods
      ratio of overal G prior: C[ci]->allpcalc.probg
      *ratio of GH given G prior for locus li
      migweight, slideweight from genealogy proposal
      beta
  if accept:
    update C[ci]->allpcalc
    update G->pdg
  else [reject update]:
    restoreedges()
 *maskgenealogy()     [restore the old genealogy, given the old hg]
    storehgstats()  [copy stats back]
    copy_probcalc()
    copy_treeinfo() [locus li, chain ci]

*/
int  update_hidden_genealogy(int ci,  int li, int *topolchange, int *tmrcachange)
{
  int ai;
  int edge, oldsis, newsis, freededge, accp;
  double newpdg, newpdg_a[MAXLINKED];
  double migweight, metropolishastingsratio;
  double likelihoodratio,priorratio,proposalratio;
  double Aterm[MAXLINKED], Atermsum;
  double slidedist;
  double slideweight, holdslidedist, slidestdv;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int autoreject;
  double like;
  int joinstatus;
  double oldmhg;


#ifdef TURNONCHECKS
checkgenealogyweights(ci);


//gtreeprint(ci,li);
//checkgenealogy(ci,li,0);
//  printgenealogyweights(ci,li);
checkprobs(ci,li);
#endif //TURNONCHECKS

  //C[ci]->G[li].mhg =  mprior * MPRIORFRACFORHG; //  this is  done in initialize

// initialize and make copies structures that hold quantities for calculating prob of genealogy
  copy_treeinfo (&holdgweight_updategenealogy, &G->gweight);
  copy_treeinfo (&holdallgweight_updategenealogy, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc_updategenealogy, &C[ci]->allpcalc); /* save old values */

  // store summary stats of the genealogy
  storegenealogystatshg (ci, li, 0);
  *tmrcachange = 0;
  *topolchange = 0;
  Atermsum = 0;   // Atermsum only used for Stepwise mutation model

/* pick an edge, identify freedup edge (the down edge) and the sister edge */
  do
  {
    edge = randposint (L[li].numlines);
  } while (gtree[edge].down == -1);
  freededge = gtree[edge].down;
  if ((oldsis = gtree[freededge].up[0]) == edge)
    oldsis = gtree[freededge].up[1];
  /* copy information on the edge,  and if it connects to the root, then the sister edge as well */
  if (gtree[edge].down == G->root)
    fillmiginfohg (ci, li, gtree, edge, oldsis);
  else
    fillmiginfohg (ci, li, gtree, edge, -1);
  /* copy information on the genealogy before changing it */
  storeoldedgeshg (ci, li, edge, oldsis, freededge);

  /* slide edge, pick a distance and slide it  */
/* PROBLEM it is possible for this to generate too big a slide distance if the sample size is large and there are a lot of very short edges
  this causes the recursion in slider to crash*/
  slidestdv = DMIN (SLIDESTDVMAX, G->roottime/3 );
//slidestdv = slidetemp[slidei]* T[numsplittimes-1].pr.max;
  holdslidedist = slidedist = normdev (0.0, slidestdv);
/* join the sister and the down branches at the point where edge used to connect, this frees up the down branch */
  joinstatus = joinsisdown (ci, li, oldsis, tmrcachange);  // no hg version yet
  if (joinstatus < 0)
    return -1;   // awkward way to get out of here. this happens if too large a value of cmm or cmmhg happens because of the join
/* do the slide and identify the new sister branch and where new connection point for the edge is */
  newsis = oldsis;
  slider (ci, li, edge, &newsis, &(gtree[edge].time), &slidedist);  // same as used by regular genealogy update in update_gtree.cpp
  *topolchange += (oldsis != newsis);
/* now separate the new sister branch into a shorter sis branch and a down branch  */
  splitsisdown (ci, li, edge, freededge, newsis);  // same as used by regular genealogy update in update_gtree.cpp
//if (0)  use with slidetemp and slidei
  if (rootmove)  // have to calculate slideweight as slidestdv may change
  {
    slideweight = -log (normprob (0.0, slidestdv, holdslidedist));
    slidestdv = DMIN (SLIDESTDVMAX, G->roottime / 3);
    slideweight += log (normprob (0.0, slidestdv, holdslidedist));
  }
  else
  {
    slideweight = 0;
  }

/* add migration events */
  oldmhg = G->mhg;
  migweight = addmigrationhg (ci, li);
  /* can force rejection if:
      addmigration tries do exceed the migration arrays, in which case migweight == forcereject
      or
      treeweight() fails - there are some rare cases where the sequence of events handled in treeweight just does not make a valid tree,  hard to debug
  */
  if (migweight > forcereject)  // adding migration went ok
  {
  /*  copy the migration info in newedgemig and newsismig  to the genealogy */
  /* determine all the weights needed for calculating the probability of the genealogy */
    copynewmighg_to_gtree (ci, li);
    makegenealogy_from_hiddengenealogy(ci,li);
#ifdef  TURNONCHECKS
//  gtreeprint(ci,li);
#endif //TURNONCHECKS
    setzero_genealogy_weights (&G->gweight);
    int treeweightcallcode = 2;  // a debugging code,  if treeweight has an error results are written to an output file with this code
    int tw = treeweight (ci, li,treeweightcallcode);  //calculate and set C[ci]->G[li].gweight
    if (tw < 0) // some failure in treeweight()
    {
      autoreject = 1;
      accp = 0;
    }
    else
      autoreject = 0;
  }
  else // migweight too small,  should only get here if too many migration events were added in addmigrationhg()
  {
    autoreject = 1;
    accp = 0;
  }
  if (autoreject == 0)
  {
    sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight,&holdgweight_updategenealogy);// after update, add the new weights to allgweight and subtract the old ones
    newpdg = 0;  /* will hold the new value for p(data|genealogy) at locus li.  The current one is in G->pdg  */
    switch (L[li].model)
    {
    case HKY:
        newpdg = newpdg_a[0] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, edge,
                         freededge, oldsis, newsis);
      break;
    case INFINITESITES:
      newpdg = newpdg_a[0] = like = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      {
        for (ai = 0, newpdg = 0; ai < L[li].nlinked; ai++)
        {
          newpdg_a[ai] =
            G->pdg_a[ai] + oldfinishSWupdateA (ci, li, ai, edge, freededge,oldsis, newsis, G->uvals[ai], &Aterm[ai]);
          newpdg += newpdg_a[ai];
          Atermsum += Aterm[ai];
        }
  //            checklikelihoodSW(ci, li,G->u[ai].mcinf.val);
        break;
      }
    case JOINT_IS_SW:
      newpdg = newpdg_a[0] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        newpdg_a[ai] =
          G->pdg_a[ai] + oldfinishSWupdateA (ci, li, ai, edge, freededge,
                                          oldsis, newsis,
                                          G->uvals[ai], &Aterm[ai]);
        newpdg += newpdg_a[ai];
        Atermsum += Aterm[ai];
      }
      //checklikelihoodSW(ci, li,Q[ci]->us[li]);
      break;
    }

#ifdef TURNONCHECKS
    //gtreeprint(ci,li);


#endif //TURNONCHECKS

/* metropolis-hastings weight calculation */
    /*C[ci]->allpcalc.probhgg is a simple sum of hgprob values,  one of which is being updated */
    C[ci]->allpcalc.probhgg -= C[ci]->G[li].hgprob; // subtract out old value
    C[ci]->G[li].hgprob = prob_hg_given_g(ci,li);   // calculate new value
    C[ci]->allpcalc.probhgg += C[ci]->G[li].hgprob;   // add in new value
    priorratio = C[ci]->allpcalc.probhgg - holdallpcalc_updategenealogy.probhgg;
    priorratio -= C[ci]->allpcalc.probg;
    integrate_tree_prob (ci, &C[ci]->allgweight,
                         &holdallgweight_updategenealogy, &C[ci]->allpcalc,
                         &holdallpcalc_updategenealogy);
    priorratio += C[ci]->allpcalc.probg;
    /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
    likelihoodratio =  (newpdg - G->pdg);
    proposalratio = migweight + slideweight + Atermsum;
    if (calcoptions[CALCMARGINALLIKELIHOOD])
    	{
      metropolishastingsratio = beta[ci] * likelihoodratio + priorratio + proposalratio;
    }
    else
    {
       metropolishastingsratio = beta[ci] * (likelihoodratio + priorratio)+ proposalratio;
    }
    if (metropolishastingsdecide(metropolishastingsratio,1))
    {
      /* accept the update */

      /* change pdg for the locus, and in allpcalc for the entire chain */
      C[ci]->allpcalc.pdg -= G->pdg;
      C[ci]->allpcalc.pdg += newpdg;
      G->pdg = newpdg;

      for (ai = 0; ai < L[li].nlinked; ai++)
        G->pdg_a[ai] = newpdg_a[ai];
      if (L[li].model == HKY)
      {
        copyfraclike (ci, li);
        storescalefactors (ci, li);
      }
      accp = 1;
#ifdef TURNONCHECKS
    //gtreeprint(ci,li);


  checkgenealogy(ci,li,0);
  checkgenealogyweights(ci);
  checktreeweight(ci,li);
#endif //TURNONCHECKS
    }
    else
      accp = 0;
  }
  /* reject the update */
  if (accp == 0)
  {
    // put the edges back
    restoreedgeshg (ci, li, edge, oldsis, freededge, newsis);
    // copy summary stats back
    storegenealogystatshg (ci, li, 1);
    // reset HKY terms
    if (L[li].model == HKY)
      restorescalefactors (ci, li);
    // copy back all the weights and results associated with calculating the probability of the genealogy
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_updategenealogy);
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_updategenealogy);
    copy_treeinfo (&G->gweight, &holdgweight_updategenealogy);
    *topolchange = 0;
    *tmrcachange = 0;
    C[ci]->G[li].mhg = oldmhg;

#ifdef TURNONCHECKS
//  printgenealogyweights(ci,li);
//  gtreeprint(ci,li);

  checkgenealogy(ci,li,0);
  checkgenealogyweights(ci);
    checktreeweight(ci,li);
#endif //TURNONCHECKS
  }
#ifdef TURNONCHECKS
   checkprobs(ci,li);
#endif //TURNONCHECKS
  return accp;
}   /* update_hidden_genealogy */
