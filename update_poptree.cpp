/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS

#include "ima.hpp"

extern int numdistinctpopulationpairs[]; // = {0,0,1,6,25,90,301,966,3025}; /* number of distinct pairs of populations that could engage in gene flow (don't share any descendant pops) */

extern void fillplist (int ci);
extern void fillancplist (int ci);
extern void poptreewrite (int ci, char *buildstr);
extern void orderupnodes(struct popedge *poptree);

/*********** local to this file  ***************/
struct popedge poptreehold[2*MAXPOPS-1];
#define SLIDESTDVMAX 0.2   // not sure what's best, copied from update_ptree.cpp
#define SLIDESTDVMIN 0.001  // not sure what's best, copied from update_ptree.cpp
#define POPSLIDESTDV 0.05
#define POPSLIDEDIVISOR 10

static struct genealogy_weights *holdgweight_array;
static struct genealogy_weights holdallgweight;
static struct probcalc holdallpcalc;
static double **topologypoppairpriors; // holds log of prior for each possible pair of sampled sister populations


/******* local prototypes ********/
static void storegenealogystats_all_loci_hg (int ci, int mode);
static void copy_poptree (struct popedge *frompoptree,struct popedge *topoptree);
static void renumberpopnodes(int ci,double *tvals);
static void poptreejoinsisdown (struct popedge *ptree,int sis, int *roottimechanges);
static void poptreesplitsisdown (struct popedge *ptree, int slidingedge, int down, int newsis,int *roottimechanges);
static void poptreeslider (struct popedge *ptree, int slidingedge, int *sis,
                    double *timepoint, double *slidedist);
static void poptreesliderlimited (struct popedge *ptree, int slidingedge, int *sis,
                    double *timepoint, double *slidedist,double *upperiodtimelimit,double *dnperiodtimelimit);
static void init_holdgweight_array_for_hg(void);
static void free_holdgweight_array_for_hg (void);
static void copyimigpriors(int ci, int mode);

/******* local functions ********/

/* there is also one of these in update_t_RYhg.cpp */ 
void
storegenealogystats_all_loci_hg (int ci, int mode)
{
  static double holdlength[MAXLOCI], holdtlength[MAXLOCI];
  static double holdroottime[MAXLOCI];
  static int holdroot[MAXLOCI];
  static int holdmig[MAXLOCI];
  static double holdhgprob[MAXLOCI];
  int li;
  if (mode == 0)
  {
    for (li = 0; li < nloci; li++)
    {
      holdlength[li] = C[ci]->G[li].length;
      holdtlength[li] = C[ci]->G[li].tlength;
      holdroottime[li] = C[ci]->G[li].roottime;
      holdroot[li] = C[ci]->G[li].root;
      holdmig[li] = C[ci]->G[li].mignum;
      holdhgprob[li] = C[ci]->G[li].hgprob;
    }
  }
  else
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].length = holdlength[li];
      C[ci]->G[li].tlength = holdtlength[li];
      C[ci]->G[li].mignum = holdmig[li];
      C[ci]->G[li].roottime = holdroottime[li];
      C[ci]->G[li].root = holdroot[li];
      C[ci]->G[li].hgprob = holdhgprob[li];
    }
  }
  return;
}                               // storegenealogystats_all_loci_hg  

void
copy_poptree (struct popedge *frompoptree,struct popedge *topoptree)
{
  int i;
  for (i=0;i<numtreepops;i++)
    topoptree[i] = frompoptree[i];
  return;
}

/* the system requires that the ancestral node numbers be ordered in time.
   Must get the current sequence of pop nodes in time
   and then relabel them in time order  and rebuild the poptree.
   This is confusing. 
   */
/* this is old code,  hard to follow, checked that it gives same results as new code */
void renumberpopnodes_hold(int ci,double *tvals)
{
  double tlast,upt,dnt;
  struct popedge  *ptree;
  struct pnodetime  pnodetimes[MAXPOPS];
  struct popedge newpoptree[2*MAXPOPS];
  int i,j,needsort,e,b;
  int  rejigger[MAXTREEPOPS];// what was node i in the old tree will become node rejigger[i] 
  //int numtreepops = 2*npops - 1; is global
  tlast = 0.0;
  needsort = 0;
  ptree = C[ci]->poptree;
  // go through ancestral pops in order by their number and get their up times (when they began)
  // if their times are not in order then it needs sorting 
  for (i=npops;i<numtreepops;i++)
  {
    j = i-npops;
    pnodetimes[j].ptreepos = i;
    pnodetimes[j].time = ptree[ptree[i].up[0]].time; // the up time  
    if (pnodetimes[j].time < tlast)
      needsort = 1;
    tlast= pnodetimes[j].time;
  }
  // set values for population that connects to root
  pnodetimes[npops-1].ptreepos = -1;
  pnodetimes[npops-1].time = TIMEMAX;
  if (needsort)
  {
    shellpnodetimes (&pnodetimes[0],npops); // sort by time 
    // fill newpoptree[] and rejigger[]
    for (i=0;i<npops;i++)
    {
      newpoptree[i].up[0] = newpoptree[i].up[1] = -1;
      rejigger[i] = i; // sampled pops always have the same number
    }
    for (;i<numtreepops;i++)
    {
      newpoptree[i].numup = 2;
      newpoptree[i].up[0] = newpoptree[i].up[1] = -1;
      // set rejigger for ancestral pop i,  gets the current id of the population with ordered uptime in that slot
      // i.e. the id of the node with an uptime at the bottom of period npops-i should be i,  but it is  pnodetimes[i-npops].ptreepos
      rejigger[i] = pnodetimes[i-npops].ptreepos; 
    }
    for (i=0;i<numtreepops;i++)
    {
      newpoptree[i].up[0] = newpoptree[i].up[1] = -1; // already did this,  could comment out
    }
    for (i=0;i<numtreepops;i++)
    {
      j = 0; // identify 
      while (ptree[rejigger[i]].time > pnodetimes[j].time)
        j++;
      newpoptree[i].time = pnodetimes[j].time;
      if (j == npops-1)
        newpoptree[i].down = -1;
      else
        newpoptree[i].down = j+npops;
    }
    for (i=0;i<numtreepops;i++)
    {
      if (newpoptree[i].down  != -1)
      {
        if (newpoptree[newpoptree[i].down].up[0]  == -1)
          newpoptree[newpoptree[i].down].up[0] = i;
        else
          newpoptree[newpoptree[i].down].up[1] = i;
      }
    }
#ifdef TURNONCHECKS
    //poptreeprint_frompointer(ci,&newpoptree[0]);
#endif //TURNONCHECKS
    /* can comment out from here to x  if want to print the tree without affecting ptree */
    for (i=0;i<numtreepops;i++)
    {
      ptree[i].up[0] = newpoptree[i].up[0];
      ptree[i].up[1] = newpoptree[i].up[1];
      ptree[i].time = newpoptree[i].time;
      ptree[i].down = newpoptree[i].down;
    } /* x*/
  }  
  /* can comment out from here to xx  if want to print the tree without affecting ptree */
  for (i=0;i<numtreepops;i++) // set b and e regardless of whether node number order was out of wack. 
  {
    if (i<npops)
    {
      upt = 0.0;
      b = 0;
    }
    else
    {
      upt = ptree[ptree[i].up[0]].time;
      j = 0;
      while (upt > pnodetimes[j].time)
        j++;
      b = j+1;
    }
    if (ptree[i].down == -1)
    {
      e = -1;
    }
    else
    {
      dnt = ptree[i].time;
      j = 0;
      while (dnt > pnodetimes[j].time)
        j++;
      e = j+1;
    }
    ptree[i].b = b;
    ptree[i].e = e;
  }
  for (i=0;i<npops;i++)
    tvals[i] = pnodetimes[i].time; /*xx */
} // renumberpopnodes_hold

/* this is cleaner and clearer than the original code (now called renumberpopnodes_hold()) */
void renumberpopnodes(int ci,double *tvals)
{
  double tlast,upt,dnt;
  struct popedge  *ptree;
  struct pnodetime  pnodetimes[MAXPOPS];
  struct popedge newpoptree[2*MAXPOPS];
  int i,j,numtreepops,needsort,e,b;
  int  rejigger[MAXTREEPOPS];// what was node i in the old tree will become node rejigger[i] 
  numtreepops = 2*npops - 1;
  tlast = 0.0;
  needsort = 0;

  ptree = C[ci]->poptree;
  // go through ancestral pops in order by their number and get their up times (when they began)
  // if their times are not in order then it needs sorting 
  for (i=npops;i<numtreepops;i++)
  {
    j = i-npops;
    pnodetimes[j].ptreepos = i;  // current position in the ptree array
    pnodetimes[j].time = ptree[ptree[i].up[0]].time; // the up time  
    if (pnodetimes[j].time < tlast)
      needsort = 1;
    tlast= pnodetimes[j].time;
  }
  // set values for population that connects to root
  pnodetimes[npops-1].ptreepos = -1;  
  pnodetimes[npops-1].time = TIMEMAX;
  if (needsort)
  {
    // sort by time 
    // ancestral populations have id numbers that are ordered with respect to  time, 
    // so the position in the sorted list (+ npops) is the node number that should be held for an ancestral population at that time 
    shellpnodetimes (&pnodetimes[0],npops); 
    // fill rejigger[]
    for (i=0;i<npops;i++)
    {
      rejigger[i] = i; // sampled pops always have the same number
    }
    // rejigger[i] should hold the position in the new tree that i has in the old tree 
    for (;i<numtreepops;i++)
      for (j=0;j<npops-1;j++)
        if (pnodetimes[j].ptreepos == i)
          rejigger[i] = j + npops;
    // go through and fill newpoptree with rejigger values 
    for (i=0;i<npops;i++)
    {
      if (ptree[i].up[0] < npops)
        newpoptree[i].up[0] =  ptree[i].up[0];
      else
        newpoptree[i].up[0] =  rejigger[ptree[i].up[0]];
      if (ptree[i].up[1] < npops)
        newpoptree[i].up[1] =  ptree[i].up[1];
      else
        newpoptree[i].up[1] =  rejigger[ptree[i].up[1]];
      newpoptree[i].down =  rejigger[ptree[i].down];
      if (newpoptree[i].down == -1)
        newpoptree[i].time = TIMEMAX;
      else
        newpoptree[i].time = pnodetimes[newpoptree[i].down - npops].time;
    }
    for (;i<numtreepops;i++)
    {
      newpoptree[rejigger[i]].time =  ptree[i].time;
      if (ptree[i].up[0] < npops)
        newpoptree[rejigger[i]].up[0] =  ptree[i].up[0];
      else
        newpoptree[rejigger[i]].up[0] =  rejigger[ptree[i].up[0]];
      if (ptree[i].up[1] < npops)
        newpoptree[rejigger[i]].up[1] =  ptree[i].up[1];
      else
        newpoptree[rejigger[i]].up[1] =  rejigger[ptree[i].up[1]];
      if (ptree[i].down >= 0)
        newpoptree[rejigger[i]].down =  rejigger[ptree[i].down];
    }
    newpoptree[i].down = -1;
  
#ifdef TURNONCHECKS
    //poptreeprint_frompointer(ci,&newpoptree[0]);
#endif //TURNONCHECKS
    // copy back to ptree 
    for (i=0;i<numtreepops;i++)
    {
      ptree[i].up[0] = newpoptree[i].up[0];
      ptree[i].up[1] = newpoptree[i].up[1];
      ptree[i].time = newpoptree[i].time;
      ptree[i].down = newpoptree[i].down;
    }
    orderupnodes(C[ci]->poptree);
  }
  ptree[numtreepops-1].down = -1;
  for (i=0;i<numtreepops;i++) // set b and e regardless of whether node number order was out of wack. // why?? 
  {
    if (i<npops)
    {
      upt = 0.0;
      b = 0;
    }
    else
    {
      upt = ptree[ptree[i].up[0]].time;
      j = 0;
      while (upt > pnodetimes[j].time)
        j++;
      b = j+1;
    }
    if (ptree[i].down == -1)
    {
      e = -1;
    }
    else
    {
      dnt = ptree[i].time;
      j = 0;
      while (dnt > pnodetimes[j].time)
        j++;
      e = j+1;
    }
    ptree[i].b = b;
    ptree[i].e = e;
  }
  
  for (i=0;i<npops;i++)
    tvals[i] = pnodetimes[i].time;
} // renumberpopnodes

void
poptreejoinsisdown (struct popedge *ptree,int sis, int *roottimechanges)
{
  /* remove an edge
   * extend sis, and free up the down edge 
   the sliding edge (sister of sis) will start its slide on sis
   if the sliding edge and sis came together at the root node 
   then disconnecting the sliding edge and rejoining it means that 
   the time of the root node will necessarily change,  even 
   if the topology does not end up changing
  */
  int  downdown, down;
  down = ptree[sis].down;
  ptree[sis].time = ptree[down].time;

  /* set the up to which sis now connects */
  ptree[sis].down = ptree[down].down;
  downdown = ptree[sis].down;
  if (downdown != -1)
  {
    if (ptree[downdown].up[0] == down)
      ptree[downdown].up[0] = sis;
    else
      ptree[downdown].up[1] = sis;
    *roottimechanges = 0; // might change to 1 after the slide 
  }
  else
  {
    //sis now connects to the root, root time must change 
    ptree[sis].down = -1;
    ptree[sis].time = TIMEMAX;
    *roottimechanges = 1;
  }
}                               /* poptreejoinsisdown */

void
poptreesplitsisdown (struct popedge *ptree,int slidingedge, int down, int newsis,int *roottimechanges)
{
  /* split newsis into two parts, and make a new down edge out of the lower part 
    if the edge that is being split was previously connected to the root node, 
    then the root will move to the new time and roottimechanges will be set to 1. 
    This happens regardless of what roottimechanges was set to by poptreejoinsisdown().
  */
  
  int downdown;
  double curt;
  curt = ptree[slidingedge].time;
  ptree[down].time = ptree[newsis].time;
  ptree[newsis].time = curt;

  /* set the up  of the edge to which down now connects, depends on whether newsis is the root */
  downdown = ptree[newsis].down;
  if (downdown != -1)
  {
    if (ptree[downdown].up[0] == newsis)
      ptree[downdown].up[0] = down;
    else
      ptree[downdown].up[1] = down;
  }
  else
  {
    /* newsis is the current root so the root must move down */
    *roottimechanges = 1;
  }
  ptree[down].down = downdown;
  ptree[newsis].down = ptree[slidingedge].down = down;
  ptree[down].up[0] = newsis;
  ptree[down].up[1] = slidingedge;
  return;
}                               /* poptreesplitsisdown */

void
poptreeslider (struct popedge *ptree, int slidingedge, int *sis,
                    double *timepoint, double *slidedist)
/* this is copied from the version in update_gtree.cpp    */
/* timepoint points at ptree[*slidingedge].time and is the current position of the sliding point, 
slidedist is the distance it must move 
do not restructure the ptree. just figure out when and on which branch timepoint ends up on, 
this will be the new sisterbranch
use recursion */
{
  double uplimit,dnlimit;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;  // set negative value to positive 
    if (slidingedge < npops)
      uplimit = 0;
    else
      uplimit = ptree[ptree[slidingedge].up[0]].time;  // up limit is the top of the sliding edge 
    assert (*timepoint >= uplimit);

    /* if uplimit >= sis up time  - slidingedge cannot reach a node.  if sis is sampled pop,  then up time is 0 and uplimit is obviously greater than this */
    if (ptree[*sis].up[0] == -1 ||  uplimit >= ptree[ptree[*sis].up[0]].time)
    {
      if (*slidedist < (*timepoint - uplimit))
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > uplimit);
        return;
      }
      else
      {
        /* slide up and reflect, sis remains the same, leave slidedist positive so slidingedge goes down with next call to poptreeslider */
        *slidedist -= (*timepoint - uplimit);
        *timepoint = uplimit;
        assert (*slidedist > 0);
        poptreeslider (ptree, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
    else
    {
      /* uplimit is less than sis up time, and thus slidingedge can reach a node */
      assert(*timepoint >= ptree[ptree[*sis].up[0]].time);
      if (*slidedist < *timepoint - ptree[ptree[*sis].up[0]].time)
        /* slide up and stop,  sis remains the same */
      {
        *timepoint -= *slidedist;
        *slidedist = 0;
        assert (*timepoint > ptree[ptree[*sis].up[0]].time);
        return;
      }
      else
      {
        /* slide up and reach a node, pick one side at random and recurse */
        *slidedist -= *timepoint - ptree[ptree[*sis].up[0]].time;
        *timepoint = ptree[ptree[*sis].up[0]].time;
        if (bitran () /*uniform() < 0.5 */ )
        {
          *sis = ptree[*sis].up[0];
        }
        else
        {
          *sis = ptree[*sis].up[1];
        }
        /* reset slidedist to negative, so slidingedge continues up the ptree in next call to poptreeslider */
        *slidedist = -*slidedist;
        assert (*slidedist < 0);
        poptreeslider (ptree, slidingedge, sis, timepoint, slidedist);
        return;
      }
    }
  }
  else
  {
    /* go down */
    if (ptree[*sis].down == -1 )// sis edge is the root, slide heads down and moves the root
    {
      dnlimit = T[numsplittimes-1].pr.max;   // time upper bound is the maximum time of the basal split in the tree. No point in considering a branch point older than this. 
      assert(dnlimit > *timepoint);
      /* if slide distance would take edge past the maximal time, slide down and then start back up */
      if (*timepoint + *slidedist > dnlimit)
      {
        /* slide down the sis and then backup */
        *slidedist -= (dnlimit - *timepoint);
        *slidedist = -*slidedist;
        *timepoint = dnlimit;
        poptreeslider (ptree, slidingedge, sis, timepoint, slidedist);
        return;
      }
      else
      {
        /* just slide down and done */
       *timepoint += *slidedist;
       *slidedist = 0;
       return;
      }
    }
    else 
    {
      if (*timepoint + *slidedist < ptree[*sis].time) // cannot reach the next down node,  
      {
         *timepoint += *slidedist;
         *slidedist = 0;
         return;
      }
      else
      {
        /* a down node is reached */
        assert(ptree[*sis].time > *timepoint);
        *slidedist -= (ptree[*sis].time - *timepoint);
        *timepoint = ptree[*sis].time;
        if (bitran ())
        {
          /* begin to slide down the down node, slidedist stays negative */
          *sis = ptree[*sis].down;
          poptreeslider (ptree, slidingedge, sis, timepoint, slidedist);
          return;
        }
        else
        {
          /* begin to slide up the sister of sis,  identify sis and make slidedist positive */
          if (ptree[ptree[*sis].down].up[0] == *sis)
          {
            *sis = ptree[ptree[*sis].down].up[1];
          }
          else
          {
            *sis = ptree[ptree[*sis].down].up[0];
          }
          *slidedist = -*slidedist;
          poptreeslider (ptree, slidingedge, sis, timepoint, slidedist);
          return;
        }
      }
    }
  }
}                               /* poptreeslider */

void
poptreesliderlimited (struct popedge *ptree, int slidingedge, int *sis,
                    double *timepoint, double *slidedist,double *upperiodtimelimit,double *dnperiodtimelimit)
/*
move up and down, with recursion, along a single sister edge between two time points
*/
{
  double uplimit,dnlimit;
  if (*slidedist < 0)
  {
    /* go up */
    *slidedist = -*slidedist;
    uplimit = *upperiodtimelimit;
    assert (*timepoint >= uplimit);
    if (*slidedist < *timepoint - uplimit)
      /* slide up and stop,  sis remains the same */
    {
      *timepoint -= *slidedist;
      *slidedist = 0;
      assert (*timepoint > uplimit);
      return;
    }
    else
    {
      /* slide up and reflect at uplimit, sis remains the same, leave slidedist positive so slidingedge goes down with next call to poptreesliderlimited */
      *slidedist -= *timepoint - uplimit;
      *timepoint = uplimit;
      assert (*slidedist > 0);
      poptreesliderlimited (ptree, slidingedge, sis, timepoint, slidedist,upperiodtimelimit,dnperiodtimelimit);
      return;
    }
  }
  else
  {
    /* go down */
    dnlimit = DMIN(T[numsplittimes-1].pr.max, *dnperiodtimelimit);
    if (*timepoint + *slidedist < dnlimit) // slide down and done
    {
        *timepoint += *slidedist;
        *slidedist = 0;
        return;
    }
    else /* dnlimit is reached,  slide back up sis */
    {
      *slidedist -= dnlimit - *timepoint;
      assert (*slidedist > 0);
      *timepoint = dnlimit;
      *slidedist = -*slidedist;
      poptreesliderlimited (ptree, slidingedge, sis, timepoint, slidedist,upperiodtimelimit,dnperiodtimelimit);
      return;
    }
  }
}                               /* poptreesliderlimited */

void
init_holdgweight_array_for_hg(void)
{
  int li;
  holdgweight_array = static_cast<struct genealogy_weights *>             //points to an array of chains
          (malloc (nloci * sizeof (struct genealogy_weights))); 
  for (li = 0; li < nloci; li++)
  {
    init_genealogy_weights (&holdgweight_array[li]);
  }
  init_genealogy_weights (&holdallgweight);
  init_probcalc (&holdallpcalc);
}

void free_holdgweight_array_for_hg (void)
{
  int li;
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_array[li]);
  }
  XFREE(holdgweight_array);
  free_genealogy_weights (&holdallgweight);
  free_probcalc(&holdallpcalc);
}


void copyqpriorinfo(SET from[],SET to[])
{
  for (int i=0;i<numtreepops;i++)
    to[i] = from[i];
} //copyqpriorinfo


/* just copy prior and from/to info */ 
void copyimigpriors(int ci, int mode)
{
  int i;
  for (i=0;i<nummigrateparams;i++)
  {
    if (mode == 0)
    {
      holdimig[i].pr = C[ci]->imig[i].pr;
      //holdimig[i].md = C[ci]->imig[i].md;
    }
    else
    {
        C[ci]->imig[i].pr = holdimig[i].pr;
//        C[ci]->imig[i].md = holdimig[i].md;   should not need these 
    }
  }
} /*copyimigpriors */



/******* global functions ********/

void
init_change_poptree(char topologypriorinfostring[])
{
  int ci,i,j,k,m,npopsa;
  double p;
  char *c,*stringvalpos[(MAXPOPS_PHYLOGENYESTIMATION-1)*MAXPOPS_PHYLOGENYESTIMATION/2];
  init_holdgweight_array_for_hg();
  for (ci=0;ci<numchainspp;ci++)
  {
  /* set slide distance terms,  these are values to be multiple times the T prior to get the standard deviation of the slide
    target update rate is 0.2,  adjustor is 1.1,  
    starting valuue is 0.01,  with min of 0.00005 and max of 0.5
    these just seem to work
    did compare results with and without constant changing of these values,  and they seem about the same*/
    /*for (i=0;i<numchainstotal;i++) // used to put unique value in for chain, for debugging
    {
      if (beta[ci]==allbetas[i])
        break;
    }
    setupdatescalarinfo(&C[ci]->branchslideinfo,0.01,1.1,0.0001,0.5,beta[ci],100); */
    // updatescalarinfo, width, adjustval,min,max,target,numattemptscheck
    //setupdatescalarinfo(&C[ci]->branchslideinfo,0.1,1.1,0.0001,0.5,0.4,100); // lowered savmin because update rate is sometimes too low 8/29/2016
    //setupdatescalarinfo(&C[ci]->branchslideinfo,0.5,1.1,0.001,0.5,0.35,100); // lowered savmin because update rate is sometimes too low 8/29/2016
    
    //IMa3I
    //setupdatescalarinfo(&C[ci]->branchslideinfo,0.01,1.0,0.01,0.1,0.4,100);// should not update 
    //IMa3M
    //setupdatescalarinfo(&C[ci]->branchslideinfo,0.001,1.0,0.01,0.1,0.4,100);// should not update fix to 0.001

    //  11/15/2017  was not updating the branchslidescalar for some reason // notes from 5/3/2017 said it was not a good way 
    //so just using a setting that will not update because the adjustval is fixed at 1.0 
    setupdatescalarinfo(&C[ci]->branchslideinfo,0.01,1.0,0.01,0.1,0.4,100);// should not update 
    
    if ((C[ci]->ancplist = static_cast<int **> (malloc (numtreepops * sizeof (*C[ci]->ancplist)))) == NULL)
      IM_err (IMERR_MEM, "  ancplist malloc did not work.   numtreepops %d, step %d",  numtreepops, step);
    for (i = 0; i < numtreepops; i++)
    {
      if ((C[ci]->ancplist[i] =
          static_cast<int *> (malloc (npops * sizeof (*C[ci]->ancplist[i])))) == NULL)
        IM_err (IMERR_MEM, "  ancplist malloc did not work.   npops - i  %d, step %d",     npops - i, step);
    }
  }
  if (strlen(topologypriorinfostring) > 0)
  {
    usetopologypriors = 1;
    topologypoppairpriors = static_cast<double **> (malloc (npops * sizeof(double *)));
    for (i = 0;  i < npops; i++) 
    {
      topologypoppairpriors[i] = static_cast<double *> (malloc (npops* sizeof (double)));
    }
    topologypriors = static_cast<double *> (malloc (numpoptopologies* sizeof (double)));
    for (i=0;i<npops;i++) for (j=0;j<npops;j++) 
      topologypoppairpriors[i][j] = 0.0;
    c = &topologypriorinfostring[0];
    i = k = 0;
    while (*(c+i) != '\0')
    {
      if (isspace(*(c+i)))
        i += 1;
      else
      {
        stringvalpos[k] = (c+i);
        k += 1;
        while (isspace(*(c+i)) == 0)
          i += 1;
      }
    }
    m = 0;
    do
    {
      scanfval = sscanf (stringvalpos[m], "%d", &i);
      scanfval = sscanf (stringvalpos[m+1], "%d", &j);
      scanfval = sscanf (stringvalpos[m+2], "%lf", &p);
      m += 3;
      p = log(p);
      topologypoppairpriors[i][j] = p;
      topologypoppairpriors[j][i] = p;

    } while (m<k);
  }
  else
    usetopologypriors = 0;
  /*if (modeloptions[ADDGHOSTPOP]==1)
    npopsa = npops-1;
  else
    npopsa = npops; */

  npopsa = (modeloptions[ADDGHOSTPOP]==1) ? npops-1 : npops;
  if (usetopologypriors != 0)
  {
    for (k=0;k<numpoptopologies;k++)
    {
      p = 0.0;
      c = &alltreestrings[k][0];
      while (*c)
      {
        if (isdigit(*c))
        {
          i = atoi(c);  // atoi() works on c and everything after that looks like an integer up to the first char that does not 
          if (i== 2*(npopsa-1))  // at root node
            goto atend;
          if (isdigit(*(c + 1))) // for #'s with 2 digits
            c += 1;
          if (*(c+1)==',' && isdigit(*(c+2)))
          {
            c += 2;
            j = atoi(c); 
            p += topologypoppairpriors[i][j];
            if (isdigit(*(c + 1))) // for #'s with 2 digits
              c += 1;
          }
        }
        c += 1; // move one space
      }
    atend:
      topologypriors[k] = p;
    }
  }
  /*if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    holdpriors = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
  } */
    
} /* init_change_poptree */

void
free_change_poptree(void)
{
  int j,ci;
  free_holdgweight_array_for_hg ();
  for (ci=0;ci<numchainspp;ci++)
  {

    if (C[ci]->ancplist != NULL)  // hgstuff 
    {
      for (j = 0; j < numtreepops; j++)
      {
        XFREE (C[ci]->ancplist[j]);
      }
      XFREE (C[ci]->ancplist);
    }
  }
  if (usetopologypriors)
  {
    XFREE(topologypriors);
    for (j = 0; j < npops; j++) 
    {
      XFREE(topologypoppairpriors[j]);
    }
    XFREE(topologypoppairpriors);
  }
  /*if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    XFREE(holdpriors);

  } */
}

#define SCALESLIDEDISTBYPRIOR 1  /* uses prior for picking a slidedistance,  alternative is to use the roottime.  as of 3/29/2017 roottime has a bug and does not work */
int
change_poptree (int ci,int *trytopolchange, int *topolchange, int *trytmrcachange, int *tmrcachange, int topologychangeallowed)
{
  int i,li;
	int edge;
	int freededge;
	int oldsis,newsis;
	double holdt[MAXPERIODS],newt[MAXPERIODS];
	double roottime;
	double oldslidestdv,newslidestdv,slidedist,holdslidedist;
  double topologypriorratio = 0.0;
  double metropolishastingsratio;
  double priorratio,proposalratio;
  int accp;
  char holdpoptreestring[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  double upperiodtimelimit,dnperiodtimelimit;
  double temptime,temptmax;
  int tempb,oldtreenum;
  int ghostmoved = 0;
  
  double newpriorp, oldpriorp, propose_old_given_new, propose_new_given_old;
  double slideweightnum = 0.0, slideweightdenom = 0.0;
  int roottimechanges = 0;
  double mpriorratio,qpriorratio;
  double scaletemp[7] = {0.0003,0.001,0.003,0.01,0.03,0.1,0.3};// series of update widths, roughly uniform on log scale  
  static int scalei;
  SET holddescendantpops[MAXTREEPOPS];

  if (step==0)
    scalei = 0;
  else
  {
    scalei += 1;
    if (scalei > 6)
      scalei = 0;
  }
  /* change slide terms if needed */
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)
    oldtreenum = C[ci]->poptreenum;

#ifdef TURNONCHECKS
  //poptreeprint(ci);
    //gtreeprint(ci,0);
checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  checkprobs(ci,-1);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    }
checkgenealogyweights(ci);
#endif //TURNONCHECKS
  *topolchange = 0;
	numtreepops = 2*npops - 1;
// initialize and make copies structures that hold quantities for calculating prob of genealogy
  for (li=0;li <nloci;li++)
  {
    copy_treeinfo (&holdgweight_array[li],&C[ci]->G[li].gweight);
  }
  copy_treeinfo (&holdallgweight, &C[ci]->allgweight);
  copy_probcalc (&holdallpcalc, &C[ci]->allpcalc);
 
  storegenealogystats_all_loci_hg (ci, 0);
	copy_poptree(C[ci]->poptree,&poptreehold[0]);
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    copyimigpriors(ci,0);
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      copyqpriorinfo(C[ci]->descendantpops,holddescendantpops);
  }
  strcpy(&holdpoptreestring[0],C[ci]->chainpoptreestring);
  //checkupdatescalarer(&C[ci]->branchslideinfo);  as of 5/3/2017  not updating this 
	for (i = 0; i < lastperiodnumber; i++)
	  holdt[i] = C[ci]->tvals[i];
	roottime = C[ci]->tvals[npops-2];
  do
  {
    edge = randposint (numtreepops);
  } while (C[ci]->poptree[edge].down == -1);
  freededge = C[ci]->poptree[edge].down;
  if ((oldsis = C[ci]->poptree[freededge].up[0]) == edge)
    oldsis = C[ci]->poptree[freededge].up[1];
  if (topologychangeallowed)
  {
    if (SCALESLIDEDISTBYPRIOR)
      //oldslidestdv = C[ci]->branchslideinfo.updatescalarval * T[numsplittimes-1].pr.max;
      oldslidestdv = scaletemp[scalei] * T[numsplittimes-1].pr.max;
      //T[npops-2].pr.max; // period npops-2 is the last period with a splitting time 
    else
      oldslidestdv = C[ci]->branchslideinfo.updatescalarval * roottime; 
  }
  else
  {
    assert(C[ci]->poptree[edge].e != 0);
    if (C[ci]->poptree[edge].e > 0)
      tempb = C[ci]->poptree[edge].e - 1;
    else
      tempb = numsplittimes;
    assert (tempb >= 0);
    if (tempb == 0)
      upperiodtimelimit = 0.0;
    else
      upperiodtimelimit = C[ci]-> tvals[tempb-1];
    if (tempb < numsplittimes-1) 
      dnperiodtimelimit = C[ci]-> tvals[tempb+1];
    else
      dnperiodtimelimit = T[numsplittimes-1].pr.max;
    oldslidestdv = C[ci]->branchslideinfo.updatescalarval * (dnperiodtimelimit - upperiodtimelimit); // using possibly updated updatescalarval
  }
  holdslidedist = slidedist = normdev (0.0, oldslidestdv);// pick a slidedistance 
  /* update slidedistance scaler */
  
  /*if (roottimechanges)
  {
    slideweightdenom = -log (normprob (0.0, oldslidestdv, holdslidedist)); // term for updated value goes in denominator // not needed with this kind of slide update, where slidestdv is the same both for the update and the reverse
  }
  else
  {
    slideweightdenom = 0.0;
  } */
  // join the sister and the down branches at the point where edge used to connect, this frees up the down branch
  poptreejoinsisdown (C[ci]->poptree, oldsis, &roottimechanges);
  newsis = oldsis;
#ifdef TURNONCHECKS
  //poptreeprint(ci);
#endif //TURNONCHECKS
  if (topologychangeallowed)
  {
    /*
      if a priorfile is used,  then the subtree with ingroups might have a maximum that is less than the maximum for where the ghost outgroup connects 
      in this case the time of connection of an edge can come in older than the allowed maximum
      but we want to allow the simulation of the edge contact point to be able to move around on the tree (e.g. move down and then up a different edge,  not just reflect)
      so this code here allows for that
      Should only be needed when modeloptions[ADDGHOSTPOP]==1 and when calcoptions[LOADPRIORSFROMFILE]==1 and when modeloptions[POPTREETOPOLOGYUPDATE]==1
    */
    if (modeloptions[ADDGHOSTPOP]==1 && calcoptions[LOADPRIORSFROMFILE]==1 && modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      do
      {
        newsis = oldsis;
        temptime = C[ci]->poptree[edge].time;
        temptmax = T[findtindex(ci,temptime)].pr.max;
        poptreeslider (C[ci]->poptree,edge, &newsis, &temptime,&slidedist);
        assert (slidedist==0.0);
      } while (temptime > temptmax);
      C[ci]->poptree[edge].time = temptime;
    }
    else
    {
      poptreeslider (C[ci]->poptree,edge, &newsis, &C[ci]->poptree[edge].time,&slidedist);
      assert (slidedist==0.0);
    }
  }
  else   // only slide within range of flanking splitting times 
  {
    poptreesliderlimited (C[ci]->poptree,edge, &newsis, &(C[ci]->poptree[edge].time),&slidedist,
      &upperiodtimelimit,&dnperiodtimelimit);
  }
  poptreesplitsisdown (C[ci]->poptree,edge, freededge, newsis,&roottimechanges);
#ifdef TURNONCHECKS
 //poptreeprint(ci);
checkpoptree(ci,2);
#endif //TURNONCHECKS

  //renumberpopnodes_hold(ci,&newt[0]);
  renumberpopnodes(ci,&newt[0]);
  
#ifdef TURNONCHECKS
//poptreeprint(ci);
checkpoptree(ci,1);
#endif //TURNONCHECKS
  
  
  if (roottimechanges)
  {
    assert (roottime != newt[npops-2]);
    roottime = newt[npops-2];
    *trytmrcachange = 1;
    if (topologychangeallowed && SCALESLIDEDISTBYPRIOR==0) // if topoogychangeallowed == 0  then slidestdv comes from the boundaries and will not change
    {
      slideweightdenom = -log (normprob (0.0, oldslidestdv, holdslidedist));
      newslidestdv = C[ci]->branchslideinfo.updatescalarval * roottime; 
      slideweightnum = log (normprob (0.0, newslidestdv,holdslidedist));
    }
    else
      slideweightnum = slideweightdenom = 0.0;

  }
  else
  {
    *trytmrcachange = 0;
    slideweightnum = slideweightdenom = 0.0;
  }
  assert(newt[npops-2] < T[npops-2].pr.max);
  C[ci]->chainpoptreestring[0] = 0;
  poptreewrite (ci, C[ci]->chainpoptreestring);
  *topolchange = (strcmp(&holdpoptreestring[0],C[ci]->chainpoptreestring) != 0);
#ifdef TURNONCHECKS
    if (topolchange)
    {
    //poptreeprint(ci);
    checkpoptree(ci,1);
    }
#endif //TURNONCHECKS



	for (i = 0; i < lastperiodnumber; i++)
	  C[ci]->tvals[i] = newt[i];

#ifdef TURNONCHECKS
//poptreeprint(ci);
 //poptreewrite (ci, C[ci]->chainpoptreestring);
#endif //TURNONCHECKS
 
  if (*topolchange) // poptree topology has changed 
  {
    
    assert (topologychangeallowed == 1);
    assert (modeloptions[POPTREETOPOLOGYUPDATE] == 1);
    C[ci]->poptreenum = getpoptreestringnum(C[ci]->chainpoptreestring);

    ghostmoved = (C[ci]->poptreenum ==-1);  // ghost is not allowed to move,  leads to rejection 
    assert(ghostmoved == 0 || (modeloptions[ADDGHOSTPOP] && ghostmoved==1));
#ifdef DEBUG
    if (C[ci]->poptreenum == oldtreenum)
    {
      //poptreeprint(ci);
      IM_err(IMERR_MISCELLANEOUS,"ci,step,newpoptreenum %d oldtreenum %d poptreestring %s",ci,step,C[ci]->poptreenum,oldtreenum,C[ci]->chainpoptreestring);
    }
#endif //DEBUG
    assert(C[ci]->poptreenum != oldtreenum);
   
    if (ghostmoved == 0)  // this saves some time if modeloptions[ADDGHOSTPOP]==1  as proposal will get rejected if ghostmoved == 1
    {
      *trytopolchange = 1;
      fillplist(ci);
      fillancplist(ci);
      set_iparam_poptreeterms (ci);
      assert (C[ci]->poptreenum >= 0 && C[ci]->poptreenum < numpoptopologies);
      poptopologyproposedlist[C[ci]->poptreenum] += 1;
      if (usetopologypriors)
        topologypriorratio = topologypriors[C[ci]->poptreenum] - topologypriors[oldtreenum];
      mpriorratio = qpriorratio = 0.0;
      if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      {
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          mpriorratio = 0.0;
          for (i=0;i<nummigrateparampairs;i++)
          {
            //mpriorratio += C[ci]->mltorhpriors[holdcurrentpairpos[i]] - C[ci]->mltorhpriors[C[ci]->currentpairpos[i]];

            /* C[ci]->imig[2*i].descstr and C[ci]->imig[2*i+1].descstr should be the same */ 
            //mpriorratio += holdimig[2*i].pr.expomean - getval(C[ci]->imig[2*i].descstr,C[ci]->mltorhpriors);
            mpriorratio += holdimig[2*i].pr.expomean - C[ci]->imig[2*i].pr.expomean;
            assert (holdimig[2*i].md.from == holdimig[2*i+1].md.to);
            assert (holdimig[2*i].md.to == holdimig[2*i+1].md.from);

            //mpriorratio += C[ci]->mrtolhpriors[holdcurrentpairpos[i]] - C[ci]->mrtolhpriors[C[ci]->currentpairpos[i]];
            //mpriorratio += holdimig[2*i+1].pr.expomean - getval(C[ci]->imig[2*i+1].descstr,C[ci]->mrtolhpriors);
            mpriorratio += holdimig[2*i+1].pr.expomean - C[ci]->imig[2*i+1].pr.expomean;

          }
          mpriorratio /= expo_m_mean;
        }
        else
        {
            mpriorratio = 1.0;
          	 for (i=0;i<nummigrateparampairs;i++)
            {
              /* C[ci]->imig[2*i].descstr and C[ci]->imig[2*i+1].descstr should be the same */ 
              //mpriorratio *= C[ci]->mltorhpriors[C[ci]->currentpairpos[i]]/C[ci]->mltorhpriors[holdcurrentpairpos[i]];
              //mpriorratio *= getval(C[ci]->imig[2*i].descstr,C[ci]->mltorhpriors)/holdimig[2*i].pr.max;
              mpriorratio *= C[ci]->imig[2*i].pr.max/holdimig[2*i].pr.max;
              //mpriorratio *= C[ci]->mrtolhpriors[C[ci]->currentpairpos[i]]/C[ci]->mrtolhpriors[holdcurrentpairpos[i]];
              //mpriorratio *= getval(C[ci]->imig[2*i+1].descstr,C[ci]->mrtolhpriors)/holdimig[2*i+1].pr.max;
              mpriorratio *= C[ci]->imig[2*i+1].pr.max/holdimig[2*i+1].pr.max;
              assert (holdimig[2*i].md.from == holdimig[2*i+1].md.to);
              assert (holdimig[2*i].md.to == holdimig[2*i+1].md.from);
            }
            mpriorratio = -log(mpriorratio); 
        }
        qpriorratio = 1.0;
        for (i=0;i<numtreepops;i++)
        {
          qpriorratio *= C[ci]->qhpriors[(int) C[ci]->descendantpops[i]]/C[ci]->qhpriors[(int) holddescendantpops[i]];
        }
        qpriorratio = -log(qpriorratio); 
      }
    }   
  }
  else
    *trytopolchange = 0;

  if (ghostmoved == 0)   // this saves some time if modeloptions[ADDGHOSTPOP]==1 as proposal will get rejected if ghostmoved == 1
  {
    setzero_genealogy_weights (&C[ci]->allgweight);
    for (li=0;li<nloci;li++)
    {
      makegenealogy_from_hiddengenealogy(ci,li);
      // determine all the weights needed for calculating the probability of the genealogy
      setzero_genealogy_weights (&C[ci]->G[li].gweight);
#ifdef TURNONCHECKS
      //poptreeprint(ci);
  checkpoptree(ci,0);
#endif //TURNONCHECKS
      int treeweightcallcode = 3;  // a debugging code,  if treeweight has an error results are written to an output file with this code 
      int tw = treeweight (ci, li, treeweightcallcode);
      if (tw < 0)
        goto rejectjump; //  very awkward,   not sure if we get here and if this is needed as of 8/11/2017
      sum_treeinfo (&C[ci]->allgweight, &C[ci]->G[li].gweight);
    }
    
    if (*topolchange) // bad bug was here,  if poptree changed, and using held weights and pcalc vals,  you get the wrong probg value sometimes 
      initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
    else
      integrate_tree_prob (ci, &C[ci]->allgweight,&holdallgweight,&C[ci]->allpcalc, &holdallpcalc);
    C[ci]->allpcalc.probhgg = 0.0;
    for (li=0;li<nloci;li++)
    {
      C[ci]->G[li].hgprob =  prob_hg_given_g(ci,li);
        C[ci]->allpcalc.probhgg +=  C[ci]->G[li].hgprob;
    }
    newpriorp = C[ci]->allpcalc.probg + C[ci]->allpcalc.probhgg;
    oldpriorp = holdallpcalc.probg + holdallpcalc.probhgg;
    propose_old_given_new = slideweightnum;
    propose_new_given_old = slideweightdenom;
    if (usetopologypriors)
    {
       newpriorp += topologypriors[C[ci]->poptreenum];
       oldpriorp += topologypriors[oldtreenum];
    }
    priorratio  = newpriorp - oldpriorp;
    proposalratio = propose_old_given_new - propose_new_given_old;
    if (*topolchange)
    {
      priorratio += mpriorratio + qpriorratio;
    }

    /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
    /* under thermodynamic integration, only ratios of likelihoods get raised to beta. 
      but in this update the ratio of likelihoods is always 1,  so beta is not used with thermodynamic integration
    */
#ifdef TURNONCHECKS
    checkdetailedbalance(0.0,0.0,newpriorp,oldpriorp,propose_old_given_new, propose_new_given_old,beta[ci]);
#endif //TURNONCHECKS
    if (calcoptions[CALCMARGINALLIKELIHOOD]) 
    {
      metropolishastingsratio = priorratio + proposalratio;
    }
    else
    {
      metropolishastingsratio = beta[ci] * priorratio + proposalratio;
    }
  }
  else
  {
    metropolishastingsratio = DBL_MAX*1.1; // make inf,  should only get here if ghostmoved == 1
    assert(ghostmoved==1 && modeloptions[ADDGHOSTPOP]==1);
  }  
  if (metropolishastingsdecide(metropolishastingsratio,(ghostmoved==0)))
  {
#ifdef TURNONCHECKS
    //gtreeprint(ci,0);
  checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
#endif
    /* accept the update */
    accp = 1;
    *tmrcachange = roottimechanges;
    if (*topolchange)
    {
      totaltopolupdates += 1;
      if (beta[ci] == 1.0)
        chain0topolupdates += 1;
    }
    C[ci]->branchslideinfo.recentaccp += 1;
    C[ci]->branchslideinfo.allaccp += 1;

#ifdef TURNONCHECKS
    //gtreeprint(ci,0);
checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  checkprobs(ci,-1);
  checkgenealogyweights(ci);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    checkprobs(ci,li);
    }
#endif //TURNONCHECKS

  }
  else
  /* reject the update */
  {
  rejectjump:
    accp = 0;
    // put everything back
    copy_poptree(&poptreehold[0],C[ci]->poptree);
	   for (i = 0; i < lastperiodnumber; i++)
	     C[ci]->tvals[i]  = holdt[i];
    copy_treeinfo (&C[ci]->allgweight,&holdallgweight);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc);
    storegenealogystats_all_loci_hg (ci, 1);
    if (*topolchange)
    {
      assert(modeloptions[POPTREETOPOLOGYUPDATE] == 1);
      strcpy(C[ci]->chainpoptreestring,&holdpoptreestring[0]);
      C[ci]->poptreenum = getpoptreestringnum(C[ci]->chainpoptreestring);
      fillplist(ci);
      set_iparam_poptreeterms (ci);
      fillancplist(ci);
      if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      {
        copyimigpriors(ci,1);
      } 
    }
    for (li=0;li<nloci;li++)
    {
      makegenealogy_from_hiddengenealogy(ci,li);
      copy_treeinfo (&C[ci]->G[li].gweight,&holdgweight_array[li]);
#ifdef TURNONCHECKS
      checkprobs(ci,li);
      checkgenealogy(ci,li,0);
#endif //TURNONCHECKS
    }
    *topolchange = 0;
    *tmrcachange = 0;

#ifdef TURNONCHECKS
    //gtreeprint(ci,0);
checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  checkprobs(ci,-1);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    checkprobs(ci,li);
    }
    checkgenealogyweights(ci);
#endif //TURNONCHECKS
  // resetupdatescalarer(&C[ci]->branchslideinfo);
  }
//printf("step %d chain %d treenum %d accp %d uni %.4lf\n",step,ci,C[ci]->poptreenum,accp,uniform());
  return accp;
}