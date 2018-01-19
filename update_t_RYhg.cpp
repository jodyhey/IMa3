/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS

#include "ima.hpp"
#include "update_gtree_common.hpp"

/*********** LOCAL STUFF **********/

static struct genealogy_weights holdallgweight_t_RY_hg;
static struct genealogy_weights holdgweight_t_RY_hg[MAXLOCI];
static struct probcalc holdallpcalc_t_RY_hg;
static int largestsamp;
static int **numhg;

struct  mighgstruct
  { 
    int a[MIGARRAYSIZE];
  } **mighgs;

static double beforesplithg (int tnode, double oldt, double newt, double tau_u, double ptime);
static double aftersplithg (int tnode, double oldt, double newt, double tau_d, double ptime);
static void storegenealogystats_all_loci_hg (int ci, int mode);

static void getmighgpositions(struct edge *gtree, int *numhgloc, int mighgs[]);
/********* LOCAL FUNCTIONS **************/

#ifdef TURNONCHECKS
void checkmigmatch(int ci)
{
  int li,i,j;
  struct edge *gtree;
  for (li=0;li<nloci;li++)
  {
    gtree = C[ci]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++) if (C[ci]->G[li].gtree[i].down != -1)
    {
      for (j=0;j<numhg[li][i];j++)
      {
        if (mighgs[li][i].a[j] >= 0)
          assert (gtree[i].mig[mighgs[li][i].a[j]].mt==gtree[i].mighg[j].mt);
      }
    }
  }

}
#endif //TURNONCHECKS

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
}                               // storegenealogystats  

double
aftersplithg (int tnode, double oldt, double newt, double tau_d, double ptime)
{
  if (tnode == lastperiodnumber - 1)
  {
     return ptime + newt - oldt;
  }
  else
  {
    return tau_d - (tau_d - newt) * (tau_d - ptime) / (tau_d - oldt);
  }
}

double
beforesplithg (int tnode, double oldt, double newt, double tau_u, double ptime)
{
  if (tnode == 0)
  {
    return ptime * newt / oldt;
  }
  else
  {
    return tau_u + (ptime - tau_u) * (newt - tau_u) / (oldt - tau_u);
  }
}
/* go over all mighg events,
  record how many mighg events there are,  and which ones have corresponding mig
  events,  and if they do,  where in the mig array the corresponding mig event is */
void getmighgpositions(struct edge *gtree, int *numhgloc, int mighgs[])
{
  int j, jhg, mdone;
  j = jhg = 0;// j index mig[], jhg indexes mighg[]
  mdone = 0;
  while (mdone == 0)
  {
    if (gtree->mig[j].mt > -0.5)
    {
      while (gtree->mighg[jhg].mt < gtree->mig[j].mt && gtree->mighg[jhg].mt > -0.5) // go thru mighg[] events without corresponding mig[] events
      {
        mighgs[jhg] = -1;
        jhg += 1;
      }
      if (gtree->mighg[jhg].mt == gtree->mig[j].mt) // should always be true at this point because we have to come up to the matching event in mighg[]
      {
        mighgs[jhg] = j;
        jhg += 1;
      }
      else
      {
        assert(0);// should not get here 
      }
      j+= 1;
    }
    else
    {
      while (gtree->mighg[jhg].mt > -0.5)
      {
        mighgs[jhg] = -1;
        jhg += 1;
      }
      mdone = 1;
    }
  }
  *numhgloc = jhg;
  assert(jhg <= MIGMAX );
}  //getmighgpositions

/*************GLOBAL FUNCTIONS ******************/


void
init_t_RYhg (void)
{
  int li, i,j;
  init_genealogy_weights (&holdallgweight_t_RY_hg);
  for (li = 0; li < nloci; li++)
    init_genealogy_weights (&holdgweight_t_RY_hg[li]);
  init_probcalc (&holdallpcalc_t_RY_hg);
  for (largestsamp = 0, j = 0; j < nloci; j++)
  {
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
  }
  numhg = alloc2Dint (nloci, largestsamp);

  if ((mighgs = static_cast<struct  mighgstruct **> (malloc (nloci * sizeof (*mighgs)))) == NULL)
    IM_err (IMERR_MEM, "  malloc did not work.  step %d",step);
  for (i = 0; i < nloci; i++)
    if ((mighgs[i] = static_cast<struct  mighgstruct *> (malloc (largestsamp * sizeof (*mighgs[i])))) == NULL)
      IM_err (IMERR_MEM,"  malloc did not work. nloci %d,  largestsamp %d, step %d",nloci, largestsamp, step);
 }                               // init_t_RYhg

void
free_t_RYhg (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_RY_hg);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_RY_hg[li]);
  }
  free_probcalc (&holdallpcalc_t_RY_hg);
  orig2d_free2D ((void **) numhg, nloci);
  orig2d_free2D ((void **) mighgs, nloci);

}                               // free_t_RYhg

/* 
Notes on changet_RYhg() (copied from changet_RY()) implements updating of Rannala and Yang (2003)   

This version is for hidden genealogies

This application is pretty much the same as theirs - changing times on a species tree that contains a gene tree. 
The big difference is that IM includes migration.   This means that we have to count migration events and change 
migration times in the same way as we change coalescent times and count how many get changed. 

R&Y also only change coalescent times in populations affected by a changed t.  But because we have migration
there is more entanglement between populations.  It seems best to change all times within an interval that is 
affected by a changing t regardless of which populations the migration events are between.

in R&Y usage
u (upper) means older, deeper in the genealogy
l (lower)  means younger, more recent in the genealogy

In R&Y 
Tu next older split time
Tl  next more recent split time
T - current split time
T*  - proposed time 
t - time of some event 
t* - time of that event after update

If  Tu> t > T
t* = Tu - (Tu-t) (Tu-t*)/(Tu-T)

If Tl< t<= T
t* = Tl + (t-tl) (T*-Tl)/(T-Tl)

nl = number of nodes with Tl< t < T
ml = number of migration events with Tl< t < T

nu = number of nodes with Tu> t > T
mu = number of migration events with Tu> t > T

MH criteria  
p(X|G*,t*)p(G*,t*) (Tu-T*)^(nu+mu) (T*-Tl)^(nl+ml)
---------------------------------------------------
p(X|G,t)p(G,t)      (Tu-T)^(nu+mu) (T-Tl)^(nl+ml)

but this causes confusion with hiddenoptions[] usage in which 
u means upper - more recent. 

so here we use  u  for upper (meaning more recent)
use d for down  (older)

tau current split time
tau*  new value 
tau_d - older node time (i.e. time of next oldest node - deeper in the genealogy)
tau_u - more recent node time (i.e. time of next youngest node - more recent in time)
tau_d  > tau > tau_u

if tau is the first node,  then tau_u = 0. 
if tau is the last node, then tau_d = infinity

for an event at time t where tau_u < t < tau_d

if t > tau  see aftersplithg()  t is older than the current split time 
t* = tau_d - (tau_d - tau*)(tau_d - t)/(tau_d-tau)  (A7)

if t <= tau  see beforesplit()
t* = tau_u + (tau* - tau_u)(t - tau_u)/(tau - tau_u) (A8)

m is the number of events moved using A7,  n is the number moved using A8

then Hastings term for the update is:
 tau_d - tau*      tau* - tau_u 
(------------)^m  (------------)^n
 tau-u - tau        tau - tau_u

 For IM,  we use the same except m and n include both includes migation and coalescent events

For edges where downtime < tau_d || uptime > tau_d  are not involved in the update
For any edge that spends zero time in either splitpop  or the ancestral pop, during the tau_u/tau_d interval
it is  possible to not update the coalescent time or migration times of 

The difficulty is that it is possible that an uptime for one edge gets moved because the times on one of its daughter edges got moved. 
This means that for checking about skipping an edge, because it is not involved in any
population associated with the splittin time update
we need to check entire groups of branches that descend from 
an edge that crosses the tau_u boundary. 

*/

/* let u refer to the more recent time  and d to the older time  */
int
changet_RYhg (int ci, int timeperiod)    // after Rannala and Yang (2003)  - rubberband method
{
  double metropolishastingsratio, newt, oldt;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum,
    temppdg;
  double t_u_hterm, t_d_hterm,t_u_htermnum,t_u_htermdenom,t_d_htermnum, t_d_htermdenom;
  double likelihoodratio,priorratio,proposalratio;
  int li, i, j, ecd, ecu, emd, emu, ai, ui;
  struct genealogy *G;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;
  double oldlikelihood,newlikelihood;
  double newpriorp, oldpriorp, propose_old_given_new, propose_new_given_old;

 checkupdatescalarer(&C[ci]->RYwidthinfo[timeperiod]);
#ifdef TURNONCHECKS
  int  cecd = 0, cecu = 0, cemd = 0, cemu = 0;
    //gtreeprint(ci,0);
   checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  checkprobs(ci,-1);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    }
#endif //TURNONCHECKS

  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];

  /* this may be a bug,  i.e. using TIMEMAX,  try using the prior max instead    restored this,  does not work to change it*/
  t_d = (timeperiod == (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  //t_d =  (timeperiod == (lastperiodnumber - 1)) ? T[timeperiod].pr.max : C[ci]->tvals[timeperiod + 1];

  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];

  newt = getnewt (t_u_prior, t_d_prior, oldt,C[ci]->RYwidthinfo[timeperiod].updatescalarval);
  t_u_htermnum = (newt - t_u);
  t_u_htermdenom = (oldt - t_u);
  t_u_hterm = t_u_htermnum /t_u_htermdenom;
  if (timeperiod == lastperiodnumber - 1) 
  {
    t_d_htermnum = t_d_htermdenom = 1.0;
    t_d_hterm = 1.0;
  }
  else 
  {
    t_d_htermnum = (t_d - newt);
    t_d_htermdenom = (t_d - oldt);
    t_d_hterm = t_d_htermnum/t_d_htermdenom;
  }  

  C[ci]->poptree[C[ci]->droppops[timeperiod+1][0]].time = newt;
  C[ci]->poptree[C[ci]->droppops[timeperiod+1][1]].time = newt;
  C[ci]->tvals[timeperiod] = newt;
  copy_treeinfo (&holdallgweight_t_RY_hg, &C[ci]->allgweight);  
#ifdef TURNONCHECKS
  compare_genealogy_weights (&holdallgweight_t_RY_hg, &C[ci]->allgweight);  
#endif
  copy_probcalc (&holdallpcalc_t_RY_hg, &C[ci]->allpcalc);
  storegenealogystats_all_loci_hg (ci, 0);
  pdgoldsum = C[ci]->allpcalc.pdg;
  setzero_genealogy_weights (&C[ci]->allgweight);
  ecd = ecu = emd = emu = 0;
  pdgnewsum = 0;
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_RY_hg[li], &G->gweight);
    for (i = 0; i < L[li].numlines; i++)
    {
      if (gtree[i].down != -1)
      {
        if (gtree[i].time <= oldt && gtree[i].time > t_u) 

        {
          gtree[i].time =
            beforesplithg (timeperiod, oldt, newt, t_u, gtree[i].time);
          assert (gtree[i].time != newt);
          ecu++;
        }
        else
        {
          if (gtree[i].time >= oldt && gtree[i].time < t_d)
          {
            gtree[i].time =
              aftersplithg (timeperiod, oldt, newt, t_d, gtree[i].time);
            assert (gtree[i].time != newt);
            ecd++;
          }
          //else  do not change the time
        }
        // identify the locations of all mig[] events in the mighhg[] array
        getmighgpositions(&gtree[i],&numhg[li][i],mighgs[li][i].a);
        // reset the migration events 
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          assert (gtree[i].mig[j].mt <TIMEMAX);
          if (gtree[i].mig[j].mt <= oldt && gtree[i].mig[j].mt > t_u) 
          {
            gtree[i].mig[j].mt =
              beforesplithg (timeperiod, oldt, newt, t_u,
                           gtree[i].mig[j].mt);
            emu++;
          }
          else
          {
            assert (oldt <TIMEMAX);
            if (gtree[i].mig[j].mt >= oldt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplithg (timeperiod, oldt, newt, t_d, 
                            gtree[i].mig[j].mt);
              emd++;
            }
            // else no need to change the time
          }
          j++;
        }
        // now do hidden genealogy migrations using mighgs to ensure we have identical times for events that are in both mig[] and mighg[]
        for (j=0;j<numhg[li][i];j++)
        {
          if (mighgs[li][i].a[j] >= 0)
          {
            gtree[i].mighg[j].mt = gtree[i].mig[mighgs[li][i].a[j]].mt;
          }
          else
          {
            assert (gtree[i].mighg[j].mt <TIMEMAX);
            if (gtree[i].mighg[j].mt <= oldt && gtree[i].mighg[j].mt > t_u) 
            {
              gtree[i].mighg[j].mt =
                beforesplithg (timeperiod, oldt, newt, t_u,
                             gtree[i].mighg[j].mt);
              emu++;
            }
            else
            {
              assert (oldt <TIMEMAX);
              if (gtree[i].mighg[j].mt >= oldt && gtree[i].mighg[j].mt < t_d)
              {
                gtree[i].mighg[j].mt = aftersplithg (timeperiod, oldt, newt, t_d, gtree[i].mighg[j].mt);
                emd++;
              }
            }
          } 
        }  // mighg work
      }

    }
    G->roottime = gtree[gtree[G->root].up[0]].time;
    // had a bug here,  roottime was not getting close enough to what it should have been 

    setzero_genealogy_weights (&G->gweight);
    int treeweightcallcode = 6;  // a debugging code,  if treeweight has an error results are written to an output file with this code 
    int tw = treeweight (ci, li,treeweightcallcode);
    if (tw < 0)
      goto rejectjump;
    sum_treeinfo (&C[ci]->allgweight, &G->gweight);
    ai = 0;
    ui = L[li].uii[ai];

    switch (L[li].model)
    {
      assert (pdgnew[ui] == 0);
    case HKY:
        temppdg = pdgnew[ui] =
          likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
      break;
    case INFINITESITES:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      break;
    case STEPWISE:
      temppdg = 0;
      for (; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    case JOINT_IS_SW:
      temppdg = pdgnew[ui] = likelihoodIS (ci, li, G->uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        ui = L[li].uii[ai];
        assert (pdgnew[ui] == 0);
        pdgnew[ui] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
        temppdg += pdgnew[ui];
      }
      break;
    }
    pdgnewsum += temppdg;
  }
#ifdef TURNONCHECKS
    //getmighgpositions(&gtree[i],&numhg[li][i],mighgs[li][i].a);
    checkmigmatch(ci);
#endif
  assert (!ODD (ecd));
  assert (!ODD (ecu));
  ecd /= 2;
  ecu /= 2;
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY_hg, &C[ci]->allpcalc, &holdallpcalc_t_RY_hg);   // try enforcing full calculation
  //initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
#ifdef TURNONCHECKS
  //  gtreeprint(ci,0);
 //   poptreeprint(ci);
  checkgenealogy(ci,0,0);
#endif //TURNONCHECKS
  C[ci]->allpcalc.probhgg = 0.0;
  for (li=0;li<nloci;li++)
  {
    C[ci]->G[li].hgprob =  prob_hg_given_g(ci,li);
      C[ci]->allpcalc.probhgg +=  C[ci]->G[li].hgprob;
  }
  oldlikelihood  = pdgoldsum;
  newlikelihood  = pdgnewsum;
  newpriorp = C[ci]->allpcalc.probhgg + C[ci]->allpcalc.probg;
  oldpriorp =  holdallpcalc_t_RY_hg.probhgg + holdallpcalc_t_RY_hg.probg;

  propose_old_given_new = (ecd + emd) * log (t_d_htermnum) + (ecu + emu) * log (t_u_htermnum);
  propose_new_given_old = (ecd + emd) * log (t_d_htermdenom) + (ecu + emu) * log (t_u_htermdenom);

  likelihoodratio = newlikelihood - oldlikelihood;
  priorratio = newpriorp - oldpriorp;
  proposalratio = propose_old_given_new - propose_new_given_old;
#ifdef TURNONCHECKS
    //propose_new_given_old = propose_old_given_new = 0.0;
    checkdetailedbalance(newlikelihood,oldlikelihood,newpriorp,oldpriorp,propose_old_given_new, propose_new_given_old,beta[ci]);
#endif //TURNONCHECKS
  //priorratio = ( C[ci]->allpcalc.probg - holdallpcalc_t_RY_hg.probg) +  (C[ci]->allpcalc.probhgg - holdallpcalc_t_RY_hg.probhgg);
  //proposalratio = (ecd + emd) * log (t_d_hterm) + (ecu + emu) * log (t_u_hterm);
/* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */

  if (calcoptions[CALCMARGINALLIKELIHOOD]) 
    {
    metropolishastingsratio = beta[ci] * likelihoodratio + priorratio + proposalratio;
  }
  else
  {
      metropolishastingsratio = beta[ci] * (likelihoodratio + priorratio) + proposalratio;
  }
  if (metropolishastingsdecide(metropolishastingsratio,1))
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].pdg = 0;
      for (ai = 0; ai < L[li].nlinked; ai++)
      {
        C[ci]->G[li].pdg_a[ai] = pdgnew[L[li].uii[ai]];
        C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
      }
      if (L[li].model == HKY)
      {
        storescalefactors (ci, li);
        copyfraclike (ci, li);
      }
#ifdef TURNONCHECKS
      checkprobs(ci,li);
#endif //TURNONCHECKS
    }
    C[ci]->allpcalc.pdg = pdgnewsum;

    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;

    C[ci]->RYwidthinfo[timeperiod].recentaccp += 1;
    C[ci]->RYwidthinfo[timeperiod].allaccp += 1;
#ifdef TURNONCHECKS
    //gtreeprint(ci,0);
checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  checkprobs(ci,-1);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    }
#endif //TURNONCHECKS
    return 1;
  }
  else
  {
  rejectjump:
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY_hg);
    C[ci]->poptree[C[ci]->droppops[timeperiod+1][0]].time = oldt;
    C[ci]->poptree[C[ci]->droppops[timeperiod+1][1]].time = oldt;
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY_hg);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storegenealogystats_all_loci_hg (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY_hg[li]);
      for (i = 0; i < L[li].numlines; i++)
      {
        if (gtree[i].down != -1)
        {
          if (gtree[i].time <= newt && gtree[i].time > t_u) 
          {
            gtree[i].time =
              beforesplithg (timeperiod, newt, oldt, t_u, gtree[i].time);
#ifdef TURNONCHECKS
            cecu++;
#endif
          }
          else
          {
            if (gtree[i].time >= newt && gtree[i].time < t_d)
            {
              gtree[i].time =
                aftersplithg (timeperiod, newt, oldt, t_d, gtree[i].time);
#ifdef TURNONCHECKS
            cecd++;
#endif

            }
          }
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            if (gtree[i].mig[j].mt <= newt && gtree[i].mig[j].mt > t_u) 
            {
              gtree[i].mig[j].mt =
                beforesplithg (timeperiod, newt, oldt, 
                             t_u, gtree[i].mig[j].mt);
#ifdef TURNONCHECKS
              cemu++;
#endif

            }
            else if (gtree[i].mig[j].mt >= newt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplithg (timeperiod, newt, oldt, t_d, 
                            gtree[i].mig[j].mt);
#ifdef TURNONCHECKS
              cemd++;
#endif

            }
            j++;
          }
          // now do mighg
          for (j=0;j<numhg[li][i];j++)
          {
            if (mighgs[li][i].a[j] >= 0)
            {
              gtree[i].mighg[j].mt = gtree[i].mig[mighgs[li][i].a[j]].mt;
            }
            else
            {
              if (gtree[i].mighg[j].mt <= newt && gtree[i].mighg[j].mt > t_u) 
              {
                gtree[i].mighg[j].mt =
                  beforesplithg (timeperiod, newt, oldt, 
                               t_u, gtree[i].mighg[j].mt);
#ifdef TURNONCHECKS
                cemu++;
#endif                
              }
              else if (gtree[i].mighg[j].mt >= newt && gtree[i].mighg[j].mt < t_d)
              {
                gtree[i].mighg[j].mt =
                  aftersplithg (timeperiod, newt, oldt, t_d, 
                              gtree[i].mighg[j].mt);
#ifdef TURNONCHECKS
                cemd++;
#endif
              }
            }
          }
        }
      }
//        assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
      G->roottime = gtree[gtree[G->root].up[0]].time;// 4/7/2017  roottime was not being reset with rejection,  a bug 
    }


#ifdef TURNONCHECKS// check to make sure that the same number of changes made putting things back after rejection
    cecu /=2;
    cecd /=2;
    assert (cecu==ecu && cecd == ecd && cemu == emu && cemd==emd);
#endif

#ifdef TURNONCHECKS
  //gtreeprint(ci,0);
  checkgenealogy(ci,0,0);
#endif //TURNONCHECKS
    for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        restorescalefactors (ci, li);
      /* have to reset the dlikeA values in the genealogies for stepwise model */
      if (L[li].model == STEPWISE)
        for (ai = 0; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      if (L[li].model == JOINT_IS_SW)
        for (ai = 1; ai < L[li].nlinked; ai++)
          likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      // assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);   


#ifdef TURNONCHECKS
    checkprobs(ci,li);
#endif //TURNONCHECKS
    }
#ifdef TURNONCHECKS
    //gtreeprint(ci,0);
    checkpoptree(ci,0);
  checkgenealogy(ci,0,0);
  //gtreeprint(0,0);

  checkprobs(ci,-1);
  check_hgprob_sums(ci);
    for (li=0;li<nloci;li++)
    {
    checkgenealogy(ci,li,0);
    }
checkgenealogyweights(ci);
#endif //TURNONCHECKS
    return 0;
  }
}                               /* changet_RYhg */

