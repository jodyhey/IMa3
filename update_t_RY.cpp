/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS

#include "ima.hpp"
#include "update_gtree_common.hpp"

/*********** LOCAL STUFF **********/

static struct genealogy_weights holdallgweight_t_RY;
static struct genealogy_weights holdgweight_t_RY[MAXLOCI];
static struct probcalc holdallpcalc_t_RY;
static int largestsamp;

static double beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime);
static double aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime);
static void storegenealogystats_all_loci (int ci, int mode);

/********* LOCAL FUNCTIONS **************/

void
storegenealogystats_all_loci (int ci, int mode)
{
  static double holdlength[MAXLOCI], holdtlength[MAXLOCI];
  static double holdroottime[MAXLOCI];
  static int holdroot[MAXLOCI];
  static int holdmig[MAXLOCI];
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
    }
  }
  return;
}                               // storegenealogystats  

double
aftersplit (int tnode, double oldt, double newt, double tau_d, double ptime)
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
beforesplit (int tnode, double oldt, double newt, double tau_u, double ptime)
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

/*************GLOBAL FUNCTIONS ******************/

void
init_t_RY (void)
{
  int li, j;
  init_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
    init_genealogy_weights (&holdgweight_t_RY[li]);
  init_probcalc (&holdallpcalc_t_RY);
  for (largestsamp = 0, j = 0; j < nloci; j++)
    if (largestsamp < L[j].numlines)
      largestsamp = L[j].numlines;
}                               // init_changet_RY

void
free_t_RY (void)
{
  int li;
  free_genealogy_weights (&holdallgweight_t_RY);
  for (li = 0; li < nloci; li++)
  {
    free_genealogy_weights (&holdgweight_t_RY[li]);
  }
  free_probcalc (&holdallpcalc_t_RY);

}                               // free_changet_RY

/* 
Notes on changet_RY()  implements updating of Rannala and Yang (2003)   

This application is pretty much the same as theirs - changing times on a species tree that contains a gene tree. 
The big difference is that IM includes migration.   This means that we have to count migration events and change 
migration times in the same was as we change coalescent times and count how many get changed. 

R&Y also only change coalescent times in populations affected by a changed t.  But because we have migration
there is more entanglement between populations.  It seems best to change all times within an interval that is 
affected by a changing t. 

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

if t > tau  see aftersplit()  t is older than the current split time 
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
changet_RY1 (int ci, int timeperiod)    // after Rannala and Yang (2003)  - rubberband method
{
  double metropolishastingsratio, newt, oldt;
  double likelihoodratio,priorratio,proposalratio;
  double pdgnew[MAXLOCI + MAXLINKED], pdgnewsum, pdgoldsum, probgnewsum,
    temppdg;
  double t_u_hterm, t_d_hterm;
  int li, i, j, ecd, ecu, emd, emu, ai, ui;
  struct genealogy *G;
  struct edge *gtree;
  double t_d, t_u, t_u_prior, t_d_prior;

  checkupdatescalarer(&C[ci]->RYwidthinfo[timeperiod]);
  t_u = (timeperiod == 0) ? 0 : C[ci]->tvals[timeperiod - 1];
  t_d =
    (timeperiod ==
     (lastperiodnumber - 1)) ? TIMEMAX : C[ci]->tvals[timeperiod + 1];
  t_d_prior = DMIN (T[timeperiod].pr.max, t_d);
  t_u_prior = DMAX (T[timeperiod].pr.min, t_u);
  oldt = C[ci]->tvals[timeperiod];
  //newt = getnewt (t_u_prior, t_d_prior, oldt, 1);
  newt = getnewt (t_u_prior, t_d_prior, oldt,C[ci]->RYwidthinfo[timeperiod].updatescalarval);
  
  t_u_hterm = (newt - t_u) / (oldt - t_u);
  if (timeperiod == lastperiodnumber - 1)
  {
    t_d_hterm = 1;
  }
  else
  {
    t_d_hterm = (t_d - newt) / (t_d - oldt);
  }

  copy_treeinfo (&holdallgweight_t_RY, &C[ci]->allgweight);  // try turning this off and forcing all recalculations
  copy_probcalc (&holdallpcalc_t_RY, &C[ci]->allpcalc);

  pdgoldsum = C[ci]->allpcalc.pdg;
  setzero_genealogy_weights (&C[ci]->allgweight);
  ecd = ecu = emd = emu = 0;
  pdgnewsum = 0;
  probgnewsum = 0;
  storegenealogystats_all_loci (ci, 0);
  C[ci]->tvals[timeperiod] = newt;
  for (i = 0; i < nurates; i++)
    pdgnew[i] = 0;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    copy_treeinfo (&holdgweight_t_RY[li], &G->gweight);
    for (i = 0; i < L[li].numlines; i++)
    {
      if (gtree[i].down != -1)
      {
        if (gtree[i].time <= oldt && gtree[i].time > t_u)

        {
          gtree[i].time =
            beforesplit (timeperiod, oldt, newt, t_u, gtree[i].time);
          assert (gtree[i].time != newt);
          ecu++;
        }
        else
        {
          if (gtree[i].time > oldt && gtree[i].time < t_d)
          {
            gtree[i].time =
              aftersplit (timeperiod, oldt, newt, t_d, gtree[i].time);
            assert (gtree[i].time != newt);
            ecd++;
          }
          //else  do not change the time
        }
        j = 0;
        while (gtree[i].mig[j].mt > -0.5)
        {
          assert (gtree[i].mig[j].mt <TIMEMAX);
          if (gtree[i].mig[j].mt <= oldt && gtree[i].mig[j].mt > t_u)
          {
            gtree[i].mig[j].mt =
              beforesplit (timeperiod, oldt, newt, t_u,
                           gtree[i].mig[j].mt);
            emu++;
          }
          else
          {
            assert (oldt <TIMEMAX);
            if (gtree[i].mig[j].mt > oldt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, oldt, newt, t_d, 
                            gtree[i].mig[j].mt);
              emd++;
            }
            // else no need to change the time
          }
          j++;
        }
      }
    }
    if (G->roottime <= oldt && G->roottime > t_u)
      G->roottime =
        beforesplit (timeperiod, oldt, newt, t_u, G->roottime);
    else if (G->roottime > oldt && G->roottime < t_d)
      G->roottime =
        aftersplit (timeperiod, oldt, newt, t_d, G->roottime);
    setzero_genealogy_weights (&G->gweight);
    int treeweightcallcode = 5;  // a debugging code,  if treeweight has an error results are written to an output file with this code 
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

  assert (!ODD (ecd));
  assert (!ODD (ecu));
  ecd /= 2;
  ecu /= 2;
  integrate_tree_prob (ci, &C[ci]->allgweight, &holdallgweight_t_RY,  &C[ci]->allpcalc, &holdallpcalc_t_RY);   // try enforcing full cacullation
//initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  likelihoodratio = (pdgnewsum - pdgoldsum);
  priorratio = (C[ci]->allpcalc.probg - holdallpcalc_t_RY.probg);
  proposalratio = (ecd + emd) * log (t_d_hterm) + (ecu + emu) * log (t_u_hterm);
/* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  if  (calcoptions[CALCMARGINALLIKELIHOOD])
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
    }
#ifdef TURNONCHECKS
    for (li = 0; li < nloci; li++)
      checkprobs(ci,li);
#endif //TURNONCHECKS
    C[ci]->allpcalc.pdg = pdgnewsum;
    C[ci]->poptree[C[ci]->droppops[timeperiod + 1][0]].time =
      C[ci]->poptree[C[ci]->droppops[timeperiod + 1][1]].time = newt;
    C[ci]->RYwidthinfo[timeperiod].recentaccp += 1;
    C[ci]->RYwidthinfo[timeperiod].allaccp += 1;

    return 1;
  }
  else
  {
  rejectjump:
    copy_treeinfo (&C[ci]->allgweight, &holdallgweight_t_RY);
    copy_probcalc (&C[ci]->allpcalc, &holdallpcalc_t_RY);
    assert (pdgoldsum == C[ci]->allpcalc.pdg);
    C[ci]->tvals[timeperiod] = oldt;
    for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      storegenealogystats_all_loci (ci, 1);
      copy_treeinfo (&G->gweight, &holdgweight_t_RY[li]);
      for (i = 0; i < L[li].numlines; i++)
      {
        if (gtree[i].down != -1)
        {
          if (gtree[i].time <= newt && gtree[i].time > t_u)
          {
            gtree[i].time =
              beforesplit (timeperiod, newt, oldt, t_u, gtree[i].time);
            //cecu++;
          }

          else
          {
            if (gtree[i].time > newt && gtree[i].time < t_d)
            {
              gtree[i].time =
                aftersplit (timeperiod, newt, oldt, t_d, gtree[i].time);
              //cecl++;
            }
          }
          j = 0;
          while (gtree[i].mig[j].mt > -0.5)
          {
            if (gtree[i].mig[j].mt <= newt && gtree[i].mig[j].mt > t_u)
            {
              gtree[i].mig[j].mt =
                beforesplit (timeperiod, newt, oldt, 
                             t_u, gtree[i].mig[j].mt);
              //cemu++;
            }
            else if (gtree[i].mig[j].mt > newt && gtree[i].mig[j].mt < t_d)
            {
              gtree[i].mig[j].mt =
                aftersplit (timeperiod, newt, oldt, t_d, 
                            gtree[i].mig[j].mt);
              //ceml++;
            }
            j++;
          }
        }
      }
//        assert(fabs(C[ci]->G[li].gtree[  C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time - C[ci]->G[li].roottime) < 1e-8);    
    }
    /*    assert(ecu==cecu/2);
       assert(ecd==cecl/2);
       assert(emu==cemu);
       assert(emd==ceml); */
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
    }
#ifdef TURNONCHECKS
    for (li = 0; li < nloci; li++)
      checkprobs(ci,li);
#endif //TURNONCHECKS
    return 0;
  }
}                               /* changet_RY1 */

