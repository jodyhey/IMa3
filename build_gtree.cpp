/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/* build starting genealogies */
/*********** LOCAL STUFF **********/

/* prototypes of local functions*/
static void addedge (int ci, int li, int newedge, int c[], int curgenes,
                     double *ptime, int nomigration);
static void addedgehg (int ci, int li, int newedge, int c[], int curgenes,
                     double *ptime);
static void addA (int ci, int li, int newedge, int *nd0, int *nd1);

/*** LOCAL FUNCTIONS ***/

/* version of addedge that builds the hidden genealogy when there is no population tree topology updating.  Requires migration 
 pick a coalescent time. If genes are in different populations
 add a migration event at half the time to the coalescent
*/ 
void
addedgehg (int ci, int li, int newedge, int c[], int curgenes, double *ptime)
{
  struct edge *gtree = C[ci]->G[li].gtree;
  int j;
  int p[2];
  //double nt[2]; not used
  double start_coalesce_rate;
  double rantime,migtime,upt;

  //start_coalesce_rate = (curgenes*(curgenes-1));
  start_coalesce_rate = curgenes / tperiodpriors[numsplittimes-1];
  for (j=0;j<2;j++)
  {
    if (c[j] < L[li].numgenes)
      upt = 0.0;
    else
      upt = gtree[gtree[c[j]].up[0]].time;
    p[j] = gtree[c[j]].pophg;
    gtree[c[j]].pop =  C[ci]->ancplist[p[j]][ findperiod(ci,upt) ];
    //nt[j] = TIMEMAX; not used
    assert (p[j] < numtreepops);
  }
  rantime = expo (start_coalesce_rate);
  
  if (p[0] != p[1])             // add 1 migration
  {
    j = 0;                      // first migration event on a branch (can get only 1 in this routine)
    migtime = *ptime + rantime/2.0;
    if (uniform () < 0.5)
    {
      gtree[c[0]].mighg[j].mt = migtime;
      gtree[c[0]].mighg[j + 1].mt = -1;
      gtree[c[0]].cmmhg = j+1;
      gtree[c[0]].mighg[j].mp = p[1];
      p[0] = p[1];
    }
    else
    {
      gtree[c[1]].mighg[j].mt = migtime;
      gtree[c[1]].mighg[j + 1].mt = -1;
      gtree[c[1]].cmmhg = j+1;
      gtree[c[1]].mighg[j].mp = p[0];
      p[1] = p[0];
    }
  }
  *ptime += rantime;
  assert (p[0] == p[1]);
  for (j=0;j<2;j++)
  {
    gtree[c[j]].fpophg = p[0];
    gtree[c[j]].time = *ptime;
    gtree[c[j]].fpop =  C[ci]->ancplist[gtree[c[j]].fpophg][ findperiod(ci,gtree[c[j]].time) ];
    gtree[newedge].up[j] = c[j];
    gtree[c[j]].down = newedge;
  }
  assert(gtree[c[0]].fpop == gtree[c[1]].fpop);
  gtree[newedge].pophg = p[0];
  gtree[newedge].pop = gtree[c[0]].fpop;
  assert (gtree[newedge].pophg >= 0 && gtree[newedge].pophg < npops);
  return;
}                               /* addedgehg */

void
addedge (int ci, int li, int newedge, int c[], int curgenes, double *ptime,
         int nomigration)
{
  struct edge *gtree = C[ci]->G[li].gtree;
  int i, j;
  int p[2];
  double nt[2];
  double start_coalesce_rate;

  // set small value for coalescent rate relative to splitting times when building initial genealogy
  if (npops == 1)
  {
    start_coalesce_rate = curgenes;
    p[0] = 0;
    p[1] = 0;
    nt[0] = TIMEMAX;
    nt[1] = TIMEMAX;
  }
  else
  {
    start_coalesce_rate = curgenes / tperiodpriors[numsplittimes-1];
    for (i = 0; i < 2; i++)  // get current ancestral populations and down times of those populations
    {
      j = gtree[c[i]].pop;
      while (*ptime > C[ci]->poptree[j].time)
        j = C[ci]->poptree[j].down;
      assert (j < numtreepops);
      assert (j >= 0);
      p[i] = j;
      nt[i] = C[ci]->poptree[p[i]].time;
    }
  }
  assert (p[0] < numtreepops && p[1] < numtreepops);
  if (p[0] != p[1])             // different populations at *ptime,  loop thru until picked time is less than next split or both in same pop
  {
    *ptime += expo (start_coalesce_rate);
    while (*ptime > nt[0] || *ptime > nt[1] || nomigration)

    {
      assert (p[0] < numtreepops && p[1] < numtreepops);
      if (C[ci]->poptree[p[0]].down == C[ci]->poptree[p[1]].down)
      {
        *ptime = nt[0];
        p[0] = p[1] = C[ci]->poptree[p[0]].down;
      }
      else
      {
        if (nt[0] < nt[1])
        {
          *ptime = nt[0];
          p[0] = C[ci]->poptree[p[0]].down;
          nt[0] = C[ci]->poptree[p[0]].time;
        }
        else
        {
          *ptime = nt[1];
          p[1] = C[ci]->poptree[p[1]].down;
          nt[1] = C[ci]->poptree[p[1]].time;
        }
      }
      if (p[0] == p[1])
        break;
      else
        *ptime += expo (start_coalesce_rate);
    }
  }
  if (p[0] != p[1])             // add 1 migration
  {
    j = 0;                      // first migration event on a branch (can get only 1 in this routine)
    if (uniform () < 0.5)
    {
      gtree[c[0]].mig[j].mt = *ptime;
      gtree[c[0]].mig[j + 1].mt = -1;
      gtree[c[0]].cmm = j+1;

      assert (gtree[c[0]].mig[j].mt <TIMEMAX);
      gtree[c[0]].mig[j].mp = p[1];
      p[0] = p[1];
    }
    else
    {
      gtree[c[1]].mig[j].mt = *ptime;
      gtree[c[1]].mig[j + 1].mt = -1;
      gtree[c[1]].mig[j].mp = p[0];
      gtree[c[1]].cmm = j+1;

      p[1] = p[0];
    }
  }
  assert (p[0] == p[1]);
  *ptime += expo (start_coalesce_rate);
  while (*ptime > C[ci]->poptree[p[0]].time)
    p[0] = C[ci]->poptree[p[0]].down;
  gtree[newedge].pop = p[0];
  gtree[c[0]].fpop = gtree[c[1]].fpop = p[0]; 
  assert (gtree[newedge].pop >= 0 && gtree[newedge].pop < numtreepops);
  gtree[c[0]].time = gtree[c[1]].time = *ptime;
  gtree[newedge].up[0] = c[0];
  gtree[newedge].up[1] = c[1];
  gtree[c[0]].down = gtree[c[1]].down = newedge;
  return;
}                               /* addedge */

/* under JOINT_IS_SW model, alleles states must be added for the SW portion of the genealogy at the beginning */
/* pick a number over the span of the interval between the descendant allele values (+1 on each side), following the bessel function distribution */
void
addA (int ci, int li, int newedge, int *nd0, int *nd1)
{
  int j, ai, aj, up[2], newA;
  struct edge *gtree = C[ci]->G[li].gtree;
  int d0, d1;
  double u, t0, t1, r;
  static double *p;
  if (L[li].model == STEPWISE)
    ai = 0;
  else
    ai = 1;
  for (; ai < L[li].nlinked; ai++)
  {
    u = C[ci]->G[li].uvals[ai];
    up[0] = gtree[newedge].up[0];
    up[1] = gtree[newedge].up[1];
    p = static_cast<double *> (malloc ((abs (gtree[up[0]].A[ai] - 
                gtree[up[1]].A[ai]) + 3) * sizeof (double)));
    if (up[0] < L[li].numgenes)
      t0 = gtree[up[0]].time;
    else
      t0 = gtree[up[0]].time - gtree[gtree[up[0]].up[0]].time;
    assert (t0 > 0);
    if (up[1] < L[li].numgenes)
      t1 = gtree[up[1]].time;
    else
      t1 = gtree[up[1]].time - gtree[gtree[up[1]].up[0]].time;
    assert (t1 > 0);
    for (j = 0, aj = IMIN (gtree[up[0]].A[ai], gtree[up[1]].A[ai]) - 1;
         aj <= 1 + IMAX (gtree[up[0]].A[ai], gtree[up[1]].A[ai]); aj++, j++)
    {
      d0 = abs (aj - gtree[up[0]].A[ai]);
      d1 = abs (aj - gtree[up[1]].A[ai]);
      p[j] = exp (-t0 * u + log (bessi (d0, t0 * u)) + -t1 * u + log (bessi (d1, t1 * u)));
      if (j > 0)
        p[j] += p[j - 1];
    }
    r = uniform () * p[j - 1];
    j = 0;
    while (r > p[j])
      j++;
    newA = IMIN (gtree[up[0]].A[ai], gtree[up[1]].A[ai]) - 1 + j;
    newA = IMAX (L[li].minA[ai], newA);
    newA = IMIN (L[li].maxA[ai], newA);
    *nd0 = abs (newA - gtree[up[0]].A[ai]);
    *nd1 = abs (newA - gtree[up[1]].A[ai]);
    XFREE (p);
    gtree[newedge].A[ai] = newA;
    assert (newA >= L[li].minA[ai]);
  }
}                               /* addA */

/*** GLOBAL FUNCTIONS ***/

/* this is not working for hidden genealogies */
void
makeHKY (int ci, int li, int nosimmigration)
{
  int c[2], newedge, lastnewedge, curgenes, i, j, k, ai, aj, *edgevec;
  double **distmat;
  int *ids;
  double ptime, distmax;
  struct edge *gtree = C[ci]->G[li].gtree;
  int count;
  edgevec = static_cast<int *> 
        (malloc ((2 * L[li].numgenes - 1) * sizeof (int)));
  for (i = 0; i < L[li].numgenes; i++)
    edgevec[i] = 1;
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
    edgevec[i] = -1;

  /* SAGNCHUL: Tue Sep 30 10:41:45 EDT 2008
   * NOTE: The initial population assignment is done 
   * by function IMA_genealogy_assignpopulation, which
   * is called before the call of make function such as
   * makeSW and makeIS. We remove the codes of the initial
   * population assignment part that used to be here. */

  C[ci]->G[li].mignum = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    gtree[i].up[0] = gtree[i].up[1] = gtree[i].down = -1;
  }
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
  {
    gtree[i].hkyi.scalefactor = static_cast<double *> 
            (malloc (L[li].numsites * sizeof (double)));
    gtree[i].hkyi.oldscalefactor = static_cast<double *> 
            (malloc (L[li].numsites * sizeof (double)));

    for (k = 0; k < L[li].numsites; k++)
    {
      gtree[i].hkyi.scalefactor[k] = 0.0;
      gtree[i].hkyi.oldscalefactor[k] = 0.0;
    }
    /* cr 110907.1 call new allocAndInit function so 2D double is initialized */
    gtree[i].hkyi.frac = allocAndInit2Ddouble (L[li].numsites,4);
    gtree[i].hkyi.newfrac = allocAndInit2Ddouble (L[li].numsites,4);
  } 
  

  C[ci]->G[li].mignum = 0;
  ptime = 0.0;
  curgenes = L[li].numgenes;
  distmat = orig2d_alloc2Ddouble (L[li].numgenes, L[li].numgenes);
  ids = static_cast<int *> (malloc ((L[li].numgenes) * (sizeof (int))));
  for (i = 0; i < L[li].numgenes; i++)
    ids[i] = i;
  newedge = L[li].numgenes;
  for (i = 0; i < L[li].numgenes; i++)
  {
    distmat[i][i] = 0;
    for (j = i + 1; j < L[li].numgenes; j++)
    {
      distmat[i][j] = 0;
      for (k = 0; k < L[li].numsites; k++)
      {
        if (!(L[li].model == INFINITESITES && L[li].badsite[k] == 0))
        {
          distmat[i][j] += (L[li].seq[ids[j]][k] != L[li].seq[ids[i]][k]);
        }
      }
      distmat[j][i] = distmat[i][j];
    }
  }
  ai = aj = -1;

  //should be upgma
  do
  {
    if (ai > -1)
    {
      for (i = 0; i < curgenes - 1; i++)
      {
        if (i == ai)
          for (j = i + 1; j < curgenes; j++)
            distmat[j][i] = distmat[i][j] =
              (distmat[ai][j] + distmat[aj][j]) / 2.0;
      }
      for (i = aj; i < curgenes - 1; i++)
      {
        for (j = i + 1; j < curgenes - 1; j++)
          distmat[j][i] = distmat[i][j] = distmat[i][j + 1];
      }
      ids[ai] = lastnewedge;
      for (i = aj; i < curgenes; i++)
        ids[i] = ids[i + 1];
    }
    distmax = -1e20;            // large value;         
    for (i = 0; i < curgenes - 1; i++)
      for (j = i + 1; j < curgenes; j++)
        if (distmax < distmat[i][j])

        {
          distmax = distmat[i][j];
        }
    for (count = 0, i = 0; i < curgenes - 1; i++)
      for (j = i + 1; j < curgenes; j++)
        count += distmax == distmat[i][j];
    count = randposint (count);
    for (k = -1, i = 0; i < curgenes - 1; i++)
      for (j = i + 1; j < curgenes; j++)
      {
        k += distmax == distmat[i][j];
        if (k == count)

        {
          ai = i;
          aj = j;
          assert (distmat[ai][aj] == distmax);
          k++;
        }
      }

    /* give new ids to positions ai and aj */
    //replace columns and rows associated with coal[0] with the average of values involving coal[0] and coal[1];
    //like a upgma contraction of matrix
    c[0] = ids[ai];
    c[1] = ids[aj];
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
      addedgehg (ci, li, newedge, c, curgenes, &ptime);
    else
      addedge (ci, li, newedge, c, curgenes, &ptime,nosimmigration);
    lastnewedge = newedge;
    newedge++;
    curgenes--;
  } while (newedge < 2 * L[li].numgenes - 1);

  C[ci]->G[li].roottime = ptime;
  C[ci]->G[li].root = L[li].numlines - 1;
  gtree[C[ci]->G[li].root].mig[0].mt = -1;
  gtree[C[ci]->G[li].root].cmm = 0;
  gtree[C[ci]->G[li].root].time = TIMEMAX;
  if (hiddenoptions[HIDDENGENEALOGY] == 1)
  {
    gtree[C[ci]->G[li].root].fpop = gtree[C[ci]->G[li].root].pop;
    gtree[C[ci]->G[li].root].fpophg = gtree[C[ci]->G[li].root].pophg;
  }
  XFREE (edgevec);
  orig2d_free2D ((void **) distmat, L[li].numgenes);
  XFREE (ids);
  C[ci]->G[li].hilike = -1e20;
}                               /* makeHKY */

void
makeIS (int ci, int li, int nosimmigration)
{
  int c[2], newedge, cursites, curgenes, num, coalnum, i, j, k, site, *curid,
    **distmat, coal[2], *singletons, **tempseq;
  double ptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  distmat = alloc2Dint (L[li].numgenes, L[li].numgenes);
  tempseq = alloc2Dint (L[li].numgenes, L[li].numsites);
  curid = static_cast<int *> (malloc ((L[li].numgenes) * (sizeof (int))));
  singletons = static_cast<int *> (malloc ((L[li].numsites) * (sizeof (int))));
  curgenes = L[li].numgenes;
  cursites = L[li].numsites;
  for (i = 0; i < L[li].numgenes; i++)
    curid[i] = i;

  /* SAGNCHUL: Tue Sep 30 10:41:45 EDT 2008
   * NOTE: The initial population assignment is done 
   * by function IMA_genealogy_assignpopulation, which
   * is called before the call of make function such as
   * makeSW and makeIS. We remove the codes of the initial
   * population assignment part that used to be here. */

  C[ci]->G[li].mignum = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    gtree[i].up[0] = gtree[i].up[1] = gtree[i].down = -1;
  }
  for (i = 0; i < L[li].numgenes; i++)
  {
    for (j = 0; j < L[li].numsites; j++)
      tempseq[i][j] = L[li].seq[i][j];
  }
  ptime = 0.0;
  while (curgenes > 1)
  {
    num = 0;
    for (i = 0; i < curgenes; i++)
    {
      for (j = i + 1; j < curgenes; j++)
      {
        distmat[i][j] = distmat[j][i] = 0;
        for (site = 0; site < cursites; site++)
        {
          if (tempseq[i][site] != tempseq[j][site])
          {
            distmat[i][j]++;
            distmat[j][i]++;
          }
        }
        if (distmat[i][j] == 0)
          num++;
      }
    }
    if (num == 0)
    {
      for (i = 0; i < cursites; i++)
      {
        num = 0;
        for (j = 0; j < curgenes; j++)
          num = num + tempseq[j][i];
        if (num == 1 || num == (curgenes - 1))
          singletons[i] = 1;

        else
          singletons[i] = 0;
      }
      num = 0;
      for (i = 0; i < curgenes; i++)
      {
        for (j = i + 1; j < curgenes; j++)
        {
          distmat[i][j] = distmat[j][i] = 0;
          for (site = 0; site < cursites; site++)
          {
            if (tempseq[i][site] != tempseq[j][site] && singletons[site] == 0)
            {
              distmat[i][j]++;
              distmat[j][i]++;
            }
          }
          if (distmat[i][j] == 0)
            num++;
        }
      }
    }
    if (num == 0)
    {
      printf ("Data not compatible with infinite sites model\n");
      for (i = 0; i < curgenes; i++)
      {
        printf ("seq %i: ", curid[i]);
        for (j = 0; j < cursites; j++)
          printf ("%i", tempseq[i][j]);
        printf ("\n");
      }
      IM_err(IMERR_INFINITESITESFAIL,"locus #: %d, name: %s,  sequence#: %d", li,L[li].name, curid[i]);
      
    }
    coalnum = (int) (uniform () * num);
    num = 0;
    for (i = 0; i < curgenes; i++)
    {
      for (j = i + 1; j < curgenes; j++)
      {
        if (distmat[i][j] == 0)
        {
          num++;
          if (num > coalnum)
          {
            if (i < j)
            {
              coal[0] = i;
              coal[1] = j;
            }
            else
            {
              coal[0] = j;
              coal[1] = i;
            }
            i = curgenes;
            j = curgenes;
          }
        }
      }
    }
    for (i = 0; i < cursites; i++)
    {
      if (tempseq[coal[0]][i] != tempseq[coal[1]][i])
      {
        for (j = 0; j < curgenes; j++)
        {
          for (k = i; k < cursites - 1; k++)
            tempseq[j][k] = tempseq[j][k + 1];
        }
        i--;
        cursites--;
      }
    }
    newedge = 2 * L[li].numgenes - curgenes;
    c[0] = curid[coal[0]];
    c[1] = curid[coal[1]];
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
      addedgehg (ci, li, newedge, c, curgenes, &ptime);
    else
      addedge (ci, li, newedge, c, curgenes, &ptime,nosimmigration);
   
    for (i = coal[1]; i < curgenes - 1; i++)
    {
      for (j = 0; j < cursites; j++)
        tempseq[i][j] = tempseq[i + 1][j];
      curid[i] = curid[i + 1];
    }
    curid[coal[0]] = newedge;
    curgenes--;
  }

  C[ci]->G[li].roottime = ptime;
  C[ci]->G[li].root = L[li].numlines - 1;
  gtree[C[ci]->G[li].root].fpop = gtree[C[ci]->G[li].root].pop;
  if (hiddenoptions[HIDDENGENEALOGY] == 1)
  {
    gtree[C[ci]->G[li].root].fpophg = gtree[C[ci]->G[li].root].pophg;
  }
  gtree[C[ci]->G[li].root].mig[0].mt = -1;
  gtree[C[ci]->G[li].root].cmm = 0;
  gtree[C[ci]->G[li].root].time = TIMEMAX;
  orig2d_free2D ((void **) tempseq, L[li].numgenes);
  orig2d_free2D ((void **) distmat, L[li].numgenes);
  XFREE (curid);
  XFREE (singletons);
  C[ci]->G[li].hilike = -1e20;
}                               /* makeIS */

void
makeJOINT_IS_SW (int ci, int li, int nosimmigration)
{
  int c[2], newedge, cursites, curgenes, num, i, j, k, site, *curid,
    **distmat, *singletons, **tempseq;
  double ptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  double **distmatA;            // for microsats if locus is JOINT_IS_SW
  double distAmax;
  int ai, aj, minA, maxA;
  int found;                    //debugging
  int nd0, nd1;
  distmat = alloc2Dint (L[li].numgenes, L[li].numgenes);
  distmatA = orig2d_alloc2Ddouble (L[li].numgenes, L[li].numgenes);

  for (ai = 1; ai < L[li].nlinked; ai++)
  {
    minA = 1000;
    maxA = 0;
    for (i = 0; i < L[li].numgenes; i++)
    {
      if (minA > gtree[i].A[ai])
        minA = gtree[i].A[ai];
      if (maxA < gtree[i].A[ai])
        maxA = gtree[i].A[ai];
    }
    L[li].minA[ai] = IMAX (1, ((maxA + minA) / 2) - (maxA - minA));
    L[li].maxA[ai] = IMIN (1000, ((maxA + minA) / 2) + (maxA - minA));
    for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
      gtree[i].A[ai] = -1;
  }
  tempseq = alloc2Dint (L[li].numgenes, L[li].numsites);
  curid = static_cast<int *> (malloc ((L[li].numgenes) * (sizeof (int))));
  singletons = static_cast<int *> (malloc ((L[li].numsites) * (sizeof (int))));
  curgenes = L[li].numgenes;
  cursites = L[li].numsites;
  for (i = 0; i < L[li].numgenes; i++)
    curid[i] = i;

  /* SAGNCHUL: Tue Sep 30 10:41:45 EDT 2008
   * NOTE: The initial population assignment is done 
   * by function IMA_genealogy_assignpopulation, which
   * is called before the call of make function such as
   * makeSW and makeIS. We remove the codes of the initial
   * population assignment part that used to be here. */

  C[ci]->G[li].mignum = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    gtree[i].up[0] = gtree[i].up[1] = gtree[i].down = -1;
  }
  for (i = 0; i < L[li].numgenes; i++)
  {
    for (j = 0; j < L[li].numsites; j++)
      tempseq[i][j] = L[li].seq[i][j];
  }
  ptime = 0.0;
  while (curgenes > 1)
  {
    num = 0;
    for (i = 0; i < curgenes; i++)
    {
      for (j = i + 1; j < curgenes; j++)
      {
        distmat[i][j] = distmat[j][i] = 0;
        for (site = 0; site < cursites; site++)
        {
          if (tempseq[i][site] != tempseq[j][site])
          {
            distmat[i][j]++;
            distmat[j][i]++;
          }
        }
        if (distmat[i][j] == 0)
          num++;
      }
    }
    if (num == 0)
    {
      for (i = 0; i < cursites; i++)
      {
        num = 0;
        for (j = 0; j < curgenes; j++)
          num = num + tempseq[j][i];
        if (num == 1 || num == (curgenes - 1))
          singletons[i] = 1;
        else
          singletons[i] = 0;
      }
      num = 0;
      for (i = 0; i < curgenes; i++)
      {
        for (j = i + 1; j < curgenes; j++)
        {
          distmat[i][j] = distmat[j][i] = 0;
          for (site = 0; site < cursites; site++)
          {
            if (tempseq[i][site] != tempseq[j][site] && singletons[site] == 0)
            {
              distmat[i][j]++;
              distmat[j][i]++;
            }
          }
          if (distmat[i][j] == 0)
            num++;
        }
      }
    }
    if (num == 0)
    {
      printf ("Data not compatible with infinite sites model\n");
      for (i = 0; i < curgenes; i++)
      {
        printf ("seq %i: ", curid[i]);
        for (j = 0; j < cursites; j++)
          printf ("%i", tempseq[i][j]);
        printf ("\n");
      }
      IM_err(IMERR_INFINITESITESFAIL,"locus %d sequence# %d", li, curid[i]);
    }
    for (i = 0; i < L[li].numgenes - 1; i++)
    {
      distmatA[i][i] = 0;
      for (j = i + 1; j < L[li].numgenes; j++)
        distmatA[j][i] = distmatA[i][j] = 0;
    }
    for (i = 0; i < curgenes - 1; i++)
    {
      distmatA[i][i] = 0;
      for (j = i + 1; j < curgenes; j++)
      {
        distmatA[i][j] = distmatA[j][i] = 0;
        for (k = 1; k < L[li].nlinked; k++)
          distmatA[i][j] +=
            SQR (abs (gtree[curid[i]].A[k] - gtree[curid[j]].A[k]));
        distmatA[j][i] = distmatA[i][j];
      }
    }
    distAmax = 1e20;            // large value;         
    num = 0;
    found = 0;
    for (i = 0; i < curgenes - 1; i++)
    {
      for (j = i + 1; j < curgenes; j++)
      {
        if (distmat[i][j] == 0)
        {
          if (distAmax > distmatA[i][j])
          {
            distAmax = distmatA[i][j];
            ai = i;
            aj = j;
            found = 1;
          }
        }
      }
    }
    for (i = 0; i < cursites; i++)
    {
      if (tempseq[ai][i] != tempseq[aj][i])
      {
        for (j = 0; j < curgenes; j++)
        {
          for (k = i; k < cursites - 1; k++)
            tempseq[j][k] = tempseq[j][k + 1];
        }
        i--;
        cursites--;
      }
    }
    c[0] = curid[ai];
    c[1] = curid[aj];
    assert (found);
    newedge = 2 * L[li].numgenes - curgenes;
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
      addedgehg (ci, li, newedge, c, curgenes, &ptime);
    else
      addedge (ci, li, newedge, c, curgenes, &ptime,nosimmigration);
    addA (ci, li, newedge, &nd0, &nd1);
    for (i = aj; i < curgenes - 1; i++)
    {
      for (j = 0; j < cursites; j++)
        tempseq[i][j] = tempseq[i + 1][j];
      curid[i] = curid[i + 1];
    }
    curid[ai] = newedge;
    curgenes--;
  }
  C[ci]->G[li].roottime = ptime;
  C[ci]->G[li].root = L[li].numlines - 1;
  gtree[C[ci]->G[li].root].fpop = gtree[C[ci]->G[li].root].pop;
  if (hiddenoptions[HIDDENGENEALOGY] == 1)
  {
    gtree[C[ci]->G[li].root].fpophg = gtree[C[ci]->G[li].root].pophg;
  }
  gtree[C[ci]->G[li].root].mig[0].mt = -1;
  gtree[C[ci]->G[li].root].cmm = 0;
  gtree[C[ci]->G[li].root].time = TIMEMAX;
  orig2d_free2D ((void **) tempseq, L[li].numgenes);
  orig2d_free2D ((void **) distmat, L[li].numgenes);
  orig2d_free2D ((void **) distmatA, L[li].numgenes);
  XFREE (curid);
  XFREE (singletons);
  C[ci]->G[li].hilike = -1e20;
}                               /* makeJOINT_IS_SW */

/*Not really sure how best to do the SW genealogies.  here is the rough procedure 
- build a distance matrix
- pick pair with lowest distance
- make a node
- pack an ancestral allele, based on distribution of alleles sizes, given the mutation rate and the branck lengths
- build a new distance matrix and repeat */
void
makeSW (int ci, int li, int nosimmigration)
{
  int i, j, k, ai, aj;
  int curgenes, newedge, c[2];
  double ptime;
  struct edge *gtree = C[ci]->G[li].gtree;
  double **distmatA;
  int *ids;
  double distAmax;
  int minA, maxA;
  int nd0, nd1;

  for (ai = 0; ai < L[li].nlinked; ai++)
  {
    minA = 1000;
    maxA = 0;
    for (i = 0; i < L[li].numgenes; i++)
    {
      if (minA > gtree[i].A[ai])
        minA = gtree[i].A[ai];
      if (maxA < gtree[i].A[ai])
        maxA = gtree[i].A[ai];
    }
    L[li].minA[ai] = IMAX (1, ((maxA + minA) / 2) - (maxA - minA));
    L[li].maxA[ai] = IMIN (1000, ((maxA + minA) / 2) + (maxA - minA));
    for (i = L[li].numgenes; i < L[li].numlines; i++)
      gtree[i].A[ai] = -1;
  }

  distmatA = orig2d_alloc2Ddouble (L[li].numgenes, L[li].numgenes);
  ids = static_cast<int *> (malloc ((L[li].numgenes) * (sizeof (int))));
  for (i = 0; i < L[li].numgenes; i++)
  {
    ids[i] = i;

    //      distmatA[i] = malloc((L[li].numgenes)*(sizeof(double)));
  }
  ptime = 0;
  curgenes = L[li].numgenes;

  /* SAGNCHUL: Tue Sep 30 10:41:45 EDT 2008
   * NOTE: The initial population assignment is done 
   * by function IMA_genealogy_assignpopulation, which
   * is called before the call of make function such as
   * makeSW and makeIS. We remove the codes of the initial
   * population assignment part that used to be here. */

  for (i = 0; i < L[li].numlines; i++)
  {
    gtree[i].up[0] = gtree[i].up[1] = gtree[i].down = -1;
    gtree[i].mig[0].mt = -1;
    gtree[i].cmm = 0;
  }
  C[ci]->G[li].mignum = 0;
  newedge = L[li].numgenes;

  do
  {
    for (i = 0; i < L[li].numgenes - 1; i++)
    {
      distmatA[i][i] = 0;
      for (j = i + 1; j < L[li].numgenes; j++)
        distmatA[j][i] = distmatA[i][j] = 0;
    }
    for (i = 0; i < curgenes - 1; i++)
      for (j = 0; j < curgenes; j++)
      {
        distmatA[i][j] = distmatA[j][i] = 0;
        for (k = 0; k < L[li].nlinked; k++)
          distmatA[i][j] +=
            SQR (abs (gtree[ids[i]].A[k] - gtree[ids[j]].A[k]));
        distmatA[j][i] = distmatA[i][j];
      }
    distAmax = 1e20;            // large value;         
    for (i = 0; i < curgenes; i++)
      for (j = i + 1; j < curgenes; j++)
      {
        if (distAmax > distmatA[i][j])

        {
          distAmax = distmatA[i][j];
          ai = i;
          aj = j;
        }
      }

    /* give new ids to positions ai and aj */

    //replace columns and rows associated with coal[0] with the average of values involving coal[0] and coal[1];
    //like a upgma contraction of matrix
    c[0] = ids[ai];
    c[1] = ids[aj];
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
      addedgehg (ci, li, newedge, c, curgenes, &ptime);
    else
      addedge (ci, li, newedge, c, curgenes, &ptime,nosimmigration);
    addA (ci, li, newedge, &nd0, &nd1);
    ids[ai] = newedge;
    for (i = aj; i < curgenes - 1; i++)
      ids[i] = ids[i + 1];
    newedge++;
    curgenes--;
  } while (newedge < 2 * L[li].numgenes - 1);

  C[ci]->G[li].roottime = ptime;
  C[ci]->G[li].root = L[li].numlines - 1;
  gtree[C[ci]->G[li].root].fpop = gtree[C[ci]->G[li].root].pop;
  if (hiddenoptions[HIDDENGENEALOGY] == 1)
  {
    gtree[C[ci]->G[li].root].fpophg = gtree[C[ci]->G[li].root].pophg;
  }
  gtree[C[ci]->G[li].root].mig[0].mt = -1;
  gtree[C[ci]->G[li].root].cmm = 0;
  gtree[C[ci]->G[li].root].time = TIMEMAX;
  orig2d_free2D ((void **) distmatA, L[li].numgenes);
  XFREE (ids);
  C[ci]->G[li].hilike = -1e20;
}                               /* makeSW */
