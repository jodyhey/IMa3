/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* for debugging - do not include in release */


/*
some useful functions to check things

the following functions create files:
  gtreeprint(ci,li) prints a genealogy in a readable form
  poptreeprint(ci) same for population tree
  void printgenealogyweights(int ci,int li)  prints genealogy weights in a readable form

  be careful because they can get VERY large if just left on

following functions check that some things match what is expected

  checkgenealogy()
  checkpoptree()
  void checkprobs(int ci, int li)  recalculate the likelihood,  and the main priors to see if results are valid numbers and if they match what is currently held
  checkgenealogyweights(int ci)  check that sums of weights equal the overall weights


*/

#undef GLOBVARS
#include "ima.hpp"
#ifdef TURNONCHECKS
extern void poptreewrite (int ci, char *buildstr);

/*********** LOCAL STUFF **********/

/* treevals only used in gtreeprint() and poptreeprint() for debugging */
struct treevals
{
  int nodenum;
  int c_or_m;
  int up1;
  int up2;
  int down;
  double utime;
  double dtime;
  double timei;
  int migcount;
  int migcounthg; //hgstuff
  int startpop;
  int finpop;
  int startpophg; //hgstuff
  int finpophg; //hgstuff
  int mut;
  int A[MAXLINKED];
  double dlikeA[MAXLINKED];
  int flag;
  int check;
} *savetree[MAXLOCI];

/* local function prototypes */
static void shellgtreevals (struct treevals *lptr, int length);

/* local functions */

void
shellgtreevals (struct treevals *lptr, int length)
{
  double aln = 1.442695022, tiny = 1.0e-5;
  struct treevals t;
  static int nn, m, lognb2, i, j, k, l;
  lognb2 = (int) floor (log (static_cast<double>(length)) * aln + tiny);
  m = length;
  for (nn = 1; nn <= lognb2; nn++)
  {
    m = m / 2;
    k = length - m;
    for (j = 0; j <= k - 1; j++)
    {
      i = j;
    reloop:l = i + m;
      if ((lptr + l)->dtime < (lptr + i)->dtime)
      {
        t = *(lptr + i);
        *(lptr + i) = *(lptr + l);
        *(lptr + l) = t;
        i = i - m;
        if (i >= 0)
          goto reloop;
      }
    }
  }
}                               /* shellfithet */

/********** GLOBAL FUNCTIONS ***********/

/* callsource not that useful,  commented out on 4/1/2008 */

/* this version prints out very long numbers - substitute gtreeprinhold for shorter numbers */
void
gtreeprint_printlongnumbers (int ci, int li/*, int step , int callsource */ )
/* use this in debugging mode to see what a genealogy looks like for a particular parameter set */
{
  int p, i, j, k, site, node, up1, up2, nm, totalm;
  int startpop, periodi;
  double upt;
  struct treevals *tvcptr;
  struct treevals *tvmptr;
  char sc[MAXPERIODS][5];
  int checkpt[MAXPERIODS] = {0}; // added initialization 8/18/2016, might have been a bug
  FILE *treeprintfile;
  char treefilename[] = "gtreeprint.out";
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int fromi, toi;
  /* hgstuff */
  int totalmhg;
  int nmhg,khg,found;
  struct treevals *tvmhgptr;
  int checkfpop;

/*	switch (callsource)
		{
		case 8 : strcpy(sourcestring,"make_JOINT_IS_SW()"); break;
		case 7 : strcpy(sourcestring,"makeSW()"); break;
		case 6 : strcpy(sourcestring,"setup_chains()"); break;
		case 5 : strcpy(sourcestring,"addmigration()"); break;
		case 3 : strcpy(sourcestring,"makeIS()"); break;
		case 4 : strcpy(sourcestring,"main()"); break;
		case 0 : strcpy(sourcestring,"updategenealogy()"); break;
		case 1 : strcpy(sourcestring,"treeweight()"); break;
		case 2 : strcpy(sourcestring,"integrate_tree_prob()"); break;
		default : strcpy(sourcestring,"UNKNOWN");
		} */
  for (i = 0; i < lastperiodnumber; i++)
  {
    sprintf (sc[i], "t%d>", i);
    checkpt[i] = 1;
  }
  totalm = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
      j++;
    totalm += j;
  }
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    totalmhg = 0;
    for (i = 0; i < 2 * L[li].numgenes - 1; i++)
    {
      j = 0;
      while (gtree[i].mighg[j].mt > -0.5)
        j++;
      totalmhg += j;
    }
  }
  tvcptr = static_cast<treevals *>
	  (malloc ((2 * L[li].numgenes - 1) * (sizeof (struct treevals))));
  tvmptr = static_cast<treevals *>
	  (malloc (totalm * sizeof (struct treevals)));
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    tvmhgptr = static_cast<treevals *>
	    (malloc (totalmhg * sizeof (struct treevals)));
  }
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
    tvcptr[i].mut = 0;
  for (site = 0; site < L[li].numsites; site++)
  {
    for (j = 0; j < L[li].numgenes; j++)
      gtree[j].mut = L[li].seq[j][site];
    for (j = L[li].numgenes; j < L[li].numlines; j++)
      gtree[j].mut = -1;
    for (j = 0; j < L[li].numgenes; j++)
      labelgtree (ci, li, j);
    for (node = L[li].numgenes; node < L[li].numlines; node++)
    {
      i = node - L[li].numgenes;
      up1 = gtree[node].up[0];
      up2 = gtree[node].up[1];
      if (gtree[node].down != -1)
        if (i != gtree[gtree[node].down].up[0]
            && i != gtree[gtree[node].down].up[1])
          tvcptr[node].flag = 2;
      if (gtree[node].up[0] != -1)
        if (i != gtree[gtree[node].up[0]].down)
          tvcptr[node].flag = 3;
      if (gtree[node].up[1] != -1)
        if (i != gtree[gtree[node].up[1]].down)
          tvcptr[node].flag = 4;
      if ((gtree[up1].mut == 0
           && gtree[up2].mut == 1)
          || (gtree[up1].mut == 1 && gtree[up2].mut == 0))
      {
        if (gtree[node].down == -1)     /* root - not clear where mutation is , put it on left */
        {
          tvcptr[up1].mut++;
        }
        else
        {
          if (gtree[gtree[node].down].mut == gtree[up1].mut)
          {
            tvcptr[up2].mut++;
          }
          else
          {
            tvcptr[up1].mut++;
          }
        }
      }
    }
  }

  for (i = 0, k = 0,khg=0; i < 2 * L[li].numgenes - 1; i++)
  {
    tvcptr[i].nodenum = i;
    tvcptr[i].c_or_m = 0;
    tvcptr[i].flag = 0;
    tvcptr[i].up1 = gtree[i].up[0];
    tvcptr[i].up2 = gtree[i].up[1];
    tvcptr[i].down = gtree[i].down;
    tvcptr[i].dtime = gtree[i].time;
    if (i >= L[li].numgenes)
      upt = gtree[gtree[i].up[0]].time;
    else
      upt = 0;
    tvcptr[i].utime = upt;
    tvcptr[i].timei = gtree[i].time - upt;
    tvcptr[i].startpop = gtree[i].pop;
    tvcptr[i].finpop =  gtree[i].fpop;
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      tvcptr[i].startpophg = gtree[i].pophg;
    j = 0;
    while (gtree[i].mig[j].mt > -1)
      j++;
    tvcptr[i].migcount = j;

    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      nmhg = 0;
      while (gtree[i].mighg[nmhg].mt > -1)
        nmhg++;
      tvcptr[i].migcounthg = nmhg;
      if (nmhg > 0)
        tvcptr[i].finpophg = gtree[i].mighg[nmhg - 1].mp;
      else
        tvcptr[i].finpophg = tvcptr[i].startpophg;
    }

    if (j > 0)
      checkfpop = gtree[i].mig[j - 1].mp;
    else
      checkfpop = tvcptr[i].startpop;
    while (tvcptr[i].dtime > C[ci]->poptree[checkfpop].time)
    {
      checkfpop = C[ci]->poptree[checkfpop].down;
    }

    if (checkfpop != tvcptr[i].finpop || tvcptr[i].finpop >= numtreepops)
      tvcptr[i].flag = 5;
    if (L[li].model == STEPWISE)
    {
      tvcptr[i].A[0] = gtree[i].A[0];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[0];
    }
    if (L[li].model == JOINT_IS_SW)
    {
      tvcptr[i].A[0] = gtree[i].A[1];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[1];
    }
    nm = j;
    if (nm > 0)
      for (startpop = gtree[i].pop, j = 0; j < nm; j++)
      {
        tvmptr[k].nodenum = i;
        tvmptr[k].c_or_m = 1;
        tvmptr[k].dtime = gtree[i].mig[j].mt;
        periodi = findperiod (ci, tvmptr[k].dtime);

        /* could change population by passing into the next period */
        while (C[ci]->poptree[startpop].e <= periodi
               && C[ci]->poptree[startpop].e != -1)
          startpop = C[ci]->poptree[startpop].down;
        tvmptr[k].startpop = startpop;
        tvmptr[k].finpop = gtree[i].mig[j].mp;
        tvmptr[k].check = 0;
        startpop = tvmptr[k].finpop;
        k++;
      }
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      if (nmhg> 0)
        for (startpop = gtree[i].pophg, j = 0; j < nmhg; j++)
        {
          tvmhgptr[khg].nodenum = i;
          tvmhgptr[khg].c_or_m = 1;
          tvmhgptr[khg].dtime = gtree[i].mighg[j].mt;
          tvmhgptr[khg].startpop = startpop;
          tvmhgptr[khg].finpop = gtree[i].mighg[j].mp;
          tvmhgptr[khg].check = 0;
          startpop = tvmhgptr[khg].finpop;
          khg++;
          assert(khg<= totalmhg);

        }
    }
  }
  shellgtreevals (tvcptr, L[li].numgenes);
  shellgtreevals (tvcptr + L[li].numgenes, L[li].numgenes - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  /*fprintf (treeprintfile,
           "chain: %d Locus: %d  Step: %d  Calling function: %s\n", ci, li,
           step, sourcestring);
           */
  fprintf (treeprintfile,  "chain: %d Locus: %d  Step: %d   Population tree: %s\n", ci, li,step,C[ci]->chainpoptreestring);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.14f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  if (hiddenoptions[HIDDENGENEALOGY]==0)
    fprintf (treeprintfile, "  current p(D|G): %.15f current p(G) (joint all loci) %.15f\n", C[ci]->G[li].pdg , C[ci]->allpcalc.probg);
  else
    fprintf (treeprintfile, "  current p(D|G): %.15f current p(G_h|G) (joint all loci) %.15f  current p(G) (joint all loci) %.15f\n", C[ci]->G[li].pdg , C[ci]->allpcalc.probhgg, C[ci]->allpcalc.probg);
  if (hiddenoptions[HIDDENGENEALOGY]==0)
  {
    fprintf (treeprintfile,
             "\tNode#\tup1\tup2\tdown\tdtime\tutime\ttimei\tpop\tfinpop\t#mig\tflag");
  }
  else
  {
    fprintf (treeprintfile,
      "\tNode#\tup1\tup2\tdown\tdtime\tutime\ttimei\tpop\tfinpop\tpophg\tfinpophg\t#mighg\t#mig\tflag");
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    fprintf (treeprintfile, "\tA\tld\n");
  else
    fprintf (treeprintfile, "\n");
  p = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    while (tvcptr[i].dtime > C[ci]->tvals[p] && checkpt[p])

    {
      if (hiddenoptions[HIDDENGENEALOGY]==0)
        fprintf (treeprintfile,
                 "%s==============================================================================================================\n",
                 sc[p]);
      else
        fprintf (treeprintfile,
                 "%s===========================================================================================================================================\n",
                 sc[p]);

      checkpt[p] = 0;
      p++;
    }
    //  assert(tvcptr[i].down == -1 || gtree[tvcptr[i].down].pop == tvcptr[i].finpop);
    if (tvcptr[i].down != -1)
      fprintf (treeprintfile,
             "\t%3d\t%3d\t%3d\t%3d\t%7.14f\t%7.14f\t%7.14f\t%3d\t%3d",
             tvcptr[i].nodenum, tvcptr[i].up1, tvcptr[i].up2,
             tvcptr[i].down, tvcptr[i].dtime, tvcptr[i].utime,
             tvcptr[i].timei, tvcptr[i].startpop, tvcptr[i].finpop);
    else
      fprintf (treeprintfile,
             "\t%3d\t%3d\t%3d\t%3d\t%7.11f\t%7.14f\t%7.11f\t%3d\t%3d",
             tvcptr[i].nodenum, tvcptr[i].up1, tvcptr[i].up2,
             tvcptr[i].down, tvcptr[i].dtime, tvcptr[i].utime,
             tvcptr[i].timei, tvcptr[i].startpop, tvcptr[i].finpop);
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      fprintf (treeprintfile,"\t%3d\t%3d\t%d",tvcptr[i].startpophg, tvcptr[i].finpophg,tvcptr[i].migcounthg);
    fprintf (treeprintfile,"\t%3d\t%d",tvcptr[i].migcount, tvcptr[i].flag);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      fprintf (treeprintfile, "\t%3d\t%7.15f\n", tvcptr[i].A[0],
               tvcptr[i].dlikeA[0]);
    else
      fprintf (treeprintfile, "\n");
    if (hiddenoptions[HIDDENGENEALOGY]==0)
    {
      for (j = 0; j < totalm; j++)
      {
        if (tvmptr[j].nodenum == tvcptr[i].nodenum)
        {
          fprintf (treeprintfile, "\t%3d\t\t\t\t%7.4f\t\t\t%3d\t%3d\n",
                   tvmptr[j].nodenum, tvmptr[j].dtime,
                   tvmptr[j].startpop, tvmptr[j].finpop);
        }
      }
    }
    else
    {
      j = 0;
      for (j = 0; j < totalmhg; j++)
      {
        if (tvmhgptr[j].nodenum == tvcptr[i].nodenum)
        {
          found = 0;
          for (nm = 0; nm < totalm; nm++)
          {
            if (tvmptr[nm].nodenum == tvcptr[i].nodenum  && tvmptr[nm].dtime==tvmhgptr[j].dtime )
            {
              tvmptr[nm].check = 1;
              fprintf (treeprintfile, "\t%3d\t\t\t\t%7.14f\t\t\t%3d\t%3d\t%3d\t%3d\n",
                  tvmptr[nm].nodenum, tvmptr[nm].dtime,
                  tvmptr[nm].startpop, tvmptr[nm].finpop,tvmhgptr[j].startpop, tvmhgptr[j].finpop);
              found = 1;
            }
          }
          if (found == 0)
          {
            fprintf (treeprintfile, "\t%3d\t\t\t\t%7.14f\t\t\t\t\t%3d\t%3d\t\t\tH\n",
                tvmhgptr[j].nodenum, tvmhgptr[j].dtime,
                tvmhgptr[j].startpop, tvmhgptr[j].finpop);
          }
        }
      }
    }
  }
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Migration counts :\n");
  fprintf (treeprintfile,"Dir\t#");
  if (hiddenoptions[HIDDENGENEALOGY]==1)
    fprintf(treeprintfile,"\t#hg");
  fprintf (treeprintfile, "\n");

  for (i = 0; i < nummigrateparams; i++)
  {
    fprintf (treeprintfile, "  %s\t", C[ci]->imig[i].str);
    fromi = atoi (&C[ci]->imig[i].str[1]);
    toi = atoi (&C[ci]->imig[i].str[3]);
    for (j = 0, k = 0; j < totalm; j++)
    {
      k += (fromi == tvmptr[j].startpop && toi == tvmptr[j].finpop);
    }
    fprintf (treeprintfile, "%d", k);
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      for (j = 0, k = 0; j < totalmhg; j++)
      {
        k += (fromi == tvmhgptr[j].startpop && toi == tvmhgptr[j].finpop);
      }
      fprintf (treeprintfile, "\t%d", k);
    }
    fprintf (treeprintfile, "\n");
  }
  fprintf (treeprintfile, "\n");
  FCLOSE (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvcptr);
  XFREE (tvmptr);
  if (hiddenoptions[HIDDENGENEALOGY]==1)
    XFREE(tvmhgptr);
}                               /* gtreeprint_printlongnumbers */

void
gtreeprint (int ci, int li/*, int step , int callsource */ )
/* use this in debugging mode to see what a genealogy looks like for a particular parameter set */
{
  int p, i, j, k, site, node, up1, up2, nm, totalm;
  int startpop, periodi;
  double upt;
  struct treevals *tvcptr;
  struct treevals *tvmptr;
  char sc[MAXPERIODS][5];
  int checkpt[MAXPERIODS] = {0}; // added initialization 8/18/2016, might have been a bug
  FILE *treeprintfile;
  char treefilename[] = "gtreeprint.out";
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  int fromi, toi;
  /* hgstuff */
  int totalmhg;
  int nmhg,khg,found;
  struct treevals *tvmhgptr;
  int checkfpop;

/*	switch (callsource)
		{
		case 8 : strcpy(sourcestring,"make_JOINT_IS_SW()"); break;
		case 7 : strcpy(sourcestring,"makeSW()"); break;
		case 6 : strcpy(sourcestring,"setup_chains()"); break;
		case 5 : strcpy(sourcestring,"addmigration()"); break;
		case 3 : strcpy(sourcestring,"makeIS()"); break;
		case 4 : strcpy(sourcestring,"main()"); break;
		case 0 : strcpy(sourcestring,"updategenealogy()"); break;
		case 1 : strcpy(sourcestring,"treeweight()"); break;
		case 2 : strcpy(sourcestring,"integrate_tree_prob()"); break;
		default : strcpy(sourcestring,"UNKNOWN");
		} */
  for (i = 0; i < lastperiodnumber; i++)
  {
    sprintf (sc[i], "t%d>", i);
    checkpt[i] = 1;
  }
  totalm = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    j = 0;
    while (gtree[i].mig[j].mt > -0.5)
      j++;
    totalm += j;
  }
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    totalmhg = 0;
    for (i = 0; i < 2 * L[li].numgenes - 1; i++)
    {
      j = 0;
      while (gtree[i].mighg[j].mt > -0.5)
        j++;
      totalmhg += j;
    }
  }
  tvcptr = static_cast<treevals *>
	  (malloc ((2 * L[li].numgenes - 1) * (sizeof (struct treevals))));
  tvmptr = static_cast<treevals *>
	  (malloc (totalm * sizeof (struct treevals)));
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    tvmhgptr = static_cast<treevals *>
	    (malloc (totalmhg * sizeof (struct treevals)));
  }
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
    tvcptr[i].mut = 0;
  for (site = 0; site < L[li].numsites; site++)
  {
    for (j = 0; j < L[li].numgenes; j++)
      gtree[j].mut = L[li].seq[j][site];
    for (j = L[li].numgenes; j < L[li].numlines; j++)
      gtree[j].mut = -1;
    for (j = 0; j < L[li].numgenes; j++)
      labelgtree (ci, li, j);
    for (node = L[li].numgenes; node < L[li].numlines; node++)
    {
      i = node - L[li].numgenes;
      up1 = gtree[node].up[0];
      up2 = gtree[node].up[1];
      if (gtree[node].down != -1)
        if (i != gtree[gtree[node].down].up[0]
            && i != gtree[gtree[node].down].up[1])
          tvcptr[node].flag = 2;
      if (gtree[node].up[0] != -1)
        if (i != gtree[gtree[node].up[0]].down)
          tvcptr[node].flag = 3;
      if (gtree[node].up[1] != -1)
        if (i != gtree[gtree[node].up[1]].down)
          tvcptr[node].flag = 4;
      if ((gtree[up1].mut == 0
           && gtree[up2].mut == 1)
          || (gtree[up1].mut == 1 && gtree[up2].mut == 0))
      {
        if (gtree[node].down == -1)     /* root - not clear where mutation is , put it on left */
        {
          tvcptr[up1].mut++;
        }
        else
        {
          if (gtree[gtree[node].down].mut == gtree[up1].mut)
          {
            tvcptr[up2].mut++;
          }
          else
          {
            tvcptr[up1].mut++;
          }
        }
      }
    }
  }

  for (i = 0, k = 0,khg=0; i < 2 * L[li].numgenes - 1; i++)
  {
    tvcptr[i].nodenum = i;
    tvcptr[i].c_or_m = 0;
    tvcptr[i].flag = 0;
    tvcptr[i].up1 = gtree[i].up[0];
    tvcptr[i].up2 = gtree[i].up[1];
    tvcptr[i].down = gtree[i].down;
    tvcptr[i].dtime = gtree[i].time;
    if (i >= L[li].numgenes)
      upt = gtree[gtree[i].up[0]].time;
    else
      upt = 0;
    tvcptr[i].utime = upt;
    tvcptr[i].timei = gtree[i].time - upt;
    tvcptr[i].startpop = gtree[i].pop;
    tvcptr[i].finpop =  gtree[i].fpop;
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      tvcptr[i].startpophg = gtree[i].pophg;
    j = 0;
    while (gtree[i].mig[j].mt > -1)
      j++;
    tvcptr[i].migcount = j;

    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      nmhg = 0;
      while (gtree[i].mighg[nmhg].mt > -1)
        nmhg++;
      tvcptr[i].migcounthg = nmhg;
      if (nmhg > 0)
        tvcptr[i].finpophg = gtree[i].mighg[nmhg - 1].mp;
      else
        tvcptr[i].finpophg = tvcptr[i].startpophg;
    }

    if (j > 0)
      checkfpop = gtree[i].mig[j - 1].mp;
    else
      checkfpop = tvcptr[i].startpop;
    while (tvcptr[i].dtime > C[ci]->poptree[checkfpop].time)
    {
      checkfpop = C[ci]->poptree[checkfpop].down;
    }

    if (checkfpop != tvcptr[i].finpop || tvcptr[i].finpop >= numtreepops)
      tvcptr[i].flag = 5;
    if (L[li].model == STEPWISE)
    {
      tvcptr[i].A[0] = gtree[i].A[0];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[0];
    }
    if (L[li].model == JOINT_IS_SW)
    {
      tvcptr[i].A[0] = gtree[i].A[1];
      tvcptr[i].dlikeA[0] = gtree[i].dlikeA[1];
    }
    nm = j;
    if (nm > 0)
      for (startpop = gtree[i].pop, j = 0; j < nm; j++)
      {
        tvmptr[k].nodenum = i;
        tvmptr[k].c_or_m = 1;
        tvmptr[k].dtime = gtree[i].mig[j].mt;
        periodi = findperiod (ci, tvmptr[k].dtime);

        /* could change population by passing into the next period */
        while (C[ci]->poptree[startpop].e <= periodi
               && C[ci]->poptree[startpop].e != -1)
          startpop = C[ci]->poptree[startpop].down;
        tvmptr[k].startpop = startpop;
        tvmptr[k].finpop = gtree[i].mig[j].mp;
        tvmptr[k].check = 0;
        startpop = tvmptr[k].finpop;
        k++;
      }
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      if (nmhg> 0)
        for (startpop = gtree[i].pophg, j = 0; j < nmhg; j++)
        {
          tvmhgptr[khg].nodenum = i;
          tvmhgptr[khg].c_or_m = 1;
          tvmhgptr[khg].dtime = gtree[i].mighg[j].mt;
          tvmhgptr[khg].startpop = startpop;
          tvmhgptr[khg].finpop = gtree[i].mighg[j].mp;
          tvmhgptr[khg].check = 0;
          startpop = tvmhgptr[khg].finpop;
          khg++;
          assert(khg<= totalmhg);

        }
    }
  }
  shellgtreevals (tvcptr, L[li].numgenes);
  shellgtreevals (tvcptr + L[li].numgenes, L[li].numgenes - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  /*fprintf (treeprintfile,
           "chain: %d Locus: %d  Step: %d  Calling function: %s\n", ci, li,
           step, sourcestring);
           */
  fprintf (treeprintfile,  "chain: %d Locus: %d  Step: %d   Population tree: %s\n", ci, li,step,C[ci]->chainpoptreestring);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.5f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  if (hiddenoptions[HIDDENGENEALOGY]==0)
    fprintf (treeprintfile, "  current p(D|G): %.6f current p(G) (joint all loci) %.6f\n", C[ci]->G[li].pdg , C[ci]->allpcalc.probg);
  else
    fprintf (treeprintfile, "  current p(D|G): %.6f current p(G_h|G) (joint all loci) %.6f  current p(G) (joint all loci) %.6f\n", C[ci]->G[li].pdg , C[ci]->allpcalc.probhgg, C[ci]->allpcalc.probg);
  if (hiddenoptions[HIDDENGENEALOGY]==0)
  {
    fprintf (treeprintfile,
             "\tNode#\tup1\tup2\tdown\tdtime\tutime\ttimei\tpop\tfinpop\t#mig\tflag");
  }
  else
  {
    fprintf (treeprintfile,
      "\tNode#\tup1\tup2\tdown\tdtime\tutime\ttimei\tpop\tfinpop\tpophg\tfinpophg\t#mighg\t#mig\tflag");
  }
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    fprintf (treeprintfile, "\tA\tld\n");
  else
    fprintf (treeprintfile, "\n");
  p = 0;
  for (i = 0; i < 2 * L[li].numgenes - 1; i++)
  {
    while (tvcptr[i].dtime > C[ci]->tvals[p] && checkpt[p])

    {
      if (hiddenoptions[HIDDENGENEALOGY]==0)
        fprintf (treeprintfile,
                 "%s==============================================================================================================\n",
                 sc[p]);
      else
        fprintf (treeprintfile,
                 "%s===========================================================================================================================================\n",
                 sc[p]);

      checkpt[p] = 0;
      p++;
    }
    //  assert(tvcptr[i].down == -1 || gtree[tvcptr[i].down].pop == tvcptr[i].finpop);
    if (tvcptr[i].down != -1)
      fprintf (treeprintfile,
             "\t%3d\t%3d\t%3d\t%3d\t%7.5f\t%7.5f\t%7.5f\t%3d\t%3d",
             tvcptr[i].nodenum, tvcptr[i].up1, tvcptr[i].up2,
             tvcptr[i].down, tvcptr[i].dtime, tvcptr[i].utime,
             tvcptr[i].timei, tvcptr[i].startpop, tvcptr[i].finpop);
    else
      fprintf (treeprintfile,
             "\t%3d\t%3d\t%3d\t%3d\t%7.1f\t%7.5f\t%7.1f\t%3d\t%3d",
             tvcptr[i].nodenum, tvcptr[i].up1, tvcptr[i].up2,
             tvcptr[i].down, tvcptr[i].dtime, tvcptr[i].utime,
             tvcptr[i].timei, tvcptr[i].startpop, tvcptr[i].finpop);
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      fprintf (treeprintfile,"\t%3d\t%3d\t%d",tvcptr[i].startpophg, tvcptr[i].finpophg,tvcptr[i].migcounthg);
    fprintf (treeprintfile,"\t%3d\t%d",tvcptr[i].migcount, tvcptr[i].flag);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      fprintf (treeprintfile, "\t%3d\t%7.5f\n", tvcptr[i].A[0],
               tvcptr[i].dlikeA[0]);
    else
      fprintf (treeprintfile, "\n");
    if (hiddenoptions[HIDDENGENEALOGY]==0)
    {
      for (j = 0; j < totalm; j++)
      {
        if (tvmptr[j].nodenum == tvcptr[i].nodenum)
        {
          fprintf (treeprintfile, "\t%3d\t\t\t\t%7.5f\t\t\t%3d\t%3d\n",
                   tvmptr[j].nodenum, tvmptr[j].dtime,
                   tvmptr[j].startpop, tvmptr[j].finpop);
        }
      }
    }
    else
    {
      j = 0;
      for (j = 0; j < totalmhg; j++)
      {
        if (tvmhgptr[j].nodenum == tvcptr[i].nodenum)
        {
          found = 0;
          for (nm = 0; nm < totalm; nm++)
          {
            if (tvmptr[nm].nodenum == tvcptr[i].nodenum  && tvmptr[nm].dtime==tvmhgptr[j].dtime )
            {
              tvmptr[nm].check = 1;
              fprintf (treeprintfile, "\t%3d\t\t\t\t%7.5f\t\t\t%3d\t%3d\t%3d\t%3d\n",
                  tvmptr[nm].nodenum, tvmptr[nm].dtime,
                  tvmptr[nm].startpop, tvmptr[nm].finpop,tvmhgptr[j].startpop, tvmhgptr[j].finpop);
              found = 1;
            }
          }
          if (found == 0)
          {
            fprintf (treeprintfile, "\t%3d\t\t\t\t%7.5f\t\t\t\t\t%3d\t%3d\t\t\tH\n",
                tvmhgptr[j].nodenum, tvmhgptr[j].dtime,
                tvmhgptr[j].startpop, tvmhgptr[j].finpop);
          }
        }
      }
    }
  }
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Migration counts :\n");
  fprintf (treeprintfile,"Dir\t#");
  if (hiddenoptions[HIDDENGENEALOGY]==1)
    fprintf(treeprintfile,"\t#hg");
  fprintf (treeprintfile, "\n");

  for (i = 0; i < nummigrateparams; i++)
  {
    fprintf (treeprintfile, "  %s\t", C[ci]->imig[i].str);
    fromi = atoi (&C[ci]->imig[i].str[1]);
    toi = atoi (&C[ci]->imig[i].str[3]);
    for (j = 0, k = 0; j < totalm; j++)
    {
      k += (fromi == tvmptr[j].startpop && toi == tvmptr[j].finpop);
    }
    fprintf (treeprintfile, "%d", k);
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      for (j = 0, k = 0; j < totalmhg; j++)
      {
        k += (fromi == tvmhgptr[j].startpop && toi == tvmhgptr[j].finpop);
      }
      fprintf (treeprintfile, "\t%d", k);
    }
    fprintf (treeprintfile, "\n");
  }
  fprintf (treeprintfile, "\n");
  FCLOSE (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvcptr);
  XFREE (tvmptr);
  if (hiddenoptions[HIDDENGENEALOGY]==1)
    XFREE(tvmhgptr);
}                               /* gtreeprint*/

void
poptreeprint (int ci)
/* use this in debugging mode to see what a population tree looks like */
{
  int i;
  double upt;
  struct treevals *tvptr;
  FILE *treeprintfile;
  char treefilename[] = "poptreeprint.out";
  tvptr = static_cast<treevals *>
	  (malloc ((2 * npops - 1) * (sizeof (struct treevals))));

  for (i = 0; i < 2 * npops - 1; i++)
  {
    tvptr[i].nodenum = i;
    tvptr[i].up1 = C[ci]->poptree[i].up[0];
    tvptr[i].up2 = C[ci]->poptree[i].up[1];
    tvptr[i].down = C[ci]->poptree[i].down;
    tvptr[i].dtime = C[ci]->poptree[i].time;
    if (i >= npops)
      upt = C[ci]->poptree[C[ci]->poptree[i].up[0]].time;
    else
      upt = 0;
    tvptr[i].utime = upt;
    tvptr[i].timei = C[ci]->poptree[i].time - upt;
  }
  shellgtreevals (tvptr, 2 * npops - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  fprintf (treeprintfile, "Chain: %d  Step: %d  Population tree: %s\n", ci, step,C[ci]->chainpoptreestring);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.4f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Node#\tup1\tup2\tdown\tdtime\tutime\ttimei\n");
  for (i = 0; i < 2 * npops - 1; i++)
  {
    if (tvptr[i].down != -1)
      fprintf (treeprintfile, "%3d\t%3d\t%3d\t%3d\t%7.5f\t%7.4f\t%7.4f\n",
             tvptr[i].nodenum, tvptr[i].up1, tvptr[i].up2,
             tvptr[i].down, tvptr[i].dtime, tvptr[i].utime, tvptr[i].timei);
    else
      fprintf (treeprintfile, "%3d\t%3d\t%3d\t%3d\t%7.1f\t%7.4f\t%7.1f\n",
             tvptr[i].nodenum, tvptr[i].up1, tvptr[i].up2,
             tvptr[i].down, tvptr[i].dtime, tvptr[i].utime, tvptr[i].timei);

  }
  fprintf (treeprintfile, "\n");
  FCLOSE (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvptr);
}                               /* poptreeprint */

void
poptreeprint_frompointer (int ci, struct popedge *poptree/*, int step*/)
/* use this in debugging mode to see what a genealogy looks like for a particular parameter set */
{
  int i;
  double upt;
  struct treevals *tvptr;
  FILE *treeprintfile;
  char treefilename[] = "poptreeprint.out";
  tvptr = static_cast<treevals *>
	  (malloc ((2 * npops - 1) * (sizeof (struct treevals))));
  for (i = 0; i < 2 * npops - 1; i++)
  {
    tvptr[i].nodenum = i;
    tvptr[i].up1 = poptree[i].up[0];
    tvptr[i].up2 = poptree[i].up[1];
    tvptr[i].down = poptree[i].down;
    tvptr[i].dtime = poptree[i].time;
    if (i >= npops)
      upt = poptree[poptree[i].up[0]].time;
    else
      upt = 0;
    tvptr[i].utime = upt;
    tvptr[i].timei = poptree[i].time - upt;
  }
  shellgtreevals (tvptr, 2 * npops - 1);
  if (!(treeprintfile = fopen (treefilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening tree file for appending");
  }
  fprintf (treeprintfile, "Chain: %d  Step: %d  Population tree: %s\n", ci, step,C[ci]->chainpoptreestring);
  fprintf (treeprintfile, "   split times:");
  for (i = 0; i < lastperiodnumber; i++)
    fprintf (treeprintfile, "  %8.4f", C[ci]->tvals[i]);
  fprintf (treeprintfile, "\n");
  fprintf (treeprintfile, "Node#\tup1\tup2\tdown\tdtime\tutime\ttimei\n");
  for (i = 0; i < 2 * npops - 1; i++)
  {
    if (tvptr[i].down != -1)
      fprintf (treeprintfile, "%3d\t%3d\t%3d\t%3d\t%7.4f\t%7.4f\t%7.4f\n",
             tvptr[i].nodenum, tvptr[i].up1, tvptr[i].up2,
             tvptr[i].down, tvptr[i].dtime, tvptr[i].utime, tvptr[i].timei);
    else
      fprintf (treeprintfile, "%3d\t%3d\t%3d\t%3d\t%7.1f\t%7.4f\t%7.1f\n",
             tvptr[i].nodenum, tvptr[i].up1, tvptr[i].up2,
             tvptr[i].down, tvptr[i].dtime, tvptr[i].utime, tvptr[i].timei);

  }
  fprintf (treeprintfile, "\n");
  FCLOSE (treeprintfile);
  treeprintfile = NULL;
  XFREE (tvptr);
}                               /* poptreeprint_frompointer */

/* 1 check migs, but not length , 2 check mighg but not  mig, 3 check neither mig or mighg, 4 check it all except mignum, 0 check everything including lengths*/
#define epsilon 0.0000000001  //1e-10
void checkgenealogy(int ci, int li, int mode )
{

  int i,j,jj,k;
  struct genealogy *G = &(C[ci]->G[li]);
  struct edge *gtree = G->gtree;
  double sumlength = 0.0;
  int mignum = 0,mignumhg = 0;
  for (i= 0;i< L[li].numlines;i++)
  {
    if (i < L[li].numgenes)
    {
      assert(gtree[i].up[0] == -1);
      assert(gtree[i].up[1] == -1);
      assert (gtree[i].pop>= 0 &&  gtree[i].pop < npops);
      assert (ISELEMENT (gtree[i].pop, C[ci]->periodset[0]));
      if (hiddenoptions[HIDDENGENEALOGY]==1)
        assert(gtree[i].pophg ==gtree[i].pop );
    }
    else
    {
      assert (ISELEMENT (gtree[i].pop, C[ci]->periodset[findperiod(ci,gtree[gtree[i].up[0]].time)]));
    }
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      assert (gtree[i].pophg>= 0 &&  gtree[i].pophg < npops);
      assert (gtree[i].fpophg>= 0 &&  gtree[i].pophg < npops);
      //assert (gtree[i].cmm <= gtree[i].cmmhg);
    }
    assert (gtree[i].fpop < numtreepops);
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      assert (gtree[i].fpophg < npops);

    if (gtree[i].down != -1)
    {
      assert(gtree[i].down >= L[li].numgenes);
      assert(i == gtree[gtree[i].down].up[0] || i==gtree[gtree[i].down].up[1]);
      assert(gtree[i].time < gtree[gtree[i].down].time);
      assert (ISELEMENT (gtree[i].fpop, C[ci]->periodset[findperiod(ci,gtree[i].time)]));
      if (i < L[li].numgenes)
        sumlength += gtree[i].time;
      else
        sumlength += gtree[i].time - gtree[gtree[i].up[0]].time;
    }
    else
    {
      assert(gtree[i].mig[0].mt < -0.5);
      if (hiddenoptions[HIDDENGENEALOGY]==1)
        assert(gtree[i].mighg[0].mt < -0.5);
      assert (checkfloatsclose(C[ci]->G[li].roottime,gtree[gtree[i].up[0]].time));
    }
    assert(gtree[i].time > 0.0);
    if (mode != 2  && mode != 3)
    {
      j = 0;
      while (gtree[i].mig[j].mt > -0.5)
      {
        assert(j == 0 || gtree[i].mig[j].mt > gtree[i].mig[j-1].mt);
        j++;
      }
      assert (gtree[i].cmm == j);
      j=0;
      k = gtree[i].pop;
      while (gtree[i].mig[j].mt > -0.5)
      {
        assert(gtree[i].mig[j].mp != k);

        k = gtree[i].mig[j].mp;
        assert (ISELEMENT (gtree[i].mig[j].mp, C[ci]->periodset[findperiod(ci,gtree[i].mig[j].mt)]));
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          assert (gtree[i].mighg[j].mt > -0.5);
          jj = 0;
          while (gtree[i].mighg[jj].mt > -0.5)
          {
            assert(jj==0 || gtree[i].mighg[jj].mt > gtree[i].mighg[jj-1].mt);
            jj++;
          }
          assert (gtree[i].cmmhg == jj);
          /* check to see if mig and mighg match up for populations */
          jj = 0;
          while (gtree[i].mighg[jj].mt > -0.5 && gtree[i].mighg[jj].mt != gtree[i].mig[j].mt )
            jj++;
          assert(gtree[i].mighg[jj].mt > -0.5 &&  gtree[i].mighg[jj].mt == gtree[i].mig[j].mt  );
          assert(gtree[i].mig[j].mp == C[ci]->ancplist[gtree[i].mighg[jj].mp][ findperiod(ci,gtree[i].mig[j].mt)]);
        }
        j += 1;
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          jj = 0;
          while (gtree[i].mighg[jj].mt > -0.5 && gtree[i].mighg[jj].mt != gtree[i].mig[j].mt )
            jj++;
          assert(jj>= j);
        }
      }
      mignum += j;
    }
    if (hiddenoptions[HIDDENGENEALOGY]==1 && mode != 3 && mode != 4)
    {
      j=0;
      k = gtree[i].pophg;
      while (gtree[i].mighg[j].mt > -0.5)
      {
        assert(gtree[i].mighg[j].mp != k);
        k = gtree[i].mighg[j].mp;
        assert (ISELEMENT (gtree[i].mighg[j].mp, C[ci]->periodset[0]));
        j += 1;
      }
      mignumhg +=j;
      if (gtree[i].down != -1)
        assert(gtree[gtree[i].down].pophg == gtree[i].fpophg);
    }
  }
  if (mode == 0)
  {
    assert(G->mignum == mignum);
    if (calcoptions[DONTCALCGENEALOGYPRIOR]==0)
      assert((fabs(G->length-sumlength)/sumlength ) < epsilon);  // this will fail because length is not set in treeweight with  calcoptions[DONTCALCGENEALOGYPRIOR]==1
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      assert(mignumhg >= mignum);
    }
  }

}    /* checkgenealogy */

void checkpoptree(int ci,int mode) // mode: 0 do all checks, 1 skip period check, 2 also skip node order checks
{
  int i,period;
  double upt;
  char temppoptreestring[POPTREESTRINGLENGTHMAX];
  int  temppoptreenum;
  if (mode==0 && modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    poptreewrite (ci, temppoptreestring);
    temppoptreenum = getpoptreestringnum(temppoptreestring);
    assert(strcmp(C[ci]->chainpoptreestring,temppoptreestring)==0);
    assert (C[ci]->poptreenum ==temppoptreenum);
  }

  for (i=0;i<numtreepops;i++)
  {
    if (i<npops)
    {
      assert(C[ci]->poptree[i].up[0]==-1);
      assert(C[ci]->poptree[i].up[1]==-1);
      upt = 0.0;
    }
    else
    {
      assert(C[ci]->poptree[C[ci]->poptree[i].up[0]].down == i);
      assert(C[ci]->poptree[C[ci]->poptree[i].up[1]].down == i);
      upt = C[ci]->poptree[C[ci]->poptree[i].up[0]].time;

    }
    if (mode==0)
    {
      assert (C[ci]->poptree[i].b == findperiod(ci,upt));
     // poptreeprint(ci);
      assert ((C[ci]->poptree[i].e == -1 && C[ci]->poptree[i].time == TIMEMAX) || (C[ci]->poptree[i].e == findperiod(ci,C[ci]->poptree[i].time)));
      assert (C[ci]->poptree[i].e == -1 || C[ci]->poptree[i].e == findperiod(ci,C[ci]->poptree[i].time));
      for (period = 0;period < lastperiodnumber;period++)  // identify which time point in tvals
        if (C[ci]->poptree[i].time==C[ci]->tvals[period])
          break;
      assert(C[ci]->poptree[i].time == TIMEMAX || period < lastperiodnumber);
      assert(C[ci]->poptree[i].time == TIMEMAX || (C[ci]->poptree[i].time <= T[period].pr.max && C[ci]->poptree[i].time > T[period].pr.min));
    }

    if (C[ci]->poptree[i].down != -1)
    {
      if (mode < 2)
        assert (i < C[ci]->poptree[i].down);
      assert(C[ci]->poptree[i].time < C[ci]->poptree[C[ci]->poptree[i].down].time);
      assert(C[ci]->poptree[C[ci]->poptree[i].down].up[0] == i || C[ci]->poptree[C[ci]->poptree[i].down].up[1] == i );
    }
    if (mode == 0)
    {
      period = findperiod(ci,upt);
      if (i<npops)
        assert(period == 0);
      else
        assert(period +npops-1 == i);
    }
  }

} /* checkpoptree */

#undef epsilon


/* recalculate the likelihood,  and the main priors
see if results are valid numbers and if they match what is currently held */
void checkprobs(int ci, int locusi)
{
  int li,lbi,lui,ai = 0;
  double pdg, phggg;
  struct probcalc pcalc;
  struct genealogy_weights holdgweight;
  double sumpdg = 0.0;
  double sumphggg = 0.0;
  int tw;
  assert(ci < numchainspp);
  assert(locusi < nloci);
  if (locusi < 0)
  {
    lbi = 0;
    lui = nloci;
  }
  else
  {
    lbi = locusi;
    lui = lbi + 1;
  }

  init_probcalc (&pcalc);
  init_genealogy_weights (&holdgweight);
  for (li=lbi;li<lui;li++)
  {

    switch (L[li].model)
    {
    case HKY:
        pdg  = likelihoodHKY (ci, li, C[ci]->G[li].uvals[0], C[ci]->G[li].kappaval, -1, -1, -1, -1);
      break;
    case INFINITESITES:
        pdg = likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
      break;
    case STEPWISE:
      pdg = 0;
      for (ai=0; ai < L[li].nlinked; ai++)
      {
        pdg += likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      }
      break;
    case JOINT_IS_SW:
      pdg = likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
      for (ai = 1; ai < L[li].nlinked; ai++)
      {
        pdg += likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
      }
      break;
    }
    assert( isnan_(pdg)==0); // chec if is nan
    assert( isnotinf_DBL(pdg));
    if (pdg == 0.0)
      assert(C[ci]->G[li].pdg==0.0);
    else
      assert (checkfloatsclose(C[ci]->G[li].pdg,pdg));
    sumpdg += pdg;

    /* to check probg have to check the overall probg while changing just locus li */

    copy_treeinfo (&holdgweight,  &C[ci]->G[li].gweight);
    setzero_genealogy_weights (&C[ci]->G[li].gweight);
    int treeweightcallcode = 9;
    tw = treeweight (ci, li,treeweightcallcode);
    sum_subtract_treeinfo (&C[ci]->allgweight, &C[ci]->G[li].gweight,&holdgweight);
    initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &pcalc);

    assert( isnan_(pcalc.probg)==0); // chec if is nan
    assert( isnotinf_DBL(pcalc.probg));
    if (calcoptions[DONTCALCGENEALOGYPRIOR]==0)
      assert( pcalc.probg != 0.0);


    assert(checkfloatsclose(C[ci]->allpcalc.probg,pcalc.probg));
   // assert(fclose05(C[ci]->allpcalc.probg,pcalc.probg));// checks within 5%

    if (hiddenoptions[HIDDENGENEALOGY])
    {
      phggg = prob_hg_given_g(ci,li);
      assert( isnan_(phggg)==0); // chec if is nan
      assert( isnotinf_DBL(phggg));
      if (phggg == 0.0)
        assert(C[ci]->G[li].hgprob==0.0);
      else
        assert (checkfloatsclose(C[ci]->G[li].hgprob,phggg));

      sumphggg += phggg;
    }
  }
  if (locusi < 0)   // multiple loci, check summed values
  {
    assert( isnan_(sumpdg)==0); // chec if is nan
    assert( isnotinf_DBL(sumpdg));
    if (calcoptions[DONTCALCLIKELIHOODMUTATION] != 1)
    {
      assert(sumpdg != 0.0);
      assert(checkfloatsclose(C[ci]->allpcalc.pdg,sumpdg));
    }
    if (hiddenoptions[HIDDENGENEALOGY])
    {
      assert( isnan_(sumphggg)==0); // chec if is nan
      assert( isnotinf_DBL(sumphggg));
      if (sumphggg == 0.0)
        assert(C[ci]->allpcalc.probhgg==0.0);
      else
        assert (checkfloatsclose(C[ci]->allpcalc.probhgg,sumphggg));
    }
  }
  free_probcalc (&pcalc);
  free_genealogy_weights (&holdgweight);
}    /* checkprobs() */


/* print to a file the current genealogy weights used for calculated probg */
#define   fpg  fprintf(gweightprintfile,
void printgenealogyweights(int ci,int li)
{
  struct genealogy_weights *gweight;
  FILE *gweightprintfile;
  char gweightfilename[] = "gweightprint.out";
  int i,j,c,csum;
  double f,hc,fsum,hcsum;
  int isall = (li<0);

  if (!(gweightprintfile = fopen (gweightfilename, "a")))
  {
    IM_err(IMERR_APPENDFILEFAIL,"Error opening gweight file for appending");
  }
  if (isall)
  {
    fpg "Step: %d Chain: %d Genealogy weights, summed over loci\n",step,ci);
    gweight = &C[ci]->allgweight;
  }
  else
  {
    fpg "Step: %d  Chain: %d  Locus: %d  Genealogy weights\n", step, ci, li);
    gweight = &C[ci]->G[li].gweight;
  }
  fpg "Population size parameters\n");
  for (i = 0; i < numpopsizeparams; i++)
  {
    csum = 0;
    fsum = hcsum = 0.0;
    fpg " pop %d #periods %d\n",i,C[ci]->itheta[i].wp.n);
    for (j = 0; j < C[ci]->itheta[i].wp.n; j++)
    {
      c = (int) gweight->cc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      csum += c;
      f = (float) gweight->fc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      fsum += f;
      hc = (float) gweight->hcc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      hcsum += hc;
      fpg "   j %d C[ci]->itheta[i].wp.p[j] %d C[ci]->itheta[i].wp.r[j] %d c %d    f %.5lf   hc %.5lf\n",j, C[ci]->itheta[i].wp.p[j],C[ci]->itheta[i].wp.r[j],c,f,hc);

    }
      fpg "   sum   c %d    f %.5lf   hc %.5lf",csum,fsum,hcsum);
      if (isall)
        fpg "  pcalc %.5lf\n",C[ci]->allpcalc.qintegrate[i]);
      else
        fpg "\n");
  }
  if (!modeloptions[NOMIGRATION])
  {
    fpg "Migration parameters\n");
    for (i = 0; i < nummigrateparams; i++)
    {
      csum = 0;
      fsum = 0.0;
      fpg " mig param# %d #periods %d\n",i,C[ci]->imig[i].wp.n);
      for (j = 0; j < C[ci]->imig[i].wp.n; j++)
      {
        c = (int) gweight->mc[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].wp.c[j]];
        csum += c;
        f =  (float) gweight->fm[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].wp.c[j]];
        fsum += f;
        fpg "   j %d C[ci]->imig[i].wp.p[j] %d C[ci]->imig[i].wp.r[j] %d C[ci]->imig[i].wp.c[j] %d c %d  f %.5lf\n",j,C[ci]->imig[i].wp.p[j],C[ci]->imig[i].wp.r[j],C[ci]->imig[i].wp.c[j],c,f);

      }
      fpg "   sum   c %d  f %.5lf",csum,fsum);
      if (isall)
        fpg "  pcalc %.5lf\n",C[ci]->allpcalc.mintegrate[i]);
      else
        fpg "\n");
    }
  }
  fpg "===========================================\n\n");
  FCLOSE (gweightprintfile);
  gweightprintfile = NULL;
}  //printgenealogyweights

#undef fpg

/* check that weights add up.  Also check some stuff involving migration hyperpriors  */
void checkgenealogyweights(int ci) // check that sums of weights equal the overall weights
{
  int li, i,j,k;
  int c, ca;
  double f,hc,fa,hca;


  if (calcoptions[DONTCALCGENEALOGYPRIOR]==1)
    return;
  // first check that the total count of coalescent events for each locus is correct
  for (li=0;li<nloci;li++)
  {
    c=0;
    for (i = 0; i < numpopsizeparams; i++)
    {
      for (j = 0; j < C[ci]->itheta[i].wp.n; j++)
      {
        c +=  C[ci]->G[li].gweight.cc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      }
    }
    assert (c== L[li].numgenes - 1);
  }


  for (i = 0; i < numpopsizeparams; i++)
  {
    c = 0;
    f = hc = 0.0;
    for (li=0;li<nloci;li++)
    {
      for (j = 0; j < C[ci]->itheta[i].wp.n; j++)
      {
        c +=  C[ci]->G[li].gweight.cc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
        f +=  C[ci]->G[li].gweight.fc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
        hc +=  C[ci]->G[li].gweight.hcc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      }
    }
    ca = 0;
    fa = hca = 0.0;
    for (j = 0; j < C[ci]->itheta[i].wp.n; j++)
    {
      ca += (int) C[ci]->allgweight.cc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      fa += C[ci]->allgweight.fc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];
      hca += C[ci]->allgweight.hcc[C[ci]->itheta[i].wp.p[j]][C[ci]->itheta[i].wp.r[j]];

    }
    //assert ((C[ci]->allpcalc.qintegrate[i] != 0.0 && (ca > 0 || fa > 0.0)) || (C[ci]->allpcalc.qintegrate[i] == 0.0 && ca == 0 && fa == 0.0));
    // this would only work with prior of 1  - see qintegrate()
    assert (checkfloatsclose( (double) c, (double) ca));
    assert (checkfloatsclose(f,fa));
    assert (checkfloatsclose(hc,hca));
  }
  if (!modeloptions[NOMIGRATION])
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      c = 0;
      f = 0.0;
      for (li=0;li<nloci;li++)
      {
        for (j = 0; j < C[ci]->imig[i].wp.n; j++)
        {
          c +=
            (int) C[ci]->G[li].gweight.mc[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].
                                                                wp.c[j]];
          f += (float)
            C[ci]->G[li].gweight.fm[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].wp.c[j]];
        }
      }
      ca = 0;
      fa = 0.0;
      for (j = 0; j < C[ci]->imig[i].wp.n; j++)
      {
        ca +=
          (int) C[ci]->allgweight.mc[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].
                                                              wp.c[j]];
        fa +=
          C[ci]->allgweight.fm[C[ci]->imig[i].wp.p[j]][C[ci]->imig[i].wp.r[j]][C[ci]->imig[i].wp.c[j]];

      }
      //assert ((C[ci]->allpcalc.mintegrate[i] != 0.0 && (ca > 0 || fa > 0.0)) || (C[ci]->allpcalc.mintegrate[i] == 0.0 && ca == 0 && fa == 0.0));
      // this does not actually does not work  because mintegrate uses a cutoff
      // also would only work with prior of 1  - see qintegrate()
      assert (checkfloatsclose( (double) c, (double) ca));
      assert (checkfloatsclose(f,fa));
    }
    // check some hyperor prior stuff
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]==1) for (i = 0; i < nummigrateparams; i++)
    {

      assert(C[ci]->imig[i].dir == makepairstring(C[ci]->descendantpops[C[ci]->imig[i].md.from],C[ci]->descendantpops[C[ci]->imig[i].md.to],C[ci]->imig[i].descstr));
      if (modeloptions[EXPOMIGRATIONPRIOR]==1)
      {
        if (C[ci]->imig[i].dir==0)
          assert(C[ci]->imig[i].pr.expomean ==getvalue(C[ci]->imig[i].descstr, C[ci]->mltorhpriors));
        else
          assert(C[ci]->imig[i].pr.expomean ==getvalue(C[ci]->imig[i].descstr, C[ci]->mrtolhpriors));
      }
      else
      {
        if (C[ci]->imig[i].dir==0)
          assert(C[ci]->imig[i].pr.max ==getvalue(C[ci]->imig[i].descstr, C[ci]->mltorhpriors));
        else
          assert(C[ci]->imig[i].pr.max ==getvalue(C[ci]->imig[i].descstr, C[ci]->mrtolhpriors));
      }
    }
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]==1) for (i = 0; i < nummigrateparampairs; i++)
    {
      assert (strcmp(C[ci]->imig[2*i].descstr, C[ci]->imig[2*i+1].descstr)==0);
    }
  }
  //check plist and wp.n

  for (int ii = 0; ii < numpopsizeparams; ii++)
  {
    int n = 0;
    for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
    {
      for (j = 0; j < k; j++)
      {
        if (ii == C[ci]->plist[i][j])
          n += 1;
      }
    }
    assert (n==C[ci]->itheta[ii].wp.n);
  }
  return;
}                               /*void checkgenealogyweights(int ci) */

/* simple check of whether C[ci]->allpcalc.probhgg equals the sum of the individual locus values */
void check_hgprob_sums(int ci)
{
  int li;
  double temp,tempnewcalc,tempsum = 0.0;
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    for (li = 0;li<nloci;li++)
    {
      temp =  C[ci]->G[li].hgprob;
      tempnewcalc = prob_hg_given_g(ci,li);
      assert ( checkfloatsortaclose(temp,tempnewcalc));
      tempsum += temp;
    }
    assert ( checkfloatsclose(tempsum,C[ci]->allpcalc.probhgg));
  }
}



void checkdetailedbalance(double newlikelihood, double oldlikelihood, double newprior, double oldprior, double propose_old_given_new, double propose_new_given_old, double beta)
{
	double aforward,abackward;
	double temp;
	double dbforward, dbbackward;
	if (calcoptions[CALCMARGINALLIKELIHOOD])
	{
		temp = beta * (newlikelihood - oldlikelihood) + (newprior - oldprior) + (propose_old_given_new - propose_new_given_old);
		aforward = (temp >= 0.0) ? 0.0 : temp;
		temp = beta * (-(newlikelihood - oldlikelihood))  - (newprior - oldprior) - (propose_old_given_new - propose_new_given_old);
		abackward = (temp >= 0.0) ? 0.0 : temp;
	}
	else
	{
		temp = beta * (newlikelihood - oldlikelihood + newprior - oldprior) + (propose_old_given_new - propose_new_given_old);
		aforward = (temp >= 0.0) ? 0.0 : temp;
		temp = beta * (-(newlikelihood - oldlikelihood + newprior - oldprior)) - (propose_old_given_new - propose_new_given_old);
		abackward = (temp >= 0.0) ? 0.0 : temp;
	}
  if (calcoptions[CALCMARGINALLIKELIHOOD])
  {
	  dbforward = beta * oldlikelihood + oldprior + propose_new_given_old + aforward;
	  dbbackward = beta * newlikelihood + newprior + propose_old_given_new + abackward;
  }
  else
  {
	  dbforward = beta * (oldlikelihood + oldprior) + propose_new_given_old + aforward;
	  dbbackward = beta * (newlikelihood + newprior) + propose_old_given_new + abackward;
  }

	assert(checkfloatsclose(dbforward,dbbackward));
}

void checkdetailedbalance_chainswap(double likelihood_i, double likelihood_k, double prior_i, double prior_k, double beta_i, double beta_k)
{
	double aforward,abackward;
	double temp;
	double dbforward, dbbackward;
	if (calcoptions[CALCMARGINALLIKELIHOOD])
	{
		temp = (beta_i - beta_k) * (likelihood_k - likelihood_i);
		aforward =  (temp >= 0.0) ? 0.0 : temp;
    abackward = ( (-temp) >= 0.0) ? 0.0 : -temp;
	}
	else
	{
    temp = (beta_i - beta_k) * (likelihood_k + prior_k - likelihood_i -  prior_i) ;
		aforward =  (temp >= 0.0) ? 0.0 : temp;
    abackward = ( (-temp) >= 0.0) ? 0.0 : -temp;


	}
  if (calcoptions[CALCMARGINALLIKELIHOOD])
  {
	  dbforward = beta_i * likelihood_i + prior_i + beta_k * likelihood_k + prior_k +  aforward;
	  dbbackward = beta_i * likelihood_k + prior_k + beta_k * likelihood_i + prior_i + abackward;
  }
  else
  {
    dbforward = beta_i *  (likelihood_i + prior_i) + beta_k * (likelihood_k + prior_k) + aforward;
	  dbbackward = beta_i * (likelihood_k + prior_k) + beta_k * (likelihood_i + prior_i) + abackward;

  }
	assert(checkfloatsclose(dbforward,dbbackward));
}


/* check identity of two struct genealogy_weights */
void
compare_genealogy_weights (struct genealogy_weights *gw1,struct genealogy_weights *gw2)
{
  int i,j,k;

  for (i = 0; i < numsplittimes + 1; i++)for (j = 0; j < npops-i; j++)
    assert(gw1->cc[i][j] == gw2->cc[i][j]);
  for (i = 0; i < numsplittimes + 1; i++)for (j = 0; j < npops-i; j++)
    assert(gw1->hcc[i][j] == gw2->hcc[i][j]);
  for (i = 0; i < numsplittimes + 1; i++)for (j = 0; j < npops-i; j++)
    assert(gw1->fc[i][j] == gw2->fc[i][j]);
  if (modeloptions[NOMIGRATION] != 0)
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      for (j=0;j<npops-i;j++)for(k=0;k<npops-i;k++)
      {
        assert(gw1->mc[i][j][k]==gw2->mc[i][j][k]);
        assert(gw1->fm[i][j][k]==gw2->fm[i][j][k]);
      }
    }
  }
} // compare_genealogy_weights()

/* some debugging code to catch probability floating point problems
  simple,  used only after call to an updating function
  this really should not be necessary*/
void pcheck(int ci, int fromcode)
{
  if (isnan_(C[ci]->allpcalc.pdg))
    IM_err(IMERR_MISCPROBPROBLEM,"Likelihood is nan in chain %d, step %d, from function %d",ci,step,fromcode);
  if (isninf_DBL(C[ci]->allpcalc.pdg))
    IM_err(IMERR_MISCPROBPROBLEM,"Likelihood is -inf in chain %d, step %d, from function %d",ci,step,fromcode);
  if (ispinf_DBL(C[ci]->allpcalc.pdg))
    IM_err(IMERR_MISCPROBPROBLEM,"Likelihood is +inf in chain %d, step %d, from function %d",ci,step,fromcode);
  if (isnan_(C[ci]->allpcalc.probg))
    IM_err(IMERR_MISCPROBPROBLEM,"Genealogy prior is nan in chain %d, step %d, from function %d",ci,step,fromcode);
  if (isninf_DBL(C[ci]->allpcalc.probg))
    IM_err(IMERR_MISCPROBPROBLEM,"Genealogy prior is -inf in chain %d, step %d, from function %d",ci,step,fromcode);
  if (ispinf_DBL(C[ci]->allpcalc.probg))
    IM_err(IMERR_MISCPROBPROBLEM,"Genealogy prior is +inf in chain %d, step %d, from function %d",ci,step,fromcode);
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    if (isnan_(C[ci]->allpcalc.probhgg))
      IM_err(IMERR_MISCPROBPROBLEM,"Hidden Genealogy prior is nan in chain %d, step %d, from function %d",ci,step,fromcode);
    if (isninf_DBL(C[ci]->allpcalc.probhgg))
      IM_err(IMERR_MISCPROBPROBLEM,"Hidden Genealogy prior is -inf in chain %d, step %d, from function %d",ci,step,fromcode);
    if (ispinf_DBL(C[ci]->allpcalc.probhgg))
      IM_err(IMERR_MISCPROBPROBLEM,"Hidden Genealogy prior is +inf in chain %d, step %d, from function %d",ci,step,fromcode);
    }
}



#endif //TURNONCHECKS
