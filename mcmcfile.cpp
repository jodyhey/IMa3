/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/*
checkpoint file usage

runoptions[]
    2 Save the state of the Markov chain in a file - named with extension .mcf (MCMC mode only)
    3 Start run by loading a previously saved *.mcf file; requires -f (data and priors must be the same)
    6 With mcf files (-r2,-r3) record and load only the MCMC state space (useful to treat a previous run as a burnin)
    7 Load *.mcf file if present AND save to *.mcf file. Use for repeating command lines.  No burn done if *.mcf present


r2 The mcf file(s) are saved after the burnin in done  (unless -b0) and at the end of the run.
Saves state space in all chains,  with one file per cpu.
If -j0  (population topology updating), it also records the population tree for each chain and the current values of poptopologycounts
Unless -r6 is invoked it also saves update rate stats,  autocorrelation values and (if -j0) tree phylogeny counts

-r3 loads the mcf file(s) using the filename given with -f  and starts a new run.

using -r23 to overwrite previous files  (e.g. to continue a run that has been stopped)
If a previous -r2 run has been done,  then additional -r23 runs can be done using the same previous filename for -o and and using the same name for -f
This causes the previous results file and the previous mcf file(s) to be overwritten
If the burnin is all done,  then the -r23 run should be started using -b0.  This will prevent resetting of accumulated update rate statistics.

-r6 prevents reading of information from the mcf file(s) on update rate statistics,  autocorrelation values, and tree phylogeny counts
This can be useful if the mcf file being loaded is only needed to initialize the state space.
For example, if treating the previous run as just a burnin run and you want to start collecting samples from that point.

The order of things written to the mcf file is nearly arbitrary.  The mcf file can be large.
*/

/* stuff for writing and reading the state of the mcmc to a file
to write the state of the mcmc state space to a file
call writemcmc(int ci)
to load it from a file
call readmcmc(int ci)

writemcmc() and readmcmc() are basically copies of each other for the reading and writing of data.

To simplify writing them and so they can use mostly the same actual text,
a simple macro 'aa' is reset for each to call either aread() or awrite() respectively

readmcmc() also has some initiation stuff

codes for atype
0  for integer
1  for long
2  for float
3  for double
4  for character
5  for unsigned long
6 for unsigned short
*/

int netsteps; //jh 1_17_2018


/* some extern prototypes */
extern void fillplist (int ci);
extern void fillancplist (int ci);

/* some extern data */
extern double *thermosum,*thermosum_rec;
extern int numdistinctpopulationpairs[]; // = {0,0,1,6,25,90,301,966,3025}; /* number of distinct pairs of populations that could engage in gene flow (don't share any descendant pops) */

static void awrite (FILE * mcffile, const char *name, int atype, int iu, void *a);
static void aread (FILE * mcffile, const char *name, int atype, int iucheck, void *a);
static void init_p (void);
void write_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x);
void read_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x);
void write_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x);
void write_autoc_record(FILE *mcffile,char *vstr,struct autoc *ac);
void read_autoc_record(FILE *mcffile,char *vstr,struct autoc *ac);
void read_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x);

/**** LOCAL FUNCTIONS *****/

#define MCFVARNAMEWRITE  // causes names of variables to be written and read along with data, helps with debugging and  has only a small effect on file size so makes sense to leave it in
// don't use if reading in mcf files from IMa2
//#undef MCFVARNAMEWRITE

void
awrite (FILE * mcffile, const char *name, int atype, int iu, void *a)
{
  int *ip;
  long *lip;
  float *fp;
  double *dp;
  char *cp;
  unsigned long *ulp;
  unsigned short *usp;
  int i;
  int isnumerical;

#ifdef MCFVARNAMEWRITE
  fprintf (mcffile, "%s %d %d ", name, atype, iu);
#else
  fprintf (mcffile, "%d %d ", atype, iu);
#endif
  isnumerical = atype != 4;

  switch (atype)
  {
  case 0:
    for (i = 0, ip = static_cast<int *> (a); i < iu; i++)
    {
      if (isnan_(*(ip + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018, probably can't get nan int, but not sure best way to handle this situation
      fprintf (mcffile, "%d ", *(ip + i));
    }
    break;
  case 1:
    for (i = 0, lip = static_cast<long *> (a); i < iu; i++)
    {
      if (isnan_(*(lip + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018, probably can't get nan int, but not sure best way to handle this situation
      fprintf (mcffile, "%ld ", *(lip + i));
    }
    break;
  case 2:
    for (i = 0, fp = static_cast<float *> (a); i < iu; i++)
    {
      if (isnan_(*(fp + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018
      if (isninf_FLT(*(fp + i)))
        fprintf (mcffile, "%.12g ", -FLT_MAX);
      else if (ispinf_FLT(*(fp + i)))
        fprintf (mcffile, "%.12g ", FLT_MAX);
      else
        fprintf (mcffile, "%g ", *(fp + i));
    }
    break;
  case 3:
    for (i = 0, dp = static_cast<double *> (a); i < iu; i++)
    {
      if (isnan_(*(dp + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018
      if (isninf_DBL(*(dp + i)))
        fprintf (mcffile, "%.18lg ", -DBL_MAX);
      else if (ispinf_DBL(*(dp + i)))
        fprintf (mcffile, "%.18lg ", DBL_MAX);
      else
        fprintf (mcffile, "%.18lg ", *(dp + i));
    }
    break;
  case 4:
    for (i = 0, cp = static_cast<char *> (a); i < iu; i++)
      fprintf (mcffile, "%c", *(cp + i));
    if (*(cp + i) != '\0')
      fprintf (mcffile, "%c", '\0');
    break;
  case 5:
    for (i = 0, ulp = static_cast<unsigned long *> (a); i < iu; i++)
    {
      if (isnan_(*(ulp + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018, probably can't get nan int, but not sure best way to handle this situation
      fprintf (mcffile, "%lu ", *(ulp + i));
    }
    break;
  case 6:
    for (i = 0, usp = static_cast<unsigned short *> (a); i < iu; i++)
    {
      if (isnan_(*(usp + i)))
        IM_err(IMERR_MCFWRITEFAIL,"attemp to write non-numerical value to mcf file: %s"); // added 1/18/2018, probably can't get nan int, but not sure best way to handle this situation
      fprintf (mcffile, "%hu ", *(usp + i));
    }
    break;
  }
  fprintf (mcffile, "\n");
}                               /* arraywrite */

void
aread (FILE * mcffile, const char *name, int atype, int iucheck, void *a)
{
  int *ip;
  long *lip;
  float *fp;
  double *dp;
  char *cp;
  unsigned long *ulp;
  unsigned short *usp;
  int typecheck, i, iu;
  char namecheck[100];

  iucheck = 0; /* To remove warning: unused parameter. */
#ifdef MCFVARNAMEWRITE
  scanfval = fscanf (mcffile, "%99s %d %d ", namecheck, &typecheck, &iu);
  if (strcmp (name, namecheck) != 0)
  {
    IM_err(IMERR_MCFREADFAIL,"variable names do not match: %s  <> %s ", name, namecheck);
  }
#else
  scanfval = fscanf (mcffile, "%d %d ", &typecheck, &iu);
#endif

  if (typecheck != atype)
  {
    if (!(typecheck==6 && atype == 0))  // allow a 6/0 mismatch, which is ok when reading. should only be relevant for some mcfs written before 2/2018 as I changed some writes from  6 to 0
      IM_err(IMERR_MCFREADFAIL,"variable types do not match: %d  <> %d ", atype, typecheck);
  }
  switch (atype)
  {
  case 0:
    for (i = 0, ip = static_cast<int *> (a); i < iu; i++, ip++)
      scanfval = fscanf (mcffile, "%d ", ip);
    break;
  case 1:
    for (i = 0, lip = static_cast<long *> (a); i < iu; i++, lip++)
      scanfval = fscanf (mcffile, "%ld ", lip);
    break;
  case 2:
    for (i = 0, fp = static_cast<float *> (a); i < iu; i++, fp++)
      scanfval = fscanf (mcffile, "%g ", fp);
    break;
  case 3:
    for (i = 0, dp = static_cast<double *> (a); i < iu; i++, dp++)
      scanfval = fscanf (mcffile, "%lg ", dp);
    break;
  case 4:
      for (i = 0, cp = static_cast<char *> (a); i < iu; i++, cp++)
        scanfval = fscanf (mcffile, "%c", cp);
      *cp = '\0';
    break;
  case 5:
    for (i = 0, ulp = static_cast<unsigned long *> (a); i < iu; i++, ulp++)
      scanfval = fscanf (mcffile, "%lu ", ulp);
    break;
  case 6:
    for (i = 0, usp = static_cast<unsigned short *> (a); i < iu; i++, usp++)
      scanfval = fscanf (mcffile, "%hu ", usp);
    break;
  }
  return;
}                               /* arrayread */

void
init_p (void)
{
  int i,ci, li, ai;
  for (ci = 0; ci < numchainspp; ci++)
  {
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
       reset_poptree (ci, C[ci]->chainpoptreestring);
       fillplist (ci);
       fillancplist (ci);
       set_iparam_poptreeterms (ci); // 5/11/2017  this seems to be needed here to set the wp.n terms correctly// also fills descendantpops
        for (i = 0; i < 2 * npops - 1; i++)
        {
          if (C[ci]->poptree[i].e == -1)
            C[ci]->poptree[i].time = TIMEMAX;
          else
            C[ci]->poptree[i].time = C[ci]->tvals[C[ci]->poptree[i].e - 1];
        }
#ifdef TURNONCHECKS
       //poptreeprint(ci);
#endif //TURNONCHECKS
    }
    setzero_genealogy_weights (&C[ci]->allgweight);
    C[ci]->allpcalc.pdg = 0;
    if (hiddenoptions[HIDDENGENEALOGY]==1)
      C[ci]->allpcalc.probhgg  = 0.0;
    for (li = 0; li < nloci; li++)
    {
      setzero_genealogy_weights (&C[ci]->G[li].gweight);
      int treeweightcallcode = 0;  // a debugging code,  if treeweight has an error results are written to an output file with this code
      int tw = treeweight (ci, li,treeweightcallcode);
      switch (L[li].model)
      {
      case INFINITESITES:
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        break;
      case HKY:

          C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
            likelihoodHKY (ci, li, C[ci]->G[li].uvals[0],
                            C[ci]->G[li].kappaval, -1, -1, -1, -1);
          copyfraclike (ci, li);
          storescalefactors (ci, li);
        break;
      case STEPWISE:
        somestepwise = 1;
        C[ci]->G[li].pdg = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      case JOINT_IS_SW:
        somestepwise = 1;
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        for (ai = 1; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      }
      sum_treeinfo (&(C[ci]->allgweight), &(C[ci]->G[li].gweight));
      C[ci]->allpcalc.pdg += C[ci]->G[li].pdg;
      if (hiddenoptions[HIDDENGENEALOGY]==1)
      {
#ifdef TURNONCHECKS
        //gtreeprint(ci,li);
#endif //TURNONCHECKS
        C[ci]->G[li].hgprob =  prob_hg_given_g(ci,li);
        C[ci]->allpcalc.probhgg +=   C[ci]->G[li].hgprob;
      }
    }
    initialize_integrate_tree_prob (ci, &(C[ci]->allgweight),
                                    &C[ci]->allpcalc);
    /* if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]) don't need this here,  its done in setup_iparams()
    {
      init_hyperprior_arrays(ci);
      if (ci==ARBCHAIN)
        int numpossiblepoppairstrings = fillmigratepairs();
    } */
  }
}                               // init_p

int findallbetasci(double beta)
{
  int i;
  for (i=0;i<numchainstotal;i++)
  {
    if (allbetas[i] == beta)
      break;
  }
  assert(i < numchainstotal);
  return i;
}

/* write values for an instance of struct autoc */
#ifdef aa
#undef aa
#endif /* aa */
#define aa   awrite (mcffile,
void write_autoc_record(FILE *mcffile,char *vstr,struct autoc *ac)
{
  char tempstr[NAMELENGTH];
  const char *cctempstr;
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_s");
  aa cctempstr,3,1,&(ac->cov.s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_ss");
  aa cctempstr,3,1,&(ac->cov.ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_s2");
  aa cctempstr,3,1,&(ac->cov.s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_n");
  aa cctempstr,0,1,&(ac->cov.n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_s");
  aa cctempstr,3,1,&(ac->var[0].s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_ss");
  aa cctempstr,3,1,&(ac->var[0].ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_s2");
  aa cctempstr,3,1,&(ac->var[0].s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_n");
  aa cctempstr,0,1,&(ac->var[0].n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_s");
  aa cctempstr,3,1,&(ac->var[1].s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_ss");
  aa cctempstr,3,1,&(ac->var[1].ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_s2");
  aa cctempstr,3,1,&(ac->var[1].s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_n");
  aa cctempstr,0,1,&(ac->var[1].n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_vals");
  aa cctempstr,3,AUTOCNEXTARRAYLENGTH,&(ac->vals[0]));
}
#undef aa


/* write values for an instance of struct chainstate_record_updates_and_values */
#ifdef aa
#undef aa
#endif /* aa */
#define aa   awrite (mcffile,
void write_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x)
{
  int i;
  char tempstr[NAMELENGTH];
  char holdstr[NAMELENGTH];
  const char *cctempstr;

  cctempstr = strcat(strcpy(tempstr,x->str),"_num_uptypes");
  aa cctempstr,0,1,&x->num_uptypes);
  for (i=0;i<x->num_uptypes;i++)
  {
    cctempstr = strcat(strcpy(tempstr,x->str),"_");
    cctempstr = strcat(tempstr,x->upnames[i]);
    while (isspace(tempstr[strlen(tempstr)-1]))
      tempstr[strlen(tempstr)-1] = 0;
    strcpy(holdstr,tempstr);
    cctempstr = strcat(tempstr,"_accp");
    aa cctempstr,3,1,&x->upinf[i].accp);
    strcpy(tempstr,holdstr);
    cctempstr = strcat(tempstr,"_tries");
    aa cctempstr,3,1,&x->upinf[i].tries);
  }
  if (x->v != NULL)
  {
    cctempstr = strcat(strcpy(tempstr,x->str),"_do_trend");
    aa cctempstr,0,1,&x->v->do_trend);

    if (x->v->do_trend)
    {
      cctempstr = strcat(strcpy(tempstr,x->str),"_trend");
      aa cctempstr,3,TRENDDIM,x->v->trend);
    }
    cctempstr = strcat(strcpy(tempstr,x->str),"_do_autoc");
    aa cctempstr,0,1,&x->v->do_autoc);
    if (x->v->do_autoc)
    {
      for (i=0;i<AUTOCTERMS;i++)
      {
        write_autoc_record(mcffile,x->str,&(x->v->ac[i]));
      }
    }
    cctempstr = strcat(strcpy(tempstr,x->str),"_beforemin");
    aa cctempstr,0,1,&(x->v->beforemin));
    cctempstr = strcat(strcpy(tempstr,x->str),"_aftermax");
    aa cctempstr,0,1,&(x->v->aftermax));
  }
} //write_struct_chainstate_record_updates_and_values
#undef aa

/*
  read_autoc_record and read_struct_chainstate_record_updates_and_values() are pretty much copies of
  write_autoc_record() and   write_struct_chainstate_record_updates_and_values()
*/
/* read values for an instance of struct autoc*/
#ifdef aa
#undef aa
#endif /* aa */
#define aa  aread(mcffile,
void read_autoc_record(FILE *mcffile,char *vstr,struct autoc *ac)
{
  char tempstr[NAMELENGTH];
  const char *cctempstr;
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_s");
  aa cctempstr,3,1,&(ac->cov.s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_ss");
  aa cctempstr,3,1,&(ac->cov.ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_s2");
  aa cctempstr,3,1,&(ac->cov.s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_cov_n");
  aa cctempstr,0,1,&(ac->cov.n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_s");
  aa cctempstr,3,1,&(ac->var[0].s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_ss");
  aa cctempstr,3,1,&(ac->var[0].ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_s2");
  aa cctempstr,3,1,&(ac->var[0].s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[0]_n");
  aa cctempstr,0,1,&(ac->var[0].n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_s");
  aa cctempstr,3,1,&(ac->var[1].s));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_ss");
  aa cctempstr,3,1,&(ac->var[1].ss));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_s2");
  aa cctempstr,3,1,&(ac->var[1].s2));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_var[1]_n");
  aa cctempstr,0,1,&(ac->var[1].n));
  cctempstr = strcat(strcpy(tempstr,vstr),"_autoc_vals");
  aa cctempstr,3,AUTOCNEXTARRAYLENGTH,&(ac->vals[0]));
}
#undef aa

/* read values for an instance of struct chainstate_record_updates_and_values */
#ifdef aa
#undef aa
#endif /* aa */
#define aa  aread(mcffile,
void read_struct_chainstate_record_updates_and_values(FILE *mcffile,struct chainstate_record_updates_and_values *x)
{
  int i;
  char tempstr[NAMELENGTH];
  char holdstr[NAMELENGTH];
  const char *cctempstr;

  cctempstr = strcat(strcpy(tempstr,x->str),"_num_uptypes");
  aa cctempstr,0,1,&x->num_uptypes);
  for (i=0;i<x->num_uptypes;i++)
  {
    cctempstr = strcat(strcpy(tempstr,x->str),"_");
    cctempstr = strcat(tempstr,x->upnames[i]);
    while (isspace(tempstr[strlen(tempstr)-1]))
      tempstr[strlen(tempstr)-1] = 0;
    strcpy(holdstr,tempstr);
    cctempstr = strcat(tempstr,"_accp");
    aa cctempstr,3,1,&x->upinf[i].accp);
    strcpy(tempstr,holdstr);
    cctempstr = strcat(tempstr,"_tries");
    aa cctempstr,3,1,&x->upinf[i].tries);
  }
  if (x->v != NULL)
  {
    cctempstr = strcat(strcpy(tempstr,x->str),"_do_trend");
    aa cctempstr,0,1,&x->v->do_trend);
    if (x->v->do_trend)
    {
      cctempstr = strcat(strcpy(tempstr,x->str),"_trend");
      aa cctempstr,3,TRENDDIM,x->v->trend);
    }
    cctempstr = strcat(strcpy(tempstr,x->str),"_do_autoc");
    aa cctempstr,0,1,&x->v->do_autoc);
    if (x->v->do_autoc)
    {
      for (i=0;i<AUTOCTERMS;i++)
      {
        read_autoc_record(mcffile,x->str,&(x->v->ac[i]));
      }
    }
    cctempstr = strcat(strcpy(tempstr,x->str),"_beforemin");
    aa cctempstr,0,1,&(x->v->beforemin));
    cctempstr = strcat(strcpy(tempstr,x->str),"_aftermax");
    aa cctempstr,0,1,&(x->v->aftermax));
  }
} //read_struct_chainstate_record_updates_and_values
#undef aa

/***********GLOBAL FUNCTIONS **********/


#ifdef TURNONCHECKS //  readima2mcf is for debugging when we want to read in genealogies generated by IMa2

#ifdef aa
#undef aa
#endif /* aa */
#define aa  aread(mcffile,

/* copy of ima2 readmcf
  this is for debugging when we want to read in genealogies generated by IMa2
*/
/*readmcf() is very similar to writemcf(), except it includes some mallocs() and a small number of initializations at the end*/

/* this reads the mcffile
if the number of chains in the current run (call this cc for now) is different than the number of chains in the mcffile (call this cf) then:
- if cc <= cf,  then the first cc chains in the mcffile are loaded
- if cc > cf,  then the first cf chains in the mcffile are loaded into the first cf positions in C[]
	then the mcffile is closed and reopened and the chains are reloaded in to the next available positions
	this keeps getting done until all the cc positions in C[] have been loaded  */
/* need to revise this so that if cc > cf,  then all chains above cf are copied from the chain cf
  this way we are less likely to have a chain that is much better than its beta value,
  at least compared to what happens when we start copying low number chains into highly heated positions */

void readima2mcf (char ima2mcffilename[])
{
  int i, j, li, ci, lastci = -1;
  double uptime;
  FILE *mcffile;
  char checkeofc;
  int largetimeflag;

  if ((mcffile = fopen (ima2mcffilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening mcffile: %s", ima2mcffilename);
  }


  for (ci = 0; ci < numchainspp; ci++)
  {
// tvalues
    largetimeflag = 0;
    for (i = 0; i < numsplittimes; i++)
    {
      aa "tvalue", 3, 1, &(C[ci]->tvals[i]));
      largetimeflag = largetimeflag || (C[ci]->tvals[i] > T[i].pr.max);
      //assert(C[ci]->tvals[i] < T[i].pr.max);
      C[ci]->poptree[C[ci]->droppops[i + 1][0]].time =
        C[ci]->poptree[C[ci]->droppops[i + 1][1]].time = C[ci]->tvals[i];
    }

    //mutation scalars
    for (li = 0; li < nloci; li++)
    {
      if (nloci > 1 || L[li].nlinked > 0)
      {
        for (i = 0; i < L[li].nlinked; i++)
          aa "uvalue", 3, 1, &(C[ci]->G[li].uvals[i]));
      }
      if (L[li].model == HKY)
        aa "kappavalue", 3, 1, &(C[ci]->G[li].kappaval));
      aa "pi[4]", 3, 4, &(C[ci]->G[li].pi[0]));
      C[ci]->G[li].mignum = 0;
      C[ci]->G[li].tlength = 0;
      C[ci]->G[li].length = 0;
      for (i = 0; i < L[li].numlines; i++)
      {
        aa "up[2]", 0, 2, &(C[ci]->G[li].gtree[i].up[0]));
        aa "down", 0, 1, &(C[ci]->G[li].gtree[i].down));
        if (C[ci]->G[li].gtree[i].down == -1)
        {
          C[ci]->G[li].root = i;
        }

        aa "mut", 0, 1, &(C[ci]->G[li].gtree[i].mut));
        aa "pop", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          aa "A[]", 0, L[li].nlinked, &(C[ci]->G[li].gtree[i].A[0]));
          aa "dlikeA[]", 3, L[li].nlinked,
            &(C[ci]->G[li].gtree[i].dlikeA[0]));
        }
        aa "time", 3, 1, &(C[ci]->G[li].gtree[i].time));
        if (i != C[ci]->G[li].root)
        {
          if (i < L[li].numgenes)
            uptime = 0;
          else
            uptime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[i].up[0]].time;

          C[ci]->G[li].tlength += C[ci]->G[li].gtree[i].time - uptime;
          if (uptime < C[ci]->tvals[npops - 1])
            C[ci]->G[li].tlength +=
              DMIN (C[ci]->G[li].gtree[i].time - uptime,
                    C[ci]->tvals[npops - 1]) - uptime;
        }

        //aa "cmm",0,1,&(C[ci]->G[li].gtree[i].cmm));  don't read this, but reset current value if neeed
        j = 0;
        do
        {
          aa "mig[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mig[j].mt));
          if (C[ci]->G[li].gtree[i].mig[j].mt > 0)
          {
            aa "mig[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mig[j].mp));
            C[ci]->G[li].mignum++;
          }
          j++;
        } while (C[ci]->G[li].gtree[i].mig[j - 1].mt > 0);
        C[ci]->G[li].gtree[i].cmm = j-1; // cmm is handled differently in ima3
       /* if (** JH 8/12/2014 COMMENTED OUT POP ASSIGNMENT STUFF ** assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          aa "asn", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        } */

        if (L[li].model == HKY && i >= L[li].numgenes)
        {
          /* don't think need to save scalefactors, as they get recalculated by makefrac
             aa "C[ci]->G[li].gtree[i].scalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].scalefactor[0]));
             aa "C[ci]->G[li].gtree[i].oldscalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].oldscalefactor[0]));  */
          for (j = 0; j < L[li].numsites; j++)
          {
            aa "C[ci]->G[li].gtree[i].hkyi.frac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.frac[j][0]));
            aa "C[ci]->G[li].gtree[i].hkyi.newfrac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.newfrac[j][0]));
          }
        }
      }
      C[ci]->G[li].roottime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time;
      /* must set fpop,  which is not written by ima2 */
      for (i = 0; i < L[li].numlines; i++)
      {
        if (C[ci]->G[li].gtree[i].down != -1)
          C[ci]->G[li].gtree[i].fpop = C[ci]->G[li].gtree[C[ci]->G[li].gtree[i].down].pop;
        else
          C[ci]->G[li].gtree[i].fpop = -1;
      }
    }
    if ((checkeofc = (char) fgetc (mcffile)) == EOF)
    {
      if (ci == lastci)
        IM_err (IMERR_MCFSPLITTIMEPROB,
                " can't load trees, probably because of multiple instances of splittime time conflict with t prior");
      if (ci < numchainspp - 1)
      {
        FCLOSE (mcffile);
        if ((mcffile = fopen (ima2mcffilename, "r")) == NULL)
        {
           IM_err(IMERR_READFILEOPENFAIL,"Error reopening mcffile: %s", ima2mcffilename);
        }
        lastci = ci;
      }
    }
    else
    {
      ungetc (checkeofc, mcffile);
    }
    if (largetimeflag)  // skip that locus
      ci--;
  }
  fclose (mcffile);
  init_p ();
}                               // readima2mcf
#undef aa
#endif



#ifdef aa
#undef aa
#endif /* aa */
#define aa   awrite (mcffile,

void
writemcf (char mcffilename[],char commandline[],int mcmcrecords,int mcmcrecords_old,int genealogysamples_old,int burninsteps_old,int runsteps_old, double hilike,double hiprob,int currentid)
{
  int i, j, li, ci;
  FILE *mcffile;
  struct chainstate_record_updates_and_values tempstruct1;
  struct update_rate_calc tempstruct2;
  strnl tempname;

  if ((mcffile = fopen (mcffilename, "w")) == NULL)
  {
    IM_err(IMERR_CREATEFILEFAIL,"Error creating mcffile: %s", mcffilename);
  }
  aa "commandline",4,(int) strlen(commandline),commandline);
  aa "numprocesses",0,1,&numprocesses);
  for (ci = 0; ci < numchainspp; ci++)
  {
    aa "chainnumber",0,1,&ci);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      aa "poptreestring",4,(int) strlen(C[ci]->chainpoptreestring),&C[ci]->chainpoptreestring);
#ifdef TURNONCHECKS
       //poptreeprint(ci);
#endif //TURNONCHECKS
    }
    //beta
    if (numchainstotal >  1)
    {
      i = findallbetasci(beta[ci]);
      aa "allbetaci",0,1,&i);
      aa "currallbetapos",0,1,&C[ci]->currallbetapos);
      assert(i==C[ci]->currallbetapos); // i should be redundant once started using currallbetapos
    }

    //updatescalars if called for
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      aa "branchslidescaler",3,1,&(C[ci]->branchslideinfo.updatescalarval));
    }
    for (i = 0; i < numsplittimes; i++)
    {
      if (doRYupdate)
        aa "RYwidthscalar", 3, 1, &(C[ci]->RYwidthinfo[i].updatescalarval));
      if (doNWupdate)
        aa "NWwidthscalar", 3, 1, &(C[ci]->NWwidthinfo[i].updatescalarval));
    }
    // tvalues
    for (i = 0; i < numsplittimes; i++)
      aa "tvalue", 3, 1, &(C[ci]->tvals[i]));
    // hyperpriors
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      /* mltorhpriors and mrtolpriors are each a dictionary, must copy to a list of doubles */
      double *temppriors =  static_cast<double *> (malloc (numdistinctpopulationpairs[npops] * sizeof (double)));
      for (i=0;i<numdistinctpopulationpairs[npops];i++)
        temppriors[i] = getvalue(poppairs[i],C[ci]->mltorhpriors);
      aa "mltorhpriors",3,numdistinctpopulationpairs[npops],temppriors);
      for (i=0;i<numdistinctpopulationpairs[npops];i++)
        temppriors[i] = getvalue(poppairs[i],C[ci]->mrtolhpriors);
      aa "mrtolhpriors",3,numdistinctpopulationpairs[npops],temppriors);
      XFREE(temppriors);
      aa "qhpriors",3,numpopsets,C[ci]->qhpriors);
    }
    //mutation scalars
    for (li = 0; li < nloci; li++)
    {
      if (nloci > 1 || L[li].nlinked > 0)
      {
        for (i = 0; i < L[li].nlinked; i++)
          aa "uvalue", 3, 1, &(C[ci]->G[li].uvals[i]));
      }
      if (L[li].model == HKY)
        aa "kappavalue", 3, 1, &(C[ci]->G[li].kappaval));
      aa "pi[4]", 3, 4, &(C[ci]->G[li].pi[0]));
      if (hiddenoptions[HIDDENGENEALOGY]==1)
        aa "hgprob", 3, 1, &(C[ci]->G[li].hgprob));
      if (hiddenoptions[UPDATEMRATEFORHGUPDATE]==1)
        aa "mhg", 3, 1, &(C[ci]->G[li].mhg));
      for (i = 0; i < L[li].numlines; i++)
      {
        aa "up[2]", 0, 2, &(C[ci]->G[li].gtree[i].up[0]));
        aa "down", 0, 1, &(C[ci]->G[li].gtree[i].down));
        aa "mut", 0, 1, &(C[ci]->G[li].gtree[i].mut));
        aa "pop", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        aa "fpop", 0, 1, &(C[ci]->G[li].gtree[i].fpop));
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          aa "pophg", 0, 1, &(C[ci]->G[li].gtree[i].pophg));
          aa "fpophg", 0, 1, &(C[ci]->G[li].gtree[i].fpophg));
        }
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          aa "A[]", 0, L[li].nlinked, &(C[ci]->G[li].gtree[i].A[0]));
          aa "dlikeA[]", 3, L[li].nlinked, &(C[ci]->G[li].gtree[i].dlikeA[0]));
        }
        aa "time", 3, 1, &(C[ci]->G[li].gtree[i].time));
        aa "cmm",0,1,&(C[ci]->G[li].gtree[i].cmm));
        for (j=0;j<C[ci]->G[li].gtree[i].cmm;j++)
        {
          aa "mig[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mig[j].mt));
          aa "mig[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mig[j].mp));
        }

        if (hiddenoptions[HIDDENGENEALOGY])
        {
          aa "cmmhg",0,1,&(C[ci]->G[li].gtree[i].cmmhg));
          for (j=0;j<C[ci]->G[li].gtree[i].cmmhg;j++)
          {
            aa "mighg[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mighg[j].mt));
            aa "mighg[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mighg[j].mp));
          }
        }
        if (L[li].model == HKY && i >= L[li].numgenes)
        {
          /* don't think need to save scalefactors, as they get recalculated by makefrac
             aa "C[ci]->G[li].gtree[i].scalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].scalefactor[0]));
             aa "C[ci]->G[li].gtree[i].oldscalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].oldscalefactor[0]));  */
          for (j = 0; j < L[li].numsites; j++)
          {
            aa "C[ci]->G[li].gtree[i].hkyi.frac[j]", 3, 4, &(C[ci]->G[li].gtree[i].hkyi.frac[j][0]));
            aa "C[ci]->G[li].gtree[i].hkyi.newfrac[j]", 3, 4, &(C[ci]->G[li].gtree[i].hkyi.newfrac[j][0]));
          }
        }
      }
#ifdef TURNONCHECKS
      //gtreeprint(ci,li);
#endif //TURNONCHECKS
    }
  }
  if (hiddenoptions[READOLDMCFFILE]==0) //jh 1_17_2018
  {
    int temp;
    temp = burninsteps + burninsteps_old;
    aa "burninsteps",0,1,&temp);
    temp = runsteps + runsteps_old;
    aa "runsteps",0,1,&temp);
    aa "numpriormcfruns",0,1,&numpriormcfruns);
    long temptotaltime = (long) totaltime;// copy to a long so aa still works
    aa "totaltime",1,1,&temptotaltime);
  }
  // now do the various instances of struct chainstate_record_updates_and_values
  // everything from this point on goes at the end of the file
  // in the case that runoptions[MCFLOADONLYSTATESPACE] == 1 all this stuff from here to the end gets skipped when reading the mcf file
  for (i=0;i<nloci;i++)
  {
    write_struct_chainstate_record_updates_and_values(mcffile,L[i].g_rec);

    for (j=0;j<L[i].nlinked;j++)
    {
      write_struct_chainstate_record_updates_and_values(mcffile,&L[i].u_rec[j]);
      if (L[i].umodel[j] == STEPWISE)
        write_struct_chainstate_record_updates_and_values(mcffile,&L[i].A_rec[j]);
    }
    if (L[i].model == HKY)
      write_struct_chainstate_record_updates_and_values(mcffile,L[i].kappa_rec);
  }
  for (i=0;i<npops-1;i++)
    write_struct_chainstate_record_updates_and_values(mcffile,&T[i]);
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    write_struct_chainstate_record_updates_and_values(mcffile,poptreeuinfo);
  // make temporary instance of struct chainstate_record_updates_and_values to hold lpgpd_v
  strcpy(tempstruct1.str,"PopP");
  strcpy(tempname,"");
  tempstruct1.num_uptypes = 1;
  tempstruct1.upnames = &tempname;
  tempstruct2.accp = tempstruct2.tries = 0;
  tempstruct1.upinf = &tempstruct2;
  tempstruct1.v = lpgpd_v;
  write_struct_chainstate_record_updates_and_values(mcffile,&tempstruct1);

  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    for (i = 0; i <nummigrateparams; i++)
      write_struct_chainstate_record_updates_and_values(mcffile,&mh[i]);
    if (modeloptions[POPTREETOPOLOGYUPDATE])
      write_struct_chainstate_record_updates_and_values(mcffile,mhnit);

    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      for (i = 0; i <numpopsizeparams; i++)
        write_struct_chainstate_record_updates_and_values(mcffile,&qh[i]);
      if (modeloptions[POPTREETOPOLOGYUPDATE])
        write_struct_chainstate_record_updates_and_values(mcffile,qhnit);
    }
  }

  if (calcoptions[CALCMARGINALLIKELIHOOD]==1)
  {
  	aa "thermosum", 3, numchainstotal, thermosum);
  	/*
    aa "thermosum2", 3, numchainstotal, thermosum2);
	  aa "stepstone2_sumA", 3, numchainstotal, stepstone2_sumA);
	  aa "stepstone2_sumB", 3, numchainstotal, stepstone2_sumB);
	  aa "stepstone_Lmax", 3, numchainstotal, stepstone_Lmax);
    */
  }
  if (hiddenoptions[READOLDMCFFILE]) //jh 1_17_2018
  {
    aa "netsteps",0,1,&netsteps);
    aa "mcmcrecords",0,1,&mcmcrecords);
  }
  else
  {
    int temp;
    if (runoptions[DONTSAVEGENEALOGIES] == 1)  // put 0 in for genealogysamples because none will have been saved in the .ti file
      temp = 0;
    else
      temp = genealogysamples + genealogysamples_old;
    aa "genealogysamples",0,1,&temp);
    temp = mcmcrecords + mcmcrecords_old;
    aa "mcmcrecords",0,1,&temp);
  }

  aa "hilike",3,1,&hilike);
  aa "hiprob",3,1,&hiprob);

  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    aa "numpoptopologies",0,1,&numpoptopologies);
    aa "poptopologycounts",0,numpoptopologies,poptopologycounts);
    aa "totaltopolupdates",0,1,&totaltopolupdates);
    aa "chain0topolupdates",0,1,&chain0topolupdates);
    aa "chain0topolswaps",0,1,&chain0topolswaps);

    aa "poptopologysequence.currentlength",0,1,&poptopologysequence.currentlength);
    aa "poptopologysequence.maxlength",0,1,&poptopologysequence.maxlength);
    assert (poptopologysequence.maxlength >  poptopologysequence.currentlength);
    if (currentid == HEADNODE)
    {
      aa "poptopologysequence.vals",0,poptopologysequence.currentlength,poptopologysequence.vals);
      // note.  disvals is actually full of doubles,  but we don't want to take up that much space and the values are integers. so first copy into a list of unsigned ints
      unsigned short *tempa = static_cast<unsigned short *>  (malloc (poptopologysequence.currentlength * sizeof(unsigned short)));
      for (i=0;i< poptopologysequence.currentlength;i++)
        tempa[i] = (unsigned short) poptopologysequence.disvals[i];
      aa "poptopologysequence.disvals",6,poptopologysequence.currentlength,tempa);
      XFREE(tempa);
    }
  }

  fclose (mcffile);
#ifdef TURNONCHECKS
      //gtreeprint(0,0);
#endif //TURNONCHECKS
//prob_hg_given_g(0,0);
}                               /* writemcf */



#ifdef aa
#undef aa
#endif /* aa */
#define aa  aread(mcffile,
/*readmcf() is very similar to writemcf(), except it includes some mallocs() and a small number of initializations at the end*/

/* this reads the mcffile
if the number of chains in the current run (call this cc for now) is different than the number of chains in the mcffile (call this cf) then:
- if cc <= cf,  then the first cc chains in the mcffile are loaded
- if cc > cf,  then the first cf chains in the mcffile are loaded into the first cf positions in C[]
	then the mcffile is closed and reopened and the chains are reloaded in to the next available positions
	this keeps getting done until all the cc positions in C[] have been loaded  */
/* need to revise this so that if cc > cf,  then all chains above cf are copied from the chain cf
  this way we are less likely to have a chain that is much better than its beta value,
  at least compared to what happens when we start copying low number chains into highly heated positions */
void readmcf (char mcffilename[],int *mcmcrecords,double *hilike,double *hiprob,int currentid)
{
  int i, j, li, ci, lastci = -1;
  int dummy;
  double uptime;
  FILE *mcffile;
  char checkeofc;
  int largetimeflag;
  int checknumprocesses;
  char oldcommandline[COMMANDLINESTRINGLENGTHMAX+15];
  struct chainstate_record_updates_and_values tempstruct1;
  struct update_rate_calc tempstruct2;
  strnl tempname;

  if ((mcffile = fopen (mcffilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening mcffile: %s", mcffilename);
  }
  aa "commandline",4,-1,oldcommandline);
  aa "numprocesses",0,1,&checknumprocesses);
  //if (numprocesses != checknumprocesses)
    //IM_err (IMERR_MPI_CPUNUM, "Error - number of CPUs (%d) in mcf files different than current run (%d) ",checknumprocesses, numprocesses);
  for (ci = 0; ci < numchainspp; ci++)
  {
    aa "chainnumber",0,1,&dummy);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
      aa "poptreestring",4,(int) strlen(C[ci]->chainpoptreestring),&C[ci]->chainpoptreestring);
    //beta
    if (numchainstotal >  1)
    {
      aa "allbetaci",0,1,&i);
      beta[ci] = allbetas[i];
      aa "currallbetapos",0,1,&C[ci]->currallbetapos);
      assert(i==C[ci]->currallbetapos);
    }

    //updatescalars if called for
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      aa "branchslidescaler",3,1,&(C[ci]->branchslideinfo.updatescalarval));
    }
    for (i = 0; i < numsplittimes; i++)
    {
      if (doRYupdate)
        aa "RYwidthscalar", 3, 1, &(C[ci]->RYwidthinfo[i].updatescalarval));
      if (doNWupdate)
        aa "NWwidthscalar", 3, 1, &(C[ci]->NWwidthinfo[i].updatescalarval));
    }

// tvalues
    largetimeflag = 0;
    for (i = 0; i < numsplittimes; i++)
    {
      aa "tvalue", 3, 1, &(C[ci]->tvals[i]));
      largetimeflag = largetimeflag || (C[ci]->tvals[i] > T[i].pr.max);
      //assert(C[ci]->tvals[i] < T[i].pr.max);
      C[ci]->poptree[C[ci]->droppops[i + 1][0]].time =
        C[ci]->poptree[C[ci]->droppops[i + 1][1]].time = C[ci]->tvals[i];
    }
    // hyperpriors
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      /* mltorhpriors and mrtolpriors are each a dictionary, must first put priors in a list of doubles */
      double *temppriors =  static_cast<double *> (malloc (numdistinctpopulationpairs[npops] * sizeof (double)));
      aa "mltorhpriors",3,numdistinctpopulationpairs[npops],temppriors);
      struct dictionary_node_kr *temp;
      for (i=0;i<numdistinctpopulationpairs[npops];i++)
        temp = dictionary_install(poppairs[i],temppriors[i],C[ci]->mltorhpriors);
      aa "mrtolhpriors",3,numdistinctpopulationpairs[npops],temppriors);
      for (i=0;i<numdistinctpopulationpairs[npops];i++)
        temp = dictionary_install(poppairs[i],temppriors[i],C[ci]->mrtolhpriors);
      XFREE(temppriors);
      aa "qhpriors",3,numpopsets,C[ci]->qhpriors);
    }
    //mutation scalars
    for (li = 0; li < nloci; li++)
    {
      if (nloci > 1 || L[li].nlinked > 0)
      {
        for (i = 0; i < L[li].nlinked; i++)
          aa "uvalue", 3, 1, &(C[ci]->G[li].uvals[i]));
      }
      if (L[li].model == HKY)
        aa "kappavalue", 3, 1, &(C[ci]->G[li].kappaval));
      aa "pi[4]", 3, 4, &(C[ci]->G[li].pi[0]));
      if (hiddenoptions[HIDDENGENEALOGY])
        aa "hgprob", 3, 1, &(C[ci]->G[li].hgprob));
       if (hiddenoptions[UPDATEMRATEFORHGUPDATE]==1)
        aa "mhg", 3, 1, &(C[ci]->G[li].mhg));
      C[ci]->G[li].mignum = 0;
      C[ci]->G[li].tlength = 0;
      C[ci]->G[li].length = 0;
      for (i = 0; i < L[li].numlines; i++)
      {
        aa "up[2]", 0, 2, &(C[ci]->G[li].gtree[i].up[0]));
        aa "down", 0, 1, &(C[ci]->G[li].gtree[i].down));
        if (C[ci]->G[li].gtree[i].down == -1)
        {
          C[ci]->G[li].root = i;
        }

        aa "mut", 0, 1, &(C[ci]->G[li].gtree[i].mut));
        aa "pop", 0, 1, &(C[ci]->G[li].gtree[i].pop));
        aa "fpop", 0, 1, &(C[ci]->G[li].gtree[i].fpop));
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          aa "pophg", 0, 1, &(C[ci]->G[li].gtree[i].pophg));
          aa "fpophg", 0, 1, &(C[ci]->G[li].gtree[i].fpophg));
        }
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          aa "A[]", 0, L[li].nlinked, &(C[ci]->G[li].gtree[i].A[0]));
          aa "dlikeA[]", 3, L[li].nlinked,
            &(C[ci]->G[li].gtree[i].dlikeA[0]));
        }
        aa "time", 3, 1, &(C[ci]->G[li].gtree[i].time));
        if (i != C[ci]->G[li].root)
        {
          if (i < L[li].numgenes)
            uptime = 0;
          else
            uptime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[i].up[0]].time;

          C[ci]->G[li].tlength += C[ci]->G[li].gtree[i].time - uptime;
          if (uptime < C[ci]->tvals[npops - 1])
            C[ci]->G[li].tlength +=
              DMIN (C[ci]->G[li].gtree[i].time - uptime,
                    C[ci]->tvals[npops - 1]) - uptime;
        }
        aa "cmm",0,1,&(C[ci]->G[li].gtree[i].cmm));
        for (j=0;j<C[ci]->G[li].gtree[i].cmm;j++)
        {
          aa "mig[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mig[j].mt));
          aa "mig[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mig[j].mp));
        }
        C[ci]->G[li].gtree[i].mig[j].mt = -1.0;

        if (hiddenoptions[HIDDENGENEALOGY])
        {
          aa "cmmhg",0,1,&(C[ci]->G[li].gtree[i].cmmhg));
          for (j=0;j<C[ci]->G[li].gtree[i].cmmhg;j++)
          {
            aa "mighg[].mt", 3, 1, &(C[ci]->G[li].gtree[i].mighg[j].mt));
            aa "mighg[].mp", 0, 1, &(C[ci]->G[li].gtree[i].mighg[j].mp));
          }
          C[ci]->G[li].gtree[i].mighg[j].mt = -1.0;
        }
        if (L[li].model == HKY && i >= L[li].numgenes)
        {
          /* don't think need to save scalefactors, as they get recalculated by makefrac
             aa "C[ci]->G[li].gtree[i].scalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].scalefactor[0]));
             aa "C[ci]->G[li].gtree[i].oldscalefactor[0]",3,L[li].numsites,&(C[ci]->G[li].gtree[i].oldscalefactor[0]));  */
          for (j = 0; j < L[li].numsites; j++)
          {
            aa "C[ci]->G[li].gtree[i].hkyi.frac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.frac[j][0]));
            aa "C[ci]->G[li].gtree[i].hkyi.newfrac[j]", 3, 4,
              &(C[ci]->G[li].gtree[i].hkyi.newfrac[j][0]));
          }
        }
      }
      C[ci]->G[li].roottime = C[ci]->G[li].gtree[C[ci]->G[li].gtree[C[ci]->G[li].root].up[0]].time;
#ifdef TURNONCHECKS
      //gtreeprint(ci,li);
#endif //TURNONCHECKS
    }
    // not sure if this section makes sense for latest code 7/20/2016
    if ((checkeofc = (char) fgetc (mcffile)) == EOF)
    {
      if (ci == lastci)
        IM_err (IMERR_MCFSPLITTIMEPROB,
                " can't load trees, probably because of multiple instances of splittime time conflict with t prior");
      if (ci < numchainspp - 1)
      {
        FCLOSE (mcffile);
        if ((mcffile = fopen (mcffilename, "r")) == NULL)
        {
           IM_err(IMERR_READFILEOPENFAIL,"Error reopening mcffile: %s", mcffilename);
        }
        lastci = ci;
      }
    }
    else
    {
      ungetc (checkeofc, mcffile);
    }
    if (largetimeflag)  // skip that locus
      ci--;
  }
  if (hiddenoptions[READOLDMCFFILE]==0) //jh 1_17_2018
  {
    aa "burninsteps",0,1,&burninsteps);
    aa "runsteps",0,1,&runsteps);
    aa "numpriormcfruns",0,1,&numpriormcfruns);
    long temptotaltime;
    aa "totaltime",1,1,&temptotaltime);
    totaltime = (time_t) temptotaltime;// copy from a longo aa still works
  }

  init_p (); // initialize various things
  // now do the various instances of struct chainstate_record_updates_and_values
  if (runoptions[MCFLOADONLYSTATESPACE] == 0 && runoptions[LOADMCSTATE] == 0)  // do not read this part of the file if runoptions[MCFLOADONLYSTATESPACE] == 1 or runoptions[LOADMCSTATE] == 1
  {
    for (i=0;i<nloci;i++)
    {
      read_struct_chainstate_record_updates_and_values(mcffile,L[i].g_rec);

      for (j=0;j<L[i].nlinked;j++)
      {
        read_struct_chainstate_record_updates_and_values(mcffile,&L[i].u_rec[j]);
        if (L[i].umodel[j] == STEPWISE)
          read_struct_chainstate_record_updates_and_values(mcffile,&L[i].A_rec[j]);
      }
      if (L[i].model == HKY)
        read_struct_chainstate_record_updates_and_values(mcffile,L[i].kappa_rec);
    }
    for (i=0;i<npops-1;i++)
      read_struct_chainstate_record_updates_and_values(mcffile,&T[i]);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
      read_struct_chainstate_record_updates_and_values(mcffile,poptreeuinfo);
    // make temporary instance of struct chainstate_record_updates_and_values to hold lpgpd_v
    strcpy(tempstruct1.str,"PopP");
    strcpy(tempname,"");
    tempstruct1.num_uptypes = 1;
    tempstruct1.upnames = &tempname;
    tempstruct2.accp = tempstruct2.tries = 0;
    tempstruct1.upinf = &tempstruct2;
    tempstruct1.v = lpgpd_v;
    read_struct_chainstate_record_updates_and_values(mcffile,&tempstruct1);

    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      for (i = 0; i <nummigrateparams; i++)
        read_struct_chainstate_record_updates_and_values(mcffile,&mh[i]);
      if (modeloptions[POPTREETOPOLOGYUPDATE])
        read_struct_chainstate_record_updates_and_values(mcffile,mhnit);

      if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      {
        for (i = 0; i <numpopsizeparams; i++)
          read_struct_chainstate_record_updates_and_values(mcffile,&qh[i]);
        if (modeloptions[POPTREETOPOLOGYUPDATE])
          read_struct_chainstate_record_updates_and_values(mcffile,qhnit);
      }
    }

    if (calcoptions[CALCMARGINALLIKELIHOOD]==1)
    {
  	  aa "thermosum", 3, numchainstotal, thermosum);
      /*
  	  aa "thermosum2", 3, numchainstotal, thermosum2);
	    aa "stepstone2_sumA", 3, numchainstotal, stepstone2_sumA);
	    aa "stepstone2_sumB", 3, numchainstotal, stepstone2_sumB);
	    aa "stepstone_Lmax", 3, numchainstotal, stepstone_Lmax);
      */
    }
    if (hiddenoptions[READOLDMCFFILE]) //jh 1_17_2018
    {
      aa "netsteps",0,1,&netsteps);
      aa "recordsteps",0,1,&mcmcrecords);  // used to be "recordsteps"
    }
    else
    {
      aa "genealogysamples",0,1,&genealogysamples);
      aa "mcmcrecords",0,1,mcmcrecords);
    }
    aa "hilike",3,1,hilike);
    aa "hiprob",3,1,hiprob);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      aa "numpoptopologies",0,1,&numpoptopologies);
      aa "poptopologycounts",0,numpoptopologies,poptopologycounts);
      aa "totaltopolupdates",0,1,&totaltopolupdates);
      aa "chain0topolupdates",0,1,&chain0topolupdates);
      aa "chain0topolswaps",0,1,&chain0topolswaps);
      aa "poptopologysequence.currentlength",0,1,&poptopologysequence.currentlength);
      aa "poptopologysequence.maxlength",0,1,&poptopologysequence.maxlength);
      assert (poptopologysequence.maxlength >  poptopologysequence.currentlength);
      if (currentid == HEADNODE)
      {
        if (poptopologysequence.maxlength > POPTOPOLOGYSEQUENCELENGTH)
        {
          poptopologysequence.vals =  static_cast<int *>  (realloc (poptopologysequence.vals,(poptopologysequence.maxlength) * sizeof(int)));
          poptopologysequence.disvals =  static_cast<double *>  (realloc (poptopologysequence.disvals,(poptopologysequence.maxlength) * sizeof(double)));
        }
        aa "poptopologysequence.vals",0,poptopologysequence.currentlength,poptopologysequence.vals);
        // note.  disvals is actually full of doubles,  but we don't want to take up that much space and the values are actually integers. so first copy into a list of unsigned ints
        unsigned short *tempa = static_cast<unsigned short *>  (malloc (poptopologysequence.currentlength * sizeof(unsigned short)));
        aa "poptopologysequence.disvals",6,poptopologysequence.currentlength,tempa);
        for (i=0;i< poptopologysequence.currentlength;i++)
          poptopologysequence.disvals[i] = (double) tempa[i];
        XFREE(tempa);
      }
    }
  }
  fclose (mcffile);
#ifdef TURNONCHECKS
  /*for (ci = 0; ci < numchainspp; ci++)
    for (li = 0; li < nloci; li++)
      gtreeprint(ci,li); */
#endif //TURNONCHECKS
}                               // readmcf
