/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"
#include "update_gtree_common.hpp"

/* Most initializations  */

/*********** LOCAL STUFF **********/

extern double pi[MAXLOCI][4];
extern int urri[2 * MAXLOCI][2 * MAXLOCI];      // used mostly in update_mc_params.c
extern double urrlow[2 * MAXLOCI][2 * MAXLOCI], urrhi[2 * MAXLOCI][2 * MAXLOCI];        // used mostly in update_mc_params.c

static double **uvals;
//static char startpoptreestring[POPTREESTRINGLENGTHMAX]; moved to ima.hpp 
//static double geomeanvar; JH  4/14/2017  not used
static int **numsitesIS;        // maxpops * maxloci   temporarily holds the number of polymorphic sites in loci with infinite sites model
static double uval_preliminary_sum;
static int maxpossiblemigrateparams;
static int numpossiblepoppairstrings;
extern int numdistinctpopulationpairs[]; 
extern int hashvalmaxes[]; 
/* prototypes */
static int numvarHKY (int li, int b, int e);
static int checkaresis (int ci, int period, int i, int j);
static void imaAsnInitAssign (int ci);
static void set_nomigrationchecklist ();
static void setuinfo (double summut);
static double setup_uval (void);
void set_x (struct value_record *v, int isint);
void init_value_record (struct value_record *v, int isint);
static void init_g_rec (int li);
static void init_a_rec (int li);
static void init_lpgpd_v (void);
static void init_i_params (int ci);
static void getparamnums(void);
static void setup_iparams (int ci); 

/* 8/26/2011 */
void fillm_times_rec(int j, int i, int pi, int pj);
static void init_migration_counts (void);
//static void fillmrec(int j, int i, int pi, int pj);
//static void init_migration_counts_times (void);

static void init_mutation_scalar_rec (int li);
static void fixmutationratescalars(void);
static void start_setup_L (char infilename[], int *fpstri, char fpstr[], int currentid);
static void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[]);
static void start_setup_C (void);
static void finish_setup_C (int currentid);
static void set_tvalues (int ci);
static void setup_T ();
static void setup_migprior_recording();
static void setup_qprior_recording();
static void finish_setup_L (void);
static void reportparamcounts(int *fpstri, char fpstr[]);
/* changes made to include hidden genealogies */
static void start_setup_poptree(char *ps);
static void finish_setup_poptree();

/* some extern prototypes */
extern void fillplist (int ci);
extern void fillancplist (int ci);
/******** LOCAL FUNCTIONS ************/

int
numvarHKY (int li, int b, int e)
{
  int i, j, tot = 0;
  for (i = 0; i < L[li].numsites; i++)

  {
    j = b + 1;
    while (j < e && L[li].seq[b][i] == L[li].seq[j][i])
      j++;
    if (j < L[li].numgenes)
      tot++;
  }
  return tot;
}

int
checkaresis (int ci, int period, int i, int j)  // returns 1 if populations i and j are both in period, and they are sister populations
{
  int aresis, k;
  aresis = 0;

  if (ISELEMENT (C[ci]->plist[period][i], C[ci]->periodset[period])
      && ISELEMENT (C[ci]->plist[period][j], C[ci]->periodset[period]))
  {
    for (k = period + 1; k < npops; k++)
    {
      if ((C[ci]->droppops[k][0] == C[ci]->plist[period][i]
           && C[ci]->droppops[k][1] == C[ci]->plist[period][j])
          || (C[ci]->droppops[k][0] == C[ci]->plist[period][j]
              && C[ci]->droppops[k][1] == C[ci]->plist[period][i]))
        aresis = 1;
    }
  }  return aresis;
}                               //checkaresis 

void
imaAsnInitAssign (int ci)
{
  struct genealogy *G; 
  struct edge *gtree; 
  int j;
  int k;
  int li;
  int pi;

  G = NULL;
  gtree = NULL;
  for (li = 0; li < nloci; li++)
  {
    G = &(C[ci]->G[li]);
    gtree = G->gtree;
    j = 0;
    k = 0;
    for (pi = 0; pi < npops; pi++)
      {
        while (j < L[li].samppop[pi] + k)
          {
            gtree[j].pop = pi;
            if (hiddenoptions[HIDDENGENEALOGY]==1)
              gtree[j].pophg = pi;
            j++;
          }
        k += L[li].samppop[pi];
      }
  }
  return;
}

/*
  
  length of the p and r arrays were messed up
  no popsize or migration parameter can extend for more than npops-1 periods. 
*/

void init_i_params (void)
{
  int i, ci,n;
  for (ci = 0;ci<numchainspp;ci++)
  {
    C[ci]->itheta = static_cast<i_param *> 
                  (malloc (numpopsizeparams * sizeof (struct i_param)));
    for (i = 0; i < numpopsizeparams; i++)
    {
      C[ci]->itheta[i].xy = static_cast<plotpoint *> (calloc (GRIDSIZE, sizeof (struct plotpoint)));
    /* try just making these big enough, I think just need for each population which should be 2*npops - 1 */
      C[ci]->itheta[i].wp.p = static_cast<int *> (malloc ((npops-1)* sizeof (int)));
                  //(malloc (2*npops* sizeof (int)));
      C[ci]->itheta[i].wp.r = static_cast<int *> (malloc ((npops-1)* sizeof (int)));
                  //(malloc (2*npops * sizeof (int)));
    }
    C[ci]->imig = static_cast<i_param *> 
                  (malloc (nummigrateparams * sizeof (struct i_param))); 
       /* initialize the parts of wp  for imig*/
    for (i = 0; i < nummigrateparams; i++) // allow for 2 in each, hope that's enough, confusing // changed to 3
    {
      if (modeloptions[ONEMIGRATIONPARAMETER]==1) // n must be big enough to hold 
         n = ((npops *npops - 1) * npops) /3;  // sum of total number of migration types across periods (>= # of migration parameters)
      else
        n =  (npops-1); //# periods

      C[ci]->imig[i].xy = static_cast<plotpoint *> 
                  (calloc (GRIDSIZE, sizeof (struct plotpoint)));
      C[ci]->imig[i].wp.p = static_cast<int *> 
                  (malloc (n * sizeof (int)));
      C[ci]->imig[i].wp.r = static_cast<int *> 
                  (malloc (n * sizeof (int)));
      C[ci]->imig[i].wp.c = static_cast<int *> 
                  (malloc (n * sizeof (int)));
    }
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    holdimig = static_cast<i_param *> 
                (malloc (nummigrateparams * sizeof (struct i_param))); 
    init_migration_prior_update();
  }
}

/* initialize the splitting rate, population size and migration instances of struct i_param. 
   This function is complicated with lots of sections because of many user options */

void getparamnums()
{
  int i, j, k, mi, mcheck;
  
  /* set up the population size parameters */
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    numpopsizeparams = (npops * (npops + 1)) / 2;       // every population in every period gets a parameter
  }
  else
  {
    numpopsizeparams = numtreepops;     // every distinct population gets a parameter
  }
  numpopsets = 1 << npops;//# of possible subsets of sampled populations, =2^npops,  this includes null set and full set 
  if (!modeloptions[PARAMETERSBYPERIOD])
    maxpossiblemigrateparams = 2*(npops-1)*(npops-1);
  /*count how many migration parameters are needed, use mi 
    set up a standard sequence of checks that get repeated when building imig[] 
     there is a complex series of checks to ensure that the intended model is being followed
     loop through periods
     check order:
     MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
     PARAMETERSBYPERIOD  - every population size and migration parameter applies for only 1 period
     NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
     SINGLEMIGRATIONBOTHDIRECTIONS
     if ONEMIGRATIONPARAMETER applies  then there is only one migration rate,  whenever there is migration in the model
   */
  mi = 0;
  if (!modeloptions[NOMIGRATION] && npops > 1)
  {
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 1;
      goto outsidefirstloop; // get out of this big condition section
    }
    /* each chain will have the same number of migration parameters, so C[ARBCHAIN] is ok */
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ( k == 0 ||
              ( modeloptions[MIGRATIONBETWEENSAMPLED]==0 &&
               (modeloptions[PARAMETERSBYPERIOD]==1 ||
                ( !ISELEMENT (C[ARBCHAIN]->plist[k][i], C[ARBCHAIN]->periodset[k - 1]) ||
                  !ISELEMENT (C[ARBCHAIN]->plist[k][j], C[ARBCHAIN]->periodset[k - 1])  
                  ) ) ) )
          /* tricky condition.  proceed to set mcheck if:
            where are considering only sampled populations (k==0)  or
            (not considering just sampled populations  AND
            EITHER we are considering each parameter to apply only to one period
                  OR one of the two populations being considered was not in the previous (k-1) period (meaning it is a new ancestral pop) */
          {
          /* if mcheck ends up as 1, we have one more parameter to add
             if it is also the case that migration is not the same in both directions,  then we have two parameters to add
             
             mcheck can get set to 0 if we require two populations to be sisters, and they are not
             or if the prior given for that pair of populations is zero 
             */
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (ARBCHAIN,k, i, j);
            else
            {
              if (mprior > 0.0 || (calcoptions[LOADPRIORSFROMFILE] && mprior < 0.0) || (calcoptions[LOADPRIORSFROMFILE] && mprior_fromfile[C[ARBCHAIN]->plist[k][i]][C[ARBCHAIN]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
              else
                mcheck = 0;
            }
            if (mcheck)
            {
              mi++;
              mi += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]==0;
            }
          }
outsidefirstloop:  ;

  }

    /* now do another loop, much like the previous one, 
      but now start to build the imig[] structures, including str, b, e
      start to build the wp part of imig[] and determine wp.n*/
  nummigrateparams = mi;
  nummigrateparampairs = nummigrateparams / 2;
  
/*
	type     |  # values  |  cumulative total at end
    cc	       numpopsizeparams   numpopsizeparams
	fc	       numpopsizeparams   2*numpopsizeparams
	hcc	       numpopsizeparams   3*numpopsizeparams
	mc         nummigrateparams   3*numpopsizeparams + nummigrateparams
	fm         nummigrateparams   3*numpopsizeparams + 2*nummigrateparams
	qintegrate numpopsizeparams   4*numpopsizeparams + 2*nummigrateparams
	mintegrate nummigrateparams   4*numpopsizeparams + 3*nummigrateparams
	pdg             1        4*numpopsizeparams + 3*nummigrateparams +  1 
	probg           1        4*numpopsizeparams + 3*nummigrateparams +  2
	t          numsplittimes 4*numpopsizeparams + 3*nummigrateparams +  numsplittimes + 2
*/

  /* set values of position markers for gsampinf */
  
  
  gsamp_ccp = 0;
  gsamp_fcp = gsamp_ccp + numpopsizeparams;
  gsamp_hccp = gsamp_fcp + numpopsizeparams;
  gsamp_mcp = gsamp_hccp + numpopsizeparams;
  gsamp_fmp = gsamp_mcp + nummigrateparams;
  gsamp_qip = gsamp_fmp + nummigrateparams;
  gsamp_mip = gsamp_qip + numpopsizeparams;
  gsamp_pdgp = gsamp_mip + nummigrateparams;
  gsamp_probgp =  gsamp_pdgp + 1;
  gsamp_tp = gsamp_probgp + 1;

}  /* getparamnums */


void set_iparam_poptreeterms(int ci)
{
  
  int i, j, k, ii, jj, kk, pj, ni, mi, mcheck;
  char temps[PARAMSTRLEN];
  
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    filldescendantpops(ci);
  // ithetas
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    for (i = 0, ii = 0, k = npops; i <= lastperiodnumber; i++, k--)
    {
      for (j = 0; j < k; j++)
      {
        pj = C[ci]->plist[i][j];
        sprintf (C[ci]->itheta[ii].str, "q%d,%d", i, pj);
        C[ci]->itheta[ii].b = C[ci]->poptree[pj].b;
        C[ci]->itheta[ii].e = C[ci]->poptree[pj].e;
        C[ci]->itheta[ii].wp.n = 1;
        *C[ci]->itheta[ii].wp.p = i;
        *C[ci]->itheta[ii].wp.r = j;
        ii++;
      }
    }
  }
  else 
  {
    for (ii = 0; ii < numpopsizeparams; ii++)
    {
      sprintf (C[ci]->itheta[ii].str, "q%d", ii);
      C[ci]->itheta[ii].b = C[ci]->poptree[ii].b;
      C[ci]->itheta[ii].e = C[ci]->poptree[ii].e;
      C[ci]->itheta[ii].wp.n = 0;
      for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
      {
        for (j = 0; j < k; j++)
        {
          if (ii == C[ci]->plist[i][j])
            C[ci]->itheta[ii].wp.n += 1;
        }
      }
      assert((ii==numpopsizeparams-1 && C[ci]->itheta[ii].wp.n==1) || (C[ci]->itheta[ii].wp.n == C[ci]->poptree[ii].e - C[ci]->poptree[ii].b ));
      ni = 0;
      for (i = 0, k = npops; i <= lastperiodnumber; i++, k--)
        for (j = 0; j < k; j++)
        {
          if (ii == C[ci]->plist[i][j])
          {
            C[ci]->itheta[ii].wp.p[ni] = i; // the period number
            C[ci]->itheta[ii].wp.r[ni] = j; // the position in plist for pop ii in that period (i.e. period i)
            ni++;
          }
        }
      assert (ni == C[ci]->itheta[ii].wp.n);
    }
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      for (i=npops;i<numpopsizeparams;i++)
      {
        C[ci]->itheta[i].pr.max = C[ci]->qhpriors[(int) C[ci]->descendantpops[i]];
      }
    }
  }
  /* set up the migration parameters 
   two loops here
    1. build  most of the imig[] structures, including str, b, e; start to build wp part of imig[]
    2.  finish building wp parts of imig[]
     there is a complex series of checks to ensure that the intended model is being followed
     loop through periods
     check order:
     MIGRATIONBETWEENSAMPLED  - migration only between sampled populations
     PARAMETERSBYPERIOD  - every population size and migration parameter applies for only 1 period
     NOMIGBETWEENNONSISTERS  - set migration to zero between non-sister populations
     SINGLEMIGRATIONBOTHDIRECTIONS
     if ONEMIGRATIONPARAMETER applies  then there is only one migration rate,  whenever there is migration in the model
   */
    /* this loops is similar to that used for counting migration parameters
      but now start to build the imig[] structures, including str, b, e
      start to build the wp part of imig[] and determine wp.n*/
  if (!modeloptions[NOMIGRATION] && npops > 1)
  {
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {

      mi = 0;
      C[ci]->imig[mi].wp.n = 0;
      C[ci]->imig[mi].b = 0;
      C[ci]->imig[mi].e = lastperiodnumber - 1;
      sprintf (C[ci]->imig[mi].str, "m_all");
    }
    else
      mi = -1;
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ((k == 0) ||
              (!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
                ((!ISELEMENT(C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
                 ||  (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
          {
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (ci,k, i, j);
            else
            {
               if (mprior > 0.0 || (calcoptions[LOADPRIORSFROMFILE] && mprior < 0.0) || (calcoptions[LOADPRIORSFROMFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
               else
                 mcheck = 0;
            }
            if (mcheck)
            {
              if (!modeloptions[ONEMIGRATIONPARAMETER])
              {
                mi++;
                C[ci]->imig[mi].b = k;
                C[ci]->imig[mi].e = k;
                C[ci]->imig[mi].wp.n = 1;
                // for m hyperprior stuff
                if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
                {
                  C[ci]->imig[mi].md.from = C[ci]->plist[k][i];
                  C[ci]->imig[mi].md.to = C[ci]->plist[k][j];
                  //dir is 0 if 'from' pops are on the left side of descstr and 'to' pops are on the right, else dir is 1 
                  C[ci]->imig[mi].dir = makepairstring(C[ci]->descendantpops[C[ci]->imig[mi].md.from],C[ci]->descendantpops[C[ci]->imig[mi].md.to],C[ci]->imig[mi].descstr);
                  if (modeloptions[EXPOMIGRATIONPRIOR])
                  {
                    C[ci]->imig[mi].pr.min =  -1.0;
                    C[ci]->imig[mi].pr.max = -1.0;
                    if (C[ci]->imig[mi].dir == 0)
                      C[ci]->imig[mi].pr.expomean = getvalue(C[ci]->imig[mi].descstr,C[ci]->mltorhpriors);
                      //C[ci]->imig[mi].pr.expomean = C[ci]->mltorhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]];
                    else
                      C[ci]->imig[mi].pr.expomean = getvalue(C[ci]->imig[mi].descstr,C[ci]->mrtolhpriors);
                      //C[ci]->imig[mi].pr.expomean = C[ci]->mrtolhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]];
                  }
                  else
                  {
                    C[ci]->imig[mi].pr.min  = 0.0;
                    C[ci]->imig[mi].pr.expomean =  -1.0;
                    if (C[ci]->imig[mi].dir == 0)
                      C[ci]->imig[mi].pr.max = getvalue(C[ci]->imig[mi].descstr,C[ci]->mltorhpriors);
                      //C[ci]->imig[mi].pr.max = C[ci]->mltorhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]]; 
                    else
                      C[ci]->imig[mi].pr.max = getvalue(C[ci]->imig[mi].descstr,C[ci]->mrtolhpriors);
                        //C[ci]->imig[mi].pr.max = C[ci]->mrtolhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]]; 
//printf("step %d chain %d treenum %d dir %d  mi %d descstr %s prior %.4lf\n",step,ci,C[ci]->poptreenum, dir,mi,C[ci]->imig[mi].descstr,C[ci]->imig[mi].pr.max);
                  }
                }
                if (hiddenoptions[HIDDENGENEALOGY]&& hiddenoptions[GSAMPINFOEXTRA])
                {
                  //sprintf (C[ci]->imig[mi].str, "%d>%d", k, C[ci]->plist[k][j], C[ci]->plist[k][i]); // bug,  why is k here?
                  sprintf (C[ci]->imig[mi].str, "%d>%d", C[ci]->plist[k][j], C[ci]->plist[k][i]);
                }
              }
              else
                C[ci]->imig[mi].wp.n++;
              C[ci]->imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
              if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
              {
                sprintf (C[ci]->imig[mi].str, "m%d,", k);
                if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
                  sprintf (temps, "%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                else
                  sprintf (temps, "%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                strcat (C[ci]->imig[mi].str, temps);
              }
              else
              {
                kk = k + 1;
                while (kk < lastperiodnumber
                       && ISELEMENT (C[ci]->plist[k][i],
                                     C[ci]->periodset[kk])
                       && ISELEMENT (C[ci]->plist[k][j],
                                     C[ci]->periodset[kk]))
                {
                  if (!modeloptions[ONEMIGRATIONPARAMETER])
                    C[ci]->imig[mi].e = kk;
                  kk++;
                  C[ci]->imig[mi].wp.n++;
                  C[ci]->imig[mi].wp.n += modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]; // can't also have ONEMIGRATIONPARAMETER
                }
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])  // can't also have ONEMIGRATIONPARAMETER
                    sprintf (C[ci]->imig[mi].str, "m%d<>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                  else
                    sprintf (C[ci]->imig[mi].str, "m%d>%d", C[ci]->plist[k][i], C[ci]->plist[k][j]);
                }
              }
              if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
              {
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  mi++;
                  C[ci]->imig[mi].b = k;
                  C[ci]->imig[mi].e = k;
                  C[ci]->imig[mi].wp.n = 1;
                  // for m hyperprior stuff
                  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
                  {
                    C[ci]->imig[mi].md.from = C[ci]->plist[k][j];
                    C[ci]->imig[mi].md.to = C[ci]->plist[k][i];
                    //dir is 0 if 'from' pops are on the left side of descstr and 'to' pops are on the right, else dir is 1 
                    C[ci]->imig[mi].dir = makepairstring(C[ci]->descendantpops[C[ci]->imig[mi].md.from],C[ci]->descendantpops[C[ci]->imig[mi].md.to],C[ci]->imig[mi].descstr);
                    if (modeloptions[EXPOMIGRATIONPRIOR])
                    {
                      C[ci]->imig[mi].pr.min = -1.0;
                      C[ci]->imig[mi].pr.max = -1.0;
                      if (C[ci]->imig[mi].dir == 0)
                        C[ci]->imig[mi].pr.expomean = getvalue(C[ci]->imig[mi].descstr,C[ci]->mltorhpriors);
                        //C[ci]->imig[mi].pr.expomean = C[ci]->mltorhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]];
                      else
                        C[ci]->imig[mi].pr.expomean = getvalue(C[ci]->imig[mi].descstr,C[ci]->mrtolhpriors);
                        //C[ci]->imig[mi].pr.expomean = C[ci]->mrtolhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]];
                    }
                    else // uniform migration prior
                    {
                      C[ci]->imig[mi].pr.min  = 0.0;
                      C[ci]->imig[mi].pr.expomean  = -1.0;
                      if (C[ci]->imig[mi].dir == 0)
                        C[ci]->imig[mi].pr.max = getvalue(C[ci]->imig[mi].descstr,C[ci]->mltorhpriors);
                        //C[ci]->imig[mi].pr.max = C[ci]->mltorhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]]; 
                      else
                        C[ci]->imig[mi].pr.max = getvalue(C[ci]->imig[mi].descstr,C[ci]->mrtolhpriors);
                        //C[ci]->imig[mi].pr.max = C[ci]->mrtolhpriors[hashedpairpos[pairhash(C[ci]->imig[mi].descstr)]]; 
                    }
                  }
                 }
                else
                  C[ci]->imig[mi].wp.n++;
                if (modeloptions[PARAMETERSBYPERIOD] && !modeloptions[ONEMIGRATIONPARAMETER])
                {
                  sprintf (C[ci]->imig[mi].str, "m%d,%d>%d", k, C[ci]->plist[k][j], C[ci]->plist[k][i]);
                }
                else
                {
                  kk = k + 1;
                  while (kk < lastperiodnumber &&
                         ISELEMENT (C[ci]->plist[k][i],
                                    C[ci]->periodset[kk])
                         && ISELEMENT (C[ci]->plist[k][j],
                                       C[ci]->periodset[kk]))
                  {
                    if (!modeloptions[ONEMIGRATIONPARAMETER])
                      C[ci]->imig[mi].e = kk;
                    kk++;
                    C[ci]->imig[mi].wp.n++;
                  }
                  if (!modeloptions[ONEMIGRATIONPARAMETER])
                    sprintf (C[ci]->imig[mi].str, "m%d>%d", C[ci]->plist[k][j], C[ci]->plist[k][i]);
                }
              }
            }
           }
    /* now do another loop, much like the previous one. this time finish building the wp parts of C[ci]->imig[]
      using wp.n that was determined in the previous loop */
    if (modeloptions[ONEMIGRATIONPARAMETER])
    {
      mi = 0;
      ni=-1;
    }
    else
      mi = -1;
    for (k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++) if ((k == 0) || 
              (!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
                ((!ISELEMENT (C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
                 || (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
          {
            if (modeloptions[NOMIGBETWEENNONSISTERS])
              mcheck = checkaresis (ci,k, i, j);
            else
            {
               if (mprior > 0.0 || (calcoptions[LOADPRIORSFROMFILE]==1 && mprior < 0.0)  || (calcoptions[LOADPRIORSFROMFILE] && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL))
                mcheck = 1;
               else
                 mcheck = 0;
            }
            if (mcheck)
            {
               if (!modeloptions[ONEMIGRATIONPARAMETER])
               {
                mi++;
                ni = 0;
               }
               else
                 ni++;
              C[ci]->imig[mi].wp.p[ni] = k;
              C[ci]->imig[mi].wp.r[ni] = i;
              C[ci]->imig[mi].wp.c[ni] = j;
              assert (ni < C[ci]->imig[mi].wp.n);
              if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])   // can't also have ONEMIGRATIONPARAMETER
              {
                ni++;
                C[ci]->imig[mi].wp.p[ni] = k;
                C[ci]->imig[mi].wp.r[ni] = j;
                C[ci]->imig[mi].wp.c[ni] = i;
                assert (ni < C[ci]->imig[mi].wp.n);
              }
              if (!modeloptions[PARAMETERSBYPERIOD])
              {
                kk = k + 1;
                while (kk < lastperiodnumber &&
                       ISELEMENT (C[ci]->plist[k][i],
                                  C[ci]->periodset[kk])
                       && ISELEMENT (C[ci]->plist[k][j],
                                     C[ci]->periodset[kk]))
                {
                  ii = 0;
                  while (C[ci]->plist[kk][ii] != C[ci]->plist[k][i])
                    ii++;
                  assert (ii < npops - kk);
                  jj = 0;
                  while (C[ci]->plist[kk][jj] != C[ci]->plist[k][j])
                    jj++;
                  assert (jj < npops - kk);
                  ni++;
                  C[ci]->imig[mi].wp.p[ni] = kk;
                  C[ci]->imig[mi].wp.r[ni] = ii;
                  C[ci]->imig[mi].wp.c[ni] = jj;
                  assert (ni < C[ci]->imig[mi].wp.n);
                  if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])   // can't also have ONEMIGRATIONPARAMETER
                  {
                    ni++;
                    C[ci]->imig[mi].wp.p[ni] = kk;
                    C[ci]->imig[mi].wp.r[ni] = jj;
                    C[ci]->imig[mi].wp.c[ni] = ii;
                    assert (ni < C[ci]->imig[mi].wp.n);
                  }
                  kk++;
                }
              }
              if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS]) 
              {
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  mi++;
                  ni = 0;
                }
                else
                  ni++;
                C[ci]->imig[mi].wp.p[ni] = k;
                C[ci]->imig[mi].wp.r[ni] = j;
                C[ci]->imig[mi].wp.c[ni] = i;
                assert (ni < C[ci]->imig[mi].wp.n);
                if (!modeloptions[PARAMETERSBYPERIOD])
                {
                  kk = k + 1;
                  while (kk < lastperiodnumber &&
                         ISELEMENT (C[ci]->plist[k][i],
                                    C[ci]->periodset[kk])
                         && ISELEMENT (C[ci]->plist[k][j],
                                       C[ci]->periodset[kk]))
                  {
                    ii = 0;
                    while (C[ci]->plist[kk][ii] != C[ci]->plist[k][i])
                      ii++;
                    assert (ii < npops - kk);
                    jj = 0;
                    while (C[ci]->plist[kk][jj] != C[ci]->plist[k][j])
                      jj++;
                    assert (jj < npops - kk);
                    ni++;
                    C[ci]->imig[mi].wp.p[ni] = kk;
                    C[ci]->imig[mi].wp.r[ni] = jj;
                    C[ci]->imig[mi].wp.c[ni] = ii;
                    kk++;
                    assert (ni < C[ci]->imig[mi].wp.n);
                  }
                }
              }
            }
          }                           // initialize wp

  } /* set_iparam_poptreeterms */

}

/* set priors and xy plot stuff for population size and migration parameters */
void
setup_iparams (int ci)
{
  int i, j, k, ii, jj, mi, mcheck;


  // numpopsizeparams = numtreepops;  this is now done in getmaramnums()   // every distinct population gets a parameter
  if (calcoptions[LOADPRIORSFROMFILE])
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      C[ci]->itheta[i].pr.max = popsizeprior_fromfile[i];
      if (C[ci]->itheta[i].pr.max <= 0.0)
          IM_err(IMERR_PRIORFILEVALS,"prior for population %d is set less than or equal to zero",i);
      C[ci]->itheta[i].pr.min = 0;
      C[ci]->itheta[i].pr.expomean = -1.0; 
    }
  }
  else
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      C[ci]->itheta[i].pr.max = thetaprior;
      C[ci]->itheta[i].pr.min = 0;
    }
  }
  for (i = 0; i < numpopsizeparams; i++)
  {
    for (j = 0; j < GRIDSIZE; j++)
    {
      C[ci]->itheta[i].xy[j].x =
        C[ci]->itheta[i].pr.min +
        ((j + 0.5) * (C[ci]->itheta[i].pr.max - C[ci]->itheta[i].pr.min)) / GRIDSIZE;
    }
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    init_hyperprior_arrays(ci);
    if (ci==ARBCHAIN)
    {
      numpossiblepoppairstrings = fillmigratepairs();
      assert(numpossiblepoppairstrings==numdistinctpopulationpairs[npops]);
    }
  }
  /* set up the migration parameter priors if a prior file is uses
  this loop mirrors the first loop for mig parameters in set_iparam_poptreeterms()
  but here we just deal with setting the priors 
  */

  if (calcoptions[LOADPRIORSFROMFILE])
  {
    if (!modeloptions[NOMIGRATION] && npops > 1)
    {
      if (modeloptions[ONEMIGRATIONPARAMETER])
      {
        mi = 0;
      }
      else
        mi = -1;
      for (k = 0; k < lastperiodnumber; k++)
        for (i = 0; i < npops - k - 1; i++)
          for (j = i + 1; j < npops - k; j++) if ((k == 0) ||
                (!modeloptions[MIGRATIONBETWEENSAMPLED] && (modeloptions[PARAMETERSBYPERIOD] ||
                  ((!ISELEMENT(C[ci]->plist[k][i], C[ci]->periodset[k - 1]))
                    ||  (!ISELEMENT(C[ci]->plist[k][j], C[ci]->periodset[k - 1]))))))
            {
              if (modeloptions[NOMIGBETWEENNONSISTERS])
                mcheck = checkaresis (ci,k, i, j);
              else
              {
                  if (mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] > MINPARAMVAL)
                  mcheck = 1;
                  else
                    mcheck = 0;
              }
              if (mcheck)
              {
                if (!modeloptions[ONEMIGRATIONPARAMETER])
                {
                  mi++;
                  /*C[ci]->imig[mi].pr.max = mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
                  C[ci]->imig[mi].pr.min = 0.0;
                  if (C[ci]->imig[mi].pr.max <= MINPARAMVAL)
                        IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j); */
                  if (modeloptions[EXPOMIGRATIONPRIOR])
                  {
                    C[ci]->imig[mi].pr.expomean = mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
                    C[ci]->imig[mi].pr.min =  0.0;
                    C[ci]->imig[mi].pr.max = EXPOMIGPLOTSCALE * C[ci]->imig[mi].pr.expomean;

                  }
                  else
                  {
                    C[ci]->imig[mi].pr.max =mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]];
                    C[ci]->imig[mi].pr.min = 0;
                    C[ci]->imig[mi].pr.expomean = -1.0;
                  }
                  if (C[ci]->imig[mi].pr.max <= MINPARAMVAL)
                    IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
                }
                if (!modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
                {
                  if (!modeloptions[ONEMIGRATIONPARAMETER])
                  {
                    mi++;
                  /*  C[ci]->imig[mi].pr.max = mprior_fromfile[C[ci]->plist[k][j]][C[ci]->plist[k][i]];
                    C[ci]->imig[mi].pr.min = 0.0;
                    if (C[ci]->imig[mi].pr.max <= MINPARAMVAL)
                        IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j); */
                    if (modeloptions[EXPOMIGRATIONPRIOR])
                    {
                      C[ci]->imig[mi].pr.expomean = mprior_fromfile[C[ci]->plist[k][j]][C[ci]->plist[k][i]];
                      C[ci]->imig[mi].pr.min =  0.0;
                      C[ci]->imig[mi].pr.max = EXPOMIGPLOTSCALE * C[ci]->imig[mi].pr.expomean;

                    }
                    else
                    {
                      C[ci]->imig[mi].pr.max =mprior_fromfile[C[ci]->plist[k][j]][C[ci]->plist[k][i]];
                      C[ci]->imig[mi].pr.min = 0;
                      C[ci]->imig[mi].pr.expomean = -1.0;
                    }
                    if (C[ci]->imig[mi].pr.max <= MINPARAMVAL)
                      IM_err(IMERR_PRIORFILEVALS,"migration rate set too low from %d to %d",i,j);
                  }
                }
              }
            }
    }
  }
  else
  {
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      {
        for (i = 0; i < numpossiblepoppairstrings; i++)
        {
          if (modeloptions[EXPOMIGRATIONPRIOR])
          {
            //fill the dictionaries,  keys are poppair strings,  values are priors
            struct dictionary_node_kr *temp;
            temp = dictionary_install(poppairs[i],uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_expo_m_mean),C[ci]->mltorhpriors);
            temp = dictionary_install(poppairs[i],uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_expo_m_mean),C[ci]->mrtolhpriors);
            //C[ci]->mltorhpriors[i] = uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_expo_m_mean); // use hyperprior_expo_m_mean as max for now,  even if it is exponential
            //C[ci]->mrtolhpriors[i] = uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_expo_m_mean);
          }
          else
          {
            //fill the dictionaries,  keys are poppair strings,  values are priors
            struct dictionary_node_kr *temp;
            temp = dictionary_install(poppairs[i],uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_uniform_m_max),C[ci]->mltorhpriors);
            temp = dictionary_install(poppairs[i],uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_uniform_m_max),C[ci]->mrtolhpriors);
            //C[ci]->mltorhpriors[i] = uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_uniform_m_max); 
            //C[ci]->mrtolhpriors[i] = uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_uniform_m_max);
          }
        }
        C[ci]->qhpriors[0] = -1.0; // does not get used because indexes are SET values and 0 is NULL set . 
        for (i=1;i<numpopsets;i++)
          C[ci]->qhpriors[i] = uniforminterval(MINPRIORFROMHYPERPRIOR,hyperprior_uniform_q_max);
      }
    }
    else
    {
      for (ii=0;ii< nummigrateparams;ii++)
      {
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          C[ci]->imig[ii].pr.expomean = expo_m_mean;
          C[ci]->imig[ii].pr.min =  0.0;
          C[ci]->imig[ii].pr.max = EXPOMIGPLOTSCALE * expo_m_mean;
        }
        else
        {
          C[ci]->imig[ii].pr.max = DMAX (m_max, MPRIORMIN);
          C[ci]->imig[ii].pr.min = 0;
          C[ci]->imig[ii].pr.expomean = -1.0;
        }
      }
    }
  }
  for (ii=0;ii< nummigrateparams;ii++)
  {
    for (jj = 0; jj < GRIDSIZE; jj++)
    {
      if (modeloptions[EXPOMIGRATIONPRIOR])
      {
        if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]==0)
          C[ci]->imig[ii].xy[jj].x =
            ((jj + 0.5) * (C[ci]->imig[ii].pr.expomean * EXPOMIGPLOTSCALE)) / GRIDSIZE;
        else
          C[ci]->imig[ii].xy[jj].x =
            ((jj + 0.5) * (hyperprior_expo_m_mean * EXPOMIGPLOTSCALE)) / GRIDSIZE;
      }
      else
      {
        if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]==0)
          C[ci]->imig[ii].xy[jj].x = C[ci]->imig[ii].pr.min +
            ((jj + 0.5) * (C[ci]->imig[ii].pr.max - C[ci]->imig[ii].pr.min)) / GRIDSIZE;
        else
          C[ci]->imig[ii].xy[jj].x = C[ci]->imig[ii].pr.min +
            ((jj + 0.5) * (hyperprior_uniform_m_max)) / GRIDSIZE;
      }
    }
    C[ci]->imig[ii].wp.n = 0;
  }
}                               /* setup_iparams */

/* nomigrationchecklist is three arrays (p for period, r for row, and c for column), each of a set length (nomigrationchecklist.n )

nomigrationchecklist is used to check gweight->mc[][][]
gweight->mc[k][i][j] is the number of migration events in the genealogies in period k 
from population C[ci]->plist[k][i] to population C[ci]->plist[k][j]

So the ith element of the arrays in nomigrationchecklist is used to check the value of 
gweight->mc[nomigrationchecklist.p[i]][nomigrationchecklist.r[i]][nomigrationchecklist.c[i]]
(this is done in update_gree_common.c).  If that value if not zero the update gets rejected. 

building nomigrationchecklist requires identifying which pairs of populations, 
C[ci]->plist[k][i] and C[ci]->plist[k][j], in each period k,  should not be
exchanging migrants.

These pairs are those that meet any of the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 and the populations are not both sampled populations
2. modeloptions[NOMIGBETWEENNONSISTERS]==1  and the populations are not sisters
3. they have had their priors set to zero in a priorfile

These correspond to the following:
1. modeloptions[MIGRATIONBETWEENSAMPLED]==1 && C[ci]->plist[k][i] >= npops && C[ci]->plist[k][j] >= npops
2. modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (ci,k, i, j)==1
3. calcoptions[LOADPRIORSFROMFILE]==1 && mprior_fromfile[C[ci]->plist[k][i]][C[ci]->plist[k][j]] <= MINPARAMVAL

To set this up we must go through two loops,  first to cound how many (i.e. set nomigrationchecklist.n)
Then after initializing the arrays in nomigrationchecklist, we must do another similar loop to fill them up. 
*/

/* this should not be needed if modeloptions[NOMIGRATION]==1 */
void
set_nomigrationchecklist ()  // not used for population tree updating so work on C[ARBCHAIN]
{
  int n, i, j, k;

  nomigrationchecklist.n = 0;
  nomigrationchecklist.p = NULL;
  nomigrationchecklist.r = NULL;
  nomigrationchecklist.c = NULL;
  for (k = 0; k < lastperiodnumber; k++)
    for (i = 0; i < npops - k - 1; i++)
      for (j = i + 1; j < npops - k; j++)
      {
        if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ARBCHAIN]->plist[k][i] >= npops || C[ARBCHAIN]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (ARBCHAIN,k, i, j)==0) ||
            (calcoptions[LOADPRIORSFROMFILE]==1 && mprior < 0.0 && mprior_fromfile[C[ARBCHAIN]->plist[k][i]][C[ARBCHAIN]->plist[k][j]] <= MINPARAMVAL))
        {
          nomigrationchecklist.n+=2;
        }
      }
    if (nomigrationchecklist.n > 0)
    {
      nomigrationchecklist.p = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
      nomigrationchecklist.r = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
      nomigrationchecklist.c = static_cast<int *> 
                        (malloc (nomigrationchecklist.n * sizeof (int)));
    }
    for (n = -1, k = 0; k < lastperiodnumber; k++)
      for (i = 0; i < npops - k - 1; i++)
        for (j = i + 1; j < npops - k; j++)
        {
          if ((modeloptions[MIGRATIONBETWEENSAMPLED]==1 && (C[ARBCHAIN]->plist[k][i] >= npops || C[ARBCHAIN]->plist[k][j] >= npops)) ||
            (modeloptions[NOMIGBETWEENNONSISTERS]==1 && checkaresis (ARBCHAIN,k, i, j)==0) ||
            (calcoptions[LOADPRIORSFROMFILE]==1 && mprior < 0.0 && mprior_fromfile[C[ARBCHAIN]->plist[k][i]][C[ARBCHAIN]->plist[k][j]] <= MINPARAMVAL))
          {
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = i;
              nomigrationchecklist.c[n] = j;
              n++;
              nomigrationchecklist.p[n] = k;
              nomigrationchecklist.r[n] = j;
              nomigrationchecklist.c[n] = i;
          }
        }
  assert(nomigrationchecklist.n == n+1);
  return;
}                               // set_nomigrationchecklist

#define MAXPRIORSCALETRY 10000000       // max number of times to try getting starting mutation rates compatible with priors
void
setuinfo (double summut)   // called only for chain 0 because mutation rate info gets copied into other chains
{
  int i, li, ui, uj, k;
  int mutchain0 = 0;
  int numuprior = 0;
  double priorscaletry = 0;
  double U, r;
  int doneu, rcheck, numupair, *upriorpairlist1, *upriorpairlist2;
  double prodcheck, maxr, newr, d;

/*  ul is global.  It is a list of struct upairlist that contains a all of the locations of the mutation rate scalars. 
    each upairlist has a .l value (has the locus #) and .u value which is the scalar number for that locus 
    Each locus L[li] has a list uii that has the position in ul of that scalar

    uii is the position in that list
	in set_mcparam_values() the priors on mutation rate scalars are set to standard pos (max) and neg (min) values (they are on a log scale)
	   and are stored in the mcinf.pr   (e.g. C[mutchain0]->G[0].u[0].mcinf.pr)
	if the input file has prior ranges on mutation yets,  these are read into uperyear.pr  (e.g. C[ci]->G[li].u[ui].uperyear.pr)
	The ratios of these uperyear values among loci must be taken to reset the values of mcinf.pr   
*/
  for (i = 0, li = 0; li < nloci; li++)
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      ul[i].l = li;
      ul[i].u = ui;
      L[li].uii[ui] = i;
      i++;
    }




  if (hiddenoptions[NOMUTATIONSCALARUPATES] == 1) // fix them all to be 1 and then do not update them
  {
    domutationscalarupdate = 0;
    for (li=0;li<nloci;li++) if (L[li].model != INFINITESITES)
      IM_err (IMERR_MUTSCALEPROBLEM,"locus %d is not infinite sites mutation model.  This is not allowed when fixing mutation rate scalars. ",li);
    for (ui = 0; ui < nurates; ui++)
      C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u] =  1.0;
  }
  else if (hiddenoptions[FIXMUTATIONSCALARUPATES] == 1)
  {
    domutationscalarupdate = 0;
    fixmutationratescalars();
  }
  else
  {
    /* reset the mutation rate values to have geometric mean of 1 */ // why only for chain 0 ? - because later they are copied from chain 0 to the other chains 
    for (prodcheck = 1, ui = 0; ui < nurates; ui++)
    {
      C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u] =  exp (uvals[ul[ui].l][ul[ui].u] - summut / nurates);
      prodcheck *= C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u];
    }
    if (fabs (log (prodcheck)) > 1e-5)
      IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
    domutationscalarupdate = 1;

  }
  for (ui = 0; ui < nurates; ui++)
  {
    if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0)
      numuprior++;
  }

  if (calcoptions[MUTATIONPRIORRANGE])
  {
    assert (numuprior > 1);
    numupair = numuprior * (numuprior - 1) / 2;
    upriorpairlist1 = static_cast<int *> 
                    (calloc ((size_t) numupair, sizeof (int)));
    upriorpairlist2 = static_cast<int *> 
                    (calloc ((size_t) numupair, sizeof (int)));
    k = 0;
    for (ui = 0; ui < nurates - 1; ui++)
    {
      for (uj = ui + 1; uj < nurates; uj++)
      {
        if (L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)

        {
          upriorpairlist1[k] = ui;
          upriorpairlist2[k] = uj;
          k++;
        }
      }
    }

    /* urri[i][j] has a 0 if neither i nor j has a prior. 1 if i has a prior and j does not,  -1 if i does not have a pior and j does  */
    for (ui = 0; ui < nurates; ui++)
    {
      for (uj = 0; uj < nurates; uj++)
      {
        urri[ui][uj] = 0;
        if (ui != uj
            && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
            && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0)
        {
          urrlow[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].min /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].max);
          urrhi[ui][uj] =
            log (L[ul[ui].l].uperyear_prior[ul[ui].u].max /
                 L[ul[uj].l].uperyear_prior[ul[uj].u].min);
          urri[ui][uj] = urri[uj][ui] = 2;
        }
        else
        {
          if (ui != uj
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min != 0
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min == 0)
            urri[ui][uj] = 1;
          if (ui != uj
              && L[ul[uj].l].uperyear_prior[ul[uj].u].min != 0
              && L[ul[ui].l].uperyear_prior[ul[ui].u].min == 0)
            urri[ui][uj] = -1;
        }
      }
    }
    /* need to set all of the uscalers so that their ratios are in the ranges defined by the priors */
    /* it is possible that suiteable sets of scalars will not be able to be found */
    maxr = 3 * L[0].u_rec[0].pr.max;

    do
    {
      doneu = 1;
      for (i = 0; i < numupair; i++)

      {
        ui = upriorpairlist1[i];
        uj = upriorpairlist2[i];
        r =log (C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u] / C[mutchain0]->G[ul[uj].l].uvals[ul[uj].u]);
        rcheck = (r >= urrlow[ui][uj] && r <= urrhi[ui][uj]);
        doneu = doneu && rcheck;
        if (!rcheck)
        {
          do
          {
            U = uniform ();
            if (U > 0.5)
              newr = r + maxr * (2.0 * U - 1.0);

            else
              newr = r - maxr * U * 2.0;
            if (newr > maxr)
              newr = 2.0 * maxr - newr;

            else if (newr < -maxr)
              newr = 2.0 * (-maxr) - newr;
            d = exp ((newr - r) / 2);
          } while ((newr <= urrlow[ui][uj] || newr >= urrhi[ui][uj]));
          C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u] *= d;
          C[mutchain0]->G[ul[uj].l].uvals[ul[uj].u] /= d;
        }
      }
      priorscaletry++;
    } while (!doneu && priorscaletry < MAXPRIORSCALETRY);

    if (priorscaletry >= MAXPRIORSCALETRY)
      IM_err(IMERR_MUTSCALARPRIORRANGEFAIL,"More than %d failed attempts at finding a valid set of mutation scalars, prior ranges probably too restrictive",MAXPRIORSCALETRY);
    for (prodcheck = 1, ui = 0; ui < nurates; ui++)
      prodcheck *= C[mutchain0]->G[ul[ui].l].uvals[ul[ui].u];

    if (fabs (log (prodcheck)) > 1e-5)
    {
      IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
    }
    XFREE (upriorpairlist1);
    XFREE (upriorpairlist2);
  }
}                               /* setuinfo */

double
setup_uval (void)
  /* set mutation rate scalars */
  /* get a relative mutation rate for a locus by summing up estimates of 4Nu for all the populations
     the geometric mean of these relative values are then used to set the starting mutation rates */
{
  int j, k, li, ui;
  int b, e;
  double w, w2;
  double temp, uval_temp_sum, sumq;
  int npopstemp;
  

  npopstemp = npops;
  uval_temp_sum = 0;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].model == INFINITESITES
        || L[li].model == JOINT_IS_SW 
        || L[li].model == HKY)
    {
      sumq = 0;
      b = 0;
      e = L[li].samppop[0] - 1;
      for (k = 0; k < npopstemp; k++)
      {
        if (e >= b)
        {
          for (j = 1, w = 0.0; j <= e; j++)
            w += 1 / (double) j;
          if (w == 0)
            w = 1;
          if (L[li].model == HKY)
            temp = numvarHKY (li, b, e) / w;
          else
            temp = numsitesIS[li][k] / w;
        }
        else
        {
          temp = 0;
        }
        b = e + 1;
        if (k + 1 == npopstemp)
        {
          e = e + L[li].numgenesunknown;
        }
        else
        {
          e = e + L[li].samppop[k + 1];
        }
        sumq += temp;

      }
      sumq = DMAX (0.1, sumq);  // nominal low values in case of zero variation 
      uvals[li][0] = log (sumq);
      uval_temp_sum += uvals[li][0];
    }

    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      somestepwise = 1;
      if (L[li].model == JOINT_IS_SW)
        ui = 1;
      else
        ui = 0;

      for (; ui < L[li].nlinked; ui++)
      {
        b = 0;
        e = L[li].samppop[0] - 1;
        sumq = 0;
        for (k = 0; k < npopstemp; k++)
        {
          if (e >= b)
          {
            for (j = b, w = 0, w2 = 0; j <= e; j++)
            {
              w += L[li].A[ui][j];
              w2 += SQR ((double) L[li].A[ui][j]);
            }
            if (k == npops)
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].numgenesunknown) /
                (L[li].numgenesunknown);
            }
            else
            {
              temp =
                (double) 2 *(w2 -
                             SQR (w) / L[li].samppop[k]) / (L[li].samppop[k]);
            }
            sumq += temp;
          }
          else
          {
            temp = 0;
          }
          b = e + 1;
          if (k + 1 == npopstemp)
          {
            e = e + L[li].numgenesunknown;
          }
          else
          {
            e = e + L[li].samppop[k + 1];
          }
        }
        sumq = DMAX (0.1, sumq);        // nominal low values in case of zero variation 
        uvals[li][ui] = log (sumq);
        uval_temp_sum += uvals[li][ui];
      }
    }
  }

  // geomeanvar is used to set the prior on thetas, in setup_iparams(), depending on command line options 
  /* sum is the sum of the logs, across loci, of the maximal amount of variation found among the sampled populations */
  /* the scalar of thetaprior is multiplied times the geometric mean of the maximal estimates of variation across loci */
  // geomeanvar = exp (uval_temp_sum / nurates);  JH 4/14/2017  noticed that this was not used at all  only uval_temp_sum (the sum ov uvals) gets used later 
  return uval_temp_sum;
}                               /*setup_uval */

void
set_x (struct value_record *v, int isint)
{
  int i;
  for (i = 0; i < GRIDSIZE; i++)
  {
    if (v->do_logplot)
    {
      v->xy[i].x = exp (v->plotrange.min + (i + 0.5) * v->plotrescale); // the same scale is used for all mutation rate scalars 
    }
    else
    {
      if (isint)
        v->xy[i].x = (double) i;
      else
        v->xy[i].x = v->plotrange.min + ((i + 0.5) * (v->plotrange.max * v->plotrescale - v->plotrange.min)) / GRIDSIZE;
    }
  }
}                               // set_x 

void
init_value_record (struct value_record *v, int isint)
{
  if (v->do_xyplot)
  {
    v->xy = static_cast<plotpoint *> 
                (calloc (GRIDSIZE, sizeof (struct plotpoint)));
    set_x (v, isint);
  }

  if (v->do_trend)
    v->trend = static_cast<double *> (calloc (TRENDDIM, sizeof (double)));

  v->beforemin = v->aftermax = 0;
}                               // init_value_Record

void
init_g_rec (int li)
{
  int num_g_update_types = IM_UPDATE_GENEALOGY_NUMBER;
  L[li].g_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (sizeof (struct chainstate_record_updates_and_values)));
  sprintf (L[li].g_rec->str, "gtree_%d", li);
  L[li].g_rec->num_uptypes = num_g_update_types;
  L[li].g_rec->upnames = static_cast<strnl *> 
                (malloc (L[li].g_rec->num_uptypes * sizeof (strnl)));
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_ANY],      "branch     ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TOPOLOGY], "topology   ");
  sprintf (L[li].g_rec->upnames[IM_UPDATE_GENEALOGY_TMRCA],    "tmrca      ");
  L[li].g_rec->upinf = static_cast<update_rate_calc *> 
        (calloc ((size_t) L[li].g_rec->num_uptypes, 
        sizeof (struct update_rate_calc)));
  L[li].g_rec->num_vals = 1;
  L[li].g_rec->v = static_cast<value_record *> 
            (malloc (L[li].g_rec->num_vals * sizeof (struct value_record)));
  sprintf (L[li].g_rec->v->str, "g%d_tmrca", li);
  sprintf (L[li].g_rec->v->strshort, "tmrca%d", li);
  L[li].g_rec->v->plotrange.min = 0;
  L[li].g_rec->v->do_autoc = 1;
  L[li].g_rec->v->do_xyplot = outputoptions[PRINTTMRCA];
  L[li].g_rec->v->do_trend = 1;
  L[li].g_rec->v->plotrescale = 1.0;
  L[li].g_rec->v->do_logplot = 0;
  

}                               // init_g_rec

void
init_lpgpd_v (void)
{
  lpgpd_v = static_cast<value_record *> (malloc (sizeof (struct value_record)));
  sprintf (lpgpd_v->str, "Log[P(G)+P(D|G)]");
  sprintf (lpgpd_v->strshort, "Log[P]");
  lpgpd_v->plotrange.min = lpgpd_v->plotrange.max = 0;
  lpgpd_v->do_xyplot = 0;
  lpgpd_v->do_logplot = 0;
  lpgpd_v->do_trend = 1;
  lpgpd_v->do_autoc = 1;
  lpgpd_v->plotrescale = 1.0;
  init_value_record (lpgpd_v, 0);
  return;
}                               // init_lpgpd_rec
/* 8/26/2011 */
void fillm_times_rec(int j, int i, int pi, int pj)
{
  char stemp[6];
  sprintf(stemp,"m%d>%d",pi,pj);
  strcpy (migration_counts[j][i].str,stemp);
  migration_counts[j][i].plotrange.max = (double) GRIDSIZE - 1;       // used for counts so the grid position is the count #
  migration_counts[j][i].plotrange.min = 0;
}

/*void fillmrec(int j, int i, int pi, int pj)
{
  char stemp[6];
  sprintf(stemp,"%d>%d",pi,pj);
  strcpy (migration_counts_times[j][i].str, stemp);
  strcat (migration_counts_times[j][i].str, "_t");
  migration_counts_times[j][i].plotrange = T[numsplittimes - 1].pr; // same range as for splitting times
  strcpy (migration_counts_times[j][i+1].str,stemp);
  strcat (migration_counts_times[j][i+1].str, "_#");
  migration_counts_times[j][i+1].plotrange.max = (double) GRIDSIZE - 1;       // used for counts so the grid position is the count #
  migration_counts_times[j][i+1].plotrange.min = 0;
}*/

void
init_migration_counts (void)  
{
  int i, j, k, pi, pj, ji, jj, mrows;//, nummigdirs;
  mrows = nloci + (nloci > 1); 
  nummigdirs = 2*(npops-1)*(npops-1);
  migration_counts = static_cast<value_record **> (malloc (mrows * sizeof (struct value_record *)));
  for (j = 0; j < mrows; j++)
    migration_counts[j] = static_cast<value_record *>
      (malloc (nummigdirs * sizeof (struct value_record)));
  for (j = 0; j < mrows; j++)
  {
    i=0;
    for (k= 0; k < lastperiodnumber; k++)
    {
      for (ji=0; ji< npops-k;ji++) 
      {
          pi = C[ARBCHAIN]->plist[k][ji];
          if (k==0)
          {
            for (jj=ji+1; jj< npops-k;jj++) 
            {
              pj = C[ARBCHAIN]->plist[k][jj];
              if (pi != pj)
              {
                fillm_times_rec(j,i,pi,pj);
                fillm_times_rec(j,i+1,pj,pi);
                i+=2;
}
            }
          }
          else
          {
            pj = C[ARBCHAIN]->addpop[k];
            if (pi != pj)
            {
              fillm_times_rec(j,i,pi,pj);
              fillm_times_rec(j,i+1,pj,pi);
              i+=2;
            }
          }
      }
    }
    assert(i==nummigdirs);
  }
  for (j = 0; j < mrows; j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migration_counts[j][i].do_xyplot = 1;
      migration_counts[j][i].do_logplot = 0;
      migration_counts[j][i].do_trend = 0;
      migration_counts[j][i].do_autoc = 0;
      migration_counts[j][i].plotrescale = 1.0;
      init_value_record (&migration_counts[j][i], 1);
    }
}                               //init_migration_counts


void
init_migration_counts_hold (void)   // only for cold chain
{
  int i, j, k, pi, pj, ji, jj, mrows;//, nummigdirs;
  int z = whichiscoldchain();
  mrows = nloci + (nloci > 1); 
  nummigdirs = 2*(npops-1)*(npops-1);
  /*if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
    nummigdirs = 2*nummigrateparams;  
  else
    nummigdirs = nummigrateparams;   */ 
  migration_counts = static_cast<value_record **> (malloc (mrows * sizeof (struct value_record *)));
  for (j = 0; j < mrows; j++)
    migration_counts[j] = static_cast<value_record *>
      (malloc (nummigdirs * sizeof (struct value_record)));
  for (j = 0; j < mrows; j++)
  {
    i=0;
    for (k= 0; k < lastperiodnumber; k++)
    {
      for (ji=0; ji< npops-k;ji++) 
      {
          pi = C[z]->plist[k][ji];
          if (k==0)
          {
            for (jj=ji+1; jj< npops-k;jj++) 
            {
              pj = C[z]->plist[k][jj];
              if (pi != pj)
              {
                fillm_times_rec(j,i,pi,pj);
                fillm_times_rec(j,i+1,pj,pi);
                i+=2;
}
            }
          }
          else
          {
            pj = C[z]->addpop[k];
            if (pi != pj)
            {
              fillm_times_rec(j,i,pi,pj);
              fillm_times_rec(j,i+1,pj,pi);
              i+=2;
            }
          }
      }
    }
    assert(i==nummigdirs);
  }
  for (j = 0; j < mrows; j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migration_counts[j][i].do_xyplot = 1;
      migration_counts[j][i].do_logplot = 0;
      migration_counts[j][i].do_trend = 0;
      migration_counts[j][i].do_autoc = 0;
      migration_counts[j][i].plotrescale = 1.0;
      init_value_record (&migration_counts[j][i], 1);
    }
}                               //init_migration_counts_hold

void
init_mutation_scalar_rec (int li)
{
  int i, ui, ai;
  int num_u_update_types = 1;   // number of different types of updates for mutation rate scalars and kappa values
  int num_A_update_types = 1;   // number of different types of updates for A states 
  double uscale;

  // initialize u_rec
  L[li].u_rec = static_cast<chainstate_record_updates_and_values *> 
            (malloc (L[li].nlinked *
                     sizeof (struct chainstate_record_updates_and_values)));
  // set priors and windows
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].pr.max = log (UMAX);        // 10000
    L[li].u_rec[ui].pr.min = -log (UMAX);
    L[li].u_rec[ui].win = L[li].u_rec[ui].pr.max / nloci;
  }
  uscale = (2 * log (UMAX)) / GRIDSIZE; // 10000
  // name the mutation rate scalars     
  if (L[li].model == STEPWISE)
  {
    for (ai = 0; ai < L[li].nlinked; ai++)
      sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  else
  {
    sprintf (L[li].u_rec[0].str, "%du", li);
    if (L[li].model == JOINT_IS_SW)
      for (ai = 1; ai < L[li].nlinked; ai++)
        sprintf (L[li].u_rec[ai].str, "%dSW%d", li, ai);
  }
  for (ui = 0; ui < L[li].nlinked; ui++)
  {
    L[li].u_rec[ui].num_uptypes = num_u_update_types;
    L[li].u_rec[ui].upnames = static_cast<strnl *> 
                    (malloc (L[li].u_rec[ui].num_uptypes * sizeof (strnl)));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].u_rec[ui].upnames[i], "scalar_update");
    }
    L[li].u_rec[ui].upinf = static_cast<update_rate_calc *> 
            (calloc ((size_t) L[li].u_rec[ui].num_uptypes, 
                     sizeof (struct update_rate_calc)));
    L[li].u_rec[ui].num_vals = 1;
    L[li].u_rec[ui].v = static_cast<value_record *> 
            (malloc (L[li].u_rec[ui].num_vals * sizeof (struct value_record)));
    for (i = 0; i < L[li].u_rec[ui].num_vals; i++)
    {
      strcpy (L[li].u_rec[ui].v[i].str, L[li].u_rec[ui].str);   // will this ever get used ? 
      L[li].u_rec[ui].v[i].do_xyplot = 1;
      L[li].u_rec[ui].v[i].plotrescale = uscale;
      L[li].u_rec[ui].v[i].do_logplot = 1;
      L[li].u_rec[ui].v[i].do_autoc = 0;
      L[li].u_rec[ui].v[i].do_trend = 1;
      L[li].u_rec[ui].v[i].plotrange = L[li].u_rec[ui].pr;
      init_value_record (&L[li].u_rec[ui].v[i], 0);
    }
  }

  // do kappa_rec 
  if (L[li].model == HKY)
  {
    L[li].kappa_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (sizeof (struct chainstate_record_updates_and_values)));
    L[li].kappa_rec->pr.max = KAPPAMAX;
    L[li].kappa_rec->pr.min = 0.0;
    L[li].kappa_rec->win = 2.0;
    sprintf (L[li].kappa_rec->str, "%d_Ka", li);
    L[li].kappa_rec->num_uptypes = num_u_update_types;
    L[li].kappa_rec->upnames = static_cast<strnl *> 
                            (malloc (num_u_update_types * sizeof (strnl)));
    for (i = 0; i < num_u_update_types; i++)
    {
      sprintf (L[li].kappa_rec->upnames[i], "kappa_update");
    }
    L[li].kappa_rec->upinf = static_cast<update_rate_calc *> 
                    (calloc ((size_t) num_u_update_types, 
                             sizeof (struct update_rate_calc)));
    L[li].kappa_rec->num_vals = 1;
    L[li].kappa_rec->v = static_cast<value_record *> 
            (malloc (L[li].kappa_rec->num_vals * sizeof (struct value_record)));
    strcpy (L[li].kappa_rec->v->str, L[li].kappa_rec->str);
    strncpy (L[li].kappa_rec->v->strshort, L[li].kappa_rec->str,
             PARAMSTRLENSHORT - 1);
    L[li].kappa_rec->v->do_xyplot = 1;
    L[li].kappa_rec->v->do_logplot = 0;
    L[li].kappa_rec->v->do_trend = 1;
    L[li].kappa_rec->v->plotrange = L[li].kappa_rec->pr;
    L[li].kappa_rec->v->plotrescale = 1.0;
    L[li].kappa_rec->v->do_autoc = 0;

    init_value_record (L[li].kappa_rec->v, 0);
  }

  // do A_rec 
  // A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
  {
    L[li].A_rec = static_cast<chainstate_record_updates_and_values *> 
                (malloc (L[li].nlinked * 
                         sizeof (struct chainstate_record_updates_and_values)));
    for (ai = 0; ai < L[li].nlinked; ai++)
    {
      if (L[li].umodel[ai] == STEPWISE)
        sprintf (L[li].A_rec[ai].str, "%dSW%d", li, ai);
      else
        sprintf (L[li].A_rec[ai].str, "NOTSTEPWISE");
      L[li].A_rec[ai].num_uptypes = num_A_update_types;
      L[li].A_rec[ai].upnames = static_cast<strnl *> 
                            (malloc (num_A_update_types * sizeof (strnl)));
      for (i = 0; i < num_A_update_types; i++)
      {
        sprintf (L[li].A_rec[ai].upnames[i], "STR_update");
      }
      L[li].A_rec[ai].upinf = static_cast<update_rate_calc *> 
                        (calloc ((size_t) num_A_update_types, 
                         sizeof (struct update_rate_calc)));
      L[li].A_rec[ai].num_vals = 0;
      L[li].A_rec[ai].v  = NULL;
    }
  }
}                               // init_mutation_scalar_rec

void fixmutationratescalars()
{
  int li,ui;
  double gm,prodcheck,lprod = 0.0; 
  for (li=0;li<nloci;li++)
  {
    if (L[li].model != INFINITESITES)
      IM_err (IMERR_MUTSCALEPROBLEM,"locus %d is not infinite sites mutation model.  This is not allowed when fixing mutation rate scalars. ",li);
    if (L[li].nlinked != 1)
      IM_err (IMERR_MUTSCALEPROBLEM,"locus %d has multiple components, each with a scalar.  This is not allowed when fixing mutation rate scalars. ",li);
    if (L[li].uperyear_vals[0] <= 0.0)
      IM_err (IMERR_MUTSCALEPROBLEM,"locus %d has a missing or negative mutation rate per year. This is not allowed when fixing mutation rate scalars. ",li);
    lprod += log(L[li].uperyear_vals[0]);
  }
  gm = exp(lprod/nloci);

  for (prodcheck = 1.0,ui = 0; ui < nurates; ui++)
  {
    C[0]->G[ul[ui].l].uvals[ul[ui].u] = L[ul[ui].l].uperyear_vals[ul[ui].u]/gm ; // why only for chain 0 ? - because later they are copied from chain 0 to the other chains 
    prodcheck *= C[0]->G[ul[ui].l].uvals[ul[ui].u];
  }
  if (fabs (log (prodcheck)) > 1e-5)
   IM_err(IMERR_MUTSCALEPRODUCTFAIL,"product of mutation scalars not close to or equal to 1: %lf",prodcheck); 
  
} //fixmutationscalars()

void add_priorinfo_to_output(char priorfilename[],int *fpstri, char fpstr[])
{
  int i;
  SP "\nParameter Priors\n");
  SP "-----------------\n");
  if (!calcoptions[LOADPRIORSFROMFILE])
  {
    SP "  Population size parameters maximum value : %.4lf \n", thetaprior);
    if (modeloptions[EXPOMIGRATIONPRIOR])
      SP "  Migration rate parameters exponential distribution mean : %.4lf \n",
        mprior);
    else
      SP "  Migration rate parameters maximum value: %.4lf \n", mprior);
    SP "  Splitting time : %.4lf\n", tprior);
    SP "\n");
  }
  else
  {
    SP"  Prior distribution terms given in file: %s\n\n",priorfilename);
    SP"  Splitting time maximum values \n");
    SP"\tPeriod\tPrior maximum value\n");
    for (i=0;i<numsplittimes;i++)
      SP"\t%s\t%.3lf\n",T[i].str,T[i].pr.max);
    SP"  Population size parameters maximum values \n");
    SP"\tPopulation\tPrior maximum value\n");
    for (i=0;i<numtreepops;i++)
      SP"\t%s\t%.3lf\n",C[ARBCHAIN]->itheta[i].str,C[ARBCHAIN]->itheta[i].pr.max);
    if (modeloptions[EXPOMIGRATIONPRIOR]==0)
    {
      SP"  Migration parameters maximum values \n");
      SP"\tMigration rate\tMaximum value\n");
      for (i=0;i<nummigrateparams;i++)
        SP"\t%s\t%.3lf\n",C[ARBCHAIN]->imig[i].str,C[ARBCHAIN]->imig[i].pr.max);
    }
    else
    {
      SP"  Migration parameter prior means (exponential priors)\n");
      SP"\tMigration rate\tPrior mean value\n");
      for (i=0;i<nummigrateparams;i++)
        SP"\t%s\t%.3lf\n",C[ARBCHAIN]->imig[i].str,C[ARBCHAIN]->imig[i].pr.expomean);
    }

  }
}
void
start_setup_L (char infilename[], int *fpstri, char fpstr[], int currentid)
{
  int li, i;

  /* get the number of loci and the number of populations from the top of the datafile */
  read_datafile_top_lines (infilename, fpstri, fpstr);
  /* setup a temporary struture to record how much variation there is,  used for picking starting values of mutation scalars */
  numsitesIS = static_cast<int **> (malloc (nloci * sizeof (int *)));
                                // rows in this matrix are malloced in setup_L
  uvals = static_cast<double **> (malloc (nloci * sizeof (double *)));   
  for (i = 0; i < nloci; i++)
  {
    numsitesIS[i] = static_cast<int *> 
                        (calloc ((size_t) npops + 1, sizeof (int)));
  }
  readdata (infilename, fpstri, fpstr, numsitesIS, currentid);

  for (li = 0; li < nloci; li++)
  {
    init_mutation_scalar_rec (li);
    init_g_rec (li);
    uvals[li] = static_cast<double *> 
                    (malloc (L[li].nlinked * sizeof (double)));
  }
  uval_preliminary_sum = setup_uval ();
  if (modeloptions[ADDGHOSTPOP])
  {
    if (npops + 1 > MAXPOPS)
    {
      IM_err (IMERR_INPUTFILEINVALID, 
              "ghost population makes number of population (%d) greater than MAXPOPS [%d]", 
              npops, MAXPOPS);
    }
     for (li = 0; li < nloci; li++)
     {
      L[li].samppop[npops] = 0;
     }
    npops++;
    numtreepops += 2;
    lastperiodnumber++;
    numsplittimes++;
   }
  /* no need for this stuff,  not used anywhere 
  gi_largestngenes = 0;
  gi_largestnumsites = 0;
  for (li = 0; li < nloci; li++)
  {
    if (gi_largestngenes < L[li].numgenes)
    {
      gi_largestngenes = L[li].numgenes;
    }
    if (gi_largestnumsites < L[li].numsites)
    {
      gi_largestnumsites = L[li].numsites;
    }
  } */

 return;
}                               //start_setup_L

void
start_setup_C (void)
{
  int ci, li;

  C = static_cast<chain **>             //points to an array of chains
            (malloc (numchainspp * sizeof (struct chain *)));     
  for (ci = 0; ci < numchainspp; ci++)
  {
    C[ci] = static_cast<chain *> (malloc (sizeof (struct chain)));
    sprintf (C[ci]->name, "%d", ci);
  }
  for (ci = 0; ci < numchainspp; ci++)
  {
    init_genealogy_weights (&C[ci]->allgweight);

    /* CR: 110516.2
     * change malloc to calloc so that memory would be initialized to
     * a known state before use.
     */
    C[ci]->G = static_cast<genealogy *> 
                (calloc ((size_t) nloci, sizeof (struct genealogy)));
    for (li = 0; li < nloci; li++)
    {
      init_genealogy_weights (&(C[ci]->G[li].gweight));
    }

  }
}                               // start_setup_C

void start_setup_poptree(char *ps)
{
  int ci,i;
  for (ci=0;ci<numchainspp;ci++)
  {
    if (modeloptions[POPTREETOPOLOGYUPDATE] == 0)
      setup_poptree (ci, ps);
    else
    {
      i = randposint(numpoptopologies); // pick a random tree 
  /* JH note - tried having an update in which topology was never updated within a chain,  but only updated by swapping in chains
  and then at the beginning all tree numbers were used evenly across chains.   Tried it to see if bugs were in the tree 
  updating code.  Dumb idea. Did not work.  This forces a uniform distribution of trees among chains,  and the different 
  chains no longer share the same state space.  */
      poptopologyproposedlist[i] += 1;
      setup_poptree (ci, alltreestrings[i]); 
        
    }
    if ((C[ci]->plist = static_cast<int **> (malloc (npops * sizeof (*C[ci]->plist)))) == NULL)
      IM_err (IMERR_MEM, "  plist malloc did not work.   npops %d, step %d",
              npops, step);
    for (i = 0; i < npops; i++)
    {
      if ((C[ci]->plist[i] =
          static_cast<int *> (malloc ((npops - i) * sizeof (*C[ci]->plist[i])))) == NULL)
        IM_err (IMERR_MEM,
                "  plist malloc did not work.   npops - i  %d, step %d",
                npops - i, step);
    }
    //assert (strlen(C[ci]->chainpoptreestring) >= 40);  
    //printf("%d %s\n",ci,C[ci]->chainpoptreestring);
    fillplist (ci);
    set_tvalues (ci);
  }

  //if (modeloptions[POPTREETOPOLOGYUPDATE]==1)  
    if (hiddenoptions[HIDDENGENEALOGY]==1)  // fixed 5_3_2017
    set_poptree_update_record();
}

void finish_setup_poptree(char topologypriorinfostring[])
{
  int ci;

  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    init_change_poptree(topologypriorinfostring);
  }
  for (ci=0;ci<numchainspp;ci++)
  {
    if (hiddenoptions[HIDDENGENEALOGY]==1)
    {
      fillancplist (ci);
    }
  }
}

void
finish_setup_C (int currentid)
{
  int ci, li, i, ai;
  int largestsamp;  
  int nosimmigration;
  int copy_u_chainnum = 0;
  int treeweightcallcode = 7;
  init_treeweight ();
  for (i = 0, li = 0; li < nloci; li++)
    if (i < L[li].numgenes)
      i = L[li].numgenes;
  largestsamp = i;  

  for (i = 0, li = 0; li < nloci; li++)
    if (i < L[li].numgenes)
      i = L[li].numgenes;
  largestsamp = i;
  nnminus1 = static_cast<double *> 
             (malloc ((largestsamp + 1) * sizeof (double)));
  nnminus1[0] = 0;
  for (i = 1; i <= largestsamp; i++)
  {
    nnminus1[i] = (double) (i) * ((double) i - 1);
  }


  if (hiddenoptions[HIDDENGENEALOGY]==0)
    init_gtreecommon ();
  else
  {
    init_gtreecommonhg ();
    setcopyedge();
  }
  
  if (modeloptions[NOMIGRATION]==0)
  {
    set_nomigrationchecklist ();
  }
  for (ci = 0; ci < numchainspp; ci++)
  {
    for (li = 0; li < nloci; li++)
    {
      C[ci]->G[li].gtree = static_cast<edge *> 
                        (malloc (L[li].numlines * sizeof (struct edge)));
      C[ci]->G[li].uvals = static_cast<double *> 
                        (malloc (L[li].nlinked * sizeof (double)));
      C[ci]->G[li].pdg_a = static_cast<double *> 
                        (malloc (L[li].nlinked * sizeof (double)));
    }
  }

  for (ci = 0; ci < numchainspp; ci++)
  {
    init_probcalc (&(C[ci]->allpcalc));
    for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].numlines; i++)
      {
        C[ci]->G[li].gtree[i].cmm = 0;
        C[ci]->G[li].gtree[i].mig[0].mt = -1;
        C[ci]->G[li].gtree[i].up[0] = -1;
        C[ci]->G[li].gtree[i].up[1] = -1;
        C[ci]->G[li].gtree[i].down = -1;
        C[ci]->G[li].gtree[i].time = 0;
        C[ci]->G[li].gtree[i].mut = -1;
        C[ci]->G[li].gtree[i].pop = -1;
        C[ci]->G[li].gtree[i].ei = i;
      //  C[ci]->G[li].gtree[i].exist = 'T';  does not appear to be used 5/12/2016
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          C[ci]->G[li].gtree[i].cmmhg = 0;
          C[ci]->G[li].gtree[i].mighg[0].mt = -1;
          C[ci]->G[li].gtree[i].pophg = -1;
        }
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          /* CR: 110516.1
           * change malloc to calloc so that memory would be initialized to
           * a known state before use.
           */
          C[ci]->G[li].gtree[i].A = static_cast<int *> 
                            (calloc ((size_t) L[li].nlinked, sizeof (int)));
          C[ci]->G[li].gtree[i].dlikeA = static_cast<double *> 
                            (calloc ((size_t) L[li].nlinked, sizeof (double)));
        }

      }
    }

    imaAsnInitAssign (ci);
    for (i = 0; i < 2 * npops - 1; i++)
    {
      if (C[ci]->poptree[i].e == -1)
        C[ci]->poptree[i].time = TIMEMAX;
      else
        C[ci]->poptree[i].time = C[ci]->tvals[C[ci]->poptree[i].e - 1];
    }

    if (ci == 0)                //&& nurates > 1)
    {
      setuinfo (uval_preliminary_sum);
    }
    copy_u_chainnum = 0;  // makes clear the uvals get copied from chain 0 
    for (li = 0; li < nloci; li++)
    {
      if (ci > 0 && nurates > 1)
      {
        for (ai = 0; ai < L[li].nlinked; ai++)
          C[ci]->G[li].uvals[ai] = C[copy_u_chainnum]->G[li].uvals[ai];
      }
      if (nurates == 1)
      {
        C[ci]->G[li].uvals[0] = 1.0;
      }

      nosimmigration = !(modeloptions[PARAMETERSBYPERIOD] || nummigrateparams == maxpossiblemigrateparams);
      int tw;
      switch (L[li].model)
      {
      case INFINITESITES:
        {
        makeIS (ci, li,nosimmigration);
        if (hiddenoptions[HIDDENGENEALOGY] == 1)
          makegenealogy_from_hiddengenealogy(ci,li);
#ifdef TURNONCHECKS
      checkgenealogy(ci,li,1);
#endif //TURNONCHECKS
        tw = treeweight (ci, li,treeweightcallcode);
#ifdef TURNONCHECKS
   checkgenealogy(ci,li,0);
#endif //TURNONCHECKS
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] = likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        break;
        }
      case HKY:
        for (i = 0; i < 4; i++)
        {
          C[ci]->G[li].pi[i] = pi[li][i];
        }
        C[ci]->G[li].kappaval = 2.0;    // starting kappa value
        makeHKY (ci, li,nosimmigration);
        if (hiddenoptions[HIDDENGENEALOGY] == 1)
          makegenealogy_from_hiddengenealogy(ci,li);
#ifdef TURNONCHECKS
    checkgenealogy(ci,li,1);
#endif //TURNONCHECKS
        tw = treeweight (ci, li,treeweightcallcode);
        C[ci]->G[li].pdg_a[0] = likelihoodHKY (ci, li, C[ci]->G[li].uvals[0], 
                                                C[ci]->G[li].kappaval, 
                                                -1, -1, -1, -1);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0];
        copyfraclike (ci, li);
        storescalefactors (ci, li);
        break;
      case STEPWISE:
        somestepwise = 1;

        for (ai = 0; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        makeSW (ci, li,nosimmigration);
        if (hiddenoptions[HIDDENGENEALOGY] == 1)
          makegenealogy_from_hiddengenealogy(ci,li);
        tw = treeweight (ci, li,treeweightcallcode);
        C[ci]->G[li].pdg = 0;
        for (ai = 0; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }

        break;
      case JOINT_IS_SW:
        for (ai = 1; ai < L[li].nlinked; ai++)
          for (i = 0; i < L[li].numgenes; i++)
            C[ci]->G[li].gtree[i].A[ai] = L[li].A[ai][i];
        somestepwise = 1;
        makeJOINT_IS_SW (ci, li,nosimmigration);
        if (hiddenoptions[HIDDENGENEALOGY] == 1)
          makegenealogy_from_hiddengenealogy(ci,li);
        tw = treeweight (ci, li,treeweightcallcode);
        C[ci]->G[li].pdg = C[ci]->G[li].pdg_a[0] =
          likelihoodIS (ci, li, C[ci]->G[li].uvals[0]);
        for (ai = 1; ai < L[li].nlinked; ai++)
        {
          C[ci]->G[li].pdg_a[ai] = likelihoodSW (ci, li, ai, C[ci]->G[li].uvals[ai], 1.0);
          C[ci]->G[li].pdg += C[ci]->G[li].pdg_a[ai];
        }
        break;
      }
#ifdef DEBUG
      C[ci]->G[li].hiprob = -1e20; // debugging 8/18/2016
#endif

      C[ci]->allpcalc.pdg += C[ci]->G[li].pdg;
      sum_treeinfo (&(C[ci]->allgweight), &(C[ci]->G[li].gweight));
      if (hiddenoptions[HIDDENGENEALOGY]==1)  // not sure if any of this is useful,  for now fix mhg in update_hg_2_16_2017.cpp
      {
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          if (hiddenoptions[UPDATEMRATEFORHGUPDATE])
            C[ci]->G[li].mhg = expo(1.0/mprior);
          else
            C[ci]->G[li].mhg = (1.0/mprior) * MPRIORFRACFORHG;
        }
        else
        {
          if (hiddenoptions[UPDATEMRATEFORHGUPDATE])
            C[ci]->G[li].mhg = uniform() * mprior;
          else
            C[ci]->G[li].mhg =  mprior * MPRIORFRACFORHG;
        }
        C[ci]->G[li].hgprob =  prob_hg_given_g(ci,li);
        C[ci]->allpcalc.probhgg +=   C[ci]->G[li].hgprob;
      }
    }

    initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
    if (doRYupdate)
    {
       C[ci]->RYwidthinfo = static_cast<updatescalarinfo *>  (malloc (numsplittimes * sizeof (struct updatescalarinfo)));
         /* set width for updates,  these are values to be multiple times the T prior to get the standard deviation of the slide
      target update rate is 0.25,  adjustor is 1.1,  
      starting valuue is 0.05,  with min of 0.00005 and max of 1.0*/
   
      for (i=0;i<numsplittimes;i++)
      {
        // updatescalarinfo, startwidth, adjustval,min,max,target,numattemptscheck
        // set the adjustval to 1 if you don't want to update it 
        setupdatescalarinfo(&C[ci]->RYwidthinfo[i],0.5,1.1,0.001,1.0,0.3,1000);
      }
     }

     if (doNWupdate)  // only do NW updates when not using hidden genealogies 
     {
        C[ci]->NWwidthinfo = static_cast<updatescalarinfo *>  (malloc (numsplittimes * sizeof (struct updatescalarinfo)));
        for (i=0;i<numsplittimes;i++)
        {
          // updatescalarinfo, width, adjustval,min,max,target,numattemptscheck
          //setupdatescalarinfo(&C[ci]->NWwidthinfo[i],0.05,1.0,0.00005,1.0,0.25,100); // no change
          setupdatescalarinfo(&C[ci]->NWwidthinfo[i],0.05,1.1,0.00005,1.0,0.3,1000);
        }
     }

  }

  for (ci = 0; ci < numchainspp; ci++)
    C[ci]->currallbetapos = currentid*numchainspp + ci;
  return;
}                               // finish_setup_C

/* set_tvalues() set random times for poptree */
  /* first pick trandom times over interval from 0 to tprior 
     use broken stick model - simple dirichlet distribution
   */
/* to fix values for debugging:
  in ima.h  turn on these undefs
  #undef DO_RY1UPDATE
  #undef DO_RY1UPDATE
  and at the end of set_tvalues()
  use code like 
  C[ci]->tvals[0]=0.1;C[ci]->tvals[1]=0.2;C[ci]->tvals[2]=0.4;C[ci]->tvals[3]=0.8; 
*/

void
set_tvalues (int ci)
{
  int i;
  double sum;
  double times[MAXPOPS];

  if (npops == 1)
  {
    C[ci]->tvals = static_cast<double *> (malloc (sizeof (double)));
    C[ci]->tvals[0] = TIMEMAX;
  }
  else 
  {
    C[ci]->tvals = static_cast<double *> 
                    (malloc ((lastperiodnumber + 1) * sizeof (double)));
    for (sum = 0, i = 0; i < npops; i++)
    {
      times[i] = expo (1.0);
      sum += times[i];
    }

    for (i = 0; i < lastperiodnumber; i++)
    {
      times[i] *= (T[i].pr.max - T[i].pr.min) / sum;
      times[i] += T[i].pr.min;
      if (i > 0)
        times[i] += times[i - 1];
      C[ci]->tvals[i] = times[i];
      assert (C[ci]->tvals[i] < T[i].pr.max && C[ci]->tvals[i] > T[i].pr.min);
    }
    C[ci]->tvals[i] = TIMEMAX;
  }
  return;
}                               /* set_tvalues */

/* JH 2/18/10 inserted this */
#define MAXTRANGE_IF_ISLANDMODEL  100.0  // have to pick something 
void
setup_T ()
{
  int i, li;
  int tui;


  time_update_type_count = doRYupdate + doNWupdate;
  tui = 0;
  if (doNWupdate)
  {
    update_type_NW = tui;
    tui += 1;
  }
  else
    update_type_NW = -1;
  if (doRYupdate)
    update_type_RY = tui;
  else
    update_type_RY = -1;

  if (npops == 1)
  {
    T = NULL;
    tprior = 0.0;
    for (li=0;li<nloci;li++)
    {
      L[li].g_rec->v->plotrange.max = MAXTRANGE_IF_ISLANDMODEL;
      init_value_record (L[li].g_rec->v, 0);
    }
  }
  else
  {
    T = static_cast<chainstate_record_updates_and_values *> 
            (calloc ((size_t) (npops - 1), 
             sizeof (struct chainstate_record_updates_and_values)));
    for (i = 0; i < lastperiodnumber; i++)
    {
      if (tprior > 0.0 )
        tperiodpriors[i] = tprior;
      T[i].pr.max = tperiodpriors[i];
      T[i].pr.min = 0;
    }
    //T[lastperiodnumber].pr.max = TIMEMAX;
    //T[lastperiodnumber].pr.min = 0.0;
    for (i = 0; i < lastperiodnumber; i++)
    {
      sprintf (T[i].str, "t%d", i);
      T[i].num_uptypes = time_update_type_count;
      T[i].upnames = static_cast<strnl *> 
                    (malloc (T[i].num_uptypes * sizeof (strnl)));
      if (doNWupdate)
      {
        sprintf (T[i].upnames[update_type_NW], "NielsenWakeley");
      }
      //else
        //sprintf (T[i].upnames[IM_UPDATE_TIME_NW], "NW(NOT_in_USE)");
      if (doRYupdate)
      {
        sprintf (T[i].upnames[update_type_RY], "RannalaYang");
      }
      T[i].upinf = static_cast<update_rate_calc *> 
                (calloc ((size_t) T[i].num_uptypes, 
                         sizeof (struct update_rate_calc)));
      T[i].num_vals = 1;
      T[i].v = static_cast<value_record *> 
                (malloc (T[i].num_vals * sizeof (struct value_record)));
      strcpy (T[i].v->str, T[i].str);
      strncpy (T[i].v->strshort, T[i].v->str, PARAMSTRLENSHORT - 1);
      T[i].v->do_xyplot = 1;
      T[i].v->do_trend = 1;
      T[i].v->do_logplot = 0;
      T[i].v->do_autoc = 1;
      T[i].v->plotrange.max = tperiodpriors[i];
      T[i].v->plotrange.min = 0;
      T[i].v->plotrescale = 1.0;
      init_value_record (T[i].v, 0);
      /* 2/18/10 JH  inserted this,  needed for plotting TMRCAs */
      if (i==0)
        for (li=0;li<nloci;li++)
        {
          L[li].g_rec->v->plotrange.max = tperiodpriors[lastperiodnumber-1] * TIMEPRIORMULTIPLIER;
          init_value_record (L[li].g_rec->v, 0);
        }
    }
  }

}                               //setup_T

void setup_migprior_recording()
{
  int i,j,k;
  int nm0terms=0;
  strn tempstr;
  mh = static_cast<chainstate_record_updates_and_values *> 
          (calloc ((size_t) (nummigrateparams), 
            sizeof (struct chainstate_record_updates_and_values)));
  for (i = 0; i <nummigrateparams; i++)
  {
    if (modeloptions[EXPOMIGRATIONPRIOR])
      mh[i].pr.max = EXPOMIGPLOTSCALE  * hyperprior_expo_m_mean;
    else
      mh[i].pr.max = hyperprior_uniform_m_max;
    mh[i].pr.min = 0;
  }
  /* put strings in for period 0 as these do not change */

  
  if (runmode == GHYPERPRIORmode2 || runmode == HGHYPERPRIORmode5)
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      sprintf (mh[i].str, "%s_h",C[ARBCHAIN]->imig[i].str);
    }
    nm0terms = 2*(npops*(npops-1));
  }
  else
  {
    k = 0;
    for (i = 0; i < npops - 1; i++)
      for (j = i + 1; j < npops ; j++)
      {
        int dir = makepairstring(C[ARBCHAIN]->descendantpops[C[ARBCHAIN]->plist[0][i]],C[ARBCHAIN]->descendantpops[C[ARBCHAIN]->plist[0][j]],tempstr);
        if (dir==0)
        {
          sprintf (mh[2*k].str, "m%d>%d_h", C[ARBCHAIN]->plist[0][i], C[ARBCHAIN]->plist[0][j]);
          sprintf (mh[2*k+1].str, "m%d<%d_h", C[ARBCHAIN]->plist[0][i], C[ARBCHAIN]->plist[0][j]);
        }
        else
        {
          sprintf (mh[2*k].str, "m%d<%d_h", C[ARBCHAIN]->plist[0][i], C[ARBCHAIN]->plist[0][j]);
          sprintf (mh[2*k+1].str, "m%d>%d_h", C[ARBCHAIN]->plist[0][i], C[ARBCHAIN]->plist[0][j]);
        }
        k += 1;
      }
    for (i = nm0terms; i < nummigrateparams; i++)  // these values can change
    {
      if (i>=nm0terms)
      {
        if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
        {
          nm0terms = 2*k;
          if (i % 2 ==0)
            sprintf (mh[i].str, "m%d/%d_h", i/2,i/2+1);
          else
            sprintf (mh[i].str, "m%d/%d_h", i/2+1,i/2);
        }
      }
    }
  }
  for (i = 0; i < nummigrateparams; i++)
  {
    mh[i].num_uptypes = NUM_PRIOR_UPDATE_RECORD_TYPES;
    mh[i].upnames = static_cast<strnl *> (malloc (mh[i].num_uptypes * sizeof (strnl)));
    sprintf (mh[i].upnames[PRIOR_UPDATE], "InPoptree");
    mh[i].upinf = static_cast<update_rate_calc *> 
              (calloc ((size_t) mh[i].num_uptypes, 
                        sizeof (struct update_rate_calc)));
    mh[i].num_vals = 1;
    mh[i].v = static_cast<value_record *> 
              (malloc (mh[i].num_vals * sizeof (struct value_record)));
    strcpy (mh[i].v->str, mh[i].str);
    strncpy (mh[i].v->strshort, mh[i].v->str, PARAMSTRLENSHORT - 1);
    mh[i].v->do_xyplot = 1;
    mh[i].v->do_trend = 0;
    mh[i].v->do_logplot = 0;
    //mh[i].v->do_autoc = 1;
    mh[i].v->plotrange.max = mh[i].pr.max ;
    mh[i].v->plotrange.min = 0;
    mh[i].v->plotrescale = 1.0;
    mh[i].v->do_autoc = 0;
    init_value_record (mh[i].v, 0);
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)
  {
    mhnit = static_cast<chainstate_record_updates_and_values *>  (malloc (sizeof (struct chainstate_record_updates_and_values)));
    if (modeloptions[EXPOMIGRATIONPRIOR])
      mhnit->pr.max = EXPOMIGPLOTSCALE  * hyperprior_expo_m_mean;
    else
      mhnit->pr.max = hyperprior_uniform_m_max;
    mhnit->pr.min = 0;
    sprintf (mhnit->str, "allm");
    mhnit->num_uptypes= NUM_PRIOR_UPDATE_RECORD_TYPES;
    mhnit->upinf = static_cast<update_rate_calc *>  (calloc ((size_t) mhnit->num_uptypes,sizeof (struct update_rate_calc)));
    mhnit->upnames = static_cast<strnl *> (malloc (mhnit->num_uptypes * sizeof (strnl)));
    sprintf (mhnit->upnames[PRIOR_UPDATE], "NotInPoptree");
    mhnit->num_vals = 0;
    mhnit->v  = NULL;
  }
} /* setup_migprior_recording */


void setup_qprior_recording()
{
  int i;
  qh = static_cast<chainstate_record_updates_and_values *> 
          (calloc ((size_t) (numpopsizeparams), 
            sizeof (struct chainstate_record_updates_and_values)));
  for (i = 0; i <numpopsizeparams; i++)
  {
    qh[i].pr.max = hyperprior_uniform_q_max;
    qh[i].pr.min = 0;
  }
  /* put strings in for period 0 as these do not change */
  for (i = 0; i <numpopsizeparams; i++)
  {
    sprintf (qh[i].str, "q%d_h", i);
    qh[i].num_uptypes = NUM_PRIOR_UPDATE_RECORD_TYPES;
    qh[i].upnames = static_cast<strnl *> (malloc (qh[i].num_uptypes * sizeof (strnl)));
    sprintf (qh[i].upnames[PRIOR_UPDATE], "InPoptree");
    qh[i].upinf = static_cast<update_rate_calc *> 
              (calloc ((size_t) mh[i].num_uptypes, 
                        sizeof (struct update_rate_calc)));
    qh[i].num_vals = 1;
    qh[i].v = static_cast<value_record *> 
              (malloc (qh[i].num_vals * sizeof (struct value_record)));
    strcpy (qh[i].v->str, qh[i].str);
    strncpy (qh[i].v->strshort, qh[i].v->str, PARAMSTRLENSHORT - 1);
    qh[i].v->do_xyplot = 1;
    qh[i].v->do_trend = 0;
    qh[i].v->do_logplot = 0;
    //qh[i].v->do_autoc = 1;
    qh[i].v->plotrange.max = qh[i].pr.max ;
    qh[i].v->plotrange.min = 0;
    qh[i].v->plotrescale = 1.0;
    qh[i].v->do_autoc = 0;
    init_value_record (qh[i].v, 0);
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)
  {
    qhnit = static_cast<chainstate_record_updates_and_values *>  (malloc (sizeof (struct chainstate_record_updates_and_values)));
    qhnit->pr.max = hyperprior_uniform_q_max;
    qhnit->pr.min = 0;
    sprintf (qhnit->str, "allq");
    qhnit->num_uptypes= NUM_PRIOR_UPDATE_RECORD_TYPES;
    qhnit->upinf = static_cast<update_rate_calc *>  (calloc ((size_t) qhnit->num_uptypes,sizeof (struct update_rate_calc)));
    qhnit->upnames = static_cast<strnl *> (malloc (qhnit->num_uptypes * sizeof (strnl)));
    sprintf (qhnit->upnames[PRIOR_UPDATE], "NotInPoptree");
    qhnit->num_vals = 0;
    qhnit->v  = NULL;
  }
} /* setup_qprior_recording */

  void
finish_setup_L (void)
{
  int li;
  for (li = 0; li < nloci; li++)
  {
    XFREE (numsitesIS[li]);
    XFREE (uvals[li]);
  }
  XFREE (numsitesIS);
  XFREE (uvals);
}                               // finish_setup_L

void reportparamcounts(int *fpstri, char fpstr[])
{
  SP"\nParameter Counts\n");
  SP"----------------\n");
  SP"  Population sizes : %d\n",numpopsizeparams);
  SP"  Migration rates  : %d\n",nummigrateparams);
  SP"  Parameters in the MCMC simulation\n");
  SP"      Splitting times : %d\n",numsplittimes);
  SP"      Mutation scalars: %d\n",nurates); 
  SP"      HKY Kappa (ti/tv) ratios: %d\n",nkappas); 
} // reportparamcounts(void)

/**********  GLOBAL FUNCTIONS  *******/

void 
setup (char infilename[], int *fpstri, char fpstr[], char priorfilename[],char topologypriorinfostring[],int currentid)
{
  int ci;
  start_setup_L (infilename, fpstri, fpstr, currentid);
  if (calcoptions[LOADPRIORSFROMFILE])
  {
    if (thetaprior < 0.0)
      popsizeprior_fromfile = static_cast<double *>  (malloc (numtreepops * sizeof(double)));
    else
      popsizeprior_fromfile = NULL;

    if (mprior < 0.0)// when is this true? 
      mprior_fromfile = alt2d_alloc2Ddouble(numtreepops,numtreepops);
    else
      mprior_fromfile = NULL;
    readpriorfile(priorfilename,popsizeprior_fromfile,mprior_fromfile);
  }
  
  setup_T ();
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      alltreestrings = allocalltreestrings();
      if (modeloptions[ADDGHOSTPOP]==1)
      {
        alltreestrings_noghost = allocalltreestrings();
      }
      numpoptopologies = buildpoptreestringarray();
      init_RF_nodeinfo();
      poptopologycounts = static_cast<int *> (calloc ((size_t) numpoptopologies, sizeof (int)));
      // poptopologysequence is only used on cpu 0,  but go ahead an initialize on all cpus
      poptopologysequence.vals =  static_cast<int *>  (malloc (POPTOPOLOGYSEQUENCELENGTH * sizeof(int)));
      poptopologysequence.disvals =  static_cast<double *>  (malloc (POPTOPOLOGYSEQUENCELENGTH * sizeof(double)));
      poptopologysequence.currentlength = 0;
      poptopologysequence.maxlength = POPTOPOLOGYSEQUENCELENGTH;
    }
  start_setup_C ();
  start_setup_poptree(&startpoptreestring[0]);
  getparamnums();
  if (nummigrateparams == 0 && hiddenoptions[HIDDENGENEALOGY]==1)
    IM_err (IMERR_MIGRATIONPRIOR0, " Migration must be included in model when using hidden genealogies (e.g. with population topology updating)");
  init_i_params();
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    hashsize = hashvalmaxes[npops];
  for (ci=0;ci<numchainspp;ci++)
  {
    setup_iparams (ci);
    set_iparam_poptreeterms (ci);
  }

  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    setup_migprior_recording();
    setup_qprior_recording();
  }
  finish_setup_poptree(topologypriorinfostring);
  finish_setup_C (currentid);
  reportparamcounts(fpstri, fpstr);
  finish_setup_L ();            // somethings in L need info from T , free up numsitesIS and uvals
  if (npops > 1  && hiddenoptions[HIDDENGENEALOGY]==0)
  {
    init_t_NW ();
    init_t_RY ();
  }
  if(hiddenoptions[HIDDENGENEALOGY]==1)
    init_t_RYhg();
  init_updategenealogy ();
  init_lpgpd_v ();
  if (outputoptions[MIGRATEHIST]) 
    init_migration_counts ();
  add_priorinfo_to_output(priorfilename,fpstri, fpstr);
  init_autoc_pointers ();

  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  if (calcoptions[CALCMARGINALLIKELIHOOD])
    initmarginlikecalc();
#ifdef TURNONCHECKS
  //for (ci=0;ci<numchainspp;ci++)  checkpoptree(ci,0);
for (ci=0;ci<numchainspp;ci++)
  checkgenealogyweights(ci);
#endif //TURNONCHECKS
   
}                               // setup() 
