/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/* 

Calculating autocorrelations and ESS values:
---------------------------------------------
the general idea is to initialize things so that autocorrelations and ESS values can be calculated for a bunch 
of diverse things. 
At the start a set of three arrays of pointers are set up.  These point to the information needed, 
then they are passed to the function that does the calculation and the printing of values. 

void init_autoc_pointers(void):  initializes three arrays:
static struct autoc **autoc_pointer - an array of pointers to struct autoc which holds the information needed to do the calculations
static char **autoc_str_pointer - an array of pointers to strings that hold the names of things for which calculationsa 
are being done 
static double *autoc_vals - an array of numerical values that are being autocorrelated

See also 
static struct autoc **autoc_a_pointer = NULL;
static char **autoc_a_str_pointer = NULL;
static double *autoc_a_vals = NULL;
These are the same kinds of things, but for assignment related work

How to have something included in autocorrelation and ESS calculations and printouts:
-------------------------------------------------------------------------------------
The only things that should be needed, to add something new to the list of things to have 
autocorrelations and ESS values calculated is to have them included in the list of pointers in init_autoc_pointers() 
and to have the values to be included in the autocorrelation identifed in set_autoc_vals() 
Once things are set up,  the rest of the work should be done. 
So the only place changes should be needed is in init_autoc_pointers() and set_autoc_vals()

*/

static struct autoc **autoc_pointer = NULL;
static char **autoc_str_pointer = NULL;
static double *autoc_vals;
static int num_autoc = 0;

/* For ASSIGNMENT */
static struct autoc **autoc_a_pointer = NULL;
static char **autoc_a_str_pointer = NULL;
static double *autoc_a_vals = NULL;
static int num_autoc_a = 0;

/* these are the lag intervals overwhich autocorrelations are measured
 on 2/9/09  added AUTOCSTEPSCALAR as a multiplier of these to make 
 measurements over a much wider interval */

static int autoc_checkstep[AUTOCTERMS] =
  { 1, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };

/****** LOCAL FUNCTIONS ***********/
static void fillautoc (struct autoc *a, int n, double v);
static double calcautoc (struct autoc a);
static double printautocvalue (FILE * outto, struct autoc *ta);
static void integrate_autoc (FILE * outto, struct autoc ta[], double ac[]/*,int step */);
void printautoctable (FILE * outto, int numautoc, struct autoc **ac, char **ac_str);
void iautoc (struct autoc *a);  // set autoc values to zero
void set_autoc_vals (double *pautoc_vals, int currentid);
// struct autoc  used for recording values for measuring autocorrelations

/* fillauto() adds new values to the cov and var terms that are needed to calculate an autocorrelation
these are calculated using the value v, which is the current value of the variable being autocorrelated and a value
that was saved previously in position n of the the vals[] array */
void
fillautoc (struct autoc *a, int n, double v)
{
  a->cov.n++;
  a->cov.s += a->vals[n];
  a->cov.ss += v;
  a->cov.s2 += (a->vals[n] * v);
  a->var[0].s += a->vals[n];
  a->var[1].s += v;
  a->var[0].s2 += SQR (a->vals[n]);
  a->var[1].s2 += SQR (v);
}

/* calculated a correlation using the cov and var terms in a struct autoc */
double
calcautoc (struct autoc a)
{
  double ac;
  /* there is a bug that sometimes causes a negative value to go to SQR and generate Nan*/
  ac =  (a.cov.s2 - (a.cov.s * a.cov.ss) / a.cov.n) / sqrt ((a.var[0].s2 - SQR (a.var[0].s) / a.cov.n) * (a.var[1].s2 - SQR (a.var[1].s) / a.cov.n));
  if (isnan_(ac))
  {
    return 0.0;
    //IM_err(IMERR_MISCELLANEOUS," autocorrelation calculation returned NaN. var0s %.4lf var0s2 %.4lf var1s %.4lf var1s2 %.4lf covn %d",a.var[0].s,a.var[0].s2,a.var[1].s,a.var[0].s2,a.cov.n);
  }
  return ac;
}

double
printautocvalue (FILE * outto, struct autoc *ta)
{
  double ac;
  ac = calcautoc (*ta);
  if (ac > 1)
    ac = 1;
  if (isnan_(ac))
  {  // should not get here if calcautoc traps Nan
    fprintf (outto, "\tNan");
  }
  else
    fprintf (outto, "\t% .4f", ac);
  return ac;
}

/* integrate over autocorrelations values for values > AUTOCMIN */
#define AUTOCMIN  0.03
void
integrate_autoc (FILE * outto, struct autoc ta[], double ac[]/*, int step*/)
{
  double ess, temp1, temp2;
  long last_checkstep;
  int i, stopsum;
  ess = 0;
  stopsum = 0;
/* if AUTOCSTEPSCALAR is > 1 then we use the autocorrelation at step 0 and assume that it is equal to 1.
otherwise we use calculated autocorrelation measured for a lag of 1 */

  if (AUTOCSTEPSCALAR > 1)
    i = 0;
  else
    i = 1;
  for (; i < AUTOCTERMS; i++)
  {
    if (ta[i].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))
    {
      if (ac[i] < AUTOCMIN)
        temp1 = 0;
      else
        temp1 = ac[i];
      if (AUTOCSTEPSCALAR > 1 && i == 0)
      {
        temp2 = 1.0;
        last_checkstep = 0;
      }
      else
      {
        if (ac[i - 1] < AUTOCMIN)
          temp2 = 0;
        else
          temp2 = ac[i - 1];
        last_checkstep = autoc_checkstep[i - 1];
      }
      if (stopsum == 0)
        ess += AUTOCSTEPSCALAR * (autoc_checkstep[i] - last_checkstep) * (temp1 + temp2) / 2;
      stopsum = (stopsum || temp1 == 0);
    }
  }
  if (stopsum)
  {
    if (outto != NULL)  fprintf(outto, "\t%.0f", step / (1 + 2 * ess));
  }
  else
  {
    if (outto != NULL)  fprintf(outto, "\t< %.0f", step / (1 + 2 * ess));
  }
}                               //integrate_autoc 

#define IFSTDOUTMAXSHOW 10
void
printautoctable (FILE * outto, int numautoc, struct autoc **ac, char **ac_str)
{
  int i, j;
  double **acvals;
  int numacprint;
  char numstr[20];
  numacprint = (outto == stdout) ? IMIN (numautoc, IFSTDOUTMAXSHOW) : numautoc;

  acvals = orig2d_alloc2Ddouble (numacprint, AUTOCTERMS);
  if (outto != NULL)  
  {
    fprintf(outto,"\nAutocorrelations and Effective Sample Size Estimates\n");
    fprintf(outto,  "----------------------------------------------------\n");

  }
  if (outto != NULL)  fprintf(outto,
           "   # Steps Between Values and Autocorrelation Estimates \n");
  if (outto != NULL)  fprintf(outto, "\tSteps ");
  for (j = 0; j < numacprint; j++)
    if (outto != NULL)  fprintf(outto, "\t %s", ac_str[j]);
  if (outto != NULL)  fprintf(outto, "\n");
  for (i = 0; i < AUTOCTERMS; i++)
  {
    if (ac[0][i].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))
    {
      if (((float) autoc_checkstep[i] * AUTOCSTEPSCALAR) > 1e4)
      {
        sprintf (&numstr[0], "%.1e",
                 (float) autoc_checkstep[i] * AUTOCSTEPSCALAR);
        if (outto != NULL)  fprintf(outto, "\t%s", shorten_e_num (&numstr[0]));
        //fprintf (outto, "\t%.1e", (float) autoc_checkstep[i]*AUTOCSTEPSCALAR);
      }
      else
      {
        if (outto != NULL)  fprintf(outto, "\t%d", autoc_checkstep[i] * AUTOCSTEPSCALAR);
      }
      for (j = 0; j < numacprint; j++)
      {
        acvals[j][i] = printautocvalue (outto, &ac[j][i]);
      }
      if (outto != NULL)  fprintf(outto, "\n");
    }
  }
  if (outto != NULL)  fprintf(outto, "\tESS");

/* integrate over autocorrelations values for values > 0.03 */

  for (j = 0; j < numacprint; j++)
    integrate_autoc (outto, ac[j], acvals[j]/*, step*/);
  if (outto != NULL)  fprintf(outto, "\n");
  orig2d_free2D ((void **) acvals, numacprint);
}                               //printautoctable

void
iautoc (struct autoc *a)        // initialize structure
{
  int i;
  for (i = 0; i < AUTOCTERMS; i++)
  {
    ieevent (&a[i].cov);
    ieevent (&a[i].var[0]);
    ieevent (&a[i].var[1]);
    for (int k = 0;k < AUTOCNEXTARRAYLENGTH;k++) // added zeroing of the vals array,  4/18/2018
          a[i].vals[k] = 0.0;
  }

}                               //iautoc 

/* this just puts the current values that are needed for the autocorrelations all in one array,  on the head node, temporarily
have to do a sends and receives if not on head node
pautoc_vals holds the values,  in order, for each thing that has autocorrelation calculated 
there are num_autoc values in pautoc_vals

pautoc_vals  points to autoc_vals in calling function
the order of things in this array should be the same as for autoc_pointer
autoc_pointer: order of values
	poptreeuinfo->v
	lpgpd_v
	T[j].v
	g_rec->v

 set_autoc_vals  
  for a value for which the autocorrelation is calculated
    there will be only 1 value across the run at any point in time, and that is the one associated with the cold chain

    this function returns a pointer but the values are meaningless unless it is the Headnode 

  if no MPI and #cpu > 1, then always true that z >= 0 && currentid == HEADNODE

  if numprocesses > 0:
    if z>=0 && currentid == HEADNODE: copy value to pautoc_vals
    if z >= 0 && currentid != HEADNODE: send value to pautocrec
    if z < 0 && currentid == HEADNODE: receive value to pautocrec
    else pautoc_vals is not relevant and holds nothing of value 

  call this kind of operation PASS_TO_HEADNODE_TO_SAVE
*/
void
set_autoc_vals (double *pautoc_vals, int currentid)
{
  int j, li, i;
  int z = whichiscoldchain();
#ifdef MPI_ENABLED
   if (z < 0 && currentid != HEADNODE) 
     return;  // nothing to do in this case
   int rc = 0; 
   MPI_Status status;
   double pautocrec;
#endif
    // follow same sequence as used with autoc_pointer and autoc_vals arrays
    //  poptreeuinfo (if present)
    //  lpgpd
    //  T (npops - 1 vals)
    //  g_rec  (nloci vals) 
  i=0;
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1 && poptreeuinfo->v->do_autoc == 1)
  {
    if (z >= 0 && currentid == HEADNODE) 
    {	
      pautoc_vals[i] =  (double) RFtreedis[C[z]->poptreenum];
    }
#ifdef MPI_ENABLED
    if (z < 0 && currentid == HEADNODE) 
    {
      rc = MPI_Recv(&pautocrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 6561, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
      pautoc_vals[i] = pautocrec;
    }
    if (z >= 0 && currentid != HEADNODE) 
    {
      pautocrec = (double) RFtreedis[C[z]->poptreenum];
      rc = MPI_Send(&pautocrec, 1, MPI_DOUBLE, 0, 6561, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif
    i++;
  }
  

  if (lpgpd_v->do_autoc == 1)
  {
    if (z >= 0 && currentid == HEADNODE) 
    {	
      if (hiddenoptions[HIDDENGENEALOGY]==0)
        pautoc_vals[i] = C[z]->allpcalc.pdg  + C[z]->allpcalc.probg;
      else
        pautoc_vals[i] = C[z]->allpcalc.pdg  + C[z]->allpcalc.probg + C[z]->allpcalc.probhgg;
    }
#ifdef MPI_ENABLED
    if (z < 0 && currentid == HEADNODE) 
    {
	    rc = MPI_Recv(&pautocrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 5656, MPI_COMM_WORLD, &status);
	    if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	    pautoc_vals[i] = pautocrec;
    }
    if (z >= 0 && currentid != HEADNODE) 
    {
      if (hiddenoptions[HIDDENGENEALOGY] == 0)
        pautocrec = C[z]->allpcalc.pdg + C[z]->allpcalc.probg;
      else
        pautocrec = C[z]->allpcalc.pdg + C[z]->allpcalc.probg + C[z]->allpcalc.probhgg;
    	rc = MPI_Send(&pautocrec, 1, MPI_DOUBLE, 0, 5656, MPI_COMM_WORLD);
	    if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
    }
    #endif
    i++;
  }
  

  for (j = 0; j < numsplittimes; j++)if (T[j].v->do_autoc == 1)
  {
    if (z >= 0 && currentid == HEADNODE) 
    {
	     pautoc_vals[i] = C[z]->tvals[j];
	   }
#ifdef MPI_ENABLED
    if (z < 0 && currentid == HEADNODE) 
    {
		    rc = MPI_Recv(&pautocrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 6767, MPI_COMM_WORLD, &status);
		    if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		    pautoc_vals[i] = pautocrec;
	   }
	   if (z >= 0 && currentid != HEADNODE) 
    {
		    pautocrec = C[z]->tvals[j];
		    rc = MPI_Send(&pautocrec, 1, MPI_DOUBLE, 0, 6767, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
	   }
#endif
    i++;
  }

  for (li = 0; li < nloci; li++) if (L[li].g_rec->v->do_autoc == 1)
  {
	   if (z >= 0 && currentid == HEADNODE) 
    {
	     pautoc_vals[i] = C[z]->G[li].roottime;
	   }
#ifdef MPI_ENABLED
	   if (z < 0 && currentid == HEADNODE) 
    {
		    rc = MPI_Recv(&pautocrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 7878, MPI_COMM_WORLD, &status);
		    if (rc != MPI_SUCCESS)MPI_Abort(MPI_COMM_WORLD, rc);
		    pautoc_vals[i] = pautocrec;
	   }
	   if (z >= 0 && currentid != HEADNODE) 
    {
		    pautocrec = C[z]->G[li].roottime;
		    rc = MPI_Send(&pautocrec, 1, MPI_DOUBLE, 0, 7878, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);		
	   }
#endif
    i++;
  }
}                               //set_autoc_vals


/*********** GLOBAL FUNCTIONS *************/

void
free_autoc_pointers (void)
{

  XFREE (autoc_pointer);
  XFREE (autoc_str_pointer);
  XFREE (autoc_vals);
}                               // free_autoc_pointers 

/* num_autoc  counts up how many things have autocorrelation calculations
   any struct value_record has an int do_autoc
   if do_autoc == 1  then autocorrelation is done 
   num_autoc counts these up

   those for which do_autoc == 1:
    poptreeuinfo->v->do_autoc
    (lpgpd_v->do_autoc == 1);
    T[j].v->do_autoc
    (L[li].g_rec->v->do_autoc == 1)


   some value_records for which do_autoc == 0:
    mh[i].v
    qh[i].v
    L[li].u_rec[ui].v[i]
    migration_counts[j][i]
    L[li].kappa_rec->v
*/

/* 
autoc_pointer: order of values
	poptreeuinfo->v
	lpgpd_v
	T[j].v
	g_rec->v

*/ 
void
init_autoc_pointers (void)
{
  int i, j, li;
  if (autoc_pointer == NULL)
  {
    // count up the number of things to record ESS values 
    // follow same sequence as used with autoc_pointer and autoc_vals arrays
    //  poptreeuinfo (if present)
    //  lpgpd
    //  T (npops - 1 vals)
    //  g_rec  (nloci vals) 
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
      num_autoc += poptreeuinfo->v->do_autoc == 1;
    num_autoc += (lpgpd_v->do_autoc == 1);      // lpgpd
    for (j = 0; j < numsplittimes; j++)
      num_autoc += (T[j].v->do_autoc == 1);
    for (li = 0; li < nloci; li++)
      num_autoc += (L[li].g_rec->v->do_autoc == 1);     // tmrca values
    autoc_pointer = static_cast<autoc **> (malloc (num_autoc * sizeof (struct autoc *)));
    autoc_str_pointer = 
            static_cast<char **> (malloc (num_autoc * sizeof (char *)));
    autoc_vals = static_cast<double *> (malloc (num_autoc * sizeof (double)));
    i = 0;
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1 && poptreeuinfo->v->do_autoc == 1)
    {
      //autoc_pointer[i] = &(poptreeuinfo->v->ac[0]);
      autoc_pointer[i] = poptreeuinfo->v->ac;
      autoc_str_pointer[i] = poptreeuinfo->v->strshort;
      i++;
    }
    if (lpgpd_v->do_autoc == 1)
    {
      //autoc_pointer[i] = &(lpgpd_v->ac[0]);
      autoc_pointer[i] = lpgpd_v->ac;
      autoc_str_pointer[i] = lpgpd_v->strshort;
      i++;
    }
    for (j = 0; j < numsplittimes; j++)
    {
      if (T[j].v->do_autoc == 1)
      {
        //autoc_pointer[i] = &(T[j].v->ac[0]);
        autoc_pointer[i] = T[j].v->ac;
        autoc_str_pointer[i] = T[j].v->strshort;
        i++;
      }
    }
    for (li = 0; li < nloci; li++)
      if (L[li].g_rec->v->do_autoc == 1)
      {
        //autoc_pointer[i] = &(L[li].g_rec->v->ac[0]);
        autoc_pointer[i] = L[li].g_rec->v->ac;
        autoc_str_pointer[i] = L[li].g_rec->v->strshort;
        i++;
      }
    assert (i == num_autoc);
    for (i = 0;i< num_autoc;i++)
      iautoc (autoc_pointer[i]) ; 
 
  }
  return;
}                               // init_autoc_pointers

/*
records things for calculating autocorrelations  - tricky code 
called with start_autocorrelations == 1 at beginning and after burn 
AUTOCTERMS is a list of the different lag values in steps for which autocorrelations are recorded
nextstepcalc[AUTOCTERMS] for each lag value, the next step number at which another term is accumulated for the autocorrelation term for that lag
nextpossave[AUTOCTERMS] for each lag value, the position in the autoc_pointer[][].vals array that should get the next value to be recorded
nextposcalc[AUTOCTERMS] for each lag value, the position in the autoc_pointer[][].vals array that is involved in the next autocorrelation calculation
maxpos[AUTOCTERMS] for each lag value, the length of the autoc_pointer[][].vals array 

Much of this function is not used on nodes that are not the head node.
With multiple cpus  we only need to accumulate autoc stuff on the head node
but values will have to be passed to the head node from other nodes when they host the cold chain
those values are put into an array by set_autoc_vals (autoc_vals, currentid) 
so set_autoc_vals (autoc_vals, currentid) needs to be called from within this code
even though those value are only used when this function is running on the headnode 
*/

void
checkautoc (int start_autocorrelations, int burndone, int burninsteps, int currentid)
{
  int i, j;
  int dofillautoc,z;
  int autoc_vals_recorded = 0;
  static int nextstepcalc[AUTOCTERMS],
    nextpossave[AUTOCTERMS], nextposcalc[AUTOCTERMS], maxpos[AUTOCTERMS];
#ifdef MPI_ENABLED
  z = whichiscoldchain();
  if (currentid != HEADNODE)
    dofillautoc = 0;
  else
    dofillautoc = 1;
#else 
  dofillautoc = 1;
#endif

  if (start_autocorrelations == 1)
  {
    // if not loading autoc values from mcf files,  initialize them
    if (runoptions[LOADMCSTATE]==0) for (i = 0; i < num_autoc; i++)
      iautoc (autoc_pointer[i]); 
    for (i = 0; i < AUTOCTERMS; i++)
    {
      nextstepcalc[i] = CHECKAUTOCWAIT + AUTOCINT * AUTOCSTEPSCALAR + autoc_checkstep[i] * AUTOCSTEPSCALAR + (burndone * burninsteps);
      nextpossave[i] = 0;
      nextposcalc[i] = 0;
      if (autoc_checkstep[i] <= AUTOCINT)
        maxpos[i] = 0;
      else
        maxpos[i] = (autoc_checkstep[i] / AUTOCINT) - 1;
      assert(maxpos[i] >= 0 && maxpos[i] < AUTOCNEXTARRAYLENGTH);
    }
  }
  else
  {
    for (i = 0; i < AUTOCTERMS; i++)
    {
      if (step == nextstepcalc[i])
      {
        if (!autoc_vals_recorded)
        {
         set_autoc_vals (autoc_vals, currentid);
         autoc_vals_recorded = 1;
        }
        if (currentid == HEADNODE) for (j = 0; j < num_autoc; j++)
        {
          fillautoc (&autoc_pointer[j][i], nextposcalc[i], autoc_vals[j]);  // record values 
        }
        nextposcalc[i]++;
        if (nextposcalc[i] > maxpos[i])
          nextposcalc[i] = 0;
        nextstepcalc[i] += AUTOCINT * AUTOCSTEPSCALAR;
      }
    }
    if ((long) AUTOCINT * AUTOCSTEPSCALAR * (long) (step / (AUTOCINT * AUTOCSTEPSCALAR)) == step)
    {
      if (!autoc_vals_recorded)
      {
        set_autoc_vals (autoc_vals, currentid);
      }

      for (i = 0; i < AUTOCTERMS; i++)
      {
        assert(nextpossave[i]>= 0 && nextpossave[i] < AUTOCNEXTARRAYLENGTH);
        for (j = 0; j < num_autoc; j++)
        {
          autoc_pointer[j][i].vals[nextpossave[i]] = autoc_vals[j];
        }
        nextpossave[i]++;
        if (nextpossave[i] > maxpos[i])
          nextpossave[i] = 0;
      }
    }
  }
}                               /*checkautoc */

/* this is called from both intervaloutput() and printoutput() */
void
callprintautoctable (FILE * outto/*, int step */)
{

  if (autoc_pointer[0][0].cov.n > (AUTOCCUTOFF / AUTOCSTEPSCALAR))      // check one of the autocorrelation accumulators to see if enough counts have been made
  {
    printautoctable (outto, num_autoc, autoc_pointer, autoc_str_pointer);
  }
}                               // callprintautoctable 
