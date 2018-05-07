/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/* JH added 9/17/09
based on IMadeopt program
uses differential evolution to find joint posterior density for full and nested models */

/*

Differential evolution algorithm
--------------------------------

This is a conventional differential evolution algorithm.
Many things were tried to see if it could be improved (e.g. multiple populations with migration,
seeding population with diverse runs using smaller samples of genealogies,
lots of playing with the weight and recombination terms of the algorithm,
(0.8 and 0.9 seem best, respectively)
lots of playing with the size of the population (100 seems best).
All to no avail.  Don't mess with it.

Models for joint parameter estimation:
--------------------------------------
for 2 population models all tests are done against the full 5 parameter set
All nested models are tested against this full model.
A nested model may include both p and m terms.

For models with more than two populations,  two models are developed
1)The full model for population size parameters includes all of the population sizes
2)The full model  for the migration parameters includes only all of the migration rate parameters.
A nested model may include either p terms or m terms, but not both.
Depending on whether it includes p or m terms it is tested against the corresponding full model.

Define model types as of three types:
0)  both m and p
1)  only p
2)  only m

Formatting of the nested model file
-----------------------------------
the model file consists of lines each of which begins with one of three words ('model' 'equal' or 'constant'):

the first line is an integer with the number of models to be read in

model   - indicates that a new model is being specified,  the string 'model ' is followed by the name of the model, a user specified string
  -all text that follows this up to the end of the line is treated as the name of the model, each time this word appears, a new model is started
equal   - a set of parameters that are always the same in value
  - followed by either an 'm' or a 'p'  for migration or population
  - then followed by two or more integers for the parameters that are equal to each other
  - each integer refers to a column position in the ti file whether or not they are population
  - if 'p' the integers are the ids for the populations
  - if 'm' the integers are the ids for the migration parameters, in the same order as they appear in full IMa2 analysis.
  - if 's'  the it is a constant term for the split time parameter
constant -  set of parameters that are equal to a constant
  - followed by either an 'm' or a 'p'  for migration or population
  - then followed by a number, that is the value that the parameters are fixed at
  - then followed by one or more integers that indicate which parameters are fixed at that value
  - if 'p' the integers are the ids for the populations
  - if 'm' the integers are the ids for the migration parameters, in the same order as they appear in full IMa2 analysis.

If there are two populations,  then having one or more nested models invokes a full model of all 5 parameters as the first model to be run.

If there are three or more populations, then a nested model can have only p terms or m terms but not both.
If a nested model has m terms then a full model of all migration parameters is also run.
If a nested model has p terms then a full model of all population size parameters is also run

e.g. for a two population case this is a nested model with both m and p terms in the nested model
1
model  sizes of populations  0, 1 and 2 equal to each other, no migration between 0 and 1
equal p 0  1  2    // meaning that 1st, 2nd and third populations are set to always have the same parameter values
constant m  0.0  0 1  // meaning that the first and second migration parameters in the model are set to a constant value of 0.0

understanding xmap
------------------
xmap contains information on how to copy one array for a nested model into a larger array for the full model
the length of xmap is as many parameters as there are in the largest (non-nested) model.
if xmap[i] is a negative value,  then the absolute value of this is plugged into the full model array at position i as a constant
if xmap[i] is a positive value, then that value is treated as an integer and position i in the full model array is copied over
  from position xmap[i] in the nested model array

the differential evolution algorithm works on an array of varying parameters
but if we are using a nested model the length of this list will be shorter than
the total number of parameters (because some are shared or identical).

parameters are counted:  first population size,  then migration, and finally a splitting rate if it is in the model.
total number of parameters in a full model is nparams

Using xmap to generate a full list for function calculation, using a shorter list of varying numbers:
  given a list f of varying numbers,  corresponding to the things that can vary in a nested model
  if we want to generate a full length list x,  that has numbers in all parameter positions:
	for (i=0;i<nparams;i++)
	{
		if (xmap[i] < 0)
			x[i] = -xmap[i];
		else
			x[i] = f[(int) xmap[i]];
	}
Using xmap to make a shorter list of varying numbers under a nested model,  given a full list of nparam numbers
  the full list is x and shorter list is f:
	for (i = 0,toi = 0;i<nparams;i++)
	{
		{
			if (toi == (int) xmap[i])
			{
				f[toi] = x[i];
				toi++;
			}
		}
	}

p: 0 1 2  m: 3  4
e.g. 2 pop model with both migration rates zero

[0,1,2,-MINPARAMVAL,-MINPARAMVAL]

e.g. 2 pops with all population sizes equal
[0,0,0,1,2]

e.g. 2 pop model with 3 popsize, 2 migration and one splitrate term

p: 0 1 2  m: 3  4   s:  5

xmap for a nested model with equal population sizes, freely varying migration and splitrate
[0,0,0,1,2,3]

xmap for a nested model with two equal migration rates,  other things vary
[0,1,2,3,3,4]

xmap for a nested model with a constant for population 1 of 2.5 and a constant for splitrate of 0.7
[0,-2.5,1,2,3,-0.7]

*/

#define STARTHI  1E200
#define STARTLOW  -1E200
#define SPREADTOL 1.0e-7 // criteria for stopping search for peak
#define PLOOPTOL  1.0e-6  // criteria for whether or not current peak is the same as best previously found
#define LOOPMATCHCRITERIA 2  // number of times the best peak must be found before stopping
#define MAXRESTART 10     // maximum number of restarts in search for peak
#define MAXJOINTMODELS 100  // maximum number of models in nested model file
#define MAXPLENGTH 200   // largest number of parameters in the model, should be good for models with up to 10 sampled populations
#define MAXMODELTEXTLINE 500 // max length of name of model in nestedmodelfile
#define DEFAULTRECRATE 0.9  // recombination rate in differential evolution algorithm.  played wth lots of values this seems best compromise
#define DEFAULTFWEIGHT 0.8  // mutation weight term in differential evolution algorithm.  played wth lots of values this seems best compromise
#define DEFAULTPOPSIZEMULTIPLIER 100  // population size multiplier in differential evolution algorithm.  should be 200
#define MAXPTERMS 20  // maximum number of 'constant' and 'equal' lines in a nested model in the nested model file

/* SANGCHUL: Thu Nov 12 22:17:34 EST 2009
 * There are too many global variables. Could these be file static by prepending
 * static?
 */
static double *bestvals, *tempvals;
static double *parameter_lower_bound;
static double *parameter_upper_bound;
static double **pop, **trialpop;
static double *popnest;
static int nparams;
static int popsizemultiplier = DEFAULTPOPSIZEMULTIPLIER;
static int num_nestedmodels;
static int depopsize;
static double holdxmaps[MAXJOINTMODELS][MAXPLENGTH];
static int paramsnotused[MAXJOINTMODELS][MAXPLENGTH]= {{0}};
static int paramsallused[MAXPLENGTH]= {0};
static int numpopsizeparams_nested[MAXJOINTMODELS];
static int nummigrateparams_nested[MAXJOINTMODELS];
static int nparams_nested[MAXJOINTMODELS];
static int nestedmodel_boundary_in_model[MAXJOINTMODELS] ={0};
static int atleastone_boundary_in_model=0;
static char modelnamelines[MAXJOINTMODELS][MAXMODELTEXTLINE];
static double recrate = DEFAULTRECRATE;
static double fweight = DEFAULTFWEIGHT;
static char logpfstr[40];
static int calculate_ess = 0;
static double effective_sample_size;
static double *logx, *divx, *log2diffx;//used by jointp
static double **double_gsamp_fcp;
static double *probgp_popsize, *probgp_migrate;
static int fullmodeltype[MAXJOINTMODELS] = {-1};
static int modeltypecounts[3] = {0}; // 0, both m and p;  1, only p;  2, only m
static char modelstartstr[3][40] = {"FULL","ALL POPULATION SIZE PARAMETERS", "ALL MIGRATION PARAMETERS"};
static int nparamrange[2];
static int nowmodeltype;
static double parameter_lower_bound_mapped[MAXPLENGTH];
static double parameter_upper_bound_mapped[MAXPLENGTH];


// listelement is for a linked list
struct listelement
{
  double v;
  struct listelement *l;
  struct listelement *r;
}**list, *first, *last, *now;

/* prototypes */
char *logpstrformat (double pval);
int fillparamlist(char *c, int *paramlist, int startparamcount, int notused[]);
int reduce_mappos(int mappos[], int k, int len);
void setup_mapping(char *nestedmodelfilename);
void reversemapvals(double *from, double *to, double *xmap);
void mapvals(double *from, double *to, double *xmap);
void setbounds();
void nextgen(double *lowpd, double *hipd, int modelparams, int modelnum);
void startpop(int modelparams,int modelnum);
double difeloop(int modelparams,int modelnum);
void modelloop(int modelnum);
void copybest(/* int modelparams */);
void printjointpeakvals(FILE *outfile,int notused[], double *printvals);
void startjointpeakouttable(FILE *outfile,char *fname);
double jointp (double *x,int calc_ess, double *effective_n);
void freejointpmem();
void addfirst(double val, int k);
void addlast(double val, int k);
void addright(double val, int k);
void addleft(double val, int k);
void listput(double val, int k);
void setuplist();

/****** LOCAL FUNCTIONS *********/

/* functions for linkedlist */
/* simple linkedlist,  with the biggest value at the right end, where last always points */
void addfirst(double val, int k)
{
  list[k]->v = val;
  list[k]->r = first;
  list[k]->l = NULL;
  first->l = list[k];
  first = list[k];
}

void addlast(double val, int k)
{
  list[k]->v = val;
  list[k]->r = NULL;
  list[k]->l = last;
  last->r = list[k];
  last = list[k];
}

void addright(double val, int k)
{
  list[k]->v = val;
  list[k]->r = now->r;
  list[k]->l = now;
  now->r->l = list[k];
  now->r = list[k];
  now = list[k];
}

void addleft(double val, int k)
{
  list[k]->v = val;
  list[k]->l = now->l;
  list[k]->r = now;
  now->l->r = list[k];
  now->l = list[k];
  now = list[k];
}

void listput(double val, int k)
{
  if (val <= now->v)
  {
    while (val <= now->v && now != first)
      now = now->l;
    if (now == first && val <= now->v)
      addfirst(val,k);
    else
      addright(val,k);
  }
  else
  {
    while (val > now->v && now != last)
      now = now->r;
    if (now == last && val > now->v)
      addlast(val,k);
    else
      addleft (val, k);
  }
} // listput

/* make enough of these listelements at the beginning rather than malloc and free them as we go along */
void setuplist()
{
  int i;
  list = (struct listelement **) malloc(genealogysamples * sizeof (struct listelement *));
  for (i=0;i<genealogysamples;i++)
    list[i] = (struct listelement *) malloc(sizeof (struct listelement));
}

/* for formatting a floating point number */
char *logpstrformat (double pval)
{
  logpfstr[0] = '\0';
  if (fabs (pval) < 1e-2)
    sprintf (logpfstr, "%.6lf", pval);
  else if (fabs (pval) < 1e-1)
    sprintf (logpfstr, "%.5lf", pval);
  else if (fabs (pval) < 1e-0)
    sprintf (logpfstr, "%.4lf", pval);
  else if (fabs (pval) < 1e1)
    sprintf (logpfstr, "%.3lf", pval);
  else if (fabs (pval) < 1e2)
    sprintf (logpfstr, "%.2lf", pval);
  else if (fabs (pval) < 1e3)
    sprintf (logpfstr, "%.1lf", pval);
  else if (fabs (pval) < 1e4)
    sprintf (logpfstr, "%.0lf", pval);
  return &logpfstr[0];
}

int fillparamlist(char *c, int *paramlist, int startparamcount, int notused[])
{
	int i;
	i = 0;
	if (isspace(*c))
		c = nextnonspaceafterspace(c);
	while (c && sscanf (c, "%d",&paramlist[i]) >0)
	{
      paramlist[i] += startparamcount;
      notused[paramlist[i]] = 1;
	  c = nextnonspaceafterspace(c);
	  i++;
	}
	paramlist[i] = -1;
	return  i;
}  // fill paramlist

int reduce_mappos(int mappos[], int k, int len)
{
	int i=0, j;
	while (i < nparams && mappos[i] != k) i++;
	assert(i < nparams);
	for (j=i;j< nparams;j++)
		mappos[j] = mappos[j+1];
	return len - 1;
}

void fillxmap(double xmap[], int pos[], int pl[][MAXPLENGTH], double terms[], int len, int type[], int num)
{
  int i, j;
  /* put positions into xmap of parameters that are still in the model  */
  for (i=0;i<nparams;i++)
  {
    j = 0;
    /*CR 110912.1 change compound test order to eliminate invalid mem access */
    while ( j < len && pos[j] != i ) j++;
    if (j<len)
      xmap[i] = (double) j;
  }
  /* now put constants and other positions into xmap  (use negative values for consants ) */
  for (i=0;i<num;i++)
  {
    if (type[i])
    {
      j = 0;
      while (pl[i][j] != -1)
      {
        xmap[pl[i][j]] = -terms[i];  // negative of constant is stored
        j++;
      }
    }
    else
    {
      j = 0;
      while (pl[i][j] != -1)
      {
        xmap[pl[i][j]] = (double) (int) xmap[(int) terms[i]];
        j++;
      }
    }
  }
} //fillxmap

/* reads the model file and sets up the arrays that determine how a full model maps onto a nested model */
/* really horrible code in this function */
void setup_mapping(char *fname)// JH fixed signficiant bug in this functino 4/5/2010
{
  FILE *nestedmodelfile;
  char *modeltextline;
  char keyword[12];
  char *c, m_or_p;
  int  mappos[MAXPLENGTH], mappos_len;
  int paramlist[MAXPTERMS][MAXPLENGTH];
  double pterms[MAXPTERMS]; // first number given on a constant or equal line,  if constant line it is a constant, else it is the first parameter number
  int ptype[MAXPTERMS];
  int i, n,numm, inmodel, thismodeltype;

  if ((nestedmodelfile = fopen (fname, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening nested model file: %s", fname);
  }
  modeltextline = static_cast<char *> (malloc(MAXMODELTEXTLINE*sizeof(char)));
  do
  {
    fgetval = fgets(modeltextline,MAXMODELTEXTLINE,nestedmodelfile);
  } while (!isdigit(modeltextline[0]));
  scanfval = sscanf (modeltextline, "%d", &num_nestedmodels);
  inmodel = -1;
  numm = -1;   // JH fixed bug not sure if it caused problems
  while (fgets(modeltextline,MAXMODELTEXTLINE,nestedmodelfile)!= NULL) if (isalpha(modeltextline[0]))
  {
    for (i=0;i< (int) strlen(modeltextline);i++) // make sure text beginning line is lowercase
    {
      if (isspace(modeltextline[i]))   break;
      modeltextline[i] = tolower(modeltextline[i]);
    }
    c = modeltextline;
    scanfval = sscanf (c, "%s", keyword);
    c = nextnonspaceafterspace(c);
    if (strcmp(keyword,"model")==0)
    {
      if (inmodel >= 0) // finish up the previous model
      {
        if (numm < 1)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d does not include any 'constant' or 'equal' specifications",inmodel);
        if (npops ==2 && fullmodeltype[inmodel]!= 0)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has m or p type but should be both when there are only two sampled populations",inmodel);
        if (npops > 2 && fullmodeltype[inmodel]== 0)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has type 0 but should be m or p type",inmodel);
        modeltypecounts[fullmodeltype[inmodel]]++;
        fillxmap(holdxmaps[inmodel], mappos,paramlist,pterms,mappos_len,/*nparams_nested[inmodel],*/ptype,numm);
      }
      inmodel++;
      strcpy(modelnamelines[inmodel],c);
      if (strlen(modelnamelines[inmodel]) && modelnamelines[inmodel][strlen(modelnamelines[inmodel])-1] == '\n')
      {
        modelnamelines[inmodel][strlen(modelnamelines[inmodel])-1] = 0;
      }
      numm=0;
      numpopsizeparams_nested[inmodel] =numpopsizeparams;
      nummigrateparams_nested[inmodel] = nummigrateparams;
      mappos_len = nparams;
      for (i=0;i<nparams;i++)
      {
	    mappos[i] = i;
	    holdxmaps[inmodel][i] = -1;
      }
      if (npops == 2)
      {
        fullmodeltype[inmodel] = 0;
        thismodeltype = 0;  //JH  fixed a bug 5/10/2010 that popped up under cygwin/linuxv
        nparams_nested[inmodel] = nparams;
      }
      else
      {
        thismodeltype = -1;
      }
    }
    if (strcmp(keyword,"equal")==0 || (strcmp(keyword,"constant")==0))
    {
      scanfval = sscanf (c, "%c",&m_or_p);
      if (thismodeltype == -1)
      {
        if (m_or_p == 'm')
        {
          thismodeltype = 2; //jh fixed bug 5/27/2010
          fullmodeltype[inmodel] = 2;
          nparams_nested[inmodel] = nummigrateparams;
        }
        if (m_or_p == 'p')
        {
          thismodeltype = 1; //jh fixed bug 5/27/2010
          fullmodeltype[inmodel] = 1;
          nparams_nested[inmodel] = numpopsizeparams;
        }
      }
      else
      {
        if (m_or_p == 'm' && thismodeltype == 1)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d specified as both p and m types",inmodel);
        if (m_or_p == 'p' && thismodeltype == 2)
            IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d specified as both p and m types",inmodel);
      }
      c = nextnonspaceafterspace(c);
      scanfval = sscanf (c, "%lf",&pterms[numm]); // first #, can be either a constant or a parameter number
      if (m_or_p == 'm' && strcmp(keyword,"equal")==0) //then first element of pterms is a parameter number, else its a constant and not changed
        pterms[numm] += numpopsizeparams;
      c = nextwhite(c);
      if (strcmp(keyword,"equal")==0) // found a match
      {
        ptype[numm] = 0;
      }
      else  // "constant"
      {
        ptype[numm] = 1;
        if (pterms[numm] <= MINPARAMVAL )
        {
	        pterms[numm] = (double) MINPARAMVAL;
            nestedmodel_boundary_in_model[inmodel] =1 ;
        }
        else
        {
          if (!calcoptions[LOADPRIORSFROMFILE]) //unlikely situation of setting constant values greater than the prior, but having mig priors read from file
          {
            if (m_or_p == 'm' && pterms[numm] >= mprior)
            {
              pterms[numm] = mprior;
              nestedmodel_boundary_in_model[inmodel] = 1;
            }
            if (m_or_p == 'p' && pterms[numm] >= thetaprior)
            {
              pterms[numm] = thetaprior;
              nestedmodel_boundary_in_model[inmodel] = 1;
            }
          }
        }
        atleastone_boundary_in_model |= nestedmodel_boundary_in_model[inmodel];
      }
      if (m_or_p == 'p')
        n = fillparamlist(c,&paramlist[numm][0], 0, paramsnotused[inmodel]);
      if (m_or_p == 'm')
        n = fillparamlist(c,&paramlist[numm][0],numpopsizeparams,paramsnotused[inmodel]);
      for (i=0;i<n;i++)
      {
        if (m_or_p == 'p')
          numpopsizeparams_nested[inmodel]--;
        else
          nummigrateparams_nested[inmodel]--;
        nparams_nested[inmodel] =  reduce_mappos(mappos,paramlist[numm][i],nparams_nested[inmodel]);
        mappos_len--;
      }
      numm++;
    }
  }
  if (inmodel >= 0) // finish up the last model
  {
    if (numm < 1)
      IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d does not include any 'constant' or 'equal' specifications",inmodel);
    if (npops ==2 && fullmodeltype[inmodel]!= 0)
        IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has m or p type but should be both",inmodel);
    if (npops > 2 && fullmodeltype[inmodel]== 0)
        IM_err (IMERR_NESTEDMODELLSPECIFYLFAIL,"nested model %d has type 0 but should be m or p type",inmodel);
    modeltypecounts[fullmodeltype[inmodel]]++;
    fillxmap(holdxmaps[inmodel], mappos,paramlist,pterms,mappos_len,/*nparams_nested[inmodel],*/ptype,numm);
  }

  XFREE(modeltextline);
  FCLOSE(nestedmodelfile);
}// setup_mapping;

//from a full parameter set down to a nested one
void reversemapvals(double *from, double *to, double *xmap)
{
  int i, toi;

  /*  */
  for (i = 0,toi = 0;i<nparams;i++)
  {
    if (toi == (int) xmap[i])
    {
      to[toi] = from[i];
      toi++;
    }
  }
} //reversemapvals

//from a nested parameter set to a full one
void mapvals(double *from, double *to, double *xmap)
{
  int i = 0;
  for (i=0;i<nparams;i++)
  {
    if (xmap[i] < 0)
      to[i] = -xmap[i];
    else
      to[i] = from[(int) xmap[i]];
  }
} //mapvals */

void setbounds()
{
  int i;
  parameter_lower_bound = static_cast<double *>
                        (malloc(nparams * sizeof(double)));
  parameter_upper_bound = static_cast<double *>
                        (malloc(nparams * sizeof(double)));
  for (i = 0; i < nparams; ++i)
  {
    parameter_lower_bound[i] =  MINPARAMVAL;
    if (i < numpopsizeparams)
    {
      parameter_upper_bound[i] = C[ARBCHAIN]->itheta[i].pr.max;
    }
    else
    {
      parameter_upper_bound[i] = C[ARBCHAIN]->imig[i-numpopsizeparams].pr.max;
    }
  }
} /* setbounds */


/* applies the differential evolution algorithm */
void nextgen(double *lowpd, double *hipd, int modelparams,int modelnum)
{
  int A,B,Cvalue, i,j, jmax;
  double temp;
  double u;

  *lowpd = STARTHI;
  *hipd = STARTLOW;

  jmax = nparamrange[0] + modelparams;
  assert(jmax <= nparamrange[1]);
  /*   CR 110929.2 from JH  9/26/2011 note that for nested models this may
   *   not use the right bounds, if bounds vary among parameters in
   *   a class - this is abug
   */
  for (i=0;i<depopsize;i++)
  {
    A = randposint(depopsize);
    B = randposint(depopsize);
    Cvalue = randposint(depopsize);
    //for (j=nparamrange[0];j<nparamrange[1];j++)
    for (j=nparamrange[0];j<jmax;j++)
    {
      u = uniform();
      if (u < recrate)
      {
      temp = pop[Cvalue][j] + fweight* (pop[A][j] - pop[B][j]);
      if (temp < parameter_lower_bound_mapped[j]) // move only part of the way towards the lower bound
        temp = pop[Cvalue][j] -  uniform() * (pop[Cvalue][j] - parameter_lower_bound_mapped[j]);
      assert ( temp >= parameter_lower_bound_mapped[j]);
      if (temp > parameter_upper_bound_mapped[j])// move only part of the way towards the upper bound
        temp = pop[Cvalue][j] +  uniform() * (parameter_upper_bound_mapped[j]- pop[Cvalue][j]);
      assert(temp <= parameter_upper_bound_mapped[j]);
      trialpop[i][j] = temp;
      }
      else
        trialpop[i][j] = pop[i][j];
    }
    if (modelnum >= 0)
    {
      mapvals(trialpop[i],popnest,holdxmaps[modelnum]);
      trialpop[i][nparams] = jointp(popnest, calculate_ess,&effective_sample_size);
    }
    else
      trialpop[i][nparams] = jointp(trialpop[i], calculate_ess,&effective_sample_size);
  }
  for (i=0;i<depopsize;i++)
  {
    if (trialpop[i][nparams] < pop[i][nparams])
      for (j=0;j<=nparams;j++)
        pop[i][j] = trialpop[i][j];
    if (pop[i][nparams] < *lowpd)
      *lowpd = pop[i][nparams];
    if (pop[i][nparams] > *hipd)
       *hipd = pop[i][nparams];
  }

} //nextgen

/* create a starting population by sampling uniformally over the priors */
void startpop(int modelparams,int modelnum)
{
  int i,j, jmax;
  jmax = nparamrange[0] + modelparams;
  if (modelnum < 0) for (j=0;j<nparams;j++)
  {
    parameter_lower_bound_mapped[j]=parameter_lower_bound[j];
    parameter_upper_bound_mapped[j]=parameter_upper_bound[j];
  }
  else
  {
    reversemapvals(parameter_lower_bound,parameter_lower_bound_mapped,holdxmaps[modelnum]);
    reversemapvals(parameter_upper_bound,parameter_upper_bound_mapped,holdxmaps[modelnum]);
  }

  for (i=0;i< depopsize;i++)
  {
    for (j=0;j<nparams;j++)
    {
                                /* CR 110929.2 from JH 9/26/2011 this
                                 * turned off because  parameter vals will be
                                 * mapped into popnest */
      //if (j>=nparamrange[0] && j<nparamrange[1])
      if (j>=nparamrange[0] && j<jmax)   /* CR 110929.2 from JH 9/26/2011 but
                                          * for nested models this may not
                                          * use the right bounds, if bounds
                                          * vary among parameters in a
                                          * class - this is abug */
        pop[i][j] = parameter_lower_bound_mapped[j] + uniform()*(parameter_upper_bound_mapped[j] - parameter_lower_bound_mapped[j]);
      else
        pop[i][j] = -1;
    }
    if (modelnum >= 0)
    {
      mapvals(pop[i],popnest,holdxmaps[modelnum]);
      pop[i][nparams] = jointp(popnest, calculate_ess,&effective_sample_size);
    }
    else
      pop[i][nparams] = jointp(pop[i], calculate_ess,&effective_sample_size);
  }
}// startpop

/* keeps calling nextgen() until a peak has been found */
double difeloop(int modelparams, int modelnum)
{
  double spread, lowpd, hipd;
  /* static double pd; Never used */

  do
  {
    nextgen(&lowpd, &hipd,modelparams,modelnum);
    spread = hipd - lowpd;
  }
  while (spread > SPREADTOL);
  return  lowpd;
} // difeloop

/* find the best individual in the population and copy it into bestvals[] */
void copybest(/* int modelparams */)
{
  int i,k;
  double temp = STARTHI;
  for (i=0;i<depopsize;i++)
  {
    if (pop[i][nparams] < temp)
    {
      temp = pop[i][nparams];
      k = i;
    }
  }
  for (i=0;i<=nparams;i++)
    bestvals[i] = pop[k][i];
  bestvals[nparams] = temp;
}//copybest

/*
  find a peak,  record location,  then begin a loop that repeats the search
 - when a better peak is found it is saved
 - if the same best peak has been found more than LOOPMATCHCRITERIA times
 - or if restarts have been done MAXRESTART times,
    it finishes.
*/
void modelloop(int modelnum)
{
  int i,newstartloop;
  int modelparams;
  int countloop;
  double local_pd;
  double global_pd = STARTHI;
  double jointProb;  /* CR 110929.5  variable name changed */

  if (modelnum < 0)
  {
      switch (nowmodeltype)
      {
        case 0 : modelparams =  nparams; break;
        case 1 : modelparams =  numpopsizeparams;break;
        case 2 : modelparams = nummigrateparams;break;
      }
  }
  else
  {
    modelparams = nparams_nested[modelnum];
  }
  depopsize = modelparams * popsizemultiplier;
  startpop(modelparams,modelnum);
  for (i=0;i< depopsize;i++)
  {
    if (pop[i][nparams] < global_pd)
      global_pd = pop[i][nparams];
  }
  newstartloop = 0;
  countloop = 0;
  do
  {
    if (newstartloop > 0)
    {
      startpop(modelparams,modelnum);
      printf("         restart #%d (out of %d maximum)\n",newstartloop,MAXRESTART);
    }
    local_pd = difeloop(modelparams,modelnum);
    if (fabs(local_pd - global_pd) < PLOOPTOL)
      countloop++;
    else
    {
      if (local_pd < global_pd)
        countloop = 0;
    }
    if (local_pd < global_pd)
    {
      copybest(/* modelparams Does nothing */);
      global_pd  = local_pd;
    }
    newstartloop ++;
  }
  while (countloop < LOOPMATCHCRITERIA && newstartloop < MAXRESTART);
  calculate_ess = 1;

  /*  CR 110929.5  although return val from jointP() call is not used
   *  locally, it may be useful during debuggin, also the return
   *  variable name was changed
   */
  if (modelnum >= 0)
  {
     mapvals(bestvals,popnest,holdxmaps[modelnum]);
     jointProb = jointp(popnest,calculate_ess, &effective_sample_size);
  }
  else
     jointProb = jointp(bestvals, calculate_ess,&effective_sample_size);
  calculate_ess = 0;
} // modelloop

void printjointpeakvals(FILE *outfile,int notused[], double *printvals)
{
  int i;
  for (i=0;i<nparams;i++)
  {
    if (i>= nparamrange[0] && i < nparamrange[1])
    {
      if (notused[i])
      {
        if (printvals[i] < 0.001)
          FP"\t[%.5lf]",printvals[i]);
        else
          FP"\t[%.4lf]",printvals[i]);
      }
      else
      {
        if (printvals[i] < 0.001)
          FP"\t%.5lf",printvals[i]);
        else
          FP"\t%.4lf",printvals[i]);
      }
    }
    else
      FP"\t-");
  }
  FP"\n");
} //printjointpeakvals

void startjointpeakouttable(FILE *outfile,char *fname)
{
  int i, mi, typeloopstart, typeloopstop,nowmodeltype_local,j=0;
  FP "%s",outputbanner("Joint Peak Locations and Posterior Probabilities"));

  FP"  estimates based on %d sampled genealogies\n",genealogysamples);
  if (strlen(fname) > 0)
    FP"  nested model filename:%s\n",fname);
  FP"\nModel#  Model Description\n");
  if (npops == 2)
  {
    typeloopstart = 0;
    typeloopstop = 0;
  }
  else
  {
    typeloopstart = 1;
    typeloopstop = 2;
  }
  for (nowmodeltype_local = typeloopstart; nowmodeltype_local <= typeloopstop;nowmodeltype_local++)
  {
    j++;
    FP"%2d     %s\n",j,modelstartstr[nowmodeltype_local]);
    for (mi = 0;mi<num_nestedmodels;mi++) if (fullmodeltype[mi] == nowmodeltype_local)
    {
      j++;
      FP"%2d   %s\n",j,modelnamelines[mi]);
    }
  }

  FP"\nModel#\tlog(P)\t#terms\tdf\t2LLR\tESS");
  for (i = 0; i < numpopsizeparams; i++)
    FP"\t%s",C[ARBCHAIN]->itheta[i].str);
  for (i = 0; i < nummigrateparams; i++)
    if (C[ARBCHAIN]->imig[i].pr.max > MPRIORMIN)
      FP"\t%s",C[ARBCHAIN]->imig[i].str);
  FP"\n");
}

/* JH new jointp brought in 9/16/09 */
#define OCUTOFF  10
#define POW10I(a) ((a)+308)
#define PRANGELOG 10
double jointp (double *x, int calc_ess,double *effective_n)
/* calculate the joint likelihhood function  */
/* the point at which to calculate the value of the function is in x */
/* return the negative of the logarithm of the joint posterior probability at x */
/* if *effective_n == 1.0  calculate the effective number of samples */
/* this has been tweaked a lot for speed
  only saves the best genealogies (uses a linked list)
  precalculates as much stuff as possible that will get reused later */
{
  int gi, i, i1;
  int gin, iin, ii;
  double sum, p, q;
  static int ccp, fcp, hccp, mcp, fmp, qip, mip, probgp, pdgp;
  static double log_numtrees;
  static int init = 0;
  int zadj, maxz = -10000000;
  double acumm = 0;
  float *g;
  double *dg;
  double acumm_sqr = 0;


  static double pow10[617];

  if (init == 0)
  {
    init = 1;
    log_numtrees = log ((double) genealogysamples);
    /********* this section copied from IMa2  initialize.c setup_iparams()  on 8/24/09 */
    ccp = 0;
    fcp = ccp + numpopsizeparams;
    hccp = fcp + numpopsizeparams;
    mcp = hccp + numpopsizeparams;
    fmp = mcp + nummigrateparams;
    qip = fmp + nummigrateparams;
    mip = qip + numpopsizeparams;
    pdgp = mip + nummigrateparams;
    /* CR 110830.1
     * probgp was not correct, the old variable plgp in jointp function
     * (jointfind.c, line 912) was no longer used.
     */
    probgp = pdgp + 1;
    /**********/
    /* pow10[],logx, divx, log2diffx and doubling gsampinf[gi][fcp+i] are done to speed things up*/
    for (i=-308;i<=308;i++)
      pow10[POW10I(i)] = pow(static_cast<double>(10),i);
    logx = (double *) malloc(nparams * sizeof (double));
    divx = (double *) malloc(nparams * sizeof (double));
    log2diffx = (double *) malloc(numpopsizeparams * sizeof (double));
    double_gsamp_fcp = alt2d_alloc2Ddouble(genealogysamples, numpopsizeparams);
    if (npops > 2)
    {
      probgp_popsize = (double *) malloc(genealogysamples * sizeof(double));
      probgp_migrate = (double *) malloc(genealogysamples * sizeof(double));
    }
    for (gi = 0; gi < genealogysamples; gi++)
    {
	  for (i = 0; i < numpopsizeparams; i++)
        double_gsamp_fcp[gi][i] = 2.0 * gsampinf[gi][fcp + i];
      if (npops > 2)
      {
        for (i = 0,probgp_popsize[gi] = 0; i < numpopsizeparams; i++)
          probgp_popsize[gi] += gsampinf[gi][qip + i];
        for (i = 0,probgp_migrate[gi] = 0; i < nummigrateparams; i++)
          probgp_migrate[gi] += gsampinf[gi][mip + i];
      }
    }
  }
  sum = 0;
  /* CR 110929.2 from JH 9/26/2011 fixed a bug  the initialization of
   * log2diffx was being done even when no population size parameters
   * are in the model being optimized (i.e. nowmodeltype==2)
   */
  if (nowmodeltype == 0 || nowmodeltype == 1)
  for (i=0;i<numpopsizeparams;i++)
    log2diffx[i] = LOG2 - log (x[i]);
  /* CR 110929.2 from JH 9/26/2011 fixed a bug  logx and divx were being
   * set for parameters that are not in the model being optimized
   * (i.e. nowmodeltype==2)
   */
  //for (i=0;i<nparams;i++)
  for (i=nparamrange[0];i< nparamrange[1];i++)
  {
    logx[i] = log(x[i]);
    divx[i] = 1.0/x[i];
  }
  for (gi = 0, ii=0; gi < genealogysamples; gi++)
  {
    g = gsampinf[gi];
    dg = double_gsamp_fcp[gi];
    if (npops == 2)
      p =  - g[probgp];
    else
    {
      if (nowmodeltype == 1)
        p = -probgp_popsize[gi];
      else
        p = -probgp_migrate[gi];
    }
    for (i = nparamrange[0]; i < nparamrange[1]; i++)
    {
      if (i< numpopsizeparams)
        p +=	g[ccp + i] *log2diffx[i] - g[hccp + i] -  dg[i] * divx[i];
      else
      {
        i1 = i - numpopsizeparams;
        p += g[mcp + i1] * logx[i] - g[fmp +i1] * x[i];
      }

    }
    /* add the value to the linked list,  but only if it is within PRANGELOG of the biggest value found so far */
    if (gi == 0)
    {
      list[0]->v = p;
      list[0]->r = list[0]->l = NULL;
      last = first = now = list[0];
    }
    else
    {
      if (last->v - p < PRANGELOG) // only save values that are close to the biggest found so far
      {
        ii++;
        listput(p, ii);
      }
    }
  }
  iin = ii;
  gi = 0;
  now = last;
  do
  {
    eexp (now->v, &eexpsum[gi].m, &eexpsum[gi].z);
    if (eexpsum[gi].z > maxz)
      maxz = eexpsum[gi].z;
    now = now->l;
    gi++;
  }
  while (gi < iin && last->v - now->v < PRANGELOG);
  gin = gi;
  maxz -= OCUTOFF;
  for (gi = 0; gi < gin /*genealogysamples */; gi++)
    {
	  zadj = eexpsum[gi].z - maxz;
      if (zadj > -308 && zadj < 308)
	    eexpsum[gi].m *= pow10[POW10I(zadj)];
      else
      {
        if (zadj <= -308)
          eexpsum[gi].m = 0.0;
        else
          eexpsum[gi].m = DBL_MAX;//trap overflow, shouldn't get here in this program
      }
	  acumm += eexpsum[gi].m;
      if (calc_ess)
        acumm_sqr += eexpsum[gi].m * eexpsum[gi].m;
    }
    if (calc_ess)
      *effective_n = acumm*acumm/acumm_sqr;
    sum = log (acumm) + maxz * LOG10;
    q = log_numtrees - sum;  // negative of log of posterior probability

    return (double) q;
} /* jointp*/

#undef OCUTOFF

void freejointpmem()
{
  int i;
  alt2d_free2D(trialpop);
  free(bestvals);
  free(tempvals);
  free(parameter_lower_bound);
  free(parameter_upper_bound);
  free(logx);
  free(divx);
  free(log2diffx);
  alt2d_free2D(double_gsamp_fcp);
  alt2d_free2D(pop);
  if (npops > 2)
  {
    free(probgp_popsize);
    free(probgp_migrate);
  }
  free(popnest);
  for (i=0;i<genealogysamples;i++)
    free(list[i]);
  free(list);
}

/***********GLOBAL FUNCTIONS **********/

/* this function can only be called for 2- and 3- population analyses
  for 3 population analyses there are two different 'full' models
    - a model with all population size parameters
    - a model with all migration parameters
  unfortunately cannot do a 3 pop model with both population size and migration parameters
  just too many parameters */
/* JH  renamed holdparams  to largestmodel_nparams */
/* CR 110921.1  Change type declaration of outfile parameter so that this
 * function can pass in the proper parameter type when calling closeopenout()
 */

void findjointpeaks(FILE **outfile,char *outfilename, char *nestfname,int number_of_parameters)
{
  int mi,j=0, typeloopstart, typeloopstop;
  double holdml;
  int largestmodel_nparams;

  nparams = number_of_parameters;
  bestvals = static_cast<double *> (malloc((nparams+1) * sizeof(double)));
  setbounds();
  setuplist();
  depopsize = nparams * popsizemultiplier;
  pop = alt2d_alloc2Ddouble(depopsize, nparams+1);
  trialpop = alt2d_alloc2Ddouble(depopsize, nparams+1);
  popnest = (double *) malloc(nparams * sizeof(double));
  if (strlen(nestfname)> 0)
    setup_mapping(nestfname);
  else
    num_nestedmodels = 0;
  /* CR 110921.1  Change type of outfile parameter */
  startjointpeakouttable(*outfile, nestfname);
  if (npops == 2)
  {
    typeloopstart = 0;
    typeloopstop = 0;
    nparamrange[0] = 0;
    nparamrange[1] = nparams;
  }
  else
  {
    typeloopstart = 1;
    typeloopstop = 2;
  }
  for (nowmodeltype = typeloopstart; nowmodeltype <= typeloopstop;nowmodeltype++)
  {
    j++;
   // if (modeltypecounts[nowmodeltype] > 0)
    {
      printf(" starting model: %s\n", modelstartstr[nowmodeltype]);
      switch (nowmodeltype)
      {
        case 0 :
              largestmodel_nparams = nparams;
              nparamrange[0] = 0;
              nparamrange[1] = nparams;
              break;
        case 1 :
              largestmodel_nparams = numpopsizeparams;
              nparamrange[0] = 0;
              nparamrange[1] = numpopsizeparams;
          break;
        case 2 :
              largestmodel_nparams = nummigrateparams;
              nparamrange[0] = numpopsizeparams;
              nparamrange[1] = numpopsizeparams+nummigrateparams;
          break;
      }
      modelloop(-1);
      holdml = bestvals[nparams];
      /* CR 110921.1  Change type of outfile parameter */
      fprintf(*outfile, "%d\t%s\t%d\t-\t-",j,
                         logpstrformat(-bestvals[nparams]), largestmodel_nparams);
      fprintf(*outfile, "\t%s",logpstrformat(effective_sample_size));
      printjointpeakvals(*outfile,paramsallused, bestvals);
//should this be &outfile ??
      closeopenout (outfile, outfilename);
      printf ("done model : %s\n", modelstartstr[nowmodeltype]);
      printf("joint density: %f\n",-bestvals[nparams]);
    }
    for (mi = 0;mi<num_nestedmodels;mi++) if (fullmodeltype[mi] == nowmodeltype)
    {
      j++;
      printf(" starting model: %s\n", modelnamelines[mi]);// to stdout
      modelloop(mi);
      mapvals(bestvals,popnest,holdxmaps[mi]);
      fprintf(*outfile, "%d\t%s\t%d",j,logpstrformat(-bestvals[nparams]),
                         nparams_nested[mi]);
      fprintf(*outfile,"\t%d",largestmodel_nparams-nparams_nested[mi]);
      if (nestedmodel_boundary_in_model[mi])
        fprintf(*outfile,"*");
      fprintf(*outfile, "\t%s", logpstrformat(2*(bestvals[nparams]-holdml)));
      fprintf(*outfile, "\t%s",logpstrformat(effective_sample_size));
      printjointpeakvals(*outfile,paramsnotused[mi], popnest);
      closeopenout (outfile, outfilename);
      // print some info to stdout
      printf ("done model : %s\n", modelnamelines[mi]);
      printf("joint density: %f\n",-bestvals[nparams]);
    }
  }
  if (atleastone_boundary_in_model)
    fprintf(*outfile,"    * test distribution of 2LLR is a mixture\n");
  fprintf(*outfile,"\n");
  freejointpmem();
}
