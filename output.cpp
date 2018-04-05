/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

#undef GLOBVARS
#include "ima.hpp"

extern int numtreesarray[];   /* number of possible ordered trees,  for up to 7 populations */


/*********** LOCAL STUFF **********/

static double calcx (int ei, int pnum, int mode);
extern int numuniquenodes[]; // in alltreestrings.cpp
void fill_2Nm_vals (struct plotpoint **popmigxy, char  **popmigstr);


/******* LOCAL FUNCTIONS ***********/

/* used when calculating means variances and correlations */
double
calcx (int ei, int pnum, int mode)
{
  int cc, mc,p;
  double fc, fm, hval, denom, max;
  double tempval;

  p = pnum;
  if (p < numpopsizeparams)
  {
    if (C[ARBCHAIN]->itheta[p].pr.max <= MINPARAMVAL)
      return -1;
    cc = (int) gsampinf[ei][gsamp_ccp + p];
    fc = gsampinf[ei][gsamp_fcp + p];
    assert(fc >= 0.0);
    hval = gsampinf[ei][gsamp_hccp + p];
    denom = gsampinf[ei][gsamp_qip + p];
    max = C[ARBCHAIN]->itheta[p].pr.max;
    if (mode == 0)
    {
      if (cc == 0 && fc == 0)
      {
        tempval = (SQR (max) / 2) / exp (denom);
      }
      else
      {
        if (cc > 1)
        {
          tempval = exp (2 * LOG2 - hval + (2 - cc) * log (fc) + uppergamma (cc - 2, 2 * fc / max) - denom);
        }
        else
        {
          if (cc == 1)
          {
            tempval = exp (LOG2 - hval + log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
          else                  //cc==0
          {
            assert (cc == 0);
            tempval = exp (log ((max / 2) * (max - 2 * fc) * exp (-2 * fc / max) + 2 * SQR (fc) * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
        }
      }
    }
    else                        //mode == 1
    {
      if (cc == 0 && fc == 0)
        tempval = (max * SQR (max) / 3) / exp (denom);
      else
      {
        if (cc > 2)
        {
          tempval = exp (uppergamma (cc - 3, 2 * fc / max) + 3 * LOG2 - hval + (3 - cc) * log (fc) - denom);
        }
        else
        {
          if (cc == 2)
          {
            tempval = exp (2 * LOG2 - hval + log (max * exp (-2 * fc / max) - 2 * fc * exp (uppergamma (0, 2 * fc / max))) - denom);
          }
          else
          {
            if (cc == 1)
            {
              tempval = exp (-hval + log (max * (max - 2 * fc) * exp (-2 * fc / max) + 4 * SQR (fc) * exp (uppergamma (0, 2 * fc / max))) - denom);
            }
            else                //cc==0
            {
              assert (cc == 0);
              tempval = exp (-log (3.0) + log (max * (2 * SQR (fc) - fc * max + SQR (max)) * exp (-2 * fc / max) - 4 * pow ((double) fc, 3.0) * exp (uppergamma (0, 2 * fc / max))) - denom);
            }
          }
        }
      }
    }
  }
  else
  {
    p -= numpopsizeparams;
    if (C[ARBCHAIN]->imig[p].pr.max <= MINPARAMVAL)
      return -1;
    mc = (int) gsampinf[ei][gsamp_mcp + p];
    fm = gsampinf[ei][gsamp_fmp + p];
    assert(fm >= 0);
    if (fm < 0.0)
      {
        IM_errloc (AT, "fm cannot be negative.");
      }
    denom = gsampinf[ei][gsamp_mip + p];
    max = C[ARBCHAIN]->imig[p].pr.max;
    if (mode == 0)
    {
      if (mc == 0 && fm == 0)
      {
        tempval = (SQR (max) / 2) / exp (denom);
      }
      else
      {
        if (mc > 0)
          tempval = exp (lowergamma ((int) mc + 2, fm * max) - (mc + 2) * log (fm) - denom);
        else
          tempval = (1 - (1 + fm * max) * exp (-fm * max)) / SQR (fm) / exp (denom);
      }
    }
    else                      //mode == 1
    {
      if (mc == 0 && fm == 0)
      {
        tempval = (pow (max, 3.0) / 3) / exp (denom);
      }
      else
      {
        tempval = exp (lowergamma ((int) mc + 3, fm * max) - (mc + 3) * log (fm) - denom);
      }
    }
    
  }
  /* Jh changed this from 0 to 0.0, some compilers return a false because tempval is a double */ 
  assert (tempval >= 0.0);         // it should be greater than or equal to 0   (not sure really,  is it ok to be 0 ? 
  return tempval;
}                               //calcx

/* fill strings and xy values for printing 2Nm ascii curves */
void fill_2Nm_vals (struct plotpoint **popmigxy, char  **popmigstr)
{
  int i, j, k, hpi, mpop, thetai, mi, found;
  double pmmax, tempy, tempx, maxxfind;

  hpi = 0;
  if (modeloptions[PARAMETERSBYPERIOD])
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        mpop = C[ARBCHAIN]->plist[k][i];
        thetai = 0;
        found = 0;
        while (!found && thetai < numpopsizeparams)
        {
          found = (k == atoi (&C[ARBCHAIN]->itheta[thetai].str[1])
                   && mpop == atoi (&C[ARBCHAIN]->itheta[thetai].str[3]));
          if (!found)
            thetai++;
        }
        for (mi = 0; mi < nummigrateparams; mi++)
        if ((modeloptions[EXPOMIGRATIONPRIOR]==1 && C[ARBCHAIN]->imig[mi].pr.expomean > MPRIORMIN) ||
          (modeloptions[EXPOMIGRATIONPRIOR]==0 && C[ARBCHAIN]->imig[mi].pr.max > MPRIORMIN))
        {
          found = 0;
          if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
          {
            found = (k == atoi (&C[ARBCHAIN]->imig[mi].str[1])
                     && (mpop == atoi (&C[ARBCHAIN]->imig[mi].str[3])
                         || mpop == atoi (&C[ARBCHAIN]->imig[mi].str[6])));
          }
          else
          {
            found = (k == atoi (&C[ARBCHAIN]->imig[mi].str[1])
                     && mpop == atoi (&C[ARBCHAIN]->imig[mi].str[3]));
          }
          if (found)
          {
            sprintf (popmigstr[hpi], "%d,2N%d", k, mpop);
            if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
            {
              strcat (popmigstr[hpi], "m");
              strcat (popmigstr[hpi], &C[ARBCHAIN]->imig[mi].str[3]);
            }
            else
            {
              strcat (popmigstr[i], C[ARBCHAIN]->imig[mi].str);
            }
            if (modeloptions[EXPOMIGRATIONPRIOR])
            {
              pmmax = EXPOMIGPLOTSCALE * C[ARBCHAIN]->imig[mi].pr.expomean;
            }
            else
            {
              pmmax = C[ARBCHAIN]->itheta[thetai].pr.max * C[ARBCHAIN]->imig[mi].pr.max / 2.0;
            }
            j = GRIDSIZE-1;
            maxxfind = 0;
            while (maxxfind==0 && j >= 0)// see if the maximum x value can be reduced from pmmax.  
            {
              tempx = (j + 0.5) * pmmax / GRIDSIZE;
              if (modeloptions[EXPOMIGRATIONPRIOR])
                tempy = calc_pop_expomig (thetai, mi,tempx , 0);
              else
                tempy = calc_popmig (thetai, mi,tempx , 0);
              if (tempy > 1e-6)  // 1e-6 fairly arbitrary cutoff for finding upper bound 
                maxxfind = tempx;
              else 
                j--;
            }
            if (j < 0)
              maxxfind = pmmax;
            for (j = 0; j < GRIDSIZE; j++)
            {
              popmigxy[hpi][j].x = (j + 0.5) * maxxfind / GRIDSIZE;
              if (modeloptions[EXPOMIGRATIONPRIOR])
              {
                popmigxy[hpi][j].y = calc_pop_expomig (thetai, mi, popmigxy[hpi][j].x, 0);
              }
              else
              {
                popmigxy[hpi][j].y = calc_popmig (thetai, mi, popmigxy[hpi][j].x, 0);
              }
            }
            hpi++;
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < numtreepops - 1; i++)
    {
      thetai = i;
      for (mi = 0; mi < nummigrateparams; mi++)
      if ((modeloptions[EXPOMIGRATIONPRIOR]==1 && C[ARBCHAIN]->imig[mi].pr.expomean > MPRIORMIN) ||
        (modeloptions[EXPOMIGRATIONPRIOR]==0 && C[ARBCHAIN]->imig[mi].pr.max > MPRIORMIN))
      {
        if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
        {
          found = ((thetai == atoi (&C[ARBCHAIN]->imig[mi].str[1])) || (thetai == atoi (&C[ARBCHAIN]->imig[mi].str[4])));
        }
        else
        {
          found = thetai == atoi (&C[ARBCHAIN]->imig[mi].str[1]);
        }
        if (found)
        {
          sprintf (popmigstr[hpi], "2N%d", thetai);
          strcat (popmigstr[hpi], C[ARBCHAIN]->imig[mi].str);
          if (modeloptions[EXPOMIGRATIONPRIOR])
          {
            pmmax = EXPOMIGPLOTSCALE * C[ARBCHAIN]->imig[mi].pr.expomean;
          }
          else
          {
            pmmax = C[ARBCHAIN]->itheta[thetai].pr.max * C[ARBCHAIN]->imig[mi].pr.max / 2.0;
          }
          j = GRIDSIZE-1;
          maxxfind = 0;
          while (maxxfind==0 && j >= 0)
          {
            tempx = (j + 0.5) * pmmax / GRIDSIZE;
            if (modeloptions[EXPOMIGRATIONPRIOR])
              tempy = calc_pop_expomig (thetai, mi,tempx , 0);
            else
              tempy = calc_popmig (thetai, mi,tempx , 0);
            if (tempy > 1e-9)
              maxxfind = tempx;
            else 
              j--;
          }
          if (j < 0)
            maxxfind = pmmax;
          for (j = 0; j < GRIDSIZE; j++)
          {
            popmigxy[hpi][j].x = (j + 0.5) * maxxfind / GRIDSIZE;
            if (modeloptions[EXPOMIGRATIONPRIOR])
              popmigxy[hpi][j].y = calc_pop_expomig (thetai, mi, popmigxy[hpi][j].x, 0);
            else
              popmigxy[hpi][j].y = calc_popmig (thetai, mi, popmigxy[hpi][j].x, 0);
          }
          hpi++;
        }
      }
    }
  }
}                               //void fill_2Nm_vals



/********** GLOBAL functions ********/

void closeopenout (FILE ** p_to_file, char fname[]) // why a pointer to a pointer? so the address can be passed back
{
  FCLOSE (* p_to_file);
  if ((* p_to_file = fopen (fname, "a+")) == NULL)
  {
    IM_err(IMERR_OUTPUTFILECHECK,"Error opening text file for writing");
  }
}
/*takes a pointer to a piece of text (length less than BANNERMAXLENGTH)
  puts a capitalized version in between two rows of equal signs '=' into the holder array,  bannerall 
  returns a pointer to the banner to be printed 
  bannerall[] is static so should be accessible to whatever calls it upon return 
*/
char* outputbanner(const char *bannertext)
{
  int i,l;
  static char bannerall[3*BANNERMAXLENGTH]; // needs to be static so it exists upon return so it can be printed (always happens immediately )
  char holdbanner[BANNERMAXLENGTH];
  strcpy(holdbanner,bannertext);
  char *b;
  l = strlen(holdbanner);
  assert (l < BANNERMAXLENGTH - 2);
  b = bannerall;
  *b = 0;
  strcat(b,"\n\n");
  for (i=0;i<l;i++)
    strcat(b,"=");
  strcat(b,"\n");
  for(char* c=&holdbanner[0],i = 0;i<l;c++,i++)
    *c=toupper(*c);
  strcat(b,holdbanner);
  strcat(b,"\n");
  for (i=0;i<l;i++)
    strcat(b,"=");
  strcat(b,"\n\n");
  return b;
}  // outputbanner


/*
 * This function is used for debugging only.  CR 110825.1
 checks to see if a file can be opened
 if already open, then error

 there was a problem that if the file does not exist then
 a file of size zero was being created
 so code was added 8/24/2016 to check the file size
 and delete it if it was zero
 */
void
checkoutfileclosed (FILE ** outfile, char outfilename[])
{
  long int fsize;

  if ((*outfile = fopen (outfilename, "a")) == NULL)
  {
    IM_err(IMERR_OUTPUTFILECHECK,"Error  - file for results output is not closed ");
  }
  fseek(*outfile,0,SEEK_END);
  fsize = ftell(*outfile);
  FCLOSE (*outfile);
  if (fsize==0)
    remove(outfilename);
}                               // checkoutfileclosed

/* prints to outfile it its not NULL, which should only be true if it is the head node - cpu 0
  also gathers topologycounts to cpu 0 from other cpus */ 
void
printrunbasics (FILE * outfile, int loadrun, char fpstr[], int burninsteps,int burninsteps_old,int runsteps_old, int mcmcrecords_old,int genealogysamples_old,
                int recordint, int mcmcrecords, int savegenealogyint,double hilike,double hiprob)
{
  //double mcalc;
  int totaltopol_rec, chain0topol_rec, chain0topolswaps_rec;
  double ml_estimate;//, ml_stdev1, ml_stdev2;
  if (calcoptions[CALCMARGINALLIKELIHOOD])  // make this call for all processors (i.e. regardless of outfile value) 
    move_calcmarglike_vals();
  if (outfile != NULL)
  {
 //printf("here 1 %ld",(long) strlen(fpstr));
  FP "%s", fpstr);
//printf("here 2");
  if (!loadrun == 1)
    {
      FP "%s",outputbanner("Sampling summaries"));
      FP "\nIntervals\n---------\n");
      FP "  Number of steps between recording values : %d\n", recordint);
      if (modeloptions[POPTREETOPOLOGYUPDATE]==0)
        FP "  Number of steps between saving genealogy information: %d\n", savegenealogyint);
      FP "\nBurnin Period Steps\n-------------------\n");
      if (burninsteps  > 0)
        FP "  Number of steps in burnin: %10d\n", burninsteps);
      if (burninsteps_old  > 0)
      {
        FP "  Number of steps in previous burnin: %10d\n", burninsteps_old);
        FP "  Total number of burnin steps: %10d\n", burninsteps_old+burninsteps);
      }
      FP "\nSampling Period Steps\n---------------------\n");
      if (runsteps > 0)
        FP "  Number of steps in chain: %10d \n", runsteps);
      if (runsteps_old > 0)
      {
        FP "  Number of steps in previous runs: %10d \n",runsteps_old);
        FP "  Total number of steps in all runs: %10d \n",runsteps_old + runsteps);
      }
      if (modeloptions[POPTREETOPOLOGYUPDATE]==0)
      {
        FP "\nGenealogy Sampling\n------------------\n");
       if (genealogysamples > 0)
        FP "  Number of genealogies sampled (per locus): %d \n", genealogysamples);
       if (genealogysamples_old > 0)
       {
          FP "  Number of genealogies sampled in prior runs: %d \n", genealogysamples_old);
          FP "  Total number of genealogies sampled: %d \n", genealogysamples_old + genealogysamples);
       }
      }
      FP "\nMCMC Sampling (other than genealogies)\n--------------------------------------\n");
      if (mcmcrecords > 0)
        FP "  Number of samples: %d \n", mcmcrecords);
      if (mcmcrecords_old > 0)
      {
        FP "  Number of samples from prior runs: %d \n", mcmcrecords_old);
        FP "  Total number of samples: %d \n", mcmcrecords_old+mcmcrecords);
      }
      FP "\nHighest Likelihood and Genealogy Probability (checked at time of sampling)\n--------------------------------------------------------------------------\n");
      FP "  Highest Joint P(G) (log): %10.3f\n", hiprob);
      FP "  Highest Joint P(D|G) (log): %10.3f \n", hilike);
      if (calcoptions[CALCMARGINALLIKELIHOOD])
      {
        thermo_marginlike_calc(mcmcrecords, &ml_estimate);
        FP "  Marginal Likelihood P(D) (log) estimate: %10.3lf\n",ml_estimate);
        /* 8/22/2016 turned off most marginal likelihood stuff as it was not working well when intervals between beta values were not uniform 
          need to get back to this to figure out how to calculate the variance of the estimate */
        /*FP "Marginal Likelihood P(D) (log) estimates\n");
        if (modeloptions[POPTREETOPOLOGYUPDATE]==0)
        {
          mcalc = harmonicmarginlikecalc(mcmcrecords);
          FP "     harmonic mean : %10.3f\n",mcalc);
        }
        thermo_marginlike_calc(mcmcrecords, &ml_estimate);
        //FP "     thermodynamic integration : %10.3f    discretization error : %10.3f   sampling error : %10.3f\n",ml_estimate,ml_stdev1,ml_stdev2); standard error not working it seems ??
        FP "     thermodynamic integration : %10.3f \n",ml_estimate);
        steppingstone_marginlike_calc2(mcmcrecords,&ml_estimate, &ml_stdev1);
        //FP "     steppingstone sampling (Aude Grelaud): %10.3f   sampling error : %10.3f\n",ml_estimate, ml_stdev1); standard error not working it seems ??
        FP "     steppingstone sampling: %10.3f\n",ml_estimate);
        */
      }
/*#ifdef DEBUG  5/15/2017 just not useful 
      // turn this on for debugging 8/18/2016,though rarely useful. only works when there is just one cpu, have not done mpi for this
      if (numprocesses == 1)
      {
        int li;
        FP "\n");
        FP "Highest P(D|G) (log) for each Locus \n");
        FP "\tLocus\tP(D|G)\tP(G|Theta)\n");
        for (li = 0; li < nloci; li++)
        {
          FP "\t%d\t%.3f\t%.3lf\n", li, C[0]->G[li].hilike,C[0]->G[li].hiprob);
        }
        FP "\n");
      }
#endif */

    }
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    //AS: Wed Apr 20 16:46:50 EDT 2016 totaltopolupdates, chain0topolupdates,
    //chaintopolswaps have to be reduced before printing
    int z;
    int poptreenum_rec;
 #ifdef MPI_ENABLED
    int rc;
    MPI_Status status;
    if (numprocesses > 1) // if multiple processes reduce sum 
    {
	     rc = MPI_Reduce(&totaltopolupdates, &totaltopol_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
			    MPI_Abort(MPI_COMM_WORLD, rc);
			   rc = MPI_Reduce(&chain0topolupdates, &chain0topol_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			   if (rc != MPI_SUCCESS)
				    MPI_Abort(MPI_COMM_WORLD, rc);
			   rc = MPI_Reduce(&chain0topolswaps, &chain0topolswaps_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			   if (rc != MPI_SUCCESS)
				    MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif	
    z = whichiscoldchain();
    if (z >=0 && outfile != NULL) 
      poptreenum_rec = C[z]->poptreenum;
#ifdef MPI_ENABLED
    else
    {
      if (z >=0 && outfile == NULL) 
      {
         rc = MPI_Send(&C[z]->poptreenum, 1, MPI_INT, 0, 7313, MPI_COMM_WORLD);
	    if (rc != MPI_SUCCESS)
		    MPI_Abort(MPI_COMM_WORLD, rc);
      }
      if (z < 0 && outfile != NULL) 
      {
         rc = MPI_Recv(&poptreenum_rec, 1,MPI_INT, MPI_ANY_SOURCE, 7313, MPI_COMM_WORLD, &status);
	    if (rc != MPI_SUCCESS)
		    MPI_Abort(MPI_COMM_WORLD, rc);
      }
    }
#endif
    if(outfile != NULL) 
    {
      if (numchainstotal > 1) 
      {
        //FP "Last Sampled Population Topology #: %d   String: %s\n\n",poptreenum_rec,alltreestrings[poptreenum_rec]);   this just not very useful in output file
        FP "\nPopulation Topology Sampling\n----------------------------\n");
        if (mcmcrecords > 0)
          FP "  Number of sampled topologies: %d\n",mcmcrecords);
        if (mcmcrecords_old > 0)
        {
          FP "  Number of topologies sampled in previous runs: %d\n",mcmcrecords_old);
          FP "  Total number of sampled topologies: %d\n",mcmcrecords+mcmcrecords_old);
        }
        if (numprocesses > 1)
        {
  		    FP "  Total Number of Accepted Topology Updates Across all Chains: %d\n",totaltopol_rec);
	        FP "  Number of Accepted Topology Updates to the Cold Chain: %d\n",chain0topol_rec);
          FP "  Number of Accepted Swaps that Changed Topology for the Cold Chain: %d\n",chain0topolswaps_rec);
        }
        else
        {
          FP "  Total Number of Topology Updates Across all Chains: %d\n",totaltopolupdates);
          FP "  Number of Accepted Topology Updates to the Cold Chain: %d\n",chain0topolupdates);
	        FP "  Number of Accepted Swaps that Changed Topology for the Cold Chain: %d\n",chain0topolswaps);
        }
        FP  "\n");
      }
      else  // 1 chain
      {
        FP "  Total Number of Topology Updates Across: %d\n",totaltopolupdates);
      }
    }
  }
  /* interval output of detailed chain info - used for debugging 
  int j;
  int ci,li,i,temp;
  FP  "\nchain#\tallbeta\tbetachain#\t#Swaps\tswaprate\tpdg\tprobhgg\tprobg\tsumP\t#mig\ttree#\n");

  for (ci=0;ci<numchainstotal-1;ci++)
  {
    for (i=0;i<numchainstotal;i++)
      if (allbetas[ci] == beta[i])
        break;
    for (j=0;j<numchainstotal;j++)
      if (allbetas[ci+1] == beta[j])
        break;
    temp = 0;
    for (li=0;li<nloci;li++)
      temp += C[i]->G[li].mignum;
    FP  "%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n",ci,allbetas[ci],i,tempbasedswapcount[ci][ci + 1],tempbasedswapcount[ci][ci + 1] / (float) tempbasedswapcount[ci + 1][ci], C[i]->allpcalc.pdg,C[i]->allpcalc.probhgg,C[i]->allpcalc.probg,(C[i]->allpcalc.pdg+C[i]->allpcalc.probhgg+C[i]->allpcalc.probg),temp,C[i]->poptreenum);
  }*/
}                               /* printrunbasics */


#define MAXLENX  200
#define MAXLENY  104
#define DEFAULTLENX 150;
#define DEFAULTLENY  40;

/* 4/27/09  JH modified so that the full x axis is used even if there are < TRENDDIM points in the array */ 

void asciitrend (FILE * outfile, struct value_record *v, int trenddoublepoint,int trendspot)
/* structure of graph 
0 name
1 ytop    |
40 ybot   |
41         -------
42        xbot    xtop

trenddoublepoint is the position in the array below which the points have twice the density 
with respect to the points abouve the trenddoublepoint 

To plot the array so that the plotted positions are proportional to time,  we adjust the plotting positions as a function 
of the density difference and the position of trenddoublepoint. 
If a is where trenddoublepoint would be plotted without adjustment and x is the value by which to multiply positions from 0 up to a,  to get 
their adjusted plotposition  then 

Solve[x a + (x/2) (d - a) == d, x]

This is because we know that the total distance must be d (i.e. available length of x axis)
then {x -> (2 d)/(a + d)}}

let f(a,d) = (2 d)/(a + d)
Then for a point at position c in the array and that is lower than position a,  the plot position is c*f(a, d)
position a itself is plotted at a * f(a,d)
For a point at position c in the array that is higher than position a,  the plot position is 
 a * f(a,d) + (c - a)*f(a, d)/2
Thus we rescale points up to point trendoublepoint by a factor of x 
and we rescale points above that by a factor of x/2 
*/
{
  char graph[MAXLENY][MAXLENX];
  char tc[20];
  double yt, ymax = -1e200, ymin = 1e200;
  double rescalerbelow, rescalerabove;
  double *y;
  double logscale;
  int i, xspot, xtemp, xmax, yspot, notplot = 0, rescaledoublepoint;
  int xlen, ylen, dim, localdoublepoint;
  y = v->trend;
  if (v->do_logplot)
    logscale = UMAX;
  xlen = DEFAULTLENX;
  ylen = DEFAULTLENY;
  if (trendspot < TRENDDIM)
  {
    dim = trendspot;
    localdoublepoint = trendspot;
  }
  else
  {
    dim = TRENDDIM;
    localdoublepoint = trenddoublepoint;
  }

  for (i = 0; i < dim; i++)
  {
    yt = y[i];
    if (ymax < yt)
      ymax = yt;
    if (ymin > yt)
      ymin = yt;
  }
  if (v->do_logplot)
  {
    if (ymax > (double) logscale)
      ymax = (double) logscale;
    if (ymin < 1 / (double) logscale)
      ymin = 1 / (double) logscale;
  }
  strcpy (graph[0], v->str);
  strcat (graph[0], " trend");
  if (v->do_logplot)
    strcat (graph[0], " - log scale");
  sprintf (tc, "%8.4g", ymax);
  strcpy (graph[1], tc);
  while ((int) strlen (graph[1]) < xlen - 1)
    strcat (graph[1], " ");
  sprintf (tc, "%8.4g", ymin);
  strcpy (graph[ylen - 2], tc);
  while ((int) strlen (graph[ylen - 2]) < xlen - 1)
    strcat (graph[ylen - 2], " ");
  strcpy (graph[ylen - 1], "          ");
  while ((int) strlen (graph[ylen - 1]) < xlen - 1)
    strcat (graph[ylen - 1], "-");
  for (i = 2; i < ylen - 2; i++)
  {
    graph[i][0] = '\0';
    while ((int) strlen (graph[i]) < xlen - 1)
      strcat (graph[i], " ");
  } 
  
  for (i = 1; i < ylen - 1; i++)
    graph[i][10] = '|';

  xmax = (int) (xlen - 12);
  rescaledoublepoint =  INTEGERROUND (((float) (xmax * localdoublepoint) / dim));
  rescalerbelow = (2 * (float) xmax) / (float) (rescaledoublepoint + xmax);
  rescalerabove = rescalerbelow / 2;

  for (i = 0; i < dim; i++)
  {
    if (y[i] >= ymin && y[i] <= ymax)
    {
      if (v->do_logplot)
        yspot = (int) 1 + (ylen - 2) - (int) ((ylen - 2) * (log (y[i]) - log (ymin)) / (log (ymax) - log (ymin)));

      else
        yspot = (int) 1 + (ylen - 2) - (int) ((ylen - 2) * (y[i] - ymin) / (ymax - ymin));
    }
    else
    {
      yspot = -1;
      notplot = 1;
    }
    xtemp = INTEGERROUND (((float) (xmax * i) / dim));
    if (i <= localdoublepoint)
      xspot = (int) (xtemp * rescalerbelow);
    else
      xspot = (int) (rescaledoublepoint * rescalerbelow + (xtemp - rescaledoublepoint) * rescalerabove);
    xspot += (int) 11;
    if (xspot < xlen - 1 && yspot >= 1)
      graph[yspot][xspot] = '*';
  }
  for (i = 0; i < ylen; i++)
    FP "%s \n", graph[i]);
  if (notplot)
    FP "    - not all values plotted, some exceed bounds: %5.4f - %5.0f \n", ymin, ymax);
  FP "\n");
}                               /* asciitrend */

// Constants for asciicurve ACXPLOT and AXYPLOT  determine the area of the plot
#define  ACXLEFTSPACE 10
#define  ACXPLOT  75
#define ACXMAX (ACXLEFTSPACE + ACXPLOT)
#define ACYPLOT 50
#define  ACYBOTSPACE 3
#define ACYMAX (ACYPLOT + ACYBOTSPACE)
#define  ASCIICURVEMINVAL 1e-6

/*for logscale,  any integer not 0  makes the plot use a log scale */
/* does not set the correct scale for t */
/* position of rows of graph 
0 name
1               |
ACYPLOT         |
ACYPLOT + 1     -------
ACYPLOT + 2    
*/

void asciicurve (FILE * outfile, struct plotpoint *a, char *qlabel,
                 int logscale, int mcmcrecords)
{
  char graph[ACYMAX][ACXMAX + 1];
  char tc[13];
  double ymax = -1e10, ymin = 1e10;
  int xmax;
  int i, xspot, yspot;

  for (i = 0; i < ACYMAX; i++)
    graph[i][0] = '\0';
  for (i = 0; i < GRIDSIZE; i++)
  {
    if (ymax < a[i].y)
      ymax = a[i].y;
    if (ymin > a[i].y)
      ymin = a[i].y;
  }
  ymax /= mcmcrecords;
  ymin = 0;                     // don't shift plot on y axis 
  if (logscale)
  {
    xmax = GRIDSIZE - 1;
  }
  else
  {
    xmax = -1;
    i = GRIDSIZE - 1;
    while (xmax == -1)
    {
      if (fabs (a[i].y) > ASCIICURVEMINVAL)
        xmax = i;
      i--;
    }
    if (xmax < 0)
      xmax = GRIDSIZE - 1;
  }
// set up plot name line 
  strcat (graph[0], qlabel);
  strcat (graph[0], " curve");
  //sprintf(graph[0],"%s curve %d pts",qlabel,mcmcrecords);
  while (strlen (graph[0]) < ACXMAX)
    strcat (graph[0], " ");
//set up the upper y axis label
  sprintf (tc, "%8.4g", ymax);
  strcpy (graph[1], tc);
  while (strlen (graph[1]) < ACXMAX)
    strcat (graph[1], " ");
//set up the lower y axis label
  sprintf (tc, "%8.4f", ymin);
  strcpy (graph[ACYPLOT], tc);
  while (strlen (graph[ACYPLOT]) < ACXMAX)
    strcat (graph[ACYPLOT], " ");
//set up the horizontal line on the x axis
  for (i = 0; i < ACXLEFTSPACE; i++)
    strcat (graph[ACYPLOT + 1], " ");
  while (strlen (graph[ACYPLOT + 1]) < ACXMAX)
    strcat (graph[ACYPLOT + 1], "-");
// set up the x axis label line
  for (i = 0; i < ACXLEFTSPACE; i++)
    strcat (graph[ACYPLOT + 2], " ");
  sprintf (tc, "%8.4f", a[0].x);
  strcat (graph[ACYPLOT + 2], tc);
  sprintf (tc, "%8.4f", a[xmax].x);
  if (logscale)
    strcat (graph[ACYPLOT + 2], "           Log Scale");
  while (strlen (graph[ACYPLOT + 2]) < ACXMAX - strlen (tc) - 1)
    strcat (graph[ACYPLOT + 2], " ");
  strcat (graph[ACYPLOT + 2], tc);
  while (strlen (graph[ACYPLOT + 2]) < ACXMAX)
    strcat (graph[ACYPLOT + 2], " ");
// fill the plot with empty space
  for (i = 2; i < ACYPLOT; i++)
    while (strlen (graph[i]) < ACXMAX)
      strcat (graph[i], " ");
  // set up the vertical line on the Y axis
  for (i = 1; i < ACYPLOT + 1; i++)
    graph[i][ACXLEFTSPACE] = '|';
  i = 0;
  while (i <= xmax)
  {
    yspot = (int) 1 + ACYPLOT - (int) (ACYPLOT * (a[i].y / mcmcrecords - ymin) / (ymax - ymin));
    if (logscale)
      xspot = (int) ACXLEFTSPACE + 1 + (int) ((ACXPLOT - 2) * (log (a[i].x) - log (a[0].x)) / (2 * log (a[xmax].x)));
    else
      xspot = (int) ACXLEFTSPACE + 1 + (int) ((ACXPLOT - 2) * (a[i].x - a[0].x) / (a[xmax].x - a[0].x));
//    assert (xspot < ACXMAX);  some bug here cauess xspot to be greater than ACXMAX ?? 
    assert (yspot < ACYMAX);
    if (xspot < ACXMAX && xspot >= (ACXLEFTSPACE + 1) && yspot < ACYMAX && yspot >= 1)
      graph[yspot][xspot] = '*';
    i++;
  }
  for (i = 0; i < ACYMAX; i++)
    FP "%s \n", graph[i]);
  fprintf (outfile, "\n");
}                               /* asciicurve */

// print acceptance rates for each the numrec elements pointed to in the array of pointers rec 
// assume that all elements have the same value for num_uptypes
#undef MAXLENX
#undef MAXLENY
void printacceptancerates (FILE * outto, int numrec,
                           struct chainstate_record_updates_and_values *rec[],
                           const char *printstring)
{

  int i, j;
  char numstr[20];

  if (outto != NULL)  fprintf(outto, "\n%s\n", printstring);
  for (i = 0; i < (int) strlen (printstring); i++)
    if (outto != NULL)  fprintf(outto, "-");
  if (outto != NULL)  fprintf(outto, "\n");

  if (outto != NULL)  fprintf(outto, " Update Type:");
  for (i = 0; i < rec[0]->num_uptypes; i++)
    if (outto != NULL)  fprintf(outto, "\t%s\t", rec[0]->upnames[i]);
  if (outto != NULL)  fprintf(outto, "\n");
  if (outto != NULL)  fprintf(outto, "            ");
  for (i = 0; i < rec[0]->num_uptypes; i++)
    if (outto != NULL)  fprintf(outto, "\t#Tries\t#Accp\t%%");
  if (outto != NULL)  fprintf(outto, "\n");
  for (j = 0; j < numrec; j++)
  {
    if (outto != NULL)  fprintf(outto, " %s", rec[j]->str);
    i = (int) strlen (rec[j]->str);
    while (i < 13)
    {
      if (outto != NULL)  fprintf(outto, " ");
      i++;
    }
    for (i = 0; i < rec[j]->num_uptypes; i++)
    {
      sprintf (&numstr[0], "%.3e", (float) rec[j]->upinf[i].tries);
      if (outto != NULL)  fprintf(outto, "\t%s", shorten_e_num (&numstr[0]));
      sprintf (&numstr[0], "%.3e", (float) rec[j]->upinf[i].accp);
      if (outto != NULL)  fprintf(outto, "\t%s", shorten_e_num (&numstr[0]));
      if (rec[j]->upinf[i].tries >= 0.5)  // tries is a double,  but accumulates as if an int
        {
          if (outto != NULL)  fprintf(outto, "\t%.3f", (float) 100 * rec[j]->upinf[i].accp / rec[j]->upinf[i].tries);
      }
      else
      {
        if (outto != NULL)  fprintf(outto, "\tna");
      }
    }
    if (outto != NULL)  fprintf(outto, "\n");
  }
}                               /* printacceptancerates() */

                           
/* To print acceptance rates:
	reclist[] is an array of pointers to struct chainstate_record_updates_and_values
	set values of reclist[] to those structures for which you want to print acceptance rates
	call printacceptancerates ()
	*/

void callprintacceptancerates (FILE * outto, int currentid)
{
  int i, j, li;
  // length of this array must be fairly long, although it is technically possible to have MAXLOCI * MAXLINKED records,  but very unlikely
  struct chainstate_record_updates_and_values *reclist[MAXLOCI + MAXLINKED];
  struct chainstate_record_updates_and_values *reclist_saved[MAXLOCI + MAXLINKED];
  int rc = 0; //AS: return code for MPI C bindings
  double *y1 = static_cast<double *> (malloc (numsplittimes * sizeof (double)));
  double *y_rec = static_cast<double *> (malloc (numsplittimes * sizeof (double)));
  double *z1 = static_cast<double *> (malloc (numsplittimes * sizeof (double)));
  double *z_rec = static_cast<double *> (malloc (numsplittimes * sizeof (double)));


/* SPLITTING TIME UPDATE RATES 
   ---------------------------*/

/* RY update MPI_REDUCE */
#ifdef MPI_ENABLED
  if (numprocesses > 1 && doRYupdate) 
  {
    for (int x = 0; x < numsplittimes; x++) 
    {
	    y1[x] = T[x].upinf[update_type_RY].tries;
	    z1[x] = T[x].upinf[update_type_RY].accp;
    }
		rc = MPI_Reduce(y1, y_rec, numsplittimes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
		rc = MPI_Reduce(z1, z_rec, numsplittimes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS)
	    MPI_Abort(MPI_COMM_WORLD, rc);
	  if (currentid == HEADNODE) 
    {
	    for (int x = 0; x < numsplittimes; x++) 
      {	
		    T[x].upinf[update_type_RY].tries = y_rec[x];
		    T[x].upinf[update_type_RY].accp = z_rec[x];
	    }
    }
  }
#endif //MPI_ENABLED

  double *y2 = static_cast<double *> (malloc (numsplittimes * sizeof (double)));
  double *z2 = static_cast<double *> (malloc (numsplittimes * sizeof (double)));

/* NW update MPI_REDUCE */

#ifdef MPI_ENABLED
  if (doNWupdate && numprocesses > 1) 
  {
    for (int x = 0; x < numsplittimes; x++) 
    {
      y2[x] = T[x].upinf[update_type_NW].tries;
      z2[x] = T[x].upinf[update_type_NW].accp;
    }
    rc = MPI_Reduce(y2, y_rec, numsplittimes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
	
    rc = MPI_Reduce(z2, z_rec, numsplittimes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
	  if (currentid == HEADNODE) 
    {
      for (int x = 0; x < numsplittimes; x++) 
      {
        T[x].upinf[update_type_NW].tries = y_rec[x];
        T[x].upinf[update_type_NW].accp = z_rec[x];
      }
    }
  }
#endif //MPI_ENABLED


/* print RY and NW update rates  */
  for (i = 0; i < numsplittimes; i++)
    reclist[i] = (T + i);

  if (currentid == HEADNODE)
    if (numsplittimes > 0) printacceptancerates (outto, numsplittimes, reclist,
                          "Update Rates (chain 0) -- Population Splitting Times");
  XFREE(y_rec);
  XFREE(z_rec);


#ifdef MPI_ENABLED
  if (currentid == HEADNODE && numprocesses > 1 && doRYupdate) 
  {
    for (int j = 0; j < numsplittimes; j++) 
    {
      T[j].upinf[update_type_RY].tries = y1[j];
      T[j].upinf[update_type_RY].accp = z1[j];
    }
  }
#endif
	 XFREE(y1);
	 XFREE(z1);

#ifdef MPI_ENABLED
  if (doNWupdate && currentid == HEADNODE && numprocesses > 1)
  {
    for (int j  = 0; j < numsplittimes; j++) 
    {
			
      T[j].upinf[update_type_NW].tries = y2[j];
      T[j].upinf[update_type_NW].accp = z2[j];
    }
  }
#endif
	 XFREE(y2);
	 XFREE(z2);

/* DONE RY AND NW SPLITTING TIME UPDATE ACCEPTANCE RATE STUFF*/

/* HYPERPRIOR UPDATE RATES 
   -----------------------*/

/* update migration hyperprior stuff, mh */
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    double *y1 = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
    double *y_rec = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
    double *z1 = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
    double *z_rec = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));
#ifdef MPI_ENABLED
    if (numprocesses > 1) 
    {
      for (int x = 0; x < nummigrateparams; x++) 
      {
	      y1[x] = mh[x].upinf[PRIOR_UPDATE].tries;
	      z1[x] = mh[x].upinf[PRIOR_UPDATE].accp;
      }
		    rc = MPI_Reduce(y1, y_rec, nummigrateparams, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
		    rc = MPI_Reduce(z1, z_rec, nummigrateparams, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
	        MPI_Abort(MPI_COMM_WORLD, rc);
	      if (currentid == HEADNODE) 
      {
	      for (int x = 0; x < nummigrateparams; x++) 
        {	
		      mh[x].upinf[PRIOR_UPDATE].tries = y_rec[x];
		      mh[x].upinf[PRIOR_UPDATE].accp = z_rec[x];
	      }
      }
    }
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      if (numprocesses > 1) 
      {
        y1[0] = mhnit->upinf[PRIOR_UPDATE].tries;
	       z1[0] = mhnit->upinf[PRIOR_UPDATE].accp;
		      rc = MPI_Reduce(y1, y_rec,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		      rc = MPI_Reduce(z1, z_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	       if (currentid == HEADNODE) 
        {
		        mhnit->upinf[PRIOR_UPDATE].tries = y_rec[0];
		        mhnit->upinf[PRIOR_UPDATE].accp = z_rec[0];
	        }
      }
    }
#endif //MPI_ENABLED
/* print migration hyperprior updates */
    for (i = 0; i < nummigrateparams; i++)
      reclist[i] = (mh + i);

    if (currentid == HEADNODE)
      if (nummigrateparams > 0) printacceptancerates (outto, nummigrateparams, reclist,
                            "Update Rates (chain 0) -- Migration Priors in Poptree");
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      reclist[0] = mhnit;
      if (currentid == HEADNODE)
        printacceptancerates (outto, 1, reclist, "Update Rates (chain 0) -- Migration Priors Not in Poptree");
    }
    XFREE(y_rec);
    XFREE(z_rec);
#ifdef MPI_ENABLED
    if (currentid == HEADNODE && numprocesses > 1) 
    {
      for (int j = 0; j < nummigrateparams; j++) 
      {
        mh[j].upinf[PRIOR_UPDATE].tries = y1[j];
        mh[j].upinf[PRIOR_UPDATE].accp = z1[j];
      }
    }
#endif
	  XFREE(y1);
	  XFREE(z1);
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])   // now do popsize hyperprior stuff 
  {
    double *y1 = static_cast<double *> (malloc (numpopsizeparams * sizeof (double)));
    double *y_rec = static_cast<double *> (malloc (numpopsizeparams * sizeof (double)));
    double *z1 = static_cast<double *> (malloc (numpopsizeparams * sizeof (double)));
    double *z_rec = static_cast<double *> (malloc (numpopsizeparams * sizeof (double)));
#ifdef MPI_ENABLED
    if (numprocesses > 1) 
    {
      for (int x = 0; x < numpopsizeparams; x++) 
      {
	      y1[x] = qh[x].upinf[PRIOR_UPDATE].tries;
	      z1[x] = qh[x].upinf[PRIOR_UPDATE].accp;
      }
		    rc = MPI_Reduce(y1, y_rec, numpopsizeparams, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
		    rc = MPI_Reduce(z1, z_rec, numpopsizeparams, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
	        MPI_Abort(MPI_COMM_WORLD, rc);
	      if (currentid == HEADNODE) 
      {
	      for (int x = 0; x < numpopsizeparams; x++) 
        {	
		      qh[x].upinf[PRIOR_UPDATE].tries = y_rec[x];
		      qh[x].upinf[PRIOR_UPDATE].accp = z_rec[x];
	      }
      }
    }
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      if (numprocesses > 1) 
      {
        y1[0] = qhnit->upinf[PRIOR_UPDATE].tries;
	       z1[0] = qhnit->upinf[PRIOR_UPDATE].accp;
		      rc = MPI_Reduce(y1, y_rec,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		      rc = MPI_Reduce(z1, z_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	       if (currentid == HEADNODE) 
        {
		        qhnit->upinf[PRIOR_UPDATE].tries = y_rec[0];
		        qhnit->upinf[PRIOR_UPDATE].accp = z_rec[0];
	        }
      }
    }
#endif //MPI_ENABLED
/* print migration hyperprior updates */
    for (i = 0; i < numpopsizeparams; i++)
      reclist[i] = (qh + i);
    if (currentid == HEADNODE)
      if (numpopsizeparams > 0) printacceptancerates (outto, numpopsizeparams, reclist,
                            "Update Rates (chain 0) -- Population Size Priors in Poptree");
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      reclist[0] = qhnit;
      if (currentid == HEADNODE)
        printacceptancerates (outto, 1, reclist, "Update Rates (chain 0) -- Population Priors Not in Poptree");
    }
    XFREE(y_rec);
    XFREE(z_rec);
#ifdef MPI_ENABLED
    if (currentid == HEADNODE && numprocesses > 1) 
    {
      for (int j = 0; j < numpopsizeparams; j++) 
      {
        qh[j].upinf[PRIOR_UPDATE].tries = y1[j];
        qh[j].upinf[PRIOR_UPDATE].accp = z1[j];
      }
    }
#endif
	  XFREE(y1);
	  XFREE(z1);
  } // popsize hyper prior updates

/* DONE HYPERPRIOR UPDATE ACCEPTANCE RATE STUFF*/

/* POPULATION TOPOLOGY UPDATE RATES 
   --------------------------------*/
  double a1,a2,a3, a_rec, b1,b2,b3, b_rec;
  //if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)   // fixed 5/3/2017
  if (hiddenoptions[HIDDENGENEALOGY] == 1)
  {
    i = 0;
    reclist_saved[i] = poptreeuinfo;
#ifdef MPI_ENABLED
	//AS: 4/13/2016 Changing this part of the function
//    int x, y, y_rec, z, z_rec;
    //if (numprocesses > 1)  // 3/29/2017  commented this out,  was causing no results printed when only 1 cpu 
    {
      a1 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries;
      b1 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp;
    }
    rc = MPI_Reduce(&a1, &a_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
		  MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(&b1, &b_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
  		MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries = a_rec;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp = b_rec;
	  }
    a2 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries;
    b2 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp;
    rc = MPI_Reduce(&a2, &a_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  	if (rc != MPI_SUCCESS)
	  	MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(&b2, &b_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  	if (rc != MPI_SUCCESS)
	  	MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries = a_rec;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp = b_rec;
	  }
    a3 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries;
    b3 = poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp;
    rc = MPI_Reduce(&a3, &a_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD,rc);
    rc = MPI_Reduce(&a3, &a_rec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  if (rc != MPI_SUCCESS)
		  MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries = a_rec;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp = b_rec;
	  }
#endif //MPI_ENABLED
    i = 0;
    reclist[0] = poptreeuinfo;
    if (currentid == HEADNODE)
      printacceptancerates (outto, 1, reclist,
                          "Update Rates (chain 0) -- Branch Slide Updates to Population Trees");
#ifdef MPI_ENABLED
    if (numprocesses > 1 && currentid == HEADNODE) 
    {
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries = a1;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries = a2;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries = a3;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp = b1;

      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp = b2;
      poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp = b3;
    }
#endif //MPI_ENABLED
  } // hiddengenealogy==1
   
  double *c1 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *c2 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *c3 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *c_rec = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *d1 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *d2 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *d3 = static_cast<double *> (malloc (nloci * sizeof (double)));
  double *d_rec = static_cast<double *> (malloc (nloci * sizeof (double)));

#ifdef MPI_ENABLED

  if (numprocesses > 1) 
  {
    for (int x = 0; x < nloci; x++) 
    {
      c1[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries;
      d1[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp;
    }
    rc = MPI_Reduce(c1, c_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(d1, d_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      for (int x = 0; x < nloci; x++) 
      {
        L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries = c_rec[x];
        L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp = d_rec[x];
      }
    }
    for (int x = 0; x < nloci; x++) 
    {
      c2[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries;
      d2[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp;
    }
    rc = MPI_Reduce(c2, c_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(d2, d_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      for (int x = 0; x < nloci; x++) 
      {
        L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries = c_rec[x];
        L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp = d_rec[x];
      }
    }
    for (int x = 0; x < nloci; x++) 
    {
      c3[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries;
      d3[x] = L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp;
    }
    rc = MPI_Reduce(c3, c_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(d3, d_rec, nloci, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    if (currentid == HEADNODE) 
    {
      for (int x = 0; x < nloci; x++) 
      {
      L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries = c_rec[x];
      L[x].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp = d_rec[x];
      }
    }
  }

#endif  //MPI_ENABLED
  if (hiddenoptions[SKIPMOSTUSCALAROUTPUT] == 0)
  {
    for (li = 0; li < nloci; li++)
      reclist[li] = L[li].g_rec;
    if (currentid == HEADNODE)
	     printacceptancerates (outto, nloci, reclist, "Update Rates (chain 0) -- Genealogies");
  }

#ifdef MPI_ENABLED
  if (numprocesses > 1 && currentid == HEADNODE) 
  {
    for (li = 0; li < nloci; li++) 
    {
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries = c1[li];
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries = c2[li];
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries = c3[li];
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp = d1[li];
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp = d2[li];
      L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp = d3[li];
    }
  }
#endif  //MPI_ENABLED
  XFREE(c1);
  XFREE(c2);
  XFREE(c3);
  XFREE(d1);
  XFREE(d2);
  XFREE(d3);
  XFREE(c_rec);
  XFREE(d_rec);

/* DONE TOPOLOGY SLIDING  UPDATE ACCEPTANCE RATE STUFF*/

/* MUTATION SCALAR UPDATE RATES 
   ----------------------------*/
  if (domutationscalarupdate)
    {
    double **e = static_cast<double **> (malloc (nloci * sizeof(double *)));
    double **e_rec = static_cast<double **> (malloc (nloci * sizeof(double *)));
    double **f = static_cast<double **> (malloc (nloci * sizeof(double *)));
    double **f_rec = static_cast<double **> (malloc (nloci * sizeof(double *)));

    for (li = 0;  li < nloci; li++) 
    {
      e[li] = static_cast<double *> (malloc (L[li].nlinked * sizeof (double)));
      e_rec[li] = static_cast<double *> (malloc (L[li].nlinked * sizeof (double)));
      f[li] = static_cast<double *> (malloc (L[li].nlinked * sizeof (double)));
      f_rec[li] = static_cast<double *> (malloc (L[li].nlinked * sizeof (double)));
    }
  // mutation rate scalars
    if (nurates > 1  && (runoptions[PRINTMUTATIONUPDATESTOSCREEN] || outto != stdout))
    {
#ifdef MPI_ENABLED
      if (numprocesses > 1) 
      {
        for (li = 0; li < nloci; li++) 
        {
          for (j = 0; j < L[li].nlinked; j++) 
          {
            e[li][j] = L[li].u_rec[j].upinf->tries;
            f[li][j] = L[li].u_rec[j].upinf->accp;
          }
          rc = MPI_Reduce(e[li], e_rec[li], L[li].nlinked, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
          rc = MPI_Reduce(f[li], f_rec[li], L[li].nlinked, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
          if (currentid == HEADNODE) 
          {
            for (j = 0; j < L[li].nlinked; j++) 
            {
	            L[li].u_rec[j].upinf->tries = e_rec[li][j];
              L[li].u_rec[j].upinf->accp = f_rec[li][j];
            }
          }
        }
      }
#endif  //MPI_ENABLED
      if (hiddenoptions[SKIPMOSTUSCALAROUTPUT]==0)
      {
      for (i = 0, li = 0; li < nloci; li++)
        for (j = 0; j < L[li].nlinked; j++)
        {
          reclist[i] = &L[li].u_rec[j];
          i++;
        }
  // kappa values for HKY model
        if (currentid == HEADNODE  && domutationscalarupdate) 
          printacceptancerates (outto, i, reclist,"Update Rates (chain 0) -- Mutation Rate Scalars");
      }
#ifdef MPI_ENABLED
      if (numprocesses > 1 && currentid == HEADNODE) 
      {
        for (li = 0; li < nloci; li++) 
        {
          for (j = 0; j < L[li].nlinked; j++) 
          {
            L[li].u_rec[j].upinf->tries = e[li][j];
            L[li].u_rec[j].upinf->accp = f[li][j];
          }
        }
      }
#endif  //MPI_ENABLED

#ifdef MPI_ENABLED
      if (numprocesses > 1) 
      {
        for (li = 0;  li < nloci; li++) 
        {
          for (j = 0; j < L[li].nlinked; j++) 
          {
            if (L[li].umodel[j] == HKY) 
            {
              e[li][j] = L[li].kappa_rec[j].upinf->tries;
              f[li][j] = L[li].kappa_rec[j].upinf->accp;
            }
          }
          rc = MPI_Reduce(e[li], e_rec[li], L[li].nlinked, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
          rc = MPI_Reduce(f[li], f_rec[li], L[li].nlinked, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          if (rc != MPI_SUCCESS)
          MPI_Abort(MPI_COMM_WORLD, rc);
          if (currentid == HEADNODE) 
          {
            for (j = 0; j < L[li].nlinked; j++) 
            {
              if (L[li].umodel[j] == HKY) 
              {
                L[li].kappa_rec[j].upinf->tries = e_rec[li][j];
                L[li].kappa_rec[j].upinf->accp = f_rec[li][j];
              }
            }
          }
        }
      }
#endif  //MPI_ENABLED
      if (hiddenoptions[SKIPMOSTUSCALAROUTPUT]==0)
      {
        for (i = 0, li = 0; li < nloci; li++)
        for (j = 0; j < L[li].nlinked; j++)
        {
          if (L[li].umodel[j] == HKY)
          {
            reclist[i] = L[li].kappa_rec;
            i++;
          }
        }
        if (i > 0 && currentid == HEADNODE)
          printacceptancerates (outto, i, reclist,
                              "Update Rates (chain 0) -- HKY Model Kappa parameter");
      }
#ifdef MPI_ENABLED
      if (numprocesses > 1 && currentid == HEADNODE) 
      {
        for (li = 0; li < nloci; li++) 
        {
          for (j = 0; j < L[li].nlinked; j++) 
          {
            if (L[li].umodel[j] == HKY) 
            {
              L[li].kappa_rec[j].upinf->tries = e[li][j];
              L[li].kappa_rec[j].upinf->accp = f[li][j];
            }
          }
        }
      }
#endif //MPI_ENABLED
    }
    for (li = 0; li < nloci; li++) 
    {
      XFREE(e[li]);
      XFREE(f[li]);
      XFREE(e_rec[li]);
      XFREE(f_rec[li]);
    }
    XFREE(e);
    XFREE(f);
    XFREE(e_rec);
    XFREE(f_rec);
  }  
/* DONE mutation rate scalar and kappa update acceptatance rate stuff */


// STR ancestral allele states 
// A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
//8/26/2011  turn this printing section off, as it only ever prints zeros when A updating is not used 
    /*
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == STEPWISE)
        {
          reclist[i] = &L[li].A_rec[j];
          i++;
        }
      }
      
    if (i > 0)
    {
      printacceptancerates (outto, i, reclist,
                            "Update Rates (chain 0) -- STR Genealogy Allele States");
    } */

  return;
}

/* this code does not deal with MPI and multiple processors, so the values that are getting output  could be from any chain */ 
/* 4/19/2017  edited it so nothing is output if multiple cpus are in use  */

void printcurrentvals (FILE * outto)
{
  int li, i;
  if (numprocesses >  1)
    return;
  if (numsplittimes > 0)
  {
    if (outto != NULL)  fprintf(outto, "\nCurrent Splitting Time Values\n");
    if (outto != NULL)  fprintf(outto, "------------------------------\n");
    if (npops > 1)
    {
      if (outto != NULL)  fprintf(outto, "Split times");
      for (i = 0; i < lastperiodnumber; i++)
        if (outto != NULL)  fprintf(outto, "   t%i: %.3f", i, C[ARBCHAIN]->tvals[i]);
      if (outto != NULL)  fprintf(outto, "\n\n");
    }
  }
  if (runoptions[PRINTMUTATIONUPDATESTOSCREEN] && nloci > 1)
  {
    if (outto != NULL)  fprintf(outto, "Current Locus Likelihoods and Mutation Rate Scalars\n");
    if (outto != NULL)  fprintf(outto, "---------------------------------------------------\n");
    if (outto != NULL)  fprintf(outto, "Locus#    p(D|G)     u  ");

    if (outto != NULL)  fprintf(outto, "\n");
    for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].nlinked; i++)
      {
        if (outto != NULL)  fprintf(outto, "%2i", li);
        if (L[li].nlinked > 1)
        {
          if (outto != NULL)  fprintf(outto, "_%-2d", i);
        }
        else
        {
          if (outto != NULL)  fprintf(outto, "   ");
        }
        if (outto != NULL)  fprintf(outto, " %9.3lg %9.4lg", C[ARBCHAIN]->G[li].pdg_a[i], C[ARBCHAIN]->G[li].uvals[i]);
        if (outto != NULL)  fprintf(outto, "\n");
      }
    }
  }
  if (outto == stdout)
    fflush(stdout);
}                               /* printcurrentvals */

/* used for when there are multiple cpus,  only prints tvalues */ 
void printcurrent_tvals (FILE * outto,int  currentid)
{
  int i;
  double t_rec[MAXPOPS];
  int z = whichiscoldchain();
#ifdef MPI_ENABLED
  MPI_Status status;
#endif
  if (numsplittimes > 0)
  {
	      if (z >= 0 && currentid == HEADNODE) 
       {
         for (i=0;i<numsplittimes;i++)
           t_rec[i] = C[z]->tvals[i];
	      }
	#ifdef MPI_ENABLED
	      if (z < 0 && currentid == HEADNODE) 
        {
		        int rc = MPI_Recv(&t_rec, numsplittimes, MPI_DOUBLE, MPI_ANY_SOURCE, 1397, MPI_COMM_WORLD, &status);
		        if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
        }
	      if (z >=0 && currentid != 0) 
       {
         for (i=0;i<numsplittimes;i++)
           t_rec[i] = C[z]->tvals[i];
 		      int rc = MPI_Send(&t_rec, numsplittimes, MPI_DOUBLE, 0, 1397, MPI_COMM_WORLD);
	 	      if (rc != MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD, rc);
	      }
	#endif
  if (currentid == HEADNODE)
    {
      if (outto != NULL)  fprintf(outto, "\nCurrent Splitting Time Values\n");
      if (outto != NULL)  fprintf(outto, "------------------------------\n");
      if (outto != NULL)  fprintf(outto, "Split times");
      for (i = 0; i < lastperiodnumber; i++)
        if (outto != NULL)  fprintf(outto, "   t%i: %.3f", i, t_rec[i]);
      if (outto != NULL)  fprintf(outto, "\n\n");
      if (outto == stdout)
        fflush(stdout);
    }
  }
}                               /* printcurrentvals */


void savegenealogyfile (char genealogyinfosavefilename[], FILE * genealogyinfosavefile,
                   int *lastgenealogysavedvalue, int gsampinflength)
{
  int j, i;
  if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "a")) == NULL)
  {
    //printf ("Error opening treeinfosave file for writing\n");
    IM_err(IMERR_CREATEFILEFAIL,"Error opening treeinfosave file for writing");
  }
  for (j = *lastgenealogysavedvalue; j < genealogysamples; j++)
  {
    // save everything as a float  - round integers later, when they get used
    //std::cout << "What is going on here???";
//	for (i = 0; i < gsampinflength; i++) {
//		std::cout << gsampinf[0][i] << "\n";
//	}
    for (i = 0; i < gsampinflength; i++)
      fprintf (genealogyinfosavefile, "%.6f\t", (float) gsampinf[j][i]);
    if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
      fprintf(genealogyinfosavefile,"%s",debug_ti_addinfo[j]);
    fprintf (genealogyinfosavefile, "\n");
  } 
  *lastgenealogysavedvalue = genealogysamples - 1;
  FCLOSE (genealogyinfosavefile);
  return;
}                               /* savegenealogyfile */

void print_means_variances_correlations (FILE * outfile)
{
  double *means, *variances, **correlations;
  int i, p, q, np;

  np = numpopsizeparams + nummigrateparams;

  means = static_cast<double *> (calloc ((size_t) np, sizeof (double)));
  variances = static_cast<double *> (calloc ((size_t) np, sizeof (double)));
  if (npops > 1)
  {
    correlations = orig2d_alloc2Ddouble (np, np);
  }

  for (i = 0; i < genealogysamples; i++)
  {
    for (p = 0; p < np; p++)
    {
      means[p] += calcx (i, p, 0);
      variances[p] += calcx (i, p, 1);
    }
  }
  for (p = 0; p < np; p++)
  {
    if (means[p] >= 0.0)
      means[p] /= genealogysamples;
    if (variances[p] >= 0.0)
    {
      variances[p] /= genealogysamples;
      variances[p] -= SQR (means[p]);
    }
  }
  if (npops > 1)
  {
    for (p = 0; p < np; p++)
      for (q = 0; q < np; q++)
        correlations[p][q] = 0;
    for (i = 0; i < genealogysamples; i++)
    {
      for (p = 0; p < np - 1; p++)
        for (q = p + 1; q < np; q++)
          correlations[p][q] += calcx (i, p, 0) * calcx (i, q, 0);
    }
    for (p = 0; p < np - 1; p++)
      for (q = p + 1; q < np; q++)
      {
        if (correlations[p][q] >= 0.0)
        {
          correlations[p][q] /= genealogysamples;
          correlations[p][q] -= (means[p] * means[q]);
          correlations[p][q] /= sqrt (variances[p] * variances[q]);
        }
        else
          correlations[p][q] = -1.0;
      }
  }
  FP "%s",outputbanner("MEANS, VARIANCES and CORRELATIONS OF PARAMETERS ('$' r > 0.4  '*' r > 0.75)"));
  FP "Param:");
  for (p = 0; p < np; p++)
  {
    if (p < numpopsizeparams)
    {
      FP "\t%s", C[ARBCHAIN]->itheta[p].str);
    }
    else
    {
      if (C[ARBCHAIN]->imig[p - numpopsizeparams].pr.max > MPRIORMIN)
          FP "\t%s", C[ARBCHAIN]->imig[p - numpopsizeparams].str);
    }
  }
  FP "\nMean:");
  for (p = 0; p < np; p++)
    FP "\t%-.3lf", means[p]);
  FP "\nStdv:");
  for (p = 0; p < np; p++)
    FP "\t%-.3lf", sqrt (variances[p]));
  if (npops > 1)
  {
    FP "\n\nCorrelations\n");
    for (p = 0; p < np; p++)
    {
      if (p < numpopsizeparams)
      {
        FP "\t%s", C[ARBCHAIN]->itheta[p].str);
      }
      else  //p < numpopsizeparams + nummigrateparams
      {
        if (C[ARBCHAIN]->imig[p - numpopsizeparams].pr.max > MPRIORMIN)  
            FP "\t%s", C[ARBCHAIN]->imig[p - numpopsizeparams].str);
      }
    }
    FP "\n");
    for (p = 0; p < np; p++)
    {
      if (p < numpopsizeparams)
        FP "%s", C[ARBCHAIN]->itheta[p].str);
      else //p < numpopsizeparams + nummigrateparams
      {
        FP "%s", C[ARBCHAIN]->imig[p - numpopsizeparams].str);
      }
      for (q = 0; q < np; q++)
      {
        if (q == p)
          FP "\t  - ");
        else
        {
          if (q > p)
          {
            if (fabs (correlations[p][q]) < 0.4)
            {
              FP "\t%-.3lf", correlations[p][q]);
            }
            else
            {
              if (fabs (correlations[p][q]) < 0.75)
                FP "\t%-.3lf%s", correlations[p][q], "$");
              else
                FP "\t%-.3lf%s", correlations[p][q], "*");
            }
          }
          else
          {
            if (fabs (correlations[q][p]) < 0.4)
            {
              FP "\t%-.3lf", correlations[q][p]);
            }
            else
            {
              if (fabs (correlations[q][p]) < 0.75)
                FP "\t%-.3lf%s", correlations[q][p], "$");
              else
                FP "\t%-.3lf%s", correlations[q][p], "*");
            }
          }

        }
      }
      FP "\n");
    }
  }
  FP "\n");

  XFREE (means);
  XFREE (variances);
  if (npops > 1)
  {
    orig2d_free2D ((void **) correlations, np);
  }
  return;
}                               // print_means_variances_correlations


/* for printing population tree posteriors,  sorts and prints 
  A fairly complex function:
    -puts info about sampled trees in fatssarray[]
    - sorts this on basis of observed frequency
    - identifies all unique internal nodes (based on descendant sampled populations)
    - counts all unique internal nodes
    - uses these counts to calculate the product of the posterior clade probabilities (ppcp)
    - prints summaries and trees with non-zero probability (sorted) 
  Calls functions in alltreestrings.cpp 
*/
void sort_and_print_alltreestrings(FILE * outfile, int *poptopologycounts, int *poptopologyproposedlist_rec, char *topologypriorinfostring)
{
  foralltreestringsort *fatssarray;
  char **uniquenodes;
  int *nodecounts;
  int nunique = 0;
  int i,j;
  int totaltreecount = 0, numnonzero = 0;
  double temphi = -1.0;
  int hii;
  int counthipcp;
  int numnotproposed = 0;
  int uniformprior = (strlen(topologypriorinfostring) == 0);
  double temp1 = 1.0;
  int *set1counts,*set2counts, setl;

  for (i = 0; i < numpoptopologies; i++)
  {
    numnotproposed += (poptopologyproposedlist_rec[i] == 0);
  }

  fatssarray = static_cast<foralltreestringsort *> (malloc (numpoptopologies * sizeof (foralltreestringsort)));
  set1counts = static_cast<int *> (calloc (numpoptopologies, sizeof (int)));
  set2counts = static_cast<int *> (calloc (numpoptopologies, sizeof (int)));

  for (i = 0; i < numpoptopologies; i++)
  {
    totaltreecount += poptopologycounts[i];
    if (poptopologycounts[i] > 0)
      numnonzero += 1;
  }
  setl = poptopologysequence.currentlength/2;
  for (i=0,j=setl;i<setl;i++,j++)
  {
    set1counts[(int) (poptopologysequence.vals[i])]++;
    set2counts[(int) (poptopologysequence.vals[j])]++;
  }


  for (i = 0; i < numpoptopologies; i++)
  {
    strcpy(fatssarray[i].treestr,alltreestrings[i]);
    if (modeloptions[ADDGHOSTPOP]==1)
      strcpy(fatssarray[i].treestrnoghost,alltreestrings_noghost[i]);
    fatssarray[i].count = poptopologycounts[i];
    fatssarray[i].freqall =  poptopologycounts[i]/ (double) totaltreecount;
    fatssarray[i].origi = i;
    fatssarray[i].freqset1 = set1counts[i]/(double) setl;
    fatssarray[i].freqset2 = set2counts[i]/(double) setl;
  }
  
  qsort(fatssarray,numpoptopologies,sizeof(foralltreestringsort),foralltreestringsort_comp);
  uniquenodes = static_cast<char **> (malloc ((numuniquenodes[(npops - modeloptions[ADDGHOSTPOP])]) * sizeof (char *)));
  
  for (i = 0; i < numuniquenodes[(npops - modeloptions[ADDGHOSTPOP])]; i++)
  {
     uniquenodes[i] = static_cast<char *> (malloc (POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION * sizeof (char)));
  }
  nodecounts = static_cast<int *> (malloc ((numuniquenodes[(npops - modeloptions[ADDGHOSTPOP])]+1) * sizeof (int)));
  nodecounts[0] = 0;
  for (i=0;i<numpoptopologies;i++) if (fatssarray[i].count > 0)
  {
    if (modeloptions[ADDGHOSTPOP]==0)
      getnodecounts(uniquenodes,fatssarray[i].treestr,fatssarray[i].count,nodecounts, &nunique);
    else
      getnodecounts(uniquenodes,fatssarray[i].treestrnoghost,fatssarray[i].count,nodecounts, &nunique);
  }

  
//for (i=0;i<nunique;i++) printf("%s ",uniquenodes[i]);
#ifdef DEBUG
  int tempsum = 0;
  for (i=0;i<nunique;i++)
    tempsum += nodecounts[i];
  assert (tempsum == totaltreecount * (npops - (modeloptions[ADDGHOSTPOP]==1) - 2));
#endif
  
  counthipcp=0;
  for (i=0;i<numpoptopologies;i++) //if (fatssarray[i].count > 0)
  {
    if (modeloptions[ADDGHOSTPOP]==0)
      fatssarray[i].ppcp = calcppcp(uniquenodes,fatssarray[i].treestr,nodecounts,totaltreecount,nunique);
    else
      fatssarray[i].ppcp = calcppcp(uniquenodes,fatssarray[i].treestrnoghost,nodecounts,totaltreecount,nunique);
    if (fatssarray[i].ppcp > temphi)
    {
      temphi = fatssarray[i].ppcp;
      counthipcp = 1;
      hii = i;
    }
    if (fatssarray[i].ppcp == temphi)
      counthipcp += 1;
  }
  FP "%s",outputbanner("Estimated Posterior Probabilities of Population Tree Topologies"));
  FP "Population Names\n");
  FP "----------------\n");
  for (i = 0; i < npops; i++)
  {
    if (i==npops-1 && modeloptions[ADDGHOSTPOP]==1)
      FP " Population %d : ghost \n", i);
    else
      FP " Population %d : %s \n", i, popnames[i]);
  }
  FP "\n");
  FP "Priors\n");
  FP "------\n");
  if (uniformprior)
    FP " Uniform prior assumed for population tree topologies\n");
  else
    FP " Sister population priors specified on command line: %s\n",topologypriorinfostring);
  FP "\n");
  FP "Summary\n");
  FP "-------\n");
  if (uniformprior)
  {
    FP " Uniform prior assumed for population tree topologies\n");
  }
  FP " Total number of tree topologies recorded: %d\n",totaltreecount); 

  FP " Number of possible unique tree topologies for %d sampled populations: %d\n",(npops - modeloptions[ADDGHOSTPOP]),numpoptopologies);
  FP " Number of unique tree topologies recorded at least once: %d\n",numnonzero);
  FP " Number of unique tree topologies not observed: %d\n",numpoptopologies-numnonzero);
  FP " Number of unique tree topologies never proposed in MCMCMC: %d\n",numnotproposed);
  if (modeloptions[ADDGHOSTPOP]==0)
  {
    FP " Tree topology with highest estimated posterior probability: %s\n",fatssarray[0].treestr);
    FP "     -Newick string (unit branch lengths):\n       ");
    printnewickstring(outfile,fatssarray[0].treestr,NULL,0);
    FP "\n");
  }
  else
  {
    FP " Tree topology with highest estimated posterior probability: %s\n",fatssarray[0].treestrnoghost);
    FP "     -estimated posterior probability: %.6f\n",fatssarray[0].freqall);
    if (uniformprior == 0)
      FP "     -prior probability: %.6f\n",exp(topologypriors[fatssarray[0].origi]));
    FP "     -Newick string without ghost (unit branch lengths):\n       ");
    printnewickstring(outfile,fatssarray[0].treestrnoghost,NULL,0);
    FP "\n");
    FP "     -topology with ghost outgroup: %s\n",fatssarray[0].treestr);
    FP "     -Newick string with ghost outgroup (unit branch lengths):\n       ");
    printnewickstring(outfile,fatssarray[0].treestr,NULL,1);
    FP "\n");
  }
  FP " Product of the posterior clade probabilities (ppcp), example tree topology with highest value : %s\n",fatssarray[hii].treestr);
  if (modeloptions[ADDGHOSTPOP]==1)
    FP "     -example topology without ghost outgroup included: %s\n",fatssarray[hii].treestrnoghost);
  FP "     -estimated product of the posterior clade probabilities (ppcp): %.7f\n",temphi);
  FP "     -total number of topologies with this same ppcp value: %d\n",counthipcp);

  FP "\nTree Topology Recorded Frequencies - only trees with nonzero counts listed,  sorted high to low\n");
  FP "------------------------------------------------------------------------------------------------\n");
  if (modeloptions[ADDGHOSTPOP]==0)
  {
    FP "Tree_string\tpriorprob\tcount\tfreq_ALL\tfreq_set1\tfreq_set2\tppcp\n");
    for (i=0;i<numpoptopologies;i++) if (fatssarray[i].count > 0)
    {
     if (uniformprior == 0)
        temp1 = exp(topologypriors[fatssarray[i].origi]);
      FP "%s\t%.4f\t%d\t%.6f\t%.6f\t%.6f\t%.7f\n",fatssarray[i].treestr,temp1,fatssarray[i].count,fatssarray[i].freqall,fatssarray[i].freqset1,fatssarray[i].freqset2,fatssarray[i].ppcp);
    }
    FP "\n");
  }
  else
  {
    FP "Tree\tTree(w/ghost)\tpriorprob\tcount\tfreq_ALL\tfreq_set1\tfreq_set2\tppcp\n");
    for (i=0;i<numpoptopologies;i++) if (fatssarray[i].count > 0)
    {
      if (uniformprior == 0)
        temp1 = exp(topologypriors[fatssarray[i].origi]);
      FP "%s\t%s\t%.4f\t%d\t%.6f\t%.6f\t%.6f\t%.7f\n",fatssarray[i].treestrnoghost,fatssarray[i].treestr,temp1,fatssarray[i].count,fatssarray[i].freqall,fatssarray[i].freqset1,fatssarray[i].freqset2,fatssarray[i].ppcp);
    }
    FP "\n");
  }
  if (npops - modeloptions[ADDGHOSTPOP] <= 6)
    printallpoptreesamples(outfile,poptopologycounts,fatssarray,poptopologyproposedlist_rec,uniformprior);
  else
    FP "Full Population Tree Sample Counts and Frequencies (%d entries) - not printed\n\n",numpoptopologies);


/*  FP" Roberson-Foulds distances\n");
for (i=0;i<numpoptopologies;i++)
  FP"%d\n",RFtreedis[i]);
FP"\n"); */


  XFREE(fatssarray);
  for (i = 0; i < numuniquenodes[(npops - modeloptions[ADDGHOSTPOP])]; i++)
    XFREE(uniquenodes[i]);
  XFREE(uniquenodes);
  XFREE(nodecounts);
  XFREE(set1counts);
  XFREE(set2counts);
  
}//sort_and_print_alltreestrings

                               //callprintacceptancerates

// makes calls to asciitrend
void callasciitrend (FILE * outtofile,int trenddoublepoint,int trendspot)
{
  int i, li, ui;
  asciitrend (outtofile, lpgpd_v, trenddoublepoint, trendspot);

  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    asciitrend (outtofile, poptreeuinfo->v, trenddoublepoint, trendspot);
  
  for (i = 0; i < lastperiodnumber; i++)
    asciitrend (outtofile, T[i].v, trenddoublepoint, trendspot);
  if (hiddenoptions[SKIPMOSTUSCALAROUTPUT]==0)
  {
    if (nurates > 1 && runoptions[LOADRUN] == 0 && domutationscalarupdate)
    {
      for (li = 0; li < nloci; li++)
        for (ui = 0; ui < L[li].nlinked; ui++)
          asciitrend (outtofile, L[li].u_rec[ui].v, trenddoublepoint, trendspot);
    }
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY && runoptions[LOADRUN] == 0  && domutationscalarupdate)
        asciitrend (outtofile, L[li].kappa_rec->v, trenddoublepoint, trendspot);
  }

  
  return;
}                               // callasciitrend 


/* set up arrays pointing to information to put in curve plots,  then call asciicurve
   some things that are plotted are based on struct value_record and others on struct i_param
   this is why we cannot simply call asciicurve() with a pointer to a single type of structure  */
/* this is called only for the head node */
void callasciicurves (FILE *outfile,int mcmcrecords)
{
  struct plotpoint **curvexy;
  char **curvestr;
  int *curve_do_logplot;
  int numcurve = 0, nummigcurve = 0;
  int i, j, li, ui;
  int *nrecstep;
  int dotcurves = 0,dodemogcurves = 0, doucurves = 0, dohypercurves = 0,dohkycurves=0,do2Nmcurves = 0;
  struct plotpoint **tempxyarrays;
  char  **tempparamstr;

  // find out how many curves

/* use C[ARBCHAIN] 
  all chains have itheta and imig  and each element of each of these arrays has xy
  we need one to use for printing,  and it does not matter because the curves will be based on the saved genealogies and these came from the cold chain
*/ 
  //if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  numcurve = numsplittimes; // all run modes generate splittime curves
  dotcurves = 1;
  if (runmode != 4  &&  domutationscalarupdate && hiddenoptions[SKIPMOSTUSCALAROUTPUT]==0)
  {
    if (nurates > 1)
    {
      numcurve += nurates;
      doucurves = 1;
    }
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
      {
        numcurve++;
        dohkycurves = 1;
      }
  }
  if (runmode == POPTREEHYPERPRIORmode0 || runmode == GHYPERPRIORmode2 || runmode == HGHYPERPRIORmode5)
  {
    numcurve += nummigrateparams;
    numcurve += numpopsizeparams;
    dohypercurves = 1;
  }
  //if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])// 
   // numcurve += numpopsizeparams;
  if (runmode == Gmode3 || runmode == LOADGmode4 || runmode == HGmode6)
  {
    if (outputoptions[NOPOPMIGPARAMHIST] == 0 && hiddenoptions[PRINT2NMASCIICURVES] == 1)
      do2Nmcurves = 1;
    dodemogcurves = 1;
    numcurve += numpopsizeparams;
    for (i = 0,nummigcurve = 0; i < nummigrateparams; i++)
    {
      if ((modeloptions[EXPOMIGRATIONPRIOR]==1 && C[ARBCHAIN]->imig[i].pr.expomean > MPRIORMIN) ||
        (modeloptions[EXPOMIGRATIONPRIOR]==0 && C[ARBCHAIN]->imig[i].pr.max > MPRIORMIN))
      {
        numcurve++;
        nummigcurve++;
        if (do2Nmcurves)
          numcurve++;
      }
    }
  }
// allocate
  curvexy = static_cast<plotpoint **> 
                (malloc (numcurve * sizeof (struct plotpoint *)));
  curvestr = static_cast<char **> (malloc (numcurve * sizeof (char *)));
  nrecstep = static_cast<int *> (malloc (numcurve * sizeof (int)));
  curve_do_logplot = static_cast<int *> (malloc (numcurve * sizeof (int)));
// assign
  j = 0;
  if (dotcurves)
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      curvexy[j] = T[i].v->xy;
      curvestr[j] = &T[i].v->str[0];
      curve_do_logplot[j] = T[i].v->do_logplot;
      nrecstep[j] = mcmcrecords;
      j++;
    }
  }
  if (dodemogcurves)
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      /* curvexy for itheta and imig for C[ARBCHAIN] should have been filled by fillvec() from the call to printhistograms() when the histograms were made */
      curvexy[j] = C[ARBCHAIN]->itheta[i].xy;
      curvestr[j] = &C[ARBCHAIN]->itheta[i].str[0];
      curve_do_logplot[j] = 0;
      nrecstep[j] = 1;
      j++;
    }
    for (i = 0; i < nummigrateparams; i++)
      if ((modeloptions[EXPOMIGRATIONPRIOR]==1 && C[ARBCHAIN]->imig[i].pr.expomean > MPRIORMIN) ||
        (modeloptions[EXPOMIGRATIONPRIOR]==0 && C[ARBCHAIN]->imig[i].pr.max > MPRIORMIN))
      {
        /* curvexy for itheta and imig for C[ARBCHAIN] should have been filled by fillvec() from the call to printhistograms() when the histograms were made */
        curvexy[j] = C[ARBCHAIN]->imig[i].xy;
        curvestr[j] = &C[ARBCHAIN]->imig[i].str[0];
        curve_do_logplot[j] = 0;
        nrecstep[j] = 1;
        j++;
      }
  }
  if (do2Nmcurves)  
    /* added 1/11/2018,  bit of a kludge.  histograms have already been printed,  including those for 2Nm,  but 
      we want to print the curves for 2Nm. But we have not saved the values in the histograms.  Even though those values
      were generated and saved for q and m  by fillvec(). Anyway,  here we set up some new temporary 
      arrays and calculate the values for 2Nm again with a call to the new function fill_2Nm_vals()
    */ 
  {
    tempxyarrays = static_cast<plotpoint **> (malloc (nummigcurve * sizeof (struct plotpoint*)));
    tempparamstr = static_cast<char **> (malloc(nummigcurve *sizeof(char*)));
    for (i = 0; i < nummigcurve; i++)
    {
      tempxyarrays[i]=static_cast<plotpoint *> (calloc (GRIDSIZE, sizeof (struct plotpoint)));
      tempparamstr[i] = static_cast<char *> (malloc(PARAMSTRLEN *sizeof(char)));
    }
    fill_2Nm_vals (tempxyarrays,tempparamstr);
    for (i = 0; i < nummigcurve; i++)
    {
      curvexy[j] = tempxyarrays[i];
      curvestr[j] = &tempparamstr[i][0];
      curve_do_logplot[j] = 0;
      nrecstep[j] = 1;
      j++;
    }
  }
  if (doucurves)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
      {
        curvexy[j] = L[li].u_rec[ui].v->xy;
        curvestr[j] = &L[li].u_rec[ui].v->str[0];
        curve_do_logplot[j] = L[li].u_rec[ui].v->do_logplot;
        nrecstep[j] = mcmcrecords;
        j++;
      }
  }
  if (dohkycurves)
  {
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
      {
        curvexy[j] = L[li].kappa_rec->v->xy;
        curvestr[j] = &L[li].kappa_rec->v->str[0];
        curve_do_logplot[j] = L[li].kappa_rec->v->do_logplot;
        nrecstep[j] = mcmcrecords;
        j++;
      }
  }
  if (dohypercurves)
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      curvexy[j] = qh[i].v->xy;
      curvestr[j] = &qh[i].v->str[0];
      curve_do_logplot[j] = qh[i].v->do_logplot;
      nrecstep[j] = mcmcrecords;
      j++;
    }
    for (i = 0; i < nummigrateparams; i++)
    {
      curvexy[j] = mh[i].v->xy;
      curvestr[j] = &mh[i].v->str[0];
      curve_do_logplot[j] = mh[i].v->do_logplot;
      nrecstep[j] = mcmcrecords;
      j++;
    }
  }
  assert (numcurve == j);
  for (j = 0; j < numcurve; j++)
    asciicurve (outfile, curvexy[j], curvestr[j], curve_do_logplot[j],
                nrecstep[j]);
//free
  XFREE (curvexy);
  XFREE (curvestr);
  XFREE (curve_do_logplot);
  XFREE (nrecstep);
  if (do2Nmcurves) /* added 1/11/2018,  free the temporary arrays used for 2Nm values */ 
  {
    for (i = 0; i < nummigrateparams; i++) 
    {
      XFREE(tempxyarrays[i]);
      XFREE(tempparamstr[i]);
    }
    XFREE(tempxyarrays);
    XFREE(tempparamstr);
  }
}                               //callasciicurve 

void printsteps (FILE * outto, double like, double probg,int burndone, int burninsteps)
{
  if (!burndone)
  {
    if (outto != NULL)  fprintf(outto, "=BURNIN-PERIOD===============================\n");
    if (outto != NULL)  fprintf(outto, "STEP # %d  p(D|G): %.3lf p(G): %.3lf\n", step,
             like, probg);
  }
  else
  {
    if (outto != NULL)  fprintf(outto, "=============================================\n");
    if (runmode == POPTREEHYPERPRIORmode0 || runmode == POPTREEmode1)
    {
      if (outto != NULL)  fprintf(outto, "STEP # %d  Topologies Sampled: %d p(D|G): %.3lf p(G): %.3lf\n",runsteps,poptopologiessampled, like, probg);
    }
    else
    {
      if (genealogysamples > 0)
      {
        if (outto != NULL)  fprintf(outto,
                 "STEP # %d # Genealogies Saved: %d p(D|G): %.1lf p(G): %.1f\n",
                 runsteps, genealogysamples, like, probg);
      }
      else
      {
        if (outto != NULL)  fprintf(outto, "STEP # %d  p(D|G): %.3lf p(G): %.3lf\n",
                 runsteps, like, probg);
      }
    }
  }
  return;
}  //printsteps

