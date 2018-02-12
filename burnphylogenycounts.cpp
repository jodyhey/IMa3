/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* added this file on 1/31/2017
  it records and generates output on the phylogeny distribution during the run 
  the idea is that buy watching that distribution you can tell if the chain has reached stationarity */

#undef GLOBVARS
#include "ima.hpp"


int *burnpoptopologycounts,*totalburnpoptopologycounts;
int lastburnarray;



/* the full 2D array has NUMPCOUNTARRAYS columns and (numpoptopologies+1) rows
  we handle this as a single block of ints 

  if it were a 2D array we'd refer to position in row i, col j as a[i][j]
  but instead we essentially have one long list 
  let the first NUMPCOUNTARRAYS values be in row 0 (i.e. for topology #0)
  and the next NUMPCOUNTARRAYS values be in row 1 (i.e. for topology #1)

  so to get  topology i, array j  (row i and col j) we go to position i*NUMPCOUNTARRAYS  + j

  If we want to move a value from column j to j + 1 but keep it in the same row
  then it would go from i*NUMPCOUNTARRAYS  + j  to   i*NUMPCOUNTARRAYS  + j + 1
  */

void init_burn_phylogeny_counts()
{
  int i, asize = NUMPCOUNTARRAYS * (numpoptopologies+1);
  lastburnarray = 0;
  burnpoptopologycounts = static_cast<int *> (malloc ((size_t) asize * sizeof (int)));
  for (i=0;i<=numpoptopologies;i++)  // initialize column 0 
    burnpoptopologycounts[i * NUMPCOUNTARRAYS ] = 0;
if (numprocesses > 1)
  totalburnpoptopologycounts = static_cast<int *> (malloc ((size_t) asize * sizeof (int))); 
}

void free_burn_phylogeny_counts()
{
  XFREE (burnpoptopologycounts);
if (numprocesses > 1)
  XFREE (totalburnpoptopologycounts);

}

void recordburntopology(void)
{
  int z = whichiscoldchain();
  if (z >= 0) 
  {
    burnpoptopologycounts[C[z]->poptreenum * NUMPCOUNTARRAYS ] += 1;
    burnpoptopologycounts[numpoptopologies * NUMPCOUNTARRAYS ] += 1;
    //*(burnpoptopologycounts  + C[z]->poptreenum) += 1;
    //*(burnpoptopologycounts  + numpoptopologies) += 1;
  }
}

void outputburntopologycounts(FILE *burnfile, int currentid)  
{
  int i,j,asize = NUMPCOUNTARRAYS * (numpoptopologies+1);
  int rc;
  if (numprocesses > 1)
  {
#ifdef MPI_ENABLED
    rc = MPI_Reduce(burnpoptopologycounts, totalburnpoptopologycounts, asize,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS)
		    MPI_Abort(MPI_COMM_WORLD, rc);
#endif
  }
  else
  {
    // if only 1 process totalburnpoptopologycounts has no allocation but is just a pointer
    totalburnpoptopologycounts = burnpoptopologycounts; 
  }
  if (currentid == HEADNODE) 
  {
    fprintf (burnfile, "========================================\n\n");

    fprintf (burnfile,"\nPopulation Names\n");
    fprintf (burnfile,"----------------\n");
    for (i = 0; i < npops; i++)
    {
      if (i==npops-1 && modeloptions[ADDGHOSTPOP]==1)
        fprintf (burnfile,"Population %d : ghost \n", i);
      else
        fprintf (burnfile,"Population %d : %s \n", i, popnames[i]);
    }
    fprintf (burnfile,"\n");
    if (numpoptopologies <= 180)
    {
      fprintf (burnfile, "Phylogeny Distribution (%d unique topologies with %d populations)\n   First column = values from the most recent burnin interval, second column = previous interval, etc etc.\n",numpoptopologies,npops);
      for (i=0;i<numpoptopologies;i++)
      {
        fprintf(burnfile,"%d\t%s\t",i,alltreestrings[i]);
        for (j=0;j<=lastburnarray;j++)
        {
          //to get  topology i, array j  (row i and col j) we go to position i*NUMPCOUNTARRAYS  + j
          fprintf(burnfile,"%.4lf\t",totalburnpoptopologycounts[ i * NUMPCOUNTARRAYS + j]/(double) totalburnpoptopologycounts[ NUMPCOUNTARRAYS * numpoptopologies + j]);
        }
        fprintf(burnfile,"\n");
      }
    }
    else
    {
      fprintf (burnfile, "Phylogeny Distribution (%d unique topologies with %d populations, unobserved trees not listed)\n   First column = values from the most recent burnin interval, second column = previous interval, etc etc.\n",numpoptopologies,npops);
      char tempstr[200] = {'\0'};   
      for (i=0;i<numpoptopologies;i++)
      {
        tempstr[0] = '\0';
        int tsi = 0;
        int nonzero = 0;
        for (j=0;j<=lastburnarray;j++)
        {
          //to get  topology i, array j  (row i and col j) we go to position i*NUMPCOUNTARRAYS  + j
          tsi += sprintf(&tempstr[tsi],"%.4lf\t",totalburnpoptopologycounts[ i * NUMPCOUNTARRAYS + j]/(double) totalburnpoptopologycounts[ NUMPCOUNTARRAYS * numpoptopologies + j]);
          if (totalburnpoptopologycounts[ i * NUMPCOUNTARRAYS + j] > 0)
            nonzero = 1;
        }
        if (nonzero==1)
        {
          fprintf(burnfile,"%d\t%s\t",i,alltreestrings[i]);
          fprintf(burnfile,"%s\n",tempstr);
        }
      }
    }
    fprintf(burnfile,"Total\t");
    for (j=0;j<=lastburnarray;j++)
      fprintf(burnfile,"%6d\t",totalburnpoptopologycounts[ NUMPCOUNTARRAYS * numpoptopologies + j]);
    fprintf(burnfile,"\n\n");

  }
  for (j=IMIN(lastburnarray,NUMPCOUNTARRAYS-2);j>=0;j--)
    for (i=0;i<=numpoptopologies;i++)
      burnpoptopologycounts[ i* NUMPCOUNTARRAYS + j + 1] = burnpoptopologycounts[ i * NUMPCOUNTARRAYS + j];
  if (lastburnarray < NUMPCOUNTARRAYS-1)
    lastburnarray += 1;
  for (i=0;i<=numpoptopologies;i++)
    burnpoptopologycounts[i * NUMPCOUNTARRAYS ] = 0;
}

