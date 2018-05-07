/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* for debugging - do not include in release */

/*

void chaininfo_print(int currentid,int recordint)

used for debugging swapping

*/


#undef GLOBVARS
#include "ima.hpp"

#ifdef TURNONCHECKS

/*********** LOCAL STUFF **********/

#define LOCALMAXCHAINS  1000
struct chaindebugstruct{
  double betaval;
  int poptreenum;
//  int processornum;
//  int chainnum;
} cdarray[LOCALMAXCHAINS];

static void cdarraysort ();
static void getchaininfo_debug(int currentid);
static void output_chaininfo_poptreenum (int currentid);
static void output_chaininfo_betaval (int currentid);


/******* LOCAL FUNCTIONS ***********/

void
cdarraysort ()
{
  unsigned long i, ir, j, l;
  struct chaindebugstruct ct;
  int n = numchainstotal;
  struct chaindebugstruct *lptr = &cdarray[0] - 1;
  if (n < 2)
    return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
    {
      ct = *(lptr + --l);
    }
    else
    {
      ct = *(lptr + ir);
      *(lptr + ir) = *(lptr + 1);
      if (--ir == 1)

      {
        *(lptr + 1) = ct;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir)
    {
      if (j < ir && (lptr + j)->betaval < (lptr + (j + 1))->betaval)
        j++;
      if (ct.betaval < (lptr + j)->betaval)
      {
        *(lptr + i) = *(lptr + j);
        i = j;
        j <<= 1;
      }
      else
      {
        j = ir + 1;
      }
    }
    *(lptr + i) = ct;
  }
}                               /*cdarraysort */

void getchaininfo_debug(int currentid)
{
  int cipp,cit;
  for (cipp=0;cipp<numchainspp;cipp++)
  {
    cit = currentid*numchainspp + cipp;
    cdarray[cit].betaval = beta[cipp];
    cdarray[cit].poptreenum = C[cipp]->poptreenum;
  }
}
void output_chaininfo_poptreenum (int currentid)
{
  /* send and receive cdarray values
  sort the cdarray
  print the cdarray */
  int rc;
  int ci,cit,np;
  FILE *chaininfo_poptreenum_file;
  static char chaininfo_poptreenum_filename[32] = "debug_chaininfo_poptreenum.out";
#ifdef MPI_ENABLED
	MPI_Status status;
#endif

#ifdef MPI_ENABLED
  if (currentid != 0)
  {
    for (ci=0;ci<numchainspp;ci++)
    {
      cit = currentid*numchainspp + ci;
      rc = MPI_Send(&cdarray[cit].betaval, 1, MPI_DOUBLE, 0, cit*1237, MPI_COMM_WORLD);
	    if (rc != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, rc);
      rc = MPI_Send(&cdarray[cit].poptreenum, 1, MPI_INT, 0, cit*7891, MPI_COMM_WORLD);
	    if (rc != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
  }
  else
  {
    for (np = 1; np < numprocesses;np++)
    {
      for (ci=0;ci<numchainspp;ci++)
      {
        cit = np*numchainspp + ci;
        rc = MPI_Recv(&cdarray[cit].betaval, 1, MPI_DOUBLE,np /* MPI_ANY_SOURCE */, cit*1237, MPI_COMM_WORLD, &status);
  	    if (rc != MPI_SUCCESS)
	  	    MPI_Abort(MPI_COMM_WORLD, rc);
        rc = MPI_Recv(&cdarray[cit].poptreenum, 1, MPI_INT,np /* MPI_ANY_SOURCE */, cit*7891, MPI_COMM_WORLD, &status);
	      if (rc != MPI_SUCCESS)
  		    MPI_Abort(MPI_COMM_WORLD, rc);
      }
    }
  }
#endif

  if (currentid == HEADNODE)
  {
    chaininfo_poptreenum_file = fopen (chaininfo_poptreenum_filename, "a");
    cdarraysort();
    fprintf(chaininfo_poptreenum_file,"%d\t",step);
    for (ci=numchainstotal-1;ci>=0;ci--)
      fprintf(chaininfo_poptreenum_file,"%d\t",cdarray[ci].poptreenum);
    fprintf(chaininfo_poptreenum_file,"\n");
    FCLOSE(chaininfo_poptreenum_file);
  }
#ifdef MPI_ENABLED
	MPI_Barrier(MPI_COMM_WORLD);// is this needed here?
#endif

  return;
 } /* output_chaininfo_poptreenum */


void output_chaininfo_betaval (int currentid)
{
  /* send and receive cdarray values
  sort the cdarray
  print the cdarray */
  int rc;
  int ci;
  FILE *chaininfo_betaval_file;
  static char chaininfo_betaval_filename[32] = "debug_chaininfo_betaval.out";
  double betavals_rec[LOCALMAXCHAINS],*brp;
  int pos;
#ifdef MPI_ENABLED
	//MPI_Status status;
#endif

#ifdef MPI_ENABLED
  if (numprocesses > 1)
  {
    pos = currentid *numchainspp;
    brp = &betavals_rec[pos];
    //if (currentid != 0)
    //{
      rc = MPI_Gather (beta,numchainspp,MPI_DOUBLE,brp,numchainspp,MPI_DOUBLE,0, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, rc);
    /*}
    else
    {
      for (ci = 0;ci<numchainspp;ci++)
        betavals_rec[pos+ci] = beta[ci];
    } */
  	//MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  if (currentid == HEADNODE)
  {
    chaininfo_betaval_file = fopen (chaininfo_betaval_filename, "a");
    fprintf(chaininfo_betaval_file,"%d\t",step);
    for (ci=0;ci<numchainstotal;ci++)
    {
      if (numprocesses > 1)
        fprintf(chaininfo_betaval_file,"%.5f\t",betavals_rec[ci]);
      else
        fprintf(chaininfo_betaval_file,"%.5f\t",beta[ci]);
    }
    fprintf(chaininfo_betaval_file,"\n");
    FCLOSE(chaininfo_betaval_file);
  }
  return;
 } /* output_chaininfo_betaval */



/********** GLOBAL functions ********/
/*
writes to two files

"debug_chaininfo_poptreenum.out"
for every time the program records,  it prints a line to that file
the line contains the current poptree number that is held by each chain,  order by high beta value to low


"debug_chaininfo_betaval.out"
for every time the program records, writes a line to this file
each line has the beta value for each chain, in order. useful to see how the beta values move around

useful for debugging swapping, particularly with multiple processors
*/

void chaininfo_print(int currentid,int recordint)
{
 static int debuginit = 0;
 static int debugi = 0;
  if (debuginit == 0)
  {
    debugi = recordint;
    debuginit = 1;
  }
  if (debugi == recordint)
  {
    getchaininfo_debug(currentid);
    output_chaininfo_poptreenum (currentid);
    output_chaininfo_betaval (currentid);
    debugi = 1;
  }
  else
  {
    debugi++;
  }
}
#endif //TURNONCHECKS
