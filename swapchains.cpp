/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/*********** LOCAL STUFF **********/

/*
consider implementing the following to cut down on loops over allbetas[]
allbetapos[ci]  for the current node,  for chain ci  (ci < numchainspp)
allbetapos[ci] is the position in allbetas[]  of beta[ci]

ciforallbeta[i] for the current node, for position i in allbetas[], ciforallbeta[i] is the corresponding ci value,  the index for the matching beta on the current node
	if that beta is not on the current node,  the ciforallbeta[i] == -1
  */

static int **swapcount;

// these are only used if numprocesses > 1
static int swapdist;
static int **swapcount_rec; 

static int sizeofswapmatrix;
void initswapstuff(void);

/* local functions */
static double swapweight (int ci, int cj);
static double calcpartialswapweight(int c);
static double swapweight_bwprocesses(double sumi, double sumj, double betai, double betaj);
void clearswapinfo(void);
void swapbetas (int ci, int cj);


/* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */


/*
 3/5/2019  fixed bugs that arose when prior ratio heating was off 
 
 apparently there was confusion about the use of the prior under different circumstances when swapping chains

 for mc^3 the full mh ratio includes the priors for both chains in both the numerator and the denominator

 so if the priors are each raised to a power,  they completely cancel. 

 if they are raised to powers,  then they do not

 so unlike other update, in which the ratio of the priors is always present (raised to beta or not, but still present)
 for the chain swapping update,  the priors are completely missing if the betas are not being applied to the priors
 */
// beta[ci] - beta[cj] ;   likelihood[cj]-likelihood[ci]  and prior[cj] - prior[ci] if it is used 
double
swapweight (int ci, int cj)  
{
  int i;
  double sumi, sumj, w;
  double likelihoodratio,priorratio;
  double tempi,tempj;

  #ifdef TURNONCHECKS
    checkdetailedbalance_chainswap(C[ci]->allpcalc.pdg, C[cj]->allpcalc.pdg, 0.0,0.0, beta[ci],beta[cj]);
  #endif //TURNONCHECKS
  likelihoodratio = C[cj]->allpcalc.pdg - C[ci]->allpcalc.pdg;  // ci in the denominator 
  if (hiddenoptions[PRIORRATIOHEATINGOFF])  // prior ratios cancel out in this case,  even with hyperparameters
     w = (beta[ci] - beta[cj]) * likelihoodratio ;
  else
  {
    priorratio = 0.0;
   // if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]) JH 4/26/2018  not entirely clear what to do in this case,  as far as calculating marginal likelihood 
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR] == 1)
    {
      if (modeloptions[EXPOMIGRATIONPRIOR])
      {
        for (i=0,tempi = 0.0,tempj = 0.0;i<nummigrateparams;i++)
        {
          tempi += C[ci]->imig[i].pr.expomean;
          tempj += C[cj]->imig[i].pr.expomean;
        }
        priorratio  =  (tempi-tempj)/expo_m_mean; // ratio of the product of the exponential densities for the priors
      }
      else
      {
        for (i=0,tempi = 1.0,tempj = 1.0;i<nummigrateparams;i++)
        {
          tempi *= C[ci]->imig[i].pr.max;
          tempj *= C[cj]->imig[i].pr.max;
        }
        priorratio = log(tempi/tempj); // i.e. ratio of  product of the inverse of migration priors // i.e. log [ (1/max_cj)/(1/max_ci) ]
      }
      // now do popsize terms
      for (i=0,tempi = 1.0,tempj = 1.0;i<numpopsizeparams;i++)
      {
        tempi *= C[ci]->qhyperparams[(int) C[ci]->descendantpops[i]];
        tempj *= C[cj]->qhyperparams[(int) C[cj]->descendantpops[i]];
      }
      priorratio += log(tempi/tempj); 
    }
    else
      priorratio = 0.0;   // demographic priors if hyperpriors not being used 

    sumi = C[ci]->allpcalc.probg;
    sumj = C[cj]->allpcalc.probg;
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
    {
      sumi += C[ci]->allpcalc.probhgg;
      sumj += C[cj]->allpcalc.probhgg;
      if (usetopologypriors)
      {
        sumi += topologypriors[C[ci]->poptreenum];
        sumj += topologypriors[C[cj]->poptreenum];
      } 
    }
    priorratio += sumj - sumi;
    w =  ((beta[ci] - beta[cj]) * (likelihoodratio + priorratio));
#ifdef TURNONCHECKS
    checkdetailedbalance_chainswap(C[ci]->allpcalc.pdg, C[cj]->allpcalc.pdg, sumi,sumj, beta[ci],beta[cj]);
#endif //TURNONCHECKS
  }
  return (w);
} //swapweight

///AS: Adding a function to only calculate the sums for a particular chain
//this would then be shared with the swapper process/chain
// also see notes on swapweight() above 
void
calcpartialswapweight (int ci, double *likelihood, double *prior)
{
 int i;
	double priorsum,tempp;

  *likelihood = C[ci]->allpcalc.pdg;
  if (hiddenoptions[PRIORRATIOHEATINGOFF])  // do not include the prior ratio when the priors are not heated, in this case prior for each chain cancels itself out
    *prior = 0.0; 
  else  // prior ratios are heated and do not cancel out
  {
    priorsum = C[ci]->allpcalc.probg;
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
      priorsum += C[ci]->allpcalc.probhgg;
    if (usetopologypriors)
      priorsum += topologypriors[C[ci]->poptreenum];
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])  // need to include ratio of demographic priors,  but if hyperpriors are not being used then demographic priors cancel out
    {
      if (modeloptions[EXPOMIGRATIONPRIOR])
      {
        for (i=0,tempp = 0.0;i<nummigrateparams;i++)
        {
          tempp +=  C[ci]->imig[i].pr.expomean;
        }
        priorsum -= tempp/mprior; // multiple by the product of the exponential densities for the priors
      }
      else
      {
        for (i=0,tempp = 1.0;i<nummigrateparams;i++)
        {
          tempp *=  C[ci]->imig[i].pr.max;
        }
        priorsum -= log(tempp); // i.e. divide by the product of the migration priors 
      }
      // now do popsize terms
      for (i=0,tempp = 1.0;i<numpopsizeparams;i++)
      {
        tempp *= C[ci]->qhyperparams[(int) C[ci]->descendantpops[i]];
      }
      priorsum -= log(tempp);  // divide by the product of the theta priors 
    }
    *prior = priorsum;
  }
}   /* calcpartialswapweight() */

double
swapweight_bwprocesses(double likelihoodratio, double priorratio,double betai, double betaj)
{
  double w;
  if (hiddenoptions[PRIORRATIOHEATINGOFF])
    w = ((betai - betaj) * likelihoodratio); // prior ratios cancel out in this case,  even with hyperparameters
  else
    w = ((betai - betaj) * (likelihoodratio + priorratio));
	return(w);
} /* End of function swapweight_bwprocesses() */


void
swapbetas (int ci, int cj)   // used for swaps between chains on the same processor
{
	 double btemp;
  int ctemp = 0;
	 btemp = beta[ci];
	 beta[ci] = beta[cj];
	 beta[cj] = btemp;
  ctemp = C[ci]->currallbetapos;
  C[ci]->currallbetapos = C[cj]->currallbetapos;
  C[cj]->currallbetapos = ctemp;
	return;
}				/* swapbetas */


/************** global functions *****************/

#define SWAPDISTANCEDEFAULT 7 // default maximum distance in chain array of two chains with betas being swapped, not sure how much this matters
void initswapstuff(void)
{
  beta = static_cast<double *> (malloc ((numchainspp) * sizeof (double)));
  allbetas = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  sizeofswapmatrix = numchainstotal*numchainstotal;
  if (numchainstotal > 1)
  {

    swapcount = alt2d_alloc2Dint (numchainstotal,numchainstotal);
    if (numprocesses > 1)
    {
      swapcount_rec = alt2d_alloc2Dint (numchainstotal,numchainstotal); 
    }
    clearswapinfo();
  }
  swapdist = IMAX((int) SWAPDISTANCEDEFAULT,(int) numchainstotal*0.1);  // JH 4/30/2018,  changed from fixed value of SWAPDISTANCEDEFAULT  (7),  allowed for a wider range of chains 
}

void clearswapinfo(void)  // does not appear to be used
{
  int ci,cj;

  for (ci=0;ci<numchainstotal;ci++) for (cj=0;cj<numchainstotal;cj++)
  {
    swapcount[ci][cj] = 0;
    if (numprocesses > 1)
    {
      swapcount_rec[ci][cj] = 0;
    }
  }
}


/************ GLOBAL FUNCTIONS ******************/

void freeswapstuff(void)
{
  XFREE (allbetas);
  XFREE (beta);
  if (numchainstotal > 1)
  {
    alt2d_free2Dint( swapcount); 
    if (numprocesses > 1)
    {
      alt2d_free2Dint( swapcount_rec); 
    }
  }
}

int setswaptries(void)
{
   /* 
   4/17/2018

   keep in mind that the # of swaptries should not be a function of the # of cpus 
   did some playing around to see the affect of # of swaptries per step
   did not see a huge effect. 

   Assume that the number of swaptries go up with # of chains. So 
   we settle on going for one swap attempt per chain every 5 steps,  or 20% chance of swapping for each chain in each step. 

   When there are n chains,  the probability that a chain gets picked on a try is 2/n.  
   Let m be the # of swaptries

   so approximately  m x 2/n = 0.2  

   m = 0.1 x n  

   But this will be super high
   so set m = 0.1 x n  if n <= 400  and m = 0.05 x n for n > 400 

   
   
   previous old discussion,  mostly wrong because it worries about # of cpus:
   --------------------------------------------------------------------------
    Consider the # of swap attemps as a function of # of chains and # of nodes

     each chain can swap with any other that is within + or- swapdist positions from it
     if # chains is less than 2*swapdist then pretty much all chains can swap 
     the number of possible swaps is approx  numchainstotal * (2*swapdist)/2 = numchainstotal*swapdist

     target swap rate is to have each chain attempt a swap with every other that is within reach of it once every 10 steps 
      this is pretty arbitrary  - hope its a good happy medium between mixing and time spent swapping 

     then the total number of swap attempts per step would be 
      (numchainstotal * swapdist) / 10
     since SWAPDISTANCEDEFAULT == 10,  this gives numchainstotal swap attempts per step 
  
  BUT - turns out that this approach to putting a number on  numchainstotal swap attempts slows things down REALLY badly 
  
  4/19/2017  
  recalled that with multiple nodes, when we pick 2 chains, if neither selected chain is on the current node then no attempt is made 
  
  So a different approach is to have on average one acceptable swap attemp per chain  every m steps
    which means that we must calculate the probability that a swap attempt is valid (i.e. at least one of the picked chains is on the current node)
      so with n chains  and k nodes  (n/k per node) the probability p that no swap attempt is made will be 
      p = ((n - n/k)/n) ( (n - (n/k) - 1)/(n - 1))  
      As n->infinity this goes to ((k-1)/k)^2 which make sense 
      for various values of k:
      k  p
      2  0.25
      5  0.64
      10 0.81
      20 0.9025
      30 0.934
      40 0.9506
      50 0.9604
      100 0.9801 
      
      Suppose we attempt 1 swap per node 
      then the # of valid attempts per node is (1-p) and the total # of valid attempts is k(1-p) 

      Each valid attempt affects 2 chains. 

      so the number of chains engaged at a rate of 1 swap per node is 2k(1-p)

      so the proportion of chains  engaged per step is about (2k(1-p))/n 

      suppose we want each chain to be engaged on average once every m steps
      what is the scalar x that sets the proportion of chains engaged per step to be 1/m ? 
      x ((2k(1-p))/n) = 1/m
      x = n/(2 k m (1-p))
      so x would be the number of swap attempts per node 
      This turns out to be consistely close to n/(m 4),  almost regardless of k 

      e.g.  for m=1
      if n = 200 and k=40, then  1-p = 0.05  and 2 k (1-p) = 4 and x = n/4
      if n = 200 and k=10,  (1-p) approx 0.2  and 2 k (1-p) = 4  and x = n/4 
      
   
      if only 1 node then swaptries would be about  n/(2 m)

      We could consider swapdist effect:
        using swapdist,  there are range of possible swap partners 2*swapdist possible partners  out of n-1 
        this increases the chances of accepting a swap by a lot,  and we could increase m 
        so to keep the level of swapping the same,  when using swapdist, as compared to without
        we would divide the # of tries by (n-1)/(2*swapdist) 

        but the effect of swapdist depends strongly on how many chains and the heating parameters
        don't include it for now 
     
      no matter what, set a min of 2 and a max of 40 
      the upper bound of 40 is there because swaps take a lot of time
      if m=1 and n=500  then x=125,  so having a max of 40 has the effect of setting m to about 3.  
      if m=1 and n=1000, then x=250  and having a max of 40 has the effect of setting m to about 6
    */
  int st;
  int m = 1;  // number of steps on average between times each chain is involved in an attempted swap
  // 4/17/2018 JH redid setting of swaptries 
  /*double temp;
  
  if (numprocesses > 1)
  {
    temp = numchainstotal/(4.0 * m);
  }
  else
  {
    temp = numchainstotal/(2.0 * m);
  }
  st = IMIN(40,IMAX(2,INTEGERROUND(temp))); */
  if (numchainstotal <= 400)
    st = INTEGERROUND(numchainstotal * 0.1);
  else
    st = INTEGERROUND(numchainstotal * 0.05);
  return st;
}

#define CHAIN0FRAC 0.1  // 0.1 is best it seems.  deal with difficulty of using beta=0 in lowest chain when doing thermodynamic integration 
/* this use of CHAIN0FRAC is an awkward kludge,  but could not find a way to 
  have beta = 0.0 and get thermodynamic integration to work. 
  with CHAIN0FRAC == 0.1  the marginal likelihood does come out to be pretty consistent
  even as chain numbers vary,  so at least there is that.  */


void
setheat (double hval1, double hval2, int heatmode, int currentid)
{
  int ci;
  int x = 0;
  double h = 0.0;
  int no2;

  initswapstuff();

  if (numprocesses == 1 && numchainspp == 1) 
  {
      allbetas[0] = 1.0;
      beta[0] = 1.0;
      return;
  }
	for (int i = 0; i < numchainstotal; i++) 
 {
			switch (heatmode)
			{
			case HLINEAR:
				allbetas[x] = 1.0 / (1.0 + hval1 * i);
				x++;
				break;
			case HGEOMETRIC:
				allbetas[x] = 1 - (1 - hval2) * (i) * 
				pow (hval1, (double) (numchainstotal - 1 - (i))) / (double) (numchainstotal - 1);
				x++;
				break;
   case HFULL:   //JH added to deal with hidden genealogies and topology updating
    no2 = numchainstotal/2;
    if (i<= no2)
    {
	  			allbetas[x] = 1 - (0.5) * (i) * pow (hval1, (double) (no2 - i)) / (double) no2;
    }
    else
    {
      allbetas[x] = (0.5) * (numchainstotal - i - 1) * pow (hval1, (double) (i-no2)) / ((double) no2 - 1.0);
    }
    x++;
				break;
			case HEVEN: // default for calculating marginal likelihood // no longer default
      //h = 1.0 / (numchainstotal - 1);
      h = (1.0  - hval2) / (numchainstotal - 1); // 7/26/2019 changed so a lower bound can be used with EVEN 
						allbetas[x] = 1.0 - i * h;	
					x++;
					break;
			}
	}
  if (calcoptions[CALCMARGINALLIKELIHOOD]==1)
  {
    assert (nearlyequaldouble(allbetas[numchainstotal-1],0.0,1e-10)); //should be very close to 0
    allbetas[numchainstotal-1] = allbetas[numchainstotal-2]* CHAIN0FRAC;  // for some reason the last chain causes problems with doing thermodynamic integration
  }
  int i = 0;
  for (ci = currentid * numchainspp; ci < currentid * numchainspp + numchainspp; ci++)
  {
    beta[i] = allbetas[ci];
    if (beta[i] < 0.0 || beta[i] > 1.0)
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "command line heating terms have caused a heating value out of range. chain %d beta %lf", ci, beta[i]);
    i += 1;
  }

}                               /* setheat */


/*  ======= swapchains_bwprocesses  ========== 
  explanation: 
    only used with mpi
    complex function with a lot of tricky mpi work
    
    main outer loop is over swaptries
    for each attempt all nodes are tied up with implementing that attempt 
      - first the head node picks the betas
      - then the node(s) that have those betas do some work,  while all other nodes wait. 

  FOR EACH SWAP ATTEMPT:
  STEP 1. HEADNODE picks 2 heating values randomly but with indices that are within swapdist positions
          does this by picking two positions in allbetas[], which every cpus has and which never changes 
  STEP 2. Bcast those values to all the other nodes.
          After step 2,  all nodes know the indices in allbetas[]
  STEP 3. each node: check to see if either beta belongs to a chain on the current node, set doISswap  (1 if either is on current node. else 0) 
            keep in mind several possibilities:
              both betas are on the same node
                that node could be the current node, or not  
              else betas are different nodes
                one of those could be the current node
                or neither could be the current node 
        if so,  send the node number to HEADNODE 
        After step 3 the HEADNODE knows which two nodes have the chains corresponding to the picked indices in allbetas[]
        These passages of node numbers to the Headnode are PASS_TO_HEADNODE_TO_SAVE  operations
        
  STEP 4. broadcast from HEADNODE to others the node numbers for the beta values 
        After step 4,  all nodes know the node numbers corresponding to the picked beta values 

  if doISwap:  (one of the beta values is on the current node), then attempt a swap 
    STEP 5.  
		 	 send/receive between nodes for A, B   poptreestrings,  put in holdpoptreestring vars
			  set treematch (1 if A tree is same as B tree, else 0) 
					if nodes for A and B are the same (i.e. both on the current node): (STEP 5.1)
				    swapchains without MPI 
			  else  
				    if current node is the node for A (STEP 5.2.A )
					     calc likelihood and prior using calcpartialswapweight() for chain  A 
					     send/receive from node for B info on RY and NW updating 
					     send/receive from node for B info on beta vals and allbetapos and likelihoods and priors
					     increment swapcount below diagonal
          Do a series of SEND_RECV_SWAP operations
				    if current node is the node for B (STEP 5.2.B )
					     calc likelihood and prior using calcpartialswapweight() for chain  B 
					     send/receive from node for A info on RY and NW updating 
					     send/receive from node for A info on beta vals and allbetapos and likelihoods and priors			
          Do a series of SEND_RECV_SWAP operations
			
			     if current node is the node for A (STEP 5.3 )  (use node for A for actual metropolis hastings)
				      increment attempts
				      calculate metropolishastings ratio 
				      if success
            change beta, betapos
					       increment successes
					       increment swapcount above diagonal
					       increment chain0 swap count
				      set swapvar indicating update accepted or not
				      send swapvar to B 
			     else  (A was on another node, where swap attempt is made) (STEP 5.4 )
          receive swapvar from A
				      if swapvar 
					       change beta, betapos
					       set RY NW  update info
  else (i.e. doISwap == 0)  do nothing 	

 */

#ifdef MPI_ENABLED	 // only use this function with mpi, otherwise just use swapchains
void 
swapchains_bwprocesses(int currentid, int swaptries,int *numattemptwithin,int *numattemptbetween,int *numsuccesswithin,int *numsuccessbetween)
{
  int swapvar = 0;
	 MPI_Status status;

  int rc = 0; 
	 double metropolishastingsratio;
	 int whichElementA, whichElementB;
	 int procIdForA, procIdForB;
	 int doISwap, areWeA;
  double likelihoodi = 0.0;
  double likelihoodj = 0.0;
  double priori = 0.0;
  double priorj = 0.0;
  double abeta = 0.0;
  double bbeta = 0.0;
  int p = 0;
  int q = 0;
  int sa = 0;
  int sb = 0;
  int x = 0;
  int y = 0;

  /* poptreestrings involved in the swaps are compared so that we can count how many swaps update the topology */
  char holdpoptreestringA[POPTREESTRINGLENGTHMAX],holdpoptreestringB[POPTREESTRINGLENGTHMAX]; 
  int treematch;
  int a_allbetapos,b_allbetapos;
  extern MPI_Datatype MPI_updatescalar;
  struct updatescalarinfo RYsendinfo[MAXPOPS-1],RYrecvinfo[MAXPOPS-1],NWsendinfo[MAXPOPS-1],NWrecvinfo[MAXPOPS-1];
  int ti;

  *numattemptwithin = 0; *numattemptbetween = 0;*numsuccesswithin = 0;*numsuccessbetween = 0;
	 for (x = 0; x < swaptries; x++) //## swaptries 
  {
		  swapvar = 0;
    holdpoptreestringA[0] = 0;
    holdpoptreestringB[0] = 0;
    treematch = INT_MAX;
    
/**** STEP 1 *****/
      /*  only the head node picks the chains that swap */ 
		  if (currentid == HEADNODE)  //## pick sa and sb, the positions in allbetas of the beta values and send them to all the other processors
    {
			   sa = (int) (uniform() * numchainstotal);
      /* ## pick swap chain indices sa and sb  using the swapdist range,  sa and sb are indices in allbetas[] (not in betas[]) */
      int sbmin;
      int sbrange;
      if (numchainstotal < 2*swapdist + 3)
      {
        sbmin = 0;
        sbrange = numchainstotal;
      }
      else
      {
        sbmin = IMAX(0,sa-swapdist);
        sbrange = IMIN(numchainstotal, sa+swapdist) -sbmin;
      } 
      do
      {
        sb = sbmin + (int) (uniform () * sbrange);
      } while (sa == sb || sb < 0 || sb >= numchainstotal);

		  } //## end loop pick the positions in allbetas of the beta values 
/**** END STEP 1 *****/

/**** STEP 2 *****/
     /* use Bcast to send sa and sb,  should be faster than sends and receives  */
				rc = MPI_Bcast(&sa, 1, MPI_INT, HEADNODE, MPI_COMM_WORLD);
				if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
				rc = MPI_Bcast(&sb, 1, MPI_INT, HEADNODE, MPI_COMM_WORLD);
				if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc); 
/**** END STEP 2 *****/    

/**** STEP 3 *****/
		/* at this point all processors know which chain numbers were picked  */
  /* now identify whether the current node  has either of these chains,  and then send that information to HEADNODE */

  /* identify the beta values */
    abeta = allbetas[sa];
    bbeta = allbetas[sb];

    /* identify which processor the corresponding chains are on, and what chain numbers on those nodes they are  */
    whichElementA = 0; // will be the chain index (position in C[])  on whatever processor sa applies to 
    whichElementB = 0; // will be the chain index (position in C[])  on whatever processor sb applies to 
    doISwap = 0;
    areWeA = 0; // 1 if current node has A 
    procIdForA = UNDEFINEDINT;   //  is sa applies to the current node,  then procIdForA=currentid else procIdForA remains = -1
    procIdForB = UNDEFINEDINT;  //  is sb applies to the current node,  then procIdForB=currentid else procIdForB remains = -1

    /*  loop over chains on current node, see if abeta is the beta for one of the chains in the current processor
     if so, and this is not head node, send procIdForA to the head node */
    for (y = 0; y < numchainspp; y++)  
    {  
			   if (beta[y] == abeta) // can only be true if chain sa is on the current processor 
      {
        /* at this point y is the chain number on the current node for abeta, */
				    procIdForA = currentid;
				    doISwap = 1;
				    areWeA = 1;
				    whichElementA = y;  // the position in C[] on the current processor of chain sa
        strcpy(holdpoptreestringA,C[whichElementA]->chainpoptreestring);
				    if (currentid != HEADNODE)  // send the processor number for chain sa  to the head node 
        {
          rc = MPI_Send(&procIdForA, 1, MPI_INT, HEADNODE, 12345, MPI_COMM_WORLD);
          if (rc != MPI_SUCCESS)			MPI_Abort(MPI_COMM_WORLD, rc);
        }
			   }
		  }
    /* now do the same thing for bbeta */
		  for (y = 0; y < numchainspp; y++)  // see if bbeta is the beta for one of the chains in the current processor, if so and this is not head node, send procIdForB to the head node
    { 
			   if (beta[y] == bbeta)   // true if chain sb is on the current processor 
      {
        /* at this point y is the chain number on the current node for bbeta, */
				    procIdForB = currentid;
				    whichElementB = y;    // the position in C[] on the current processor of chain sb
				    doISwap = 1;
        strcpy(holdpoptreestringB,C[whichElementB]->chainpoptreestring);
				    if (currentid != HEADNODE) // send the processor number for chain sb  to processor 0 
        {
					      rc = MPI_Send(&procIdForB, 1, MPI_INT, 0, 1973, MPI_COMM_WORLD);  
					      if (rc != MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD, rc);
				    }
			   }
		  }
    /* corresponding Recv's for 12345,1973 */
    /* at this point doISwap == 1 if either sa or sb is on the current node */
		  if (currentid == HEADNODE && procIdForA == UNDEFINEDINT)  // receive from some other processor the processor number for which abeta is  beta[y]
    {  
      rc = MPI_Recv(&procIdForA, 1, MPI_INT, MPI_ANY_SOURCE, 12345, MPI_COMM_WORLD, &status);
			   if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
		  }
		  if (currentid == HEADNODE && procIdForB == UNDEFINEDINT)   // receive from some other processor the processor number for which bbeta is  beta[y]
    {
      rc = MPI_Recv(&procIdForB, 1, MPI_INT, MPI_ANY_SOURCE, 1973, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
		  }
/**** END STEP 3 *****/

/**** STEP 4 *****/
    /* at this point only HEADNODE and the processors holding chains sa and sb knows which processors hold chains sa and sb 
    use Bcast to broadcast the procid vals to all the nodes  - use bcast rather than sends and receives   */
				rc = MPI_Bcast(&procIdForA, 1, MPI_INT, HEADNODE, MPI_COMM_WORLD);
				if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Bcast(&procIdForB, 1, MPI_INT, HEADNODE, MPI_COMM_WORLD);
				if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
/**** END STEP 4 *****/    

	   if (doISwap == 1) // if either selected chain, A or B,  is on the current node,  we enter this loop
    {
/**** STEP 5 *****/
      if (procIdForA  != procIdForB) /* need to share poptreestrings between nodes for A and B  */
      {
         if (areWeA) // 1_26_2018  break these send/recv up 
        {
           rc = MPI_Send(C[whichElementA]->chainpoptreestring, POPTREESTRINGLENGTHMAX, MPI_CHAR, procIdForB, 26123, MPI_COMM_WORLD);// the corresponding MPI_Receive puts this in holdpoptreestringA[]
          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
          rc = MPI_Recv(holdpoptreestringB, POPTREESTRINGLENGTHMAX, MPI_CHAR, procIdForB, 26123, MPI_COMM_WORLD, &status);// the corresponding MPI_Receive puts this in holdpoptreestringA[]
          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
        }
        else
        {
          rc = MPI_Recv(holdpoptreestringA, POPTREESTRINGLENGTHMAX, MPI_CHAR, procIdForA, 26123, MPI_COMM_WORLD, &status);// the corresponding MPI_Receive puts this in holdpoptreestringA[]
          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
          rc = MPI_Send(C[whichElementB]->chainpoptreestring, POPTREESTRINGLENGTHMAX, MPI_CHAR, procIdForA, 26123, MPI_COMM_WORLD);// the corresponding MPI_Receive puts this in holdpoptreestringA[]
          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
        } 
      }
      assert (holdpoptreestringA != 0 && holdpoptreestringB != 0);
      treematch = strcmp(holdpoptreestringA,holdpoptreestringB)==0; // 1 if they match.  

	     if (procIdForA == procIdForB )//&& procIdForA == currentid) // same node,  must be current node, no need to swap between nodes
      {
        /* ==== STEP 5.1 ====*/
	    //	do within process swap(),  whichElementB and whichElementA are the chain numbers on the current node for which betas are swapped 
        assert (procIdForA == currentid);  // doISwap is only true if the current node is has A or B 
        *numattemptwithin += 1;
		      swapvar = swapchains(3,1,whichElementB,whichElementA);
        numsuccesswithin += swapvar;
        /* ==== END STEP 5.1 ====*/
	     } 
      else    // must be different nodes and one of them is the current node 
      {
        assert (procIdForA != procIdForB);
		      if (areWeA == 1) // currentnode is node A 
        {
          /* ==== STEP 5.2.A ====*/
			       abeta = beta[whichElementA];
          a_allbetapos = C[whichElementA]->currallbetapos;  // the position of abeta in the allbeta list 
          assert (allbetas[a_allbetapos] == abeta);
			       calcpartialswapweight(whichElementA,&likelihoodi,&priori);
          for (ti=0;ti<numsplittimes;ti++)  // have to also swap the current updating info RYupdate and NWupdate 
          {
            if (doRYupdate)
            {
              // send from A to B 
              RYsendinfo[ti] = C[whichElementA]->RYwidthinfo[ti];
              rc = MPI_Send(&RYsendinfo[ti],1,MPI_updatescalar,procIdForB,123+(ti+1), MPI_COMM_WORLD);
				          if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
              // receive from B into A 
              rc = MPI_Recv(&RYrecvinfo[ti],1,MPI_updatescalar,procIdForB,123+(ti+1), MPI_COMM_WORLD, &status);
				          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
            }
           if (doNWupdate)  // only do NW updates when not using hidden genealogies 
           {
              NWsendinfo[ti] = C[whichElementA]->NWwidthinfo[ti];
              rc = MPI_Send(&NWsendinfo[ti],1,MPI_updatescalar,procIdForB,323+(ti+1), MPI_COMM_WORLD);
				          if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
              rc = MPI_Recv(&NWrecvinfo[ti],1,MPI_updatescalar,procIdForB,323+(ti+1), MPI_COMM_WORLD, &status);
				          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
           }
          } 
          // send abeta to node for B  
				      rc = MPI_Send(&abeta, 1, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
          // receive bbeta from node for B 
				      rc = MPI_Recv(&bbeta, 1, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
           // send a_allbetapos to node for B
				      rc = MPI_Send(&a_allbetapos, 1, MPI_INT, procIdForB, 2, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);   
          // receive b_allbetapos from node for B 
				      rc = MPI_Recv(&b_allbetapos, 1, MPI_INT, procIdForB, 2, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
            // node for A sends lratio for A to node for B
				      rc = MPI_Send(&likelihoodi, 1, MPI_DOUBLE, procIdForB, 4, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
            //node for A receives lratio for B from node for B 
				      rc = MPI_Recv(&likelihoodj, 1, MPI_DOUBLE, procIdForB, 4, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)				MPI_Abort(MPI_COMM_WORLD, rc);
          // node for A sends priorratio for A to node for B
				      rc = MPI_Send(&priori, 1, MPI_DOUBLE, procIdForB, 6, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS)				MPI_Abort(MPI_COMM_WORLD, rc);
          // node for A receives priorratio for B  from node for B
				      rc = MPI_Recv(&priorj, 1, MPI_DOUBLE, procIdForB, 6, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)				MPI_Abort(MPI_COMM_WORLD, rc);

#ifdef TURNONCHECKS
			    for (p = 0; p < numchainstotal; p++) {
				    if (allbetas[p] == abeta) 
            {
              assert(p == a_allbetapos);
					    break;
				    }
			    }
			    for (q = 0; q < numchainstotal; q++) 
          {
				    if (allbetas[q] == bbeta) 
            {
              assert(q == b_allbetapos);
					    break;
				    }
			    } 
#endif //TURNONCHECKS 
          // only keep track of swapcount[][] on node A (it will get summed across nodes later )
			       p = a_allbetapos;
          q = b_allbetapos;
			       if (p < q) 
          {
				        swapcount[q][p]++;
			       } else 
          {
				        swapcount[p][q]++;
			       }
        /* ==== END STEP 5.2.A ====*/
		      }
		      if (areWeA != 1) // then currentid must be node for B 
        {
          /* ==== STEP 5.2.B ====*/
          assert (procIdForB == currentid);
			       bbeta = beta[whichElementB];
          b_allbetapos = C[whichElementB]->currallbetapos;
          assert (allbetas[b_allbetapos] == bbeta);
          calcpartialswapweight(whichElementB,&likelihoodj,&priorj);
          for (ti=0;ti<numsplittimes;ti++)
          {
            if (doRYupdate)
            {
              RYsendinfo[ti] = C[whichElementB]->RYwidthinfo[ti];
              rc = MPI_Recv(&RYrecvinfo[ti], 1,MPI_updatescalar,procIdForA,123+(ti+1), MPI_COMM_WORLD, &status);
				          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
              rc = MPI_Send(&RYsendinfo[ti],1,MPI_updatescalar,procIdForA,123+(ti+1), MPI_COMM_WORLD);
				          if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);

            }
            if (doNWupdate)  // only do NW updates when not using hidden genealogies 
            {
              NWsendinfo[ti] = C[whichElementB]->NWwidthinfo[ti];
              rc = MPI_Recv(&NWrecvinfo[ti], 1,MPI_updatescalar,procIdForA,323+(ti+1), MPI_COMM_WORLD, &status);
				          if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
              rc = MPI_Send(&NWsendinfo[ti],1,MPI_updatescalar,procIdForA,323+(ti+1), MPI_COMM_WORLD);
				          if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
            }
          }
				      rc = MPI_Recv(&abeta, 1, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)				MPI_Abort(MPI_COMM_WORLD, rc);
				      rc =  MPI_Send(&bbeta, 1, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
				      rc = MPI_Recv(&a_allbetapos, 1, MPI_INT, procIdForA, 2, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)			MPI_Abort(MPI_COMM_WORLD, rc);
				      rc =  MPI_Send(&b_allbetapos, 1, MPI_INT, procIdForA, 2, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS)			MPI_Abort(MPI_COMM_WORLD, rc);
				      rc = MPI_Recv(&likelihoodi, 1, MPI_DOUBLE, procIdForA, 4, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)			MPI_Abort(MPI_COMM_WORLD, rc);
				      rc = MPI_Send(&likelihoodj, 1, MPI_DOUBLE, procIdForA, 4, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS)			MPI_Abort(MPI_COMM_WORLD, rc);
						    rc = MPI_Recv(&priori, 1, MPI_DOUBLE, procIdForA, 6, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
				      rc = MPI_Send(&priorj, 1, MPI_DOUBLE, procIdForA, 6, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
          /* ==== END STEP 5.2.B ====*/
		      }
		      swapvar = 0;
		      if (areWeA == 1) // only keep track on the node for A 
        {
          /* ==== STEP 5.3 ====*/
          *numattemptbetween += 1;
#ifdef TURNONCHECKS
          checkdetailedbalance_chainswap(likelihoodi, likelihoodj, priori, priorj, abeta,bbeta);
#endif //TURNONCHECKS
          metropolishastingsratio = metropolishastingsratio = swapweight_bwprocesses(likelihoodj - likelihoodi, priorj - priori, abeta, bbeta);
          if (metropolishastingsdecide(metropolishastingsratio,1))
          {
            *numsuccessbetween  += 1;
				        beta[whichElementA] = bbeta;
            C[whichElementA]->currallbetapos = b_allbetapos;
				        if (p < q) 
            {
					         swapcount[p][q]++;
				        } 
            else 
            {
					         swapcount[q][p]++;
				        }
				        swapvar = 1;
            if ((abeta == 1.0 || bbeta == 1.0) && treematch ==0)  
            {
              chain0topolswaps +=1;  // only count this once per swap,  so count it when current node is for A chain 
              //printf("1 ");
            }
            for (ti=0;ti<numsplittimes;ti++)
            {
              if (doRYupdate)
                C[whichElementA]->RYwidthinfo[ti] = RYrecvinfo[ti];
              if (doNWupdate)
                C[whichElementA]->NWwidthinfo[ti] = NWrecvinfo[ti];
            } 
			       } 
          else 
          {
				        swapvar = 0;
			       }
			    // send swapvar to procIdForB
				      rc = MPI_Send(&swapvar, 1, MPI_INT, procIdForB, 1, MPI_COMM_WORLD);
				      if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
          /* ==== END STEP 5.3 ====*/
		      }
		      if (areWeA != 1) 
        {
          /* ==== STEP 5.4 ====*/
            // receive swapvar from procIdForA
				      rc = MPI_Recv(&swapvar, 1, MPI_INT, procIdForA, 1, MPI_COMM_WORLD, &status);
				      if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
			       if (swapvar != 0) 
          {
				        beta[whichElementB] = abeta;
            C[whichElementB]->currallbetapos = a_allbetapos;
            // added this,  why was it not here ?  5/3/2017
            for (ti=0;ti<numsplittimes;ti++)
            {
              if (doRYupdate)
                C[whichElementB]->RYwidthinfo[ti] = RYrecvinfo[ti];
              if (doNWupdate)
                C[whichElementB]->NWwidthinfo[ti] = NWrecvinfo[ti];
            } 
			       }
          /* ==== END STEP 5.4 ====*/
		      }
	     } //if (procIdForA != procIdForB) 
/**** END STEP 5 *****/
	   }


    #ifdef TURNONCHECKS
    if (areWeA == 1)
      pcheck(whichElementA,15);
    if (procIdForB == currentid)
      pcheck(whichElementB,16);
    #endif
	 } //## end loop swaptries 
}  //swapchains_bwprocesses
#endif  // MPI_ENABLED



/* 
  swaps between chains on the same processor  
  called either from qupdate or from swapchains_bwprocesses() 
  first value after nargs must be the # of swaptries
    after that there could be two integers for the specific chains to be swapped

  if called from qupdate, it will attempt multiple swaps, each time picking two chains within swapdist positions of each other
  if called from swapchains_bwprocesses()  the two chain indices will be passed into this function 
    in this case ci and cj are the chain numbers on the current node for which betas are swapped 
    and there will be just 1 attempt
*/

int
swapchains (int nargs, ...)    // used for swaps between chains on the same processor 
{
  int betai, betaj, i, swap0ok;
  int ci,cj;
  double metropolishastingsratio;
  
  
  //struct updatescalarinfo tempforswap;
  //int ti;
  int swapok;
  int swapattempts;
  int vai;
  va_list ap;
  va_start(ap, nargs);
  for (i=1;i<=nargs;i++)
  {
    vai = va_arg(ap,int);
    if (i==1)
      swapattempts = vai;
    if (i==2)
      ci = vai;
    if (i==3)
      cj = vai;
  }
  for (i = 0, swap0ok = 0,swapok = 0; i < swapattempts; i++)
  {
//printf("step %d attempt %d uni %.4lf",step,i,uniform());
    if (nargs < 2)  //ci and cj are not specified 
    {
      assert (numchainspp == numchainstotal);
			   betai = (int) (uniform() * numchainstotal);  // random chain index
      int cmin;
      int crange;
      if (numchainstotal < 2*swapdist + 3)
      {
        cmin = 0;
        crange = numchainstotal;
      }
      else
      {
        cmin = IMAX(0,betai-swapdist);
        crange = IMIN(numchainstotal, betai+swapdist) -cmin;
      } 
      do
      {
        betaj = cmin + (int) (uniform () * crange);
      } while (betai == betaj || betaj < 0 || betaj >= numchainstotal);  // 2nd random chain index different from first,  but not too far away */
	     for (ci = 0; ci < numchainstotal; ci++) 
      {
		      if (allbetas[betai] == beta[ci])
        {
			       break;
		      }
	     }
	     for (cj = 0; cj < numchainstotal; cj++) 
      { 
        if (allbetas[betaj] == beta[cj])
        {
          break;
        }
      }
      assert(betai == C[ci]->currallbetapos);
      assert(betaj == C[cj]->currallbetapos); 
	   }
    else
    {
      betai = C[ci]->currallbetapos;
      betaj = C[cj]->currallbetapos;
    }
	// record swaps between beta values, as they are listed in allbetas
  //   below the diagonal for attempts,  above the diagonal for successes 
    // record attempts 
    if (betai < betaj) 
    {
      swapcount[betaj][betai]++;
    } 
    else 
    {
      swapcount[betai][betaj]++;
    }	
    metropolishastingsratio = swapweight (ci, cj);
    if (metropolishastingsdecide(metropolishastingsratio,1))
    {
		    swapbetas(ci, cj);
  //TODO: These slideinfoscalers have to be swapped across processors as well
      /* turn these off   // why? 5/3/2017
      if (hiddenoptions[HIDDENGENEALOGY]==1) // also swap slide adjusters 
      {
        tempforswap = C[ci]->branchslideinfo;
        C[ci]->branchslideinfo = C[cj]->branchslideinfo;
        C[cj]->branchslideinfo = tempforswap;
      }
      for (ti=0;ti<numsplittimes;ti++)
      {
  #ifdef DO_RY1UPDATE
        tempforswap = C[ci]->RYwidthinfo[ti];
        C[ci]->RYwidthinfo[ti] = C[cj]->RYwidthinfo[ti];
        C[cj]->RYwidthinfo[ti] = tempforswap;
  #endif
  #ifdef DO_NWUPDATE  
        if (hiddenoptions[HIDDENGENEALOGY] == 0)  // only do NW updates when not using hidden genealogies 
        {
          tempforswap = C[ci]->NWwidthinfo[ti];
          C[ci]->NWwidthinfo[ti] = C[cj]->NWwidthinfo[ti];
          C[cj]->NWwidthinfo[ti] = tempforswap;
        }
  #endif
      } */
      if (C[ci]->poptreenum != C[cj]->poptreenum)
      {
        assert(strcmp(C[ci]->chainpoptreestring,C[cj]->chainpoptreestring)!=0);
        if (beta[ci] == 1.0 || beta[cj] == 1.0)
        {
          chain0topolswaps +=1;
        }
      }
      else
        assert(strcmp(C[ci]->chainpoptreestring,C[cj]->chainpoptreestring)==0);
	    if (betai < betaj) 
     {
		    swapcount[betai][betaj]++;
	    } else 
     {
		    swapcount[betaj][betai]++;
	    }
      if (ci == 0 || cj == 0)
        swap0ok |= 1;
      swapok += 1;
    }
    #ifdef TURNONCHECKS
    pcheck(ci,13);
    pcheck(cj,14);
    #endif
 }
  //return swap0ok;
#ifdef TURNONCHECKS
 check_hgprob_sums(ci);
 check_hgprob_sums(cj);
#endif //TURNONCHECKS
  va_end(ap);
  return swapok;
}                               /* swapchains */



/**********  writes to file  - this is the only function in swapchains.cpp  that does any file I/O ****************/
void
printchaininfo (FILE * outto, int heatmode, double hval1,
                double hval2, int currentid)
{
  int i;
#ifdef MPI_ENABLED
  int rc = 0; //AS: return code for MPI C bindings
  //MPI_Status status;
#endif
  if (currentid == HEADNODE && outto != NULL)
  {

    fprintf(outto, "\nHeated Chain Swapping\n");
    fprintf(outto,   "---------------------\n");
  }

#ifdef MPI_ENABLED
	if (numprocesses > 1) 
  {
    MPI_Barrier(MPI_COMM_WORLD);  // would sometimes get hangs here without putting in the barrier 
    rc = MPI_Reduce(swapcount[0],swapcount_rec[0], sizeofswapmatrix, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  // should be the same sd  rc = MPI_Reduce(&(swapcount[0][0]),&(swapcount_rec[0][0]), sizeofswapmatrix, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	}
#endif
  if (outto == stdout)  // print to screen if currentid == HEADNODE
  {
    if (currentid == HEADNODE) 
    {
      if (outto != NULL)  fprintf(outto, "Temp1     Temp2    #Swaps    #Attempts   Rate\n");
	     for (i = 0; i < numchainstotal - 1; i++) 
      {
#ifdef MPI_ENABLED
        if (numprocesses > 1) 
        {
  			     if (/*i < j && */swapcount_rec[i+1][i] > 0) 
          {
				        if (outto != NULL)  fprintf(outto, " %.6f    %.6f   %5d   %5d   %7.4f\n", allbetas[i], allbetas[i+1], swapcount_rec[i][i+1], swapcount_rec[i+1][i], (float) swapcount_rec[i][i+1]/ (float) swapcount_rec[i+1][i]);
			       }
			       else if (/*i < j && */swapcount_rec[i+1][i] == 0) 
          {
				        if (outto != NULL)  fprintf(outto, " %.6f    %.6f   %5d   %5d   na\n", allbetas[i], allbetas[i+1], swapcount_rec[i][i+1], swapcount_rec[i+1][i]);
			       }
	       } 
#endif
        if (numprocesses == 1) 
        {
    	     if (/*i < j && */swapcount[i+1][i] > 0) 
          {
				        if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   %7.4f\n", allbetas[i], allbetas[i+1], swapcount[i][i+1], swapcount[i+1][i], (float) swapcount[i][i+1]/ (float) swapcount[i+1][i]);
			       }
			       else if (/*i < j && */swapcount[i+1][i] == 0) 
          {
				        if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   na\n", allbetas[i], allbetas[i+1], swapcount[i][i+1], swapcount[i+1][i]);
			       }
        }
	     }
	     if (outto != NULL)  fprintf(outto, "\n\n");
	   }
  }
  else  // printing to a file 
  {
    if (numprocesses == 1 && currentid == HEADNODE)  // only 1 node 
    {
	     if (outto != NULL)  fprintf(outto, "Temp1     Temp2    #Swaps    #Attempts   Rate\n");	
      for (i = 0; i < numchainspp - 1; i++) 
      {
    	   if (/*i < j && */swapcount[i+1][i] > 0) 
        {
				      if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   %7.4f\n", allbetas[i], allbetas[i+1], swapcount[i][i+1], swapcount[i+1][i], (float) swapcount[i][i+1]/ (float) swapcount[i+1][i]);
			     }
			     else if (/*i < j && */swapcount[i+1][i] == 0) 
        {
				      if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   na\n", allbetas[i], allbetas[i+1], swapcount[i][i+1], swapcount[i+1][i]);
			     }
	     }
    }
    if (numprocesses > 1 && currentid == HEADNODE) 
    {
	     if (outto != NULL)  fprintf(outto, "Temp1     Temp2    #Swaps    #Attempts   Rate\n");
	     for (i = 0; i < numchainstotal - 1; i++) 
      {
        if (/*i < j && */swapcount_rec[i+1][i] > 0) 
        {
	         if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   %7.4f\n", allbetas[i], allbetas[i+1], swapcount_rec[i][i+1], swapcount_rec[i+1][i], (float) swapcount_rec[i][i+1]/ (float) swapcount_rec[i+1][i]);
        }
        else if (/*i < j && */swapcount_rec[i+1][i] == 0) 
        {
	         if (outto != NULL)  fprintf(outto, " %7.4f    %7.4f   %5d   %5d   na\n", allbetas[i], allbetas[i+1], swapcount_rec[i][i+1], swapcount_rec[i+1][i]);
        }
	     }
	     if (outto != NULL)  
        fprintf(outto, "\n\n");
	   }  
	 }
}                               /* printchaininfo */
