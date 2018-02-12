/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/* code for updating migration rate priors from the hyper prior */

extern int numdistinctpopulationpairs[]; 


/*********** local to this file  ***************/

static struct probcalc localholdallpcalc;


/******* local prototypes ********/
double getnewpriorval(double hyperparam, double oldval,int morq);  //morq == 1 if migration parameter, 0 otherwise
void freekrlinkedlist(struct dictionary_node_kr* head);


/******* local functions ********/

double getnewpriorval(double hyperparam, double oldval,int morq)  //morq == 1 if migration parameter, 0 otherwise
{
  double U, newval, windowsize,priormax;
 
  if (morq==1 && modeloptions[EXPOMIGRATIONPRIOR])
  {
    windowsize = hyperparam;//10.0;
    priormax = EXPOMIGPLOTSCALE*hyperparam;
  }
  else
  {
    windowsize = hyperparam/5.0;//10.0;
    priormax = hyperparam;
  }

  U=uniform();
  newval = (oldval - windowsize/2.0) + U*windowsize;
  while ((newval > priormax) || (newval < 0.0))
	 {
    if (newval > priormax)
			   newval = 2.0 * priormax - newval;
		  else if (newval < 0.0)
				  newval =  - newval;
	} 
  if (newval < MINPRIORFROMHYPERPRIOR)
    return MINPRIORFROMHYPERPRIOR; 
  assert (newval > 0);
  return newval;  
}  /* getnewpriorval */


/* added this to free up the linked lists that are created in the K&R dictionary 
need to free individual nodes,  and the linked lists descending from each node
also for each node that is not null, must free name 
*/ 
void freekrlinkedlist(struct dictionary_node_kr* head)
{
  struct dictionary_node_kr *tmphead,*tmp;
  if (head != NULL)
  {
    XFREE(head->name);
    tmphead = head->next;
    while (tmphead != NULL)
    {
       tmp = tmphead;
       tmphead = tmphead->next;
       XFREE(tmp->name);
       XFREE(tmp);
    } 
  }
  XFREE(head); // must be freed even if null because it was malloced and initialized null when it was created 
}




/******* global functions ********/

int update_migration_prior_intree(int ci, int mi) // 0<mi<nummigrateparams
{
  struct priorvalues holdpr;
  double newmpriorval;
  int accp;
  double metropolishastingsratio;
  double priorratio,proposalratio;

  copy_probcalc (&localholdallpcalc, &C[ci]->allpcalc);
  
  if (modeloptions[EXPOMIGRATIONPRIOR])
    // hastings ratio when hyperprior_expo_m_mean is the mean of the exponential hyperprior 
  { // log of the hastings ratio of two exponentials from the same density
    holdpr.expomean = C[ci]->imig[mi].pr.expomean;
    newmpriorval = getnewpriorval(expo_m_mean,holdpr.expomean,1);
    proposalratio = 0.0;
    priorratio = (holdpr.expomean - newmpriorval)/expo_m_mean;
    //priorratio =  -2*proposalratio;// yujin's bug fix -   count it twice because we are updating two migration rate priors
    C[ci]->imig[mi].pr.expomean = newmpriorval;
  }
  else  // uniform
  {
    holdpr.max =C[ci]->imig[mi].pr.max;
    holdpr.min = 0.0;
    newmpriorval = getnewpriorval(m_max,holdpr.max,1);
    C[ci]->imig[mi].pr.max = newmpriorval;
    
    priorratio = -log(newmpriorval/holdpr.max); 
    proposalratio = 0.0;
  }
  // initialize_integrate_tree_prob is overkill because most priors have not changed
  // could write a new function that just does integrations for the terms with new priors 
  initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  priorratio += C[ci]->allpcalc.probg - localholdallpcalc.probg;
  
  if (calcoptions[CALCMARGINALLIKELIHOOD]) 
  {
    metropolishastingsratio = priorratio + proposalratio;
  }
  else
  {
    metropolishastingsratio = beta[ci] * priorratio + proposalratio;
  }
  if (metropolishastingsdecide(metropolishastingsratio,1))
  {
    accp = 1;
    if (modeloptions[EXPOMIGRATIONPRIOR])
    {
      struct dictionary_node_kr *temp;
      if (C[ci]->imig[mi].dir==0)
        temp = dictionary_install(C[ci]->imig[mi].descstr,C[ci]->imig[mi].pr.expomean,C[ci]->mltorhpriors);
      else
        temp = dictionary_install(C[ci]->imig[mi].descstr,C[ci]->imig[mi].pr.expomean,C[ci]->mrtolhpriors);
    }
    else
    {
      struct dictionary_node_kr *temp;
      if (C[ci]->imig[mi].dir==0)
        temp = dictionary_install(C[ci]->imig[mi].descstr,C[ci]->imig[mi].pr.max,C[ci]->mltorhpriors);
      else
        temp = dictionary_install(C[ci]->imig[mi].descstr,C[ci]->imig[mi].pr.max,C[ci]->mrtolhpriors);
    }
  }
  else
  {
//if (mii+dir == 0)   
  //printf("reject \n");
    if (modeloptions[EXPOMIGRATIONPRIOR])
      C[ci]->imig[mi].pr.expomean = holdpr.expomean;
    else
      C[ci]->imig[mi].pr.max = holdpr.max;
    copy_probcalc (&C[ci]->allpcalc, &localholdallpcalc);
    accp = 0;
  }
  return accp;
}  //update_migration_prior_intree


void update_migration_prior_not_intree(int ci, int *attempted, int *accepted)
{
  double newpriorval,oldpriorval;
  double metropolishastingsratio;
  double priorratio,proposalratio;
  static int numattempts;
  static int currenti = -1;
  int i,j,found,accp;
  static int numnotintree[]={0,0,0,5,22, 84, 291, 951,9302};
  static int dir = -1;

  if (dir==-1)
    dir = 0;
  else
    dir = 1-dir;  //switch back and forth with alternate calls. 

 /* If we update 1 term for 3 pops  each time we come through,  and we want the same rate of update
  for more npops  then we numnotintree/5 attempts each time we come thru  */

  if (currenti == -1)
  {
    numattempts = numnotintree[npops]/5;
    currenti = 0;
  }
  i = 0;
  accp = 0;
  while (i < numattempts)
  {
    if (strlen(poppairs[currenti]) > 3) // only check pairs involving ancestral populations 
    {
      found = 0;
      for (j=0;j<nummigrateparampairs;j++)
      {
        if (strcmp(C[ci]->imig[2*j].descstr, poppairs[currenti])==0)  // 2*j because there are two imig[] values that have the same descstr
        {
          found = 1;
          break;
        } 
      }
      if (found == 0) // attempt update, poppairs[currenti] is not in the current tree 
      {
        if (dir==0)
          //oldpriorval = C[ci]->mltorhpriors[currenti];
          oldpriorval = getvalue(poppairs[currenti],C[ci]->mltorhpriors);  // mltorhpriors are for dir==0
        else
          //oldpriorval = C[ci]->mrtolhpriors[currenti];
          oldpriorval = getvalue(poppairs[currenti],C[ci]->mrtolhpriors); // mrtolhpriors are for dir==1

        if (modeloptions[EXPOMIGRATIONPRIOR])
          // hastings ratio when hyperprior_expo_m_mean is the mean of the exponential hyperprior 
        { // log of the hastings ratio of two exponentials from the same density
          newpriorval = getnewpriorval(expo_m_mean,oldpriorval,1);
          proposalratio = 0.0;
          priorratio = (oldpriorval-newpriorval)/expo_m_mean;
          /*if (calcoptions[CALCMARGINALLIKELIHOOD])   // calculating marginal likelihood does not seem to work when updating population tree  - not sure why 9/13/2016
          {
            metropolishastingsratio = exp(priorratio + proposalratio);
          }
          else
          {
            metropolishastingsratio = exp(beta[ci] * priorratio + proposalratio);
          }
          U = uniform (); 
          //U = 2;
          if (U <= DMIN(1.0, metropolishastingsratio)) */
          if (calcoptions[CALCMARGINALLIKELIHOOD]) 
          {
            metropolishastingsratio = priorratio + proposalratio;
          }
          else
          {
            metropolishastingsratio = beta[ci] * priorratio + proposalratio;
          }
          if (metropolishastingsdecide(metropolishastingsratio,1))
          {
            accp += 1;
            struct dictionary_node_kr *temp;
            if (dir==0)
              temp = dictionary_install(poppairs[currenti],newpriorval,C[ci]->mltorhpriors); // mltorhpriors are for dir==0
              //C[ci]->mltorhpriors[currenti] = newpriorval;
            else
              temp = dictionary_install(poppairs[currenti],newpriorval,C[ci]->mrtolhpriors); // mrtolhpriors are for dir==1
              //C[ci]->mrtolhpriors[currenti] = newpriorval;
          }
        }
        else  // uniform, accept all updates 
        {
          struct dictionary_node_kr *temp;
          if (dir == 0)
            temp = dictionary_install(poppairs[currenti],getnewpriorval(m_max,oldpriorval,1),C[ci]->mltorhpriors);
            //C[ci]->mltorhpriors[currenti] = getnewpriorval(m_max,oldpriorval,1);
          else
            temp = dictionary_install(poppairs[currenti],getnewpriorval(m_max,oldpriorval,1),C[ci]->mrtolhpriors);
            //C[ci]->mrtolhpriors[currenti] = getnewpriorval(m_max,oldpriorval,1);
          accp += 1; 
        }
        i += 1;
      }
      currenti += 1;
      if (currenti >= numdistinctpopulationpairs[npops])
        currenti = 0;
    }
    else
      currenti += 1;
  }
  *accepted = accp;
  *attempted = i;

}  //update_migration_prior_not_intree


int update_popsize_prior_intree(int ci, int qi)
{
  struct priorvalues holdpr;
  double newmpriorval;
  int accp;
  double metropolishastingsratio;
  double priorratio; //proposalratio is zero

  assert((int) C[ci]->descendantpops[qi] > 0); // no population should have a Null set 
  copy_probcalc (&localholdallpcalc, &C[ci]->allpcalc);
 
  holdpr.max =C[ci]->itheta[qi].pr.max;
  holdpr.min = 0.0;
  newmpriorval = getnewpriorval(q_max,holdpr.max,0);
  C[ci]->itheta[qi].pr.max = newmpriorval;
  priorratio = log(holdpr.max/newmpriorval);  // i.e. log((1/newpriorval)/(1/holdpr.max))
 // proposalratio = 0.0; // with uniform hyperprior 

    // initialize_integrate_tree_prob is overkill because most priors have not changed
  // could write a new function that just does integrations for the terms with new priors 
  initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  priorratio += C[ci]->allpcalc.probg - localholdallpcalc.probg;

  if (calcoptions[CALCMARGINALLIKELIHOOD]) 
  {
    metropolishastingsratio = priorratio;
  }
  else
  {
    metropolishastingsratio = beta[ci] * priorratio;
  }
  if (metropolishastingsdecide(metropolishastingsratio,1))
  {
    accp = 1;
    C[ci]->qhpriors[C[ci]->descendantpops[qi]]= C[ci]->itheta[qi].pr.max;
  }
  else
  {
    C[ci]->itheta[qi].pr.max= holdpr.max;
    copy_probcalc (&C[ci]->allpcalc, &localholdallpcalc);
    accp = 0;
  }
  return accp;

}  //update_popsize_prior_intree


void update_popsize_prior_not_intree(int ci, int *attempted, int *accepted)
{
  double oldpriorval;
  static int numattempts;
  static int currenti = -1;
  int i,j,found,accp;
  static int numnotintree[]={0, 0, 0, 2, 8, 22, 52, 114, 240};

 /* 
with npops,  there are 2*npops-1 populations in the tree
and there are 2^npops - 1  distinct subsets (populations) not counting null set
so in any one tree,  there are 2^npops-1 - (2*npops-1) = 2^npops - 2*npops populations not present
these are values in numnotintreep[]
we can try to update half of these each time we come through
*/
  if (currenti == -1)
  {
    numattempts = numnotintree[npops]/2;
    currenti = 1;  // currenti is cast as a SET and so 0 is the NULL set and is not allowed
  }
  i = 0;
  accp = 0;
  while (i < numattempts)
  {
    found = 0;
    for (j=0;j<numtreepops;j++)
    {
      if (C[ci]->descendantpops[j] == (SET) currenti)
      {
        found = 1;
        break;
      }
    }
    if (found == 0) // attempt update
    {
      oldpriorval = C[ci]->qhpriors[currenti];
      C[ci]->qhpriors[currenti] = getnewpriorval(q_max,oldpriorval,0);
      accp += 1;
      i += 1;
    }
    currenti += 1;
    if (currenti >= numpopsets)
      currenti = 1;
  }
  *accepted = accp;
  *attempted = i;

}  //update_popsize_prior_not_intree



void init_migration_prior_update()
{
  init_probcalc (&localholdallpcalc);
}

void free_migration_prior_update()
{
  free_probcalc(&localholdallpcalc);
}

void init_hyperprior_arrays(int ci)
{
  int i;
  if (ci==ARBCHAIN)
  {
    poppairs = static_cast<char **> (malloc (numdistinctpopulationpairs[npops] * sizeof (char *)));
    for (i = 0; i < numdistinctpopulationpairs[npops]; i++)
    {
      poppairs[i] = static_cast<char *> (malloc ((MAXPOPS_PHYLOGENYESTIMATION +2) * sizeof (char)));
    }
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    //C[ci]->mltorhpriors = static_cast<double *> (malloc (numdistinctpopulationpairs[npops] * sizeof (double)));
    C[ci]->mltorhpriors = static_cast<struct dictionary_node_kr **> (malloc (hashsize* sizeof (struct dictionary_node_kr *)));
    //C[ci]->mrtolhpriors = static_cast<double *> (malloc (numdistinctpopulationpairs[npops] * sizeof (double)));
    C[ci]->mrtolhpriors = static_cast<struct dictionary_node_kr **> (malloc (hashsize* sizeof (struct dictionary_node_kr *)));
    for (i=0;i<hashsize;i++)
    {
      C[ci]->mltorhpriors[i] = 0;
      C[ci]->mrtolhpriors[i] = 0;
    }

    C[ci]->qhpriors = static_cast<double *> (malloc (numpopsets * sizeof (double)));
  }
}

void free_hyperprior_arrays(int ci)
{
  int i;
  if (ci==ARBCHAIN)
  {
    for (i = 0; i < numdistinctpopulationpairs[npops]; i++)
      XFREE(poppairs[i] );
    XFREE(poppairs);
    //XFREE(hashedpairpos);
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    for (i=0;i<hashsize;i++)
    {
      freekrlinkedlist(C[ci]->mrtolhpriors[i]);
      freekrlinkedlist(C[ci]->mltorhpriors[i]);
    }
    XFREE(C[ci]->mltorhpriors);
    XFREE(C[ci]->mrtolhpriors);
    XFREE(C[ci]->qhpriors);
  }
}
