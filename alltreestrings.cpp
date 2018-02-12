/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
/*
build a list of all possible population tree strings

MAXPOPS_PHYLOGENYESTIMATION  should be 7   (as of 5_31_2016) 
works for up to npops=7,  or if there is a ghost population, then npops=8
for larger numbers the list of possbilities just gets too large. 

uses a lot of strings of fixed length to keep dynamic memory stuff down. 
*/

#undef GLOBVARS
#include "ima.hpp"

int numuniquenodes[] = {0,0,0,3,10,25,56,119,246,501};  /* number of unique internal branches,  
     internal branches that occur at least once among all the possible trees (not counting entire tree or leaf nodes , equal to 2^n- (n+2)   
     also this is extern in output.cpp*/
extern int numtreesarray[];  // =  {0,1,1,3,18,180,2700,56700,1587600};  /* number of possible ordered trees,  for up to 7 populations */
extern int numdistinctpopulationpairs[]; // = {0,0,1,6,25,90,301,966,3025}; /* number of distinct pairs of populations that could engage in gene flow (don't share any descendant pops) */

/*********** local to this file  ***************/

#define localmaxpairs  ((MAXPOPS_PHYLOGENYESTIMATION  * (MAXPOPS_PHYLOGENYESTIMATION -1))/2)
typedef char nodearray[MAXPOPS_PHYLOGENYESTIMATION +1][POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];  /* holds the nodes that are being built into a tree string */

/******* local prototypes ********/
int cardinality(SET x);
int makesubsets(SET nset, SET *subsets);
void stringfromset(SET set,char *str);
void findreplace(char *ss,char *oldtext,char *newtext);
void addoutgroup(char s[]);  // also extern in readdata.cpp
void compstrings(char a[], char b[]);
void makenewnode(nodearray nodelist, char *left, char * right, int nextnodenum);
void addnode(char **a,nodearray nodes,int nodeslength,int nextnodenum, int *treecount,int addghost);
int compnodes (const void * a, const void * b);
void findpopsinnode(char *ns,char *popsinnodestr);
void getinternalnodes(char *tstr, char nodes[][20]);


/******* local functions ********/


int cardinality(SET x)  // returns count of items in the set
{
  int count = 0;
	 while (x != EMPTYSET) 
  {
		  x ^= (x & -x); ++count;
		}
	   return(count);
} 
/* fills up subsets with all possible subsets of nset */
int makesubsets(SET nset, SET *subsets)
{
  int i,j;
  SET newset;
  int ssize,last,lslen;
  int memlist[MAXPOPS];

  subsets[0] = EMPTYSET;
  ssize = cardinality(nset);
  j = 0;
  FORALL(i,nset)
  {
    memlist[j] = i;
    j++;
  }
  last = 1;
  for (i=0;i<ssize;i++)
  {
    lslen = last;
    for (j=0;j<lslen;j++)
    {
      newset = SETADD(subsets[j],memlist[i]);
      subsets[lslen+j] = newset;
    }
    last += j;
  }
  return last;

}
/* take a SET and return a string containing the values in order from low to high
   this version assumes no element higher than 9 so there are no separators between
   the elements */ 
void stringfromset(SET set,char *str)
{
  int i;
  char temps[3];
  str[0] = '\0';
  FORALL(i,set)
  {
    sprintf(temps,"%d",i);
    strcat(str,temps);
  }
} //stringfromset



/* find oldtext in ss and replace it with newtext
   works for strings of max fixed length of POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION 
   replaces the old string with the new one */
void findreplace(char *ss,char *oldtext,char *newtext)
{
  char ret[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION],tempss[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  char *holdss;
  int i; //, count = 0; why here
  size_t newlen = strlen(newtext);
  size_t oldlen = strlen(oldtext);

  holdss = &tempss[0];
  strcpy(holdss,ss);
  for (i = 0; holdss[i] != '\0'; i++) 
  {
    if (strstr(&holdss[i], oldtext) == &holdss[i]) 
    {
      //count++;  why was this here? 
      i += (int) (oldlen - 1);
    }
  }
  
  i = 0;
  while (*holdss) 
  {
    if (strstr(holdss, oldtext) == holdss) 
    {
      strcpy(&ret[i], newtext);
      i += (int) newlen;
      holdss += oldlen;
    } 
    else
      ret[i++] = *holdss++;
  }
  ret[i] = '\0';
  strcpy(ss,ret);
}

/* add the outgroup to the newick string,  always goes on the end */ 
void addoutgroup(char s[])
{
  int j,k,outgroupnum;
  int firsta,lasta;
  char temps[3],olds[3],news[3];
  char builds[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];

  j = (int) strlen(s);
  
  while (s[j] != ':')
  {
    j-= 1;
  }
  strcpy(temps,&s[j+1]);
  k = atoi(temps);
  outgroupnum = (k+2)/2;
  lasta = k;
  firsta = outgroupnum;
  k = lasta;
  while (k >= firsta) /* increase all the internal node #'s by 1 */ 
  {
    sprintf(olds,"%d",k);
    sprintf(news,"%d",k+1);
    findreplace(&s[0],olds,news);
    k -= 1;
  }
  sprintf(builds,"(%s,%d):%d",s,outgroupnum,2*outgroupnum);
  strcpy(s,builds);
}

/* compares two node strings and swap if needed */ 
void compstrings(char a[], char b[])
{
  int ia,ib;
  char hold[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  ia = 0;
  //while (a[ia] != '/0' && isdigit(a[ia])==0)
  while (a[ia] != 0 && isdigit(a[ia])==0)
    ia += 1;
  ib = 0;
  //while (b[ib] != '/0' && isdigit(b[ib])==0)
  while (b[ib] != 0 && isdigit(b[ib])==0)
    ib += 1;
  if (atoi(&a[ia])> atoi(&b[ib]))
  {
    strcpy(hold,b);
    strcpy(b,a);
    strcpy(a,hold);
  }
}

/* build a new node from two smaller nodes
removes left and right from nodelist and puts back the new 
joined node into nodelist*/
void makenewnode(nodearray nodelist, char *left, char * right, int nextnodenum)
{
  int i;
  char newnode[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  i=0;
  while (strcmp(nodelist[i],left) != 0)
    i += 1;
  do
  {
    strcpy(nodelist[i],nodelist[i+1]);
    i+=1;
  } while (nodelist[i][0] != '\0');
  i=0;
  while (strcmp(nodelist[i],right) != 0)
    i += 1;
  do
  {
    strcpy(nodelist[i],nodelist[i+1]);
    i+=1;
  } while (nodelist[i][0] != '\0');
  sprintf(newnode,"(%s,%s):%d",left,right,nextnodenum);
  strcpy(nodelist[i-1],newnode);

  /*int temp[MAXTREEPOPS];
  findpopsinnode(left,temp);
  findpopsinnode(right,temp);
  findpopsinnode(newnode,temp); */
}

/* the main recursive part of the function 
makes all possible pairs of nodes  and makes a treestring for each one*/
void addnode(char **a,nodearray nodes,int nodeslength,int nextnodenum, int *treecount,int addghost)
{
  int i,j,k,numpairs;
  nodearray holdnodes;

  char pairs[localmaxpairs][2][POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  k = 0;
  for (i=0;i<nodeslength-1;i++)
    for (j=i+1;j<nodeslength;j++)
    {
      strcpy(pairs[k][0],nodes[i]);
      strcpy(pairs[k][1],nodes[j]);
      compstrings(pairs[k][0],pairs[k][1]);
      k += 1;
    }
  numpairs = k;

  nextnodenum += 1;
  for (i=0;i<numpairs;i++)
  {
    nodeslength = 0;
    do
    {
      strcpy(holdnodes[nodeslength],nodes[nodeslength]);
      nodeslength += 1;
    }
    while (nodes[nodeslength][0] != 0);
    for (;j<=MAXPOPS_PHYLOGENYESTIMATION ;j++)
      holdnodes[j][0] = 0;       //had been sprintf(holdnodes[j],"");
    makenewnode(holdnodes,pairs[i][0],pairs[i][1],nextnodenum);
    nodeslength -= 1;
    if (nodeslength == 1)
    {
      if (addghost)
      {
        addoutgroup(holdnodes[0]);
      }
      assert( a[*treecount]);
      strcpy(a[*treecount],holdnodes[0]);
      *treecount += 1;
    }
    else
    {
      addnode(a,holdnodes,nodeslength,nextnodenum,treecount,addghost);
    }
  }
}

int compnodes (const void * a, const void * b)  // used for qsort in findpopsinnode, sort from low to high 
{
	return ( *(int*)a - *(int*)b );
}

/* take a string containing a node from a tree and return a pointer to a string of sorted (low to hi) array of integers for each 
of the populations in the string separated by commas*/
void findpopsinnode(char *ns,char *popsinnodestr)
{
  char *c;
  int i,j;
  int popsinnode[MAXTREEPOPS];
  char temps[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  int ts;
  i = 0;
  c = ns;
  while (*c != '\0')
  {
    if (isdigit(*c))
    {
      j = atoi(c);  // atoi() works on c and everything after that looks like an integer up to the first char that does not 
      if (j < (npops - modeloptions[ADDGHOSTPOP])) // don't use ancestral pops
      {
        popsinnode[i] = j;
        i+= 1;
      }
      if (j>= 100)
        c += 3;  // move three spaces
      else if (j>= 10)
        c += 2;  // move two spaces
      else
        c += 1;  // move one space 
    }
    else
      c += 1; // move one space

  }
  qsort(popsinnode,i, sizeof(int), compnodes);
  popsinnode[i] = -1;
  ts = 0;
  i = 0;
  do
  {
    ts += sprintf(&temps[ts],"%d,",popsinnode[i]);
    i += 1;
  }while (popsinnode[i] >= 0);
  strcpy(popsinnodestr,temps);
}  //findpopsinnode

// return a list of internal nodes (each represented as an ordered sequence of external node numbers in a string)
void getinternalnodes(char *tstr, char nodes[][20])
{
  int i,ii,j;
  int b,e,cp;
  char tempc[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  char tempcsorted[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  j = 0;
  i = 1; // start at 1 to skip the root node 
  do
  {
    if ( tstr[i]=='(')
    {
      b = i;
      cp = 1;
      ii = i;
      do
      {
        ii += 1;
        if (tstr[ii]=='(')
          cp += 1;
        if (tstr[ii]==')')
          cp -= 1;
      }while (!(tstr[ii] == ')' && cp == 0));
      e = ii; 
      strncpy(&tempc[0],tstr+b,e-b+1);
      tempc[e-b+1] = '\0';
      findpopsinnode(tempc,tempcsorted);
      strcpy(nodes[j],tempcsorted);
      j+= 1;
    }
    i += 1; 
  }while (i < strlen(tstr));
} //getinternalnodes


/******* global functions ********/
/* makes a string out of two sets
   returns 0 if d0 ends up on the left
   returns 1 if d1 ends up on the left */
int makepairstring(SET d0,SET d1, char s[])
{
  char psl[MAXPOPS],psr[MAXPOPS];
  stringfromset(d0,psl);
  stringfromset(d1,psr);
  if (atoi(psl) < atoi(psr))
  {
    sprintf(s,"%s|%s",psl,psr);
    return 0;
  }
  else
  {
    sprintf(s,"%s|%s",psr,psl);
    return 1;
  }
}
/* fill descendantpops[] array for the poptree in chain ci*/
void filldescendantpops(int ci)
{
  int i,j;
  for (i=0;i<numtreepops;i++)
    C[ci]->descendantpops[i] = EMPTYSET;
  for (i=0;i<npops;i++)
  {
    C[ci]->descendantpops[i] = SINGLESET(i);
    j = i;
    while  (C[ci]->poptree[j].down != -1)
    {
      C[ci]->descendantpops[C[ci]->poptree[j].down] = SETADD(C[ci]->descendantpops[C[ci]->poptree[j].down],i);
      j = C[ci]->poptree[j].down; // an ancestral pop # that has i as a descendant 
    }
  }
}  // filldescendantpops


/*
  fillmigratepairs() creates a list of strings, each which contains a pair of populations (including sampled and ancestral populations).
  Each string is a pair of possible populations represented as a string, with a '|' in between two strings of integers
  Each population is represented as a string of the sampled populations that are descendant from it. 
  Each pair string is unique. 
  makes extensive use of SET 

  only called once per cpu,  in alltreestrings.cpp just because it does a lot of string operations and calls other functions in here.  could be in initialize.cpp
*/
int fillmigratepairs(void)
{
  int i,ii,j,k,numsubsets;
  SET *subsets,**subsetsbyk,fullset,seta,comp;
  int countbyk[MAXPOPS] = {0};
  char checkstr[MAXPOPS+1],strseta[MAXPOPS],strcomp[MAXPOPS]; // 1_19_2018 fixed bug added 1 to checkstr length,  is this enough? 
  int numcompsubsets,countpairstrings,found;

  numsubsets = 1 << npops;  // 2^npops
  subsets = static_cast<SET *> (malloc (numsubsets * sizeof (SET)));
  fullset = EMPTYSET;
  for (i=0;i<npops;i++)
    fullset = SETADD(fullset,i);
  int checknum = makesubsets(fullset,subsets);
  assert(checknum==numsubsets);
  subsetsbyk = static_cast<SET **> (malloc ((npops) * sizeof (SET *)));
  for (k=0;k<npops;k++)
    subsetsbyk[k] = static_cast<SET *> (malloc (numsubsets * sizeof (SET )));
  for (i=1;i<numsubsets-1;i++)
  {
    k = cardinality(subsets[i]);
    assert(k<npops);
    subsetsbyk[k][countbyk[k]] = subsets[i];
    countbyk[k] += 1;
  }
  countpairstrings = 0;
  for (k=1;k<npops;k++)
  {
    for (i=0;i<countbyk[k];i++)
    {
      seta = subsetsbyk[k][i];
      stringfromset(seta,strseta);
      comp = SETDIFF(fullset,seta);
      numcompsubsets = makesubsets(comp,subsets);
      for (j=1;j<numcompsubsets;j++)
      {
        stringfromset(subsets[j],strcomp);
        if (atoi(strcomp) < atoi(strseta))
          sprintf(checkstr,"%s|%s",strcomp,strseta);
        else
          sprintf(checkstr,"%s|%s",strseta,strcomp);
        found = 0;
        for (ii=0;ii<countpairstrings;ii++)
          if (strcmp(checkstr,poppairs[ii])==0)
          {
            found = 1;
            break;
          }
        if (found == 0)
        {
          sprintf(poppairs[countpairstrings],"%s",checkstr);
          //hashedpairpos[pairhash(poppairs[countpairstrings])] = countpairstrings;

          countpairstrings++;
        }
      }
    }
  }
  for (k=0;k<npops;k++)
  {
    XFREE(subsetsbyk[k]);
  }
  XFREE(subsetsbyk);
  XFREE(subsets);
  return countpairstrings;
} 



char **allocalltreestrings(void)
{
  char **ats;
  int i;
  int npopsa;
  if (modeloptions[ADDGHOSTPOP]==1)
    npopsa = npops-1;
  else
    npopsa = npops;
  assert (npopsa <= MAXPOPS_PHYLOGENYESTIMATION );
  ats = static_cast<char **> (malloc (numtreesarray[npopsa] * sizeof (char *)));
  for (i = 0; i < numtreesarray[npopsa]; i++)
  {
      ats[i] = static_cast<char *> (malloc (POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION * sizeof (char)));
  }
  return ats;
}

/* alltreestrings is a global char** */ 
void freepoptreestringarrays(void)
{
  int i;
  for (i = 0; i < numpoptopologies; i++)
  {
    //printf("%s\n",alltreestrings[i]);
    XFREE(alltreestrings[i]);
  }
  XFREE(alltreestrings);

  if (modeloptions[ADDGHOSTPOP]==1)
  {
    for (i = 0; i < numpoptopologies; i++)
      XFREE(alltreestrings_noghost[i]);
    XFREE(alltreestrings_noghost);
  }

  XFREE(poptopologycounts);
}

int buildpoptreestringarray(void)
{
  int numtrees, numtreesng;
  nodearray nodes;
  int i,numnodes;
  int npopsa;

  if (modeloptions[ADDGHOSTPOP]==1)
    npopsa = npops-1;
  else
    npopsa = npops;
  for (i=0;i<npopsa;i++)
  {
    sprintf(nodes[i],"%d",i);
  }
  numnodes = i;
  for (;i<=MAXPOPS_PHYLOGENYESTIMATION ;i++)
    nodes[i][0] = 0; // had been sprintf(nodes[i],"");
  numtrees = 0;
  addnode(alltreestrings,nodes,numnodes,npopsa-1, &numtrees,modeloptions[ADDGHOSTPOP]==1);
  for (i=0;i<numtrees;i++)
    rewrite(alltreestrings[i]);
  if (modeloptions[ADDGHOSTPOP]==1)  // redo tree string making, but without the ghost.  should get the same trees but without ghost as outgroup
  {
    numtreesng = 0;
    addnode(alltreestrings_noghost,nodes,numnodes,npops-2, &numtreesng,0);
    for (i=0;i<numtreesng;i++)
      rewrite(alltreestrings_noghost[i]);
    assert(numtrees==numtreesng);
  }
  assert(numtrees == numtreesarray[npops - modeloptions[ADDGHOSTPOP]]);

  
  return numtrees;
} /* buildpoptreestringarray */

/* 
  void printnewickstring(FILE * outto, char *ps, double *tvals, int ghostintree):
    -print a newick string for a poptreestring 
    -tree must come out ultrametric or else a bug 
    -format should correspond to http://evolution.genetics.washington.edu/phylip/newicktree.html
    -ps is the pointer to the poptreestring
    -if *tvals is NULL  just use unit values for each successive time period 
      otherwise use the times in tvals 
    - if ghostintree != 0,  then use "ghostpop" as the name of the ghost population
*/
/* may be a bug in here that does not get the unit branch lengths for deepest nodes in longer trees  
  initialize branchlengths with 1's to deal with this,  but not well. */
void printnewickstring(FILE * outto, char *ps, double *tvals, int ghostintree)
{
  int i,j,k;
  int neededlength = 0;
  double branchlengths[MAXTREEPOPS] = {1.0},sumbl = 0.0;
  double unitlengths[MAXPOPS-1];
  char *c,*cp,*ns;
  int tinc,pcount;
  double t;
  int checklength = 0;
  char tempstr[NAMELENGTH]; // temp holder string, NAMELENGTH should be plenty
  int npopsa;
  int useunitlengths = 0;

  if (modeloptions[ADDGHOSTPOP]==1 && ghostintree == 0 )
    npopsa = npops-1;
  else
    npopsa = npops;
  // estimate needed length of newick string,  should be overkill
  for (i=0;i < (npops-(modeloptions[ADDGHOSTPOP]==1));i++)
    neededlength += (int) strlen(popnames[i]);
  neededlength += (int) strlen(ps) + (numtreepops * 7); // room for branch lengths
  ns = static_cast<char *> (malloc(neededlength*sizeof(char)));  // ns is the pointer to the newick string
  ns[0] = '\0';
  if (tvals==NULL) // set unit length vector for tvals to point to
  {
    for (i=0;i< (npopsa-1);i++)
      unitlengths[i] = 1;
    tvals = unitlengths;
    checklength = 1;
    useunitlengths = 1;
  } 
  tinc = npopsa;  // what to subtract from internal node number to get position in tvals
  /*
  tricky code for getting branch lengths
  for a given node # j move right to find the ancestor k
  get the time t of that ancestor k from its number
  then if the node j is itself is an ancestor
  subtract from t the time of the node j
  */
  c = ps;
  while (*c)
  {
    if (isdigit(*c))
    {
      j = atoi(c);  // atoi() works on c and everything after that looks like an integer up to the first char that does not 
      if (j== 2*(npopsa-1))  // at root node
        goto done;
      if (isdigit(*(c + 1))) // for #'s with 2 digits
        c += 1;
      cp = c;
      pcount = 0;
      while (pcount != 1)
      {
        pcount -= (cp[0] == '(');
        if (cp[0] == ')')
        {
          pcount += 1;
        }
        cp += 1;
      }
      cp +=1; // skip past ':' to node number
      k = atoi(cp) - tinc;  // the position of the t value of the down node
      t = 0.0;
      for (i=0;i<= k;i++)
        t += tvals[i];
      if (j >= npopsa)
      {
        for (i=0;i<= (j-tinc);i++)
          t -= tvals[i];
      }
      branchlengths[j] = t;
      sumbl += t; // for unit t's sumbl = (npopsa (npopsa+1) / 2) -1
    }
    c += 1; // move one space
  }
 done:
  if (checklength)  // if using unit lengths,  then we can calculate the expected length of the tree easily and check it
    assert (INTEGERROUND(sumbl) == ((npopsa * (npopsa+1))/2 - 1) );
  c = ps;
  while (*c)
  {
    if (c[0]==')' || c[0] == '(')
    {
      tempstr[0] = c[0];
      tempstr[1] = '\0';
      strcat(ns,tempstr);
    }
    if (isdigit(*c))
    {
      j = atoi(c); 
      if (isdigit(*(c + 1))) // for #'s with 2 digits
        c += 1;
      if (j<npopsa)
      {
        if (ghostintree && j==npopsa-1)
        {
          sprintf(tempstr,"%s","ghostpop");
        } 
        else
          sprintf(tempstr,"%s",popnames[j]);
        strcat(ns,tempstr);
      }
      tempstr[0] = ':';
      tempstr[1] = '\0';
      strcat(ns,tempstr);
      if (useunitlengths)
        sprintf(tempstr,"%.1f",branchlengths[j]);
      else
        sprintf(tempstr,"%.5f",branchlengths[j]);
      strcat(ns,tempstr);
      if (*(c + 1) == ',')
      {
        tempstr[0] = ',';
        tempstr[1] = '\0';
        strcat(ns,tempstr);
      }
    }
    c += 1;
  }
  tempstr[0] = ';';
  tempstr[1] = '\0';
  strcat(ns,tempstr);
  fprintf(outto,"%s",ns);
  XFREE(ns);
}


/*
  For a population tree strint tstr
  Find each internal node (not root) get the string of sampled populations 
  check to see if that is in the list of nodes (na) that has been found before
  if not add it to the list
  add the count to nc 
*/
//getnodecounts(uniquenodes,fatssarray[i].treestr,fatssarray[i].count,nodecounts, &nunique);
/* na is a list of unique internal nodes  (each represented as an ordered sequence of external node numbes in a string)
   nc is a corresponding set of counts 
   nunique is the number of nodes in na
   tstr is a treestring 
*/
void getnodecounts(char **na,char *tstr,int count,int *nc, int *nunique)
{
  int i,j,ni,found;
  char tree_nodes[MAXPOPS_PHYLOGENYESTIMATION-2][20];
  getinternalnodes(tstr,tree_nodes);
  if (modeloptions[ADDGHOSTPOP]==1)
    ni = npops - 3;
  else
    ni = npops-2;
  for (i=0;i<ni;i++)
  {
    found = 0;
    for (j=0;j<*nunique;j++)
    {
      if (strcmp(na[j],tree_nodes[i])==0)
      {
        found = 1;
        break;
      }
    }
    if (found)
    {
      nc[j] += count;
    }
    else
    {
      nc[j] = count;
      strcpy(na[j],tree_nodes[i]);
      *nunique += 1;
    }
  }
} //  getnodecounts

/*
  calculate the product of the posterior clade probabilities
   na is a list of unique internal nodes  (each represented as an ordered sequence of external node numbes in a string)
   nc is a corresponding set of counts 
   nunique is the number of nodes in na
   tstr is a treestring  */
double calcppcp(char **na,char *tstr, int *nc, int totaltreecount, int nunique)
{
  int i,j,ni,found;
  double ppcp = 1.0;
  char tree_nodes[MAXPOPS_PHYLOGENYESTIMATION-2][20];
  getinternalnodes(tstr,tree_nodes);
  if (modeloptions[ADDGHOSTPOP]==1)
    ni = npops - 3;
  else
    ni = npops-2;
  for (i=0;i<ni;i++)
  {
    found = 0;
    for (j=0;j< nunique;j++)
    {
      if (strcmp(na[j],tree_nodes[i])==0)
      {
        found = 1;
        break;
      }
    }
    if (found)
    {
      if (j < nunique)
        ppcp *= nc[j]/ (double) totaltreecount;
      else
        ppcp = 0.0;
    }
  }
  return ppcp;
} //calcppcp

int foralltreestringsort_comp (const void * a, const void * b)  // used for qsort, sort from high to low // used in  sort_and_print_alltreestrings()
{
  foralltreestringsort *fA = (foralltreestringsort *)a;
  foralltreestringsort *fB = (foralltreestringsort *)b;
  return ( fB->count - fA->count);
}

/* for printing full table of population tree posteriors,  simple */
void printallpoptreesamples (FILE * outfile, int *poptopologycounts,foralltreestringsort *fa, int *poptreeproposed,int uniformprior )
{
  int i,j, totaltreecount = 0;
  double ppcp;
  char poptreetabletitle[100];
  double temp1 = 1.0;

  for (i=0;i<numpoptopologies;i++)
    totaltreecount += poptopologycounts[i];
  /* format table title */
  sprintf(poptreetabletitle,"Full Population Tree Sample Counts and Frequencies (%d entries) - not sorted",numpoptopologies);
  FP "\n%s\n",poptreetabletitle);
  for (i = 0; i< (signed) strlen(poptreetabletitle);i++)
    FP "%c",'-');
  FP "\n");

  if (modeloptions[ADDGHOSTPOP]==0)
  {
    FP "Tree_string\tpriorprob\tproposals\tsampcount\tfrequency\tppcp\n");
    for (i=0;i<numpoptopologies;i++)
    {
      j = 0;
      while (fa[j].count > 0 && fa[j].origi != i)
        j++;
      if (fa[j].count == 0)
        ppcp = 0.0;
      else
        ppcp = fa[j].ppcp;
      if (uniformprior==0)
        temp1 = exp(topologypriors[i]);
      FP "%s\t%.6f\t%d\t%d\t%.6f\t%.7f\n",alltreestrings[i],temp1,poptreeproposed[i],poptopologycounts[i],poptopologycounts[i]/ (double) totaltreecount,ppcp);
    }
    FP "\n");
  }
  else
  {
    FP "Tree\tTree(w/ghost)\tpriorprob\tproposals\tsampcount\tfrequency\tppcp\n");
    for (i=0;i<numpoptopologies;i++)
    {
      j = 0;
      while (fa[j].count > 0 && fa[j].origi != i)
        j++;
      if (fa[j].count == 0)
        ppcp = 0.0;
      else
        ppcp = fa[j].ppcp;
      if (uniformprior==0)
        temp1 = exp(topologypriors[i]);
      FP "%s\t%s\t%.6f\t%d\t%d\t%.6f\t%.7f\n",alltreestrings_noghost[i],alltreestrings[i],temp1,poptreeproposed[i],poptopologycounts[i],poptopologycounts[i]/ (double) totaltreecount,ppcp);
    }
    FP "\n");
  }
}   /*printallpoptreesamples */

/* make RFtreedis array to hold Robinson Foulds distance
  - RF is simply the number of nodes not held in common,  i.e. an integer
  - we are using the distrance from tree 0,  so these values go into the RFtreedis
  - set alltreestrings[0] as the reference and record the list of nodes for that strings 
    -  for all other nodes,  count how many nodes match a node in the list of nodes for alltreestrings[0]
  - store this as the RF distance for that tree (i.e. tree 0 has RF distance of 0)
 */

void init_RF_nodeinfo(void)
{
  char tstr[POPTREESTRINGLENGTHMAX_PHYLOGENYESTIMATION];
  int i,j,k;
  int found,ni,d;
  char tree0nodes[MAXPOPS_PHYLOGENYESTIMATION-2][20];// big enough to hold strings of lists of internal nodes (represetned as strings of external nodes)
  char tree_nodes[MAXPOPS_PHYLOGENYESTIMATION-2][20];
   
  RFtreedis = static_cast<unsigned short *> (malloc (numtreesarray[npops - modeloptions[ADDGHOSTPOP]] * sizeof (unsigned short)));
  if (npops - modeloptions[ADDGHOSTPOP] == 3) // just use treenumber as all distances are 0
  {
    for (k=0;k<numpoptopologies;k++) 
      RFtreedis[k] = k;
    
  }
  else
  {
    if (modeloptions[ADDGHOSTPOP]==1)
    {
      strcpy(tstr,alltreestrings_noghost[0]);
      ni = npops - 3;
    }
    else
    {
      strcpy(tstr,alltreestrings[0]);
      ni = npops - 2;
    }
    getinternalnodes(tstr,tree0nodes);
    RFtreedis[0] = 0;
    for (k=1;k<numpoptopologies;k++) 
    {
      if (modeloptions[ADDGHOSTPOP]==1)
        strcpy(tstr,alltreestrings_noghost[k]);
      else
        strcpy(tstr,alltreestrings[k]);
      getinternalnodes(tstr,tree_nodes);
      d = 0;
      for (i=0;i<ni;i++)
      {
        found = 0;
        for (j=0;j<ni;j++)
        {
          if (strcmp(tree0nodes[i],tree_nodes[j])==0)
            found = 1;
        }
        d += (found==0);
      }
      RFtreedis[k] = d;
      //printf("%d ",d);
    }
  }
}   //init_RF_nodeinfo(void)