/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

extern void set_x (struct value_record *v, int isint); // in initialize.cpp
extern void init_value_record (struct value_record *v, int isint); // in initialize.cpp

/*********** LOCAL STUFF but declared extern in other files**********/
void fillplist (int ci);
void fillancplist (int ci);
void poptreewrite (int ci, char *buildstr);
void orderupnodes(struct popedge *poptree);

/*********** LOCAL STUFF **********/
static int numpopnodes;
static char *treestringspot;
static int pos;
static int ancestralpopnums[2 * MAXPOPS - 1];
static void parenth0 ();
static void parenth (int ci, int tempcurrent, int startparenth);
static void simpoptree (int ci);
static void poptreeread (int ci, char *poptreestring);
static void makepoptreestring (int ci, int curpop, char *buildstr);
//static void set_x_local (struct value_record *v, int isint); don't need 
//static void init_value_record_local (struct value_record *v, int isint);  don't need, can extern to the other one

/* popstring primary format:
for npops,  the populations are numbered 0 thru npops-1
format includes branching pattern and node order in time
every node is flanked by parentheses, and every pair of items within a node is separated by a comma
every closed parentheses is followed by a colon and a node sequence number
node sequence numbers are ordered from most recent (npops) oldest,  2*(npops - 1)

In other words, the external nodes (sampled populations) are numbered
from 0 to npops - 1,  and the internal nodes are numbered from npops to 2*(npops - 1)
*/

/* parenth0() identifies the order of nodes as they will be reached by parenth()
then it associates these with the correct period and ancestral population size numbers
and puts these into ancestralpopnums[] to be used by parenth()  

How it works:
come in on first opening parenthesis
the ith opening parenthesis is associated with a number after its corresponding closing parenthesis

parenth0() will go through parentheses from left to right,
when it comes to a close parenthesis it records the ancestral node number for that pairing

ancestral node numbers range from numpops to 2*numpops-2 and proceed from lowest
to highest in order of what time they occured

the ancestral node number is recorded in the array ancestralpopnums[]

the position in the array is the count of which parenthesis pair has just closed, plus numpops

in other words if it is the 0'th parenthesis pair (i.e. the first one that opened, meaning it is the
outermost pair),  then the ancestral node number is recorded in ancestralpopnums[numpops]

If it is the ith pair that has closed, it is recorded in ancestralpopnums[numpops + i]

This seemed to be the only way to get the correct labeling of internal nodes. 

Then when the function parenth() is called,  the correct times and populatino sizes can be associated 
with these ancestral populations. 
*/

/**** LOCAL FUNCTIONS *****/
/* don't need these copies of set_x and init_value_record
void
set_x_local (struct value_record *v, int isint) // same as in initialize,  but needed it in this file as well, so made a local copy
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
}                               // set_x_local 

void
init_value_record_local (struct value_record *v, int isint)   // same as in initialize,  but needed it in this file as well, so made a local copy
{
  if (v->do_xyplot)
  {
    v->xy = static_cast<plotpoint *> 
                (calloc (GRIDSIZE, sizeof (struct plotpoint)));
    set_x_local (v, isint);
  }

  if (v->do_trend)
    v->trend = static_cast<double *> (calloc (TRENDDIM, sizeof (double)));

  v->beforemin = v->aftermax = 0;
}                               // init_value_Record_local */

void
parenth0 ()
{
  int itemp;
  char *ne;
  int psetlist[MAXPOPS], nextlistspot, popennum;
  nextlistspot = 0;
  popennum = 0;
  ne = treestringspot;
  while (*ne != '\0')
  {
    if (*ne == '(')
    {
      psetlist[nextlistspot] = popennum;
      nextlistspot++;
      popennum++;
      ne++;
    }
    else
    {
      if (*ne == ')')
      {
        ne += 2;
        itemp = strtol (ne, &ne, 10);
        ancestralpopnums[npops + psetlist[nextlistspot - 1]] = itemp;
        nextlistspot--;
      }
      else
      {
        ne++;
      }
    }
  }
}                               /* parenth0 */

/* parenth()  reads the tree topology and puts it into poptree structure */

void
parenth (int ci, int tempcurrent, int startparenth)
/* recursive - reenters for every open parentheses */
{
  int itemp, current, i;
  static int nextnode;
  static int periodi;

  if (startparenth == 1)
  {
    nextnode = -1;
    periodi = 0;
  }
  current = ancestralpopnums[tempcurrent];

  treestringspot++;
  while (isspace (*(treestringspot + 1))) // why  + 1 ? 
    treestringspot++;

  /* next could be:
     - a number  (a simple upnode) read it and get ':' and float number after it
     - an open parenthesis (a complex upnode - call parenth)
     - a comma  skip it
     - a close parenthesis - close the node 
   */
  do
  {
    if (isdigit (*treestringspot))
    {
      itemp = atoi (treestringspot);
      treestringspot++;
      C[ci]->poptree[current].up[C[ci]->poptree[current].numup] = itemp;
      C[ci]->poptree[current].numup++;
      C[ci]->poptree[current].time = 0;
      C[ci]->poptree[itemp].down = current;
    }
    if (*treestringspot == ',')
      treestringspot++;
    if (*treestringspot == '(')
    {
      if (nextnode == -1)
        nextnode = npops + 1;
      else
        nextnode++;
      C[ci]->poptree[ancestralpopnums[nextnode]].down = current;
      C[ci]->poptree[current].up[C[ci]->poptree[current].numup] =
        ancestralpopnums[nextnode];
      C[ci]->poptree[current].numup++;
      C[ci]->poptree[current].time = 0;
      parenth (ci, nextnode, 0);
    }
  } while (*treestringspot != ')');
  treestringspot++;             /* skip parentheses */
  if (*treestringspot == ':')
  {
    treestringspot++;
    i = atoi (treestringspot);
    if (i < npops)
      IM_err (IMERR_POPTREESTRINGFAIL,
              " wrong number of ancestral populations indicated. string %s, step %d",
              treestringspot, step);
    assert (i >= npops);
    periodi = i - npops;
    C[ci]->poptree[current].b = periodi + 1;
    C[ci]->poptree[C[ci]->poptree[current].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[current].up[1]].e = periodi + 1;
    if (i >= 10)
      treestringspot += 2;
    else
      treestringspot++;
  }
  else
  { // is it possible to get here? 
    C[ci]->poptree[current].b = periodi + 1;
    C[ci]->poptree[C[ci]->poptree[current].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[current].up[1]].e = periodi + 1;
    periodi++;
  }
  if (C[ci]->poptree[current].down != -1)
  {
    numpopnodes++;
    current = C[ci]->poptree[current].down;
  }
  else
  {
    periodi++;
    C[ci]->poptree[current].e = -1;
  }
}                               /* parenth */

/* find root,  move up, recurse, build the string, makes strings 
 * that contain node order information 
 */
void
makepoptreestring (int ci, int curpop, char *buildstr)
{
  int i;
  char ss[6];
  if (curpop == -1)
  {
    i = 0;
    while (C[ci]->poptree[i].down != -1)
      i++;
    curpop = i;
    pos = 0;
  }
  sprintf (ss, "(");
  strinsert (buildstr, ss, pos);
  pos += (int) strlen(ss);
  for (i = 0; i < C[ci]->poptree[curpop].numup; i++)
  {
    if (C[ci]->poptree[curpop].up[i] < npops)
    {
      sprintf (ss, "%d", C[ci]->poptree[curpop].up[i]);
      strinsert (buildstr, ss, pos);
      pos += (int) strlen(ss);
    }
    else
    {
      makepoptreestring (ci, C[ci]->poptree[curpop].up[i],
                         buildstr /*, pos */ );
    }
    if (i < C[ci]->poptree[curpop].numup - 1)
    {
      sprintf (ss, ",");
      strinsert (buildstr, ss, pos);
      pos += (int) strlen(ss);
    }
    else
    {
      sprintf (ss, "):%d", curpop);
      strinsert (buildstr, ss, pos);
      pos += (int) strlen(ss);
    }
  }
}                               /* makepoptreestring */

/* rewrite() rewrites the treestring in a standard order
    swivels nodes,  if both have node sequence values, the one with the lower node sequence value (periodi[]) goes on the left
    if only one has a node sequence value,  it goes on the right
	when neither has a node sequence value, the one with the lowest node number go on the left */
/* this is simply sorting for a pair.  To handle multifurcations, must put in proper sorting */
/* it actually works,  checked on simulated random trees on 4/17/07 */
/*rewrite if for trees that have node sequence values, it is based on an older version rewrite for trees without node sequence values */
/* works recursively */

void rewritecheckchar(char c)
{
  if ((isdigit(c) || c=='(' || c==',' || c==')' || c==':') == 0)
    IM_err (IMERR_POPTREESTRINGFAIL," something wrong in formatting of population string in input file");
}

void rewrite (char *substr)
{
  int slengths[MAXPOPS];
  int pcount, subpos, subcount;
  char holdsubs[MAXPOPS][POPTREESTRINGLENGTHMAX];
  int firstint[MAXPOPS];
  int i, j, k;
  int periodi[MAXPOPS];
  pos = 1;
  subpos = pos;
  subcount = 0;
  pcount = 0;
  slengths[subcount] = 0;

  do
  {
    if (substr[pos] == '(')
      pcount++;
    if (substr[pos] == ')')
      pcount--;
    pos++;
    slengths[subcount]++;
    if (pcount == 0)
    {
      if (slengths[subcount] > 1)
      {
        pos++;
        i = atoi (&substr[pos]);
        periodi[subcount] = i;
        if (i >= 10)
        {
          pos += 2;
          slengths[subcount] += 3;
        }
        else
        {
          pos++;
          slengths[subcount] += 2;
        }
      }
      else
      {
        periodi[subcount] = -1;
      }
      assert (!(pos < subpos));
      strncpy (holdsubs[subcount], &substr[subpos], (size_t) (pos - subpos));
      holdsubs[subcount][slengths[subcount]] = '\0';
      i = 0;
      while (!isdigit (holdsubs[subcount][i]))
        i++;
      firstint[subcount] = atoi (&holdsubs[subcount][i]);
      subcount++;
      slengths[subcount] = 0;
      if (substr[pos] == ',')
      {
        pos++;
      }
      subpos = pos;
    }
  } while (pos < (int) strlen (substr));
  if ((periodi[0] > periodi[1] && periodi[0] >= 0 && periodi[1] >= 0)
      || (periodi[0] >= 0 && periodi[1] < 0))
  {
    substr[0] = '(';
    j = slengths[1];
    for (i = 1, k = 0; i <= j; i++, k++)
    { 
      rewritecheckchar(holdsubs[1][k]);
      substr[i] = holdsubs[1][k];
    }
    subpos = 1;
    substr[i] = '\0';
    if (slengths[1] > 2)
      rewrite (&substr[subpos]);
    substr[i] = ',';
    i++;
    subpos = i;
    j += 1 + slengths[0];
    for (k = 0; i <= j; i++, k++)
    {
      rewritecheckchar(holdsubs[0][k]);
      substr[i] = holdsubs[0][k];
    }
    substr[i] = '\0';
    if (slengths[0] > 2)
      rewrite (&substr[subpos]);
    substr[i] = ')';
  }
  else
  {
    if (firstint[0] > firstint[1] && periodi[0] < 0 && periodi[1] < 0)
    {
      substr[0] = '(';
      j = slengths[1];
      for (i = 1, k = 0; i <= j; i++, k++)
      {
        rewritecheckchar(holdsubs[1][k]);
        substr[i] = holdsubs[1][k];
      }
      subpos = 1;
      if (slengths[1] > 2)
        rewrite (&substr[subpos]);
      substr[i] = ',';
      i++;
      subpos = i;
      j += 1 + slengths[0];
      for (k = 0; i <= j; i++, k++)
      {
        rewritecheckchar(holdsubs[0][k]);
        substr[i] = holdsubs[0][k];
      }
      if (slengths[0] > 2)
        rewrite (&substr[subpos]);
      substr[i] = ')';
    }
    else
    {
      substr[0] = '(';
      subpos = 1;
      substr[slengths[0] + 1] = '\0';
      if (slengths[0] > 2)
        rewrite (&substr[subpos]);
      substr[slengths[0] + 1] = ',';
      subpos = slengths[0] + 2;
      substr[slengths[0] + slengths[1] + 2] = '\0';
      if (slengths[1] > 2)
        rewrite (&substr[subpos]);
      substr[slengths[0] + slengths[1] + 2] = ')';
    }
  }
} /* rewrite */ ;

/* create a 2d list ancplist[i][j] = k   means that sampled population i, during period j,   is in ancestral population k 
  create this after plist is created
  can always find the next down time for any sampled population  by first finding its current ancestral population
    and then looking up the down time for that population in poptree
 e.g. for population i in period j    dntime = C[ci]->poptree[C[ci]->ancplist[i][j]].time
 to get the ancestral population that a given sampled population is in at some time t,  use findperiod(ci,t)  [look in update_gree_common.cpp]
 e.g.  for population i at time t    
  per = findperiod (ci, t)
  apop = C[ci]->ancplist[i][per]
*/
void
fillancplist (int ci)
{
  int i, j;
  int numperiods;
  numperiods = npops;


  for (i = 0; i < npops; i++)
  {
    C[ci]->ancplist[i][0] = i;
  }
  for (i = 0; i < npops; i++)
  {
    for (j=1;j<npops; j++)
      if (C[ci]->ancplist[i][j-1] == C[ci]->droppops[j][0] || C[ci]->ancplist[i][j-1] == C[ci]->droppops[j][1])
      {
        C[ci]->ancplist[i][j] = C[ci]->addpop[j];
      }
      else
        C[ci]->ancplist[i][j] = C[ci]->ancplist[i][j-1];
  }
  for (i = npops; i < numtreepops; i++)
  {
    for (j=0;j<=i-npops;j++)
      C[ci]->ancplist[i][j] = -1;
    C[ci]->ancplist[i][j] = i;
  }
  for (i = npops; i < numtreepops; i++)
  {
    for (j=(i-npops)+2;j<numperiods; j++)
      if (C[ci]->ancplist[i][j-1] == C[ci]->droppops[j][0] || C[ci]->ancplist[i][j-1] == C[ci]->droppops[j][1])
      {
        C[ci]->ancplist[i][j] = C[ci]->addpop[j];
      }
      else
        C[ci]->ancplist[i][j] = C[ci]->ancplist[i][j-1];
  }
}                               /* fillancplist */


/*
 fillplist: fills 
	C[ci]->periodset[i]  the set of populations that exist in period i
	C[ci]->addpop[i] the population that appears at the beginning of period i
	C[ci]->droppops[i][2] the two populations that join at the beginning 
		of period i (i.e. that are dropped  at the begining of period i)
	C[ci]->plist[i][j] is a 2d array, one row for each period. In row i 
		the elements are the population numbers that exist during that period,  in 
		order from low to high.  These are the same population numbers that are in 
		periodset[i], with the elements in order from lowest to highest. 
*/
void
fillplist (int ci)
{
  int i, j, k;
  SET tempset;

  C[ci]->periodset[0] = EMPTYSET;
  for (i = 0; i < npops; i++)
    C[ci]->periodset[0] = UNION (C[ci]->periodset[0], SINGLESET (i));
  tempset = C[ci]->periodset[0];
  C[ci]->addpop[0] = C[ci]->droppops[0][0] = C[ci]->droppops[0][1] = -1;
  C[ci]->addpop[npops] = C[ci]->droppops[npops][0] = C[ci]->droppops[npops][1] = 0;
  for (i = 1; i < npops; i++)   // loop over periods
  {
    k = 0;
    for (j = 0; j < numtreepops; j++)
    {
      if (C[ci]->poptree[j].e == i)
      {
        if (k >= 2)
        {
//#ifdef DEUBG
      //    poptreeprint(ci);
//#endif
          IM_err (IMERR_POPTREESTRINGFAIL,  " wrong number ancestral populations indicated. step %d treestring %s",
                  step,C[ci]->chainpoptreestring);
        }

        assert (ISELEMENT (j, tempset));
        tempset = SETREMOVE (tempset, j);       // remove j from the set
        C[ci]->droppops[i][k] = j;
        k++;
      }
      if (C[ci]->poptree[j].b == i)
      {
        tempset = SETADD (tempset, j);  // add j to the set
        C[ci]->addpop[i] = j;
      }
    }
    C[ci]->periodset[i] = tempset;
  }
  for (i = 0; i < npops; i++)
    C[ci]->plist[0][i] = i;
  for (i = 0; i < npops; i++)
  {
    j = 0;
    FORALL (k, C[ci]->periodset[i])
    {
      C[ci]->plist[i][j] = k;
      j++;
    }
  }
}                               /* fillplist */
/* 1/4/2017  jh added this.  */
void orderupnodes(struct popedge *poptree)
{
  int i,temp;
  //int numtreepops = 2*npops - 1; is global
  for (i=npops;i<numtreepops;i++)
  {
    assert (poptree[i].up[0] >= 0 && poptree[i].up[0] < numtreepops);
    assert (poptree[i].up[1] >= 0 && poptree[i].up[1] < numtreepops);
    assert (poptree[i].up[0] != poptree[i].up[1]);
    if (poptree[i].up[0] > poptree[i].up[1]) 
    {
      temp = poptree[i].up[0];
      poptree[i].up[0] = poptree[i].up[1];
      poptree[i].up[1] = temp;
    }
  }
}

/*
 * This function is used for debugging only.  CR 110825.1  
 */
void
simpoptree (int ci)
{
  int i, k1, k2, newpop, n, periodi;
  int list[2 * MAXPOPS - 1];
  for (i = 0; i < numtreepops; i++)
    C[ci]->poptree[i].up[0] = C[ci]->poptree[i].up[1] =
      C[ci]->poptree[i].down = -1;
  for (i = 0; i < npops; i++)
  {
    C[ci]->poptree[i].b = 0;
    list[i] = i;
  }
  for (periodi = 1, newpop = npops, n = npops; newpop < numtreepops;
       newpop++, n--, periodi++)
  {
    C[ci]->poptree[newpop].numup = 2;
    k1 = randposint (n);

    do
    {
      k2 = randposint (n);
    } while (k2 != k1);
    C[ci]->poptree[newpop].up[0] = list[k1];
    C[ci]->poptree[newpop].up[1] = list[k2];
    C[ci]->poptree[newpop].b = periodi;
    for (i = k1; i < n; i++)
      list[i] = list[i + 1];
    for (i = k2; i < n - 1; i++)
      list[i] = list[i + 1];
    list[n - 2] = newpop;
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].time =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].time =
      C[ci]->tvals[periodi];
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].e =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].e = periodi;
    C[ci]->poptree[C[ci]->poptree[newpop].up[0]].down =
      C[ci]->poptree[C[ci]->poptree[newpop].up[1]].down = newpop;
  }
  C[ci]->poptree[2 * npops - 2].down = -1;
  C[ci]->poptree[2 * npops - 2].time = TIMEMAX;
}                               /* simpoptree */

void
poptreeread (int ci, char *poptreestring)
{
  int i, j;

  /* read in the tree string until enough parentheses are found */
  /* pcount counts parentheses '(' is +1 ')' is -1 repeat until 0 */
  numpopnodes = 0;
  for (i = 0; i < npops; i++)
  {
    C[ci]->poptree[i].b = 0;
    C[ci]->poptree[i].numup = 0;
    /*C[ci]->poptree[i].up = static_cast<int *> (malloc (2 * sizeof (int)));  */  /*changed up to a fixed array for hidden genealogy stuff */
    for (j = 0; j < 2; j++)
      C[ci]->poptree[i].up[j] = -1;
    C[ci]->poptree[i].down = -1;
  }
  for (; i < numtreepops; i++)
  {
    C[ci]->poptree[i].numup = 0;
    /*C[ci]->poptree[i].up = static_cast<int *> (malloc (2 * sizeof (int))); */  /*changed up to a fixed array for hidden genealogy stuff */
    for (j = 0; j < 2; j++)
      C[ci]->poptree[i].up[j] = -1;
    C[ci]->poptree[i].down = -1;
  }
  rewrite (poptreestring);
  C[ci]->poptree[npops].down = -1;
  treestringspot = poptreestring;
  if (ci == 0|| modeloptions[POPTREETOPOLOGYUPDATE] == 1)
    parenth0 ();
  parenth (ci, npops, 1);
  rewrite(poptreestring);
  return;
}                               /* end treeread */

/* write the treestring for chain ci  to buildstr */
void
poptreewrite (int ci, char *buildstr)
{
#ifdef TURNONCHECKS
//  poptreeprint(ci);
  checkpoptree(ci,1);
#endif
  buildstr[0] = 0;
  makepoptreestring (ci, -1, buildstr);
  rewrite (buildstr);
}


/********** GLOBAL FUNCTIONS ***********/

int getpoptreestringnum(char *s)
{
  int i=0;
  int dcheck,gcheck;
  if (modeloptions[ADDGHOSTPOP])  // check to see if ghost population number is in position 2, if not return -1
  {
    dcheck = isdigit(*(s+1));
    gcheck = atoi(s+1);
    if (!(dcheck && gcheck==(npops-1)))
      return -1;
  }
  while (i < numpoptopologies && strcmp(s,alltreestrings[i]) != 0)
    i += 1;
  if (i>= numpoptopologies)
  {
    IM_err (IMERR_POPTREESTRINGFAIL, " poptreestring not found in list: %s",s);
  }
  return i;
}

void
add_ghost_to_popstring (char poptreestring[])
{
  size_t i;
  int n;
  char stringBuf[POPTREESTRINGLENGTHMAX]; /* temp buffer to build newstring */
  char newstring[POPTREESTRINGLENGTHMAX];

  strcpy (newstring, "(");
  for (i=0;i<strlen(poptreestring);i++)
  {
    if (poptreestring[i] == ':')
    {
      n = atoi(&poptreestring[i+1]);
      sprintf(stringBuf,"%s%c%d",newstring,':',n+1);
      i++;
      if (n>=10)
        i++;
      strcpy(newstring, stringBuf);
    }
    else
    {
      strncat(newstring,&poptreestring[i],1);
    }
  }
  sprintf (poptreestring, "%s,%d):%d", newstring, npops, 2*npops);
}

void
setup_poptree (int ci, char poptreestring[])
{
  if (npops == 1 )
  {
    C[ci]->poptree = static_cast<popedge *> (malloc (sizeof (struct popedge)));
    C[ci]->poptree[0].numup = 2;
    /*C[ci]->poptree[0].up = static_cast<int *> (malloc (2 * sizeof (int))); */  /*changed up to a fixed array for hidden genealogy stuff */
    C[ci]->poptree[0].down = -1;
    C[ci]->poptree[0].time = TIMEMAX;
    C[ci]->poptree[0].b = 0;
    C[ci]->poptree[0].e = -1;
    C[ci]->poptree[0].up[0] = -1;
    C[ci]->poptree[0].up[1] = -1;
  }
  else   
  {
    checktreestring (poptreestring);
    strcpy (C[ci]->chainpoptreestring, poptreestring);

    
    C[ci]->poptree = static_cast<popedge *> 
            (malloc (numtreepops * sizeof (struct popedge)));
    if (strlen (C[ci]->chainpoptreestring) > 0)
    {
      poptreeread (ci, C[ci]->chainpoptreestring);
    }
    else
    {
      // this was used for debugging poptreestring code 
#ifdef TURNONCHECKS
      assert(0); // need to see if we get here
      simpoptree (ci);
      poptreewrite (ci, C[ci]->chainpoptreestring);
#endif //TURNONCHECKS

    }
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      C[ci]->poptreenum = getpoptreestringnum(C[ci]->chainpoptreestring);

    }
  }
  orderupnodes(C[ci]->poptree);
  return;
}                               /* setup_poptree */

/* when reading in a tree from mcf file, the tree has already been initialized, but it must now be reset
with the poptreestring that was read in */
void
reset_poptree (int ci, char poptreestring[])
{

  if (npops > 1 )
  {
    checktreestring (poptreestring);
    strcpy (C[ci]->chainpoptreestring, poptreestring);
    poptreeread (ci, C[ci]->chainpoptreestring);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      C[ci]->poptreenum = getpoptreestringnum(C[ci]->chainpoptreestring);

    }
  }
  return;
}                               /* reset_poptree */


void set_poptree_update_record(void)
{
  int i;
  poptreeuinfo = static_cast<chainstate_record_updates_and_values *>  (malloc (sizeof (struct chainstate_record_updates_and_values)));
  sprintf (poptreeuinfo->str, "PopTree");
  poptreeuinfo->num_uptypes = IM_UPDATE_POPTREE_NUMBER;
  poptreeuinfo->upnames = static_cast<strnl *>  (malloc (poptreeuinfo->num_uptypes * sizeof (strnl)));
  sprintf (poptreeuinfo->upnames[IM_UPDATE_POPTREE_ANY],      "branch     ");
  sprintf (poptreeuinfo->upnames[IM_UPDATE_POPTREE_TOPOLOGY], "topology   ");
  sprintf (poptreeuinfo->upnames[IM_UPDATE_POPTREE_TMRCA],    "tmrca      ");
  poptreeuinfo->upinf = static_cast<update_rate_calc *> 
        (calloc ((size_t) poptreeuinfo->num_uptypes,sizeof (struct update_rate_calc)));
  for (i=0;i<poptreeuinfo->num_uptypes;i++)
  {
    poptreeuinfo->upinf[i].accp = 0;
    poptreeuinfo->upinf[i].tries = 0;
  }
  //poptreeuinfo->num_vals = 0;
  poptreeuinfo->num_vals = 1;
  poptreeuinfo->v = static_cast<value_record *> 
            (malloc (poptreeuinfo->num_vals * sizeof (struct value_record)));
  if (npops - modeloptions[ADDGHOSTPOP] < MINPOPSFORDISTANCE)
  {
    sprintf (poptreeuinfo->v->str, "topology");
    sprintf (poptreeuinfo->v->strshort, "Topol");
  }
  else
  {
    sprintf (poptreeuinfo->v->str, "topology distance");
    sprintf (poptreeuinfo->v->strshort, "TopolD");
  }
  poptreeuinfo->v->plotrange.min = 0;
  poptreeuinfo->v->plotrange.max = numpoptopologies-1;
  poptreeuinfo->v->do_autoc = 1;
  poptreeuinfo->v->do_xyplot = 1;
  poptreeuinfo->v->do_trend = 1;
  poptreeuinfo->v->plotrescale = 1.0;
  poptreeuinfo->v->do_logplot = 0;
  //init_value_record_local (poptreeuinfo->v, 0); 
  init_value_record (poptreeuinfo->v, 0); 
  return;
}