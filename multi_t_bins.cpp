/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

#undef GLOBVARS
#include "ima.hpp"

/****** LOCAL STUFF **********/

#define MAX_FOR_MULTI_T_ARRAY  4

static int **t2array;
static int ***t3array;
static int ****t4array;
static double tarraybinvals[MAX_FOR_MULTI_T_ARRAY][NUMTARRAYBINS];
static int tvalarraypos[MAX_FOR_MULTI_T_ARRAY] = { 0 };
static int highpos[MAX_FOR_MULTI_T_ARRAY];
static int numt;
static int init = 0;

/* set up multidimensional arrays for t values */

/******* GLOBAL FUNCTIONS ***********/

void
setup_multi_t_arrays (int z)
{

  int i, j,k;
  static double tpriormax;
  if (init == 0)
  {
    numt = npops - 1;
    tpriormax = T[numt-1].pr.max;
    for (i = 0; i < numt; i++)
      for (j = 0; j < NUMTARRAYBINS; j++)
        tarraybinvals[i][j] = T[i].pr.min + ((j + 0.5) * (T[i].pr.max  - T[i].pr.min)) / NUMTARRAYBINS;

    if (numt == 2)
    {
      if ((t2array = static_cast<int **> 
                    (malloc (NUMTARRAYBINS * sizeof (*t2array)))) == NULL)
        IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
      for (i = 0; i < NUMTARRAYBINS; i++)
        if ((t2array[i] = static_cast<int *> (calloc ((size_t) (NUMTARRAYBINS - i), sizeof (*t2array[i])))) == NULL)
          IM_err (IMERR_MEM, "  t2array problem in setup_multi_t_arrays()");
    }

    if (numt == 3)
    {
      if ((t3array = static_cast<int ***> 
                (malloc (NUMTARRAYBINS * sizeof (*t3array)))) == NULL)
        IM_err (IMERR_MEM, "  t3array problem in setup_multi_t_arrays()");
      for (i = 0; i < NUMTARRAYBINS; i++)
      {
        if ((t3array[i] = static_cast<int **> (malloc ((NUMTARRAYBINS - i) * sizeof (*t3array[i])))) == NULL)
          IM_err (IMERR_MEM, "  t3array problem in setup_multi_t_arrays()");
        for (j = 0; j < (NUMTARRAYBINS - i); j++)
          if ((t3array[i][j] = static_cast<int *> (calloc ((size_t) (NUMTARRAYBINS - j), sizeof (*t3array[i][j])))) == NULL)
            IM_err (IMERR_MEM, "  t3array problem in setup_multi_t_arrays()");
      }
    }
    if (numt == 4)
    {
      if ((t4array = static_cast<int ****> 
                    (malloc (NUMTARRAYBINS * sizeof (*t4array)))) == NULL)
        IM_err (IMERR_MEM, "  t4array problem in setup_multi_t_arrays()");
      for (k = 0; k < NUMTARRAYBINS; k++)
      {
        if ((t4array[k] = static_cast<int ***> 
                (malloc ((NUMTARRAYBINS-k) * sizeof (*t4array[k])))) == NULL)
          IM_err (IMERR_MEM, "  t4array problem in setup_multi_t_arrays()");
        for (i = 0; i < (NUMTARRAYBINS - k); i++)
        {
          if ((t4array[k][i] = static_cast<int **> (malloc ((NUMTARRAYBINS - i) * sizeof (*t4array[k][i])))) == NULL)
            IM_err (IMERR_MEM, "  t4array problem in setup_multi_t_arrays()");
          for (j = 0; j < (NUMTARRAYBINS - i); j++)
            if ((t4array[k][i][j] = static_cast<int *> (calloc ((size_t) (NUMTARRAYBINS - j), sizeof (*t4array[k][i][j])))) == NULL)
              IM_err (IMERR_MEM, "  t4array problem in setup_multi_t_arrays()");
        }
      }
    }
    for (i = 0; i < numt; i++)
      highpos[i] = 0;
  }
  init++;

	if (z >= 0) {
  for (i = 0; i < numt; i++)
  {
    if (C[z]->tvals[i] > T[i].pr.max)
      return;
    if (i == 0)
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((C[z]->tvals[i] - T[i].pr.min) / (tpriormax  - T[i].pr.min)));
    else
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((C[z]->tvals[i] - T[i].pr.min) / (tpriormax  - T[i].pr.min))) - tvalarraypos[i - 1];
  }
	}
  if (numt == 2)
  {
    t2array[tvalarraypos[0]][tvalarraypos[1]]++;
    if (t2array[tvalarraypos[0]][tvalarraypos[1]] >
        t2array[highpos[0]][highpos[1]])
    {
      highpos[0] = tvalarraypos[0];
      highpos[1] = tvalarraypos[1];
    }
  }
  if (numt == 3)
  {
    t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]]++;
    if (t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]] >
        t3array[highpos[0]][highpos[1]][highpos[2]])
    {
      highpos[0] = tvalarraypos[0];
      highpos[1] = tvalarraypos[1];
      highpos[2] = tvalarraypos[2];
    }
  }
  if (numt == 4)
  {
    t4array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]][tvalarraypos[3]]++;
    if (t4array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]][tvalarraypos[3]] >
        t4array[highpos[0]][highpos[1]][highpos[2]][highpos[3]])
    {
      highpos[0] = tvalarraypos[0];
      highpos[1] = tvalarraypos[1];
      highpos[2] = tvalarraypos[2];
      highpos[3] = tvalarraypos[3];
    }
  }
}                               //setup_multi_t_arrays()

void
free_multi_t_arrays ()
{
  int i, j,k;

  if (numt == 2)
  {
    orig2d_free2D ((void **) t2array, NUMTARRAYBINS);
  }

  if (numt == 3)
  {
    for (i = 0; i < NUMTARRAYBINS; i++)
    {
      for (j = 0; j < (NUMTARRAYBINS - i); j++)
        XFREE (t3array[i][j]);
      XFREE (t3array[i]);
    }
    XFREE (t3array);
  }
  if (numt == 4)
  {
    for (k = 0; k < NUMTARRAYBINS; k++)
    {
      for (i = 0; i < NUMTARRAYBINS - k; i++)
      {
        for (j = 0; j < (NUMTARRAYBINS - i); j++)
          XFREE (t4array[k][i][j]);
        XFREE (t4array[k][i]);
      }
      XFREE (t4array[k]);
    }
    XFREE (t4array);
  }
  if (numt < 2 || numt > 4)
  {
    // should not get here
    return;
  }
}                               //free_multi_t_arrays()

void
return_joint_t (double tvals[])
{
  int i;
  for (i = 0; i < numt; i++)
  {
    if (i == 0)
      tvals[i] = tarraybinvals[i][highpos[i]];
    else
      tvals[i] = tarraybinvals[i][highpos[i] + highpos[i - 1]];
  }
}

// estimate joint posterior probability of a t value 
//  not sure what this would be used for 
double
joint_t_prob (double *tvals)
{
  int i;
  for (i = 0; i < numt; i++)
  {
    if (i == 0)
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((tvals[i] - T[i].pr.min) / (T[i].pr.max - T[i].pr.min)));
    else
      tvalarraypos[i] = (int) (NUMTARRAYBINS * ((tvals[i] - T[i].pr.min) / (T[i].pr.max  - T[i].pr.min))) - tvalarraypos[i - 1];
  }
  if (numt == 2)
  {
    return (double) t2array[tvalarraypos[0]][tvalarraypos[1]] / (double) init;
  }
  if (numt == 3)
  {
    return (double) t3array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]] / (double) init;
  }
  if (numt == 4)
  {
    return (double) t4array[tvalarraypos[0]][tvalarraypos[1]][tvalarraypos[2]][tvalarraypos[3]] / (double) init;
  }
  return -1;                    // should not get here

}                               /* joint_t_prob */
