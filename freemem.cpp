/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"
#include "update_gtree_common.hpp"

/* misc functions for freeing memory at the end */

/********LOCAL PROTOTYPES ********/
static void free_value_record (struct value_record *v);
//static void free_iparam (struct i_param *ip, int n, int m_or_p_or_s); obsolete with hg stuff
static void free_i_params (void);
static void free_locus ();
static void free_T (void);
static void free_mh (void);
static void free_lpgpd_v (struct value_record *l);
static void free_migration_counts (struct value_record **l);
//static void free_migration_counts_times (struct value_record **l);

/**** extern prototypes ******/
extern void free_gtreecommonhg (void);   

/********LOCAL FUNCTIONS ********/
void
free_value_record (struct value_record *v)
{
  if (v->do_xyplot)
    XFREE (v->xy);
  if (v->do_trend)
    XFREE (v->trend);
  return;
}

void
free_lpgpd_v (struct value_record *l)
{
  free_value_record (l);
  XFREE (l);
}

void
free_migration_counts (struct value_record **l)
{
  int i, j;

  for (j = 0; j < nloci + (nloci > 1); j++)
  {
    for (i = 0; i < nummigdirs; i++)
      free_value_record (&l[j][i]);
    XFREE (l[j]);
  }
  XFREE (l);
}


void
free_T (void)
{
  int i;
  for (i = 0; i < lastperiodnumber; i++)
  {
    XFREE (T[i].upnames);
    XFREE (T[i].upinf);
    free_value_record (T[i].v);
    XFREE (T[i].v);
  }
  XFREE (T);
}                               // free_T


void
free_mh (void)
{
  int i;
  for (i = 0; i < nummigrateparams; i++)
  {
    XFREE (mh[i].upnames);
    XFREE (mh[i].upinf);
    free_value_record (mh[i].v);
    XFREE (mh[i].v);
  }
  XFREE (mh);
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    XFREE (mhnit->upnames);
    XFREE (mhnit->upinf);
    XFREE (mhnit->upnames);
    XFREE (mhnit);
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])  // also do popsize terms
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      XFREE (qh[i].upnames);
      XFREE (qh[i].upinf);
      free_value_record (qh[i].v);
      XFREE (qh[i].v);
    }
    XFREE (qh);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      XFREE (qhnit->upnames);
      XFREE (qhnit->upinf);
      XFREE (qhnit->upnames);
      XFREE (qhnit);
    }
  }
}                               // free_mh

void
free_poptreeupinfo (void)
{
  XFREE (poptreeuinfo->upnames);
  XFREE (poptreeuinfo->upinf);
  int j;
  for (j = 0; j < poptreeuinfo->num_vals; j++)
    free_value_record (&(poptreeuinfo->v[j]));
  if (poptreeuinfo->num_vals)
    XFREE (poptreeuinfo->v);
  XFREE (poptreeuinfo);
}                               // free_poptreeupinfo

//m_or_p_or_s refers to migration (1) or population size (0) or split (-1)
//obsolete with hg stuff
/*
void
free_iparam (struct i_param *ip, int n, int m_or_p_or_s)
{
  int i;
  for (i = 0; i < n; i++)
  {
    XFREE (ip[i].xy);
    if (ip[i].wp.n > 0)
    {
      XFREE (ip[i].wp.p);
      XFREE (ip[i].wp.r);
      if (m_or_p_or_s == 1)
        XFREE (ip[i].wp.c);
    }
  }
  if (m_or_p_or_s >= 0)
    XFREE (ip);
  ip = NULL;
  return;
} */                              //free_iparam 

void
free_chainstate_record_updates_and_values (struct
                                           chainstate_record_updates_and_values
                                           *rec, int nrec)
{
  int i, j;
  for (i = 0; i < nrec; i++)
  {
    XFREE ((rec + i)->upnames);
    XFREE ((rec + i)->upinf);
    for (j = 0; j < (rec + i)->num_vals; j++)
      free_value_record (&((rec + i)->v[j]));
    if ((rec + i)->num_vals)
      XFREE ((rec + i)->v);
  }
  XFREE (rec);
}                               //free_chainstate_record_updates_and_values

void
free_locus ()
{
  int li, ai, i;
  for (li = 0; li < nloci; li++)
  {
    if (L[li].model == INFINITESITES || L[li].model == HKY
        || L[li].model == JOINT_IS_SW)
    {
      for (i = 0; i < L[li].numgenes; i++)
        XFREE (L[li].seq[i]);
      XFREE (L[li].seq);
    }
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {
      if (L[li].model == STEPWISE)
        ai = 0;
      else
        ai = 1;
      for (; ai < L[li].nlinked; ai++)
      {
        XFREE (L[li].A[ai]);
      }

      XFREE (L[li].A);
    }
    if (L[li].model == INFINITESITES || L[li].model == JOINT_IS_SW)
      XFREE (L[li].badsite);
    free_chainstate_record_updates_and_values (L[li].u_rec, L[li].nlinked);
    if (L[li].model == HKY)
      free_chainstate_record_updates_and_values (L[li].kappa_rec, 1);
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      free_chainstate_record_updates_and_values (L[li].A_rec, L[li].nlinked);
    free_chainstate_record_updates_and_values (L[li].g_rec, 1);
  }
}                               // free_locus

void free_i_params (void)
{
  int ci,i;
  for (ci=0;ci<numchainspp;ci++)
  {
    for (i = 0; i < numpopsizeparams; i++)
    {
      XFREE(C[ci]->itheta[i].xy);
      XFREE(C[ci]->itheta[i].wp.p);
      XFREE(C[ci]->itheta[i].wp.r);
     
    }
    XFREE(C[ci]->itheta);
    for (i = 0; i < nummigrateparams; i++) 
    {
      XFREE(C[ci]->imig[i].wp.p);
      XFREE(C[ci]->imig[i].wp.r);
      XFREE(C[ci]->imig[i].wp.c);
      XFREE(C[ci]->imig[i].xy); 
    }
    XFREE(C[ci]->imig);
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      XFREE(holdimig);
  }
}  /* free_i_params */

/****** GLOBAL FUNCTIONS  *****/

void
freeanymemory (void)
{
  int ci;
  int i;
  int j;
  int li;
  //int npnodes; no nlonger used
  int nperiods;

  unsetseeds ();
 if (npops == 1)
  {
    nperiods = 1;
  }
  else
  {
    nperiods = npops;
    for (i = 0; i < lastperiodnumber; i++)
    {
      free_value_record (T[i].v);   // is this getting done more than once see free_T()
    }
  }
  
  if (hiddenoptions[HIDDENGENEALOGY]==1)
  {
    free_change_poptree();
/*for (int ix = 0; ix < numpoptopologies; ix++)
{
  assert (strcmp(alltreestrings[ix]+46,"14")==0);
  assert (strlen(alltreestrings[ix])==48);
}     */
		  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)
    {
      freepoptreestringarrays();
      free_poptreeupinfo ();
      XFREE(poptopologyproposedlist); 
      XFREE(poptopologysequence.vals);
      XFREE(poptopologysequence.disvals);
      XFREE(RFtreedis);
    }
    else
    {
      if (hiddenoptions[HIDDENGENEALOGY])
      {
        free_poptreeupinfo ();// this is used for pop branch slide 
      }
    }
  }
/*for (int ix = 0; ix < numpoptopologies; ix++)
{
  assert (strcmp(alltreestrings[ix]+46,"14")==0);
  assert (strlen(alltreestrings[ix])==48);
} */
  free_i_params();  // replace free_iparam for hg stuff
  for (ci = 0; ci < numchainspp; ci++)
  {
    XFREE (C[ci]->tvals);
    XFREE (C[ci]->poptree);
    free_genealogy_weights (&(C[ci]->allgweight));
    free_probcalc (&(C[ci]->allpcalc));
    if (C[ci]->plist != NULL)
    {
      for (j = 0; j < nperiods; j++)
      {
        XFREE (C[ci]->plist[j]);
      }
      XFREE (C[ci]->plist);
    }
    for (li = 0; li < nloci; li++)
    {
      free_genealogy_weights (&(C[ci]->G[li].gweight));
      XFREE (C[ci]->G[li].uvals);
      XFREE (C[ci]->G[li].pdg_a);
      if (L[li].model == HKY)
      {
        XFREE (L[li].mult);
        for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
        {
          XFREE (C[ci]->G[li].gtree[i].hkyi.scalefactor);
          XFREE (C[ci]->G[li].gtree[i].hkyi.oldscalefactor);
          for (j = 0; j < L[li].numsites; j++)
          {
            XFREE (C[ci]->G[li].gtree[i].hkyi.frac[j]);
            XFREE (C[ci]->G[li].gtree[i].hkyi.newfrac[j]);
          }
          XFREE (C[ci]->G[li].gtree[i].hkyi.frac);
          XFREE (C[ci]->G[li].gtree[i].hkyi.newfrac);
        }
      }

      for (i = 0; i < L[li].numlines; i++)
      {
        if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {
          XFREE (C[ci]->G[li].gtree[i].A);
          XFREE (C[ci]->G[li].gtree[i].dlikeA);
        }
      }
      XFREE (C[ci]->G[li].gtree);
    }
    if (doRYupdate)
      XFREE(C[ci]->RYwidthinfo);
    if (doNWupdate)
       XFREE(C[ci]->NWwidthinfo);
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      free_hyperprior_arrays(ci);
    XFREE (C[ci]->G);
    XFREE (C[ci]);
  }
  free_locus ();

/* don't need to do this allocation when hiiden genealogies are used */
  if (hiddenoptions[HIDDENGENEALOGY]==0)
  {
    XFREE (oldedgemig.mtimeavail);
    XFREE (oldedgemig.mp);

    XFREE (oldsismig.mtimeavail);
    XFREE (oldsismig.mp);

    XFREE (newedgemig.mtimeavail);
    XFREE (newedgemig.mp);

    XFREE (newsismig.mtimeavail);
    XFREE (newsismig.mp);
  }
  XFREE (C);
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    free_mh();
    free_migration_prior_update();
  }
  if (npops > 2 && npops <= 5 && outputoptions[PRINTJOINTTEST])
    free_multi_t_arrays ();

  if (npops > 1  && hiddenoptions[HIDDENGENEALOGY]==0)
  {
    free_t_NW ();
    free_t_RY ();
  }
  if(hiddenoptions[HIDDENGENEALOGY]==1)
    free_t_RYhg();
  free_updategenealogy ();
  free_treeweight ();

  if (hiddenoptions[HIDDENGENEALOGY]==0)
    free_gtreecommon ();
  else
    free_gtreecommonhg ();
   XFREE (nnminus1);  
  free_sumlogk ();
  for (i = 0; i < nomigrationchecklist.n; i++)
  {
    XFREE (nomigrationchecklist.p);
    XFREE (nomigrationchecklist.r);
    XFREE (nomigrationchecklist.c);
  }
  if (npops > 1)
    free_T(); /* else no split time */
  free_lpgpd_v (lpgpd_v);
  if (outputoptions[MIGRATEHIST])
    /* 8/26/2011 */
    free_migration_counts (migration_counts);
    //free_migration_counts_times (migration_counts_times);
  free_autoc_pointers ();

  XFREE (L);
  
  //XFREE (nnminus1);
  freeswapstuff();
  if (calcoptions[CALCMARGINALLIKELIHOOD])
    freemarginlikecalc();
  if (calcoptions[LOADPRIORSFROMFILE])
  {
    if (popsizeprior_fromfile != NULL)
      XFREE(popsizeprior_fromfile);
    if (mprior_fromfile != NULL)
      alt2d_free2D(mprior_fromfile);
  }  
  if (runoptions[PRINTBURNTREND] && modeloptions[POPTREETOPOLOGYUPDATE])
    free_burn_phylogeny_counts();

  return;
}                               //freeanymemory
