/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"

/* functions for basic operations on struct genealogy_weights and struct probcalc data structures */

/*********** LOCAL STUFF **********/

/********** GLOBAL FUNCTIONS ***********/

void
init_genealogy_weights (struct genealogy_weights *gweight)
{
  int i;
  gweight->cc = static_cast<int **> (malloc ((numsplittimes + 1) * sizeof (int *)));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->cc[i] = static_cast<int *> (malloc ((npops - i) * sizeof (int)));
  gweight->hcc = static_cast<double **> (malloc ((numsplittimes + 1) * sizeof (double *)));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->hcc[i] = static_cast<double *> (malloc ((npops - i) * sizeof (double)));
  gweight->fc = static_cast<double **> (malloc ((numsplittimes + 1) * sizeof (double *)));
  for (i = 0; i < numsplittimes + 1; i++)
    gweight->fc[i] = static_cast<double *> (malloc ((npops - i) * sizeof (double)));
  
  if (modeloptions[NOMIGRATION] == 0)
  {
    gweight->mc = static_cast<int ***> (malloc (lastperiodnumber * sizeof (int **)));
    gweight->fm = static_cast<double ***> (malloc (lastperiodnumber * sizeof (double **)));

    for (i = 0; i < lastperiodnumber; i++)
    {
      gweight->mc[i] = alloc2Dint (npops - i, npops - i);
      gweight->fm[i] = orig2d_alloc2Ddouble (npops - i, npops - i);
    }
  }
  setzero_genealogy_weights (gweight);
}

void
setzero_genealogy_weights (struct genealogy_weights *gweight)
{
  int i, j, k;
  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      gweight->cc[i][j] = 0;
      gweight->fc[i][j] = 0;
      gweight->hcc[i][j] = 0;
    }
  }
  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)
        {
          gweight->mc[k][i][j] = 0;
          gweight->fm[k][i][j] = 0;
        }
      }
    }
  }
  return;
}

void
free_genealogy_weights (struct genealogy_weights *gweight)
{
  int i;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    XFREE (gweight->cc[i]);
    XFREE (gweight->fc[i]);
    XFREE (gweight->hcc[i]);
  }
  XFREE (gweight->cc);
  XFREE (gweight->fc);
  XFREE (gweight->hcc);

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      orig2d_free2D ((void **) gweight->mc[i], npops - i);     /* orig2d_free2D ((void **) gweight->mc[i], npops - i); */
      orig2d_free2D ((void **) gweight->fm[i], npops - i);     /* orig2d_free2D ((void **) gweight->fm[i], npops - i); */
    }
    XFREE (gweight->mc);
    XFREE (gweight->fm);
  }
  return;
}

void
init_probcalc (struct probcalc *pcalc)
{
  int i;
  pcalc->qintegrate = static_cast<double *> (malloc (numpopsizeparams * sizeof (double)));

  if (nummigrateparams > 0)
  {
    pcalc->mintegrate = static_cast<double *> (malloc (nummigrateparams * sizeof (double)));

  }

  for (i = 0; i < numpopsizeparams; i++)
    pcalc->qintegrate[i] = 0;
  for (i = 0; i < nummigrateparams; i++)
    pcalc->mintegrate[i] = 0;
  pcalc->pdg = 0;
  pcalc->probg = 0;
  pcalc->probhgg = 0;  /*hgstuff */
}

void
free_probcalc (struct probcalc *pcalc)
{
  XFREE (pcalc->qintegrate);
  if (nummigrateparams > 0)
  {
    XFREE (pcalc->mintegrate);
  }
}

// compared memcpy with loops and memcpy did ok
void
copy_treeinfo (struct genealogy_weights *dest, struct genealogy_weights *srce)
{

  int k, i, j;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      dest->cc[i][j] = srce->cc[i][j];
      dest->hcc[i][j] = srce->hcc[i][j];
      dest->fc[i][j] = srce->fc[i][j];
    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        memcpy (dest->mc[k][i], srce->mc[k][i], (npops - k) * sizeof (int));
        memcpy (dest->fm[k][i], srce->fm[k][i], (npops - k) * sizeof (double));
      }
    }
  }
  return;
}

// copies over the things in probcalc  
// used for genealogy as well as RY and NW updates 
// (NOTE: does not do pdg, i.e. prob(data|genealogy), 
// because this needs to get handled separately 
// in different ways depending on context in update_genealogy(),  changet()  
// and changeu() 
// compared memcpy with loops and memcpy did ok

void
copy_probcalc (struct probcalc *dest, struct probcalc *srce)
{

  memcpy (dest->qintegrate, srce->qintegrate, numpopsizeparams * sizeof (double));
  if (nummigrateparams > 0)
  {
    memcpy (dest->mintegrate, srce->mintegrate, nummigrateparams * sizeof (double));
  }
  /* Note that pdg is not copied!!! */
  dest->probg = srce->probg;
  if (hiddenoptions[HIDDENGENEALOGY]==1)
    dest->probhgg = srce->probhgg;
  return;
}

void
sum_treeinfo (struct genealogy_weights *addup,
              struct genealogy_weights *addto)
{

  int i, j, k;
  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      addup->cc[i][j] += addto->cc[i][j];
      addup->fc[i][j] += addto->fc[i][j];
      assert (addup->cc[i][j] == 0 || addup->fc[i][j] > 0);
      assert (addup->fc[i][j] < DBL_MAX);
      addup->hcc[i][j] += addto->hcc[i][j];

    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)

        {
          addup->mc[k][i][j] += addto->mc[k][i][j];
          addup->fm[k][i][j] += addto->fm[k][i][j];
        }
      }
    }
  }

  return;
}

void
sum_subtract_treeinfo (struct genealogy_weights *addup,
                       struct genealogy_weights *addtoplus,
                       struct genealogy_weights *addtominus)
{
  int i, j, k;

  for (i = 0; i < numsplittimes + 1; i++)
  {
    for (j = 0; j < npops - i; j++)
    {
      addup->cc[i][j] += addtoplus->cc[i][j] - addtominus->cc[i][j];
      addup->fc[i][j] -= addtominus->fc[i][j];
      addup->fc[i][j] += addtoplus->fc[i][j];
      addup->fc[i][j] = DMAX(0,addup->fc[i][j]);
      addup->hcc[i][j] -= addtominus->hcc[i][j];
      addup->hcc[i][j] += addtoplus->hcc[i][j];
    }
  }

  if (modeloptions[NOMIGRATION] == 0)
  {
    for (k = 0; k < lastperiodnumber; k++)
    {
      for (i = 0; i < npops - k; i++)
      {
        for (j = 0; j < npops - k; j++)
        {
          addup->mc[k][i][j] += addtoplus->mc[k][i][j] - addtominus->mc[k][i][j];
          addup->fm[k][i][j] -= addtominus->fm[k][i][j];
          addup->fm[k][i][j] += addtoplus->fm[k][i][j];
          addup->fm[k][i][j] = DMAX(0,addup->fm[k][i][j]);
        }
      }
    }
  }
  return;
}                               /*sum_subtract_treeinfo */

/* gsampinf is an array of floats that will hold all of the stuff in allgweight  and pcalc for chain 0
nummigrateparams is determined in setup_iparams()
the sequence in the gsampinf array:
the sequence in this array:

	type     |  # values  |  cumulative total at end
    cc	       numpopsizeparams   numpopsizeparams
	fc	       numpopsizeparams   2*numpopsizeparams
	hcc	       numpopsizeparams   3*numpopsizeparams
	mc         nummigrateparams   3*numpopsizeparams + nummigrateparams
	fm         nummigrateparams   3*numpopsizeparams + 2*nummigrateparams
	qintegrate numpopsizeparams   4*numpopsizeparams + 2*nummigrateparams
	mintegrate nummigrateparams   4*numpopsizeparams + 3*nummigrateparams
	pdg             1        4*numpopsizeparams + 3*nummigrateparams +  1 
	probg           1        4*numpopsizeparams + 3*nummigrateparams +  2
	t          numsplittimes 4*numpopsizeparams + 3*nummigrateparams +  numsplittimes + 2
*/

int
calc_gsampinf_length (void)
{
  int i;
  i = 4 * numpopsizeparams;     //cc, hcc, fc,qintegrate 
  if (!modeloptions[NOMIGRATION]) 
    i += 3 * nummigrateparams;  //mc, fm, mintegrate
  i += numsplittimes;           // times
  i += 2;        // pdg and probg 
  if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
  {
    i += 1; // probhgg
    i += 1; // cmmhg
    if (modeloptions[POPTREETOPOLOGYUPDATE])
      i += 1;// treenumber
  }
  return i;
}

void
savegsampinf (float *g, int z)  // *g is gsampinf, an array for holding genealogy stats
{

  int i, j, c;
  float f, hc;

  // positions where these types (mc,fc,mintegrate, qintegrate) begin in array g
  // AS: Have to figure out which chain is cold before I can do this. So C[0] has to be replaced with C[z]
  // where z is the index of the cold chain, i.e. where beta = 1.0

 // *fpstri += sprintf(&fpstr[*fpstri], 
   // debug_ti_addj0info

  struct genealogy_weights *gweight = &C[z]->allgweight;
  for (i = 0; i < numpopsizeparams; i++)
  {
    c = 0;
    f = hc = 0.0;
    for (j = 0; j < C[z]->itheta[i].wp.n; j++)
    {
      c += (int) gweight->cc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];
      f += (float) gweight->fc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];
      hc += (float) gweight->hcc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];

    }
    g[gsamp_ccp + i] = (float) c;
    assert(g[gsamp_ccp + i] >= 0);
    g[gsamp_fcp + i] = f;
    assert(g[gsamp_fcp + i] >= 0);
    g[gsamp_hccp + i] = hc;
  }
  if (!modeloptions[NOMIGRATION])
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      c = 0;
      f = 0.0;
      for (j = 0; j < C[z]->imig[i].wp.n; j++)
      {
        c +=
          (int) gweight->mc[C[z]->imig[i].wp.p[j]][C[z]->imig[i].wp.r[j]][C[z]->imig[i].wp.c[j]];
        f += (float)
          gweight->fm[C[z]->imig[i].wp.p[j]][C[z]->imig[i].wp.r[j]][C[z]->imig[i].wp.c[j]];
      }
      g[gsamp_mcp + i] = (float) c;
      assert(g[gsamp_mcp + i] >= 0);
      g[gsamp_fmp + i] = f;
      assert(g[gsamp_fmp + i] >= 0);
    }
  }
  for (i = 0; i < numpopsizeparams; i++)
  {
    g[gsamp_qip + i] = (float) C[z]->allpcalc.qintegrate[i];
  }
  for (i = 0; i < nummigrateparams; i++)
  {
    g[gsamp_mip + i] = (float) C[z]->allpcalc.mintegrate[i];
  }
  g[gsamp_pdgp] = (float) C[z]->allpcalc.pdg;
  g[gsamp_probgp] = (float) C[z]->allpcalc.probg;
  for (i = 0; i < numsplittimes; i++)
    g[gsamp_tp + i] = (float) C[z]->tvals[i];
}                               /* savesampinf */


/* used when if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1) */
void
savegsampinf_debug_ti_addinfo (float *g, int z, char *a)
{

  int i, j, c;
  float f, hc;
  int gsamp_hg;
  int gsamp_cmmhg,li;

  int debugtii = 0; // debugging genealogies when sampling trees

  // positions where these types (mc,fc,mintegrate, qintegrate) begin in array g
  // AS: Have to figure out which chain is cold before I can do this. So C[0] has to be replaced with C[z]
  // where z is the index of the cold chain, i.e. where beta = 1.0

 // *fpstri += sprintf(&fpstr[*fpstri], 
   // debug_ti_addj0info

    /*if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
    {
      a[0] = 0;
      debugtii += sprintf(a + debugtii, "\t%.3f\t%d\t",C[z]->allpcalc.probhgg,C[z]->poptreenum);
    } */


  struct genealogy_weights *gweight = &C[z]->allgweight;
  for (i = 0; i < numpopsizeparams; i++)
  {
    c = 0;
    f = hc = 0.0;
    for (j = 0; j < C[z]->itheta[i].wp.n; j++)
    {
      c += (int) gweight->cc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];
      f += (float) gweight->fc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];
      hc += (float) gweight->hcc[C[z]->itheta[i].wp.p[j]][C[z]->itheta[i].wp.r[j]];
    }
    if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
    {
      debugtii += sprintf(a + debugtii,"%s\t",C[z]->itheta[i].str);
    }
    g[gsamp_ccp + i] = (float) c;
    assert(g[gsamp_ccp + i] >= 0);
    g[gsamp_fcp + i] = f;
    assert(g[gsamp_fcp + i] >= 0);
    g[gsamp_hccp + i] = hc;
  }
  if (!modeloptions[NOMIGRATION])
  {
    for (i = 0; i < nummigrateparams; i++)
    {
      c = 0;
      f = 0.0;
      for (j = 0; j < C[z]->imig[i].wp.n; j++)
      {
        c +=
          (int) gweight->mc[C[z]->imig[i].wp.p[j]][C[z]->imig[i].wp.r[j]][C[z]->imig[i].
                                                              wp.c[j]];
        f += (float)
          gweight->fm[C[z]->imig[i].wp.p[j]][C[z]->imig[i].wp.r[j]][C[z]->imig[i].wp.c[j]];
      }

      g[gsamp_mcp + i] = (float) c;
      assert(g[gsamp_mcp + i] >= 0);
      g[gsamp_fmp + i] = f;
      assert(g[gsamp_fmp + i] >= 0);
      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
      {
        debugtii += sprintf(a + debugtii,"%s\t",C[z]->imig[i].str);
      }
    }
  }
  for (i = 0; i < numpopsizeparams; i++)
  {
    g[gsamp_qip + i] = (float) C[z]->allpcalc.qintegrate[i];
  }
  for (i = 0; i < nummigrateparams; i++)
  {
    g[gsamp_mip + i] = (float) C[z]->allpcalc.mintegrate[i];
  }
  g[gsamp_pdgp] = (float) C[z]->allpcalc.pdg;
  g[gsamp_probgp] = (float) C[z]->allpcalc.probg;
  for (i = 0; i < numsplittimes; i++)
    g[gsamp_tp + i] = (float) C[z]->tvals[i];
  if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  
  {
    gsamp_hg = gsamp_tp + numsplittimes;
    g[gsamp_hg] = (float) C[z]->allpcalc.probhgg;
    //g[gsamp_hg] = (float) C[z]->G[0].hgprob;
    gsamp_cmmhg = gsamp_hg +1;
    g[gsamp_cmmhg] = 0.0;
    for (li=0;li<nloci;li++)
      g[gsamp_cmmhg] += (float) C[z]->G[li].gtree->cmmhg;
    if (modeloptions[POPTREETOPOLOGYUPDATE])
      g[gsamp_cmmhg+1]  = (float) C[z]->poptreenum;
    debugtii += sprintf(a + debugtii,"\t");
    for (li=0;li<nloci;li++)
      debugtii += sprintf(a + debugtii,"%.4f\t",C[z]->G[li].hgprob);


  }

}                               /* savegsampinf_debug_ti_addinfo*/

