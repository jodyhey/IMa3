/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"
/* calculate marginal likelihood */

/*
implement thermodynamic integration over a large number of intervals

we want the marginal likelihood under the model

p(D) 

for p_Bi(D|G)

where Bi is a heating value i,  for a total of  j values   0<i<j-1

Let L_i  be the mean of p_Bi(D|G) sampled over the course of the run

L_i = Sum[p_Bi(D|G)]/k    for k samples

Then 

p(d) = Sum[ (Bi-B(i-1)) (L_i + L_(i-1))/2, {i, 1, j-1}] 

this is trapezoidal rule 

Also record the harmonic mean  - have to use eexp()

p(d) = 1/ [ SUM[ 1/(p(D|G)]/k ]

*/
#define LOG_10_2  0.30102999566398119521
#define OCUTOFF  10

//double thermosum[MAXCHAINS]; 


double *thermosum,*thermosum_rec;
//double *thermosum2,*thermosum2_rec;
//double *stepstone1_sum,*stepstone1_sum_rec;
//double *stepstone1_sum2,*stepstone1_sum2_rec;
//double *stepstone2_sumA,*stepstone2_sumA_rec;
//double *stepstone2_sumB,*stepstone2_sumB_rec;
//double *stepstone_Lmax,*stepstone_Lmax_rec;
double betawidth;
#define  ifromk(k,K)  ((K) - (k))  // chains are ordered from high beta to low,  but some marginal likelihood estimators use the reverse,  so this little macro makes it easier to work with those


void initmarginlikecalc()
{
  int i;
  thermosum = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  thermosum_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //thermosum2 = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //thermosum2_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
//  stepstone1_sum = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
//  stepstone1_sum_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
//  stepstone1_sum2 = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
//  stepstone1_sum2_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone2_sumA = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone2_sumA_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone2_sumB = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone2_sumB_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone_Lmax = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  //stepstone_Lmax_rec = static_cast<double *> (malloc ((numchainstotal) * sizeof (double)));
  for (i=0;i<numchainstotal;i++)
  {
	  thermosum[i] = 0.0;
	  thermosum_rec[i] = 0.0;
//	  thermosum2[i] = 0.0;
//	  thermosum2_rec[i] = 0.0;
//	  stepstone1_sum[i] = 0.0;
//	  stepstone1_sum_rec[i] = 0.0;
//	  stepstone1_sum2[i] = 0.0;
//	  stepstone1_sum2_rec[i] = 0.0;
//	  stepstone2_sumA[i] = 0.0;
//	  stepstone2_sumA_rec[i] = 0.0;
//	  stepstone2_sumB[i] = 0.0;
//	  stepstone2_sumB_rec[i] = 0.0;
//	  stepstone_Lmax[i] = -DBL_MAX;
//	  stepstone_Lmax_rec[i] = -DBL_MAX;
  }
  betawidth = 1.0 / (numchainstotal-1);

}

void freemarginlikecalc(void)
{
	XFREE(thermosum);
	XFREE(thermosum_rec);
	//XFREE(thermosum2);
	//XFREE(thermosum2_rec);
	//XFREE(stepstone1_sum);
	//XFREE(stepstone1_sum_rec);
	//XFREE(stepstone1_sum2);
	//XFREE(stepstone1_sum2_rec);
	//XFREE(stepstone2_sumA);
	//XFREE(stepstone2_sumA_rec);
	//XFREE(stepstone2_sumB);
	//XFREE(stepstone2_sumB_rec);
	//XFREE(stepstone_Lmax);
	//XFREE(stepstone_Lmax_rec);

}

void move_calcmarglike_vals(void)
{
  int i;
  int rc;
// move values from thermosum into thermosum_rec
  if (numprocesses > 1)
  {
#ifdef MPI_ENABLED
    rc = MPI_Reduce(thermosum,thermosum_rec,numchainstotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
/*
    rc = MPI_Reduce(thermosum2,thermosum2_rec,numchainstotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(stepstone_Lmax,stepstone_Lmax_rec,numchainstotal, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(stepstone2_sumA,stepstone2_sumA_rec,numchainstotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    rc = MPI_Reduce(stepstone2_sumB,stepstone2_sumB_rec,numchainstotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
 */
#endif  // MPI_ENABLED
  }
  else
  {
    for (i=0;i<numchainstotal;i++)
    {
      thermosum_rec[i] = thermosum[i];
/*      thermosum2_rec[i] = thermosum2[i];
      stepstone_Lmax_rec[i] = stepstone_Lmax[i];
      stepstone2_sumA_rec[i] =   stepstone2_sumA[i];
      stepstone2_sumB_rec[i] =   stepstone2_sumB[i];
 */
    }
  }
  /* old code that was here sorted thermosum,  but now that using currallbetapos for indexing thermosum should have right values in right places */
}


 /* JH 1/5/2012  added stepstone_Lmax initialization for steppingstone sampling for calculation of marginal likelihood 
 this is used to record the largest (in this case most negative ??) value of the likelihood for each chain during the burnin
void stepstone_get_Lmax()
{
  int ci;
  for (ci=0;ci<numchainspp;ci++)
  {
    if (stepstone_Lmax[C[ci]->currallbetapos] < C[ci]->allpcalc.pdg)
      stepstone_Lmax[C[ci]->currallbetapos] = C[ci]->allpcalc.pdg;
  }
} */


void summarginlikecalc()
{
  int ci,cc;
  //double betawidth;
  //double temp;
  //int rc;
  // add the likelihood from chain ci to its appropriate accumulator in thermosum    
#ifdef MPI_ENABLED
   /*rc = MPI_Allreduce(stepstone_Lmax,stepstone_Lmax_rec,numchainstotal, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   if (rc != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
    */
#endif

  for (ci = 0;ci<numchainspp;ci++)
  {
    cc = C[ci]->currallbetapos;
    thermosum[cc] += C[ci]->allpcalc.pdg;
    //thermosum2[cc] += SQR(C[ci]->allpcalc.pdg);
    //betawidth = allbetas[cc] - allbetas[cc+1]; 
    //assert(betawidth >= 0.0);
    //stepstone1_sum[i_k] += exp(betawidth * (C[i_kminus1]->allpcalc.pdg - stepstone_Lmax[i_k]));
    //stepstone1_sum2[i_k]  +=  exp(2*betawidth * (C[i_kminus1]->allpcalc.pdg - stepstone_Lmax[i_k]));
    //temp = exp(-betawidth * (C[ci]->allpcalc.pdg - stepstone_Lmax_rec[cc]));
    //if (isnotinf_DBL(temp))
      //stepstone2_sumA[cc] += temp;
    //temp =  exp(-2*betawidth * (C[ci]->allpcalc.pdg - stepstone_Lmax_rec[cc]));
    //if (isnotinf_DBL(temp))
      //stepstone2_sumB[cc] += temp;
  } 

} 

/* should not bother with this as it does not work well */
double harmonicmarginlikecalc(void)
{
  double hmlog;
  int gi, zadj;
  int pdgp;
  int tempz;
  double tempm;
  struct extendnum *harmonicsump;
  int maxz = 0;
  double  harmonicsum_eexp = 0.0;

  harmonicsump = (struct extendnum *) malloc ((size_t) ((genealogysamples + 1) * sizeof (struct extendnum)));

  pdgp = 4*numpopsizeparams + 3* nummigrateparams;  // position in gsampinf[gi] that holds pdg
  for (gi = 0; gi < genealogysamples; gi++)
  {
    eexp(gsampinf[gi][pdgp],&tempm,&tempz);
    harmonicsump[gi].m = 1.0/tempm;
    harmonicsump[gi].z = -tempz;
    if (harmonicsump[gi].z > maxz)
      maxz = harmonicsump[gi].z;
  }
  for (gi = 0; gi < genealogysamples; gi++)
  {
    zadj = harmonicsump[gi].z - (maxz - OCUTOFF);
    harmonicsum_eexp += harmonicsump[gi].m * pow (10.0, (double) zadj);
  }
  hmlog = -(log (harmonicsum_eexp) - log( (double) genealogysamples) + (maxz - OCUTOFF) * LOG10);
  XFREE(harmonicsump);
  return hmlog;
}

/* does simpson's rule, but this requires even interval widths 
void thermo_marginlike_calc_hold(int n, double *estimator, double *stdev_d, double *stdev_s)
{
  int i;
  double sum, varsum; 
 

  sum = 0.0;
  varsum = 0.0;
 
  for (i = 0;i<numchainstotal-1;i++)  //do not use the last chain 
  {
    thermosum_rec[i] /= n;
    thermosum2_rec[i]  /= n;
    if (i==0 ) // || i == numchainstotal-1  do not use the last chain )
      varsum +=  (thermosum2_rec[i] - SQR(thermosum_rec[i]))/4.0;
    else
      varsum +=  thermosum2_rec[i] - SQR(thermosum_rec[i]);
  }
  *stdev_d = fabs(thermosum_rec[0] - thermosum_rec[numchainstotal-2])/ (double) numchainstotal;  // use second to last chain, do not use the last chain 
  *stdev_s = sqrt(varsum/SQR(numchainstotal));
  // Simpson's rule  - must ensure previously that the number of intervals is even
  // so the number of chains must  be odd.   but since we are counting from zero,  the value of numchainstotal must  even 
  for (i = 0;i<=numchainstotal-1;i+=2)
  {
    if (i != (numchainstotal - 1))  // for chain[numchainstotal-1] use a value of 0.0 because under thermodynamic integration the mean likelihood at Beta=0 is 0.0 
      sum += 4.0 * thermosum_rec[i];
  }
  for (i = 1;i<=numchainstotal-2;i+=2)
  {
    sum += 2.0 * thermosum_rec[i];
  }
  *estimator = betawidth * sum/3.0; 
} */

/* trapezoid rule,  allows variable width  
this can cause problems with the last chain with beta = 0
see setheat(),  where this last chain is checked to have beta > 0 
with doing thermodynamic integration
  */
void thermo_marginlike_calc(int n, double *estimator)
{
  int i;
  double sum; 
  double betawidth;

  sum = 0.0;
  for (i = 0;i<=numchainstotal-1;i++) 
  {
    thermosum_rec[i] /= n;
  }
  /* trapezoid rule*/
  for (i = 0;i<=numchainstotal-2;i+=1)
  {
    // had bad bug here,  was using abs instead of fabs. mvc++ handled it but not gcc 
    betawidth = fabs(allbetas[i] - allbetas[i+1]);  //fabs is here because some bug turns up negative values for betawidth
    sum += betawidth * (thermosum_rec[i] + thermosum_rec[i+1]);
  }
  *estimator = sum/2.0;
}

void thermo_marginlike_calc_hold(int n, double *estimator, double *stdev_d, double *stdev_s)
{
  int i;
  double sum; 
  double betawidth;
  //double varsum;

  sum = 0.0;
  //varsum = 0.0;
 
  for (i = 0;i<=numchainstotal-2;i++)  /*do not use the last chain */
  {
    thermosum_rec[i] /= n;
    /*
    thermosum2_rec[i]  /= n;
    if (i==0 )
      varsum +=  (thermosum2_rec[i] - SQR(thermosum_rec[i]))/4.0;
    else
      varsum +=  thermosum2_rec[i] - SQR(thermosum_rec[i]);
     */
  }
  //resuse the second to last chain as if it were last
  i = numchainstotal-2;
  /*varsum +=  (thermosum2_rec[i] - SQR(thermosum_rec[i]))/4.0;
  *stdev_d = fabs(thermosum_rec[0] - thermosum_rec[numchainstotal-2])/ (double) numchainstotal;  // use second to last chain, do not use the last chain 
  *stdev_s = sqrt(varsum/SQR(numchainstotal)); */
  /* trapezoid rule, but for chain with beta = 0 use value for second to last chain */
  for (i = 0;i<=numchainstotal-2;i+=2)
  {
    // had bad bug here,  was using abs instead of fabs. mvc++ handled it but not gcc 
    betawidth = fabs(allbetas[i] - allbetas[i+1]);  //fabs is here because some bug turns up negative values for betawidth
    if (i== numchainstotal-2)
      sum += betawidth * 2 * thermosum_rec[i];
    else
    sum += betawidth * (thermosum_rec[i] + thermosum_rec[i+1]);

  }
  *estimator = sum/2.0;
}  // thermo_marginlike_calc_hold

/* works using the alternative stepping stone method that Aude develoepd*/
/* not being used,  as not fully debugged
difficult to use because it needs a value for L_max  which can't be known until after a long run

void steppingstone_marginlike_calc2(int n,double *estimator2, double *stdev2)
{
  int K,k,i_k;
  double logA,logB, sum,varsum; 
  double betawidth;

  // use a loop over k  where the unheated chain is chain K 
  //  this helps keep things in line with notation in the Marginal_likelihood_estimation.nb  notebook
  K = numchainstotal - 1;
  sum = 0.0;
  varsum = 0.0;
  for (k=1;k<= K; k++)
  {
    i_k = ifromk(k,K);  // goes from numchainstotal-2  to 0 
    betawidth = allbetas[i_k] - allbetas[i_k+1]; 
    assert(betawidth >= 0.0);
    logA = betawidth * stepstone_Lmax_rec[i_k] - log(stepstone2_sumA_rec[i_k]);
    logB = 2*betawidth * stepstone_Lmax_rec[i_k] - log(stepstone2_sumB_rec[i_k]);
    sum += logA;
//    printf("logA  %lf  logB  %lf   exp(loga-2logB)  %lf\n",logA,logB,exp(logA-2*logB));
    varsum += exp(logA-2*logB);
  }
  *estimator2 = sum + (double) numchainstotal* log((double) n);
  *stdev2 = sqrt( -((double) numchainstotal/ (double) n) + varsum);
} */