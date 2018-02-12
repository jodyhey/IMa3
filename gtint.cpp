/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

/* calculations for assessing if one parameter is greater than another 
print matrices of results */ 

/* JH 4/23/2010  extensively revised this file */ 

#undef GLOBVARS
#include "ima.hpp"

static int cci, ccj, wi, wj;
static int treeinc, hitreenum, numtreesused;
static double fci, fcj, hval, denom, qmax, fmi,fmj, mmax;
static double pgt_fcj_gt_0 (double qi);
static double mgt_wj_gt_0 (double mi);
static double  trapzd(double  (*func)(double ), double  a, double  b, int n);
static double  qtrap(double  (*func)(double ), double  a, double  b);
static double gtmig(int mi, int mj);
static double gtpops(int pi, int pj);

#define SWITCH_TO_LOWERGAMMA_CRIT  1e-15   // when the difference between a gamma and uppergamma is smaller than this, use the lowergamma
#define USETREESMAX  20000   // only use this many sampled genealogies for calculations of probabilities that one param is greater than another

double mgt_wj_gt_0 (double mi)
{
  double temp1, temp2, a,b;
  if (mi < MINPARAMVAL)
    return 0.0;
  a = logfact[wj];
  b = uppergamma(wj+1,fmj*mi);
  if (a<=b)
    //IM_err (IMERR_LOGDIFF, " in mgt_wj_gt_0() mi %d a %lf  b %lf",mi,a,b);
    return 0.0;  // sometimes get here for very low mi values, floating point thing 
  if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
  {
    temp1 = lowergamma (wj+1,fmj*mi);
  }
  else
  {
    LogDiff(temp1,a,b);
  }
  temp2 =  wi*log(mi) - fmi*mi -(wj+1)*log(fmj) + temp1;
  temp2 -= denom;
  return exp(temp2);
}

double pgt_fcj_gt_0 (double qi)
{
  double fcj2, fci2,temp1, temp2, temp3, a,b;

  if (qi < MINPARAMVAL)
    return 0.0;
  fcj2 = 2*fcj;
  fci2 = 2*fci;
  if (ccj == 0)
  {
    a = log(qi)-fcj2/qi;
    b = log(fcj2)+ uppergamma(0,fcj2/qi);
    if (a>b)
    {
      LogDiff(temp1,a,b);
      temp2 =  - fci2/qi + cci*log(2/qi) ; 
      temp3  = temp1 + temp2 - hval - denom; 
      return exp(temp3);
    }
    else 
      return 0.0;
  }
  else
  {
    temp1 =  uppergamma(ccj-1,fcj2/qi);
    temp2 = LOG2 + cci*log(2/qi) + (1-ccj)*log(fcj) -fci2/qi; 
    temp3 = temp2 + temp1 - hval - denom;
    return exp(temp3);
  }
}

#define FUNC(x) ((*func)(x))

double  trapzd(double  (*func)(double ), double  a, double  b, int n)
{
   double  x,tnm,sum,del;
   static double  s;
   int it,j;

   if (n == 1) {
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   } else {
      for (it=1,j=1;j<n-1;j++) 
        it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) 
        sum += FUNC(x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

#define EPS 1.0e-4
#define JMAX 20

double  qtrap(double  (*func)(double ), double  a, double  b)
{
   /* void nrerror(char error_text[]); */
   int j;
   double  s,olds;
   olds = -1.0e100;
   for (j=1;j<=JMAX;j++) {
      s=trapzd(func,a,b,j);
      if (j > 5) //Avoid spurious early convergence.
      {
        if (fabs(s-olds) < EPS*fabs(olds) ||(s == 0.0 && olds == 0.0)) 
          return s;
      }
      olds=s;
   }
   return s;
   //nrerror("Too many steps in routine qtrap");
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

double gtmig(int mi, int mj)
{
  int ei;
  double sum, temp, temp1, temp2, temp3, temp4, a,b,c;
  double  (*func)(double );

  func = mgt_wj_gt_0;
  mmax = C[ARBCHAIN]->imig[mi].pr.max;
  sum = 0;
  for (ei = 0; ei < hitreenum; ei+= treeinc) 
  {
    fmi = gsampinf[ei][gsamp_fmp + mi];
    fmj = gsampinf[ei][gsamp_fmp + mj];
    wi = (int) gsampinf[ei][gsamp_mcp + mi];
    wj = (int) gsampinf[ei][gsamp_mcp + mj];
    denom = gsampinf[ei][gsamp_mip + mi]+gsampinf[ei][gsamp_mip + mj];
    if (wj==0) 
    {
      if (fmj>0.0) //
      {
        if (wi > 0) //fmi > 0.0
        {
          a = logfact[wi];
          b = uppergamma(wi+1,(fmi + fmj)*mmax);
          c = uppergamma(wi+1,fmi*mmax);
          if (a<=b || a <= c)
            temp = 0.0;
            //IM_err (IMERR_LOGDIFF, " pos 1 in gtmig() mi %d  mj %d a %lf  b %lf",mi,mj,a,b);
          else
          {
            if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp1 = lowergamma (wi+1,(fmi + fmj)*mmax);
            }
            else
            {
              LogDiff(temp1,a,b);
            }
            temp1 += -(wi+1)*log(fmi + fmj);
            if ((a-c) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp2 = lowergamma (wi+1,fmi*mmax);
            }
            else
            {
              LogDiff(temp2,a,c);
            }
            temp2 += -(wi+1)*log(fmi);
            if (temp2<=temp1)
              temp = 0.0;
              //IM_err (IMERR_LOGDIFF, " pos 3 in gtmig() mi %d  mj %d temp2 %lf  temp1 %lf",mi,mj,temp2,temp1);
            else
            {
              LogDiff(temp3,temp2,temp1);
              temp4 = temp3 - log(fmj) - denom;
              temp = exp(temp4);
            }
          }
        }
        else
        {
          if (fmi > 0.0)
          {
            temp1 = fmi*(exp(-(fmi+fmj)*mmax)- exp(-fmi*mmax) );
            temp2 = fmj*(1-exp(-fmi*mmax));
            temp3 = (temp1 + temp2) / (fmi*fmj*(fmi+fmj));
            temp4 = log(temp3) - denom;
            temp = exp(temp4);
          }
          else // fmi == 0.0 
          {
            temp1 = (fmj*mmax - 1.0 + exp(-fmj*mmax))/(fmj*fmj);
            temp4 = log(temp1) - denom;
            temp = exp(temp4);
          }
        }
      }
      else  //fmj == 0.0 && wj==0
      {
        if (wi > 0) //fmi > 0.0
        {
          a = logfact[wi+1];
          b = uppergamma(wi+2,fmi*mmax);
          if (a<=b)
            temp = 0.0;
            //IM_err (IMERR_LOGDIFF, " pos 4 in gtmig() mi %d  mj %d a %lf  b %lf",mi,mj,a,b);
          else
          {
            if ((a-b) < SWITCH_TO_LOWERGAMMA_CRIT)
            {
              temp1 = lowergamma (wi+2,fmi*mmax);
            }
            else
            {
              LogDiff(temp1,a,b);
            }
            temp2 = -(wi+2) * log(fmi);
            temp3 = temp2 + temp1 - denom;
            temp = exp(temp3);
          }
        }
        else
        {
          if (fmi > 0.0)
          {
            temp1 = (1.0 - exp(-fmi * mmax)*(fmi*mmax + 1.0))/(fmi*fmi);
            temp2 = log(temp1)-denom;
            temp = exp(temp2);
          }
          else // fmi == 0.0 
          {
            temp1 = log(mmax*mmax/2.0) - denom;
            temp = exp(temp1); 
          } 
        }
      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, mmax);
    }
    temp = DMIN(1.0, temp);  // JH 4/22/2010 values should not be greater than 1 , except for small fp issues
    sum += temp;
  }
  sum /= numtreesused;
  return sum;

} // gtmig

/* return the probability that pram pi >param pj */
double gtpops(int pi, int pj)
{
  int ei;
  double sum, temp, temp1, temp2, temp3;
  double  (*func)(double );

  func = pgt_fcj_gt_0;

  qmax = C[ARBCHAIN]->itheta[pi].pr.max;

  sum = 0;
  for (ei = 0; ei < hitreenum; ei+= treeinc)
  {
    cci = (int) gsampinf[ei][gsamp_ccp + pi];
    ccj = (int) gsampinf[ei][gsamp_ccp + pj];
    fci = gsampinf[ei][gsamp_fcp + pi];
    fcj = gsampinf[ei][gsamp_fcp + pj];
    hval = gsampinf[ei][gsamp_hccp + pi] + gsampinf[ei][gsamp_hccp + pj];
    denom = gsampinf[ei][gsamp_qip + pi]+gsampinf[ei][gsamp_qip + pj];

    if (ccj==0)
    {
      if (fcj == 0)
      {
        if (fci==0)
        {
          temp = exp(2.0 * log(qmax)- LOG2 - hval - denom);
        }
        else
        {
          if (cci>= 2)
          {
            temp1 = 2*LOG2 + (2-cci)* log(fci);
            temp2 = uppergamma(cci-2,2*fci/qmax);
            temp3 = temp1 + temp2 - hval-denom;
            temp = exp(temp3);
          }
          else
          {
            if (cci == 1)
            {
              temp1 = 4*fci*exp(uppergamma(0,2*fci/qmax));
              temp2 = 2*qmax*exp(-2*fci/qmax) - temp1;
              temp3 = log(temp2) - hval - denom;
              temp=exp(temp3); // JH added 4/22/2010 
            }
            else  // cci==0
            {
              temp1 = exp(uppergamma(0,2*fci/qmax));  // use for ei() of a negative value 
              temp2 = (qmax/2) * (qmax-2*fci) * exp(-2*fci/qmax) + 2*fci*fci*temp1;
              temp3 = log(temp2) - hval -denom;
              temp=exp(temp3);  // JH added 4/22/2010 
            }
          }
        }
      }
      else //ccj=0 fcj > 0 
      {
        temp = qtrap(func, MINPARAMVAL, qmax);

      }
    }
    else
    {
      temp = qtrap(func, MINPARAMVAL, qmax);
    }
    temp = DMIN(1.0, temp);  // JH 4/22/2010 values should not be greater than 1 except for small fp issues
    sum += temp;
  }
  sum /= numtreesused;
  return sum;
} // gtpops 

#define CHECKGREATERTHAN 0.02  //JH 4/22/2010  to see if reciprocal values sum approximately to 1
void print_greater_than_tests (FILE * outfile)
{
  double **gt_popsize, **gt_mig;
  int i, j;
  int printwarning = 0;
  if (genealogysamples > USETREESMAX)
  {
    treeinc = (int) genealogysamples/(int) USETREESMAX;
    numtreesused = USETREESMAX;
    hitreenum = treeinc * USETREESMAX;
  }
  else
  {
    treeinc = 1;
    hitreenum = numtreesused = genealogysamples;
  }
  gt_popsize= orig2d_alloc2Ddouble (numpopsizeparams, numpopsizeparams);
  for (i=0;i < numpopsizeparams ; i++) //JH 4/22/2010  calculate full matrix
    for (j = 0; j < numpopsizeparams; j++)
    {
      if (i != j && C[ARBCHAIN]->itheta[i].pr.max == C[ARBCHAIN]->itheta[j].pr.max)
        gt_popsize[i][j] = gtpops(i,j); 
      else
        gt_popsize[i][j] = -1;
    } 
  if (!modeloptions[EXPOMIGRATIONPRIOR])
  {
    gt_mig = orig2d_alloc2Ddouble (nummigrateparams, nummigrateparams);
    for (i=0;i < nummigrateparams ; i++) //JH 4/22/2010  calculate full matrix
      for (j = 0; j < nummigrateparams; j++)
      {
        if (i!=j && C[ARBCHAIN]->imig[i].pr.max > MINPARAMVAL && C[ARBCHAIN]->imig[j].pr.max > MINPARAMVAL && C[ARBCHAIN]->imig[i].pr.max ==C[ARBCHAIN]->imig[j].pr.max)
          gt_mig[i][j] = gtmig(i,j);
        else 
          gt_mig[i][j] = -1.0;
      }
  }
  FP "\nPARAMETER COMPARISONS, PROBABILITY THAT ROW PARAMETER IS GREATER THAN COLUMN PARAMETER\n");
  FP "========================================================================================\n");
  if (calcoptions[LOADPRIORSFROMFILE]) 
    FP"    Comparisons only done for parameters with identical prior distributions \n");
  FP "Population Sizes\n");
  for (i=0;i<numpopsizeparams;i++)
    FP "\t%s", C[ARBCHAIN]->itheta[i].str);
  FP "\n");
  for (i=0;i<numpopsizeparams;i++)
  {
    FP "%s", C[ARBCHAIN]->itheta[i].str);
    for (j=0;j<numpopsizeparams;j++)
    {
      if (i == j)
        FP "\t  - ");
      else
      {
        if (gt_popsize[i][j] < 0.0)  
          FP"\tna");
        else
        {
          if ( fabs(1.0 - (gt_popsize[i][j] + gt_popsize[j][i])) >= CHECKGREATERTHAN) 
          {
            FP "\t%.3lf?", gt_popsize[i][j]);  //JH 4/22/2010  calculate full matrix
            printwarning = 1;
          }
          else
            FP "\t%.3lf", gt_popsize[i][j]);  //JH 4/22/2010  calculate full matrix
        }
      }
    }
    FP "\n");
  }
  orig2d_free2D ((void **) gt_popsize, numpopsizeparams);
  FP "\nMigration Rates\n");
  if (modeloptions[EXPOMIGRATIONPRIOR])
    FP"  NOT IMPLEMENTED FOR MIGRATION RATES WITH EXPONENTIAL PRIORS \n");
  else
  {
    for (i=0;i<nummigrateparams;i++)
    {
      FP "\t%s", C[ARBCHAIN]->imig[i].str);
    }
    FP "\n");
    for (i=0;i<nummigrateparams;i++)
    {
      FP "%s", C[ARBCHAIN]->imig[i].str);
      for (j=0;j<nummigrateparams;j++)
      {
        if (i == j)
          FP "\t  - ");
        else
        {
          if (gt_mig[i][j] < 0.0)
            FP"\tna");
          else
          {
            if ( fabs(1.0 - (gt_mig[i][j] + gt_mig[j][i])) >= CHECKGREATERTHAN)
            {
              FP "\t%.3lf?", gt_mig[i][j]);  //JH 4/22/2010  calculate full matrix
              printwarning = 1;
            }
            else
              FP "\t%.3lf", gt_mig[i][j]);  //JH 4/22/2010  calculate full matrix
          }
        }
      }
      FP "\n");
    }
    orig2d_free2D ((void **) gt_mig, nummigrateparams);
  }
  if (printwarning)
    FP"  ""?"" indicates that reciprocal values do not sum to approximately 1, possibly due to a small sample of genealogies\n");
  FP "\n\n");
  return;
}                               // print_greater_than_tests
#undef CHECKGREATERTHAN 

#undef USETREESMAX
#undef SWITCH_TO_LOWERGAMMA_CRIT
