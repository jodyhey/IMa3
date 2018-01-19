/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#undef GLOBVARS
#include "ima.hpp"
/* #include "imagsl.h" */

/*********** LOCAL STUFF **********/

static double *sumlogk[MAXLOCI];

static int makefrac (int ci, int li, int node, double mutrate, double kappa,
                     int e1, int e2, int e3, int e4);
static double getstandfactor (int ci, int li, double kappa);

/******** LOCAL FUNCTIONS ***********/

#if 0
/*  CR 110825.1
 *  pijt() & makefrac_rasmus() are no longer compiled into code, 
 *  they are being kept for historical reasons only, the
 *  functions were replaced by Jim Long's makefrac()
 */

// Rasmus's original makefrac and pijt

double pijt(int ci, int li, double mutrate, double t, double kappa, int from, int to)
    {      /*ACGT*/
    double A, PIj;
    double *pi;
    pi = C[ci]->G[li].pi;

    if  (to==0 || to==2) PIj=pi[0]+pi[2];
    else PIj=pi[1]+pi[3];
    /*A=kappa+1.0; */ /*F84*/
    A=1.0+PIj*(kappa-1.0); /*HKY85*/
    if (from==to)
            return pi[to]+pi[to]*exp(-mutrate*t)*(1.0/PIj-1.0)+exp(-mutrate*t*A)*((PIj-pi[to])/PIj);
    else if (from+to==2 || from+to==4)
            return pi[to]+pi[to]*(1.0/PIj-1.0)*exp(-mutrate*t)-(pi[to]/PIj)*exp(-mutrate*t*A);
    else return pi[to]*(1.0-exp(-mutrate*t));
    }

int makefrac_rasmus(int ci, int li, int node, double mutrate, double kappa, int e1, int e2, int e3, int e4)
    {
     int i, j, k, up1, up2, ret=0;
    double **fracpoint1, **fracpoint2, sum[4];
    struct edge *gtree;

    gtree = C[ci]->G[li].gtree;
    if (node == e1 || node ==  e2 || node ==  e3 || node ==  e4 || e1 == -1) ret=1; 
    up1=gtree[node].up[0];
    if (up1!=-1)
        {
        up2=gtree[node].up[1];
        i=makefrac(ci,li,up1, mutrate, kappa, e1, e2, e3, e4);
        if (i==0) fracpoint1 = gtree[up1].hkyi.frac;
        else {fracpoint1 = gtree[up1].hkyi.newfrac; ret=1;}
        i=makefrac(ci,li,up2, mutrate, kappa, e1, e2, e3, e4);
        if (i==0) fracpoint2 = gtree[up2].hkyi.frac;
        else {fracpoint2 = gtree[up2].hkyi.newfrac; ret=1;}
        if (ret>0)
            {
            if (gtree[up1].up[0]==-1 && gtree[up2].up[0]==-1){
                    for (i=0; i<L[li].numsites; i++)
                            for (j=0; j<4; j++)
                                   gtree[node].hkyi.newfrac[i][j]=pijt(ci,li,mutrate, gtree[up1].time, kappa, j, L[li].seq[up1][i])*pijt(ci, li,mutrate, gtree[up2].time, kappa, j, L[li].seq[up2][i]);
                    }
           else if (gtree[up1].up[0]==-1)
               {
                for (i=0; i<L[li].numsites; i++){
                    for (j=0; j<4; j++)
                        {
                        gtree[node].hkyi.newfrac[i][j]=0;
                        for (k=0; k<4; k++)
                              gtree[node].hkyi.newfrac[i][j]+=pijt(ci, li,mutrate, gtree[up2].time-gtree[gtree[up2].up[0]].time, kappa, j, k)*fracpoint2[i][k];
                        gtree[node].hkyi.newfrac[i][j]=gtree[node].hkyi.newfrac[i][j]*pijt(ci,li,mutrate, gtree[up1].time, kappa, j, L[li].seq[up1][i]);
                        }
                    }
                }
           else 
               if (gtree[up2].up[0]==-1)
                   {
                   for (i=0; i<L[li].numsites; i++)
                       {
                        for (j=0; j<4; j++)
                            {
                            gtree[node].hkyi.newfrac[i][j]=0;
                            for (k=0; k<4; k++)
                                  gtree[node].hkyi.newfrac[i][j]+=pijt(ci, li,mutrate, gtree[up1].time-gtree[gtree[up1].up[0]].time, kappa, j, k)*fracpoint1[i][k];
                            gtree[node].hkyi.newfrac[i][j]=gtree[node].hkyi.newfrac[i][j]*pijt(ci,li,mutrate, gtree[up2].time, kappa, j, L[li].seq[up2][i]);
                            }
                        }
                    }
           else
               {
                for (i=0; i<L[li].numsites; i++){
                    for (j=0; j<4; j++)
                        {
                        sum[j]=0;
                        for (k=0; k<4; k++)
                              sum[j]+=pijt(ci, li,mutrate, gtree[up1].time-gtree[gtree[up1].up[0]].time, kappa, j, k)*fracpoint1[i][k];
                        gtree[node].hkyi.newfrac[i][j]=0;
                        for (k=0; k<4; k++)
                              gtree[node].hkyi.newfrac[i][j]+=pijt(ci,li,mutrate, gtree[up2].time-gtree[gtree[up2].up[0]].time, kappa, j, k)*fracpoint2[i][k];
                        gtree[node].hkyi.newfrac[i][j]=gtree[node].hkyi.newfrac[i][j]*sum[j];
                        }
                    }
                }
            }
        }
    return ret;
    } /* Rasmus's original makefrac */

#endif  // #if 0

// Jim Long's unwrapped makefrac() - much faster than Rasmus's original makefrac() that used recursion 
int
makefrac (int ci, int li, int node, double mutrate, double kappa, int e1,
          int e2, int e3, int e4)
{
  int i, j, up1, up2, ret = 0;
  double **fracpoint1, **fracpoint2;

  /* variables to precompute values */
  int stmp[2];
  double expp[2], kpa, A_j[2], L_to[4], PI_j[2], sum[4], temp[4];
  double minus_t_mutrate1, minus_t_mutrate2, exp1[2], exp2[2], onem[2];
  double A[2], PIj[2], divm[2], mdiv[2], pijt1, pijt2, ptmp[4];
  double *pi;
  struct edge *gtree;
  double max;

  gtree = C[ci]->G[li].gtree;
  pi = C[ci]->G[li].pi;
  kpa = kappa - 1.0;
  if (node == e1 || node == e2 || node == e3 || node == e4 || e1 == -1)
    ret = 1;
  up1 = gtree[node].up[0];
  if (up1 != -1)

  {
    up2 = gtree[node].up[1];
    i = makefrac (ci, li, up1, mutrate, kappa, e1, e2, e3, e4);
    if (i == 0)
    {
      fracpoint1 = gtree[up1].hkyi.frac;
    }
    else
    {
      fracpoint1 = gtree[up1].hkyi.newfrac;
      ret = 1;
    }
    i = makefrac (ci, li, up2, mutrate, kappa, e1, e2, e3, e4);
    if (i == 0)
    {
      fracpoint2 = gtree[up2].hkyi.frac;
    }
    else
    {
      fracpoint2 = gtree[up2].hkyi.newfrac;
      ret = 1;
    }
    if (ret > 0)
    {
      if (gtree[up1].up[0] == -1 && gtree[up2].up[0] == -1)
      {
        minus_t_mutrate1 = -mutrate * gtree[up1].time;
        minus_t_mutrate2 = -mutrate * gtree[up2].time;
        expp[0] = exp (minus_t_mutrate1);
        expp[1] = exp (minus_t_mutrate2);
        for (i = 0; i < L[li].numsites; i++)
        {
          stmp[0] = L[li].seq[up1][i];
          stmp[1] = L[li].seq[up2][i];
          L_to[0] = pi[stmp[0]];
          L_to[1] = pi[stmp[1]];
          if (stmp[0] == 0 || stmp[0] == 2)
            PIj[0] = pi[0] + pi[2];
          else
            PIj[0] = pi[1] + pi[3];

          if (stmp[1] == 0 || stmp[1] == 2)
            PIj[1] = pi[0] + pi[2];
          else
            PIj[1] = pi[1] + pi[3];
          A[0] = 1.0 + PIj[0] * kpa;
          divm[0] = 1.0 / PIj[0] - 1.0;
          mdiv[0] = (PIj[0] - L_to[0]) / PIj[0];
          A[1] = 1.0 + PIj[1] * kpa;
          divm[1] = 1.0 / PIj[1] - 1.0;
          mdiv[1] = (PIj[1] - L_to[1]) / PIj[1];
          max = 0.0;
          for (j = 0; j < 4; j++)
          {
            if (stmp[0] == j)
              pijt1 = L_to[0] + L_to[0] * expp[0] * divm[0] + exp (minus_t_mutrate1 * A[0]) * mdiv[0];
            else if (stmp[0] + j == 2 || stmp[0] + j == 4)
              pijt1 = L_to[0] + L_to[0] * expp[0] * divm[0] - (L_to[0] / PIj[0]) * exp (minus_t_mutrate1 * A[0]);
            else
              pijt1 = L_to[0] * (1.0 - expp[0]);

            if (stmp[1] == j)
              pijt2 = L_to[1] + L_to[1] * expp[1] * divm[1] + exp (minus_t_mutrate2 * A[1]) * mdiv[1];
            else if (stmp[1] + j == 2 || stmp[1] + j == 4)
              pijt2 = L_to[1] + L_to[1] * expp[1] * divm[1] - (L_to[1] / PIj[1]) * exp (minus_t_mutrate2 * A[1]);
            else
              pijt2 = L_to[1] * (1.0 - expp[1]);

            gtree[node].hkyi.newfrac[i][j] = pijt1 * pijt2;
            if (gtree[node].hkyi.newfrac[i][j] > max)
              max = gtree[node].hkyi.newfrac[i][j];
          }
          for (j = 0; j < 4; j++)
            gtree[node].hkyi.newfrac[i][j] =
              gtree[node].hkyi.newfrac[i][j] / max;
          gtree[node].hkyi.scalefactor[i] = log (max);
        }
      }
      else if (gtree[up1].up[0] == -1)
      {
        minus_t_mutrate1 = -mutrate * (gtree[up2].time - gtree[gtree[up2].up[0]].time);
        minus_t_mutrate2 = -mutrate * gtree[up1].time;
        expp[0] = exp (minus_t_mutrate1);
        onem[0] = 1.0 - expp[0];
        expp[1] = exp (minus_t_mutrate2);
        PI_j[0] = pi[0] + pi[2];
        PI_j[1] = pi[1] + pi[3];
        A_j[0] = 1.0 + PI_j[0] * kpa;
        A_j[1] = 1.0 + PI_j[1] * kpa;
        divm[0] = 1.0 / PI_j[0] - 1.0;
        divm[1] = 1.0 / PI_j[1] - 1.0;
        exp1[0] = exp (minus_t_mutrate1 * A_j[0]);
        exp1[1] = exp (minus_t_mutrate1 * A_j[1]);
        exp2[0] = exp (minus_t_mutrate2 * A_j[0]);
        exp2[1] = exp (minus_t_mutrate2 * A_j[1]);
        for (i = 0; i < L[li].numsites; i++)
        {
          max = 0.0;
          gtree[node].hkyi.newfrac[i][0] = 0;
          gtree[node].hkyi.newfrac[i][1] = 0;
          gtree[node].hkyi.newfrac[i][2] = 0;
          gtree[node].hkyi.newfrac[i][3] = 0;
          stmp[0] = L[li].seq[up1][i];
          if (stmp[0] == 0 || stmp[0] == 2)     /* use PI_j[0] & A_j[0] */
          {
            for (j = 0; j < 4; j++)
            {
              if (stmp[0] == j)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[0] + exp2[0] * ((PI_j[0] - pi[stmp[0]]) / PI_j[0]);
              else if (stmp[0] + j == 2 || stmp[0] + j == 4)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[0] - (pi[stmp[0]] / PI_j[0]) * exp2[0];
              else
                ptmp[j] = pi[stmp[0]] * (1.0 - expp[1]);
            }
          }
          else /* use PI_j[1] & A_j[1] */
          {
            for (j = 0; j < 4; j++)
            {
              if (stmp[0] == j)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[1] + exp2[1] * ((PI_j[1] - pi[stmp[0]]) / PI_j[1]);
              else if (stmp[0] + j == 2 || stmp[0] + j == 4)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[1] - (pi[stmp[0]] / PI_j[1]) * exp2[1];
              else
                ptmp[j] = pi[stmp[0]] * (1.0 - expp[1]);
            }
          }
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[0]) / PI_j[0])) * fracpoint2[i][0];
          assert (fracpoint2 != 0x00000000);
          temp[1] = (pi[1] * onem[0]) * fracpoint2[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] - exp1[0] * (pi[2] / PI_j[0])) * fracpoint2[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][0] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][0] *= ptmp[0];
          temp[0] = (pi[0] * onem[0]) * fracpoint2[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] + exp1[1] * ((PI_j[1] - pi[1]) / PI_j[1])) * fracpoint2[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint2[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] - exp1[1] * (pi[3] / PI_j[1])) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][1] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][1] *= ptmp[1];
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] - exp1[0] * (pi[0] / PI_j[0])) * fracpoint2[i][0];
          temp[1] = (pi[1] * onem[0]) * fracpoint2[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[2]) / PI_j[0])) * fracpoint2[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][2] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][2] *= ptmp[2];
          temp[0] = (pi[0] * onem[0]) * fracpoint2[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] - exp1[1] * (pi[1] / PI_j[1])) * fracpoint2[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint2[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] + exp1[1] * ((PI_j[1] - pi[3]) / PI_j[1])) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][3] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][3] *= ptmp[3];
          for (j = 0; j < 4; j++)
            if (gtree[node].hkyi.newfrac[i][j] > max)
              max = gtree[node].hkyi.newfrac[i][j];
          for (j = 0; j < 4; j++)
            gtree[node].hkyi.newfrac[i][j] =
              gtree[node].hkyi.newfrac[i][j] / max;
          gtree[node].hkyi.scalefactor[i] =
            gtree[up2].hkyi.scalefactor[i] + log (max);
        }
      }
      else if (gtree[up2].up[0] == -1)
      {
        minus_t_mutrate1 = -mutrate * (gtree[up1].time - gtree[gtree[up1].up[0]].time);
        minus_t_mutrate2 = -mutrate * gtree[up2].time;
        expp[0] = exp (minus_t_mutrate1);
        onem[0] = 1.0 - expp[0];
        expp[1] = exp (minus_t_mutrate2);
        PI_j[0] = pi[0] + pi[2];
        PI_j[1] = pi[1] + pi[3];
        A_j[0] = 1.0 + PI_j[0] * kpa;
        A_j[1] = 1.0 + PI_j[1] * kpa;
        divm[0] = 1.0 / PI_j[0] - 1.0;
        divm[1] = 1.0 / PI_j[1] - 1.0;
        exp1[0] = exp (minus_t_mutrate1 * A_j[0]);
        exp1[1] = exp (minus_t_mutrate1 * A_j[1]);
        exp2[0] = exp (minus_t_mutrate2 * A_j[0]);
        exp2[1] = exp (minus_t_mutrate2 * A_j[1]);
        for (i = 0; i < L[li].numsites; i++)
        {
          max = 0.0;
          gtree[node].hkyi.newfrac[i][0] = 0;
          gtree[node].hkyi.newfrac[i][1] = 0;
          gtree[node].hkyi.newfrac[i][2] = 0;
          gtree[node].hkyi.newfrac[i][3] = 0;
          stmp[0] = L[li].seq[up2][i];
          if (stmp[0] == 0 || stmp[0] == 2)     /* use PI_j[0] & A_j[0] */
          {
            for (j = 0; j < 4; j++)
            {
              if (stmp[0] == j)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[0] + exp2[0] * ((PI_j[0] - pi[stmp[0]]) / PI_j[0]);
              else if (stmp[0] + j == 2 || stmp[0] + j == 4)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[0] - (pi[stmp[0]] / PI_j[0]) * exp2[0];
              else
                ptmp[j] = pi[stmp[0]] * (1.0 - expp[1]);
            }
          }
          else /* use PI_j[1] & A_j[1] */
          {
            for (j = 0; j < 4; j++)
            {
              if (stmp[0] == j)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[1] + exp2[1] * ((PI_j[1] - pi[stmp[0]]) / PI_j[1]);
              else if (stmp[0] + j == 2 || stmp[0] + j == 4)
                ptmp[j] = pi[stmp[0]] + pi[stmp[0]] * expp[1] * divm[1] - (pi[stmp[0]] / PI_j[1]) * exp2[1];
              else
                ptmp[j] = pi[stmp[0]] * (1.0 - expp[1]);
            }
          }
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[0]) / PI_j[0])) * fracpoint1[i][0];
          temp[1] = (pi[1] * onem[0]) * fracpoint1[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] - exp1[0] * (pi[2] / PI_j[0])) * fracpoint1[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint1[i][3];
          gtree[node].hkyi.newfrac[i][0] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][0] *= ptmp[0];
          temp[0] = (pi[0] * onem[0]) * fracpoint1[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] + exp1[1] * ((PI_j[1] - pi[1]) / PI_j[1])) * fracpoint1[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint1[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] - exp1[1] * (pi[3] / PI_j[1])) * fracpoint1[i][3];
          gtree[node].hkyi.newfrac[i][1] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][1] *= ptmp[1];
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] - exp1[0] * (pi[0] / PI_j[0])) * fracpoint1[i][0];
          temp[1] = (pi[1] * onem[0]) * fracpoint1[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[2]) / PI_j[0])) * fracpoint1[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint1[i][3];
          gtree[node].hkyi.newfrac[i][2] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][2] *= ptmp[2];
          temp[0] = (pi[0] * onem[0]) * fracpoint1[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] - exp1[1] * (pi[1] / PI_j[1])) * fracpoint1[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint1[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] + exp1[1] * ((PI_j[1] - pi[3]) / PI_j[1])) * fracpoint1[i][3];
          gtree[node].hkyi.newfrac[i][3] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][3] *= ptmp[3];
          for (j = 0; j < 4; j++)
            if (gtree[node].hkyi.newfrac[i][j] > max)
              max = gtree[node].hkyi.newfrac[i][j];
          for (j = 0; j < 4; j++)
            gtree[node].hkyi.newfrac[i][j] = gtree[node].hkyi.newfrac[i][j] / max;
          gtree[node].hkyi.scalefactor[i] = gtree[up1].hkyi.scalefactor[i] + log (max);
        }
      }
      else
      {
        minus_t_mutrate1 = -mutrate * (gtree[up1].time - gtree[gtree[up1].up[0]].time);
        minus_t_mutrate2 = -mutrate * (gtree[up2].time - gtree[gtree[up2].up[0]].time);
        expp[0] = exp (minus_t_mutrate1);
        onem[0] = 1.0 - expp[0];
        expp[1] = exp (minus_t_mutrate2);
        onem[1] = 1.0 - expp[1];
        PI_j[0] = pi[0] + pi[2];
        PI_j[1] = pi[1] + pi[3];
        A_j[0] = 1.0 + PI_j[0] * kpa;
        A_j[1] = 1.0 + PI_j[1] * kpa;
        divm[0] = 1.0 / PI_j[0] - 1.0;
        divm[1] = 1.0 / PI_j[1] - 1.0;
        exp1[0] = exp (minus_t_mutrate1 * A_j[0]);
        exp1[1] = exp (minus_t_mutrate1 * A_j[1]);
        exp2[0] = exp (minus_t_mutrate2 * A_j[0]);
        exp2[1] = exp (minus_t_mutrate2 * A_j[1]);
        for (i = 0; i < L[li].numsites; i++)
        {
          max = 0.0;
          sum[0] = 0;
          sum[1] = 0;
          sum[2] = 0;
          sum[3] = 0;
          gtree[node].hkyi.newfrac[i][0] = 0;
          gtree[node].hkyi.newfrac[i][1] = 0;
          gtree[node].hkyi.newfrac[i][2] = 0;
          gtree[node].hkyi.newfrac[i][3] = 0;
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[0]) / PI_j[0])) * fracpoint1[i][0];
          temp[1] = (pi[1] * onem[0]) * fracpoint1[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] - exp1[0] *
                     (pi[2] / PI_j[0])) * fracpoint1[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint1[i][3];
          sum[0] += temp[0] + temp[1] + temp[2] + temp[3];
          temp[0] = (pi[0] + pi[0] * expp[1] * divm[0] + exp2[0] * ((PI_j[0] - pi[0]) / PI_j[0])) * fracpoint2[i][0];
          temp[1] = (pi[1] * onem[1]) * fracpoint2[i][1];
          temp[2] = (pi[2] + pi[2] * expp[1] * divm[0] - exp2[0] * (pi[2] / PI_j[0])) * fracpoint2[i][2];
          temp[3] = (pi[3] * onem[1]) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][0] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][0] *= sum[0];
          temp[0] = (pi[0] * onem[0]) * fracpoint1[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] + exp1[1] * ((PI_j[1] - pi[1]) / PI_j[1])) * fracpoint1[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint1[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] - exp1[1] * (pi[3] / PI_j[1])) * fracpoint1[i][3];
          sum[1] += temp[0] + temp[1] + temp[2] + temp[3];
          temp[0] = (pi[0] * onem[1]) * fracpoint2[i][0];
          temp[1] = (pi[1] + pi[1] * expp[1] * divm[1] + exp2[1] * ((PI_j[1] - pi[1]) / PI_j[1])) * fracpoint2[i][1];
          temp[2] = (pi[2] * onem[1]) * fracpoint2[i][2];
          temp[3] = (pi[3] + pi[3] * expp[1] * divm[1] - exp2[1] * (pi[3] / PI_j[1])) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][1] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][1] *= sum[1];
          temp[0] = (pi[0] + pi[0] * expp[0] * divm[0] - exp1[0] * (pi[0] / PI_j[0])) * fracpoint1[i][0];
          temp[1] = (pi[1] * onem[0]) * fracpoint1[i][1];
          temp[2] = (pi[2] + pi[2] * expp[0] * divm[0] + exp1[0] * ((PI_j[0] - pi[2]) / PI_j[0])) * fracpoint1[i][2];
          temp[3] = (pi[3] * onem[0]) * fracpoint1[i][3];
          sum[2] += temp[0] + temp[1] + temp[2] + temp[3];
          temp[0] = (pi[0] + pi[0] * expp[1] * divm[0] - exp2[0] * (pi[0] / PI_j[0])) * fracpoint2[i][0];
          temp[1] = (pi[1] * onem[1]) * fracpoint2[i][1];
          temp[2] = (pi[2] + pi[2] * expp[1] * divm[0] + exp2[0] * ((PI_j[0] - pi[2]) / PI_j[0])) * fracpoint2[i][2];
          temp[3] = (pi[3] * onem[1]) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][2] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][2] *= sum[2];
          temp[0] = (pi[0] * onem[0]) * fracpoint1[i][0];
          temp[1] = (pi[1] + pi[1] * expp[0] * divm[1] - exp1[1] * (pi[1] / PI_j[1])) * fracpoint1[i][1];
          temp[2] = (pi[2] * onem[0]) * fracpoint1[i][2];
          temp[3] = (pi[3] + pi[3] * expp[0] * divm[1] + exp1[1] *
                     ((PI_j[1] - pi[3]) / PI_j[1])) * fracpoint1[i][3];
          sum[3] += temp[0] + temp[1] + temp[2] + temp[3];
          temp[0] = (pi[0] * onem[1]) * fracpoint2[i][0];
          temp[1] = (pi[1] + pi[1] * expp[1] * divm[1] - exp2[1] * (pi[1] / PI_j[1])) * fracpoint2[i][1];
          temp[2] = (pi[2] * onem[1]) * fracpoint2[i][2];
          temp[3] = (pi[3] + pi[3] * expp[1] * divm[1] + exp2[1] * ((PI_j[1] - pi[3]) / PI_j[1])) * fracpoint2[i][3];
          gtree[node].hkyi.newfrac[i][3] += temp[0] + temp[1] + temp[2] + temp[3];
          gtree[node].hkyi.newfrac[i][3] *= sum[3];
          for (j = 0; j < 4; j++)
            if (gtree[node].hkyi.newfrac[i][j] > max)
              max = gtree[node].hkyi.newfrac[i][j];
          for (j = 0; j < 4; j++)
            gtree[node].hkyi.newfrac[i][j] = gtree[node].hkyi.newfrac[i][j] / max;
          gtree[node].hkyi.scalefactor[i] = gtree[up1].hkyi.scalefactor[i] + gtree[up2].hkyi.scalefactor[i] + log (max);
        }
      }
    }
  }
  return ret;
}                               // Jim Long's unwrapped makefrac() 

double
getstandfactor (int ci, int li, double kappa)
{
  int i, j;
  double p = 0;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      if (i != j)
      {
        if (i + j == 2 || i + j == 4)
          p += C[ci]->G[li].pi[i] * C[ci]->G[li].pi[j] * kappa;
         /*HKY*/
        else
          p += C[ci]->G[li].pi[i] * C[ci]->G[li].pi[j];
      }
    }
  }
  return p;
}

/************ GLOBAL FUNCTIONS ******************/

/* older version with recursion  
  it is a little slower than the one without recursion

void labelgtree( int ci, int li, int edge)
	{
	int dow1, sis;
	struct edge *gtree = C[ci]->G[li].gtree;

	dow1 = gtree[edge].down;
	if (dow1 != -1)
		{
		
		if ((sis = gtree[dow1].up[0]) == edge)
				sis = gtree[dow1].up[1];
		if (gtree[dow1].mut != gtree[edge].mut)
			{
			if (gtree[edge].mut == 2)
				{
				if (gtree[sis].mut != -1)
					gtree[dow1].mut = gtree[sis].mut;
				else 
					gtree[dow1].mut = 2;
				}
			else if ((gtree[sis].mut != -1 && gtree[sis].mut != 2) && (gtree[sis].mut != gtree[edge].mut))
				gtree[dow1].mut = 2;
			else
				gtree[dow1].mut = gtree[edge].mut;
			labelgtree(ci, li, dow1);
			}
		else if (gtree[dow1].mut == 2)
			{
			if (gtree[sis].mut != -1 && gtree[sis].mut != 2)
				{
				gtree[dow1].mut = gtree[sis].mut;	
				labelgtree(ci, li, dow1);
				}
			}
		}
	}	older version labelgtree */

/* version without recursion, suggested by Jim Long as it allows inlining of the function */
//__forceinline void
// inline    just commented this out, was causing problems compiling under linux, not sure if it was very very helpful JH 6/8/2016
// tried profiling with and without inline,  saw no difference 
void labelgtree (int ci, int li, int edge)
{
  int dow1, sis, flag;
  struct edge *gtree = C[ci]->G[li].gtree;
  flag = 1;
  while (flag)
  {
    flag = 0;
    dow1 = gtree[edge].down;
    if (dow1 == -1)
      break;
//    if ((sis = gtree[dow1].up[0]) == edge)// broke this up,  profiling said it was slow 
  //    sis = gtree[dow1].up[1];
    if (gtree[dow1].up[0] == edge)
      sis = gtree[dow1].up[1];
    else
      sis = gtree[dow1].up[0];

    if (gtree[dow1].mut != gtree[edge].mut)
    {
      if (gtree[edge].mut == 2)
      {
        if (gtree[sis].mut != -1)
          gtree[dow1].mut = gtree[sis].mut;
        else
          gtree[dow1].mut = 2;
      }
      else if ((gtree[sis].mut != -1 && gtree[sis].mut != 2)
               && (gtree[sis].mut != gtree[edge].mut))
      {
        gtree[dow1].mut = 2;
      }
      else
      {
        gtree[dow1].mut = gtree[edge].mut;
      }
      flag = 1;
    }
    else if (gtree[dow1].mut == 2)
    {
      if (gtree[sis].mut != -1 && gtree[sis].mut != 2)
      {
        gtree[dow1].mut = gtree[sis].mut;
        flag = 1;
      }
    }
    edge = dow1;
  }
}                               /* labelgtree non-recursive */

double
likelihoodHKY (int ci, int li, double mutrate, double kappa, int e1, int e2,
               int e3, int e4)
{
  int i, j;
  double fracp, p = 0;
  struct edge *gtree = C[ci]->G[li].gtree;
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
    return 0;
  /*JH  5/15/09  fixed bad HKY model bug replaced numsites with totsites */
  mutrate = mutrate / (L[li].totsites * getstandfactor (ci, li, kappa));
  for (i = L[li].numgenes; i < 2 * L[li].numgenes - 1; i++)
    gtree[i].hkyi.newfrac[0][0] = -1;
  makefrac (ci, li, C[ci]->G[li].root, mutrate, kappa, e1, e2, e3, e4);
  for (i = 0; i < L[li].numsites; i++)
  {
    fracp = 0;
    for (j = 0; j < 4; j++)
    {
      fracp += C[ci]->G[li].pi[j] * gtree[C[ci]->G[li].root].hkyi.newfrac[i][j];
    }
    p += L[li].mult[i] * (log (fracp) + gtree[C[ci]->G[li].root].hkyi.scalefactor[i]);
  }
  return p;
}

void
calc_sumlogk (int ci, int li, double *psumlogk)
/* basically a copy of the likelihood function, and probably overkill for what it does */
{
  int ret, node, site, up1, up2 /*, upup */;

#if 0
/*  
 *  CR 110929.3 from JH 9/27/2011
 *  no longer compiled into code, kept for historical reasons only
 */
/* double t;    stopped using.  it was left over from the original 
 *              likelihood function that Rasmus wrote, but it is needed 
 *              to get the sumlogk constant 
 */
#endif  // #if 0
  //int summ = 0, *mutcount, i, j; summ was not being used 
  int *mutcount, i, j;
  double fact, sum;
  struct edge *gtree = C[ci]->G[li].gtree;
  int ng = L[li].numgenes;
  // shouldn't this be sizeof(int) ?? 
  //mutcount = static_cast<int *> (calloc ((size_t) ng, (sizeof (int *)))); bug here,  should not be sizeof(int *) 
  mutcount = static_cast<int *> (calloc ((size_t) ng, (sizeof (int ))));
  for (site = 0; site < L[li].numsites; site++)
  {
    for (j = 0; j < ng; j++)
      gtree[j].mut = L[li].seq[j][site];
    for (j = ng; j < L[li].numlines; j++)
      gtree[j].mut = -1;
    for (j = 0; j < ng; j++)
      labelgtree (ci, li, j);
    ret = -1;
#if 0
/*  
 *  CR 110929.3 from JH 9/27/2011
 *  no longer compiled into code, kept for historical reasons only
 */
    t = 0;
#endif  // #if 0
    for (node = ng; node < L[li].numlines; node++)
    {
      i = node - ng;
      up1 = gtree[node].up[0];
      up2 = gtree[node].up[1];
      if ((gtree[up1].mut == 0 && gtree[up2].mut == 1)
          || (gtree[up1].mut == 1 && gtree[up2].mut == 0))
      {
        if (node != ret && ret != -1)
        {
          XFREE (mutcount);
          return;
        }
        else
        {
#if 0
/*  
 *  CR 110929.3 from JH 9/27/2011
 *  no longer compiled into code, kept for historical reasons only
 */
          /*  stopped using.  this block of code was left 
           *  over from the original likelihood function that Rasmus wrote,
           *  but it is not needed to get the sumlogk constant 
           
          if (gtree[node].down == -1)
          {
            if ((upup = gtree[up1].up[0]) == -1)
              t = gtree[up1].time;
            else
              t = gtree[up1].time - gtree[upup].time;
            if ((upup = gtree[up2].up[0]) == -1)
              t = t + gtree[up2].time;
            else
              t = t + gtree[up2].time - gtree[upup].time;
          }
          else
          {
            if (gtree[gtree[node].down].mut == gtree[up1].mut)
            {
              if ((upup = gtree[up2].up[0]) == -1)
                t = gtree[up2].time;
              else
                t = gtree[up2].time - gtree[upup].time;
            }
            else
            {
              if ((upup = gtree[up1].up[0]) == -1)
                t = gtree[up1].time;
              else
                t = gtree[up1].time - gtree[upup].time;
            }
          }  */
          #endif  // #if 0
          ret = node;
          mutcount[i]++;
          //summ++; not used
        }
      }
    }
    if (ret == -1)
    {
      IM_err(IMERR_INFINITESITESFAIL,"Data incompatible with current genealogy, likelihood calculation failed, chain %d, locus %d, site %d, step %d",ci,li,site, step);
    }
  }
  for (sum = 0.0, i = 0; i < L[li].numlines - ng; i++)
  {
    fact = 1.0;
    for (j = 1; j <= mutcount[i]; j++)
      fact *= (double) j;
    sum += log (fact);
  } 
  *psumlogk = sum;
  XFREE (mutcount);
} //calc_sumlogk 

void
free_sumlogk (void)
{
  int li;

  for (li = 0; li < nloci; li++)
    if (sumlogk[li])
      XFREE (sumlogk[li]);
}

double
likelihoodIS (int ci, int li, double mutrate)
/* written by Rasmus to return the following value
 - (total length of gtree in time) * theta + SUM(log(t_i * theta)*k_i, over all branches)
  where k_i is the number of mutations on branch i. 
  This is equal to the logarithm of the product of a bunch of poisson variables, one for each branch, except that 
  one term is dropped out. [ -SUM(log(k_i !)  ]  is not included.  However this will be the same for two infinite 
  sites gtrees that fit the same data set.
 */
 /* under new parameterization
    - (total length of gtree in time) * mutrate + SUM(log(t_i * mutrate)*k_i, over all branches) - SUM(log(k_i !) 
    where k_i is the number of mutations on branch i. 
    as before, do not need to include  the term : -SUM(log(k_i !) 

    Can also do it as 
    - (total length of gtree in time) * mutrate + log(mutrate) * SUM(k_i) + Sum(log(t_i)*k_i, over all branches) - SUM(log(k_i !) 
    check out the method of getting gtree length - seems ok.
  */
{
  int ret, j, node, site, up1, up2, upup;
  double ptime, p;
  static int done_sumlogk[MAXLOCI] = { 0 };

  int ng = L[li].numgenes;
  struct edge *gtree;
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    return 0;
  }
  if (done_sumlogk[li] == 0)
  {
    sumlogk[li] = static_cast< double *> (calloc (1, sizeof (double)));
    calc_sumlogk (ci, li, sumlogk[li]);
    done_sumlogk[li] = 1;
  }
  p = -C[ci]->G[li].length * mutrate;
  gtree = C[ci]->G[li].gtree;
  for (site = 0; site < L[li].numsites; site++)
  {
    for (j = 0; j < ng; j++)
      gtree[j].mut = L[li].seq[j][site];
    for (j = ng; j < L[li].numlines; j++)
      gtree[j].mut = -1;
    for (j = 0; j < ng; j++)
      labelgtree (ci, li, j);
    ret = -1;
    for (node = ng; node < L[li].numlines; node++)
    {
      up1 = gtree[node].up[0];
      up2 = gtree[node].up[1];
      if ((gtree[up1].mut == 0 && gtree[up2].mut == 1)
          || (gtree[up1].mut == 1 && gtree[up2].mut == 0))
      {
        if (node != ret && ret != -1)
        {
          return FORCEREJECTIONCONSTANT;
        }
        else
        {
          /* i is the root node add together both branch times */
          if (gtree[node].down == -1)
          {
            if ((upup = gtree[up1].up[0]) == -1)
              ptime = gtree[up1].time;
            else
              ptime = gtree[up1].time - gtree[upup].time;
            if ((upup = gtree[up2].up[0]) == -1)
              ptime = ptime + gtree[up2].time;
            else
              ptime = ptime + gtree[up2].time - gtree[upup].time;
          }
          else /* i is not the root node  */
          {
            if (gtree[gtree[node].down].mut == gtree[up1].mut)
            {
              if ((upup = gtree[up2].up[0]) == -1)
                ptime = gtree[up2].time;
              else
                ptime = gtree[up2].time - gtree[upup].time;
            }
            else
            {
              if ((upup = gtree[up1].up[0]) == -1)
                ptime = gtree[up1].time;
              else
                ptime = gtree[up1].time - gtree[upup].time;
            }
          }
          ret = node;
        }
      }
    }
    if (ret == -1)
    {
      IM_err(IMERR_INFINITESITESFAIL,"Data incompatible with current genealogy, likelihood calculation failed, chain %d, locus %d, site %d, step %d",ci,li,site, step);
    }
    else
    {
      assert (ptime * mutrate != 0);
      p = p + log (ptime * mutrate);
    }
  }
  p -= *sumlogk[li];
//std::cout << "likelihoodIS returned is " << p << "\n";
  return p;
}                               /* likelihood */

/* likelihoodSW() calculates and  sums the dlikeA values on each branch */
/* if a locus is joint then ai goes from 0 to nlinked - 1  otherwise ai goes from 0 to nlinked */
double
likelihoodSW (int ci, int li, int ai, double u, double tr)
{
  int i, up, d;
  double t, like;
  struct edge *gtree; 
  double bessiv;
  int iszero; 

  gtree = C[ci]->G[li].gtree;
  like = 0.0;
  iszero = 0;
  for (i = 0; i < L[li].numlines; i++)
  {
    if (gtree[i].down != -1)
    {
      assert (gtree[i].A[ai] >= L[li].minA[ai]
              && gtree[i].A[ai] <= L[li].maxA[ai]);
      d = gtree[i].A[ai] - gtree[gtree[i].down].A[ai];
      if (gtree[i].up[0] == -1)
      {
        t = gtree[i].time * tr;
      }
      else
      {
        up = gtree[i].up[0];
        t = tr * (gtree[i].time - gtree[up].time);
      }

      /* SANGCHUL: Wed Apr 15 17:40:24 EDT 2009 
       * When d is positive,
       * u is very small, bessi(d,t*u) would return 0, which makes
       * log infinity. If t * u is too small then, we have zero probability
       * of likelihood and we should reject any proposal. */
      bessiv = bessi (d, t * u);
      if (!(bessiv > 0.0))
      {
        iszero = 1;
        gtree[i].dlikeA[ai] = IM_BESSI_MIN;
      }
      else
      {
        gtree[i].dlikeA[ai] = -(t * u) + log (bessiv);
      }
      like += gtree[i].dlikeA[ai];
    }
    else
    {
      gtree[i].dlikeA[ai] = 0.0;
    }
  }

  /* SANGCHUL: We return the smallest value for the likelihood */
  if (iszero == 1)
  {
    like = -DBL_MAX;
  }

  /* SANGCHUL: Wed Apr 15 17:42:29 EDT 2009
   * Can this statement be called before the for-loop above? */
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])  
  {
    return 0.0;
  }
  else
  {
    return like;
  }
}                               /* likelihoodSW */

#ifdef TURNONCHECKS
/* use this for debugging, to make sure likelihoodSW has set the dlikeA values properly */
void
checklikelihoodSW (int ci, int li, int ai, double u)
{
  int i, up, d;
  double t, dlikeA;
  struct edge *gtree;
  double bessiv;

  gtree = C[ci]->G[li].gtree;
  for (i = 0; i < L[li].numlines; i++)
  {
    if (gtree[i].down != -1)
    {
      d = gtree[i].A[ai] - gtree[gtree[i].down].A[ai];
      if (gtree[i].up[0] == -1)
      {
        t = gtree[i].time;
      }
      else
      {
        up = gtree[i].up[0];
        t = gtree[i].time - gtree[up].time;
      }
      bessiv = bessi (d, t * u);
      if (!(bessiv > 0.0))
      {
        IM_err (IMERR_SWCHECK, 
                "Negative or zero value (%lf) of Bessel function", bessiv);
      }
      dlikeA = -(t * u) + log (bessi (d, t * u));
      if (gtree[i].dlikeA[ai] != dlikeA)
      {
        IM_err (IMERR_SWCHECK, "Chain %d, Locus %d", ci, li);
      }
    }
  }
  return;
}                               /* checklikelihoodSW */
#endif //TURNONCHECKS