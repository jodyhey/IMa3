/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
/* modified from standalone code Yujin wrote */ 

#include "ima.hpp"
//#include <iostream>  
//#include <math.h>       // log10, floor, isnan
#include <algorithm>    // std::max
#include <vector>
//#include <fstream>
//#include "geweke.hpp"


class Eureka
{
private:
  std::vector<std::vector<double> > coefs;
  std::vector<double> var;
  
public:
  std::vector<std::vector<double> > get_coefs(){return coefs;}
  std::vector<double> get_var(){return var;}
  void eureka(int lr, std::vector<double> r, std::vector<double> g); //, std::vector<std::vector<double> > f, std::vector<double> var, std::vector<double> a);
};


class AR
{
private:
  double var_pred;
  std::vector<double> ar_coef;
  int ar_nan;
public:
  int get_ar_nan(){return ar_nan;}
  double get_var_pred(){return var_pred;}
  std::vector<double> get_ar_coef(){return ar_coef;}
  void ar(std::vector<double> vars);
};

// Copied from R STATS package, filter.c
/* now allows missing values */
static void
acf0(double *x, int n, int ns, int nl, int correlation, double *acf)
{
    int d1 = nl+1, d2 = ns*d1;

    for(int u = 0; u < ns; u++)
        for(int v = 0; v < ns; v++)
            for(int lag = 0; lag <= nl; lag++) {
                double sum = 0.0; int nu = 0;
                for(int i = 0; i < n-lag; i++)
                    if(!isnan_(x[i + lag + n*u]) && !isnan_(x[i + n*v]))
		      {
                        nu++;
                        sum += x[i + lag + n*u] * x[i + n*v];
		      }
                //acf[lag + d1*u + d2*v] = (nu > 0) ? sum/(nu + lag) : double.Nan;//nan("");
                if (nu>0)
                      acf[lag + d1*u + d2*v] = sum/(nu + lag);
                    else
                    {
                      double b0 = 0.0;
                      acf[lag + d1*u + d2*v] = 0.0/b0;
                    }
            }
    /*
    if(correlation) {
        if(n == 1) {
            for(int u = 0; u < ns; u++)
                acf[0 + d1*u + d2*u] = 1.0;
        } else {
            double *se = (double *) R_alloc(ns, sizeof(double));
            for(int u = 0; u < ns; u++)
                se[u] = sqrt(acf[0 + d1*u + d2*u]);
            for(int u = 0; u < ns; u++)
                for(int v = 0; v < ns; v++)
                    for(int lag = 0; lag <= nl; lag++) { // ensure correlations remain in  [-1,1] :
                        double a = acf[lag + d1*u + d2*v] / (se[u]*se[v]);
                        acf[lag + d1*u + d2*v] = (a > 1.) ? 1. : ((a < -1.) ? -1. : a);
                    }
        }
    }
    */
    return;
}

// acf <- .Call(C_acf, x, lag.max, type == "correlation")
std::vector<double> SEXP_acf(std::vector<double> x, double lmax) // , SEXP sCor)
{
  int nx = x.size(), ns = 1, lagmax = (int) lmax;
  //        cor = asLogical(sCor);
  //  x = PROTECT(coerceVector(x, REALSXP));
  //  SEXP ans = PROTECT(allocVector(REALSXP, (lagmax + 1)*ns*ns));
  double *ans, *xp;
  ans  = (double *)malloc(sizeof(double)*(lagmax + 1)*ns*ns);
  xp = (double *)malloc(sizeof(double)*x.size());
  for(unsigned int i=0; i< x.size(); i++)
    xp[i] = x.at(i);
  acf0(xp, nx, ns, lagmax, 1, ans);
  /*
    SEXP d = PROTECT(allocVector(INTSXP, 3));
    INTEGER(d)[0] = lagmax + 1;
    INTEGER(d)[1] = INTEGER(d)[2] = ns;
    setAttrib(ans, R_DimSymbol, d);
    UNPROTECT(3);
  */

  std::vector<double> res;
  for(unsigned int i=0; i< (lagmax + 1)*ns*ns; i++)
    res.push_back(ans[i]);
  XFREE(ans);
  XFREE(xp);
  
  return res;
}

// autocovariance or autocorrelation function
// returns An array with the same dimensions as ‘lag’ containing the
//          estimated acf
// acf(x, type = "covariance", lag.max = order.max, plot = FALSE, demean = demean)$acf
std::vector<double> acf(std::vector<double> x, double order_max)
{
  std::vector <double> acf;

  
  // acf <- .Call(C_acf, x, lag.max, type == "correlation")
  acf = SEXP_acf(x,order_max);

  /*
  for(int i=0; i<acf.size(); i++)
    std::cout <<"i = " <<i <<" acf.at(i) = " << acf.at(i) <<"\n";
  */
  return acf;
}


/*
z <- .Fortran(C_eureka, as.integer(order.max), r, r,
                      coefs = double(order.max^2),
                      vars = double(order.max),
                      double(order.max))
*/
void Eureka::eureka(int lr, std::vector<double> r, std::vector<double> g) // , std::vector<std::vector<double> > f, std::vector<double> var, std::vector<double> a)
{
/*
  f: coefs
  vars: variance matrix
 */
//c
//c      solves Toeplitz matrix equation toep(r)f=g(1+.)
//c      by Levinson's algorithm
//c      a is a workspace of size lr, the number
//c      of equations
//c
// double r(lr + 1), g(lr + 1), f(lr, lr), a(lr), var(lr);

  /*
  std::cout <<"\nIn eureka()\n";
  std::cout <<" lr = "<< lr <<"\n";
  for(int i=0; i< r.size();i++)
    std::cout <<"i=" <<i <<" r.at(i) ="<< r.at(i)<<"\n";
  */  

  double v = r.at(0); //r(1);
  double d = r.at(1); // r(2);

  std::vector<double> a; a.resize(lr);
  var.resize(lr);
  coefs.resize(lr);
  for(int i=0; i<lr; i++)
    coefs.at(i).resize(lr);  
  
  a.at(0)=1; // a(1) = 1.0d0;
  coefs.at(0).at(0)=g.at(1)/v; // f(1, 1) = g(2)/v;
  double q = coefs.at(0).at(0)*r.at(1); //q = f(1, 1)*r(2);
  var.at(0) = (1-coefs.at(0).at(0)*coefs.at(0).at(0))*r.at(0); //var(1) = (1 - f(1, 1)*f(1, 1))*r(1);
  // std::cout <<"var.at(0) = "<< var.at(0)<<"\n";
  
  if (lr  ==  1) // return;
    return;
    // std::cout <<"STOP: lr =1\n";
  for(int l=2; l<=lr; l++) //for(l=2; l<=lr; l++)
    {
      a.at(l-1) = -d/v; // a(l) = -d/v;
      if (l  >  2) 
	{
	  int l1 = (l - 2)/2;
	  int  l2 = l1 + 1;
	  for(int j=2; j<=l2; j++)
	    {
	      double hold = a.at(j-1); //hold = a(j);
	      int k = l - j + 1;
	      a.at(j-1) = a.at(j-1) + a.at(l-1)*a.at(k-1);// a(j) = a(j) + a(l)*a(k);
	      a.at(k-1) = a.at(k-1) + a.at(l-1)*hold;//  a(k) = a(k) + a(l)*hold;
            
	    }
	  if (2*l1  !=  l - 2)
	    a.at(l2) = a.at(l2 )*(1 + a.at(l-1)); //a(l2 + 1) = a(l2 + 1)*(1.0d0 + a(l));
	}
      v = v + a.at(l-1)*d;// v = v + a(l)*d;
      coefs.at(l-1).at(l-1) = (g.at(l) - q)/v; //f(l, l) = (g(l + 1) - q)/v;
      for(int j = 1; j <= l - 1; j++)// for(40 j = 1, l - 1)\n{
	{
	  coefs.at(l-1).at(j-1) = coefs.at(l - 2).at(j-1) + coefs.at(l-1).at(l-1)*a.at(l - j);     // f(l, j) = f(l - 1, j) + f(l, l)*a(l - j + 1);        
	}
      //c  estimate the innovations variance
      var.at(l-1) = var.at(l - 2)*(1 - coefs.at(l-1).at(l-1)*coefs.at(l-1).at(l-1)); // var(l) = var(l - 1)*(1 - f(l, l)*f(l, l));
      // std::cout <<"l="<<l <<" var.at(l-1) = "<<var.at(l-1) <<"\n";
      if (l  ==  lr) // return;
	return;
	//	std::cout <<"STOP: l=lr\n";
      d = 0;// 0.0d0;
      q = 0;// 0.0d0;
      for(int i=1; i<=l; i++)
	{
	  int k = l - i + 2;
	  d = d + a.at(i-1)*r.at(k-1); // d = d + a(i)*r(k);
	  q = q + coefs.at(l-1).at(i-1)*r.at(k-1); // q = q + f(l, i)*r(k);
	  
	}
      
    }
  return;
}




// ar.yw.default in R package STATS
// ar.out <- ar(x[,i], aic=TRUE)
// v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
void AR::ar(std::vector<double> vars)
{
  int n_used =vars.size();
  double order_max = std::min((double) n_used - 1, floor(10* log10((double) n_used)) );
  if(order_max <1)
    std::cout <<"'order_max' must be >= 1";
  else if(order_max >= n_used)
    std::cout <<"'order_max' must be >= n_used";

  double xm = 0;
  for(int i=0; i<vars.size(); i++)
    xm += vars.at(i);
  xm /= n_used;
    
  std::vector<double> newX;
  for(int i=0; i<vars.size(); i++)
    newX.push_back(vars.at(i)-xm);

  std::vector<double> xacf = acf(newX, order_max);
  //   xacf <- acf(x, type = "covariance", lag.max = order.max, plot = FALSE,
  //              demean = demean)$acf

  //nser = 1 (univariate case)
  if(xacf.at(0) ==0)
    std::cout <<"zero-variance series\n";

  Eureka Z;
  Z.eureka((int) order_max, xacf, xacf);
  std::vector<std::vector<double> > coefs_mat = Z.get_coefs();

  /*
  for(int i=0; i<coefs_mat.size(); i++)
    for(int j=0; j<coefs_mat.at(i).size(); j++)
      {
	std::cout <<"i="<<i<<", j="<<j<<", coefs_mat.at(i).at(j) = "<<coefs_mat.at(i).at(j)<<"\n";
      }
  */
 

  std::vector<double> partialacf; 
  for(int i=0; i<coefs_mat.size(); i++)
    {
      partialacf.push_back(coefs_mat.at(i).at(i));
      // std::cout <<"i="<<i <<" partialacf.at(i) = "<<partialacf.at(i)<<"\n";
    }

 
  std::vector<double> var_pred_vector, xaic;
  var_pred_vector.push_back(xacf.at(0));
  xaic.push_back( (double)n_used*log(xacf.at(0)) + 2);
  std::vector<int> order; order.push_back(0);
  double min_xaic = xaic.at(0);
  for(int i=0; i<Z.get_var().size(); i++)
    {
      double tmp = Z.get_var().at(i);
      var_pred_vector.push_back(tmp);
      xaic.push_back( (double)n_used*log(tmp)+2*(i+1)+2);
      if(xaic.at(i+1) < min_xaic)
	{
	  min_xaic = xaic.at(i+1);
	  order.resize(0); order.push_back(i+1);
	}
      else if(xaic.at(i+1) == min_xaic)
	order.push_back(i+1);
	    
    }

  /*
  for(int i=0; i<xaic.size(); i++)
    {
      std::cout <<"i="<<i 
		<<" xaic.at(i)="<< xaic.at(i) <<"\n";
    }
  std::cout <<"min_xaic = "<<min_xaic << " order.at(0) = "<< order.at(0)<< "\n";
  */
  
  // std::vector<double>::const_iterator it;
  // it = std::min_element(xaic.begin(), xaic.end());
  // double min_xaic=*it;
  for(int i=0; i<xaic.size(); i++)
    xaic.at(i) -= min_xaic;

  /*
  for(int i=0; i<var_pred_vector.size(); i++)
    {
      std::cout <<"i="<<i <<" var_pred_vector.at(i) = "<<var_pred_vector.at(i)
		<<" xaic.at(i)="<< xaic.at(i) <<"\n";
    }
  */
  
   if(order.size()>1)
    {
      std::cout <<"order.size() >1\n";
    }
  if(order.at(0)>=1)
    {
      ar_nan = 0;
      for(unsigned int i=0; i<order.at(0); i++)
	       ar_coef.push_back(coefs_mat.at(order.at(0)-1).at(i));
      // ar_coef.push_back(coefs_mat.at(order.at(0)-1).at(order.size()-1));
    }
  else
    {
      ar_nan = 1;
    }

  
  var_pred = var_pred_vector.at(order.at(0))*(double)n_used/((double) n_used -(order.at(0)+1)) ;
  
 /*    
  std::cout <<"order.at(0) ="<< order.at(0) <<"\n";
  std::cout <<"var_pred = "<<var_pred <<"\n";
 
  for(int i=0; i<ar_coef.size(); i++)
    {
      std::cout <<"i= "<< i<<" ar_coef.at(i) = "<< ar_coef.at(i)<<"\n";
    }
  */
  
  return ;
}

// spectrum0.ar(y)$spec
double spectrum0_ar(std::vector<double> vars)
  {  
  AR ar_out;
  ar_out.ar(vars); //  ar.out <- ar(x[,i], aic=TRUE)
  double var_pred =ar_out.get_var_pred();
  
  double sum_coef=0;
  if(!ar_out.get_ar_nan())
    {
      std::vector<double> coef = ar_out.get_ar_coef();
      for(int i=0; i<coef.size(); i++)
	sum_coef += coef.at(i);
    }

  //  std::cout <<"var_pred = "<<var_pred <<" sum_coef="<< sum_coef <<"\n";
  
  double spec =  var_pred/((1-sum_coef)*(1-sum_coef)) ;//ar.out$var.pred/(1 - sum(ar.out$ar))^2
  
  // std::cout <<"spec = "<<spec <<"\n";
  return spec;
}


double gewekez(std::vector<double> vars1, std::vector<double> vars2)
{

  double mean1 = 0, mean2 = 0;
  for(int i=0; i<vars1.size(); i++)
    mean1 += vars1.at(i);
  for(int i=0; i<vars2.size(); i++)
    mean2 += vars2.at(i);
  mean1 /= vars1.size(); mean2 /= vars2.size();

  double variance1 =0, variance2 =0;
  variance1 = spectrum0_ar(vars1)/vars1.size();
  variance2 = spectrum0_ar(vars2)/vars2.size();

  
  double z = (mean1 - mean2)/sqrt(variance1+variance2);
  return z;
}

double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}


double gewekez(double *vars1, double *vars2, int lenvars1, int lenvars2)
{
  std::vector<double> x1, x2;
  for(int i=0; i<lenvars1; i++)
    x1.push_back(vars1[i]);
  for(int i=0; i<lenvars2; i++)
    x2.push_back(vars2[i]);
  double z = gewekez(x1,x2);
  return z;
}
/* only get here if currentid == HEADNODE */
void 
printgewekez(FILE * outto)
{
  double z,p;
  int h,b;

  if (outto != NULL && poptopologysequence.currentlength > 1000)
  {
    //fprintf(outto,"\n===============================================\n");
    //fprintf(outto,"Geweke Z statistic of Tree Topology Stationarity\n");
    //fprintf(outto,"=================================================\n");
    //fprintf(outto,"\nGeweke Z statistic of Tree Topology Stationarity:\n");
    b = 0;
    h = (poptopologysequence.currentlength-b)/2;
    z = gewekez(&poptopologysequence.disvals[0]+b,&poptopologysequence.disvals[0] + b + h,h,h);
    p = phi(fabs(z));
    fprintf(outto,"\nGeweke Z statistic of Tree Topology Stationarity\n------------------------------------------------\n");
    fprintf(outto,"   All sampled values (%d) split in two parts:   %.4lf   p(z): %.4lf\n",poptopologysequence.currentlength,fabs(z),p);
    b = h = poptopologysequence.currentlength/3;
    z = gewekez(&poptopologysequence.disvals[0]+b,&poptopologysequence.disvals[0] + b + h,h-1,h-1);  // subtract one in case of overrun ?? 
    p = phi(fabs(z));
    fprintf(outto,"   Skip first %d sampled values out of %d:  %.4lf   p(z): %.4lf\n",h,poptopologysequence.currentlength,fabs(z),p);

/*
    fprintf(outto,"\t#skip\t#values\tmean1\tmean2\tz\tp(z)\n");
    for (int i=0;i<=5;i++)
    {
      int b= INTEGERROUND(poptopologysequence.currentlength * (i/10.0));
      int h = (poptopologysequence.currentlength-b)/2;
      double s1 = 0.0;
      double s2 = 0.0;
      for (int j=0;j<h;j++)
      {
        s1 += poptopologysequence.disvals[0+b+j];
        s2 += poptopologysequence.disvals[0+b+h+j];
      }
      double z = gewekez(&poptopologysequence.disvals[0]+b,&poptopologysequence.disvals[0] + b + h,h,h);
      double p = phi(fabs(z));
      
      fprintf(outto,"\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",b,poptopologysequence.currentlength-b, s1/h,s2/h,z,p);
      */
      /*
      s1 = 0.0;
      s2 = 0.0;
      int lag = 10;
      int h2 = (int) h/lag;
      double *v1 = (double *) malloc(sizeof(double)*h2);
      double *v2 = (double *) malloc(sizeof(double)*h2);
      int count =0;
      for (int j=0;j<h& count<h2;j+=lag)
      {
        s1 +=  poptopologysequence.disvals[0+b+j];
	v1[count] =  poptopologysequence.disvals[0+b+j];
        s2 +=  poptopologysequence.disvals[0+b+h+j];
	v2[count] =  poptopologysequence.disvals[0+b+h+j];
	count++;
      }
      z = gewekez(v1,v2,h2,h2);
      p = phi(fabs(z));
      fprintf(outto,"\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n",b,h2, s1/h2,s2/h2,z,p);
      */
      
    //}
     }
}

/*
int main()
{
  std::ifstream infile;
  infile.open("toy.txt");
  std::vector<double> x;
  double tmp;
  while(infile >> tmp)
    {
      x.push_back(tmp);
    }
  infile.close();
 
  infile.open("toy2.txt");
  std::vector<double> y;
  while(infile >> tmp)
    {
      y.push_back(tmp);
    }
  infile.close();
  */
  /*
  for(int i=0; i< x.size(); i++)
    std::cout << "i="<<i<<" x.at(i)="<<x.at(i)<<"\n";
  for(int i=0; i< y.size(); i++)
    std::cout << "i="<<i<<" y.at(i)="<<y.at(i)<<"\n";
  */

/*
  double z = gewekez(x,y);
  std::cout <<"(Using 'vector') Geweke's z = " << z <<"\n";

  int len1 = x.size();
  int len2 = y.size();
  double *var1 = (double *) malloc(sizeof(double)*len1);
  double *var2 = (double *) malloc(sizeof(double)*len1);
  for(int i=0; i<len1; i++)
    var1[i] = x.at(i);
  for(int i=0; i<len2; i++)
    var2[i] = y.at(i);

    
  z = gewekez(var1,var2,len1,len2);
  std::cout <<"(Using 'double *') Geweke's z = " << z <<"\n";
  
  
  return 0;
}
*/
