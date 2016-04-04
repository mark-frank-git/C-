/********************************************************************************
 *                                                                              *
 * This subclass of object implements a class for generating samples from       *
 * different types of distributions.                                            *
 * References:                                                                  *
 * 1. G.E. Johnson, "Constructions of particular random processes," Proc. IEEE  *
 *    vol. 82, no. 2, Feb. 1994, pp. 270-285.
 *                                                                              *
 * File: /User/frank/Objc_Classes/RandomDist/RandomDist.cc                      *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 05/04/92 - Started                                                       *
 *  2. 06/10/92 - Free arrays in here                                           *
 *  3. 07/09/92 - Add gamma, log_gamma, log_normal, and correlation             *
 *                coefficient                                                   *
 *  4. 11/19/92 - Add GAUSS_MARKOV, distributionFunction -> densityFunction     *
 *  5. 12/17/97 - Add setNewSeeds().                                            *
 *  6. 01/17/00 - Went to GNU uniform generator.                                *
 *  7. 08/22/00 - Went to MT19937 random number generator, and changed name.    *
 *  8. 04/17/06 - Added CONSTANT_DIST for testing RadarScatterers.              *
 *  9. 05/03/06 - Added sqrt(2) scaling on Rayleigh samples                     *
 ********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "RandomDist.h"
#include "UniformNumber.h"

#if defined (_WIN32) && defined (NeXT_PDO)
#include <C_Libraries/constants.h>
#include <Polynomials/DoublePoly.h>
#include <Specfuns/specfuns.h>
#include <GNU/Complex.h>
#include <GSLSpecfuns/gsl_sf_result.h>
#include <GSLSpecfuns/gsl_sf_bessel.h>
#define ADD_RICIAN_OUTPUT       1
#define ADD_CHI_SQUARED_OUTPUT  1
#define ADD_NON_CENTRAL_OUTPUT  1
#else
#include "constants.h"                          // For compiling on Sun
#endif

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

#define MIN_RHO         -0.99999999
#define MAX_RHO         0.99999999
#define MIN_VARIANCE    1.e-20
#define SQRT2           1.414213562             /* sqrt(2)      */
#define SQRT3           1.732050808             /* sqrt(3)      */
#define SQRT2PI         2.506628275             /* sqrt(2*PI)   */
#define YES             1
#define NO              0

//
// The following is a polynomial fit to the Rician mean vs SNR plot as shown, e.g., in Fig. 4-6 in Whalen's book
//
#define NUMBER_RICE_MEAN        9               // 1 + order of polynomial for fit
#define MAX_RICE_MEAN           20.             // Beyond 20, the mean == SNR
#define MIN_RICE_MEAN           1.25            // Mean not defined less than 1.25 (SNR < 0)
#define MIN_RICE_ALPHA          0.              // Min value of alpha (SNR > 0)
double rice_mean[NUMBER_RICE_MEAN] = {7.98818e-009, -7.69816e-007, 3.13612e-005, -0.000701982, 0.00938983, -0.0764299,
                                     0.369008, 0.0167069, 1.24687};
// ############################# Private Function ###############################
// initMeanVariance --Initialize means and variances
//
// Input:                   None
//          
// Output:                  None
//
// Notes:
// ############################# Private Function ###############################
void RandomDist::initMeanVariance()
{
  int i;

  delete [] _mean;
  delete [] _variance;
  delete [] _stdDev;
  delete [] _ricianAlpha;
  
  _mean         = new double[_dimensions];
  _variance     = new double[_dimensions];
  _stdDev       = new double[_dimensions];
  _ricianAlpha  = new double[_dimensions];
  for(i=0; i<_dimensions; i++)
  {
    _mean[i]            = 0.;
    _variance[i]        = _stdDev[i] = 1.;
    _ricianAlpha[i]     = 0.;
  }
  return;
}

// ############################# Private Function ###############################
// findRiceAlphaFromMean -- Returns the alpha value from the mean for a Rician
//                          distribution.
//
// Input:       mean:           Input mean value
//          
// Output:                      alpha (SNR)
//
// Notes:
// 1. From the specified mean, we estimate the SNR, alpha of the underlying process
//    and thereby get the parameters for the distribution
// ############################# Private Function ###############################
double RandomDist::findRiceAlphaFromMean(double mean)
{
#if(ADD_RICIAN_OUTPUT)
  int           j;
  double        alpha, closest_root;
  DoublePoly    poly_fit;
  Complex       *poly_roots;
//
// Estimate SNR from input mean:
//
  if(mean <= MIN_RICE_MEAN )
    alpha       = MIN_RICE_ALPHA;
  else if (mean >= MAX_RICE_MEAN)
    alpha       = mean;
  else
  {
    poly_fit.assign(NUMBER_RICE_MEAN-1, rice_mean);
    poly_fit            -= mean;
    poly_roots          = poly_fit.roots();             // SNR = solution of p(x)-mean=0
    closest_root        = real(poly_roots[0]);
    for(j=1; j<NUMBER_RICE_MEAN-1; j++)                 // The closest root should be near mean
    {
      if( ABS((real(poly_roots[j])-mean)) < ABS((closest_root-mean)) )
      {
        closest_root    = real(poly_roots[j]);
      }
    }
    alpha               = closest_root;
  }
  return alpha;
#else
  return 1.+mean/100.;
#endif
}

// ############################# Private Function ###############################
// normalSample -- Returns samples from the Normal distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// 1. Uses Box Muller algorithm as given, for example in "Computer Methods for
//    Mathematical Computations," Forsythe, Malcolm and Moler.
// ############################# Private Function ###############################
double  *RandomDist::normalSample()
{
  int    i;
  double *sample;
  double u1,s,ln_s;
  double v1,v1_sq,v2_sq;
  static double v2, sqrt_lns, temp;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    if(_oddSample)
    {
      s = 2.;
      while(s > 1.)
      {
        u1              = (*_uniformGenerator)();       // returns 0<=x<=1
        v2              = u1+u1-1.;
        v2_sq           = v2*v2;
        u1              = (*_uniformGenerator)();
        v1              = u1 + u1 - 1.;
        v1_sq           = v1*v1;
        s               = v1_sq + v2_sq;
      }
      ln_s              = log(s);
      sqrt_lns          = sqrt(-(ln_s+ln_s)/s);
      _oddSample        = 0;
      sample[i]         = v1*sqrt_lns;
    }
    else
    {
      _oddSample        = 1;
      sample[i]         = v2*sqrt_lns;
    }
    if(i==1)
    {
      sample[i]         = _stdDev[i]*(_rho*temp + sample[i]*sqrt(1.-_rho*_rho));
      sample[i]         += _mean[i];
    }
    else
    {
      temp              = sample[i];
      sample[i]         *= _stdDev[i];
      sample[i]         += _mean[i];
    }
  }
  return        sample;
}

// ############################# Private Function ###############################
// uniformSample -- Returns samples from the Uniform distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// ############################# Private Function ###############################
double  *RandomDist::uniformSample()
{
  int    i;
  double *sample, range;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    sample[i]   = (*_uniformGenerator)();               // sample in [0,1]
    sample[i]   -= 0.5;
    range       = SQRT3*_stdDev[i];
    sample[i]   *= range + range;
    sample[i]   += _mean[i];
  }

  return        sample;
}

// ############################# Private Function ###############################
// exponentialSample -- Returns samples from the Exponential distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// ############################# Private Function ###############################
double  *RandomDist::exponentialSample()
{
  int    i;
  double *sample;
  double u1,alpha;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    alpha       = _mean[i] - _stdDev[i];
    u1          = (*_uniformGenerator)();
    if(u1>0.)
      sample[i] = -_stdDev[i]*log(u1);
    else
      sample[i] = 0.;
    sample[i]   += alpha;
  }
  return        sample;
}

// ############################# Private Function ###############################
// erlangSample -- Returns samples from the Erlang distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// ############################# Private Function ###############################
double  *RandomDist::erlangSample()
{
  int           i, j, k;
  double        *sample;
  double        a, tr;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    a           = _mean[i]/_variance[i];
    k           = ROUND(_mean[i]*a);
    tr          = 1.;
    for(j=0; j<k; j++)
      tr        *= (*_uniformGenerator)();
    if(tr>0.)
      sample[i] = -log(tr)/a;
    else
      sample[i] = 0.;
  }
  return        sample;
}

// ############################# Private Function ###############################
// lognormalSample -- Returns samples from the Lognormal distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// ############################# Private Function ###############################
double  *RandomDist::lognormalSample()
{
  int    i, j;
  double *sample;
  double std_devx, std_devy, term, mean_y, sum;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    std_devx            = _stdDev[i];
    term                = log(_variance[i]/_mean[i]/_mean[i] + 1.);
    std_devy            = sqrt(term);
    if(_mean[i]>0.)
      mean_y            = log(_mean[i]) - 0.5*term;
    else
      mean_y            = 0.;
    sum                 = -6.0;
    for(j=0; j<12; j++)
      sum               += (*_uniformGenerator)();
    sample[i]           = exp(mean_y + std_devy*sum);
  }

  return        sample;
}

// ############################# Private Function ###############################
// gaussMarkovSample -- Returns samples from the Gauss-Markov distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// ############################# Private Function ###############################
double  *RandomDist::gaussMarkovSample()
{
  int    i;
  double *sample;
  double u1,s,ln_s;
  double v1,v1_sq,v2_sq;
  static double v2, sqrt_lns;

  sample                = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    if(_oddSample)
    {
      s = 2.;
      while(s > 1.)
      {
        u1              = (*_uniformGenerator)();
        v2              = u1+u1-1.;
        v2_sq           = v2*v2;
        u1              = (*_uniformGenerator)();
        v1              = u1 + u1 - 1.;
        v1_sq           = v1*v1;
        s               = v1_sq + v2_sq;
      }
      ln_s              = log(s);
      sqrt_lns          = sqrt(-(ln_s+ln_s)/s);
      _oddSample        = 0;
      sample[i]         = _rho*_markovOld + sqrt(1.-_rho*_rho)*v1*sqrt_lns;
      _markovOld        = sample[i];
    }
    else
    {
      _oddSample        = 1;
      sample[i]         = _rho*_markovOld + sqrt(1.-_rho*_rho)*v2*sqrt_lns;
      _markovOld        = sample[i];
    }
  }

  return        sample;
}

// ############################# Private Function ###############################
// ricianSample -- Returns samples from the Rician distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// 1. The mean and variance for the Rician distribution are linked, therefore, we
//    use only the mean in generating the sample.
// 2. From the specified mean, we estimate the SNR, alpha of the underlying process
//    and thereby get the parameters for the distribution
// ############################# Private Function ###############################
double  *RandomDist::ricianSample()
{
  int           i;
  double        *sample, n1, n2, u1, u2, alpha;

  sample        = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    alpha       = _ricianAlpha[i];
//
// Use the formula in Table 1 of Johnson's paper.
//
    u1          = -1.;
    while (u1 <= 0.)
      u1        = (*_uniformGenerator)();
    u2          = (*_uniformGenerator)();
    n1          = sqrt(-2.*log(u1))*cos(TWOPI*u2);
    u1          = -1.;
    while (u1 <= 0.)
      u1        = (*_uniformGenerator)();
    u2          = (*_uniformGenerator)();
    n2          = sqrt(-2.*log(u1))*cos(TWOPI*u2);
    sample[i]   = sqrt(n1*n1 + (n2+alpha)*(n2+alpha));
  }

  return        sample;
}

// ############################# Private Function ###############################
// rayleighSample -- Returns samples from the Rayleigh distribution
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// 1. The mean and variance for the Rayleigh are linked, therefore, we use only
//    the mean in generating the sample.
// ############################# Private Function ###############################
double  *RandomDist::rayleighSample()
{
  int           i;
  double        *sample, two_sigma_sq, u1;

  sample                = _doubleOutputs;
  double scale          = 4/PI;
  for(i=0; i<_dimensions; i++)
  {
    two_sigma_sq        = scale*_mean[i]*_mean[i];              // Whalen 4-48
    u1                  = -1.;
    while(u1 <= 0.)
      u1                = (*_uniformGenerator)();
    sample[i]           = sqrt(-two_sigma_sq*log(u1));
  }
  return        sample;
}

// ############################# Private Function ###############################
// constantSample -- Returns the mean values (used for testing RadarScatterers)
//
// Input:                       None
//          
// Output:                      array of samples of size, _dimensions
//
// Notes:
// 1. .
// ############################# Private Function ###############################
double  *RandomDist::constantSample()
{
  int           i;
  double        *sample;

  sample                = _doubleOutputs;
  for(i=0; i<_dimensions; i++)
  {
    sample[i]           = _mean[i];
  }
  return        sample;
}

// ############################# Private Function ###############################
// normalDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// ############################# Private Function ###############################
double  RandomDist::normalDensityAt(double *x)
{
  int    i;
  double p_of_x, temp;
  
  temp          = (1-_rho*_rho);
  p_of_x        = 1.;
  if(_dimensions == 2)
  {
    p_of_x = exp(_rho*(x[0]-_mean[0])*(x[1]-_mean[1])/_stdDev[0]/
                 _stdDev[1]/temp);
  }
  temp *= 2.;
  for(i=0; i<_dimensions; i++)
  {
    p_of_x      *= exp(-(x[i]-_mean[i])*(x[i]-_mean[i])/temp/_variance[i]);
    p_of_x      /= SQRT2PI*_stdDev[i];
  }
  return p_of_x;
}

// ############################# Private Function ###############################
// uniformDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// ############################# Private Function ###############################
double  RandomDist::uniformDensityAt(double *x)
{
  int    i;
  int    in_range;
  double p_of_x, range, volume;
  
  in_range = YES;
  volume   = 1.;
  for(i=0; i<_dimensions; i++)
  {
    range   = SQRT3*_stdDev[i];
    volume *= range;
    if( (x[i]<(_mean[i]+range)) && (x[i]>=(_mean[i]-range)) )
          ;
    else
    {
      in_range = NO;
      break;
    }
  }
  if(in_range)
    p_of_x = 1./volume;
  else
    p_of_x = 0.;
  return p_of_x;
}

// ############################# Private Function ###############################
// exponentialDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// ############################# Private Function ###############################
double  RandomDist::exponentialDensityAt(double *x)
{
  int           in_range, i;
  double        alpha, p_of_x;
  in_range = YES;
  for(i=0; i<_dimensions; i++)
  {
    alpha   = _mean[i] - _stdDev[i];
    if( x[i]>=alpha )
      ;
    else
    {
      in_range = NO;
      break;
    }
  }
  if(in_range)
  {
    p_of_x = 1.;
    for(i=0; i<_dimensions; i++)
    {
      p_of_x    *= exp(-(x[i]-_mean[i])/_stdDev[i]);
      p_of_x    /= _stdDev[i];
    }
  }
  else
    p_of_x = 0.;
  return p_of_x;
}

// ############################# Private Function ###############################
// erlangDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// ############################# Private Function ###############################
double  RandomDist::erlangDensityAt(double *x)
{
  int           in_range, i, j, n;
  double        p_of_x, c, sum, temp;

  in_range      = YES;
  p_of_x        = 1.;
  for(i=0; i<_dimensions; i++)
  {
    if(x[i] <= 0.)
    {
      in_range  = NO;
      break;
    }
    c           = _mean[i]/_variance[i];
    n           = ROUND(_mean[i]*c);
    sum         = 0.;
    for(j=1; j<n; j++)
      sum       += log((double)j);
    temp        = n*log(c) + (n-1)*log(x[i]) - c*x[i] - sum;
    p_of_x      *= exp(temp);
  }
  if(!in_range)
    p_of_x      = 0.;
  return p_of_x;
}

// ############################# Private Function ###############################
// lognormalDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// ############################# Private Function ###############################
double  RandomDist::lognormalDensityAt(double *x)
{
  int   in_range, i;
  double        p_of_x, sigma_sq, mu, temp;
  
  in_range      = YES;
  p_of_x        = 1.;
  for(i=0; i<_dimensions; i++)
  {
    if(x[i] <= 0.)
    {
      in_range  = NO;
      break;
    }
    if(_mean[i]>0.)
    {
      sigma_sq  = log(_variance[i] + _mean[i]*_mean[i]) - 2.*log(_mean[i]);
      mu        = log(_mean[i]) - sigma_sq/2;
    }
    else
    mu          = sigma_sq              = 0.;
    temp        = log(x[i]) - mu;
    p_of_x      *= exp(-temp*temp/2./sigma_sq)/SQRT2PI/sqrt(sigma_sq)/x[i];
  }
  if(!in_range)
    p_of_x = 0.;
  return p_of_x;
}

// ############################# Private Function ###############################
// gaussMarkovDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. Non-functional at present time
// ############################# Private Function ###############################
double  RandomDist::gaussMarkovDensityAt(double *x)
{
  return        x[0];
}

#define MIN_EXP_ARG     -15.
#define BESSEL_I0       0
#define NO_EXP_BESSEL   0
#define EXP_BESSEL      1
// ############################# Private Function ###############################
// ricianDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. One dimensional only
// 2. Plots (4-53) in Whalen, p(v) = v exp[-(v**2+alpha**2)/2]*Io(alpha*v)
// ############################# Private Function ###############################
double  RandomDist::ricianDensityAt(double *v)
{
  double        bessel_arg, exp_arg;
  double        pdf;
  
  if(v[0] <= 0.)
    return 0.;

  bessel_arg    = v[0]*_ricianAlpha[0];
  exp_arg       = (v[0]*v[0] + _ricianAlpha[0]*_ricianAlpha[0])/2.;
  exp_arg       = -exp_arg;
  pdf           = 0.;
#if(ADD_RICIAN_OUTPUT)
  if(exp_arg < MIN_EXP_ARG)
  {
    exp_arg     += bessel_arg;
    pdf         = v[0]*exp(exp_arg)*modbes(bessel_arg, BESSEL_I0, EXP_BESSEL);
  }
  else
  {
    pdf         = v[0]*exp(exp_arg)*modbes(bessel_arg, BESSEL_I0, NO_EXP_BESSEL);
  }
#endif
  return        pdf;
}

// ############################# Private Function ###############################
// rayleighDensityAt -- Returns the probability density function at the given input
//                    value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. One dimensional only
// 2. Plots (4-47) in Whalen, p(z) = z/sigma**2 exp[-z**2/2/sigma**2]
// 3. The mean and variance are linked, we use the mean to set the variance.
// ############################# Private Function ###############################
double  RandomDist::rayleighDensityAt(double *x)
{
  double        sigma_sqd, pdf;
  sigma_sqd     = 2*_mean[0]*_mean[0]/PI;       // E{z} = sigma*sqrt(PI/2.)

  pdf           = x[0]*exp(-x[0]*x[0]/2./sigma_sqd)/sigma_sqd;
  return        pdf;
}

// ############################# Private Function ###############################
// chiSquaredDensityAt -- Returns the probability density function at the given input
//                        value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. One dimensional only
// 2. Plots (4-65) in Whalen.
// 3. The mean and variance are linked, we use the mean to set the variance.
// ############################# Private Function ###############################
double  RandomDist:: chiSquaredDensityAt (double *y)
{
//
// Check if degrees is even:
//
  double p_of_x = 0.;
#if(ADD_CHI_SQUARED_OUTPUT)
//
// degrees = mean
//
  int degrees           = ROUND(_mean[0]/2.);
  double x              = y[0];
  if( (degrees%2) == 0)
  {
    int degrees_half    = degrees/2;
    p_of_x              = ipow(x/2., degrees_half-1)*exp(-x/2.);
    p_of_x              /= 2.*gamma_n(degrees_half);
  }
  else
  {
    double n_over_2     = degrees/2.;
    int degrees_half    = degrees/2;
    p_of_x              = pow(x/2., n_over_2-1.)*exp(-x/2.);
    p_of_x              /=2.*gamma_n12(degrees_half);
  }
#endif
  return p_of_x;
}

#define MAX_EXP         15.
#define LINEAR_COEFF    0.0281386               // Coefficient for checking for using normal approx
#define LINEAR_OFFSET   -4.5                    // offset coefficient for normal approx
// ############################# Private Function ###############################
// nonCentralDensityAt -- Returns the probability density function at the given input
//                        value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. One dimensional only
// 2. Plots (4-77) in Whalen.
// 3. The mean and variance are related to the noncentral parameter and degrees
//    of freedom as, E{q} = lambda+N, var{q} = 4*lambda+N
// ############################# Private Function ###############################
double  RandomDist::nonCentralDensityAt (double *y)
{
  double        p_of_x  = 0.;
#if(ADD_NON_CENTRAL_OUTPUT)
  gsl_sf_result result;
//
// Find lambda and N from mean and variance
//
  double lambda = (_variance[0] - 2.*_mean[0])/2.;
  int degrees   = ROUND((_mean[0] - lambda)/2.);
  double x      = y[0];
//
  if(x <= 0.)
    p_of_x      = 0.;
  else
  {
    double arg  = LINEAR_COEFF*degrees + LINEAR_OFFSET;
    if(sqrt(lambda) < arg)
    {
          return normalDensityAt(y);                    // Out of range of Bessel, use normal approx
    }
    arg         = sqrt(x*lambda);
    gsl_sf_bessel_Inu_scaled_impl((degrees-1.), sqrt(x*lambda), &result);
    double temp = result.val;
    temp        *= exp( -(lambda+x-2.*arg)/2.);
    temp *= pow(x/lambda, (degrees-1.)/2.)/2.;
    p_of_x = temp;
  }
#endif
  return p_of_x;
}
// ############################# Private Function ###############################
// nonCentralDensityAt -- Returns the probability density function at the given input
//                        value.
//
// Input:       x:              input function value
//          
// Output:                      PDF
//
// Notes:
// 1. One dimensional only
// 2. Plots (4-77) in Whalen.
// 3. The mean and variance are related to the noncentral parameter and degrees
//    of freedom as, E{q} = lambda+N, var{q} = 4*lambda+N
// ############################# Private Function ###############################
double  RandomDist::oldnonCentralDensityAt (double *y)
{
  double        p_of_x = 0.;
#if(ADD_NON_CENTRAL_OUTPUT)
//
// Find lambda and N from mean and variance
//
  double lambda = (_variance[0] - 2.*_mean[0])/2.;
  int degrees   = ROUND((_mean[0] - lambda)/2.);
  double x      = y[0];
//
  double temp, arg;
  int   error;
  if(x <= 0.)
    p_of_x      = 0.;
  else
  {
    int islct   = 0;
    switch(degrees)
    {
      
      case 4:
        islct   = 1;
      case 2:
        arg     = sqrt(x*lambda);
        if(arg > MAX_EXP)
        {
          temp = modbes(arg, islct, 1);
          temp *= exp( -(lambda+x-2.*arg)/2.);
        }
        else
        {
          temp = modbes(arg, islct, 0);
          temp *= exp( -(lambda+x)/2.);
        }
        break;
      default:
        arg     = sqrt(x*lambda);
        if(arg > MAX_EXP)
        {
          temp = ModBesselILargeExp((double)(degrees/2.-1.), sqrt(x*lambda));
          temp *= exp( -(lambda+x-2.*arg)/2.);
        }
        else
        {
          temp = ModBesselI((double)(degrees/2.-1.), sqrt(x*lambda), &error);
          temp *= exp( -(lambda+x)/2.);
        }
        break;
    }
    temp *= pow(x/lambda, (degrees-2.)/4.)/2.;
    p_of_x = temp;
  }
#endif
  return p_of_x;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RandomDist class.
//
// Input:           type:           window type
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
RandomDist::RandomDist ()
{
// 
// Initialize instance variables:
//
  _mean                 = _variance = _stdDev = NULL;
  _ricianAlpha          = NULL;
  _distributionType     = NORMAL;
  _probability          = 1.;
  _rho                  = 0.;
  _doubleOutputs        = NULL;
  setDimensions(1);
  initMeanVariance();
  
  _uniformGenerator             = new UniformNumber();
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the RandomDist class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
RandomDist::~RandomDist ()
{
  delete _uniformGenerator;
  delete [] _doubleOutputs;
  delete [] _mean;
  delete [] _variance;
  delete [] _stdDev;
  delete [] _ricianAlpha;
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the noise generator to a defined state.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void RandomDist::initGenerator()
{
  _oddSample            = 1;
  _markovOld            = 0.;
  if(_uniformGenerator != NULL)
    _uniformGenerator->reset();
  return;
}

// ############################# Public Function ###############################
// setMean -- sets a new set of means
//
// Input:       newMean:        Array of means
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setMean( double *newMean)
{
  int i;
  for(i=0; i<_dimensions; i++)
  {
    _mean[i]            = newMean[i];
    _ricianAlpha[i]     = findRiceAlphaFromMean(newMean[i]);
  }
}

// ############################# Public Function ###############################
// setVariance -- sets a new set of variances
//
// Input:       newVariance:    Array of variances
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setVariance( double *newVariance)
{
  int i;
  for(i=0; i<_dimensions; i++)
  {
    _variance[i] = MAX(MIN_VARIANCE, newVariance[i]);
    _stdDev[i]   = sqrt(_variance[i]);
  }
}

// ############################# Public Function ###############################
// setCorrelationCoefficient -- sets a correlation coefficient for 2D distributions
//
// Input:       newRho:         New correlation coefficient
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setCorrelationCoefficient(double newRho)
{
  _rho = MAX(MIN_RHO,newRho);
  _rho = MIN(MAX_RHO,_rho);
}

// ############################# Public Function ###############################
// setProbability -- sets a new probability of occurrence
//
// Input:       newRandomDist:          Probability of occurring
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setProbability(double probability)
{
  _probability = probability;
}

// ############################# Public Function ###############################
// setDistributionType -- sets a new type of probability distribution
//
// Input:       newDistribution:        New probability distribution
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setDistributionType(int newDistribution)
{
  _distributionType = newDistribution;
}

// ############################# Public Function ###############################
// setDimensions -- sets a new number of dimensions for sample outputs
//
// Input:       newDimensions:          New number of dimensions
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setDimensions(int newDimensions)
{
  delete [] _doubleOutputs;
  _dimensions           = newDimensions;
  _doubleOutputs        = new double[_dimensions];
  initMeanVariance();
}

// ############################# Public Function ###############################
// setSeed -- sets a new seed for the random number generator
//
// Input:       seed:                   New seed
//          
// Output:                              None
//
// Notes:
// ############################# Public Function ###############################
void RandomDist::setSeed(long seed)
{
  _uniformGenerator->setSeed(seed);
  return;
}


// ############################# Public Function ###############################
// newSample -- Returns a set of samples from our distribution
//
// Input:                               None
//          
// Output:                              Array of samples, array size = dimensions
//
// Notes:
// ############################# Public Function ###############################
double * RandomDist::newSample()
{
  switch(_distributionType)
  {
    case UNIFORM:
      return uniformSample();
    case NORMAL:
    default:
      return normalSample();
    case EXPONENTIAL_DIST:
      return exponentialSample();
    case ERLANG_DIST:
      return erlangSample();
    case LOG_NORMAL:
      return lognormalSample();
    case GAUSS_MARKOV:
      return gaussMarkovSample();
    case RICIAN_DIST:
      return ricianSample();
    case RAYLEIGH_DIST:
      return rayleighSample();
    case CONSTANT_DIST:
      return constantSample();
  }
}

// ############################# Public Function ###############################
// densityFunctionAtX -- Returns the value of the density function at the input
//                       point.
//
// Input:       x:      point to evaluate PDF x[_dimensions].
//          
// Output:              PDF
//
// Notes:
// ############################# Public Function ###############################
double RandomDist::densityFunctionAtX( double * x)
{
  switch(_distributionType)
  {
    case UNIFORM:
      return uniformDensityAt(x);
    case NORMAL:
    default:
      return normalDensityAt(x);
    case EXPONENTIAL_DIST:
      return exponentialDensityAt(x);
      break;
    case ERLANG_DIST:       /* see Papoulis p. 77 */
      return erlangDensityAt(x);
      break;
    case LOG_NORMAL:        /* see Whalen p. 23 */
      return lognormalDensityAt(x);
      break;
    case GAUSS_MARKOV:
      return gaussMarkovDensityAt(x);
      break;
    case RICIAN_DIST:
      return ricianDensityAt(x);
      break;
    case RAYLEIGH_DIST:
      return rayleighDensityAt(x);
      break;
    case CHI_SQUARED_DIST:
      return chiSquaredDensityAt(x);
      break;
    case NON_CENTRAL_DIST:
      return nonCentralDensityAt(x);
      break;
  }
}

