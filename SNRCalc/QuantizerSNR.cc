/********************************************************************************
 *                                                                              *
 * This class calculates the output SNR from an n level quantizer.              *
 *                                                                              *
 * File: /User/frank/C++/SNRCalc/QuantizerSNR.h                                 *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 10/20/00 - Started.                                                      *
 *  2. 12/05/01 - Modify noise density when setting SNR.                        *
 ********************************************************************************/
#include "QuantizerSNR.h"
#include "NoiseCorrelation.h"
#include <math.h>
#include <Specfuns/specfuns.h>
#include <stdio.h>
#import <CumulativeDistributions/multi_variate_normal.h>

#define N_DIMENSIONS    2
#define USE_CUM         0
// ############################ Private Function ###################################
// lfnh - Calculates and returns the cumulative binomial distribution function.
// Input:       h:      first lower limit
//              k:      second lower limit
//              rho:    correlation coefficient
//
// Output:              cumulative binomial distribution function
//
// Notes:
// ############################ Private Function ###################################
double QuantizerSNR::lfnh(double h, double k, double rho)
{
  double        lfn_rtn;
#if USE_CUM
  int           error;
  double        d[N_DIMENSIONS], rhos[N_DIMENSIONS];
  d[0]          = h;
  d[1]          = k;
  rhos[0]       = rho;
  rhos[1]       = rho;
  lfn_rtn       = multiNormalCDF(d, rhos, N_DIMENSIONS, &error);
#else
  lfn_rtn       = lfn(h, k, rho);
#endif
  return lfn_rtn;
}

// ############################ Private Function ###################################
// mean2Level - Calculates and returns the mean value out of a 1 bit quantizer.
// Input:               None
// Output:              Mean value in volts
//
// Notes:
// ############################ Private Function ###################################
double QuantizerSNR::mean2Level()
{
  double        mean_value, snr_in, sqrt_snr;
//
  snr_in        = pow(10., _inputSNR/10.);      // SNR as a ratio
  sqrt_snr      = sqrt(snr_in);
  mean_value    = 1. - 2.*qx(sqrt_snr);
  return mean_value;
}

// ############################ Private Function ###################################
// varianceIJ2Level - Calculates and returns the variance term out of the quantizer when
//            I not equal to J, for the 1 bit quantizer
// Input:       tau:    Autocorrelation lag in seconds
// Output:              variance value in watts
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::varianceIJ2Level(double tau)
{
  double        rn, variance_term, snr_in, sqrt_snr, sigma_n2, mean_value;
//
  snr_in        = pow(10., _inputSNR/10.);      // SNR as a ratio
  sqrt_snr      = sqrt(snr_in);
  sigma_n2      = _noiseAuto->outputNoisePower();
  rn            = _noiseAuto->autocorrelationAt(tau);
//
// Use arcsin approximation:
//
  variance_term = 2.*asin(rn/sigma_n2)/PI;
  mean_value    = mean2Level();
  variance_term += mean_value*mean_value;
  return variance_term;
}

// ############################ Public Function ###################################
// varianceII2Level - Calculates and returns the variance term out of the quantizer when
//            I equals J, for the 1 bit quantizer
// Input:               None
// Output:              variance value in watts
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::varianceII2Level()
{
  double        mean_value, variance_term;
  
  mean_value    = mean2Level();
  variance_term = 1. + mean_value*mean_value;
  return variance_term;
}

// ############################# Class Constructor #################################
// QuantizerSNR -- Constructor for the QuantizerSNR class
// Input:       levels:         Number of quantizer levels
//              step:           type of correlation operation
//              samples:        sample/chip for input signal
//
// Output:              None
// ############################# Class Constructor #################################
QuantizerSNR::QuantizerSNR(int levels, double stepSize, double inputSNR)
{
// 
// Initialize instance variables:
//
  _noiseAuto    = new NoiseCorrelation();
  setNumberLevels(levels);
  setQuantizerStepSize(stepSize);
  setInputSNR(inputSNR);
  return;
}

// ############################# Class Destructor ###############################
// QuantizerSNR -- Destructor for the QuantizerSNR class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
QuantizerSNR::~QuantizerSNR()
{
//
// Delete space:
//
  delete _noiseAuto;
  return;
}

// ############################ Public Function ###################################
// setInputSNR - Sets a new input SNR (to quantizer).
// Input:       snr:    SNR in dB
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setInputSNR(double snr)
{
  double        noise_density, noise_power, snr_desired, snr_current;
  _inputSNR     = snr;
//
// Modify noise density to achieve given SNR assuming signal power = 1
//
  noise_power   = _noiseAuto->outputNoisePower();
  if(noise_power > 0.)
    snr_current = 1./noise_power;
  else
    snr_current = 1.;
  noise_density = _noiseAuto->noiseDensity();
//
// Scale noise density to achieve SNR:
//
  snr_desired   = pow(10., _inputSNR/10.);      // SNR as a ratio
  noise_density *= snr_current/snr_desired;
  _noiseAuto->setNoiseDensity(noise_density);
  return;
}

// ############################ Public Function ###################################
// setFilterType - Sets a new type of IF filter.
// Input:       type:   Ideal, Butterworth
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setFilterType(int type)
{
  _noiseAuto->setFilterType(type);
  setInputSNR(_inputSNR);               // Affected by filter
  return;
}

// ############################ Public Function ###################################
// setFilterOrder - Sets the order of the Butterworth filter
// Input:       order:  Order of Butterworth filter (if selected)
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setFilterOrder(int order)
{
  _noiseAuto->setFilterOrder(order);
  setInputSNR(_inputSNR);               // Affected by filter
  return;
}

// ############################ Public Function ###################################
// setFilterCutoff - Sets a new filter cutoff.
// Input:       cutoff: One sided bandwidth in Hz
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setFilterCutoff(double cutoff)
{
  _noiseAuto->setFilterCutoff(cutoff);
  setInputSNR(_inputSNR);               // Affected by filter
  return;
}

// ############################ Public Function ###################################
// setNoiseDensity - Sets the noise density in W/Hz.
// Input:       density:        Noise density in W/Hz
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setNoiseDensity(double density)
{
  _noiseAuto->setNoiseDensity(density);
  return;
}

// ############################ Public Function ###################################
// setCenterFrequency - Sets the center frequency of the IF filter in Hz.
// Input:       freq:   Center frequency in Hz
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void QuantizerSNR::setCenterFrequency(double freq)
{
  _noiseAuto->setCenterFrequency(freq);
  return;
}

// ############################ Public Function ###################################
// filterType - Returns the IF filter type.
// Input:               None
// Output:              Ideal, Butterworth
//
// Notes:
// ############################ Public Function ###################################
int QuantizerSNR::filterType()
{
  return _noiseAuto->filterType();
}

// ############################ Public Function ###################################
// filterOrder - Returns the order of the IF filter.
// Input:               None
// Output:              Order of Butterworth filter
//
// Notes:
// ############################ Public Function ###################################
int QuantizerSNR::filterOrder()
{
  return _noiseAuto->filterOrder();
}

// ############################ Public Function ###################################
// filterCutoff - Returns the 1 sided BW of IF filter in Hz.
// Input:               None
// Output:              1 sided BW in Hz
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::filterCutoff()
{
  return _noiseAuto->filterCutoff();
}

// ############################ Public Function ###################################
// noiseDensity - Returns the noise density in W/Hz.
// Input:               None
// Output:              Noise density in W/Hz
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::noiseDensity()
{
  return _noiseAuto->noiseDensity();
}

// ############################ Public Function ###################################
// noiseSigma - Returns the standard deviation of the noise.
// Input:               None
// Output:              Noise standard deviation
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::noiseSigma()
{
  double        sigma_n;
  sigma_n       = _noiseAuto->outputNoisePower();
  if(sigma_n>0.)
    sigma_n     = sqrt(sigma_n);
  return sigma_n;
}

// ############################ Public Function ###################################
// centerFrequency - Returns the center frequency of the IF filter in Hz.
// Input:               None
// Output:              center frequency in Hz
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::centerFrequency()
{
  return _noiseAuto->centerFrequency();
}

// ############################ Public Function ###################################
// meanValue - Calculates and returns the mean value out of quantizer.
// Input:               None
// Output:              Mean value in volts
//
// Notes:
// 1. This function only works for an even number of quantizer levels
// ############################ Public Function ###################################
double QuantizerSNR::meanValue()
{
  int           j, m;
  double        mean_value, snr_in, sqrt_snr;
  double        sigma_n, j1_delta, j_delta;

  if(_quantizerLevels==2)
    return mean2Level();                        // Special case
  m             = _quantizerLevels/2 -1;
  snr_in        = pow(10., _inputSNR/10.);      // SNR as a ratio
  sqrt_snr      = sqrt(snr_in);
  mean_value    = 0.;
  sigma_n       = _noiseAuto->outputNoisePower();
  if(sigma_n>0.)
    sigma_n     = sqrt(sigma_n);
  for(j=1; j<=m; j++)
  {
    j1_delta            = (j-1)*_quantizerStepSize/SQRT2/sigma_n;
    j_delta             = j*_quantizerStepSize/SQRT2/sigma_n;
    mean_value          += j*erf(j_delta-sqrt_snr/SQRT2)/2.;
    mean_value          += j*erf(j1_delta+sqrt_snr/SQRT2)/2.;
    mean_value          -= j*erf(j_delta+sqrt_snr/SQRT2)/2.;
    if(j==1)
      mean_value        += erf(sqrt_snr/SQRT2)/2.;
    else
      mean_value        -= j*erf(j1_delta-sqrt_snr/SQRT2)/2.;
  }
  mean_value    += (m+1)*(qx(m*_quantizerStepSize/sigma_n-sqrt_snr)- px(-m*_quantizerStepSize/sigma_n-sqrt_snr));
  return mean_value;
}

// ############################ Public Function ###################################
// varianceIJ - Calculates and returns the variance term out of the quantizer when
//            I not equal to J.
// Input:       tau:    Autocorrelation lag in seconds
// Output:              variance value in watts
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::varianceIJ(double tau)
{
  int           cap_m, n, m;
  double        rho, n_a, m_a, n_1a, m_1a, n_m1a, m_m1a;
  double        mn_pa, mn_1pa, mn_m1pa, cap_ma, cap_mpa;
  double        variance_term, delta, snr_in, sqrt_snr, sigma_n, sum;
//
  if(_quantizerLevels == 2)
    return varianceIJ2Level(tau);
  cap_m         = _quantizerLevels/2 -1;
  snr_in        = pow(10., _inputSNR/10.);      // SNR as a ratio
  sqrt_snr      = sqrt(snr_in);
  sigma_n       = _noiseAuto->outputNoisePower();
  rho           = 1.;                           // default
  if(sigma_n>0.)
  {
    rho         = _noiseAuto->autocorrelationAt(tau)/sigma_n;
    sigma_n     = sqrt(sigma_n);
  }
  delta         = _quantizerStepSize;
//
  variance_term = 0.;
  cap_ma        = cap_m*delta/sigma_n-sqrt_snr;
  cap_mpa       = cap_m*delta/sigma_n+sqrt_snr;
//
// First double sum:
//
  sum           = 0.;
  sum           = 0.;
  for(n=-cap_m; n<=-1; n++)
  {
    n_1a        = (n+1)*delta/sigma_n-sqrt_snr;
    n_a         = n*delta/sigma_n-sqrt_snr;
    for(m=-cap_m; m<=-1; m++)
    {
      m_a       = m*delta/sigma_n-sqrt_snr;
      m_1a      = (m+1)*delta/sigma_n-sqrt_snr;
      sum       += (lfnh(n_a,m_a,rho)+lfnh(n_1a,m_1a,rho))*n*m;
      sum       -= (lfnh(n_1a,m_a,rho)+lfnh(n_a,m_1a,rho))*n*m;
    }
  }
  variance_term = sum;
//
// Second double sum:
//
  sum           = 0.;
  for(n=1; n<=cap_m; n++)
  {
    n_a         = n*delta/sigma_n-sqrt_snr;
    n_m1a       = (n-1)*delta/sigma_n-sqrt_snr;
    for(m=1; m<=cap_m; m++)
    {
      m_m1a     = (m-1)*delta/sigma_n - sqrt_snr;
      m_a       = m*delta/sigma_n - sqrt_snr;
      sum       += (lfnh(n_m1a,m_m1a,rho)+lfnh(n_a,m_a,rho))*n*m;
      sum       -= (lfnh(n_a,m_m1a,rho)+lfnh(n_m1a,m_a,rho))*n*m;
    }
  }
  variance_term += sum;
//
// Third double sum:
//
  sum           = 0.;
  for(n=-cap_m; n<=-1; n++)
  {
    n_a         = n*delta/sigma_n-sqrt_snr;
    n_1a        = (n+1)*delta/sigma_n-sqrt_snr;
    for(m=1; m<=cap_m; m++)
    {
      m_m1a     = (m-1)*delta/sigma_n - sqrt_snr;
      m_a       = m*delta/sigma_n - sqrt_snr;
      sum       += (lfnh(n_a,m_m1a,rho)+lfnh(n_1a,m_a,rho))*n*m;
      sum       -= (lfnh(n_1a,m_m1a,rho)+lfnh(n_a,m_a,rho))*n*m;
    }
  }
  variance_term += sum + sum;
//
// First single sum:
//
  sum           = 0.;
  for(n=-cap_m; n<=-1; n++)
  {
    n_a         = n*delta/sigma_n-sqrt_snr;
    n_1a        = (n+1)*delta/sigma_n-sqrt_snr;
    mn_1pa      = -(n+1)*delta/sigma_n+sqrt_snr;
    mn_pa       = -n*delta/sigma_n+sqrt_snr;
    sum         += (lfnh(n_a, cap_ma, rho) - lfnh(n_1a, cap_ma, rho))*n*(cap_m+1);
    sum         -= (lfnh(mn_1pa, cap_mpa, rho) - lfnh(mn_pa, cap_mpa, rho))*n*(cap_m+1);
  }
  variance_term += sum + sum;
//
// Second single sum:
//
  sum           = 0.;
  for(n=1; n<=cap_m; n++)
  {
    n_m1a       = (n-1)*delta/sigma_n-sqrt_snr;
    n_a         = n*delta/sigma_n-sqrt_snr;
    mn_pa       = -n*delta/sigma_n+sqrt_snr;
    mn_m1pa     = -(n-1)*delta/sigma_n+sqrt_snr;
    sum         += (lfnh(n_m1a, cap_ma, rho) - lfnh(n_a, cap_ma, rho))*n*(cap_m+1);
    sum         -= (lfnh(mn_pa, cap_mpa, rho) - lfnh(mn_m1pa, cap_mpa, rho))*n*(cap_m+1);
  }
  variance_term += sum + sum;
//
// Additional terms:
//
  variance_term += (cap_m+1)*(cap_m+1)*(lfnh(cap_mpa,cap_mpa,rho) - 2.*lfnh(cap_mpa,cap_ma,-rho)
                                       +lfnh(cap_ma,cap_ma,rho));
  return variance_term;
}

// ############################ Public Function ###################################
// varianceII - Calculates and returns the variance term out of the quantizer when
//            I equals J.
// Input:               None
// Output:              variance value in watts
//
// Notes:
// ############################ Public Function ###################################
double QuantizerSNR::varianceII()
{
  int           n, cap_m;
  double        arg1, arg2, arg3, arg4, arg5, arg6;
  double        variance_term, delta, snr_in, sqrt_snr, sigma_n, sum;
//
  if(_quantizerLevels == 2)
    return varianceII2Level();
  cap_m         = _quantizerLevels/2 - 1;
  snr_in        = pow(10., _inputSNR/10.);      // SNR as a ratio
  sqrt_snr      = sqrt(snr_in);
  sigma_n       = _noiseAuto->outputNoisePower();
  if(sigma_n>0.)
    sigma_n     = sqrt(sigma_n);
  delta         = _quantizerStepSize;
//
// Sum the terms:
//
  arg5          = cap_m*delta/sigma_n - sqrt_snr;
  arg6          = -cap_m*delta/sigma_n - sqrt_snr;
  sum           = 0.;
  for(n=1; n<=cap_m; n++)
  {
    arg1        = n*delta/SQRT2/sigma_n + sqrt_snr/SQRT2;
    arg2        = (n-1)*delta/SQRT2/sigma_n + sqrt_snr/SQRT2;
    if(n==1)
    {
      arg3      = delta/SQRT2/sigma_n - sqrt_snr/SQRT2;
      arg4      = sqrt_snr/SQRT2;
      sum       += (erf(arg1) - erf(arg2) + erf(arg3) + erf(arg4))/2.;
    }
    else
    {
      arg3      = n*delta/SQRT2/sigma_n - sqrt_snr/SQRT2;
      arg4      = (n-1)*delta/SQRT2/sigma_n - sqrt_snr/SQRT2;
      sum       += n*n*(erf(arg1) - erf(arg2) + erf(arg3) - erf(arg4))/2.;
    }
  }
  variance_term = sum + (cap_m+1)*(cap_m+1)*(qx(arg5)+px(arg6));
  return variance_term;
}
