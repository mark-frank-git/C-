/********************************************************************************
 *                                                                              *
 * This class calculates the autocorrelation value for filtered noise.          *
 *                                                                              *
 * File: /User/frank/C++/SNRCalcu/NoiseCorrelation.h                            *
 *                                                                              *
 ********************************************************************************/
#include <math.h>
#include "NoiseCorrelation.h"
#include <C_Libraries/constants.h>

#define ABS(a)          ((a) >= 0 ? (a) : (-a))

// ############################# Class Constructor #################################
// NoiseCorrelation -- Constructor for the NoiseCorrelation class
// Input:       type:           low pass filter type
//              density:        Noise density in W/Hz
//              cutoff:         1 sided bandwidth in Hz
//
// Output:              None
// ############################# Class Constructor #################################
NoiseCorrelation::NoiseCorrelation(int type, double density, double cutoff)
{
// 
// Initialize instance variables:
//
  setFilterType(type);
  setNoiseDensity(density);
  setFilterCutoff(cutoff);
  setCenterFrequency(0.);
  setFilterOrder(1);
  return;
}

// ############################# Class Destructor ###############################
// NoiseCorrelation -- Destructor for the NoiseCorrelation class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
NoiseCorrelation::~NoiseCorrelation()
{
//
// Delete space:
//
  return;
}

// ############################ Public Function ###################################
// outputNoisePower - Calculates and returns the output noise power in watts.
// Input:               None
// Output:              Output noise power in watts (sigma**2)
//
// Notes:
// 1. We assume a real (not complex) noise source
// ############################ Public Function ###################################
double NoiseCorrelation::outputNoisePower()
{
  double        noise_bw, noise_power;

  noise_bw              = 0.;
  switch(_filterType)
  {
    case IDEAL_FILTER:
      noise_bw          = _filterCutoff;
      break;
    case BUTTER_FILTER:
      if(_filterOrder > 0)
        noise_bw        = PI*_filterCutoff/2./_filterOrder/sin(PI/2./_filterOrder);
      break;
  }
  noise_power           = noise_bw*_noiseDensity*2.;            // noise_bw is 1 sided
  return noise_power;
}

// ############################ Public Function ###################################
// autocorrelationAt - Calculates and returns the autocorrelation value.
// Input:      tau:     Autocorrelation lag in seconds
// Output:              Autocorrelation value
//
// Notes:
// 1. We assume a real (not complex) noise source
// ############################ Public Function ###################################
double NoiseCorrelation::autocorrelationAt(double tau)
{
  double        b, auto_corr, arg, noise_bw;
  double        omega;

  noise_bw              = 0.;
  tau                   = ABS(tau);                             // Autocorrelation is symmetric about 0
  if( tau < MIN_TAU)
    return outputNoisePower();                                  // R(0)
  b                     = TWOPI*_filterCutoff;                  // B is in radians
  auto_corr             = 0.;
  switch(_filterType)
  {
    case IDEAL_FILTER:
    default:
      auto_corr         = _noiseDensity*sin(b*tau)/PI/tau;
      break;
    case BUTTER_FILTER:
      switch(_filterOrder)
      {
        case    1:
        default:
          auto_corr     = _noiseDensity*b*exp(-b*tau)/2.;
          break;
        case    2:
          arg           = b*tau/SQRT2;
          auto_corr     = _noiseDensity*b*exp(-arg)*cos(arg)*sin(arg)/2/SQRT2;
          break;
      }
      break;
  }
//
// Add in effects due to shift in frequency:
//
  omega                 = _centerFrequency*TWOPI;
  auto_corr             *= cos(omega*tau);
  return auto_corr;
}
