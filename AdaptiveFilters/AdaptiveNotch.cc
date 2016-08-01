/************************************************************************
 *                                                                      *
 * This class implements an adaptive notch filter.  It can be used,     *
 * for example, for detecting sine waves in noise.                      *
 *                                                                      *
 * File:AdaptiveNotch.cc                                                *
 *                                                                      *
 * The transfer function of the filter is:                              *
 *                                                                      *
 *           1 + a*z^(-1) + z^(-2)                                      *
 *    H(z) = --------------------                                       *
 *           1 + rho*a*z^(-1) + rho^2*z^(-2)                            *
 *                                                                      *
 * Where a = -2*cos(omega), with omega equal to the notch frequency,    *
 * and rho = contraction factor to ensure that the poles stay inside    *
 * the unit circle.                                                     *
 * See: Y. Xiao and K. Shida, "New gradient-based algorithms for        *
 * adaptive IIR notch filters, Proc. ICSP '98, pp. 453-456.             *
 *                                                                      *
 ************************************************************************/

#include "AdaptiveNotch.h"                                      // Object prototypes

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )


// ############################# Private Method ###############################
// filterOutput -- Returns output of filter given filter input.
//
// Input:       inputData               : Filter input data
//
// Output:                                      : Filter output data
//
// Notes:
// ############################# Private Method ###############################
float AdaptiveNotch::filterOutput(float inputData)
{
  float         filter_output;
//
// inputs[0] = x[n-1], inputs[1] = x[n-2]
//
  filter_output = inputData + _aCoefficient*_inputs[0] + _inputs[1];
  filter_output -= _contractionFactor*(_aCoefficient*_outputs[0] + _contractionFactor*_outputs[1]);
  return filter_output;
}

// ############################# Private Method ###############################
// gradientOutput -- Returns the approximation to the gradient signal
//
// Input:               None
//          
// Output:              gradient signal approximation
//
// Notes:
// ############################# Private Method ###############################
float AdaptiveNotch:: gradientOutput ()
{
  float         gradient_output;
//
// inputs[0] = x[n-1], outputs[0] = y[n-1]
//
  gradient_output       = _inputs[0] - _contractionFactor*_outputs[0];
  return gradient_output;
}

// ############################# Private Method ###############################
// updateFilter -- Updates the filter coefficients, and past inputs and outputs
//
// Input:       output:         Current output
//              gradient:       Current gradient
//          
// Output:                      None
//
// Notes:
// 1. The instance variables, _aCoefficient, _inputs[], and _outputs[] are
//    modified.
// ############################# Private Method ###############################
void AdaptiveNotch::updateFilter(float output, float gradient, float inputData)
{
  float temp;
//
// Update filter coefficient based on _adaptType:
//
  switch(_adaptType)
  {
    case NOTCH_SG:
    default:
      _aCoefficient     = _aCoefficient - _stepSize*SGN(output)*gradient;
      break;
    case NOTCH_PG:
      _aCoefficient     = _aCoefficient - _stepSize*output*gradient;
      break;
    case NOTCH_LMP:
      _aCoefficient     = _aCoefficient - _stepSize*SGN(output)*output*output*gradient;
      break;
    case NOTCH_MNG:
      temp              = 1. + _epsilon*gradient*gradient;
      _aCoefficient     = _aCoefficient - _stepSize*output*gradient/temp;
      break;
    case NOTCH_SG_MNG:
      temp              = 1. + _epsilon*gradient*gradient;
      _aCoefficient     = _aCoefficient - _stepSize*SGN(output)*gradient/temp;
      break;
    case NOTCH_LMP_MNG:
      temp              = 1. + _epsilon*gradient*gradient;
      _aCoefficient     = _aCoefficient - _stepSize*SGN(output)*output*output*gradient/temp;
      break;
  }
//
// Update inputs and outputs:
//
//
// inputs[0] = x[n-1], inputs[1] = x[n-2]
// outputs[0] = y[n-1], outputs[1] = y[n-2]
//
  _inputs[1]    = _inputs[0];
  _inputs[0]    = inputData;
  _outputs[1]   = _outputs[0];
  _outputs[0]   = output;
  
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the AdaptiveNotch class.
//
// Input:       type:           Adaptation type, SG, etc.
//              stepSize:       Step size, mu
//              epsilon:        Step parameter for MNG, etc.
//              contraction:    Contraction factor for ensuring poles inside unit circle
//              aCoeff:         Initial value of a coefficient of filter.
//              sampling:       Sampling frequency in Hertz
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
AdaptiveNotch::AdaptiveNotch(int type, float stepSize, float epsilon, float contraction,
                             float aCoeff)
{
  setAdaptType(type);
  setStepSize(stepSize);
  setEpsilon(epsilon);
  setContractionFactor(contraction);
  setACoefficient(aCoeff);

  resetFilter();
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the AdaptiveNotch class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
AdaptiveNotch::~AdaptiveNotch()
{
  return;
}

// ############################# Public Method ###############################
// resetFilter -- Resets the filter's parameters to initial values
// Input:               None
//          
// Output:              None
// Notes:
// ############################# Public Method ###############################
void AdaptiveNotch::resetFilter()
{
  int   i;
  _aCoefficient = _storedAValue;
  for(i=0; i< NOTCH_ORDER; i++)
    _outputs[i] = _inputs[i]    = 0.;
}

// ############################# Public Method ###############################
// setAdaptType -- Sets a new type of adaptation type.
// Input:       type:           New adaptation type, PG, etc.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void AdaptiveNotch::setAdaptType(int type)
{
  _adaptType    = type;
  return;
}

// ############################# Public Method ###############################
// setStepSize -- Sets a new adaptation step size.
// Input:       size:           New adaptation step size
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void AdaptiveNotch::setStepSize(float size)
{
  _stepSize     = size;
  return;
}

// ############################# Public Method ###############################
// setEpsilon -- Sets a new adaptation parameter, epsilon.
// Input:       epsilon:        New value for epsilon
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void AdaptiveNotch::setEpsilon(float epsilon)
{
  _epsilon      = epsilon;
  return;
}

// ############################# Public Method ###############################
// setContractionFactor -- Sets a new contraction factor, rho.
// Input:       contraction:    New value for rho.
//          
// Output:                      None
// Notes:
// ############################# Public Method ###############################
void AdaptiveNotch::setContractionFactor(float factor)
{
  _contractionFactor    = factor;
  return;
}

// ############################# Public Method ###############################
// setACoeff -- Sets a new value for the filter coefficient.
// Input:       a:              New value for filter coefficient.
//          
// Output:                      None
// Notes:
// 1. Normally, this coefficient adaptively changes, this function sets the reset
//    value for this parameter.
// 2. See also, setNotchFrequency.
// ############################# Public Method ###############################
void AdaptiveNotch::setACoefficient(float a)
{
  _storedAValue = a;
  return;
}

// ############################# Public Method ###############################
// setNotchFrequency -- Sets a new value for the notch frequency in rad/s.
// Input:       notch:          New value for notch frequency from -PI to PI.
//          
// Output:                      None
// Notes:
// 1. The notch frequency set the filter coefficient using routine above.
// ############################# Public Method ###############################
void AdaptiveNotch::setNotchFrequency(float notch)
{
  float a;
  a     = -2.*cos(notch);
  setACoefficient(a);
  return;
}
// ############################# Public Method ###############################
// notchFrequency -- Returns the adapted notch frequency.
// Input:       None
//          
// Output:                      Notch frequency in digital rad/s
// Notes:
// ############################# Public Method ###############################
float AdaptiveNotch:: notchFrequency()
{
  float omega;
  omega = acos(-_aCoefficient/2.);
  return omega;
}

// ############################# Public Method ###############################
// processInput -- Filters input data, returns output, adapts filter.
// Input:       inputData:      New data to be filtered.
//          
// Output:                      Filter output
// Notes:
// ############################# Public Method ###############################
float AdaptiveNotch::processInput(float inputData)
{
  float filter_output, gradient_output;
  
  filter_output         = filterOutput(inputData);              // Output of filter
  gradient_output       = gradientOutput();                     // Gradient signal calc
  updateFilter(filter_output, gradient_output, inputData);      // Update filter coefficients
  return filter_output;
}
