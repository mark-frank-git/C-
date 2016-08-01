/********************************************************************************
 *                                                                              *
 * This class generates sine wave data.                                         *
 *                                                                              *
 * File: /User/frank/C++/Generators/SineGenerator.h                             *
 *                                                                              *
********************************************************************************/
#include <stdio.h>
#include "SineGenerator.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Class Constructor #################################
// SineGenerator -- Constructor for the SineGenerator class
// Input:       frequency:      Sine frequency in Hz
//              sampling:       Sampling frequency in Hz
// Output:              None
// ############################# Class Constructor #################################
SineGenerator::SineGenerator(double frequency, double sampling)
{
// 
// Initialize instance variables:
//
  setSineGenType(SINE_WAVE_GEN);
  setInitialPhase(0.);
  setSineFrequency(frequency);
  setSampling(sampling);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// SineGenerator -- Destructor for the SineGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
SineGenerator::~SineGenerator()
{
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the sine wave generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void SineGenerator::initGenerator()
{
  _theta                = _initialPhase*RAD_DEG;                // initial phase in degrees
  return;
}

// ############################ Public Function ###################################
// setInitialPhase - Sets a new initial phase in degrees.
// Input:       phase:          New initial phase in degrees
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SineGenerator::setInitialPhase(float phase)
{
  _initialPhase = phase;
  return;
}

// ############################ Public Function ###################################
// setSineFrequency - Sets a new sine wave frequency in Hertz.
// Input:       frequency:      New frequency in Hertz
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SineGenerator::setSineFrequency(double frequency)
{
  _omega         = TWOPI*frequency;
  return;
}

// ############################ Public Function ###################################
// setSampling - Sets a new sampling rate in Hertz.
// Input:       sampling:       New sampling frequency in Hz
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SineGenerator::setSampling(double sampling)
{
  if(sampling > 0.)
    _deltaTime  = 1./sampling;
  else
    _deltaTime  = 1.;

  return;
}

// ############################ Public Function ###################################
// newData - Gets new data from the Pulse generator.
// Input:               None
// Output:              New pulse data
//
// Notes:
//
// ############################ Public Function ###################################
Complex SineGenerator::newData()
{
  Complex       output_data;
  switch(_sineGenType)
  {
    default:
    case SINE_WAVE_GEN:                 // sin(wt)
      output_data               = Complex(sin(_theta), 0.);
      break;
    case EXP_WAVE_GEN:                  // Exp(jwt)
      output_data               = Complex(cos(_theta), sin(_theta));
      break;
    case COS_WAVE_GEN:                  // cos(wt)
      output_data               = Complex(cos(_theta), 0.);
      break;
   }
   _theta                       += _deltaTime*_omega;
   if(_theta > TWOPI)                                   // prevent overflow
     _theta                     -= TWOPI;
   else if (_theta < -(TWOPI))
     _theta                     += TWOPI;
  return output_data;
}
