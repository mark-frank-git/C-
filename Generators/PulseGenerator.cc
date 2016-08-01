/********************************************************************************
 *                                                                              *
 * This class generates pulse data.                                             *
 *                                                                              *
 * File: /User/frank/C++/Generators/PulseGenerator.h                            *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "PulseGenerator.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Class Constructor #################################
// PulseGenerator -- Constructor for the PulseGenerator class
// Input:       width:          Pulse width in seconds
//              period:         Pulse period in seconds
//              sampling:       Sampling frequency in Hz
// Output:              None
// ############################# Class Constructor #################################
PulseGenerator::PulseGenerator(float width, float period, float sampling)
{
// 
// Initialize instance variables:
//
  setPulseWidth(width);
  setPulsePeriod(period);
  setPulseGenType(RECTANGULAR_PULSE);
  setSampling(sampling);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// PulseGenerator -- Destructor for the PulseGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
PulseGenerator::~PulseGenerator()
{
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the pulse generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PulseGenerator::initGenerator()
{
  _pulseTimer   = _sampleTime;
  _pulseBits    = 1;
  _pulseOn      = YES;
  return;
}

// ############################ Public Function ###################################
// setPulseWidth - Sets a new pulse width in seconds.
// Input:       width:          New pulse width in seconds
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void PulseGenerator::setPulseWidth(float width)
{
  _pulseWidth   = MAX(MIN_WIDTH, width);
  return;
}

// ############################ Public Function ###################################
// setPulsePeriod - Sets a new pulse period in seconds.
// Input:       period:         New pulse period in seconds
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void PulseGenerator::setPulsePeriod(float period)
{
  _pulsePeriod  = MAX(MIN_PERIOD, period);
  return;
}

// ############################ Public Function ###################################
// setSampling - Sets a new sampling rate in Hertz.
// Input:       sampling:       New sampling frequency in Hz
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void PulseGenerator::setSampling(float sampling)
{
  if(sampling > 0.)
    _sampleTime = 1./sampling;
  else
    _sampleTime = 1.;

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
int PulseGenerator::newData()
{
  int   output_data;
  if(_pulseOn)
  {
    if(_pulseTimer > ((_pulseBits-1)*_pulsePeriod+_pulseWidth))
    {
      _pulseOn                  = NO;
      output_data               = 0;
    }
    else
      output_data               = 1;
  }
  else
  {
    output_data                 = 0;
    if(_pulseTimer >= _pulseBits*_pulsePeriod)
    {
      _pulseBits++;
      _pulseOn                  = YES;
      output_data               = 1;
    }
  }
  _pulseTimer                   += _sampleTime;
  return output_data;
}

// ############################ Public Function ###################################
// pulseOverflow - Returns the pulse overflow for use in linear interpolation.
// Input:               None
// Output:              Pulse overflow value
//
// Notes:
// 1. Call this after a new bit comes out.
//
// ############################ Public Function ###################################
double PulseGenerator::pulseOverflow()
{
  double        overflow;

  overflow      = (_pulseTimer-_sampleTime) - (_pulseBits-1)*_pulsePeriod;
  return overflow;
}
