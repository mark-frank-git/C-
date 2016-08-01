/********************************************************************************
 *                                                                              *
 * This class generates sawtooth data.                                          *
 *                                                                              *
 * File: /User/frank/C++/Generators/SawGenerator.h                              *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "SawGenerator.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Class Constructor #################################
// SawGenerator -- Constructor for the SawGenerator class
// Input:       width:          Saw width in seconds
//              period:         Saw period in seconds
//              slope:          Saw slope in volts/second
//              sampling:       Sampling frequency in Hz
// Output:                      None
// ############################# Class Constructor #################################
SawGenerator::SawGenerator(double width, double period, double slope, double sampling)
{
// 
// Initialize instance variables:
//
  setSawWidth(width);
  setSawPeriod(period);
  setSampling(sampling);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// SawGenerator -- Destructor for the SawGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
SawGenerator::~SawGenerator()
{
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the saw generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void SawGenerator::initGenerator()
{
  _sawTimer     = 0.;
  _sawOn        = YES;
  _sawBits      = 1;
  _sawOutput    = 0.;
  return;
}

// ############################ Public Function ###################################
// setSawWidth - Sets a new saw width in seconds.
// Input:       width:          New saw width in seconds
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SawGenerator::setSawWidth(double width)
{
  _sawWidth     = MAX(MIN_WIDTH, width);
  return;
}

// ############################ Public Function ###################################
// setSawPeriod - Sets a new saw period in seconds.
// Input:       period:         New saw period in seconds
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SawGenerator::setSawPeriod(double period)
{
  _sawPeriod    = MAX(MIN_PERIOD, period);
  return;
}

// ############################ Public Function ###################################
// setRampSlope - Sets a new saw ramp in volts/second.
// Input:       slope:          New slope in volts/second
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SawGenerator::setRampSlope(double slope)
{
  _rampSlope    = slope;
  return;
}

// ############################ Public Function ###################################
// setSampling - Sets a new sampling rate in Hertz.
// Input:       sampling:       New sampling frequency in Hz
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void SawGenerator::setSampling(double sampling)
{
  if(sampling > 0.)
    _sampleTime = 1./sampling;
  else
    _sampleTime = 1.;

  return;
}

// ############################ Public Function ###################################
// newData - Gets new data from the Saw generator.
// Input:               None
// Output:              New saw data
//
// Notes:
//
// ############################ Public Function ###################################
double SawGenerator::newData()
{
  double        output_data;
  if(_sawOn)
  {
    if(_sawTimer > ((_sawBits-1)*_sawPeriod+_sawWidth))
    {
      _sawOn                    = NO;
      output_data               = 0;
    }
    else
    {
      _sawOutput                += _sampleTime*_rampSlope;
      output_data               = _sawOutput;
    }
  }
  else
  {
    output_data                 = 0;
    if(_sawTimer > _sawBits*_sawPeriod)
    {
      _sawBits++;
      _sawOn                    = YES;
      _sawOutput                = 0.;
      output_data               = 0.;
    }
  }
  _sawTimer                     += _sampleTime;
  return output_data;
}
