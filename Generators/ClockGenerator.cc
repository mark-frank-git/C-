/********************************************************************************
 *                                                                              *
 * This class generates clock data.  The default behavior is that after reset   *
 * the clock puts out a high (1), and then a 0 after the high period (depends   *
 * on the duty cycle.                                                           *
 *                                                                              *
 * File: /User/frank/C++/Generators/ClockGenerator.h                            *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "ClockGenerator.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )


// ############################ Private Function ###################################
// driftClockPeriod - Drift the clock high and low periods.
// Input:               None
// Output:              None
//
// Notes:
// 1. _driftRate is given in Hertz/sec.
// 2. The instance variables, _highPeriod and _lowPeriod are modified if _driftRate != 0
// 3. We should probably ensure that clock period does not drift lower than sample period.
//    Still need to do this.
// ############################ Private Function ###################################
void ClockGenerator::driftClockPeriod()
{
  double        drift_freq, freq_ratio;
  double        old_freq, old_period;
  double        high_freq, low_freq;
//
// Calculate drift for current sample period:
//
  old_period            = _highPeriod + _lowPeriod;
  if( (_highPeriod > 0.) && (_lowPeriod > 0.) )
  {
    old_freq            = 1./old_period;
    drift_freq          = _driftRate*_samplePeriod;     // # of Hertz of drift
    freq_ratio          = drift_freq/old_freq;
    high_freq           = 1./_highPeriod;               // use reciprocals to avoid round off
    high_freq           += high_freq*freq_ratio;
    if(high_freq > 0.)
      _highPeriod       = 1./high_freq;
    low_freq            = 1./_lowPeriod;                // use reciprocals to avoid round off
    low_freq            += low_freq*freq_ratio;
    if(low_freq > 0.)
      _lowPeriod        = 1./low_freq;
  }
  return;
}

// ############################# Class Constructor #################################
// ClockGenerator -- Constructor for the ClockGenerator class
// Input:       period:         Full cycle period of the clock in seconds
//              sampling:       Sampling frequency in Hz
// Output:                      None
// ############################# Class Constructor #################################
ClockGenerator::ClockGenerator(double period, double samplingFreq)
{
//
// Initialize instance variables:
//
  setDriftRate(0.);
  setDutyCycle(0.5);
  setClockPeriod(period);
  setSampleFrequency(samplingFreq);
  setTruncateClockPeriod(YES);
  initClock();

  return;
}

// ############################# Class Destructor ###############################
// ClockGenerator -- Destructor for the ClockGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
ClockGenerator::~ClockGenerator()
{  
  return;
}

// ############################ Public Function ###################################
// initClock - Initializes the clock generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void ClockGenerator::initClock()
{
  _clockTimer           = -_samplePeriod;               // Makes first clock high period OK
  _clockValue           = 1;
  _highPeriod           = _dutyCycle*_clockPeriod;      // This gets drifted
  _lowPeriod            = _clockPeriod - _highPeriod;   // This also gets drifted
  return;
}

// ############################ Public Function ###################################
// setDriftRate - Sets a new clock drift rate in Hz/sec.
// Input:       drift:          New drift in Hz/sec
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void ClockGenerator::setDriftRate(float drift)
{
  _driftRate    = drift;
  return;
}

// ############################ Public Function ###################################
// setDutyCycle - Sets a new duty cycle in [0,1].
// Input:       duty:           New duty cycle as a fraction of 1
// Output:                      None
//
// Notes:
// 1. A duty cycle of 0.8 will mean that the clock stays high 80% of the time.
// ############################ Public Function ###################################
void ClockGenerator::setDutyCycle(float duty)
{
  _dutyCycle    = MAX(MIN_DUTY, duty);
  _dutyCycle    = MIN(MAX_DUTY, _dutyCycle);
  return;
}

// ############################ Public Function ###################################
// setClockPeriod - Sets a new full cycle period in seconds.
// Input:       period:         New full cycle period in seconds
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void ClockGenerator::setClockPeriod(double period)
{
  _clockPeriod  = MAX(MIN_PERIOD, period);
  return;
}

// ############################ Public Function ###################################
// setClockFrequency - Sets a new clock frequency in Hertz.
// Input:       frequency:      New clock frequency in Hertz
// Output:                      None
//
// Notes:
// 1. This is a convenience method for those who don't want to convert frequency to
//    time
// ############################ Public Function ###################################
void ClockGenerator::setClockFrequency(double frequency)
{
  if(frequency > 0.)
    setClockPeriod(1./frequency);
  return;
}

// ############################ Public Function ###################################
// setSamplePeriod - Sets the underlying sampling clock period in seconds.
// Input:       period:         New sampling clock period in seconds
// Output:                      None
//
// Notes:
// 1. There is no error checking to insure that the sampling clock period is shorter
//    then the clock generator clock period.
// ############################ Public Function ###################################
void ClockGenerator::setSamplePeriod(double period)
{
  _samplePeriod = MAX(MIN_PERIOD, period);
  return;
}

// ############################ Public Function ###################################
// setSampleFrequency - Sets the underlying sampling frequency in Hertz.
// Input:       frequency:      New sampling frequency in Hertz
// Output:                      None
//
// Notes:
// 1. This is a convenience method for those who don't want to convert frequency to
//    time
// ############################ Public Function ###################################
void ClockGenerator::setSampleFrequency(double frequency)
{
  if(frequency > 0.)
    setSamplePeriod(1./frequency);
  return;
}

// ############################ Public Function ###################################
// setTruncateClockPeriod - Sets whether or not to truncate the clock period to the
//                          next lowest sample period.
// Input:       flag:           If TRUE, truncate, otherwise round
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void ClockGenerator::setTruncateClockPeriod(LOGICAL flag)
{
  _truncNotRound        = flag;
  return;
}

// ############################ Public Function ###################################
// newData - Gets a new clock output from the generator.
// Input:               None
// Output:              New clock data: 0 or 1
//
// Notes:
//
// ############################ Public Function ###################################
int ClockGenerator::newData()
{
//
// Update the counter:
//
  _clockTimer           += _samplePeriod;
//
// Update the clock value based on counter:
//
  if(_clockValue == 1)
  {
    if(_clockTimer > _highPeriod)               // This works for rounded or trunc behavior
    {
      _clockTimer       -= _highPeriod;         // Reset timer with remainder
      _clockValue       = 0;
    }
    else if(!_truncNotRound)                    // Rounded behavior
    {
      if( (_highPeriod - _clockTimer) < _samplePeriod/2.)
      {
        _clockTimer     -= _highPeriod;
        _clockValue     = 0;
      }
    }
  }
  else                                          // Clock value == 0
  {
    if(_clockTimer > _lowPeriod)                // This works for rounded or trunc behavior
    {
      _clockTimer       -= _lowPeriod;          // Reset timer with remainder
      _clockValue       = 1;
    }
    else if(!_truncNotRound)                    // Rounded behavior
    {
      if( (_lowPeriod - _clockTimer) < _samplePeriod/2.)
      {
        _clockTimer     -= _lowPeriod;
        _clockValue     = 1;
      }
    }
  }
//
// Update the clock period based on drift:
//
  driftClockPeriod();
//
// Return the current value of the clock:
//
  return _clockValue;
}
