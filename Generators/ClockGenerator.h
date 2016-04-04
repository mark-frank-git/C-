#ifndef _CLOCK_GENERATOR_H
#define _CLOCK_GENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates clock data.  The default behavior is that after reset   *
 * the clock puts out a high (1), and then a 0 after the high period (depends   *
 * on the duty cycle.                                                           *
 *                                                                              *
 * File: /User/frank/C++/Generators/ClockGenerator.h                            *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 10/10/00 - Started.                                                      *
 ********************************************************************************/
#define MIN_PERIOD              1.e-30          // min pulse period
#define MIN_DUTY                0.0001          // Minimum duty cycle
#define MAX_DUTY                0.9999

#ifndef LOGICAL
#define LOGICAL char
#endif

class ClockGenerator
{
protected:
  int           _clockValue;                    // Current value of the clock 0 or 1
  
  float         _driftRate;                     // Drift rate of clock in Hz/second
  float         _dutyCycle;                     // duty cycle in [0,1], e.g.,  a value of 0.6
                                                // indicates that the clock spends 60% of time
                                                // in high state

  double        _clockPeriod;                   // Full cycle period of clock in seconds
  double        _samplePeriod;                  // 1/(underlying sampling freq) given in seconds
  double        _highPeriod;                    // (Drifted) amount of time the clock stays high
  double        _lowPeriod;                     // (Drifted) amount of time the clock stays low
  double        _clockTimer;                    // This is the timer for the clock

  LOGICAL       _truncNotRound;                 // If 1 truncate clock period to nearest sample period
                                                // else round to nearest sample period
//
// Private methods
//
  void          driftClockPeriod();             // Change the high and low periods based on drift
public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  ClockGenerator(double period = 1., double samplingFreq=10.);                  // Constructor
  virtual ~ClockGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initClock();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setDriftRate(float drift);
  void  setDutyCycle(float duty);
  void  setClockPeriod(double period);
  void  setClockFrequency(double frequency);    // Reciprocal of setClockPeriod()
  void  setSamplePeriod(double period);
  void  setSampleFrequency(double frequency);   // Reciprocal of setSamplePeriod()
  void  setTruncateClockPeriod(LOGICAL flag);

/*******************************
 * These methods get parameters*
 *******************************/
  int           clockValue()                    {return _clockValue;}
  float         driftRate()                     {return _driftRate;}
  float         dutyCycle()                     {return _dutyCycle;}
  double        clockPeriod()                   {return _clockPeriod;}
  double        clockFrequency()                {return 1./_clockPeriod;}
  double        samplePeriod()                  {return _samplePeriod;}
  double        sampleFrequency()               {return 1./_samplePeriod;}
  LOGICAL       truncateClockPeriod()           {return _truncNotRound;}

/********************************
 * These methods generate       *
 * new clock data.              *
 ********************************/
  int   newData();                              // Gets the next clock output, this function must be called
                                                // every sample time


};

#endif
