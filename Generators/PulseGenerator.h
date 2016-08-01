#ifndef _PULSE_GENERATOR_H
#define _PULSE_GENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates pulse data.                                             *
 *                                                                              *
 * File: /User/frank/C++/Generators/PulseGenerator.h                            *
 *                                                                              *
 ********************************************************************************/

#define RECTANGULAR_PULSE       0               // pulse types

#define MIN_WIDTH               1.e-35          // min pulse width
#define MIN_PERIOD              1.e-30          // min pulse period

class PulseGenerator
{
protected:
  int           _pulseGenType;                  // Rectangular, etc.
  int           _pulseBits;                     // Counts the # of bits that have occurred
  int           _pulseOn;                       // YES = pulse currently on
  
  float         _pulseWidth;                    // Pulse width in seconds
  float         _pulsePeriod;                   // Pulse period in second

  double        _pulseTimer;                    // Width/period counter
  double        _sampleTime;                    // 1./fs
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  PulseGenerator(float width = 1., float period = 10., float sampling=1.);              // Constructor
  virtual ~PulseGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setPulseGenType(int type)               {_pulseGenType  = type; return;}
  void  setPulseWidth(float width);
  void  setPulsePeriod(float period);
  void  setSampling(float sampling);

/*******************************
 * These methods get parameters*
 *******************************/
  int   pulseGenType()                          {return _pulseGenType;}
  float pulseWidth()                            {return _pulseWidth;}
  float pulsePeriod()                           {return _pulsePeriod;}

/********************************
 * These methods generate       *
 * new Pulse data.              *
 ********************************/
  int   newData();                              // Gets the next data bit from gen
  double pulseOverflow();                       // Returns the overflow of _pulseTimer


};

#endif
