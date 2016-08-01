#ifndef _SINE_GENERATOR_H
#define _SINE_GENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates sine wave data.                                         *
 *                                                                              *
 * File: /User/frank/C++/Generators/SineGenerator.h                             *
 *                                                                              *
 ********************************************************************************/

#if defined(WIN32)
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#else
#include "Complex.h"
#include "constants.h"
#endif

#define SINE_WAVE_GEN   0                       // _sineGenType
#define COS_WAVE_GEN    1
#define EXP_WAVE_GEN    2


class SineGenerator
{
protected:
  int           _sineGenType;                   // Type of generator.
  float         _initialPhase;                  // Initial starting phase in degrees
  double        _deltaTime;                     // 1/sampling
  double        _omega;                         // Frequency in radians/s
  double        _theta;                         // Accumulated angle
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  SineGenerator(double frequency=10., double sampling=100.);    // Constructor
  virtual ~SineGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setSineGenType(int type)                {_sineGenType   = type; return;}
  void  setInitialPhase(float phase);
  void  setSineFrequency(double frequency);
  void  setSampling(double sampling);

/*******************************
 * These methods get parameters*
 *******************************/
  int   sineGenType()                           {return _sineGenType;}
  float initialPhase()                          {return _initialPhase;}
  float sineFrequency()                         {return _omega/TWOPI;}

/********************************
 * These methods generate       *
 * new Pulse data.              *
 ********************************/
  Complex       newData();                      // Gets the next sample from generator

};

#endif
