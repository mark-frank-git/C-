#ifndef _SAW_GENERATOR_H
#define _SAW_GENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates sawtooth data.                                          *
 *                                                                              *
 * File: /User/frank/C++/Generators/SawGenerator.h                              *
 *                                                                              *
 ********************************************************************************/

#define MIN_WIDTH               1.e-35          // min saw width
#define MIN_PERIOD              1.e-30          // min saw period

class SawGenerator
{
protected:
  int           _sawOn;                         // YES = saw currently on
  int           _sawBits;                       // Number of pulse bits generated
  
  double        _sawWidth;                      // Saw width in seconds
  double        _sawPeriod;                     // Saw period in second
  double        _rampSlope;                     // Slope of saw in volts/second

  double        _sawTimer;                      // Width/period counter
  double        _sawOutput;                     // Saw output
  double        _sampleTime;                    // 1./fs
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  SawGenerator(double width = 1., double period = 10., double slope = 10., double sampling=1.);
  virtual ~SawGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setSawWidth(double width);
  void  setSawPeriod(double period);
  void  setRampSlope(double slope);
  void  setSampling(double sampling);

/*******************************
 * These methods get parameters*
 *******************************/
  double        sawWidth()                      {return _sawWidth;}
  double        sawPeriod()                     {return _sawPeriod;}
  double        rampSlope()                     {return _rampSlope;}
  double        sampleTime()                    {return _sampleTime;}

/********************************
 * These methods generate       *
 * new Saw data.                *
 ********************************/
  double        newData();              // Gets the next data bit from gen


};

#endif
