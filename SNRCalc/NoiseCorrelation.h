#ifndef _NOISE_CORRELATION_H
#define _NOISE_CORRELATION_H 1
/********************************************************************************
 *                                                                              *
 * This class calculates the autocorrelation value for filtered noise.          *
 *                                                                              *
 * File: /User/frank/C++/SNRCalcu/NoiseCorrelation.h                            *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 10/20/00 - Started.                                                      *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_CUTOFF          5.e6
#define DEFAULT_DENSITY         1.e-7
#define MIN_TAU                 1.e-30

//
// The following are the filter types:
//
#define IDEAL_FILTER            0               // Brickwall filter
#define BUTTER_FILTER           1               // Butterworth filter

class NoiseCorrelation
{
protected:
  int           _filterType;                    // Ideal, Butterworth
  int           _filterOrder;                   // Order of the Butterworth filter

  double        _noiseDensity;                  // Input noise density in Watts/Hz
  double        _filterCutoff;                  // Filter 1 sided bandwidth in Hz
  double        _centerFrequency;               // center frequency of the shifted noise in Hz
//
// Private functions:
//
public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  NoiseCorrelation(int type=IDEAL_FILTER, double density=DEFAULT_DENSITY, double cutoff=DEFAULT_CUTOFF); // Constructor
  virtual ~NoiseCorrelation();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setFilterType(int type)                 {_filterType            = type;         return;}
  void  setFilterOrder(int order)               {_filterOrder           = order;        return;}
  void  setFilterCutoff(double cutoff)          {_filterCutoff          = cutoff;       return;}
  void  setNoiseDensity(double density)         {_noiseDensity          = density;      return;}
  void  setCenterFrequency(double freq)         {_centerFrequency       = freq;         return;}

/*******************************
 * These methods get parameters:*
 *******************************/
  int           filterType()                    {return _filterType;}
  int           filterOrder()                   {return _filterOrder;}
  double        filterCutoff()                  {return _filterCutoff;}
  double        noiseDensity()                  {return _noiseDensity;}
  double        centerFrequency()               {return _centerFrequency;}

/********************************
 * These methods make the       *
 * calculations.                *
 ********************************/
  double        outputNoisePower();                                     // Returns sigma**2
  double        autocorrelationAt(double tau);                          // Returns autocorrelation at tau

};

#endif
