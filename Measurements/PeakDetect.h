#ifndef _PEAKDETECT_H
#define _PEAKDETECT_H 1
/************************************************************************
 *                                                                      *
 * This is used for finding the peak position of a set of FFT bins.     *
 * The methods used here are from DSP GURU.                             *
 *                                                                      *
 * File: /User/frank/C++/Measurements/PeakDetect.h                      *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 07/17/99 - Started.                                              *
 ************************************************************************/

#ifndef LOGICAL
#define LOGICAL    char
#endif

#define NO_INTERPOLATE          0       // peak detect types
#define QUADRATIC_METHOD        1
#define BARYCENTRIC_METHOD      2
#define QUINNS_FIRST            3
#define QUINNS_SECOND           4       // least RMS error
#define JAINS_METHOD            5

class PeakDetect
{
protected:
  int           m_detectType;                           // Number of zero crossings

//
// Private methods
//
  inline        double  apFN(float realData1, float reaData0, float imagData1, float imagData0);
  inline        double  tauFN(float x);

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  PeakDetect();
  virtual ~PeakDetect();

/*******************************
 * These methods set parameters*
 *******************************/
  void  setPeakDetectType(int type)             {m_detectType = type; return;}

/****************************************
 * These methods find the peak position *
 * from a set of FFT data.              *
 ****************************************/
  float peakPosition(float *realFFTData, float *imagFFTData, int low_index, int high_index);

};

#endif