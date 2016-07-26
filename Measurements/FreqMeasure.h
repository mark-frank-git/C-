#ifndef _FREQMEASURE_H
#define _FREQMEASURE_H 1
/************************************************************************
 *                                                                      *
 * This is used for making frequency measurements on a sine wave.       *
 *                                                                      *
 * File: /User/frank/C++/Measurements/FreqMeasure.h                     *
 *                                                                      *
 ************************************************************************/
#include <GNU/SmplStat.h>

#ifndef LOGICAL
#define LOGICAL    char
#endif

class FreqMeasure
{
protected:
  int           m_numberCrossings;                      // Number of zero crossings
  int           m_numberSamples;                        // number samples since first crossing
  int           m_lastCrossingPoint;                    // Sample # for last crossing
  int           m_currentSign;                          // Current sign of samples +/-

  float         m_samplingFrequency;                    // Sampling frequency in Hertz

  LOGICAL       m_gotFirstCrossing;                     // Did we get a zero crossing

  SampleStatistic       *m_freqStatistics;              // Statistics for carrier frequency

//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  FreqMeasure(float sampleFrequency);
  virtual ~FreqMeasure();

/********************************
 * These methods reset          *
 * statistics and counters:     *
 ********************************/
  void resetCounters();                         // Resets counters for a new trial
  void resetStatistics();                       // Resets statistics variables for a new channel

/*******************************
 * These methods set parameters*
 *******************************/
  void  setSamplingFrequency(float sampleFrequency);

/********************************
 * These methods get outputs    *
 * from the statistics          *
 * accumulators, they should be *
 * called after processing      *
 * symbols, and updating stats. *
 ********************************/
  int           numberCrossings()               {return m_numberCrossings;}
  int           numberSamples()                 {return m_lastCrossingPoint;}
  int           totalSamples()                  {return m_numberSamples;}
  double        stdDevFreq()                    {return m_freqStatistics->stdDev();}
  double        meanFreq()                      {return m_freqStatistics->mean();}

/****************************************
 * These methods process the input      *
 * sample data.                         *
 ****************************************/
  void          processNewData(double sample);
  void          processArrayData(const float *data, int numberSamples);
  void          processQPSKData(const float *iData, const float *qData, int numberSamples);
  
/****************************************
 * These methods accumulate the         *
 * statistics based on processed        *
 * samples.                             *
 ****************************************/
  void          updateStatistics();

};

#endif