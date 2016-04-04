/************************************************************************
 *                                                                      *
 * This is used for making frequency measurements on a sine wave.       *
 *                                                                      *
 * File: /User/frank/C++/Measurements/FreqMeasure.cc                    *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 03/01/99 - Started.                                              *
 ************************************************************************/

#include "FreqMeasure.h"
#include <stdio.h>
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define SGN(a)    ( ((a)>=0.) ? 1 : -1 )

#ifndef YES
#define YES 1
#define NO  0
#endif

// ############################# Class Constructor #################################
// FreqMeasure -- Constructor for the FreqMeasure class
//
// Input:       fs:     Sampling frequency
//
// Output:              None
// ############################# Class Constructor #################################
FreqMeasure::FreqMeasure(float fs)
{
// 
// Initialize instance variables:
//
  setSamplingFrequency(fs);
              
//
// Allocate statistics objects:
//
  m_freqStatistics                      = new SampleStatistic();  

// 
// Reset all counters and statistics:
//
  resetCounters();
  resetStatistics();

  return;
}


// ############################# Class Destructor ###############################
// FreqMeasure -- Destructor for the FreqMeasure class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
FreqMeasure::~FreqMeasure()
{
  delete [] m_freqStatistics;
  
  return;
}

// ############################ Public Function ###################################
// resetCounters - Resets the counters between trials.
//
// Input:                       None
// Output:                      None
//
// Notes:
// 1. This function is also called from resetStatistics().
// ############################ Public Function ###################################
void FreqMeasure::resetCounters()
{
//
// Reset counters:
//
  m_numberCrossings     = 0;
  m_lastCrossingPoint   = 0;
  m_currentSign         = 0;

  m_gotFirstCrossing    = NO;

  return;
}

// ############################ Public Function ###################################
// resetStatistics - Resets the statistics between channel changes.
//
// Input:                       None
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void FreqMeasure::resetStatistics()
{
//
// Add these:
//
  m_freqStatistics->reset();
  resetCounters();

  return;
}

// ############################ Public Function ###################################
// setSamplingFrequency - This sets a new sampling frequency.
//
// Input:       sampleFrequency:        sampling frequency in Hz
//
// Output:                              None
//
// Notes:
// ############################ Public Function ###################################
void FreqMeasure::setSamplingFrequency(float sampleFrequency)
{
  m_samplingFrequency   = sampleFrequency;
  return;
}

// ############################ Public Function ###################################
// processNewData - This method processes a new sample.
//
// Input:       sample:         New sample from a sine wave
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void FreqMeasure::processNewData(double sample)
{
  if(m_currentSign == 0)                        // First sample?
  {
    m_currentSign       = SGN(sample);
    return;
  }
//
// Check if we got our first crossing:
//
  if(!m_gotFirstCrossing)
  {
    if(m_currentSign != SGN(sample))
    {
      m_gotFirstCrossing        = YES;
      m_lastCrossingPoint       = m_numberSamples       = 0;
      m_currentSign             = SGN(sample);
      m_numberCrossings         = 1;
    }
    return;
  }
//
// Check for crossing:
//
  m_numberSamples++;
  if(m_currentSign != SGN(sample))
  {
    m_numberCrossings++;
    m_currentSign       = SGN(sample);
    m_lastCrossingPoint = m_numberSamples;
  }

  return;
}

// ############################ Public Function ###################################
// processArrayData - This method processes an array of signal samples.
//
// Input:       data:           Array of signal data
//              numberSamples:  Size of arrays
//
// Output:                      None
//
// Notes:
// 1. The input signal is raised to the fourth before checking for zero crossings
// ############################ Public Function ###################################
void FreqMeasure::processArrayData(const float *data, int numberSamples)
{
  int   n;

//
// Loop over number of samples:
//
  for(n=0; n<numberSamples; n++)
  {
//
// Check for first sample:
//
    if(m_currentSign == 0)                      // First sample?
    {
      m_currentSign     = SGN(data[n]);
    }
//
// Check if we got our first crossing:
//
    else if(!m_gotFirstCrossing)
    {
      if(m_currentSign != SGN(data[n]))
      {
        m_gotFirstCrossing      = YES;
        m_lastCrossingPoint     = m_numberSamples       = 0;
        m_currentSign           = SGN(data[n]);
        m_numberCrossings       = 1;
      }
    }
//
// Check for crossing:
//
    else
    {
      m_numberSamples++;
      if(m_currentSign != SGN(data[n]))
      {
        m_numberCrossings++;
        m_currentSign   = SGN(data[n]);
        m_lastCrossingPoint     = m_numberSamples;
      }
    }
  }  // End loop over n

  return;
}

// ############################ Public Function ###################################
// processQPSKData - This method processes an array of QPSK samples.
//
// Input:       iData:          Array of I QPSK data
//              qData:          Array of Q QPSK data
//              numberSamples:  Size of arrays
//
// Output:                      None
//
// Notes:
// 1. The input signal is raised to the fourth before checking for zero crossings
// ############################ Public Function ###################################
void FreqMeasure::processQPSKData(const float *iData, const float *qData, int numberSamples)
{
  int   n;
  float i_squared;
  float imag_fourth;

//
// Loop over number of samples:
//
  for(n=0; n<numberSamples; n++)
  {
//
// First raise the signal to the 4th power, and consider only imag part (fewer mults than real):
//
    i_squared   = iData[n]*iData[n];
    imag_fourth = 4.*iData[n]*(i_squared*qData[n] - qData[n]*qData[n]*qData[n]);
//
// Check for first sample:
//
    if(m_currentSign == 0)                      // First sample?
    {
      m_currentSign     = SGN(imag_fourth);
    }
//
// Check if we got our first crossing:
//
    else if(!m_gotFirstCrossing)
    {
      if(m_currentSign != SGN(imag_fourth))
      {
        m_gotFirstCrossing      = YES;
        m_lastCrossingPoint     = m_numberSamples       = 0;
        m_currentSign           = SGN(imag_fourth);
        m_numberCrossings       = 1;
      }
    }
//
// Check for crossing:
//
    else
    {
      m_numberSamples++;
      if(m_currentSign != SGN(imag_fourth))
      {
        m_numberCrossings++;
        m_currentSign   = SGN(imag_fourth);
        m_lastCrossingPoint     = m_numberSamples;
      }
    }
  }  // End loop over n

  return;
}

// ############################ Public Function ###################################
// updateStatistics - This method updates the statistics after a number of samples.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void FreqMeasure::updateStatistics()
{
  float freq, cycles;
//
//
  freq          = 0.;
  if(m_lastCrossingPoint > 0)
  {
    cycles      = (m_numberCrossings-1.)/2.;
    freq        = m_samplingFrequency*cycles/m_lastCrossingPoint;
  }
  (*m_freqStatistics)   += freq;
  return;
}

