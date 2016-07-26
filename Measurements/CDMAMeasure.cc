/************************************************************************
 *                                                                      *
 * This is used for making measurements on QAM I and Q symbols.  The    *
 * measurements are geared for efficiency in that variables common to   *
 * multiple measurements should only be calculated once.                *
 *                                                                      *
 *                                                                      *
 * File: /User/frank/C++/QAM/CDMAMeasure.cc                             *
 *                                                                      *
 ************************************************************************/

#include "CDMAMeasure.h"
#include <stdio.h>
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )


// ############################# Class Constructor #################################
// CDMAMeasure -- Constructor for the CDMAMeasure class
//
// Input:       bits:   # of bits per constellation point
//              scale:  qam scale factor
//
// Output:              None
// ############################# Class Constructor #################################
CDMAMeasure::CDMAMeasure(int bits, float scale)
            :IQMeasure(bits, scale)
{
// 
// Initialize instance variables:
//
  m_walshLength                         = WALSH_LENGTH;
  m_numberWalshSums                     = WALSH_SUMS;
              
//
// Allocate statistics objects:
//
  m_rhoStatistics                       = new SampleStatistic();
  m_dotProductStatistics                = new SampleStatistic();
  m_signalNormStatistics                = new SampleStatistic();
  

// 
// Reset all counters and statistics:
//
  resetCounters();
  resetStatistics();

  return;
}


// ############################# Class Destructor ###############################
// CDMAMeasure -- Destructor for the CDMAMeasure class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
CDMAMeasure::~CDMAMeasure()
{
  delete [] m_rhoStatistics;
  delete [] m_dotProductStatistics;
  delete [] m_signalNormStatistics;
  
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
void CDMAMeasure::resetCounters()
{
//
// Do everything in base class:
//
  IQMeasure::resetCounters();
//
// Add these:
//
  m_walshCounter        = 0;
  m_sumCounter          = 0;
  m_realDotProduct      = m_imagDotProduct      = m_meanDotProduct      = 0.;
  m_idealNorm           = m_meanIdealNorm       = 0.;
  m_signalNorm          = m_meanSignalNorm      = 0.;

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
void CDMAMeasure::resetStatistics()
{
//
// Do everything in base class:
//
  IQMeasure::resetStatistics();
//
// Add these:
//
  m_rhoStatistics->reset();
  m_dotProductStatistics->reset();
  m_signalNormStatistics->reset();
  resetCounters();

  return;
}

// ############################ Public Function ###################################
// setWalshSums - This sets a new number of Walsh sums.
//
// Input:       sums:           new number of Walsh sums for calculating rho
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CDMAMeasure::setWalshSums(int sums)
{
  m_numberWalshSums     = sums;
  return;
}

// ############################ Public Function ###################################
// processQPSKSymbol - This method adds a new I/Q symbol to the statistics counters.
//
// Input:       iData:          I (in-phase) value of the input symbol
//              qData:          Q (quadrature) value of the input symbol
//              iIdeal:         I value of the ideal constellation point
//              qIdeal:         Q value of the ideal constellation point
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void CDMAMeasure::processQPSKSymbol(double iData, double qData, double iIdeal, double qIdeal)
{
  double        rho;

//
// First, find dot product of r and w*:
//
  m_realDotProduct      += iData*iIdeal + qData*qIdeal;
  m_imagDotProduct      += qData*iIdeal - iData*qIdeal;
//
// Now, find dot product of w and w*
  m_idealNorm           += iIdeal*iIdeal + qIdeal*qIdeal;
//
// Now, find dot product of r and r*
//
  m_signalNorm          += iData*iData + qData*qData;
//
// Check for length of dot product:
//
  if(++m_walshCounter == m_walshLength)
  {
    m_walshCounter      = 0;
    m_meanDotProduct    += m_realDotProduct*m_realDotProduct + m_imagDotProduct*m_imagDotProduct;
    m_meanIdealNorm     += m_idealNorm;
    m_meanSignalNorm    += m_signalNorm;
//
// The following are a little redundant with the above:
//
    *m_dotProductStatistics     += m_realDotProduct*m_realDotProduct + m_imagDotProduct*m_imagDotProduct;
    *m_signalNormStatistics     += m_signalNorm;

    m_realDotProduct    = m_imagDotProduct      = 0.;
    m_idealNorm         = m_signalNorm          = 0.;
//
// Note, in calculating rho, mult by N to account for the summing of ||wo||^2
//
    if(++m_sumCounter == m_numberWalshSums)
    {
      m_sumCounter      = 0;
      if( (m_meanSignalNorm > 0.) && (m_meanIdealNorm > 0.) )
      {
        rho                     = m_numberWalshSums*m_meanDotProduct/m_meanIdealNorm/m_meanSignalNorm;
        *m_rhoStatistics        += rho;
      }
      m_meanSignalNorm          = m_meanDotProduct      = m_meanIdealNorm       = 0.;
    }
  }
  return;
}
