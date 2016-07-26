#ifndef _CDMAMEASURE_H
#define _CDMAMEASURE_H 1
/************************************************************************
 *                                                                      *
 * This is used for making measurements on CDMA QPSK symbols.           *
 *                                                                      *
 * File: /User/frank/C++/QAM/CDMAMeasure.h                              *
 *                                                                      *
 ************************************************************************/
#include "IQMeasure.h"

#define CDMA_BITS       2                               // For QPSK
#define WALSH_LENGTH    64                              // Default Walsh length
#define WALSH_SUMS      20                              // N

class CDMAMeasure: public IQMeasure
{
protected:
  int           m_walshLength;                          // Length of Walsh codes in chips
  int           m_walshCounter;                         // Counter for dot products and norms
  int           m_numberWalshSums;                      // # of sums to accumulate, N in the rho equation
  int           m_sumCounter;                           // counter for # of sums

  double        m_realDotProduct;                       // real part of dot product for rho
  double        m_imagDotProduct;                       // imag part of dot product for rho
  double        m_meanDotProduct;                       // mean value of dot product for rho

  double        m_idealNorm;                            // current accumulation of ideal norm
  double        m_meanIdealNorm;                        // mean value of norm of ideal signal for rho

  double        m_signalNorm;                           // current accumulation of signal norm
  double        m_meanSignalNorm;                       // mean value of norm of input signal for rho

  SampleStatistic       *m_rhoStatistics;               // Statistics for waveform quality measurement
  SampleStatistic       *m_dotProductStatistics;        // Statistics for dot product
  SampleStatistic       *m_signalNormStatistics;        // Statistics for signal norm

//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  CDMAMeasure(int bits=CDMA_BITS, float scale=DEFAULT_SCALE);
  virtual ~CDMAMeasure();

/********************************
 * These methods reset          *
 * statistics and counters:     *
 ********************************/
  void resetCounters();                         // Resets counters for a new trial
  void resetStatistics();                       // Resets statistics variables for a new channel

/*******************************
 * These methods get parameters*
 *******************************/
  int   walshLength()                           {return m_walshLength;}

/*******************************
 * These methods set parameters*
 *******************************/
  void  setWalshSums(int sums);

/********************************
 * These methods get outputs    *
 * from the statistics          *
 * accumulators, they should be *
 * called after processing      *
 * symbols, and updating stats. *
 ********************************/
  double        meanRho()                       {return m_rhoStatistics->mean();}
  double        stdDevRho()                     {return m_rhoStatistics->stdDev();}
  double        meanDotProduct()                {return m_dotProductStatistics->mean();}
  double        stdDevDotProduct()              {return m_dotProductStatistics->stdDev();}
  double        meanSignalNorm()                {return m_signalNormStatistics->mean();}
  double        stdDevSignalNorm()              {return m_signalNormStatistics->stdDev();}

/************************************
 * These methods process the input  *
 * I/Q symbol data.                 *
 ************************************/
  void          processQPSKSymbol(double iData, double qData, double iIdeal, double qIdeal);
  

};

#endif