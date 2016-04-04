#ifndef _IQMEASURE_H
#define _IQMEASURE_H 1
/************************************************************************
 *                                                                      *
 * This is used for making measurements on QAM I and Q symbols.  The    *
 * measurements are geared for efficiency in that variables common to   *
 * multiple measurements should only be calculated once.                *
 *                                                                      *
 *                                                                      *
 * File: /User/frank/C++/QAM/IQMeasure.h                                *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 07/28/97 - Started.                                              *
 *  2. 01/23/98 - Added processIQSymbolInRing().                        *
 *  3. 03/10/98 - Add m_targetVariance for finding normalized phase     *
 *                jitter, normalized MER.                               *
 *  4. 04/06/98 - Make # of bits settable.                              *
 *  5. 04/09/98 - Add idealSigPower().                                  *
 ************************************************************************/
#include <GNU/SmplStat.h>

#define DEFAULT_BITS            6                               // # of qam bits
#define MAX_BITS                8                               // Maximum # of bits
#define DEFAULT_SCALE           0.25                            // default QAM values are 0.25, 0.75, etc.
#define FILTER_ALPHA            0.8                             // For averaging signal power
#define RAD_TO_DEG              (57.29577951)                   // Radians to degree


class IQMeasure
{
protected:
  int           m_qamBits;                              // # of bits per constellation point
  int           m_qamLevels;                            // # of levels in I and Q dimensions, not valid for odd const
  int           m_qamPoints;                            // # of constellation points
  int           m_numberMeanPoints;                     // Counter for statistics
  int           m_numberPowerPoints;                    // Counter for average power
  int           *m_numberTargetPoints;                  // For calculating amplitude imbalance statistics

  double        m_meanI;                                // mean value of I channel
  double        m_meanQ;                                // mean value of Q channel
  double        *m_meanTargetI;                         // mean values of I points
  double        *m_meanTargetQ;                         // mean values of Q points
  double        *m_meanTargetErrorI;                    // mean values of I target error points
  double        *m_meanTargetErrorQ;                    // mean values of Q target error points
  double        *m_targetVarianceI;                     // variances of I target error points
  double        *m_targetVarianceQ;                     // variances of Q target error points
  double        *m_idealTargetI;                        // ideal I channel point for amplitude imbalance
  double        *m_idealTargetQ;                        // ideal Q channel point for amplitude imbalance
  double        *m_slopeI;                              // slopes of curve fit lines parallel to I axis
  double        *m_interceptI;                          // intercepts of curve fit lines parallel to I axis
  double        *m_slopeQ;                              // slopes of curve fit lines parallel to Q axis
  double        *m_interceptQ;                          // intercepts of curve fit lines parallel to Q axis
  double        m_meanSigPower;                         // Mean constellation power
  double        m_currentMeanPower;                     // Current accumulated signal power
  double        m_modulationError;                      // constellation error
  double        m_maxConstellationPt;                   // Maximum constellation point magnitude
  double        m_qamMinScale;                          // Minimum constellation point I or Q magnitude;
  double        m_meanErrorAngle;                       // Mean error angle for phase jitter
  double        m_varianceErrorAngle;                   // Variance of error angle for phase jitter

  SampleStatistic       *m_merStatistics;               // Statistics for MER measurement
  SampleStatistic       *m_evmStatistics;               // Statistics for EVM measurement
  SampleStatistic       *m_suppressionStatistics;       // Statistics for carrier suppression measurement
  SampleStatistic       *m_imbalanceStatistics;         // Statistics for amplitude imbalance measurement
  SampleStatistic       *m_quadErrorStatistics;         // Statistics for quadrature error
  SampleStatistic       *m_phaseJitterStatistics;       // Statistics for phase jitter
  SampleStatistic       *m_normalizedMERStatistics;     // Statistics for MER with phase jitter removed
  SampleStatistic       *m_normalizedPhaseJitterStatistics;     // Statistics for Phase jitter with 'MER' removed

//
// Private methods
//
  void          updateTargetStatistics();
  void          leastSquaresFit(double *xData, double *yData, int numberDataPoints, 
                       double *lineSlope, double *lineIntercept);
  double        findImbalance();
  double        findQuadError();
  void          findNormalizedMERAndPJ(double *mer, double *pj);
  void          processIQSymbolAtIndex(double iData, double qData, double iIdeal, double qIdeal, int qamIndex);
  void          allocateArrays();

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  IQMeasure(int bits=DEFAULT_BITS, float scale=DEFAULT_SCALE);
  virtual ~IQMeasure();

/********************************
 * These methods reset          *
 * statistics and counters:     *
 ********************************/
  void resetCounters();                         // Resets counters for a new trial
  void resetStatistics();                       // Resets statistics variables for a new channel

/*******************************
 * These methods get parameters*
 *******************************/
  int       numberMeanPoints()                  {return m_numberMeanPoints;}

  double    meanI()                             {return m_meanI;}
  double    meanQ()                             {return m_meanQ;}
  double    meanSigPower()                      {return m_meanSigPower;}
  double    modulationError()                   {return m_modulationError;}
  double    maxConstellationPt()                {return m_maxConstellationPt;}

/*******************************
 * These methods set parameters*
 *******************************/
  void          setQAMMinScale(float scale);
  void          setQAMBits(int bits);

/********************************
 * These methods get outputs    *
 * from the statistics          *
 * accumulators, they should be *
 * called after processing      *
 * symbols, and updating stats. *
 ********************************/
  double    meanMER()                           {return m_merStatistics->mean();}
  double    stdDevMER()                         {return m_merStatistics->stdDev();}
  double    meanNormalizedMER()                 {return m_normalizedMERStatistics->mean();}
  double    stdDevNormalizedMER()               {return m_normalizedMERStatistics->stdDev();}
  double    meanEVM()                           {return m_evmStatistics->mean();}
  double    stdDevEVM()                         {return m_evmStatistics->stdDev();}
  double    meanSuppression()                   {return m_suppressionStatistics->mean();}
  double    stdDevSuppression()                 {return m_suppressionStatistics->stdDev();}
  double    meanImbalance()                     {return m_imbalanceStatistics->mean();}
  double    stdDevImbalance()                   {return m_imbalanceStatistics->stdDev();}
  double    meanQuadError()                     {return m_quadErrorStatistics->mean();}
  double    stdDevQuadError()                   {return m_quadErrorStatistics->stdDev();}
  double    meanPhaseJitter()                   {return m_phaseJitterStatistics->mean();}       // In radians
  double    stdDevPhaseJitter()                 {return m_phaseJitterStatistics->stdDev();}     // In radians
  double    meanNormalizedPJ()                  {return m_normalizedPhaseJitterStatistics->mean();} // In radians
  double    stdDevNormalizedPJ()                {return m_normalizedPhaseJitterStatistics->stdDev();} // In radians

/************************************
 * These methods process the input  *
 * I/Q symbol data.                 *
 ************************************/
  void          processIQSymbol(double iData, double qData, double iIdeal, double qIdeal);
  void          processIQSymbolInRing(double iData, double qData, double iIdeal, double qIdeal, int ring);

/************************************
 * These methods accumulate the     *
 * statistics based on processed    *
 * I/Q symbols.                     *
 ************************************/
  void          updateStatistics();

/****************************************
 * Miscellaneous methods.               *
 ****************************************/
  double        idealSigPower();
  

};

#endif
