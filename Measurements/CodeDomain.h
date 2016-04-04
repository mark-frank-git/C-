#ifndef _CODE_DOMAIN_H
#define _CODE_DOMAIN_H 1
/************************************************************************
 *                                                                      *
 * This class is used for making code domain power measurements on CDMA *
 * QPSK symbols.                                                        *
 *                                                                      *
 * File: /User/frank/C++/Measurements/CodeDomain.h                      *
 *                                                                      *
 * Revision History:                                                    *
 *  1. 03/20/99 - Started.                                              *
 *  2. 05/01/99 - Add processComplexSignal().                           *
 *  3. 05/31/99 - Add calculateRho().                                   *
 ************************************************************************/

#define WALSH_LENGTH            64                      // Default Walsh length
#define DEF_SAMPLES_CHIP        2                       // Default samples/chip
#define MIN_CODE_POWER          0.01                    // 20 dB down

class   WalshCodes;                                     // Class prototype for walsh code generator

class CodeDomain
{
protected:
  int           m_samplesPerChip;                       // Samples/chip for input signal
  int           m_walshLength;                          // length of walsh codes
  int           m_accumulatedWalshPeriods;              // Number of Walsh periods over which phase jitter has been acc'd

  double        *m_codePowers;                          // Code domain powers
  double        *m_integrateIWalsh;                     // Integrated I Walsh values
  double        *m_integrateQWalsh;                     // Integrated Q Walsh values
  double        *m_crossIWalsh;                         // Integrated I Walsh values from cross-spreading
  double        *m_crossQWalsh;                         // Integrated Q Walsh values from cross-spreading

  float         m_minCodePower;                         // For EVM calculation
  float         *m_idealIWalsh;                         // Average I values over the Walsh period
  float         *m_idealQWalsh;                         // Average Q values over the Walsh period

  double        m_meanErrorAngle;                       // Average angle error for phase jitter measurement
  double        m_varianceErrorAngle;                   // Variance angle error for phase jitter measurement

  const char    **m_walshCodes;                         // Walsh code functions

  WalshCodes    *m_walshGen;                            // Walsh code generator

//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  CodeDomain(int length=WALSH_LENGTH, int samples=DEF_SAMPLES_CHIP);
  ~CodeDomain();

/********************************
 * These methods reset          *
 * accumulated measurement      *
 * statistics.                  *
 ********************************/
  void  resetStatistics();

/*******************************
 * These methods get parameters*
 *******************************/
  int   walshLength()                   {return m_walshLength;}

/*******************************
 * These methods set parameters*
 *******************************/
  void  setMinCodePower(float power);
  void  setSamplesPerChip(int samples);
  void  setWalshLength(int length);

/********************************
 * These methods get outputs    *
 ********************************/
  const double  *codePowers()           {return m_codePowers;}
  const float   *idealIWalsh()          {return m_idealIWalsh;}
  const float   *idealQWalsh()          {return m_idealQWalsh;}

/****************************************
 * These methods process the input      *
 * I/Q symbol data over a single        *
 * Walsh period to aid in code power    *
 * and EVM calculations.  Results are   *
 * stored in m_integrateIWalsh[] and    *
 * m_integrateQWalsh[].                 *
 ****************************************/
  void  processRealWalshPeriod(const float *inputI);
  void  processComplexWalshPeriod(const float *inputI, const float *inputQ);
  void  processCrossWalshPeriod(const float *inputI, const float *inputQ);
  
/****************************************
 * These methods process input arrays   *
 * to find code domain powers.          *
 ****************************************/
  void  findCodePowers(const float *inputArray, int size);
  void  findCodePowers(const float *inputI, const float *inputQ, int size);     // rotation insensitive
  void  findCodePowersEVM(const float *inputI, const float *inputQ, int size);  // rotation insensitive
  
/****************************************
 * These methods generate reference     *
 * CDMA signals for EVM and phase       *
 * jitter calculations.                 *
 ****************************************/
  void  generateWalshPeriod();
  void  generateWalshPeriodGeneral();
  
/****************************************
 * These methods calculate rotation     *
 * angles for EVM calculation.          *
 ****************************************/
  double calculateRotationAngle();

/****************************************
 * These methods accumulate             *
 * measurement statistics over a single *
 * Walsh period.                        *
 ****************************************/
  void  accumulatePhaseJitter(const float *inputI, const float *inputQ);

/****************************************
 * These methods make measurements      *
 * on the processed data over a single  *
 * Walsh period.                        *
 ****************************************/
  float calculateEVM(const float *inputI, const float *inputQ);
  float calculateEVM(const float *inputI, const float *inputQ, double theta);
  float calculateRho(const float *inputI, const float *inputQ, const char *iPN, const char *qPN);

/****************************************
 * These methods make measurements      *
 * based on accumulated statistics.     *
 ****************************************/
  float calculatePhaseJitter();


};

#endif