#ifndef _FIRFILTER_H
#define _FIRFILTER_H    1
/************************************************************************
 *                                                                      *
 * This subclass of AbstractFIR adds functionality for calculating      *
 * the filter coefficients of certain types of FIR                      *
 * filters.                                                             *
 *                                                                      *
 * File:FIRFilter.h                                                     *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[n] + b[n-1]z + ... + b[0]z^n                            *
 *         =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 01/06/98 - Derived from DigitalFilter.                           *
 *  2. 03/05/98 - Override filtering methods from super class.          *
 *  3. 06/09/98 - Add FIR_USER_DEFINED.                                 *
 *  4. 06/15/98 - Ovverride setPassBandGain, findTransferResponse.      *
 *  5. 07/27/98 - Add GMSK filter.                                      *
 *  6. 03/12/99 - Abstracted out AbstractFIR.                           *
 *  7. 09/27/00 - Remove Remez stuff (moved to RemezFIR).               *
 *  8. 10/30/01 - Add output filterDelay().                             *
 *  9. 11/06/01 - Add FIR_EDGE.                                         *
 * 10. 11/29/02 - Corrected, GMSK_CENTER from 0 to 0.5                  *
 ************************************************************************/

#include        "AbstractFIR.h"

//
// _firFilterType
//
#define FIR_WINDOW              0               // FIR window
#define FIR_RAISED_COS          1               // Raised cosine
#define FIR_ROOT_RAISED_COS     2               // Square root raised cosine
#define FIR_POWER_LAW           3               // 1/(f^alpha) response
#define FIR_USER_DEFINED        4               // User input coefficients
#define FIR_GMSK                5               // GMSK filter
#define FIR_EDGE                6               // EDGE modified GMSK filter


#define DEFAULT_ALPHA           0.15            // Root raised cosine alpha value
#define DEFAULT_SYM_SAMP        0.25            // 1 symbol per 4 samples
#define DEFAULT_POWER_ALPHA     1.              // Power law alpha
#define DEFAULT_BT              0.3             // GMSK BT product
#define GMSK_K1                 7.546           // GMSK k1 constant
#define GMSK_CENTER             0.5             // Center of GMSK window
#define MIN_ALPHA               0.01            // Min value for alpha
#define MAX_ALPHA               2.0             // Max value for alpha

#ifndef YES
#define YES     1
#define NO      0
#endif

class Complex;                                  // class prototype
class DataWindow;

class FIRFilter: public AbstractFIR
{
private:
  int           _firFilterType;                         // FIR_WINDOW, etc.

  float         _raisedCosineAlpha;                     // raised cosine alpha value
  float         _symbolsPerSample;                      // raised cosine and GMSK symbols/sample
  float         _powerLawAlpha;                         // 1/(f^alpha)
  float         _gmskBT;                                // GMSK bandwidth symbol product

  DataWindow    *_dataWindow;                           // Window function for FIR filters

//
// Private functions:
//
  double        windowImpulseResponse(int i, double alpha, double wc, double w1, double w2);
                                                        // Used by transferForFIRWindow()
  LOGICAL       transferForFIRWindow();                 // Find transfer function for FIR filter
  LOGICAL       transferForPowerLawFilter();            // Find transfer function for FIR power law filter
  LOGICAL       transferForGMSKFilter();                // Find transfer function for GMSK filter
  LOGICAL       transferForEDGEFilter();                // Find transfer function for EDGE filter
  void          initInstanceVariables();                // Initializes object's instance variables

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  FIRFilter(int pass=LOW_PASS, double centerFreq=0., double cutoffFreq=10.,
                double samplingFreq=100., int order=2);
  ~FIRFilter();

/**********************
 * Set parameters:    *
 **********************/
  void      setFIRFilterType(int type)          {_firFilterType         = type;         return;}
  void      setRaisedCosineAlpha(float alpha)   {_raisedCosineAlpha     = alpha;        return;}
  void      setSymbolsPerSample(float sPs)      {_symbolsPerSample      = sPs;          return;}
  void      setPowerLawAlpha(float alpha)       {_powerLawAlpha         = alpha;        return;}
  void      setGMSKBT(float bt)                 {_gmskBT                = bt;           return;}
  void      setWindowType(int type);

/**********************
 * Get parameters:    *
 **********************/
  int           firFilterType()                 {return _firFilterType;}
  float         raisedCosineAlpha()             {return _raisedCosineAlpha;}
  float         symbolsPerSample()              {return _symbolsPerSample;}
  float         powerLawAlpha()                 {return _powerLawAlpha;}
  float         gmskBT()                        {return _gmskBT;}
  float         filterDelay();
  int           windowType();
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an FIR structure.  These    *
 * are defined in the superclass.       *
 ****************************************/

/**********************
 * Getting Outputs:   *
 **********************/
  void          findTransferFunction();                         // Find transfer function

};
#endif
