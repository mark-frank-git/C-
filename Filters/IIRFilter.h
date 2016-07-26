#ifndef _IIRFILTER_H
#define _IIRFILTER_H    1
/************************************************************************
 *                                                                      *
 * This subclass of DigitalFilter adds functionality for calculating    *
 * the transfer function and poles and zeros of certain types of IIR    *
 * filters.  The actual filtering functions are defined in the super-   *
 * class.                                                               *
 *                                                                      *
 * File:IIRFilter.h                                                     *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *           b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                       *
 *  H(z)   = ------------------------------------                       *
 *           1    + a[1]z^(-1) + ... + a[n]z^-(n)                       *
 *                                                                      *
 *  and:                                                                *
 *           (z-zero[0]) * (z-zero[1]) ... (z-zero[n_zero])             *
 *    H(z) = ----------------------------------------------             *
 *           (z-pole[0]) * (z-pole[1]) ... (z-pole[n_pole])             *
 *                                                                      *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 * NOTE: The order of the filter is the # taps - 1.                     *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1.  01/06/01  - Derived from DigitalFilter.                         *
************************************************************************/

#include        "DigitalFilter.h"

//
// _iirFilterType
//
#define ANALOG_PROTOTYPE        0               // Bessel, Cheb, etc.
#define IIR_POWER_LAW           1               // 1/(f^alpha) response filter of above
#define CIC                     2               // cascaded integrator comb filter
#define USER_DEFINED_FILTER     3               // User defined filter coefficients

#define DEFAULT_POWER_ALPHA     1.              // Power law alpha
#define DEFAULT_CIC_DECIMATION  1               // CIC filter decimation factor
#define DEFAULT_RIPPLE          3.              // Passband ripple in dB

#define POWER_LAW_BETA          0.99            // move power law filter poles inside unit circle

#ifndef YES
#define YES     1
#define NO      0
#endif

class Complex;                                  // class prototype
class AnalogFilter;

class IIRFilter: public DigitalFilter
{
private:
  int           _iirFilterType;                         // ANALOG_PROTO, IIR_POWER_LAW, etc.

  float         _powerLawAlpha;                         // 1/(f^alpha)
  float         _passbandRipple;                        // Analog prototype passband ripple in dB

  AnalogFilter  *_analogFilter;                         // Used for designing analog prototype filters

//
// Private functions:
//
  void          findDigitalPolesZerosFromAnalog();      // Find the filter's poles and zeros
  LOGICAL       transferForIIRPowerLawFilter();         // Find transfer function for IIR power law filter
  LOGICAL       transferForCICFilter();                 // Cascaded integrator comb filter
  void          initInstanceVariables();                // Initializes object's instance variables
  Complex       residueAtPoint(DoublePoly &numerator, DoublePoly &denominator, DoublePoly &denomRev, Complex &point);

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  IIRFilter(int pass=LOW_PASS, double centerFreq=0., double cutoffFreq=10.,
                double samplingFreq=100., int order=2);
  ~IIRFilter();

  /**********************
 * Set parameters:    *
 **********************/
  void          setIIRFilterType(int type)      {_iirFilterType         = type;         return;}
  void          setPowerLawAlpha(float alpha)   {_powerLawAlpha         = alpha;        return;}
  void          setPassbandRipple(float ripple) {_passbandRipple        = ripple;       return;}

/*
 * Set functions inherited from DigitalFilter
 *   void               setFilterParameters(double centerFreq, double cutoff, double samplingFreq);
 *   void               setSamplingFrequency(double fs)
 *   void               setDecimateFactor(int factor)
 *   virtual    void    setFilterOrder(int order);
 *   virtual    void    setPassBandGain(LOGICAL usePolesZeros=NO);
 *   virtual    void    zeroOutTaps();
 *   void               setFilterACoeffs(const double *a);
 *   void               setFilterBCoeffs(const double *b);
 *   void               scaleCoefficientsByPassBandGain();
 */

/**********************
 * Get parameters:    *
 **********************/
  int           iirFilterType()                 {return _iirFilterType;}
  float         powerLawAlpha()                 {return _powerLawAlpha;}
  const         Complex *filterPoles();
  const         Complex *filterZeros();

/************************
 * Calculating filters: *
 ************************/
  void          findTransferFunction();                 // Find transfer function based on filter type
  double        warpFrequency(double omega, double fs); // warp the input frequency according to fs
  Complex       sToZ(const Complex &sPoint);            // Bilinear transform
  void          convertCoefficientsToZ();               // Convert coefficients in powers of s to z
  void          convertCoefficientsToZUsingBackwardDifference();
  double        findPower();                            // Find filter power

};
#endif
