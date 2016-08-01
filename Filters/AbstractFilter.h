#ifndef _ABSTRACTFILTER_H
#define _ABSTRACTFILTER_H 1
/************************************************************************
 *                                                                      *
 * This subclass of object implements an abstract filter class.  For    *
 * actual implementations, see the subclasses, DigitalFilter and        *
 * AnalogFilter.                                                        *
 *                                                                      *
 * File:AbstractFilter.h                                                *
 *                                                                      *
 * The filter is stored in two forms:                                   *
 *                                                                      *
 *          b[n] + b[n-1]s + ... + b[0]s**(n)                           *
 *    H(s) = ------------------------------------                       *
 *          1    + a[n-1]s + ... + a[0]s**(n)                           *
 *                                                                      *
 *           (s-zero[0]) * (s-zero[1]) ... (s-zero[n_zero])             *
 *    H(s) = ----------------------------------------------             *
 *           (s-pole[0]) * (s-pole[1]) ... (s-pole[n_pole])             *
 *                                                                      *
 * Or, equivalently for digital filters:                                *
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
 * NOTE: The order of the filter is n, but the number of coefficients   *
 *       is equal to n+1.                                               *
 * Where the a and b polynomials are stored in the polynomial objects.  *
 *                                                                      *
 ************************************************************************/

#define LOW_PASS                0                       // filterPassType
#define HIGH_PASS               1
#define BAND_PASS               2
#define BAND_STOP               3

#define TRANSFER_FUNCTION       0                       // filterStructureType
#define POLE_ZERO               1

#define BUTTERWORTH             0                       // analogType
#define CHEBYSHEV               1
#define BESSEL                  2

#define MAGNITUDE               0                       // response/plot type
#define DB_MAGNITUDE            1
#define PHASE                   2                       // In radians
#define PHASE_DELAY             3
#define GROUP_DELAY             4
#define POLE_ZERO_PLOT          5
#define SIGNAL_RESPONSE         6
#define SIGNAL_FFT              7
#define WINDOW_FUNCTION         8                       // The windowing function
#define PHASE_NORMALIZED        9
#define ALIAS_RESPONSE          10                      // For decimation filters
#define DB_ALIAS_RESPONSE       11

#ifndef LOGICAL
#define LOGICAL char
#endif

#ifndef YES
#define YES                     1
#define NO                      0
#endif

class   Complex;

class AbstractFilter
{
protected:
  int           _filterPassType;                // Type of pass band
  int           _filterStructureType;           // transfer, pole-zero, etc.
  int           _analogType;                    // type of analog filter, Butter, Cheby, etc.
                                                // this is inherited by AnalogFilter and IIRFilter
  int           _responseType;                  // magnitude, phase, etc.
  int           _filterOrder;                   // Order of the filter
  int           _numberPoles, _numberZeros;     // # of poles not always = order
  int           _oldNumberFrequencies;           // number of points for the frequency response
  LOGICAL       _phaseInDegrees;                // Yes == calculate phase in degrees
  float         *_filterResponse;               // Filter's response at input frequencies
  double        _fo;                            // Filter center frequency in Hertz
  double        _fc;                            // Filter cutoff frequency in Hertz
  double        _thetaOld;                      // used in phase response calculations
  double        _subAngle;
  double        _deltaOmega;
  double        _passBandGain;                  // Divisor for unity pass gain

  Complex       *_filterZeros;                  // Filter's zeros
  Complex       *_filterPoles;                  // Filter's poles

//
// Private methods:
//
  virtual       void initInstanceVariables();   // Init the object's instance variables
  virtual       void initPoleZeroArrays();      // Initializes poles and zeros arrays
  Complex       poleZeroResponseAt(Complex &point, Complex *poles, Complex *zeros);
  virtual float outputResponseFor(Complex complexResponse , float omega);
                                                // Convert complex response to mag, phase, etc.
//
// Public methods:
//
public:
/********************************
 * Constructors, destructor     *
 *******************************/
  AbstractFilter(int type=LOW_PASS, double centerFreq=0., double cutoffFreq=10., int order=2);
  virtual  ~AbstractFilter();

/**************************
 * Setting parameters:    *
 **************************/
  void  setPassType(int type);
  void  setFilterStructureType(int type);
  void  setAnalogType(int type);
  void  setResponseType(int type);
  void  setPhaseInDegrees(LOGICAL flag);
  void  setFilterFrequencies(double centerFreq, double cutoff);
  virtual       void setFilterOrder(int order);         // Check this for digitalFilter
  virtual       void setPassBandGain(LOGICAL usePolesZeros=NO);

/**********************
 * Get parameters:    *
 **********************/
  int           filterPassType()        {return _filterPassType;}
  int           filterStructureType()   {return _filterStructureType;}
  int           analogType()            {return _analogType;}
  int           filterOrder()           {return _filterOrder;}
  int           numberOfPoles()         {return _numberPoles;}
  int           numberOfZeros()         {return _numberZeros;}
  int           structureType()         {return _filterStructureType;}
  int           responseType()          {return _responseType;}
  double        fo()                    {return _fo;}
  double        fc()                    {return _fc;}
  double        passBandGain()          {return _passBandGain;}

};

#endif
