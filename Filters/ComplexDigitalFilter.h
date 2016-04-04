#ifndef _COMPLEX_DIGITAL_FILTER_H
#define _COMPLEX_DIGITAL_FITLER_H       1
/************************************************************************
 *                                                                      *
 * This subclass of ComplexFilter implements a digital filtering        *
 * object.  Note that this class has some capability for designing      *
 * its own coefficients, otherwise, they must be calculated externally, *
 * and loaded in.                                                       *
 *                                                                      *
 * File:ComplexDigitalFilter.cc                                         *
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
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 10/11/99  - Started.                                             *
 *  2. 01/30/02  - Add shiftCoefficientsByFrequency().                  *
 *  3. 03/19/02  - Add convertCoefficientsToZ().                        *
 ************************************************************************/

#include        "ComplexFilter.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

class Complex;                                          // class prototype

class ComplexDigitalFilter: public ComplexFilter
{
protected:
  int           _shiftSize;                             // Number of shift stages allocated


  double        _fs;                                    // sampling frequency
  Complex       *_complexShift;                         // complex shift register

//
// Private functions:
//
  void  checkShiftSize(int size);                       // Check for size of shifts
  inline Complex filterNextDouble(double input, const Complex *a, const Complex *b, int a_order,
                                     int b_order);
                                                        // Local function for filtering
  inline Complex filterNextComplex(Complex input, const Complex *a, const Complex *b, int a_order,
                                     int b_order);
                                                        // Local function for filtering
  virtual float outputResponseFor(Complex complexResponse , float omega);

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  ComplexDigitalFilter(Complex *aCoeffs, Complex *bCoeffs, int order);
  ComplexDigitalFilter(int type=LOW_PASS, double centerFreq=0., double cutoffFreq=10., double samplingFreq=100.,
                       int order=2);
  virtual ~ComplexDigitalFilter();
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an IIR structure.  A sub-   *
 * class may wish to override these     *
 * with more efficient implementations. *
 ****************************************/
  virtual       Complex filterFloatData(float input);                   // Filter float data
  virtual       Complex filterDoubleData(double input);                 // Filter double data
  virtual       LOGICAL filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  virtual       Complex filterComplexData(Complex &input);              // Filter complex data

/**********************
 * Set parameters:    *
 **********************/
  void          setSamplingFrequency(double fs)         {_fs = fs; return;}
  virtual       void    setPassBandGain(LOGICAL usePolesZeros=NO);
  virtual       void    zeroOutTaps();
  void          setFilterACoeffs(const Complex *a);
  void          setFilterBCoeffs(const Complex *b);
  void          scaleCoefficientsByPassBandGain();

/**********************
 * Get parameters:    *
 **********************/
  double                fs()                                    {return _fs;}
  virtual const         Complex *filterPoles();
  virtual const         Complex *filterZeros();

/**********************
 * Getting Outputs:   *
 **********************/
  float     *filterResponseAtFrequencies(float *frequencies, int number);
  float     filterResponseAtFrequency(float frequency);
  Complex   complexResponseAtFrequency(float frequency);
  virtual   void      findTransferFunction();                   // Find transfer function

/************************
 * Calculating filter   *
 * coefficients:        *
 ************************/
  virtual void  shiftCoefficientsByFrequency(float analogFrequency);
  void          convertCoefficientsToZ();               // Convert coefficients in powers of s to z
  void          convertCoefficientsToZUsingBackwardDifference();
  double        warpFrequency(double omega, double fs);

};
#endif
