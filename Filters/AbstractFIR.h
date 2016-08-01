#ifndef _ABSTRACTFIR_H
#define _ABSTRACTFIR_H    1
/************************************************************************
 *                                                                      *
 * This subclass of DigitalFilter adds functionality for filtering      *
 * using an FIR structure.  For calculating filter coefficients, see    *
 * the subclasses.                                                      *
 *                                                                      *
 * File:AbstractFIR.h                                                   *
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
 ************************************************************************/

#include        "DigitalFilter.h"

class Complex;                                  // class prototype

#ifndef YES
#define YES     1
#define NO      0
#endif
#define         MIN_TAPS        1
#define         MAX_TAPS        50000

class AbstractFIR: public DigitalFilter
{
protected:
  int           _numberTaps;                            // taps = order + 1
  int           _shiftPtr;                              // Shift buffer pointer


//
// Private functions:
//
  Complex       transferResponseAt(Complex &point);     // Return transfer function response at a point
  virtual float outputResponseFor(Complex complexResponse , float omega);

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  AbstractFIR(int pass=LOW_PASS, double centerFreq=0.25, double cutoffFreq=0.5,
                                  double samplingFreq=1., int order=2);
  virtual ~AbstractFIR();

/**********************
 * Set parameters:    *
 **********************/
  virtual void  setNumberTaps(int taps);
  virtual void  setFilterOrder(int order);

/**********************
 * Get parameters:    *
 **********************/
  int           numberTaps()                            {return _numberTaps;}
  const         Complex *filterPoles();
  const         Complex *filterZeros();
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an FIR structure.  These    *
 * are overridden from the superclass.  *
 ****************************************/
  void          setPassBandGain(LOGICAL usePolesZeros);
  void          zeroOutTaps();
  LOGICAL       filterFloatArray(float * input, int numberPts); // Filter float array
  double        filterFloatData(float input);                   // Filter float data
  double        filterDoubleData(double input);                 // Filter double data
  LOGICAL       filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  Complex       filterComplexData(Complex &input);              // Filter complex data
  void          complexFilterComplexArray(Complex *coeff, int taps, Complex *x, int numberPts);
                                                                // Filter complex array with complex coef_fs

/**********************
 * Getting Outputs:   *
 **********************/
  virtual       float   *filterResponseAtFrequencies(float *frequencies, int number, LOGICAL calc=YES);
  virtual       float   filterResponseAtFrequency(float frequency, LOGICAL calc=YES);
  virtual       Complex complexResponseAtFrequency(float frequency, LOGICAL calc=YES);
  virtual void  findTransferFunction();                         // Needs to be overriden in subclasses
  double        findPower();                                    // Find and return the filter power


};
#endif
