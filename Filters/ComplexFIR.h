#ifndef _COMPLEX_FIR_H
#define _COMPLEX_FIR_H    1
/************************************************************************
 *                                                                      *
 * This subclass of ComplexDigitalFilter adds functionality for         *
 * filtering using an FIR structure.  For calculating filter            *
 * coefficients, see the subclasses.                                    *
 *                                                                      *
 * File:ComplexFIR.h                                                    *
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


class Complex;                                  // class prototype
class ComplexDigitalFilter;

#ifndef YES
#define YES     1
#define NO      0
#endif
#define         MIN_TAPS        1
#define         MAX_TAPS        50000

class ComplexFIR: public ComplexDigitalFilter
{
protected:
  int           _numberTaps;                            // taps = order + 1
  int           _shiftPtr;                              // Shift buffer pointer


//
// Private functions:
//

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  ComplexFIR(Complex *bCoeffs, int numberCoeffs);
  virtual ~ComplexFIR();

/**********************
 * Set parameters:    *
 **********************/
  virtual void  setNumberTaps(int taps);
  virtual void  setFilterOrder(int order);
/* Inherited from ComplexDigitalFilter ******************
  void          setSamplingFrequency(double fs)         {_fs = fs; return;}
  virtual       void    setPassBandGain(LOGICAL usePolesZeros=NO);
  virtual       void    zeroOutTaps();
  void          setFilterACoeffs(const Complex *a);
  void          setFilterBCoeffs(const Complex *b);
  void          scaleCoefficientsByPassBandGain();
*********************************************************/
  
/**********************
 * Get parameters:    *
 **********************/
  int           numberTaps()                            {return _numberTaps;}
  const         Complex *filterPoles();
  const         Complex *filterZeros();
/* Inherited from ComplexDigitalFilter ******************
  double        fs();
  ******************************************************/

/****************************************
 * Filtering float and Complex data     *
 * assuming an FIR structure.  These    *
 * are overridden from the superclass.  *
 ****************************************/
  void          setPassBandGain(LOGICAL usePolesZeros);
  void          zeroOutTaps();
  Complex       filterFloatData(float input);                   // Filter float data
  Complex       filterDoubleData(double input);                 // Filter double data
  LOGICAL       filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  Complex       filterComplexData(const Complex &input);        // Filter complex data

/**********************
 * Getting Outputs:   *
 **********************/
  double        findPower();                                    // Find and return the filter power


};
#endif
