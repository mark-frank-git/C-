#ifndef _DIGITALFILTER_H
#define _DIGITALFILTER_H    1
/************************************************************************
 *                                                                      *
 * This subclass of RealFilter implements a digital filtering           *
 * object.  Note that this class is somewhat abstract, see the sub-     *
 * classes, FIRFilter and IIRFilter for actual implementations.         *
 * At the present time, it contains code for general purpose IIR        *
 * filtering, but can not calculate the filter coefficients.            *
 *                                                                      *
 * File:DigitalFilter.h                                                 *
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
 ************************************************************************/

#include        "RealFilter.h"

#ifndef YES
#define YES     1
#define NO      0
#endif

#define NUMBER_NBW_POINTS       200                     // Integration points for NBW calculation

class Complex;                                          // class prototype
class DataWindow;
class AnalogFilter;

class DigitalFilter: public RealFilter
{
protected:
  int           _shiftSize;                             // Number of shift stages allocated
  int           _decimateFactor;                        // CIC filter decimation factor (also used for
                                                        // calculating alias response).

  double        _fs;                                    // sampling frequency
  double        *_doubleShift;                          // Real shift register
  Complex       *_complexShift;                         // complex shift register

//
// Private functions:
//
  void  checkShiftSize(int size);                       // Check for size of shifts
  inline double filterNextDouble(double input, const double *a, const double *b, int a_order, int b_order);
                                                        // Local function for filtering
  inline Complex filterNextComplex(Complex input, const double *a, const double *b, int a_order,
                                   int b_order);
                                                        // Local function for filtering
  inline Complex complexFilterNextComplex(Complex input, const Complex *a, const Complex *b, int a_order,
                                     int b_order);      // Local function for filtering
  Complex aliasResponseAt(float frequency);             // Finds alias response at a frequency point

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  DigitalFilter(int pass=LOW_PASS, double centerFreq=0., double cutoffFreq=10.,
                double samplingFreq=100., int order=2);
  DigitalFilter(double *aCoeffs, double *bCoeffs, int order);
  virtual ~DigitalFilter();
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an IIR structure.  A sub-   *
 * class may wish to override these     *
 * with more efficient implementations. *
 ****************************************/
  virtual       LOGICAL filterFloatArray(float * input, int numberPts); // Filter float array
  virtual       double  filterFloatData(float input);                   // Filter float data
  virtual       double  filterDoubleData(double input);                 // Filter double data
  virtual       LOGICAL filterComplexArray(Complex *x, int numberPts);  // Filter complex array
  virtual       LOGICAL complexFilterComplexArray(Complex *a, Complex *b, int n, Complex *x, int numberPts);
                                                                        // Filter complex array with complex coef_fs
  virtual       Complex filterComplexData(Complex &input);              // Filter complex data

/**********************
 * Set parameters:    *
 **********************/
  void          setFilterParameters(double centerFreq, double cutoff, double samplingFreq);
  void          setSamplingFrequency(double fs)         {_fs = fs; return;}
  void          setDecimateFactor(int factor)           {_decimateFactor        = factor;       return;}
  virtual       void    setFilterOrder(int order);
  virtual       void    setPassBandGain(LOGICAL usePolesZeros=NO);
  virtual       void    zeroOutTaps();
  void          setFilterACoeffs(const double *a);
  void          setFilterBCoeffs(const double *b);
  void          scaleCoefficientsByPassBandGain();
  
/* Set functions inherited from RealFilter
  virtual       void setFilterACoeffs(const double *a);
  virtual       void setFilterBCoeffs(const double *b);
  virtual       LOGICAL readACoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL readBCoeffsFromFile(const char *fileName, int fileType);
  virtual       LOGICAL writeACoeffsToFile(const char *fileName, int fileType);
  virtual       LOGICAL writeBCoeffsToFile(const char *fileName, int fileType);
  virtual       LOGICAL writePolesToFile(const char *fileName, int fileType);
  virtual       LOGICAL writeZerosToFile(const char *fileName, int fileType);
  virtual       void setFilterOrder(int order);
 */
  
/**********************
 * Get parameters:    *
 **********************/
  double        fs()                                    {return _fs;}
  int           decimateFactor()                        {return _decimateFactor;}
/*
// Inherited from RealFilter:
//
  int   aOrder();
  int   bOrder();
  const double *aCoeffs();
  const double *bCoeffs();
  const DoublePoly *aPoly()             {return _aPolyObject;}
  const DoublePoly *bPoly()             {return _bPolyObject;}
  virtual const Complex *filterPoles()  {return NULL;}  // override in subclass
  virtual const Complex *filterZeros()  {return NULL;}  // override in subclass
*/

/**********************
 * Getting Outputs:   *
 **********************/
  virtual       float   *filterResponseAtFrequencies(float *frequencies, int number, LOGICAL calc=YES,
                                             LOGICAL normalize=NO);
  virtual       float   filterResponseAtFrequency(float frequency, LOGICAL calc=YES);
  virtual       Complex complexResponseAtFrequency(float frequency, LOGICAL calc=YES);
  virtual       void    findTransferFunction();                 // Find transfer function
  float         findNBW();                                      // Find and return the noise bandwidth
  virtual       double  findPower();                            // Find and return the filter power
  double        findMagSquaredMax();                            // Find the maximum mag squared

};
#endif
