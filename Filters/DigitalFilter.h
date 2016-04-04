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
 *                                                                      *
 * Revision history:                                                    *
 *  1. 12/30/94  - Objective-C -> C++                                   *
 *  2. 01/17/95  - Add filterDoubleData().                              *
 *  3. 02/22/95  - Remove Complex *garb - was used for debug.           *
 *  4. 01/16/96  - typedef int BOOL for DOS.                            *
 *  5. 11/07/96  - Make windowFunction() public.                        *
 *  6. 02/12/97  - Add complexFilterComplexArray().                     *
 *  7. 02/20/97  - Add correlateVectors(), lsInverseFilter().           *
 *  8. 02/20/97  - Reversed order of numerator, denominator polys       *
 *  9. 02/24/97  - Add getWindowScale().                                *
 * 10. 03/12/97  - Add setPassBandGain().                               *
 * 11. 04/29/97  - Move windowing stuff to DataWindow class.            *
 * 12. 06/06/97  - Add remez algorithm filter.                          *
 * 13. 07/05/97  - Add raised cosine stuff.                             *
 * 14. 07/14/97  - Add root raised cosine stuff.                        *
 * 15. 12/05/97  - Override setFilterACoeff().                          *
 * 16. 01/05/98  - Make subclass of AbstractFilter.                     *
 * 17. 01/06/98  - Make superclass of IIRFilter and FIRFilter.          *
 * 18. 06/15/98  - Add setFilterACoeffs, setFilterBCoeffs.              *
 * 19. 09/27/00  - Add findNBW.                                         *
 * 20. 12/18/00  - Add aliasResponse().                                 *
 * 21. 03/27/02  - Add findPower.                                       *
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
