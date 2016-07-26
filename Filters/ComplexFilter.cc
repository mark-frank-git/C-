/************************************************************************
 *                                                                      *
 * This subclass of AbstractFilter is for filters with complex          *
 * coefficients.                                                        *
 *                                                                      *
 * File:ComplexFilter.cc                                                *
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
 * Revision history:                                                    *
 *  1. 10/08/02  - Started as subclass of AbstractFilter                *
 ************************************************************************/

#include "ComplexFilter.h"                                      // Object prototypes
#if defined(WIN32)
#include <DataIO/TestData.h>
#include <Polynomials/ComplexPoly.h>
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#else
#include "ComplexPoly.h"
#include "Complex.h"
#include "constants.h"
#include <math.h>
#endif


#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)       ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define SGN01(a)        ( ((a)>=0.) ? 1 :  0 )


// ############################# Private Method ###############################
// transferResponseAt -- Finds the response of the filter at a complex point using
//                       the transfer function representation of the filter.
//
// Input:       point:          Complex point in the s plane
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
Complex ComplexFilter::transferResponseAt(Complex &point)
{
  Complex num, den, complex_response;

//
// Evaluate numerator and denominator polynomials at the point:
//
  num = _bPolyObject->evaluateAtComplexPoint(point);
  den = _aPolyObject->evaluateAtComplexPoint(point);

  if(abs(den) > 0.)
    complex_response = num/den;
  else
    complex_response = num;
  
  return complex_response;
}

// ############################# Private Method ###############################
// readDataFile -- Reads in a data file, and stores the result in 
//                  _fileData, and _numberFileData.
//
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//              real:           YES = coefficients are real, not re im re im
//          
// Output:                      None
//
// Notes:
// 1. Data is stored in _fileData, and _numberFileData.
// ############################# Private Method ###############################
LOGICAL ComplexFilter::readDataFile(const char *fileName, int fileType, LOGICAL realCoeff)
{
  int   i, k;
  const float   *real_coeffs = NULL;
  const float   *imag_coeffs = NULL;
#if defined(WIN32)
  if(_fileReader == NULL)
    _fileReader = new TestData();
  switch(fileType)
  {
    case READ_TWO_BYTE_FILE:
      if(!_fileReader->readTwoByteFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberMagPoints();
      real_coeffs       = _fileReader->getMagnitudeData();
      break;
    case READ_MAT_FILE:
      if(!_fileReader->readMatlabFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberMagPoints();
      real_coeffs       = _fileReader->getMagnitudeData();
      break;
    case READ_WFM_FILE:
      if(!_fileReader->readWFMFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberMagPoints();
      real_coeffs       = _fileReader->getMagnitudeData();
      break;
    case READ_WAVE_FILE:
      if(!_fileReader->readWaveFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberMagPoints();
      real_coeffs       = _fileReader->getMagnitudeData();
      break;
    case READ_ASCII_FILE:
      if(!_fileReader->readASCIIFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberMagPoints();
      real_coeffs       = _fileReader->getMagnitudeData();
      break;
    case READ_TEK_STD_FILE:
      if(!_fileReader->readTekStdFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberIPoints();
      real_coeffs       = _fileReader->getIData();
      imag_coeffs       = _fileReader->getQData();
      break;
    case READ_SONY_TEK_FILE:
      if(!_fileReader->readSonyTekFile(fileName))
        return NO;
      _numberFileData   = _fileReader->getNumberIPoints();
      real_coeffs       = _fileReader->getIData();
      imag_coeffs       = _fileReader->getQData();
      break;
  }
//
// Store coefficients in _fileData;
//
  if(_oldNumberData < _numberFileData)
  {
    _oldNumberData      = _numberFileData;
    delete [] _fileData;
    _fileData           = new Complex [ _numberFileData];
  }
  if(imag_coeffs == NULL)
  {
    if(realCoeff)                                       // Coeffs are real, set imaginary parts = 0
    {
      for(i=0; i<_numberFileData; i++)
      {
        _fileData[i]    = Complex(real_coeffs[i], 0.);
      }      
    }
    else                                                // Assumed stored as re im re im ..
    {
      _numberFileData   /= 2;
      k                 = 0;
      for(i=0; i<_numberFileData; i++)
      {
        _fileData[i]    = Complex(real_coeffs[k], real_coeffs[k+1]);
        k+=2;
      }
    }
  }
  else
  {
    for(i=0; i<_numberFileData; i++)
      _fileData[i]      = Complex(real_coeffs[i], imag_coeffs[i]);
  }
#endif
  return YES;
}

// ############################# Private Method ###############################
// denominatorCenterOfMass -- Finds the center of mass of denominator coefficients.
//
// Input:               None
//          
// Output:              Center of mass (useful for delay calculations)
//
// Notes:
// ############################# Private Method ###############################
float ComplexFilter::denominatorCenterOfMass()
{
  int           i, num_coeffs;
  double        sum_coeff, sum_mass;
  const         Complex *b;
//
  b             = _bPolyObject->getCoefficients();
  num_coeffs    = bOrder() + 1;
  if(num_coeffs > 0)
  {
    sum_coeff   = sum_mass      = 0.;
    for(i=0; i<num_coeffs; i++)
    {
      sum_mass  += i*norm(b[i]);
      sum_coeff += norm(b[i]);
    }
    if(sum_coeff > 0.)
      return    sum_mass/sum_coeff;
  }
  return 0.;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexFilter::ComplexFilter(Complex *aCoeffs, Complex *bCoeffs, int order)
        :AbstractFilter(0, 0., 0., order)
{
  _aPolyObject          = new ComplexPoly();
  _bPolyObject          = new ComplexPoly();
  _fileReader           = NULL;
  _fileData             = NULL;
  _numberFileData       = 0;
  _oldNumberData        = 0;
  setFilterOrder(order);
  if(aCoeffs == NULL)
    setFilterACoeffs(bCoeffs);          // This is to handle the ComplexFIR constructor not having a coeffs
  else
    setFilterBCoeffs(bCoeffs);
  setFilterBCoeffs(bCoeffs);

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the ComplexFilter class.
//
// Input:           type:       Band pass, low pass, etc.
//                  centerFreq: Filter center frequency in Hertz
//                  cutoffFreq: Filter cutoff frequency in Hertz
//                  order:      The order of the filter
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
ComplexFilter::ComplexFilter(int type, double centerFreq, double cutoffFreq, int order)
        :AbstractFilter(type, centerFreq, cutoffFreq, order)
{
  _aPolyObject          = new ComplexPoly();
  _bPolyObject          = new ComplexPoly();
  _fileReader           = NULL;
  _fileData             = NULL;
  _numberFileData       = 0;
  _oldNumberData        = 0;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the ComplexFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
ComplexFilter::~ComplexFilter()
{
  delete _aPolyObject;
  delete _bPolyObject;
  delete [] _fileData;
#if defined(WIN32)
  delete _fileReader;
#endif
  return;
}

// ############################# Public Method ###############################
// setFilterACoeffs -- Sets filter denominator polynomial
// Input:   aCoeff:         coefficients of denominator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void ComplexFilter::setFilterACoeffs(const Complex *aCoeff)
{
  if(_aPolyObject == NULL)
    _aPolyObject = new ComplexPoly(_filterOrder, aCoeff);
  else
    _aPolyObject->assign(_filterOrder, aCoeff);
  _passBandGain = 1.;
  return;
}

// ############################# Public Method ###############################
// setFilterBCoeffs -- Sets filter numerator polynomial
// Input:   bCoeff:         coefficients of numerator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void ComplexFilter::setFilterBCoeffs(const Complex *bCoeff)
{
  if(_bPolyObject == NULL)
    _bPolyObject = new ComplexPoly(_filterOrder, bCoeff);
  else
    _bPolyObject->assign(_filterOrder, bCoeff);
  _passBandGain = 1.;
  return;
}

// ############################# Public Method ###############################
// setFilterPoles -- Sets filter denominator polynomial from pole locations
// Input:       number:         number of poles
//              poles:          Complex pole locations
//          
// Output                       None
//
// Notes:
// 1. _aPolyObject is modified.
// ############################# Public Method ###############################
void ComplexFilter::setFilterPoles(int number, const Complex *poles)
{
  
  if(_aPolyObject == NULL)
    _aPolyObject = new ComplexPoly();
  _aPolyObject->setRoots(number, poles);
  _numberPoles  = number;
  _passBandGain = 1.;
  return;
}

// ############################# Public Method ###############################
// setFilterZeros -- Sets filter numerator polynomial from zero locations
// Input:       number:         number of zeros
//              zeros:          Complex zero locations
//          
// Output                       None
//
// Notes:
// 1. _bPolyObject is modified.
// ############################# Public Method ###############################
void ComplexFilter::setFilterZeros(int number, const Complex *zeros)
{
  
  if(_bPolyObject == NULL)
    _bPolyObject = new ComplexPoly();
  _bPolyObject->setRoots(number, zeros);
  _numberZeros  = number;
  _passBandGain = 1.;
  return;
}

// ############################# Public Method ###############################
// readACoeffsFromFile -- Sets filter denominator polynomial by reading coefficients
//                        from file.
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//              realCoeff:      YES is coefficients are assumed to be real, else re im re im ...
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL ComplexFilter::readACoeffsFromFile(const char *fileName, int fileType, LOGICAL realCoeff)
{
  const Complex *a;
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType, realCoeff))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileData-1);
  setFilterACoeffs(_fileData);

  a     = _aPolyObject->getCoefficients();
  return YES;
}

// ############################# Public Method ###############################
// readBCoeffsFromFile -- Sets filter denominator polynomial by reading coefficients
//                        from file.
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//              realCoeff:      YES is coefficients are assumed to be real, else re im re im ...
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL ComplexFilter::readBCoeffsFromFile(const char *fileName, int fileType, LOGICAL realCoeff)
{
  const Complex *b;
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType, realCoeff))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileData-1);
  setFilterBCoeffs(_fileData);

  b     = _bPolyObject->getCoefficients();
  return YES;
}

// ############################# Public Method ###############################
// readPolesFromFile -- Reads in a set of poles from file, and stores poles
//                      in A coefficients.
// Input:       fileName:       name of file to get poles
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL ComplexFilter::readPolesFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType, NO))
    return NO;
//
// Now, set the poles:
//
  setFilterPoles(_numberFileData, _fileData);
  return YES;
}

// ############################# Public Method ###############################
// readZerosFromFile -- Reads in a set of zeros from file, and stores zeros
//                      in B coefficients.
// Input:       fileName:       name of file to get zeros
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL ComplexFilter::readZerosFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType, NO))
    return NO;
//
// Now, set the zeros:
//
  setFilterZeros(_numberFileData, _fileData);
  return YES;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets the order of the filter
// Input:   order:              New filter order
//          
// Output:                      None
// ############################# Public Method ###############################
void ComplexFilter::setFilterOrder(int order)
{
  AbstractFilter::setFilterOrder(order);
  return;
}

// ############################# Public Method ###############################
// aOrder -- Returns the order of the denominator polynomial
// Input:                   None
//          
// Output:              order of denominator polynomial
//
// ############################# Public Method ###############################
int ComplexFilter::aOrder()
{
  return _aPolyObject->getOrder();
}

// ############################# Public Method ###############################
// bOrder -- Returns the order of the numerator polynomial
// Input:                   None
//          
// Output:              order of numerator polynomial
//
// ############################# Public Method ###############################
int ComplexFilter::bOrder()
{
  return _bPolyObject->getOrder();
}

// ############################# Public Method ###############################
// aCoeffs -- Returns the denominator coefficients for the filter
// Input:                   None
//          
// Output:                  The complex Zeros
// ############################# Public Method ###############################
const Complex *ComplexFilter::aCoeffs()
{
  return _aPolyObject->getCoefficients();
}

// ############################# Public Method ###############################
// bCoeffs -- Returns the numerator coefficients for the filter
// Input:                   None
//          
// Output:                  The complex Zeros
// ############################# Public Method ###############################
const Complex *ComplexFilter::bCoeffs()
{
  return _bPolyObject->getCoefficients();
}
