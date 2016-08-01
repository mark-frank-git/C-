/************************************************************************
 *                                                                      *
 * This subclass of AbstractFilter is for filters with real             *
 * coefficients.                                                        *
 *                                                                      *
 * File:RealFilter.cc                                                   *
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
 * or equivalently, for digital filters:                                *
 *                                                                      *
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

#include "RealFilter.h"                                 // Object prototypes

#if defined(WIN32)
#include <Polynomials/DoublePoly.h>
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#include <DataIO/TestData.h>
#else
#include "DoublePoly.h"
#include "Complex.h"
#include "constants.h"
#endif

#include <math.h>

#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)       ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )


// ############################# Private Method ###############################
// transferFromPoles -- Find transfer function form from poles and zeros
// Input:       poles:          Complex poles
//              zeros:          Complex zeros
//              numPoles:       # of complex poles
//              numZeros:       # of complex zeros
//              realAxisPoles:  YES = poles on the real axis
//          
// Output:                      None
//
// Notes:
//  1. The instance variables, _aPolyObject and _bPolyObject are modified.
// ############################# Private Method ###############################
void RealFilter::transferFromPoles(Complex *poles, Complex *zeros, int numPoles, int numZeros,
                                       LOGICAL realAxisPoles)
{
  int           k, n, number_conj_poles;
  double        real_pole, imag_pole, real_zero, imag_zero;
  double        c[1], d[3];
  DoublePoly    x_minus_root;

  n     = numPoles;

/***************************
 * Find the denominator    *
 * from the poles:         *
 ***************************/
  c[0]              = 1.;
  _aPolyObject->assign(0, c);
  
  number_conj_poles = n/2;
  if(realAxisPoles)
    number_conj_poles--;
  for(k=0; k<number_conj_poles; k++)
  {                         
    real_pole       = real(poles[k]);
    imag_pole       = imag(poles[k]);
    d[0]            = 1.;
    d[1]            = -(real_pole + real_pole);
    d[2]            = real_pole*real_pole + imag_pole*imag_pole;
    x_minus_root.assign(2, d);
    (*_aPolyObject)  *= x_minus_root; 
  }   
/*********************************
 * If n is odd, or real axis     *
 * poles from band pass or band  *
 * stop transformations, mult    *
 * by these poles:               *
 *********************************/
  if( n%2 )
  {
    d[0]        = 1.;
    d[1]        = -real(poles[n/2]);
    x_minus_root.assign(1, d);
    (*_aPolyObject)  *= x_minus_root; 
  }
  else if(realAxisPoles)
  {
    for(k=0; k<2; k++)
    {
      d[0]      = 1.;
      d[1]      = -real(poles[n/2-1]);
      if(k == 0)
        d[1]    = -real(poles[n-1]);
      x_minus_root.assign(1, d);
      (*_aPolyObject)    *= x_minus_root; 
    }
  }
                                     
/***************************
 * Find the numerator      *
 * from the zeros:         *
 ***************************/
  n     = numZeros;
  c[0]  = 1.;
  _bPolyObject->assign(0, c);
  
  for(k=0; k<n/2; k++)
  {                         
    real_zero   = real(zeros[k]);
    imag_zero   = imag(zeros[k]);
    d[0]        = 1.;
    if(_filterPassType == BAND_PASS) /* bandpass has all real axis */
    {
      d[1]      = 0.;
      d[2]      = -real_zero*real_zero;
    }
    else
    {
      d[1]      = -(real_zero + real_zero);
      d[2]      = real_zero*real_zero + imag_zero*imag_zero;
    }
    x_minus_root.assign(2, d);
    (*_bPolyObject)  *= x_minus_root; 
  }   
/*********************************
 * If n is odd, mult by the real *
 * axis zero:                    *
 *********************************/
  if( n%2 )
  {
    d[0]        = 1.;
    d[1]        = -real(zeros[n/2]);
    x_minus_root.assign(1, d);
    (*_bPolyObject)  *= x_minus_root; 
  }
  
  return;
}

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
Complex RealFilter::transferResponseAt(Complex &point)
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
// readCoeffFile -- Reads in a coefficient file, and stores the result in 
//                  _fileCoeffs, and _numberFileCoeffs.
//
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      None
// ############################# Private Method ###############################
LOGICAL RealFilter::readCoeffFile(const char *fileName, int fileType)
{
#if defined (WIN32)
  int   i;
  const float   *coeffs = NULL;
  if(_fileReader == NULL)
    _fileReader = new TestData();
  switch(fileType)
  {
    case READ_TWO_BYTE_FILE:
      if(!_fileReader->readTwoByteFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberMagPoints();
      coeffs            = _fileReader->getMagnitudeData();
      break;
    case READ_MAT_FILE:
      if(!_fileReader->readMatlabFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberMagPoints();
      coeffs            = _fileReader->getMagnitudeData();
      break;
    case READ_WFM_FILE:
      if(!_fileReader->readWFMFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberMagPoints();
      coeffs            = _fileReader->getMagnitudeData();
      break;
    case READ_WAVE_FILE:
      if(!_fileReader->readWaveFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberMagPoints();
      coeffs            = _fileReader->getMagnitudeData();
      break;
    case READ_ASCII_FILE:
      if(!_fileReader->readASCIIFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberMagPoints();
      coeffs            = _fileReader->getMagnitudeData();
      break;
    case READ_TEK_STD_FILE:
      if(!_fileReader->readTekStdFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberIPoints();
      coeffs            = _fileReader->getIData();
      break;
    case READ_SONY_TEK_FILE:
      if(!_fileReader->readSonyTekFile(fileName))
        return NO;
      _numberFileCoeffs = _fileReader->getNumberIPoints();
      coeffs            = _fileReader->getIData();
      break;
  }
//
// Store coefficients in _fileCoeffs;
//
  if(_oldNumberCoeffs < _numberFileCoeffs)
  {
    _oldNumberCoeffs    = _numberFileCoeffs;
    delete [] _fileCoeffs;
    _fileCoeffs         = new double [ _numberFileCoeffs];
  }
  for(i=0; i<_numberFileCoeffs; i++)
    _fileCoeffs[i]      = coeffs[i];
  return YES;
#else
  printf("TestData not implemented in RealFilter.cc\n");
  return NO;
#endif
}

// ############################# Private Method ###############################
// writeCoeffFile -- Writes a set of coefficients to file
//
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//              coeffs:         the coefficients to write to file
//              numberCoeffs:   number of coefficients
//          
// Output:                      None
// ############################# Private Method ###############################
LOGICAL RealFilter::writeCoeffFile(const char *fileName, int fileType, const double *coeffs, int numberCoeffs)
{
  LOGICAL       rtn_value;
#if defined(WIN32)
  int           i;
  float         *i_coeffs, *q_coeffs;

  i_coeffs      = new float[numberCoeffs];
  q_coeffs      = new float[numberCoeffs];
  for(i=0; i<numberCoeffs; i++)
  {
    i_coeffs[i] = (float)coeffs[i];
    q_coeffs[i] = 0.;
  }
  if(_fileReader == NULL)
    _fileReader = new TestData();               // Also is a file writer
  
  rtn_value     = YES;
  switch(fileType)
  {
    case WRITE_TWO_BYTE_FILE:
      if(!_fileReader->writeTwoByteFile(fileName, i_coeffs, numberCoeffs))
        rtn_value       = NO;
      break;
    case WRITE_ASCII_FILE:
      if(!_fileReader->writeASCIIFile(fileName, i_coeffs, numberCoeffs))
        rtn_value       = NO;
      break;
    case WRITE_MAT_FILE:
    case WRITE_WFM_FILE:
    case WRITE_WAVE_FILE:
    case WRITE_SONY_TEK_FILE:
    default:
      rtn_value         = NO;                                   // not implemented
      break;
    case WRITE_TEK_STD_FILE:
      if(!_fileReader->writeTekComplexFile(fileName, i_coeffs, q_coeffs, numberCoeffs))
        rtn_value       = NO;
      break;
  }
// delete allocated arrays:
  delete [] i_coeffs;
  delete [] q_coeffs;
#else
  printf("TestData not implemented in RealFilter.cc\n");
  rtn_value     = NO;
#endif
  return rtn_value;
}

// ############################# Private Method ###############################
// writePoleZeroToFile -- Writes a set of poles or zeros to file
//
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//              roots:          the poles or zeros to write to file
//              numberRoots:    number of poles or zeros
//          
// Output:                      None
// ############################# Private Method ###############################
LOGICAL RealFilter::writePoleZeroToFile(const char *fileName, int fileType, const Complex *roots, int numberRoots)
{
  LOGICAL       rtn_value;
#if defined(WIN32)
  int           i, j;
  float         *i_roots, *q_roots, *iq_roots;

  i_roots       = new float[numberRoots];
  q_roots       = new float[numberRoots];
  iq_roots      = new float[2*numberRoots];
  for(i=0; i<numberRoots; i++)
  {
    i_roots[i]  = real(roots[i]);
    q_roots[i]  = imag(roots[i]);
  }
  j             = 0;
  for(i=0; i<numberRoots; i++)
  {
    iq_roots[j++]       = real(roots[i]);
    iq_roots[j++]       = imag(roots[i]);
  }

  if(_fileReader == NULL)
    _fileReader = new TestData();               // Also is a file writer
  
  rtn_value     = YES;
  switch(fileType)
  {
    case WRITE_TWO_BYTE_FILE:
      if(!_fileReader->writeTwoByteFile(fileName, iq_roots, 2*numberRoots))
        rtn_value       = NO;
      break;
    case WRITE_ASCII_FILE:
      if(!_fileReader->writeASCIIFile(fileName, iq_roots, 2*numberRoots))
        rtn_value       = NO;
      break;
    case WRITE_MAT_FILE:
    case WRITE_WFM_FILE:
    case WRITE_WAVE_FILE:
    case WRITE_SONY_TEK_FILE:
    default:
      rtn_value         = NO;                                   // not implemented
      break;
    case WRITE_TEK_STD_FILE:
      if(!_fileReader->writeTekComplexFile(fileName, i_roots, q_roots, numberRoots))
        rtn_value       = NO;
      break;
  }
// delete allocated arrays:
  delete [] i_roots;
  delete [] q_roots;
  delete [] iq_roots;
#else
  printf("TestData not implemented in RealFilter.cc\n");
  rtn_value     = NO;
#endif
  return rtn_value;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RealFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
RealFilter::RealFilter(double *aCoeffs, double *bCoeffs, int order)
        :AbstractFilter(0, 0., 0., order)
{
  _aPolyObject          = new DoublePoly();
  _bPolyObject          = new DoublePoly();
  _fileReader           = NULL;
  _fileCoeffs           = NULL;
  _numberFileCoeffs     = 0;
  _oldNumberCoeffs      = 0;
  setFilterOrder(order);
  setFilterACoeffs(aCoeffs);
  setFilterBCoeffs(bCoeffs);

  return;
}
  
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the RealFilter class.
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
RealFilter::RealFilter(int type, double centerFreq, double cutoffFreq, int order)
        :AbstractFilter(type, centerFreq, cutoffFreq, order)
{
  _aPolyObject          = new DoublePoly();
  _bPolyObject          = new DoublePoly();
  _fileReader           = NULL;
  _fileCoeffs           = NULL;
  _numberFileCoeffs     = 0;
  _oldNumberCoeffs      = 0;
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the RealFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
RealFilter::~RealFilter()
{
  delete _aPolyObject;
  delete _bPolyObject;
  delete _fileReader;
  delete [] _fileCoeffs;
  return;
}

// ############################# Public Method ###############################
// setFilterACoeffs -- Sets filter denominator polynomial
// Input:   aCoeff:         coefficients of denominator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void RealFilter::setFilterACoeffs(const double *aCoeff)
{
  
  if(_aPolyObject == NULL)
    _aPolyObject = new DoublePoly(_filterOrder, aCoeff);
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
void RealFilter::setFilterBCoeffs(const double *bCoeff)
{
  
  if(_bPolyObject == NULL)
    _bPolyObject = new DoublePoly(_filterOrder, bCoeff);
  else
    _bPolyObject->assign(_filterOrder, bCoeff);
  _passBandGain = 1.;
  return;
}

// ############################# Public Method ###############################
// readACoeffsFromFile -- Sets filter denominator polynomial by reading coefficients
//                        from file.
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL RealFilter::readACoeffsFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileCoeffs and _numberFileCoeffs
//
  if(!readCoeffFile(fileName, fileType))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileCoeffs-1);
  setFilterACoeffs(_fileCoeffs);
  return YES;
}

// ############################# Public Method ###############################
// readBCoeffsFromFile -- Sets filter denominator polynomial by reading coefficients
//                        from file.
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if OK
// ############################# Public Method ###############################
LOGICAL RealFilter::readBCoeffsFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileCoeffs and _numberFileCoeffs
//
  if(!readCoeffFile(fileName, fileType))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileCoeffs-1);
  setFilterBCoeffs(_fileCoeffs);
  return YES;
}

// ############################# Public Method ###############################
// writeACoeffsToFile -- Writes filter coefficients to file.
//
// Input:       fileName:       name of file to write coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if write OK
// ############################# Public Method ###############################
LOGICAL RealFilter::writeACoeffsToFile(const char *fileName, int fileType)
{
  int   number_coeffs;
  const double *coeffs;
//
// Get the coefficients:
//
  number_coeffs = aOrder()+1;
  coeffs        = aCoeffs();
//
// Write the coefficients to file.
//
  return writeCoeffFile(fileName, fileType, coeffs, number_coeffs);
}

// ############################# Public Method ###############################
// writeBCoeffsToFile -- Writes filter coefficients to file.
//
// Input:       fileName:       name of file to write coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if OK
// ############################# Public Method ###############################
LOGICAL RealFilter::writeBCoeffsToFile(const char *fileName, int fileType)
{
  int   number_coeffs;
  const double *coeffs;
//
// Get the coefficients:
//
  number_coeffs = bOrder()+1;
  coeffs        = bCoeffs();
//
// Write the coefficients to file.
//
  return writeCoeffFile(fileName, fileType, coeffs, number_coeffs);
}

// ############################# Public Method ###############################
// writePolesToFile -- Writes filter poles to file.
//
// Input:       fileName:       name of file to write poles
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if OK
// ############################# Public Method ###############################
LOGICAL RealFilter:: writePolesToFile (const char *fileName, int fileType)
{
  int   number_poles;
  const Complex *poles;
//
// Get the coefficients:
//
  number_poles  = aOrder();
  poles         = filterPoles();
//
// Write the poles to file.
//
  return writePoleZeroToFile(fileName, fileType, poles, number_poles);
}

// ############################# Public Method ###############################
// writeZerosToFile -- Writes filter zeros to file.
//
// Input:       fileName:       name of file to write zeros
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if OK
// ############################# Public Method ###############################
LOGICAL RealFilter:: writeZerosToFile (const char *fileName, int fileType)
{
  int   number_zeros;
  const Complex *zeros;
//
// Get the coefficients:
//
  number_zeros  = bOrder();
  zeros         = filterZeros();
//
// Write the poles to file.
//
  return writePoleZeroToFile(fileName, fileType, zeros, number_zeros);
}

// ############################# Public Method ###############################
// setAbstract_filterOrder -- Sets the order of the filter
// Input:   order:              New filter order
//          
// Output:                      None
// ############################# Public Method ###############################
void RealFilter::setFilterOrder(int order)
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
int RealFilter::aOrder()
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
int RealFilter::bOrder()
{
  return _bPolyObject->getOrder();
}

// ############################# Public Method ###############################
// aCoeffs -- Returns the denominator coefficients for the filter
// Input:                   None
//          
// Output:                  The complex Zeros
// ############################# Public Method ###############################
const double *RealFilter::aCoeffs()
{
  return _aPolyObject->getCoefficients();
}

// ############################# Public Method ###############################
// bCoeffs -- Returns the numerator coefficients for the filter
// Input:                   None
//          
// Output:                  The complex Zeros
// ############################# Public Method ###############################
const double *RealFilter::bCoeffs()
{
  return _bPolyObject->getCoefficients();
}
