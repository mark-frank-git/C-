/************************************************************************
 *                                                                      *
 * This class implements a generalized feedforward feedback filter      *
 * for use in delta sigma modulators.  See 5.6.4 in Delta-Sigma Data    *
 * Converters: Theory, Design, and Simulation.                          *
 *                                                                      *
 * File: FFFBFilter.cc                                                  *
 *                                                                      *
 * The structure of the filter is as follows:                           *
 *                                                                      *
 * x -----------------------------                                      *
 *    |            |              |                                     *
 *    b0           b1             bn                                    *
 *    |     1      |       1      |     1                               *
 *    +  -------   +    -------   +  -------   -----> Quant() ----> y   *
 *   -|  z - p0   -|    z - p1   -|  z - pn                   |         *
 *    |            |              |                           |         *
 *    a0           a1             an                          |         *
 *    |            |              |                           |         *
 *    ---------------------------------------------------------         *
 *                                                                      *
 * where the pi's set the transfer function poles, the bi's set the     *
 * signal transfer function zeros, and the ci's set the noise transfer  *
 * function zeros. Also, note that the coefficients and poles are       *
 * complex in order to effect a bandpass filter.                        *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/24/00  - Started.                                             *
 *  2. 01/05/01  - Modified (corrected) filter operation.               *
 ************************************************************************/

#include "FFFBFilter.h"                                 // Object prototypes
#include <GNU/Complex.h>
#include <C_Libraries/constants.h>
#include <math.h>

#define ABS(a)          ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b)       ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define SGN01(a)        ( ((a)>=0.) ? 1 :  0 )
#define SGN(a)          ( ((a)>=0.) ? 1. : -1.)


// ############################# Private Method ###############################
// readDataFile -- Reads in a data file, and stores the result in 
//                  _fileData, and _numberFileData.
//
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      None
//
// Notes:
// 1. Data is stored in _fileData, and _numberFileData.
// ############################# Private Method ###############################
LOGICAL FFFBFilter::readDataFile(const char *fileName, int fileType)
{
  int   i, k;
  const float   *real_coeffs = NULL;
  const float   *imag_coeffs = NULL;
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
    _numberFileData     /= 2;
    k                   = 0;
    for(i=0; i<_numberFileData; i++)                    // Assumed stored as re im re im ...
    {
      _fileData[i]      = Complex(real_coeffs[k], real_coeffs[k+1]);
      k+=2;
    }
  }
  else
  {
    for(i=0; i<_numberFileData; i++)
      _fileData[i]      = Complex(real_coeffs[i], imag_coeffs[i]);
  }
  return YES;
}

// ############################# Private Method ###############################
// checkShiftSize -- Allocates the arrays for the digital shift registers.
//
// Input:           size:       New shift size needed
//          
// Output:                      None
//
// Notes:
// ############################# Private Method ###############################
void FFFBFilter::checkShiftSize(int size)
{  
//
// Check if we are big enough:
//
  if(_shiftSize < size)
  {
    delete [] _inputShift;
    delete [] _backShift;

    _shiftSize          = size;
    _inputShift         = new Complex[_shiftSize+1];
    _backShift          = new Complex[_shiftSize+1];
    zeroOutTaps();                              // Clear the filter's taps
  }
  return;
}

// ############################# Private Method ###############################
// filterNextDouble -- Filter a single double input variable.
//                     See the block diagram above.
// Input:   input:          Next double input
//          
// Output:                  Filtered and quantized output.
// ############################# Private Method ###############################
Complex FFFBFilter::filterNextDouble(double input)
{
  int           i;
  double        real_out, imag_out;
  Complex       yout, complex_input, old_shift;
  Complex       quant_out, stage_out, last_stage;

  complex_input = Complex(input, 0.);
//
// Get quantized output for feedback:
//
  last_stage    = _inputShift[_filterOrder-1] +
                  _backShift[_filterOrder-1]*_filterPoles[_filterOrder-1];
  real_out      = real(last_stage);
  imag_out      = imag(last_stage);
  if(_quantLevels == 2)
  {
   real_out     = SGN(real_out);
   imag_out     = SGN(imag_out);
  }
  quant_out     = Complex(_feedbackGain*real_out, _feedbackGain*imag_out);
//
// Loop through all the shift register stages starting with the first one:
//
  last_stage    = Complex(0.,0.);
  for(i=0; i<_filterOrder; i++)
  {
    stage_out           = _inputShift[i] + _backShift[i]*_filterPoles[i];
    _inputShift[i]      = complex_input*_bCoeffs[i] - quant_out*_aCoeffs[i] + last_stage;
    _backShift[i]       = stage_out;
    last_stage          = stage_out;
  }
//
// Get the quantized output:
//
  real_out      = real(stage_out);
  imag_out      = imag(stage_out);
  if(_quantLevels == 2)
  {
    real_out    = SGN(real_out);
    imag_out    = SGN(imag_out);
  }
  quant_out     = Complex(real_out, imag_out);

  return quant_out;
}

// ############################# Private Method ###############################
// filterNextComplex -- Filter a single complex input variable.
//                     See the block diagram above.
// Input:   input:          Next complex input
//          
// Output:                  Filtered and quantized output.
// ############################# Private Method ###############################
Complex FFFBFilter::filterNextComplex(Complex complex_input)
{
  int           i;
  double        real_out, imag_out;
  Complex       yout, old_shift, last_stage;
  Complex       quant_out, stage_out;

//
// Get quantized output for feedback:
//
  last_stage    = _inputShift[_filterOrder-1] +
                  _backShift[_filterOrder-1]*_filterPoles[_filterOrder-1];
  real_out      = real(last_stage);
  imag_out      = imag(last_stage);
  if(_quantLevels == 2)
  {
   real_out     = SGN(real_out);
   imag_out     = SGN(imag_out);
  }
  quant_out     = Complex(_feedbackGain*real_out, _feedbackGain*imag_out);
//
// Loop through all the shift register stages starting with the first one:
//
 last_stage     = Complex(0.,0.);
 for(i=0; i<_filterOrder; i++)
 {
   stage_out            = _inputShift[i] + _backShift[i]*_filterPoles[i];
   _inputShift[i]       = complex_input*_bCoeffs[i] - quant_out*_aCoeffs[i] + last_stage;
   _backShift[i]        = stage_out;
   last_stage           = stage_out;
 }
//
// Get the quantized output:
//
  real_out      = real(stage_out);
  imag_out      = imag(stage_out);
  if(_quantLevels == 2)
  {
    real_out    = SGN(real_out);
    imag_out    = SGN(imag_out);
  }
 quant_out      = Complex(real_out, imag_out);

 return quant_out;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the FFFBFilter class.
//
// Input:           aCoeffs:    The coefficients of the denominator polynomial
//                  bCoeffs:    The coefficients of the numerator polynomial
//                  order:      The order of the polynomials
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
FFFBFilter::FFFBFilter(Complex *aCoeffs, Complex *bCoeffs, Complex *poles, int order)
{
  _numberFileData       = 0;
  _oldNumberData        = 0;
  _shiftSize            = -1;
  _fileReader           = NULL;
  _fileData             = NULL;
  _inputShift           = NULL;
  _backShift            = NULL;
  _aCoeffs              = NULL;
  _bCoeffs              = NULL;
  _filterPoles          = NULL;
  setQuantLevels(2);
  setQuantScale(1.);
  setFeedbackGain(1.);
  setFilterOrder(order);
  setFilterACoeffs(aCoeffs);
  setFilterBCoeffs(bCoeffs);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FFFBFilter class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FFFBFilter::~FFFBFilter()
{
  delete _fileReader;
  delete [] _fileData;
  delete [] _inputShift;
  delete [] _backShift;
  delete [] _aCoeffs;
  delete [] _bCoeffs;
  delete [] _filterPoles;
  return;
}

// ############################# Public Method ###############################
// setFilterACoeffs -- Sets filter denominator polynomial
// Input:   aCoeff:         coefficients of denominator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void FFFBFilter::setFilterACoeffs(const Complex *aCoeff)
{
  int   i;
  delete [] _aCoeffs;
  _aCoeffs      = new Complex[_filterOrder];

  if(aCoeff != NULL)
  {
    for(i=0; i<_filterOrder; i++)
    {
      _aCoeffs[i]       = aCoeff[i];
      printf("a[%d] = %g %g\n", i, real(aCoeff[i]), imag(aCoeff[i]));
    }
  }
  return;
}

// ############################# Public Method ###############################
// setFilterBCoeffs -- Sets filter numerator polynomial
// Input:   bCoeff:         coefficients of numerator polynomial
//          
// Output:                  None
// ############################# Public Method ###############################
void FFFBFilter::setFilterBCoeffs(const Complex *bCoeff)
{
  int   i;
  delete [] _bCoeffs;
  _bCoeffs      = new Complex[_filterOrder];

  if(bCoeff != NULL)
  {
    for(i=0; i<_filterOrder; i++)
    {
      _bCoeffs[i]       = bCoeff[i];
      printf("b[%d] = %g %g\n", i, real(bCoeff[i]), imag(bCoeff[i]));
    }
  }
  return;
}

// ############################# Public Method ###############################
// setFilterPoles -- Sets filter denominator polynomial from pole locations
// Input:       poles:          Complex pole locations
//          
// Output                       None
//
// Notes:
// ############################# Public Method ###############################
void FFFBFilter::setFilterPoles(const Complex *poles)
{
  int   i;
  delete [] _filterPoles;
  _filterPoles  = new Complex[_filterOrder];

  if(poles != NULL)
  {
    for(i=0; i<_filterOrder; i++)
    {
      _filterPoles[i]   = poles[i];
      printf("p[%d] = %g %g\n", i, real(poles[i]), imag(poles[i]));
    }
  }
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
LOGICAL FFFBFilter::readACoeffsFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileData);
  setFilterACoeffs(_fileData);

  return YES;
}

// ############################# Public Method ###############################
// readBCoeffsFromFile -- Sets filter denominator polynomial by reading coefficients
//                        from file.
// Input:       fileName:       name of file to get coefficients
//              fileType:       type of file (see TestData.h)
//          
// Output:                      YES if read OK
// ############################# Public Method ###############################
LOGICAL FFFBFilter::readBCoeffsFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType))
    return NO;
//
// Now, set the coefficients:
//
  setFilterOrder(_numberFileData);
  setFilterBCoeffs(_fileData);

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
LOGICAL FFFBFilter::readPolesFromFile(const char *fileName, int fileType)
{
//
// Read in the file, data is stored in _fileData and _numberFileData
//
  if(!readDataFile(fileName, fileType))
    return NO;
//
// Now, set the poles:
//
  setFilterOrder(_numberFileData);
  setFilterPoles(_fileData);
  return YES;
}

// ############################# Public Method ###############################
// setFilterOrder -- Sets the order of the filter
// Input:   order:              New filter order
//          
// Output:                      None
// ############################# Public Method ###############################
void FFFBFilter::setFilterOrder(int order)
{
  _filterOrder  = MAX(1, order);
  _filterOrder  = MIN(_filterOrder, MAX_FILTER_ORDER);
  checkShiftSize(_filterOrder+1);
  
  return;
}

// ############################# Public Method ###############################
// filterFloatData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex FFFBFilter::filterFloatData(float input)
{
  Complex       yout;
//
// Filter the data:
//
  yout      = filterNextDouble((double)input);
  return yout;
}

// ############################# Public Method ###############################
// filterDoubleData -- Filter an input data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   input:          An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex FFFBFilter::filterDoubleData(double input)
{
  Complex       yout;
//
// Filter the data:
//
  yout      = filterNextDouble(input);
  return yout;
}

// ############################# Public Method ###############################
// filterComplexArray -- Filter the input array and overwrite it with output.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   x:              An array of floating point data, modified on output
//          numberPts:      Number of points in the input array.
//          
// Output:                  YES if went OK, No if _aPolyObject doesn't exist.
// ############################# Public Method ###############################
LOGICAL FFFBFilter::filterComplexArray(Complex *x, int numberPts)
{
  int       j;

  for(j=0; j<numberPts; j++)
     x[j]   = filterNextComplex(x[j]);
  return YES;
}

// ############################# Public Method ###############################
// filterComplexData -- Filter a complex data sample, and return the sampled value.
//                     It uses the direct form II network, p. 151, Oppenheim &
//                     Schafer.
// Input:   xin:            An input data sample
//          
// Output:                  The filtered data sample
// ############################# Public Method ###############################
Complex  FFFBFilter::filterComplexData(Complex &input)
{
  Complex   yout;

//
// Filter the data:
//
  yout      = filterNextComplex(input);
  return yout;
}

// ############################# Public Method ###############################
// zeroOutTaps -- checks the sizer of the shift register, and zeros out the
//                digital filter's taps.
// Input:                       None
//          
// Output:                      None
// ############################# Public Method ###############################
void FFFBFilter::zeroOutTaps()
{
  int i;
//
// Check if we are big enough:
//
  if(_shiftSize < (_filterOrder+1))
  {
    delete [] _inputShift;
    delete [] _backShift;

    _shiftSize          = _filterOrder+1;
    _inputShift         = new Complex[_shiftSize+1];
    _backShift          = new Complex[_shiftSize+1];
  }

  for(i=0; i<_shiftSize; i++)
  {
    _inputShift[i]      = Complex(0., 0.);
    _backShift[i]       = Complex(0., 0.);
  }
  return;
}
