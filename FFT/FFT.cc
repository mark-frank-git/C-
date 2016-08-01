/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for computing fast      *
 * Fourier transforms, as well as DFTs and DCTs                         *
 *                                                                      *
 * File: FFT.cc                                                         *
 *                                                                      *
 ************************************************************************/

/**************************
 * Include files:         *
 **************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "FFT.h"


// ############################# Class Constructor #################################
// FFT -- Constructor for the FFT class
// Input:       sampling:   Sampling Frequency in Hertz
// Output:                  None
//
// ############################# Class Constructor #################################
FFT::FFT(float sampling)
{
  _samplingFreq      = sampling;
  _numberPoints      = 0;
  _operationType     = FORWARD_FFT;

  _frequencyPoints   = _twoSidedFrequency = NULL;
  _floatRealOutput   = _floatImagOutput   = NULL;
  _floatMagOutput    = NULL;

  _doubleRealOutput  = _doubleImagOutput  = NULL;
  _doubleMagOutput   = NULL;

  _complexArray      = NULL;

  _oldFrequencyPoints    = _oldComplexPoints      = _oldRealPoints         = _oldImagPoints         = 0;
  _oldFloatMagPoints     = _oldDoubleRealPoints   = _oldDoubleImagPoints   = _oldDoubleMagPoints    = 0;
  _oldTwoSidedPoints     = 0;

  return;
}


// ############################# Class Destructor #################################
// FFT -- Destructor for the FFT class
// Input:               None
// Output:              None
//
// ############################# Class Destructor #################################
FFT::~FFT()
{
  delete [] _frequencyPoints;
  delete [] _twoSidedFrequency;
  delete [] _floatRealOutput;
  delete [] _floatImagOutput;
  delete [] _floatMagOutput;
  delete [] _doubleRealOutput;
  delete [] _doubleMagOutput;
  delete [] _complexArray;

  return;
}

// ############################ Public Function ###################################
// setSamplingFrequency - Set the new sampling freq.
//
// Input:       frequency:  Sampling frequency in Hertz
// Output:                  None
//
// ############################ Public Function ###################################
void FFT::setSamplingFrequency(float frequency)
{
  _samplingFreq = frequency;
  return;
}

// ############################ Public Function ###################################
// setOperationType - Sets a new operation type: forward or inverse.
//
// Input:       type:       FORWARD_FFT or INVERSE_FFT
// Output:                  None
//
// ############################ Public Function ###################################
void FFT::setOperationType(int type)
{
  _operationType = type;
  return;
}

// ############################ Public Function ###################################
// getFrequencyPoints - Returns the frequency points at which the FFT was calculated.
//
// Input:                   None
// Output:                  The frequency points in Hertz
//
// ############################ Public Function ###################################
float *FFT::frequencyPoints()
{
  int    i;
  double frequency, frequency_step;

  if( _oldFrequencyPoints < _numberPoints)
  {
    delete [] _frequencyPoints;
    _frequencyPoints     = new float[_numberPoints];
    _oldFrequencyPoints  = _numberPoints;
  }
  frequency_step    = 0.;
  if(_numberPoints > 0)
    frequency_step = _samplingFreq/_numberPoints;
  frequency         = 0.;
  for(i=0; i<_numberPoints; i++)
  {
    _frequencyPoints[i]  = frequency;
    frequency           += frequency_step;
  }
  return _frequencyPoints;  
} 

// ############################ Public Function ###################################
// getTwoSidedFrequency - Returns the two-sided frequency points at which the FFT was calculated.
//
// Input:                   None
// Output:                  The frequency points in Hertz
//
// ############################ Public Function ###################################
float *FFT::twoSidedFrequency()
{
  int    i;
  double frequency, frequency_step;

  if(_oldTwoSidedPoints < _numberPoints)
  {
    delete [] _twoSidedFrequency;
    _twoSidedFrequency       = new float[_numberPoints];
    _oldTwoSidedPoints       = _numberPoints;
  }
  frequency_step = 0.;
  if(_numberPoints > 0)
    frequency_step          = _samplingFreq/_numberPoints;
  frequency                 = -_numberPoints*frequency_step/2.;
  for(i=0; i<_numberPoints; i++)
  {
    _twoSidedFrequency[i]    = frequency;
    frequency               += frequency_step;
  }
  return _twoSidedFrequency;  
} 

// ############################ Public Function ###################################
// deltaFrequency - Returns the spacing between frequency pts at which the FFT was calculated.
//
// Input:       number:         The number of points in FFT
//
// Output:                      The delta frequency in Hertz
//
// ############################ Public Function ###################################
float FFT::deltaFrequency(int number)
{
  float frequency_step;

  frequency_step        = 0.;
  if(number > 0)
    frequency_step      = _samplingFreq/number;
  else if(_numberPoints > 0)
    frequency_step      = _samplingFreq/_numberPoints;
  return frequency_step;  
}

// ############################ Public Function ###################################
// isPowerOfTwo - Returns whether or not the input integer is a power of 2.
//
// Input:           n:      integer to be tested
// Output:                  YES if n is a power of 2
//
// ############################ Public Function ###################################
LOGICAL FFT:: isPowerOfTwo(int n)
{
  int old_n, two_to_n;
  _powerOfTwo = -1;

  old_n = n;
  while(1)
  {
    if(n==0)
      break;
    n >>= 1;            // Divide by 2
    _powerOfTwo++;
  }
  if(_powerOfTwo < 1)
    return NO;
  two_to_n  = 1<<_powerOfTwo;
  if(old_n == two_to_n)
    return YES;
  return NO;
}

// ############################ Public Function ###################################
// fft - Calculates the FFT of the input arrays.
//
// Input:       realInput:      real input array
//              imagInput:      imaginary input array
//              n:              number of points in the above array.
//
// Output:                      None
//
// Notes:
// 1. The number of points is not required to be power of two.  If not a power of 2,
//    the DFT is used, see the function isPowerOfTwo() for checking.
// 2. To get the output of the FFT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT:: fft(const float *realInput, const float *imagInput, int n)
{
  int i;
  Complex *dft_input;

  _numberPoints  = n;
  if(_oldComplexPoints   < n)
  {
    _oldComplexPoints    = n;
    delete [] _complexArray;
    _complexArray        = new Complex[n];
  }

  if(isPowerOfTwo(n))
  {
    for(i=0; i<_numberPoints; i++)
      _complexArray[i]   = Complex(realInput[i], imagInput[i]);
    complex_fft(_complexArray, _powerOfTwo, _operationType);
  }
  else
  {
    dft_input   = new Complex[_numberPoints];
    for(i=0; i<_numberPoints; i++)
      dft_input[i]  = Complex(realInput[i], imagInput[i]);
    complex_dft(dft_input, _complexArray, n, _operationType);
    delete [] dft_input;
  }
  return;
}

// ############################ Public Function ###################################
// fft - Calculates the FFT of the input arrays.
//
// Input:       realInput:      real input array
//              imagInput:      imaginary input array
//              n:              number of points in the above array.
//
// Output:                      None
//
// Notes:
// 1. The number of points is not required to be power of two.  If not a power of 2,
//    the DFT is used, see the function isPowerOfTwo() for checking.
// 2. To get the output of the FFT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT::fft(const double *realInput, const double *imagInput, int n)
{
  int i;
  Complex *dft_input;

  _numberPoints  = n;
  if(_oldComplexPoints   < n)
  {
    delete [] _complexArray;
    _oldComplexPoints    = n;
    _complexArray        = new Complex[n];
  }

  if(isPowerOfTwo(n))
  {
    for(i=0; i<_numberPoints; i++)
      _complexArray[i]   = Complex(realInput[i], imagInput[i]);
    complex_fft(_complexArray, _powerOfTwo, _operationType);
  }
  else
  {
    dft_input   = new Complex[_numberPoints];
    for(i=0; i<_numberPoints; i++)
      dft_input[i]  = Complex(realInput[i], imagInput[i]);
    complex_dft(dft_input, _complexArray, n, _operationType);
    delete [] dft_input;
  }
  return;
}

// ############################ Public Function ###################################
// fft - Calculates the FFT of the input complex array.
//
// Input:       complexInput:   complex input array
//              n:              number of points in the above array.
//
// Output:                      None
//
// Notes:
// 1. The number of points is not required to be power of two.  If not a power of 2,
//    the DFT is used, see the function isPowerOfTwo() for checking.
// 2. To get the output of the FFT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT:: fft(const Complex *complexInput, int n)
{
  int i;
  Complex *dft_input;

  _numberPoints  = n;
  if(_oldComplexPoints   < n)
  {
    delete [] _complexArray;
    _oldComplexPoints    = n;
    _complexArray        = new Complex[n];
  }

  if(isPowerOfTwo(n))
  {
    for(i=0; i<_numberPoints; i++)
      _complexArray[i]   = complexInput[i];
    complex_fft(_complexArray, _powerOfTwo, _operationType);
  }
  else
  {
    dft_input   = new Complex[_numberPoints];
    for(i=0; i<_numberPoints; i++)
      dft_input[i]  = complexInput[i];
    complex_dft(dft_input, _complexArray, n, _operationType);
    delete [] dft_input;
  }
  return;
}


// ############################ Public Function ###################################
// czt - Calculates the chirp z transform of the real and imaginary input arrays.
//
// Input:       realInput:      real input array
//              imagInput       imag input array
//              inputPoints:    number of points in the above arrays.
//              omegaStart:     starting point of evaluation [0, 2*PI]
//              omegaEnd        ending point of evaluation [omega_start, 2*PI]
//              numberOutputPoints # of output Chirp Z points
//
// Output:                      None
//
// Notes:
// 1. To get the output of the CZT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT:: czt(const float *realInput, const float *imagInput, int inputPoints, float omegaStart,
               float omegaEnd, int numberOutputPoints)
{

  _numberPoints         = numberOutputPoints;
  if(_oldComplexPoints   < numberOutputPoints)
  {
    _oldComplexPoints    = numberOutputPoints;
    delete [] _complexArray;
    _complexArray        = new Complex[numberOutputPoints];
  }
//
// Send input to chirp z:
//
  complex_czt(realInput, imagInput, inputPoints, omegaStart, omegaEnd, _complexArray, numberOutputPoints);

  return;
}

// ############################ Public Function ###################################
// realFFT - Calculates the FFT of the real input array.
//
// Input:       realInput:      real input array
//              n:              number of points in the above array.
//
// Output:                      None
//
// Notes:
// 1. The number of points is not required to be power of two.  If not a power of 2,
//    the DFT is used, see the function isPowerOfTwo() for checking.
// 2. To get the output of the FFT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT:: realFFT(const float *realInput, int n)
{
  int i;
  Complex *dft_input;

  _numberPoints  = n;
  if(_oldComplexPoints   < n)
  {
    delete [] _complexArray;
    _oldComplexPoints    = n;
    _complexArray        = new Complex[n];
  }
  if(isPowerOfTwo(n))
  {
    for(i=0; i<_numberPoints; i++)
      _complexArray[i]   = realInput[i];
    complex_fft(_complexArray, _powerOfTwo, _operationType);
  }
  else
  {
    dft_input   = new Complex[_numberPoints];
    for(i=0; i<_numberPoints; i++)
      dft_input[i]  = realInput[i];
    complex_dft(dft_input, _complexArray, n, _operationType);
    delete [] dft_input;
  }
  return;
}

// ############################ Public Function ###################################
// realFFT - Calculates the FFT of the real input array.
//
// Input:       realInput:      real input array
//              n:              number of points in the above array.
//
// Output:                      None
//
// Notes:
// 1. The number of points is not required to be power of two.  If not a power of 2,
//    the DFT is used, see the function isPowerOfTwo() for checking.
// 2. To get the output of the FFT calculation, call realFloatOutput(), etc.
// 3. This hasn't yet been optimized to use faster real_fft() routine.
// ############################ Public Function ###################################
void FFT:: realFFT(const double *realInput, int n)
{
  int i;
  Complex *dft_input;

  _numberPoints  = n;
  if(_oldComplexPoints   < n)
  {
    delete [] _complexArray;
    _oldComplexPoints    = n;
    _complexArray        = new Complex[n];
  }

  if(isPowerOfTwo(n))
  {
    for(i=0; i<_numberPoints; i++)
      _complexArray[i]   = realInput[i];
    complex_fft(_complexArray, _powerOfTwo, _operationType);
  }
  else
  {
    dft_input   = new Complex[_numberPoints];
    for(i=0; i<_numberPoints; i++)
      dft_input[i]  = realInput[i];
    complex_dft(dft_input, _complexArray, n, _operationType);
    delete [] dft_input;
  }
  return;
}

// ############################ Public Function ###################################
// realCZT - Calculates the chirp z transform of the real input array.
//
// Input:       realInput:      real input array
//              n:              number of points in the above array.
//              omegaStart:     starting point of evaluation [0, 2*PI]
//              omegaEnd        ending point of evaluation [omega_start, 2*PI]
//              numberOutputPoints # of output Chirp Z points
//
// Output:                      None
//
// Notes:
// 1. To get the output of the CZT calculation, call realFloatOutput(), etc.
// ############################ Public Function ###################################
void FFT:: realCZT(const float *realInput, int numberInputPoints, float omegaStart, float omegaEnd,
                   int numberOutputPoints)
{
  _numberPoints          = numberOutputPoints;
  if(_oldComplexPoints   < numberOutputPoints)
  {
    delete [] _complexArray;
    _oldComplexPoints    = numberOutputPoints;
    _complexArray        = new Complex[numberOutputPoints];
  }

  real_czt(realInput, numberInputPoints, omegaStart, omegaEnd, _complexArray, numberOutputPoints);
  return;
}

// ############################ Public Function ###################################
// realFloatOutput - Returns the real part of a previous FFT calculation.
//
// Input:                       None
//
// Output:                      The real part of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const float *FFT:: realFloatOutput()
{
  int i;

  if(_oldRealPoints  < _numberPoints)
  {
    delete [] _floatRealOutput;
    _oldRealPoints   = _numberPoints;
    _floatRealOutput = new float[_numberPoints];
  }
  for(i=0; i<_numberPoints; i++)
    _floatRealOutput[i] = real(_complexArray[i]);
  return _floatRealOutput;
}

// ############################ Public Function ###################################
// imagFloatOutput - Returns the imaginary part of a previous FFT calculation.
//
// Input:                       None
//
// Output:                      The imaginary part of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const float *FFT:: imagFloatOutput()
{
  int i;

  if(_oldImagPoints  < _numberPoints)
  {
    delete [] _floatImagOutput;
    _oldImagPoints   = _numberPoints;
    _floatImagOutput = new float[_numberPoints];
  }
  for(i=0; i<_numberPoints; i++)
    _floatImagOutput[i] = imag(_complexArray[i]);
  return _floatImagOutput;
}

// ############################ Public Function ###################################
// magFloatOutput - Returns the magnitude of a previous FFT calculation.
//
// Input:       magType:        MAGNITUDE_IN_DB, MAGNITUDE_IN_WATTS, COMPLEX_PHASE, etc
//
// Output:                      The magnitude of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const float *FFT:: magFloatOutput(int magType)
{
  int i;
  double magnitude;

  if(_oldFloatMagPoints  < _numberPoints)
  {
    delete [] _floatMagOutput;
    _oldFloatMagPoints   = _numberPoints;
    _floatMagOutput      = new float[_numberPoints];
  }
  switch(magType)
  {
   case MAGNITUDE_IN_DB:
    for(i=0; i<_numberPoints; i++)
    {
      magnitude = norm(_complexArray[i]);
      if(magnitude>0.)
        _floatMagOutput[i]   = 10.*log10(magnitude);
      else
        _floatMagOutput[i]   = 0.;
    }
    break;
   case MAGNITUDE_IN_WATTS:
   default:
    for(i=0; i<_numberPoints; i++)
      _floatMagOutput[i] = norm(_complexArray[i]);
    break;
   case MAGNITUDE_IN_VOLTS:
    for(i=0; i<_numberPoints; i++)
      _floatMagOutput[i] = abs(_complexArray[i]);
    break;
   case COMPLEX_PHASE:
    for(i=0; i<_numberPoints; i++)
      _floatMagOutput[i] = arg(_complexArray[i]);
    break;
   case COMPLEX_REAL_PART:
    for(i=0; i<_numberPoints; i++)
      _floatMagOutput[i] = real(_complexArray[i]);
    break;
   case COMPLEX_IMAG_PART:
    for(i=0; i<_numberPoints; i++)
      _floatMagOutput[i] = imag(_complexArray[i]);
    break;
  }
  return _floatMagOutput;
}

// ############################ Public Function ###################################
// realDoubleOutput - Returns the real part of a previous FFT calculation.
//
// Input:                       None
//
// Output:                      The real part of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const double *FFT:: realDoubleOutput()
{
  int i;

  if(_oldDoubleRealPoints    < _numberPoints)
  {
    delete [] _doubleRealOutput;
    _oldDoubleRealPoints = _numberPoints;
    _doubleRealOutput    = new double[_numberPoints];
  }
  for(i=0; i<_numberPoints; i++)
    _doubleRealOutput[i] = real(_complexArray[i]);
  return _doubleRealOutput;
}

// ############################ Public Function ###################################
// imagDoubleOutput - Returns the imaginary part of a previous FFT calculation.
//
// Input:                       None
//
// Output:                      The imaginary part of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const double *FFT:: imagDoubleOutput()
{
  int i;

  if(_oldDoubleImagPoints    < _numberPoints)
  {
    delete [] _doubleImagOutput;
    _oldDoubleImagPoints = _numberPoints;
    _doubleImagOutput    = new double[_numberPoints];
  }
  for(i=0; i<_numberPoints; i++)
    _doubleImagOutput[i] = imag(_complexArray[i]);
  return _doubleImagOutput;
}

// ############################ Public Function ###################################
// magDoubleOutput - Returns the magnitude of a previous FFT calculation.
//
// Input:       magType:        MAGNITUDE_IN_DB, MAGNITUDE_IN_WATTS, COMPLEX_REAL_PART, etc.
//
// Output:                      The magnitude of the FFT output
//
// Notes:
// 1. You must first call one of the FFT routines before getting the output from here.
// ############################ Public Function ###################################
const double *FFT:: magDoubleOutput(int magType)
{
  int i;
  double magnitude;

  if(_oldDoubleMagPoints < _numberPoints)
  {
    delete [] _doubleMagOutput;
    _oldDoubleMagPoints  = _numberPoints;
    _doubleMagOutput     = new double[_numberPoints];
  }
  switch(magType)
  {
   case MAGNITUDE_IN_DB:
    for(i=0; i<_numberPoints; i++)
    {
      magnitude = abs(_complexArray[i]);
      if(magnitude>0.)
        _doubleMagOutput[i]  = 10.*log10(magnitude);
      else
        _doubleMagOutput[i]  = 0.;
    }
    break;
   case MAGNITUDE_IN_WATTS:
   default:
    for(i=0; i<_numberPoints; i++)
      _doubleMagOutput[i]    = abs(_complexArray[i]);
    break;
   case MAGNITUDE_IN_VOLTS:
    for(i=0; i<_numberPoints; i++)
      _doubleMagOutput[i]    = norm(_complexArray[i]);
    break;
   case COMPLEX_PHASE:
    for(i=0; i<_numberPoints; i++)
      _doubleMagOutput[i]    = arg(_complexArray[i]);
    break;
   case COMPLEX_REAL_PART:
    for(i=0; i<_numberPoints; i++)
      _doubleMagOutput[i]    = real(_complexArray[i]);
    break;
   case COMPLEX_IMAG_PART:
    for(i=0; i<_numberPoints; i++)
      _doubleMagOutput[i]    = imag(_complexArray[i]);
    break;
  }
  return _doubleMagOutput;
}

