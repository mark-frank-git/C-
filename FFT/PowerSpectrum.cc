/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for computing PSDs,     *
 * autocorrelations, etc.                                               *
 *                                                                      *
 * File: PowerSpectrum.cc                                               *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 10/06/95  - Started                                              *
 *  2. 12/28/95  - Fixed error using atan2(), in autocorrelation.       *
 *  3. 02/07/97  - Add avgPowerSpectrumFromReal().                      *
 *  4. 02/11/97  - Add smoothSpectrum()                                 *
 *  5. 02/13/97  - Add shiftRight for Complex data.                     *
 *  6. 04/07/97  - Add calibrationMode.                                 *
 *  7. 04/29/97  - Move windowing stuff to _dataWindow class.           *
 *  8. 04/30/97  - Add resolutionBandwidth(), numberDataSamples().      *
 *  9. 05/01/97  - Add chirp z for avgPowerSpectrum().                  *
 * 10. 05/19/97  - Add windowSize to avgPowerSpectrumFromReal().        *
 * 11. 05/28/97  - Add maxBinFrequency().                               *
 * 12. 08/25/99  - Compensate scaling for zero padding, and coherent    *
 *                 versus non-coherent gain of window.                  *
 * 13. 09/02/99  - Fixed power spectrum scale factor -> 1/N.            *
 * 14. 09/10/99  - Use 3 dB BW in RBW calculation.                      *
 * 15. 09/22/99  - Add averagingType in avgPowerSpectrum().             *
 * 16. 11/06/00  - Add averagingType to avgPowerSpectrum with Chirp Z.  *
 ************************************************************************/

/**************************
 * Include files:         *
 **************************/
#include "FFT.h"
#include "PowerSpectrum.h"
#include "fft_routines.h"
#if defined (_WIN32) && defined (NeXT_PDO)
#include <Filters/DataWindow.h>
#else
#include "DataWindow.h"
#endif
#include <math.h>
#include <stdlib.h>


#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define TWOPI   6.283185307                     // 2*PI

// ############################# Class Constructor #################################
// PowerSpectrum -- Constructor for the PowerSpectrum class
//
// Input:       sampling:   Sampling Frequency in Hertz
//
// Output:                  None
//
// ############################# Class Constructor #################################
PowerSpectrum::PowerSpectrum(float sampling)
{
  _fft                   = new FFT(sampling);
  _dataWindow            = new DataWindow();
  _numberSpectrumPoints  = _numberAutoPoints  = _numberFloatPoints = 0;
  _numberPhasePoints    = 0;
  _oldSpectrumPoints     = _oldAutoPoints     = _oldPhasePoints = 0;
  _meanIndexType         = MEAN_FROM_SCM;
  _calibrateMode         = CALIBRATE_FOR_NOISE;
  _twoSidedBandwidth     = 0.;
  _outputSpectrum        = _autocorrelation   = _outputFloat   = NULL;
  _outputPhase          = NULL;


  return;
}

// ############################# Class Destructor #################################
// ~PowerSpectrum -- Destructor for the PowerSpectrum class
//
// Input:       sampling:   Sampling Frequency in Hertz
//
// Output:                  None
//
// ############################# Class Destructor #################################
PowerSpectrum::~PowerSpectrum()
{
  delete    _fft;
  delete    _dataWindow;
  delete [] _outputSpectrum;
  delete [] _autocorrelation;
  delete [] _outputFloat;
  delete [] _outputPhase;
  return;
}

// ############################# Public Method ###############################
// setWindowType -- Set a new windowing type for the data window.
//
// Input:   type:           E.g., Hamming
//          
// Output:                  None
//
// Notes:
// ############################# Public Method ###############################
void PowerSpectrum::setWindowType(int type)
{
  _dataWindow->setWindowType(type);
  return;
}

// ############################# Public Method ###############################
// resolutionBandwidth -- Returns the resolution BW of the PSD given the number
//                        of data samples.
//
// Input:   numberSamples:  Number of data samples
//          
// Output:                  Resolution BW in Hertz
//
// Notes:
// 1. The resolution bandwidth is a function of the sampling frequency, the
//    number of data points, and the window function type.
// 2. Uses 3 dB RBW.
// ############################# Public Method ###############################
float PowerSpectrum::resolutionBandwidth(int numberSamples)
{
  float res_bw = 0.;

  if(numberSamples > 0)
    res_bw      = _fft->samplingFreq()/numberSamples;
  res_bw        *= _dataWindow->bandwidth3dB();
  return  res_bw;
}

// ############################# Public Method ###############################
// frequencyFromIndex -- Returns the frequency corresponding to FFT array index.
//
// Input:   index:              (interpolated) array index
//          
// Output:                      Corresponding frequency in Hertz
//
// Notes:
// 1. The index is a floating point number to account for interpolated index
// ############################# Public Method ###############################
float PowerSpectrum::frequencyFromIndex(float index)
{
  return  index*_fft->deltaFrequency();
}

// ############################# Public Method ###############################
// numberDataSamples -- Returns the number of data samples needed to achieve
//                      the given resolution bw
//
// Input:   resolutionBandwidth:    PSD resolution BW in Hertz
//          
// Output:                          number of data samples needed
//
// Notes:
// 1. The resolution bandwidth is a function of the sampling frequency, the
//    number of data points, and the window function type.
// 2. Uses 3 dB RBW.
// ############################# Public Method ###############################
int PowerSpectrum::numberDataSamples(float resolutionBandwidth)
{
  int   number_samples  = 1;

  if(resolutionBandwidth > 0.)
    number_samples  = ROUND(_dataWindow->bandwidth3dB()*_fft->samplingFreq()/resolutionBandwidth);
  return  number_samples;
}

// ############################# Public Method ###############################
// nextFFTSize -- Returns the next power of 2 greater than or equal to the input
//                number of samples.
//
// Input:   numberSamples:  The number of data samples
//          
// Output:                  The next higher FFT size
//
// Notes:
// ############################# Public Method ###############################
int PowerSpectrum::nextFFTSize(int numberSamples)
{
  int   size    = 2;

  while(size < numberSamples)
    size    <<= 1;                  // Multiply by 2
  return  size;
}

// ############################# Public Method ###############################
// indexFromFrequency -- Returns the FFT array index corresponding to the input
//                      frequency
//
// Input:   frequency:          frequency in Hertz
//          
// Output:                      index of FFT array
//
// Notes:
// ############################# Public Method ###############################
int PowerSpectrum::indexFromFrequency(float frequency)
{
  int   index = 0;
  float delta;
  delta = _fft->deltaFrequency();
  if(delta > 0.)
    index       = ROUND((frequency/delta));
  return  index;
}

// ############################# Public Method ###############################
// windowType -- Gets the current windowing type for the data window.
//
// Input:                   None
//          
// Output:                  window type, e.g., Hamming
//
// Notes:
// ############################# Public Method ###############################
int PowerSpectrum::windowType()
{
  return _dataWindow->windowType();
}

// ############################# Public Method #################################
// powerSpectrum -- Calculates the power spectrum of real input data.
//
// Input:       input:          Real time series data
//              numberPoints:   Number of points in above array
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// ############################# Public Method #################################
float *PowerSpectrum::powerSpectrum(float *input, int numberPoints)
{
  int   i;
  float scale;
  const float *fft_magnitude;
//
// Allocate space
//
  _numberSpectrumPoints  = numberPoints;
  if(_oldSpectrumPoints < numberPoints)
  {
    _oldSpectrumPoints   = numberPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[numberPoints];
  }
//
// First find the FFT, and then scale.
//
  _fft->realFFT(input, numberPoints);
  fft_magnitude = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
  scale = 1.;
//
// Calculate scale factor
//
  if(numberPoints)
  {
    scale               = 1./numberPoints;
    if(_calibrateMode == CALIBRATE_FOR_NOISE && (_fft->samplingFreq() > 0.))
      scale             /= _fft->samplingFreq();        // given in dB/Hz
  }
//
// Scale FFT output to get power spectrum
//
  for(i=0; i<numberPoints; i++)
    _outputSpectrum[i] = fft_magnitude[i]*scale;
  return _outputSpectrum;
}

// ############################# Public Method #################################
// psdFromAuto -- Calculates the power spectrum from autocorrelation data.
//
// Input:       input:          Real autocorrelation data
//              numberPoints:   Number of points in above array
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. The power spectrum data is given in Watts/Hz.
// ############################# Public Method #################################
float *PowerSpectrum::psdFromAuto(float *input, int number)
{
  int       n, k;
  const     double *window;
  static    float *coefficients = NULL;
  static    int old_coeff_points;
  double    arg, window_scale;
//
// Allocate space
//
  _numberSpectrumPoints  = 2*number-1;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints   = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[_numberSpectrumPoints];
  }
  if(old_coeff_points != _numberSpectrumPoints)
  {
    old_coeff_points = _numberSpectrumPoints;
    delete [] coefficients;
    coefficients    = new float[_numberSpectrumPoints];
    arg             = TWOPI/_numberSpectrumPoints;
    for(n=0; n<_numberSpectrumPoints; n++)
       coefficients[n] = cos(arg*n);                    // Used for taking Inverse DFT
  }
  window        = _dataWindow->windowFunction(_numberSpectrumPoints);
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
  for(k=0; k<_numberSpectrumPoints; k++)
  {
    _outputSpectrum[k] = window[number-1]*input[0];
    for(n=1; n<number; n++)
      _outputSpectrum[k] += window[number-1+n]*input[n]*
        (coefficients[ (k*n)%_numberSpectrumPoints] +
         coefficients[ (k*(_numberSpectrumPoints-1))%_numberSpectrumPoints])/window_scale;
  }
  return _outputSpectrum;
}

// ############################# Public Method #################################
// avgPowerSpectrum -- Calculates the average power spectrum of complex input data.
//
// Input:       realInput:      Real part of time series data
//              imagInput:      Imaginary part of time series data
//              fftSize:        Size of each FFT to take
//              overlapSize:    Size of overlap
//              numberAvgs:     Number of FFTs to average
//              type:           type of spectrum, e.g., MAGNITUDE_IN_WATTS
//              windowSize:     Size of data to window, zero pad up to fftSize.  The
//                              default value is zero.  If windowSize=0, no zero pad.
//              initialize:     YES = initialize avg to zero, else add current FFT to running avg
//              averageType:    PSD_AVG, PSD_MIN, PSD_MAX
//
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// 3. The size of the input arrays have to be (fftSize-overlapSize)*numberAvgs
// 4. The time series data is windowed prior to each FFT, see windowType.
// ############################# Public Method #################################
float *PowerSpectrum::avgPowerSpectrum(const float *realInput, const float *imagInput, int fftSize,
                                       int overlapSize, int numberAvgs, int type, int windowSize,
                                       LOGICAL initialize, int averageType)
{
  int           i, n, shift;
  float         scale;
  const float   *fft_magnitude;
  const float   *real_ptr, *imag_ptr;
  float         *windowed_real, *windowed_imag;
  const         double *window;
  double        window_scale;
//
// Allocate space
//
  _numberSpectrumPoints  = fftSize;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints   = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[_numberSpectrumPoints];
  }
//
// Calculate shift from overlap and FFT size, check for errors
//
  shift         = fftSize - overlapSize;
  shift         = MAX(0, shift);
  shift         = MIN(fftSize, shift);
//
// Get data window
//
  windowSize    = MIN(windowSize, _numberSpectrumPoints);               // error check
  if(windowSize < 1)
    windowSize  = _numberSpectrumPoints;
  windowed_real = new float[_numberSpectrumPoints];
  windowed_imag = new float[_numberSpectrumPoints];
  window        = _dataWindow->windowFunction(windowSize);
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
//
// Initialize the averages
//
  if(initialize)
  {
    for(i=0; i<_numberSpectrumPoints; i++)
      _outputSpectrum[i] = 0.;
  }
  real_ptr  = realInput;
  imag_ptr  = imagInput;
//
// Set up scale factor for resolution bw, window, and zero padding:
//
  scale                 = 1./window_scale;
  if( (windowSize>0) && (numberAvgs>0) )
  {
    scale               /= (float)numberAvgs*windowSize;
    if(_calibrateMode == CALIBRATE_FOR_NOISE && (_fft->samplingFreq() > 0.))
      scale             /= _fft->samplingFreq();        // given in dB/Hz
  }
  for(n=0; n<numberAvgs; n++)
  {
    for(i=0; i<windowSize; i++)                         // Window the input arrays
    {
      windowed_real[i]  = window[i]*real_ptr[i];
      windowed_imag[i]  = window[i]*imag_ptr[i];
    }
    for(; i<_numberSpectrumPoints; i++)                 // Zero pad if window size < fft size
      windowed_real[i]  = windowed_imag[i]      = 0.;
    _fft->fft(windowed_real, windowed_imag, _numberSpectrumPoints);
    fft_magnitude       = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
    switch(averageType)
    {
      case PSD_AVG:
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    += fft_magnitude[i]*scale;
        break;
      case PSD_MIN:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
        {
          if( _outputSpectrum[i] > 0.)
            _outputSpectrum[i]  = MIN(_outputSpectrum[i], (fft_magnitude[i]*scale));
          else
            _outputSpectrum[i]  = fft_magnitude[i]*scale;
        }
        break;
      case PSD_MAX:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    = MAX(_outputSpectrum[i], (fft_magnitude[i]*scale));
        break;
    }
    real_ptr += shift;
    imag_ptr += shift;
  }
  switch(type)
  {
    default:
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i]>0.)
          _outputSpectrum[i] = 10.*log10(_outputSpectrum[i]);
        else
          _outputSpectrum[i] = 0.;
      }
      break;
   case MAGNITUDE_IN_VOLTS:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i] > 0.)
          _outputSpectrum[i] = sqrt(_outputSpectrum[i]);
      }
      break;
  }
  delete [] windowed_real;
  delete [] windowed_imag;
  return _outputSpectrum;
}

// ############################# Public Method #################################
// avgPowerSpectrum -- Calculates the average power spectrum of complex input data.
//
// Input:       complexInput:   Real and imaginary parts of time series data
//              fftSize:        Size of each FFT to take
//              overlapSize:    Size of overlap
//              numberAvgs:     Number of FFTs to average
//              type:           type of spectrum, e.g., MAGNITUDE_IN_WATTS
//              windowSize:     Size of data to window, zero pad up to fftSize.  The
//                              default value is zero.  If windowSize=0, no zero pad.
//              initialize:     YES = initialize avg to zero, else add current FFT to running avg
//              averageType:    PSD_AVG, PSD_MIN, PSD_MAX
//
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// 3. The size of the input arrays have to be (fftSize-overlapSize)*numberAvgs
// 4. The time series data is windowed prior to each FFT, see windowType.
// ############################# Public Method #################################
float *PowerSpectrum::avgPowerSpectrum(const Complex *complexInput, int fftSize,
                                       int overlapSize, int numberAvgs, int type, int windowSize,
                                       LOGICAL initialize, int averageType)
{
  int           i, n, shift;
  float         scale;
  const float   *fft_magnitude;
  const Complex *input_ptr;
  float         *windowed_real, *windowed_imag;
  const         double *window;
  double        window_scale;
//
// Allocate space
//
  _numberSpectrumPoints  = fftSize;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints   = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[_numberSpectrumPoints];
  }
//
// Calculate shift from overlap and FFT size, check for errors
//
  shift         = fftSize - overlapSize;
  shift         = MAX(0, shift);
  shift         = MIN(fftSize, shift);
//
// Get data window
//
  windowSize    = MIN(windowSize, _numberSpectrumPoints);               // error check
  if(windowSize < 1)
    windowSize  = _numberSpectrumPoints;
  windowed_real = new float[_numberSpectrumPoints];
  windowed_imag = new float[_numberSpectrumPoints];
  window        = _dataWindow->windowFunction(windowSize);
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
//
// Initialize the averages
//
  if(initialize)
  {
    for(i=0; i<_numberSpectrumPoints; i++)
      _outputSpectrum[i] = 0.;
  }
  input_ptr             = complexInput;
//
// Set up scale factor for resolution bw, window, and zero padding:
//
  scale                 = 1./window_scale;
  if( (windowSize>0) && (numberAvgs>0) )
  {
    scale               /= (float)numberAvgs*windowSize;
    if(_calibrateMode == CALIBRATE_FOR_NOISE && (_fft->samplingFreq() > 0.))
      scale             /= _fft->samplingFreq();        // given in dB/Hz
  }
  for(n=0; n<numberAvgs; n++)
  {
    for(i=0; i<windowSize; i++)                         // Window the input arrays
    {
      Complex windowed  = window[i]*input_ptr[i];
      windowed_real[i]  = real(windowed);
      windowed_imag[i]  = imag(windowed);
    }
    for(; i<_numberSpectrumPoints; i++)                 // Zero pad if window size < fft size
      windowed_real[i]  = windowed_imag[i]      = 0.;
    _fft->fft(windowed_real, windowed_imag, _numberSpectrumPoints);
    fft_magnitude       = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
    switch(averageType)
    {
      case PSD_AVG:
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    += fft_magnitude[i]*scale;
        break;
      case PSD_MIN:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
        {
          if( _outputSpectrum[i] > 0.)
            _outputSpectrum[i]  = MIN(_outputSpectrum[i], (fft_magnitude[i]*scale));
          else
            _outputSpectrum[i]  = fft_magnitude[i]*scale;
        }
        break;
      case PSD_MAX:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    = MAX(_outputSpectrum[i], (fft_magnitude[i]*scale));
        break;
    }
    input_ptr   += shift;
  }
  switch(type)
  {
    default:
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i]>0.)
          _outputSpectrum[i] = 10.*log10(_outputSpectrum[i]);
        else
          _outputSpectrum[i] = 0.;
      }
      break;
   case MAGNITUDE_IN_VOLTS:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i] > 0.)
          _outputSpectrum[i] = sqrt(_outputSpectrum[i]);
      }
      break;
  }
  delete [] windowed_real;
  delete [] windowed_imag;
  return _outputSpectrum;
}
// ############################# Public Method #################################
// avgPowerSpectrum -- Calculates the average power spectrum of complex input
//                      data using the Chirp Z algorithm.
//
// Input:       realInput:      Real part of time series data
//              imagInput:      Imag part of time series data
//              inputSize:      Size of each data window to take
//              overlapSize:    Size of overlap
//              numberAvgs:     Number of FFTs to average
//              numberOutputPts: Number of output points for chirp z
//              startFreq:      Start frequency in Hertz for chirp z
//              endFreq:        End frequency in Hertz for chirp z
//              type:           type of spectrum, e.g., MAGNITUDE_IN_WATTS
//              initialize:     YES = initialize avg to zero, else add current FFT to running avg
//              averageType:    PSD_AVG, PSD_MIN, PSD_MAX
//
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// 3. The size of the input arrays have to be (fftSize-overlapSize)*numberAvgs
// 4. The time series data is windowed prior to each FFT, see windowType.
// ############################# Public Method #################################
float *PowerSpectrum::avgPowerSpectrum(const float *realInput, const float *imagInput, int inputSize, int overlapSize,
                                               int numberAvgs, int numberOutputPts,
                                               float startFreq, float endFreq, int type, LOGICAL initialize, int averageType)
{
  int           i, n, shift;
  float         scale;
  const float   *fft_magnitude;
  const float   *real_ptr, *imag_ptr;
  float         digital_start, digital_end;
  float         *windowed_real, *windowed_imag;
  const         double *window;
  double        window_scale;
//
// Allocate space
//
  _numberSpectrumPoints = numberOutputPts;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints  = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum     = new float[_numberSpectrumPoints];
  }
//
// Calculate shift from overlap and FFT size, check for errors
//
  shift         = inputSize - overlapSize;
  shift         = MAX(0, shift);
  shift         = MIN(inputSize, shift);
//
// Get data window
//
  window        = _dataWindow->windowFunction(inputSize);
  windowed_real = new float[inputSize];
  windowed_imag = new float[inputSize];
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
//
// Get the start and stop frequencies:
//
  digital_start = TWOPI*startFreq/_fft->samplingFreq();
  digital_end   = TWOPI*endFreq/_fft->samplingFreq();
//
// Initialize the averages
//
  if(initialize)
  {
    for(i=0; i<_numberSpectrumPoints; i++)
      _outputSpectrum[i] = 0.;
  }
  real_ptr      = realInput;
  imag_ptr      = imagInput;
  scale         = 1./window_scale;
  if( (inputSize>0) && (numberAvgs>0) )
  {
    scale               /= (float)numberAvgs*inputSize;
    if((_calibrateMode == CALIBRATE_FOR_NOISE) && (_fft->samplingFreq() > 0.))
      scale             /= _fft->samplingFreq();        // given in dB/Hz
  }
  for(n=0; n<numberAvgs; n++)
  {
    for(i=0; i<inputSize; i++)                          // Window the input arrays
    {
      windowed_real[i]  = window[i]*real_ptr[i];
      windowed_imag[i]  = window[i]*imag_ptr[i];
    }
    _fft->czt(windowed_real, windowed_imag, inputSize, digital_start, digital_end, _numberSpectrumPoints);
    fft_magnitude   = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
    switch(averageType)
    {
      case PSD_AVG:
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    += fft_magnitude[i]*scale;
        break;
      case PSD_MIN:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
        {
          if( _outputSpectrum[i] > 0.)
            _outputSpectrum[i]  = MIN(_outputSpectrum[i], (fft_magnitude[i]*scale));
          else
            _outputSpectrum[i]  = fft_magnitude[i]*scale;
        }
        break;
      case PSD_MAX:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    = MAX(_outputSpectrum[i], (fft_magnitude[i]*scale));
        break;
    }
    real_ptr    += shift;
    imag_ptr    += shift;
  }
  switch(type)
  {
    default:
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i]>0.)
          _outputSpectrum[i] = 10.*log10(_outputSpectrum[i]);
        else
          _outputSpectrum[i] = 0.;
      }
      break;
    case MAGNITUDE_IN_VOLTS:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i] > 0.)
          _outputSpectrum[i] = sqrt(_outputSpectrum[i]);
      }
      break;
  }
  delete [] windowed_real;
  delete [] windowed_imag;
  return _outputSpectrum;
}

// ############################# Public Method #################################
// avgPowerSpectrumFromReal -- Calculates the average power spectrum of real input data.
//
// Input:       realInput:      Real part of time series data
//              fftSize:        Size of each FFT to take
//              overlapSize:    Size of overlap
//              numberAvgs:     Number of FFTs to average
//              type:           type of spectrum, e.g., MAGNITUDE_IN_WATTS
//              windowSize:     Size of data to window, zero pad up to fftSize.  The
//                              default value is zero.  If windowSize=0, no zero pad.
//              initialize:     YES = initialize avg to zero, else add current FFT to running avg
//              averageType:    PSD_AVG, PSD_MIN, PSD_MAX
//
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// 3. The size of the input arrays have to be (fftSize-overlapSize)*numberAvgs
// 4. The time series data is windowed prior to each FFT, see windowType.
// ############################# Public Method #################################
float *PowerSpectrum::avgPowerSpectrumFromReal(const float *realInput, int fftSize, int overlapSize,
                                               int numberAvgs, int type, int windowSize, LOGICAL initialize,
                                               int averageType)
{
  int           i, n, shift;
  float         scale;
  const float   *fft_magnitude;
  float         *windowed_data;
  const         double *window;
  double        window_scale;
//
// Allocate space
//
  _numberSpectrumPoints  = fftSize;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints   = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[_numberSpectrumPoints];
  }
//
// Calculate shift from overlap and FFT size, check for errors
//
  shift         = fftSize - overlapSize;
  shift         = MAX(0, shift);
  shift         = MIN(fftSize, shift);
//
// Get data window
//
  windowSize    = MIN(windowSize, _numberSpectrumPoints);               // error check
  if(windowSize < 1)
    windowSize  = _numberSpectrumPoints;
  windowed_data = new float[_numberSpectrumPoints];
  window        = _dataWindow->windowFunction(windowSize);
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
//
// Initialize the averages
//
  if(initialize)
  {
    for(i=0; i<_numberSpectrumPoints; i++)
      _outputSpectrum[i] = 0.;
  }
//
// scale for res bw, window and zero padding:
//
  scale     = 1./window_scale;
  if( (windowSize>0) && (numberAvgs>0) )
  {
    scale       /= (float)numberAvgs*windowSize;
    if(_calibrateMode == CALIBRATE_FOR_NOISE && (_fft->samplingFreq() > 0.))
      scale     /= _fft->samplingFreq();        // given in dB/Hz
  }
  for(n=0; n<numberAvgs; n++)
  {
    for(i=0; i<windowSize; i++)                     // Window the input arrays
      windowed_data[i]  = window[i]*realInput[i];
    for(; i<_numberSpectrumPoints; i++)              // Zero pad if window size < fft size
      windowed_data[i]  = 0.;
    _fft->realFFT(windowed_data, _numberSpectrumPoints);
    fft_magnitude       = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
    switch(averageType)
    {
      case PSD_AVG:
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i] += fft_magnitude[i]*scale;
        break;
      case PSD_MIN:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
        {
          if( _outputSpectrum[i] > 0.)
            _outputSpectrum[i]  = MIN(_outputSpectrum[i], (fft_magnitude[i]*scale));
          else
            _outputSpectrum[i]  = fft_magnitude[i]*scale;
        }
        break;
      case PSD_MAX:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    = MAX(_outputSpectrum[i], (fft_magnitude[i]*scale));
        break;
    }
    realInput           += shift;
  }
  switch(type)
  {
    default:
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i]>0.)
          _outputSpectrum[i] = 10.*log10(_outputSpectrum[i]);
        else
          _outputSpectrum[i] = 0.;
      }
      break;
   case MAGNITUDE_IN_VOLTS:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i] > 0.)
          _outputSpectrum[i] = sqrt(_outputSpectrum[i]);
      }
      break;
  }
  delete [] windowed_data;
  return _outputSpectrum;
}

// ############################# Public Method #################################
// avgPowerSpectrumFromReal -- Calculates the average power spectrum of real input
//                              data using the Chirp Z algorithm.
//
// Input:       realInput:      Real part of time series data
//              inputSize:      Size of each data window to take
//              overlapSize:    Size of overlap
//              numberAvgs:     Number of FFTs to average
//              numberOutputPts: Number of output points for chirp z
//              startFreq:      Start frequency in Hertz for chirp z
//              endFreq:        End frequency in Hertz for chirp z
//              type:           type of spectrum, e.g., MAGNITUDE_IN_WATTS
//              initialize:     YES = initialize avg to zero, else add current FFT to running avg
//              averageType:    PSD_AVG, PSD_MIN, PSD_MAX
//
//              
//
// Output:                      Pointer to the output power spectrum data
//
// Notes:
// 1. Set _calibrateMode to get noise PSDs right, or get sinusoid magnitudes right.
// 2. The power spectrum data is given in Watts/Hz.
// 3. The size of the input arrays have to be (fftSize-overlapSize)*numberAvgs
// 4. The time series data is windowed prior to each FFT, see windowType.
// ############################# Public Method #################################
float *PowerSpectrum::avgPowerSpectrumFromReal(const float *realInput, int inputSize, int overlapSize,
                                               int numberAvgs, int numberOutputPts,
                                               float startFreq, float endFreq, int type, LOGICAL initialize,
                                               int averageType)
{
  int           i, n, shift;
  float         scale;
  float         *windowed_data;
  const float   *fft_magnitude;
  const float   *real_ptr;
  float         digital_start, digital_end;
  const         double *window;
  double        window_scale;
//
// Allocate space
//
  _numberSpectrumPoints  = numberOutputPts;
  if(_oldSpectrumPoints < _numberSpectrumPoints)
  {
    _oldSpectrumPoints   = _numberSpectrumPoints;
    delete [] _outputSpectrum;
    _outputSpectrum      = new float[_numberSpectrumPoints];
  }
//
// Calculate shift from overlap and FFT size, check for errors
//
  shift         = inputSize - overlapSize;
  shift         = MAX(0, shift);
  shift         = MIN(inputSize, shift);
//
// Get data window
//
  windowed_data = new float[inputSize];
  window        = _dataWindow->windowFunction(inputSize);
  if(_calibrateMode == CALIBRATE_FOR_NOISE)
    window_scale  = _dataWindow->noncoherentGain();
  else
    window_scale  = _dataWindow->coherentGain();
//
// Get the start and stop frequencies:
//
  digital_start = TWOPI*startFreq/_fft->samplingFreq();
  digital_end   = TWOPI*endFreq/_fft->samplingFreq();
//
// Initialize the averages
//
  if(initialize)
  {
    for(i=0; i<_numberSpectrumPoints; i++)
      _outputSpectrum[i] = 0.;
  }
  real_ptr      = realInput;
  scale         = 1./window_scale;
  if( (inputSize>0) && (numberAvgs>0) )
  {
    scale       /= (float)numberAvgs*inputSize;
    if(_calibrateMode == CALIBRATE_FOR_NOISE && (_fft->samplingFreq() > 0.))
      scale     /= _fft->samplingFreq();        // given in dB/Hz
  }
  for(n=0; n<numberAvgs; n++)
  {
    for(i=0; i<inputSize; i++)                  // Window the input arrays
    {
      windowed_data[i]  = real_ptr[i]*window[i];
    }
    _fft->realCZT(windowed_data, inputSize, digital_start, digital_end, _numberSpectrumPoints);
    fft_magnitude   = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);
    switch(averageType)
    {
      case PSD_AVG:
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i] += fft_magnitude[i]*scale;
        break;
      case PSD_MIN:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
        {
          if( _outputSpectrum[i] > 0.)
            _outputSpectrum[i]  = MIN(_outputSpectrum[i], (fft_magnitude[i]*scale));
          else
            _outputSpectrum[i]  = fft_magnitude[i]*scale;
        }
        break;
      case PSD_MAX:
        scale                   *= numberAvgs;
        for(i=0; i<_numberSpectrumPoints; i++)
          _outputSpectrum[i]    = MAX(_outputSpectrum[i], (fft_magnitude[i]*scale));
        break;
    }
    real_ptr += shift;
  }
  switch(type)
  {
    default:
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i]>0.)
          _outputSpectrum[i] = 10.*log10(_outputSpectrum[i]);
        else
          _outputSpectrum[i] = 0.;
      }
      break;
   case MAGNITUDE_IN_VOLTS:
      for(i=0; i<_numberSpectrumPoints; i++)
      {
        if(_outputSpectrum[i] > 0.)
          _outputSpectrum[i] = sqrt(_outputSpectrum[i]);
      }
      break;
  }
  delete [] windowed_data;
  return _outputSpectrum;
}

// ############################# Public Method #################################
// phaseSpectrum -- Calculates the phase spectrum of real input data.
//
// Input:       realInput:      Real time series data
//              imagInput:      Imaginary time series data
//              numberPoints:   Number of points in above array
//              
//
// Output:                      Pointer to the phase spectrum data in radians
//
// Notes:
// ############################# Public Method #################################
float *PowerSpectrum::phaseSpectrum(const float *realInput, const float *imagInput, int numberPoints)
{
  int   i;
  const  float *fft_magnitude;
//
// Allocate space
//
  _numberPhasePoints    = numberPoints;
  if(_oldPhasePoints < numberPoints)
  {
    _oldPhasePoints     = numberPoints;
    delete [] _outputPhase;
    _outputPhase                = new float[numberPoints];
  }
//
// First find the FFT, and then convert to phase.
//
  _fft->fft(realInput, imagInput, numberPoints);
  fft_magnitude = _fft->magFloatOutput(COMPLEX_PHASE);
//
// Save output
//
  for(i=0; i<numberPoints; i++)
    _outputPhase[i] = fft_magnitude[i];
  return _outputPhase;
}

// ############################# Public Method #################################
// complexSpectrum -- Calculates the complex spectrum of complex input data.
//
// Input:       realInput:      Real time series data (over-written on output)
//              imagInput:      Imaginary time series data (over-written on output)
//              numberPoints:   Number of points in above array
//              windowData:     YES = window data before taking spectrum
//              
//
// Output:                      Pointer to the phase spectrum data in radians
//
// Notes:
// ############################# Public Method #################################
void PowerSpectrum::complexSpectrum(float *realInput, float *imagInput, int numberPoints, LOGICAL windowData)
{
  int   i;
  const float *fft_real, *fft_imag;
  float *windowed_real, *windowed_imag;
  const double *window;
//
//  Window data if requested:
//
  if(windowData)
  {
    windowed_real       = new float[numberPoints];
    windowed_imag       = new float[numberPoints];
    window              = _dataWindow->windowFunction(numberPoints);
    for(i=0; i<numberPoints; i++)
    {
      windowed_real[i]  = window[i]*realInput[i];
      windowed_imag[i]  = window[i]*imagInput[i];
    }
    _fft->fft(windowed_real, windowed_imag, numberPoints);
    delete [] windowed_real;
    delete [] windowed_imag;
  }
//
// Else, no window, so just FFT input
//
  else
    _fft->fft(realInput, imagInput, numberPoints);
//
// Save output
//
  fft_real              = _fft->realFloatOutput();
  fft_imag              = _fft->imagFloatOutput();
  for(i=0; i<numberPoints; i++)
  {
    realInput[i]        = fft_real[i];
    imagInput[i]        = fft_imag[i];
  }
  return;
}

// ############################# Public Method #################################
// autoFromPSD -- Calculates the autocovariance from the averaged power spectrum
//
//
// Input:       inputSpectrum:  Real power spectrum data
//              numberPoints:   Number of points in above array
//              type:           Complex to real conversion type, see defines in FFT.h file
//              
//
// Output:      scale:          Magnitude of autocorrelation at lag = 0
//                              Pointer to the output power spectrum data
//
// Notes:
// ############################# Public Method #################################
float *PowerSpectrum::autoFromPSD(const float *inputSpectrum, int numberPoints, int type, float *scale)
{
  int   i;
  const float *auto_magnitude;

  _numberAutoPoints  = numberPoints;
  if(_oldAutoPoints < numberPoints)
  {
    _oldAutoPoints   = numberPoints;
    delete [] _autocorrelation;
    _autocorrelation = new float[numberPoints];
  }
//
// Now take inverse FFT to find autocorrelation
//
  _fft->setOperationType(INVERSE_FFT);
  _fft->realFFT(inputSpectrum, numberPoints);
  auto_magnitude    = _fft->magFloatOutput(MAGNITUDE_IN_WATTS);  // Use to get scale
  *scale            = auto_magnitude[0];
  auto_magnitude    = _fft->magFloatOutput(type);
  for(i=0; i<numberPoints; i++)
    _autocorrelation[i] = auto_magnitude[i];
  _fft->setOperationType(FORWARD_FFT);           // Restore for future FFTs
  return _autocorrelation;
}

// ############################# Public Method #################################
// autocorrelationOf -- This  method finds the modified short-time autocorrelation
//                      function from (4.30), page 143, Rabiner and Schafer. A pointer
//                      to the output autocorrelation is returned. 
//
// Input:       realInput:      Real part of time series data
//              imagInput:      Imaginary part of time series data
//              numberPoints:   Size of each FFT to take
//              windowLength:   Window to calculate autocorrelation over
//              lagLength:      Number of lags of autocorrelation to calculate
//              type:           Complex to real conversion type, see defines in FFT.h file
//
// Output:      scale:          Magnitude of autocorrelation at lag = 0
//                              Pointer to the output power spectrum data
//
// Notes:
// 1. Result is divided by N.
// ############################# Public Method #################################
float *PowerSpectrum::autocorrelationOf(const float *realInput, const float *imagInput, int numberPoints,
                                    int windowLength, int lagLength, int type, float *scale)
{
  int       m,k;   
  float     real_part, imag_part;
  const     double *window;
        
// Check # of points: 
  if( (numberPoints<windowLength) || (windowLength<=lagLength) )
    return NULL;                                             

  _numberAutoPoints  = lagLength;
  if(_oldAutoPoints < lagLength)
  {
    _oldAutoPoints = lagLength;
    delete [] _autocorrelation;
    _autocorrelation = new float[lagLength];
  }
  window        = _dataWindow->windowFunction(windowLength);
/****************************
 * Now, calculate it:       *  
 * choose n = 0  in (4.36)  *
 ****************************/
  switch(type)
  {
    case MAGNITUDE_IN_WATTS:
    default:
      for(k=0; k<lagLength; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<windowLength-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        _autocorrelation[k] = sqrt(real_part*real_part + imag_part*imag_part);
        _autocorrelation[k] /= (double) (windowLength-k);
      }
      *scale = _autocorrelation[0];
      break;
    case COMPLEX_PHASE:
      for(k=0; k<lagLength; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<windowLength-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        _autocorrelation[k] = atan2(imag_part, real_part);
      }
      *scale = 1;
      break;
    case COMPLEX_REAL_PART:
      for(k=0; k<lagLength; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<windowLength-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        _autocorrelation[k] = real_part/(double) (windowLength-k);
        if(k==0)
          *scale = sqrt(real_part*real_part + imag_part*imag_part)/(double) (windowLength-k);
      }
      break;
    case COMPLEX_IMAG_PART:
      for(k=0; k<lagLength; k++)
      { 
        real_part = imag_part = 0.;
        for(m=0; m<windowLength-k; m++)
        {
          real_part += realInput[m]*window[m]*realInput[m+k]*window[m+k] + 
                      imagInput[m]*window[m]*imagInput[m+k]*window[m+k];
          imag_part += imagInput[m]*window[m]*realInput[m+k]*window[m+k] - 
                      realInput[m]*window[m]*imagInput[m+k]*window[m+k];
        }
        _autocorrelation[k] = imag_part/ (windowLength-k);
        if(k==0)
          *scale = sqrt(real_part*real_part + imag_part*imag_part)/(double) (windowLength-k);
      }
      break;
  }
  return _autocorrelation;
}

// ############################# Public Method #################################
// meanIndexOf -- Calculates the mean index of the power spectrum
//
// Input:       inputSpectrum:  Real power spectrum data
//              numberPoints:   Number of points in above array
//              
// Output:                      index of the center of the spectrum (float value)
//
// Notes:
// 1. The center of mass is calculated according to _meanIndexType, which should
//    be preset.  If _meanIndexType==MEAN_FROM_SCM, the spectral center of mass is
//    returned.  If _meanIndexType==MEAN_FROM_THRESHOLD, the mean index is found
//    as the mean value of the crossing points of the spectrum from DB_THRESHOLD.
// 2. This function also finds the _twoSidedBandwidth using the crossing points,
//    and can be returned by post-calling getTwoSidedBandwidth().
// ############################# Public Method #################################
float PowerSpectrum::meanIndexOf(const float *inputSpectrum, int numberPoints)
{
  int       i, max_index, flo, fhi;
  float     max_spect, bandwidth_threshold;
  double    mean_sum, spectrum_sum;
//
// First find the two sided bandwidth by checking when
// spectrum falls below threshold:
//
  max_spect = inputSpectrum[0];
  max_index = 0;
  for(i=1; i<numberPoints; i++)
  {
    if(inputSpectrum[i]>max_spect)
    {
      max_index = i;
      max_spect = inputSpectrum[i];
    }
  }
  bandwidth_threshold = max_spect/pow(10., DB_THRESHOLD/10.);
  flo = fhi = 0;
  for(i=max_index; i>0; i--)
    if(inputSpectrum[i] < bandwidth_threshold)
    {
      flo = i;
      break;
    }
  for(i=max_index; i<numberPoints; i++)
    if(inputSpectrum[i] < bandwidth_threshold)
    {
      fhi = i;
      break;
    }
  _twoSidedBandwidth = (fhi - flo)*_fft->deltaFrequency();
  if(_meanIndexType  == MEAN_FROM_THRESHOLD)
    mean_sum        = (float)(fhi+flo)/2.;
// Else, find spectral center of mass
  else
  {
    mean_sum        = spectrum_sum = 0.;
    for(i=1; i<numberPoints; i++)
    {
      mean_sum      += (double)i*inputSpectrum[i];
      spectrum_sum  += inputSpectrum[i];
    }
    if(spectrum_sum>0.)
     mean_sum       /= spectrum_sum;
  }
  return (float)mean_sum;
}

// ############################# Public Method #################################
// maxBinFrequency -- Calculates the frequency of the bin having the max value
//
// Input:       inputSpectrum:  Real power spectrum data
//              numberPoints:   Number of points in above array
//              flo:            Low frequency search point
//              fhi:            high frequency search point
//              
// Output:                      bin having the max value
//
// Notes:
// ############################# Public Method #################################
float PowerSpectrum::maxBinFrequency(const float *inputSpectrum, int numberPoints, float flo, float fhi)
{
  int       i, start, end, max_bin;
  float     max_spect, delta_f, max_freq;
//
// First, get the freq spacing, and the search range:
//
  delta_f   = _fft->deltaFrequency();
  start     = ROUND(flo/delta_f);
  start     = MAX(0, start);
  start     = MIN(start, numberPoints-1);

  end       = ROUND(fhi/delta_f);
  end       = MIN(end, numberPoints-1);
//
//  Now, search for max
//
  max_spect = inputSpectrum[start];
  max_bin   = start;
  for(i=start+1; i<end; i++)
  {
    if(max_spect < inputSpectrum[i])
    {
      max_spect = inputSpectrum[i];
      max_bin   = i;
    }
  }
  max_freq  = max_bin * delta_f;
  return max_freq;
}

// ############################# Public Method #################################
// smoothSpectrum -- Smooths the input spectrum by passing a averaging window
//                   over the spectrum array.
//
// Input:       inputSpectrum:  Real power spectrum data
//              numberPoints:   Number of points in above array
//              windowSize:     Size of averaging window
//              
// Output:      inputSpectrum:  smoothed spectrum
//
// Notes:
// ############################# Public Method #################################
void PowerSpectrum::smoothSpectrum(float *inputSpectrum, int numberPoints, int windowSize)
{
  int       i, j, start, end;
  static    int old_number_points=0;
  float     scale;
  double    sum;
  static    float *temp_array=NULL;
//
// Allocate the temp array and copy the spectrum into the temp array
//
  if(old_number_points < numberPoints)
  {
    delete [] temp_array;
    old_number_points   = numberPoints;
    temp_array          = new float[numberPoints];
  }
  for(i=0; i<numberPoints; i++)
    temp_array[i]       = inputSpectrum[i];
//
// Now, pass the windowing array over the temp array:
//
  for(i=0; i<numberPoints; i++)
  {
    start   = MAX(0, (i-windowSize/2));
    end     = MIN((numberPoints-1), (i+windowSize/2));
    sum     = 0.;
    scale   = 1.;
    if(end >= start)
      scale = 1./(end-start+1);
    for(j=start; j<=end; j++)
      sum   += temp_array[j];
    inputSpectrum[i]    = sum*scale;
  }
  return;
}

// ############################# Public Method #################################
// shiftRight -- Circular shift the input array (e.g., spectrum) to the right.
//
// Input:       inputArray:     Real power spectrum data
//              shift:          Number of right shifts to do
//              arraySize:      Number of points in above array
//              
// Output:      None
//
// Notes:
// 1. This is useful for shifting a two sided spectrum for display.
// ############################# Public Method #################################
void PowerSpectrum::shiftRight(float *inputArray, int shift, int arraySize)
{
  int   i, j;
  float *temp;

  shift %= arraySize;
  temp  = new float[arraySize];
  j     = arraySize - shift;
  j     %= arraySize;
  for(i=0; i<arraySize; i++)
  {
    temp[i] = inputArray[j++];
    j       %= arraySize;
  }
  for(i=0; i<arraySize; i++)
    inputArray[i]   = temp[i];
  delete [] temp;
  return;
}

// ############################# Public Method #################################
// shiftRight -- Circular shift the input array (e.g., spectrum) to the right.
//
// Input:       inputArray:     Real power spectrum data
//              shift:          Number of right shifts to do
//              arraySize:      Number of points in above array
//              
// Output:      None
//
// Notes:
// 1. This is useful for shifting a two sided spectrum for display.
// ############################# Public Method #################################
void PowerSpectrum::shiftRight(Complex *inputArray, int shift, int arraySize)
{
  int       i, j;
  Complex   *temp;

  shift %= arraySize;
  temp  = new Complex[arraySize];
  j     = arraySize - shift;
  j     %= arraySize;
  for(i=0; i<arraySize; i++)
  {
    temp[i] = inputArray[j++];
    j       %= arraySize;
  }
  for(i=0; i<arraySize; i++)
    inputArray[i]   = temp[i];
  delete [] temp;
  return;
}

// ############################# Public Method #################################
// magnitudeToDB -- Convert the input array from watts to dB.
//
// Input:       inputSpectrum:  Real power spectrum data in watts
//              numberPoints:   Number of points in above array
//              
// Output:                      inputSpectrum is modified
//
// Notes:
// ############################# Public Method #################################
void PowerSpectrum::magnitudeToDB(float *inputSpectrum, int numberPoints)
{
  int i;

  for(i=0; i<numberPoints; i++)
  {
    if(inputSpectrum[i]>0.)
      inputSpectrum[i] = 10.*log10(inputSpectrum[i]);
  }
  return;
}

// ############################# Public Method #################################
// normalizeAuto -- Normalize the autocorrelation by a scale factor, e.g., its dc
//                  value.
//
// Input:       autoc:          Input autocorrelation factor
//              numberPoints:   Number of points in above array
//              scale:          Scale factor to divide through by
//              
// Output:                      autoc is modified
//
// Notes:
// ############################# Public Method #################################
void PowerSpectrum::normalizeAuto(float *autoc, int numberPoints, float scale)
{
  int i;

  if(scale != 0.)
    for(i=0; i<numberPoints; i++)
      autoc[i] /= scale;
  return;
}


// ############################# Public Method #################################
// complexToMag -- Convert complex array to magnitude output
//
// Input:       inputArray:     Real power spectrum data
//              arraySize:      Size of the array
//              type:           Complex to float conversion
//              
// Output:                      index of the center of the spectrum (float value)
//
// Notes:
// 1. This is useful for shifting a two sided spectrum for display.
// ############################# Public Method #################################
float *PowerSpectrum::complexToFloat(const Complex *inputArray, int arraySize, int type)
{
  int       i;
  double    magnitude;

  if(_numberFloatPoints < arraySize)
  {
    _numberFloatPoints   = arraySize;
    delete [] _outputFloat;
    _outputFloat     = new float[_numberFloatPoints];
  }
  switch(type)
  {
    case MAGNITUDE_IN_WATTS:
    default:
      for(i=0; i<_numberFloatPoints; i++)
      {
        _outputFloat[i]  = abs(inputArray[i]);
      }
      break;
    case MAGNITUDE_IN_DB:
      for(i=0; i<_numberFloatPoints; i++)
      {
        magnitude           = abs(inputArray[i]);
        if(magnitude>0.)
          _outputFloat[i]    = 10.*log10(magnitude);
        else
          _outputFloat[i]    = 0.;
      }
      break;
    case COMPLEX_PHASE:
      for(i=0; i<_numberFloatPoints; i++)
      {
        _outputFloat[i]  = arg(inputArray[i]);
      }
      break;
    case COMPLEX_REAL_PART:
      for(i=0; i<_numberFloatPoints; i++)
      {
        _outputFloat[i]  = real(inputArray[i]);
      }
      break;
    case COMPLEX_IMAG_PART:
      for(i=0; i<_numberFloatPoints; i++)
      {
        _outputFloat[i]  = imag(inputArray[i]);
      }
      break;
  }

  return _outputFloat;
}

// ############################# Public Method #################################
// oneSidedToTwoSided -- Converts the one sided spectrum to a two-sided spectrum.
//
// Input:                       None
//              
// Output:                      The two sided spectrum
//
// Notes:
// ############################# Public Method #################################
float *PowerSpectrum::oneSidedToTwoSided()
{
  int       i, j, old_number_points;
  float     *temp_array;
//
// Allocate the temp array and copy the spectrum into the temp array
//
  old_number_points     = _numberSpectrumPoints;
  _numberSpectrumPoints  = 2*_numberSpectrumPoints - 1;   // Two sided has twice the number of points
  _oldSpectrumPoints     = _numberSpectrumPoints;
  temp_array            = new float [_numberSpectrumPoints];
//
// Copy spectrum into temp:
//
  j     = old_number_points - 1;
  for(i=0; i<old_number_points; i++)
     temp_array[i]  = _outputSpectrum[j--];
  j     = 1;
  for(; i<_numberSpectrumPoints; i++)
     temp_array[i]  = _outputSpectrum[j++];
//
// Remove old storage, and return
//
  delete [] _outputSpectrum;
  _outputSpectrum    = temp_array;
  return _outputSpectrum;
}

