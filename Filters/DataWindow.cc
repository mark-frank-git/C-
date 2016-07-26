/************************************************************************
 *                                                                      *
 * This class implements a data windowing function.  It is based on the *
 * paper, "On the use of windows for harmonic analysis with the dis-    *
 * crete Fourier transform," by F.J. Harris, proc. IEEE vol. 66, no.1.  *
 *                                                                      *
 * File:DataWindow.cc                                                   *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/29/02  - Extracted from DigitalFilter                         *
 ************************************************************************/

#include "DataWindow.h"                                   // Object prototypes
#if defined(WIN32)
#include <Specfuns/specfuns.h>
#include <GNU/Complex.h>
#else
#include "specfuns.h"
#include "Complex.h"
#endif

#include <stdio.h>
#include <math.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)      ((a) >= 0 ? (a) : (-a))
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
//
// The following are taken from Table 1 in the reference paper:
// Note that the bandwidths will depend on the alpha values for the
// Kaiser, Gaussian and Hanning.  The numbers below reflect the BWs for
// the default alpha values
// Also, note that the BWs for the Tukey and Raised Cosine are bogus.
float bw_6dB[NUMBER_WINDOW_TYPES]           = {1.21, 1.78, 2.0,  1.81, 1.81, 2.72, 2.18, 2.20, 1.38, 1.5};
float bw_3dB[NUMBER_WINDOW_TYPES]           = {0.89, 1.28, 1.44, 1.30, 1.66, 1.90, 1.55, 1.57, 1.01, 1.5};


#define NUMBER_COEFFICIENTS 4
float blackman3_coeffs[NUMBER_COEFFICIENTS] = {0.42323, -0.49755, 0.07922, 0.};
float blackman4_coeffs[NUMBER_COEFFICIENTS] = {0.35875, -0.48829, 0.14128, -0.01168};


// ############################# Private Method ###############################
// dirichletKernel -- This routine returns the complex Dirichlet kernel
//                   at the input frequency.
// Input:       omega:          digital frequency [-PI, PI]
//              length:         window length
//          
// Output:                      magnitude of the DFT at omega
//
// Notes:
// ############################# Private Method ###############################
Complex DataWindow::dirichletKernel(double omega, int length)
{
  Complex       carg;
  double        sin_omega_2;
  if(omega==0.)
    return Complex(length, 0.);
  sin_omega_2   = sin(omega/2.);
  if(sin_omega_2 != 0.)
  {
    carg        = Complex(0., -omega*(length-1)/2.);
    return exp(carg)*sin(length*omega/2.)/sin_omega_2;
  }
  return Complex(1., 0.);                       // undefined?
}

// ############################# Private Method ###############################
// dirichletMagnitude -- This routine returns the magnitude of the Dirichlet kernel
//                   at the input frequency.
// Input:       omega:          digital frequency [-PI, PI]
//              length:         window length
//          
// Output:                      magnitude of the DFT at omega
//
// Notes:
// ############################# Private Method ###############################
double DataWindow::dirichletMagnitude(double omega, int length)
{
  double        sin_omega_2;
  if(omega==0.)
    return length;
  else if( (omega<(-PI)) || (omega>PI) )
    return 0.;
  sin_omega_2   = sin(omega/2.);
  if(sin_omega_2 != 0.)
    return sin(length*omega/2.)/sin_omega_2;
  else
    return 1.;                  // undefined?
}
// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the DataWindow class.
//
// Input:           type:           window type
//
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
DataWindow::DataWindow(int type)
{

  _windowType       = type;
//
// Set default values:
//
  setHanningAlpha(HANNING_ALPHA);
  setGaussianAlpha(GAUSSIAN_ALPHA);
  setKaiserBeta(KAISER_ALPHA);
  setTukeyAlpha(TUKEY_ALPHA);
  setRaisedCosineSamplesPerSymbol(RAISED_SAMPLES);

  _window           = NULL;
  
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the DataWindow class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
DataWindow::~DataWindow()
{
  delete [] _window;
  return;
}

// ############################# Public Method ###############################
// bandwidth3dB -- This routine returns the 3 dB bandwidth of the window in bins
//                 as taken from Table 1 in the reference.
// Input:                           None
//          
// Output:                          3 dB BW in bins
//
// Notes:
// 1. We should probably calculate this bandwidth on the fly rather than use
//    a table look up.
// ############################# Public Method ###############################
float DataWindow::bandwidth3dB()
{
  return bw_3dB[_windowType];
}

// ############################# Public Method ###############################
// bandwidth6dB -- This routine returns the 6 dB bandwidth of the window in bins
//                 as taken from Table 1 in the reference.
// Input:                           None
//          
// Output:                          6 dB BW in bins
//
// Notes:
// 1. We should probably calculate this bandwidth on the fly rather than use
//    a table look up.
// ############################# Public Method ###############################
float DataWindow::bandwidth6dB()
{
  return bw_6dB[_windowType];
}

// ############################# Public Method ###############################
// windowFunction -- This routine calculates the FIR filter windowing function:
// Input:           length:         Length of the window function
//          
// Output:                          Pointer to the window array
//
// Notes:
// 1. The _coherentGain and _noncoherentGain variables are also calculated here.
//    They are equal to the dc signal gain Harris Eqn. (13), and the noise power
//    gain, Harris Eqn. (9)
// ############################# Public Method ###############################
const double *DataWindow::windowFunction(int length)
{
  int           i, j;
  double        arg;
  double        sum_sq, sum;
  float         *coefficients;
  
  if(_window != NULL)
    delete [] _window;
  _window = new double[MAX(2, length)];

  if(length < 2)                                        // Check for illegal length
    return _window;
 
  coefficients  = blackman4_coeffs;                     // Default for BLACMAN_HARRIS4
  arg           = TWOPI/(length-1);                     // Default value for arg, note Harris uses TWOPI/N?
  switch(_windowType)
  {
    case RECTANGULAR:
    default:      
      for(i=0; i<length; i++)
        _window[i]  = 1.;
      break;
    case TRIANGULAR:
      findTriangularWindow(length, _window);
      break;
    case HANNING:
      findHanningWindow(length, _hanningAlpha, _window);
      break;
    case HAMMING:
     findHammingWindow(length, _window);
     break;
    case BLACKMAN_HARRIS3:
      coefficients  = blackman3_coeffs;
    case BLACKMAN_HARRIS4:
      for(i=0; i<length; i++)
      {
        _window[i]  = coefficients[0];
        for(j=1; j<NUMBER_COEFFICIENTS; j++)
          _window[i]    += coefficients[j]*cos(i*j*arg);
      }
      break;
    case KAISER:
      findKaiserWindow(length, _kaiserAlpha, _window);
      break;
    case GAUSSIAN_WINDOW:
      findGaussianWindow(length, _gaussianAlpha, _window);
      break;
    case TUKEY_WINDOW:
      findTukeyWindow(length, _tukeyAlpha, _window);
      break;
    case RAISED_COSINE:
      findRaisedCosineWindow(length, _raisedSamples, _window);
      break;
  }
//
// Calculate the normalization factors:
//
  sum_sq        = sum           = 0.;
  for(i=0; i<length; i++)
  {
    sum         += _window[i];
    sum_sq      += _window[i]*_window[i];
  }
  if(length > 0)
  {
    _coherentGain       = sum*sum/length;
    _noncoherentGain    = sum_sq/length;
  }
  
  return _window;
}

#define CONVOLUTION_POINTS      500
// ############################# Public Method ###############################
// dftMagnitudeAt -- This routine calculates magnitude of the DFT of the window
//                   function at the given frequency
// Input:       omega:          digital frequency [-PI, PI]
//              length:         length of the window
//          
// Output:                      magnitude of the DFT at omega
//
// Notes:
// ############################# Public Method ###############################
double DataWindow::dftMagnitudeSquaredAt(double omega, int length)
{
  int           m, sum_length;
  double        window_mag;
  double        arg, io_beta, sqrt_arg;
  const         double *gaussian_window;
  float         *coefficients;
  Complex       blackman_sum, gaussian_sum, exp_arg;

  if(length < 2)                // error check
    return 1.;
//
// Calculate based on window type:
//
  coefficients  = blackman4_coeffs;                     // Default for BLACMAN_HARRIS4
  blackman_sum  = Complex(0., 0.);
  gaussian_sum  = Complex(0., 0.);
  sum_length    = 4;
  switch(_windowType)
  {
    case RECTANGULAR:
    default:      
      window_mag        = dirichletMagnitude(omega, length);
      break;
    case TRIANGULAR:
      window_mag        = 2*dirichletMagnitude(omega, length/2)/length;
      break;
    case HANNING:
      window_mag        = 0.5*dirichletMagnitude(omega, length) +
                          0.25*(dirichletMagnitude(omega-TWOPI/length, length) +
                                dirichletMagnitude(omega+TWOPI/length, length));
      break;
    case HAMMING:
      window_mag        = 0.54*dirichletMagnitude(omega, length) +
                          0.5*(1-0.54)*(dirichletMagnitude(omega-TWOPI/length, length) +
                                       dirichletMagnitude(omega+TWOPI/length, length));
      break;
    case BLACKMAN_HARRIS3:
      coefficients      = blackman3_coeffs;
      sum_length        = 3;
    case BLACKMAN_HARRIS4:
      for(m=0; m<sum_length; m++)
        blackman_sum    += coefficients[m]*(dirichletKernel(omega-TWOPI*m/length, length) +
                                             dirichletKernel(omega+TWOPI*m/length, length))/2.;
      window_mag        = abs(blackman_sum);
      break;
    case KAISER:
      arg               = PI*_kaiserAlpha;
      io_beta           = modbes(arg, 0, 0);
      sqrt_arg          = (arg*arg - (length*length*omega*omega/4.));
      if(sqrt_arg > 0.)
        window_mag      = length*sinh(sqrt(sqrt_arg))/sqrt(sqrt_arg)/io_beta;
      else if (sqrt_arg < 0.)
        window_mag      = length*sin(sqrt(-sqrt_arg))/sqrt(-sqrt_arg)/io_beta;
      else
        window_mag      = 1.;
      break;
    case GAUSSIAN_WINDOW:                       // Perform explicit FT, couldn't get Harris' formula to work
      gaussian_window   = windowFunction(length);
      for(m=0; m<length; m++)
      {
        exp_arg         = Complex(0., -omega*m);
        gaussian_sum    += gaussian_window[m]*exp(exp_arg);
      }
      window_mag        = abs(gaussian_sum);
      break;
  }
  return        window_mag*window_mag;
}

// ############################# Private Method ###############################
// findTriangularWindow -- This routine calculates triangular window:
// Input:           length:         Length of the window function
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow::findTriangularWindow(
                            const int length,
                            double    *window) const
{
  int i;
  double arg            = 2.0/(length-1);
  for(i=0; i<=(length-1)/2; i++)
    window[i]   = i*arg;
  for(   ; i<length; i++)
    window[i]    = 2. - i*arg;
  return;
}
// ############################# Private Method ###############################
// findHanningWindow -- This routine calculates Hanning window:
// Input:           length:         Length of the window function
//                                                                      hanningAlpha:           Alpha
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow:: findHanningWindow (
                            const int   length,
                            const float hanningAlpha,
                            double       *window) const
{
  int i;
  double arg    = TWOPI/(length-1);
  for(i=0; i<length; i++)
  {
    if(hanningAlpha    == 1)
      window[i]    = sin(0.5*arg*i);
    else
      window[i]    = 0.5*(1. - cos(arg*i));            // Normal Hann window
  }
  return;
}


// ############################# Private Method ###############################
// findHammingWindow -- This routine calculates Hamming window:
// Input:           length:         Length of the window function
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow:: findHammingWindow (
                            const int   length,
                            double      *window) const
{
  int i;
  double arg    = TWOPI/(length-1);
  for(i=0; i<length; i++)
  {
    window[i]  = 0.54 - 0.46*cos(arg*i);
  }
  return;
}


// ############################# Private Method ###############################
// findKaiserWindow -- This routine calculates Kaiser window:
// Input:           length:         Length of the window function
//                                  kaiserAlpha:                Alpha
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow:: findKaiserWindow (
                            const int   length,
                                                                                                          const float kaiserAlpha,
                            double      *window) const
{
  int i;
  double arg      =  PI*kaiserAlpha;
  double io_beta  = modbes(arg, 0, 0);
  int half_length  = length/2;
  int j            = half_length;
  for(i=0; i<=half_length; i++)
  {
    double arg_bessel  = arg*sqrt(1. - ((double)j*j/half_length/half_length));
    window[i]          = modbes(arg_bessel, 0, 0)/io_beta;
    j--;
  }
  if(length%2)
    i           = half_length+1;                    // length is odd
  else
    i           = half_length;                      // length is even
  j             = half_length-1;
  for(; i<length; i++)
    window[i]   = window[j--];
  return;
}

// ############################# Private Method ###############################
// findGaussianWindow -- This routine calculates Gaussian window:
// Input:           length:         Length of the window function
//                  gaussianAlpha:    Alpha
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow:: findGaussianWindow (
                            const int   length,
                            const float gaussianAlpha,
                            double      *window) const
{
  int i;
  double arg        = TWOPI/(length-1);
  double mid_point  = (length - 1.)/2.;
  for(i=0; i<length; i++)
  {
    arg         = (i-mid_point)* gaussianAlpha/length;
    arg         = 2.*arg*arg;
    window[i]  = exp(-arg);
  }
  return;
}


// ############################# Private Method ###############################
// findTukeyWindow -- This routine calculates Tukey window:
// Input:           length:         Length of the window function
//                  tukeyAlpha:      Alpha
//          
// Output:          window:         Pointer to the window array
//
// Notes:
//  None
// ############################# Private Method ###############################
void DataWindow:: findTukeyWindow (
                            const int   length,
                            const float tukeyAlpha,
                            double      *window) const
{
  int i              = 0;
  double arg        = TWOPI/(length-1);
  double den;
  int half_length    = length/2;
  int lower          = -half_length;
  int upper          = half_length;
  if(length%2 == 0)        // Even length
  {
    lower    = -half_length;
    upper    = half_length-1;
  }
  for(int n=lower; n<=upper; n++)
  {
    if(ABS(n) <= ROUND(((1.-_tukeyAlpha)*half_length)) )
      window[i]  = 1.;
    else
    {
          if(n > 0)
          {
            arg   = ABS(n) - (1.-tukeyAlpha)*upper;
            den  = upper * tukeyAlpha;
          }
          else
          {
            arg   = ABS(n) - (1.-tukeyAlpha)*lower;
            den  = lower * tukeyAlpha;
          }
          if(den != 0.)
            arg  /= den;
          window[i]  = 0.5*(1+cos(PI*arg));
    }
    i++;
  }
  return;
}


// ############################# Private Method ###############################
// findRaisedCosineWindow -- This routine calculates the EDGE raised cosine window:
// Input:           length:            Length of the window function
//                  samplesPerSymbol:  Alpha
//
// Output:          window:         Pointer to the window array
//
// Notes:
//  See GSM 05.05, 8-PSK modulation
// ############################# Private Method ###############################
void DataWindow:: findRaisedCosineWindow (
                            const int   length,
                            const float samplesPerSymbol,
                            double      *window) const
{
  int limit1        = ROUND(2.25*samplesPerSymbol);
  limit1            = MIN(limit1, length);
  int i;
  for(i=0; i<=limit1; i++)
  {
    window[i]        = 0.5*(1+cos(PI*(i/samplesPerSymbol/2.25-1.)));
  }
  int limit2        = ROUND(5.25*samplesPerSymbol);
  limit2            = MIN(limit2, length);
  for( ;i<limit2; i++)
  {
    window[i]        = 1.;
  }
  int limit3        = ROUND(7.5*samplesPerSymbol);
  limit3            = MIN(limit3, length);
  for( ;i<limit3; i++)
  {
    window[i]        = 0.5*(1+cos(PI*(i/samplesPerSymbol-5.25)/2.25));
  }
  for( ;i<length; i++)
  {
    window[i]        = 0.;
  }
  return;
}

