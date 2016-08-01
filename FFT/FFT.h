#ifndef _FFT_H
#define _FFT_H  1
/************************************************************************
 *                                                                      *
 * This subclass of object implements an object for computing fast      *
 * Fourier transforms, as well as DFTs and DCTs                         *
 *                                                                      *
 * File: FFT.h                                                          *
 *                                                                      *
 ************************************************************************/
#if defined (_WIN32) && defined (NeXT_PDO)
#include <GNU/Complex.h>
#else
#include "Complex.h"
#endif
#include "fft_routines.h"


#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif

//
// These routines are located in fft_routines.cc
//

    
//! FFT output types enum
/*! The output types enum categories the different types of output from the FFT */
enum output_type_enum
{
  MAGNITUDE_IN_DB,                                      //!< return magnitude in dB
  MAGNITUDE_IN_WATTS,                                   //!< return magnitude in watts
  MAGNITUDE_IN_VOLTS,                                   //!< return magnitude in volts
  COMPLEX_PHASE,                                        //!< return phase in radians
  COMPLEX_REAL_PART,                                    //!< return real part of result
  COMPLEX_IMAG_PART                                     //!< return imaginary part of result
};

//! FFT Class
/*!
   The FFT class allows an easy interface to the FFT routines.
*/
class FFT
{
private:
  int           _numberPoints;                          //!< Number of FFT points
  int           _operationType;                         //!< Forward or Inverse FFT
  int           _powerOfTwo;                            //!< Power of two for number of points
  int           _oldFrequencyPoints;                    //!< Used for allocating memory
  int           _oldTwoSidedPoints;                     //!< Used for allocating memory
  int           _oldComplexPoints;                      //!< Used for allocating memory
  int           _oldRealPoints;                         //!< Used for allocating memory
  int           _oldImagPoints;                         //!< Used for allocating memory
  int           _oldFloatMagPoints;                     //!< Used for allocating memory
  int           _oldDoubleRealPoints;                   //!< Used for allocating memory
  int           _oldDoubleImagPoints;                   //!< Used for allocating memory
  int           _oldDoubleMagPoints;                    //!< Used for allocating memory
  float         _samplingFreq;                          //!< Sampling frequency in Hertz
  float         *_frequencyPoints;                      //!< frequencies of FFT samples
  float         *_twoSidedFrequency;                    //!< Negative and positive freqs
  float         *_floatRealOutput;                      //!< Real output from FFT
  float         *_floatImagOutput;                      //!< Imag output from FFT
  float         *_floatMagOutput;                       //!< Magnitude output from FFT
  double        *_doubleRealOutput;                     //!< Real output from FFT
  double        *_doubleImagOutput;                     //!< Imag output from FFT
  double        *_doubleMagOutput;                      //!< Magnitude output from FFT
  Complex       *_complexArray;                         //!< Complex input/output to FFT

public:
//
// Public methods:
//
//! Class constructor, the sampling frequency is used for frequency axis scaling
  FFT(float sampling=1000.);
//!< Class destructor
  ~FFT();

/**********************
 * Set parameters:    *
 **********************/
//! setSamplingFrequency() allows the setting of a new sampling frequency
/**
  * a normal member taking one argument with no return value
  * @param frequency the new sampling frequency in Hertz
  * @see frequencyPoints()
  * @see twoSidedFrequency()
  * @return void
  */
  void setSamplingFrequency(float frequency);
  
//! setOperationType() allows the setting of a new operation type (forward or inverse)
/**
  * a normal member taking one argument with no return value
  * @param type -1 = inverse FFT, 1 = forward FFT
  * @return void
  */
  void setOperationType(int type);

/**********************
 * Get parameters:    *
 **********************/
  float samplingFreq()          { return _samplingFreq;}        //!< returns the sampling frequency in Hz
  float *frequencyPoints();                                     //!< returns a set of frequency points
                                                                //!< for plotting FFT results
  float *twoSidedFrequency();                                   //!< returns a set of 2 sided frequency points
                                                                //!< for plotting FFT results
  float deltaFrequency(int number=0);                           //!< Returns the spacing between frequency pts at which
                                                                //!< the FFT was calculated, if number is passed in,
                                                                //!< it is used for calculation, otherwise size of
                                                                //!< last FFT calculation is used

//
// The following functions calculate FFTs, call these before getting the outputs
//
//! isPowerOfTwo() returns YES if the input integer is a power of two
/**
  * a normal member taking one arguments with a LOGICAL return value
  * @param n the input integer
  * @return YES or NO
  */
  LOGICAL isPowerOfTwo(int n);
  
//! fft() calculates the FFT of an array of float real and imaginary data, the input
//! arrays are not modified
/**
  * a normal member taking three arguments with no return value
  * @param realInput an array of real data
  * @param imagInput an array of imag data
  * @param n the number of points in the array
  * @see realFloatOutput(), imagFloatOutput(), etc for getting results
  * @return void
  */
  void fft(const float *realInput, const float *imagInput, int n);
  
//! fft() calculates the FFT of an array of double real and imaginary data, the input
//! arrays are not modified
/**
  * a normal member taking three arguments with no return value
  * @param realInput an array of real double data
  * @param imagInput an array of imag double data
  * @param n the number of points in the array
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void fft(const double *realInput, const double *imagInput, int n);

//! fft() calculates the FFT of an array of complex data, the input
//! array is not modified
/**
  * a normal member taking three arguments with no return value
  * @param complexInput an array of real data
  * @param n the number of points in the array
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void fft(const Complex *complexInput, int n);

//! czt() calculates the chirp Z transform of an array of real/imag data, the input
//! arrays are not modified. This function is useful, when the number of points is
//! not a power of two, or a custom frequency spacing is desired.
/**
  * a normal member taking six arguments with no return value
  * @param realInput an array of real data
  * @param imagInput an array of imag data
  * @param n the number of points in the input arrays
  * @param omegaStart the starting output frequency in [0 2*PI]
  * @param omegaEnd the ending output frequency in [0 2*PI]
  * @param numberOutputPoints the number of output points to calculate the Fourier transform
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void czt(const float *realInput, const float *imagInput, int n, float omegaStart, float omegaEnd,
               int numberOutputPoints);
  

//! fft() calculates the FFT of an array of real data, the input
//! array is not modified
/**
  * a normal member taking two arguments with no return value
  * @param realInput an array of real data
  * @param n the number of points in the array
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void realFFT(const float *realInput, int n);

//! fft() calculates the FFT of an array of real double data, the input
//! array is not modified
/**
  * a normal member taking two arguments with no return value
  * @param realInput an array of real data
  * @param n the number of points in the array
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void realFFT(const double *realInput, int n);

//! realCZT() calculates the chirp Z transform of an array of real data, the input
//! array is not modified. This function is useful, when the number of points is
//! not a power of two, or a custom frequency spacing is desired.
/**
  * a normal member taking six arguments with no return value
  * @param realInput an array of real data
  * @param n the number of points in the input arrays
  * @param omegaStart the starting output frequency in [0 2*PI]
  * @param omegaEnd the ending output frequency in [0 2*PI]
  * @param numberOutputPoints the number of output points to calculate the Fourier transform
  * @see realDoubleOutput(), imagDoubleOutput(), etc for getting results
  * @return void
  */
  void realCZT(const float *realInput, int numberInputPoints, float omegaStart, float omegaEnd,
               int numberOutputPoints);
//
// The following functions return the output from the FFT
//
  const float *realFloatOutput();                       //!< returns the float real part of the FT
  const float *imagFloatOutput();                       //!< returns the float imag part of the FT
  const float *magFloatOutput(int magType);             //!< returns the float magnitude of the FT
  const double *realDoubleOutput();                     //!< returns the double real part of the FT
  const double *imagDoubleOutput();                     //!< returns the double imag part of the FT
  const double *magDoubleOutput(int magType);           //!< returns the double magnitude of the FT
  const Complex *complexOutput() { return _complexArray;} //!< returns the FT as a complex array

};
#endif
