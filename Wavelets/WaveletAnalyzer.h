#ifndef _WAVELET_ANALYZER_H
#define _WAVELET_ANALYZER_H     1
/************************************************************************
 *                                                                      *
 * This class implements analyzes a set of data using quadrature mirror *
 * filter banks.                                                        *
 *                                                                      *
 * File:WaveletAnalyzer.h                                               *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/23/05 - Started.                                              *
 *                                                                      *
 * References:                                                          *
 *  1. P.E. Pace, Low Probability of Intercept Radar.                   *
 *                                                                      *
 ************************************************************************/

#ifndef         YES
#define         YES             1
#define         NO              0
#define         LOGICAL         char
#endif

const   int DEFAULT_FILTER_LENGTH       = 64;

class   RealFIR;                                // The H and G filters
class   WaveletFilterCoefficients;              // Generates coefficients
class   Complex;

//! WaveletAnalyzer Class
/*!
   The WaveletAnalyzer class can be used for analyzing a set of
   sampled data.
*/
class WaveletAnalyzer                   //!< WaveletAnalyzer class
{
protected:
  int           _filterLength;                  //!< Number of coefficients in H and G
  int           _oldInputSize;                  //!< saves on calls to delete/new

  RealFIR       *_lowPassFilter;                //!< "H" filter
  RealFIR       *_highPassFilter;               //!< "G" filter
  WaveletFilterCoefficients
                *_coefficientGenerator;         //!< Generates coefficients for H and G

  Complex       *_analyzerOutput;               //!< Output data from wavelet analyzer
  Complex       *_analyzerInput;                //!< Input data to analyzer

//
// Private functions:
//
  void  filterData(const Complex *inputData, Complex *outputData, int length, RealFIR *filter);
                                               //!< Filter a set of data, and store in output
  int   updateFilterPointer(char *currentFilter, char *nextFilter, int dataSize);

public:

// 
// Public functions:
//
//! Class constructor
  WaveletAnalyzer();
  virtual ~WaveletAnalyzer();

/****************************************
 * These methods set parameters:        *
 ****************************************/  
//! setFilterLength() allows the setting of new filter length
/**
  * a normal member taking one arguments with no return value
  * @param type the new filter length
  * @return void
  */
  void  setFilterLength(int length);

/**********************
 * Get parameters:    *
 **********************/
  int   filterLength()          {return _filterLength;}                 //!< returns the length of the filter
  const Complex *analyzerOutput() {return _analyzerOutput;}             //!< returns the result of waveletAnalyze
  int   numberStages(int length, int desiredStages);                    //!< returns number of actual stages

/************************
 * Analyze using QMF.   *
 ************************/
//! waveletAnalyze() this method analyzes the input signal using a wavelet filter
//! bank.  The output data is returned as a single array, with the first part of the
//! array containing the output of the lowest frequency filter, and the last part of
//! the array containing the output of the highest frequency filter.
/**
  * a normal member taking three arguments with a return value
  * @param inputData the input data
  * @param inputLength the length of the input data
  * @param numberStages the number of stages of filter banks to perform
  * @return the actual number of stages performed (if inputLength does not divide by 2, numberStages times)
  * @see analyzerOutput() to get the results
  */
  int   waveletAnalyze(const Complex *inputData, int inputLength, int numberStages);

//! waveletAnalyze() this method analyzes the input signal using a wavelet filter
//! bank.  The output data is returned as a single array, with the first part of the
//! array containing the output of the lowest frequency filter, and the last part of
//! the array containing the output of the highest frequency filter.
/**
  * a normal member taking four arguments with a return value
  * @param inputReal the input data (real part)
  * @param inputImag the input data (imag part)
  * @param inputLength the length of the input data
  * @param numberStages the number of stages of filter banks to perform
  * @return the actual number of stages performed (if inputLength does not divide by 2, numberStages times)
  * @see analyzerOutput() to get the results
  */
  int waveletAnalyze(const float *inputReal, const float *inputImag, int inputLength, int numberStages);

};
#endif
