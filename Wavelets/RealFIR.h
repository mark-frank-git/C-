#ifndef _REAL_FIR_H
#define _REAL_FIR_H     1
/************************************************************************
 *                                                                      *
 * This class implements a digital FIR filter having real coefficients  *
 *                                                                      *
 * File:RealFIR.h                                                       *
 *                                                                      *
 * The filter is stored in the forms:                                   *
 *                                                                      *
 *    H(z) =  b[0] + b[1]z^(-1) + ... + b[n]z^-(n)                      *
 *                                                                      *
 * Where the b coefficients are calculated outside of this class.       *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/09/05 - Started, modified from AbstractFilter.                *
 ************************************************************************/


class Complex;                                  // class prototype

#ifndef         YES
#define         YES             1
#define         NO              0
#endif
#define         MIN_TAPS        2
#define         MAX_TAPS        50000

//! RealFIR Class
/*!
   The RealFIR class can be used for implementing an FIR filter
   with real taps.  The filter coefficients need to be calculated
   outside of this class.
*/
class RealFIR                                   //!< RealFIR class
{
protected:
  int           _numberTaps;                    //!< Number of filter coefficients
  int           _shiftPtr;                      //!< Shift buffer pointer
  int           _shiftSize;                     //!< Size of shift register = numberTaps-1
  int           _numberRealOutput;              //!< Size of output array
  int           _oldRealOutput;                 //!< Old size to save on calls to delete/new
  int           _numberComplexOutput;           //!< Size of complex output array
  int           _oldComplexOutput;              //!< Old size to save on calls to delete/new

  float         *_filterCoefficients;           //!< The filter coefficients
  float         *_realShiftBuffer;              //!< Real circular shift register
  Complex       *_complexShiftBuffer;           //!< Complex circular shift register
  float         *_realOutputData;               //!< Filtered real data
  Complex       *_complexOutputData;            //!< Filtered complex data

//
// Private functions:
//

// 
// Public functions:
//
public:
//! Class constructor, the length is the filter length
  RealFIR(const float *coeffs=NULL, int length = 0);
  virtual ~RealFIR();

/****************************************
 * These methods set parameters:        *
 ****************************************/
//! setFilterCoefficients() allows the setting of new filter coefficients
/**
  * a normal member taking two arguments with no return value
  * @param coeffs the new filter coefficients
  * @param length the number of coefficients
  * @return void
  */
  void  setFilterCoefficients(const float *coeffs, int length);
  
//! reset() resets the filter by zeroing out the taps
/**
  * a normal member taking no arguments with no return value
  * @return void
  */
  void  reset();

/**********************
 * Get parameters:    *
 **********************/
  int           numberTaps()                    {return _numberTaps;}           //!< returns the length of the filter
  const         float *filterCoefficients()     {return _filterCoefficients;}   //!< returns the filter coefficients
  const         float *realOutputData()         {return _realOutputData;}       //!< returns the filtered data
  const         Complex *complexOutputData()    {return _complexOutputData;}    //!< returns the filtered complex data
        
/****************************************
 * Filtering float and Complex data     *
 * assuming an FIR structure.           *
 ****************************************/
//! filterFloatArray() filters the input array, and stores result in _realOutputData
/**
  * a normal member taking two arguments with no return value
  * @param input the input data
  * @param length the number of data points
  * @see realOutputData() for getting results
  * @return void
  */
  void          filterFloatArray(const float * input, int numberPts);
  
//! filterFloatData() filters a single input data, and returns the filtered output
/**
  * a normal member taking one argument with a return value
  * @param input the input data
  * @return filtered datum
  */
  float         filterFloatData(float input);                   // Filter float data

//! filterComplexArray() filters the input array, and stores result in _complexOutputData
/**
  * a normal member taking two arguments with no return value
  * @param input the input data
  * @param length the number of data points
  * @see complexOutputData() for getting results
  * @return void
  */
  void          filterComplexArray(const Complex *x, int numberPts);
  
//! filterComplexData() filters a single complex input datum, and returns the filtered output
/**
  * a normal member taking one argument with a return value
  * @param input the input data
  * @return filtered datum
  */
  Complex       filterComplexData(const Complex &input);                // Filter complex data


};
#endif
