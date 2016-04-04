#ifndef _WAVELET_FILTER_COEFFICIENTS_H
#define _WAVELET_FILTER_COEFFICIENTS_H  1
/************************************************************************
 *                                                                      *
 * This class implements generates coefficients for wavelet filters     *
 *                                                                      *
 * File:WaveletFilterCoefficients.h                                     *
 *                                                                      *
 *                                                                      *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/09/05 - Started.                                              *
 *                                                                      *
 ************************************************************************/

#ifndef         YES
#define         YES             1
#define         NO              0
#define         LOGICAL         char
#endif


//! WaveletFilterCoefficients Class
/*!
   The WaveletFilterCoefficients class can be used for designing coefficients
   for wavelet filters.
*/
class WaveletFilterCoefficients                 //!< WaveletFilterCoefficients class
{
protected:
  int           _filterType;                    //!< Type of wavelet filter
  int           _windowType;                    //!< Type of windowing for sinc type filter
  int           _filterLength;                  //!< Number of coefficients

  double        _compressionVariable;           //!< Sinc filter parameter
  double        _scalingVariable;               //!< Sinc filter parameter
  float         *_lowPassCoefficients;          //!< The low pass filter coefficients
  float         *_highPassCoefficients;         //!< The high pass filter coefficients

//
// Private functions:
//
  double        sinc(double x);                 //!< Calculates the sinc function
  
public:
//! Filter types enum
/*! The filter_type_enum categories the different type of wavelet filter types */
enum filter_type_enum
{
  SINC_WAVELET                          //!< Sinc type of wavelet 
};

// 
// Public functions:
//
//! Class constructor
  WaveletFilterCoefficients();
  virtual ~WaveletFilterCoefficients();

/****************************************
 * These methods set parameters:        *
 ****************************************/
//! setFilterType() allows the setting of new filter type
/**
  * a normal member taking one arguments with no return value
  * @param type the new filter type
  * @return void
  */
  void  setFilterType(int type);

//! setWindowType() allows the setting of new window type for sinc filter
/**
  * a normal member taking one arguments with no return value
  * @param type the new window type
  * @return void
  */
  void  setWindowType(int type);
  
//! setFilterLength() allows the setting of new filter length
/**
  * a normal member taking one arguments with no return value
  * @param type the new filter length
  * @return void
  */
  void  setFilterLength(int length);
  
//! setCompressionVariable() allows the setting of new compression variable
/**
  * a normal member taking one arguments with no return value
  * @param compression the new compression value
  * @return void
  */
  void  setCompressionVariable(double compression);
  
//! setScalingVariable() allows the setting of new scaling variable
/**
  * a normal member taking one arguments with no return value
  * @param scaling the new scaling value
  * @return void
  */
  void  setScalingVariable(double scaling);

/**********************
 * Get parameters:    *
 **********************/
  int           filterType()                    {return _filterType;}           //!< returns the type of the filter
  int           windowType()                    {return _windowType;}           //!< returns the type of the window
  int           filterLength()                  {return _filterLength;}         //!< returns the length of the filter

/************************
 * Calculate and return *
 * filter coefficients: *
 ************************/
  const         float *calculateLowPassCoefficients();                          //!< calculates and returns LPF coeffs
  const         float *calculateHighPassCoefficients();                         //!< calculates and returns LPF coeffs


};
#endif
