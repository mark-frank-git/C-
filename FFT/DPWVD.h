#ifndef _DPWVD_H
#define _DPWVD_H        1
/********************************************************************************
 *                                                                              *
 * This implements an object for computing the discrete pseudo                  *
 * Wigner Ville distribution.                                                   *
 *                                                                              *
 * File: DPWVD.h                                                                *
 *                                                                              *
 * REFERENCES:                                                                  *
 *  1. [PAC04]:P.E. Pace, Detecting and Classifying Low Probability of          *
 *  Intercept Radar, Artech House, 2004.                                        *
 *                                                                              *
 *                                                                              *
 * Revision history:                                                            *
 *  1. 08/04/05  - Started.                                                     *
 *********************************************************************************/

#ifndef LOGICAL
#define LOGICAL                 char
#endif
#ifndef YES
#define YES                     1
#define NO                      0
#endif


//! DPWVD Class
/*!
   The DPWVD class allows an easy interface to the Wigner Ville routines.
*/
class DPWVD
{
private:
  int           _kernelSize;                            //!< kernel size = N
  int           _inputBufferSize;                       //!< should be equal to 2*N-1
  int           _windowType;                            //!< data windowing type, in dataWindow.h
  int           _readPointer;                           //!< Input buffer read pointer
  int           _writePointer;                          //!< Input buffer write pointer
  double        *_windowArray;                          //!< w(n)*w(-n)
  Complex       *_inputBuffer;                          //!< Ring buffer for input data
  Complex       *_kernelRegister;                       //!< Kernel register buffer

public:
//
// Public methods:
//
//! Class constructor, the sampling frequency is used for frequency axis scaling
  DPWVD(int kernelSize=8, int windowType=0);
//!< Class destructor
  ~DPWVD();

/**********************
 * Set parameters:    *
 **********************/
//! setKernelSize() allows the setting of a new kernel size (should be a power of 2)
/**
  * a normal member taking one argument with no return value
  * @param kernelSize the new kernel size, N
  * @return void
  */
  void setKernelSize(int size);
  
//! setWindowType() allows the setting of a new window type (see DataWindow.h)
/**
  * a normal member taking one argument with no return value
  * @param type RECTANGULAR, TRIANGULAR, etc.
  * @return void
  */
  void setWindowType(int type);

  //! initialize() zeros out the input buffer
/**
  * a normal member taking no argument with no return value
  * @return void
  */
  void initialize();

/**********************
 * Get parameters:    *
 **********************/
  int   kernelSize()            {return _kernelSize;}           //!< returns the kernel size, N
  int   outputSize()            {return _inputBufferSize;}      //!< returns size of spectral data
  int   windowType()            {return _windowType;}           //!< returns the window type

  
//! processInput() puts an input datum into the input buffer, and calculates a new kernel
//! register
/**
  * a normal member taking one argument with no return value
  * @param input a single complex input data
  * @see spectralOutput() for getting an array of spectral values for current input
  * @return void
  */
  void processInput(const Complex &input);
  void goodprocessInput(const Complex &input);
  
//! spectralOutput() calculates and returns the spectral output for the current input
/**
  * a normal member taking no arguments with Complex * return
  * @param input a single complex input data
  * @see processInput() for inputting data into input buffer
  * @return const Complex*
  */
  Complex *spectralOutput();

};
#endif
