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

/**************************
 * Include files:         *
 **************************/
#include <math.h>
#include <stdio.h>
#include <GNU/Complex.h>
#include <Filters/DataWindow.h>
#include "DPWVD.h"
#include "fft_routines.h"


// ############################# Class Constructor #################################
// DPWVD -- Constructor for the DPWVD class
// Input:       kernelSize:     size of kernel, N
//              windowType:     type of window function, see DataWindow.h
// Output:                      None
//
// ############################# Class Constructor #################################
DPWVD::DPWVD(int kernelSize, int windowType)
{
//
// Initialize arrays to NULL
//
  _windowArray          = NULL;
  _inputBuffer          = NULL;
  _kernelRegister       = NULL;
//
// Set parameters:
//
  _windowType           = windowType;
  setKernelSize(kernelSize);
//
// reset
//
  initialize();

  return;
}

// ############################# Class Destructor #################################
// DPWVD -- Destructor for the DPWVD class
// Input:               None
// Output:              None
//
// ############################# Class Destructor #################################
DPWVD::~DPWVD()
{
  delete [] _windowArray;
  delete [] _inputBuffer;
  delete [] _kernelRegister;

  return;
}

// ############################ Public Function ###################################
// setKernelSize - Sets a new kernel size, N (should be a power of 2)
//
// Input:       kernelSize:     new kernel size
// Output:                      None
//
// ############################ Public Function ###################################
void DPWVD:: setKernelSize(int kernelSize)
{
//
// Make kernel size a power of two
//
  _kernelSize           = kernelSize;
  while(!powerOfTwo(_kernelSize))
    _kernelSize++;
  _inputBufferSize      = 2*_kernelSize;
//
// Allocate arrays for new sizes
// Note that the kernel register only needs to be kernelSize
// long, but we make it as long as the input buffer size, to
// facilitate taking the FFT
//
  delete [] _inputBuffer;
  delete [] _kernelRegister;

  _inputBuffer                  = new Complex[_inputBufferSize];
  _kernelRegister               = new Complex[_inputBufferSize];
//
// Initialize the new window array:

  setWindowType(_windowType);

  return;
}

// ############################ Public Function ###################################
// setWindowType - Sets a new window type (see DataWindow.h).
//
// Input:       type:       new window type
// Output:                  None
//
// Notes:
// 1. The _windowArray, w(n)*w(n) is modified.
// ############################ Public Function ###################################
void DPWVD:: setWindowType(int type)
{
//
// Allocate array:
//
  delete []     _windowArray;
  _windowArray  = new double[_kernelSize];
//
// Get window
//
  DataWindow    data_window(type);
  const double *window  = data_window.windowFunction(_inputBufferSize);
//
// Create window array = w(n)*w(-n)
//
  _windowArray[0]       = 1.;
  int minus_n           = _inputBufferSize/2 - 1;
  int plus_n            = _inputBufferSize/2 + 1;
  for(int i=1; i<_kernelSize; i++)
  {
    _windowArray[i]     = window[minus_n--]*window[plus_n++];
  }
  return;
}

// ############################ Public Function ###################################
// initialize - Initialize the buffers for new data.
//
// Input:                   None
// Output:                  None
//
// Notes:
// 1. _readPointer, _writePointer, and _inputBuffer are modified.
// ############################ Public Function ###################################
void DPWVD::initialize()
{
  _readPointer          = _writePointer = 0;
  for(int i=0; i<_inputBufferSize; i++)
  {
    _inputBuffer[i]     = Complex(0., 0.);
    _kernelRegister[i]  = Complex(0., 0.);
  }
  return;  
} 

// ############################ Public Function ###################################
// processInput - Processes a new input datum: 1) adds to input buffer, 2) updates
//                kernel register.
//
// Input:       input           New complex input datum
// Output:                      None
//
// Notes:
// 1. Call spectralOutput after each input data is processed.
// 2. This function implements the computational structure in Figure 9.1 of Pace'
//    book.  Note that we actual generate fl' in order to simplify call to FFT.
// ############################ Public Function ###################################
void DPWVD::processInput(const Complex &input)
{
//
// Add new datum to buffer:
//
  _inputBuffer[_writePointer]   = input;
//
// Process data into kernel register.  Start with
// the oldest and newest data, and work inward:
//
  int   recent_input_ptr        = _writePointer;
  int   past_input_ptr          = _writePointer + 1;
  if(past_input_ptr == _inputBufferSize)
    past_input_ptr              = 0;
  int   lower_kernel_ptr        = _kernelSize - 1;
  int   upper_kernel_ptr        = _kernelSize + 1;
  for(int i=0; i<_kernelSize-1; i++)
  {
    int                 j       = lower_kernel_ptr;
    _kernelRegister[j]          = _windowArray[j]*_inputBuffer[recent_input_ptr] *
                                   conj(_inputBuffer[past_input_ptr]);
    int                 k       = upper_kernel_ptr;
    _kernelRegister[k]          = conj(_kernelRegister[j]);
//
// update the indexes
//
    lower_kernel_ptr--;
    upper_kernel_ptr++;
    past_input_ptr++;
    recent_input_ptr--;
    if(past_input_ptr == _inputBufferSize)
      past_input_ptr    = 0;
    if(recent_input_ptr < 0)
      recent_input_ptr  = _inputBufferSize - 1;
  }
//
// Fill in fl(0)
//
  _kernelRegister[0]    = _inputBuffer[past_input_ptr]*_windowArray[lower_kernel_ptr];
//
// Set middle value of kernel register to zero for FFT:
//
  _kernelRegister[_kernelSize]  = Complex(0., 0.);
//
// update write pointer:
//
  _writePointer++;
  if(_writePointer == _inputBufferSize)
    _writePointer               = 0;
  
  return;  
} 


// ############################ Public Function ###################################
// spectralOutput - Returns the spectral output of the current input by running DFT
//                  on _kernelRegister.
//
// Input:                       None
// Output:                      None
//
// Notes:
// 1. Call processInput before calling this routine.
// ############################ Public Function ###################################
Complex *DPWVD::spectralOutput()
{
//
// Find the FFT of the kernel register
//
  int   m = log2(_inputBufferSize);
  complex_fft(_kernelRegister, m, FORWARD_FFT);
//
// Return the result:
//
  
  return _kernelRegister;  
} 
