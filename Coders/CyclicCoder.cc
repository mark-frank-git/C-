/********************************************************************************
 * This subclass of Object implements a cyclic coder.                           *
 * The generator polynomial for the degree r code is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/22/02 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "CyclicCoder.h"
#include "LFSR.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#ifndef YES
#define YES 1
#define NO  0
#endif


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the CyclicCoder class.
// Input:               degree:         Degree of generator polynomial
//                      poly:           Generator poly, poly[0] = g0, poly[1] = g1
//          
// Output:                          None
//
// Notes:
// g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.
// ############################# Class Constructor ###############################
CyclicCoder::CyclicCoder(int degree, BIT *poly)

{
  _oldOutputSize        = 0;
  _outputSize           = 0;
  _coderOutput          = NULL;
  _lfsr                 = NULL;
  setGeneratorPoly(degree, poly);
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the CyclicCoder class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
CyclicCoder::~CyclicCoder()
{
  delete [] _coderOutput;
  delete _lfsr;
  return;
}

// ############################# Public Function ###############################
// initCoder -- Initializes the coder to start accepting new input.
//
// Input:       registerLoad:   Initial load into shift register
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void CyclicCoder::initCoder(BIT registerLoad)
{
  _lfsr->initLFSR(registerLoad);
  return;
}

// ############################# Public Function ###############################
// setGeneratorPoly -- Sets a new generator polynomial.
//
// Input:       degree:         degree of new polynomial.
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void CyclicCoder::setGeneratorPoly(int degree, const BIT *poly)
{
  if(_lfsr == NULL)
    _lfsr       = new LFSR(degree, poly);
  else
    _lfsr->setGeneratorPoly(degree, poly);
  return;
}

// ############################# Public Function ###############################
// blockOutput -- Returns a block of coded data using systematic coding.
//
// Input:       inputBits       set of input bits
//              inputLength     length of the input bits array
//              flipParity      flip the parity bits (seems to be used in GSM)
//          
// Output:                      pointer to the output array
//
// Notes:
// 1. systematic coder
// 2. Call outputSize() to get length of output array
// ############################# Public Function ###############################
const BIT *CyclicCoder::blockOutput(const BIT *inputBits, int inputLength, LOGICAL flipParity)
{
  int   i, k, degree;
  const BIT *shift_register;
//
// Allocate output array:
//
  degree        = _lfsr->degree();
  _outputSize   = inputLength + degree;         // input plus parity
  if(_oldOutputSize < _outputSize)
  {
    _oldOutputSize      = _outputSize;
    delete [] _coderOutput;
    _coderOutput        = new BIT[_outputSize];
  }
//
// Send the input through the LFSR:
//
  k                     = 0;
  for(i=0; i<inputLength; i++)
  {
    _lfsr->coderOutput(inputBits[i]);
    _coderOutput[k++]   = inputBits[i];
  }
//
// Get the parity bits:
//
  shift_register        = _lfsr->shiftRegister();
  for(i=0; i< degree; i++)
  {
    if(flipParity)
    {
      _coderOutput[k++] = (shift_register [degree-i-1] == 1) ? 0 : 1;
    }
    else
    {
      _coderOutput[k++] = shift_register [degree-i-1];
    }
  }
  return _coderOutput;
}

// ############################# Public Function ###############################
// bitErrors -- Returns YES if there are bit errors in the received code word.
//
// Input:       inputBits       set of input bits
//              inputLength     length of the input bits array
//              flipParity      flip the parity bits (seems to be used in GSM)
//          
// Output:                      YES if errors, else NO
//
// Notes:
// 1. systematic coder
// ############################# Public Function ###############################
LOGICAL CyclicCoder::bitErrors(const BIT *inputBits, int inputLength, LOGICAL flipParity)
{
  int   i, k, degree;
  const BIT *shift_register;
//
// Send the input through the LFSR:
//
  k                     = 0;
  for(i=0; i<inputLength; i++)
  {
    _lfsr->syndromeOutput(inputBits[i]);
  }
//
// Check the syndrome:
//
  shift_register        = _lfsr->shiftRegister();
  degree                = _lfsr->degree();
  if(flipParity)
  {
    for(i=0; i< degree; i++)
    {
      if(shift_register[i] == 0)
      {
        return YES;
      }
    }
  }
  else
  {
    for(i=0; i< degree; i++)
    {
      if(shift_register[i] == 1)
      {
        return YES;
      }
    }
  }
  return NO;
}
