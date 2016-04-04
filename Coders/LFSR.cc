/********************************************************************************
 * This subclass of Object implements a linear feedback shift register.         *
 * The generator polynomial for the degree r LFSR is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/02 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
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
// Class Constructor -- Constructor for the LFSR class.
// Input:               degree:         Degree of generator polynomial
//                      poly:           Generator poly, poly[0] = g0, poly[1] = g1
//          
// Output:                          None
//
// Notes:
// g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.
// ############################# Class Constructor ###############################
LFSR::LFSR(int degree, const BIT *poly)

{
  _oldDegree            = 0;
  _oldOutputSize        = 0;
  _oldShiftStages       = 0;
  _degree               = 0;
  _outputSize           = 0;
  _generatorPoly        = NULL;
  _lfsrOutput           = NULL;
  _shiftRegister        = NULL;
  if(poly != NULL)
    setGeneratorPoly(degree, poly);
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the LFSR class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
LFSR::~LFSR()
{
  delete [] _generatorPoly;
  delete [] _lfsrOutput;
  delete [] _shiftRegister;
  return;
}

// ############################# Public Function ###############################
// initLFSR -- Initializes the LFSR to start accepting new input.
//
// Input:       registerLoad:   Initial load into shift register
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void LFSR::initLFSR(BIT registerLoad)
{
//
// Add this:
//
  if(_oldShiftStages < _degree)
  {
    _oldShiftStages     = _degree;
    delete [] _shiftRegister;
    _shiftRegister      = new BIT[_degree];
  }
  for(int i=0; i<_degree; i++)
    _shiftRegister[i]   = registerLoad;
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
void LFSR::setGeneratorPoly(int degree, const BIT *poly)
{
  int   i;

  _degree                       = MIN(degree, MAX_POLY_DEGREE);
  if(_oldDegree < _degree)
  {
    _oldDegree                  = _degree;
    delete [] _generatorPoly;
    _generatorPoly              = new BIT[_degree];
  }
  for(i=0; i< _degree; i++)
  {
    _generatorPoly[i]   = poly[i];
  }
  return;
}

// ############################# Public Function ###############################
// coderOutput -- Returns data out of the LFSR from a new data input, assuming
//                a systematic coder with the data fed into an exclusive-or at
//                the output of the last register
//
// Input:       inputBit        new input bit
//          
// Output:                      output bit
//
// Notes:
// 1. systematic LFSR, as in Figure 5-12 in Wicker's 1995 book on Error Control
//    Systems.
// ############################# Public Function ###############################
BIT LFSR:: coderOutput(BIT inputBit)
{
  int   i;
  BIT   output, shift_input, shift_output;
//
// Get the output
//
  output        = _shiftRegister[_degree-1] ^ inputBit;
//
// Feed back the output
//
  shift_output  = 0;
  for(i=0; i<_degree; i++)
  {
    if(_generatorPoly[i] == 1)
    {
      shift_input       = output;
    }
    else
    {
      shift_input       = 0;
    }
    shift_input         ^= shift_output;
    shift_output        = _shiftRegister[i];
    _shiftRegister[i]   = shift_input;
  }
  return        output;
}

// ############################# Public Function ###############################
// syndromeOutput -- Returns data out of the syndrome circuit from a new data input.
//                   The data is fed into an exclusive or circuit which XORs the
//                   output from the shift register.  The output of the XOR then
//                   feeds into the first stage of the shift register.
//
// Input:       inputBit        new input bit
//          
// Output:                      output bit
//
// Notes:
// 1. systematic LFSR
// ############################# Public Function ###############################
BIT LFSR:: syndromeOutput(BIT inputBit)
{
  int   i;
  BIT   output, shift_input, shift_output;
//
// Get the output
//
  output        = _shiftRegister[_degree-1];
//
// Feed back the output
//
  shift_output  = inputBit;
  for(i=0; i<_degree; i++)
  {
    if(_generatorPoly[i] == 1)
    {
      shift_input       = output;
    }
    else
    {
      shift_input       = 0;
    }
    shift_input         ^= shift_output;
    shift_output        = _shiftRegister[i];
    _shiftRegister[i]   = shift_input;
  }
  return        output;
}
