/********************************************************************************
 *                                                                              *
 * This class generates PN data using PN polynomial as described in Ziemer      *
 * and Peterson.                                                                *
 * Notes:                                                                       *
 * 1. Both SSRG and MSRG configurations as in Fig. 8-5,6 in Ziemer and          *
 * Peterson are implemented.                                                    *
 * 2. diagram:                                                                  *
 *              SSRG (W-CDMA):                                                  *
 *    ------- + <----- + <--------                                              *
 *   |        ^        ^          |                                             *
 *   |        | g[n]   | g[n-1]     |g[0]                                       *
 *   -> a[n] -> a[n-1] -> ... -> a[0]                                           *
 *                                                                              *
 *              MSRG (IS-95):                                                   *
 *    ---------------------------------------                                   *
 *   |        |          |         |        |                                   *
 *   |        | g[1]     | g[n-2]  |g[n-1]  |                                   *
 *   -> a[n]->+ a[n-1]-> +-> ... -> + ->a[0]---->                               *
 *                                                                              *
 * NOTE: This class uses the PNGenerator class to do the work, and then saves   *
 *       the output in a buffer for fast storage.                               *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNBuffered.h                                *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "PNGenerator.h"
#include "PNBuffered.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Private Method #################################
// generateNewData -- Generates a new set of PN generator data, and stores in buffer.
// Input:               None
// Output:              None
//
// Notes:
// 1. New PN data is stored in _pnData.
// ############################# Private Method #################################
void PNBuffered::generateNewData()
{
  int   i, pn_length;

  _outputLength = pn_length     = _pnGenerator->numberBits();
  if(_addTrailingZero)
    _outputLength++;
  if(_oldLength < _outputLength)
  {
    delete [] _pnData;
    _oldLength  = _outputLength;
    _pnData     = new char[_outputLength];
  }
  _pnGenerator->initGenerator();
  for(i=0; i<pn_length; i++)
    _pnData[i]          = _pnGenerator->newData();
  if(_addTrailingZero)
    _pnData[pn_length]  = 0;
  return;
}


// ############################# Class Constructor #################################
// PNBuffered -- Constructor for the PNBuffered class
// Input:       poly:   New primitive polynomial given in octal notation
// Output:              None
// ############################# Class Constructor #################################
PNBuffered::PNBuffered(int poly)
{
// 
// Initialize instance variables:
//
  _pnGenerator  = new PNGenerator(poly);
  _pnData       = NULL;
  _oldLength    = 0;
  setPNInitialLoad(DEFAULT_PN_LOAD);
  setPNOffset(0);
  setGeneratorType(MSRG_TYPE);
  setOutputTapNumber(0);
  setAddTrailingZero(0);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// PNBuffered -- Destructor for the PNBuffered class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
PNBuffered::~PNBuffered()
{
//
// Delete space:
//
  delete _pnGenerator;
  delete [] _pnData;
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the shift register stages of the PN generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNBuffered::initGenerator()
{
  _outputPointer        = 0;
  return;
}

// ############################ Public Function ###################################
// setPNPolynomial - Sets a new polynomial in octal, and initializes the shift registers.
// Input:       poly:   New primitive polynomial given in octal notation
// Output:              None
//
// Notes:
// 1. The shift register is initialized to a[0] = 1, a[i] = 0;
// ############################ Public Function ###################################
void PNBuffered::setPNPolynomial(int poly)
{
  _pnGenerator->setPNPolynomial(poly);
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// setPNOffset - Sets the offset of the PN polynomial.
// Input:       offset: PN offset
// Output:              None
//
// Notes:
// 1. The offset is used during initPNBuffered to run the PNBuffered for a number
//    of samples equal to the offset
// ############################ Public Function ###################################
void PNBuffered::setPNOffset(int offset)
{
  _pnGenerator->setPNOffset(offset);
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// setPNInitialLoad - Sets the initial load for the shift registers.
// Input:       load:   Initial PN register load
// Output:              None
//
// Notes:
// 1. The shift register is initialized to a[0] = load&1, a[1] = load&2, etc.;
// ############################ Public Function ###################################
void PNBuffered::setPNInitialLoad(int load)
{
  _pnGenerator->setPNInitialLoad(load);
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// setGeneratorType - Sets the PN generator type, MSRG or SSRG.
// Input:       type:   PN generator type
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNBuffered::setGeneratorType(int type)
{
  _pnGenerator->setGeneratorType(type);
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// setOutputTapNumber - Sets a new shift register tap to take output.
// Input:       number: tap number
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNBuffered::setOutputTapNumber(int number)
{
  _pnGenerator->setOutputTapNumber(number);
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// setAddTrailingZero - Sets the flag for adding a trailing zero.
// Input:       flag:   YES = add trailing zero
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNBuffered::setAddTrailingZero(LOGICAL flag)
{
  _addTrailingZero      = flag;
  generateNewData();
  return;
}

// ############################ Public Function ###################################
// pnPolynomial - Returns the PN polynomial in octal.
// Input:               None
// Output:              PN polynomial
//
// Notes:
// ############################ Public Function ###################################
int PNBuffered::pnPolynomial()
{
  return (_pnGenerator->pnPolynomial());
}

// ############################ Public Function ###################################
// numberBits - Returns number of bits in sequence.
// Input:               None
// Output:              length of sequence
//
// Notes:
// 1. Length of sequence = 2**n - 1, where n = degree of the polynomial.
// ############################ Public Function ###################################
int PNBuffered::numberBits()
{
  if(_addTrailingZero)
    return (_pnGenerator->numberBits() + 1);
  else
    return (_pnGenerator->numberBits());
}

// ############################ Public Function ###################################
// newData - Gets new data from the PN generator.
// Input:               None
// Output:              New data bit (0,1)
//
// Notes:
// 1. The generator type (MSRG or SSRG) needs to be set in setGeneratorType
//
// ############################ Public Function ###################################
int PNBuffered::newData()
{
  int   output_data;

  output_data           = _pnData[_outputPointer];
  if(++_outputPointer >= _outputLength)
    _outputPointer      = 0;
  return output_data;
}
