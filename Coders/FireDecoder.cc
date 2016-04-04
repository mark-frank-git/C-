/********************************************************************************
 * This subclass of Object implements a decoder for a Fire code.                *
 * The generator polynomial for the degree r code is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/02 - Started                                                      *
 *                                                                              *
 * Notes:                                                                       *
 * 1. Based on the algorithm in Wicker's channel coding book.                   *
 ********************************************************************************/
#include "FireDecoder.h"
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
// Class Constructor -- Constructor for the FireDecoder class.
// Input:       burstLength:            Length of the correctable burst, b
//              generatorPoly:          LFSR2 generator poly, poly[0] = g0, poly[1] = g1
//              degree:                 Degree of generatorPoly
//          
// Output:                              None
//
// Notes:
// g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.
// ############################# Class Constructor ###############################
FireDecoder::FireDecoder(int burstLength, const BIT *generatorPoly, int degree)

{
  _oldOutputSize        = 0;
  _outputSize           = 0;
  _decoderOutput        = NULL;
  burstLength           = MAX(1, burstLength);
//
// Allocate the LFSRs
//
  _lfsr1                = new LFSR();
  _lfsr2                = new LFSR();
  if(generatorPoly != NULL)
    setGeneratorPoly2(degree, generatorPoly);
  setBurstLength(burstLength);
  initDecoder();
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the FireDecoder class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
FireDecoder::~FireDecoder()
{
  delete _lfsr1;
  delete _lfsr2;
  delete [] _decoderOutput;
  return;
}

// ############################# Public Function ###############################
// initDecoder -- Initializes the coder to start accepting new input.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void FireDecoder::initDecoder()
{
  if(_lfsr1 != NULL)
    _lfsr1->initLFSR();
  if(_lfsr2 != NULL)
    _lfsr2->initLFSR();
  return;
}

// ############################# Public Function ###############################
// setBurstLength -- Sets the length of the burst of errors that can be corrected.
//
// Input:       burstLength:    length of correctable burst.
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void FireDecoder::setBurstLength(int burstLength)
{
  int   i, poly_degree;
  BIT   *polynomial;

  _burstLength  = MAX(1, burstLength);
  _burstLength  = MIN(MAX_BURST_LENGTH, _burstLength);

//
// Set the single feedback LFSR:
//
  poly_degree   = 2*_burstLength - 1;
  polynomial    = new BIT[poly_degree+1];
  for(i=1; i<poly_degree; i++)
  {
    polynomial[i]       = 0;
  }
  polynomial[0]         = polynomial[poly_degree]       = 1;
  _lfsr1->setGeneratorPoly(poly_degree, polynomial);

  delete [] polynomial;
  return;
}

// ############################# Public Function ###############################
// setGeneratorPoly2 -- Sets a new generator polynomial for LFSR2
//
// Input:       degree:         degree of new polynomial.
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void FireDecoder::setGeneratorPoly2(int degree, const BIT *poly)
{
  _lfsr2->setGeneratorPoly(degree, poly);
  return;
}

// ############################# Public Function ###############################
// shiftInBit -- Shift a received bit into the LFSRs.
//
// Input:       inputBit        new input bit
//          
// Output:                      None
//
// Notes:
// 1. systematic coder
// ############################# Public Function ###############################
void FireDecoder:: shiftInBit(BIT inputBit)
{
  _lfsr1->syndromeOutput(inputBit);
  _lfsr2->syndromeOutput(inputBit);
  return;
}

// ############################# Public Function ###############################
// decodeBits -- Decodes a set of received bits.
//
// Input:       inputBits       set of input bits
//              inputLength     length of the input bits array
//          
// Output:                      pointer to the output array
//
// Notes:
// 1. Call decodeStatus(), decoderOutput(), and outputSize() to get the decoder
//    output.
// ############################# Public Function ###############################
void FireDecoder::decodeBits(const BIT *inputBits, int inputLength)
{
    int i;
//
// Allocate output array:
//
  _outputSize   = inputLength;                  // input plus parity
  if(_oldOutputSize < _outputSize)
  {
    _oldOutputSize      = _outputSize;
    delete [] _decoderOutput;
    _decoderOutput      = new BIT[_outputSize];
  }
//
// Send the data through the LFSR
//
  for(i=0; i<inputLength; i++)
  {
    shiftInBit(inputBits[i]);
  }
  return;
}










  
  return _decoderOutput;
}
