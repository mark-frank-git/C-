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
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNGenerator.cc                              *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "PNGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################ Private Function ###################################
// newMSRGData - Gets new data from the PN generator.
// Input:               None
// Output:              New data bit (0,1)
//
// Notes:
//              MSRG:
//    ---------------------------------------
//   |        |          |         |        |
//   |        | g[1]     | g[n-2]  |g[n-1]  |
//   -> a[n]->+ a[n-1]-> +-> ... -> + ->a[0]---->
//
// ############################ Private Function ###################################
int PNGenerator::newMSRGData()
{
  int   i, output;
  int   rtn_value;
//
// Linear feedback shift register with MSRG configuration.  For example, with
// PN = 45(octal), g[0]=1, g[1]=0, g[2]=1, g[3]=0, g[4]=0, g[5]=1
//
  output        = _a[0];
  rtn_value     = _a[_outputTapNumber];
  for(i=1; i<_polynomialDegree; i++)
  {
    if(_g[_polynomialDegree-i])
      _a[i-1]   = _a[i] ^ output;
    else
      _a[i-1]   = _a[i];
  }
  _a[_polynomialDegree-1]       = output;
  return rtn_value;
}

// ############################ Private Function ###################################
// newSSRGData - Gets new data from the PN generator.
// Input:               None
// Output:              New data bit (0,1)
//
// Notes:
// 1. SSRG configuration is in Fig. 8-6 in Ziemer and Peterson
// 2. diagram:
//
//    ------- + <----- + <--------
//   |        ^        ^          |
//   |        | g[n]   | g[n-1]     |g[0]
//   -> a[n] -> a[n-1] -> ... -> a[0]
//
// 3. OUTPUT is taken from a[0] as in W-CDMA, NOT LIKE Ziemer and Peterson.
// 4. g[] order is BACKWARDS from Ziemer and Peterson
// ############################ Private Function ###################################
int PNGenerator::newSSRGData()
{
  int   i, sum, output;
//
// Linear feedback shift register with SSRG configuration.  For example, with
// PN = 45(octal), g[0]=1, g[1]=0, g[2]=1, g[3]=0, g[4]=0, g[5]=1
//
  sum           = _a[0];                                        // g[polynomialDegree] always == 1
  output        = _a[_outputTapNumber];
  for(i=1; i<_polynomialDegree; i++)
  {
    if(_a[i])
      sum       ^= _g[i];
    _a[i-1]     = _a[i];
  }
  _a[_polynomialDegree-1]       = sum;
  return output;

}
// ############################# Class Constructor #################################
// PNGenerator -- Constructor for the PNGenerator class
// Input:       poly:   New primitive polynomial given in octal notation
// Output:              None
// ############################# Class Constructor #################################
PNGenerator::PNGenerator(int poly)
{
// 
// Initialize instance variables:
//
  _a                    = _g                    = NULL;
  setPNInitialLoad(DEFAULT_PN_LOAD);
  setPNPolynomial(poly);
  setPNOffset(0);
  setGeneratorType(MSRG_TYPE);
  setOutputTapNumber(0);

  return;
}

// ############################# Class Destructor ###############################
// PNGenerator -- Destructor for the PNGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
PNGenerator::~PNGenerator()
{
//
// Delete space:
//
  delete [] _a;
  delete [] _g;
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the shift register stages of the PN generator.
// Input:               None
// Output:              None
//
// Notes:
// 1. The shift register is initialized to a[0] = _pnInitialLoad&1,
//    a[1]= _pnInitialLoad&2, a[2] = _pnInitialLoad&4, etc.
// ############################ Public Function ###################################
void PNGenerator::initGenerator()
{
  int   i, load;
//
// Set the initial load:
//
  load          = _pnInitialLoad;
  for(i=0; i <= _polynomialDegree; i++)
  {
    _a[i]       = 1 & load;
    load        >>= 1;
  }
//
// run until offset samples:
//
  for(i=0; i<_pnOffset; i++)
    newData();
  return;
}

#define OCTAL_BITS  3
// ############################ Public Function ###################################
// setPNPolynomial - Sets a new polynomial in octal, and initializes the shift registers.
// Input:       poly:   New primitive polynomial given in octal notation
// Output:              None
//
// Notes:
// 1. The shift register is initialized to a[0] = 1, a[i] = 0;
// ############################ Public Function ###################################
void PNGenerator::setPNPolynomial(int poly)
{
  int i, pn_binary, number_bits;
  int pn_poly, octal_power, octal_digit;

// First convert octal representation to binary
  pn_poly       = _pnPolynomial = poly;
  pn_binary     = number_bits   = 0;
  octal_power   = 1;
  while(pn_poly/10)
  {
    octal_digit = pn_poly % 10;
    pn_binary   += octal_digit*octal_power;
    octal_power *= 8;
    pn_poly     /= 10;
    number_bits += OCTAL_BITS;
  }
  pn_binary     += pn_poly*octal_power;  
  while(pn_poly)
  {
    number_bits++;
    pn_poly >>= 1;
  }
// Save binary representation in the g polynomial:
  _polynomialDegree     = MAX(1, number_bits-1);
  delete [] _g;
  delete [] _a;
  _g                    = new short[(_polynomialDegree+1)];
  _a                    = new short[(_polynomialDegree+1)];
  for(i=0; i<=_polynomialDegree; i++)
  {
    _g[i]               = pn_binary & 1;
    _a[i]               = 0;
    pn_binary >>= 1;
  }
  _a[0] = 1;

//
// Error check _outputTapNumber:
//
  _outputTapNumber      = MIN(_polynomialDegree, _outputTapNumber);

  return;
}

// ############################ Public Function ###################################
// setPNOffset - Sets the offset of the PN polynomial.
// Input:       offset: PN offset
// Output:              None
//
// Notes:
// 1. The offset is used during initPNGenerator to run the pnGenerator for a number
//    of samples equal to the offset
// ############################ Public Function ###################################
void PNGenerator::setPNOffset(int offset)
{
  _pnOffset     = MAX(0, offset);
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
void PNGenerator::setPNInitialLoad(int load)
{
  _pnInitialLoad        = MAX(1, load);
  return;
}

// ############################ Public Function ###################################
// setGeneratorType - Sets the PN generator type, MSRG or SSRG.
// Input:       type:   PN generator type
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNGenerator::setGeneratorType(int type)
{
  _generatorType        = type;
  return;
}

// ############################ Public Function ###################################
// setOutputTapNumber - Sets a new shift register tap to take output.
// Input:       number: tap number
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void PNGenerator::setOutputTapNumber(int number)
{
  _outputTapNumber      = MAX(0, number);
  _outputTapNumber      = MIN(_polynomialDegree, _outputTapNumber);
  return;
}

// ############################ Public Function ###################################
// numberBits - Returns number of bits in sequence.
// Input:               None
// Output:              length of sequence
//
// Notes:
// 1. Length of sequence = 2**n - 1, where n = degree of the polynomial.
// ############################ Public Function ###################################
int PNGenerator::numberBits()
{
  long  bits, stages;

  bits          = 2;
  stages        = _polynomialDegree;
  while(--stages > 0)
    bits        <<= 1;
  return bits - 1;
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
int PNGenerator::newData()
{
  if(_generatorType == MSRG_TYPE)
    return newMSRGData();
  else
    return newSSRGData();
}
