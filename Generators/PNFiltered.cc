/********************************************************************************
 *                                                                              *
 * This class generates PN data using PN polynomial as described in Ziemer      *
 * and Peterson.                                                                *
 * Notes:                                                                       *
 * 1. SSRG configuration is in Fig. 8-6 in Ziemer and Peterson                  *
 * 2. diagram:                                                                  *
 *                                                                              *
 *    ------- + <----- + <--------                                              *
 *   |        ^        ^          |                                             *
 *   |        | g[1]   | g[2]     |                                             *
 *   -> a[n] -> a[n-1] -> ... -> a[0]                                           *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNFiltered.cc                               *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 02/06/01 - Started.                                                      *
 ********************************************************************************/
#include <stdio.h>
#include <C_Libraries/constants.h>
#include "PNGenerator.h"
#include "PNFiltered.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define BIT_MAP(a)      ( ((a)==0 ) ? 1 : -1 )          // (0,1) -> (1, -1)

// ############################# Private Function #################################
// generateData -- Generates a new set of filtered bits data.  The bits just indicate
//                 whether the data is +1, -1, or sloping in between.  It is up to
//                 the getEarlyLateOutput() function to interpret this data.
//
// Input:               None
// Output:              None
//
// Notes:
// 1. The instance variables, _numberBits and _filteredBits are modified.
// 2. The filtered bits are for the late arm of the early/late gate DLL.
// ############################# Private Function #################################
void PNFiltered::generateData()
{
  int   i, bit, old_bit, first_bit;
  
  _numberBits           = _pnGenerator->numberBits();
  delete [] _filteredBits;
  _filteredBits         = new char[_numberBits];
  old_bit               = first_bit     = _pnGenerator->newData();
  bit                   = 0;                            // keep compiler happy
  for(i=1; i<_numberBits; i++)
  {
    bit                 = _pnGenerator->newData();
    if(bit == old_bit)
    {
      if(bit == 0)
        _filteredBits[i]        = PLUS1_BIT;
      else
        _filteredBits[i]        = MINUS1_BIT;
    }
    else
    {
      if(old_bit == 0)
        _filteredBits[i]        = MINUS_SLOPE_BIT;
      else
        _filteredBits[i]        = PLUS_SLOPE_BIT;
    }
    old_bit                     = bit;
  }
//
// Take care of the last bit, wrap around to the first:
//
  i     = 0;
  if(first_bit == bit)
  {
    if(bit == 0)
      _filteredBits[i]          = PLUS1_BIT;
    else
      _filteredBits[i]          = MINUS1_BIT;
  }
  else
  {
    if(first_bit == 0)
      _filteredBits[i]          = MINUS_SLOPE_BIT;
    else
      _filteredBits[i]          = PLUS_SLOPE_BIT;
  }
  return;
}

// ############################# Private Function #################################
// generateData -- Generates a new set of filtered bits data.  The bits just indicate
//                 whether the data is +1, -1, or sloping in between.  It is up to
//                 the getEarlyLateOutput() function to interpret this data.
//
// Input:       pnSequence:     PN sequence mapped to (-1, 1)
//              length:         Length of sequence
// Output:                      None
//
// Notes:
// 1. The instance variables, _numberBits and _filteredBits are modified.
// 2. The filtered bits are for the late arm of the early/late gate DLL.
// ############################# Private Function #################################
void PNFiltered::generateData(const char *pnSequence, int length)
{
  int   i, bit, old_bit, first_bit;
  
  _numberBits           = length;
  delete [] _filteredBits;
  _filteredBits         = new char[_numberBits];
  old_bit               = first_bit     = pnSequence[0];
  bit                   = 0;                            // keep compiler happy
  for(i=1; i<_numberBits; i++)
  {
    bit                 = pnSequence[i];
    if(bit == old_bit)
    {
      if(bit == 1)
        _filteredBits[i]        = PLUS1_BIT;
      else
        _filteredBits[i]        = MINUS1_BIT;
    }
    else
    {
      if(old_bit == 1)
        _filteredBits[i]        = MINUS_SLOPE_BIT;
      else
        _filteredBits[i]        = PLUS_SLOPE_BIT;
    }
    old_bit                     = bit;
  }
//
// Take care of the last bit, wrap around to the first:
//
  i     = 0;
  if(first_bit == bit)
  {
    if(bit == 1)
      _filteredBits[i]          = PLUS1_BIT;
    else
      _filteredBits[i]          = MINUS1_BIT;
  }
  else
  {
    if(first_bit == 1)
      _filteredBits[i]          = MINUS_SLOPE_BIT;
    else
      _filteredBits[i]          = PLUS_SLOPE_BIT;
  }
  return;
}

// ############################# Class Constructor #################################
// PNFiltered -- Constructor for the PNFiltered class
// Input:       poly:   New primitive polynomial given in octal notation
// Output:              None
// ############################# Class Constructor #################################
PNFiltered::PNFiltered(int poly)
{
// 
// Allocate objects
//
  _pnGenerator          = new PNGenerator(poly);
  _filteredBits         = NULL;
//
// Generate bits based on PN generator:
//
  generateData();
  return;
}

// ############################# Class Constructor #################################
// PNFiltered -- Constructor for the PNFiltered class
//
// Input:       pnSequence:     PN sequence mapped to (-1, 1)
//              length:         Length of sequence
// Output:                      None
// ############################# Class Constructor #################################
PNFiltered::PNFiltered(const char *pnSequence, int length)
{
// 
// initialize objects
//
  _pnGenerator          = NULL;
  _filteredBits         = NULL;
//
// Generate bits based on input PN sequence:
//
  generateData(pnSequence, length);
  return;
}

// ############################# Class Destructor ###############################
// PNFiltered -- Destructor for the PNFiltered class
// Input:               None
// Output:              None
//
// NOTES:
// ############################# Class Destructor ###############################
PNFiltered::~PNFiltered()
{
//
// Delete allocated space:
//
  delete _pnGenerator;
  delete [] _filteredBits;
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the shift register stages of the PN generator.
// Input:               None
// Output:              None
//
// Notes:
// 1. The shift register is initialized to a[0] = 1, a[i] = 0;
// ############################ Public Function ###################################
void PNFiltered::initGenerator()
{
  if(_pnGenerator != NULL)
  {
    _pnGenerator->initGenerator();
  }
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
void PNFiltered::setPNPolynomial(int poly)
{
  if(_pnGenerator != NULL)
  {
    _pnGenerator->setPNPolynomial(poly);
    generateData();
  }
  return;
}

// ############################ Public Function ###################################
// getEarlyLateOutput - Gets new data from the PN generator for the early/late PNs.
//
// Input:       index:  Index into _filteredBits array
//              offset: Fractional offset [0, 2*PI]
//
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
float PNFiltered::getEarlyLateOutput(int index, float offset)
{
  index         = MAX(0, index);
  index         %= _numberBits;
  switch (_filteredBits[index])
  {
    case PLUS1_BIT:
    default:
      return 1.;
    case MINUS1_BIT:
      return -1.;
    case MINUS_SLOPE_BIT:
      return (1. - offset/PI);
    case PLUS_SLOPE_BIT:
      return (offset/PI - 1.);
  }
}

// ############################ Public Function ###################################
// getOnTimeOutput - Gets new data from the PN generator.
//
// Input:       index:  Index into _filteredBits array
//
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
float PNFiltered::getOnTimeOutput(int index)
{
  index         = MAX(0, index);
  index         %= _numberBits;
  switch (_filteredBits[index])
  {
    case PLUS_SLOPE_BIT:
    case MINUS1_BIT:
    default:
      return -1.;
    case MINUS_SLOPE_BIT:
    case PLUS1_BIT:
      return 1.;
  }
}