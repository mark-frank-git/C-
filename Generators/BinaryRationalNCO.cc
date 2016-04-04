/********************************************************************************
 * This class implements a binary-rational NCO.                                 *
 * It consists of two accumulators with the second accumulator taking           *
 * the carry out of the first accumulator.                                      *
 *                                                                              *
 * File: /User/frank/C++/Generators/BinaryRationalNCO.h                         *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 01/18/01 - Started.                                                      *
 ********************************************************************************/
#include <stdio.h>
#include "BinaryRationalNCO.h"
#include <C_Libraries/constants.h>

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )


// ############################# Class Constructor #################################
// BinaryRationalNCO -- Constructor for the BinaryRationalNCO class
// Input:       phases:         number of phases in the output
//              r0:             input accumulator increment
//              r1:             input accumulator alternate increment
//              r2:             output accumulator increment
// Output:                      None
// ############################# Class Constructor #################################
BinaryRationalNCO::BinaryRationalNCO(int phases, int r0, int r1, int r2)
{
//
// Initialize instance variables:
//
  setNumberPhases(phases);
  setR0(r0);
  setR1(r1);
  setR2(r2);
  initNCO();

  return;
}

// ############################# Class Destructor ###############################
// BinaryRationalNCO -- Destructor for the BinaryRationalNCO class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
BinaryRationalNCO::~BinaryRationalNCO()
{  
  return;
}

// ############################ Public Function ###################################
// initNCO - Initializes the clock generator.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void BinaryRationalNCO::initNCO()
{
  _inputAccumulator     = _outputAccumulator    = 0;
  _inputCarryOut        = NO;
  return;
}

// ############################ Public Function ###################################
// setNumberPhases - Sets a new number of output phases, this is also equal to the
//                   output accumulator modulus.
// Input:       phases:         Number of phases
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void BinaryRationalNCO::setNumberPhases(int phases)
{
  _numberPhases = MAX(1, phases);
  return;
}

// ############################ Public Function ###################################
// setR0 - Sets a new increment for the input accumulator.
// Input:       r0:             Default increment for input accumulator
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void BinaryRationalNCO::setR0(int r0)
{
  _r0                   = MAX(1, r0);
  _inputOverflow        = _r0 + _r1;
  return;
}

// ############################ Public Function ###################################
// setR1 - Sets a new increment for the input accumulator.
// Input:       r1:             Alternate increment for input accumulator
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void BinaryRationalNCO::setR1(int r1)
{
  _r1                   = MAX(1, r1);
  _inputOverflow        = _r0 + _r1;
  return;
}

// ############################ Public Function ###################################
// setR2 - Sets a new increment for the output accumulator.
// Input:       r2:             Increment for output accumulator
//
// Output:                      None
//
// Notes:
// ############################ Public Function ###################################
void BinaryRationalNCO::setR2(int r2)
{
  _r2                   = MAX(1, r2);
  return;
}

// ############################ Public Function ###################################
// newPhase - Gets a new phase output from the generator.
// Input:               None
//
// Output:              New phase: [0, _numberPhases-1]
//
// Notes:
// 1. This routine needs to be called once per sample time.
// ############################ Public Function ###################################
int BinaryRationalNCO::newPhase()
{
  int   output_increment;
//
// Update the input accumulator:
//
  if(_inputCarryOut)
    _inputAccumulator   -= _r1;
  else
    _inputAccumulator   += _r0;
  if(_inputAccumulator >= _inputOverflow)
    _inputCarryOut      = YES;
  else
    _inputCarryOut      = NO;
//
// Update the output accumulator:
//
  output_increment      = _r2;
  if(_inputCarryOut)
     output_increment++;
  _outputAccumulator    += output_increment;
  if(_outputAccumulator >= _numberPhases)
    _outputAccumulator  -= _numberPhases;
//
// Return the current value of the clock:
//
  return _outputAccumulator;
}

// ############################ Public Function ###################################
// newRadianPhase - Gets a new phase output from the generator.  The phase is given
//                  in radians
// Input:               None
//
// Output:              New phase: [0, TWOPI]
//
// Notes:
// 1. This routine needs to be called once per sample time.
// ############################ Public Function ###################################
double BinaryRationalNCO::newRadianPhase()
{
  int           integer_phase;
  double        output_phase;
//
// Get the phase from the routine above:
//
  integer_phase         = newPhase();
//
// Convert to radians:
//
  output_phase          = integer_phase*TWOPI;
  output_phase          /= _numberPhases;
//
// Return the current value of the clock:
//
  return output_phase;
}
