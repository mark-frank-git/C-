/********************************************************************************
 *                                                                              *
 * This class generates EDGE data according to GSM 05.02.                       *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/EDGEGenerator.cc                            *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "EDGEGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define EDGE_NB_TRAINING_BITS           26
#define MAX_TRAINING_SEQ                8
int edge_training_sequences[MAX_TRAINING_SEQ*EDGE_NB_TRAINING_BITS*EIGHT_PSK_BITS_PER_SYMBOL] =
{
  1,1,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1,
  1,1,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1,

  1,1,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1,
  0,0,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1,

  1,1,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1,
  1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1,

  1,1,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1,
  0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1,

  1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 1,1,1,
  0,0,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1,

  1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 1,1,1,
  1,1,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 1,1,1,

  0,0,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1,
  1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1, 0,0,1, 1,1,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1,

  0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 1,1,1, 1,1,1, 0,0,1, 1,1,1,
  1,1,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,1,1, 1,1,1,
};


// ############################# Class Constructor #################################
// EDGEGenerator -- Constructor for the EDGEGenerator class
// Input:       tsc:    New Training Sequence Code number
// Output:              None
// ############################# Class Constructor #################################
EDGEGenerator::EDGEGenerator(int tsc)
{
// 
// Initialize instance variables:
//
  setTSCNumber(tsc);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// EDGEGenerator -- Destructor for the EDGEGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
EDGEGenerator::~EDGEGenerator()
{
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the EDGE generator.
// Input:               None
// Output:              None
//
// Notes:
// 1. None
// ############################ Public Function ###################################
void EDGEGenerator::initGenerator()
{
  _bitCount     = 0;
  return;
}

#define OCTAL_BITS  3
// ############################ Public Function ###################################
// setTSCNumber - Sets a new TSC number.
// Input:       tsc:    New training sequence code
// Output:              None
//
// Notes:
// 1. None;
// ############################ Public Function ###################################
void EDGEGenerator:: setTSCNumber(int tsc)
{
//
  _tscNumber    = MAX(0, tsc);
  _tscNumber    = MIN(tsc, (MAX_TRAINING_SEQ-1));

  return;
}

// ############################ Public Function ###################################
// newData - Gets new data from the EDGE generator.
// Input:               None
// Output:              3 new data bits in (0,1)
//
// Notes:
// 1. None
//
// ############################ Public Function ###################################
const int *EDGEGenerator::newData()
{
  int index     = _tscNumber*EDGE_NB_TRAINING_BITS* EIGHT_PSK_BITS_PER_SYMBOL + _bitCount*EIGHT_PSK_BITS_PER_SYMBOL;
  int k         = EIGHT_PSK_BITS_PER_SYMBOL - 1;
  for(int i=0; i< EIGHT_PSK_BITS_PER_SYMBOL; i++)
  {
    _edgeBits[k]        = edge_training_sequences[index++];
    k--;
  }
  if(_bitCount++ == EDGE_NB_TRAINING_BITS)
    _bitCount   = 0;
  return _edgeBits;
}
