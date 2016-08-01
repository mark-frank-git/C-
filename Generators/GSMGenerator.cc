/********************************************************************************
 *                                                                              *
 * This class generates GSM data according to GSM 05.02.                        *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMGenerator.cc                             *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>
#include "GSMGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define GSM_NB_TRAINING_BITS            26
#define MAX_TRAINING_SEQ                8
int training_sequences[MAX_TRAINING_SEQ*GSM_NB_TRAINING_BITS] =
{
  0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1,
  0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1,
  0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0,
  0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0,
  0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1,
  0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0,
  1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1,
  1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0
};


// ############################# Class Constructor #################################
// GSMGenerator -- Constructor for the GSMGenerator class
// Input:       tsc:    New Training Sequence Code number
// Output:              None
// ############################# Class Constructor #################################
GSMGenerator::GSMGenerator(int tsc)
{
// 
// Initialize instance variables:
//
  setTSCNumber(tsc);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// GSMGenerator -- Destructor for the GSMGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
GSMGenerator::~GSMGenerator()
{
  
  return;
}

// ############################ Public Function ###################################
// initGenerator - Initializes the GSM generator.
// Input:               None
// Output:              None
//
// Notes:
// 1. None
// ############################ Public Function ###################################
void GSMGenerator::initGenerator()
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
void GSMGenerator:: setTSCNumber(int tsc)
{
//
  _tscNumber    = MAX(0, tsc);
  _tscNumber    = MIN(tsc, (MAX_TRAINING_SEQ-1));

  return;
}

// ############################ Public Function ###################################
// newData - Gets new data from the GSM generator.
// Input:               None
// Output:              New data bit (0,1)
//
// Notes:
// 1. None
//
// ############################ Public Function ###################################
int GSMGenerator::newData()
{
  int index     = _tscNumber*GSM_NB_TRAINING_BITS + _bitCount;
  int data      = training_sequences[index];
  if(_bitCount++ == GSM_NB_TRAINING_BITS)
    _bitCount   = 0;
  return data;
}
