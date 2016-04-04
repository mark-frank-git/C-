/********************************************************************************
 *                                                                              *
 * This class generates GSM Dummy burst data according to GSM 05.02.            *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMDummyGenerator.cc                        *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 04/19/2012 - Started.                                                    *
 ********************************************************************************/
#include <stdio.h>
#include "GSMDummyGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define GSM_DUMMY_BITS                  142
#define MAX_TRAINING_SEQ                3
int dummy_training_sequence[GSM_DUMMY_BITS] =
{
 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0,
0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1,
0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0,
1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0,
1, 0, 1, 0
};


// ############################# Class Constructor #################################
// GSMDummyGenerator -- Constructor for the GSMDummyGenerator class
// Input:               None
// Output:              None
// ############################# Class Constructor #################################
GSMDummyGenerator::GSMDummyGenerator()
{
// 
// Initialize instance variables:
//
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// GSMDummyGenerator -- Destructor for the GSMDummyGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
GSMDummyGenerator::~GSMDummyGenerator()
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
void GSMDummyGenerator::initGenerator()
{
  _bitCount     = 0;
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
int GSMDummyGenerator::newData()
{
  int index     = _bitCount;
  int data      = dummy_training_sequence[index];
  if(_bitCount++ == GSM_DUMMY_BITS)
    _bitCount   = 0;
  return data;
}
