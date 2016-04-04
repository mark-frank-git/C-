/********************************************************************************
 *                                                                              *
 * This class generates GSM Sync data according to GSM 05.02.                   *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMSyncGenerator.cc                         *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 04/19/2012 - Started.                                                    *
 ********************************************************************************/
#include <stdio.h>
#include "GSMSyncGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define GSM_SYNC_TRAINING_BITS          64
#define MAX_TRAINING_SEQ                3
int sync_training_sequences[MAX_TRAINING_SEQ*GSM_SYNC_TRAINING_BITS] =
{
  1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
  1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1,
  1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0,
  0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0,
  1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1,
  1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0
};


// ############################# Class Constructor #################################
// GSMSyncGenerator -- Constructor for the GSMSyncGenerator class
// Input:       tsc:    New Training Sequence Code number
// Output:              None
// ############################# Class Constructor #################################
GSMSyncGenerator::GSMSyncGenerator(int tsc)
{
// 
// Initialize instance variables:
//
  setTSCNumber(tsc);
  initGenerator();

  return;
}

// ############################# Class Destructor ###############################
// GSMSyncGenerator -- Destructor for the GSMSyncGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
GSMSyncGenerator::~GSMSyncGenerator()
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
void GSMSyncGenerator::initGenerator()
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
void GSMSyncGenerator:: setTSCNumber(int tsc)
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
int GSMSyncGenerator::newData()
{
  int index     = _tscNumber*GSM_SYNC_TRAINING_BITS + _bitCount;
  int data      = sync_training_sequences [index];
  if(_bitCount++ == GSM_SYNC_TRAINING_BITS)
    _bitCount   = 0;
  return data;
}
