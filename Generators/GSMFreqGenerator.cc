/********************************************************************************
 *                                                                              *
 * This class generates GSM Dummy burst data according to GSM 05.02.            *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMFreqGenerator.cc                         *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 04/19/2012 - Started.                                                    *
 ********************************************************************************/
#include <stdio.h>
#include "GSMFreqGenerator.h"


// ############################# Class Constructor #################################
// GSMFreqGenerator -- Constructor for the GSMFreqGenerator class
// Input:               None
// Output:              None
// ############################# Class Constructor #################################
GSMFreqGenerator::GSMFreqGenerator()
{
  return;
}

// ############################# Class Destructor ###############################
// GSMFreqGenerator -- Destructor for the GSMFreqGenerator class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
GSMFreqGenerator::~GSMFreqGenerator()
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
void GSMFreqGenerator::initGenerator()
{
  return;
}

// ############################ Public Function ###################################
// newData - Gets new data from the GSM generator.
// Input:               None
// Output:              New data bit, always zero for freq burst
//
// Notes:
// 1. None
//
// ############################ Public Function ###################################
int GSMFreqGenerator::newData()
{
  return 0;
}
