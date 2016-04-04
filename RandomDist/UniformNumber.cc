/********************************************************************************
 *                                                                              *
 * This subclass of object implements a class for generating random numbers     *
 * It uses the MT19937 generator in the GSL library.                            *
 *                                                                              *
 * File: /User/frank/Objc_Classes/UniformNumber/UniformNumber.h                 *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 08/22/00 - Started.                                                      *
 ********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "UniformNumber.h"
#include "mtrand.h"

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a) ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the UniformNumber class.
//
// Input:           seed:       random generator seed
//
// Output:                      None
// Notes:
// ############################# Class Constructor ###############################
UniformNumber::UniformNumber (long seed)
{
// 
// Initialize instance variables:
//
  _seed         = seed;
  _randomNumber = new MTRand(_seed);
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the UniformNumber class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
UniformNumber::~UniformNumber ()
{
  delete _randomNumber;
  return;
}

// ############################ Public Function ###################################
// reset - Initializes the noise generator to a defined state.
// Input:               None
// Output:              None
//
// Notes:
// ############################ Public Function ###################################
void UniformNumber::reset()
{
  _randomNumber->seed(_seed);                   // Set old seed to initialize
  return;
}

// ############################ Public Function ###################################
// setSeed - Sets a new seed for generator.
// Input:       seed:           new seed
//
// Output:                      None
//
// Notes:
// 1. We need to make sure that seed is odd, so we multiply by 2 and add 1.  This
//    ensures if seed is set in a loop for(i=0; i<n) that each seed will be unique
// ############################ Public Function ###################################
void UniformNumber::setSeed(long seed)
{
  _seed         = seed+seed+1;
  _randomNumber->seed(_seed);                   // Set old seed to initialize
  return;
}

// ############################ Public Function ###################################
// () - Gets a new number in [0,1].
// Input:                       None
//
// Output:                      New random number in [0,1]
//
// Notes:
// ############################ Public Function ###################################
double UniformNumber::operator ()()
{
  double        output;
  output        = (*_randomNumber)();
  return (output/RANDOM_NUMBER_MAX_DBL);
}
