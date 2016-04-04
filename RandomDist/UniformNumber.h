#ifndef _UNIFORM_NUMBER_H
#define _UNIFORM_NUMBER_H 1
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

#define SEED_OFFSET     1167                    // Default seed

class   MTRand;

class UniformNumber
{
private:
  MTRand        *_randomNumber;                 // Random number generator
  long          _seed;                          // Random number seed
//
// Private functions:
//

public:
/*********************************
 * Constructors/destructors:     *
 *********************************/
  UniformNumber(long seed=SEED_OFFSET);         // Class Constructor
  ~UniformNumber();                             // Class Destructor

/*******************************
 * These methods initialize:    *
 *******************************/
  void  reset();

/*******************************
 * These methods set parameters:*
 *******************************/
  void setSeed(long seed);

/********************************
 * These methods get outputs    *
 * from the generator.          *
 ********************************/
  double operator ()();                         // Returns random number in [0, 1];
  
};

#endif
