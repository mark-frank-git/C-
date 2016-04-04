#ifndef _EDGEGENERATOR_H
#define _EDGEGENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates EDGE data according to GSM 05.02.                       *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/EDGEGenerator.h                             *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 04/19/2012 - Started.                                                    *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif
#define EIGHT_PSK_BITS_PER_SYMBOL       3

class EDGEGenerator
{
protected:
  int           _tscNumber;                     // Training sequence code number (0 to 7)
  int           _bitCount;                      // Bit number count
  int           _edgeBits[EIGHT_PSK_BITS_PER_SYMBOL];   //Output bits
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  EDGEGenerator(int tsc);                       // Constructor with input tsc number
  virtual ~EDGEGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setTSCNumber(int tsc);                  // Sets a new TSC

/*******************************
 * These methods get parameters*
 *******************************/
  int   getTSCNumber () {return _tscNumber;}

/********************************
 * These methods generate       *
 * new EDGE data.               *
 ********************************/
  const int     *newData();                     // Gets the next set of data bits from gen


};

#endif
