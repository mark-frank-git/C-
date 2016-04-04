#ifndef _GSMGENERATOR_H
#define _GSMGENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates GSM data according to GSM 05.02.                        *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMGenerator.h                              *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 04/19/2012 - Started.                                                    *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif


class GSMGenerator
{
protected:
  int           _tscNumber;                     // Training sequence code number (0 to 7)
  int           _bitCount;                      // Bit number count
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  GSMGenerator(int tsc);                        // Constructor with input tsc number
  virtual ~GSMGenerator();

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
 * new GSM data.                *
 ********************************/
  int   newData();                              // Gets the next data bit from gen


};

#endif
