#ifndef _GSMSYNCGENERATOR_H
#define _GSMSYNCGENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates GSM Sync Burst data according to GSM 05.02.             *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMSyncGenerator.h                          *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 06/06/2012 - Started.                                                    *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif


class GSMSyncGenerator
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
  GSMSyncGenerator(int tsc);                    // Constructor with input tsc number
  virtual ~GSMSyncGenerator();

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
