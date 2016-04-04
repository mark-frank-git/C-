#ifndef _GSMDUMMYGENERATOR_H
#define _GSMDUMMYGENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates GSM Sync Burst data according to GSM 05.02.             *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMDummyGenerator.h                         *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 06/06/2012 - Started.                                                    *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif


class GSMDummyGenerator
{
protected:
  int           _bitCount;                      // Bit number count
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  GSMDummyGenerator();                          // Constructor with input tsc number
  virtual ~GSMDummyGenerator();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/

/*******************************
 * These methods get parameters*
 *******************************/

/********************************
 * These methods generate       *
 * new GSM data.                *
 ********************************/
  int   newData();                              // Gets the next data bit from gen


};

#endif
