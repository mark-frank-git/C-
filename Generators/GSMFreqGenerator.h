#ifndef _GSMFREQGENERATOR_H
#define _GSMFREQGENERATOR_H 1
/********************************************************************************
 *                                                                              *
 * This class generates GSM Freq Burst data according to GSM 05.02.             *
 * Notes:                                                                       *
 * 1. None                                                                      *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/GSMFreqGenerator.h                          *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 06/06/2012 - Started.                                                    *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif


class GSMFreqGenerator
{
protected:
//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  GSMFreqGenerator();                           // Constructor with input tsc number
  virtual ~GSMFreqGenerator();

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
