#ifndef _DOUBLEBUFFER_H
#define _DOUBLEBUFFER_H  1
/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Misc/DoubleBuffer.h                                    *
 *                                                                              *
 *                                                                              *
 ********************************************************************************/
#include "AbstractBuffer.h"

class DoubleBuffer: public AbstractBuffer
{
private:
  double     **_buffer;                      // The bufferSize x bufferLength array
//
// private functions:
//

//
// public functions:
//
public:
  DoubleBuffer(int length, int size=1);              // Constructor
  ~DoubleBuffer();                                   // Destructor

/********************************
 * Initializing buffer          *
 ********************************/
  void  resetBuffer();
  
/********************************
 * Writing into buffer.         *
 ********************************/
  void writeToBuffer(double *data);                  // Write to buffer, increment write pointer

/********************************
 * Reading from buffer.         *
 ********************************/
  const double *readFromBuffer();                    // Read from buffer, increment read pointer
//
// This function zeros out the buffer:
//

};

#endif
