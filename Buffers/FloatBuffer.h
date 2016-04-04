#ifndef _FLOATBUFFER_H
#define _FLOATBUFFER_H  1
/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Misc/FloatBuffer.h                                     *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 08/15/96 - Started                                                       *
 *  2. 06/09/98 - Made subclass of AbstractBuffer.                              *
 *                                                                              *
 ********************************************************************************/
#include "AbstractBuffer.h"

class FloatBuffer: public AbstractBuffer
{
private:
  float     **_buffer;                      // The bufferSize x bufferLength array
//
// private functions:
//

//
// public functions:
//
public:
  FloatBuffer(int length, int size=1);              // Constructor
  ~FloatBuffer();                                   // Destructor
//
// These functions write into the buffer:
//
  void writeToBuffer(float *data);                  // Write to buffer, increment write pointer
//
// These functions get from the buffer:
//
  const float *readFromBuffer();                    // Read from buffer, increment read pointer
//
// This function zeros out the buffer:
//
  void  resetBuffer();

};

#endif
