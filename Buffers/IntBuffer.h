#ifndef _INT_BUFFER_H
#define _INT_BUFFER_H  1
/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Misc/IntBuffer.h                                       *
 *                                                                              *
 *                                                                              *
 ********************************************************************************/
#include "AbstractBuffer.h"

class IntBuffer: public AbstractBuffer
{
private:
  int           **_buffer;                      // The bufferSize x bufferLength array
//
// private functions:
//

//
// public functions:
//
public:
  IntBuffer(int length, int size=1);              // Constructor
  ~IntBuffer();                                   // Destructor
//
// These functions write into the buffer:
//
  void writeToBuffer(int *data);                  // Write to buffer, increment write pointer
//
// These functions get from the buffer:
//
  const int *readFromBuffer();                    // Read from buffer, increment read pointer
//
// This function zeros out the buffer:
//
  void  resetBuffer();

};

#endif
