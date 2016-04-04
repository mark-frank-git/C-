#ifndef _COMPLEXBUFFER_H
#define _COMPLEXBUFFER_H  1
/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Buffers/ComplexBuffer.h                                *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 10/24/01 - Started.                                                      *
 *                                                                              *
 ********************************************************************************/
#include "AbstractBuffer.h"
class   Complex;

class ComplexBuffer: public AbstractBuffer
{
private:
  Complex       **_buffer;                              // The bufferSize x bufferLength array
//
// private functions:
//

//
// public functions:
//
public:
  ComplexBuffer(int length, int size=1);                // Constructor
  ~ComplexBuffer();                                     // Destructor
//
// These functions write into the buffer:
//
  void writeToBuffer(Complex *data);                    // Write to buffer, increment write pointer
//
// These functions get from the buffer:
//
  const Complex *readFromBuffer();                      // Read from buffer, increment read pointer
  const Complex *readOldSampleAt(int index);
//
// This function zeros out the buffer:
//
  void  resetBuffer();

};

#endif
