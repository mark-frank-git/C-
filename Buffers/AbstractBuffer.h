#ifndef _ABSTRACTBUFFER_H
#define _ABSTRACTBUFFER_H  1
/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer.    *
 * This is an abstract base class, see the subclasses, ComplexBuffer,           *
 * DoubleBuffer, FloatBuffer, and IntBuffer for actual implementations.         *
 *                                                                              *
 * File: /User/frank/C++/Misc/AbstractBuffer.h                                  *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 06/09/98 - Started                                                       *
 *                                                                              *
 ********************************************************************************/
class AbstractBuffer
{
protected:
  int       _bufferLength;                  // Number of stages in the buffer
  int       _bufferSize;                    // The width of the buffer
  int       _writePointer, _readPointer;    // Ring buffer pointers
//
// private functions:
//

//
// public functions:
//
public:
  AbstractBuffer(int length, int size=1);               // Constructor
  virtual ~AbstractBuffer();                            // Destructor
//
// These functions manipulate the pointers:
//
  void  setReadPointerToOldest();                       // Set the read buffer pointer to oldest data
  void  setReadPointerToWritePointer();
//
// These functions return parameters:
//
  int bufferLength()            {return _bufferLength;}
  int bufferSize()              {return _bufferSize;}
//
// This function zeros out the buffer:
//
  virtual void  resetBuffer();

};

#endif
