/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer.    *
 * This is an abstract base class, see the subclasses, FloatBuffer, and         *
 * DoubleBuffer for actual implementations.                                     *
 *                                                                              *
 * File: /User/frank/C++/Misc/AbstractBuffer.h                                  *
 *                                                                              *
 *                                                                              *
 ********************************************************************************/
#include    "AbstractBuffer.h"

#define MAX_LENGTH  2048
#define MAX_SIZE    128
#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )

// ############################ Class Constructor ################################
// AbstractBuffer -- Constructor for the AbstractBuffer class
// Input:       length:         Number of stages in the buffer
//              size:           The width of the buffer
//
// Output:                          None
// ############################ Class Constructor ################################
AbstractBuffer::AbstractBuffer(int length, int size)
{
//
// Init local variables and buffer:
//
  _bufferLength = MIN(MAX_LENGTH, length);
  _bufferLength = MAX(1, _bufferLength);
  _bufferSize   = MIN(MAX_SIZE, size);
  _readPointer  = _writePointer         = 0;

  return;
}

// ############################ Class Destructor #################################
// ~AbstractBuffer -- Destructor for the ISU class
// Input:     None
// Output:    None
// Delete all the allocated objects:
// ############################ Class Destructor #################################
AbstractBuffer::~AbstractBuffer()
{
  return;
}

// ############################ Public Function ###################################
// setReadPointerToOldest -- Sets the read pointer to point to oldest data
//
// Input:                   None
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
void AbstractBuffer::setReadPointerToOldest()
{
  _readPointer  = _writePointer;
  return;
}


// ############################ Public Function ###################################
// setReadPointerToWritePointer -- Sets the read pointer to write pointer in order
//                                 to get oldest data first.
//
// Input:                   None
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
void AbstractBuffer::setReadPointerToWritePointer()
{
  _readPointer  = _writePointer;
  return;
}


// ############################ Public Function ###################################
// resetBuffer -- Sets the buffer to all zeros, should be overridden in subclasses.
//
// Input:                   None
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
void AbstractBuffer::resetBuffer()
{
  _readPointer  = _writePointer = 0;
  return;
}
