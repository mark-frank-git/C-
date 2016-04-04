/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Misc/DoubleBuffer.cc                                   *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 08/15/96 - Started                                                       *
 *                                                                              *
 ********************************************************************************/
#include    "DoubleBuffer.h"

// ############################ Class Constructor ################################
// DoubleBuffer -- Constructor for the DoubleBuffer class
// Input:       length:         Number of stages in the buffer
//              size:           The width of the buffer
//
// Output:                          None
// ############################ Class Constructor ################################
DoubleBuffer::DoubleBuffer(int length, int size)
            :AbstractBuffer(length, size)
{
  int i;
//
// Init local variables and buffer:
//
  _buffer       = new double *[_bufferLength];
  for(i=0; i<_bufferLength; i++)
    _buffer[i]  = new double[_bufferSize];

  return;
}

// ############################ Class Destructor #################################
// ~DoubleBuffer -- Destructor for the DoubleBuffer class
// Input:     None
// Output:    None
// Delete all the allocated objects:
// ############################ Class Destructor #################################
DoubleBuffer::~DoubleBuffer()
{
  int i;
  for(i=0; i<_bufferLength; i++)
    delete [] _buffer[i];
  delete _buffer;
  return;
}

// ############################ Public Function ###################################
// resetBuffer -- Sets the buffer to all zeros.
//
// Input:                   None
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
void DoubleBuffer::resetBuffer()
{
  int   i, j;

  _readPointer  = _writePointer = 0;
  for(i=0; i<_bufferLength; i++)
  {
    for(j=0; j<_bufferSize; j++)
      _buffer[i][j]     = 0.;
  }
  return;
}

// ############################ Public Function ###################################
// writeToBuffer -- Writes a new set of values into the buffer
//
// Input:       data:       The array which needs to be of size, _bufferSize
// Output:                  None
//
// NOTES:
// 1. This function writes the data array into one tap of the buffer
// ############################ Public Function ###################################
void DoubleBuffer::writeToBuffer(double *data)
{
  int i;
  for(i=0; i<_bufferSize; i++)
    _buffer[_writePointer][i]   = data[i];
  _writePointer++;
  if(_writePointer > _bufferLength-1)
    _writePointer               = 0;
  return;
}

// ############################ Public Function ###################################
// readFromBuffer -- Returns a set of values from the buffer
//
// Input:       data:       The array which needs to be of size, _bufferSize
// Output:                  None
//
// NOTES:
// 1. This function post increments the read pointer, see also, setReadPointerToOldest()
// ############################ Public Function ###################################
const double *DoubleBuffer::readFromBuffer()
{
  const double *output;

  output                = _buffer[_readPointer];
  _readPointer++;
  if(_readPointer > _bufferLength-1)
    _readPointer        = 0;
  return output;
}
