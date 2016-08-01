/********************************************************************************
 *                                                                              *
 * This class implements the functionality of a variable length ring buffer     *
 *                                                                              *
 * File: /User/frank/C++/Buffers/ComplexBuffer.cc                               *
 *                                                                              *
 *                                                                              *
 ********************************************************************************/
#include        "ComplexBuffer.h"
#include        <GNU/Complex.h>

// ############################ Class Constructor ################################
// ComplexBuffer -- Constructor for the ComplexBuffer class
// Input:       length:         Number of stages in the buffer
//              size:           The width of the buffer
//
// Output:                          None
// ############################ Class Constructor ################################
ComplexBuffer::ComplexBuffer(int length, int size)
            :AbstractBuffer(length, size)
{
  int i;
//
// Init local variables and buffer:
//
  _buffer       = new Complex *[_bufferLength];
  for(i=0; i<_bufferLength; i++)
    _buffer[i]  = new Complex[_bufferSize];

  return;
}

// ############################ Class Destructor #################################
// ~ComplexBuffer -- Destructor for the ComplexBuffer class
// Input:     None
// Output:    None
// Delete all the allocated objects:
// ############################ Class Destructor #################################
ComplexBuffer::~ComplexBuffer()
{
  int i;
  for(i=0; i<_bufferLength; i++)
    delete [] _buffer[i];
  delete _buffer;
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
void ComplexBuffer::writeToBuffer(Complex *data)
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
const Complex *ComplexBuffer::readFromBuffer()
{
  const Complex *output;

  output                = _buffer[_readPointer];
  _readPointer++;
  if(_readPointer > _bufferLength-1)
    _readPointer        = 0;
  return output;
}


// ############################ Public Function ###################################
// readOldSampleAt -- Reads a previously written sample at an index
//
// Input:       index:  Past sample location, if index == 0, return last written sample
//                      if index == 1, return previous sample, etc.
// Output:              Past sample
//
// NOTES:
// 1. This function does not affect the read or write pointers
// ############################ Public Function ###################################
const Complex *ComplexBuffer::readOldSampleAt(int index)
{
  int   read_pointer;
  const Complex *output;

  read_pointer          = _writePointer - index - 1;
  while(read_pointer < 0)
    read_pointer        += _bufferLength;
  output                = _buffer[read_pointer];

  return output;
}

// ############################ Public Function ###################################
// resetBuffer -- Sets the buffer to all zeros.
//
// Input:                   None
// Output:                  None
//
// NOTES:
// ############################ Public Function ###################################
void ComplexBuffer::resetBuffer()
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
