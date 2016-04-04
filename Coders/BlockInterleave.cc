/********************************************************************************
 * This class implements a block interleaver.  The interleaver writes in by     *
 * columns and reads by out by rows.                                            *
 * File: BlockInterleave.h                                                      *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "BlockInterleave.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#ifndef YES
#define YES 1
#define NO  0
#endif

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the BlockInterleave class.
// Input:               rows:           number of rows in interleaver
//                      columns:        number of columns interleaver
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
BlockInterleave::BlockInterleave(int rows, int columns)
             :BlockLeaveDeleave(rows, columns)
{
  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the BlockInterleave class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
BlockInterleave::~BlockInterleave()
{
  return;
}

// ############################# Public Function ###############################
// writeData -- Puts a new floating point data into the interleaver.
//
// Input:       data:           New data point to write into interleaver.
//          
// Output:                      None
//
// Notes: Interleaver writes in by columns.
// ############################# Public Function ###############################
void BlockInterleave:: writeData(float data)
{
  if(_blockWritePtr == 0)
  {
    _writePtr           = _rowPtr*_blockColumns + _columnPtr;
    _block0[_writePtr]  = data;
    _rowPtr++;
    if(_rowPtr >= _blockRows)
    {
      _rowPtr           = 0;
      _columnPtr++;
      if(_columnPtr >= _blockColumns)
      {
        _columnPtr      = 0;
        _blockWritePtr++;
      }
    }
  }
  else
  {
    _writePtr           = _rowPtr*_blockColumns + _columnPtr;
    _block1[_writePtr]  = data;
    _rowPtr++;
    if(_rowPtr >= _blockRows)
    {
      _rowPtr           = 0;
      _columnPtr++;
      if(_columnPtr >= _blockColumns)
      {
        _columnPtr      = 0;
        _blockWritePtr  = 0;
      }
    }
  }
  return;
}

// ############################# Public Function ###############################
// readData -- Gets a new floating point data from the interleaver.
//
// Input:                       None
//          
// Output:                      New floating point data
//
// Notes:
// 1) Deleaver reads out by columns
// 2) call readyToRead() before calling this function.
// ############################# Public Function ###############################
float BlockInterleave:: readData()
{
  float output_data;
  if(_blockReadPtr == 0)
  {
    output_data         = _block0[_readPtr];
    _readPtr++;
    if(_readPtr == _blockSize)
    {
      _readPtr          = 0;
      _blockReadPtr++;
    }
  }
  else
  {
    output_data         = _block1[_readPtr];
    _readPtr++;
    if(_readPtr == _blockSize)
    {
      _readPtr          = 0;
      _blockReadPtr--;
    }
  }
  return output_data;
}

// ############################# Public Function ###############################
// readyToRead -- Determines whether the interleaver is ready to output data.
//
// Input:                       None
//          
// Output:                      YES if data ready to be read
//
// NOTES:
// 1. readyToRead is TRUE when the writePtr is such that when we start reading
//    out we won't run out of new data.
// ############################# Public Function ###############################
LOGICAL BlockInterleave:: readyToRead()
{
  int                   number_written;
  if(_readyToRead)
    return YES;
  number_written        = _columnPtr*_blockRows + _rowPtr;
  if(number_written > (_blockSize - _blockRows - _blockColumns + 1) )
    _readyToRead        = YES;

  return _readyToRead;
}
