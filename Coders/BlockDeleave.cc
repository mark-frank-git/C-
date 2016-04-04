/********************************************************************************
 * This class implements a block deleaver.  The deleaver writes in by rows,     *
 * and reads by out by columns.                                                 *
 * File: BlockDeleave.h                                                         *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "BlockDeleave.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#ifndef YES
#define YES 1
#define NO  0
#endif

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the BlockDeleave class.
// Input:               rows:           number of rows in deleaver
//                      columns:        number of columns deleaver
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
BlockDeleave::BlockDeleave(int rows, int columns)
             :BlockLeaveDeleave(rows, columns)
{
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the BlockDeleave class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
BlockDeleave::~BlockDeleave()
{
  return;
}

// ############################# Public Function ###############################
// writeData -- Puts a new floating point data into the deleaver.
//
// Input:       data:           New data point to write into deleaver.
//          
// Output:                      None
//
// Notes: Deleaver writes in by rows.
// ############################# Public Function ###############################
void BlockDeleave:: writeData(float data)
{
  if(_blockWritePtr == 0)
  {
    _block0[_writePtr]  = data;
    _writePtr++;
    if(_writePtr == _blockSize)
    {
      _writePtr         = 0;
      _blockWritePtr++;
    }
  }
  else
  {
    _block1[_writePtr]  = data;
    _writePtr++;
    if(_writePtr == _blockSize)
    {
      _writePtr         = 0;
      _blockWritePtr--;
    }
  }
  return;
}

// ############################# Public Function ###############################
// readData -- Gets a new floating point data from the deleaver.
//
// Input:                       None
//          
// Output:                      New floating point data
//
// Notes:
// 1) Deleaver reads out by columns
// 2) call readyToRead() before calling this function.
// ############################# Public Function ###############################
float BlockDeleave:: readData()
{
  float output_data;
  
  if(_blockReadPtr == 0)
  {
    _readPtr            = _rowPtr*_blockColumns + _columnPtr;
    output_data         = _block0[_readPtr];
    _rowPtr++;
    if(_rowPtr >= _blockRows)
    {
      _rowPtr           = 0;
      _columnPtr++;
      if(_columnPtr >= _blockColumns)
      {
        _columnPtr      = 0;
        _blockReadPtr++;
      }
    }
  }
  else
  {
    _readPtr            = _rowPtr*_blockColumns + _columnPtr;
    output_data         = _block1[_readPtr];
    _rowPtr++;
    if(_rowPtr >= _blockRows)
    {
      _rowPtr           = 0;
      _columnPtr++;
      if(_columnPtr >= _blockColumns)
      {
        _columnPtr      = 0;
        _blockReadPtr--;
      }
    }
  }
  return output_data;
}

// ############################# Public Function ###############################
// readyToRead -- Determines whether the deleaver is ready to output data.
//
// Input:                       None
//          
// Output:                      YES if data ready to be read
//
// NOTES:
// 1. readyToRead is TRUE when the writePtr is such that when we start reading
//    reading out we won't run out of new data.
// ############################# Public Function ###############################
LOGICAL BlockDeleave:: readyToRead()
{
  if(_readyToRead == YES)
    return YES;
  if(_writePtr > (_blockSize - _blockRows - _blockColumns + 1) )
    _readyToRead        = YES;
  return _readyToRead;
}
