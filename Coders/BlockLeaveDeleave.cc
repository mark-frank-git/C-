/********************************************************************************
 * This class is the super class for the block interleaver and deleaver.        *
 * See the subclasses, BlockDeleave and BlockInterleave for real implementations*
 * File: BlockLeaveDeleave.cc                                                   *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "BlockLeaveDeleave.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#ifndef YES
#define YES 1
#define NO  0
#endif


// ############################# Public Function ###############################
// allocateBlocks -- Allocate the blocks to be used in leaver/deleaver.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void BlockLeaveDeleave::allocateBlocks()
{

  _blockSize    = _blockRows * _blockColumns;
  if(_oldSize < _blockSize)
  {
    _oldSize    = _blockSize;
    delete [] _block0;
    delete [] _block1;
    _block0     = new float[_blockSize];
    _block1     = new float[_blockSize];
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the BlockLeaveDeleave class.
// Input:               rows:           number of rows in inter/deleaver
//                      columns:        number of columns inter/deleaver
//          
// Output:                          None
//
// Notes:
// ############################# Class Constructor ###############################
BlockLeaveDeleave::BlockLeaveDeleave(int rows, int columns)
{
  _block0               = _block1       = NULL;
  _oldSize              = 0;
  _blockRows            = _blockColumns = 0;
  setRows(rows);
  setColumns(columns);
  reset();
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the BlockLeaveDeleave class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
BlockLeaveDeleave::~BlockLeaveDeleave()
{
  delete [] _block0;
  delete [] _block1;
  return;
}

// ############################# Public Function ###############################
// reset -- Reset the leaver/deleaver.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void BlockLeaveDeleave::reset()
{
  _readPtr              = _writePtr     = 0;
  _blockWritePtr        = _blockReadPtr = 0;
  _rowPtr               = _columnPtr    = 0;
  _readyToRead          = NO;
  return;
}

// ############################# Public Function ###############################
// setRows -- Sets the number of rows.
//
// Input:       rows:           Number of rows in leaver/deleaver
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void BlockLeaveDeleave::setRows (int rows)
{
  _blockRows            = MIN(rows, MAX_ROWS);
  _blockRows            = MAX(_blockRows, MIN_ROWS);

  allocateBlocks();
  return;
}

// ############################# Public Function ###############################
// setColumns -- Sets the number of columns.
//
// Input:       columns:        Number of columns in leaver/deleaver
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void BlockLeaveDeleave:: setColumns (int cols)
{
  _blockColumns         = MIN(cols, MAX_COLUMNS);
  _blockColumns         = MAX(_blockColumns, MIN_COLUMNS);

  allocateBlocks();
  return;
}
