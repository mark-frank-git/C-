#ifndef _BLOCK_LEAVE_DELEAVE_H
#define _BLOCK_LEAVE_DELEAVE_H 1
/********************************************************************************
 * This class is the super class for the block interleaver and deleaver.        *
 * See the subclasses, BlockDeleave and BlockInterleave for real implementations*
 * File: BlockLeaveDeleave.h                                                    *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif

#define MIN_ROWS        1
#define MIN_COLUMNS     1
#define MAX_ROWS        10000
#define MAX_COLUMNS     10000

class BlockLeaveDeleave
{
protected:
  int           _blockRows;             // # of rows in inter/deleaver
  int           _blockColumns;          // # of columns in inter/deleaver
  int           _blockSize;             // rows * columns
  int           _writePtr;              // For writing into inter/deleaver
  int           _readPtr;               // For reading out of inter/deleaver
  int           _rowPtr;                // Row pointer for reading/writing
  int           _columnPtr;             // Column pointer for reading/writing
  int           _blockWritePtr;         // Block 0 or block 1 pointer for writes
  int           _blockReadPtr;          // Block 0 or block 1 pointer for reads
  int           _oldSize;               // Old block sizes for news and deletes

  float         *_block0;               // One of the blocks for saving/reading out
  float         *_block1;               // The alternative block for saving/reading out

  LOGICAL       _readyToRead;           // YES = ready to read data out of inter/deleaver

//
// Private functions:
//
  void          allocateBlocks();       // Allocate blocks for inter/deleavers

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  BlockLeaveDeleave (int rows=2, int columns=2);
  virtual ~ BlockLeaveDeleave ();

/********************************
 * Initializing inter/deleaver: *
 ********************************/
  virtual void  reset();
  
/********************************
 * Setting inter/deleave        *
 * parameters:                  *
 ********************************/
  void          setRows(int rows);
  void          setColumns(int cols);

/********************************
 * Getting inter/deleave        *
 * parameters:                  *
 ********************************/
  int           rows()                  {return _blockRows;}
  int           columns()               {return _blockColumns;}

};
#endif


