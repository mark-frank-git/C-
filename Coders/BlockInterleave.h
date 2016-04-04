#ifndef _BLOCK_INTERLEAVE_H
#define _BLOCK_INTERLEAVE_H 1
/********************************************************************************
 * This class implements a block interleaver.  The interleaver writes in by     *
 * columns and reads by out by rows.                                            *
 * File: BlockInterleave.h                                                      *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "BlockLeaveDeleave.h"

#ifndef LOGICAL
#define LOGICAL char
#endif

class BlockInterleave: public BlockLeaveDeleave
{
protected:
//
// Private functions:
//

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  BlockInterleave (int rows=2, int columns=2);
  ~ BlockInterleave ();

/********************************
 * Initializing interleaver:    *
 ********************************/
  
/********************************
 * Setting interleaver          *
 * parameters:                  *
 ********************************/

/********************************
 * Getting interleaver          *
 * parameters:                  *
 ********************************/

/********************************
 * Writing to and reading       *
 * from interleaver.            *
 ********************************/
  void          writeData(float data);
  float         readData();
  LOGICAL       readyToRead();

};
#endif


