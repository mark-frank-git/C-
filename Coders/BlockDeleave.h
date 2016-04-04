#ifndef _BLOCK_DELEAVE_H
#define _BLOCK_DELEAVE_H 1
/********************************************************************************
 * This class implements a block deleaver.  The deleaver writes in by rows,     *
 * and reads by out by columns.                                                 *
 * File: BlockDeleave.h                                                         *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/24/01 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "BlockLeaveDeleave.h"

#ifndef LOGICAL
#define LOGICAL char
#endif

class BlockDeleave: public BlockLeaveDeleave
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
  BlockDeleave (int rows=2, int columns=2);
  ~ BlockDeleave ();

/********************************
 * Initializing deleaver:       *
 ********************************/
  
/********************************
 * Setting deleaver             *
 * parameters:                  *
 ********************************/

/********************************
 * Getting deleaver             *
 * parameters:                  *
 ********************************/

/********************************
 * Writing to and reading       *
 * from deleaver.               *
 ********************************/
  void          writeData(float data);
  float         readData();
  LOGICAL       readyToRead();

};
#endif


