/********************************************************************************
 *                                                                              *
 * This class generates the I and Q short code PN sequences for IS-95.  The     *
 * outputs from the shift register are mapped according to (0,1)->(1,-1).       *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/IS95PNGen.h                                 *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 03/16/01 - Started.                                                      *
 ********************************************************************************/
#include <stdio.h>
#include "IS95PNGen.h"
#include "PNGenerator.h"

#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define BIT_MAP(a)      ( ((a)==0 ) ? 1 : -1 )          // (0,1) -> (1, -1)

// ############################# Class Constructor #################################
// IS95PNGen -- Constructor for the IS95PNGen class
// Input:               None
// Output:              None
//
// Notes:
// 1. In the constructor, the data is generated and stored in the local arrays
// ############################# Class Constructor #################################
IS95PNGen::IS95PNGen()
{
  int   i;
  PNGenerator   *pn_generator;
// 
// Allocate the local storage:
//
  _iData        = new char[PN_LENGTH];
  _qData        = new char[PN_LENGTH];
//
// Generate the I data:
//
  pn_generator  = new PNGenerator(I_PN_POLY);
  pn_generator->initGenerator();
  for(i=0; i<(PN_LENGTH-1); i++)
    _iData[i]   = BIT_MAP((pn_generator->newData()));
  _iData[PN_LENGTH-1]   = BIT_MAP(0);
//
// Generate the Q data:
//
  pn_generator->setPNPolynomial(Q_PN_POLY);
  pn_generator->initGenerator();
  for(i=0; i<(PN_LENGTH-1); i++)
    _qData[i]   = BIT_MAP((pn_generator->newData()));
  _qData[PN_LENGTH-1]   = BIT_MAP(0);

  delete pn_generator;
  return;
}

// ############################# Class Destructor ###############################
// IS95PNGen -- Destructor for the IS95PNGen class
// Input:               None
// Output:              None
//
// NOTES:
// 1. Deletes are done in the subclasses
// ############################# Class Destructor ###############################
IS95PNGen::~IS95PNGen()
{
//
// Delete space:
//
  delete [] _iData;
  delete [] _qData;

  return;
}

// ############################ Public Function ###################################
// iData - Gets new data from the I PN generator.
//
// Input:       offset: Offset into PN array
//
// Output:              Array of PN data mapped from (0,1) -> (1,-1)
//
// Notes:
// ############################ Public Function ###################################
const char *IS95PNGen::iData(int offset)
{
  offset        = MIN(offset, (PN_LENGTH-1));
  return &_iData[offset];
}

// ############################ Public Function ###################################
// qData - Gets new data from the Q PN generator.
//
// Input:       offset: Offset into PN array
//
// Output:              Array of PN data mapped from (0,1) -> (1,-1)
//
// Notes:
// ############################ Public Function ###################################
const char *IS95PNGen::qData(int offset)
{
  offset        = MIN(offset, (PN_LENGTH-1));
  return &_qData[offset];
}
