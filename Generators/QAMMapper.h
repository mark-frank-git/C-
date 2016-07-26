#ifndef _QAM_MAPPER_H
#define _QAM_MAPPER_H 1
/********************************************************************************
 * This subclass of Object implements a QAM constellation mapper.               *
 * File: QAMMapper.h                                                            *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 03/28/03 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif


//
// QAM types:
//
#define QAM_DVBC                0
#define QAM_J83                 1

#define MAX_BITS                8

class QAMMapper
{
protected:
  int   _qamType;                       // DVB-C, J.83, etc.
  int   _qamBits;                       // # of bits per constellation point
  int   _qamLevels;                     // # of levels in I and Q dimensions, not valid for odd const
  int   _qamPoints;                     // # of constellation points

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
  QAMMapper(int type=QAM_DVBC, int bits=6);
  ~QAMMapper();

/********************************
 * Setting parameters:          *
 ********************************/
  void  setQAMType(int type);
  void  setQAMBits(int bits);

/********************************
 *  Getting parameters:         *
 ********************************/
  int   qamType()       {return _qamType;}
  int   qamBits()       {return _qamBits;}
  int   qamLevels()     {return _qamLevels;}

/********************************
 * Getting QAM mappings:        *
 * The index referred to below  *
 * starts at 0 in the lower     *
 * left corner of the constel-  *
 * lation and increases by ones *
 * along the I direction.       *
 ********************************/
  void  iAndQAtIndex(int index, float *iData, float *qData);    // return the I and Q constellation values
  void  iAndQForBitMap(int bitMap, float *iData, float *qData); // return the I and Q for bit map
  int   bitMappingAtIndex(int index);                           // returns the bit mapping for the given point
  int   indexAtIAndQ(float iData, float qData);                 // returns the index for the I and Q data
  int   indexForBitMap(int bitMap);                             // This should only be used for table gen.


};
#endif


