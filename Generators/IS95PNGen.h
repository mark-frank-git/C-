#ifndef _IS95PNGEN_H
#define _IS95PNGEN_H 1
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

#ifndef LOGICAL
#define LOGICAL char
#endif

#define PN_LENGTH       32768                   // Length of sequences
#define I_PN_POLY       121641                  // I channel polynomial, octal representation
#define Q_PN_POLY       116171                  // Q channel polynomial, octal representation

class   PNGenerator;                            // class prototype

class IS95PNGen
{
protected:
  char          *_iData;                        // Stored I PN sequence data
  char          *_qData;                        // Stored Q PN sequence data

//
// Private methods
//

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  IS95PNGen();                                  // Constructor for class
  ~IS95PNGen();

/*******************************
 * These methods get parameters*
 *******************************/
  int           pnLength()                      {return PN_LENGTH;}

/********************************
 * These methods generate       *
 * new PN data.                 *
 ********************************/
  const char    *iData(int offset=0);           // Returns an array of data from generator
  const char    *qData(int offset=0);


};

#endif
