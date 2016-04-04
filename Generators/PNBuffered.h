#ifndef _PNBUFFERED_H
#define _PNBUFFERED_H 1
/********************************************************************************
 *                                                                              *
 * This class generates PN data using PN polynomial as described in Ziemer      *
 * and Peterson.                                                                *
 * Notes:                                                                       *
 * 1. Both SSRG and MSRG configurations as in Fig. 8-5,6 in Ziemer and          *
 * Peterson are implemented.                                                    *
 * 2. diagram:                                                                  *
 *              SSRG (W-CDMA):                                                  *
 *    ------- + <----- + <--------                                              *
 *   |        ^        ^          |                                             *
 *   |        | g[n]   | g[n-1]     |g[0]                                       *
 *   -> a[n] -> a[n-1] -> ... -> a[0]                                           *
 *                                                                              *
 *              MSRG (IS-95):                                                   *
 *    ---------------------------------------                                   *
 *   |        |          |         |        |                                   *
 *   |        | g[1]     | g[n-2]  |g[n-1]  |                                   *
 *   -> a[n]->+ a[n-1]-> +-> ... -> + ->a[0]---->                               *
 *                                                                              *
 * NOTE: This class uses the PNGenerator class to do the work, and then saves   *
 *       the output in a buffer for fast storage.                               *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNBuffered.h                                *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 01/21/00 - Started.                                                      *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif

class   PNGenerator;

class PNBuffered
{
protected:
//
// Private methods
//
  PNGenerator   *_pnGenerator;                  // The actual pn generator
  char          *_pnData;                       // Stored PN sequence data
  int           _outputLength;                  // Length of output sequence
  int           _oldLength;                     // Old length of PN sequence
  int           _outputPointer;                 // Index into _pnData array
  LOGICAL       _addTrailingZero;               // YES = add zero at end of PN sequence

public:
  void          generateNewData();              // Generate a new set of PN data

/*********************************
 * Constructors/destructors:     *
 *********************************/
  PNBuffered(int poly);                         // Constructor with input polynomial in octal
  virtual ~PNBuffered();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setPNPolynomial(int poly);              // Sets a new PN polynomial in octal
  void  setPNOffset(int offset);                // Sets the offset of the PN generator
  void  setPNInitialLoad(int load);             // Sets the initial load into PN generator
  void  setGeneratorType(int type);             // MSRG or SSRG
  void  setOutputTapNumber(int number);         // Set tap number to take output
  void  setAddTrailingZero(LOGICAL flag);       // YES = add trailing zero

/*******************************
 * These methods get parameters*
 *******************************/
  int   pnPolynomial();                         // returns the PN polynomial in octal
  int   numberBits();                           // returns number of bits in generator

/********************************
 * These methods generate       *
 * new PN data.                 *
 ********************************/
  int   newData();                              // Gets the next data bit from gen

};

#endif
