#ifndef _PNGENERATOR_H
#define _PNGENERATOR_H 1
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
 *   |        | g[n]   | g[n-1]   |g[0]                                         *
 *   -> a[n] -> a[n-1] -> ... -> a[0]                                           *
 *                                                                              *
 *              MSRG (IS-95):                                                   *
 *    ---------------------------------------                                   *
 *   |        |          |         |        |                                   *
 *   |        | g[1]     | g[n-2]  |g[n-1]  |                                   *
 *   -> a[n]->+ a[n-1]-> +-> ... -> + ->a[0]---->                               *
 *                                                                              *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNGenerator.h                               *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 02/06/00 - Started.                                                      *
 ********************************************************************************/
#define MSRG_TYPE       0
#define SSRG_TYPE       1

#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_PN_LOAD 1

class PNGenerator
{
protected:
  short         *_g, *_a;                       // PN generator polynomials
  int           _polynomialDegree;              // PN polynomial degree.
  int           _pnPolynomial;                  // PN generator polynomial, in octal digits
  int           _pnInitialLoad;                 // Initial load of registers
  int           _pnOffset;                      // Offset of PN generator from 0
  int           _generatorType;                 // MSRG or SSRG
  int           _outputTapNumber;               // The shift register tap to take output
//
// Private methods
//
  int   newMSRGData();                          // Gets the next data bit from gen (using MSRG)
  int   newSSRGData();                          // Gets next data bit using SSRG.

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  PNGenerator(int poly);                        // Constructor with input polynomial in octal
  virtual ~PNGenerator();

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

/*******************************
 * These methods get parameters*
 *******************************/
  int   pnPolynomial()  {return _pnPolynomial;}
  int   numberBits();                           // returns number of bits in generator

/********************************
 * These methods generate       *
 * new PN data.                 *
 ********************************/
  int   newData();                              // Gets the next data bit from gen


};

#endif
