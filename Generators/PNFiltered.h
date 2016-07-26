#ifndef _PNFILTERED_H
#define _PNFILTERED_H 1
/********************************************************************************
 *                                                                              *
 * This class generates PN data using PN polynomial as described in Ziemer      *
 * and Peterson.  However, this class generates filtered PN data in the shape   *
 * of triangular waves.
 * Notes:                                                                       *
 * 1. SSRG configuration is in Fig. 8-6 in Ziemer and Peterson                  *
 * 2. diagram:                                                                  *
 *                                                                              *
 *    ------- + <----- + <--------                                              *
 *   |        ^        ^          |                                             *
 *   |        | g[1]   | g[2]     |                                             *
 *   -> a[n] -> a[n-1] -> ... -> a[0]                                           *
 *                                                                              *
 *                                                                              *
 * File: /User/frank/C++/Generators/PNFiltered.h                                *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 02/06/01 - Started.                                                      *
 ********************************************************************************/

#ifndef LOGICAL
#define LOGICAL char
#endif

#define PLUS1_BIT       0               // defines the filtered bits
#define MINUS1_BIT      1
#define MINUS_SLOPE_BIT 2
#define PLUS_SLOPE_BIT  3

class   PNGenerator;

class PNFiltered
{
protected:
  PNGenerator   *_pnGenerator;                  // PN generator
  int           _numberBits;
  char          *_filteredBits;                 // The stored bits
//
// Private methods
//
  void          generateData();                 // Pre-generate filtered data
  void          generateData(const char *pnSequence, int length);

public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  PNFiltered(int poly);
  PNFiltered(const char *pnData, int length);
  virtual ~PNFiltered();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initGenerator();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setPNPolynomial(int poly);

/********************************
 * These methods get output from*
 * the generator.               *
 ********************************/
  float getEarlyLateOutput(int index, float offset);
  float getOnTimeOutput(int index);


};

#endif
