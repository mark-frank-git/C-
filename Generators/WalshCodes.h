#ifndef _WALSH_CODES_H
#define _WALSH_CODES_H 1
/********************************************************************************
 * This subclass of Object generates Walsh functions (codes) and stores them    *
 * in arrays for output.                                                        *
 *                                                                              *
 * File: WalshCodes.h                                                           *
 *                                                                              *
 * Notes:                                                                       *
 * 1. See Lewis Franks, "Signal Theory" for Walsh function generation formulae. *
 * 2. The Walsh code functions have values of 0 or 1.                           *
 * 3. These are really indexed according to the Hadamard matrix as in IS-95.    *
 *                                                                              *
 ********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_NUMBER_CODES    64              // # of codes for IS-95
#define MAX_CODES               1024            // Max value for _numberCodes

class WalshCodes
{
protected:
  int           _numberCodes;                   // Number of codes to be generated
  int           _codeLength;                    // Length of the codes
  char          **_walshCodes;                  // The generated codes

//
// Private functions:
//
  int   log2(int n);
  int   ipow(int i, int n);
  int   *generateHadamardMatrix(int *hn, int n);
  void  generateCodes();

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  WalshCodes(int numberCodes=DEFAULT_NUMBER_CODES);
  ~WalshCodes();

/********************************
 * Getting a code:              *
 ********************************/
  const char    *codeForIndex(int index);

/********************************
 * Printing the codes:          *
 ********************************/
  void          printCodes();

};
#endif

