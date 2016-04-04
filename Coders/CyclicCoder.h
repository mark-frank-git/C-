#ifndef _CYCLIC_CODER_H
#define _CYCLIC_CODER_H 1
/********************************************************************************
 * This subclass of Object implements a cyclic coder.                           *
 * The generator polynomial for the degree r code is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/22/02 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include <stdio.h>

#define MAX_POLY_DEGREE         200

#ifndef LOGICAL
#define LOGICAL char
#endif

#ifndef BIT
#define BIT int
#endif

#ifndef NO
#define NO      0
#endif

class   LFSR;                                   // Does the work

class CyclicCoder
{
protected:
  int   _oldOutputSize;                         // Old size for allocs and deletes
  int   _outputSize;                            // Size of output array
  int   _cycleCounts;                           // Cycle count estimator
  BIT   *_coderOutput;                          // Outputs from coder

  LFSR  *_lfsr;                                 // linear feedback shift register
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
  CyclicCoder(int degree=0, BIT *poly=NULL);
  ~CyclicCoder();

/********************************
 * Initializing coder:          *
 ********************************/
  void          initCoder(BIT registerLoad=0);

/********************************
 * Setting coder parameters:    *
 ********************************/
  void          setGeneratorPoly(int degree, const BIT *poly);

/********************************
 * Getting coder parameters:    *
 ********************************/

/********************************
 * Getting outputs from coder:  *
 ********************************/
  const BIT     *blockOutput(const BIT *inputBits, int inputLength, LOGICAL flipParity=NO);
                                                        // Returns a block of output
  int           outputSize()            {return _outputSize;}
  int           cycleCounts()           {return _cycleCounts;}

/********************************
 * Error detection:             *
 ********************************/
  LOGICAL       bitErrors(const BIT *inputBits, int inputLength, LOGICAL flipParity=NO);

};
#endif


