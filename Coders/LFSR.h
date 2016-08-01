#ifndef _LFSR_H
#define _LFSR_H 1
/********************************************************************************
 * This subclass of Object implements a linear feedback shift register.         *
 * The generator polynomial for the degree r LFSR is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
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

class LFSR
{
protected:
  int   _degree;                                // Degree of generator poly
  int   _oldDegree;                             // Old size for allocs and deletes
  int   _oldShiftStages;                        // Old size for allocs and deletes
  int   _oldOutputSize;                         // Old size for allocs and deletes
  int   _outputSize;                            // Size of output array
  BIT   *_generatorPoly;                        // _generatorPoly[0] = g0, _generatorPoly[1] = g1
  BIT   *_lfsrOutput;                           // Outputs from LFSR
  BIT   *_shiftRegister;                        // Shift register

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
  LFSR(int degree=0, const BIT *poly=NULL);
  ~LFSR();

/********************************
 * Initializing LFSR:           *
 ********************************/
  void          initLFSR(BIT registerLoad=0);

/********************************
 * Setting LFSR parameters:    *
 ********************************/
  void          setGeneratorPoly(int degree, const BIT *poly);

/********************************
 * Getting LFSR parameters:    *
 ********************************/
  int           degree()                {return _degree;}
  const         BIT *shiftRegister()    {return _shiftRegister;}

/********************************
 * Getting outputs from LFSR:   *
 ********************************/
  BIT           coderOutput(BIT inputBit);              // Returns coder output given new data
  BIT           syndromeOutput(BIT inputBit);

};
#endif


