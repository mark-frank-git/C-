#ifndef _FIRE_DECODER_H
#define _FIRE_DECODER_H 1
/********************************************************************************
 * This subclass of Object implements a decoder for a Fire code.                *
 * The generator polynomial for the degree r code is assumed to be of the form: *
 * g(x) = go + g1x + g2*x^2 + g3*x^3 + .. x^r.                                  *
 *                                                                              *
 *                                                                              *
 * Notes:                                                                       *
 * 1. Based on the algorithm in Wicker's channel coding book.                   *
 ********************************************************************************/
#include <stdio.h>

#define MAX_BURST_LENGTH        80

#ifndef LOGICAL
#define LOGICAL char
#endif

#ifndef BIT
#define BIT int
#endif

#ifndef NO
#define NO      0
#endif

class   LFSR;                                   // Class prototype

class FireDecoder
{
protected:
  int   _burstLength;                           // b, degree of lfsr1 is 2b-1
  int   _oldDegree;                             // Old size for allocs and deletes
  int   _oldShiftStages;                        // Old size for allocs and deletes
  int   _oldOutputSize;                         // Old size for allocs and deletes
  int   _outputSize;                            // Size of output array
  BIT   *_generatorPoly2;                       // _generatorPoly[0] = g0, _generatorPoly[1] = g1
  BIT   *_decoderOutput;                        // Outputs from decoder
  BIT   *_shiftRegister;                        // Shift register
  LFSR  *_lfsr1;                                // The syndrome calculator for LFSR 1
  LFSR  *_lfsr2;                                // The syndrome calculator for LFSR 2

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
   FireDecoder(int burstLength=0, const BIT *generatorPoly=NULL, int degree=0);
  ~FireDecoder();

/********************************
 * Initializing coder:          *
 ********************************/
  void          initDecoder();

/********************************
 * Setting coder parameters:    *
 ********************************/
  void          setBurstLength(int burstLength);
  void          setGeneratorPoly2(int degree, const BIT *poly);

/********************************
 * Getting coder parameters:    *
 ********************************/
  int           burstLength()           {return _burstLength;}

/********************************
 * Getting outputs from coder:  *
 ********************************/
  void          shiftInBit (BIT inputBit);              // Shifts a bit into the decoder
  void          decodeBits(const BIT *inputBits, int inputLength);
                                                        // Returns a block of output
  const         BIT *decoderOutput()    {return _decoderOutput; }
  int           outputSize()            {return _outputSize;    }


};
#endif


