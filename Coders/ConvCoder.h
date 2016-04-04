#ifndef _CONV_CODER_H
#define _CONV_CODER_H 1
/***************************************************************************//*!
**  \file         ConvCoder.h
**  \brief        Convolutional coder
**  \description  implements a (n, k, m) convolutional coder.
*                 Where: k = number inputs, n = number outputs, and m = memory size.
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        The convolutional coders and decoders are defined by the output path
**                polynomials in octal.  So, for example, for the GSM Sync channel,
**                G0 = 1 + D^3 + D^4, so that the octal rep of the polynomial would
**                be, 23. Similarly for G1 = 1 + D + D^3 + D^4 -> 33 (octal)
********************************************************************************/

#define PUNCTURED_OUTPUT        -1              // Indicates a punctured bit

#include        "ConvCodeDecode.h"              //!< Base class

class ConvCoder: public ConvCodeDecode
{
protected:
  int   m_coderOutput[MAX_CODER_OUTPUTS];        //!< Outputs from coder
  int   *m_shiftRegister;                        //!< Shift register

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
  ConvCoder(int numberOutputs=2,
            int numberInputs=1,
            int numberStages=5,
            int g1=23,                            //Default to GSM Sync channel
            int g2=33,
            int g3=0,
            int p1=37,
            int p2=37,
            int p3=0);
  ~ConvCoder();

/********************************
 * Initializing coder:          *
 ********************************/
  void          initCoder();

/********************************
 * Setting coder parameters:    *
 ********************************/
  void          setState(const int state);

/********************************
 * Getting coder parameters:    *
 ********************************/
  const int     state() const;

/********************************
 * Getting outputs from coder:  *
 ********************************/
  const int     *oldDataOutput() const;                   //!< Returns old output if above fn returns YES
  const int     *newDataOutput(const int inputData);      //!< Returns output given new data
  const int     *unpuncturedOutput(const int inputData);  //!< Returns (unpunctured) output given new data

};
#endif



