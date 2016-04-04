#ifndef _CONV_DECODER_H
#define _CONV_DECODER_H 1
/***************************************************************************//*!
**  \file         ConvDecoder.h
**  \brief        Convolutional decoder
**  \description  implements a (n, k, m) convolutional decoder.
*                 Where: k = number inputs, n = number outputs, and m = memory size.
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        The convolutional coders and decoders are defined by the output path
**                polynomials in octal.  So, for example, for the GSM Sync channel,
**                G0 = 1 + D^3 + D^4, so that the octal rep of the polynomial would
**                be, 23. Similarly for G1 = 1 + D + D^3 + D^4 -> 33 (octal)
********************************************************************************/
#define DEFAULT_DECODE_DEPTH    5               // Back trace depth

#define HAMMING_DIST            0               // distance types
#define SQUARE_DIST             1
#define ABS_DIST                2

#include                "ConvCodeDecode.h"

struct StateTransition                          //!< State transition information for Trellis
{
  int   nextState[MAX_CODER_TRANSITIONS];       //!< Next state for given input
  int   stateOutputs[MAX_CODER_TRANSITIONS];    //!< Coder output for given input, the outputs are packed
                                                //!< into an int using stateOutput = output[0] + output[1] << 1 +..
};

class                   TrellisState;           //!< Class prototypes

class ConvDecoder: public ConvCodeDecode
{
protected:
  int                   m_numberStates;          //!< Number of states in the trellis
  int                   m_numberTransitions;     //!< Number of transitions between states in the trellis
  int                   m_oldNumberStates;       //!< Old number of states for reducing calls to new
  int                   m_oldDecodeDepth;        //!< Old decode depth for reducing calls to new
  int                   m_decodeDepth;           //!< Depth of the decoder (trace back history)
  int                   m_distanceType;          //!< HAMMING_DIST (hard), SQUARE_DIST, ABS_DIST, etc.
  float                 m_minDistance;           //!< Minimum distance calculated in getSymbolAt (for debug)
  StateTransition       *m_transitions;          //!< State transition information from coder
  TrellisState          *m_states;               //!< The states for the trellis
  TrellisState          *m_oldStates;            //!< The old states for updating
//
// Private functions:
//
  const float  findDistance(int packedBits,      //!< Find distance between input symbol and trellis
                            const float *floatData)
                            const;

//
// Public functions:
//
public:

/********************************
 * Constructors, destructors    *
 ********************************/
  ConvDecoder(const int numberOutputs=2,
              const int numberInputs=1,
              const int numberStages=5,
              const int g1=23,                            //Default to GSM Sync channel
              const int g2=33,
              const int g3=0,
              const int p1=37,
              const int p2=37,
              const int p3=0);
  ~ConvDecoder();

/********************************
 * Initializing decoder:        *
 ********************************/
  void          initCoder();

/********************************
 * Setting decoder parameters:  *
 ********************************/
  void          setDecodeDepth(const int depth);
  void          setDistanceType(const int type);

/********************************
 * Getting decoder parameters:  *
 ********************************/
  const int     decodeDepth()   const {return m_decodeDepth;}
  const int     distanceType()  const {return m_distanceType;}

/********************************
 * Updating the trellis:        *
 ********************************/
  void          newInput(const float *inputData);                // Inputs new symbol to trellis

/********************************
 * Getting output from the      *
 * trellis:                     *
 ********************************/
  const int     getSymbolAt(const int history);                 // Returns best symbol at 'history' samples ago

};

#endif

