/***************************************************************************//*!
**  \file         ConvDecoder.cpp
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
#include "ConvDecoder.h"
#include "ConvCoder.h"
#include "TrellisState.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

#define BIT_MAP(a)      ( ((a)==0 ) ? (float)1. : (float)-1. )          // (0,1) -> (1, -1)
#define NEG_SGN01(a)    ( ((a)>=0.) ? 0 :  1 )
#define ABS(a)          ((a) >= 0 ? (a) : (-a))

/********************** Private Function **********************************//**
**  \brief              Return the distance between the state transition output and
**                      the input data
**  \description        Constructor
**  \param[in]          packedBits  = state transition output packed into an int
**  \param[in]          floatData   = input array data to the decoder of size m_numberOutputs
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**            1. The LSB of packedBits corresponds to floatData[0].
**            2. The mapping for the bits is as defined above for BIT_MAP()
**            3. The distance function used depends on the variable, m_distanceType
**
**  \TODO         None
*******************************************************************************/
const float ConvDecoder::findDistance(int packedBits,
                                      const float *floatData)
                                      const
{
  int           i, transition_bit, input_bit;
  float         mapped_bit, distance, difference;
  distance      = 0.;
  switch(m_distanceType)
  {
    case HAMMING_DIST:
      distance  = 0.;
      for(i=0; i<m_numberOutputs; i++)
      {
        transition_bit  = packedBits & 1;
        input_bit       = NEG_SGN01(floatData[i]);
        if(transition_bit != input_bit)
        {
          distance     += 1.;
        }
        packedBits     /= 2;
      }
      break;
    case SQUARE_DIST:
      distance  = 0.;
      for(i=0; i<m_numberOutputs; i++)
      {
        transition_bit  = packedBits & 1;
        mapped_bit      = BIT_MAP(transition_bit);
        difference      = mapped_bit - floatData[i];
        distance       += difference*difference;
        packedBits     /= 2;
      }
      break;
    case ABS_DIST:
      distance  = 0.;
      for(i=0; i<m_numberOutputs; i++)
      {
        transition_bit  = packedBits & 1;
        mapped_bit      = BIT_MAP(transition_bit);
        difference      = mapped_bit - floatData[i];
        distance       += ABS(difference);
        packedBits     /= 2;
      }
      break;
  }
  return distance;
}

/********************** Class Constructor **********************************//**
**  \brief              Constructor
**  \description        Constructor
**  \param[in]          numberOutputs  = number of output bits per input bit
**  \param[in]          numberInputs   = number of input shift registers
**  \param[in]          numberStages   = number of stages in coder
**  \param[in]          g1->g3         = generator polynomials in octal
**  \param[in]          p1->p3         = puncture polynomials in octal
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
ConvDecoder::ConvDecoder(const int numberOutputs,
                         const int numberInputs,
                         const int numberStages,
                         const int g1,
                         const int g2,
                         const int g3,
                         const int p1,
                         const int p2,
                         const int p3)
            :ConvCodeDecode(numberOutputs, numberInputs, numberStages, g1, g2, g3, p1, p2, p3)
{
  setDecodeDepth(DEFAULT_DECODE_DEPTH);
  setDistanceType(SQUARE_DIST);
  m_transitions          = NULL;
  m_states               = NULL;
  m_oldStates            = NULL;
  m_numberStates         = 0;
  m_numberTransitions    = 0;
  m_oldNumberStates      = m_oldDecodeDepth       = -1;
  initCoder();
  return;
}
  
/********************** Class Destructor **********************************//**
**  \brief              Class Destructor
**  \description        Class Destructor
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
ConvDecoder::~ConvDecoder()
{
  delete [] m_transitions;
  delete [] m_states;
  delete [] m_oldStates;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Initializes the decoder to start accepting new input
**  \description        Initializes the decoder to start accepting new input
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_transitions, m_states, m_oldStates, m_oldNumberStates 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvDecoder::initCoder()
{
  int           n, m, i;
  ConvCoder     *conv_coder;            // The coder for calculating the transitions
  const int     *coder_output;
//
// Call super's function:
//
  ConvCodeDecode::initCoder();
//
// Add this:
//
// Allocate the coder for calculating the transitions:
//
  conv_coder            = new ConvCoder(m_numberOutputs,
                                        m_numberInputs,
                                        m_numberStages,
                                        m_g[0],
                                        m_g[1],
                                        m_g[2],
                                        m_p[0],
                                        m_p[1],
                                        m_p[2]);
  conv_coder->initCoder();
//
// Allocate the transitions and trellis states:
//
  m_numberStates         = conv_coder->numberStates();
  if(m_oldNumberStates < m_numberStates)
  {
    delete [] m_transitions;
    delete [] m_states;
    delete [] m_oldStates;
    m_transitions        = new StateTransition[m_numberStates];
    m_states             = new TrellisState[m_numberStates];
    m_oldStates          = new TrellisState[m_numberStates];
    m_oldNumberStates    = m_numberStates;
  }
//
// Initialize the trellis states:
//
  for(n=0; n<m_numberStates; n++)
  {
    m_states->setPathLength(m_decodeDepth);
    m_states->initState();
  }
//
// Get the state transition information:
//
  m_numberTransitions                    = 1 << m_numberInputs;
  for(n=0; n<m_numberStates; n++)
  {
    for(m=0; m<m_numberTransitions; m++)         // This only works for m_numberInputs == 1
    {
      conv_coder->setState(n);                  // Initialize the state to current state
      coder_output                       = conv_coder->unpuncturedOutput(m);
      m_transitions[n].nextState[m]      = conv_coder->state();
      unsigned int shift                 = m_numberOutputs-1;
      m_transitions[n].stateOutputs[m]   = 0;
      for(i=0; i<m_numberOutputs; i++)           // pack output into an int
      {
        int mult                          = 1 << shift;
        m_transitions[n].stateOutputs[m] += coder_output[i] * mult;
        shift--;
      }
    }
  }
  delete conv_coder;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets a new size of the back trace buffer
**  \description        Sets a new size of the back trace buffer
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_decodeDepth 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvDecoder::setDecodeDepth(const int depth)
{
  m_decodeDepth          = MAX(depth, MIN_PATH_LENGTH);
  m_decodeDepth          = MIN(m_decodeDepth, MAX_PATH_LENGTH);
  return;
}

/********************** Public Function **********************************//**
**  \brief              Sets the type of distance type to use
**  \description        Sets the type of distance type to use
**  \param[in]          type    = HAMMING_DIST, SQUARE_DIST, etc
**  \return             None
**   \post         
**  <b>Modified:</b>    m_decodeDepth 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvDecoder::setDistanceType(const int type)
{
  m_distanceType         = type;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Update the trellis given the new constellation input
**  \description        Update the trellis given the new constellation input
**  \param[in]          data    = Array of symbol values, the size of the array needs
**                                to be equal to the number of outputs in the coder.  E.g.,
**                                for a rate 1/2 coder, the size of the array is 2.
**  \return             None
**   \post         
**  <b>Modified:</b>    m_oldStates, m_states 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void ConvDecoder::newInput(const float *data)
{
  int           i, n, output_statem_number;
  int           transition_symbol;
  float         distance;
  TrellisState  *input_state, *output_state;
//
// Save the current value of the trellis states by copying:
//
  for(i=0; i<m_numberStates; i++)
  {
    m_oldStates[i]       = m_states[i];
    m_states[i].initForNewStep();
  }
//
// Loop over all the states updating the distance into each trellis state:
// NOTE: By passing 'i' into updateState, we save the transition bit in the
// back trace.
//
  for(n=0; n<m_numberStates; n++)
  {
    input_state                 = &m_oldStates[n];
    for(i=0; i<m_numberTransitions; i++)
    {
      transition_symbol         = m_transitions[n].stateOutputs[i];
      distance                  = findDistance(transition_symbol, data);
      output_statem_number      = m_transitions[n].nextState[i];
      output_state              = &m_states[output_statem_number];
      output_state->updateState(input_state, distance, i);
    }
  }
//
// After all the new distances have been updated, find the new best input path
// for each state:
//
  for(n=0; n<m_numberStates; n++)
  {
    m_states[n].findBestInput();
  }
  return;
}

/********************** Public Function **********************************//**
**  \brief              getSymbolAt()
**  \description        Returns the symbol in the best path at 'history' samples in the past.
**  \param[in]          history    = # of samples in past
**  \return             None
**   \post         
**  <b>Modified:</b>    m_oldStates, m_states 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int ConvDecoder::getSymbolAt(const int history)
{
  int           i, best_state, bit_output;
  float         distance;
//
// First find the minimum distance state:
//
  m_minDistance  = m_states[0].accumulatedDistance();
  best_state    = 0;
  for(i=1; i<m_numberStates; i++)
  {
    distance    = m_states[i].accumulatedDistance();
    if(distance < m_minDistance)
    {
      m_minDistance      = distance;
      best_state        = i;
    }
  }
//
// We need to get transition corresponding to symbol, this was stored in the back trace
// when updateState() was called above.
//
  bit_output            = m_states[best_state].symbolInPast(history);
//
// return the symbol corresponding to the back trace of the best state:
//
  return bit_output;
}

