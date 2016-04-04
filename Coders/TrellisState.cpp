/***************************************************************************//*!
**  \file         TrellisState.cpp
**  \brief        Trellis for viterbi back trace
**  \description  implements a single state for a trellis demodulator.
**                Should be able to be used for any type of Viterbi decoding.
**                Each state keeps the back trace path, along with the accumulated distance
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        None
********************************************************************************/
#include "TrellisState.h"
#include <stdio.h>

#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

/********************** Class Constructor **********************************//**
**  \brief              Constructor
**  \description        Constructor
**  \param[in]          pathLength      = back trace path length
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
TrellisState::TrellisState(const int pathLength)
{
  m_backTracePath        = NULL;
  m_oldLength            = 0;
  setPathLength(pathLength);
  initState();
  initForNewStep();
  
  return;
}

/********************** Class Constructor **********************************//**
**  \brief              Constructor
**  \description        Constructor
**  \param[in]          trellis         = A previously allocated TrellisState
**  \return             None
**   \post         
**  <b>Modified:</b>    None
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
TrellisState::TrellisState(const TrellisState& trellis)
{
  int   i;
  const int *back_trace;
//
// copy instance variables:
//
  m_backTracePath        = NULL;
  m_oldLength            = 0;
  setPathLength(trellis.pathLength());
  back_trace            = trellis.backTracePath();
  for(i=0; i<m_pathLength; i++)
  {
    m_backTracePath[i]   = back_trace[i];
  }
  m_accumulatedDistance  = trellis.accumulatedDistance();
  m_pathBufferPtr        = trellis.pathBufferPtr();

  initForNewStep();
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
TrellisState::~TrellisState()
{
  delete [] m_backTracePath;
  return;
}

/********************** Public Function **********************************//**
**  \brief              = operator
**  \description        Sets the TrellisState equal to the input TrellisState
**  \param[in]          trellis    = input TrellisState
**  \return             the result TrellisState
**   \post         
**  <b>Modified:</b>    this
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
TrellisState& TrellisState::operator = (const TrellisState& trellis)
{
  int   i, newm_path;
  const int *back_trace;

  if( this == &trellis )                        // Check for x=x
    return *this;
//
// copy instance variables:
//
  newm_path              = trellis.pathLength();
  setPathLength(newm_path);
  back_trace            = trellis.backTracePath();
  for(i=0; i<m_pathLength; i++)
    m_backTracePath[i]   = back_trace[i];
  m_accumulatedDistance  = trellis.accumulatedDistance();
  m_pathBufferPtr        = trellis.pathBufferPtr();

  initForNewStep();

  return *this;
}

/********************** Public Function **********************************//**
**  \brief              Init State
**  \description        Initializes the state
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_pathBufferPtr, m_accumulatedDistance, m_backTracePath 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void TrellisState::initState()
{
  int   i;
  
  m_pathBufferPtr        = 0;
  m_accumulatedDistance  = 0.;
  for(i=0; i<m_pathLength; i++)
    m_backTracePath[i]           = 0;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Init State
**  \description        Initializes the state for a new Viterbi step
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_minimumDistance, m_bestInputSymbol, m_bestSecondSymbol, m_bestInputState 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void TrellisState::initForNewStep()
{
  m_minimumDistance      = BIG_DISTANCE;
  m_bestInputSymbol      = -1;
  m_bestSecondSymbol     = -1;
  m_bestInputState       = NULL;
  return;
}

/********************** Public Function **********************************//**
**  \brief              Set path length
**  \description        Sets a new size of the back trace buffer
**  \param[in]          length    = new pathLength size
**  \return             None
**   \post         
**  <b>Modified:</b>    m_pathLength, m_backTracePath, m_oldLength 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void TrellisState::setPathLength(const int length)
{
  m_pathLength           = MAX(length, MIN_PATH_LENGTH);
  m_pathLength           = MIN(m_pathLength, MAX_PATH_LENGTH);
  if(m_oldLength < m_pathLength)
  {
    delete [] m_backTracePath;
    m_backTracePath      = new int[m_pathLength];
    m_oldLength          = m_pathLength;
  }
  return;
}

/********************** Public Function **********************************//**
**  \brief              Symbol from the past
**  \description        Returns a symbol from the back trace buffer
**  \param[in]          pastLength    = time in the past
**  \return             None
**   \post         
**  <b>Modified:</b>    None 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
const int TrellisState::symbolInPast(int pastLength)
                                    const
{
  int           read_ptr;
  
  pastLength    = MAX(0, pastLength);                   // error check
  read_ptr      = m_pathBufferPtr - pastLength;
  while(read_ptr < 0)
    read_ptr    += m_pathLength;
  return m_backTracePath[read_ptr];
}

/********************** Public Function **********************************//**
**  \brief              Normalize distance
**  \description        Normalize the accumulated distance by the input subtractor
**  \param[in]          normalization    = Normalization factor
**  \return             None
**   \post         
**  <b>Modified:</b>    m_accumulatedDistance 
**  <b>Notes:</b>       This function should be periodically called to avoid over/under flows
**
**  \TODO         None
*******************************************************************************/
void TrellisState::normalizeDistanceBy(const float normalization)
{
  m_accumulatedDistance  -= normalization;
  m_accumulatedDistance  = MIN(m_accumulatedDistance, BIG_DISTANCE);
  return;
}

/********************** Public Function **********************************//**
**  \brief              Update state
**  \description        Update the state based on the minimum distance path into this
**                      state
**  \param[in]          inputStates         = An array of input states
**  \param[in]          transitionDistances = The corresponding state transition distances
**  \param[in]          transitionSymbols   = The symbols corresponding to the transitions
**  \param[in]          numberInputs        = The size of the above array
**  \return             None
**   \post         
**  <b>Modified:</b>    m_minimumDistance, m_bestInputSymbol, m_bestInputState 
**  <b>Notes:</b>       1. After calling this function for all inputs, call findBestInput()
**                      2. For convolutional decoder, it might be better to set transitionSymbol to
**                         input bit
**
**  \TODO         None
*******************************************************************************/
void TrellisState::updateState(TrellisState *inputStates,
                               const float *transitionDistances,
                               const int *transitionSymbols,
                               const int numberInputs)
{
  int           i;
  float distance;
//
// Update the minimum distance input:
//
  for(i=0; i<numberInputs; i++)
  {
    distance            = inputStates[i].accumulatedDistance() + transitionDistances[i];
    if(distance < m_minimumDistance)
    {
      m_minimumDistance  = distance;
      m_bestInputSymbol  = transitionSymbols[i];                 // Save the best input
      m_bestInputState   = &inputStates[i];
    }
  }
  return;
}

/********************** Public Function **********************************//**
**  \brief              Update state
**  \description        Update the state based on the minimum distance path into this
**                      state
**  \param[in]          inputState          = pointer to an input state
**  \param[in]          transitionDistance  = The corresponding state transition distance
**  \param[in]          transitionSymbol    = The symbol corresponding to the transition
**  \return             None
**   \post         
**  <b>Modified:</b>    m_minimumDistance, m_bestInputSymbol, m_bestInputState 
**  <b>Notes:</b>       1. After calling this function for all inputs, call findBestInput()
**                      2. For convolutional decoder, it might be better to set transitionSymbol to
**                         input bit
**
**  \TODO         None
*******************************************************************************/
void TrellisState::updateState(TrellisState *inputState,
                               const float transitionDistance,
                               const int transitionSymbol)
{
  float distance;
  if(inputState == NULL)
  {
    return;
  }
//
// Update the minimum distance input:
//
  distance      = inputState->accumulatedDistance() + transitionDistance;
  if(distance < m_minimumDistance)
  {
      m_minimumDistance  = distance;
      m_bestInputSymbol  = transitionSymbol;                     // Save the best input
      m_bestInputState   = inputState;
  }
  return;
}

/********************** Public Function **********************************//**
**  \brief              Update state
**  \description        Update the state based on the minimum distance path into this
**                      state
**  \param[in]          inputState          = pointer to an input state
**  \param[in]          transitionDistance  = The corresponding state transition distance
**  \param[in]          transitionSymbol    = The symbol corresponding to the transition
**  \return             None
**   \post         
**  <b>Modified:</b>    m_minimumDistance, m_bestInputSymbol, m_bestInputState 
**  <b>Notes:</b>       1. After calling this function for all inputs, call findBestInput()
**                      2. This function is similar to the above, except it handles the case of two output
**                         symbols per transition
**
**  \TODO         None
*******************************************************************************/
void TrellisState::updateState(TrellisState *inputState,
                               const float transitionDistance,
                               const int firstSymbol,
                               const int secondSymbol)
{
  float distance;

  if(inputState == NULL)
  {
    return;
  }
//
// Update the minimum distance input:
//
  distance      = inputState->accumulatedDistance() + transitionDistance;
  if(distance < m_minimumDistance)
  {
      m_minimumDistance  = distance;
      m_bestInputSymbol  = firstSymbol;                  // Save the best input
      m_bestSecondSymbol = secondSymbol;
      m_bestInputState   = inputState;
  }
  return;
}


/********************** Public Function **********************************//**
**  \brief              Find best input
**  \description        Update the state based on the current minimum distance input
**  \param[in]          None
**  \return             None
**   \post         
**  <b>Modified:</b>    m_pathBufferPtr, m_backTracePath, m_accumulatedDistance 
**  <b>Notes:</b>       None
**
**  \TODO         None
*******************************************************************************/
void TrellisState::findBestInput()
{
  int   i;
  const int *back_trace;
  if(m_bestInputState != NULL)
  {
//
// Copy the best input state to ourself:
//
    setPathLength(m_bestInputState->pathLength());
    m_pathBufferPtr      = m_bestInputState->pathBufferPtr();
    back_trace           = m_bestInputState->backTracePath();
    for(i=0; i<m_pathLength; i++)
      m_backTracePath[i] = back_trace[i];
//
// Update the accumulated distance, and the transition symbol:
//
    m_accumulatedDistance                = m_minimumDistance;
    m_pathBufferPtr++;
    if(m_pathBufferPtr == m_pathLength)
    {
      m_pathBufferPtr                    = 0;
    }
    m_backTracePath[m_pathBufferPtr]     = m_bestInputSymbol;
//
// Check for a second symbol:
//
    if(m_bestSecondSymbol > -1)
    {
      m_pathBufferPtr++;
      if(m_pathBufferPtr == m_pathLength)
      {
        m_pathBufferPtr                  = 0;
      }
      m_backTracePath[m_pathBufferPtr]   = m_bestSecondSymbol;
    }
  }
  else
  {
    m_accumulatedDistance                = BIG_DISTANCE;
  }
  
  return;
}

