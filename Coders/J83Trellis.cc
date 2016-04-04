/********************************************************************************
 * This subclass of Object implements the trellis coded demodulator for J.83.   *
 * 64 QAM                                                                       *
 * File: J83Trellis.cc                                                          *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 03/29/98 - Started                                                      *
 *                                                                              *
 ********************************************************************************/
#include "J83Trellis.h"
#include <stdio.h>

#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Private Function ###############################
// mapSymbolFromArray -- Find the map symbol given an array of bits.
//
// Input:       c:                      array of bits
//          
// Output:                              map symbol
//
// Notes:
// ############################# Public Function ###############################
int J83Trellis::mapSymbolFromArray(int *c)
{
  int   i;
  int   power_of_2      = 1;
  int   output          = 0;
  for(i=0; i<_qamBits; i++)
  {
    if(c[i])
      output    |= power_of_2;
    power_of_2  <<= 1;
  }
  return output;
}

// ############################# Private Function ###############################
// updateParallelPaths -- Calculate the distances over the parallel paths, and update
//                        the transition trellis state.
//
// Input:       c0:                     Output bit from first convolutional coder
//              c3:                     Output bit from second convolutional coder
//              transitionState:        Pointer to transition trellis state object
//              fromState:              Pointer to from trellis state object
//          
// Output:                              None
//
// Notes:
// 1. transitionState is updated on output.
// ############################# Public Function ###############################
void J83Trellis::updateParallelPaths(int c0, int c3, TrellisState *transitionState, TrellisState *fromState)
{
  register      int     k;
  static        int     b[QAM_BITS], c[QAM_BITS];
  int           map_symbol, source_symbol;
  float         distance, min_distance;
//
// Find the uncoded bits:
//
  b[0]          = c[0]  = c0;
  b[3]          = c[3]  = c3;
  min_distance  = BIG_DISTANCE;
  source_symbol = 0;
  for(k=0; k<_numberParallel; k++)
  {
    b[1]        = c[1]          = k & 1;                // These really should be shifted right, but
    b[2]        = c[2]          = k & 2;                // mapSymbolFromArray() only checks for non-zero
    b[4]        = c[4]          = k & 4;
    b[5]        = c[5]          = k & 8;
//
// Find the QAM symbol corresponding to C0->C5, the corresponding I and Q components, and
// then the distance^2 to the input symbol:
//
    map_symbol  = mapSymbolFromArray(c);
    distance    = _distances[map_symbol];
    if(distance < min_distance)
    {
      min_distance      = distance;
      if(DEBUG)
        source_symbol   = map_symbol;
      else
        source_symbol   = mapSymbolFromArray(b);
      if(min_distance<PRUNE_DISTANCE)
        break;
    }
  } // end for(k=0; k<number_parallel; ...)
//
// Finally, update transition state with input:
//
  transitionState->updateState(fromState, min_distance, source_symbol);

  return;
}

// ############################# Private Function ###############################
// updateParallelPaths -- Calculate the distances over the parallel paths, and update
//                        the transition trellis state.
//
// Input:       c0:                     Output bit from first convolutional coder 1st posn
//              c02:                    Output from first coder 2nd commutator posn
//              c3:                     Output bit from second convolutional coder
//              c02:                    Output from first coder 2nd commutator posn
//              transitionState:        Pointer to transition trellis state object
//              fromState:              Pointer to from trellis state object
//          
// Output:                              None
//
// Notes:
// 1. transitionState is updated on output.
// ############################# Public Function ###############################
void J83Trellis::updateParallelPaths(int c0, int c02, int c3, int c32, TrellisState *transitionState,
                                     TrellisState *fromState)
{
  register      int     k;
  static        int     b[QAM_BITS], c[QAM_BITS];
  int           map_symbol;
  int           min_first_source, min_second_source;
  float         distance, min_1st_distance, min_2nd_distance;
  
  min_1st_distance      = min_2nd_distance      = BIG_DISTANCE;
  min_first_source      = min_second_source     = 0;
//
// Find the uncoded bits:
//
  for(k=0; k<_numberParallel; k++)
  {
    b[1]        = c[1]          = k & 1;                // These really should be shifted right, but
    b[2]        = c[2]          = k & 2;                // mapSymbolFromArray() only checks for non-zero
    b[4]        = c[4]          = k & 4;
    b[5]        = c[5]          = k & 8;
//
// Find the 1st QAM symbol corresponding to C0->C5, the corresponding I and Q components, and
// then the distance^2 to the input symbol:
//
    b[0]                        = c[0]  = c0;
    b[3]                        = c[3]  = c3;
    map_symbol                  = mapSymbolFromArray(c);
    distance                    = _oldDistances[map_symbol];

    if(distance < min_1st_distance)
    {
      min_1st_distance          = distance;
      if(DEBUG)
        min_first_source        = map_symbol;
      else
        min_first_source        = mapSymbolFromArray(b);
    }
//
// Find the 2nd QAM symbol corresponding to C0->C5, the corresponding I and Q components, and
// then the distance^2 to the input symbol:
//
    b[0]                        = c[0]  = c02;
    b[3]                        = c[3]  = c32;
    map_symbol                  = mapSymbolFromArray(c);
    distance                    = _distances[map_symbol];
    if(distance < min_2nd_distance)
    {
      min_2nd_distance          = distance;
      if(DEBUG)
        min_second_source       = map_symbol;
      else
        min_second_source       = mapSymbolFromArray(b);
    }
//
// Find the cumulative distance to both symbols:
//
  } // end for(k=0; k<number_parallel; ...)
  distance              = min_1st_distance + min_2nd_distance;
  transitionState->updateState(fromState, distance, min_first_source, min_second_source);

  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the J83Trellis class.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
J83Trellis::J83Trellis()
{
//
// We may want to make the following variable:
//
  _numberStates         = NUMBER_STATES;
  _statesPerCoder       = STATES_PER_CODER;
  _qamBits              = QAM_BITS;
  _qamPoints            = 1 << _qamBits;
  _numberParallel       = 1 << UNCODED_BITS;
//
// Allocate objects and arrays:
//
  _distances            = new float[_qamPoints];
  _oldDistances         = new float[_qamPoints];
  _coder[0]             = new ConvCoder();
  _coder[1]             = new ConvCoder();
  _qamMapper            = new QAMMapper(QAM_J83, _qamBits);
  _states               = new TrellisState[_numberStates];
  _oldStates            = new TrellisState[_numberStates];
//
// Reset trellis:
//
  initTrellis();
  
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the J83Trellis class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
J83Trellis::~J83Trellis()
{
  delete [] _distances;
  delete [] _oldDistances;
  delete _coder[0];
  delete _coder[1];
  delete _qamMapper;
  delete [] _states;
  delete [] _oldStates;
  
  return;
}

// ############################# Public Function ###############################
// initTrellis -- Initializes the trellis.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void J83Trellis::initTrellis()
{
  int   i;

  _punctureCounter      = _resetCounter = 0;
  for(i=0; i<_numberStates; i++)
    _states[i].initState();
  return;
}

// ############################# Public Function ###############################
// normalizeDistances -- Normalize the accumulated distances of all of the states.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// 1. This function should be periodically called to avoid over/under flows.
// ############################# Public Function ###############################
void J83Trellis::normalizeDistances()
{
  int   i;
  float min_distance, distance;

  min_distance  = _states[0].accumulatedDistance();
  for(i=1; i<_numberStates; i++)
  {
    distance    = _states[i].accumulatedDistance();
    if(distance < min_distance)
      min_distance      = distance;
  }
  if(min_distance > 0.)
  {
    for(i=0; i<_numberStates; i++)
    {
      _states[i].normalizeDistanceBy(min_distance);
    }
  }

  return;
}

// ############################# Public Function ###############################
// newInput -- Update the trellis given the new constellation input.
//
// Input:       iData:                  I value of constellation point
//              qData:                  Q value of constellation point
//          
// Output:                      None
//
// Notes:
// 1. There are two convolutional encoders, so we need to set their states
//    separately.
// ############################# Public Function ###############################
void J83Trellis::newInput(float iData, float qData)
{
  int           i, j, number_inputs;
  int           first_coder_state, second_coder_state;
  int           input_state, transition_state;
  const int     *coder_output;
  int           c0, c3, c0_second, c3_second;                   // QAM mapping bits
  int           b0, b3;                                         // source data bits
  float         i_map, q_map, delta_i, delta_q;
  TrellisState  *old_state_object, *transition_state_object;

//
// Calculate and save the distances to the input constellation point:
//
  for(i=0; i<_qamPoints; i++)
  {
    _qamMapper->iAndQForBitMap(i, &i_map, &q_map);
    delta_i             = i_map - iData;
    delta_q             = q_map - qData;
    _distances[i]       = delta_i*delta_i + delta_q*delta_q;
  }

//
// Check to see if we need to store the distances for next call due to unpunctured output:
// If so, return for next call to this function
//
  if(_punctureCounter == FIRST_COMMUTATOR_CNT)
  {
    _punctureCounter++;
    _resetCounter++;
    for(i=0; i<_qamPoints; i++)
      _oldDistances[i]  = _distances[i];
    return;
  }
//
// Now, save the states by copying:
//
  for(i=0; i<_numberStates; i++)
  {
    _oldStates[i]       = _states[i];
    _states[i].initForNewStep();
  }

//
// Loop over all the states finding outputs from each state:
//
  number_inputs         = 1<<CONV_CODER_INPUTS;
  input_state           = 0;
  for(second_coder_state=0; second_coder_state<_statesPerCoder; second_coder_state++)
  {
    for(first_coder_state=0; first_coder_state<_statesPerCoder; first_coder_state++)
    {
      old_state_object  = &_oldStates[input_state];             // Save the old state for efficiency
       for(j=0; j<number_inputs; j++)
       {
//
// Set the states for the coder:
//
        _coder[0]->setState(first_coder_state);
        _coder[1]->setState(second_coder_state);

//
// Get the outputs from the coder, for C0 and C3 bits:
//
        b0              = j&1;
        b3              = (j&2)>>1;
        coder_output    = _coder[0]->unpuncturedOutput(b0);
        c0              = coder_output[FIRST_COMMUTATOR];
        c0_second       = coder_output[SECOND_COMMUTATOR];
        coder_output    = _coder[1]->unpuncturedOutput(b3);
        c3              = coder_output[FIRST_COMMUTATOR];
        c3_second       = coder_output[SECOND_COMMUTATOR];
//
// Get the transition (to) state:
//
        transition_state        = _coder[1]->state();
        transition_state        <<= _coder[1]->numberStages();
        transition_state        |= _coder[0]->state();
        transition_state_object = &_states[transition_state];   // Save the transition state for efficiency

//
// Calculate the distances for parallel paths, and update the transition state:
// When _puncturePtr == SECOND_COMMUTATOR_CNT, we have two outputs (both positions of commutator)
// for each input.  Therefore we need to send in c0_second and c3_second
//
        if(_punctureCounter==SECOND_COMMUTATOR_CNT)
          updateParallelPaths(c0, c0_second, c3, c3_second, transition_state_object, old_state_object);
        else
          updateParallelPaths(c0_second, c3_second, transition_state_object, old_state_object);
      } // End for (j=0...)
      input_state++;
      
    }   // End for (first_coder_state=0;...)
  }     // End for (second_coder_state=0; ...)
//
// Now, go through all the states, and select the best input path:
//
  for(i=0; i<_numberStates; i++)
    _states[i].findBestInput();
//
// Update counters:
//
  if(++_punctureCounter>SECOND_COMMUTATOR_CNT)
    _punctureCounter    = 0;
  if(++_resetCounter > NORMALIZATION_COUNT)
  {
    _resetCounter       = 0;
    normalizeDistances();
  }

  return;
}

// ############################# Public Function ###############################
// getSymbolAt -- Returns the symbol in the best path at 'history' samples in the
//                past.
//
// Input:       history:                sample time in past
//          
// Output:                              best symbol
//
// Notes:
// ############################# Public Function ###############################
int J83Trellis::getSymbolAt(int history)
{
  int           i, best_state;
  float         distance;
//
// First find the minimum distance state:
//
  _minDistance  = _states[0].accumulatedDistance();
  best_state    = 0;
  for(i=1; i<_numberStates; i++)
  {
    distance    = _states[i].accumulatedDistance();
    if(distance < _minDistance)
    {
      _minDistance      = distance;
      best_state        = i;
    }
  }
//
// We modify history according to whether we had a non-punctured output:
//
  if(_punctureCounter == SECOND_COMMUTATOR_CNT)
    history--;
//
// return the symbol corresponding to the back trace of the best state:
//
  return (_states[best_state].symbolInPast(history));
}

// ############################# Public Function ###############################
// getSymbolAt -- Returns the symbol in the best path at 'history' samples in the
//                past.
//
// Input:       history:                sample time in past
//          
// Output:      iData:                  I component of best symbol
//              qData:                  Q component of best symbol
//
// Notes:
// 1. This function is used for debug.
// ############################# Public Function ###############################
void J83Trellis::getIQAt(int history, float *iData, float *qData)
{
  int           i, best_state, best_symbol;
  float         distance;
//
// First find the minimum distance state:
//
  _minDistance  = _states[0].accumulatedDistance();
  best_state    = 0;
  for(i=1; i<_numberStates; i++)
  {
    distance    = _states[i].accumulatedDistance();

    
    if(distance < _minDistance)
    {
      _minDistance      = distance;
      best_state        = i;
    }
  }
//
// return the symbol corresponding to the back trace of the best state:
//
// We modify history according to whether we had a non-punctured output:
//
  if(_punctureCounter == SECOND_COMMUTATOR_CNT)
    history--;
  best_symbol   = _states[best_state].symbolInPast(history);
  _qamMapper->iAndQForBitMap(best_symbol, iData, qData);
  return;
}
