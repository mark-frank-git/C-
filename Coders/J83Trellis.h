#ifndef _J83_TRELLIS_H
#define _J83_TRELLIS_H 1
/********************************************************************************
 * This subclass of Object implements the trellis coded demodulator for J.83.   *
 * 64 QAM                                                                       *
 * File: J83Trellis.h                                                           *
 *                                                                              *
 *                                                                              *
 ********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif

#define         STATES_PER_CODER        16                      // 16 states per coder
#define         NUMBER_STATES           (STATES_PER_CODER * STATES_PER_CODER)
#define         QAM_BITS                6                       // 64 QAM
#define         UNCODED_BITS            4                       // 4 of 6 are uncoded
#define         NORMALIZATION_COUNT     30                      // Normalize distance count
#define         FIRST_COMMUTATOR_CNT    3                       // puncture counter to take 1st out
#define         SECOND_COMMUTATOR_CNT   4                       // puncture counter to take 2nd out
#define         FIRST_COMMUTATOR        0                       // 1st commutator position
#define         SECOND_COMMUTATOR       1                       // 2nd commutator position
#define         CONV_CODER_INPUTS       2                       // 2 bits input to conv coder
#define         PRUNE_DISTANCE          4.                      // Distance for pruning states
#define         TWO_PRUNE_DISTANCE      8.                      // Twice the above distance
#define         DEBUG                   1                       // 1 to output map symbols, else input symbols

#include        <Generators/QAMMapper.h>
#include        "ConvCoder.h"
#include        "TrellisState.h"

class J83Trellis
{
protected:
  int           _numberStates;          // Number of states in the trellis, normally = 256
  int           _statesPerCoder;        // Number of states per coder, normally = 16
  int           _qamBits;               // Number of QAM bits, normally = 6
  int           _qamPoints;             // 2**_qamBits;
  int           _punctureCounter;       // Counter for puncturing of the convolutional coders
  int           _resetCounter;          // Counter for normalizing distances
  int           _numberParallel;        // Number of parallel paths for each transition
  float         _minDistance;           // Minimum distance calculated in getSymbolAt (for debug)
  float         *_distances;            // Calculated distances
  float         *_oldDistances;         // Previous symbol calculated distances
  ConvCoder     *_coder[2];             // Convolutional coder
  QAMMapper     *_qamMapper;            // QAM symbol mapper
  TrellisState  *_states;               // Trellis states
  TrellisState  *_oldStates;            // Copy of trellis states

//
// Private functions:
//
  inline        int     mapSymbolFromArray(int *c);     // Find the map symbol from the C0->C5 bits
  inline        void    updateParallelPaths(int c0, int c3, TrellisState *transitionState, TrellisState *fromState);
                                                        // Find the parallels path distances and update transition state
  inline        void    updateParallelPaths(int c0, int c02, int c3, int c32, TrellisState *transitionState,
                                     TrellisState *fromState);
                                                        // Similar to above, except for unpunctured (double) outputs
//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  J83Trellis();
  ~J83Trellis();

/********************************
 * Initializing the trellis:    *
 ********************************/
  void          initTrellis();
  void          normalizeDistances();

/********************************
 * Updating the trellis:        *
 ********************************/
  void          newInput(float iData, float qData);                     // Inputs new symbol to trellis

/********************************
 * Getting output from the      *
 * trellis:                     *
 ********************************/
  int           getSymbolAt(int history);                               // Returns best symbol at 'history' samples ago
  void          getIQAt(int history, float *iData, float *qData);       // Best I/Q point for debug
  float         getMinDistance()                {return _minDistance;}


};
#endif


