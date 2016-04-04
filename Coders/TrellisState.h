#ifndef _TRELLIS_STATE_H
#define _TRELLIS_STATE_H 1
/***************************************************************************//*!
**  \file         TrellisState.h
**  \brief        Trellis for viterbi back trace
**  \description  implements a single state for a trellis demodulator.
**                Should be able to be used for any type of Viterbi decoding.
**                Each state keeps the back trace path, along with the accumulated distance
**  \project      EE 598 class
**  \author       Mark Frank
**  \Notes        None
********************************************************************************/
#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_PATH_LENGTH     24
#define MIN_PATH_LENGTH         1
#define MAX_PATH_LENGTH         128

#define BIG_DISTANCE            (float)1.e20     // For initializing m_minimumDistance

class TrellisState
{
protected:
  int           m_pathBufferPtr;                 //!< The circular buffer pointer for the back trace
  int           m_pathLength;                    //!< The size of the back trace buffer
  int           m_oldLength;                     //!< Old size to reduce calls to new
  int           m_bestInputSymbol;               //!< Best input symbol for current update
  int           m_bestSecondSymbol;              //!< Best input second symbol when there are two outputs/transition
  int           *m_backTracePath;                //!< The back trace of the path
  float         m_minimumDistance;               //!< Minimum distance for current update
  float         m_accumulatedDistance;           //!< The accumulated distance for the min path to this state
  TrellisState  *m_bestInputState;               //!< The best input state to this one

//
// Public functions:
//
public:
    
/********************************
 * Constructors, destructors    *
 ********************************/
  TrellisState(const int pathLength=DEFAULT_PATH_LENGTH);
  TrellisState(const TrellisState &trellis);
  ~TrellisState();

/********************************
 * Define the equals operator   *
 ********************************/
  TrellisState&   operator =  (const TrellisState& y);      //!< copy constructor

/********************************
 * Initializing the state:      *
 ********************************/
  void          initState();                                //!< Initialize state
  void          initForNewStep();                           //!< Initialize for new step
  
/********************************
 * Setting parameters:          *
 ********************************/
  void          setPathLength(const int length);            //!< Sets a new path length
  
/********************************
 * Getting outputs from state:  *
 ********************************/
  const float   accumulatedDistance() const     {return m_accumulatedDistance;}
  const int     pathLength()          const     {return m_pathLength;}
  const int     pathBufferPtr()       const     {return m_pathBufferPtr;}
  const int     *backTracePath()      const     {return m_backTracePath;}
  const int     bestSymbol()                    {return m_bestInputSymbol;}
  const int     bestSecondSymbol()              {return m_bestSecondSymbol;}
  const int     symbolInPast(int pastLength) const;
  
/********************************
 * Normalizing the distance:    *
 ********************************/
  void          normalizeDistanceBy(const float normalization);

/********************************
 * Updating the state with      *
 * input states:                *
 ********************************/
  void          updateState(TrellisState *inputStates,
                            const float *transitionDistances,
                            const int *transitionSymbols,
                            const int numberInputs);
  void          updateState(TrellisState *inputState,
                            const float transitionDistance,
                            const int transitionSymbol);
  void          updateState(TrellisState *inputState,
                            const float transitionDistance,
                            const int firstSymbol,
                            const int secondSymbol);
  

/********************************
 * Finding the best input path  *
 ********************************/
  void          findBestInput();
};
#endif



