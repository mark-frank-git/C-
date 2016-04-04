#ifndef _BINARY_RATIONAL_NCO_H
#define _BINARY_RATIONAL_NCO_H  1
/********************************************************************************
 * This class implements a binary-rational NCO.                                 *
 * It consists of two accumulators with the output accumulator taking           *
 * the carry out of the input accumulator.                                      *
 *                                                                              *
 * File: /User/frank/C++/Generators/BinaryRationalNCO.h                         *
 *                                                                              *
 * Revision History:                                                            *
 *  1. 01/18/01 - Started.                                                      *
 ********************************************************************************/
#define DEFAULT_PHASES          128             // # of phases in output accumulator
#define DEFAULT_R2              8               // output accumulator increment
#define DEFAULT_R0              23              // first accumulator increment
#define DEFAULT_R1              102             // first accumulator alternate increment

#ifndef LOGICAL
#define LOGICAL char
#endif

class BinaryRationalNCO
{
protected:
  int           _numberPhases;                  // Number of output phases
  int           _r0;                            // input accumulator increment
  int           _r1;                            // input accumulator increment on overflow
  int           _r2;                            // output accumulator increment
  int           _inputAccumulator;              // Input accumulation
  int           _outputAccumulator;             // Output accumulation
  int           _inputOverflow;                 // Modulus of input accumulator

  LOGICAL       _inputCarryOut;                 // YES = first accumulator had a carry out
//
// Private methods
//
public:

/*********************************
 * Constructors/destructors:     *
 *********************************/
  BinaryRationalNCO(int phases=DEFAULT_PHASES, int r0=DEFAULT_R0, int r1=DEFAULT_R1,
                    int r2=DEFAULT_R2);         // Constructor
  virtual ~BinaryRationalNCO();

/*******************************
 * These methods initialize:    *
 *******************************/
  void  initNCO();

/*******************************
 * These methods set parameters:*
 *******************************/
  void  setNumberPhases(int phases);
  void  setR0(int r0);
  void  setR1(int r1);
  void  setR2(int r2);

/*******************************
 * These methods get parameters*
 *******************************/
  int           numberPhases()                  {return _numberPhases;}
  int           r0()                            {return _r0;}
  int           r1()                            {return _r1;}
  int           r2()                            {return _r2;}

/********************************
 * These methods generate       *
 * new clock data.              *
 ********************************/
  int   newPhase();                             // Gets the next NCO phase, this function must
                                                // be called every sample time
  double newRadianPhase();                      // Alternative to the above, returns phase in radians


};

#endif
