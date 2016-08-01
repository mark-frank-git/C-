#ifndef _CARRIER_LOCK_H
#define _CARRIER_LOCK_H  1
/************************************************************************
 *                                                                      *
 * This subclass of PhaseLocked supports locking to a carrier.          *
 *                                                                      *
 * File:CarrierLock.h                                                   *
 *                                                                      *
 ************************************************************************/
#include "PhaseLocked.h"

//
// The following are the types of outputs that are available:
//
#define PD_OUTPUT               0               // output phase detector
#define LOOP_FILTER_OUTPUT      1               // output loop filter output
#define INTEGRATED_LOOP_FILTER  2               // output integrated loop filter output
#define VCO_OUTPUT              3               // output VCO output
#define OSCILLATOR_OUTPUT       4               // Input to phase detector
#define DEMOD_OUTPUT            5               // Equal to zn for data directed
#define VCO_SLOPE               6               // VCO slope for idle operation
#define LOCK_DETECTOR           7               // Lock detector variable
#define LOCK_STATUS             8
#define MAX_TYPE                LOCK_STATUS

#define LOCK_DETECT_COUNT       1000            // # of counts after lock detected for estimating carrier freq

class CarrierLock: public PhaseLocked
{
private:
  Complex       _zn;                            // Data directed demod output
  int           _loopLockCounter;               // Counter for after loop locked
  LOGICAL       _loopRelocked;                  // Allows a carrier freq estimate to be made after 1st lock
//
// Private functions:
//
  inline        LOGICAL loopLocked(Complex &input);             // Returns YES if loop is locked
  inline        void    findPDOutput(Complex &input);           // Processes input through phase detector
  inline        void    findPDOutput(Complex &input, Complex &sliceData);
                                                                // Processes input through phase detector
  inline        void    processPDOutput();                      // process PD output through loop

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  CarrierLock(int loopFilter=PHASE_LEAD, int output=PD_OUTPUT, float input=1., float zeta=0.707, float fn=1000.,
              float sampling=1000., float tune=0.);
  ~CarrierLock();
        
/**********************
 * Set parameters:    *
 **********************/

/**********************
 * Get parameters:    *
 **********************/

/************************
 * Initializing the     *
 * PLL.                 *
 ************************/
  
/**********************
 * Getting Outputs:   *
 **********************/
  Complex       pllOutputFor(Complex &input);                   // Returns a single output
  Complex       *pllOutputFor(Complex *input, int number);      // Returns an array of outputs
  Complex       *pllOutputFor(const float *realInput, const float *imagInput, int number);
  Complex       pllOutputForSlice(Complex &input, Complex &sliceData);

};

#endif
