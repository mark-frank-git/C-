#ifndef _PN_LOCK_NON_COH_H
#define _PN_LOCK_NON_COH_H  1
/************************************************************************
 *                                                                      *
 * This subclass of PNLock supports locking to a PN code, assuming      *
 * non-coherent operation (i.e., carrier frequency offset).             *
 *                                                                      *
 * File:PNLockNonCoh.h                                                  *
 *                                                                      *
 ************************************************************************/
#include "PNLock.h"
class   IIRFilter;

#define FILTER_ORDER            4                       // BPF order
#define BW_FACTOR               2.                      // BPF factor


class PNLockNonCoh: public PNLock
{
private:
  IIRFilter     *_earlyFilter;                  // early channel BPF
  IIRFilter     *_lateFilter;                   // late channel BPF
  IIRFilter     *_onTimeFilter;                 // on time filter for lock detector

  float         _bpfFrequency;                  // BPF center frequency for non-coherent loop.  Should
                                                // be set to carrier frequency in Hz
  
//
// Private functions:
//
  inline        LOGICAL loopLocked(Complex &input);             // Returns YES if loop is locked
  inline        void    findPDOutput(Complex &input);           // Processes input through phase detector
  float         findPDVCOGain(float inputAmplitude);            // Find the gain throught the phase det/
                                                                // VCO combination
// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  PNLockNonCoh(int loopFilter=PHASE_LEAD, int output=PD_OUTPUT, float input=1., float zeta=0.707, float fn=1000.,
              float sampling=DEFAULT_SAMPLING, float tune=DEFAULT_TUNE, int pnPolyI=IS95_I_PN,
               IS95PNGen *gen=NULL);
  ~PNLockNonCoh();
        
/**********************
 * Set parameters:    *
 **********************/
  void  setBPFFrequency(float fo);
  void  setSamplingAndTune(float sampling, float tune);

/**********************
 * Get parameters:    *
 **********************/


/************************
 * Initializing the     *
 * PLL.                 *
 ************************/
  void  initializePLL();
  
/**********************
 * Getting Outputs:   *
 **********************/

};

#endif
