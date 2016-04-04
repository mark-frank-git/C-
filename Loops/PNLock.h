#ifndef _PN_LOCK_H
#define _PN_LOCK_H  1
/************************************************************************
 *                                                                      *
 * This subclass of PhaseLocked supports locking to a PN code.          *
 *                                                                      *
 * File:PNLock.h                                                        *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 11/23/98  - Subclassed from PhaseLocked.                         *
 *  2. 03/17/99  - Removed non coherent stuff.  Moved it to a subclass. *
 *  3. 03/27/99  - Add _samplesPerChip, _sampleCount, etc.              *
 *  4. 06/15/99  - Add _timingVariance.                                 *
 *  5. 01/16/01  - Implement pllOutput from super class, add            *
 ************************************************************************/
#include "PhaseLocked.h"
#include <stdio.h>

class   PNGenerator;                                    // Class prototypes
class   PNFiltered;
class   IS95PNGen;

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
#define PN_EARLY_OUTPUT         9               // _iPNEarly, _qPNEarly
#define PN_LATE_OUTPUT          10              // _iPNLate, _qPNLate
#define PN_ON_OUTPUT            11              // _iPNOnTime, _qPNOnTime
#define DESPREAD_OUTPUT         12              // _iPNOnTime*inputI, _qPNOnTime*inputQ
#define MAX_TYPE                DESPREAD_OUTPUT

#define IS95_I_PN               121641          // PN polynomial in octal
#define IS95_Q_PN               116171

#ifndef LOGICAL
#define LOGICAL char
#endif

#define DEFAULT_SAMPLING        4.9152e6                // Sampling rate in Hertz
#define DEFAULT_TUNE            (DEFAULT_SAMPLING/4.)   // 1.2288 MHz
#define DEFAULT_INTEGRATE       10                      // I/D default size 
#define DEFAULT_OFFSET          0                       // PN offset
#define MAX_OFFSET              33000                   // size 2**15
#define WALSH_CHIPS             64                      // # of chips in Walsh sequence
#define NUMBER_QUADS            4                       // Number of quadrants

class PNLock: public PhaseLocked
{
protected:
  PNFiltered    *_iPN;                          // I channel PN generator
  PNFiltered    *_qPN;

  int           _integrateSize;                 // Integrate and dump samples
  int           _integrateSamples;              // Count of integration samples
  int           _pnOffset;                      // PN offset in samples
  int           _numberChips;                   // This could overflow for large sequences!!!
  int           _samplesPerChip;                // # of samples per each chip
  int           _numberTimingAvgs;              // Counter for statistics
  int           _quadrantCount[NUMBER_QUADS];   // For rotation estimation
  
  float         _pdOutputReal;                  // Phase detector output I channel
  float         _pdOutputImag;                  // Q channel
  float         _iPNEarly;                      // Local I PN generator outputs
  float         _iPNLate;
  float         _iPNOnTime;
  float         _qPNEarly;                      // Local Q PN generator outputs
  float         _qPNLate;
  float         _qPNOnTime;
//
// The following only need to be instance variables if ROTATION_INSENSITIVE is used
//
  float         _iEarlyOut;                     // Early arm I channel output
  float         _qEarlyOut;                     // Early arm Q channel output
  float         _iLateOut;                      // Late arm I channel output
  float         _qLateOut;                      // Late arm Q channel output
  float         _iOnTimeOut;                    // Needed for rotation estimation
  float         _qOnTimeOut;                    // Needed for rotation estimation
  double        _timingAverage;                 // Average of integrated filter output
  double        _timingVariance;                // Variance of integrated filter output
  double        _qOverIAverage[NUMBER_QUADS];   // Average of Q/I rotation output
  LOGICAL       _rotationInsensitive;           // YES = loop insensitive to rotations in input constell.
//
// Private functions:
//
  virtual inline LOGICAL        loopLocked(Complex &input);     // Returns YES if loop is locked
  virtual inline void   findPDOutput(Complex &input);           // Processes input through phase detector
  float                 findSamplingTime();                     // Returns filter sampling time in secs
  virtual float         findPDVCOGain(float inputAmplitude);    // Find the gain throught the phase det/
                                                                // VCO combination
  inline void           processPDOutput();                      // process PD output through loop

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  PNLock(int loopFilter=PHASE_LEAD, int output=PD_OUTPUT, float input=1., float zeta=0.707, float fn=1000.,
         float sampling=DEFAULT_SAMPLING, float tune=DEFAULT_TUNE, int pnPolyI=IS95_I_PN,
         IS95PNGen *gen=NULL);
  virtual ~PNLock();
        
/**********************
 * Set parameters:    *
 **********************/
  void  setIPoly(int poly);                             // PN generator polynomials
  void  setQPoly(int poly);
  void  setIntegrateSize(int size);
  void  setRotationInsensitive(LOGICAL flag);           // Rotation insensitive or not
  void  setPNOffset(int offset);                        // pre-load shift register
  void  setSamplingAndTune(float sampling, float tune); // Set new sampling and tune frequencies

/********************************
 * The following methods        *
 * return calculate values:     *
 ********************************/
  float timingEstimate();
  float timingVariance();
  float rotationEstimate();


/********************************
 * The following method         *
 * is used to reset the DLL     *
 * before running a new set of  *
 * samples.                     *
 ********************************/
  virtual void  initializePLL();
  
/**********************
 * Getting Outputs:   *
 **********************/
  Complex       pllOutputFor(Complex &input);                   // Returns a single output
  Complex       *pllOutputFor(Complex *input, int number);      // Returns an array of outputs
  Complex       *pllOutputFor(const float *realInput, const float *imagInput, int number);

};

#endif
