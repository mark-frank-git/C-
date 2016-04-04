#ifndef _PHASE_LOCKED_H
#define _PHASE_LOCKED_H  1
/************************************************************************
 *                                                                      *
 * This class implements a phase-locked loop.  It designs and imple-    *
 * ments different types of phase-locked loops.                         *
 *                                                                      *
 * File:PhaseLocked.h                                                   *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 12/19/97  - Started.                                             *
 *  2. 12/31/97  - Added _tuneFrequency.                                *
 *  3. 01/27/98  - Went to phaseDetectorType, data directed PLL.        *
 *  4. 02/09/98  - Add _unModuloVCOOutput.                              *
 *  5. 11/23/98  - Made this an abstract class, see CarrierLock, Code-  *
 *                 lock subclasses.                                     *
 *  6. 01/14/99  - Add lockedLoopBandwidth, etc.                        *
 *  7. 02/25/99  - Add more lock detector stuff.                        *
 *  8. 03/27/99  - Add _lockDetectorSum[MAX_SUMS]                       *
 *  9. 06/14/99  - Make _lockDetectorSum a scalar.                      *
 * 10. 11/17/99  - Output _loopLocked.                                  *
 * 11. 01/16/01  - Move pllOutputFor() to subclasses, stub here.        *
 ************************************************************************/
#include <GNU/Complex.h>

#ifndef LOGICAL
#define LOGICAL char
#endif
//
// The following are the types of loop filters available
//
#define NO_FILTER               0               // F(s) = 1
#define ONE_POLE                1               // F(s) = 1/(1+sTau)
#define PHASE_LEAD              2               // F(s) = (1+sTau)/s
#define ONE_POLE_PHASE_LEAD     3               // F(s) = (1+sT2)/(1+sT1)
//
// The following are the types of phase detectors available
//
#define COMPLEX_PHASE_DETECTOR  0               // input*exp(-j*phi)
#define REAL_PHASE_DETECTOR     1               // input*cos(phi)
#define DATA_DIRECTED_PD        2               // data directed phase detector
#define EXTERNAL_PD             3               // external PD and VCO

//
// The following are the types of lock detectors available
//
#define NO_LOCK_DETECTOR        0               // don't implement lock detector
#define Q_LOCK_DETECTOR         1               // lock detector against imag output of phase detector
#define I_LOCK_DETECTOR         2               // lock detector against real output of phase detector
#define MAX_DET_TYPE            I_LOCK_DETECTOR

#define MIN_ZETA                0.1             // Minimum damping factor
#define MAX_ORDER               3               // Maximum loop order
#define PEAK_SCALE              0.02            // Also used in lock detect circuit
#define LOCK_BW_SCALE           10.             // Lock BW is 1/10 of normal loop BW
#define DEFAULT_LOCK_THRESH     0.9             // Default lock threshold for PN DLL
#define DEFAULT_LOCK_COUNT      100             // Lock count

class PhaseLocked
{
protected:
  int           _loopOrder;                     // Order of the loop
  int           _loopType;                      // Type 1, etc., # of integrators in GH
  int           _loopFilterType;                // loop filter type, e.g., ONE_POLE
  int           _lockDetectorType;              // type of lock detector
  int           _outputType;                    // PD_OUTPUT, etc.
  int           _numberNumCoeffs;               // # of coefficients in numerator of loop filter
  int           _numberDenCoeffs;               // # of coefficients in denominator of loop filter
  int           _oldNumberOutput;               // To save on mallocs
  int           _numberSlopePoints;             // Used for idle mode
  int           _lockDetectCounter;             // Used for detecting loop in lock
  int           _lockDetectorCount;             // Count to compare lock detector counter against
  int           _phaseDetectorType;             // Complex, real, data directed

  float         _inputAmplitude;                // Input sine wave amplitude for calc'ing phase det gain
  float         _samplingFrequency;             // Sampling frequency in Hertz
  float         _dampingFactor;                 // Loop damping factor, zeta
  float         _loopBandwidth;                 // fn in Hertz
  float         _lockedLoopBandwidth;           // fn in Hertz after lock is achieved
  float         _tuneFrequency;                 // initial offset frequency of VCO in Hertz
  float         _vcoOffset;                     // derived from _tuneFrequency
  float         _lockThreshold;                 // Lock detect threshold
  float         _lockDetectorSum;               // Lock detector variable to be checked against thresh.

  double        _loopFilterOutput;              // Output of loop filter
  double        _vcoOutput;                     // VCO accumulator output
  double        _unModuloVCOOutput;             // Not modulo by TWOPI
  double        _vcoSlope, _oldVCOOutput;       // Needed for open loop operation
  double        _peakDetect;                    // Also needed for lock detect/ open loop operation
  double        _integratedFilterOutput;        // Integrated filter output for phase jitter analysis
  double        _phaseRotationOutput;           // phase rotation estimation from DLL
  double        *_numCoeffs;                    // Loop filter numerator coefficients
  double        *_denCoeffs;                    // Loop filter denominator coefficients
  double        *_shiftArray;                   // Shift register values for loop filter

  Complex       *_outputArray;                  // Array of output values
  Complex       _pdOutput;                      // Output of phase detector
  Complex       _oscillatorOutput;              // Output of VCO
  Complex       _pllOutput;                     // Selected output of PLL

  LOGICAL       _openLoopOperation;             // Free running loop
  LOGICAL       _loopLocked;                    // YES = loop lock was detected
//
// Private functions:
//
  inline        double  loopFilter(double input);               // loop filter next input
  virtual       inline  LOGICAL loopLocked(Complex &input);     // Returns YES if loop is locked
  virtual       inline  void findPDOutput(Complex &input);      // Processes input through phase detector
  virtual       float   findPDVCOGain(float inputAmplitude);    // Find the gain throught the phase det
                                                                // VCO combination
  virtual       float   findSamplingTime();                     // Find the loop filter sampling time
  virtual       void    processPDOutput();                      // Process the output of phase detector

  void          initializeLoopFilterFor(int type, float pdGain, LOGICAL initializeTaps=1);

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  PhaseLocked(int loopFilter=PHASE_LEAD, int output=0, float input=1., float zeta=0.707, float fn=1000.,
              float sampling=1000., float tune=0.);
  virtual ~PhaseLocked();
        
/**********************
 * Set parameters:    *
 **********************/
  void          setOpenLoopOperation(LOGICAL flag)      {_openLoopOperation     = flag; return;}
  void          setPhaseDetectorType(int type)          {_phaseDetectorType     = type; return;}
  void          setLoopFilterType(int type);
  void          setOutputType(int type);
  void          setLockDetectorType(int type);
  void          setInputAmplitude(float amplitude);
  void          setDampingFactor(float zeta);
  void          setLoopBandwidth(float fn);
  void          setLockedLoopBandwidth(float fn);
  void          setLockThreshold(float thresh);
  void          setLockDetectorCount(int count);
  virtual       void    setSamplingFrequency(float sampling);
  virtual       void    setTuneFrequency(float tune);

/**********************
 * Get parameters:    *
 **********************/
  int           loopOrder()                             {return _loopOrder;}
  int           loopType()                              {return _loopType;}
  int           loopFilterType()                        {return _loopFilterType;}
  int           outputType()                            {return _outputType;}
  LOGICAL       loopIsLocked()                          {return _loopLocked;}
  float         inputAmplitude()                        {return _inputAmplitude;}
  float         samplingFrequency()                     {return _samplingFrequency;}
  float         dampingFactor()                         {return _dampingFactor;}
  float         loopBandwidth()                         {return _loopBandwidth;}
  float         tuneFrequency()                         {return _tuneFrequency;}

/************************
 * Initializing the     *
 * PLL.                 *
 ************************/
  virtual       void    initializePLL();
  
/**********************
 * Getting Outputs:   *
 **********************/
  virtual Complex       pllOutputFor(Complex &input);                   // Returns a single output
  virtual Complex       *pllOutputFor(Complex *input, int number);      // Returns an array of outputs
  virtual Complex       *pllOutputFor(const float *realInput, const float *imagInput, int number);

};

#endif
