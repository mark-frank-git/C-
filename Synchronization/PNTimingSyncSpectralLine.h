#ifndef _SPECTRAL_LINE_H
#define _SPECTRAL_LINE  1
/************************************************************************
 *                                                                      *
 *                                                                      *
 * This subclass of object estimates sampling time in synchronization   *
 * with PN chip timing using spectral line method                       *
 *                                                                      *
 * File:SpectralLine.h                                                  *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 10/06/99  - Started.                                             *
 ************************************************************************/
#define SPLN_OUTPUT_SIZE        2000
#define POST_NL                 0
#define TIME_TONE               1
#define TIME_UP_ZC              2
#define TIME_DOWN_ZC            3

#ifndef LOGICAL
#define LOGICAL char
#endif

class   TestData;

class PNTimingSyncSpectralLine
{
protected:

  int           _samplesPerChip;                // # of samples per each chip of the input signal
  int           _upSampleFactorOnInputSignal;   // required up sample ratio for the input signal
  int           _totalUpSampleRatio;            // total upsample ratio of the signal in this processing

  int           _inputDataLength; 
  int           _dataLength;                    // to keep trace of the valid data length
  int           _oldDataLength;                 // To keep track of news
  int           _numChips;                      // Number of chips for zero crossing calculation
  int           _oldNumChips;                   // For allocating _timeUpZC, _timeDownZC
  int           _preBPFLength;
  int           _postBPFLength;
  int           _outputSize;                    // Size of output array for plotting

  float         *_iUpSampledData;
  float         *_qUpSampledData;
  float         *_iPostPreBPFData;
  float         *_qPostPreBPFData;
  float         *_postNLData;
  float         *_postBPFData;
  float         *_timeUpZC;
  float         *_timeDownZC;
 
  float         _timingPhaseInUpSampledSpace;

  class AbstractFIR     *preBPF;
  class AbstractFIR     *postBPF;
  TestData      *_fileReader;                   // Reads coefficients from file
//
// Private function:
//
  LOGICAL       writeArrayToFile(const char *fileName, float *iArray, float *qArray, int size);

// 
// Public functions:
//
public:
/********************************
 * Constructors, destructors    *
 ********************************/
  PNTimingSyncSpectralLine(int upSampleRatioIn,int dataLen);
  ~PNTimingSyncSpectralLine();

/********************************
 * functions to calculate       *
 * the timing offset value      *
 ********************************/
  void zeroPaddingUpSampling(float *iData, float *qData);
  void preBPFiltering();
  void absNonlinearProcess();
  void postBPFiltering();
  float zeroCrossTimingPhaseEstimate(float *tmTone=NULL, int len=2000, int upSampleRatioTone=8);
  float getPNtimingOffsetWithNLbiasAdjusted(float ZCtimingPhase);

//*********************************
//* functions to set parameters   *
//*********************************
 void           setInputDataLength(int Len);
 void           setSamplesPerChip(int samples);
 LOGICAL        readPreBPFCoeffsFromFile(const char *fileName, int fileType);
 LOGICAL        readPostBPFCoeffsFromFile(const char *fileName, int fileType);

/********************************
 * Functions to get parameters: *
 ********************************/
 int    totalUpSampleRatio()    {return _totalUpSampleRatio;}   // Determines coeffs of BPF

/********************************
 * Functions to get outputs.    *
 ********************************/
  const float *getOutputData(int type);
  int   getOutputSize()         {return _outputSize;}

};

#endif
