/************************************************************************
 *                                                                      *
 * This subclass contains functions of estimating sampling time in      *
 * synchronization with PN chip timing using spectral line method       *
 *                                                                      *
 * File:PNTimingSyncSpectralLine.cc                                     *
 *                                                                      *
 ************************************************************************/
#include <stdio.h>
#include <math.h>

#include <C_Libraries/constants.h>
#include <Filters/AbstractFIR.h>
#include <DataIO/TestData.h>

#include "PNTimingSyncSpectralLine.h"                                   // Object prototypes


#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)        ( ((a)>(b)) ? (a) : (b) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )
#define ABS(a)          ( (a) >= 0 ? (a) : (-a) )


// ######################## Private Function ################################
// writeArrayToFile -- Upsample the input baseband signal to 8~9 
//                          samples per chip by using zero padding method
// Input:       fileName:       name of file
//              iArray:         real data array
//              qArray:         imaginary data array 
//              size:           size of array
//              
//
// Output:                      YES if OK
//
// Notes: 
// 1. program starts on 10/12/99
//
//###########################################################################
LOGICAL PNTimingSyncSpectralLine::writeArrayToFile(const char *fileName, float *iArray,float *qArray, int size)
{
  int   rtn_val;

  if(_fileReader == NULL)
    _fileReader = new TestData();

  rtn_val       = _fileReader->writeTekComplexFile(fileName, iArray, qArray, size);
  return rtn_val;
}


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the SpectralLine class.
//
// Input:               upSampleRatioIn         samples per chip of the input signal
//                      dataLen                 number of input data
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
 PNTimingSyncSpectralLine::PNTimingSyncSpectralLine(int upSampleRatioIn,int dataLen)
{
   _oldDataLength               = _dataLength   = 0;
   _oldNumChips                 = 0;
   _iUpSampledData              = NULL;
   _qUpSampledData              = NULL;
   _timeUpZC                    = NULL;
   _timeDownZC                  = NULL;
   setSamplesPerChip(upSampleRatioIn);
   setInputDataLength(dataLen);

  preBPF                        = new AbstractFIR();
  postBPF                       = new AbstractFIR();

  _fileReader                   = NULL;

 return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the SpectralLine class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
PNTimingSyncSpectralLine::~PNTimingSyncSpectralLine()
{
  delete        preBPF;
  delete        postBPF;
  delete        _fileReader;
  delete[]      _iUpSampledData;
  delete[]      _qUpSampledData;
  delete[]      _timeUpZC;
  delete[]      _timeDownZC;
  return;
}

// ######################## Public Function ################################
// zeroPaddingUpSampling -- Upsample the input baseband signal to 8~9 
//                          samples per chip by using zero padding method
// Input:       
//              I and Q baseband data
//              
//
// Output:      none
//
// Notes: 
// 1. program starts on 10/12/99
//
//###########################################################################
void PNTimingSyncSpectralLine::zeroPaddingUpSampling(float *iData, float *qData)
{
  int n, m, i;
writeArrayToFile("d:/OriginalSamples.bin",iData,qData,_inputDataLength);
  i     = 0;
  for (n=0; n<_inputDataLength; n++)
  {
    _iUpSampledData[i]          = iData[n];
    _qUpSampledData[i]          = qData[n];
    i++;
    for (m=1; m<_upSampleFactorOnInputSignal; m++)
    {
      _iUpSampledData[i]                = 0;
      _qUpSampledData[i]                = 0;
      i++;
     }
   }
}

// ######################## Public Function ################################
// preBPFiltering -- perform bandpass filtering on the input signal using 
//                   a specially designed BPF which is centered at half of 
//                   the chip rate
//
// Input:       upsampled I and Q baseband data (implicitly)
//              pre_BPF coefficients(implicitly)
//
// Output:      bandpassed complex data(implicitly)
//
// Notes: 
// 1. program starts on 10/12/99
//
//###########################################################################
void PNTimingSyncSpectralLine::preBPFiltering()
{
  int n;
  

  _preBPFLength         = preBPF->numberTaps();
    
  preBPF->zeroOutTaps();
  preBPF->filterFloatArray(_iUpSampledData, _dataLength);
  preBPF->zeroOutTaps();
  preBPF->filterFloatArray(_qUpSampledData, _dataLength);

  n                     = _preBPFLength-1;
  _iPostPreBPFData      = _iUpSampledData + n;
  _qPostPreBPFData      = _qUpSampledData + n;
  _dataLength           -= n;

}

// ######################## Public Function #################################
// nonlinearProcess -- perform nonlinear operation on the upsampled signal 
//
// Input:       upsampled I and Q baseband data (implicitly)
//
// Output:      
//
// Notes: 
// 1. program starts on 10/13/99
//
//###########################################################################
void PNTimingSyncSpectralLine::absNonlinearProcess()
{
  int i;

  _postNLData           = _iUpSampledData;
  for (i=0;i<_dataLength;i++)
     _postNLData[i]=_iPostPreBPFData[i]*_iPostPreBPFData[i]+_qPostPreBPFData[i]*_qPostPreBPFData[i]; 
}


// ######################## Public Function ################################
// postBPFiltering -- perform bandpass filtering on the nonlinearly distorted  
//                    signal using a BPF centered at chip frequency
//                  
// Input:       the output of nonlinear process(implicitly)
//              BPF coefficients (implicitly)
//
// Output:      timing tone data(implicitly)
//
// Notes: 
// 1. program starts on 10/13/99
//
//###########################################################################
void PNTimingSyncSpectralLine::postBPFiltering()
{
  int n;

  _postBPFLength                = postBPF->numberTaps();

  postBPF->zeroOutTaps();
  postBPF->filterFloatArray(_postNLData, _dataLength);

  n                     = _postBPFLength-1;
  _postBPFData          = _postNLData + n;
  _dataLength           -= n;
 
}

// ######################## Public Function ################################
// zeroCrossTimingPhaseEstimate -- estimate timing-tone's phase using
//                                 zero-crossing method
//                  
// Input:       timing-tone
//              length of the timing tone
//              the up sample ratio of the timing tone (w.r.t chip rate)
//
// Output:      timing-tone's phase in upsampled space
//
// Notes: 
// 1. program starts on 10/13/99
//
//###########################################################################
float PNTimingSyncSpectralLine::zeroCrossTimingPhaseEstimate(float *tmTone=NULL, int len, int upSampleRatioTone)
{
 int    countUpZC, countDownZC, countZC, i, ku, kd;
 float  avgUpZCPhase, avgDownZCPhase, meanUp, meanDn, er, WrapTH, VarTH;

 // get source and parameters
 if (tmTone==NULL)
 {
   tmTone               = _postBPFData;
   len                  = _dataLength;
   upSampleRatioTone    = _totalUpSampleRatio;
  }
//
// create buffers for recording zero crossing time:
//
  _numChips             = (int) (len/upSampleRatioTone*(1+0.25));
  if(_oldNumChips < _numChips)
  {
    _oldNumChips        = _numChips;
    _timeUpZC           = new float[_numChips]; 
    _timeDownZC         = new float[_numChips];
    for(i=0; i<_numChips; i++)
      _timeUpZC[i]      = _timeDownZC[i]        = 0.;           // Only needed for plotting output
  }

// detect zeros crossings
  countUpZC             = 0;
  countDownZC           = 0;
  for (i=0; i<len;i++)
  { 
   if (tmTone[i]<=0 && tmTone[i+1]>=0)
     _timeUpZC[countUpZC++]             = i-tmTone[i]/(tmTone[i+1]-tmTone[i]);
   else if (tmTone[i]>=0 && tmTone[i+1]<=0)
     _timeDownZC[countDownZC++]         = i+tmTone[i]/(tmTone[i]-tmTone[i+1]);
  }
  
// remove periods to get initial time phase
  countZC               = MIN(countUpZC, countDownZC);
  for (i=0; i<countZC; i++)
  {
    ku                  = (int) (_timeUpZC[i]/upSampleRatioTone);
    _timeUpZC[i]        -= ku*upSampleRatioTone;
    kd                  = (int) (_timeDownZC[i]/upSampleRatioTone);
    _timeDownZC[i]              -= kd*upSampleRatioTone;
  }
// detect phase wrapping and corrected
  meanUp                = _timeUpZC[0];
  meanDn                = _timeDownZC[0];
  WrapTH                = upSampleRatioTone/2.0;
  for (i=0; i<countZC; i++)
  {
    er                  = _timeUpZC[i]-meanUp;
    if (er>WrapTH)
      _timeUpZC[i]      -= upSampleRatioTone;
    else if (er<-WrapTH)
       _timeUpZC[i]     += upSampleRatioTone;
    meanUp              = (meanUp*i+_timeUpZC[i])/(i+1);
    
    er                  = _timeDownZC[i]-meanDn;
    if (er>WrapTH)
      _timeDownZC[i]    -= upSampleRatioTone;
    else if (er<-WrapTH)
       _timeDownZC[i]   += upSampleRatioTone;
    meanDn              = (meanDn*i+_timeDownZC[i])/(i+1);
   }

 // remove largely deviated points and measure the average time phases
  avgUpZCPhase          = 0;
  avgDownZCPhase        = 0;
  countUpZC             = 0;
  countDownZC           = 0;
  VarTH                 = upSampleRatioTone/10.0;
  for (i=0; i<countZC; i++)
  {
    er                  = _timeUpZC[i]-meanUp;
    if (ABS(er)<VarTH)
    {
      avgUpZCPhase      += _timeUpZC[i];
      countUpZC++;
     }
    
    er                  = _timeDownZC[i]-meanDn;
    if (ABS(er)<VarTH)
    {
      avgDownZCPhase    += _timeDownZC[i];
      countDownZC++;
     }
   }

  avgUpZCPhase          /= countUpZC;
  avgDownZCPhase        /= countDownZC;

 // measure timing phase in terms of delay
  ku                    = (int) (avgUpZCPhase/upSampleRatioTone);
  avgUpZCPhase          -= ku*upSampleRatioTone;
  kd                    = (int) (avgDownZCPhase/upSampleRatioTone);
  avgDownZCPhase        -= kd*upSampleRatioTone;

// avgerage upward and downward measures
 if (avgUpZCPhase > avgDownZCPhase)
   _timingPhaseInUpSampledSpace         = (avgUpZCPhase - upSampleRatioTone/2.0 + avgDownZCPhase)/2.0;
 else
   _timingPhaseInUpSampledSpace         = (avgUpZCPhase + upSampleRatioTone/2.0 + avgDownZCPhase)/2.0;

 printf("timing phase in Tone =%f,  avgU=%f, avgD=%f \n",_timingPhaseInUpSampledSpace, avgUpZCPhase,avgDownZCPhase);

 return  _timingPhaseInUpSampledSpace;
}

// ######################## Public Functions #####################################
// getPNtimingOffsetWithNLbiasAdjusted: adjust phase shift by nonlinear process,
//                                      and calulate PN timing offset based on 
//                                      original sample unit
//
// Input:               ZCtimingPhase:  timing phase of the zero-crossing tone
//
// Output:              PNtimingOffset: PN timing offset in fractions of a chip
//
// Notes:
// ############################# Class Constructor ###############################
#define PHASE_SHIFT_BY_AbsNL            0.25
float PNTimingSyncSpectralLine::getPNtimingOffsetWithNLbiasAdjusted(float ZCtimingPhase)
{
  float tau;
  
  //adjust phase shift by absolute-nonlinearity
  tau           = (ZCtimingPhase - _totalUpSampleRatio*PHASE_SHIFT_BY_AbsNL)/_totalUpSampleRatio;

  //relate timing position in delay dimension
  tau           -= floor((double)tau);

  //the amount of delay needed for adjustment
  tau           = 1-tau;

  //measure PN timing offset in terms of original sample unit
  //tau         = tau*_samplesPerChip;
  
  return tau;
}

// ######################## Public Functions #####################################
// setInputDataLength: Sets a new input data length, and allocates arrays
//                      appropriately
//
// Input:       len:    new data length
//
// Output:              none
//
// Notes:
// ############################# Class Constructor ###############################
 void PNTimingSyncSpectralLine::setInputDataLength(int Len)
 {
   _inputDataLength     = Len;
   _dataLength          = _upSampleFactorOnInputSignal*_inputDataLength;
   if(_dataLength > _oldDataLength)
   {
     _oldDataLength     = _dataLength;
     delete     []_iUpSampledData;
     delete     []_qUpSampledData;
     _iUpSampledData    = new float[_dataLength];
     _qUpSampledData    = new float[_dataLength];
   }
   return;
 }


// ######################## Public Functions #####################################
// setSamplesPerChip: Sets a new number of samples per chip.
//
// Input:       upSampleRatioIn:        new samples/chip
//
// Output:                              none
//
// Notes:
// 1. This should be called 
// ############################# Class Constructor ###############################
void PNTimingSyncSpectralLine::setSamplesPerChip(int upSampleRatioIn)
{
  _samplesPerChip               = upSampleRatioIn;
  _upSampleFactorOnInputSignal  = ROUND(8.0/_samplesPerChip);
  _totalUpSampleRatio           = _upSampleFactorOnInputSignal*upSampleRatioIn;
  setInputDataLength(_dataLength);
  return;
}

// ######################## Public Functions #####################################
// readPreBPFCoeffsFromFile: Reads in the pre-BPF filter coefficients from file.
//
// Input:       fileName:               name of filter file
//              type:                   file type 
//
// Output:                              YES if file read succeeded
//
// Notes:
// ############################# Class Constructor ###############################
LOGICAL PNTimingSyncSpectralLine::readPreBPFCoeffsFromFile(const char *fileName, int fileType)
{
  return preBPF->readBCoeffsFromFile(fileName, fileType);
}


// ######################## Public Functions #####################################
// readPostBPFCoeffsFromFile: Reads in the post-BPF filter coefficients from file.
//
// Input:       fileName:               name of filter file
//              type:                   file type 
//
// Output:                              YES if file read succeeded
//
// Notes:
// ############################# Class Constructor ###############################
LOGICAL PNTimingSyncSpectralLine::readPostBPFCoeffsFromFile(const char *fileName, int fileType)
{
  return postBPF->readBCoeffsFromFile(fileName, fileType);
}

// ######################## Public Functions #####################################
// getOutputData: Get output data for display/debug
//
// Input:       type:   Type of data output
//
// Output:              Output data array
//
// Notes:
// ############################# Class Constructor ###############################
const float* PNTimingSyncSpectralLine::getOutputData(int type)
{
   float        *output = NULL;

   switch(type)
   {
     default:
     case POST_NL:
       output           = _postNLData;
       _outputSize      = SPLN_OUTPUT_SIZE;
       break;
     case TIME_TONE:
       output           = _postBPFData;
       _outputSize      = SPLN_OUTPUT_SIZE;
       break;
     case TIME_UP_ZC:
       output           = _timeUpZC;
       _outputSize      = _numChips;
       break;
     case TIME_DOWN_ZC:
       output           = _timeDownZC;
       _outputSize      = _numChips;
       break;
    }
   
   return output;
 }
