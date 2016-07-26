/************************************************************************
 *                                                                      *
 * This class implements analyzes a set of data using quadrature mirror *
 * filter banks.                                                        *
 *                                                                      *
 * File:WaveletAnalyzer.h                                               *
 *                                                                      *
 *                                                                      *
 *                                                                      *
 * References:                                                          *
 *  1. P.E. Pace, Low Probability of Intercept Radar.                   *
 *                                                                      *
 ************************************************************************/
#include <stdio.h>
#include "WaveletAnalyzer.h"                                    // Object prototypes
#include "RealFIR.h"
#include "WaveletFilterCoefficients.h"
#if defined(WIN32)
#include <GNU/Complex.h>
#else
#include "Complex.h"
#endif

#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ( ((a)>(b)) ? (b) : (a) )
#define MAX(a,b)  ( ((a)>(b)) ? (a) : (b) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

const int DECIMATION_SIZE       = 2;

// ############################# Private Method ###############################
// filterData -- Filters and decimates by two the input data through the given filter.
// Input:       inputData:      data to be filtered
//              length:         number of data points
//              filter:         filtering object
//          
// Output:      outputData:     filtered data
// Notes:
// ############################# Private Method ###############################
void WaveletAnalyzer::filterData(const Complex *inputData, Complex *outputData, int length, RealFIR *filter)
{
//
// Initialize filter:
//
  filter->reset();
  int filter_length     = filter->numberTaps();
  int half_filter       = filter_length/2;
//
// If data length < half filter, then we have to send in the data, send
// in zeros, and continue sending in zeros to get output.
// Otherwise, we just send in data, and flush at the end.
//
  int fill_size         = MIN(half_filter, length);
//
// Fill up half the filter:
//
  Complex zero_data     = Complex(0., 0.);
  Complex filter_output;
  int i;
  for(i=0; i<fill_size; i++)
  {
    filter->filterComplexData(inputData[i]);
  }
  for(;i<half_filter; i++)                      // This is only occurs for length < half_filter
  {
    filter->filterComplexData(zero_data);
  }
  int output_index      = 0;
  for(; i<length; i++)                          // Now continue filtering while extracting the output
  {
    filter_output       = filter->filterComplexData(inputData[i]);
    if(output_index%DECIMATION_SIZE==0)         // decimation
    {
      outputData[output_index/DECIMATION_SIZE]  = filter_output;
    }
    output_index++;
  }
  for(i=0; i<fill_size; i++)
  {
    filter_output       = filter->filterComplexData(inputData[i]);
    if(output_index%DECIMATION_SIZE==0)         // decimation
    {
      outputData[output_index/DECIMATION_SIZE]  = filter_output;
    }
    output_index++;
  }

  return;
}

// ############################# Private Method ###############################
// updateFilterPointer -- Updates the filter and data pointers.
// Input:       currentFilter:  current filter 'H' or 'G', modified on output
//              nextFilter:     next filter 'H' or 'G', modified on output
//              dataSize:       size of data for incrementing pointers
//          
// Output:                      0, or data size depending on whether pointers
//                              need to be incremented
// Notes:
// 1. Refer to Figure 10.12 in Pace' book, and note that we start from the
//    bottom of the figure.
// ############################# Private Method ###############################
int WaveletAnalyzer::updateFilterPointer(char *currentFilter, char *nextFilter, int dataSize)
{
//
// Update pointers based on our current state:
//
  int   increment_size  = 0;
  if( *currentFilter == 'H' )
  {
    if( *nextFilter == 'G' )                    // HG
    {
      *currentFilter    = 'G';
      *nextFilter       = 'G';
    }
    else                                        // HH
    {
      *nextFilter       = 'G';
      increment_size    = dataSize;             // move to next block of data
    }
  }
  else                                          // current == G
  {
    if( *nextFilter == 'G' )                    // GG
    {
      *nextFilter       = 'H';
      increment_size    = dataSize;             // move to next block of data
    }
    else                                        // GH
    {
      *currentFilter    = 'H';
      *nextFilter       = 'H';
    }
  }
  
  return increment_size;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the WaveletAnalyzer class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
WaveletAnalyzer::WaveletAnalyzer()
{
//
// Initialize instance variables:
//
  _analyzerOutput       = NULL;
  _analyzerInput        = NULL;
  _lowPassFilter        = NULL;
  _highPassFilter       = NULL;
  _oldInputSize         = 0;
  _coefficientGenerator = new WaveletFilterCoefficients();
  
  setFilterLength(DEFAULT_FILTER_LENGTH);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the WaveletAnalyzer class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
WaveletAnalyzer::~WaveletAnalyzer()
{
  delete [] _analyzerOutput;
  delete [] _analyzerInput;
  delete    _lowPassFilter;
  delete    _highPassFilter;
  delete    _coefficientGenerator;
  return;
}

// ############################# Public Method ###############################
// setFilterLength -- Sets a new wavelet filter coefficient size/length.
// Input:       length:         new coefficient length
//          
// Output:                      none
// Notes:
// 1. _lowPassFilter, _highPassFilter, and _coefficientGenerator
//    are modified.
// ############################# Public Method ###############################
void WaveletAnalyzer::setFilterLength(int length)
{
  delete _lowPassFilter;
  delete _highPassFilter;
  _coefficientGenerator->setFilterLength(length);

  _lowPassFilter  = new RealFIR(_coefficientGenerator->calculateLowPassCoefficients(), length);
  _highPassFilter = new RealFIR(_coefficientGenerator->calculateHighPassCoefficients(), length);

  return;
}


// ############################# Public Method ###############################
// numberStages -- Calculates the number of actual stages that can be implemented
//                 given the desired number of stages.
//
// Input:       length:         the input data length
//              desiredStages:  desired number of stages
//          
// Output:                      actual number of stages if input length does
//                              not divide evenly by two, the desired-number
//                              of-stages times.
// Notes:
// ############################# Public Method ###############################
int WaveletAnalyzer::numberStages(int length, int desiredStages)
{
  int actual_stages = 0;
  while(length > 1)
  {
    if(length%2 == 0)           // still divisible by 2
    {
      actual_stages++;
      length    /= 2;
    }
    else
      break;
  }
  return MIN(desiredStages, actual_stages);
}

// ############################# Public Method ###############################
// waveletAnalyze -- Analyze the input data using wavelet filter banks.
//
// Input:       inputData:      set of input samples
//              inputLength:    length of the above
//              numberStages:   desired number of stages of analyzer
//          
// Output:                      actual number of stages if input length does
//                              not divide evenly by two, the desired-number
//                              of-stages times.
// Notes:
// ############################# Public Method ###############################
int WaveletAnalyzer::waveletAnalyze(const Complex *inputData, int inputLength, int numberStages)
{
//
// Allocate the output arrays:
//
  delete [] _analyzerOutput;
  _analyzerOutput       = new Complex [inputLength];
  Complex *temp_output  = new Complex [inputLength];
//
// See Figure 10.12 in Pace' book.  We start at the bottom with the
// H filter, and then filter as H GG HH GG HH ..
//
  int   number_filters          = 2;                    // # of filters at each stage
  int   data_size               = inputLength;
  Complex *output_ptr           = _analyzerOutput;
  const Complex *input_ptr      = inputData;
  char  output_array            = 'A';                  // Analyzer output
  int   actual_stages;
  for(actual_stages=0; actual_stages<numberStages; actual_stages++)
  {
    char current_filter         = 'H';
    char next_filter            = 'G';
    /*
    if(actual_stages == (numberStages-1))
    {
      int index                 = 0;
      for(int j=0; j<number_filters; j++)
      {
//
// Filter the data, and store it in output
//

        for(int k=0; k<data_size/2; k++)
        {
          _analyzerOutput[index]        = temp_output[index] = Complex(j, 0.);
          index++;
        }
//
// Update the filter and data pointers
//
        output_ptr      += data_size/DECIMATION_SIZE;
      }
      break;
    }
    */
    for(int j=0; j<number_filters; j++)
    {
//
// Filter the data, and store it in output
//
      if(current_filter == 'H')
        filterData(input_ptr, output_ptr, data_size, _lowPassFilter);
      else
        filterData(input_ptr, output_ptr, data_size, _highPassFilter);
//
// Update the filter and data pointers
//
      int ptr_inc       = updateFilterPointer(&current_filter, &next_filter, data_size);
      input_ptr         += ptr_inc;
      output_ptr        += data_size/DECIMATION_SIZE;
    }
//
// change the data size, number of filters, and input and output pointers:
//
    if(data_size%2 != 0)                        // Keep analyzing until numberStages
      break;
    data_size                   /= DECIMATION_SIZE;
    number_filters              *= DECIMATION_SIZE;
    if(output_array == 'A')
    {
      input_ptr                 = _analyzerOutput;
      output_ptr                = temp_output;
      output_array              = 'T';
    }
    else
    {
      input_ptr                 = temp_output;
      output_ptr                = _analyzerOutput;
      output_array              = 'A';
    }
  }
//
// Copy temp into output, if temp contains output
//
  if(output_array == 'A')
  {
    for(int i=0; i<inputLength; i++)
    {
      _analyzerOutput[i]        = temp_output[i];
    }
  }

  delete [] temp_output;
  return actual_stages;
}


// ############################# Public Method ###############################
// waveletAnalyze -- Analyze the input data using wavelet filter banks.
//
// Input:       inputReal:      set of input samples real part
//              inputImag:      set of input samples imag part
//              inputLength:    length of the above
//              numberStages:   desired number of stages of analyzer
//          
// Output:                      actual number of stages if input length does
//                              not divide evenly by two, the desired-number
//                              of-stages times.
// Notes:
// 1. Just copies real and imaginary to a complex array, and calls the above
//    routine.
// ############################# Public Method ###############################
int WaveletAnalyzer::waveletAnalyze(const float *inputReal, const float *inputImag, int inputLength, int numberStages)
{
//
// Allocate array for copying input
//
  if(_oldInputSize < inputLength)
  {
    _oldInputSize       = inputLength;
    delete [] _analyzerInput;
    _analyzerInput      = new Complex[inputLength];
  }
  for(int i=0; i<inputLength; i++)
  {
    _analyzerInput[i]   = Complex(inputReal[i], inputImag[i]);
  }
  return waveletAnalyze(_analyzerInput, inputLength, numberStages);
}