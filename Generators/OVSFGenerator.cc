/********************************************************************************
 * This subclass of Object generates orthogonal variable spreading factor (VSF) *
 * codes for WCDMA and stores them in arrays for output.                        *
 *                                                                              *
 * File: OVSFGenerator.cc                                                       *
 *                                                                              *
 *                                                                              *
 * Notes:                                                                       *
 * 1. The VSF code functions have values of 0 or 1.                             *
 * 2. Reference is "Spreading codes for direct sequence CDMA and wideband       *
 *    CDMA cellular networks," Dinan and Jabbari, IEEE Comm. Mag. Sept. 1998.   *
 *                                                                              *
 ********************************************************************************/
#include "OVSFGenerator.h"
#include <stdio.h>

#ifndef YES
#define YES     1
#define NO      0
#endif

#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define BIT_MAP(a)      ( ((a)==0 ) ? 1 : -1 )          // (0,1) -> (1, -1)

// ############################ Public Function ###################################
// isPowerOfTwo - Returns whether or not the input integer is a power of 2.
//
// Input:           n:      integer to be tested
// Output:                  YES if n is a power of 2
//
// ############################ Public Function ###################################
LOGICAL OVSFGenerator:: isPowerOfTwo(int n)
{
  int   old_n, two_to_n;
  int   powerOfTwo;
  
  powerOfTwo = -1;

  old_n = n;
  while(1)
  {
    if(n==0)
      break;
    n >>= 1;            // Divide by 2
    powerOfTwo++;
  }
  if(powerOfTwo < 1)
    return NO;
  two_to_n  = 1<<powerOfTwo;
  if(old_n == two_to_n)
    return YES;
  return NO;
}

// ############################# Private Method ###############################
// generateVSFCode  Generates a VSF code recursively
//
// Input:       vsf:            code array
//              length:         length of code
//              row:            row of code
//          
// Output:                      vsf is modified on output
//
// Notes:
// ############################# Private Method ###############################
void OVSFGenerator::generateVSFCode(char *vsf, int length, int row)
{
  int   i, new_row, half_length;
//
  if(length == 1)
  {
    vsf[0]      = 0;                                    // reached end of recursion
    return;
  }
  else
  {
    if( (row%2) == 0)
      new_row   = row/2;                                // even row
    else
      new_row   = (row+1)/2;
    half_length = length/2;
    generateVSFCode(vsf, half_length, new_row);         // recursive call
    for(i=0; i<length/2; i++)                           // fill in length/2
    {
      if( (row%2) == 0 )
      {                                                 // even row
        if(vsf[i] == 0)
          vsf[i+half_length]    = 1;                    // complement
        else
          vsf[i+half_length]    = 0;
      }
      else
        vsf[i+half_length]      = vsf[i];               // no complement
    }
  }
  return;
}

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the OVSFGenerator class.
//
// Input:       length:         Length of the code
//              row:            Row of the tree
//          
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
OVSFGenerator::OVSFGenerator(int length, int row)
{
//
// Error check inputs:
//
  _vsfCode      = NULL;
  _codeLength   = MAX(1, length);
  if(!isPowerOfTwo(_codeLength))
  {
    printf("Length = %d, not power of two in OVSFGenerator\n", _codeLength);
    _codeLength = 2;
  }
  _codeRow      = MIN(row, _codeLength);
//
// Allocate and generate code:
//
  _vsfCode      = new char[_codeLength];
  generateVSFCode(_vsfCode, _codeLength, _codeRow);

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the OVSFGenerator class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
OVSFGenerator::~OVSFGenerator()
{
  delete [] _vsfCode;
  return;
}
