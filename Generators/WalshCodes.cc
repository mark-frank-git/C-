/********************************************************************************
 * This subclass of Object generates Walsh functions (codes) and stores them    *
 * in arrays for output.                                                        *
 *                                                                              *
 * File: WalshCodes.c                                                           *
 *                                                                              *
 * Revision history:                                                            *
 *   1. 04/15/01 - Started                                                      *
 * Notes:                                                                       *
 * 1. See Lewis Franks, "Signal Theory" for Walsh function generation formulae. *
 * 2. The Walsh code functions have values of 0 or 1.                           *
 * 3. These are really indexed according to the Hadamard matrix as in IS-95.    *
 *                                                                              *
 ********************************************************************************/
#include "WalshCodes.h"
#include <stdio.h>

#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )

// ############################# Private Method ###############################
// log2  base 2 logarithm
//
// Input:       n:              integer to take log2() of
//          
// Output:      log2(n)
//
// ############################# Private Method ###############################
int WalshCodes::log2(int n)
{
  unsigned int mask,i;

  if(n == 0) 
    return -1;                /* zero is an error, return -1 */

  n--;                        /* get the max index, x-1     */

  for(mask = 1 , i = 0 ; ; mask *= 2 , i++)
  {
    if(n == 0) 
      return i;             /* return log2 if all zero   */
    n &= (~mask);           /* AND off a bit             */
  }
  return n;
}

// ############################# Private Method ###############################
// ipow  returns i^n
//
// Input:       i:              integer to be raised to power
//              n:              power n
//          
// Output:                      i^n
//
// ############################# Private Method ###############################
int WalshCodes::ipow(int i, int n)
{
  int   pow;
  
  pow   = 1;
  while(n>0)
  {
    pow *= i;
    n--;
  }
  return pow;
}

// ############################# Private Method ###############################
// generateHadamardMatrix  Generates a new Hadamard matrix given an old one
//
// Input:       hn:             Hn matrix
//              n:              size of matrix
//          
// Output:                      H2n matrix
//
// Notes:
// 1. It is up to the calling routine to dealloc the generated matrix
// 2. The matrix is generated row-wise.
//
//        Hn   Hn
// H2n =
//        Hn  !Hn
// ############################# Private Method ###############################
int *WalshCodes::generateHadamardMatrix(int *hn, int n)
{
  int   i, j, index;
  int   *h2n;
//
// allocate new matrix:
//
  h2n   = new int[2*n*2*n];
  index = 0;
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
      h2n[index++]      = hn[i*n+j];
    for(j=0; j<n; j++)
      h2n[index++]      = hn[i*n+j];
  }
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
      h2n[index++]      = hn[i*n+j];
    for(j=0; j<n; j++)
    {
      if(hn[i*n+j] == 1)
        h2n[index]      = 0;
      else
        h2n[index]      = 1;
      index++;
    }
  }
  return h2n;
}

// ############################# Private Method ###############################
// generateCodes  Generates all of the Walsh codes
//
// Input:                       None
//          
// Output:                      None
//
// ############################# Private Method ###############################
void WalshCodes::generateCodes()
{
  int   i, j, size;
  int   *h2n, *hn;
  int   log2_code, index;
//
// iteratively generate the hadamard matrices:
//
  hn    = new int[4];
  hn[0] = 0;
  hn[1] = 0;
  hn[2] = 0;
  hn[3] = 1;
  log2_code     = log2(_numberCodes);
  size  = 2;
  for(i=1; i<log2_code; i++)
  {
    h2n         = generateHadamardMatrix(hn, size);
    delete [] hn;
    hn          = h2n;
    size        *= 2;
  }
//
// Copy the hadamard matrix into the code vectors:
//
  index         = 0;
  for(i=0; i<_numberCodes; i++)
  {
    for(j=0; j<_codeLength; j++)
      _walshCodes[i][j] = hn[index++];
  }
  delete [] hn;
  return;
}
      

// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the WalshCodes class.
//
// Input:       numberCodes:    number of codes to generate, should be power of 2
//          
// Output:                      None
//
// Notes:
// ############################# Class Constructor ###############################
WalshCodes::WalshCodes(int numberCodes)
{
  int   i;
  _numberCodes  = MAX(2, numberCodes);
  _numberCodes  = MIN(MAX_CODES, _numberCodes);
  _codeLength   = _numberCodes;
  _walshCodes   = new char * [_numberCodes];
  for(i=0; i<_numberCodes; i++)
    _walshCodes[i]      = new char[_codeLength];
//
// Generate the codes:
//
  generateCodes();

  return;
}

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the WalshCodes class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
WalshCodes::~WalshCodes()
{
  int   i;
  for(i=0; i<_numberCodes; i++)
    delete [] _walshCodes[i];
  delete [] _walshCodes;
  return;
}

// ############################# Public Function ###############################
// codeForIndex -- Returns the code at the given index.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
const char *WalshCodes::codeForIndex(int index)
{
  index = MAX(0, index);
  index = MIN(index, (_numberCodes-1));
  return _walshCodes[index];
}

// ############################# Public Function ###############################
// printCodes -- Prints out all of the codes.
//
// Input:                       None
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void WalshCodes::printCodes()
{
  int   i, j;

  for(i=0; i<_numberCodes; i++)
  {
    if(i<10)
     printf("i =  %d: ", i);
    else
     printf("i = %d: ", i);
    for(j=0; j<_codeLength; j++)
    {
      if(_walshCodes[i][j] < 0)
        printf(" %d", _walshCodes[i][j]);
      else
        printf("  %d", _walshCodes[i][j]);
    }
    printf("\n");
  }
  return;
}
