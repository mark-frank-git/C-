/********************************************************************************
 * This subclass of Object implements a QAM constellation mapper.               *
 * File: QAMMapper.h                                                            *
 *                                                                              *
 ********************************************************************************/
#include "QAMMapper.h"
#include <stdio.h>


#define MAX(a,b)        ( ((a)<(b)) ? (b) : (a) )
#define MIN(a,b)        ( ((a)>(b)) ? (b) : (a) )
#define ROUND(a)        ( ((a)>0.) ? (int)((a)+0.5) : (int)((a)-0.5) )

#ifndef YES
#define YES 1
#define NO  0
#endif

//
// The offset for converting from index to I and Q values:
//
int     qam_offset[MAX_BITS+1]  = {0, 0, 1, 1, 3, 5, 7, 9, 15};
//
// The bit mappings for the indexing (put in an index, get out a bit map):
//
int     dvbc_bit_16[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
int     dvbc_bit_32[32] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int     dvbc_bit_64[64] = {60, 61, 57, 56, 20, 22, 30, 28, 62, 63, 59, 58, 21, 23, 31, 29,
                           54, 55, 51, 50, 17, 19, 27, 25, 52, 53, 49, 48, 16, 18, 26, 24,
                           40, 42, 34, 32,  0,  1,  5,  4, 41, 43, 35, 33,  2,  3,  7,  6,
                           45, 47, 39, 37, 10, 11, 15, 14, 44, 46, 38, 36,  8,  9, 13, 12};
int     dvbc_bit_128[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127};
int     dvbc_bit_256[256] = {240, 241, 245, 244, 228, 229, 225, 224, 80, 82, 90, 88, 120, 122, 114, 112,
                             242, 243, 247, 246, 230, 231, 227, 226, 81, 83, 91, 89, 121, 123, 115, 113,
                             250, 251, 255, 254, 238, 239, 235, 234, 85, 87, 95, 93, 125, 127, 119, 117,
                             248, 249, 253, 252, 236, 237, 233, 232, 84, 86, 94, 92, 124, 126, 118, 116,
                             216, 217, 221, 220, 204, 205, 201, 200, 68, 70, 78, 76, 108, 110, 102, 100,
                             218, 219, 223, 222, 206, 207, 203, 202, 69, 71, 79, 77, 109, 111, 103, 101,
                             210, 211, 215, 214, 198, 199, 195, 194, 65, 67, 75, 73, 105, 107, 99, 97,
                             208, 209, 213, 212, 196, 197, 193, 192, 64, 66, 74, 72, 104, 106, 98, 96,
                             160, 162, 170, 168, 136, 138, 130, 128,  0,  1,  5,  4,  20,  21, 17, 16,
                             161, 163, 171, 169, 137, 139, 131, 129,  2,  3,  7,  6,  22,  23, 19, 18,
                             165, 167, 175, 173, 141, 143, 135, 133, 10, 11, 15, 14,  30,  31, 27, 26,
                             164, 166, 174, 172, 140, 142, 134, 132,  8,  9, 13, 12,  28,  29, 25, 24,
                             180, 182, 190, 188, 156, 158, 150, 148, 40, 41, 45, 44,  60,  61, 57, 56,
                             181, 183, 191, 189, 157, 159, 151, 149, 42, 43, 47, 46,  62,  63, 59, 58,
                             177, 179, 187, 185, 153, 155, 147, 145, 34, 35, 39, 38,  54,  55, 51, 50,
                             176, 178, 186, 184, 152, 154, 146, 144, 32, 33, 37, 36,  52,  53, 49, 48};
int     j83_bit_16[16]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
int     j83_bit_32[32]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int     j83_bit_64[64]  = {54, 60, 38, 44, 18, 30, 50, 62, 51, 57, 35, 41, 17, 29, 49, 61,
                           22, 28,  6, 12,  2, 14, 34, 46, 19, 25,  3,  9,  1, 13, 33, 45,
                           36, 40,  4,  8,  0, 10, 16, 26, 39, 43,  7, 11,  5, 15, 21, 31,
                           52, 56, 20, 24, 32, 42, 48, 58, 55, 59, 23, 27, 37, 47, 53, 63};
int     j83_bit_128[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127};
int     j83_bit_256[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
                           141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                           156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
                           171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                           186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
                           201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
                           211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
                           226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
                           241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

//
// The indexes for the bit maps (put in a bit map # get out an index)
//
int     dvbc_index_16[16]       = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
int     dvbc_index_32[32]       = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int     dvbc_index_64[64]       = {36, 37, 44, 45, 39, 38, 47, 46, 60, 61, 52, 53, 63, 62, 55, 54,
                                   28, 20, 29, 21,  4, 12,  5, 13, 31, 23, 30, 22,  7, 15,  6, 14,
                                   35, 43, 34, 42, 59, 51, 58, 50, 32, 40, 33, 41, 56, 48, 57, 49,
                                   27, 26, 19, 18, 24, 25, 16, 17,  3,  2, 11, 10,  0,  1,  8,  9};
int     dvbc_index_128[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127};
int     dvbc_index_256[256] = {136, 137, 152, 153, 139, 138, 155, 154, 184, 185, 168, 169, 187, 186, 171, 170, 
                               143, 142, 159, 158, 140, 141, 156, 157, 191, 190, 175, 174, 188, 189, 172, 173, 
                               248, 249, 232, 233, 251, 250, 235, 234, 200, 201, 216, 217, 203, 202, 219, 218, 
                               255, 254, 239, 238, 252, 253, 236, 237, 207, 206, 223, 222, 204, 205, 220, 221, 
                               120, 104, 121, 105, 72, 88, 73, 89, 123, 107, 122, 106, 75, 91, 74, 90, 
                               8, 24, 9, 25, 56, 40, 57, 41, 11, 27, 10, 26, 59, 43, 58, 42, 
                               127, 111, 126, 110, 79, 95, 78, 94, 124, 108, 125, 109, 76, 92, 77, 93, 
                               15, 31, 14, 30, 63, 47, 62, 46, 12, 28, 13, 29, 60, 44, 61, 45, 
                               135, 151, 134, 150, 183, 167, 182, 166, 132, 148, 133, 149, 180, 164, 181, 165, 
                               247, 231, 246, 230, 199, 215, 198, 214, 244, 228, 245, 229, 196, 212, 197, 213, 
                               128, 144, 129, 145, 176, 160, 177, 161, 131, 147, 130, 146, 179, 163, 178, 162, 
                               240, 224, 241, 225, 192, 208, 193, 209, 243, 227, 242, 226, 195, 211, 194, 210, 
                               119, 118, 103, 102, 116, 117, 100, 101, 71, 70, 87, 86, 68, 69, 84, 85, 
                               112, 113, 96, 97, 115, 114, 99, 98, 64, 65, 80, 81, 67, 66, 83, 82, 
                               7, 6, 23, 22, 4, 5, 20, 21, 55, 54, 39, 38, 52, 53, 36, 37, 
                               0, 1, 16, 17, 3, 2, 19, 18, 48, 49, 32, 33, 51, 50, 35, 34};
int     j83_index_16[16]        = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
int     j83_index_32[32]        = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
int     j83_index_64[64]        = {36, 28, 20, 26, 34, 44, 18, 42, 35, 27, 37, 43, 19, 29, 21, 45,
                                   38, 12,  4, 24, 50, 46, 16, 58, 51, 25, 39, 59, 17, 13,  5, 47,
                                   52, 30, 22, 10, 32, 60,  2, 40, 33, 11, 53, 41,  3, 31, 23, 61,
                                   54, 14,  6,  8, 48, 62,  0, 56, 49,  9, 55, 57,  1, 15,  7, 63};
int     j83_index_128[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127};
int     j83_index_256[256] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                           32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                           48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
                           64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                           80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                           96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                           111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                           126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
                           141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                           156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
                           171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                           186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
                           201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
                           211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
                           226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
                           241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};


// ############################# Class Constructor ###############################
// Class Constructor -- Constructor for the QAMMapper class.
// Input:               type:           QAM type, DVB-C, etc.
//                      bits:           # of bits per constellation point
//          
// Output:                              None
//
// Notes:
// ############################# Class Constructor ###############################
QAMMapper::QAMMapper(int type, int bits)
{
  setQAMType(type);
  setQAMBits(bits);
  return;
}
  

// ############################# Class Destructor ###############################
// Class Destructor -- Destructor for the QAMMapper class.
//
// Input:                       None
//
// Output:                      None
//
// Notes:
// ############################# Class Destructor ###############################
QAMMapper::~QAMMapper()
{
  return;
}

// ############################# Public Function ###############################
// setQAMType -- Sets a new QAM type.
//
// Input:       type:           new QAM type
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void QAMMapper::setQAMType(int type)
{
  _qamType      = type;
  return;
}

// ############################# Public Function ###############################
// setQAMBits -- Sets a new number of QAM bits.
//
// Input:       bits:           new number of QAM bits
//          
// Output:                      None
//
// Notes:
// ############################# Public Function ###############################
void QAMMapper::setQAMBits(int bits)
{
  _qamBits      = MAX(2, bits);
  _qamBits      = MIN(_qamBits, MAX_BITS);

  _qamPoints    = 1 << _qamBits;
  _qamLevels    = 1 << (_qamBits/2);                    // doesn't work for odd _qamBits
  return;
}

// ############################# Public Function ###############################
// iAndQAtIndex -- Returns the I and Q constellation values for the given index.
//
// Input:       index:          constellation index
//          
// Output:      iData:          I component of constellation point
//              qData:          Q component of constellation point
//
// Notes:
// 1. The index starts at 0 in the lower left corner of the constellation and
//    increases by 1s along the I direction.
// 2. The constellation is normalized so that the I and Q values are at +/- 1,
//    +/- 3, +/- 5, etc.
// ############################# Public Function ###############################
void    QAMMapper::iAndQAtIndex(int index, float *iData, float *qData)
{
  int           i_offset, q_offset;
//
// Get I and Q offsets from index
//
  i_offset      = index % _qamLevels;
  q_offset      = index / _qamLevels;
//
// Now, map to I and Q levels:
//
  i_offset      -= _qamLevels/2;
  if(i_offset >= 0)
    *iData      = 1. + i_offset*2;
  else
    *iData      = -1. + (i_offset+1)*2;
//
  q_offset      -= _qamLevels/2;
  if(q_offset >= 0)
    *qData      = 1. + q_offset*2;
  else
    *qData      = -1. + (q_offset+1)*2;
  return;
}

// ############################# Public Function ###############################
// iAndQForBitMap -- Returns the I and Q constellation values for the bit pattern.
//
// Input:       bitMap:         constellation bit map
//          
// Output:      iData:          I component of constellation point
//              qData:          Q component of constellation point
//
// Notes:
// 1. The constellation is normalized so that the I and Q values are at +/- 1,
//    +/- 3, +/- 5, etc.
// ############################# Public Function ###############################
void    QAMMapper::iAndQForBitMap(int bitMap, float *iData, float *qData)
{
  int   index;
//
// Get the bits out of the mapping table
//
  switch(_qamBits)
  {
    default:
    case 4:
      bitMap    = MAX(0, bitMap);
      bitMap    = MIN(bitMap, 15);
      if(_qamType == QAM_DVBC)
        index   = dvbc_index_16[bitMap];
      else
        index   = j83_index_16[bitMap];
      break;
    case 5:
      bitMap    = MAX(0, bitMap);
      bitMap    = MIN(bitMap, 31);
      if(_qamType == QAM_DVBC)
        index   = dvbc_index_32[bitMap];
      else
        index   = j83_index_32[bitMap];
      break;
    case 6:
      bitMap    = MAX(0, bitMap);
      bitMap    = MIN(bitMap, 63);
      if(_qamType == QAM_DVBC)
        index   = dvbc_index_64[bitMap];
      else
        index   = j83_index_64[bitMap];
      break;
    case 7:
      bitMap    = MAX(0, bitMap);
      bitMap    = MIN(bitMap, 127);
      if(_qamType == QAM_DVBC)
        index   = dvbc_index_128[bitMap];
      else
        index   = j83_index_128[bitMap];
      break;
    case 8:
      bitMap    = MAX(0, bitMap);
      bitMap    = MIN(bitMap, 255);
      if(_qamType == QAM_DVBC)
        index   = dvbc_index_256[bitMap];
      else
        index   = j83_index_256[bitMap];
      break;
  }
  iAndQAtIndex(index, iData, qData);
  return;
}


// ############################# Public Function ###############################
// bitMappingAtIndex -- Returns the constellation point bit mapping.
//
// Input:       index:          constellation point index
//          
// Output:      iData:          I component of constellation point
//              qData:          Q component of constellation point
//
// Notes:
// 1. The index starts at 0 in the lower left corner of the constellation and
//    increases by 1s along the I direction.
// 2. The bit mapping is stored in an integer for return.
// ############################# Public Function ###############################
int     QAMMapper::bitMappingAtIndex(int index)
{
  int   bit_map;
//
// Get the bits out of the mapping table
//
  switch(_qamBits)
  {
    default:
    case 4:
      index     = MAX(0, index);
      index     = MIN(index, 15);
      if(_qamType == QAM_DVBC)
        bit_map = dvbc_bit_16[index];
      else
        bit_map = j83_bit_16[index];
      break;
    case 5:
      index     = MAX(0, index);
      index     = MIN(index, 31);
      if(_qamType == QAM_DVBC)
        bit_map = dvbc_bit_32[index];
      else
        bit_map = j83_bit_32[index];
      break;
    case 6:
      index     = MAX(0, index);
      index     = MIN(index, 63);
      if(_qamType == QAM_DVBC)
        bit_map = dvbc_bit_64[index];
      else
        bit_map = j83_bit_64[index];
      break;
    case 7:
      index     = MAX(0, index);
      index     = MIN(index, 127);
      if(_qamType == QAM_DVBC)
        bit_map = dvbc_bit_128[index];
      else
        bit_map = j83_bit_128[index];
      break;
    case 8:
      index     = MAX(0, index);
      index     = MIN(index, 255);
      if(_qamType == QAM_DVBC)
        bit_map = dvbc_bit_256[index];
      else
        bit_map = j83_bit_256[index];
      break;
  }
  return bit_map;
}

// ############################# Public Function ###############################
// getIndexAtIAndQ -- Returns the constellation point given the I and Q values.
//
// Input:       iData:          I component of constellation point
//              qData:          Q component of constellation point
//          
// Output:      index:          constellation point index
//
// Notes:
// 1. The index starts at 0 in the lower left corner of the constellation and
//    increases by 1s along the I direction.
// ############################# Public Function ###############################
int     QAMMapper::indexAtIAndQ(float iData, float qData)
{
  int   offset, i_index, q_index, index;
//
// Get constellation point index for
//
  offset        = qam_offset[_qamBits];
  i_index       = (ROUND( (iData + offset)) )/2;
  q_index       = (ROUND( (qData + offset)) )/2;
  index         = i_index + q_index*_qamLevels;
  return index;
}


// ############################# Public Function ###############################
// indexForBitMap -- Returns the constellation table given the index.  This function
//                   generates the index table given the bit map table.
//
// Input:       bitMap:         QAM bit mapping
//          
// Output:      index:          constellation point index
//
// Notes:
// 1. The index starts at 0 in the lower left corner of the constellation and
//    increases by 1s along the I direction.
// 2. This function is slow, and should only be used for generation of the index
//    table.
// ############################# Public Function ###############################
int     QAMMapper::indexForBitMap(int bitMap)
{
  int   i;
  int   *table, index;
//
// Get the bits out of the mapping table
//
  switch(_qamBits)
  {
    default:
    case 4:
      if(_qamType == QAM_DVBC)
        table   = dvbc_bit_16;
      else
        table   = j83_bit_16;
      break;
    case 5:
      if(_qamType == QAM_DVBC)
        table   = dvbc_bit_32;
      else
        table   = j83_bit_32;
      break;
    case 6:
      if(_qamType == QAM_DVBC)
        table   = dvbc_bit_64;
      else
        table   = j83_bit_64;
      break;
    case 7:
      if(_qamType == QAM_DVBC)
        table   = dvbc_bit_128;
      else
        table   = j83_bit_128;
      break;
    case 8:
      if(_qamType == QAM_DVBC)
        table   = dvbc_bit_256;
      else
        table   = j83_bit_256;
      break;
  }
//
// Now, search for the index:
//
  index         = 0;
  for(i=0; i< _qamPoints; i++)
  {
    if(table[i] == bitMap)
    {
      index     = i;
      break;
    }
  }
  return index;
}
