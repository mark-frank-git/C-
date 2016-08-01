#ifndef _ENVIRONMENT_TYPES_H
#define _ENVIRONMENT_TYPES_H    1
/************************************************************************
 *                                                                      *
 * This header file defines the environment types specified in COST 207 *
 *                                                                      *
 *                                                                      *
 * File:environment_types.h                                             *
 *                                                                      *
 *                                                                      *
 * See: GSM 05.05 V8.2.0, Annex C.                                      *
 *                                                                      *
 ************************************************************************/

//
// _environmentType
//
#define RA_TYPE                 0               // Rural area
#define HT_TYPE                 1               // Hilly terrain
#define TU_TYPE                 2               // Typical urban
#define EQ_TYPE                 3               // Equalization test
#define TI_TYPE                 4               // Small cells
#define STATIC_TYPE             5               // Static (no fading) channel
#define TEST1                   6               // Test channel
#define TEST2                   7               // Test channel
#define TEST3                   8
#define TEST4                   9
#define NUMBER_ENV_TYPES        (TEST4+1)       // Total number of environment types

//
// The filter types define the Doppler spectrum filters:
//
#define CLASSICAL_TYPE  0                                               // Classical U shaped
#define RICIAN_TYPE     1                                               // Contains direct path
#define GAUSSIAN_TYPE   2
#define NO_FADING_TYPE  3                                               // No fading for static case
#endif
