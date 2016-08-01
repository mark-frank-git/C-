/***************************************************************************//*!
**  \file         SymbolEstimate.h
**  \brief        symbol estimator
**  \description  This class takes input I/Q data and outputs symbol estimates
**  \project      BlockDraw
**  \author       Mark Frank
*/ /*
*******************************************************************************/
#ifndef _SYMBOL_ESTIMATE_H
#define _SYMBOL_ESTIMATE_H    1

////////////////////////////////////////////////////////////////////////////////
// enum
typedef enum
{
  eMod_BPSK,
  eMod_QPSK
} eMod_ModType;


////////////////////////////////////////////////////////////////////////////////
// Structures

// Symbol structure
typedef struct
{
  double        iAverage;               //!< Average I value 
  double        qAverage;               //!< Average Q value
  int           numberSymbols;          //!< Number of symbols of this type
} SymbolEstimateStruct;



/***************************************************************************//**
**  \class          SymbolEstimate    
**  \brief          symbol estimator
**  \description    This class takes input I/Q data and outputs symbol estimates
*******************************************************************************/
class SymbolEstimate
{
  public:
  // LIFECYCLE
  SymbolEstimate (eMod_ModType modTYpe);                //!< Default constructor
  virtual ~ SymbolEstimate ();                          //!< Default destructor
  
  // PUBLIC MEMBER FUNCTIONS
  void ResetEstimates();                                //!< Reset all of the averages
  
  void AccumulateSymbols(                               //!< Accumulate estimate statistics for an array of symbols
                  const double *iData,                  //!< Input I data
                  const double *qData,                  //!< Input Q data
                  const int    numberSymbols);          //!< Size of above arrays

  void AccumulateSymbol(                                //!< Accumulate estimate statistics for a single symbol
                  const double  iData,                  //!< Input I data
                  const double  qData);                 //!< Input Q data
  
  void EstimateSymbols();                               //!< Find the symbols from the accumulated statistics

/*******************************
 * These methods set parameters*
 *******************************/
  void SetMaxNumberSymbols(                             //!< Sets the max number of symbols to accumulate
                          const int maxNumber);         //!< New max number of symbols

/*******************************
 * These methods get parameters*
 *******************************/
  const double *iData() const                           //!< Return an array of I values for all estimated symbols
  {
    return m_iData;
  }
  const double *qData() const                           //!< Return an array of Q values for all estimated symbols
  {
    return m_qData;
  }

protected:
  // PROTECTED MEMBER FUNCTIONS
  //
  void SetModType(                                      //!< Sets the modulation type
             const eMod_ModType modType);               //!< New modulation type
  
  const int GetSymbolType(                              //!< Gets the symbol type (number)
                          const eMod_ModType modType,   //!< modulation type
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const;          //!< Q data for symbol type to estimate
  
  const int GetSymbolTypeForBPSK(                       //!< Gets the symbol type (number) for BPSK
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const;          //!< Q data for symbol type to estimate
  
  const int GetSymbolTypeForQPSK(                       //!< Gets the symbol type (number) for QPSK
                          double iData,                 //!< I data for symbol type to estimate
                          double qData) const;          //!< Q data for symbol type to estimate
    
  // PROTECTED MEMBER DATA
  //
  eMod_ModType          m_modulationType;               //!< Modulation type
  int                   *m_symbolTypes;                 //!< Array of estimated symbol types
  int                   m_numberSymbols;                //!< Number of symbols estimated
  int                   m_numberSymbolTypes;            //!< Number of symbol types (depends on mod type)
  int                   m_maxNumberSymbols;             //!< Max # of symbols to estimate
  SymbolEstimateStruct  *m_estimatedSymbols;            //!< Estimated symbols
  double                *m_iData;                       //!< I values of estimated symbols
  double                *m_qData;                       //!< Q values of estimated symbols

 private:
  // PRIVATE MEMBER FUNCTIONS
  //

  // These are not implemented, just listed here so they are not implicitly used:
  // lint !e1704
  SymbolEstimate & operator = (const SymbolEstimate & rhs);           //!< Disabled assignment operator
  // lint !e1704
  SymbolEstimate (const SymbolEstimate & rhs);                              //!< Disabled copy constructor
    
  // PRIVATE MEMBER DATA


};
#endif
