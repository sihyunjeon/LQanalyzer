#ifndef ConversionValidation_h
#define ConversionValidation_h

#include "AnalyzerCore.h"
class ConversionValidation : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ConversionValidation();
  ~ConversionValidation();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  void DrawHistColl( TString this_string, snu::KElectron electron, double this_weight );

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ConversionValidation, 1);
};
#endif
