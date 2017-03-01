#ifndef DataValidation_EMu_h
#define DataValidation_EMu_h

#include "AnalyzerCore.h"
class DataValidation_EMu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  DataValidation_EMu();
  ~DataValidation_EMu();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( DataValidation_EMu, 1);
};
#endif
