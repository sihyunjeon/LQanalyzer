#ifndef _SKTree_KEvent_H__
#define _SKTree_KEvent_H__

/// Local includes
#include "KParticle.h"

#include <iosfwd>
#include <string>
#include <vector>
#include "TObject.h"

namespace snu {
  
  class KEvent : public TObject{
  public:
    KEvent();
    
    ///Copy constructor
    KEvent(const KEvent& el);
    
    ///Destructor    
    ~KEvent() ;

    KEvent& operator= (const KEvent& obj);
    
    /// MET variabkes
    void SetPFMET(double met); ///defailt MET for analysis
    void SetPFMETphi(double metphi);
    void SetPFSumET(double sumet);
    
    void SetPuppiMET(double met); /// puppi MET for analysis
    void SetPuppiMETphi(double metphi);
    void SetPuppiSumET(double sumet);

    void SetNoHFMET(double met); /// puppi MET for analysis
    void SetNoHFMETphi(double metphi);
    void SetNoHFSumET(double sumet);


    void SetPfMvaMET(double met); /// puppi MET for analysis
    void SetPfMvaMETphi(double metphi);
    void SetPfMvaSumET(double sumet);


    
    /// PDF
    //void SetPDFCTEQWeight( std::vector<double> pdf);
    //void SetPDFMSTWWeight( std::vector<double> pdf);
    //void SetPDFNNPDFWeight( std::vector<double> pdf);
    
    /// Process ID
    void SetLumiSection(int ls);
    /// Vertex
    void SetNVertices(int nvert);
    
    void SetNGoodVertices(int nvert);

    //    void SetVertexIsFake(bool isfakeV);
    void SetIsGoodEvent(int isgood);
    void SetVertexX(double vX);
    void SetVertexY(double vY);
    void SetVertexZ(double vZ);    
    void SetVertexNDOF(double ndof);    

    /// Event
    void SetWeight(double weight);
    void SetGenId1(int  id1);
    void SetGenId2(int  id2);
    void SetLHEWeight(double  lhe);
    void SetGenX1(double  x1);
    void SetGenX2(double  x2);
    void SetGenQ(double  Q);

    void SetEventNumber(int ev);
    void SetRunNumber(int rn);
    void SetIsData(bool isdata);
    
    // Jet rhos
    //void SetJetRho(double jetrho);    
    
    /// MET filters

    void SetPassCSCHaloFilterTight(bool pass);
    void SetPassEcalDeadCellTriggerPrimitiveFilter(bool pass);
    void SetPassHBHENoiseFilter(bool pass);
    void SetPassBadEESupercrystalFilter(bool pass);

    /// PileUp reweighting (only in MC)
    void SetPileUpInteractionsTrue(double npu);
				   
    
    //// Functions to call in analysis code/selection code
    
    /// MET variables
    inline Double_t PFMET() const {return k_PF_MET;}
    inline Double_t PFMETphi() const {return k_PF_METphi;}

    inline Double_t PFSumET() const {return k_PF_SumET;}
    

    inline Double_t PuppiMET() const {return k_Puppi_MET;}
    inline Double_t PuppiMETphi() const {return k_Puppi_METphi;}
    inline Double_t PuppiSumET() const {return k_Puppi_SumET;}

    inline Double_t NoHFMET() const {return k_NoHF_MET;}
    inline Double_t NoHFMETphi() const {return k_NoHF_METphi;}
    inline Double_t NoHFSumET() const {return k_NoHF_SumET;}

    inline Double_t PfMvaMET() const {return k_PfMva_MET;}
    inline Double_t PfMvaMETphi() const {return k_PfMva_METphi;}
    inline Double_t PfMvaSumET() const {return k_PfMva_SumET;}


    //inline std::vector<Double_t> PDFCTEQWeight() const {return k_pdf_cteq;}
    //inline std::vector<Double_t> PDFMSTWWeight() const {return k_pdf_mstw;}
    //inline std::vector<Double_t> PDFNNPDFWeight() const {return k_pdf_nnpdf;}
    
    
    /// ProcessID

    inline Int_t LumiSection()  const {return k_lumisec;}

    /// Event specific
    inline Double_t MCWeight() const{return k_mcweight;}
    inline Int_t IsData() const {return k_isData;}
    inline Int_t EventNumber() const {return k_EventNumber;}
    inline Int_t RunNumber() const {return k_RunNumber;}
    
    // Vertex
    inline Int_t nVertices()  const {return k_nvertices;}
    inline Int_t nGoodVertices()  const {return k_ngoodvertices;}

    inline Bool_t HasGoodPrimaryVertex() const {return k_isgoodevent;}
    inline Double_t VertexX() const {return k_vertexX;}
    inline Double_t VertexY() const {return k_vertexY;}
    inline Double_t VertexZ() const {return k_vertexZ;}
    inline Double_t VertexNDOF() const {return k_vertexNDOF;}


    
    /// MET filters
    inline Bool_t   PassCSCHaloFilterTight() const {return   k_passCSCHaloFilterTight;}
    inline Bool_t   PassEcalDeadCellTriggerPrimitiveFilter() const {return   k_passEcalDeadCellTriggerPrimitiveFilter;}
    inline Bool_t   PassHBHENoiseFilter() const {return  k_passHBHENoiseFilter ;}
    inline Bool_t   PassBadEESupercrystalFilter() const {return  k_passBadEESupercrystalFilter ;}

    //Pileup reweighting
    inline Double_t PileUpInteractionsTrue() const{ return k_PileUpInteractionsTrue;}
    
    inline Double_t LHEWeight() const{return k_lheweight;}
    inline Int_t Id1() const {return k_pdf_id1;}
    inline Int_t Id2() const {return k_pdf_id2;}
    inline Double_t Q() const {return k_pdf_q;}
    inline Double_t x1() const {return k_pdf_x1;}
    inline Double_t x2() const {return k_pdf_x2;}
			
    
    virtual void Reset();    
  protected:
    /// Reset function.                                                                  
    
  private:
    /// decalre private functions
    
    Int_t    k_EventNumber, k_RunNumber,k_nvertices,  k_lumisec, k_ngoodvertices,k_pdf_id1, k_pdf_id2;
    Double_t k_vertexX,k_vertexY,k_vertexZ, k_vertexNDOF,  k_mcweight, k_lheweight, k_pdf_q, k_pdf_x1, k_pdf_x2;
    //std::vector<Double_t> k_pdf_cteq,k_pdf_mstw, k_pdf_nnpdf;
    Double_t k_PF_MET, k_PF_METphi, k_PF_SumET ;
    Double_t k_Puppi_MET, k_Puppi_METphi, k_Puppi_SumET ;
    Double_t k_NoHF_MET, k_NoHF_METphi, k_NoHF_SumET ;
    Double_t k_PfMva_MET, k_PfMva_METphi, k_PfMva_SumET ;

    Bool_t k_isData, k_isgoodevent;
    Bool_t k_passBadEESupercrystalFilter,k_passCSCHaloFilterTight,k_passEcalDeadCellTriggerPrimitiveFilter,  k_passHBHENoiseFilter;
    Double_t  k_PileUpInteractionsTrue;
    
    
    ClassDef(KEvent,16)
  }; 
  
}//namespace snu

#endif
