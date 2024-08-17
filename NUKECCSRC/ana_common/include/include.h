//========================================================================
// Base class for an un-systematically shifted (i.e. CV) universe.  Implement
// your own base class with the functions you need. I've implemented GetEnu(),
// GetMuonE() and GetHadronE() as examples: you'll have other variables you
// want.
//
// To add a systematic, inherit from this class and override whichever
// functions you need to. For a "vertical" error, this will mean overriding the
// GetWeight() function to modify the event weight. For a "lateral" error, this
// will mean overriding the function that calculates the quantity that is being
// shifted (muon energy, or hadronic energy or whatever).
//
// For examples of each of those two cases, see ./LateralSystematics.h and
// PlotUtils/GenieSystematics.h. For an example of how to put the whole thing
// together and actually *use* the classes, see the runEventLoop.C macro in
// this directory. `root -l -b load.C+ runEventLoop.C+`
// ========================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

//#include "NukeCCsetbranchaddrs.h"
#include "CommonIncludes.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/CaloCorrection.h"
//#include "PlotUtils/RecoilEnergyFunctions.h" // for new particle response systematic

#include <iostream>
//#include "CCQENuUtilsNSF.h"
//using namespace globalV;

//#include "NukeCCvars.h"
namespace NUKECC_ANA
{
class CVUniverse : public PlotUtils::MinervaUniverse {
  public:
    #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
    #include "PlotUtils/TruthFunctions.h" //Getq3True
    #include "PlotUtils/WeightFunctions.h"
    #include "PlotUtils/RecoilEnergyFunctions.h"
    // Constructor
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma=0);
    CVUniverse();

    // Destructor
    virtual ~CVUniverse();
    // All functions we write here should be 'virtual', so that the universes
    // that inherit from CVUniverse can override them.

    // ========================================================================
    // Get Weight
    //
    // We override the various weight getting functions herein in different
    // vertical systematic universe classes.
    // ========================================================================
  

  virtual double GetWeight() const;
  virtual double GetWeightEmu() const;
  virtual double GetWeightQ2() const;
  virtual double GetWeighty() const;
  virtual double GetTruthWeight() const;
        virtual Long64_t GetMaxEntries(); 
    virtual double GetMuonE() const { return GetVecElem("MasterAnaDev_leptonE",3); }

//      virtual void  Show(Long64_t entry = -1);


virtual double calcq3(const double Q2,const double Enu,const double Emu	)const;
virtual double calcq0(const double Enu,const double Emu	)const;
      //! Turn off branches that I don't use
      //virtual bool SetBranches( bool useWeights = true );

virtual double GetMuonCurve()  const;
virtual double GetFiducial()  const;
virtual double GetTdead()  const;
virtual double GetTargetID()  const;
virtual double GetHelicity()  const;
virtual double GetTrueHelicity()  const;

//      virtual double GetRecoilEnergy() const; //commented out for new particle response systematic
 virtual double GetCalRecoilEnergy() const; // added for new particle response sys.
 virtual double GetNonCalRecoilEnergy() const; //added for new particle response sys.
 virtual double ApplyCaloTuning(double calRecoilE) const;

      virtual double GetEhadGeV() const;
      virtual double GetMuonEGeV() const;
      virtual double GetMuonETrueGeV() const;
 virtual double GetThetamuTrueDeg()  const ;
// virtual double GetThetamu()  const ;
 virtual double GetMuonPt()  const ;
 virtual double GetMuonPz()  const ;
 virtual double GetlepPtTrue()  const ;
 virtual double GetQ2Reco()     const ;
 virtual double GetWReco()     const;
 virtual double GetxReco()     const;
 virtual double GetyReco()     const;

// from Anezka, added May 13 for production work
virtual int GetMCRunN() const;
virtual double GetlepPzTrue()  const ;
double virtual GetVertexZMy() const;
double virtual GetVertexZMyML() const;
double virtual GetVertexZMyTBV() const;
double virtual GetVertexZTrueMy() const;
virtual double GetplaneDNNTrue()   const;

// end of added May 13 

virtual double GetVertexTBVminusML() const;
virtual double GetVertexTBVminusTrue() const;
virtual double GetVertexMLminusTrue() const;


 virtual double GetQ2RecoGeV()     const ;
 virtual double GetWRecoGeV()     const;
virtual double GetEnu() const;  
virtual double GetEnuGeV() const;  
virtual double GetEnuTrueGeV() const;  

virtual double GetMuonPT() const;
virtual double GetMuonPTTrue() const;
virtual double GetMuonPZ() const;
virtual double GetMuonPZTrue() const;
//virtual double GetVertexZ() const;
//virtual double GetVertexZTrue() const;


// for true variables
virtual double GetTruthNuPDG() const;
virtual double GetCurrent()    const; 
virtual double GetEhadTrue() const;
virtual double GetEhadTrueGeV() const;
virtual double GetThetamuTrue( )const;

 virtual double GetQ2IncTrue()     const ;
 virtual double GetWTrue()     const;
 virtual double GetyTrue()     const;
 virtual double GetxTrue()     const;
 
 virtual double GetQ2TrueGeV()     const ;
 virtual double GetWTrueGeV()     const;
virtual double GetThetamuDeg()  const;
//virtual double Getq3True()const; 
//virtual double Getq0True()const; 
//virtual double GetxVertexTrue()const; //added June 17/2021
//virtual double GetxVertexReco()const; //added June 17/2021
   
virtual double calcRecoQ2(const double Enu,const double Emu,const double Thetamu)const;	
virtual double calcWReco(const double Q2,const double Ehad) const;

virtual double calcTrueQ2(const double EnuTrue,const double EmuTrue,const double ThetamuTrue)const;	
virtual double calcWTrue(const double Q2True,const double EhadTrue) const;


virtual double calcXTrue(const double Q2True,const double EnuTrue, const double EmuTrue)const;

virtual double calcXReco(const double Q2,const double Enu, const double Emu)
const ;


virtual double calcYTrue(const double EnuTrue,const double EhadTrue)const;	

virtual double calcYReco(const double Enu,const double Ehad)const;	

       
       virtual  int GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane );
       virtual double GetVtxEnergy(  unsigned int shift = VTX_BLOB_SHIFT );  
	
       virtual double GetCCQERecoil( unsigned int shift = VTX_BLOB_SHIFT );
	
       virtual double GetMuonP() const;
	

};
}

#endif
