#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "../../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "../../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"
#include "PlotUtils/ParticleResponseDefaults.h"
using namespace NUKECC_ANA;

typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;

// Error Bands
std::map<std::string, std::vector<CVUniverse*> > GetVertErrorBands(PlotUtils::ChainWrapper* chain) {
  
  // return map
  SystMap error_bands;

  // CV 
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics) {
  
    //========================================================================
    // FLUX
    //========================================================================
    int n_flux_universes = 100 ;
    SystMap flux_systematics = 
    PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,n_flux_universes);
  
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());

    //========================================================================
    // GENIE
    //========================================================================
    SystMap genie_systematics = PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain);
    // change that true to a switch on do_nonrespi_tune
    error_bands.insert(genie_systematics.begin(), genie_systematics.end());
    // includes standard Genie + Rvx1pi (final state pion normalization) + FaCCQE

    //========================================================================
    // MnvTunes
    //======================================================================== 
  
    // RPA
    SystMap RPA_systematics = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
    error_bands.insert(RPA_systematics.begin(),RPA_systematics.end());

    // 2p2h
    SystMap _2p2h_systematics = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
    error_bands.insert(_2p2h_systematics.begin(),_2p2h_systematics.end());

  
    //========================================================================
    //  Muons
    //========================================================================

    // Muon momentum
    // MINOS
    SystMap muonminosP_systematics = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminosP_systematics.begin(), muonminosP_systematics.end());
  
    // Muon momentum
    // MINERvA
    SystMap muonminervaP_systematics = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonminervaP_systematics.begin(), muonminervaP_systematics.end());

    //Muon Momentum Resolution
    SystMap muonP_resolutions = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
    error_bands.insert(muonP_resolutions.begin(),muonP_resolutions.end());
  
    //MINOS Efficiency
    SystMap minos_efficiency = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
    error_bands.insert(minos_efficiency.begin(),minos_efficiency.end());

    //========================================================================
    // Detector
    //========================================================================
  
   // SystMap angle_systematics = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
    ///////////////////////NSFDefaults::beamThetaX_Err,NSFDefaults::beamThetaY_Err);
   // error_bands.insert(angle_systematics.begin(), angle_systematics.end());


    // Target Mass Systematics
    SystMap targetMass_systematics = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
    error_bands.insert(targetMass_systematics.begin(), targetMass_systematics.end());

    //========================================================================
    // Particle Response Systematics (Detector)
    //========================================================================

/*
    const bool NEUTRON = false;
    const bool use_new_part_resp = false;
    const bool PROTON = true;
    SystMap Response_systematics = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, NEUTRON, use_new_part_resp, PROTON);
    error_bands.insert(Response_systematics.begin(), Response_systematics.end());
*/ //was for old implementation, p3 doesn't have legacy branches needed for that

///////////


    const bool use_ID = true;
    const bool use_OD = true;
    std::string name_tag = "allNonMuonClusters";
    const bool use_neutron = false;
    const bool use_new = false;
    const bool use_proton = true;
    SystMap response_systematics = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_ID, use_OD, name_tag, use_neutron, use_new, use_proton);
    error_bands.insert(response_systematics.begin(), response_systematics.end());


    //========================================================================
    // GEANT hadrons with MnvHadronReweight
    //========================================================================
    SystMap GeantHadron_systematics = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
    error_bands.insert(GeantHadron_systematics.begin(), GeantHadron_systematics.end());
  }

  return error_bands;
}

std::map<std::string, std::vector<CVUniverse*> > GetLatErrorBands(PlotUtils::ChainWrapper* chain) {

   SystMap error_bands;
   error_bands[std::string("CV")].push_back(new CVUniverse(chain));
   if(RunCodeWithSystematics) {
      SystMap angle_systematics = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
      error_bands.insert(angle_systematics.begin(), angle_systematics.end());
   }
   return error_bands;
}

#endif  // Systematics_h
