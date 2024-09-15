//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

#include "../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/VariableEvent.h"
//#include "../include/Variable.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/FluxSystematics.h"
#include <iostream>
#include <stdlib.h>
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "TParameter.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we 
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h" 
//using namespace globalV;
using namespace NUKECC_ANA;

//======================================================================

typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;

bool isMC = false;
//void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
//helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
//std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
//int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

//=============================================================================
//=============================================================================
// MAIN FUNCTION
//=============================================================================
//=============================================================================

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

int inputpdg;

int main(int argc, char *argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "MACROS HELP:\n\n" <<
    "\t-./runEventLoop Path_to_Output_file Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
    "\t-Playlist\t \t =\t eg. minervame1A  \n" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }

  TString dir(argv[1]);
  const string playlist= argv[2];
  int targetID = 99;
  int targetZ = 99;
   
  bool doDIS=false;
  //bool RunCodeWithSystematics = false;
  //bool RunCodeWithSysFluxUniverses = false;

  const std::string plist_string(playlist);

  inputpdg = atoi(argv[3]);  // 14 for neutrino mode (wanting to look at neutrinos in the beam), -14 for anti neutrino mode

  const string name_of_file = argv[4];
//  const string fewer_files_opt = argv[5];
   const string fewer_files_opt = "";  

  //const std::string mc_file_list(Form("../include/playlists/DansMuonKludged/MasterAnaDev_MC_%s_MuonKludged%s.txt", plist_string.c_str(),fewer_files_opt.c_str()));
  const std::string mc_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_mc_%s%s.txt", plist_string.c_str(), fewer_files_opt.c_str()));

  const std::string reco_tree_name("MasterAnaDev");
  
  const bool wants_truth = true; // read in truth!

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, plist_string, wants_truth);

  util.PrintMacroConfiguration("main");

  // SYSTEMATICS

  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);

  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  //PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  //PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);  
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);


  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
  NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
  //PlotUtils::ChainWrapper* chainData = util.m_data;
  //PlotUtils::ChainWrapper* chainMC = util.m_mc;
  PlotUtils::ChainWrapper* chainTruth = util.m_truth; // reading in the Truth tree
  //HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

  HelicityType::t_HelicityType helicity;
  if(inputpdg == 14) helicity = NUKECC_ANA::HelicityType::kNeutrino;
  else helicity = NUKECC_ANA::HelicityType::kAntiNeutrino;

  //double DataPot=  util.m_data_pot; 
  double MCPot=  util.m_mc_pot;  
  //double MCscale=DataPot/MCPot;

  //std::cout << "MC Scale = " << MCscale << std::endl; 
  //std::cout << "Data POT: " << DataPot << std::endl;
  std::cout << "MC POT: " << MCPot << std::endl;

  std::vector<Var*> variablesMC;//,variablesData; 
  std::vector<Var2D*> variables2DMC;//,variables2DData; 

  TString histFileName;
  if(RunCodeWithSystematics){
    histFileName = utils->GetHistFileName( Form("EventSelection_BeamCoord_%s_%s_%ssys", name_of_file.c_str(),playlist.c_str(),fewer_files_opt.c_str()),  FileType::kAny, targetID, targetZ, helicity );
  }

  else{
    histFileName = utils->GetHistFileName( Form("EventSelection_BeamCoord_%s_%s_%snosys", name_of_file.c_str(),playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );
  }

  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
  
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainTruth, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,false,targetID, targetZ, plist_string,doDIS);
    
  for (auto v : variablesMC) {
    v->m_selected_truth_reco.SyncCVHistos();
    }


  for (auto v : variables2DMC) {
    v->m_selected_truth_reco.SyncCVHistos();
  }

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFileEff(fout, false);
  }

  // 2D Variables
  for (auto v : variables2DMC) {
    v->WriteAllHistogramsToFileEff(fout,false);
  }
  
  //Writing POT to the HistFile
  fout.cd();
  //auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  //dataPOTOut->Write();
  mcPOTOut->Write(); 
  fout.Close();

  std::cout << "DONE " << plist_string << std::endl;


}//End Main


//=============================================================================
//=============================================================================
// OTHER FUNCTIONS
//=============================================================================
//=============================================================================

std::map<std::string, std::vector<CVUniverse*> > GetErrorBands(PlotUtils::ChainWrapper* chain) {
  typedef std::map<std::string, std::vector<CVUniverse*> > SystMap;
  SystMap error_bands;
  // CV
  error_bands[std::string("CV")].push_back(new CVUniverse(chain));

  if(RunCodeWithSystematics){
    SystMap flux_systematics = PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain,CVUniverse::GetNFluxUniverses());
    error_bands.insert(flux_systematics.begin(), flux_systematics.end());
  }
  return error_bands;
}
// Based on Jeffrey Kleykamp's method
// Do not need other systematics, just the flux
// measuring the number of neutrinos which doesn't change
//  as a function of interaction model since you're dividing by the xsec in the end

// Group universes that change by vertical shift == flux universes
// Meaning: if 500 flux universes, can just apply one cut rather than run the vent 500 times
/*
std::vector<std::vector<CVUniverse*>> groupCompatibleUniverses(const std::map<std::string, std::vector<CVUniverse*>> bands) {
  std::vector<std::vector<CVUniverse*>> groupedUnivs;
  std::vector<CVUniverse*> vertical;
  for(const auto& band: bands){
    if(band.first == "cv") vertical.insert(vertical.begin(), band.second.begin(), band.second.end());
    else{
      for(const auto univ: band.second){
        if(univ->IsVerticalOnly()) vertical.push_back(univ);
        else groupedUnivs.push_back(std::vector<CVUniverse*>{univ});
      }
    }
  }
  groupedUnivs.insert(groupedUnivs.begin(), vertical);
  return groupedUnivs;
}
*/

// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  for(auto band :error_bands ){std::cout<<"Checking Universe this is universe with name : " << band.first<<std::endl;}
  //std::vector<std::vector<CVUniverse*>> error_bands_GROUP = groupCompatibleUniverses(error_bands);
  std::cout<< "error_bands.size() = " << error_bands.size()<<std::endl;
  //std::cout<< "Error_Band_GROUPS.size() = " << error_bands_GROUP.size()<<std::endl;
  std::cout<<"Number of Universes set is = "<<    MinervaUniverse::GetNFluxUniverses()<<std::endl;
  
  std::vector<double> Enubin, ringbin;
  Enubin = binsDef->GetEnergyBins("EnuFlux"); // 0.1 GeV binning from 0 to 100 GeV
  ringbin = binsDef->GetEnergyBins("ring");

  // 1D Variables
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* ring = new Var("Ring", "Ring", ringbin, &CVUniverse::GetTrueRingBeamCoord, &CVUniverse::GetTrueRingBeamCoord);

  variables = {enu, ring};//, vtxz, vtxx, vtxy};

  Var2D* enu_ring = new Var2D(*enu, *ring);
  Var2D* ring_enu = new Var2D(*ring, *enu);
  variables2d = {enu_ring, ring_enu};

  for (auto v : variables) v->InitializeAllHistograms(error_bands);
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);


  int all =0 ;
  int truth = 0 ;
  int exclude2p2h = 0;
  int trackerC = 0;
    
  //CVUniverse *dataverse = new CVUniverse(chain,0);

  //=========================================
  // Entry Loop
  //=========================================

  std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;
  for(int i=0; i<chain->GetEntries(); ++i){
    if(i%500000==0) std::cout << (i/1000) << "k " << std::endl;
      //=========================================
      // For every systematic, loop over the universes, and fill the
      // appropriate histogram in the MnvH1D
      //=========================================
      if(!isMC){     
    //    for (auto band_GROUP : error_bands_GROUP){
          for (auto band : error_bands){

          // Tell the Event which entry in the TChain it's looking at
          //=========================================
          // CUTS in each universe
          //========================================
    //      band_GROUP.front()->SetEntry(i);
          std::vector<CVUniverse*> error_band_universes = band.second;
          for (auto universe : error_band_universes){
          universe->SetEntry(i);
          all++;

          if( ! cutter->PassTrueCC(universe,helicity)) continue;
          if( ! cutter->InHexagonTrue(universe, 850.) ) continue;
          truth++;

          // exclude if true 2p2h
          if (8 == universe->GetInt("mc_intType") ) continue;
          exclude2p2h++;

          // TRACKER
          if (cutter->TrackerOnlyTrue(universe)){
            //5890.0, 8467.0
              trackerC++;
           //   for (auto universe : band_GROUP){
                for (auto v : variables) v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeightFlux());
                for (auto v : variables2d) v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeightFlux());
           //   }
            }    
          }
        }// End Band loop
      
      
}
      else{
        std::cout << "Not MC!!" << std::endl;
      }       
    }//End entries loop



  for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 

 // int univ_norm = error_bands_GROUP.size();
}
//=============================================================================

