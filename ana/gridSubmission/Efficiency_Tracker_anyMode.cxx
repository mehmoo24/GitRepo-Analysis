//==============================================================================
// Loop entries, make cuts, fill histograms.
// * Uses the New Systematics Framework and "Universe" objects.
// * loop universes, make cuts and fill histograms with the correct lateral
// shifts and weights for each universe.
// * TChain --> PlotUtils::ChainWrapper.
// * MnvHXD --> PlotUtils::HistWrapper.
// * Genie, flux, non-resonant pion, and some detector systematics calculated.
//==============================================================================

//#include "include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "../include/Variable_3rings.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
//#include "/exp/minerva/app/users/mmehmood/MAT_AL9/NSFNukeCCInclusive/ana/include/systematics/Systematics.h"
#include "/exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/include/systematics/Systematics.h"
//#include "../../NUKECCSRC/ana_common/include/LateralSystematics.h"


#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "TParameter.h"



#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/ParticleResponseDefaults.h"

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


void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A", bool doDIS=true);
 
int inputpdg;

int main(int argc, char *argv[]){
//   ROOT::Cintex::Cintex::Enable();
   TH1::AddDirectory(false);

 if(argc==1){
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     std::cout<<"MACROS HELP:\n\n"<<
       "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<< std::endl;
     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
     return 0;
   }

   TString dir(argv[1]);
   int targetID = 99;
   int targetZ = 99;
//   int targetID = atoi(argv[2]);
//   int targetZ = atoi(argv[3]);

//   TString option(argv[2]);
//  const std::string mc_file_list("mc"+option+".txt");
//  const std::string data_file_list("data"+option+".txt");

   const string playlist= argv[2];

    const std::string plist_string(playlist);
    inputpdg = atoi(argv[3]);  // 14 for neutrino mode (wanting to look at neutrinos in the beam), -14 for anti neutrino mode

    const string name_of_file = argv[4];  
//   const string fewer_files_opt = argv[5];
   const string fewer_files_opt = "";

    const std::string mc_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_mc_%s%s.txt", plist_string.c_str(), fewer_files_opt.c_str()));
    const std::string data_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_data_%s%s.txt",plist_string.c_str(), fewer_files_opt.c_str()));

	
     bool doDIS=false;
//  const std::string plist_string("minervame"+option);
//  const std::string plist_string("minervame6A");
  const std::string reco_tree_name("MasterAnaDev");
  const bool wants_truth = true;
  const bool is_grid = false;


  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, plist_string, wants_truth);
  PlotUtils::MacroUtil util_data(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);

    util.PrintMacroConfiguration("main");

  //=========================================
  // Systematics
  //=========================================
 // std::map<std::string, std::vector<CVUniverse*> > error_bands =
   //   GetErrorBands(util.m_mc);

   PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
   PlotUtils::MinervaUniverse::RPAMaterials(false);
   PlotUtils::MinervaUniverse::SetNuEConstraint(true);
   //PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
   PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
   //PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
   PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
   PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
   // MnvHadronReweighter (GEANT Hadron sytematics)
   PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
   //Neutron CV reweight is on by default (recommended you keep this on) 
   PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight( true );
   //Elastics are on by default (recommended you keep this on)  
   PlotUtils::MinervaUniverse::SetMHRWeightElastics( true );
    
   NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
   NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
   NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
    PlotUtils::ChainWrapper* chainTruth = util.m_truth;
    PlotUtils::ChainWrapper* chainMC = util.m_mc;
//    HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

  HelicityType::t_HelicityType helicity;
  if(inputpdg == 14) helicity = NUKECC_ANA::HelicityType::kNeutrino;
  else helicity = NUKECC_ANA::HelicityType::kAntiNeutrino;

    cout<<"Helicity"<<helicity<<endl;

    double DataPot=  util_data.m_data_pot; 
    double MCPot=  util.m_mc_pot; 
 //  double total_pot_data,total_pot_mc;
  // utils->getPOT(total_pot_data,total_pot_mc);  
   double  MCscale=DataPot/MCPot;
  // double  MCscale=1.0;
 
      std::cout<<" MCScale= "<<MCscale<<std::endl; 
   std::vector<Var*> variablesMC,variablesTruth; 
   std::vector<Var2D*> variables2DMC,variables2DTruth; 

   TString histFileName;
   if(RunCodeWithSystematics){
     histFileName = utils->GetHistFileName(Form("Efficiency_%s_%s_%ssys", name_of_file.c_str(),playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );
   }
   else
     histFileName = utils->GetHistFileName(Form("Efficiency_%s_%s_%snosys",name_of_file.c_str(), playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );

   //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );
   
   TFile fout(dir.Append(histFileName),"RECREATE");	
   
   // For 1D variables 
   FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true, targetID, targetZ, plist_string,doDIS);
       
   for (auto v : variablesMC) v->m_selected_mc_reco.SyncCVHistos();
   for (auto v : variables2DMC){
     v->m_selected_mc_reco.SyncCVHistos();
     v->reco_QE.SyncCVHistos();
     v->reco_2p2h.SyncCVHistos();
     v->reco_trueDIS.SyncCVHistos();
     v->reco_softDIS.SyncCVHistos();
     v->reco_other.SyncCVHistos();
     v->reco_resonant.SyncCVHistos();     
     for(int ring=0;ring<3;ring++){
        v->ring_m_selected_mc_reco[ring].SyncCVHistos();
        v->ring_reco_QE[ring].SyncCVHistos();
        v->ring_reco_2p2h[ring].SyncCVHistos();
        v->ring_reco_trueDIS[ring].SyncCVHistos();
        v->ring_reco_softDIS[ring].SyncCVHistos();
        v->ring_reco_other[ring].SyncCVHistos();
        v->ring_reco_resonant[ring].SyncCVHistos();
     }
   }
 
   FillVariable(chainTruth, helicity, utils, cutter,binsDef,variablesTruth,variables2DTruth,false, targetID, targetZ, plist_string,doDIS);
   
   for (auto v : variables2DTruth){
     v->m_selected_truth_reco.SyncCVHistos();
     v->true_QE.SyncCVHistos();
     v->true_2p2h.SyncCVHistos();
     v->true_trueDIS.SyncCVHistos();
     v->true_softDIS.SyncCVHistos();
     v->true_other.SyncCVHistos();
     v->true_resonant.SyncCVHistos(); 
     for(int ring=0;ring<3;ring++){
        v->ring_m_selected_truth_reco[ring].SyncCVHistos();
        v->ring_true_QE[ring].SyncCVHistos();
        v->ring_true_2p2h[ring].SyncCVHistos();
        v->ring_true_trueDIS[ring].SyncCVHistos();
        v->ring_true_softDIS[ring].SyncCVHistos();
        v->ring_true_other[ring].SyncCVHistos();
        v->ring_true_resonant[ring].SyncCVHistos();
     }

   }

 //  for (auto v : variablesTruth) v->data_reco_sb.SyncCVHistos();
   for (auto v : variablesTruth) v->m_selected_truth_reco.SyncCVHistos();
 
   for (auto v : variablesMC) {
     v->WriteAllHistogramsToFileEff(fout, true);
   }


   for (auto v : variablesTruth) {
     v->WriteAllHistogramsToFileEff(fout, false);
   }

   // Plotting If you want for 1D
   /*   for(int i=0; i< variablesMC.size();i++){
     PlotCVAndError(variablesData[i]->data_reco.hist,variablesMC[i]->mc_reco.hist,variablesMC[i]->GetName(),MCscale);
       
     PlotErrorSummary(variablesMC[i]->mc_reco.hist, variablesMC[i]->GetName());
     PlotStacked(variablesData[i]->data_reco_sb.hist,variablesMC[i]->mc_sb.GetHistArray(),MCscale, variablesMC[i]->mc_sb.GetName(), variablesMC[i]->mc_sb.GetName());
   }//End 1D plotting 
   */
   for (auto v : variables2DMC) {
     v->WriteAllHistogramsToFileEff(fout,true);
   }

   for (auto v : variables2DTruth) {
     v->WriteAllHistogramsToFileEff(fout,false);
   }

  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 


 
   //Plotting in 2D
   
   for(int i=0; i< variables2DMC.size();i++){
     Plot2D(variables2DMC[i]->m_selected_mc_reco.hist, variables2DMC[i]->GetName(), variables2DMC[i]->GetNameX(), variables2DMC[i]->GetNameY()); //Plotting line that I somehow cannot delete without producing memory errors, but no one else can reproduce. --ANF 2020.4.6
     //Plot2D(variables2DData[i]->data_reco.hist, variables2DData[i]->GetName(), variables2DData[i]->GetNameX(),variables2DData[i]->GetNameY());
     
   }//End 2D plotting

}//End Main

   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID, int targetZ, const string playlist, bool doDIS){
 // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain);
    
  std::map<std::string, std::vector<CVUniverse*> > error_bands =GetErrorBands(chain);
   std::vector<double> Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,dansPTBins,dansPZBins;
   if (doDIS){
     Enubin = binsDef->GetDISBins("Enu"); 
     Emubin = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     Q2bin = binsDef->GetDISBins("Q2");
     Wbin = binsDef->GetDISBins("W");
     xbin = binsDef->GetDISBins("x");
     ybin = binsDef->GetDISBins("y");
     }
   else{
     Enubin = binsDef->GetEnergyBins("Enu"); 
     Emubin = binsDef->GetEnergyBins("Emu"); 
     Ehadbin = binsDef->GetEnergyBins("Ehad");
     Q2bin = binsDef->GetEnergyBins("Q2");
     Wbin = binsDef->GetEnergyBins("W");
     xbin = binsDef->GetEnergyBins("x");
     ybin = binsDef->GetEnergyBins("y");
    dansPTBins = binsDef->GetEnergyBins("pTmu");
    dansPZBins = binsDef->GetEnergyBins("pZmu");
   }
   //Q2bin = binsDef->GetSidebandBins("Q2");
   //Wbin = binsDef->GetSidebandBins("W");

   Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
   Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
   Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
   Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
   Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
   Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
   Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
  Var* pTmu = new Var("pTmu", "pTmu", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue);
  Var* pZmu = new Var("pZmu", "pZmu", dansPZBins, &CVUniverse::GetMuonPZ, &CVUniverse::GetMuonPZTrue);

   //std::vector<Var*> variables = {enu,ehad}; 
//   variables = {emu};//{enu,ehad}; 
  variables = {};
   
   Var2D* W_Q2 = new Var2D(*W, *Q2);
   Var2D* x_Q2 = new Var2D(*x, *Q2);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_enu = new Var2D(*emu, *enu);
   Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
   Var2D* x_y = new Var2D(*x, *y);  // y var
   Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu);
   Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);

//   variables2d = {emu_ehad, W_Q2, x_Q2, x_y};//{enu_ehad, Q2_W};
   variables2d = {pZmu_pTmu}; 
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
   for (auto v : variables) v->InitializeAllHistograms(error_bands);

   int mc0=0;
   int mc1=0;
   int mc2=0;
   int mc3=0;
   int mc4=0; 
   int mc5=0; 
   int mc6=0; 
   int mc7=0; 
   int mc8=0; 
   int mc9=0; 
   int mc10=0; 
   int mc11=0; 
   int mc12=0; 
   int mc13=0; 
   int mc14=0; 
   int mc15=0; 
   int mc16=0; 
   int mc17=0; 
     int mc_truth0=0.0;
     int mc_truth1=0.0;
     int mc_truth2=0.0; 
     int mc_truth3=0.0; 
     int mc_truth4=0.0; 
     int mc_truth5=0.0; 
     int mc_truth6=0.0; 
     int mc_truth7=0.0; 
   int allcuts=0;
   
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
     for (auto band : error_bands){
       std::vector<CVUniverse*> error_band_universes = band.second;
       for (auto universe : error_band_universes){
	 // Tell the Event which entry in the TChain it's looking at
	 universe->SetEntry(i);
	 //=========================================
	 // CUTS in each universe
	 //=========================================
	 
	if(isMC){ 
/////////////////////////////////////////////////////////////////
//              universe->SetEntry(i);
              if ((universe->GetInt("has_interaction_vertex") != 1)) continue;
  
              if(!cutter->PassReco(universe,helicity)) continue;

              if(!cutter->TrackerOnly(universe)) continue;

              if(!cutter->InHexagon(universe, 850)) continue;

              if(!cutter->PassThetaCut(universe)) continue;
              int ring = universe->GetRingBeamCoord();
	
	   for (auto v : variables2d){
                    if ((universe->GetInt("mc_incoming") == inputpdg) && (universe->GetInt("mc_current") == 1)){
             v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight()); 
             v->ring_m_selected_mc_reco[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
            
            // mc breakdown
/*
            if (universe->GetInt("mc_intType") == 1){
                v->reco_QE.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                v->ring_reco_QE[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
            }
            else if (universe->GetInt("mc_intType") == 8){
                v->reco_2p2h.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                v->ring_reco_2p2h[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
            }
            else if (universe->GetInt("mc_intType") == 3){
              if (cutter->PassTrueDISCut(universe)){
                  v->reco_trueDIS.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                  v->ring_reco_trueDIS[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
              }
              else {
                  v->reco_softDIS.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                  v->ring_reco_softDIS[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
              }
            }
            else if (universe->GetInt("mc_intType") == 2){
                  v->reco_resonant.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                  v->ring_reco_resonant[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
            }
            else {
                  v->reco_other.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
                  v->ring_reco_other[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());
            }
*/
           }
            } //end variables2d for loop
	   for (auto v : variables){
	            if ((universe->GetInt("mc_incoming") == inputpdg) && (universe->GetInt("mc_current") == 1)){ 
              v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
           }
	     //v->m_selected_mc_sb.GetComponentHist("MC")->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
	   }
	 } else if(!isMC ){

         if(!cutter->TrackerOnlyTrue(universe)) continue;

         if(!cutter->InHexagonTrue(universe, 850)) continue;

         if(!cutter->PassTrueThetaCut(universe)) continue;
 
	 if( 1 != universe->GetInt("mc_current") ) continue;

         if( inputpdg != universe->GetInt("mc_incoming") ) continue;

         int ring = universe->GetTrueRingBeamCoord();
	   for (auto v : variables2d){
	     v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
             v->ring_m_selected_truth_reco[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());

             if (universe->GetInt("mc_intType") == 1){ 
                v->true_QE.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                v->ring_true_QE[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
             } 
             else if (universe->GetInt("mc_intType") == 8){
                v->true_2p2h.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                v->ring_true_2p2h[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
             }
             else if (universe->GetInt("mc_intType") == 3){
               if (cutter->PassTrueDISCut(universe)){
                  v->true_trueDIS.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                  v->ring_true_trueDIS[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
               }
               else {
                  v->true_softDIS.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                  v->ring_true_softDIS[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
               }
             }
             else if (universe->GetInt("mc_intType") == 2){
                  v->true_resonant.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                  v->ring_true_resonant[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
             }
             else {
                  v->true_other.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
                  v->ring_true_other[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetTruthWeight());
             }
      
	   } // end variables2d for loop
	   for (auto v : variables){
	       v->m_selected_truth_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetTruthWeight());
	   }
	 }
       } // End band's universe loop
     }// End Band loop
   }//End entries loop


}
//============================================================================================================================
// Main
