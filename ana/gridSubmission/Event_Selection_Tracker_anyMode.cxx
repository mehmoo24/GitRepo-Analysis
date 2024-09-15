// ./Event_Selection_Tracker . minervame1A short 1run

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
#include "../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
//#include "../include/VariableRun.h" //trying to run from /minerva/data
//#include "/minerva/app/users/mmehmood/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/include/Variable.h" //Zubair confirmed using Variable.h instead of VariableRun.h (July 11/21)
#include "../include/Variable_3rings.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
//#include "../../NUKECCSRC/ana_common/include/LateralSystematics.h"
#include <iostream>
#include <stdlib.h>
//#include "Cintex/Cintex.h"
#include "/exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/include/systematics/Systematics.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "PlotUtils/MnvPlotter.h" 
#include "TParameter.h"

//#include "PlotUtils/MinosEfficiencySystematics.h"
//#include "PlotUtils/MnvHadronReweight.h" 
//#include "PlotUtils/FluxReweighter.h"
//#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/AngleSystematics.h"
//#include "PlotUtils/MuonSystematics.h"
//#include "PlotUtils/MuonResolutionSystematics.h"
//#include "PlotUtils/MnvTuneSystematics.h"
//#include "PlotUtils/TargetMassSystematics.h"
//#include "PlotUtils/GeantHadronSystematics.h"
//#include "PlotUtils/ResponseSystematics.h"
//#include "PlotUtils/ParticleResponseDefaults.h"

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

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType
helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef,
std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1,
int targetZ=26, const string playlist="minervame1A", bool doDIS=true);

//=============================================================================
//=============================================================================
// MAIN FUNCTION
//=============================================================================
//=============================================================================

// global variables
int data0;
int data1;
int data2;
int data3;
int data4;
int data5;

int mcreco0;
int mcreco1;
int mcreco2;
int mcreco3;
int mcreco4;
int mcreco5;

int mcsignal;
int mcbkg;

int num_2Dvars;

int inputpdg;

int main(int argc, char *argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  if(argc==1){
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "MACROS HELP:\n\n" <<
    "\t-./runEventLoop Path_to_Output_file Target_number Material_atomic_number Playlist\n\n" <<
    "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<<
 //   "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n" <<
//    "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n"
  std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    return 0;
  }

  TString dir(argv[1]);
//  int targetID = atoi(argv[2]);
//  int targetZ = atoi(argv[3]);


  int targetID = 99; int targetZ = 99;  
 
  bool doDIS=false;

    const string playlist= argv[2];
    const std::string plist_string(playlist);

    inputpdg = atoi(argv[3]);  // 14 for neutrino mode (wanting to look at neutrinos in the beam), -14 for anti neutrino mode

    const string name_of_file = argv[4]; // name of file (ideally include the date)
//    const string fewer_files_opt = argv[5];
    const string fewer_files_opt = "";

    const std::string mc_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_mc_%s%s.txt", plist_string.c_str(), fewer_files_opt.c_str()));
    const std::string data_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_data_%s%s.txt",plist_string.c_str(), fewer_files_opt.c_str()));

  const std::string reco_tree_name("MasterAnaDev");
  const bool wants_truth = false;
  const bool is_grid = false;

//  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth, is_grid);

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth); //is_grid option removed when did git pull MAT May 9/22

  util.PrintMacroConfiguration("main");

  // SYSTEMATICS

  //std::map<std::string, std::vector<CVUniverse*> > error_bands =
  //GetErrorBands(util.m_mc);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::RPAMaterials(false);
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
//  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(pdg);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(true); 
//  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false); 
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
  // MnvHadronReweighter (GEANT Hadron sytematics)
  PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
 // PlotUtils::MinervaUniverse::SetReadoutVolume(5980.0,8422.0,850);
  //Neutron CV reweight is on by default (recommended you keep this on)
  PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
  //Elastics are on by default (recommended you keep this on) 
  PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
  NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
  NukeCC_Binning  *binsDef = new NukeCC_Binning();
  
  PlotUtils::ChainWrapper* chainData = util.m_data;
  PlotUtils::ChainWrapper* chainMC = util.m_mc;
  //HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

  HelicityType::t_HelicityType helicity; 
  if(inputpdg == 14) helicity = NUKECC_ANA::HelicityType::kNeutrino;
  else helicity = NUKECC_ANA::HelicityType::kAntiNeutrino; 

  double DataPot=  util.m_data_pot; 
  double MCPot=  util.m_mc_pot;  
  double MCscale=DataPot/MCPot;
 
  std::cout << "MC Scale = " << MCscale << std::endl; 
  std::cout << "Data POT: " << DataPot << std::endl;
  std::cout << "MC POT: " << MCPot << std::endl;

  std::vector<Var*> variablesMC,variablesData; 
  std::vector<Var2D*> variables2DMC,variables2DData; 

  TString histFileName; 
  if(RunCodeWithSystematics){
    histFileName = utils->GetHistFileName( Form("EventSelection_%s_%s_%ssys", name_of_file.c_str(),playlist.c_str(),fewer_files_opt.c_str()),  FileType::kAny, targetID, targetZ, helicity );
  }

  else{
    histFileName = utils->GetHistFileName( Form("EventSelection_%s_%s_%snosys", name_of_file.c_str(),playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );
  }  


  //Works good for the grid submission
  //TString histFileName = utils->GetHistFileName( "EventSelection", FileType::kAny, targetID, targetZ );

  TFile fout(dir.Append(histFileName),"RECREATE");	
   
  // MC 
  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);     
  for (auto v : variablesMC) {v->m_selected_mc_reco.SyncCVHistos();
                                v->signal_purityNum.SyncCVHistos();
                                v->bkg_total.SyncCVHistos();
  }
  for (auto v : variables2DMC) {v->m_selected_mc_reco.SyncCVHistos();
                                v->signal_purityNum.SyncCVHistos();
                                v->bkg_total.SyncCVHistos();

                                v->reco_QE.SyncCVHistos();
                                v->reco_2p2h.SyncCVHistos();
                                v->reco_trueDIS.SyncCVHistos();
                                v->reco_softDIS.SyncCVHistos();
                                v->reco_other.SyncCVHistos();
                                v->reco_resonant.SyncCVHistos();

				v->ANCC.SyncCVHistos(); v->ANNC.SyncCVHistos(); v->NCC.SyncCVHistos(); v->NNC.SyncCVHistos(); v->ANNCCNCother.SyncCVHistos();

				for(int ring=0;ring<3;ring++){
                                   v->ring_m_selected_mc_reco[ring].SyncCVHistos();
                                   v->ring_signal_purityNum[ring].SyncCVHistos();
                                   v->ring_bkg_total[ring].SyncCVHistos();
                                }
  }
   
  // DATA
  std::cout << "Processing Data and filling histograms" << std::endl;

  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();
  for (auto v : variablesData) v->m_selected_data_reco_sb.SyncCVHistos();
  for (auto v : variables2DData) {
     v->m_selected_data_reco.SyncCVHistos();
     for(int ring=0;ring<3;ring++){
        v->ring_m_selected_data_reco[ring].SyncCVHistos();
     }
  }
  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  for (auto v : variablesData) {
    v->WriteAllHistogramsToFile(fout, false);
  }

 
  // 2D Variables
  for (auto v : variables2DMC) {
    v->WriteAllHistogramsToFile(fout,true);
  }

  for (auto v : variables2DData) {
    v->WriteAllHistogramsToFile(fout,false);
  }
 
	
  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 

  auto data0Cut = new TParameter<int>("DataNoCuts", data0);
  auto data1Cut = new TParameter<int>("DataRecoCuts", data1);
  auto data2Cut = new TParameter<int>("DataTrackerOnly", data2);
  auto data3Cut = new TParameter<int>("DataInHexagon", data3);
  auto data4Cut = new TParameter<int>("DataPassThetaCut", data4);

  auto mcreco0Cut = new TParameter<int>("MCRecoNoCuts", mcreco0);
  auto mcreco1Cut = new TParameter<int>("MCReco_RecoCuts", mcreco1);
  auto mcreco2Cut = new TParameter<int>("MCRecoTrackerOnly", mcreco2);
  auto mcreco3Cut = new TParameter<int>("MCRecoInHexagon", mcreco3);
  auto mcreco4Cut = new TParameter<int>("MCRecoPassThetaCut", mcreco4);

  mcsignal = mcsignal/num_2Dvars;
  mcbkg = mcbkg/num_2Dvars;
  auto mcsignalCut = new TParameter<int>("MCSignal", mcsignal);
  auto mcbkgCut = new TParameter<int>("MCBkg", mcbkg);
 

  data0Cut->Write();
  data1Cut->Write();
  data2Cut->Write();
  data3Cut->Write();
  data4Cut->Write();
//  data5Cut->Write();
  mcreco0Cut->Write();
  mcreco1Cut->Write();
  mcreco2Cut->Write();
  mcreco3Cut->Write();
  mcreco4Cut->Write();
//  mcreco5Cut->Write();
  mcsignalCut->Write();
  mcbkgCut->Write();

  std::cout << "This is data0!!: " << data0 << std::endl;


}//End Main


//=============================================================================
//=============================================================================
// OTHER FUNCTIONS
//=============================================================================
//=============================================================================


// Fill Variables
   
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
  
  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
  
  std::vector<double> ThetaMuBin,Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,dansPTBins,dansPZBins,minosRBins;
  
  if (doDIS){
    Enubin = binsDef->GetDISBins("Enu"); 
    Emubin = binsDef->GetDISBins("Emu"); 
    Ehadbin = binsDef->GetDISBins("Ehad");
    Q2bin = binsDef->GetDISBins("Q2");
    Wbin = binsDef->GetDISBins("W");
    xbin = binsDef->GetDISBins("x");
    ybin = binsDef->GetDISBins("y");
    ThetaMuBin = binsDef->GetDISBins("ThetaMu");
  }
  else{
    Enubin = binsDef->GetEnergyBins("Enu"); 
    Emubin = binsDef->GetEnergyBins("Emu"); 
    Ehadbin = binsDef->GetEnergyBins("Ehad");
    Q2bin = binsDef->GetEnergyBins("Q2");
    Wbin = binsDef->GetEnergyBins("W");
    xbin = binsDef->GetEnergyBins("x");
    ybin = binsDef->GetEnergyBins("y");
    ThetaMuBin = binsDef->GetEnergyBins("ThetaMu");
    dansPTBins = binsDef->GetEnergyBins("pTmu");
    dansPZBins = binsDef->GetEnergyBins("pZmu");
    minosRBins = binsDef->GetEnergyBins("minosRBins");
  }
  //Q2bin = binsDef->GetSidebandBins("Q2");
  //Wbin = binsDef->GetSidebandBins("W");

  // 1D Variables
  Var* thetaMu = new Var("GetThetamuDeg", "GetThetamuDeg (Degree)", ThetaMuBin, &CVUniverse::GetThetamuDeg, &CVUniverse::GetThetamuTrueDeg);
  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  Var* ehad = new Var("Ehad", "Ehad (GeV)", Ehadbin, &CVUniverse::GetEhadGeV, &CVUniverse::GetEhadTrueGeV);
  Var* Q2 = new Var("Q2", "Q2 (GeV^2)", Q2bin, &CVUniverse::GetQ2RecoGeV, &CVUniverse::GetQ2TrueGeV);
  Var* W = new Var("W", "W (GeV)", Wbin, &CVUniverse::GetWRecoGeV, &CVUniverse::GetWTrueGeV);
  Var* emu = new Var("Emu", "Emu (GeV)", Emubin, &CVUniverse::GetMuonEGeV, &CVUniverse::GetMuonETrueGeV);
  Var* x = new Var("x", "x", xbin, &CVUniverse::GetxReco, &CVUniverse::GetxTrue);
  Var* y = new Var("y", "y", ybin, &CVUniverse::GetyReco, &CVUniverse::GetyTrue);
  Var* pTmu = new Var("pTmu", "pTmu", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue);
  Var* pZmu = new Var("pZmu", "pZmu", dansPZBins, &CVUniverse::GetMuonPZ, &CVUniverse::GetMuonPZTrue);
//  Var* xVertexTrue = new Var("xVertexTrue", "xVertexTrue",xVertexTrueBin, &CVUniverse::GetxVertexReco, &CVUniverse::GetxVertexTrue); //added June 17/21
  Var* minosR = new Var("minosR", "minosR", minosRBins, &CVUniverse::GetMinosR, &CVUniverse::GetMinosRTrue);

//  variables = {emu, ehad, enu, thetaMu, x, y, Q2,pTmu, pZmu}; //{enu,ehad}; 
  variables = {};
  //xVertexTrue->WriteAllHistogramsToFile(fout, true); //added June 17/21

  // 2D Variables 
  Var2D* W_Q2 = new Var2D(*W, *Q2);
  Var2D* enu_ehad = new Var2D(*enu, *ehad);
  Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
  Var2D* x_y = new Var2D(*x, *y);  // y var
  Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
  Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu); 
  Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu); 
 
//  variables2d = {emu_ehad,enu_ehad, x_y, W_Q2, pZmu_pTmu, pTmu_pZmu};//{enu_ehad, Q2_W};
  variables2d = {pZmu_pTmu}; 
  num_2Dvars = variables2d.size();
  for (auto v : variables2d) v->InitializeAllHistograms(error_bands);
  for (auto v : variables) v->InitializeAllHistograms(error_bands);


  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0; 


  
  CVUniverse *dataverse = new CVUniverse(chain,0);
    
  //=========================================
  // Targets combining Loop
  //=========================================

  ////for(int t = 0; t < targetIDs.size(); t++){

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
        if(isMC){     
          for (auto band : error_bands){
            std::vector<CVUniverse*> error_band_universes = band.second;
            for (auto universe : error_band_universes){
	            // Tell the Event which entry in the TChain it's looking at
	            //=========================================
	            // CUTS in each universe
	            //========================================
	            universe->SetEntry(i);
              reco0++;
              mcreco0++;

              if ((universe->GetInt("has_interaction_vertex") != 1)) continue;

              if(!cutter->PassReco(universe,helicity)) continue;
              reco1++;
              mcreco1++;

              if(!cutter->TrackerOnly(universe)) continue;
              reco2++;
              mcreco2++;

              if(!cutter->InHexagon(universe, 850)) continue;
              reco3++;
              mcreco3++;
         
              if(!cutter->PassThetaCut(universe)) continue;
              reco4++;
              mcreco4++;
     
              reco5++;
              mcreco5++;


              int ring = universe->GetRingBeamCoord();
              for (auto v : variables2d){
                    // filling reco distributions before signal defn cuts
	            v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                    v->ring_m_selected_mc_reco[ring].univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
	            
                     if ((universe->GetInt("mc_incoming") == inputpdg) && (universe->GetInt("mc_current") == 1)){
                           v->signal_purityNum.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                           v->ring_signal_purityNum[ring].univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                           if (universe->GetInt("mc_intType") == 1) v->reco_QE.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                    else if (universe->GetInt("mc_intType") == 8) v->reco_2p2h.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                    else if (universe->GetInt("mc_intType") == 3){
                            if (cutter->PassDISCut(universe)) v->reco_trueDIS.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            else {v->reco_softDIS.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
                    }
                    else if (universe->GetInt("mc_intType") == 2) v->reco_resonant.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                    else {v->reco_other.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}

			   mcsignal++;

                   } // end signal if statement  
                   else{
                            v->bkg_total.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            // breakdown
			    if ((universe->GetInt("mc_incoming") == -14) && (universe->GetInt("mc_current") == 1)) v->ANCC.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            else if ((universe->GetInt("mc_incoming") == -14) && (universe->GetInt("mc_current") != 1)) v->ANNC.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            else if ((universe->GetInt("mc_incoming") == 14) && (universe->GetInt("mc_current") == 1)) v->NCC.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            else if ((universe->GetInt("mc_incoming") == 14) && (universe->GetInt("mc_current") != 1)) {v->NNC.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}
                            else {v->ANNCCNCother.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());}

			    v->ring_bkg_total[ring].univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight());
                            mcbkg++;

                   } //end background else statement

}	
            } // End band's universe loop
          }// End Band loop
        }

        else{

        dataverse->SetEntry(i);
        reco0++;
        data0++;

        if ((dataverse->GetInt("has_interaction_vertex") != 1)) continue;

	 if(!cutter->PassReco(dataverse,helicity)) continue;
	 reco1++;
         data1++;  

         if(!cutter->TrackerOnly(dataverse)) continue;
         reco2++;
         data2++;

         if(!cutter->InHexagon(dataverse, 850)) continue;
         reco3++;
         data3++;

         if(!cutter->PassThetaCut(dataverse)) continue;
         reco4++;
         data4++;        

              reco5++;
              data5++;

              int ring = dataverse->GetRingBeamCoord();
 
 	      for (auto v : variables2d){
	        v->m_selected_data_reco.hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
                v->ring_m_selected_data_reco[ring].hist->Fill(v->GetRecoValueX(*dataverse), v->GetRecoValueY(*dataverse));
	      }
	   
      }

    }//End entries loop
  //}//Target loop closed


  for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
    delete band_universes[i_universe];
  } 
    
  delete dataverse;

  // Printing summary
  std::cout << "**********************************" << std::endl;
  std::cout << "Printing the ";
    isMC? std::cout << "MC ": std::cout << "Data ";
  std::cout << "Summary " << std::endl;
  std::cout << "No cuts = " << reco0 << std::endl;
  std::cout << "Reco = " << reco1 << std::endl;
  std::cout << "Tracker Only = " << reco2 << std::endl;
  std::cout << "In Hexagon = " << reco3 << std::endl;
  std::cout<<"Muon theta < 20 degrees  = " << reco4  << std::endl;
 // std::cout<<"Minus 9999 Cut (if applied)  = " << reco5  << std::endl;
  std::cout << "**********************************" << std::endl;
  
  //return variables;
}
//=============================================================================

