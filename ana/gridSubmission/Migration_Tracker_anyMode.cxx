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
//#include "../include/Variable_v1.h"
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
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "/exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/include/systematics/Systematics.h"
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

// global variable
int inputpdg;

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC, int targetID=1, int targetZ=26, const string playlist="minervame1A", bool doDIS=true);
 
//============================================================================================================================
// FillVariable 
//============================================================================================================================
    
void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){
   // std::map< std::string, std::vector<CVUniverse*> > error_bands = utils->GetErrorBands(chain)
   
   std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);
   
   std::map<std::string, std::vector<CVUniverse*> >::iterator itr;
   
   std::map<const std::string, int> error_name;
   for(itr = error_bands.begin(); itr != error_bands.end(); ++itr) error_name.insert(pair<std::string, const int>((itr->second)[0]->ShortName(), (itr->second).size())); 

   std::map<std::string, const int>::iterator itr_m;
   
   std::vector<double> ThetaMuBin, Enubin,Emubin,Ehadbin,xbin,ybin,Q2bin,Wbin,dansPTBins,dansPZBins;

   if (doDIS){
     Enubin  = binsDef->GetDISBins("Enu"); 
     Emubin  = binsDef->GetDISBins("Emu"); 
     Ehadbin = binsDef->GetDISBins("Ehad");
     Q2bin = binsDef->GetDISBins("Q2");
     Wbin = binsDef->GetDISBins("W");
     xbin    = binsDef->GetDISBins("x");
     ybin    = binsDef->GetDISBins("y");
     ThetaMuBin = binsDef->GetDISBins("ThetaMu");
     }
   else{
     Enubin  = binsDef->GetEnergyBins("Enu"); 
     Emubin  = binsDef->GetEnergyBins("Emu"); 
     Ehadbin = binsDef->GetEnergyBins("Ehad");
     Q2bin = binsDef->GetEnergyBins("Q2");
     Wbin = binsDef->GetEnergyBins("W");
     xbin    = binsDef->GetEnergyBins("x");
     ybin    = binsDef->GetEnergyBins("y");
     dansPTBins = binsDef->GetEnergyBins("pTmu");
    dansPZBins = binsDef->GetEnergyBins("pZmu");
   }
  
   //Q2bin = binsDef->GetSidebandBins("Q2");
   //Wbin  = binsDef->GetSidebandBins("W");
   //For 1D varriable

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

   
   //std::vector<Var*> variables = {enu,ehad}; 
   variables = {};//{enu,ehad}; 
   
   //For 2D variable

   Var2D* W_Q2     = new Var2D(*W, *Q2);
   Var2D* enu_ehad = new Var2D(*enu, *ehad);
   Var2D* emu_ehad = new Var2D(*emu, *ehad);  // y var
   Var2D* x_Q2 = new Var2D(*x, *Q2);  // y var
   Var2D* x_y = new Var2D(*x, *y);  // y var
   Var2D* pZmu_pTmu = new Var2D(*pZmu, *pTmu);
  Var2D* pTmu_pZmu = new Var2D(*pTmu, *pZmu);

   variables2d = {pZmu_pTmu};//{enu_ehad, Q2_W};
   //variables2d = {emu_ehad};//{enu_ehad, Q2_W};
   
   for (auto v : variables2d) v->InitializeAllHistograms(error_bands);

   for (auto v : variables) v->InitializeAllHistograms(error_bands);

   //Migration starts here!
   for (auto v : variables2d) {
      v->SetupResponse(error_name);
      v->SetupResponse_ring0(error_name);
      v->SetupResponse_ring1(error_name);
      v->SetupResponse_ring2(error_name);
   }
   for (auto v : variables) v->SetupResponse1D(error_name);
    
  int reco0=0;
  int reco1=0;
  int reco2=0;
  int reco3=0;
  int reco4=0; 
  int reco5=0; 
  int reco6=0;  
  int reco7 = 0;
  int reco8=0;
   
   CVUniverse *dataverse = new CVUniverse(chain,0);

   
    
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
	       int unv_count = 0;
               std::vector<CVUniverse*> error_band_universes = band.second;
	       for (auto universe : error_band_universes){
		 // Tell the Event which entry in the TChain it's looking at
		 //=========================================
		 // CUTS in each universe
		 //========================================
	 
		 universe->SetEntry(i);
//     reco0++;
     reco0=reco0+universe->GetWeight();

              if ((universe->GetInt("has_interaction_vertex") != 1)) continue;

              if(!cutter->PassReco(universe,helicity)) continue;
              reco1=reco1+universe->GetWeight();

              if(!cutter->TrackerOnly(universe)) continue;
              reco2=reco2+universe->GetWeight();

              if(!cutter->InHexagon(universe, 850)) continue;
              reco3=reco3+universe->GetWeight();

              if(!cutter->PassThetaCut(universe)) continue;
              reco4=reco4+universe->GetWeight();
     // reco cuts
//     if(!cutter->PassReco(universe,helicity)) continue;
//     reco1++;

//	   if(!cutter->IsInMaterial(universe,targetID,targetZ, false)) continue;
//     reco2++;

//	   if(targetID<10 && universe->GetInt("MasterAnaDev_targetID") != targetID) continue;
//	   reco3++;

//     if( universe->GetVecElem("ANN_plane_probs",0) < MIN_PROB_PLANE_CUT ) continue;
//     reco4++;

//     if (!cutter->PassTruth(universe, helicity)) continue; // True fiducial, true CC, true antinu
//     reco5++;

//	   if(!cutter->IsInTrueMaterial(universe,targetID, targetZ,false)) continue; // true target + material
//     reco6++;

           if ((universe->GetInt("mc_incoming") != inputpdg) || (universe->GetInt("mc_current") != 1)) continue;  
           int ring = universe->GetRingBeamCoord();
		 for (auto v : variables2d){
      
  //    v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe), universe->GetWeight()); 
      v->m_selected_Migration.univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight()); 
      v->ring_m_selected_Migration[ring].univHist(universe)->Fill(v->GetTrueValueX(*universe), v->GetTrueValueY(*universe), universe->GetWeight());  
      
      //Migration stuff
      v->FillResponse(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count); 
      if(ring==0) v->FillResponse_ring0(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count);
      else if(ring==1) v->FillResponse_ring1(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count);
      else if(ring==2) v->FillResponse_ring2(v->GetRecoValueX(*universe), v->GetRecoValueY(*universe),v->GetTrueValueX(*universe), v->GetTrueValueY(*universe),universe->ShortName(),universe->GetWeight(),unv_count);

 		 }
 //     std::cout << universe->ShortName() << " " << unv_count << std::endl;

	   for (auto v : variables){

	     //v->m_selected_mc_reco.univHist(universe)->Fill(v->GetTrueValue(*universe, 0), universe->GetWeight());
             v->m_selected_Migration.univHist(universe)->Fill(v->GetRecoValue(*universe), v->GetTrueValue(*universe), universe->GetWeight());
       //1D response
       v->FillResponse1D(v->GetRecoValue(*universe),v->GetTrueValue(*universe),universe->ShortName(),universe->GetWeight(),unv_count); 
	   	}
		 unv_count++;
    } // End band's universe loop
     }
  }//End entries loop
 
 for (auto v : variables2d) {
    v->getResponseObjects(error_bands);
    v->getResponseObjects_ring0(error_bands);
    v->getResponseObjects_ring1(error_bands);
    v->getResponseObjects_ring2(error_bands);
 }
 for (auto v : variables) v->getResponseObjects1D(error_bands);
 
 for(auto band : error_bands){
    std::vector<CVUniverse*> band_universes = band.second;
    for(unsigned int i_universe = 0; i_universe < band_universes.size(); ++i_universe)
      delete band_universes[i_universe];
 } 
    
 delete dataverse;

   std::cout<<"**********************************"<<std::endl;
      std::cout<<" Summary "<<std::endl;
   std::cout<<"**********************************"<<std::endl;
     std::cout<<"Migration Matrix "<<std::endl;
//    std::cout << "No cuts = " << reco0 << std::endl;
//    std::cout << "Reco Cut = " << reco1 << std::endl;
//    std::cout << "Material Cut = " << reco2 << std::endl;
//    std::cout << "TargetID Cuts = " << reco3 << std::endl;
//    std::cout << "Plane prob. cut = " << reco4 << std::endl;
//    std::cout << "Truth cut (fiducial, CC, antinu) = " << reco5 << std::endl;
//    std::cout << "True Material cut  = "<< reco6 << std::endl;
//    std::cout << "Muon kinematics cut  = " << reco7 << std::endl;
//    std::cout << "True muon kinematics cut  = " << reco8 << std::endl;
}

//============================================================================================================================
// Main
//============================================================================================================================

int main(int argc, char *argv[]){
//	 ROOT::Cintex::Cintex::Enable();
	 TH1::AddDirectory(false);

	 if(argc==1){
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     std::cout<<"MACROS HELP:\n\n"<<
	       "\t-./runEventLoop Path_to_Output_file\n\n"<<
	       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<< std::endl;
//	       "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
//	       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n" << std::endl;
	     std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
	     return 0;
	 }

	 TString dir(argv[1]);
//	 int targetID = atoi(argv[2]);
//	 int targetZ  = atoi(argv[3]);

         int targetID = 99;
         int targetZ = 99;

         const string playlist= argv[2];

    const std::string plist_string(playlist);

    inputpdg = atoi(argv[3]);  // 14 for neutrino mode (wanting to look at neutrinos in the beam), -14 for anti neutrino mode

    const string name_of_file = argv[4]; 
//    const string fewer_files_opt = argv[5]; //1run or ""
    const string fewer_files_opt = "";

    const std::string mc_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_mc_%s%s.txt", plist_string.c_str(), fewer_files_opt.c_str()));
    const std::string data_file_list(Form("../include/playlists/MAD_tuples_p4_all/mad_data_%s%s.txt",plist_string.c_str(), fewer_files_opt.c_str()));


		bool doDIS = false; 


    const std::string reco_tree_name("MasterAnaDev");
//    const std::string plist_string("minervame"+option);
 
//    const std::string plist_string("minervame6A");
    const bool wants_truth = false;
    const bool is_grid = false;

	 PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list, plist_string, wants_truth);

	 util.PrintMacroConfiguration("main");

	 //=========================================
	 // Systematics
	 //=========================================
	 //std::map<std::string, std::vector<CVUniverse*> > error_bands =
	 //  GetErrorBands(util.m_mc);

	 PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
         PlotUtils::MinervaUniverse::RPAMaterials(false);
	 PlotUtils::MinervaUniverse::SetNuEConstraint(true);
//	 PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
	 PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
//         PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
         PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
         PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
         // MnvHadronReweighter (GEANT Hadron sytematics)
         PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
         //Neutron CV reweight is on by default (recommended you keep this on)
         PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
         //Elastics are on by default (recommended you keep this on)  
         PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);
  
	    
	 NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
	 NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
	 NukeCC_Binning  *binsDef = new NukeCC_Binning();
	  
	 PlotUtils::ChainWrapper* chainData = util.m_data;
	 PlotUtils::ChainWrapper* chainMC = util.m_mc;
//	 HelicityType::t_HelicityType helicity = utils->GetHelicityFromPlaylist(plist_string);

         HelicityType::t_HelicityType helicity;
         if(inputpdg == 14) helicity = NUKECC_ANA::HelicityType::kNeutrino;
         else helicity = NUKECC_ANA::HelicityType::kAntiNeutrino;

	 double DataPot=  util.m_data_pot; 
	 double MCPot=  util.m_mc_pot; 
	 //double total_pot_data,total_pot_mc;
	 //utils->getPOT(total_pot_data,total_pot_mc);  
	 double  MCscale=DataPot/MCPot;
	 //double  MCscale=1.0;
	 
	 std::cout<<" MCScale= "<<MCscale<<std::endl; 
	 std::vector<Var*> variablesMC,variablesData; 
	 std::vector<Var2D*> variables2DMC,variables2DData; 


         TString histFileName;
        if(RunCodeWithSystematics){
           histFileName = utils->GetHistFileName(Form("Migration_%s_%s_%ssys",name_of_file.c_str(), playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );
         }
         else{
           histFileName = utils->GetHistFileName(Form("Migration_%s_%s_%snosys",name_of_file.c_str(), playlist.c_str(),fewer_files_opt.c_str()), FileType::kAny, targetID, targetZ, helicity );
         }
	   
	 TFile fout(dir.Append(histFileName),"RECREATE");	
	   
	 // For 1D variables 
	 FillVariable(chainMC, helicity, utils, cutter, binsDef, variablesMC, variables2DMC, true, targetID, targetZ, plist_string, doDIS);
	     
         	 
	 for (auto v : variablesMC){ v-> mresp1D.SyncCVHistos();
                                     v-> m_selected_Migration.SyncCVHistos();
                                   }
	 for (auto v : variables2DMC) {
           v-> mresp.SyncCVHistos();
          // v-> m_selected_mc_reco.SyncCVHistos();
           v-> m_selected_Migration.SyncCVHistos();
           for(int ring=0;ring<3;ring++){
              v-> ring_m_selected_Migration[ring].SyncCVHistos();
              v-> ring_mresp[ring].SyncCVHistos();
           }
         }
	 
	 
	 for (auto v : variablesMC) {
	   v->WriteAllHistogramsToFileMig(fout, true);
	 }


	 //For 2D variable

	 for (auto v : variables2DMC) {
	   v->WriteAllHistogramsToFileMig(fout,true);
	 }

  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write(); 
	 

}//End Main

