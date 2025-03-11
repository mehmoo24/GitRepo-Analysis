#include "../../NUKECCSRC/ana_common/include/CommonIncludes.h"
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"

#include "../include/Variable_dEdxStudy.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"

#include "../../NUKECCSRC/ana_common/include/NukeCC_Binning.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"

#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"

#include <iostream>
#include <stdlib.h>
#include "/exp/minerva/app/users/mmehmood/MAT_AL9/GitRepo-Analysis/ana/include/systematics/Systematics.h"
#include "../../NUKECCSRC/ana_common/include/NukeCCUtilsNSF.h"
#include "../../NUKECCSRC/ana_common/include/NukeCC_Cuts.h"
#include "PlotUtils/MnvPlotter.h"
#include "TParameter.h"

// ROOT's interpreter, CINT, doesn't understand some legitimate c++ code so we
// shield it.
#ifndef __CINT__
#include "../include/plotting_functions.h"
#endif
#include "PlotUtils/MacroUtil.h"
//using namespace globalV;
using namespace NUKECC_ANA;

typedef VarLoop::Variable Var;
typedef Var2DLoop::Variable2D Var2D;

// global variables
int inputpdg;

void FillVariable( PlotUtils::ChainWrapper* chain, HelicityType::t_HelicityType helicity, NukeCCUtilsNSF *utils , NukeCC_Cuts *cutter ,NukeCC_Binning  *binsDef ,std::vector<Var*>& variables,std::vector<Var2D*>& variables2d,bool isMC,int targetID, int targetZ, const string playlist, bool doDIS){

  std::map<std::string, std::vector<CVUniverse*> > error_bands = GetErrorBands(chain);

  std::vector<double> Enubin, gamma_E, gamma_phi, gamma_dEdx;
  Enubin = binsDef->GetEnergyBins("Enu");
  gamma_E = binsDef->GetEnergyBins("gamma_E");
  gamma_phi = binsDef->GetEnergyBins("gamma_phi");
  gamma_dEdx = binsDef->GetEnergyBins("gamma_dEdx");

  Var* enu = new Var("Enu", "Enu (GeV)", Enubin, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV);
  // be careful! Put in reco func for truth as well below!!!
  Var* gamma1_E = new Var("gamma1_E", "gamma1_E (GeV)", gamma_E, &CVUniverse::GetGamma1E_GeV, &CVUniverse::GetGamma1E_GeV); // be careful! Put in reco func for truth as well!!!
  Var* gamma2_E = new Var("gamma2_E", "gamma2_E (GeV)", gamma_E, &CVUniverse::GetGamma2E_GeV, &CVUniverse::GetGamma2E_GeV); // be careful! Put in reco func for truth as well!!!
  Var* gamma1_phi = new Var("gamma1_phi", "gamma1_phi (rad)", gamma_phi, &CVUniverse::GetGamma1phi_radIthink, &CVUniverse::GetGamma1phi_radIthink);
  Var* gamma2_phi = new Var("gamma2_phi", "gamma2_phi (rad)", gamma_phi, &CVUniverse::GetGamma2phi_radIthink, &CVUniverse::GetGamma2phi_radIthink);
  Var* gamma1_dEdx = new Var("gamma1_dEdx", "gamma1_dEdx (GeV/cm)", gamma_dEdx, &CVUniverse::GetGamma1dEdx_GeVcm, &CVUniverse::GetGamma1dEdx_GeVcm);
  Var* gamma2_dEdx = new Var("gamma2_dEdx", "gamma2_dEdx (GeV/cm)", gamma_dEdx, &CVUniverse::GetGamma2dEdx_GeVcm, &CVUniverse::GetGamma2dEdx_GeVcm);


  variables = {gamma1_E, gamma2_E, gamma1_phi, gamma2_phi, gamma1_dEdx, gamma2_dEdx};
  

  for (auto v : variables) v->InitializeAllHistograms(error_bands);

  CVUniverse *dataverse = new CVUniverse(chain,0);
  std::cout<<"# of entries = "<<chain->GetEntries()<<std::endl;
    for(int i=0; i<chain->GetEntries(); ++i){
      if(i%500000==0) std::cout << (i/1000) << "k " << std::endl;
      if(isMC){
          for (auto band : error_bands){
		  std::vector<CVUniverse*> error_band_universes = band.second;
            for (auto universe : error_band_universes){
		    universe->SetEntry(i);
		   // if ((universe->GetInt("has_interaction_vertex") != 1)) continue;
		    for (auto v : variables){
			    v->m_selected_mc_reco.univHist(universe)->Fill(v->GetRecoValue(*universe), universe->GetWeight());
		    } // end variables for loop
	    }
	   }
	} // end isMC if statement

      else{

        dataverse->SetEntry(i);
	//if ((dataverse->GetInt("has_interaction_vertex") != 1)) continue;
	for (auto v : variables){
                v->m_selected_data_reco.hist->Fill(v->GetRecoValue(*dataverse));
        } // end variables for loop
      }
     } // end entries loop
  

} // end FillVariable

int main(int argc, char *argv[]){
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);

  TString dir(argv[1]);
  int targetID = 99; int targetZ = 99;

  bool doDIS=false;

  const string playlist= argv[2];
  const std::string plist_string(playlist);

  inputpdg = atoi(argv[3]); // 14 for neutrino mode (wanting to look at neutrinos in the beam), -14 for anti neutrino mode

  const string name_of_file = argv[4]; // name of file (ideally include the date)
//  const string fewer_files_opt = argv[5];
  const string fewer_files_opt = "";

  const std::string mc_file_list(Form("../include/playlists/Noes_Tuples/mad_mc_%s%s.txt", plist_string.c_str(), fewer_files_opt.c_str()));
    const std::string data_file_list(Form("../include/playlists/Noes_Tuples/mad_data_%s%s.txt",plist_string.c_str(), fewer_files_opt.c_str()));

  const std::string reco_tree_name("MasterAnaDev");
  const bool wants_truth = false;
  const bool is_grid = false;

  cout << "check 1" << endl;

  PlotUtils::MacroUtil util("MasterAnaDev", mc_file_list, data_file_list, plist_string, wants_truth); //is_grid option removed when did git pull MAT May 9/22

  util.PrintMacroConfiguration("main");



  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::RPAMaterials(false);
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);
  PlotUtils::MinervaUniverse::SetReadoutVolume("Tracker");
  PlotUtils::MinervaUniverse::SetMHRWeightNeutronCVReweight(true);
  PlotUtils::MinervaUniverse::SetMHRWeightElastics(true);

  NukeCCUtilsNSF  *utils   = new NukeCCUtilsNSF(plist_string);
  NukeCC_Cuts     *cutter  = new NukeCC_Cuts();
  NukeCC_Binning  *binsDef = new NukeCC_Binning();

  PlotUtils::ChainWrapper* chainData = util.m_data;
  PlotUtils::ChainWrapper* chainMC = util.m_mc;

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

  TFile fout(dir.Append(histFileName),"RECREATE");

  std::cout << "Processing MC and filling histograms" << std::endl;

  FillVariable(chainMC, helicity, utils, cutter,binsDef,variablesMC,variables2DMC,true,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesMC) {v->m_selected_mc_reco.SyncCVHistos();
                               // v->signal_purityNum.SyncCVHistos();
                               // v->bkg_total.SyncCVHistos();
  }

  std::cout << "Processing Data and filling histograms" << std::endl;
  FillVariable(chainData, helicity, utils, cutter,binsDef,variablesData,variables2DData,false,targetID, targetZ, plist_string,doDIS);
  for (auto v : variablesData) v->m_selected_data_reco.SyncCVHistos();

  // WRITE HISTOGRAMS TO FILE

  // 1D variables
  for (auto v : variablesMC) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  for (auto v : variablesData) {
    v->WriteAllHistogramsToFile(fout, false);
  }

  //Writing POT to the HistFile
  fout.cd();
  auto dataPOTOut = new TParameter<double>("DataPOT", DataPot);
  auto mcPOTOut = new TParameter<double>("MCPOT", MCPot);
  dataPOTOut->Write();
  mcPOTOut->Write();

} // end main function
