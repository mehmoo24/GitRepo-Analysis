#ifndef VARIABLE_H
#define VARIABLE_H

#include <iterator>
#include "../../NUKECCSRC/ana_common/include/CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"
#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#include "PlotUtils/Variable2DBase.h"
#include "MinervaUnfold/MnvResponse.h"
#include <bitset>
#include "PlotUtils/AnaBinning.h"

#endif  // __CINT__

namespace VarLoop {

class Variable : public PlotUtils::VariableBase<NUKECC_ANA::CVUniverse> {
 private:
  typedef PlotUtils::HistWrapper<NUKECC_ANA::CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable(ARGS... args) : PlotUtils::VariableBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  // selected mc reco histwrapper
  HW m_selected_mc_reco,  m_selected_data_reco, m_selected_truth_reco, m_selected_data_reco_sb, m_selected_mc_true,reco_QE,
       reco_2p2h,
       reco_trueDIS,
       reco_softDIS,
       reco_other,
       reco_resonant,
       signal_purityNum,
       bkg_total,
       ANCC, ANNC, NCC, NNC, ANNCCNCother,
       m_ringone_mc_true, m_ringtwo_mc_true, m_ringthree_mc_true,
       m_ringfour_mc_true, m_ringfive_mc_true, m_ringsix_mc_true,
       m_petal0_mc_true, m_petal1_mc_true, m_petal2_mc_true, m_petal3_mc_true,
       m_petal4_mc_true, m_petal5_mc_true, m_petal6_mc_true, m_petal7_mc_true,
       m_petal8_mc_true, m_petal9_mc_true, m_petal10_mc_true, m_petal11_mc_true;

  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  HW2D mresp1D;
  typedef PlotUtils::MnvH2D MH2D;
  HW2D m_selected_Migration;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response1D;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2; 

  MnvH2D *migrationH2D = NULL;
  MnvH1D *h_reco1D = NULL;
  MnvH1D *h_truth1D = NULL;
 // HISTFOLIO
  // selected mc reco - signal background histfolio
  PlotUtils::HistFolio<MH1D> m_selected_mc_sb;
  // PlotUtils::MH1D* m_selected_data_sb;
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    const bool clear_bands = true;  // we want empty histograms

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH1D* dummy_selected_mc_reco = new MH1D(Form("h_mc_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_reco = HW(dummy_selected_mc_reco, univs, clear_bands);

    MH1D* dummy_reco_QE =
        new MH1D(Form("reco_QE_%s", name), name, GetNBins(), bins.data());
    reco_QE = HW(dummy_reco_QE, univs, clear_bands);
    MH1D* dummy_reco_2p2h =
        new MH1D(Form("reco_2p2h_%s", name), name, GetNBins(), bins.data());
    reco_2p2h = HW(dummy_reco_2p2h, univs, clear_bands);
    MH1D* dummy_reco_trueDIS =
        new MH1D(Form("reco_trueDIS_%s", name), name, GetNBins(), bins.data());
    reco_trueDIS = HW(dummy_reco_trueDIS, univs, clear_bands);
    MH1D* dummy_reco_softDIS =
        new MH1D(Form("reco_softDIS_%s", name), name, GetNBins(), bins.data());
    reco_softDIS = HW(dummy_reco_softDIS, univs, clear_bands);
    MH1D* dummy_reco_other =
        new MH1D(Form("reco_other_%s", name), name, GetNBins(), bins.data());
    reco_other = HW(dummy_reco_other, univs, clear_bands);
    MH1D* dummy_reco_resonant =
        new MH1D(Form("reco_resonant_%s", name), name, GetNBins(), bins.data());
    reco_resonant = HW(dummy_reco_resonant, univs, clear_bands);
    MH1D* dummy_signal_purityNum =
        new MH1D(Form("signal_purityNum_%s", name), name, GetNBins(), bins.data());
    signal_purityNum = HW(dummy_signal_purityNum, univs, clear_bands);
    MH1D* dummy_bkg_total =
        new MH1D(Form("bkg_total_%s", name), name, GetNBins(), bins.data());
    bkg_total = HW(dummy_bkg_total, univs, clear_bands);  
   
    MH1D* dummy_ANCC =
        new MH1D(Form("AntiNeutrino_CC_%s", name), name, GetNBins(), bins.data());
    ANCC = HW(dummy_ANCC, univs, clear_bands);
    MH1D* dummy_ANNC =
        new MH1D(Form("AntiNeutrino_NC_%s", name), name, GetNBins(), bins.data());
    ANNC = HW(dummy_ANNC, univs, clear_bands);
    MH1D* dummy_NCC =
        new MH1D(Form("Neutrino_CC_%s", name), name, GetNBins(), bins.data());
    NCC = HW(dummy_NCC, univs, clear_bands);
    MH1D* dummy_NNC =
        new MH1D(Form("Neutrino_NC_%s", name), name, GetNBins(), bins.data());
    NNC = HW(dummy_NNC, univs, clear_bands);
    MH1D* dummy_ANNCCNCother =
        new MH1D(Form("ANNCCNCother_%s", name), name, GetNBins(), bins.data());
    ANNCCNCother = HW(dummy_ANNCCNCother, univs, clear_bands);



    //for 1D analysis
    MH1D* dummy_selected_mc_true = new MH1D(Form("h_mc_true_%s", name), name,GetNBins(), bins.data());
    m_selected_mc_true = HW(dummy_selected_mc_true, univs, clear_bands);
    
    //////////For Efficiency 
    MH1D* dummy_selected_truth_reco = new MH1D(Form("h_truth_%s", name), name, GetNBins(), bins.data());
    m_selected_truth_reco = HW(dummy_selected_truth_reco, univs, clear_bands);

    MH1D* dummy_ringone = new MH1D(Form("h_truth_ringone_%s", name), name, GetNBins(), bins.data());
    m_ringone_mc_true = HW(dummy_ringone, univs, clear_bands);
    MH1D* dummy_ringtwo = new MH1D(Form("h_truth_ringtwo_%s", name), name, GetNBins(), bins.data());
    m_ringtwo_mc_true = HW(dummy_ringtwo, univs, clear_bands);
    MH1D* dummy_ringthree = new MH1D(Form("h_truth_ringthree_%s", name), name, GetNBins(), bins.data());
    m_ringthree_mc_true = HW(dummy_ringthree, univs, clear_bands);
    MH1D* dummy_ringfour = new MH1D(Form("h_truth_ringfour_%s", name), name, GetNBins(), bins.data());
    m_ringfour_mc_true = HW(dummy_ringfour, univs, clear_bands);
    MH1D* dummy_ringfive = new MH1D(Form("h_truth_ringfive_%s", name), name, GetNBins(), bins.data());
    m_ringfive_mc_true = HW(dummy_ringfive, univs, clear_bands);
    MH1D* dummy_ringsix = new MH1D(Form("h_truth_ringsix_%s", name), name, GetNBins(), bins.data());
    m_ringsix_mc_true = HW(dummy_ringsix, univs, clear_bands);


    MH1D* dummy_petal0 = new MH1D(Form("h_truth_petal0_%s", name), name, GetNBins(), bins.data());
    m_petal0_mc_true = HW(dummy_petal0, univs, clear_bands);
    MH1D* dummy_petal1 = new MH1D(Form("h_truth_petal1_%s", name), name, GetNBins(), bins.data());
    m_petal1_mc_true = HW(dummy_petal1, univs, clear_bands);
    MH1D* dummy_petal2 = new MH1D(Form("h_truth_petal2_%s", name), name, GetNBins(), bins.data());
    m_petal2_mc_true = HW(dummy_petal2, univs, clear_bands);
    MH1D* dummy_petal3 = new MH1D(Form("h_truth_petal3_%s", name), name, GetNBins(), bins.data());
    m_petal3_mc_true = HW(dummy_petal3, univs, clear_bands);
    MH1D* dummy_petal4 = new MH1D(Form("h_truth_petal4_%s", name), name, GetNBins(), bins.data());
    m_petal4_mc_true = HW(dummy_petal4, univs, clear_bands);
    MH1D* dummy_petal5 = new MH1D(Form("h_truth_petal5_%s", name), name, GetNBins(), bins.data());
    m_petal5_mc_true = HW(dummy_petal5, univs, clear_bands);
    MH1D* dummy_petal6 = new MH1D(Form("h_truth_petal6_%s", name), name, GetNBins(), bins.data());
    m_petal6_mc_true = HW(dummy_petal6, univs, clear_bands);
    MH1D* dummy_petal7 = new MH1D(Form("h_truth_petal7_%s", name), name, GetNBins(), bins.data());
    m_petal7_mc_true = HW(dummy_petal7, univs, clear_bands);
    MH1D* dummy_petal8 = new MH1D(Form("h_truth_petal8_%s", name), name, GetNBins(), bins.data());
    m_petal8_mc_true = HW(dummy_petal8, univs, clear_bands);
    MH1D* dummy_petal9 = new MH1D(Form("h_truth_petal9_%s", name), name, GetNBins(), bins.data());
    m_petal9_mc_true = HW(dummy_petal9, univs, clear_bands);
    MH1D* dummy_petal10 = new MH1D(Form("h_truth_petal10_%s", name), name, GetNBins(), bins.data());
    m_petal10_mc_true = HW(dummy_petal10, univs, clear_bands);
    MH1D* dummy_petal11 = new MH1D(Form("h_truth_petal11_%s", name), name, GetNBins(), bins.data());
    m_petal11_mc_true = HW(dummy_petal11, univs, clear_bands);



    MH2D* dummy_selected_Migration = new MH2D(Form("selected_Migration_%s", name), name, GetNBins(),GetBinVec().data(), GetNBins(), GetBinVec().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);


    //For Data
    MH1D* dummy_selected_data_reco = new MH1D(Form("h_data_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco = HW(dummy_selected_data_reco, univs, clear_bands);


    //For Data in sidebands
    MH1D* selected_data_reco_sb = new MH1D(Form("h_data_sb_%s", name), name, GetNBins(), bins.data());
    m_selected_data_reco_sb = HW(selected_data_reco_sb, univs, clear_bands);
  
    //HISTFOLIO
    //selected mc reco - signal background histfolio
    m_selected_mc_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(Form("selected_mc_sb_%s", name), name, GetNBins(), bins.data());
    
 //PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
 //   "dummy", "dummy", plotting::nbins, plotting::xmin, plotting::xmax);
 //  m_selected_data_sb = PlotUtils::HistFolio<PlotUtils::MnvH1D>(
 //    Form("selected_data_sb_%s", name), name, GetNBins(), bins.data());
   
    m_selected_mc_sb.AddComponentHist("DIS");
    m_selected_mc_sb.AddComponentHist("MC");
   // m_selected_data_sb.AddComponentHist("Data");
delete dummy_selected_mc_reco;
delete dummy_selected_truth_reco;
delete dummy_selected_mc_true;
delete dummy_selected_data_reco;
  }


//=======================================================================================
  // When plotting, don't make Cuts on this Variable's axes
  //=======================================================================================
  //TODO: Move into PlotUtils?
  #define cut_map_t std::bitset<64>
  std::bitset<64> ignoreTheseCuts;
  //Populate ignoreTheseCuts.
  template <class CUTS_T> //Kill two bird templates with one stone template
  cut_map_t MatchCutsToVars(const CUTS_T& cuts)
  {
     const char* name = GetName().c_str();
    //Stupid compiler!  Let me use generic lambdas!
    using cut_t = decltype(*cuts.begin());
    ignoreTheseCuts.reset(); //Make sure I start off ignoring NO Cuts.
                             //Only cuts with 1 will be ignored.
    const auto found = std::find_if(cuts.begin(), cuts.end(),
                                    [this](const cut_t cut)
                                    {
                                      return cut->getName() == this->GetName();
                                    });
  
    if(found != cuts.end()) //If we found a match
    {
      ignoreTheseCuts.flip(std::distance(cuts.begin(), found));
    }
    else std::cerr << "Failed to find ANY Cut that matches " << GetName() << ".  You should double-check your Cut names.\n";
    return ignoreTheseCuts; //Just here in case it helps with debugging
  }
  //Use ignoreTheseCuts
  inline cut_map_t IgnoreMyVars(const cut_map_t allCuts) const
  {
    return allCuts | ignoreTheseCuts;
  }




void SetupResponse1D(std::map<const std::string, int> systematics){
//void SetupResponse(T univs){
	   std::vector<double> bins = GetBinVec();
	   const char* name = GetName().c_str();

	   axis_binning bin_x;
	   bin_x.uniform=false;
	
     	   vector<double> vx;
     	   
	   for(int i=0; i<=GetNBins(); i++){vx.push_back(GetBinVec().data()[i]);}
     	   
	   bin_x.bin_edges = vx;
	   bin_x.nbins	   = GetNBins();
	   bin_x.min 	   = GetBinVec().data()[0];
	   bin_x.max       = GetBinVec().data()[GetNBins()];
           
           cout<<"bins min = "<<bin_x.min<<"bins max = "<<bin_x.max<<"Number of bins = "<<bin_x.nbins<<endl;
            
	   //Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response1d_%s", name), name, bin_x, bin_x, systematics))); 
	   Response1D.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("response1d_%s", name), name, bin_x, bin_x, systematics))); 
}


//===================================================================================
void FillResponse1D(double x_reco, double x_true, const std::string name, double w, int unv){
	
	//std::cout<<name<<std::endl;
 	for(mnv_itr = Response1D.begin(); mnv_itr != Response1D.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,x_true,name,unv, w);
	}		
}

//=====================================

template <typename T>
void getResponseObjects1D(T univs)
{
//  bool status = false;
  for(mnv_itr2 = Response1D.begin(); mnv_itr2 != Response1D.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH2D, h_reco1D, h_truth1D );
        }

//  const bool clear_bands = true; //aug28 
  const bool clear_bands = false;
  migrationH2D->PopVertErrorBand("cv");
  mresp1D = HW2D(migrationH2D, univs, clear_bands);
}
	
  //=======================================================================================
  // WRITE ALL HISTOGRAMS FOR EVENTLOOP
  //=======================================================================================
 void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
               reco_QE.hist->Write();
                 reco_2p2h.hist->Write();
                 reco_trueDIS.hist->Write();
                 reco_softDIS.hist->Write();
                 reco_other.hist->Write();
                 reco_resonant.hist->Write();
                 signal_purityNum.hist->Write();
                 bkg_total.hist->Write();
                 ANCC.hist->Write(); ANNC.hist->Write(); NCC.hist->Write(); NNC.hist->Write(); ANNCCNCother.hist->Write();
}
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
  }
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileEff(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { m_selected_mc_reco.hist->Write();
    }
    else { m_selected_truth_reco.hist->Write();
           /*m_ringone_mc_true.hist->Write();
           m_ringtwo_mc_true.hist->Write();
           m_ringthree_mc_true.hist->Write();
           m_ringfour_mc_true.hist->Write(); 
           m_ringfive_mc_true.hist->Write();
           m_ringsix_mc_true.hist->Write();
           */
          /*
           m_petal0_mc_true.hist->Write();
           m_petal1_mc_true.hist->Write();
           m_petal2_mc_true.hist->Write();
           m_petal3_mc_true.hist->Write();
           m_petal4_mc_true.hist->Write();
           m_petal5_mc_true.hist->Write();
           m_petal6_mc_true.hist->Write();
           m_petal7_mc_true.hist->Write();
           m_petal8_mc_true.hist->Write();
           m_petal9_mc_true.hist->Write();
           m_petal10_mc_true.hist->Write();
           m_petal11_mc_true.hist->Write();
          */
    }
}
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
 void WriteAllHistogramsToFileMig(TFile& f, bool isMC) const {
    f.cd();

    // selected mc reco
    if(isMC) { //m_selected_mc_reco.hist->Write();
               //m_selected_mc_true.hist->Write();
               m_selected_Migration.hist->Write();
               mresp1D.hist->Write();
}
    else {m_selected_data_reco.hist->Write();
}
    // selected mc  histfolio fir Hist Stacking
   //if(isMC) m_selected_mc_sb.WriteToFile(f);
   //else m_selected_data_reco_sb.hist->Write();
  }
};



}  // namespace VarLoop



namespace Var2DLoop {
class Variable2D : public PlotUtils::Variable2DBase<NUKECC_ANA::CVUniverse> {
 private:
  //=======================================================================================
  // TYPEDEFS CONVENIENCE
  //=======================================================================================
  typedef PlotUtils::Hist2DWrapper<NUKECC_ANA::CVUniverse> HW2D;
  //typedef Histogram<NUKECC_ANA::CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable2D(ARGS... args) : Variable2DBase<NUKECC_ANA::CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_selected_mc_reco,
       reco_QE,
       reco_2p2h,
       reco_trueDIS,
       reco_softDIS,
       reco_other,
       reco_resonant,
       true_QE,
       true_2p2h,
       true_trueDIS,
       true_softDIS,
       true_other,
       true_resonant,
       signal_purityNum,
       bkg_total,
       bkg_wrongSign,
       bkg_NC,
       bkg_QE,
       bkg_2p2h,
       bkg_trueDIS,
       bkg_softDIS,
       bkg_other,
       bkg_resonant,  
       ANCC, ANNC, NCC, NNC, ANNCCNCother, 
       h_background_subtracted_mc,
       h_background_mc_scale,
       h_background_subtracted_data,
       m_selected_data_reco,m_selected_truth_reco,m_selected_Migration, mresp;
  //MinervaUnfold::MnvResponse* Response;
  std::map<std::string,MinervaUnfold::MnvResponse*> Response;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr;
  std::map<std::string,MinervaUnfold::MnvResponse*>::iterator mnv_itr2;

  MnvH2D *migrationH = NULL;
  MnvH2D *h_reco = NULL;
  MnvH2D *h_truth = NULL;
  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("h_mc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);
//start of new histwrappers added in
    MH2D* dummy_reco_QE =
        new MH2D(Form("reco_QE_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_QE = HW2D(dummy_reco_QE, univs, clear_bands);
    MH2D* dummy_true_QE =
        new MH2D(Form("true_QE_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_QE = HW2D(dummy_true_QE, univs, clear_bands);

    MH2D* dummy_reco_2p2h =
        new MH2D(Form("reco_2p2h_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_2p2h = HW2D(dummy_reco_2p2h, univs, clear_bands);
    MH2D* dummy_true_2p2h =
        new MH2D(Form("true_2p2h_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_2p2h = HW2D(dummy_true_2p2h, univs, clear_bands);

    MH2D* dummy_reco_trueDIS =
        new MH2D(Form("reco_trueDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_trueDIS = HW2D(dummy_reco_trueDIS, univs, clear_bands);
    MH2D* dummy_true_trueDIS =
        new MH2D(Form("true_trueDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_trueDIS = HW2D(dummy_true_trueDIS, univs, clear_bands);

    MH2D* dummy_reco_softDIS =
        new MH2D(Form("reco_softDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_softDIS = HW2D(dummy_reco_softDIS, univs, clear_bands);
    MH2D* dummy_true_softDIS =
        new MH2D(Form("true_softDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_softDIS = HW2D(dummy_true_softDIS, univs, clear_bands);

    MH2D* dummy_reco_other =
        new MH2D(Form("reco_other_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_other = HW2D(dummy_reco_other, univs, clear_bands);
    MH2D* dummy_true_other =
        new MH2D(Form("true_other_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_other = HW2D(dummy_true_other, univs, clear_bands);

    MH2D* dummy_reco_resonant =
        new MH2D(Form("reco_resonant_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    reco_resonant = HW2D(dummy_reco_resonant, univs, clear_bands);
    MH2D* dummy_true_resonant =
        new MH2D(Form("true_resonant_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    true_resonant = HW2D(dummy_true_resonant, univs, clear_bands);


    MH2D* dummy_signal_purityNum =
        new MH2D(Form("signal_purityNum_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    signal_purityNum = HW2D(dummy_signal_purityNum, univs, clear_bands);
    MH2D* dummy_bkg_total =
        new MH2D(Form("bkg_total_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_total = HW2D(dummy_bkg_total, univs, clear_bands);
    MH2D* dummy_bkg_wrongSign =
        new MH2D(Form("bkg_wrongSign_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_wrongSign = HW2D(dummy_bkg_wrongSign, univs, clear_bands);
    MH2D* dummy_bkg_NC =
        new MH2D(Form("bkg_NC_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_NC = HW2D(dummy_bkg_NC, univs, clear_bands);
    MH2D* dummy_bkg_QE =
        new MH2D(Form("bkg_QE_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_QE = HW2D(dummy_bkg_QE, univs, clear_bands);
    MH2D* dummy_bkg_2p2h =
        new MH2D(Form("bkg_2p2h_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_2p2h = HW2D(dummy_bkg_2p2h, univs, clear_bands);
    MH2D* dummy_bkg_trueDIS =
        new MH2D(Form("bkg_trueDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_trueDIS = HW2D(dummy_bkg_trueDIS, univs, clear_bands);
    MH2D* dummy_bkg_softDIS =
        new MH2D(Form("bkg_softDIS_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_softDIS = HW2D(dummy_bkg_softDIS, univs, clear_bands);
    MH2D* dummy_bkg_other =
        new MH2D(Form("bkg_other_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_other = HW2D(dummy_bkg_other, univs, clear_bands);
    MH2D* dummy_bkg_resonant =
        new MH2D(Form("bkg_resonant_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    bkg_resonant = HW2D(dummy_bkg_resonant, univs, clear_bands);

    MH2D* dummy_ANCC =
        new MH2D(Form("AntiNeutrino_CC_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    ANCC = HW2D(dummy_ANCC, univs, clear_bands);
    MH2D* dummy_ANNC =
        new MH2D(Form("AntiNeutrino_NC_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    ANNC = HW2D(dummy_ANNC, univs, clear_bands);
    MH2D* dummy_NCC =
        new MH2D(Form("Neutrino_CC_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    NCC = HW2D(dummy_NCC, univs, clear_bands);
    MH2D* dummy_NNC =
        new MH2D(Form("Neutrino_NC_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    NNC = HW2D(dummy_NNC, univs, clear_bands);
    MH2D* dummy_ANNCCNCother =
        new MH2D(Form("ANNCCNCother_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    ANNCCNCother = HW2D(dummy_ANNCCNCother, univs, clear_bands);

    MH2D* dummy_h_background_subtracted_mc =
        new MH2D(Form("h_background_subtracted_mc_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    h_background_subtracted_mc = HW2D(dummy_h_background_subtracted_mc, univs, clear_bands);
    MH2D* dummy_h_background_mc_scale =
        new MH2D(Form("h_background_mc_scale_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    h_background_mc_scale = HW2D(dummy_h_background_mc_scale, univs, clear_bands);
    MH2D* dummy_h_background_subtracted_data =
        new MH2D(Form("h_background_subtracted_data_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    h_background_subtracted_data = HW2D(dummy_h_background_subtracted_data, univs, clear_bands);

//end of new histwrappers added in     
    MH2D* dummy_selected_truth_reco =
        new MH2D(Form("h_truth_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);
    
   // MH2D* dummy_selected_truth_reco =
   //     new MH2D(Form("selected_truth2d_reco_%s", name), name, setBinLogX(), GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   //              GetBinVecX().data(), setBinLogY(), GetBinVecY().data());
   // m_selected_truth_reco = HW2D(dummy_selected_truth_reco, univs, clear_bands);

    MH2D* dummy_selected_data_reco =
        new MH2D(Form("h_data_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_data_reco = HW2D(dummy_selected_data_reco, univs, clear_bands);
    

    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/
    

    /////For 2D analysis
    MH2D* dummy_selected_Migration =
        new MH2D(Form("Migration_%s", name), name, GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);///For 2D Migration
   

    ///For 1D analysis
    /*MH2D* dummy_selected_Migration =
        new MH2D(Form("selected_Migration_%s", name), name, GetNBinsX(),GetBinVecX().data(), GetNBinsX(), GetBinVecX().data());
    m_selected_Migration = HW2D(dummy_selected_Migration, univs, clear_bands);*/

    //delete Response; 
    delete dummy_selected_mc_reco;
    delete dummy_selected_truth_reco;
    delete dummy_selected_Migration;
    delete dummy_selected_data_reco;
  }
//=====



void SetupResponse(std::map<const std::string,  int> systematics){
//void SetupResponse(T univs){

	   const char* name = GetName().c_str();
	   axis_binning bin_x, bin_y;
	   bin_x.uniform=false;
	
     	   vector<double> vx;
	   
	   for(int i=0; i<=GetNBinsX(); i++){vx.push_back(GetBinVecX().data()[i]);}
     	   
	   vector<double> vy;
	   for(int j=0; j<=GetNBinsY(); j++){vy.push_back(GetBinVecY().data()[j]);}
	   bin_x.bin_edges = vx;
	   bin_x.nbins	    = GetNBinsX();
	   bin_x.min 	    = GetBinVecX().data()[0];
	   bin_x.max       = GetBinVecX().data()[GetNBinsX()];  
	   bin_y.bin_edges = vy;
	   bin_y.nbins	    = GetNBinsY();
	   bin_y.min 	    = GetBinVecY().data()[0];
	   bin_y.max       = GetBinVecY().data()[GetNBinsY()]; 

	   //Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("selected_mc_response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics))); 
	   Response.insert(pair<const std::string, MinervaUnfold::MnvResponse*>(name, new MinervaUnfold::MnvResponse(Form("response2d_%s", name), name, bin_x, bin_y, bin_x, bin_y, systematics))); 
}

//===================================================================================
//
//===================================================================================
void FillResponse(double x_reco, double y_reco, double x_true, double y_true, const std::string name, double w, int unv){
 	for(mnv_itr = Response.begin(); mnv_itr != Response.end(); ++mnv_itr){
		(mnv_itr->second)->Fill(x_reco,y_reco,x_true,y_true,name,unv, w);
	}		
	
	
}
//===================================================================================
//
//===================================================================================
template <typename T>
void getResponseObjects(T univs)
{
 // bool status = false;
  for(mnv_itr2 = Response.begin(); mnv_itr2 != Response.end(); ++mnv_itr2){
                (mnv_itr2->second)->GetMigrationObjects( migrationH, h_reco, h_truth );;
        }
//  const bool clear_bands = true; // this was causing the unfolding matrix to not have any systematics 
  const bool clear_bands = false;
  migrationH->PopVertErrorBand("cv"); //causing problems b/c this error band isn't on the background subtracted data and actually turns out to be null-pointer
  mresp = HW2D(migrationH, univs, clear_bands);
}


//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFile(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
//                 reco_QE.hist->Write();
//                 reco_2p2h.hist->Write();
//                 reco_trueDIS.hist->Write();
//                 reco_softDIS.hist->Write();
//                 reco_other.hist->Write();
//                 reco_resonant.hist->Write();
                 signal_purityNum.hist->Write();
                 bkg_total.hist->Write();
//                 bkg_wrongSign.hist->Write();
//                 bkg_NC.hist->Write();
//                 bkg_QE.hist->Write();
//                 bkg_2p2h.hist->Write();
 //                bkg_trueDIS.hist->Write();
 //                bkg_softDIS.hist->Write();
 //                bkg_other.hist->Write();
 //                bkg_resonant.hist->Write();
 //                ANCC.hist->Write(); ANNC.hist->Write(); NCC.hist->Write(); NNC.hist->Write(); ANNCCNCother.hist->Write();
//                 h_background_subtracted_mc.hist->Write();
//                 h_background_mc_scale.hist->Write();
   }//  Response->GetMigrationMatrix()->Write();}   
    else { m_selected_data_reco.hist->Write();
//           h_background_subtracted_data.hist->Write();
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileEff(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
//                 reco_QE.hist->Write();
//                 reco_2p2h.hist->Write();
//                 reco_trueDIS.hist->Write();
//                 reco_softDIS.hist->Write();
//                 reco_other.hist->Write();
//                 reco_resonant.hist->Write();
   }  
    else { m_selected_truth_reco.hist->Write();
//           true_QE.hist->Write();
//           true_2p2h.hist->Write();
//           true_trueDIS.hist->Write();
//           true_softDIS.hist->Write();
//           true_other.hist->Write();
//           true_resonant.hist->Write();
     }
    // selected mc reco
  }
//=======================================================================================
// WRITE ALL HISTOGRAMS
//=======================================================================================
void WriteAllHistogramsToFileMig(TFile& f,bool isMC) const {
    f.cd();
       if(isMC){ m_selected_mc_reco.hist->Write();
                 m_selected_Migration.hist->Write();
//       for (auto responses:Response){
//         responses.second->GetMigrationMatrix()->Write();
//         responses.second->GetReco2D()->Write();
//         responses.second->GetTruth2D()->Write();
//       }   
        mresp.hist->Write();
//        Response->GetMigrationMatrix()->Write();
//        Response->GetReco2D()->Write();
//        Response->GetTruth2D()->Write(); 
   }//  Response->GetMigrationMatrix()->Write();}   
    else {m_selected_data_reco.hist->Write();
     }
  }
};
}  // namespace Var2DLoop


#endif  // VARIABLE_H
