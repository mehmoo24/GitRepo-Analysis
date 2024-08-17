#ifndef VARIABLE_H
#define VARIABLE_H

#include "CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"

#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#endif  // __CINT__

namespace Ben {

class Variable : public PlotUtils::VariableBase<CVUniverse> {
 private:
  typedef PlotUtils::HistWrapper<CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;
  typedef PlotUtils::HistFolio<PlotUtils::MnvH1D> FOLIO;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable(ARGS... args) : PlotUtils::VariableBase<CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // Histwrappers -- selected mc, selected data
  HW m_selected_mc_reco;
  HW m_selected_data;

  // Histofolio to categorize MC by interaction channel
  FOLIO m_selected_mc_by_channel;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeMCHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    m_selected_mc_reco = HW(Form("selected_mc_reco_%s", name), name, GetNBins(), bins, univs);

    //Interaction types from https://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/classMinerva_1_1GenMinInteraction.html
    std::vector<PlotUtils::NamedCategory<int>> genieCategories = {{1, "QE"}, {2, "RES"}, {3, "DIS"}, {4, "Coherent Pi"}, {8, "MEC"}};
    m_selected_mc_by_channel = FOLIO(genieCategories, std::string("selected_mc_in_") + name + "by_channel", GetNBins(), bins.data());
  }

  template <typename T>
  void InitializeDataHistograms(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    m_selected_data = HW(Form("selected_data_%s", name), name, GetNBins(), bins, univs);
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f) const {
    f.cd();

    // selected mc reco
    if(m_selected_mc_reco.hist) m_selected_mc_reco.hist->Write();

    // selected mc reco - signal background histfolio
    if(m_selected_data.hist) m_selected_data.hist->Write();

    // selected mc reco broken down by GENIE category
    m_selected_mc_by_channel.WriteToFile(f);
  }

  //=======================================================================================
  // SYNC ALL HISTOGRAMS
  //=======================================================================================
  void SyncAllHists() {
    m_selected_mc_reco.SyncCVHistos();
    m_selected_data.SyncCVHistos();
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
};

}  // namespace Ben

#endif  // VARIABLE_H
