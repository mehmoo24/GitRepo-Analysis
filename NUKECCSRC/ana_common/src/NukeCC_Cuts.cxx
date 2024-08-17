//#ifndef MNV_NUKECC_cxx
//#define MNV_NUKECC_cxx 1

#ifndef MNV_NUKECC_CUTS_CXX
#define MNV_NUKECC_CUTS_CXX 1

//#include "include/NukeCC_Cuts.h"
#include "../include/NukeCC_Cuts.h"
#include "../include/NukeCCvars.h"
#include "PlotUtils/TargetUtils.h"
#include "../include/CVUniverse.h"

#include "../include/GlobalIncludes.h" 

#include "../include/CommonIncludes.h"
//#include "../include/NukeCCvars.h"
//#include "CCQENuUtilsNSF.h"
//#include "include/NukeUtils.h"
//#include "Acceptance/TAcceptanceTable.h"
#include <PlotUtils/MnvNormalization.h>
#include <PlotUtils/NuclModUtils.h>
#include <PlotUtils/FluxReweighter.h>
//#include <PlotUtils/FluxReweighterWithWiggleFit.h>
//#include "include/CondorInput.h"
#include <PlotUtils/MnvNuclearModelWeight.h>
//#include <PlotUtils/MnvNormalizerME.h>
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "NukeCCsetbranchaddrs.h"

#include "TFileCollection.h"


#include "TVector3.h"

using namespace std;

using namespace NUKECC_ANA;

NukeCC_Cuts::NukeCC_Cuts()
{
}

NukeCC_Cuts::~NukeCC_Cuts()
{
}

bool NukeCC_Cuts::PassMuCurveCut(CVUniverse *cv, double minCut,HelicityType::t_HelicityType h  )
{
    //! If measured by curvature, make a cut on the significance.
    double Curve = (double)( 1/cv->GetDouble((cv->GetAnaToolName() + "_minos_trk_eqp_qp").c_str()) );
  //  return ( Curve >= minCut );
    if (h==HelicityType::kNeutrino)return ( Curve <= -minCut ); //curve< -5
    if (h==HelicityType::kAntiNeutrino)return ( Curve >= minCut ); //curve>5
    else return (abs(Curve)>=minCut);//|curve|>5
}

bool NukeCC_Cuts::PassMuCurveCut(CVUniverse *cv ,HelicityType::t_HelicityType h){

    //Don't make a significance cut if we reconstructed by range

/* // commented out following for running test	
    if(cv->GetDouble((cv->GetAnaToolName() + "_minos_used_curvature").c_str()) != 1)
        return true;
*/
    return PassMuCurveCut(cv, MIN_MINOS_CURVE, h);
    
}

bool NukeCC_Cuts::PassMuCoilCut(CVUniverse *cv )
{
    
    const double coilXPos = 1219.0;
    const double coilYPos = 393.0;
    const double minos_x = cv->GetDouble((cv->GetAnaToolName() + "_minos_trk_end_x").c_str()) + coilXPos;
    const double minos_y = cv->GetDouble((cv->GetAnaToolName() + "_minos_trk_end_y").c_str()) + coilYPos;
    double minosR = sqrt(pow(minos_x,2) + pow(minos_y,2) );
    
    /*if( !((pow(minos_x,2) + pow(minos_y,2) )>= pow(coilR, 2)) )
     cout << minos_x << " " << minos_y << " " << coilR << endl;
     */
    return (minosR > MINOS_COIL_RADIUS && minosR < MAX_MINOS_RADIUS );
}


bool NukeCC_Cuts::PassMuQualityCut( CVUniverse* cv,int qual /* = 2 */ )
{
    //! If the minos track quality is unknown or not set, return false
    if( cv->GetDouble("MasterAnaDev_minos_trk_quality") <= 0 )
        return false;
    //! If the minos track quality is equal to or better than requested, return true
    return ( cv->GetDouble("MasterAnaDev_minos_trk_quality") <= qual );
}

bool NukeCC_Cuts::PassTrueMuEnergyCut(CVUniverse* cv )
{
    //data always passes
    if( ! isMC )
        return true;
    
    //same as reco limits
    return ( MIN_RECO_E_MU <  cv->GetElepTrue() &&  cv->GetElepTrue() < MAX_RECO_E_MU );
}


    

bool NukeCC_Cuts::PassMuEnergyCut(CVUniverse* cv )
{
    return ( MIN_RECO_E_MU < cv->GetEmu() && cv->GetEmu() < MAX_RECO_E_MU );
}

/*bool CVUniverse::PassMuEnergyCut(CVUniverse* cv )
{
	//cout<<"umaaa"<<endl;
    return PassMuEnergyCut( cv->GetVecElem("MasterAnaDev_leptonE",3)>MIN_RECO_E_MU );
}
*/

bool NukeCC_Cuts::PassTrueThetaCut(CVUniverse* cv)
{
//  return ( //0. <= cv->GetThetalepTrue()*rad_to_deg && 
//cv->GetThetalepTrue()*rad_to_deg < MAX_THETA_MU );
    if (cv->GetThetalepTrue() > 0.349 || TMath::IsNaN(cv->GetThetalepTrue())) return false;
  return true;
}


bool NukeCC_Cuts::PassThetaCut(CVUniverse* cv)
{
// initially I had: 
//  return ( 0. <= cv->GetThetamu()*rad_to_deg && cv->GetThetamu()*rad_to_deg < MAX_RECO_THETA_MU );

//    return ( 0. <= (cv->GetDouble("CCQENu_muon_theta"))*rad_to_deg && (cv->GetDouble("CCQENu_muon_theta"))*rad_to_deg < MAX_RECO_THETA_MU );

// Dan's version of my cut:
//   return ( 0. <= cv->GetDouble("CCQENu_muon_theta") && cv->GetDouble("CCQENu_muon_theta") < 0.349);

  // implement it the way Dan is implementing it
//  if( cv->GetDouble("CCQENu_muon_theta") > 0.349 || TMath::IsNaN(cv->GetDouble("CCQENu_muon_theta")) ) return false;
//  return true;

// UNCOMMENT: 
  if (cv->GetThetamu() > 0.349 || TMath::IsNaN(cv->GetThetamu())) return false;
  return true;

     
    //if (cv->GetThetamu() > 0.349) return false;
    //else return true;
//    return ( 0. <= cv->GetThetamu() && cv->GetThetamu() <= 0.349); //49973
//    return ( 0. <= cv->GetThetamu() && cv->GetThetamu() < 0.349);  //49973

    //return ( 0. <= (cv->GetDouble("CCQENu_muon_theta"))*rad_to_deg && (cv->GetDouble("CCQENu_muon_theta"))*rad_to_deg < MAX_RECO_THETA_MU );

}

bool NukeCC_Cuts::PassGoodTrackingCut(CVUniverse* cv)
{
    //const double MaxUpstreamE = 6.5;
    const double MaxUpstreamE   = 10000000.5;
    const double MinUpstreamE   =  1.5;
    //const int minUpstreamPlanes = 1;
    //const double MinUpstreamE = 0.5;
    const int minUpstreamPlanes = 6;
    if( minUpstreamPlanes <= cv->GetInt("usact_n_planes_tiny") && MinUpstreamE <= cv->GetDouble("usact_avg_E_tiny") && cv->GetDouble("usact_avg_E_tiny") <=MaxUpstreamE ) return false;
    //if( minUpstreamPlanes <= usact_n_planes_tiny && 1.5 <= usact_avg_E_tiny ) return false;
    if( 6 <= cv->GetInt("usact_n_planes_tiny") && 1.5 <= cv->GetDouble("usact_avg_E_tiny") && 6 <= cv->GetInt("muon_n_USclusters") ) return false;
    //muon_n_USclusters == number of clusters TrackAddClusters adds to the track
    
    return true;
}

bool NukeCC_Cuts::passDansCuts(CVUniverse* cv){
    // pass interaction vertex cut
    if (cv->GetInt("has_interaction_vertex") != 1) return false;

    // passDeadTimeCuts
    int tdead_max = 1;
    if( cv->GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj") > tdead_max ) return false;

    // pass NuHelicity Cut
    // nu_helicity_cut(1),// 1 if neutrino, 2 if anti neutrino
    if( cv->GetInt("CCQENu_nuHelicity") != 1 ) return false; //neutrino mode

    // passTrackAngleCut
    if( cv->GetDouble("CCQENu_muon_theta") > 0.349 || TMath::IsNaN(cv->GetDouble("CCQENu_muon_theta")) ) return false;

    return true;
} 


bool NukeCC_Cuts::PassZDistCut(CVUniverse* cv,double  lowZ /*= 1001. */, double highZ /*= 1001.*/ )
{
  double minZ = 4290.; // z center of mod -5 plane 1 is 4293.04mm, while plane 2 is 4313.68mm
  double maxZ = 6000.;
  if (useDNN){
  if ( cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) == 0 && ( minZ <= cv->GetVecElem("ANN_vtx",2) && cv->GetVecElem("ANN_vtx",2) <= maxZ ) ){
    return true;
  }
  if( 1 <= cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) && cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) <= 5 )
  {
     // if we only want fitted vtx events with vertex really in target
     //if( IsMultiTrack() && 1.0E-4 < fabs(CVUniverse_target_zDist) )
     //return false;
    return true;
  }
  if( FIRST_TRACKER_MOD <= cv->GetInt("ANN_vtx_modules") && cv->GetInt("ANN_vtx_modules") <= LAST_TRACKER_MOD )
    return true;

  } 
  else {
    if ( cv->GetInt((cv->GetAnaToolName() + "_targetID").c_str()) == 0 && ( minZ <= cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) && cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) <= maxZ ) ){
      return true;
    }
    if( 1 <= cv->GetInt((cv->GetAnaToolName() + "_targetID").c_str()) && cv->GetInt((cv->GetAnaToolName() + "_targetID").c_str()) <= 5 )
    {
       // if we only want fitted vtx events with vertex really in target
       // if( IsMultiTrack() && 1.0E-4 < fabs(CVUniverse_target_zDist) )
       // return false;
       return true;
    }
    if( FIRST_TRACKER_MOD <= cv->GetInt((cv->GetAnaToolName() + "_vtx_module").c_str()) && cv->GetInt((cv->GetAnaToolName() + "_vtx_module").c_str()) <= LAST_TRACKER_MOD )
        return true;

    }
    // passive target ZDist cut is made with nPlanes cut in framework
/*   if( 1 <= cv->GetInt("NukeCC_ANN_targetID") && cv->GetInt("NukeCC_ANN_targetID") <= 5 )
    {
        // if we only want fitted vtx events with vertex really in target
        //if( IsMultiTrack() && 1.0E-4 < fabs(CVUniverse_target_zDist) )
        // return false;
        return true;
    }
  
    // all events in the tracker pass for faux targets
    if (useDNN){
    if( FIRST_TRACKER_MOD <= cv->GetInt("ANN_vtx_modules") && cv->GetInt("ANN_vtx_modules") <= LAST_TRACKER_MOD )
        return true;
    }
    else {
    if( FIRST_TRACKER_MOD <= cv->GetInt("NukeCC_vtx_module") && cv->GetInt("NukeCC_vtx_module") <= LAST_TRACKER_MOD )
        return true;
    }*/

  return false;
} // end PassZDisCut

bool NukeCC_Cuts::PassTruePTCut(CVUniverse* cv)
{
    if ((cv->GetMuonPTTrue() >= 0.85) && (cv->GetMuonPTTrue() <= 4.5)) return true;
    else return false;
}
bool NukeCC_Cuts::PassRecoPTCut(CVUniverse* cv)
{
    if ((cv->GetMuonPT() >= 0.85) && (cv->GetMuonPT() <= 4.5)) return true;
    else return false;
}


bool NukeCC_Cuts::PassDistToDivisionCut(CVUniverse* cv, double xySep /* = 25.*/ )
{
//NukeCC_targetID=cv->GetInt("NukeCC_targetID"); 
    //only relevant for passive targets 1235
    //if( 0 < cv->GetInt("NukeCC_targetID") && cv->GetInt("NukeCC_targetID") < 10 && 4 != cv->GetInt("NukeCC_targetID") )
    if (useDNN){
    if( 0 < cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) && cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) < 10 && 4 != cv->GetInt((cv->GetAnaToolName() + "_ANN_targetID").c_str()) )      
        return ( xySep < cv->GetDouble((cv->GetAnaToolName() + "_ANN_target_dist_to_division" ).c_str()));   
    }
    else{ 
    if( 0 < cv->GetInt((cv->GetAnaToolName() + "_targetID").c_str()) && cv->GetInt((cv->GetAnaToolName() + "_targetID").c_str()) < 10 && 4 != cv->GetInt((cv->GetAnaToolName() +  "_targetID").c_str()) )
        return ( xySep < cv->GetDouble((cv->GetAnaToolName() +  "_target_dist_to_division" ).c_str()));
    }
    return true;
}

bool NukeCC_Cuts::PassTrueDistToDivisionCut(CVUniverse* cv, double xySep /* = 25. */ )
{
    //only relevant for passive targets 1235
    if( 0 < cv->GetInt("truth_targetID") && cv->GetInt("truth_targetID") < 10 && 4 != cv->GetInt("truth_targetID"))
        return ( xySep < cv->GetDouble("truth_target_dist_to_division") );
    
    return true;
}

bool NukeCC_Cuts::TrackerOnly(CVUniverse* cv)
{
    if (cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) >= 5980 && cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) <= 8422)
        return true;
    else
        return false;
} // added July 6/21

bool NukeCC_Cuts::TrackerOnlyTrue(CVUniverse* cv)
{
    if (cv->GetVecElem("mc_vtx",2) >= 5980 && cv->GetVecElem("mc_vtx",2) <= 8422)
        return true;
    else
        return false;
} // added July 6/21 




/*
void NukeCC_Cuts::Show(CVUniverse *cv,Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
} 
*/
//=======================================================
// Analysis Helpers/Cuts
//=======================================================
bool NukeCC_Cuts::IsInMaterial(CVUniverse* cv,int i_targetID,  int i_targetZ, bool anyTrackerMod /* = false */)
{
    
    if( i_targetID < 0)
    {
        // if targetID < 0, then we want any event in the nuclear target region mods -5-26
        double z = cv->GetVecElem("MasterAnaDev_vtx",2);
        double minZ = 4290.; // z center of mod -5 plane 1 is 4293.04mm, while plane 2 is 4313.68mm
        double maxZ = 6000.;
        if ( !( minZ <= z && z <= maxZ ) )
            return false;
        
        // -targetID is the reference target
        // THIS DOESN'T WORK FOR LOCAL PLASTIC SIDEBAND
        // Instead check if the x,y position is in target material
        /*int refTarg = -i_targetID;
        if( i_targetZ > 0 && CVUniverse_ref_targZ[refTarg-1] != i_targetZ )
          return false;*/

        if( i_targetZ > 0 && ! IsInTargetSection(i_targetID, i_targetZ, cv->GetVecElem("MasterAnaDev_vtx",0),cv->GetVecElem("MasterAnaDev_vtx",1) ) )
          return false;

        return true;
    }
    
    if(i_targetID < 10 )
    {
       //THIS IS A SPECIAL CASE FOR PLASTIC BACKGROUND
      if(cv-> GetInt("MasterAnaDev_targetID") == 0 && cv->GetDouble("MasterAnaDev_vtx_module") < 27 )
      {
	if( i_targetZ > 0 && IsInTargetSection( i_targetID, i_targetZ, cv->GetVecElem("MasterAnaDev_vtx",0),cv->GetVecElem("MasterAnaDev_vtx",1) ) )
	  return true;
      }
      //cout << "targetID < 10" << endl;
      // If targetID < 10, then you are looking for a passive target event.
      // Require that the event has the same targetID and targetZ.
      if( cv->GetInt("MasterAnaDev_targetID") == i_targetID )
      {
	if( i_targetZ > 0 )
	  return cv->GetInt("MasterAnaDev_targetZ") == i_targetZ;
	else
	  return true;
      }//targetID< 10 looking as passive target. Event doesn't have right targetID
    }
    else if( i_targetID < 100 )
    {
        // If 10 < targetID < 100, then we are looking for an event in a plastic reference target.
        // Say targetID = AT, then the event must be in the Ath active target group and the reference target is T.
        int refTarg = i_targetID % 10;
        
        // The starting module of the 4-module reference target is 6*A + 21
        int refModStart = 6*((i_targetID - refTarg )/10) + 21;
        
        int refTargID = refTarg*10000 + 6*1000 + refModStart;
        return IsInMaterial(cv,refTargID, i_targetZ, anyTrackerMod );
    }
    else
    {
        int refTarg = (i_targetID - i_targetID % 10000 ) / 10000;
        if( i_targetZ > 0 && cv->GetVecElem("MasterAnaDev_ref_targZ",(refTarg-1)) != i_targetZ )
            return false;
        
        int refModStart = i_targetID % 1000;
        int refNMod = ( ( i_targetID - refModStart ) / 1000 % 10 );
        int firstMod = refModStart;
        int lastMod  = refModStart + refNMod - 1;
        
        if( anyTrackerMod )
        {
            firstMod = FIRST_TRACKER_MOD;
            lastMod  = LAST_TRACKER_MOD;
        }
        
        // OK if the vertex module is within the range specified
        if( firstMod <= cv->GetDouble("MasterAnaDev_vtx_module") && cv->GetDouble("MasterAnaDev_vtx_module") <= lastMod )
            return true;
    }
    return false;
}

// Passive Target Region Cut (including targets and plastic sidebands) // from Anezka - added May 13 for production work
bool NukeCC_Cuts::PassiveTargetRegion(CVUniverse* cv){
    if (cv->GetVecElem("ANN_vtx",2) >= PlotUtils::TargetProp::NukeRegion::Face && cv->GetVecElem("ANN_vtx",2) <= PlotUtils::TargetProp::NukeRegion::Back)
        return true;
    else
        return false;
}

bool NukeCC_Cuts::PassiveTargetRegionTBV(CVUniverse* cv){
    if (cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) >= PlotUtils::TargetProp::NukeRegion::Face && cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) <= PlotUtils::TargetProp::NukeRegion::Back)
        return true;
    else
        return false;
}

bool NukeCC_Cuts::PassiveTargetRegionTrue(CVUniverse* cv){
    if (cv->GetVecElem("mc_vtx",2) >= PlotUtils::TargetProp::NukeRegion::Face && cv->GetVecElem("mc_vtx",2) <= PlotUtils::TargetProp::NukeRegion::Back)

        return true;
    else
        return false;
} // end of added May 13

bool NukeCC_Cuts::TrackerRegionML(CVUniverse* cv){
    if (cv->GetVecElem("ANN_vtx",2) >= 5980 && cv->GetVecElem("ANN_vtx",2) <= 8422)
        return true;
    else
        return false;
}

bool NukeCC_Cuts::TrackerRegionTBV(CVUniverse* cv){
    if (cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) >= 5980 && cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),2) <= 8422)
        return true;
    else
        return false;
}

bool NukeCC_Cuts::TrackerRegionTrue(CVUniverse* cv){
    if (cv->GetVecElem("mc_vtx",2) >= 5980 && cv->GetVecElem("mc_vtx",2) <= 8422)
        return true;
    else
        return false;
}

bool NukeCC_Cuts::IsInTargetSection( int targetID, int targetZ, double x, double y ) const
{
    // refTarg is targetID for passives
    //     ... targetID%10 for shorthand faux
    int refTarg = targetID % 10;
    
    //     ... and this for scint chunks
    if( targetID > 100 )
        refTarg = (targetID - targetID % 10000 ) / 10000;
    
    //everyone's a winner in target 4!
    if( 4 == refTarg ) return true;
    
    if( 1 == refTarg || 5 == refTarg )
    {
        double u = GetU( x, y );
        //iron is u < 205mm
        if( 26 == targetZ ) return ( u < 205 );
        //lead is 205mm < u
        if( 82 == targetZ ) return ( 205 < u );
    }
    if( 2 == refTarg )
    {
        double d = GetD( x, y );
        //iron is d < 205mm
        if( 26 == targetZ ) return ( d < 205 );
        //lead is 205mm < d
        if( 82 == targetZ ) return ( 205 < d );
    }
    if( 3 == refTarg )
    {
        double c = GetC(x,y);
        //carbon is 0mm < c
        if( 6 == targetZ ) return ( 0 < c );
        //iron is c < 0mm and x < 0
        if( 26 == targetZ ) return ( c < 0 && x < 0 );
        //lead is c < 0mm and 0 < x
        if( 82 == targetZ ) return ( c < 0 && 0 < x );
    }
    
    Error( "NukeCC_Cuts::::IsInTargetSection", Form("Given invalid targetID %d and targetZ %d", targetID, targetZ ) );
    return false;
}

double NukeCC_Cuts::GetU(double x, double y) const 
{
    return -x*cos(PI/6) + y*sin(PI/6); // perp to divide in target 1/5.  up and to the left.
}

double NukeCC_Cuts::GetD(double x, double y) const 
{
    return x*cos(PI/6) + y*sin(PI/6); //points to MINERvA racks (negative x), perp to divide in target 2. down and to the right.
}

double NukeCC_Cuts::GetC(double x, double y) const 
{
    return x*sin(PI/6) + y*cos(PI/6); // points to MINERvA racks (negative x), perp to carbon of target 3
}


int NukeCC_Cuts::GetTargetPlane( int targetID ) const{
    if( targetID == 1 ) return 9;
    else if( targetID == 2 ) return 19;
    else if( targetID == 3 ) return 30;
    else if( targetID == 4 ) return 49;
    else if( targetID == 5 ) return 55;
    else{
        std::cerr << "Bad targetID choice in NukeCC_Cuts::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
    }
    return -999;
}
//bool CVUniverse::IsInMaterial( const MatPair& mat, bool anyTrackerMod /* = false */ )
/*8{
    return IsInMaterial( mat.first, mat.second, anyTrackerMod );
}*/

bool NukeCC_Cuts::IsInTrueMaterial(CVUniverse* cv,const int i_targetID, const int i_targetZ, bool anyTrackerMod /* = false */ )
{
  int truth_targetID=cv->GetInt("truth_targetID");
  int truth_vtx_module=cv->GetInt("truth_vtx_module");
  int truth_targetZ=cv->GetInt("truth_targetZ");
    if(i_targetID < 10 )
    {
        // If targetID < 10, then you are looking for a passive target event.
        // Require that the event has the same targetID and targetZ.
        if( truth_targetID == i_targetID )
        {
            if( i_targetZ > 0 )
                return truth_targetZ == i_targetZ;
            else
                return true;
        }
    }
    else if( i_targetID < 100 )
    {
        // If 10 < targetID < 100, then we are looking for an event in a plastic reference target.
        // Say targetID = AT, then the event must be in the Ath active target group and the reference target is T.
        int refTarg = i_targetID % 10;
        
        // The starting module of the 4-module reference target is 6*A + 21
        int refModStart = 6*((i_targetID - refTarg )/10) + 21;
        
        int refTargID = refTarg*10000 + 6*1000 + refModStart;
        return IsInTrueMaterial( cv,refTargID, i_targetZ, anyTrackerMod );
    }
    else
    {
        int refTarg = (i_targetID - i_targetID % 10000 ) / 10000;
        if( i_targetZ > 0 && truth_ref_targZ[refTarg-1] != i_targetZ )
            return false;
        
        int refModStart = i_targetID % 1000;
        int refNMod = ( ( i_targetID - refModStart ) / 1000 % 10 );
        int firstMod = refModStart;
        int lastMod  = refModStart + refNMod - 1;
        
        if( anyTrackerMod )
        {
            firstMod = FIRST_TRACKER_MOD;
            lastMod  = LAST_TRACKER_MOD;
        }
        
        // OK if the vertex module is within the range specified
        if( firstMod <= truth_vtx_module && truth_vtx_module <= lastMod )
            return true;
        
    }
    return false;
}

/*bool CVUniverse::IsInTrueMaterial( const MatPair& mat, bool anyTrackerMod )
{
    return IsInTrueMaterial( mat.first, mat.second, anyTrackerMod );
}
*/



bool NukeCC_Cuts::PassCCQERecoilCut( CVUniverse* cv,const double  minE /* = MIN_CCQE_RECOIL_E */ )
{
    return minE < cv->GetCCQERecoil();
}


bool NukeCC_Cuts::PassInelasticCut(CVUniverse* cv)
{
    return PassInelasticCut( cv->GetQ2RecoGeV(), cv->GetCCQERecoil() );
}


bool NukeCC_Cuts::PassInelasticCut( double q2, double ccqeRecoil )
{  
    if( ! USE_INEL_CUT )
        return true;
    
    // cut for low_q2 region : [0-1] GeV^2
    if ( q2 < 1000.*1000. )
    {
        if ( ! TRIANGLE_INELCUT )
        {
            if( ccqeRecoil < LOW_Q2_CCQE_RECOIL_E_CUT )
                return false;
        }
        else
        {
            // Eextra = -1.5 q2 + 2.0 (GeV): (q2=0, Eextra=2 GeV), (q2=1, Eextra=0.5)
            if ( ccqeRecoil < ( 2.0*1000.*1000. - 1.5*q2 ) )
                return false;
        }
    }
    else
    {   // cut for high_q2 region >= 1 GeV^2
        if ( ccqeRecoil < HIGH_Q2_CCQE_RECOIL_E_CUT )
            return false;
    }
    
    return true;
}

bool NukeCC_Cuts::PassDISCut(CVUniverse* cv)
{

    return PassDISCut( cv->GetQ2RecoGeV(), cv->GetWRecoGeV() );
}

bool NukeCC_Cuts::PassDISCut(double q2, double W)
{
 if(q2 >= 1 && W >= 2){
        return true;
    }
    else{
        return false;
    }
}

bool NukeCC_Cuts::PassTrueDISCut(CVUniverse* cv)
{
    return  PassTrueDISCut(cv->GetQ2TrueGeV(), cv->GetWTrueGeV());

}


bool NukeCC_Cuts::PassTrueDISCut(double q2, double W)
{
  // if(q2 >= 1.0 && 1.3 <=W &&  W < 4.0){
  if(q2 >= 1 && W >= 2){
        return true;
    }
    else{
        return false;
    }
}
//low W
bool NukeCC_Cuts::PassLowWCut(CVUniverse* cv)
{
    return PassLowWCut(cv->GetWRecoGeV() );}


bool NukeCC_Cuts::PassLowWCut( double W )
{
 if(W < 1.5)
        return true;
    else
        return false;
}


//low W
bool NukeCC_Cuts::PassLowWCutTrue(CVUniverse* cv)
{
    return PassLowWCutTrue(cv->GetWTrueGeV());}


bool NukeCC_Cuts::PassLowWCutTrue( double W )
{
 if(W < 1.5)
        return true;
    else
        return false;
}
//soft DIS
bool NukeCC_Cuts::PassLowQ2Cut(CVUniverse* cv)
{
 return PassLowQ2Cut( cv->GetQ2RecoGeV(), cv->GetWRecoGeV() );
}

bool NukeCC_Cuts::PassLowQ2Cut( double q2, double W )
{

//if(q2 >= MIN_DIS_Q2 && W >= MIN_DIS_W)
if(q2 < 1 && W >= 2)
        return true;
    else
        return false;
}

// True
bool NukeCC_Cuts::PassLowQ2CutTrue(CVUniverse* cv)
{
    return PassLowQ2Cut(cv->GetQ2TrueGeV(), cv->GetWTrueGeV()  );
}

bool NukeCC_Cuts::PassLowQ2CutTrue( double q2, double W )
{

if(q2 < 1 && W >= 2.0)
        return true;
    else
        return false;
}
// lowQ2Trans
bool NukeCC_Cuts::PassLowQ2Trans(CVUniverse* cv)
{
    return PassLowQ2Trans(cv->GetQ2RecoGeV(), cv->GetWRecoGeV());
}
bool NukeCC_Cuts::PassLowQ2Trans( double q2, double W )
{
if(q2 < 1.0 && 1.5 <= W && W < 2.0)
        return true;
    else
        return false;
}


// lowQ2Trans True
bool NukeCC_Cuts::PassLowQ2TransTrue(CVUniverse* cv)
{

    return PassLowQ2TransTrue(cv->GetQ2TrueGeV(), cv->GetWTrueGeV());
}
bool NukeCC_Cuts::PassLowQ2TransTrue( double q2, double W )
{
if(q2 < 1.0 && 1.5 <= W && W < 2.0)
        return true;
    else
        return false;
}


//trans Reco
bool NukeCC_Cuts::Passtrans(CVUniverse* cv)
{
    return Passtrans(cv->GetQ2RecoGeV(), cv->GetWRecoGeV() );
}
bool NukeCC_Cuts::Passtrans( double q2, double W )
{
if(q2 > 1.0 && 1.5 <= W && W < 2.0)
        return true;
    else
        return false;
}
// True
bool NukeCC_Cuts::PasstransTrue(CVUniverse* cv)
{
    return PasstransTrue(cv->GetQ2TrueGeV(), cv->GetWTrueGeV());

}
bool NukeCC_Cuts::PasstransTrue( double q2, double W )
{
if(q2 >= 1.0 && 1.5 <= W && W < 2.0)
        return true;
    else
        return false;
}
bool NukeCC_Cuts::InHexagon(CVUniverse* cv, double apothem /*= 850.*/ ){
  
    //vector<double> newVertex = GetMLVertex();
   double x, y;
   if (useDNN){
    x = cv->GetVecElem("ANN_vtx",0);
    y = cv->GetVecElem("ANN_vtx",1);
   }
   else {
    x = cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),0);
    y = cv->GetVecElem((cv->GetAnaToolName() + "_vtx").c_str(),1);
   }
    if( pow(x,2) + pow(y,2) < pow(apothem,2) )
        return true;
   
    //Hexagon is symmetric about its x and y
    x = fabs(x);
    y = fabs(y);

    double lenOfSide = apothem * ( 2 / sqrt(3) );

    if( x > apothem )
        return false;

    if( y < lenOfSide/2.0 )
        return true;

    double slope = (lenOfSide / 2.0) / apothem;
    if( y < lenOfSide - x*slope )
        return true;

    return false; 
}



bool NukeCC_Cuts::PassHelicityCut(CVUniverse* cv,HelicityType::t_HelicityType h )
{

   int helicity = cv->GetInt((cv->GetAnaToolName() + "_nuHelicity").c_str());
// commented out below lines
 //  if(h==0) return true;
//   else
     return HelicityType::t_HelicityType(helicity) == h;
}

bool NukeCC_Cuts::PassReco(CVUniverse* cv,HelicityType::t_HelicityType h)
{
    if( ! PassHelicityCut( cv,h) )   return false; //! Is this event the right helicity?
//    if( ! cv->GetInt("MasterAnaDev_in_fiducial_area") ) return false; //! Is the event in the fiducial area?
//    if( ! PassZDistCut( cv) )   return false; //! Is the event vertex close enough in Z to the NuclearTarget?
//    if( ! PassDistToDivisionCut(cv) ) return false; //! Is the event vertex far enough away from the separation of target sections?
//    if( 120. < cv->GetEnu() * mev_to_gev ) return false; //! Anne's neutrino Energy cut
    if( 1 < cv->GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj") ) return false; //! Gabe's tdead cut

// put this cut back in
    if( ! PassMuCurveCut(cv, h) ) return false; //! Is the curvature significant?
    //if( ! PassMuCoilCut(cv) ) return false; //! Does the muon track end outside of the coil, or inside the area of MINOS which gives us good energy reco?    

    // if( !IsGoodTarget3Event( ) ) return false; //! Is this a target 3 event which clearly originated in target 1 or 2?
   
   // added this in (Dec 31,2023)
   if (cv->GetInt("isMinosMatchTrack") != 1 ) return false; // it is just a integer that is equal to 1 when the muon match in minos

   // fail minos match cut:
   //if (cv->GetInt("isMinosMatchTrack") == 1 ) return false; 

    return true;
}

//cuts Zubair is using
bool NukeCC_Cuts::PassRecoZ(CVUniverse* cv,HelicityType::t_HelicityType h)
{
  if(! PassHelicityCut(cv,h))   
    return false; 
  if (useDNN){
    if(! cv->GetInt((cv->GetAnaToolName() + "_ANN_in_fiducial_area").c_str()))
      return false; // event in fiducial area?
  }
  else{
    if(! cv->GetInt((cv->GetAnaToolName() + "_in_fiducial_area").c_str()))
      return false; // event in fiducial area?
  }
  if(! InHexagon(cv, 850.))  
    return false; // event in fiducial area?
  if(! PassZDistCut(cv))   
    return false; // event close enough in z to nuclear target?
  if(! PassDistToDivisionCut(cv)) 
    return false; // event vtx far enough away from the separation of target sections?
  if(120. < cv->GetEnu() * mev_to_gev) 
    return false; // Anne's neutrino energy cut
  if(1 < cv->GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj")) 
    return false; // Gabe's tdead cut
  if(! PassMuCurveCut(cv, h)) 
    return false; // curvature significant?
  if(! PassMuCoilCut(cv)) 
    return false; // Does the muon track end outside of the coil, or inside the area of MINOS which gives us good energy reco? 
  return true;
} // end PassRecoZ

bool NukeCC_Cuts::PassTruth(CVUniverse* cv, HelicityType::t_HelicityType h)
{
    if( ! PassTrueCC(cv,h) ) return false; //true nu_mu (bar) CC event?
    if( ! PassTrueFiducial(cv ) ) return false; //true fiducial volume cut
    
    return true;
}

bool NukeCC_Cuts::PassTrueCC(CVUniverse* cv, HelicityType::t_HelicityType h /* = HelicityType::kNeutrino */ )
{
    if( 1 != cv->GetInt("mc_current") ) return false;
    if( ! PassTrueHelicityCut(cv, h ) ) return false;
    
    return true;
}

bool NukeCC_Cuts::PassTrueFiducial(CVUniverse* cv )
{
    if( ! InHexagonTrue(cv, 850.) ) return false;
    if( ! PassTrueDistToDivisionCut(cv ) ) return false;
    
    return true;
}

bool NukeCC_Cuts::PassTrueInelasticCut(CVUniverse* cv)
{
    if( 1 == cv->GetInt("mc_intType") ) return false;
    if( 2 == cv->GetInt("mc_intType") ) return false;
	 
     return true;
}

bool NukeCC_Cuts::PassTrueHelicityCut(CVUniverse* cv,HelicityType::t_HelicityType h )
{
 
  int mc_incoming=cv->GetInt("mc_incoming");
  if(HelicityType::kAny == h)
     return true;
  else if (HelicityType::kNeutrino == h)
     return mc_incoming ==14;
  else if (HelicityType::kAntiNeutrino == h)
     return mc_incoming == -14;
  else{
     Warning( "NukeCC_Cuts::PassTrueHelicityCut", "Helicity type unknown, should not have compiled.  Returning false..." );
     return false;
  }
}

bool NukeCC_Cuts::InCircleTrue(CVUniverse* cv, double radius){
   double x = cv->GetVecElem("mc_vtx",0); 
   double y = cv->GetVecElem("mc_vtx",1);
  
   if( pow(x,2) + pow(y,2) < pow(radius,2) )
        return true;
   return false;
}


bool NukeCC_Cuts::InHexagonTrue(CVUniverse* cv, double apothem /*= 850.*/ ){
    
    //vector<double> newVertex = GetMLVertex();
    double x = cv->GetVecElem("mc_vtx",0);
    double y = cv->GetVecElem("mc_vtx",1);
    
    if( pow(x,2) + pow(y,2) < pow(apothem,2) )
        return true;
    
    //Hexagon is symmetric about its x and y
    x = fabs(x);
    y = fabs(y);
    
    double lenOfSide = apothem * ( 2 / sqrt(3) );
    
    if( x > apothem )
        return false;
    
    if( y < lenOfSide/2.0 )
        return true;
    
    double slope = (lenOfSide / 2.0) / apothem;
    if( y < lenOfSide - x*slope )
        return true;
    
    return false;
}

// hexagonal rings for flux study, sept 15/23
bool NukeCC_Cuts::PassTrueRing(CVUniverse* cv,double apothem_outer, double apothem_inner ){
   bool inBigHex = InHexagonTrue(cv,apothem_outer);
   bool inLittleHex = InHexagonTrue(cv,apothem_inner);
   if( inBigHex && !inLittleHex)
      return true;
   
   return false;
}

bool NukeCC_Cuts::PassRing(CVUniverse* cv,double apothem_outer, double apothem_inner ){
   bool inBigHex = InHexagon(cv,apothem_outer);
   bool inLittleHex = InHexagon(cv,apothem_inner);
   if( inBigHex && !inLittleHex)
      return true;

   return false;
}

bool NukeCC_Cuts::passTrueCCQE( CVUniverse* cv ){
  return passTrueCCQE( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetInt("mc_charm") );
}

bool NukeCC_Cuts::passTrueCCQE( int pdg, int current, int type, bool charm ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 1 && !charm );
  else return ( pdg == -14 && current == 1 && type == 1 && !charm );
}


bool NukeCC_Cuts::passTrueCCRES(CVUniverse* cv){
  return passTrueCCRES(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"));

}

bool NukeCC_Cuts::passTrueCCRES( int pdg, int current, int type ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 2 );
  else  return ( pdg == -14 && current == 1 && type == 2 );
}

bool NukeCC_Cuts::PassEhadCut(CVUniverse* cv){
 if (cv->GetEhadGeV() < 0.2) return true;
 else return false; 
}

bool NukeCC_Cuts::passTrueMEC(CVUniverse* cv){
  return passTrueMEC(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"),cv->GetBool("mc_charm"));
}


bool NukeCC_Cuts::passTrueMEC( int pdg, int current, int type, bool charm ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 8 && !charm );
  else  return ( pdg == -14 && current == 1 && type == 8 && !charm );
}

bool NukeCC_Cuts::passTrueCoh(CVUniverse* cv){
  return passTrueCoh(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetBool("mc_intType"));

}

bool NukeCC_Cuts::passTrueCoh( int pdg, int current, int type ){ 
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 4 ); 
  else return ( pdg == -14 && current == 1 && type == 4 ); //antinu 
}   

bool NukeCC_Cuts::passTrueCCDIS(CVUniverse *cv){
  return passTrueCCDIS(cv->GetInt("mc_incoming"),cv->GetInt("mc_current"),cv->GetInt("mc_intType"));

}

bool NukeCC_Cuts::passTrueCCDIS( int pdg, int current, int type ){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 );
  else  return ( pdg == -14 && current == 1 && type == 3 );
}


bool NukeCC_Cuts::passTrueCCTrueDIS( CVUniverse* cv ){
  return passTrueCCTrueDIS( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetQ2True(), cv->GetWTrue() );
}

bool NukeCC_Cuts::passTrueCCTrueDIS( int pdg, int current, int type, double Q2, double W){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 && Q2>=1000 && W >=2000);
  else  return ( pdg == -14 && current == 1 && type==3 && Q2>=1000 && W >=2000);
}

bool NukeCC_Cuts::passTrueCCTrueSIS( CVUniverse* cv ){
  return passTrueCCTrueSIS( cv->GetInt("mc_incoming"), cv->GetInt("mc_current"), cv->GetInt("mc_intType"), cv->GetQ2True(), cv->GetWTrue() );
}

bool NukeCC_Cuts::passTrueCCTrueSIS( int pdg, int current, int type, double Q2, double W){
  if(neutrinoMode) return ( pdg == 14 && current == 1 && type == 3 && !(Q2>=1000 && W >=2000));
  else  return ( pdg == -14 && current == 1 && type == 3 && !(Q2>=1000 && W >=2000));
}
#endif 

