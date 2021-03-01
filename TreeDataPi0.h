#pragma once


#include "BaseOutEvent.h"

#pragma link C++ class pauln::TreeDataPi0;

namespace pauln{

  class TreeDataPi0 : public chanser::BaseOutEvent{

  public:
    TreeDataPi0(){SetName("Pi0");}
    ~TreeDataPi0() final =default;
      

    Double_t no_neut_e_magP = -1;
    Double_t no_neut_phot1_magP = -1;
    Double_t no_neut_phot2_magP = -1;

    Double_t no_neut_e_theta = -1;
    Double_t no_neut_phot1_theta = -1;
    Double_t no_neut_phot2_theta = -1;

    Short_t rec_egg_PID = -1;


    //data member for tree branches below here
    // Int_t helicity = 0;
    // //Int_t helicityonline = 0;

    // Double_t W2=0;
    // Double_t Q2=0;
    // Double_t xB=0;
    // Double_t tneg=0;
    // Double_t tneg_pi0=0;
    // Double_t dtneg=0;
    
    // Double_t IM_g1g2=0;
    // Double_t MM2_total=0;
    // Double_t MP_total=0;
    // Double_t ME_total=0;
    // Double_t MPt_total=0;

    // // Double_t MM2_totalfixed=0;
    // // Double_t MP_totalfixed=0;
    // // Double_t ME_totalfixed=0;
    // // Double_t MPt_totalfixed=0;

    // Double_t pi0coneangle=0;
    // Double_t eg1coneangle=0;
    // Double_t eg2coneangle=0;
    // Double_t recoilconeangle = 0;

    // Double_t MP_rec_recoil   = 0;
    // Double_t MPt_rec_recoil  = 0;
    // Double_t MM_rec_recoil   = 0;
    // Double_t MM2_rec_recoil  = 0;

    // Double_t MP_rec_spectator  = 0;
    // Double_t MPt_rec_spectator = 0;
    // Double_t MM_rec_spectator  = 0;
    // Double_t MM2_rec_spectator = 0;

    // Double_t phi_Nvg = 0;
    // Double_t phi_Nnew = 0;
    // Double_t phi_vgnew = 0;
    // Double_t cop_Nvg_vgnew = 0;
    // Double_t cop_Nvg_Nnew = 0;
    // Double_t cop_Nnew_vgnew = 0;

    // Double_t recon_recoil_T = 0;
    // Double_t recoil_T = 0;
    // Double_t dneutT = 0;


    // ///FS components
    // Double_t rec_px  = 0; 
    // Double_t rec_py  = 0; 
    // Double_t rec_pz  = 0; 
    // Double_t rec_E   = 0; 
    // Double_t rec_magP = 0;
    // Double_t rec_pT  = 0; 
    // Double_t rec_theta  = 0;
    // Double_t rec_phi  = 0;
    Short_t rec_status = 0;
    Short_t rec_PID = 0;

    // Double_t e_px  = 0; 
    // Double_t e_py  = 0; 
    // Double_t e_pz  = 0; 
    // Double_t e_E   = 0; 
    // Double_t e_magP = 0;
    // Double_t e_pT  = 0; 
    // Double_t e_theta  = 0;
    // Double_t e_phi  = 0;
    // Short_t e_status = 0;

    // Double_t phot1_px  = 0;
    // Double_t phot1_py  = 0;
    // Double_t phot1_pz  = 0;
    // Double_t phot1_E   = 0;
    // Double_t phot1_magP = 0;
    // Double_t phot1_pT  = 0;
    // Double_t phot1_theta = 0;
    // Double_t phot1_phi = 0;
    // Short_t phot1_status = 0;

    // Double_t phot2_px  = 0;
    // Double_t phot2_py  = 0;
    // Double_t phot2_pz  = 0;
    // Double_t phot2_E   = 0;
    // Double_t phot2_magP = 0;
    // Double_t phot2_pT  = 0;
    // Double_t phot2_theta = 0;
    // Double_t phot2_phi = 0;
    // Short_t phot2_status = 0;

    // Double_t rec_Beta = 0;
    // Double_t phot1_Beta = 0;
    // Double_t phot2_Beta = 0;

    // Double_t e_sampfrac = 0;
    // Double_t rcdb_Ebeam = 0;


    // //flags
    // Short_t flag_photon1_FT = -1;
    // Short_t flag_photon1_FD = -1;
    // Short_t flag_photon2_FT = -1;
    // Short_t flag_photon2_FD = -1;      
    // Short_t flag_cut_3sigPi0IM = -1;
    // Short_t flag_MC_neutrec = -1;

    ///////////////////////////////////////////////////////////
    //LEAVE THE FOLLOWING FUNCTIONS
    //Function required to set tree branches
    void Branches(TTree* tree) final{
      BaseOutEvent::Branches(tree,Class()->GetListOfDataMembers());
    }
    void Hipo(hipo::ntuple_writer* writer) final{
      BaseOutEvent::Hipo(writer,Class()->GetListOfDataMembers());
    }
      
    ClassDefOverride(TreeDataPi0,1);
  };
}
