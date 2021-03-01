//////////////////////////////////////////////////////////////
///
///Class:		Pi0
///Description:
///
#pragma once


#include "CLAS12FinalState.h"
#include "CLAS12Particle.h"

#include "TreeDataPi0.h"


namespace pauln{

  using Particle=chanser::CLAS12Particle;

  class Pi0 : public chanser::CLAS12FinalState{

       
  public :
    Pi0()=default;
      
   TString GetUSER() final {return _USER;};
 
    //create an instance of the class
    static std::unique_ptr<Pi0> Make(TString ch,TString inc) {
      return std::unique_ptr<Pi0>{new Pi0{ch,inc}};
    }
    //create an instance of the treedata, should be used to init unique_ptr
    chanser::base_outevt_uptr TreeDataFactory() final {
      return chanser::base_outevt_uptr{new TreeDataPi0{}};
    }
    void SetOutEvent(BaseOutEvent* out) final{
      TD=static_cast<TreeDataPi0*>(out);
    }
  
    ~Pi0() final =default;

    void Define() final;
      
    BaseOutEvent* GetOutEvent() noexcept final{return TD;}
    
    void DerivedChangeRun() final {
      //If databases are implemented you can
      //set the beam energy here
      auto ebeam=GetRunInfo()->_BeamEnergy;
      auto mele = 0.00051099891;
      std::cout<<"Change beam energy to :"<<ebeam<<std::endl;
      _beam.SetXYZT(0,0,ebeam,TMath::Sqrt(ebeam*ebeam + mele*mele));
       
    }  

    void calc_angles(TVector3 Ebeam_vect, TVector3 Electron_vect, TVector3 Recoil_vect, TVector3 Newpart_vect);
    TVector3 p3(chanser::HSLorentzVector v);
    Short_t Pi0mass_3sigma_check();
    Double_t Reconstruct_Recoil_KinE();
    Double_t Calc_SamplingFract(Particle& part);

  protected :
    void Kinematics() final;
    void UserProcess() final;
      
      
   
  private:
    //constructor private so only create unique_ptr
    //using Pi0::Make(...)
    //auto fs = pauln::Pi0::Make("NONE","ALL");
  Pi0(TString ch,TString inc) : chanser::CLAS12FinalState(std::move(ch),std::move(inc)){
      //Give object class name - namespace
      //Used for compiling and loading
      SetName(chanser::Archive::BareClassName(ClassName()));
      Define();
    }

    
    //Final Particles Detected
    Particle   _electron = Particle{"e-"};//!
    // Particle   _neutron = Particle{"neutron"};//!
    Particle   _gamma1 = Particle{"gamma"};//!
    Particle   _gamma2 = Particle{"gamma"};//!
    //chanser::CLAS12Particle _PARTICLE=BaseParticle("PDG");//!
    
    //Final Parents
    Particle _pi0 = Particle{"pi0"};//!


    //Initial state
    HSLorentzVector _beam{0,0,10.6,10.6};//!
    HSLorentzVector _target{0,0,0,0.939565420};//!

    //Tree Output Data
    TreeDataPi0* TD{nullptr};//!;

   
    
    const TString _USER="pauln";
    ClassDefOverride(pauln::Pi0,1); //class Pi0
  }; //end Pi0
  
  //utility functions:

  inline Bool_t region_check(const Particle& part, short region) {
    // Checks whether a particle has a hit associate in the specified region
    // and returns true or false.
    auto c12=part.CLAS12(); //if you require other DST data
    return c12->getRegion() == region ?  kTRUE :  kFALSE;
  }

  inline Double_t calc_Beta(const Particle& part){
    Double_t p = part.P4().P();
    Double_t E = part.P4().E();

    return p/E;
  }

}
