#include "Pi0.h"
 
namespace pauln{

  //Bool_t region_check(const Particle& part, short region);
  
  ///////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$//////////////////////  
  void Pi0::Define(){
    //Set final state detected particles
    //Note if particle is added to final with a valid genID it will be used
    //to determine the correct permutation of the simulated event
    
                                                // GenePi LUND looks like:
    AddParticle("Electron",&_electron,kTRUE,0); // 1 -1 1    11 (scattered E)
    AddParticle("Neutron",&_neutron,kTRUE,2);   // 2 -1 1  2112 (spectator)
    AddParticle("Gamma1",&_gamma1,kTRUE,3);     // 3 -1 1  2212 (recoil)
    AddParticle("Gamma2",&_gamma2,kTRUE,4);     // 4 -1 1    22 (photon)
                                                // 5 -1 1    22 (photon)
                                                
    //AddParticle("Name",particle,true/false you want to write in final vector, genID for linking to generated truth value)
    // AddParticle("PARTICLE",&_PARTICLE,kTRUE,-1);

    //Set final state parents
    AddParticle("Pi0",&_pi0,kTRUE,-1);
    ConfigParent(&_pi0,&_gamma1);
    ConfigParent(&_pi0,&_gamma2);


    //_doToTopo["TOPOLOGY"]=[&](){
      //TOPOLOGY Define your topology dedendent code in here
      ///////+++++++++++++++++++++++++++++++++++///////
    ////auto miss= _beam + _target - _electron.P4() - _proton.P4()
    ////  -_pip1.P4()-_pim1.P4();
    ////TD->MissMass=miss.M();
    ////TD->MissMass2=miss.M2();
    
      ///////------------------------------------///////
    // };

    //Set Possible Topologies
    _doToTopo["Electron:Neutron:Gamma1:Gamma2"]=[&](){
      //TOPOLOGY Define your topology dedendent code in here
      ///////+++++++++++++++++++++++++++++++++++///////

      // Currently no topology specific code, as only one topology!
      // -- maybe implement something like:
      //    Particle recoil = _neutron
      //    where another Topology could use  Particle recoil = _proton for pHEMP?

       
    }; // close _doToTopo
      
  } // close Define
 
  ///////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$//////////////////////  
  void Pi0::Kinematics(){
    //Define reaction specific kinematic calculations here
    //Assign to tree data TD.var=

    auto gg_p4 = _gamma1.P4() + _gamma2.P4(); //reconstructed pi0
    TD->IM_g1g2 = gg_p4.M();

    //Thin out ouput: remove any events with very broad (~6sigma) pi0-mass cut
    //if (TD->IM_g1g2 < 0.466 || TD->IM_g1g2 > 0.975) {RejectEvent(); return;}

    //Use Kinematics to calculate electron variables
    //Note this assumes you called your electron "electron" or "Electron"
    _kinCalc.SetElecsTarget(_beam,_electron.P4(),_target);
    TD->W2 = (_neutron.P4() + (_gamma1.P4() + _gamma2.P4())).M2(); 
    TD->Q2 = _kinCalc.Q2();
    TD->xB = TD->Q2 / (2 * _target.Dot(_kinCalc.Gamma()));

    //TD->tneg=-1*_kinCalc.t(_neutron().P4(), _target); <-- DIFFERENT T!!
    // Double_t tneg = -1 * (_neutron.P4() - _target).M2();
    // Double_t tneg_pi0 = -1 * (_beam - _electron.P4() - _gamma1.P4() - _gamma2.P4()).M2();
    TD->tneg = -1 * (_neutron.P4() - _target).M2();
    TD->tneg_pi0 = -1 * (_beam - _electron.P4() - _gamma1.P4() - _gamma2.P4()).M2();
    TD->dtneg = TD->tneg - TD->tneg_pi0;

    using ROOT::Math::VectorUtil::Angle;

    Double_t DEG = TMath::RadToDeg(); //rad->deg conversion

    TD->rcdb_Ebeam = _beam.E();

    TD->helicity =  GetEventInfo()->_BeamHel;

    TD->rec_px = _neutron.P4().Px();
    TD->rec_py = _neutron.P4().Py();
    TD->rec_pz = _neutron.P4().Pz();
    TD->rec_E  = _neutron.P4().E();
    TD->rec_magP = _neutron.P4().P();
    TD->rec_pT = TMath::Sqrt(_neutron.P4().Perp2());
    TD->rec_theta = DEG*_neutron.P4().Theta();
    TD->rec_phi   = DEG*_neutron.P4().Phi();
    TD->rec_status = _neutron.CLAS12()->par()->getStatus();
    TD->rec_PID = _neutron.Truth()->_pdgCode;

    TD->e_px = _electron.P4().Px();
    TD->e_py = _electron.P4().Py();
    TD->e_pz = _electron.P4().Pz();
    TD->e_E  = _electron.P4().E();
    TD->e_magP = _electron.P4().P();
    TD->e_pT = TMath::Sqrt(_electron.P4().Perp2());
    TD->e_theta = DEG*_electron.P4().Theta();
    TD->e_phi   = DEG*_electron.P4().Phi();
    TD->e_status = _electron.CLAS12()->par()->getStatus();

    TD->phot1_px = _gamma1.P4().Px();
    TD->phot1_py = _gamma1.P4().Py();
    TD->phot1_pz = _gamma1.P4().Pz();
    TD->phot1_E  = _gamma1.P4().E();
    TD->phot1_magP = _gamma1.P4().P();
    TD->phot1_pT = TMath::Sqrt(_gamma1.P4().Perp2());
    TD->phot1_theta = DEG*_gamma1.P4().Theta();
    TD->phot1_phi   = DEG*_gamma1.P4().Phi();
    TD->phot1_status = _gamma1.CLAS12()->par()->getStatus();

    TD->phot2_px = _gamma2.P4().Px();
    TD->phot2_py = _gamma2.P4().Py();
    TD->phot2_pz = _gamma2.P4().Pz();
    TD->phot2_E  = _gamma2.P4().E();
    TD->phot2_magP = _gamma2.P4().P();
    TD->phot2_pT = TMath::Sqrt(_gamma2.P4().Perp2());
    TD->phot2_theta = DEG*_gamma2.P4().Theta();
    TD->phot2_phi   = DEG*_gamma2.P4().Phi();
    TD->phot2_status = _gamma2.CLAS12()->par()->getStatus();

    TD->rec_Beta = calc_Beta(_neutron);
    TD->phot1_Beta = calc_Beta(_gamma1);
    TD->phot2_Beta = calc_Beta(_gamma2);


    auto total_p4 = _beam + _target - (_electron.P4() + _neutron.P4() + gg_p4);
    TD->MM2_total =  total_p4.M2();
    TD->MP_total  =  total_p4.P();
    TD->ME_total  =  total_p4.E();
    TD->MPt_total =  TMath::Sqrt(total_p4.Perp2());

    // _pi0.FixP4(gg_p4); //recalc pi0-p4 with fixed PDG pi0 mass
    // auto totalfixed_p4 = _beam + _target - (_electron.P4() + _neutron.P4() + _pi0.P4());

    // TD->MM2_totalfixed  = totalfixed_p4.M2();
    // TD->MP_totalfixed   = totalfixed_p4.P();
    // TD->ME_totalfixed   = totalfixed_p4.E();
    // TD->MPt_totalfixed  = TMath::Sqrt(totalfixed_p4.Perp2());

    auto rec_pi0 = (_beam + _target) - (_electron.P4() + _neutron.P4());
    TD->pi0coneangle = Angle(rec_pi0.Vect(), gg_p4.Vect()) * DEG;

    TD->eg1coneangle = Angle(_electron.P4().Vect(), _gamma1.P4().Vect()) * DEG;
    TD->eg2coneangle = Angle(_electron.P4().Vect(), _gamma2.P4().Vect()) * DEG;

    auto rec_recoil = (_beam + _target) - (_electron.P4() + gg_p4);
    TD->MP_rec_recoil     = rec_recoil.P();
    TD->MPt_rec_recoil    = TMath::Sqrt(rec_recoil.Perp2());
    TD->MM_rec_recoil     = rec_recoil.M();
    TD->MM2_rec_recoil    = rec_recoil.M2();
    TD->recoilconeangle   = (Angle(rec_recoil.Vect(), _neutron.P4().Vect())) * DEG;

    HSLorentzVector deut{0, 0, 0, 1.876};
    auto rec_spectator = (_beam + deut) - (_electron.P4() +_neutron.P4() + gg_p4);
    TD->MP_rec_spectator  = rec_spectator.P();
    TD->MPt_rec_spectator = TMath::Sqrt(rec_spectator.Perp2());
    TD->MM_rec_spectator  = rec_spectator.M();
    TD->MM2_rec_spectator = rec_spectator.M2();


    /* This should all be tweaked into a class analogous to kinematics.   */
    /* "Set" required vectors to pass into clas, and then return phi/cops */
    TVector3 beamp3 = p3(_beam);
    TVector3 electronp3 = p3(_electron.P4());
    TVector3 recoilp3 = p3(_neutron.P4());
    TVector3 gg_p3 = p3(gg_p4);

    calc_angles(beamp3, electronp3, recoilp3, gg_p3);
    /* ------------------------------------------------------------------ */

    //****----------------****//
    //****  Various Flags ****//
    //****----------------****//

    //Photon detection region:
    //region_flags(_gamma1, TD->flag_photon1_FD, TD->flag_photon2_FT);
    TD->flag_photon1_FD = region_check(_gamma1, clas12::FD);
    TD->flag_photon2_FD = region_check(_gamma2, clas12::FD);
    TD->flag_photon1_FT = region_check(_gamma1, clas12::FT);
    TD->flag_photon2_FT = region_check(_gamma2, clas12::FT);

    ///////------------------------------------///////

    TD->flag_cut_3sigPi0IM = Pi0mass_3sigma_check();

    TD->recoil_T = _neutron.P4().E() - 0.939565420; //target assumed at rest (E = m)
    TD->recon_recoil_T = Reconstruct_Recoil_KinE();
    TD->dneutT = TD->recon_recoil_T - TD->recoil_T;

    TD->e_sampfrac = Calc_SamplingFract(_electron);

    if (_neutron.Truth()->_pdgCode==2212){
      TD->flag_MC_neutrec = 0;
    }
    else if (_neutron.Truth()->_pdgCode==2112){
      TD->flag_MC_neutrec = 1;
    }
  }
    
  ///////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$//////////////////////  
  void Pi0::UserProcess(){
    //Optional additional steps
    //This is called after the Topo function
    //and before the kinematics function
    _counter++;


    //Must call the following to fill Final trees etc.
    //Fill Final data output at the end
    FinalState::UserProcess();

  }

  void Pi0::calc_angles(TVector3 Ebeam_vect, TVector3 Electron_vect, TVector3 Recoil_vect, TVector3 Newpart_vect)
  {
    
    Double_t DEG = TMath::RadToDeg(); //rad->deg conversion

    TVector3 Virtual_photon; // virtual photon vector
    TVector3 vect_ee;		     // vector in the leptonic plane
    TVector3 vect_Nnew;		   // three vectors in the hadronic plane
    TVector3 vect_Nvg;
    TVector3 vect_vgnew;

    Virtual_photon = Ebeam_vect - Electron_vect; // vector of the virtual photon

    // Leptonic plane:
    vect_ee = Ebeam_vect.Cross(Electron_vect); // cross product of e and e', gives normal to leptonic plane

    // Hadronic plane:
    vect_Nnew = Recoil_vect.Cross(Newpart_vect);	 // normal to hadronic plane formed by recoiling nucleon and new particle (photon, meson)
    vect_Nvg = Recoil_vect.Cross(Virtual_photon);	 // normal to hadronic plane formed by recoiling nucleon and virtual photon
    vect_vgnew = Virtual_photon.Cross(Newpart_vect); // gives normal to hadronic plane formed by virtual photon and new particle (photon, meson)

    // Phi angles (between leptonic and hadronic planes):
    // Phi convention here:
    //- e from bottom left to center, e' from center to top left (so vect_ee is up out of page)
    //- gamma* from center left to center right
    //- N' from center to bottom right, new particle (gamma or pi0) from center to top right (so vect_Nnew, vect_Nvg, vect_vgnew all up out of page)
    //
    //- phi is 0 when both leptopnic and hadronic planes are parallel
    //- phi increases as you rotate hadronic plane clockwise, if looking along trajectory of virtual photon
    //- phi is 180 deg when the two planes are again parallel but the N' and the new particle (photon or meson) have swapped directions
    //- phi continues to increase to 360 deg as you continue to rotate hadronic plane clockwise, from point of view of gamma*

    Double_t phi_Nvg = (vect_ee.Angle(vect_Nvg)) * DEG; // angle between (e,e') and (N,gamma*) planes, in degrees.
    if ((vect_ee.Dot(Recoil_vect)) > 0.)
      TD->phi_Nvg = 360. - phi_Nvg; // To make sure it runs from 0 - 360 degrees.
    else 
      TD->phi_Nvg = phi_Nvg;

    Double_t phi_Nnew = (vect_ee.Angle(vect_Nnew)) * DEG; // angle between (e,e') and (N,gamma) or (N,pi0) planes, in degrees.
    if ((vect_ee.Dot(Newpart_vect)) < 0.)
      TD->phi_Nnew = 360. - phi_Nnew; // Note: the pi0 or gamma are called "New" for simplicity.
    else 
      TD->phi_Nnew = phi_Nnew;

    Double_t phi_vgnew = (vect_ee.Angle(vect_vgnew)) * DEG; // angle between (e,e') and (gamma*,gamma) or (gamma*,pi0) planes, in degrees.
    if ((vect_ee.Dot(Newpart_vect)) < 0.)
      TD->phi_vgnew = 360. - phi_vgnew; // Note: the pi0 or gamma are called "New" for simplicity.
    else 
      TD->phi_vgnew = phi_vgnew;
    // Coplanarity angles (to make sure the recoil nucleon, the virtual photon and the photon/meson are in the same plane):
    // In the naming convention, "new" refers to produced gamma or pi0.

    Double_t cop_Nvg_vgnew = (vect_Nvg.Angle(vect_vgnew)) * DEG; // angle between (N,gamma*) and (gamma*,gamma) or (gamma*,pi0) planes, in degrees.
    if ((vect_Nvg.Dot(Newpart_vect)) < 0)
      TD->cop_Nvg_vgnew = -1 * cop_Nvg_vgnew; // sort of arbitrary, but makes sure there's a nice peak at zero
    else
      TD->cop_Nvg_vgnew = cop_Nvg_vgnew;

    Double_t cop_Nvg_Nnew = (vect_Nvg.Angle(vect_Nnew)) * DEG; // angle between (N,gamma*) and (N,gamma) or (N,pi0) planes, in degrees.
    if ((vect_Nvg.Dot(Newpart_vect)) < 0)
      TD->cop_Nvg_Nnew = -1 * cop_Nvg_Nnew; // sort of arbitrary, but makes sure there's a nice peak at zero
    else
      TD->cop_Nvg_Nnew = cop_Nvg_Nnew;

    Double_t cop_Nnew_vgnew = (vect_Nnew.Angle(vect_vgnew)) * DEG; // angle between (N,gamma) and (gamma*,gamma) or between (N,pi0) and (gamma*,pi0) planes.
    if ((vect_vgnew.Dot(Recoil_vect)) < 0)
      TD->cop_Nnew_vgnew = -1 * cop_Nnew_vgnew; // sort of arbitrary, but makes sure there's a nice peak at zero
    else
      TD->cop_Nnew_vgnew = cop_Nnew_vgnew;

  } //close calc_angles 

  TVector3 Pi0::p3(chanser::HSLorentzVector v){
    TVector3 v3;
    v3.SetX(v.X());
    v3.SetY(v.Y());
    v3.SetZ(v.Z());

    return v3;
  }

  Short_t Pi0::Pi0mass_3sigma_check(){
    //###################################################################################
    // Set 3sigma IM_gg cut flags.  (FD/FT split)
    // [IM_gg fitted with cuts: MM2 +/- 0.5GeV && pi0coneangle < 20deg (somewhat arb.)]
    //###################################################################################

    /* FD/FT split (linear BG) 3 sigma cut gives the following:
            // ---- both FD: Mean: 0.13100  Sig: 0.01406 =>  Lower: 0.08882   Upper: 0.17318
            // ---- both FT: Mean: 0.13110  Sig: 0.00642 =>  Lower: 0.11184   Upper: 0.15036
            // ---- 1FD/1FT: Mean: 0.12610  Sig: 0.01483 =>  Lower: 0.08161   Upper: 0.17059

            Defined using only 2-photon events:  <<-- CURRENTLY USING
            // ---- both FD: Mean: 0.130879  Sig: 0.0125521  Lower: 0.0932222  Upper: 0.168535
            // ---- both FT: Mean: 0.13331  Sig: 0.00468502  Lower: 0.119254  Upper: 0.147365
            // ---- 1FD/1FT: Mean: 0.128562  Sig: 0.0103054  Lower: 0.097646  Upper: 0.159478
    */
    const Double_t mass = TD->IM_g1g2;
    Double_t mean = 0;
    Double_t sigma = 0;

    if ((TD->flag_photon1_FD == 1) && (TD->flag_photon2_FD == 1)){
      mean = 0.130879;
      sigma = 0.01406;
    } 
    else if ((TD->flag_photon1_FT == 1) && (TD->flag_photon2_FT == 1)){
      mean = 0.13331;
      sigma = 0.00642;
    }
    else if (((TD->flag_photon1_FD == 1) && (TD->flag_photon2_FT == 1)) 
            || ((TD->flag_photon1_FT == 1) && (TD->flag_photon2_FD == 1))){
      mean = 0.128562;
      sigma = 0.0103054;
    }

    Double_t upper = mean + 3 * sigma;
    Double_t lower = mean - 3 * sigma;

    return ((mass <= upper) && (mass >= lower)) ? 1 : 0;

  } // close Pi0mass_3sigma_check

  Double_t Pi0::Calc_SamplingFract(Particle& part){

    auto c12=part.CLAS12();
    if(c12->getRegion()!=clas12::FD) return true; //cut only applies to FD

    double part_p = part.P4().P();
    double ECIN_en = c12->cal(clas12::ECIN)->getEnergy();
    double ECOUT_en = c12->cal(clas12::ECOUT)->getEnergy();
    double PCAL_en = c12->cal(clas12::PCAL)->getEnergy();
    double total_CAL_edep = ECIN_en + ECOUT_en + PCAL_en;

    return (total_CAL_edep/part_p);

  }

  Double_t Pi0::Reconstruct_Recoil_KinE(/*Double_t beamE, HSLorentzVector product_p4, HSLorentzVector recoil_p4,
              Double_t mTarget, Double_t mRecoil, Double_t mSpectator*/)
  {
    // -----------------------------------------------------------------------------------------
    // Calculate the kinetic energy of the recoil/participant with mass 'mRecoil'
    // in photoproduction off a target (Deuteron) with mass 'mTarget' with beam energy
    // 'beamE' in a quasi-free three-body-decay:
    //
    // electron + Target(Recoil + Spectator) -> Recoil + Spectator + Product
    // (where Product could be e.g. a produced meson, or photon.)
    //
    // Use the reconstructed theta and phi angles in the 4-vector of the
    // participant 'recoil_p4' and the fully reconstructed 4-vector
    // of the mesor/photon 'product_p4' to calculate the kinetic energy of the recoil
    // participant. The spectator mass is 'mSpectator'.
    //
    // Returns Kinetic energy of recoil.
    //
    // full derivation can be found at:
    //   https://jazz.physik.unibas.ch/site/publications/project/project_dieterle.pdf (pg. 29)
    //
    // -----------------------------------------------------------------------------------------

    Double_t beamE = _beam.E();
    Double_t mTarget = 1.876; //Deuteron mass
    Double_t mRecoil = 0.939565420; //neutron mass
    Double_t mSpectator = 0.938272088; //proton mass

    auto product_p4 = _gamma1.P4() + _gamma2.P4();
    auto recoil_p4 = _neutron.P4();

    // set input kinematics variables
    Double_t prodPx = product_p4.Px();
    Double_t prodPy = product_p4.Py();
    Double_t prodPz = product_p4.Pz();
    Double_t prodE  = product_p4.E();

    Double_t partTheta = recoil_p4.Theta();
    Double_t partPhi   = recoil_p4.Phi();

    // calculate terms
    Double_t a =   prodPx * TMath::Sin(partTheta) * TMath::Cos(partPhi)
      + prodPy * TMath::Sin(partTheta) * TMath::Sin(partPhi)
      + (prodPz - beamE ) * TMath::Cos(partTheta);
    Double_t b = prodE - beamE - mTarget;
    Double_t c = (prodE + mRecoil - beamE - mTarget) * (prodE + mRecoil - beamE - mTarget)
      - mSpectator * mSpectator
      - prodPx * prodPx
      - prodPy * prodPy
      - prodPz * prodPz
      - beamE * beamE
      + 2 * beamE * prodPz;

    Double_t quadEquA = b*b - a*a;
    Double_t quadEquB = b*c - 2*a*a*mRecoil;
    Double_t quadEquC = c*c / 4;

  return (-quadEquB + TMath::Sqrt(quadEquB*quadEquB - 4*quadEquA*quadEquC)) / (2*quadEquA);
  } // close Reconstruct_Recoil_KinE

}
