{
  auto FS = pauln::Pi0::Make("ALL","ALL");
  //FS->AddTopology("Electron:Neutron:Gamma1:Gamma2");
  FS->AddTopology("Electron:Gamma1:Gamma2");
  // FS->AddTopology(OTHER_TOPOLOGY);
  
  // try to solve extra neutrals
  // checks how close together clusters are and discount things too close (5cm)
  //FS->MaskParticles(new MaskCalorSplitOffs(50,50,50,1) );//currently only works with inclusive =="ALL" 

  ////Save TreeDataPi0
  FS->UseOutputRootTree();  // removed, as being written out into PDM tree.
  //FS->UseOutputHipoNtuple();

  
  // Remove "garbage" events with p4 = (0,0,0,0) [artefact of recon.]
  ParticleCutsManager pcm_zk{"ZeroKins",1};
  pcm_zk.AddParticleCut("e-",new ZeroKins());
  pcm_zk.AddParticleCut("neutron",new ZeroKins());
  pcm_zk.AddParticleCut("gamma",new ZeroKins());
  FS->RegisterPostTopoAction(pcm_zk);  //before pdm so these events are not included there either.

  // // Perform truth-matching for simulated files
  GenePiTruthAction ev_truth("EventTruth");
  FS->RegisterPostKinAction(ev_truth); //PostKin


  ///Make particle trees first in case want to add cut flags
  // (in this case after removing zero-value parts.)
  // ParticleDataManager pdm{"particle",0};
  // pdm.SetParticleOut(new ParticleKinematics);
  // FS->RegisterPostKinAction(pdm);

  ////

  // Flag events with photons which hit ECAL but miss PCAL (suspicious).
  ParticleCutsManager pcm_PCAL{"hitPCAL",0};
  pcm_PCAL.AddParticleCut("gamma", new hitPCAL());
  FS->RegisterPostTopoAction(pcm_PCAL);

  ////Write to file for later processing
  FS->WriteToFile("finalstates/eff.root");
  //FS->WriteToFile("finalstates/PID_fullcomb.root");
  // FS->WriteToFile("finalstates/PID_photcomb.root");
  // FS->WriteToFile("finalstates/PID_exclcomb.root");


  FS->Print();


  //Delete the final state rather than let ROOT try
  FS.reset();
}
