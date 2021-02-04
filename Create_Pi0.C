{
  auto FS = pauln::Pi0::Make("ALL","gamma");
  FS->AddTopology("Electron:Neutron:Gamma1:Gamma2");
  // FS->AddTopology(OTHER_TOPOLOGY);
  
  // try to solve extra neutrals
  // checks how close together clusters are and discount things too close (5cm)
  //FS->MaskParticles(new MaskCalorSplitOffs(50,50,50,1) );//currently only works with inclusive =="ALL" 

  ////Save TreeDataPi0
  FS->UseOutputRootTree();
  //FS->UseOutputHipoNtuple();

  ///Make particle trees first in case want to add cut flags
  // ParticleDataManager pdm{"particle",1};
  // pdm.SetParticleOut(new CLAS12ParticleOutEvent0);
  // FS->RegisterPostTopoAction(pdm);

  ////
  ParticleCutsManager pcm_zk{"ZeroKins",1};
  pcm_zk.AddParticleCut("e-",new ZeroKins());
  pcm_zk.AddParticleCut("neutron",new ZeroKins());
  pcm_zk.AddParticleCut("gamma",new ZeroKins());
  FS->RegisterPostTopoAction(pcm_zk);
 
  ////Write to file for later processing
  FS->WriteToFile("finalstates/PID_photcomb.root");

  FS->Print();


  //Delete the final state rather than let ROOT try
  FS.reset();
}
