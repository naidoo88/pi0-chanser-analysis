{
  ////Set hipo file to be analysed

  clas12databases::SetRCDBRootConnection("/home/pauln/code/pi0analysis/rcdb_6336.root");

  HipoData hdata;
  //hdata.SetFile("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_14117.hipo"); //pi0 MC sample file
  hdata.SetFile("/home/pauln/work/data/pass1v0/ndvcs_006336.hipo");

  ////create FinalStateManager
  ////we can load as many saved final states as we like
  FinalStateManager fsm;
  fsm.SetBaseOutDir("output/6336_chanser.root");
  ////Connect the data to the manager
  fsm.LoadData(&hdata);

  ////load one or more FinalStates
  fsm.LoadFinalState("Pi0", "finalstates/PID_fullcomb.root");
  fsm.LoadFinalState("Pi0", "finalstates/PID_photcomb.root"); 
  // fsm.LoadFinalState("Pi0", "finalstates/PID_fullcomb_masked.root"); 
  // fsm.LoadFinalState("Pi0", "finalstates/PID_photcomb_masked.root"); 

  //Max number of particles of any 1 species
  //Whole event disgarded if this not met.
  fsm.GetEventParticles().SetMaxParticles(6);

  ////Run through all events
  fsm.ProcessAll();
  ////Run through N events
  //fsm.ProcessAll(N);


}
