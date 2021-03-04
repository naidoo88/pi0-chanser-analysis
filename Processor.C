
{
  // clas12databases::SetRCDBRootConnection("/w/work3/home/pauln/analysis/processed/rcdb/rcdb_allP1.root");
  clas12databases::SetRCDBRootConnection("/w/work3/home/pauln/analysis/processed/rcdb/rcdb_skimmedP1.root");
  // clas12databases::SetRCDBRootConnection("/w/work3/home/pauln/analysis/processed/rcdb/rcdb_tau_skimmedP1.root");

  //clas12databases::SetRCDBRootConnection("/home/pauln/code/pi0analysis/rcdb_6336.root");

  
  clas12root::HipoChain chain;
  chain.SetReaderTags({0});
  //chain.Add("/w/work3/home/pauln/data/pass1v0/ndvcs_006*.hipo");
  // chain.Add("/w/work3/home/pauln/data/skimmed_pass1v0/ndvcs_006*preskim.hipo");

  //Subset of Data:
  // chain.Add("/w/work3/home/pauln/data/skimmed_pass1v0/ndvcs_0061*preskim.hipo"); //subset of ~80GB files
  // chain.Add("/w/work3/home/pauln/data/skimmed_pass1v0/ndvcs_0062*preskim.hipo"); //--------- ^^

  //chain.Add("/w/work3/home/pauln/data/skimmed_pass1v0/ndvcs_0064*preskim.hipo"); //subset of 58 files chanser_MCsub
  // chain.Add("/scratch/pauln/data/ndvcs_006*.hipo");
  //chain.Add("/home/pauln/work/data/skimmed_pass1v0/ndvcs_006336_preskim.hipo");

  // Batches for MC data:
  //chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_14117.hipo");
  //chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_1*.hipo"); //chanser_MCbig
  chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_11*.hipo"); //*
  chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_12*.hipo"); //*
  chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_13*.hipo"); //*
  chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_14*.hipo"); //*
  chain.Add("/w/work3/home/pauln/sim/pi0/cooked_out_rgb_pi0_15*.hipo"); //*
  
  // chain.GetNRecords();
  chanser::HipoProcessor processor(&chain,"finalstates/eff.txt", "output/chanser_MCsub_eff");
  //chanser::HipoProcessor processor(&chain,"finalstates/FullPhotonExclCombs.txt", "output/chanser_MCbig");
  //chanser::HipoProcessor processor(&chain,"finalstates/FullAndPhotonCombs.txt", "output/chanser_MCsub");
  // chanser::HipoProcessor processor(&chain,"finalstates/FullPhotonExclCombs.txt", "output/chanser_ANAsub");

  //chanser::HipoProcessor processor(&chain, "finalstates/FullAndPhotonCombs.txt", "output/chanser_PIDtest");

  processor.AddOption("HIPOPROCESSOR_MAXPARTICLES","6");
  //processor.AddOption("HIPOPROCESSOR_ANADB","$CHANSER/rga_actions/anadb/RGA_ACTIONS_PASS1.db");
  processor.AddOption("HIPOPROCESSOR_ANADB","$CHANSER/anadbs/RunPeriodPass1.db"); //for MC target shift
  processor.AddOption("HIPOPROCESSOR_RUNPERIOD","spring_2019");
  gBenchmark->Start("hd");

  gProof->Process(&processor,chain.GetNRecords());

  gBenchmark->Stop("hd");
  gBenchmark->Print("hd");

}