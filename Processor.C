
{
  
 clas12root::HipoChain chain;
 chain.SetReaderTags({0});
 //chain.Add("/work/jlab/clas12data/pass1_test/skim8_005051.hipo");
 chain.Add("/home/pauln/work/data/pass1v0/ndvcs_006*preskim.hipo");
 // chain.GetNRecords();
 chanser::HipoProcessor processor(&chain,"finalstates/FullAndPhotonCombs.txt");
 processor.AddOption("HIPOPROCESSOR_MAXPARTICLES","6");

 gBenchmark->Start("hd");

 gProof->Process(&processor,chain.GetNRecords());
 
 gBenchmark->Stop("hd");
 gBenchmark->Print("hd");

}