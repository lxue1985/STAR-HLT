void hltMiniChain() {

  // Load the shared libs
  gROOT->Macro("Load.C");
  
  Info(__FILE__, "Enabling Logger...");
  gROOT->LoadMacro("LoadLogger.C");
  LoadLogger();


  gSystem->Load("RTS");
  gSystem->Load("StDaqLib"); 
  gSystem->Load("StDAQMaker"); 

  gSystem->Load("StIOMaker");
  gSystem->Load("StHltMaker");

  // Create Chain
  //  StChain *chain = new StChain();
  chain = new StBFChain(); 
  cout << "Create chain " << chain->GetName() << endl;

  chain->SetFlags("in");
  new StIOMaker("IO","r","/star/data05/scratch/yipkin/st_physics_10183035_raw_4030001.daq");
  new StHltMaker(); // The Hlt maker is to be appended to the chain
  // Execute the chain
  chain->Init();    // start the chain

  //  chain->Make();    // produce one event
  int counter = 10;
  while ( counter &&( chain->MakeEvent() == kStOK ) )   { counter--;}

  chain->Finish();  // close chain

  delete chain;
}

