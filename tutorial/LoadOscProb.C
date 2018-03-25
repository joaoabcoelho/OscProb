
void LoadOscProb(TString logLevel = "info") { 

  // This macro will try to find the OscProb library
  // and load it in your session. You may consider
  // running it in your rootlogon.C file to be loaded
  // by default when you start ROOT.
  
  bool verbose = logLevel.Contains("v");
  bool quiet   = logLevel.Contains("q") && !verbose;
  
  // ... may be needed for ROOT 5, if "undefined reference to gsl_*" ...
  gSystem->Load("libMathMore"); // automatically loads "gsl" and "cblas"

  // The library name
  TString s("libOscProb.so");

  // Search for the library here or in your library paths
  gSystem->FindFile(TString::Format("./:../:%s", gSystem->GetDynamicPath()), s);

  // Check if library file was found
  if(s.Length() == 0){
    if(!quiet) std::cout << "libOscProb not found!" << std::endl;
    return;
  } 
  
  // Try to load the library and complain if failed
  int errCode = gSystem->Load(s.Data());
  
  // If loading failed, print message
  if(errCode != 0 && !quiet){

    switch(errCode){

      case 1: // Library already loaded
        std::cout << "libOscProb already loaded at: ";
        TObjArray* libs = TString(gSystem->GetLibraries()).Tokenize(" ");
  
        for(int i =0; i<libs->GetEntries(); i++){
  
          TString lib = libs->At(i)->GetName();
    
          if(lib.Contains("libOscProb.so")) std::cout << lib << std::endl;
  
        }
        break;
        
      default: std::cout << "libOscProb could not be loaded. Error Code: " << errCode << std::endl;
      
    }

  }
  else if(verbose){
  
    std::cout << "Loaded library " << s << std::endl;

  }
    
}
