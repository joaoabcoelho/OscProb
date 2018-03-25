void LoadOscProb(void) { 

  // This macro will try to find the OscProb library
  // and load it in your session. You may consider
  // running it in your rootlogon.C file to be loaded
  // by default when you start ROOT.
  
  // ... may be needed for ROOT 5, if "undefined reference to gsl_*" ...
  gSystem->Load("libMathMore"); // automatically loads "gsl" and "cblas"

  // The library name
  TString s("libOscProb.so");

  // Search for the library here or in your library paths
  gSystem->FindFile(TString::Format("./:../:%s", gSystem->GetDynamicPath()), s);
  
  // Try to load the library and complain if failed
  if ((s.Length() != 0) && (gSystem->Load(s.Data()) == 0));
  else std::cout << "... libOscProb could not be loaded ..." << std::endl;

}
