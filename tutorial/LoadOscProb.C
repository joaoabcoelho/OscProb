#include <iostream>
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"

void LoadOscProb(TString logLevel = "info") {

  // This macro will try to find the OscProb library
  // and load it in your session. You may consider
  // running it in your rootlogon.C file to be loaded
  // by default when you start ROOT.

  TString currentDir = gSystem->DirName(gInterpreter->GetCurrentMacroName());

  bool verbose = logLevel.Contains("v");
  bool quiet   = logLevel.Contains("q") && !verbose;

  // The library name
  TString s("libOscProb");

  // Search for the library here or in your library paths
  gSystem->AddDynamicPath(currentDir+"/lib");
  gSystem->AddDynamicPath(currentDir+"/../lib");
  gSystem->AddDynamicPath(currentDir+"/lib64");
  gSystem->AddDynamicPath(currentDir+"/../lib64");
  // gSystem->FindFile(TString::Format("./lib:../lib:%s", gSystem->GetDynamicPath()), s);
  s = gSystem->FindDynamicLibrary(s);

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

      case 1: { // Library already loaded
        std::cout << "libOscProb already loaded at: ";
        TObjArray* libs = TString(gSystem->GetLibraries()).Tokenize(" ");

        for(int i =0; i<libs->GetEntries(); i++){

          TString lib = libs->At(i)->GetName();

          if(lib.Contains("libOscProb.so")){
            std::cout << lib << std::endl;
            s = lib;
          }

        }
        break;
      }
      default: std::cout << "libOscProb could not be loaded. Error Code: " << errCode << std::endl;

    }

  }
  else if(verbose){

    std::cout << "Loaded library " << s << std::endl;

  }

  TString dirname = gSystem->DirName(s);
  dirname += "/../inc";

  if(verbose) std::cout << "Adding " << dirname << " to include path" << std::endl;
  gInterpreter->AddIncludePath(dirname);

}
