#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace OscProb;

#pragma link C++ struct OscProb::NuPath+;
#pragma link C++ struct OscProb::PremLayer+;
#pragma link C++ struct OscProb::EigenPoint+;
#pragma link C++ class vector<OscProb::NuPath>;
#pragma link C++ class vector<OscProb::PremLayer>;
#pragma link C++ class vector<complexD>;
#pragma link C++ class set<OscProb::EigenPoint>;


#pragma link C++ class OscProb::PMNS_Base+;
#pragma link C++ class OscProb::PMNS_Fast+;
#pragma link C++ class OscProb::PMNS_NSI+;
#pragma link C++ class OscProb::PMNS_Sterile+;
#pragma link C++ class OscProb::PMNS_Deco+;
#pragma link C++ class OscProb::PMNS_Decay+;

#pragma link C++ class OscProb::PremModel+;

#endif
