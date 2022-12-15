#include"libs.hh"
#include"definitions.hh"
#include"invariant_mass.hh"
#include <TMath.h>
#include "mass.h"
//#include"R3BMDFWrapper.h"
//#include"R3BMDFWrapper.cxx"
//class R3BMDFWrapper;

using namespace std;

Double_t VoV(Double_t* V, Double_t* W){  //scalar product
  Double_t res = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];
  return res;
}

Double_t Estar1n(Double_t mn, Double_t mf, Double_t* Pnn1, Double_t* Pnn3){ //Calculate Estar with only one neutron
  Double_t res;
  res = TMath::Sqrt(mf*mf + mn*mn + 2*Pnn1[3]*Pnn3[3] - 2*VoV(Pnn1,Pnn3))- mf - mn;
  return res;
}

Double_t Estar2n(Double_t mn, Double_t mf, Double_t* Pnn1, Double_t* Pnn2, Double_t* Pnn3){
  Double_t res;
  res = TMath::Sqrt(mf*mf + 2*mn*mn + 2*Pnn1[3]*Pnn3[3] - 2*VoV(Pnn1,Pnn3) + 2*Pnn2[3]*Pnn3[3] - 2*VoV(Pnn2,Pnn3) + 2*Pnn2[3]*Pnn1[3] - 2*VoV(Pnn2,Pnn1))
					   - mf - 2*mn;
  return res;
}

void invariant_mass(Char_t* InputFile, Bool_t MakeCalMap, Char_t* WhichNeutronWall){

  cout << "\n\n\t*************************************************" << endl;
  cout << "\t*                                               *" << endl;
  cout << "\t*        Invariant-Mass macros in R3BRoot       *" << endl;
  cout << "\t*                                               *" << endl;
  cout << "\t*************************************************" << endl;

  Double_t Res_Pos = 0.;
  Double_t Res_Time = 0.;
  //Double_t NeutronWall_ZOffset = 1100.; //11m
  //Double_t NeutronWall_ZOffset = 800.; //8m
  Double_t NeutronWall_ZOffset = 1500.; //15m
  TString NWall = "";

  if((TString)WhichNeutronWall=="Neuland"){
    cout << "Using Neuland Resolution !! " << endl; 
    Res_Pos = 3.; //cm
    Res_Time = 0.250; //ns
    NWall = "Neuland";
  }else if((TString)WhichNeutronWall=="Mona"){
    cout << "Using Mona Resolution !! " << endl; 
    Res_Pos = 7.; //cm
    Res_Time = 0.250; //ns
    NWall = "Mona";
  }else if((TString)WhichNeutronWall=="Land"){
    cout << "Using Land Resolution !! " << endl; 
    Res_Pos = 10.; //cm
    Res_Time = 0.350; //ns
    NWall = "Land";
  }
  
  gRandom = new TRandom3();
  gRandom->SetSeed(1);

  Double_t mn = 939.565;
  
  TRandom3 *RandX = new TRandom3(0);
  
  TString inputName = (TString)"input/"+InputFile;
  TFile *f = new TFile(inputName);
  TTree *tree = (TTree*)f->Get("evt");
  if(tree==NULL){ cout << "\nERROR! No specified tree in the file " << InputFile  << endl; return; }
  TString outputName = (TString)"output/"+InputFile;
  TFile* output = new TFile(outputName, "RECREATE");
  TTree *OutTree = new TTree("OutTree","OutTree");
  
  //--------- GetData from sim tree --------
  TClonesArray * MCTrackArray   = new TClonesArray("R3BMCTrack");
  TClonesArray * NeutronWallPointArray;
  if(NWall!="Land"){
    NeutronWallPointArray = new TClonesArray("R3BNeulandPoint");
    tree->SetBranchAddress("NeulandPoints",	&NeutronWallPointArray);
  }else{
    NeutronWallPointArray = new TClonesArray("R3BLandPoint");
    tree->SetBranchAddress("LandPoint",	&NeutronWallPointArray);
  }
  tree->SetBranchAddress("MCTrack",	&MCTrackArray);


  cout << "\n-- Input file: " << InputFile;
  cout << "\n-- Number of events in the file: " << tree->GetEntries() << endl;

  const Int_t Nentries=tree->GetEntries();

  R3BMCTrack *MCTrack;
  R3BNeulandPoint *NeutronWallPoint;
  R3BLandPoint *NeutronWallPointLand;
  TVector3 NeutronWall_Pos;
  Double_t NeutronWall_Time;
  Int_t NeutronWall_ID;
  vector <Double_t> *NEUTRONWALL_X = new vector <Double_t>;
  vector <Double_t> *NEUTRONWALL_Y = new vector <Double_t>;
  vector <Double_t> *NEUTRONWALL_Z = new vector <Double_t>;
  vector <Double_t> *NEUTRONWALL_Time = new vector <Double_t>;
  vector <Int_t> *NEUTRONWALL_ID = new vector <Int_t>;
  vector <Double_t> *NEUTRON_X = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Y = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Z = new vector <Double_t>;
  vector <Double_t> *NEUTRON_P = new vector <Double_t>;
  vector <Double_t> *NEUTRON_PX = new vector <Double_t>;
  vector <Double_t> *NEUTRON_PY = new vector <Double_t>;
  vector <Double_t> *NEUTRON_PZ = new vector <Double_t>;
  vector <Double_t> *NEUTRON_E = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Time = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Beta = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Gamma = new vector <Double_t>;
  vector <Double_t> *NEUTRON_Dist = new vector <Double_t>;
  Int_t NEUTRON_Mult;
  
  //MCTrack Info
  Double_t MCTrack_X, MCTrack_Y, MCTrack_Z;
  Int_t MCTrack_Mult;
  Double_t MCTrack1_Theta, MCTrack2_Theta, MCTrack3_Theta;
  Double_t MCTrack1_Phi, MCTrack2_Phi, MCTrack3_Phi;
  Double_t MCTrack1_E, MCTrack2_E, MCTrack3_E;
  Double_t MCTrack1_Ek, MCTrack2_Ek, MCTrack3_Ek;
  Double_t MCTrack1_Px, MCTrack2_Px, MCTrack3_Px;
  Double_t MCTrack1_Py, MCTrack2_Py, MCTrack3_Py;
  Double_t MCTrack1_Pz, MCTrack2_Pz, MCTrack3_Pz;
  Double_t MCTrack1_Mass, MCTrack2_Mass, MCTrack3_Mass;
  
  OutTree->Branch("MCTrack_X",&MCTrack_X);
  OutTree->Branch("MCTrack_Y",&MCTrack_Y);
  OutTree->Branch("MCTrack_Z",&MCTrack_Z);
  OutTree->Branch("MCTrack1_Theta",&MCTrack1_Theta);
  OutTree->Branch("MCTrack2_Theta",&MCTrack2_Theta);
  OutTree->Branch("MCTrack1_Phi",&MCTrack1_Phi);
  OutTree->Branch("MCTrack2_Phi",&MCTrack2_Phi);
  OutTree->Branch("MCTrack1_E",&MCTrack1_E);
  OutTree->Branch("MCTrack2_E",&MCTrack2_E);
  OutTree->Branch("MCTrack1_Ek",&MCTrack1_Ek);
  OutTree->Branch("MCTrack2_Ek",&MCTrack2_Ek);
  OutTree->Branch("MCTrack1_Px",&MCTrack1_Px);
  OutTree->Branch("MCTrack2_Px",&MCTrack2_Px);
  OutTree->Branch("MCTrack1_Py",&MCTrack1_Py);
  OutTree->Branch("MCTrack2_Py",&MCTrack2_Py);
  OutTree->Branch("MCTrack1_Pz",&MCTrack1_Pz);
  OutTree->Branch("MCTrack2_Pz",&MCTrack2_Pz);
  OutTree->Branch("MCTrack1_Mass",&MCTrack1_Mass);
  OutTree->Branch("MCTrack2_Mass",&MCTrack2_Mass);
  OutTree->Branch("MCTrack3_Theta",&MCTrack3_Theta);
  OutTree->Branch("MCTrack3_Phi",&MCTrack3_Phi);
  OutTree->Branch("MCTrack3_E",&MCTrack3_E);
  OutTree->Branch("MCTrack3_Ek",&MCTrack3_Ek);
  OutTree->Branch("MCTrack3_Px",&MCTrack3_Px);
  OutTree->Branch("MCTrack3_Py",&MCTrack3_Py);
  OutTree->Branch("MCTrack3_Pz",&MCTrack3_Pz);
  OutTree->Branch("MCTrack3_Mass",&MCTrack3_Mass);
  
  //Invariant-mass
  Double_t Erel1n_in;
  Double_t Erel2n_in;
  Double_t Erel1n;
  Double_t Erel2n;

  OutTree->Branch("Erel1n_in",&Erel1n_in);
  OutTree->Branch("Erel1n",&Erel1n);
  OutTree->Branch("Erel2n_in",&Erel2n_in);
  OutTree->Branch("Erel2n",&Erel2n);
  
  //NEUTRONWALL
  OutTree->Branch("NEUTRONWALL_X",&NEUTRONWALL_X);
  OutTree->Branch("NEUTRONWALL_Y",&NEUTRONWALL_Y);
  OutTree->Branch("NEUTRONWALL_Z",&NEUTRONWALL_Z);
  OutTree->Branch("NEUTRONWALL_Time",&NEUTRONWALL_Time);
  OutTree->Branch("NEUTRONWALL_ID",&NEUTRONWALL_ID);
  OutTree->Branch("NEUTRON_X",&NEUTRON_X);
  OutTree->Branch("NEUTRON_Y",&NEUTRON_Y);
  OutTree->Branch("NEUTRON_Z",&NEUTRON_Z);
  OutTree->Branch("NEUTRON_P",&NEUTRON_P);
  OutTree->Branch("NEUTRON_PX",&NEUTRON_PX);
  OutTree->Branch("NEUTRON_PY",&NEUTRON_PY);
  OutTree->Branch("NEUTRON_PZ",&NEUTRON_PZ);
  OutTree->Branch("NEUTRON_E",&NEUTRON_E);
  OutTree->Branch("NEUTRON_Time",&NEUTRON_Time);
  OutTree->Branch("NEUTRON_Beta",&NEUTRON_Beta);
  OutTree->Branch("NEUTRON_Gamma",&NEUTRON_Gamma);
  OutTree->Branch("NEUTRON_Dist",&NEUTRON_Dist);
  OutTree->Branch("NEUTRON_Mult",&NEUTRON_Mult);
  
  cout << "\n-- Collecting track data.....\n";
  Int_t counter = 0;
  
  for(Int_t ev=0; ev<Nentries; ev++){

    if (ev%100==0) cout << "\r-- Working on entry " << ev << flush;

    MCTrackArray->Clear();
    NeutronWallPointArray->Clear();

    MCTrack1_Theta = -9999.;
    MCTrack2_Theta = -9999.;
    MCTrack1_Phi = -9999.;
    MCTrack2_Phi = -9999.;
    MCTrack1_E = -9999;
    MCTrack2_E = -9999;
    MCTrack1_Ek = -9999;
    MCTrack2_Ek = -9999;
    MCTrack1_Px = -9999;
    MCTrack2_Px = -9999;
    MCTrack1_Py = -9999;
    MCTrack2_Py = -9999;
    MCTrack1_Pz = -9999;
    MCTrack2_Pz = -9999;
    MCTrack1_Mass = -9999;
    MCTrack2_Mass = -9999;
    MCTrack3_Theta = -9999.;
    MCTrack3_Phi = -9999.;
    MCTrack3_E = -9999;
    MCTrack3_Ek = -9999;
    MCTrack3_Px = -9999;
    MCTrack3_Py = -9999;
    MCTrack3_Pz = -9999;
    MCTrack3_Mass = -9999;
    
    Erel1n_in = -9999.;
    Erel2n_in = -9999.;
    Erel1n = -9999.;
    Erel2n = -9999.;

    NEUTRONWALL_X->clear();
    NEUTRONWALL_Y->clear();
    NEUTRONWALL_Z->clear();
    NEUTRONWALL_Time->clear();
    NEUTRONWALL_ID->clear();
    NEUTRON_X->clear();
    NEUTRON_Y->clear();
    NEUTRON_Z->clear();
    NEUTRON_Time->clear();
    NEUTRON_Beta->clear();
    NEUTRON_Gamma->clear();
    NEUTRON_Dist->clear();
    NEUTRON_P->clear();
    NEUTRON_PX->clear();
    NEUTRON_PY->clear();
    NEUTRON_PZ->clear();
    NEUTRON_E->clear();
    NEUTRON_Mult=0;

    
    tree->GetEntry(ev);

    //MCTrack
    MCTrack_Mult = MCTrackArray->GetEntriesFast();
    Int_t TrackCounter = 0;
    for(Int_t i=0; i<MCTrackArray->GetEntriesFast(); i++){

      MCTrack = (R3BMCTrack*) MCTrackArray->At(i); 
      
      if(MCTrack->GetMotherId()==-1){

	MCTrack_X = MCTrack->GetStartX();
	MCTrack_Y = MCTrack->GetStartY();
	MCTrack_Z = MCTrack->GetStartZ();

	TVector3 Track_Mom;
	Track_Mom.SetXYZ(MCTrack->GetPx(),MCTrack->GetPy(),MCTrack->GetPz());
	
	if(TrackCounter==0){
	  MCTrack1_E = MCTrack->GetEnergy()*1000.;
	  MCTrack1_Px = MCTrack->GetPx()*1000.;
	  MCTrack1_Py = MCTrack->GetPy()*1000.;
	  MCTrack1_Pz = MCTrack->GetPz()*1000.;
	  MCTrack1_Mass = MCTrack->GetMass()*1000.;
	  Double_t MCTrack_Gamma = MCTrack->GetEnergy()/MCTrack->GetMass();
	  MCTrack1_Ek = (MCTrack_Gamma-1.)*MCTrack->GetMass()*1000.;
	  MCTrack1_Theta = Track_Mom.Theta();//*TMath::RadToDeg();
	  MCTrack1_Phi = Track_Mom.Phi();//*TMath::RadToDeg();
	  TrackCounter++;
	}else if(TrackCounter==1){
	  MCTrack2_E = MCTrack->GetEnergy()*1000.;
	  MCTrack2_Px = MCTrack->GetPx()*1000.;
	  MCTrack2_Py = MCTrack->GetPy()*1000.;
	  MCTrack2_Pz = MCTrack->GetPz()*1000.;
	  MCTrack2_Mass = MCTrack->GetMass()*1000.;
	  Double_t MCTrack_Gamma = MCTrack->GetEnergy()/MCTrack->GetMass();
	  MCTrack2_Ek = (MCTrack_Gamma-1.)*MCTrack->GetMass()*1000.;
	  MCTrack2_Theta = Track_Mom.Theta();//*TMath::RadToDeg();
	  MCTrack2_Phi = Track_Mom.Phi();//*TMath::RadToDeg();
	  TrackCounter++;
	}else{
	  MCTrack3_E = MCTrack->GetEnergy()*1000.;
	  MCTrack3_Px = MCTrack->GetPx()*1000.;
	  MCTrack3_Py = MCTrack->GetPy()*1000.;
	  MCTrack3_Pz = MCTrack->GetPz()*1000.;
	  MCTrack3_Mass = MCTrack->GetMass()*1000.;
	  Double_t MCTrack_Gamma = MCTrack->GetEnergy()/MCTrack->GetMass();
	  MCTrack3_Ek = (MCTrack_Gamma-1.)*MCTrack->GetMass()*1000.;
	  MCTrack3_Theta = Track_Mom.Theta();//*TMath::RadToDeg();
	  MCTrack3_Phi = Track_Mom.Phi();//*TMath::RadToDeg();
	  TrackCounter++;
	}
      }
    }

    Double_t Pn1[4];
    Pn1[0] = MCTrack2_Px;
    Pn1[1] = MCTrack2_Py;
    Pn1[2] = MCTrack2_Pz;
    Pn1[3] = MCTrack2_E;
    Double_t Pn2[4];
    Pn2[0] = MCTrack3_Px;
    Pn2[1] = MCTrack3_Py;
    Pn2[2] = MCTrack3_Pz;
    Pn2[3] = MCTrack3_E;
    Double_t Pf[4];
    Pf[0] = MCTrack1_Px;
    Pf[1] = MCTrack1_Py;
    Pf[2] = MCTrack1_Pz;
    Pf[3] = MCTrack1_E;
    
    Erel1n_in = Estar1n(MCTrack2_Mass,MCTrack1_Mass,Pn1,Pf);
    Erel2n_in = Estar2n(MCTrack2_Mass,MCTrack1_Mass,Pn1,Pn2,Pf);

    //NeutronWall Analysis
    Int_t NeutronWall_Mult = NeutronWallPointArray->GetEntriesFast();

    for(Int_t i=0; i<NeutronWallPointArray->GetEntriesFast(); i++){
      if(NWall!="Land"){
	NeutronWallPoint = (R3BNeulandPoint*) NeutronWallPointArray->At(i);
	NeutronWallPoint->Position(NeutronWall_Pos);
	NeutronWall_ID = NeutronWallPoint->GetPaddle();
	NeutronWall_Time = NeutronWallPoint->GetTime();
      }else{
	NeutronWallPointLand = (R3BLandPoint*) NeutronWallPointArray->At(i);
	NeutronWallPointLand->Position(NeutronWall_Pos);
	NeutronWall_ID = NeutronWallPointLand->GetSector();
	NeutronWall_Time = NeutronWallPointLand->GetTime();
      }

      Bool_t Is_New_ID = true;
      
      for(Int_t j=0 ; j<NEUTRONWALL_X->size() ; j++){
	if(NEUTRONWALL_ID->at(j)==NeutronWall_ID){ //Merge IDs
	  if(NEUTRONWALL_Time->at(j)>NeutronWall_Time){
	    NEUTRONWALL_X->at(j) = NeutronWall_Pos.X();
	    NEUTRONWALL_Y->at(j) = NeutronWall_Pos.Y();
	    NEUTRONWALL_Z->at(j) = NeutronWall_Pos.Z();
	    NEUTRONWALL_Time->at(j) = NeutronWall_Time;
	    NEUTRONWALL_ID->at(j) = NeutronWall_ID;
	  }
	  Is_New_ID = false;
	}
      }
      if(Is_New_ID){
	NEUTRONWALL_X->push_back(NeutronWall_Pos.X());
	NEUTRONWALL_Y->push_back(NeutronWall_Pos.Y());
	NEUTRONWALL_Z->push_back(NeutronWall_Pos.Z());
	NEUTRONWALL_Time->push_back(NeutronWall_Time);
	NEUTRONWALL_ID->push_back(NeutronWall_ID);
      }
    }

    //Time ordering
    if(NEUTRONWALL_X->size()>0){
      for(Int_t i=0; i<NEUTRONWALL_X->size()-1; i++){
	for(Int_t j=i+1; j<NEUTRONWALL_X->size(); j++){
	  if(NEUTRONWALL_Time->at(j)<NEUTRONWALL_Time->at(i)){
	    Double_t tempD;
	    Int_t tempI;

	    tempD = NEUTRONWALL_X->at(i);
	    NEUTRONWALL_X->at(i) = NEUTRONWALL_X->at(j);
	    NEUTRONWALL_X->at(j) = tempD;

	    tempD = NEUTRONWALL_Y->at(i);
	    NEUTRONWALL_Y->at(i) = NEUTRONWALL_Y->at(j);
	    NEUTRONWALL_Y->at(j) = tempD;

	    tempD = NEUTRONWALL_Z->at(i);
	    NEUTRONWALL_Z->at(i) = NEUTRONWALL_Z->at(j);
	    NEUTRONWALL_Z->at(j) = tempD;

	    tempD = NEUTRONWALL_Time->at(i);
	    NEUTRONWALL_Time->at(i) = NEUTRONWALL_Time->at(j);
	    NEUTRONWALL_Time->at(j) = tempD;

	    tempI = NEUTRONWALL_ID->at(i);
	    NEUTRONWALL_ID->at(i) = NEUTRONWALL_ID->at(j);
	    NEUTRONWALL_ID->at(j) = tempI;
	    
	  }
	}
      }
    }
	
    //NeutronWall Hit to Neutron -> Discretize and add Resolution
    for(Int_t i=0; i<NEUTRONWALL_X->size(); i++){
      if(NWall=="Mona"){ //Mona
	//From Time Difference
	NEUTRON_X->push_back(RandX->Gaus(NEUTRONWALL_X->at(i),Res_Pos));
	//From Bar Thickness
	Int_t WallNbr = 1+(Int_t)((NEUTRONWALL_ID->at(i)-1)/16.);
	Int_t BarInWall = NEUTRONWALL_ID->at(i) - 16*(WallNbr-1);
	Int_t PosY = -75.+10.*(BarInWall-1);
	NEUTRON_Y->push_back(PosY);
	Int_t PosZ = 5.+10.*(WallNbr-1);
	NEUTRON_Z->push_back(PosZ+NeutronWall_ZOffset);
	//Time resolution for ToF
	NEUTRON_Time->push_back(RandX->Gaus(NEUTRONWALL_Time->at(i),Res_Time));
      }
      
      if(NWall=="Neuland"){ //Neuland
	Int_t WallNbr = 1+(Int_t)((NEUTRONWALL_ID->at(i)-1)/50.);
	Int_t BarInWall = NEUTRONWALL_ID->at(i) - 50*(WallNbr-1);
	if(WallNbr%2==0){ //vertical bars
	  Int_t PosX = -122.5+5.*(BarInWall-1);
	  NEUTRON_Y->push_back(RandX->Gaus(NEUTRONWALL_Y->at(i),Res_Pos));
	  NEUTRON_X->push_back(PosX);
	}else{
	  Int_t PosY = -122.5+5.*(BarInWall-1);
	  NEUTRON_X->push_back(RandX->Gaus(NEUTRONWALL_X->at(i),Res_Pos));
	  NEUTRON_Y->push_back(PosY);
	}
	Int_t PosZ = 2.5+5.*(WallNbr-1);
	NEUTRON_Z->push_back(PosZ+NeutronWall_ZOffset);
	//Time resolution for ToF
	NEUTRON_Time->push_back(RandX->Gaus(NEUTRONWALL_Time->at(i),Res_Time));
      }

      if(NWall=="Land"){ //Land
	Int_t WallNbr = 1+(Int_t)((NEUTRONWALL_ID->at(i)-1)/20.);
	Int_t BarInWall = NEUTRONWALL_ID->at(i) - 20*(WallNbr-1);
	if(WallNbr%2==0){ //vertical bars
	  Int_t PosX = -95.+10.*(BarInWall-1);
	  NEUTRON_Y->push_back(RandX->Gaus(NEUTRONWALL_Y->at(i),Res_Pos));
	  NEUTRON_X->push_back(PosX);
	}else{
	  Int_t PosY = -95.+10.*(BarInWall-1);
	  NEUTRON_X->push_back(RandX->Gaus(NEUTRONWALL_X->at(i),Res_Pos));
	  NEUTRON_Y->push_back(PosY);
	}
	Int_t PosZ = 5.+10.*(WallNbr-1);
	NEUTRON_Z->push_back(PosZ+NeutronWall_ZOffset);
	//Time resolution for ToF
	NEUTRON_Time->push_back(RandX->Gaus(NEUTRONWALL_Time->at(i),Res_Time));
      }
    }

    //Remove cross-talk
    Double_t Dist_Cut = 0.;
    if(NWall=="Neuland"){
      Dist_Cut = 12.5*2.;
    }
    if(NWall=="Mona"){
      Dist_Cut = 20.*2.;
    }
    if(NWall=="Land"){
      Dist_Cut = 20.*2.;
    }

    if(NEUTRON_X->size()>0){
      for(Int_t i=0; i<NEUTRON_X->size()-1; i++){
	for(Int_t j=i+1; j<NEUTRON_X->size(); j++){

	  Double_t Dist_2n = TMath::Sqrt((NEUTRON_X->at(j)-NEUTRON_X->at(i))*(NEUTRON_X->at(j)-NEUTRON_X->at(i))+(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))*(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))+(NEUTRON_Z->at(j)-NEUTRON_Z->at(i))*(NEUTRON_Z->at(j)-NEUTRON_Z->at(i)));

	  if(Dist_2n<Dist_Cut){
	    NEUTRON_X->erase(NEUTRON_X->begin() + j);
	    NEUTRON_Y->erase(NEUTRON_Y->begin() + j);
	    NEUTRON_Z->erase(NEUTRON_Z->begin() + j);
	    NEUTRON_Time->erase(NEUTRON_Time->begin() + j);
	    j--;
	  }
	  
	}
      }
    }

    //Causality Cut
    if(NEUTRON_X->size()>0){
      for(Int_t i=0; i<NEUTRON_X->size()-1; i++){
	for(Int_t j=i+1; j<NEUTRON_X->size(); j++){

	  Double_t Dist_2n = TMath::Sqrt((NEUTRON_X->at(j)-NEUTRON_X->at(i))*(NEUTRON_X->at(j)-NEUTRON_X->at(i))+(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))*(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))+(NEUTRON_Z->at(j)-NEUTRON_Z->at(i))*(NEUTRON_Z->at(j)-NEUTRON_Z->at(i)));

	  if(Dist_2n<Dist_Cut){
	    NEUTRON_X->erase(NEUTRON_X->begin() + j);
	    NEUTRON_Y->erase(NEUTRON_Y->begin() + j);
	    NEUTRON_Z->erase(NEUTRON_Z->begin() + j);
	    NEUTRON_Time->erase(NEUTRON_Time->begin() + j);
	    j--;
	  }
	  
	}
      }
    }
    
    if(NEUTRON_X->size()>1){
      for(Int_t i=0; i<NEUTRON_X->size()-1; i++){
	for(Int_t j=i+1; j<NEUTRON_X->size(); j++){

	  Double_t Dist = TMath::Sqrt(NEUTRON_X->at(i)*NEUTRON_X->at(i)+NEUTRON_Y->at(i)*NEUTRON_Y->at(i)+NEUTRON_Z->at(i)*NEUTRON_Z->at(i));
	  Double_t Beta1 = Dist/NEUTRON_Time->at(i)/29.9792458;
	  Double_t Dist_2n = TMath::Sqrt((NEUTRON_X->at(j)-NEUTRON_X->at(i))*(NEUTRON_X->at(j)-NEUTRON_X->at(i))+(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))*(NEUTRON_Y->at(j)-NEUTRON_Y->at(i))+(NEUTRON_Z->at(j)-NEUTRON_Z->at(i))*(NEUTRON_Z->at(j)-NEUTRON_Z->at(i)));
	  Double_t Beta2 = Dist_2n/(NEUTRON_Time->at(j)-NEUTRON_Time->at(i))/29.9792458;

	  if(Beta1>Beta2){
	    NEUTRON_X->erase(NEUTRON_X->begin() + j);
	    NEUTRON_Y->erase(NEUTRON_Y->begin() + j);
	    NEUTRON_Z->erase(NEUTRON_Z->begin() + j);
	    NEUTRON_Time->erase(NEUTRON_Time->begin() + j);
	    j--;
	  }
	  
	}
      }
    }

    for(Int_t i=0; i<NEUTRON_X->size(); i++){
      NEUTRON_Dist->push_back(TMath::Sqrt(NEUTRON_X->at(i)*NEUTRON_X->at(i)+NEUTRON_Y->at(i)*NEUTRON_Y->at(i)+NEUTRON_Z->at(i)*NEUTRON_Z->at(i)));
      NEUTRON_Beta->push_back(NEUTRON_Dist->at(i)/NEUTRON_Time->at(i)/29.9792458);
      NEUTRON_Gamma->push_back(1./TMath::Sqrt(1-NEUTRON_Beta->at(i)*NEUTRON_Beta->at(i)));
      NEUTRON_P->push_back(NEUTRON_Gamma->at(i)*NEUTRON_Beta->at(i)*mn);
      NEUTRON_PX->push_back(NEUTRON_P->at(0)*NEUTRON_X->at(i)/NEUTRON_Dist->at(i));
      NEUTRON_PY->push_back(NEUTRON_P->at(0)*NEUTRON_Y->at(i)/NEUTRON_Dist->at(i));
      NEUTRON_PZ->push_back(NEUTRON_P->at(0)*NEUTRON_Z->at(i)/NEUTRON_Dist->at(i));
      NEUTRON_E->push_back(TMath::Sqrt(mn*mn + NEUTRON_PX->at(i)*NEUTRON_PX->at(i) + NEUTRON_PY->at(i)*NEUTRON_PY->at(i) + NEUTRON_PZ->at(i)*NEUTRON_PZ->at(i)));
    }

    Double_t MOM_Res = 0.; //33.;

    NEUTRON_Mult = NEUTRON_X->size();
    
    if(NEUTRON_X->size()>0){
      Double_t PNeutron[4];
      PNeutron[0] = NEUTRON_PX->at(0);
      PNeutron[1] = NEUTRON_PY->at(0);
      PNeutron[2] = NEUTRON_PZ->at(0);
      PNeutron[3] = NEUTRON_E->at(0);

      Double_t Pf_Res[4];
      Double_t Pf_Tot = TMath::Sqrt(Pf[0]*Pf[0] + Pf[1]*Pf[1] + Pf[2]*Pf[2]);
      Double_t Pf_Tot_Res = RandX->Gaus(Pf_Tot,MOM_Res);
      Pf_Res[0] = Pf_Tot_Res*Pf[0]/Pf_Tot;
      Pf_Res[1] = Pf_Tot_Res*Pf[1]/Pf_Tot;
      Pf_Res[2] = Pf_Tot_Res*Pf[2]/Pf_Tot;
      Pf_Res[3] = TMath::Sqrt(MCTrack1_Mass*MCTrack1_Mass + Pf_Res[0]*Pf_Res[0] + Pf_Res[1]*Pf_Res[1] + Pf_Res[2]*Pf_Res[2]);
            
      Erel1n = Estar1n(MCTrack2_Mass,MCTrack1_Mass,PNeutron,Pf_Res);

    }

    if(NEUTRON_X->size()>1){
      Double_t PNeutron[4];
      PNeutron[0] = NEUTRON_PX->at(0);
      PNeutron[1] = NEUTRON_PY->at(0);
      PNeutron[2] = NEUTRON_PZ->at(0);
      PNeutron[3] = NEUTRON_E->at(0);
      Double_t PNeutron1[4];
      PNeutron1[0] = NEUTRON_PX->at(1);
      PNeutron1[1] = NEUTRON_PY->at(1);
      PNeutron1[2] = NEUTRON_PZ->at(1);
      PNeutron1[3] = NEUTRON_E->at(1);
      Double_t Pf_Res[4];
      Double_t Pf_Tot = TMath::Sqrt(Pf[0]*Pf[0] + Pf[1]*Pf[1] + Pf[2]*Pf[2]);
      Double_t Pf_Tot_Res = RandX->Gaus(Pf_Tot,MOM_Res);
      Pf_Res[0] = Pf_Tot_Res*Pf[0]/Pf_Tot;
      Pf_Res[1] = Pf_Tot_Res*Pf[1]/Pf_Tot;
      Pf_Res[2] = Pf_Tot_Res*Pf[2]/Pf_Tot;
      Pf_Res[3] = TMath::Sqrt(MCTrack1_Mass*MCTrack1_Mass + Pf_Res[0]*Pf_Res[0] + Pf_Res[1]*Pf_Res[1] + Pf_Res[2]*Pf_Res[2]);
            
      Erel2n = Estar2n(MCTrack2_Mass,MCTrack1_Mass,PNeutron,PNeutron1,Pf_Res);
    }
    
    OutTree->Fill();
  }//End of Event Loop
  
  cout << " " << endl;

  OutTree->Write();
  output->Close();

}

//Main function
int main(Int_t argc, Char_t* argv[]){
  Char_t *in_filename=0;
  Bool_t NeedHelp = kTRUE;
  Bool_t MakeCalMap = kFALSE;
  Char_t *WhichNeutronWall=0;
  
  if (argc > 1){
    for (Int_t i = 0; i < argc; i++){
      if (strncmp(argv[i],"--file=",7) == 0){
	in_filename = argv[i]+7;
	//NeedHelp = kFALSE;
      }
      if (strncmp(argv[i],"--CalMap",8) == 0){
	MakeCalMap = kTRUE;
	//NeedHelp = kFALSE;
      }
      if (strncmp(argv[i],"--NWall=",8) == 0){
	WhichNeutronWall = argv[i]+8;
	NeedHelp = kFALSE;
      }
    }
  }
  if (NeedHelp){
    cout << "\nUse the following flags: ";
    cout << "\n--file=/path/to/your/file/filename.root   : (mandatory) simulation input file";
    cout << "\n--NWall=Neuland/Mona/Land   : (mandatory) Choice of Neutron Wall";
    cout << "  " << endl;
    return 0;
  }
  
  invariant_mass(in_filename,MakeCalMap,WhichNeutronWall);
  return 0;
}
