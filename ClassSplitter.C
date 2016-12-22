void ClassSplitter(){
   
TFile *oldfile = new TFile("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/Xicc2D+pK_TMVAresults.root");
TTree *oldtree = (TTree*)oldfile->Get("TrainTree");
Double_t nEntries = oldtree->GetEntries();
 
TFile *Sig_File = new TFile("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/Xicc2D+pK_TMVAresults_Signal_MagUp.root","recreate");
TTree *newtree1 = oldtree->CloneTree(0);

Int_t classID;
Int_t Polarity;

Float_t BDT;
Float_t BDTG;
Float_t MLP;
Float_t MLPBNN;
  
oldtree->SetBranchAddress("classID",&classID);
oldtree->SetBranchAddress("Polarity",&Polarity);

oldtree->SetBranchAddress("BDT",&BDT);
oldtree->SetBranchAddress("BDTG",&BDTG);
oldtree->SetBranchAddress("MLP",&MLP);
oldtree->SetBranchAddress("MLPBNN",&MLPBNN);

int cnt1 = 0;

for (Long64_t i=0;i<nEntries; i++) {
    if(i%10000==0) cout << "Processing event " << i << "..." << endl;
    oldtree->GetEntry(i);
	if(classID ==0 && MLPBNN > 0.08 && Polarity == 1){cnt1++ && newtree1->Fill();}
}

newtree1->Print();
newtree1->AutoSave();
cout << "-------> Number of events are: " << cnt1 << endl;

delete oldfile;
delete Sig_File;

return;
}