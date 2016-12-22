
// General script to calculate optimimum cut on MVA response variable

////////////////////////////////////////////
/// Authur: Murdo Traill 
/// Date: 3/11/16
///////////////////////////////////////////

#include <TLegend.h>
#include <TLatex.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>

#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <numeric>

using namespace std;

void printUsage()
{
  cout << "Usage: <Path/to/TFile> <Path/to/Tuple> <MVANumber> <#Steps> <Sigma> <BkgRootFile>(Maybe)"  << endl;
  cout << "0 = BDT, 1 = BDTG, 2 = MLP, 3 = MLPBNN"  << endl;
}
 
int main(int argc, char *argv[])
{
  int nargs = argc - 1;
  if (nargs!=5)
  {
    printUsage();
    return 1;
  }

  ///////////// DEFINE INPUTS ///////////

  string FileName = argv[1];
  string TupleName = argv[2];
  //string BkgName = argv[6];

  string MVAnumber = argv[3];
  string nsteps = argv[4];
  string sigma = argv[5]; 

  double nSteps = atof(nsteps.c_str());
  double Sigma = atof(sigma.c_str());
  int MVANumber = atoi(MVAnumber.c_str());
  
  /////////////////// DATA INPUT HERE //////////////////

  TFile myFile(FileName.c_str(), "READ");                         
  TTree* myTree = (TTree*) myFile.Get(TupleName.c_str());  

  //////////////////////////////////////////////////////////

  int nEntries = myTree->GetEntries(); 

  //////////////////////////////////////////////////////////

  int nMVA = 4;

  float bkg_mass_min = 3200.;                  
  float bkg_mass_max = 3900.;

  float eff_denominator = 4007970.;               
  //float eff_denominator = 1.;  

  vector<double> vector1;
  vector<double> vector2;

  vector<double> vector3;
  vector<double> vector4;

  Float_t Xicc_M;
  Int_t classID;

  Float_t pplus_PIDCalibEffWeight;
  Float_t Ksoft_PIDCalibEffWeight;
  Float_t Kminus_PIDCalibEffWeight;
  Float_t piplus1_PIDCalibEffWeight;
  Float_t piplus2_PIDCalibEffWeight;

  Float_t Event_PIDCalibEffWeight;

  Float_t MVA_Array[nMVA-1]; 

  vector<string> MVAMethods;

  MVAMethods.push_back("BDT");
  MVAMethods.push_back("BDTG");
  MVAMethods.push_back("MLP");
  MVAMethods.push_back("MLPBNN");

  myTree->SetBranchAddress("Xicc_M",&Xicc_M);
  myTree->SetBranchAddress("classID",&classID);

  myTree->SetBranchAddress("pplus_PIDCalibEffWeight", &pplus_PIDCalibEffWeight);
  myTree->SetBranchAddress("Ksoft_PIDCalibEffWeight", &Ksoft_PIDCalibEffWeight);
  myTree->SetBranchAddress("Kminus_PIDCalibEffWeight", &Kminus_PIDCalibEffWeight);
  myTree->SetBranchAddress("piplus1_PIDCalibEffWeight", &piplus1_PIDCalibEffWeight);
  myTree->SetBranchAddress("piplus2_PIDCalibEffWeight", &piplus2_PIDCalibEffWeight);

  myTree->SetBranchAddress("Event_PIDCalibEffWeight", &Event_PIDCalibEffWeight);

  for (int i = 0; i < (nMVA); ++i) 
  {
  myTree->SetBranchAddress(MVAMethods[i].c_str(),&MVA_Array[i]);   
  }

  for (int i = 0; i<nEntries; i++) { 
       myTree->GetEntry(i);
       vector1.push_back(MVA_Array[MVANumber]);
       vector2.push_back(Xicc_M);
  }

  double MVA_Max = *max_element(vector1.begin(), vector1.end());
  //double MVA_Max = 1.0;

  //cout << MVA_Max << endl;

  double MVA_Min = *min_element(vector1.begin(), vector1.end());

  double MVA_Range = MVA_Max - MVA_Min;
  double Step_Value_MVA = MVA_Range/nSteps;

  double Mass_Max = *max_element(vector2.begin(), vector2.end());
  double Mass_Min = *min_element(vector2.begin(), vector2.end());

  cout << "" << endl;

  cout << "The minimum " + MVAMethods[MVANumber] + " value is: " << MVA_Min <<endl;
  cout << "The maximum " + MVAMethods[MVANumber] + " value is: " << MVA_Max <<endl;

  cout << "" << endl;

  cout << "The minimum Mass value is: " << Mass_Min << " MeV" << endl;
  cout << "The maximum Mass value is: " << Mass_Max << " MeV" << endl;

  cout << "" << endl;

  cout << "----> Starting optimisation for " << MVAMethods[MVANumber] << "!" << endl;
  cout << "" << endl;
  cout << "The number of events in the training tree is: " << nEntries << endl;
  
  int signal_counter;
  int background_counter;

  vector<float> weights;

  for (int j = 0; j<(nSteps); j++){

    signal_counter = 0;
    background_counter =0;
    weights.clear();

    double current_MVA_value = MVA_Min + (j * Step_Value_MVA);

    int loop_number = j+1;

    cout << "" << endl;
    cout << "Beginning loop " << loop_number << " of " << nSteps << " with MVA cut at " << current_MVA_value << endl;
    cout << "Looking at Signal Data" << endl;

    for (Int_t k = 0; k < nEntries; k++) {
        myTree->GetEntry(k);

        if(k%10000 == 0){cout << "Processing event " << k << "..." << endl;}
        if(classID == 0 && MVA_Array[MVANumber] > current_MVA_value ){
          //signal_counter++; <---------For no weighted method
          weights.push_back(Event_PIDCalibEffWeight);
        }
        if(classID == 1 && MVA_Array[MVANumber] > current_MVA_value && Xicc_M > bkg_mass_min && Xicc_M < bkg_mass_max){background_counter++;}
    }

    float sum_of_elems = std::accumulate(weights.begin(), weights.end(), 0.0);
    double current_signal = sum_of_elems;

    //double current_signal = signal_counter; <---------For no weighted method
    
    double current_background = background_counter;

    ///////////// DEFINE PUNZI MERIT /////////////////////

    double punzi_numerator = current_signal/eff_denominator;
    double punzi_denominator = Sigma/2. + sqrt(current_background);                     
    double punzi_value = punzi_numerator/punzi_denominator;

    ///////////////////////////////////////////////////////

    if(punzi_value > 0){
    vector3.push_back(current_MVA_value);
    vector4.push_back(punzi_value);
    }
    else{cout << "...WARNING: No data passed this cut!" << endl;}
  }   

  ////////////////////////////////////////////////////////////

  cout << "" << endl;

  double punzi_max = *max_element(vector4.begin(), vector4.end());

  int entry_max = distance(vector4.begin(), max_element(vector4.begin(), vector4.end())) + 1;

  for (int t=0; t < vector3.size(); t++)
  {
  int entry_number1 = t+1;
  cout << " Entry " << entry_number1 << ": " << vector3[t] << "   " << vector4[t] << endl;
  }

  cout << "" << endl;

  double entry_max_MVA = vector3[entry_max -1];

  cout << "" << endl;
  cout << "Punzi Maximum is " << punzi_max << " at entry " << entry_max << ": " << entry_max_MVA  << endl;

  //////////////// BEGIN FoM PLOT ///////////////////////////////
  ///////////////////////////////////////////////////////////////

  cout << "" << endl;
  cout << "---> Beginning FoM Plotting!" << endl;

  gStyle->SetTitleY(0.97);  
  gStyle->SetTitleX(0.5);
  gStyle->SetOptStat(1110);
  gStyle->SetStatX(0.35);
  gStyle->SetStatY(0.6);
  gStyle->SetCanvasColor(0);
  gStyle->SetFuncWidth(2);

  TCanvas* c1 = new TCanvas("c1","c1",900,600); 
  c1->cd(1);

  TGraph * g1 = new TGraph(vector3.size(), &vector3[0], &vector4[0]);
 
  g1->SetLineColor(4);
  g1->SetLineWidth(2);
  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(21);

  string titlename = MVAMethods[MVANumber] + " with " + nsteps + " cuts";
  string xtitle = MVAMethods[MVANumber] + " Response Value";
  string ytitle = "Punzi FoM (#sigma = " + sigma + ")";  

  g1->SetTitle(titlename.c_str());
  g1->GetXaxis()->SetTitleOffset(1.2);
  //g1->GetXaxis()->CenterTitle();
  g1->GetXaxis()->SetTitle(xtitle.c_str());
  g1->GetYaxis()->SetTitle(ytitle.c_str());
  g1->GetYaxis()->SetTitleOffset(1.12);
  //g1->GetYaxis()->CenterTitle();

  g1->Draw("AC");

  TLatex Tl2;
  Tl2.SetNDC();
  Tl2.SetTextAlign(12);
  Tl2.SetTextSize(0.03);
  Tl2.DrawLatex(0.15,0.87,Form("Optimum cut > ~ %0.2f with FoM ~ %0.3e", entry_max_MVA, punzi_max));

  ///////////// BEGIN OPTIMISED MASS PLOT ///////////////////////
  ////////////////////////////////////////////////////////////////                       

  cout << "---> Beginning Optimised Mass Plot!" << endl;

  TCanvas* c2 = new TCanvas("c2","c2",900,600); 
  c2->cd(1);

  float nSig = 0;
  float nBkg = 0;

  for (int m = 0; m<nEntries; m++) { 
       myTree->GetEntry(m);
       if(classID == 0){nSig++;}
       if(classID == 1){nBkg++;}
  }

  cout << "" << endl;

  cout << "Total Number of Signal events in ntuple is: " << nSig << endl;
  cout << "Total number of Backgroud events in ntuple is: " << nBkg << endl;

  //////////////////////////////////////////////////////

  float cs_min = 100e-9;  // 1000nb  
  float cs_max = 1000e-9;   // 100nb

  /*
  TH1F * h2 = new TH1F("h2","h2",100,0,20);

  TFile myBkgFile(BkgName.c_str(), "READ");                         
  TTree* LumiTree = (TTree*) myBkgFile.Get("GetIntegratedLuminosity/LumiTuple");  

  Double_t IntegratedLuminosity;
  LumiTree->SetBranchAddress("IntegratedLuminosity",&IntegratedLuminosity);

  int nEntries2 = LumiTree->GetEntries(); 

  for (int i = 0; i <nEntries2 ; ++i) 
  {
      LumiTree->GetEntry(i);
      h2->Fill(IntegratedLuminosity);
  }
  h2->Draw();
  double integral = h2->Integral();
  double mean = h2->GetMean();

  double lumi_int = (integral*mean)/1000.;
  */
  // Put in manually

  double lumi_int = 1.436;

  float Br_min = 0.01;  // 1%              
  float Br_max = 0.05;  // 5%
  
  float Total_Comb_Bkg = 2.10236e7;   // From WC data

  ////////////////////////////////////////////////////////

  TH1F * h1 = new TH1F("h1","h1",100,3,4);

  cout << "" << endl;

  int cnt1 = 0;
  int cnt2 = 0;

  for (Int_t p = 0; p < nEntries; p++) {
        myTree->GetEntry(p);

      if(p%100000 == 0){cout << "---->Processing"<<endl;}
       
      if (MVA_Array[MVANumber] > entry_max_MVA){
        if(classID == 0){cnt1++;}
        }
      }

  for (Int_t p = 0; p < nEntries; p++) {
        myTree->GetEntry(p);

      if(p%100000 == 0){cout << "---->Processing"<<endl;}
       
      if (MVA_Array[MVANumber] > entry_max_MVA){
        if(classID == 1){cnt2++;}
        }
      }   

  double sig_eff = cnt1*100/nSig;
  double bkg_rej = 100.0-(cnt2*100)/nBkg;

  cout << "" << endl;

  cout << "Signal numbers after the count is: " << cnt1 << endl;
  cout << "Background numbers after the count is: " << cnt2 <<endl;

  cout << "" << endl;

  printf ("Signal efficiency is: %0.2f %% \n", sig_eff);
  printf ("Background rejection is: %0.2f %% \n", bkg_rej);

  float N_Xicc_min = (sig_eff/100.) * lumi_int * cs_min * Br_min * 1e9;
  int N_Xicc_min_rounded = N_Xicc_min;

  float N_Xicc_max = (sig_eff/100.) * lumi_int * cs_max * Br_max * 1e9;
  int N_Xicc_max_rounded = N_Xicc_max;

  cout << "" << endl;
  printf("Number of expected signal Xicc events in %0.2f fb-1 of WC data is between %0.2f and %0.2f events \n", lumi_int,  N_Xicc_min,N_Xicc_max); 
  cout << "" << endl;

   for (Int_t p = 0; p < nEntries; p++) {
        myTree->GetEntry(p);

      //if(p%100000 == 0){cout << "---->Processing"<<endl;}
       
      if (MVA_Array[MVANumber] > entry_max_MVA){
         float Mass =Xicc_M/1000.;
        if(classID == 0){cnt1++ && h1->Fill(Mass);}
        }
      }

  double weight = (Total_Comb_Bkg/N_Xicc_max_rounded)*cnt1;

  for (Int_t p = 0; p < nEntries; p++) {
        myTree->GetEntry(p);

      //if(p%100000 == 0){cout << "---->Processing"<<endl;}
       
      if (MVA_Array[MVANumber] > entry_max_MVA){
         float Mass =Xicc_M/1000.;
        if(classID == 1){cnt2++ && h1->Fill(Mass,weight);}
        }
      }   

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  h1->SetTitle("Optimised with Punzi FoM");
  
  h1->GetXaxis()->SetTitleOffset(1.2);
  //h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetTitle("m(D^{+}K^{-}p) [GeV]");

  h1->GetYaxis()->SetTitle("Number of Events / 10 MeV/c^{2}");
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->GetYaxis()->CenterTitle();
  
  h1->SetLineColor(1);   
  h1->SetLineWidth(2);

  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(19); 

  h1->SetFillStyle(3001);
  h1->SetFillColor(kRed +1);

  h1->Draw("HIST B");

  ///////////////////////////////////////////////////

  TLatex Tl1;
  Tl1.SetNDC();
  Tl1.SetTextAlign(12);
  Tl1.SetTextSize(0.05);
  Tl1.DrawLatex(0.15,0.7,"#scale[1.1]{#Xi^{+}_{cc}} #rightarrow D^{+} K^{-} p");
  Tl1.DrawLatex(0.15,0.8,"LHCb Unoffical");

  /////////////// SAVING FILES AND OBJECTS ////////////////////

  string output1 =  "./RootFiles/PunziFoM_" + MVAMethods[MVANumber] + "_" + sigma + "Sigma_" + nsteps + "steps.root";
  TFile *output = new TFile(output1.c_str(), "RECREATE");

  cout << "Created root file: " << output1 << endl;

  g1->Write(); 
  h1->Write();

  string fig1 = "./Figures/FoM_Plots/PunziFoM_" + MVAMethods[MVANumber] + "_" + sigma + "Sigma_" + nsteps + "steps.pdf"; 
  string fig2 = "./Figures/Mass_Plots/CleanedUpMass_" + MVAMethods[MVANumber] + "_" + sigma + "Sigma_" + nsteps + "steps.pdf"; 

  c1->Print(fig1.c_str());  
  c2->Print(fig2.c_str());

  myFile.Close();  

  return 0;

} 

