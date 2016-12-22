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

void AllPunzi(){

TFile BDT_File("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/RootFiles/PunziFoM_BDT_5Sigma_100steps.root");
TFile BDTG_File("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/RootFiles/PunziFoM_BDTG_5Sigma_100steps.root");
TFile MLP_File("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/RootFiles/PunziFoM_MLP_5Sigma_100steps.root");
TFile MLPBNN_File("/Users/murdotraill/Desktop/Scripts/TMVAStudies/Xicc+2D+K-p/Attempt4/RootFiles/PunziFoM_MLPBNN_5Sigma_100steps.root");

gStyle->SetTitleY(0.97);  
gStyle->SetTitleX(0.53);
gStyle->SetCanvasColor(0);
gStyle->SetOptStat(0);    
gStyle->SetStatX(0.85);
gStyle->SetStatY(0.59);                 
gStyle->SetLegendTextSize(0.03);

TCanvas* c1 = new TCanvas("c1","c1",900,600);
  
c1->cd(1);

TGraph * g1 = new TGraph(); 
g1 = (TGraph*)BDT_File.Get("Graph"); 

TGraph * g2 = new TGraph(); 
g2 = (TGraph*)BDTG_File.Get("Graph"); 

TGraph * g3 = new TGraph(); 
g3 = (TGraph*)MLP_File.Get("Graph"); 

TGraph * g4 = new TGraph(); 
g4 = (TGraph*)MLPBNN_File.Get("Graph"); 

g2->GetYaxis()->SetRangeUser(0,5.2e-6);

g2->Draw("AC");
g2->GetXaxis()->SetLimits(0.,1.);
g2->SetTitle("");

g2->GetXaxis()->SetTitleFont(132);
g2->GetXaxis()->SetLabelFont(132);
g2->GetXaxis()->SetTitleOffset(1.2);

g2->GetYaxis()->SetTitleFont(132);
g2->GetYaxis()->SetLabelFont(132);
g2->GetYaxis()->SetTitleOffset(0.7);

g2->GetYaxis()->SetTitle("Punzi FoM (#sigma = 5)");
g2->GetXaxis()->SetTitle("MVA response variable (x)");

c1->Update();
g2->SetLineColor(kRed);

g1->Draw("same");
g1->SetLineColor(kBlue);

g3->Draw("same");
g3->SetLineColor(kGreen+1);

g4->Draw("same");
g4->SetLineColor(kBlack);

TLegend *leg1 = new TLegend(0.62,0.67,0.9,0.9);
leg1->AddEntry(g1,"BDT; 0.11; 23.39%","L");
leg1->AddEntry(g2,"BDTG; 0.51; 39.80%","L");
leg1->AddEntry(g3,"MLP; 0.15; 17.39%","L");
leg1->AddEntry(g4,"MLPBNN; 0.08; 27.92%","L");
leg1->Draw();

TLatex Tl1;
Tl1.SetNDC();
Tl1.SetTextAlign(12);
Tl1.SetTextSize(0.035);
Tl1.DrawLatex(0.4,0.85,"#font[12]{FoM(x) = #frac{#epsilon(x)}{#sigma /2 + #sqrt{B(x)}}}");

c1->Update();
c1->Print("./AllPunzi.png");

return;
}