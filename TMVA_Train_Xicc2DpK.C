// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : Selection of XiccDpK                                               *
 * Package   : TMVA                                                               *
 * Root Macro: TMVA_Train_Xicc2DpK                                                *
 *                                                                                *
 * Input: /Users/austin/work/ntuple/data/Background.root                          *
 *        /Users/austin/work/ntuple/mc/DecayTree.root                             *
 *                                                                                *
 * Specify method: root -l ./'TMVAClassification.C("BDT, BDTG")'                  *
 *                                                                                *
 * Output: TMVAoutput_XiccDpK.root                                                *
 * Date:   31/10/2016                                                             *
 *********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVA_Train_Xicc2DpK( TString myMethodList = "" )
{
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 1; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ------------------------------------------------------------------------------- //

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "Xicc2D+pK_TMVAresults.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   // The second argument is the output file for the training results
   
   TMVA::Factory *factory = new TMVA::Factory( "myTMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   ////////////// ADDING TRAINING VARIABLES ////////////

   factory->AddVariable("Xicc_MAXDOCA","#Xi^{+}_{cc} MaxDOCA","",'F'); // 1
   factory->AddVariable("Xicc_ConsPV_chi2 := Xicc_ConsPV_chi2[0]/Xicc_ConsPV_nDOF[0]","#Xi^{+}_{cc} PV #chi^{2}/ndof (Fixed PV)","",'F'); // 2 
   factory->AddVariable("Xicc_ENDVERTEX_CHI2_NDOF := Xicc_ENDVERTEX_CHI2/Xicc_ENDVERTEX_NDOF ","#Xi^{+}_{cc} DecayVertex #chi^{2}/ndof","",'F'); // 3

   factory->AddVariable("pplus_PT","proton PT", "MeV/c", 'F'); //4
   factory->AddVariable("Ksoft_PT","K^{-} PT", "MeV/c", 'F'); //5
   //factory->AddVariable("Xicc_Daughter_Min_PT := min(Dplus_PT,min(pplus_PT,Ksoft_PT))","#Xi^{+}_{cc} Daughter MinPT","MeV/c",'F'); 

   factory->AddVariable("Dplus_MAXDOCA","D^{+} MaxDOCA","",'F'); // 6

   factory->AddVariable("Xicc_IPCHI2_OWNPV","#Xi^{+}_{cc} IP #chi^{2}","",'F'); // 7 
   factory->AddVariable("Dplus_IPCHI2_OWNPV","D^{+} IP #chi^{2}","",'F'); // 8
   
   //factory->AddVariable("Dplus_ENDVERTEX_CHI2","Dplus_ENDVERTEX_CHI2","",'F'); 
   //factory->AddVariable("log(1-Xicc_DIRA_OWNPV):= log(1-Xicc_DIRA_OWNPV)","log(1-Xicc_DIRA_OWNPV)", "", 'F');
   //factory->AddVariable("Dplus_FDCHI2_OWNPV","Dplus_FDCHI2_OWNPV","",'F'); 

   ////////////// ADDING SPECTATORS ///////////
   
   // These are not used to train the MVA methods but are added added as brnaches in the final n-tuple 

   factory->AddSpectator("Xicc_PT", "Xicc_PT", "MeV", 'F' );
   factory->AddSpectator("Xicc_P",  "Xicc_P", "MeV", 'F' );
   factory->AddSpectator("Xicc_M",  "Xicc_M", "MeV", 'F' );
   factory->AddSpectator("Xicc_MMERR",  "Xicc_MMERR", "MeV", 'F' );
   factory->AddSpectator("Xicc_TAU", "Xicc_TAU", "ns", 'F');
   factory->AddSpectator("Xicc_TAUERR", "Xicc_TAUERR", "", 'I');

   factory->AddSpectator("Dplus_PT", "Dplus_PT", "MeV", 'F' );
   factory->AddSpectator("Dplus_P",  "Dplus_P", "MeV", 'F' );
   factory->AddSpectator("Dplus_M",  "Dplus_M", "MeV", 'F' );
   factory->AddSpectator("Dplus_MMERR",  "Dplus_MMERR", "MeV", 'F' );
   factory->AddSpectator("Dplus_TAU", "Dplus_TAU", "ns", 'F');
   factory->AddSpectator("Dplus_TAUERR", "Dplus_TAUERR", "ns", 'F');

   factory->AddSpectator("Ksoft_PIDK", "Ksoft_PIDK", "", 'F');
   factory->AddSpectator("Ksoft_PIDp", "Ksoft_PIDp", "", 'F');
   factory->AddSpectator("Ksoft_PIDe", "Ksoft_PIDe", "", 'F');
   factory->AddSpectator("Ksoft_PIDmu", "Ksoft_PIDmu", "", 'F');
   factory->AddSpectator("Ksoft_P", "Ksoft_P", "MeV", 'F');
   factory->AddSpectator("Ksoft_PT", "Ksoft_PT", "MeV", 'F');
   factory->AddSpectator("Ksoft_PX", "Ksoft_PX", "MeV", 'F');
   factory->AddSpectator("Ksoft_PY", "Ksoft_PY", "MeV", 'F');
   factory->AddSpectator("Ksoft_PZ", "Ksoft_PZ", "MeV", 'F');

   factory->AddSpectator("pplus_PIDK","pplus_PIDK", "", 'F');
   factory->AddSpectator("pplus_PIDp", "pplus_PIDp", "", 'F');
   factory->AddSpectator("pplus_PIDe", "pplus_PIDe", "", 'F');
   factory->AddSpectator("pplus_PIDmu", "pplus_PIDmu", "", 'F');
   factory->AddSpectator("pplus_P", "pplus_P", "MeV", 'F');
   factory->AddSpectator("pplus_PT", "pplus_PT", "MeV", 'F');
   factory->AddSpectator("pplus_PX", "pplus_PX", "MeV", 'F');
   factory->AddSpectator("pplus_PY", "pplus_PY", "MeV", 'F');
   factory->AddSpectator("pplus_PZ", "pplus_PZ", "MeV", 'F');

   factory->AddSpectator("pplus_ProbNNp", "pplus_ProbNNp", "", 'F');
   factory->AddSpectator("pplus_ProbNNk", "pplus_ProbNNk", "", 'F');
   factory->AddSpectator("pplus_ProbNNpi", "pplus_ProbNNpi", "", 'F');
   factory->AddSpectator("pplus_ProbNNghost", "pplus_ProbNNghost", "", 'F');

   factory->AddSpectator("Ksoft_ProbNNp", "Ksoft_ProbNNp", "", 'F');
   factory->AddSpectator("Ksoft_ProbNNk", "Ksoft_ProbNNk", "", 'F');
   factory->AddSpectator("Ksoft_ProbNNpi", "Ksoft_ProbNNpi", "", 'F');
   factory->AddSpectator("Ksoft_ProbNNghost", "Ksoft_ProbNNghost", "", 'F');

   factory->AddSpectator("Kminus_ProbNNp", "Kminus_ProbNNp", "", 'F');
   factory->AddSpectator("Kminus_ProbNNk", "Kminus_ProbNNk", "", 'F');
   factory->AddSpectator("Kminus_ProbNNpi", "Kminus_ProbNNpi", "", 'F');
   factory->AddSpectator("Kminus_ProbNNghost", "Kminus_ProbNNghost", "", 'F');

   factory->AddSpectator("piplus1_ProbNNp", "piplus1_ProbNNp", "", 'F');
   factory->AddSpectator("piplus1_ProbNNk", "piplus1_ProbNNk", "", 'F');
   factory->AddSpectator("piplus1_ProbNNpi", "piplus1_ProbNNpi", "", 'F');
   factory->AddSpectator("piplus1_ProbNNghost", "piplus1_ProbNNghost", "", 'F');

   factory->AddSpectator("piplus2_ProbNNp", "piplus2_ProbNNp", "", 'F');
   factory->AddSpectator("piplus2_ProbNNk", "piplus2_ProbNNk", "", 'F');
   factory->AddSpectator("piplus2_ProbNNpi", "piplus2_ProbNNpi", "", 'F');
   factory->AddSpectator("piplus2_ProbNNghost", "piplus2_ProbNNghost", "", 'F');

   factory->AddSpectator("pplus_PIDCalibEffWeight", "pplus_PIDCalibEffWeight", "", 'F' );
   factory->AddSpectator("Ksoft_PIDCalibEffWeight", "Ksoft_PIDCalibEffWeight", "", 'F' );
   factory->AddSpectator("Kminus_PIDCalibEffWeight", "Kminus_PIDCalibEffWeight", "", 'F' );
   factory->AddSpectator("piplus1_PIDCalibEffWeight", "piplus1_PIDCalibEffWeight", "", 'F' );
   factory->AddSpectator("piplus2_PIDCalibEffWeight", "piplus2_PIDCalibEffWeight", "", 'F' );

   factory->AddSpectator("Event_PIDCalibEffWeight", "Event_PIDCalibEffWeight","", 'F');

   factory->AddSpectator("Polarity", "Polarity","", 'I');

   TString sfname("/Users/murdotraill/Desktop/Data/MonteCarlo/GenXicc/Xicc+2D+K-p/latest3/Xicc2DpK_BothMag_TMVASignal.root");
   TString bfname("/Users/murdotraill/Desktop/Data/Real/WrongChargeLines/Xicc+2D+K-p/latest2/Xicc+2DKp_MagDown_BKG_Preselection_withEta_cutonEta_Eff.root");    
   
   TFile *sinput = TFile::Open( sfname );
   TFile *binput = TFile::Open( bfname );
   
   std::cout << "--- TMVAClassification       : Using input signal file: " << sinput->GetName() << std::endl;
   std::cout << "--- TMVAClassification       : Using input background file: " << binput->GetName() << std::endl;
   
   // --- Register the training and test trees

   TTree *signal     = (TTree*)sinput->Get("DecayTree");
   TTree *background = (TTree*)binput->Get("DecayTree");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   
   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //     
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //   
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration 

   // <------------------------------------ Look at this on monday!!

   // Set individual event weights (the variables must exist in the original TTree)
   // So need to create a branch with branch names  weight1 and weight2
  
   factory->SetSignalWeightExpression("Event_PIDCalibEffWeight"); //<-------- NEW!!!

   //factory->SetBackgroundWeightExpression("weight1*weight2");

   // Apply additional cuts on the signal and background samples (can be different)
   
   TCut mycuts = ""; 
   
   TCut mycutb = "Xicc_M > 3000 && Xicc_M < 4000"; 

   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   factory->PrepareTrainingAndTestTree(mycuts, mycutb,"SplitMode=Random:NormMode=NumEvents:!V");

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NegWeightTreatment=NoNegWeightsInTraining:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList; 
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(","); 
      methodList += regMethod;
   }
   return TMVA_Train_Xicc2DpK(methodList); 
}
