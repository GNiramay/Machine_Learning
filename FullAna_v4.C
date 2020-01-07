// Program compare following ML models with weight rechit energy distribtuions
// 1. simple chisquare minimization.
// 2. chisquare and std. dev. minimization for each iteration
// 3. chisquare and std. dev. minimization one traninng after the other

// List of plots made:
// 1. Weights vs layer no. -> 8 plots for 8 energies
// 2. Reconstructed energy plots (without scaling)
// 3. Reconstructed energy plots (scaling mean to beam energy)
// 3. Gaussian fit plots
// 4. Energy resolution vs. beam energy -> 2 plots for sigma/mu and truncated RMS/Mean (truncated RMS/Mean not yet done)

#include"Common_Resources.C"
// #include"Results.h"
#include"Results_v2.h"

TString NOW = "02_01_20";	// Today's date
TString InDate = "03_12_19";	// Date of input files->ML-layerwise
TString DetwiseDate = "28_11_19";	// Date of input files->ML-detectorwise
TString cf = "Chi2";		// cost function
float ENRES[3][4][8];		// 3->EnRes methods, 4->ML models, 8-> energy samples

void WeightsPlot();
void RecoEn(bool);
void EnRes_All(bool);		// Plot EnRes using all the methods
void FitPlots(TH1F*,TString,TString,int);
// void Fit_EnRes(TGraph*);	// Eneregy resolution parametrization
void Fit_EnRes(TGraph*, TString , TString, TString );

void Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2);
double Find_EnRes(TH1F* HH);
void GetEnRes(TH1F*,int,int); // to fill energy resolution array
void EnResponse();	      // To plot energy response

bool RelWt = true;		// whether to plot relative weights (like alpha, beta, gamma etc.)

void FullAna_v4()
{

  // RecoEn(false);		// No mean scaling
  // RecoEn(true);			// mean scaled to beam energy
  // EnRes_All(true);
  EnResponse();

  // EnRes_All("Chi2","#chi^{2}", true);
  // // EnRes_All("RMSMean","#frac{Std. dev}{Mean}", true);
  // // EnRes_All("OnlyRMS","Std. dev", true);
  // // WeightsPlot();
}

void WeightsPlot()
{
  float ML_Wt[8][79];
  float LayerNo[79];  for(int i = 0;i<79;i++) LayerNo[i]=i+1;
  gStyle->SetTextFont(122);
  // gStyle->SetTextSize(11);

  ifstream iif("WeightsInfo.txt");
  float E[8],A[8],B[8],GE[8],GF[8];
  for(int i=0;i<8;i++) iif>>E[i]>>A[i]>>B[i]>>GE[i]>>GF[i];

  for(int i=0;i<8;i++)
    {
      TString fname = "../OUT_TXT/Weights_"+cf+"_"+InDate+"_E"+SEnergy[i]+".txt";
      ifstream iif(fname);
      float temp2;
      for(int j=0;j<79;j++){iif>>temp2; ML_Wt[i][j]=temp2;}

      if(RelWt){
	//////////// Relative weights
	for(int k=0;k<79;k++)if(k!=19) ML_Wt[i][k]=ML_Wt[i][k]/ML_Wt[i][19]; // Dividing everything by wt. of 20th layer
	ML_Wt[i][19]=1;
	/////////// End
      }

      TCanvas* tc = new TCanvas("aa","bb",800,600);
      TLegend* tl = new TLegend(0.1,0.7,0.4,0.9);
      TLine* ShWt[3];

      TGraph* tg = new TGraph(79,LayerNo,ML_Wt[i]);
      tg_Beauty(tg,"sim Pion | "+SEnergy[i]+" GeV | FTFP_BERT_EMN", "E"+SEnergy[i],4,21);
      tg->GetXaxis()->SetTitle("Layer index");
      tg->GetYaxis()->SetTitle("Weight");
      tl->AddEntry(tg,"ML weights","p");

      if(RelWt){
	//////////// Relative weights
	ShWt[0] = new TLine(1,1.0,28,1); ShWt[0]->SetLineColor(1);ShWt[0]->SetLineWidth(4);
	ShWt[1] = new TLine(29,B[i],40,B[i]); ShWt[1]->SetLineColor(2);ShWt[1]->SetLineWidth(4);
	ShWt[2] = new TLine(41,A[i]*B[i],79,A[i]*B[i]); ShWt[2]->SetLineColor(5);ShWt[2]->SetLineWidth(4);
	tl->AddEntry(ShWt[0],"reference. 1","l");
	tl->AddEntry(ShWt[1],"Shubham's #beta","l");
	tl->AddEntry(ShWt[2],"Shubham's #beta#alpha","l");

	tg->GetYaxis()->SetTitle("Relative weight(#frac{W_{i}}{W_{20}})");
	tg->GetYaxis()->SetRangeUser(-1.0,10.0);
	//////////// End
      }
      else{
	ShWt[0] = new TLine(1,1/GE[i],28,1/GE[i]); ShWt[0]->SetLineColor(1);ShWt[0]->SetLineWidth(4);
	ShWt[1] = new TLine(29,B[i]/GE[i],40,B[i]/GE[i]); ShWt[1]->SetLineColor(2);ShWt[1]->SetLineWidth(4);
	ShWt[2] = new TLine(41,A[i]*B[i]/GE[i],79,A[i]*B[i]/GE[i]); ShWt[2]->SetLineColor(5);ShWt[2]->SetLineWidth(4);
	tl->AddEntry(ShWt[0],"Shubham's #gamma","l");
	tl->AddEntry(ShWt[1],"Shubham's #gamma#beta","l");
	tl->AddEntry(ShWt[2],"Shubham's #gamma#beta#alpha","l");
      }

      tg->Draw("apl");
      tc->SetGrid();
      tc->SetTicky(2);
      ShWt[0]->Draw("same");
      ShWt[1]->Draw("same");
      ShWt[2]->Draw("same");
      tl->Draw("same");

      if(RelWt){
	// SaveMe(tc,"Relative_Weights_"+NOW+"_CF_"+cf+"_E"+SEnergy[i],false);
	delete tc;
      }
      else{
	// SaveMe(tc,"Weights_"+NOW+"_CF_"+cf+"_E"+SEnergy[i],false);
	delete tc;      
      }

    }
  
}

void RecoEn(bool DoScale)
{
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.3);
  double lyr[8],det[8],shb[8];

  for(int i=0;i<8;i++){
    TCanvas* tc = new TCanvas("aa","bb",800,600);

    TFile* tempf = new TFile("../Assignment_3_UsingSS/IN_ROOT/Pion_"+InDate+"_Sim_"+SEnergy[i]+"_ForML.root");
    Results* Rs_EE[3];
    Results* Rs_FH[3];
    Results* Rs[3];
    TH1F* h_Shb;

    Rs_EE[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_EE.root"); // One-step
    Rs_EE[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_EE_TwoStep.root"); // Two-step ML
    Rs_EE[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_EE_Twice_ML.root"); // Twice ML

    Rs_FH[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_FH.root"); // One-step
    Rs_FH[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_FH_TwoStep.root"); // Two-step ML
    Rs_FH[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_FH_Twice_ML.root"); // Twice ML

    for(int j=0;j<3;j++) Rs[j] = Rs_EE[j]->CombineIt(Rs_FH[j]);

    float EMIN[3] = {Rs[0]->emin,Rs[1]->emin,Rs[2]->emin};
    float EMAX[3] = {Rs[0]->emax,Rs[1]->emax,Rs[2]->emax};
    
    float xmin = 0;
    float xmax = 2*FEnergy[i];
    if(!DoScale){
      for(int j=0;j<3;j++){
	if(EMIN[j]<xmin) xmin = EMIN[j];
	if(EMAX[j]>xmax) xmax = EMAX[j];
      }
    }// cout<<xmin<<"\t"<<xmax<<endl;

    // TH1F* h_EE = Rs[0]->GetHist(100,0,4*FEnergy[i],DoScale);
    // TH1F* h_FH = Rs[1]->GetHist(100,0,4*FEnergy[i],DoScale);
    // TH1F* h_All = Rs[2]->GetHist(100,0,4*FEnergy[i],DoScale);

    float tempN = 50*(xmax-xmin)/FEnergy[i];
    // float tempN = 100*(xmax-xmin)/FEnergy[i];
    // // xmin = -100;
    // // float tempN = 100;

    TH1F* h_One = Rs[0]->GetHist((int)tempN,xmin,2*xmax,DoScale);
    TH1F* h_TwoStep = Rs[1]->GetHist((int)tempN,xmin,2*xmax,DoScale);
    TH1F* h_Twice = Rs[2]->GetHist((int)tempN,xmin,2*xmax,DoScale);

    TString ttl,imtitle;

    if(DoScale){h_Shb = (TH1F*)tempf->Get("ScaleTotRecEn");ttl="(means scaled to beam energy)";imtitle="Scaled";}
    else{h_Shb = (TH1F*)tempf->Get("TotRecEn");ttl="(No scaling)";imtitle="NoScale";}

    Beauty(h_One,"Beam Energy "+SEnergy[i]+" GeV"+ttl, "One loss function",2,1,0.75,0.9);
    Beauty(h_TwoStep,"Beam Energy "+SEnergy[i]+" GeV"+ttl, "two-step training",4,1,0.6,0.75);
    Beauty(h_Twice,"Beam Energy "+SEnergy[i]+" GeV"+ttl, "Twice training",6,1,0.45,0.6);
    Beauty(h_Shb,"Beam Energy "+SEnergy[i]+" GeV"+ttl, "Weighted rechit energy",1,1,0.3,0.45);

    h_Shb->Rebin();
    h_One->DrawNormalized("",1);
    h_TwoStep->DrawNormalized("sames hist",1);
    h_Twice->DrawNormalized("sames hist",1);
    h_Shb->DrawNormalized("sames hist",1);

    SaveMe(tc,"RecoEn_E"+SEnergy[i]+"_"+NOW+"_"+imtitle,true);
    delete tc;
  }
}

// void h1_Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms)
void Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2)
{
  h->SetName(hname);
  h->SetMarkerColor(mc);
  h->SetLineColor(mc);
  h->SetMarkerSize(1);
  h->SetMarkerStyle(ms);
  h->SetLineWidth(2);

  h->SetTitle(htitle);
  h->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetYaxis()->SetTitle("No. of events");

  h->Draw();
  gPad->Update();
  TPaveStats* st = (TPaveStats*)h->FindObject("stats");

  st->SetTextColor(mc);
  st->SetY1NDC(Y1);
  st->SetY2NDC(Y2);
}

double Find_EnRes(TH1F* HH)
{
  Double_t hmean = HH->GetMean();
  Double_t hsgm = HH->GetStdDev();
  // Double_t rng = 2;
  Double_t rng = 1.5;

  // TFitResultPtr rr2 = HH->Fit("gaus","SQNR","",hmean-rng*hsgm,hmean+rng*hsgm);
  TFitResultPtr rr2 = HH->Fit("gaus","SQR","",hmean-rng*hsgm,hmean+rng*hsgm);

  double MEAN=rr2->Parameter(1);
  double SIGMA=rr2->Parameter(2);
  return(SIGMA/MEAN);
}

void EnRes_All(bool DoScale)
{
  TString meth[3]={"Sigma/Mean","90% truncated RMS/Mean", "RMS/Mean"};
  TString mtitle[3]={"Gauss","Trunc","RMS"};

  for(int i=0;i<8;i++){
    TString htitle = "Sim Pion | Beam Energy "+SEnergy[i]+" GeV";
    // if(DoScale) htitle = htitle+" | (rescaled)";

    TFile* tempf = new TFile("../Assignment_3_UsingSS/IN_ROOT/Pion_"+InDate+"_Sim_"+SEnergy[i]+"_ForML.root");
    Results* Rs_EE[3];
    Results* Rs_FH[3];
    Results* Rs[3];
    TH1F* h_Shb;

    Rs_EE[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_EE.root"); // One-step
    Rs_EE[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_EE_TwoStep.root"); // Two-step ML
    Rs_EE[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_EE_Twice_ML.root"); // Twice ML

    Rs_FH[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_FH.root"); // One-step
    Rs_FH[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_FH_TwoStep.root"); // Two-step ML
    Rs_FH[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_FH_Twice_ML.root"); // Twice ML
    for(int j=0;j<3;j++) Rs[j] = Rs_EE[j]->CombineIt(Rs_FH[j]);

    TH1F* h_One = Rs[0]->GetHist(100,0,4*FEnergy[i],DoScale);
    TH1F* h_TwoStep = Rs[1]->GetHist(100,0,4*FEnergy[i],DoScale);
    TH1F* h_Twice = Rs[2]->GetHist(100,0,4*FEnergy[i],DoScale);
    h_Shb = (TH1F*)tempf->Get("TotRecEn");

    GetEnRes(h_One,0,i);
    GetEnRes(h_TwoStep,1,i);
    GetEnRes(h_Twice,2,i);
    GetEnRes(h_Shb,3,i);

    // FitPlots(h_One,"Loos function: #chi^{2}","OneStep",i);
    // FitPlots(h_TwoStep,"two-step training","TwoStep",i);
    // FitPlots(h_Twice,"Twice training","Twice",i);
    // FitPlots(h_Shb,"weighted rechit energy","Shb",i);
  }

  for(int i=0;i<3;i++){
    TCanvas* tc = new TCanvas("aa","bb",800,600);
    TGraph* tg[4];
    TLegend* tl = new TLegend(0.6,0.6,0.9,0.9);
    TMultiGraph* tm = new TMultiGraph();
    tm->SetTitle("HGCAL Pion energy resolution | Sim (FTFP_BERT_EMN);Beam Energy [GeV];"+meth[i]);
    tm->SetMaximum(0.4);
    tm->SetMinimum(0);

    tg[0] = new TGraph(8,FEnergy,ENRES[i][0]);
    tg[1] = new TGraph(8,FEnergy,ENRES[i][1]);
    tg[2] = new TGraph(8,FEnergy,ENRES[i][2]);
    tg[3] = new TGraph(8,FEnergy,ENRES[i][3]);

    tg_Beauty(tg[0],"Sim Pion | FTFP_BERT_EMN","One loss function",2,21); tl->AddEntry(tg[0],"One loss function","pl"); tm->Add(tg[0],"");
    tg_Beauty(tg[1],"Sim Pion | FTFP_BERT_EMN","two-step training",4,21); tl->AddEntry(tg[1],"two-step training","pl"); tm->Add(tg[1],"");
    tg_Beauty(tg[2],"Sim Pion | FTFP_BERT_EMN","Twice training",6,21); tl->AddEntry(tg[2],"Twice training","pl"); tm->Add(tg[2],"");
    tg_Beauty(tg[3],"Sim Pion | FTFP_BERT_EMN","Weighted rechit energy",1,21); tl->AddEntry(tg[3],"Weighted rechit energy","pl"); tm->Add(tg[3],"");

    // /////////// More Info
    // TLatex* tx = new TLatex(250.0,0.225,"Loss function: "+GoodName);
    // tx->SetTextColor(4);
    // tx->SetTextAlign(23);
    // tx->SetTextSize(0.04);
    // /////////// End

    // tc->SetGrid();
    // tm->Draw("apl"); tl->Draw(); // tx->Draw("same");
    // tm->GetYaxis()->SetTitleOffset(1.2);
    // // SaveMe(tc,"EnRes_"+mtitle[i]+"_"+NOW,false);
    // delete tc;

    // Fit_EnRes(tg[0]);
    Fit_EnRes(tg[0],"One loss function",meth[i],mtitle[i]+"_OneStep");
    Fit_EnRes(tg[1],"two-step training",meth[i],mtitle[i]+"TwoStep");
    Fit_EnRes(tg[2],"Twice training",meth[i],mtitle[i]+"TwiceML");
    Fit_EnRes(tg[3],"Weighted rechit energy",meth[i],mtitle[i]+"Shb");
  }
}

void GetEnRes(TH1F* HH, int mlmod, int energy)
{
  float* aa = TruncRMS(HH,0.9);
  float* bb = TruncRMS(HH,1.0);

  ENRES[0][mlmod][energy] = Find_EnRes(HH);
  ENRES[1][mlmod][energy] = aa[1]/aa[0];
  ENRES[2][mlmod][energy] = bb[1]/bb[0];
}

void FitPlots(TH1F* hh, TString modl, TString modl2, int ii)
{
  gStyle->SetOptStat(false);
  // gStyle->SetOptFit(true);
  gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);

  TCanvas* cc = new TCanvas("qq","ww",800,600);
  double qq = Find_EnRes(hh);
  Beauty(hh,"Beam Energy "+SEnergy[ii]+" GeV",modl,4,1,0.55,0.9);
  cc->Update();
  SaveMe(cc,"FitPlots/E"+SEnergy[ii]+"_"+NOW+"_Model_"+modl2,false);
  delete cc;
  gStyle->SetOptStat(true);
  gStyle->SetOptFit(false);
}

void Fit_EnRes(TGraph* gg, TString hname, TString yname, TString fname)
{
  MyFun->SetParNames("a","b","c");
  MyFun->SetLineStyle(7);
  MyFun->SetLineColor(4);
  TCanvas* tc = new TCanvas("aa","ss",800,600);

  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.15);

  gg->Draw("ap");
  gg->GetXaxis()->SetTitle("Pion Beam Energy(E) [GeV]");
  gg->SetTitle("Sim Pion | FTFP_BERT_EMN | "+hname);
  // gg->GetYaxis()->SetTitle("Energy resolution (#frac{#DeltaE}{E})");
  gg->GetYaxis()->SetTitle(yname);
  gg->GetYaxis()->SetTitleOffset(1.3);
  gg->GetYaxis()->SetRangeUser(0.0,0.4);
  gg->SetMarkerColor(1);
  gg->SetLineStyle(9);
  TFitResultPtr rr = gg->Fit("reso","SQ","");

  TLatex* tx = new TLatex(270,0.26,"#frac{#DeltaE}{E} = "+EnRes_TXT);
  tx->SetTextSize(0.04);
  tx->SetTextAlign(21);
  tx->SetTextColor(4);
  tx->Draw();
  tc->Update();

  SaveMe(tc,"EnRes_Fit_"+NOW+"_"+fname,false);
  delete tc;
}

void EnResponse()
{
  TCanvas* tc = new TCanvas("aa","bb",800,600);
  TMultiGraph* tm = new TMultiGraph();
  tm->SetTitle("HGCAL Pion energy response | Sim (FTFP_BERT_EMN);Beam Energy [GeV];#frac{Avg. Reco Energy}{Beam Energy}");
  tm->SetMaximum(2);
  tm->SetMinimum(0);

  TGraph* tg[4];
  TLegend* tl = new TLegend(0.6,0.6,0.9,0.9);

  Double_t OneStep[8],TwoStep[8],Twice[8],Shubham[8];

  for(int i=0;i<8;i++){
    TFile* tempf = new TFile("../Assignment_3_UsingSS/IN_ROOT/Pion_"+InDate+"_Sim_"+SEnergy[i]+"_ForML.root");
    Results* Rs_EE[3];
    Results* Rs_FH[3];
    Results* Rs[3];
    TH1F* h_Shb;

    Rs_EE[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_EE.root"); // One-step
    Rs_EE[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_EE_TwoStep.root"); // Two-step ML
    Rs_EE[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_EE_Twice_ML.root"); // Twice ML

    Rs_FH[0] = new Results("../Assignment_3_UsingSS/OUT_ROOT/CF_Chi2_18_12_19_E"+SEnergy[i]+"_FH.root"); // One-step
    Rs_FH[1] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_FH_TwoStep.root"); // Two-step ML
    Rs_FH[2] = new Results("../Assignment_5_TwoStep_ML/OUT_ROOT/CF_RMS_Chi2_26_12_19_E"+SEnergy[i]+"_temp_FH_Twice_ML.root"); // Twice ML
    for(int j=0;j<3;j++) Rs[j] = Rs_EE[j]->CombineIt(Rs_FH[j]);

    TH1F* h_One = Rs[0]->GetHist(100,0,4*FEnergy[i],false);
    TH1F* h_TwoStep = Rs[1]->GetHist(100,0,4*FEnergy[i],false);
    TH1F* h_Twice = Rs[2]->GetHist(100,0,4*FEnergy[i],false);
    h_Shb = (TH1F*)tempf->Get("TotRecEn");

    OneStep[i] = h_One->GetMean()/DEnergy[i];
    TwoStep[i] = h_TwoStep->GetMean()/DEnergy[i];
    Twice[i] = h_Twice->GetMean()/DEnergy[i];
    Shubham[i] = h_Shb->GetMean()/DEnergy[i];
  }

  tg[0] = new TGraph(8,DEnergy,OneStep);
  tg[1] = new TGraph(8,DEnergy,TwoStep);
  tg[2] = new TGraph(8,DEnergy,Twice);
  tg[3] = new TGraph(8,DEnergy,Shubham);

  tg_Beauty(tg[0],"Sim Pion | FTFP_BERT_EMN","One loss function",2,21); tl->AddEntry(tg[0],"One loss function","pl"); tm->Add(tg[0],"");
  tg_Beauty(tg[1],"Sim Pion | FTFP_BERT_EMN","two-step training",4,21); tl->AddEntry(tg[1],"two-step training","pl"); tm->Add(tg[1],"");
  tg_Beauty(tg[2],"Sim Pion | FTFP_BERT_EMN","Twice training",6,21); tl->AddEntry(tg[2],"Twice training","pl"); tm->Add(tg[2],"");
  tg_Beauty(tg[3],"Sim Pion | FTFP_BERT_EMN","Weighted rechit energy",1,21); tl->AddEntry(tg[3],"Weighted rechit energy","pl"); tm->Add(tg[3],"");

  tc->SetGrid();
  tc->SetTicky(2);
  tm->Draw("apl"); tl->Draw(); // tx->Draw("same");
  tm->GetYaxis()->SetTitleOffset(1.2);
  SaveMe(tc,"Response_"+NOW,false);
  delete tc;
}
