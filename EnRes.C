// Program to compare ML with conventional energy reconstruction method
#include"Common_Resources.C"

double Find_EnRes(TH1F* HH);
void Beauty(TH1F*,TString,TString,int,int, Double_t, Double_t);
// void Beauty(TH1F*,TString,TString,int,int, Double_t, Double_t,Tstring DrawOpt);
TString DATE ="27_11_19";

// float EnRes_All[8]={0};		// ML energy resolution - all energies combined
// float EnRes_ML[8]={0};		// ML energy resolution
// float EnRes_Wt[8]={0};		// Weighted rechit energy sum energy resolution

float SigmaMu[8];		// sigma/mu
float RMSMean[8];		// RMS/Mean
float TrEnRes[8];		// truncated RMS/Mean

void GetEnRes(TH1F* HH, int ii);

void EnRes()
{
  for(int i=0;i<8;i++){
    TCanvas* cc = new TCanvas("aq","sw",800,600);
    TFile* tf_ML = new TFile("../OUT_ROOT/Pion_Energy_"+SEnergy[i]+"_GeV_22_11_19.root");
    TH1F* h_ML = (TH1F*)tf_ML->Get("ScaleRecEn_"+SEnergy[i]);h_ML->Rebin();
    Beauty(h_ML,"Sim Pion | Beam Energy "+SEnergy[i]+" GeV", "ML separate energies",1,1,0.6,0.9);
    GetEnRes(h_ML,i);
    delete cc;
  }

  TCanvas* tc = new TCanvas("aa","bb",800,600);
  // TLegend* tl = new TLegend(0.75,0.75,0.9,0.9);
  TLegend* tl = new TLegend(0.45,0.8,0.9,0.9); tl->SetNColumns(3);
  tl->SetTextSize(0.03);
  TMultiGraph* tm = new TMultiGraph();

  TGraph* gr0 = new TGraph(8,FEnergy,SigmaMu); tm->Add(gr0); tl->AddEntry(gr0,"#frac{#sigma}{#mu}","pl");
  TGraph* gr1 = new TGraph(8,FEnergy,RMSMean); tm->Add(gr1); tl->AddEntry(gr1,"#frac{RMS}{Mean}","pl");
  TGraph* gr2 = new TGraph(8,FEnergy,TrEnRes); tm->Add(gr2); tl->AddEntry(gr2,"90% truncated #frac{RMS}{Mean}","pl");

  tg_Beauty(gr0,"","ML",4,22);	// blue
  tg_Beauty(gr1,"","ML",2,22);	// red
  tg_Beauty(gr2,"","Wt",1,21);	// black

  tm->Draw("apl"); tl->Draw();
  tm->SetTitle("Pion Energy resolution | HGCAL Test-Beam | Sim (FTFP_BERT_EMN);Pion Beam energy [GeV];Energy Resolution");
  tm->GetYaxis()->SetTitleOffset(1.2);
  SaveMe(tc,"EnRes_Compare_"+DATE,false);
}

void GetEnRes(TH1F* HH, int ii)
{
  float* temp = TruncRMS(HH,0.9);
  SigmaMu[ii] = (float)Find_EnRes(HH);
  RMSMean[ii] = HH->GetRMS()/HH->GetMean();
  TrEnRes[ii] = temp[1]/temp[0];
}

double Find_EnRes(TH1F* HH)
{
  Double_t hmean = HH->GetMean();
  Double_t hsgm = HH->GetStdDev();
  // Double_t rng = 2;
  Double_t rng = 1;

  // TFitResultPtr rr2 = HH->Fit("gaus","SQNR","",hmean-rng*hsgm,hmean+rng*hsgm);
  TFitResultPtr rr2 = HH->Fit("gaus","SQR","",hmean-rng*hsgm,hmean+rng*hsgm);

  double MEAN=rr2->Parameter(1);
  double SIGMA=rr2->Parameter(2);
  return(SIGMA/MEAN);
}

void Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2)
// void Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2,Tstring DrawOpt)
{
  gStyle->SetOptFit(1111);
  h->SetName(hname);
  h->SetMarkerColor(mc);
  h->SetLineColor(mc);
  h->SetMarkerSize(1);
  h->SetMarkerStyle(ms);
  // h->SetLineWidth(2);
  h->SetLineWidth(1);

  h->SetTitle(htitle);
  h->Draw();
  // h->Draw("sames");
  // h->Draw(DrawOpt);

  gPad->Update();
  TPaveStats* st = (TPaveStats*)h->FindObject("stats");

  st->SetTextColor(mc);
  st->SetY1NDC(Y1);
  st->SetY2NDC(Y2);
}
