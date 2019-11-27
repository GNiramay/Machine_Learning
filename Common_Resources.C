// File to store frequently used variables and functions

TString SEnergy[8] = {"300","250","200","120","100","80","50","20"}; // TString energy
Double_t DEnergy[8] = {300,250,200,120,100,80,50,20}; // double_t energy
Float_t FEnergy[8] = {300,250,200,120,100,80,50,20}; // Float_t energy
TString lyr[40]={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40"};
float lambda[40]={0.11197,0.145689,0.207199,0.240918,0.302427,0.336146,0.397655,0.431374,0.492884,0.526603,0.588112,0.621831,0.68334,0.717059,0.778569,0.812288,0.873797,0.907516,0.969025,1.00274,1.07992,1.11364,1.19081,1.22453,1.29668,1.3304,1.40402,1.64898,1.99845,2.2931,2.58776,2.88242,3.16908,3.57123,3.86589,4.16055,4.45838,4.74838,5.0352,5.13603};
TString SComp[6]={"1","2","3","4","5","6"};

void SaveMe(TCanvas* tc, TString tst, bool IsLog); // To save tcanvas in proper locations
void h1_Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms); // to beautify histograms
void tg_Beauty(TGraph* h, TString htitle, TString hname, int mc,int ms); // to beautify graphs
float* TruncRMS(TH1F*,float);	// returns RMS, Mean for truncated histogram

void SaveMe(TCanvas* tc, TString ttl, bool IsLog)
{
  TString aa;
  if(IsLog) {tc->SetLogy();aa = ttl +"_Log";}
  else aa = ttl +"_Nolog";
  tc->SaveAs("../PNG/"+aa+".png");
  tc->SaveAs("../PDF/"+aa+".pdf");
}


void h1_Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms)
// void h1_Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2)
{
  h->SetName(hname);
  h->SetMarkerColor(mc);
  h->SetLineColor(mc);
  h->SetMarkerSize(1);
  h->SetMarkerStyle(ms);

  h->SetTitle(htitle);
  // h->GetXaxis()->SetTitle("Layer location (in #lambda)");
  // h->GetYaxis()->SetTitleOffset(1.25);
  // h->GetYaxis()->SetTitle("Average total rechit energy (in MIPs)");

  // h->Draw();
  // gPad->Update();
  // TPaveStats* st = (TPaveStats*)h->FindObject("stats");

  // st->SetTextColor(mc);
  // st->SetY1NDC(Y1);
  // st->SetY2NDC(Y2);
}

void tg_Beauty(TGraph* h, TString htitle, TString hname, int mc,int ms)
// void h1_Beauty(TH1F* h, TString htitle, TString hname, int mc,int ms, Double_t Y1, Double_t Y2)
{
  h->SetName(hname);
  h->SetMarkerColor(mc);
  h->SetLineColor(mc);
  // h->SetMarkerSize(1);
  h->SetMarkerSize(1.3);
  h->SetMarkerStyle(ms);

  h->SetTitle(htitle);
  // h->GetXaxis()->SetTitle("Layer location (in #lambda)");
  // h->GetYaxis()->SetTitleOffset(1.25);
  // h->GetYaxis()->SetTitle("Average total rechit energy (in MIPs)");

  // h->Draw();
  // gPad->Update();
  // TPaveStats* st = (TPaveStats*)h->FindObject("stats");

  // st->SetTextColor(mc);
  // st->SetY1NDC(Y1);
  // st->SetY2NDC(Y2);
}

float* TruncRMS(TH1F* H,float lvl)
{
  int mx = H->GetMaximumBin();
  float content = H->GetBinContent(mx);
  float maxlim = lvl*H->Integral();// cout<<maxlim<<endl;
  float* mom = new float[2];

  mom[0] = H->GetBinContent(mx)*H->GetBinCenter(mx);
  mom[1] = H->GetBinContent(mx)*pow(H->GetBinCenter(mx),2);

  int tt=1;
  while(content<maxlim){
    content+=H->GetBinContent(mx+tt);
    mom[0]+=H->GetBinContent(mx+tt)*H->GetBinCenter(mx+tt);
    mom[1]+=H->GetBinContent(mx+tt)*pow(H->GetBinCenter(mx+tt),2);

    if(content<maxlim){
      content+=H->GetBinContent(mx-tt);
      mom[0]+=H->GetBinContent(mx-tt)*H->GetBinCenter(mx-tt);
      mom[1]+=H->GetBinContent(mx-tt)*pow(H->GetBinCenter(mx-tt),2);
    }

    tt++;
  }

  mom[0] = mom[0]/content;
  mom[1] = sqrt(mom[1]/content-pow(mom[0],2));
  return(mom);
}
