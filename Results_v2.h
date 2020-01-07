// Class to store results of ML.
class Results
{
 public:
  Results(TString fname);
  Results();
  float mean=0;
  float emin,emax;		// range of the reconstructed energy
  int nentries;			// no. of envents
  float* EPred;			// array of predicted energies
  Double_t* Etrue;			// array of beam energies
  Results* CombineIt(Results*);		// Operator overloadin
  TH1F* GetHist(int,float,float,bool);

 private:
  TTree* tr;			// Tree for given root file.
  TBranch* b_pred;		// predicted energy
  TBranch* b_tru;		// true energy
  /* Double_t t1; */
  float t1;
  float t2;

  void InitEvents();
};

Results::Results()
{
}

void Results::InitEvents()
{
  tr->SetBranchAddress("E_test_tru",&t1,&b_tru);
  tr->SetBranchAddress("E_test_pred",&t2,&b_pred);

  EPred = new float[nentries];
  Etrue = new Double_t[nentries];

  tr->GetEntry(0);
  emin=t2;
  emax=t2;

  for(int i=0;i<nentries;i++){
    tr->GetEntry(i);
    Etrue[i] = t1;
    EPred[i] = t2;
    if(emin>(float)t1) emin = (float)t1;
    if(emax<t2) emax = t2;
    mean+=t2;
  }
  mean=mean/nentries;
}

TH1F* Results::GetHist(int NBIN,float MIN,float MAX, bool DoScale)
{
  TH1F* HH = new TH1F("","",NBIN,MIN,MAX);

  float wt =1;
  if(DoScale) wt = Etrue[0]/mean;

  for(int i=0;i<nentries;i++){
    HH->Fill(wt*EPred[i]);
  }
  return(HH);
}

Results::Results(TString fname)
{
  TFile* tf = new TFile(fname);
  // tr = (TTree*)tf->Get("MainTree");
  tr = (TTree*)tf->Get("Test");
  nentries = tr->GetEntries();
  InitEvents();
}

Results* Results::CombineIt(Results* aa)
{
  Results* NU = new Results();
  NU->nentries = nentries+aa->nentries;

  NU->EPred = new float[NU->nentries];
  NU->Etrue = new Double_t[NU->nentries];

  for(int i=0;i<nentries;i++){
    NU->EPred[i] = EPred[i];
    NU->Etrue[i] = Etrue[i];
  }

  for(int i=0;i<aa->nentries;i++){
    NU->EPred[i+nentries] = aa->EPred[i];
    NU->Etrue[i+nentries] = aa->Etrue[i];
  }
  NU->emin = emin;
  NU->emax = emax;
  NU->mean = (mean*nentries+aa->mean*aa->nentries)/(nentries+aa->nentries);

  if(emin>aa->emin) NU->emin = aa->emin;
  if(emax<aa->emax) NU->emax = aa->emax;

  return(NU);
}
