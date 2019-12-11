#define AnalyzeHGCOctTB_cxx

#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>

using namespace std;
// chip 3022,44,3028

int main(int argc, char* argv[])
{

  if (argc < 6) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "configuration" 
	 <<" " << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *config          = argv[4];
  const char *energy = argv[5];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;

  hgcOctTB.EventLoop(data);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data) {

  // void MIP_Rescale_Init(vector<Double_t>* ,vector<Double_t>* ,vector<Double_t>* ,vector<Double_t>* );
  // int MyComp(int,float*);		// Sorting layers into convenient compartments
  int New_MyComp(int,float*);		// Sorting wrt layer number
  // void SS_Nevt(int,float[40],TH1F*);	// To plot normalized hist of SS location
  double alpha,beta,gammaE,gammaF;	// To apply rechit weights
  float ScalingWt = ScaleWt_Init();	// To scale mean reco energy to beam energy
  
  
  Double_t NComp[6]={0};	// # of events

  Weight_Init();		// To initialise weights. See .h for more info

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  // Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  // Long64_t nbytes3 = 0, nb3 = 0;

  Long64_t region_1_classified_events = 0;
  Long64_t region_2_classified_events = 0;
  Long64_t region_3_classified_events = 0;
  Long64_t non_classified_events = 0;

  int decade = 0;
  int ahc_zeroHit = 0;



  // The new lambda
  // Edited first, on 24.08.19
  // Edited last, on 10.09.19

  float lambda[40]={0.11197,
		    0.145689,
		    0.207199,
		    0.240918,
		    0.302427,
		    0.336146,
		    0.397655,
		    0.431374,
		    0.492884,
		    0.526603,
		    0.588112,
		    0.621831,
		    0.68334,
		    0.717059,
		    0.778569,
		    0.812288,
		    0.873797,
		    0.907516,
		    0.969025,
		    1.00274,
		    1.07992,
		    1.11364,
		    1.19081,
		    1.22453,
		    1.29668,
		    1.3304,
		    1.40402,
		    1.64898,
		    1.99845,
		    2.2931,
		    2.58776,
		    2.88242,
		    3.16908,
		    3.57123,
		    3.86589,
		    4.16055,
		    4.45838,
		    4.74838,
		    5.0352,
		    5.13603};

  //in mm
  float ahc_pos[39] = {27.45, 53.65, 79.85, 106.05, 132.25, 
		       158.45, 184.65,210.85,237.05, 263.25,
		       289.45, 315.65, 341.85, 368.05, 394.25,
		       420.45, 446.65, 472.85, 499.05, 525.25,
		       551.45, 577.65, 603.85, 630.05, 656.25, 
		       682.45, 708.65, 734.85, 761.05, 787.25,
		       813.45, 839.65, 865.85, 892.05, 918.25, 
		       944.45, 970.65, 996.85, 1075.45};
  //in cm
  float ahc_front = 169.9;


  bool DEBUG = false;
  Long64_t cut_count[28];
  Long64_t nEvents = 0;
  Long64_t MIP_pions = 0;
  int TOTAL_ACTIVE_LAYER = -1;
  int EE_LAYER = -1;
  int FH_LAYER = -1;
  if(!strcmp(conf_,"alpha") || !strcmp(conf_,"config1")) {
    TOTAL_ACTIVE_LAYER = 40;
    EE_LAYER = 28;
    FH_LAYER = 12;
  }
  else if(!strcmp(conf_,"bravo") || !strcmp(conf_,"config2")){
    TOTAL_ACTIVE_LAYER = 39;
    EE_LAYER = 28;
    FH_LAYER = 11;
  }
  else if(!strcmp(conf_,"charlie") || !strcmp(conf_,"config3")) {
    TOTAL_ACTIVE_LAYER = 20;
    EE_LAYER = 8;
    FH_LAYER = 12;
  }
  else {
    cout<<"ERROR: Unknown configuration!!!!"<<endl;
    return;
  }




  double E_beam = -1.0;




  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;
  for(int i = 0; i< 28; i++){ cut_count[i] = 0;}
  Long64_t jentry = 0;;
  for (jentry=0; jentry<nentries;jentry++) {
  // for (jentry=0; jentry<3;jentry++) {
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;
    
    // ===============read this entry == == == == == == == == == == ==

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) { break; cout<<"Breaking"<<endl;}
        // cout<<"****"<<endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    nb2 = fChain2->GetEntry(jentry); nbytes2 += nb2;
    // nb3 = fChain3->GetEntry(jentry); nbytes3 += nb3;

    if(NRechits == 0) continue;
    // if(NRechits > 50 || NRechits < 20) continue;
    if(ntracks != 1) continue;
    E_beam = beamEnergy;
    
    const int TOTAL_LAYERS = 79;
    Float_t track_x[TOTAL_LAYERS];
    Float_t track_y[TOTAL_LAYERS];

    track_x[0] = impactX_HGCal_layer_1;
    track_y[0] = impactY_HGCal_layer_1;
    track_x[1] = impactX_HGCal_layer_2;
    track_y[1] = impactY_HGCal_layer_2;
    track_x[2] = impactX_HGCal_layer_3;
    track_y[2] = impactY_HGCal_layer_3;
    track_x[3] = impactX_HGCal_layer_4;
    track_y[3] = impactY_HGCal_layer_4;
    track_x[4] = impactX_HGCal_layer_5;
    track_y[4] = impactY_HGCal_layer_5;
    track_x[5] = impactX_HGCal_layer_6;
    track_y[5] = impactY_HGCal_layer_6;
    track_x[6] = impactX_HGCal_layer_7;
    track_y[6] = impactY_HGCal_layer_7;
    track_x[7] = impactX_HGCal_layer_8;
    track_y[7] = impactY_HGCal_layer_8;
    track_x[8] = impactX_HGCal_layer_9;
    track_y[8] = impactY_HGCal_layer_9;
    track_x[9] = impactX_HGCal_layer_10;
    track_y[9] = impactY_HGCal_layer_10;
    track_x[10] = impactX_HGCal_layer_11;
    track_y[10] = impactY_HGCal_layer_11;
    track_x[11] = impactX_HGCal_layer_12;
    track_y[11] = impactY_HGCal_layer_12;
    track_x[12] = impactX_HGCal_layer_13;
    track_y[12] = impactY_HGCal_layer_13;
    track_x[13] = impactX_HGCal_layer_14;
    track_y[13] = impactY_HGCal_layer_14;
    track_x[14] = impactX_HGCal_layer_15;
    track_y[14] = impactY_HGCal_layer_15;
    track_x[15] = impactX_HGCal_layer_16;
    track_y[15] = impactY_HGCal_layer_16;
    track_x[16] = impactX_HGCal_layer_17;
    track_y[16] = impactY_HGCal_layer_17;
    track_x[17] = impactX_HGCal_layer_18;
    track_y[17] = impactY_HGCal_layer_18;
    track_x[18] = impactX_HGCal_layer_19;
    track_y[18] = impactY_HGCal_layer_19;
    track_x[19] = impactX_HGCal_layer_20;
    track_y[19] = impactY_HGCal_layer_20;
    track_x[20] = impactX_HGCal_layer_21;
    track_y[20] = impactY_HGCal_layer_21;
    track_x[21] = impactX_HGCal_layer_22;
    track_y[21] = impactY_HGCal_layer_22;
    track_x[22] = impactX_HGCal_layer_23;
    track_y[22] = impactY_HGCal_layer_23;
    track_x[23] = impactX_HGCal_layer_24;
    track_y[23] = impactY_HGCal_layer_24;
    track_x[24] = impactX_HGCal_layer_25;
    track_y[24] = impactY_HGCal_layer_25;
    track_x[25] = impactX_HGCal_layer_26;
    track_y[25] = impactY_HGCal_layer_26;
    track_x[26] = impactX_HGCal_layer_27;
    track_y[26] = impactY_HGCal_layer_27;
    track_x[27] = impactX_HGCal_layer_28;
    track_y[27] = impactY_HGCal_layer_28;
    track_x[28] = impactX_HGCal_layer_29;
    track_y[28] = impactY_HGCal_layer_29;
    track_x[29] = impactX_HGCal_layer_30;
    track_y[29] = impactY_HGCal_layer_30;
    track_x[30] = impactX_HGCal_layer_31;
    track_y[30] = impactY_HGCal_layer_31;
    track_x[31] = impactX_HGCal_layer_32;
    track_y[31] = impactY_HGCal_layer_32;
    track_x[32] = impactX_HGCal_layer_33;
    track_y[32] = impactY_HGCal_layer_33;
    track_x[33] = impactX_HGCal_layer_34;
    track_y[33] = impactY_HGCal_layer_34;
    track_x[34] = impactX_HGCal_layer_35;
    track_y[34] = impactY_HGCal_layer_35;
    track_x[35] = impactX_HGCal_layer_36;
    track_y[35] = impactY_HGCal_layer_36;
    track_x[36] = impactX_HGCal_layer_37;
    track_y[36] = impactY_HGCal_layer_37;
    track_x[37] = impactX_HGCal_layer_38;
    track_y[37] = impactY_HGCal_layer_38;
    track_x[38] = impactX_HGCal_layer_39;
    track_y[38] = impactY_HGCal_layer_39;
    track_x[39] = impactX_HGCal_layer_40;
    track_y[39] = impactY_HGCal_layer_40;

    for(int j = 0; j < 39; j++) {
      double z = ahc_front + (ahc_pos[j]/10.0) ;
      double x = (z * m_x) + b_x;
      double y = (z * m_y) + b_y;

      track_x[j+40] = x;
      track_y[j+40] = y;
    }


    Double_t rechitEnergySum = 0.0;
    // Double_t rechitEnergySum_mip_cut = 0.0;
    Double_t un_cali = 0.0;
    Long_t Nrechit_layer[40];
    // Long_t NRechits_EE[28];
    // Long_t NRechits_FH[12][7];

    Long_t NRechits_EE[28];
    Long_t NRechits_FH[12];
    Long_t NRechits_AH[39];

    Double_t RechitsEn_EE[28];
    Double_t RechitsEn_FH[12];
    Double_t RechitsEn_AH[39];
    Double_t RechitEn_layer[40];
    

    int module_part_ = -1;
    int module_layer_ = -1;
    int module_position_ = -1;

    double energy_sum_dR2_[40];
    // double cogX_[40];
    // double cogY_[40];
    
    for(int ii=0;ii<40;ii++){
      if(ii<28) {
	NRechits_EE[ii]=0;
	RechitsEn_EE[ii] = 0.0;
	// dR[ii].clear();
	// dr_min_index[ii] = -1.0;
	// dr_min[ii] = 1.e10;
      }
      else{
	NRechits_FH[ii-28]=0;
	RechitsEn_FH[ii-28] = 0.0;
      }
      Nrechit_layer[ii]=0;
      energy_sum_dR2_[ii] = 0.0;
      RechitEn_layer[ii] = 0.0;
      // cogX_[ii] = -1.0;
      // cogY_[ii] = -1.0;
      // dR[ii].clear();
      // dr_min_index[ii] = -1.0;
      // dr_min[ii] = 1.e10;
      if(ii < 39) {
	RechitsEn_AH[ii] = 0.0;
	NRechits_AH[ii] = 0;
      }

    }

    /// FIll Tree
    if(DEBUG) cout<<"DEBUG: Start Analylizing RecHits!!"<<endl;
    if(DEBUG) cout<<"DEBUG: NRechits = "<<NRechits<<endl;

    std::vector<int> temp_moduleID;
    int Nrech_L1[4] = {0,0,0,0};
    for(int i = 0 ; i < NRechits; i++){
      temp_moduleID.clear();
      int temp_layer = rechit_layer->at(i);
      int temp_chip = rechit_chip->at(i);
      int temp_channel = rechit_channel->at(i);
      int en_chan = temp_chip*1000+temp_channel;
      
      double recx = rechit_x->at(i);
      double recy = rechit_y->at(i);
      
      //channel masking
      if(en_chan == 3022 || en_chan == 3028 || en_chan == 44) continue;
      if(temp_layer==1 && temp_chip==0) continue;
      if(temp_layer==1) {
	Nrech_L1[temp_chip]++;
      }
      
      // noise cut
      std::pair<float, float> temp_mod_chip((int)rechit_module->at(i),(int)rechit_chip->at(i));
      float noise_chip = getNoise(temp_mod_chip);
      if(DEBUG) cout<<"Module, layer, chip , noise = "<<rechit_module->at(i)<<", "<<rechit_layer->at(i)<<", "<<rechit_chip->at(i)<<", "<<noise_chip<<endl;
      if(rechit_amplitudeHigh->at(i) < 3*noise_chip)
	continue;
      
      double trackx = track_x[temp_layer-1];
      double tracky = track_y[temp_layer-1];
      if(temp_layer == 2)

      if(!strcmp(data,"data")) {
	trackx = -1*trackx;
	tracky = -1*tracky;
	std::pair<float, float> dxy_al = dxy_alignment(temp_layer);
	float dx_corr = dxy_al.first;
	float dy_corr = dxy_al.second; 
	recx -= dx_corr;
	recy -= dy_corr;	
      }

      Nrechit_layer[temp_layer-1]++;
      RechitEn_layer[temp_layer-1]+=rechit_energy->at(i);
      temp_moduleID = getModuleLocation(rechit_module->at(i));



      if((temp_layer==38 && rechit_module->at(i)!=62) || (temp_layer==39 && rechit_module->at(i)!=54) || (temp_layer==40 && rechit_module->at(i)!=43)){
	cout<< jentry+1 << " = " << rechit_module->at(i)<<endl;
      }


      if(!temp_moduleID.size() || temp_moduleID.size()<3) {
	cout<<"ERROR: Could NOT locate MODULE location for module "<<rechit_module->at(i)<<endl;
	cout<<"\t more info temp_moduleID.size() = "<<temp_moduleID.size()<<endl;
	return;
      }
      module_part_ = temp_moduleID.at(0);
      module_layer_ = temp_moduleID.at(1);
      module_position_ = temp_moduleID.at(2);
      
      double dR = deltaR(recx,recy,trackx,tracky);
      if(dR < 2.0) {
	// cout<<"#######NOTICE ME SENPAI*********HELLLOOOO"<<endl;
	energy_sum_dR2_[temp_layer-1]+=rechit_energy->at(i);
      }
      
      if(module_part_ == 0) {
	// rechitEnergySum_EE+=rechit_energy->at(i);
	RechitsEn_EE[module_layer_-1]+=rechit_energy->at(i);
	NRechits_EE[module_layer_-1]++;
      }
      else if(module_part_ == 1) {
	// rechitEnergySum_FH+=rechit_energy->at(i);
	RechitsEn_FH[module_layer_-1]+=rechit_energy->at(i);
	NRechits_FH[module_layer_-1]++;
	
      }
      else {
	cout<<"ERROR: Unknown Module Part detected!!!!"<<endl;
	return;
	
      }
      rechitEnergySum+=rechit_energy->at(i);
      
    }


    ////////////////////////////////////////////
    //            AHCAL Part                  //
    ////////////////////////////////////////////
    
    // Double_t rechitEnergySumAHCAL = sim_energyAH;

    /////////////////////////////////////////////

    if(DEBUG) cout<<"DEBUG: For shower start"<<endl;
    if(DEBUG && false) {
      for (int ll = 0; ll < 40; ll++) {
	cout<<ll+1<<"\t"<<energy_sum_dR2_[ll]<<endl;
      }
    }
    
    /////////////////// shower finding algo ////////
    bool MIP = true;
    int shower_start_index = -1;
    float shower_lambda_ = -1.0;
    float shower_weight_ = 1.0;
    if(energy_sum_dR2_[0] > 20) {
      // cout<<"Shower started in layer =1"<<endl;
      shower_start_index = 0;
      shower_lambda_ = lambda[0];
      shower_weight_ = lambda[0];
      MIP = false;
    }
    
    else if(energy_sum_dR2_[1] > 20 && energy_sum_dR2_[1] > 2*energy_sum_dR2_[0]) {
      // cout<<"Shower started in layer =2"<<endl;
      shower_start_index = 1;
      shower_lambda_ = lambda[1];
      shower_weight_ = lambda[1]-lambda[0];
      MIP = false;
    }
    
    else {
      
      for(int i = 2; i < 40; i++) {
	if((energy_sum_dR2_[i] > 20) && (energy_sum_dR2_[i] > 2*energy_sum_dR2_[i-1]) && (energy_sum_dR2_[i] > 2*energy_sum_dR2_[i-2])) {

	  // if(i+1==27 || i+1==28) {
	  //   cout<<"Shower started in layer ="<<i+1<<" ,lambda_ = "<<lambda[i]<<endl;
	  // }
	  shower_start_index = i;
	  shower_lambda_ = lambda[i];
	  shower_weight_ = lambda[i]-lambda[i-1];
	  MIP = false;
	  break;
	}
      }
    }
    ///////////////////////////////////////////////
    //       Niramay - Compartmentalising        //
    ///////////////////////////////////////////////
    // int COMP = MyComp(shower_start_index,lambda);
    int COMP = New_MyComp(shower_start_index,lambda);


    ///////////////////////////////////////////////
    //                Collapsed EE               //
    ///////////////////////////////////////////////

    if(DEBUG) cout<<"DEBUG: Shower start index = "<<shower_start_index<<endl;


    ///////////////////////////////////////////////
    //          Niramay - Rechit weights         //
    ///////////////////////////////////////////////
    ifstream iff("WeightsInfo.txt");
    float bEn;
    while(iff>>bEn>>alpha>>beta>>gammaE>>gammaF){if(beamEnergy==bEn) break;}
    ///////////////////////////////////////////////

    Long_t Nrechit_EE = 0;
    Long_t Nrechit_FH = 0;

    Double_t rechitEnergySum_EE = 0.0;
    Double_t rechitEnergySum_FH = 0.0;
    Double_t rechitEnergySum_AH = 0.0;

    Double_t test = 0.0;
    bool zero_rh[40];
    
    for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){
      if(iL < EE_LAYER) {
	Nrechit_EE+=NRechits_EE[iL];
	rechitEnergySum_EE+=RechitsEn_EE[iL];
      }
      else {
	Nrechit_FH+=NRechits_FH[iL-EE_LAYER];
	rechitEnergySum_FH+=RechitsEn_FH[iL-EE_LAYER];
      }
    }
    
    if(DEBUG) cout<<"DEBUG: NRechits in EE = "<<Nrechit_EE<<"\t Rechit_energy in EE = "<<rechitEnergySum_EE<<endl;
    if(DEBUG) cout<<"DEBUG: NRechits in FH = "<<Nrechit_FH<<"\t Rechit_energy in FH = "<<rechitEnergySum_FH<<endl;
   
    double isMuonLike = false;
    bool isRegion1 = false;
    bool isRegion2 = false;
    bool isRegion3 = false;

    if(shower_start_index == -1) {  //MIP LIKE
      isRegion1 = true;
      isRegion2 = false;
      isRegion3 = false;
      region_1_classified_events++;
    }
    else if(shower_start_index > 27) { //H hadrons
      isRegion1 = false;
      isRegion2 = true;
      isRegion3 = false;
      region_2_classified_events++;
    }
    else {  //EH hadron
      isRegion1 = false;
      isRegion2 = false;
      isRegion3 = true;
      region_3_classified_events++;
    }

    /////////////////// Niramay's Plots

    // float LEAK = energyLeakTransverseEE + energyLeakTransverseFH + energyLeakTransverseAH;
    // float EABS = energyLostEE + energyLostFH + energyLostBH;
    // float invis = beamEnergy - (EABS+LEAK+energyLostBeam+energyLostOutside+energyLeakLongitudinal);

    if(shower_start_index!=-1)
      {
    	// h_Invis->Fill(invis);
	// for(int i=0;i<40;i++) Tot_En_MIP[i] = RechitEn_layer[i];
	EE_Sum = rechitEnergySum_EE;
	FH_Sum = rechitEnergySum_FH;
	AH_Sum = ahc_energySum;

	Double_t FULL_EN;

	if(isRegion2) FULL_EN = (FH_Sum+alpha*AH_Sum)/gammaF; // SS in FH
	if(isRegion3) FULL_EN = (EE_Sum+beta*FH_Sum+alpha*beta*AH_Sum)/gammaE; // SS in EE

	TotRecEn->Fill(FULL_EN);               // cout<<FULL_EN<<endl;
	ScaleTotRecEn->Fill(ScalingWt*FULL_EN);// cout<<FULL_EN<<endl;

	for(int i=0;i<79;i++) Tot_En_MIP[i]=0.0;
	for(int i=0;i<28;i++) Tot_En_MIP[i] = RechitsEn_EE[i];
	for(int i=0;i<12;i++) Tot_En_MIP[i+28] = RechitsEn_FH[i];
	for(int i=0;i<39;i++) Tot_En_MIP[i+40] = ahc_energyPerLayer->at(i);

	BEAMEN = beamEnergy;
	For_ML->Fill();
      }

    /////////////////// The END

    if(isRegion2) {
      // double a = -1.0;
      // double b = -1.0;
      // double alpha = -1.0;  // FH & AHCAL relative weight
      // double gamma = -1.0;
      
      // if(beamEnergy == 20) {
      // 	// a = 0.0611934; b = 0.0183451;
      // 	//a = 0.0630539 ;b = 0.0120932;
      // 	a = 0.0604318 ;b = 0.0307894;
      // 	alpha = 0.45; gamma = 13.205;
      // }
      // if(beamEnergy == 50) {
      // 	// a = 0.0608286 ;b = 0.0136046;
      // 	//a = 0.076075 ;b = 0.0116685;
      // 	a = 0.0679195 ;b = 0.0324864;
      // 	alpha = 0.40; gamma = 12.856;
      // }
      // if(beamEnergy == 80) {
      // 	// a = 0.0622612 ;b = 0.0152219;
      // 	//a = 0.0788491 ;b = 0.013591;
      // 	a = 0.0683456 ;b = 0.0320886;
      // 	alpha = 0.40; gamma = 13.1375;
      // }
      // if(beamEnergy == 100) {
      // 	//a = 0.0786438 ;b = 0.0151615;
      // 	a = 0.0677498 ;b = 0.031738;
      // 	alpha = 0.40; gamma = 13.38;
      // }
      // if(beamEnergy == 120) {
      // 	// a = 0.0600868 ;b = 0.018315;
      // 	//a = 0.0794755 ;b = 0.0174122;
      // 	a = 0.068558 ;b = 0.0314515;
      // 	alpha = 0.40; gamma = 13.38;
      // }
      // if(beamEnergy == 200) {
      // 	//a = 0.0802285 ;b = 0.0178579;
      // 	//a = 0.0681734 ;b = 0.031085;
      // 	a = 0.0678221 ;b = 0.0308716;
      // 	alpha = 0.40; gamma = 13.635;
      // }
      // if(beamEnergy == 250) {
      // 	// a = 0.0802285 ;b = 0.0178579;
      // 	// a = 0.0804709 ;b = 0.0182122;
      // 	a = 0.0678221 ;b = 0.0308716;
      // 	alpha = 0.40; gamma = 13.756;
      // }
      // if(beamEnergy == 300) {
      // 	//a = 0.0844364 ;b = 0.0148193;
      // 	a = 0.0703497 ;b = 0.0293021;
      // 	alpha = 0.40; gamma = 13.7266;
      // }
      
      // double FULL_EN = (rechitEnergySum_FH + alpha*rechitEnergySumAHCAL)/gamma;
      // if(FULL_EN<10){cout<<"SS "<<shower_start_index<<"\tEE "<<rechitEnergySum_EE<<"\tFH "<<rechitEnergySum_FH<<"\tAH "<<rechitEnergySumAHCAL<<"\tFULL_EN "<<FULL_EN<<endl;}


      // double beta=4.5;		// 200GeV Pion: beta=4.5

    }
    
    if(DEBUG) cout<<"Point hotel bravo"<<endl;

    //if(isRegion3 && shower_start_index>0) {
    if(isRegion3) {

      // double alpha = -1.0; // EE & FH + AHCAL relative weight
      // double beta = -1.0o; // FH & AHCAL relative weight
      // double gamma = -1.0; // FH & AHCAL relative weight
      
      // if(beamEnergy == 20) {
      // 	beta = 4.8; alpha = 0.45; gamma = 54.15;
      // }
      // if(beamEnergy == 50) {
      // 	beta = 5.3; alpha = 0.40; gamma = 62.92;
      // }
      // if(beamEnergy == 80) {
      // 	beta = 5.0; alpha = 0.40; gamma = 63.737;
      // }
      // if(beamEnergy == 100) {
      // 	beta = 4.9; alpha = 0.40; gamma = 64.41;
      // }
      // if(beamEnergy == 120) {
      // 	beta = 5.0; alpha = 0.40; gamma = 65.812;
      // }
      // if(beamEnergy == 200) {
      // 	beta = 4.5; alpha = 0.40; gamma = 63.3;
      // }
      // if(beamEnergy == 250) {
      // 	beta = 5.6; alpha = 0.40; gamma = 73.76;
      // }
      // if(beamEnergy == 300) {
      // 	beta = 5.5; alpha = 0.40; gamma = 73.567;
      // }

      // double FULL_EN = (rechitEnergySum_EE + beta*rechitEnergySum_FH + beta*alpha*rechitEnergySumAHCAL)/gamma;
      // if(FULL_EN<10){cout<<"SS "<<shower_start_index<<"\tEE "<<rechitEnergySum_EE<<"\tFH "<<rechitEnergySum_FH<<"\tAH "<<rechitEnergySumAHCAL<<"\tFULL_EN "<<FULL_EN<<endl;}
    }

    if(DEBUG) cout<<"Point hotel charlie"<<endl;
    if(DEBUG) cout<<"DEBUG: End of Event = "<<jentry+1<<endl;
    if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
    if(DEBUG && jentry > 5) return;
  } // loop over entries

  cout<<"Got Out "<<jentry<<endl;
  Long64_t total_events = (region_1_classified_events+region_2_classified_events+region_3_classified_events+non_classified_events);
  cout<<"Events with zero hits in AHCAL = "<<ahc_zeroHit<<endl;
  cout<<"MIP like events = "<<((float)region_1_classified_events*100.0)/total_events<<"%"<<endl;
  cout<<"shower start in EE = "<<((float)region_3_classified_events*100.0)/total_events<<"%"<<endl;
  cout<<"shower start in FH = "<<((float)region_2_classified_events*100.0)/total_events<<"%"<<endl;
  cout<<"Non-classified events = "<<((float)non_classified_events*100.0)/total_events<<"%"<<endl;
  //cout<<"Sum = "<<(region_1_classified_events+region_2_classified_events+region_3_classified_events+non_classified_events)<<endl;           
  cout<<"Sum = "<<total_events<<endl;

  if(E_beam < 0) {
    cout<<"E_beam negative!!!"<<endl;
    return;
  }
}


int New_MyComp(int ss , float* lambda)
{
  if(ss<0) return(-1);		// Before EE
  if(ss<6) return(0);           // l<0.35
  if(ss<14) return(1);		// l<0.78
  if(ss<22) return(2);		// l<1.12
  if(ss<28) return(3);		// l<1.80
  if(ss<31) return(4);		// l<2.70
  else return(5);		// FH end
//   if(ss<0) return(-1);
//   if(lambda[ss]<0.35) return(0);
//   if(lambda[ss]<0.75) return(1);
//   if(lambda[ss]<1.03) return(2);
//   if(lambda[ss]<1.43) return(3);
//   if(lambda[ss]<2.7) return(4);
//   else return(5);
}


// void SS_Nevt(int SS,float LMD[40], TH1F* hh)
// {
//   // Input: SS = shower siart index
//   // input: LMD = array of lambda values
//   int New_ss;			// ss after merging back-to back modules

//   if(SS<28){
//     New_ss = (int)(SS/2);
//     // if(New_ss==0) hh->Fill(LMD[0],1/LMD[0]);
//     // else hh->Fill(LMD[2*New_ss+1],1/(LMD[2*New_ss]-LMD[2*New_ss-2]));
//   }
  
//   else{
//     New_ss = SS;
//     // hh->Fill(LMD[New_ss],1/(LMD[New_ss]-LMD[New_ss-1]));
//   }
// }

void AnalyzeHGCOctTB::Weight_Init()
{
  ifstream iif("../For_All_Energies/ScalingWeights.txt");
  int enr,cmp;
  Double_t eabs,eleak,wgt,dat;
  while(iif>>enr>>cmp>>eabs>>eleak>>wgt>>dat){
    if(enr==inEnergy_){
      Wt_Eabs[cmp]=eabs;
      Wt_ELeak[cmp]=eleak;
      Wt_RecEn[cmp]=wgt;

      if(cmp==5) return;
    }
  }
}

float AnalyzeHGCOctTB::ScaleWt_Init()
{
  TString INPATH = "MeanShift_Weights.txt";           // input text file path
  float NumEn[8] = {20,50,80,100,120,200,250,300};
  cout<<"Reading scaling weights from file:"<<endl<<INPATH<<endl;

  ifstream aa(INPATH);
  float temp;
  // float ScaleWt[8];

  for(int i=7;i>=0;i--){
    // aa>>temp>>temp;Mean_All[i] = NumEn[7-i]/temp;
    // aa>>temp;Mean_Sep[i] = NumEn[7-i]/temp;
    aa>>temp>>temp>>temp>>temp;
    // ScaleWt[i]=NumEn[i]/temp;
    if(inEnergy_==(int)NumEn[i]) return(NumEn[i]/temp);
  }
  return(0.1);
}
