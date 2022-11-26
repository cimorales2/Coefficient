#include <iostream>
#include <fstream>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"

using namespace std;

void get_muon_rate(int entry_num) {

	TChain* reconT_AdSimple=new TChain("/Event/Rec/AdSimple");
    TChain* calibStatsT=new TChain("/Event/Data/CalibStats");


    //---get files
    int run_num=-1;
    int EH=-1;

    FILE* runfile=fopen("./run_list_good.txt","r");
    int iline=1;
    while(1){
        fscanf(runfile,"%d %d",&run_num,&EH);
        if(feof(runfile)) break;
        if(iline==entry_num) break;
        iline++;
    }
    fclose(runfile);

    ifstream inputfile("/global/project/projectdirs/dayabay/scratch/mkramer/p17b/release/data/v3/pub/paths.physics.good.p17b.v3.sync.txt");
    // ifstream inputfile("/global/homes/m/mkramer/projscratch/p17b/release/data/v1/pub/physics.good.txt");
    std::string line;

    while(1){  //reading the data from file
        std::getline(inputfile,line);
        if(inputfile.eof()) break;
        std::size_t found = line.find("recon.Neutrino.00");
        std::size_t found2 = line.find("recon.NoTag.00");

        if(found>600) found=found2-3;

        std::string st_runnum = line.substr(found+17,found+22);
        const char* tempchar=st_runnum.c_str();
        int temp_runnum=atoi(tempchar);

        if(temp_runnum>run_num) break;
        if(temp_runnum!=run_num) continue;

        const char* file_path=line.c_str();

        reconT_AdSimple->AddFile(file_path);
        calibStatsT->AddFile(file_path);
    }
    //---end

	int Nads_temp;
    if(EH<3) Nads_temp=2;
    else Nads_temp=4;
    const int Nads=Nads_temp;

    int nentries=reconT_AdSimple->GetEntries();
    cout<<"There are "<<nentries<<" entries"<<endl;

    //Chain---------------------------------------------------


    reconT_AdSimple->SetBranchStatus("*",0);
    calibStatsT->SetBranchStatus("*",0);

    calibStatsT->SetBranchStatus("nHit",1);

    reconT_AdSimple->SetBranchStatus("detector",1);
    reconT_AdSimple->SetBranchStatus("energy",1);
    reconT_AdSimple->SetBranchStatus("triggerTimeSec",1);
    reconT_AdSimple->SetBranchStatus("triggerTimeNanoSec",1);
    reconT_AdSimple->SetBranchStatus("x",1);
    reconT_AdSimple->SetBranchStatus("y",1);
    reconT_AdSimple->SetBranchStatus("z",1);

    reconT_AdSimple->SetMakeClass(1);
    calibStatsT->SetMakeClass(1);

    int detector,time_sec,time_nanosec;
    float energy,x,y,z;
    reconT_AdSimple->SetBranchAddress("detector",&detector);
    reconT_AdSimple->SetBranchAddress("energy",&energy);
    reconT_AdSimple->SetBranchAddress("triggerTimeSec",&time_sec);
    reconT_AdSimple->SetBranchAddress("triggerTimeNanoSec",&time_nanosec);
    reconT_AdSimple->SetBranchAddress("x",&x);
    reconT_AdSimple->SetBranchAddress("y",&y);
    reconT_AdSimple->SetBranchAddress("z",&z);

    int nHit;
    calibStatsT->SetBranchAddress("nHit",&nHit);

    //Trees---------------------------------------------------

    TTree* m_energy = new TTree("m_energy","Muon energy");
    m_energy->Branch("muon_energy",&energy,"energy/F");

    TTree* position = new TTree("position","Muon position");
    position->Branch("x",&x,"x/F");
    position->Branch("y",&y,"y/F");
    position->Branch("z",&z,"z/F");

    TTree* time = new TTree("time","Trigger time Muon");
    time->Branch("timesec",&time_sec,"time_sec/I");
    time->Branch("timenanosec",&time_nanosec,"time_nanosec/I");

    TTree* det = new TTree("detector","AD Detectior");
    det->Branch("ADdetector",&detector,"detector/I");

    //variables-----------------------------------------------

    float Ecut; //MeV

    float time_window = 2.; //us

    int last_muon_sec = 0, last_muon_nanosec = 0, last_muon_nHit = 0;
    int current_timesec = 0, current_timenanosec = 0;

    double time_diff;

    bool muon_event = false;

    //Init----------------------------------------------------

    for(int ientry=0; ientry<nentries;ientry++) {
        if(ientry%10000000==0) cout<<"Running "<<ientry<<"/"<<nentries<<endl;

    	reconT_AdSimple->GetEntry(ientry);
        calibStatsT->GetEntry(ientry);

        if(detector==5 || detector==6) {
        	last_muon_nanosec = time_nanosec;
        	last_muon_sec = time_sec;
        	last_muon_nHit = nHit;
        	continue;
        }

        Ecut = 60.;
        if(EH==3 && detector==1) Ecut=100.;
        if(energy>=Ecut) {
        	time_diff = (time_sec-last_muon_sec)*1.e6+(time_nanosec-last_muon_nanosec)*1e-3;
        	if(time_diff<time_window && last_muon_nHit>12) muon_event = true;
        	else {
        		current_timenanosec = time_nanosec;
        		current_timesec = time_sec;
        		for(int jentry=ientry;jentry<nentries;jentry++) {
        			reconT_AdSimple->GetEntry(jentry);
        			calibStatsT->GetEntry(jentry);
        			time_diff = (time_sec-current_timesec)*1.e6 + (time_nanosec-current_timenanosec)*1e-3;
        			if(time_diff>time_window) break;
        			if(detector>4) {
        				if(time_diff<time_window && nHit>12) muon_event = true;
        			}

        		}

        		reconT_AdSimple->GetEntry(ientry);
        		calibStatsT->GetEntry(ientry);

        	}

        if(muon_event) {
        	muon_event = false;
        	det->Fill();
        	m_energy->Fill();
        	time->Fill();
        	position->Fill();
        }

        }

    }

    char out_name[64];
    sprintf(out_name,"/global/cscratch1/sd/cimorale/files/run_%d_muon_EH%d.root",run_num,EH);
    TFile* output = new TFile(out_name,"RECREATE");
    output->cd();
    det->CloneTree(-1,"fast");
    m_energy->CloneTree(-1,"fast");
    time->CloneTree(-1,"fast");
    position->CloneTree(-1,"fast");
    output->Write();
    output->Close();

}
