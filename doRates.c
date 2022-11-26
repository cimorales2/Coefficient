#include <iostream>
#include <string>
#include <fstream>
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
#include "TLegend.h"
#include "TSystem.h"

using namespace std;

void do_rates(int nHall) { //name of the files fmt = run_%d_muon_EH%d.root
	int sec, nanosec, nentries, run_num;
	int det=0;
	float engy;

	const int init_run_hour = 367982, end_run_hour = 417793;
	int nro_bins = end_run_hour - init_run_hour;

	int Nads = 2;
	if(nHall==3) Nads = 4;

	char aux_name[64];
	char aux_hname[64];
	TH1D* m_rates[Nads];
	for(int i=0;i<Nads;i++) {
		sprintf(aux_name,"Muon Rate AD%d EH%d",i+1,nHall);
		sprintf(aux_hname,"m_rate_AD%d_EH%d",i+1,nHall);
		m_rates[i] = new TH1D(aux_name,aux_hname,nro_bins,init_run_hour,end_run_hour);
	}

	TH1D* h_eftv_time = new TH1D("Effective Time","effective_time",nro_bins,init_run_hour,end_run_hour);
	TH1F* energy = new TH1F("Energy","Energy",1000,0,1000);

	TFile* f_temp;
	TTree* det_temp;
	TTree* time_temp;
	TTree* en_temp;
	TTree* pos;
	int ibin;
	double current_time, init_hour, time_diff;
	double window = init_run_hour + 1.0;
	double last_time=-1., eftv_time = 0.;
	float time_limit = 6./3600.;

	char name_file[64], fileName[64];
	sprintf(name_file,"./rates_stuff/EH%d_runs.dat",nHall);
	FILE* runfile = fopen(name_file,"r");
	while (1) {
		fscanf(runfile,"%d",&run_num);
		if(feof(runfile)) break;

		sprintf(fileName,"../../beda/for_cristobal/run_%d_muon_EH%d.root",run_num,nHall);
		cout<<fileName<<endl;

		f_temp = new TFile(fileName,"READ");

		det_temp = (TTree*)f_temp->Get("detector");
		time_temp = (TTree*)f_temp->Get("time");
		en_temp = (TTree*)f_temp->Get("m_energy");
		pos = (TTree*)f_temp->Get("position");

		det_temp->SetBranchAddress("ADdetector",&det);
		time_temp->SetBranchAddress("timesec",&sec);
		time_temp->SetBranchAddress("timenanosec",&nanosec);
		en_temp->SetBranchAddress("muon_energy",&engy);

		nentries = det_temp->GetEntries();
		for(int i=0;i<nentries;i++) {
			det_temp->GetEntry(i);
			time_temp->GetEntry(i);
			en_temp->GetEntry(i);
			pos->GetEntry(i);

			if(run_num == 67633 || run_num == 67749 || run_num == 67755 || run_num == 67768) {
				if(det==1) continue;
			}

			if(engy<200 && nHall==3 && det==1) continue;
			if(engy<200 && (nHall!=3 || (nHall==3 && det>1)) )  continue;
			current_time = sec/3600.;
			m_rates[det-1]->Fill(current_time); //have to be in hours
			if(det==1) energy->Fill(engy);

			if(last_time==-1.) init_hour = current_time;

			if(current_time>window) {
				eftv_time += (last_time - init_hour);
				//cout<<eftv_time<<" "<<ibin<<" "<<current_time<<" "<<window<<endl;
				h_eftv_time->SetBinContent(ibin,eftv_time*3600);
				init_hour = current_time;
				window = (sec - sec%3600)/3600. + 1.;
				ibin = int(current_time - init_run_hour)+1;
				eftv_time = 0.;
				last_time = current_time;
				continue;
			}

			//in between check!!
			time_diff = current_time - last_time;
			if(time_diff>=time_limit && last_time!=-1.) {eftv_time -= time_diff;}

			last_time = current_time; // last thing to change
		}

		det_temp->Delete();
		time_temp->Delete();
		en_temp->Delete();
		f_temp->Close();
		f_temp->Delete();
	
	}


	// for(int i=0;i<Nads;i++) {m_rates[i]->Sumw2();}//m_rates[i]->Divide(h_eftv_time);}

	TH1F* daily_bins[Nads]; //From here is just to get daily bins
	double bin_center, content;
	float total_time,error;
	int first_day = 15332, last_day = 17409, day_bins = last_day - first_day;
	int day_temp, c_day;
	for(int i=0;i<Nads;i++) {
		sprintf(aux_name,"Muon Daily Rate AD%d EH%d",i+1,nHall);
		sprintf(aux_hname,"m_daily_rate_AD%d_EH%d",i+1,nHall);
		daily_bins[i]= new TH1F(aux_name,aux_hname,day_bins,first_day,last_day);

		content = 0.;
		total_time = 0.;
		c_day = first_day+1;
		error = 0;

		for(int jbin=1;jbin<=nro_bins;jbin++) {

			bin_center = m_rates[i]->GetBinCenter(jbin);
			day_temp = ((int)bin_center-(int)bin_center%24)/24;

			if(day_temp>c_day) {
				content /= total_time;
				if(total_time>=6.*3600.) {daily_bins[i]->SetBinContent(c_day-first_day,content);daily_bins[i]->SetBinError(c_day-first_day,sqrt(error)/total_time);}
				c_day++;
				total_time = 0.;
				content = 0;
				error = 0.;
			}

			if(m_rates[i]->GetBinContent(jbin)==0.) continue;
			content += m_rates[i]->GetBinContent(jbin); //this is the number of events
			total_time += h_eftv_time->GetBinContent(jbin);
			error += m_rates[i]->GetBinContent(jbin); //sqrt(N)
			// error += pow(m_rates[i]->GetBinError(jbin),2); //sqrt(N)

			if(jbin==nro_bins && total_time>=6.*3600.) {
				content /= total_time;
				daily_bins[i]->SetBinContent(c_day-first_day,content);
				daily_bins[i]->SetBinError(c_day-first_day,sqrt(error)/total_time);
			}

		}

	}
	// for(int i=0;i<Nads;i++) {m_rates[i]->Divide(h_eftv_time);}

	char outname[64];
	sprintf(outname,"../muon_rate_EH%d_v2.root",nHall); 
	TFile* outfile = new TFile(outname,"RECREATE");

	outfile->cd();
	for(int i=0;i<Nads;i++) {m_rates[i]->Write();daily_bins[i]->Write();}
	energy->Write();
	h_eftv_time->Write();
	outfile->Close();
	
}

