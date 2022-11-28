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
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TList.h"

using namespace std;

void coefficient() {

	TCanvas* c = new TCanvas("all","all",2000,4000); // here correlation at each ad is plotted
	c->Divide(2,4,0.0001,0.0001);

	TCanvas* asd = new TCanvas("asd","asd",3000,1000); // here correlation of combined ads is plotted
	asd->Divide(3,1,0.000001,0.01);

	for(int iEH=1; iEH<=3; iEH++) {
		int Nads;
		Nads = 2;
		if(iEH==3) Nads = 4;

		float mean_Teff,sig_Teff,aux_T, aux_R;
		float Teff[Nads][2077];
		float rate[Nads][2077],mean_rate[Nads],sig_mean_rate[Nads];
		float sigrate[Nads][2077];
		bool use[Nads][2077];
		int first_day = 15332;
		double alph[Nads][2];
		int aux4 = 0;
		TGraphErrors* gr[Nads];
		TGraphErrors* all_gr = new TGraphErrors();
		TF1* f1[Nads];
		TF1* f_all = new TF1("f_all","[0] + [1]*x");
		f_all->SetParNames("q","m");
		TH2F* Allcoeff_hist = new TH2F("All_coeff","All_coeff",30,-2,2,30,-2,2); //used for drawing countors, not needed
		TH1F* hist_rate[Nads];
		TH1F* hist_temp[Nads];
		TH2F* coeff_hist[Nads];


		char filename[64];sprintf(filename,"../muon_rate_EH%d_v3.root",iEH);
		TFile* f = new TFile(filename,"READ");

		TH1F* AD_rates[Nads];
		TH1F* eff_time = (TH1F*)f->Get("Effective Time");
		TH1F* hour_rate;
		TCanvas* canvas[Nads];
		char ADname[64],h_rate_name[64];
		for(int iAD=0;iAD<Nads;iAD++) {
			sprintf(ADname,"Muon Daily Rate AD%d EH%d",iAD+1,iEH);
			AD_rates[iAD] = (TH1F*)f->Get(ADname);
			AD_rates[iAD]->SetTitle(ADname);

			// Get Poisson Error --------------------------------------
			sprintf(h_rate_name,"Muon Rate AD%d EH%d",iAD+1,iEH);
			hour_rate = (TH1F*)f->Get(h_rate_name);
			float error = 0.,aux;
			float global_mean = 0;
			float global_tot_time = 0;
			float content,hourcontent,bin_center,total_time = 0.;
			int c_day = first_day+1,day_temp, aux_rate = 0;
			for(int ibin=1;ibin<=hour_rate->GetNbinsX();ibin++) {

				bin_center = hour_rate->GetBinCenter(ibin);
				day_temp = ((int)bin_center-(int)bin_center%24)/24;

				if(day_temp>c_day) {
					if(total_time>=6.*3600.) {AD_rates[iAD]->SetBinError(c_day-first_day,sqrt(error)/aux_rate);}
					else AD_rates[iAD]->SetBinError(c_day-first_day,0);
					c_day++;
					error = 0.;
					total_time = 0.;
					aux_rate = 0;
				}

				hourcontent = hour_rate->GetBinContent(ibin);
				if(hourcontent==0) continue;
				aux = eff_time->GetBinContent(ibin);
				if(aux!=0) error += hourcontent/aux;
				global_mean += hourcontent*aux;
				total_time += aux;
				global_tot_time += aux;
				aux_rate++;

			}
			// ---------------------------------------------------------

			sprintf(filename,"../temperature/Teff_EH%d.dat",iEH); //change here to location of the temperature data
			FILE* fp = fopen(filename,"r");
			int aux_day;

			// Retrieve temperature data - in my case there are 2078 days, so you must change it to the correspoding
			for(int iline=0;iline<2078;iline++) {
				if(iline==0) { // my case the first line was a header with mean effective temp. and its uncerainty
					fscanf(fp,"# %f %f",&mean_Teff,&sig_Teff);
					continue;
				}

				fscanf(fp,"%d %f",&aux_day,&aux_T); // entry of each day is given as: day Teff
				aux_R = AD_rates[iAD]->GetBinContent(iline);

				use[iAD][iline-1] = true;
				if((aux_day>=15342 && aux_day<=15353) && iEH==3 && iAD==2) use[iAD][iline-1] = false; //These were some outliner days found, comment to keep these days

				if(aux_R == 0) use[iAD][iline-1] = false;
				if(aux_day!=AD_rates[iAD]->GetBinLowEdge(iline) || aux_R==0) {
					use[iAD][iline-1] = false;
				}

				Teff[iAD][iline-1] = aux_T;
				rate[iAD][iline-1] = aux_R;
				sigrate[iAD][iline-1] = AD_rates[iAD]->GetBinError(iline);

			}

			int count = 0;
			sig_mean_rate[iAD] = 0;
			mean_rate[iAD] = 0;
			mean_Teff = 0;
			//Two ways of calculating the mean rate:
			//	1) "Classical" of summing all rates and divide buy the total amount of measured days
			//	2) Another is to take the total amount of muons ("global_mean") and divide by the total live time ("global_tot_time")
			//Sadly I dont remember which one we used, I think it was method 2
			for(int j=0;j<2077;j++) {if(use[iAD][j]) { mean_Teff+=Teff[iAD][j]; mean_rate[iAD]+=rate[iAD][j]; count++; sig_mean_rate[iAD]+=pow(sigrate[iAD][j],2);}} //method 1
			mean_rate[iAD] /= count; // method 1
			sig_mean_rate[iAD] = sqrt(sig_mean_rate[iAD])/count; // method 1
			// mean_rate[iAD] = global_mean/global_tot_time; // method 2
			// sig_mean_rate[iAD] = sqrt(global_mean)/global_tot_time; // method 2

			mean_Teff /=  count;

			gr[iAD] = new TGraphErrors();
			char th2_name[64];sprintf(th2_name,"coeff_AD%d",iAD+1);
			coeff_hist[iAD] = new TH2F(th2_name,th2_name,30,-2,2,30,-2,2);
			sprintf(h_rate_name,"rate_AD%d",iAD+1);
			hist_rate[iAD] = new TH1F(h_rate_name,h_rate_name,2077,15332,17409);
			sprintf(h_rate_name,"temp_AD%d",iAD+1);
			hist_temp[iAD] = new TH1F(h_rate_name,h_rate_name,2077,15332,17409);

			double ex = 100*sig_Teff/mean_Teff;
			int aux3 = 0;
			for(int ientry=0;ientry<2077;ientry++) {

				if(!use[iAD][ientry]) continue;

				double x = 100*(Teff[iAD][ientry]-mean_Teff)/mean_Teff;
				double diff_R = rate[iAD][ientry]-mean_rate[iAD];
				double y = 100*diff_R/mean_rate[iAD];
				double dR = sqrt( pow(sigrate[iAD][ientry],2) + pow(sig_mean_rate[iAD],2) );
				double ey = 100*abs(sigrate[iAD][ientry])/mean_rate[iAD];

				hist_rate[iAD]->SetBinContent(ientry,rate[iAD][ientry]);
				hist_temp[iAD]->SetBinContent(ientry,Teff[iAD][ientry]);

				gr[iAD]->TGraphErrors::SetPoint(aux3, x, y);
				gr[iAD]->SetPointError(aux3, ex, ey);
				all_gr->TGraphErrors::SetPoint(aux4, x, y);
				if(iEH==3 && iAD==0) {
					float new_err = sqrt( pow(ey,2) + pow(0.1536*y,2) )/ey; // include syst uncerainty for EH3 AD1 when calculating combined fit - must be recaulculated
					all_gr->SetPointError(aux4, ex, new_err);
				} else {
					all_gr->SetPointError(aux4, ex, ey);
				}
				coeff_hist[iAD]->Fill(x,y);
				Allcoeff_hist->Fill(x,y);
				aux3++;
				aux4++;
			}

			//Plotting options of correlation plots
			char func_name[64];sprintf(func_name,"f%d",iAD+1);
			f1[iAD] = new TF1(func_name,"[0] + x*[1]");
			f1[iAD]->SetParName(0,"q");
			f1[iAD]->SetParName(1,"m");
			if(iEH==3) {
				gr[iAD]->GetXaxis()->SetLimits(-3,2);
				gr[iAD]->GetYaxis()->SetRangeUser(-2.5,2.5);
			} else {
				gr[iAD]->GetXaxis()->SetLimits(-2.5,2);
				gr[iAD]->GetYaxis()->SetRangeUser(-1.5,2);
			}
			gr[iAD]->Fit(func_name,"Q");

			alph[iAD][0] = f1[iAD]->GetParameter(1);
			alph[iAD][1] = f1[iAD]->GetParError(1);


			c->cd(2*iEH+iAD-1);
			gr[iAD]->Draw("A P");
			gPad->Update();

			gStyle->SetOptFit();
			gStyle->SetFitFormat("4.3f");
			gStyle->SetOptTitle(0);
			c->Modified();
			c->Update();
			TPaveStats* st = (TPaveStats*)gr[iAD]->FindObject("stats");
			st->SetX1NDC(0.1);
			st->SetX2NDC(0.51);
			st->SetY1NDC(0.65);
			st->SetY2NDC(0.95);
			st->SetTextSize(0.07);

			TText* text = new TText();
            text->SetNDC();
            text->SetTextSize(0.08);
            text->SetTextAlign(22);
            text->SetTextAngle(0);
            char plot_name[64];sprintf(plot_name,"EH%d AD%d",iEH,iAD+1);
            text->DrawText(0.8, 0.25, plot_name);
            text->Draw();

			c->Modified();
			c->Update();
            gPad->SetTopMargin(0.05);
            gPad->SetBottomMargin(0.15);
			gPad->Update();
			gPad->SetGrid(1,1);
			char plotname[64];sprintf(plotname,"EH%d AD%d",iEH,iAD+1);
			gr[iAD]->SetTitle(plotname);
			gr[iAD]->SetMarkerStyle(kFullCircle);
			gr[iAD]->SetMarkerColor(4);
			gr[iAD]->GetXaxis()->SetTitle("#DeltaT [%]");
			gr[iAD]->GetYaxis()->SetTitle("#DeltaR [%]");
            gr[iAD]->GetXaxis()->SetLabelSize(0.08);
            gr[iAD]->GetYaxis()->SetLabelSize(0.08);
            gr[iAD]->GetXaxis()->SetTitleSize(0.08);
            gr[iAD]->GetYaxis()->SetTitleSize(0.08);
            gr[iAD]->GetYaxis()->SetNdivisions(4,2,1);

            gr[iAD]->GetYaxis()->SetTitleOffset(0.6);
            gr[iAD]->GetXaxis()->SetTitleOffset(1);

		}

		for(int iAD=0;iAD<Nads;iAD++) {
			float syst_err = 15.05/100.; //systematic uncertainty - change for current
			if(iEH==3 && iAD==0) syst_err = 15.35/100.; //systematic uncertainty - change for current
			cout<<"alph AD"<<iAD+1<<" = "<<alph[iAD][0]<<" +/- "<<alph[iAD][1]<<" +/- "<<alph[iAD][0]*syst_err<<" "<< sqrt( pow(alph[iAD][1],2) + pow(alph[iAD][0]*syst_err,2) ) <<endl;
		}

		// Combined results
		all_gr->Fit("f_all","Q");
		cout<<"alph EH"<<iEH<<" = "<<f_all->GetParameter(1)<< " +/- "<<f_all->GetParError(1)<<" +/- "<<f_all->GetParameter(1)*0.1505<<" "<<sqrt( pow(f_all->GetParameter(1)*0.1505,2) + pow(f_all->GetParError(1),2) )<<endl; //systematic uncertainty (0.1505) - change for current
		asd->cd(iEH);
		all_gr->Draw("A P");
		// Allcoeff_hist->Draw("SAME CONTZ");
		gStyle->SetOptFit();
		asd->Modified();
		asd->Update();
		TPaveStats* st = (TPaveStats*)all_gr->FindObject("stats");
		st->SetX1NDC(0.15);
		st->SetX2NDC(0.6);
		st->SetY1NDC(0.75);
		st->SetY2NDC(0.95);
		st->SetTextSize(0.04);
		asd->Modified();
		asd->Update();
		gPad->Update();
		gPad->SetGrid(1,1);
		char plotname[64];sprintf(plotname,"EH%d",iEH);
		all_gr->SetTitle(plotname);
		all_gr->SetMarkerStyle(kFullCircle);
		all_gr->SetMarkerColor(4);

		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		all_gr->GetXaxis()->SetTitle("#DeltaT [%]");
		all_gr->GetYaxis()->SetTitle("#DeltaR [%]");
		all_gr->GetYaxis()->SetTitleOffset(1.2);
		all_gr->GetYaxis()->SetTitleSize(0.06);
		all_gr->GetXaxis()->SetTitleSize(0.06);
		all_gr->GetYaxis()->SetLabelSize(0.06);
		all_gr->GetXaxis()->SetLabelSize(0.06);

		gPad->Update();

		TText* text = new TText();
        text->SetNDC();
        text->SetTextSize(0.06);
        text->SetTextAlign(22);
        text->SetTextAngle(0);
        char plot_name[64];sprintf(plot_name,"EH%d",iEH);
        text->DrawText(0.78, 0.25, plot_name);
        text->Draw();

		gPad->Update();

	}

	asd->SaveAs("../coeff/coeff_EH.png");
	c->SaveAs("../coeff/all_coeff.png");

}
