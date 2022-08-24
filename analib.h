#ifndef  analib_h
#define analib_h
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

Double_t FunG0(Double_t *x, Double_t *par){
//fdfdsf
	Float_t xx = x[0];
	Double_t Q0 = par[0];
	Double_t sig0 = par[1];
	Double_t w = par[4];
	Double_t alpha = par[5];
	Double_t mu = par[6];
	Double_t A = par[7];
	Double_t f = ((1-w)/(Q0*TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-( TMath::Power(xx-Q0,2) / (2* TMath::Power(sig0,2)) ) )) *TMath::Exp(-mu);
	return A*f;
}
Double_t FunNoise(Double_t *x, Double_t *par){
	Float_t xx = x[0];
	Double_t Q0 = par[0];
	Double_t sig0 = par[1];
	Double_t w = par[4];
	Double_t alpha = par[5];
	Double_t mu = par[6];
	Double_t A = par[7];
	Double_t StepFun;
	if ((xx-Q0) >= 0) StepFun = 1;
	else StepFun = 0;
	Double_t f = ( w*StepFun*alpha*TMath::Exp(-alpha*(xx-Q0))) *TMath::Exp(-mu);
	return A*f;
}

Double_t FunGn(Double_t *x, Double_t *par){
	Float_t xx = x[0];
	Double_t Q0 = par[0];
	Double_t Q1 = par[2];
	Double_t sig1 = par[3];
	Double_t w = par[4];
	Double_t alpha = par[5];
	Double_t mu = par[6];
	Double_t A = par[7];

	Double_t f = 0;
	for(Int_t n = 1; n <= 100;++n){
		f += TMath::Power(mu,n) * TMath::Exp(-mu)/(TMath::Factorial(n) * sig1 * TMath::Sqrt(TMath::TwoPi() * n)) * TMath::Exp( - TMath::Power( xx - Q0 - w/alpha - n*Q1, 2) / (2*n*TMath::Power(sig1,2))  );
	}
	return A*f;

}
Double_t Fun(Double_t *x, Double_t *par){
	Float_t xx = x[0];
	Double_t Q0 = par[0];
	Double_t sig0 = par[1];
	Double_t Q1 = par[2];
	Double_t sig1 = par[3];
	Double_t w = par[4];
	Double_t alpha = par[5];
	Double_t mu = par[6];
	Double_t A = par[7];
	
	
	Double_t StepFun;
	if ((xx-Q0) >= 0) StepFun = 1;
	else StepFun = 0;
	Double_t f = 	(  
						(1-w)/(Q0*TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-( TMath::Power(xx-Q0,2) / (2* TMath::Power(sig0,2)) ) )
						+ w*StepFun*alpha*TMath::Exp(-alpha*(xx-Q0))  
					) *TMath::Exp(-mu);
	
	for(Int_t n = 1; n <= 100;++n){
		f += TMath::Power(mu,n) * TMath::Exp(-mu)/(TMath::Factorial(n) * sig1 * TMath::Sqrt(TMath::TwoPi() * n)) * TMath::Exp( - TMath::Power( xx - Q0 - w/alpha - n*Q1, 2) / (2*n*TMath::Power(sig1,2))  );
	}	
	return A * f;
}

int GetT0(vector<float> * t){
	int i = 0;
	
	while(t->at(i) <= 0)++i;
	return i;
}
double GetQ(vector<float> * t, vector<double> * x, double max, int i_max, float & cut_t, double cut_ratio = 0.5){
	int i0 = GetT0(t);
	int i = i0;
	double Q = 0;
	double cut_thr = cut_ratio * max + (1-cut_ratio) * x->at(i0);
	double thr = x->at(i0);
	while ( i < x->size()-1 ) {
		Q += x->at(i++);
		if( x->at(i) <= cut_thr && i > i_max) break;
	}
	cut_t = t->at(i);
	return Q;	
}
	
double GetMax(vector<double> * x, int & i_max){
	double max = x->at(0);
	for(int i = 1; i < x->size();++i){
		if ( x->at(i) > max) {
			max = x->at(i);
			i_max = i;
		}
	
	}
	return max;
}
void DoHists(string DataName = ""){
	string fHistsName = "OscLecHists_"+DataName+".root";
	string fDataName = "OscLecData_"+DataName+".root";
	TFile * f = TFile::Open(fDataName.data());
	TTree * tr = (TTree*) f->Get("OscLecData");
	
	int i_Wfm;
	vector<double>  * x = 0;
	vector<float>   * t = 0;
	tr->SetBranchAddress("i_Wfm",&i_Wfm);
	tr->SetBranchAddress("x",&x);
	tr->SetBranchAddress("t",&t);
	
	tr->GetEntry(0);
	
	int n_Seg = x->size();
	//TCanvas *c = new TCanvas("c","Waveform");
	//TGraph *gr = new TGraph(n_Seg);
	//gr->SetLineColor(4);
	//gr->SetLineWidth(1);
	double top = 0.6; //max wfm amplitude
	int bins = 200;
	TFile * fOut = new TFile(fHistsName.data(),"recreate");

	
	TH2F * wfm_Stack = new TH2F("wfm_Stack","Waveform stack",t->size(),t->front(),t->back(),bins,0,top);
	TH1F * h_Cut = new TH1F( "h_Cut", "Q integration cut time", t->size(), t->front(), t->back());
	TH1D * h_Ampl = new TH1D("h_Ampl","Amplitudes",bins,0,top);
	TH1D * h_Q = new TH1D("Q","",bins,0,20);
	
	
	double max;
	int i_max;
	double Q;
	float cut_t;
	long int nentries = tr->GetEntriesFast();
	int throwed = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		tr->GetEntry(jentry);
		max = GetMax(x, i_max);
		//if (max == 0.715408 || max > 0.2) continue;
		if (max > top) {
			++throwed;
			continue;
		}
		Q = GetQ(t, x, max, i_max, cut_t,0.3);

		h_Q->Fill(Q);
		h_Cut->Fill(cut_t);
		h_Ampl->Fill(max);
		
		for (int i = 0; i < x->size(); ++i){//wfms iterator
			//gr->SetPoint(i,(double)t->at(i),x->at(i));
			wfm_Stack->Fill(t->at(i) , (float) x->at(i) );
		}
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
		
	}
	cout<<"From "<<nentries<<" wfm's, "<<throwed<<" was deleted"<<endl;
	//plot
	
	
	fOut->Write();
	
	f->Close();
	
	
}
void PlotHists(string DataName = ""){
	string fHistsName = "OscLecHists_"+DataName+".root";
	TFile * f = TFile::Open(fHistsName.data());
	TCanvas *c = new TCanvas("c","",1600,600);
	c->Divide(2,1);
	c->cd(1);
	
	TH1D * h_Ampl = (TH1D*) f->Get("h_Ampl");
	h_Ampl->SetTitle("");
	h_Ampl->GetXaxis()->SetTitle("Amplitude [V]");
	h_Ampl->GetYaxis()->SetTitle("Cnts");
	
	h_Ampl->SetStats(0);
	h_Ampl->Draw();
	
	c->cd(2);
	TH1D * h_Q = (TH1D*) f->Get("Q");
	h_Q->GetXaxis()->SetTitle("Charge");
	h_Q->GetYaxis()->SetTitle("Cnts");
	h_Q->SetStats(0);
	h_Q->Draw();
	c->Print("Hists.png");
	
	TCanvas *c1 = new TCanvas("c1","",800,1200);
	c1->Divide(1,2);
	c1->cd(1);
	TH2F * wfm_Stack = (TH2F*) f->Get("wfm_Stack");
	wfm_Stack->Draw("colz");
	gPad->SetLogz();
	c1->cd(2);
	TH1D *h_Cut = (TH1D*) f->Get("h_Cut");
	h_Cut->Draw();
}

void FitFuncDraw(TF1 *f){
	Int_t n = 100;
	TCanvas *c2 = new TCanvas("c2", "Fit function parts");
	TF1 * G0 = new TF1("G0",FunG0, f->GetXmin(), f->GetXmax(),8);
	TF1 * Noise = new TF1("Noise",FunNoise, f->GetXmin(), f->GetXmax(),8);
	TF1 * Gn = new TF1("Gn",FunGn, f->GetXmin(), f->GetXmax(),8);
	
	f->SetLineColor(1);
	G0->SetLineColor(3);
	Noise->SetLineColor(6);
	Gn->SetLineColor(2);
	G0->SetLineStyle(2);
	Noise->SetLineStyle(2);
	Gn->SetLineStyle(2);
	G0->SetParameters(f->GetParameters());
	Noise->SetParameters(f->GetParameters());
	Gn->SetParameters(f->GetParameters());
	
	f->Draw();
	
	G0->Draw("same");
	Noise->Draw("same");
	Gn->Draw("same");
	

}

void DoFit(string DataName = ""){
	string fHistsName = "OscLecHists_"+DataName+".root";
	TFile * f = TFile::Open(fHistsName.data());
	Double_t Q0 = 0.075;
	TCanvas *c_Q = new TCanvas("c_Q","Q");
	TH1D * Q = (TH1D*) f->Get("Q");
	Q->SetStats(0);
	//cout<<"Hist integral = "<<Q->GetEntries()<<endl;
	
	//Q->Scale(1/Q->GetEntries());
	//cout<<"Hist integral after normalization = "<<Q->GetEntries()<<endl;
	Q->Draw("HIST");
	
	TF1 *fitfun = new TF1("fitfun",Fun, Q0,20,8);
	fitfun->SetLineColor(1);
	fitfun->SetLineWidth(2);
	fitfun->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu","A");
	//fitfun->SetParameters(12,0.19,16.2,6,0.339,0.045,2.5,10); //12,0.19,16.2,6,0.339,0.045,2.5,10
	//fitfun->SetParameters(12,1.9,16.2,5,0.339,0.045,2.2,1);
	
	fitfun->SetParameters(Q0,0.075,1.5,0.83,0.46,0.46,2,4800);
	Q->Fit("fitfun");
	c_Q->cd();
	fitfun->Draw("same");
	
	FitFuncDraw(fitfun);
	
	
}
void FitFunc()
{
	TF1 *f = new TF1("f",Fun, 0,400,8);
	f->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu","A");
	f->SetParameters(10,0.2,35,10,0.5,0.5,1.5,1);
	cout<<"Integral: "<<f->Integral(0,300,1.e-6)<<endl;
	
	FitFuncDraw(f);

}
#endif
