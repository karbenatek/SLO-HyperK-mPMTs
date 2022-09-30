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
#include <filesystem>
#include "json.hpp"
using json = nlohmann::json;

Double_t FunG0(Double_t *x, Double_t *par){
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
						(1-w)/(sig0*TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-( TMath::Power(xx-Q0,2) / (2* TMath::Power(sig0,2)) ) )
						+ w*StepFun*alpha*TMath::Exp(-alpha*(xx-Q0))  
					) *TMath::Exp(-mu);
	
	for(Int_t n = 1; n <= 100;++n){
		f += TMath::Power(mu,n) * TMath::Exp(-mu)/(TMath::Factorial(n) * sig1 * TMath::Sqrt(TMath::TwoPi() * n)) * TMath::Exp( - TMath::Power( xx - Q0 - w/alpha - n*Q1, 2) / (2*n*TMath::Power(sig1,2))  );
	}	
	return A * f;
}

int GetT0(vector<float> * t, vector<double> * x, double thr = 0.02){
	int i = 0;
	int x_size = x->size() - 1;
	while(t->at(i) < 0) ++i;
	while(i < x_size && x->at(i) < thr) ++i;
	
	return i;
}
double GetQ(vector<float> * t, vector<double> * x, double max, int i_max, float & edge_t, float & cut_t, double cut_ratio = 0.5, double thr = 0.02){
	int i0 = GetT0(t, x, thr);
	edge_t = t->at(i0);
	int i = i0;
	double Q = 0;
	double cut_thr = cut_ratio * max + (1-cut_ratio) * x->at(i0);
	
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

void DefaultCfg(string DataName = ""){
	string fCfg = DataName + "/cfg.json";
	cout<<"Setting cfg.json to default in "<<DataName<<endl;
	fs::copy("cfg.json", fCfg, fs::copy_options::overwrite_existing);
}
void CopyCfg(string from = "", string to = ""){
	from = from + "/cfg.json";
	to = to + "/cfg.json";
	cout<<"Copying "<<from<<" to "<<to << endl;
	fs::copy(from, to, fs::copy_options::overwrite_existing);
}

void DoHists(string DataName = ""){
	
	string fHistsName = DataName + "/OscLecHists.root";
	cout<<"Making histogram at "<<fHistsName<<endl;
	string fDataName = DataName + "/OscLecData.root";
	string fCfg = DataName + "/cfg.json";
	
	TFile * f = TFile::Open(fDataName.data());
	TTree * tr = (TTree*) f->Get("OscLecData");
	
	//read parameters from config file cfg.json
	
	ifstream cfg_file(fCfg);
	json cfg = json::parse(cfg_file);
	int bins = cfg["hist"]["bins"];
	double max_ampl = cfg["hist"]["max_ampl"]; //max wfm amplitude
	double min_ampl = cfg["hist"]["min_ampl"]; //min wfm amplitude
	double cut_ratio = cfg["hist"]["cut_ratio"];
	double thr = cfg["hist"]["thr"];
	cfg_file.close();
	
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
	TFile * fOut = new TFile(fHistsName.data(),"recreate");

	
	TH2F * wfm_Stack = new TH2F("wfm_Stack","Waveform stack",t->size(),t->front(),t->back(),bins,-0.05,max_ampl);
	TH1F * h_Cut = new TH1F( "h_Cut", "Q integration cut time", t->size(), t->front(), t->back());
	TH1D * h_Ampl = new TH1D("h_Ampl","Amplitudes",bins,0,max_ampl);
	TH1D * h_Q = new TH1D("Q","",bins,0,20);
	TH1F * h_t = new TH1F("t","Pulse rising edge",t->size(), t->front(), t->back());
	
	
	double max;
	int i_max;
	double Q;
	float edge_t;
	float cut_t;
	long int nentries = tr->GetEntriesFast();
	int throwed = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		tr->GetEntry(jentry);
		max = GetMax(x, i_max);
		//if (max == 0.715408 || max > 0.2) continue;
		if (max > max_ampl || max < min_ampl) {
			++throwed; 
			continue;
		}
		Q = GetQ(t, x, max, i_max, edge_t, cut_t, cut_ratio, thr);

		h_Q->Fill(Q);
		h_Cut->Fill(cut_t);
		h_Ampl->Fill(max);
		h_t->Fill(edge_t);
		
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
	string fHistsName = DataName + "/OscLecHists.root";
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
	
	string fHistsImg;
	if (DataName.rfind('/') < DataName.length()) fHistsImg = DataName + "/HistImg_" + DataName.substr(DataName.rfind('/')+1) + ".png";
	else fHistsImg = DataName + "/HistImg_" + DataName + ".png";
	c->Print(fHistsImg.data());
	
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
	string fHistsName = DataName + "/OscLecHists.root";
	string fCfg = DataName + "/cfg.json";
	
	TFile * f = TFile::Open(fHistsName.data());
	ifstream cfg_file(fCfg);
	json cfg = json::parse(cfg_file);
	json fit = cfg["fit"];
	
	
	Double_t Q0 = fit["Q0"];
	
	TCanvas *c_Q = new TCanvas("c_Q","Q");
	TH1D * Q = (TH1D*) f->Get("Q");
	Q->SetStats(0);
	//cout<<"Hist integral = "<<Q->GetEntries()<<endl;
	
	//Q->Scale(1/Q->GetEntries());
	//cout<<"Hist integral after normalization = "<<Q->GetEntries()<<endl;
	Q->Draw("HIST");
	
	TF1 *fitfun = new TF1("fitfun",Fun, Q0,Q->GetXaxis()->GetXmax(),8);
	
	fitfun->SetLineColor(1);
	fitfun->SetLineWidth(2);
	
	fitfun->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu","A");
	
	int pNum;
	string pKey;
	double pVal, pValB, pValT;
	for (auto it = fit.begin(); it != fit.end();){
		pKey = it.key();
		pNum = fitfun->GetParNumber(pKey.data());
		pVal = it++.value();
		pValB = it++.value();
		pValT = it++.value();
		fitfun->SetParameter(pNum, pVal);
		if (pValB < pValT) fitfun->SetParLimits(pNum, pValB, pValT);
		cout<< pNum << " - " << pKey << " = " << pVal << " B: " << pValB << " T: " << pValT << endl;
	}

	//~ fitfun->SetParameters(Q0,0.075,4.1,1,0.7,0.5,2,6000);
	//~ fitfun->SetParLimits(0,-3,3);
	//~ fitfun->SetParLimits(1,0,5);
	//~ fitfun->SetParLimits(2,0,10);
	//~ fitfun->SetParLimits(3,0,5);
	//~ fitfun->SetParLimits(4,0,1);
	//~ fitfun->SetParLimits(6,0.5,20);
	
	fitfun->Draw("same");
	
	
	
	Q->Fit("fitfun","R","");//,Q0,Q->GetXaxis()->GetXmax());
	
	double chi = fitfun->GetChisquare()/( (double) fitfun->GetNDF());
	cout<<"Chi^2/ndf = "<<chi<<endl;
	
	string fFitImg;
	string fFitPars;
	if (DataName.rfind('/') < DataName.length()) {
		fFitImg = DataName + "/FitImg_" + DataName.substr(DataName.rfind('/')+1) + ".png";
		fFitPars = DataName + "/FitPars_" + DataName.substr(DataName.rfind('/')+1) + ".txt";
	}
	else {
		fFitImg = DataName + "/FitImg_" + DataName + ".png";
		fFitPars = DataName + "/FitPars_" + DataName + ".txt";
	}
	c_Q->Print(fFitImg.data());
	
	ofstream flFitPars;
	flFitPars.open(fFitPars, std::ofstream::trunc);
	for(int i = 0; i<8; ++i){
		flFitPars<<fitfun->GetParName( i )<<"\t"<<fitfun->GetParameter(i)<<"\t"<<fitfun->GetParError(i)<<endl;
	}
	flFitPars<<"chi/n\t"<<chi<<endl;
	flFitPars.close();
	//c_Q->cd();
	//fitfun->Draw("same");
	
	FitFuncDraw(fitfun);
	
	
}
void FitHist( TH1D * h , json cfg_fit){
	
	Double_t Q0 = cfg_fit["Q0"];
	
	TF1 *fitfun = new TF1("fitfun",Fun, Q0, h->GetXaxis()->GetXmax(),8);
	
	fitfun->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu","A");
	
	int pNum;
	string pKey;
	double pVal, pValB, pValT;
	for (auto it = cfg_fit.begin(); it != cfg_fit.end();){
		pKey = it.key();
		pNum = fitfun->GetParNumber(pKey.data());
		pVal = it++.value();
		pValB = it++.value();
		pValT = it++.value();
		fitfun->SetParameter(pNum, pVal);
		if (pValB < pValT) fitfun->SetParLimits(pNum, pValB, pValT);
		//~ cout<< pNum << " - " << pKey << " = " << pVal << " B: " << pValB << " T: " << pValT << endl;
	}
	//~ fitfun->SetParameter(7,1);
	cout<<endl<<endl<<"Integral = "<<fitfun->Integral(-1,5)<<endl;
	//~ fitfun->Draw();
	h->Fit("fitfun","R","");
	h->Draw();
}
void FitHistTB(){
	ifstream cfg_file("cfg.json");
	json cfg = json::parse(cfg_file);
	json fit = cfg["fit"];
	
	TFile * f = TFile::Open("OscLecHists.root");
	TH1D * h = (TH1D*) f->Get("Q");
	h->ResetStats();
	FitHist(h, fit);

	
	
}

void DoAPA(string DataName)
{
	string fDataName = DataName + "/OscLecData.root";
	
	TFile * f = TFile::Open(fDataName.data());
	TTree * tr = (TTree*) f->Get("OscLecData");

	int i_Wfm;
	vector<double>  * x = 0;
	vector<float>   * t = 0;
	tr->SetBranchAddress("x",&x);
	tr->SetBranchAddress("t",&t);
	tr->GetEntry(0);
	double thr = 0.02;
	TH1F * T = new TH1F("T", "Pulse time",t->size(), t->front(), t->back());
	
	long int nentries = tr->GetEntriesFast();
	double px = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		tr->GetEntry(jentry);
		for (int i = 0; i < x->size(); ++i){//wfms iterator
				if (x->at(i) >= thr) {
					T->Fill(t->at(i));
					while (x->at(i) >= thr) {
						++i;
					}
				}
			px = x->at(i);
		}
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
	}
}
void DoTSA(string DataName, int merge = 1)//time stability analysis
{
	string fDataName = DataName + "/OscLecData.root";
	string fTSAName = DataName + "/TSA.root";
	string fCfg = DataName + "/cfg.json";
	
	TFile * f = TFile::Open(fDataName.data());
	TTree * tr = (TTree*) f->Get("OscLecData");
	TFile * fOut = new TFile(fTSAName.data(),"recreate");

	ifstream cfg_file(fCfg);
	json cfg = json::parse(cfg_file);
	int bins = cfg["TSA"]["bins"];
	double max_ampl = cfg["hist"]["max_ampl"]; //max wfm amplitude
	double min_ampl = cfg["hist"]["min_ampl"]; //min wfm amplitude
	double cut_ratio = cfg["hist"]["cut_ratio"];
	cfg_file.close();
	
	int i_Wfm;
	int i_Seg;
	TBranch * b_i_Seg;
	char date[22];
	vector<double>  * x = 0;
	vector<float>   * t = 0;
	tr->SetBranchAddress("i_Wfm",&i_Wfm);
	tr->SetBranchAddress("i_Seg",&i_Seg, &b_i_Seg);
	tr->SetBranchAddress("Date",&date);
	tr->SetBranchAddress("x",&x);
	tr->SetBranchAddress("t",&t);
	
	double dRMS;
	
	int i_Pck = 0;
	TH1D * hQ = new TH1D("hQ", "Charges", bins, 0, 20);
	TH1D * pQ = new TH1D("pQ", "Prev Charges", bins, 0, 20);
	TH1D * dQ = new TH1D("dQ", "Charge difference", bins, 0, 20);
	
	
	
	TTree * fChain = new TTree("TSA","Histograms for TSA - Time Stabilisation Analysis");
	fChain->Branch("i_Pck", &i_Pck);
	fChain->Branch("Date", date, "Date/C", 22);
	fChain->Branch("dRMS", &dRMS);
	fChain->Branch("hQ", &hQ);
	fChain->Branch("dQ", &dQ);
	
	double max;
	int i_max;
	double Q;
	float edge_t;
	float cut_t;
	long int nentries = tr->GetEntriesFast();
	tr->GetEntry(0);
	int dnQ;
	int next_i_Seg;
	int merged = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		tr->GetEntry(jentry);
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
		max = GetMax(x, i_max);
		Q = GetQ(t, x, max, i_max, edge_t, cut_t, cut_ratio);
		//~ cout<<"Max: "<<max<<" ,Q: "<<Q<<endl;
		hQ->Fill(Q);
		b_i_Seg->GetEntry(jentry+1);
		next_i_Seg = i_Seg;
		b_i_Seg->GetEntry(jentry);
		if ( next_i_Seg < i_Seg) ++merged;
		if ( merged == merge || jentry == nentries - 1) //for end of data segment or eof 
		{
			if (merged < merge) break;
			merged = 0;
			dRMS = 0;
			if (i_Pck++ != 0)
			{
				for (int bin = 0; bin < bins; ++bin){
					dnQ = hQ->GetBinContent(bin) - pQ->GetBinContent(bin);
					dQ->SetBinContent(bin, dnQ);
					dQ->ResetStats();
					hQ->ResetStats();
					dRMS += dnQ*dnQ;
				}
				dRMS = TMath::Sqrt(dRMS)/bins;
			}
			
			fChain->Fill();
			//~ hQ->Draw();
			for (int bin = 0; bin < bins; ++bin){
				pQ->SetBinContent(bin, hQ->GetBinContent(bin));
				hQ->SetBinContent(bin,0);
			}
			
			//~ return;
		}
	
	}
	cout<<nentries<<"/"<<nentries<<endl<<"TSA Completed\n";
	
	fOut->Write();
	fOut->Close();


}
void PlotTSA(string DataName){
	string fTSAName = DataName + "/TSA.root";
	
	TFile * f = TFile::Open(fTSAName.data());
	TTree * t = (TTree*) f->Get("TSA");
	int nentries = t->GetEntries();

	TH1D *hQ = new TH1D();
	t->SetBranchAddress("hQ",&hQ);
	TH1D *hQs[nentries];

	TH1D *dQ = new TH1D();
	t->SetBranchAddress("dQ",&dQ);
	TH1D *dQs[nentries];
	
	double dRMS;
	t->SetBranchAddress("dRMS",&dRMS);
	TGraph * GdRMS = new TGraph();
	
	int h = 600;
	TCanvas *c = new TCanvas("TSA","Time Stability Analysis");
	
	
	c->Divide(2,nentries + 1);
	c->SetCanvasSize(1800,h*nentries);
	c->SetWindowSize(1800,1000);
	
	THStack *hs = new THStack("hs","");
	
	for (int jentry = 0; jentry < nentries; ++jentry){
		t->GetEntry(jentry);
		hQs[jentry] = (TH1D*) hQ->Clone( ("hQ"+ to_string(jentry+1)).data() );
		dQs[jentry] = (TH1D*) dQ->Clone( ("dQ"+ to_string(jentry+1)).data() );
		//~ hQs[jentry]->SetLineColor(800+jentry);
		hs->Add(hQs[jentry]);
		GdRMS->AddPoint(jentry,dRMS);
	}
	c->cd(1);
	hs->Draw("NOSTACK");
	c->cd(2);
	GdRMS->Draw();	
	for (int jentry = 0; jentry < nentries; ++jentry){
		c->cd( 2*jentry + 3);
		hQs[jentry]->Draw();
		c->cd( 2*jentry + 4);
		dQs[jentry]->Draw();
	}
	
	
	
}
#endif
