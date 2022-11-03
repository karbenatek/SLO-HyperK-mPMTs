#ifndef  analib_h
#define analib_h
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <iostream>
#include <TGraph.h>
#include <vector>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include <filesystem>
#include <limits>
#include "LambertW/LambertW.cc"

#include "json.hpp"
using json = nlohmann::json;


Double_t Funoise(Double_t *x, Double_t *par){
	Float_t xx = x[0];
	Double_t Q0 = par[0];
	Double_t sig0 = par[1];
	Double_t w = par[2];
	Double_t alpha = par[3];
	Double_t A = par[4];	
	Double_t StepFun;
	if ((xx-Q0) >= 0) StepFun = 1;
	else StepFun = 0;
	Double_t f = 	(  
						(1-w)/(sig0*TMath::Sqrt(TMath::TwoPi())) * TMath::Exp(-( TMath::Power(xx-Q0,2) / (2* TMath::Power(sig0,2)) ) )
						+ w*StepFun*alpha*TMath::Exp(-alpha*(xx-Q0))  
					) *A;
	return f;
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
	return f;
}

int GetT0(vector<float> * t, vector<double> * x, double thr = 0.02){
	int i = 0;
	int x_size = x->size() - 1;
	while(t->at(i) < 0) ++i;
	while(i < x_size && x->at(i) < thr) ++i;
	
	return i;
}
double GetQ(vector<float> * t, vector<double> * x, double max, int i_max, float & t0, float & t1, double cut_ratio = 0.5, double thr = 0.02){
	
	
	int i = i_max;
	double Q = 0;
	double cut_thr = cut_ratio * max + (1-cut_ratio) * x->at(GetT0(t,x,thr));
	double X;
	while ( i++ < x->size()-1 ) {
		X = x->at(i);
		Q += X;
		if( x->at(i) <= cut_thr ) break;
	}
	t1 = t->at(--i);

	i = i_max - 1;
	while ( i-- > 0 ) {
		X = x->at(i);
		Q += X;
		if( x->at(i) <= cut_thr ) break;
	}
	t0 = t->at(++i);
	
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
	double Qtop = cfg["hist"]["Qtop"];
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
	TH1F * h_t0 = new TH1F( "h_t0", "t0", t->size(), t->front(), t->back());
	TH1F * h_t1 = new TH1F( "h_t1", "t1", t->size(), t->front(), t->back());
	TH1D * h_Ampl = new TH1D("h_Ampl","Amplitudes",bins,0,max_ampl);
	TH1D * h_Q = new TH1D("Q","",bins,0,Qtop);
	TH2D * h2_AQ = new TH2D("h2_AQ","Voltage/sum",bins, 0, Qtop, bins,0,max_ampl);
	TProfile * pf_AQ = new TProfile("pf_AQ","Voltage/sum",bins, 0, Qtop,0,max_ampl);
	pf_AQ->SetMarkerColor(3);
	pf_AQ->SetLineColor(5);
	pf_AQ->SetLineStyle(1);
	pf_AQ->SetLineWidth(1);
	
	
	
	
	double max;
	int i_max;
	double Q;
	float t0;
	float t1;
	long int nentries = tr->GetEntriesFast();
	double AQ_ratio;
	int throwed = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		

		tr->GetEntry(jentry);
		max = GetMax(x, i_max);
		//if (max == 0.715408 || max > 0.2) continue;
		if (max > max_ampl || max < min_ampl) {
			++throwed;
			continue;
		}
		Q = GetQ(t, x, max, i_max, t0, t1, cut_ratio, thr);
		h2_AQ->Fill( Q, max);
		pf_AQ->Fill( Q, max);
		h_Q->Fill(Q);
		h_Ampl->Fill(max);
		h_t0->Fill(t0);
		h_t1->Fill(t1);
		AQ_ratio += max/Q;
		
		for (int i = 0; i < x->size(); ++i){//wfms iterator
			//gr->SetPoint(i,(double)t->at(i),x->at(i));
			wfm_Stack->Fill(t->at(i) , (float) x->at(i) );
		}
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
		
	}
	AQ_ratio = AQ_ratio / (nentries - throwed);
	cout<<"From "<<nentries<<" wfm's, "<<throwed<<" was deleted"<<endl;
	cout<<"A/Q ratio = "<<AQ_ratio<<endl;
	//~ TObject oQA_ratio;
	//~ oQA_ratio.SetUniqueID(QA_ratio);
	//~ oQA_ratio.Write("Q/A_ratio");
	
	TF1 *f_A = new TF1("f_A", "[0] + [1]*x", h_Q->GetXaxis()->GetXmin(),h_Q->GetXaxis()->GetXmax());
	f_A->SetParameters(cfg["hist"]["min_ampl"], AQ_ratio);
	
	pf_AQ->Fit(f_A,"N");
	TParameter<double> * A0 = new TParameter("A0",f_A->GetParameter(0));
	TParameter<double> * A1 = new TParameter("A1",f_A->GetParameter(1));
	A0->Write();
	A1->Write();
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
	
	TCanvas *c1 = new TCanvas("c1","",1600,1200);
	c1->Divide(2,2);
	c1->cd(1);
	TH2F * wfm_Stack = (TH2F*) f->Get("wfm_Stack");
	wfm_Stack->Draw("colz");
	gPad->SetLogz();
	c1->cd(3);
	TH1D *h_t0 = (TH1D*) f->Get("h_t0");
	
	TH1D *h_t1 = (TH1D*) f->Get("h_t1");
	h_t1->SetLineColor(2);
	h_t0->Draw();
	h_t1->Draw("same");
	
	c1->cd(2);
	TH2D * h2_AQ = (TH2D*) f->Get("h2_AQ");
	h2_AQ->Draw("colz");
	TProfile * pf_AQ = (TProfile*) f->Get("pf_AQ");
	pf_AQ->Draw("same L");
	TF1 * f_A = new TF1("f_A","[0] + [1]*x", h_Q->GetXaxis()->GetXmin(), h_Q->GetXaxis()->GetXmax());
	double A0 = ((TParameter<double> *) f->Get("A0"))->GetVal();
	double A1 = ((TParameter<double> *) f->Get("A1"))->GetVal();
	f_A->SetParameters(A0,A1);
	f_A->Draw("same");
	TH1D * h_Ampl_scaled = (TH1D *) h_Ampl->Clone("h_Ampl_scaled");

	h_Ampl_scaled->SetBins( h_Ampl_scaled->GetNbinsX(), - A0/A1, (h_Ampl_scaled->GetXaxis()->GetXmax() - A0) / A1);
	h_Ampl_scaled->SetLineColor(2);
	c1->cd(4);
	
	if (h_Q->GetMaximum() < h_Ampl_scaled->GetMaximum()) h_Q->SetAxisRange(0,h_Ampl_scaled->GetMaximum() +100,"Y");
	h_Q->Draw();
	h_Ampl_scaled->Draw("same");
	

	TH1D * h_QxA = (TH1D *) h_Q->Clone("h_QxA");
	double Q_center;
	for (int i = 0; i < h_QxA->GetNbinsX(); ++i) {
		Q_center = h_QxA->GetBinCenter(i);
		h_QxA->SetBinContent(i, h_QxA->GetBinContent(i) * h_Ampl_scaled->GetBinContent(h_Ampl_scaled->FindBin(Q_center)));
	
	}
	
	
	h_QxA->Scale(h_Q->GetMaximum() / h_QxA->GetMaximum());
	
	h_QxA->SetLineColor(3);
	h_QxA->SetLineStyle(1);
	h_QxA->SetLineWidth(1);
	h_QxA->ResetStats();
	h_QxA->Draw("same");
	string fGrImg;
	if (DataName.rfind('/') < DataName.length()) fGrImg = DataName + "/GrImg_" + DataName.substr(DataName.rfind('/')+1) + ".png";
	else fGrImg = DataName + "/GrImg_" + DataName + ".png";
	c1->Print(fGrImg.data());
}
int GetQ1(TH1D * h_Qs,int i){
	int i_max = 0;
	while (h_Qs->GetNbinsX() > i++) if (h_Qs->GetBinContent(i) > h_Qs->GetBinContent(i_max)) i_max = i;
	return i_max;
}
double GetQ0(TH1D * h_Q){

	int i = 0;
	while(h_Q->GetBinContent(i) < h_Q->GetBinContent(i+1) ) ++i;
	return i; 
}

//~ void FitFuncDraw(TF1 *f){
	//~ Int_t n = 100;
	//~ TCanvas *c2 = new TCanvas("c2", "Fit function parts");
	//~ TF1 * G0 = new TF1("G0",FunG0, f->GetXmin(), f->GetXmax(),8);
	//~ TF1 * Noise = new TF1("Noise",FunNoise, f->GetXmin(), f->GetXmax(),8);
	//~ TF1 * Gn = new TF1("Gn",FunGn, f->GetXmin(), f->GetXmax(),8);
	
	//~ f->SetLineColor(1);
	//~ G0->SetLineColor(3);
	//~ Noise->SetLineColor(6);
	//~ Gn->SetLineColor(2);
	//~ G0->SetLineStyle(2);
	//~ Noise->SetLineStyle(2);
	//~ Gn->SetLineStyle(2);
	//~ G0->SetParameters(f->GetParameters());
	//~ Noise->SetParameters(f->GetParameters());
	//~ Gn->SetParameters(f->GetParameters());
	
	//~ f->Draw();
	
	//~ G0->Draw("same");
	//~ Noise->Draw("same");
	//~ Gn->Draw("same");
	

//~ }
void InitFitFun(TF1 * f_Q, json fit) {
	f_Q->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu");
	
	int pNum;
	string pKey;
	double pVal, pValB, pValT;
	for (auto it = fit.begin(); it != fit.end();){
		pKey = it.key();
		pNum = f_Q->GetParNumber(pKey.data());
		pVal = it++.value();
		pValB = it++.value();
		pValT = it++.value();
		f_Q->SetParameter(pNum, pVal);
		if (pValB < pValT) f_Q->SetParLimits(pNum, pValB, pValT);
		cout<< pNum << " - " << pKey << " = " << pVal << " B: " << pValB << " T: " << pValT << endl;
	}
}
int GetValey(TH1D * h_Qs){
	int i = 2;
	while ( h_Qs->GetBinContent(i)>h_Qs->GetBinContent(i+1) ) ++i;
	
	cout<<"Valey = "<<h_Qs->GetBinCenter(i)<<":"<<h_Qs->GetBinContent(i)<<" - "<<i<<endl;
	return i;
}

void FitFunFromcfg(string DataName = ""){
	string fCfg = DataName + "/cfg.json";
	ifstream cfg_file(fCfg);
	json cfg = json::parse(cfg_file);
	json fit = cfg["fit"];
	TCanvas *c = new TCanvas("c","fitfun");
	TF1 *f_Q = new TF1("f_Q",Fun, 0 , 10,7);
	InitFitFun(f_Q,fit);
	f_Q->Draw();
	
}


void DoFit(string DataName = ""){
	string fHistsName = DataName + "/OscLecHists.root";
	string fCfg = DataName + "/cfg.json";
	
	TFile * f = TFile::Open(fHistsName.data());
	ifstream cfg_file(fCfg);
	json cfg = json::parse(cfg_file);
	json fit = cfg["fit"];
	
	
	
	int bins = cfg["hist"]["bins"];
	double A0 = ((TParameter<double> *) f->Get("A0"))->GetVal();
	double A1 = ((TParameter<double> *) f->Get("A1"))->GetVal();
	TH1D * h_Q = (TH1D*) f->Get("Q");
	h_Q->Scale(1./h_Q->GetEntries());
	TH1D * h_Qs = (TH1D*) h_Q->Clone("h_Qs");
	h_Qs->SetLineColor(2);
	TH1D * h_A = (TH1D*) f->Get("h_Ampl");
	
	h_Qs->ResetStats();
	
	TCanvas *c_s = new TCanvas("c_s","sandbox");
	c_s->SetWindowPosition(800,0);
	h_Q->SetStats(0);
	h_Qs->Smooth();
	int i_Valey = GetValey(h_Qs);
	int i_Q0 = GetQ0(h_Q);
	double Q0 = h_Qs->GetBinCenter(i_Q0);
	h_A->Smooth();	
	h_A->SetBins( h_A->GetNbinsX(), -A0/A1, (h_A->GetXaxis()->GetXmax() - A0) / A1);
	h_A->Scale(1./h_A->GetEntries());
	h_A->Draw("HIST");
	h_Qs->Draw("HIST same");
	
	TCanvas *c_Q = new TCanvas("c_Q","Q");
	h_Qs->Draw("HIST");
	h_Q->Draw("HIST same");
	
	
	
	
	
	
	
	
	//~ h_AQ->Fit(f_AQ);
	TMarker * p_Valey = new TMarker(h_Q->GetBinCenter(i_Valey),h_Q->GetBinContent(i_Valey),1);
	TF1 *f_Q = new TF1("f_Q",Fun, Q0, h_Q->GetXaxis()->GetXmax(),7);
	f_Q->SetLineColor(3);  // later change to Valey
	f_Q->SetLineWidth(2);
	
	
	InitFitFun(f_Q, fit);//loads parameters from config file

	//~ f_Q->SetParameters(Q0,0.075,4.1,1,0.7,0.5,2,6000);
	//~ f_Q->SetParLimits(0,-3,3);
	//~ f_Q->SetParLimits(1,0,5);
	//~ f_Q->SetParLimits(2,0,10);
	//~ f_Q->SetParLimits(3,0,5);spectrum</         <. .>>>>>>>>> ++00000000000000000000000000000000000
	
	
	
	
	
	
	
	
	
	
	//~ f_Q->SetParLimits(4,0,1);
	//~ f_Q->SetParLimits(6,0.5,20);
	
	
	
	f_Q->SetParameter("Q0",Q0);
	TF1 *f_N = new TF1("f_N","[0]*[1]*exp(-([1]*(x-[2])))",Q0,h_Q->GetXaxis()->GetXmax());
	f_N->SetParNames("w","alpha","Q0");
	f_N->SetParameters(fit["w"],fit["alpha"],Q0);
	f_N->FixParameter(2,Q0);
	h_Q->Fit(f_N,"B","",Q0,h_Q->GetBinCenter(i_Valey)*1.1);
	f_N->Draw("same");
	int i_Q1 = GetQ1(h_Qs,i_Valey);
	double Q1 = h_Q->GetBinCenter(i_Q1);
	double w = f_N->GetParameter(0);
	double alpha = f_N->GetParameter(1);
	double sig1;
	f_Q->SetParameter("w",w);
	f_Q->SetParameter("alpha",alpha);
	sig1 = (Q1 - h_Q->GetBinCenter(i_Valey)/2);
	Q1 = Q1 - Q0 - w/alpha;
	
	double k = h_Qs->GetBinContent(i_Q1)*sig1*1.8119051;
	double mu = - utl::LambertW<0>(-k);
	f_Q->SetParameter("mu", mu);
	
	f_Q->SetParameter("Q1",Q1);
	f_Q->SetParameter("sig1",sig1);
	f_Q->SetParLimits(3,sig1/4,sig1*1.5);
	cout<<"Q1 = "<<Q1<<endl;
	cout<<"sig1 = "<<sig1<<endl;
	cout<<"Qsh = "<<w/alpha<<endl;
	cout<<"mu = "<<mu<<endl;




	


	f_Q->Draw("same");
	//~ TCanvas *c = new TCanvas();
	//~ f_Q->Draw();
	cout<<"\n\n_______________\n";
	for(int i = 0; i<7;++i)cout<<f_Q->GetParName(i)<<" = "<<f_Q->GetParameter(i)<<endl;
	
	c_s->Update();
	//~ int i = 0;
	//~ for(double W = 0.1; W < 1  ; W += 0.05 )
	//~ for(double M = 2; M <  8; M += 0.1)
	//~ for(double A = 0.5; A <  5; A += 0.01)
	//~ {
		//~ f_Q->SetParameter("w",W);
		//~ f_Q->SetParameter("alpha", A);
		//~ f_Q->SetParameter("mu", M);
		
		//~ cout<<++i<<endl<<"w = "<<W<<"\nalpha = "<<A<<"\nmu = "<<M<<endl;
		//~ c_Q->Draw();
		//~ c_Q->Update();
		//~ cin.get();
		
	//~ }
	
	
	cout<<endl;
	//~ f_Q->FixParameter(0,Q0);
	
	h_Q->Fit(f_Q,"","", h_Q->GetBinCenter(i_Valey)*0.7, h_Q->GetXaxis()->GetXmax());//,Q0,Q->GetXaxis()->GetXmax());
	
	double chi = f_Q->GetChisquare()/( (double) f_Q->GetNDF());
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
		flFitPars<<f_Q->GetParName( i )<<"\t"<<f_Q->GetParameter(i)<<"\t"<<f_Q->GetParError(i)<<endl;
	}
	flFitPars<<"chi/n\t"<<chi<<endl;
	flFitPars.close();
	//c_Q->cd();
	//f_Q->Draw("same");	
	
}
void ParEst(TH1D *h_Q, int i_Q0, double &Q0, double &Q1, double &Qp, double &Qv, double &sig0, double &sig1, double &w, double &alpha, double &mu, vector<int> v_Qv, vector<int> v_Qp, vector<int> v_dQv, vector<int> v_dQp)
{
	int i_Qv, i_Qp;
	if (v_Qv.size() > 0 && v_Qp.size() > 0){
		i_Qv = v_Qv.at(0);
		i_Qp = v_Qp.at(0);
	}else  //here will be the magic
	{
		i_Qv = v_dQv.at(0);
		i_Qp = v_dQp.at(0);
	}
	Q0 = h_Q->GetBinCenter(i_Q0);
	Qv = h_Q->GetBinCenter(i_Qv);
	Qp = h_Q->GetBinCenter(i_Qp);
	//noise parameters
	TF1 *f_N = new TF1("f_N",Funoise,Q0,Qv,5);
	f_N->SetParameters(Q0,(Qv-Q0)/10,0.5,(Qv-Q0)/10, h_Q->Integral(i_Q0, i_Qv));
	f_N->SetParNames("Q0","sig0","w","alpha","A");
	f_N->SetParLimits(0,0,Qv);
	f_N->SetParLimits(1,0,(Qv-Q0)/2);
	f_N->SetParLimits(2,0,1);
	f_N->SetParLimits(4,0,1);
	cout<<endl<<"NOISE PARAMETER ESTIMATION FIT\n_______________________________\n";
	h_Q->Fit(f_N,"R");
	f_N->Draw("same");
	
	
	Q0 = f_N->GetParameter("Q0");
	sig0 = f_N->GetParameter("sig0");
	alpha = f_N->GetParameter("alpha");
	w = f_N->GetParameter("w");
	
	TF1 *f_Nexp = new TF1("f_N","[0]*[1]*exp(-([1]*(x-[2])))",Q0,Qv);
	f_Nexp->SetParNames("w","alpha","Q0");
	f_Nexp->SetParameters(w,alpha,Q0);
	f_Nexp->FixParameter(2,Q0);
	f_Nexp->SetLineColor(3);
	h_Q->Fit(f_Nexp,"R");
	f_Nexp->Draw("same");
	alpha = f_Nexp->GetParameter("alpha");
	w = f_Nexp->GetParameter("w");
	
	
	
	
	//gauss bunch parameters -- Q1, sig1, mu
	Q1 = Qp - Q0 - w/alpha;
	sig1 = (Qp - Qv)/2;
	mu = - utl::LambertW<0>(-h_Q->GetBinContent(i_Qv)*sig1*1.8119051);
}
	
void FitHist( TH1D * h_Q){
	double Q0,Q1,Qv,Qp,sig0,sig1,w,alpha,mu;
	int i_Q0, i_Q1, i_Qv, i_Qp, i_end, i_dQ0;
	int i = 0;
	vector<int> v_Qp, v_Qv, v_dQp, v_dQv;
	TGraph *g_Qp = new TGraph();
	TGraph *g_Qv = new TGraph();
	TGraph *g_dQp = new TGraph();
	TGraph *g_dQv = new TGraph();
	
	TH1D * h_Q0 = (TH1D*) h_Q->Clone("h_Q0");
	TH1D * h_Qs = (TH1D*) h_Q->Clone("h_Qs");
	h_Qs->Scale(1/h_Qs->GetEntries());
	h_Qs->Smooth();
	h_Qs->ResetStats();
	
	TH1D * h_dQ = (TH1D*) h_Qs->Clone("h_dQ");
	
	//go to Q0 -> pedestal gauss peak
	while ( h_Qs->GetBinContent(i-1) < h_Qs->GetBinContent(i++));
	i_Q0 = i;
	cout<<"i_Q0 = "<<i<<endl;
	//finding local extremes in spectra
	double thr = h_dQ->GetMaximum()*1e-2;
	while (h_dQ->GetBinContent(++i) > thr) {
		if (h_dQ->GetBinContent(i-1) > h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) <= h_dQ->GetBinContent(i+1)) v_Qv.push_back(i);
		if (h_dQ->GetBinContent(i-1) < h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) >= h_dQ->GetBinContent(i+1)) v_Qp.push_back(i);
	}
	i_end = i;

	
	for(int v,p,j = 0; j < v_Qv.size();++j){
		v = v_Qv.at(j);
		p = v_Qp.at(j);
		g_Qv->AddPoint(h_Qs->GetBinCenter(v), h_Qs->GetBinContent(v));
		g_Qp->AddPoint(h_Qs->GetBinCenter(p), h_Qs->GetBinContent(p));
	}
	//~ cout<<v_Qp.size()<<endl;
	
	//~ cout<<"i_Qv = "<<v_Qv<<endl;
	//~ cout<<"i_Qp = "<<v_Qp<<endl;

	double dQ;
	for (int j = 1; j < i_end; j++ ) {
		dQ = h_dQ->GetBinContent(j) - h_dQ->GetBinContent(j+1);
		h_dQ->SetBinContent(j,dQ);
		//~ if (dQ > 0) h_dQ->SetBinContent(j,dQ);
		//~ else h_dQ->SetBinContent(j,0);
	}
	h_dQ->Smooth(10);
	h_dQ->ResetStats();
	h_dQ->Scale(1/h_dQ->GetEntries());
	
	i = 0;
	while ( h_dQ->GetBinContent(i) < h_dQ->GetBinContent(++i));
	
	cout<<"dQ0 - "<<h_dQ->GetBinCenter(i)<<endl;
	i_dQ0 = i;
	
	while (++i < i_end) {
		if (h_dQ->GetBinContent(i-1) >= h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) < h_dQ->GetBinContent(i+1)) v_dQv.push_back(i);
		if (h_dQ->GetBinContent(i-1) <= h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) > h_dQ->GetBinContent(i+1)) v_dQp.push_back(i);
	}
	for(int v,p,j = 0; j < min(v_dQv.size(), v_dQp.size());++j){
		v = v_dQv.at(j);
		p = v_dQp.at(j);
		g_dQv->AddPoint(h_dQ->GetBinCenter(v), h_dQ->GetBinContent(v));
		g_dQp->AddPoint(h_dQ->GetBinCenter(p), h_dQ->GetBinContent(p));
	}
	
	TGraph *g_W = new TGraph();
	g_W->SetLineColor(1);
	g_W->SetLineWidth(1);
	g_W->SetLineStyle(2);
	
	g_W->AddPoint(h_Qs->GetBinCenter(i_Q0),h_Qs->GetBinContent(i_end));
	g_W->AddPoint(h_Qs->GetBinCenter(i_Q0),h_Qs->GetBinContent(i_Q0));
	g_W->AddPoint(h_Qs->GetBinCenter(i_end),h_Qs->GetBinContent(i_Q0));
	g_W->AddPoint(h_Qs->GetBinCenter(i_end),h_Qs->GetBinContent(i_end));
	g_W->AddPoint(h_Qs->GetBinCenter(i_Q0),h_Qs->GetBinContent(i_end));
	//~ thr = h_dQ->GetMaximum()*1e-2;
	
	

	TCanvas * c = new TCanvas("c","",1000,1400);
	
	c->Divide(1,2);
	
	g_dQp->SetMarkerStyle(8);
	g_dQp->SetMarkerColor(3);
	g_dQp->SetMarkerSize(1);
	g_dQv->SetMarkerStyle(8);
	g_dQv->SetMarkerColor(2);
	g_dQv->SetMarkerSize(1);
	g_Qp->SetMarkerStyle(8);
	g_Qp->SetMarkerColor(3);
	g_Qp->SetMarkerSize(1);
	g_Qv->SetMarkerStyle(8);
	g_Qv->SetMarkerColor(2);
	g_Qv->SetMarkerSize(1);
	
	
	c->cd(1);
	h_Q0->Scale(1/h_Q0->GetEntries());
	h_Q0->Draw("HIST");
	cout<<"vbbl "<<h_Qs->GetBinCenter(i_end)<<endl;
	h_Q0->GetXaxis()->SetRange(0, i_end*1.2);
	//~ h_Q0->GetYaxis()->SetRange(0,h_Qs->GetBinContent(i_Q0));
	g_W->Draw("same");
	
	
	h_Qs->SetLineColor(2);
	h_Qs->Draw("HIST same");
	g_Qp->Draw("P same");
	g_Qv->Draw("P same");
	c->cd(2);
	h_dQ->Draw("HIST");
	g_W->Draw("same");
	g_dQp->Draw("P same");
	g_dQv->Draw("P same");
	
	
	//case 1: spectre has valey and peak so parameters can be estimated from that
	//case 2: spectre has only negative derivative so parameters will be estimated from inflection points
	//in both cases, we can use inflection points to estimate Q1 (or try atleast)
	
	c->cd(1);
	//case 1
	ParEst(h_Qs, i_Q0, Q0, Q1, Qp, Qv, sig0, sig1, w, alpha, mu, v_Qv, v_Qp, v_dQv, v_dQp);
	TF1 *f_Q = new TF1("f_Q",Fun, Qv, h_Qs->GetBinCenter(i_end),7);
	f_Q->SetLineColor(5);
	f_Q->SetParNames("Q0","sig0","Q1","sig1","w","alpha","mu");
	f_Q->SetParameters(Q0,sig0,Q1,sig1,w,alpha,mu);
	f_Q->FixParameter(0,Q0);
	f_Q->FixParameter(1,sig0);
	f_Q->SetParLimits(2,Qv,Qp+Qv);
	f_Q->SetParLimits(3, (Qp-Qv)*0.2, (Qp-Qv)*2);
	f_Q->SetParLimits(4,0,1);
	f_Q->SetParLimits(5,0,h_Qs->GetBinCenter(i_end));
	f_Q->SetParLimits(6,0,1e2);
	TF1 *f_Qe = (TF1*) f_Q->Clone("f_Qe");
	f_Qe->SetLineColor(4);
	f_Qe->Draw("same");
	
	
	cout<<"\nParameter estimation result:\n";
	for(int i = 0; i<7;++i)cout<<f_Q->GetParName(i)<<" = "<<f_Q->GetParameter(i)<<endl;
	f_Q->Draw("same");
	h_Q0->Fit(f_Q,"R");
	
	
	//todo noise fit a odhad parametru a pak fit! wuhuu
	
	
	
	//~ h_Qs->SetLineColor(2);
	//~ TCanvas * c = new TCanvas("c","bla",1200,1000);
	//~ c->Divide(1,2);
	//~ c->cd(1);
	//~ h_Q->Draw();
	//~ h_Qs->Draw("same");
	//~ TH1D * h_dQ = (TH1D*) h_Qs->Clone("h_dQ");
	
	//~ h_dQ->SetBinContent(0, 0);
	//~ for (int i = 1; i < h_Q->GetNbinsX(); ++i) h_dQ->SetBinContent(i, h_Qs->GetBinContent(i) - h_Qs->GetBinContent(i+1));
	//~ h_dQ->Smooth();
	//~ c->cd(2);
	//~ h_dQ->Draw();
	
	//~ TH1D * h_d2Q = (TH1D*) h_dQ->Clone("h_d2Q");
	//~ h_d2Q->SetBinContent(0, 0);
	//~ for (int i = 1; i < h_Q->GetNbinsX(); ++i) h_d2Q->SetBinContent(i, h_dQ->GetBinContent(i) - h_dQ->GetBinContent(i+1));
	//~ h_d2Q->Smooth();
	
	//~ h_d2Q->SetLineColor(3);
	//~ h_d2Q->Draw("same");

}
void FitHistTB(){
	
	//~ TFile * f = TFile::Open("dr_smalwfms/TSA.root");
	//~ TTree * t = (TTree*) f->Get("TSA");
	//~ int ientry = 30;
	//~ TH1D * h_Q = new TH1D();
	//~ t->SetBranchAddress("hQ",&h_Q);
	//~ t->GetEntry(ientry);
	
	
	TFile * f = TFile::Open("dr_smalwfms/OscLecHists.root");
	TH1D * h_Q = (TH1D*) f->Get("Q");
	
	
	FitHist(h_Q);
	
		
	//~ FitHist(h);	
}

void DoAPA(string DataName)//after pulse analysis
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
