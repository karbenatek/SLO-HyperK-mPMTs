
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <cstddef>
#include <iostream>
#include <TGraph.h>
#include <vector>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TMath.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TF1.h>
#include <filesystem>
namespace fs = std::filesystem;
#include <filesystem>
#include <limits>
#include "../lib/LambertW/LambertW.cc"

#include "json.hpp"
using json = nlohmann::json;


int GetT0(vector<float> * t, vector<double> * x, double thr = 0.02){
	int i = 0;
	int x_size = x->size() - 1;
	while(t->at(i) < 0) ++i;
	while(i < x_size && x->at(i) < thr) ++i;
	
	return i;
}
double GetQ(vector<float> * t, vector<double> * x, double max, int i_max, float & T0, float & t0, float & t1, double cut_ratio = 0.5, double thr = 0.02){
	
	
	int i = i_max;
	double Q = 0;
	int i_T0 = GetT0(t,x,thr);
	T0 = t->at(i_T0);
	double cut_thr = cut_ratio * max + (1-cut_ratio) * x->at(i_T0);
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
	cout<<"\nSetting cfg.json to default in "<<DataName<<endl;
	fs::copy("cfg.json", fCfg, fs::copy_options::overwrite_existing);
}
void CopyCfg(string from = "", string to = ""){
	from = from + "/cfg.json";
	to = to + "/cfg.json";
	cout<<"\nCopying cfg.json "<<from<<" to "<<to << endl;
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
	double Amax = cfg["hist"]["max_ampl"]; //max wfm amplitude
	double Amin = cfg["hist"]["min_ampl"]; //min wfm amplitude
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

	TH2F * wfm_Stack = new TH2F("wfm_Stack","Waveform stack",t->size(),t->front(),t->back(),bins,-0.05,Amax);
	TH1F * h_t0 = new TH1F( "h_t0", "t0", t->size(), t->front(), t->back());
	TH1F * h_t1 = new TH1F( "h_t1", "t1", t->size(), t->front(), t->back());
	TH1F * h_tts = new TH1F("h_tts","Transit time distribution",t->size(), t->front(), t->back());
	
	vector<double> v_Q;
	vector<double> v_A;




	
	
	
	
	double A;
	int i_max;
	double Q, Qmax = 0;
	float t0;
	float t1;
	float T0;
	long int nentries = tr->GetEntriesFast();
	double AQ_ratio;
	int throwed = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		

		tr->GetEntry(jentry);
		A = GetMax(x, i_max);
		//if (max == 0.715408 || max > 0.2) continue;
		if (A > Amax || A < Amin) {
			++throwed;
			continue;
		}
		Q = GetQ(t, x, A, i_max, T0, t0, t1, cut_ratio, thr);
		if (Q > Qmax) Qmax = Q;
		v_Q.push_back(Q);
		v_A.push_back(A);


	
		h_tts->Fill(T0);
		h_t0->Fill(t0);
		h_t1->Fill(t1);
		AQ_ratio += A/Q;
		
		for (int i = 0; i < x->size(); ++i){//wfms iterator
			//gr->SetPoint(i,(double)t->at(i),x->at(i));
			wfm_Stack->Fill(t->at(i) , (float) x->at(i) );
		}
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
		
	}
	AQ_ratio = AQ_ratio / (nentries - throwed);
	cout<<"From "<<nentries<<" wfm's, "<<throwed<<" was deleted"<<endl;
	cout<< "A/Q ratio = "<< AQ_ratio<< endl;

	TH1D * h_Ampl = new TH1D("h_Ampl","Amplitudes", bins, 0, Amax);
	TH1D * h_Q = new TH1D("Q","Charge", bins, 0, Qmax);
	TH2D * h2_AQ = new TH2D("h2_AQ","Voltage/Charge",bins, 0, Qmax, bins,0, Amax);
	TProfile * pf_AQ = new TProfile
	("pf_AQ","Voltage/Charge",bins, 0, Qmax, 0, Amax);
	pf_AQ->SetMarkerColor(3);
	pf_AQ->SetLineColor(5);
	pf_AQ->SetLineStyle(1);
	pf_AQ->SetLineWidth(1);

	for (int i = 0; i < nentries - throwed; ++i) {
		Q = v_Q[i];
		A = v_A[i];
		
		h2_AQ->Fill( Q, A);
		pf_AQ->Fill( Q, A);
		h_Q->Fill(Q);
		h_Ampl->Fill(A);
	}
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
void PlotTTS(string DataName = ""){
	string fHistsName = DataName + "/OscLecHists.root";
	string fImgName1 = DataName + "/TTS.png";
	string fImgName2 = "img/tts/" + DataName.substr(DataName.rfind("/")+1) + "TTS.png";

	TFile * f = TFile::Open(fHistsName.data());
	TH1F * h_tts = (TH1F *) f->Get("h_tts");


	TF1 *f_tts = new TF1("f_tts","gaus",h_tts->GetXaxis()->GetXmin(), h_tts->GetXaxis()->GetXmax());
	f_tts->SetParameters(h_tts->GetMaximum(), h_tts->GetMean(), h_tts->GetStdDev());
	int i_l, i_r, i = h_tts->GetMaximumBin(), dir = 1;

	while (h_tts->GetBinContent(i) >= h_tts->GetMaximum()*0.2) ++i;
	i_r = i;
	i = h_tts->GetMaximumBin();
	while (h_tts->GetBinContent(i) >= h_tts->GetMaximum()*0.2) --i;
	i_l = i;
	
	h_tts->Fit(f_tts,"","",h_tts->GetBinCenter(i_l), h_tts->GetBinCenter(i_r));

	TCanvas *c_tts = new TCanvas("c_tts","Transit time spread");
	gStyle->SetOptFit(0111);
	h_tts->SetAxisRange(f_tts->GetParameter(1)-5*f_tts->GetParameter(2), f_tts->GetParameter(1)+20*f_tts->GetParameter(2));
	// h_tts->GetListOfFunctions()->Add(f_tts);
	h_tts->Draw();
	c_tts->Print(fImgName1.data());
	c_tts->Print(fImgName2.data());

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
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
bool overzero(vector<float> * v_t, int i, bool dir){
	
	if (sgn(i) * sgn(v_t[i-1 +2*dir]) > 0) return 0;
	else return 1;
}
void DoAPA(string MeasName)//after pulse analysis
{
	string fDataName = MeasName + "/AP/OscLecData.root";
	//how? find i , t_trig, 
	TFile * f = TFile::Open(fDataName.data());
	TTree * tr = (TTree*) f->Get("OscLecData");

	int i_Wfm;
	vector<double>  * v_x = 0;
	vector<float>   * v_t = 0;
	double t, t0;
	tr->SetBranchAddress("x",&v_x);
	tr->SetBranchAddress("t",&v_t);
	tr->SetBranchAddress("TimeSinceSeg1",&t);
	
	tr->GetEntry(0);
	// find trig time index
	int i0 = v_t->size()/2;
	if (v_t->at(i0) > 0) {
		while (v_t->at(i0) >= 0) i0--; 
		i0++;
	}
	else 
		while (v_t->at(i0) < 0) i0++;

	int trigs = 0;
	int pulses = 0;
	bool trig;
	long int nentries = tr->GetEntriesFast();
	TH1I * h_AP = new TH1I("h_AP", "After pulse counts histogram", 11, -0.5, 10.5);
	TH1F * h_t = new TH1F("h_t", "Pulse time", 100, 0,0.0012);
	TH1D * h_V = new TH1D("h_V", "", 100, -5, 0);
	for (Long64_t jentry=0; jentry<nentries; jentry++) {//events iterator
		tr->GetEntry(jentry);
		trig = v_x->at(i0) < 0;
		if (trig) {
			h_V->Fill(v_x->at(i0 - 10));
			if (trigs++ > 0) {
				h_AP->Fill( pulses );

			}
			pulses = 0;
			h_t->Fill(t - t0);

			t0 = t;
		}
		else {
			pulses++;
		}
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
	}
	cout<<"Done\n\n";
	
	h_V->Draw();
	
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
	float edge_t, cut_t, T0;
	long int nentries = tr->GetEntriesFast();
	tr->GetEntry(0);
	int dnQ;
	int next_i_Seg;
	int merged = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		tr->GetEntry(jentry);
		if( (jentry)%1000 == 0 && jentry != 0) cout<<jentry<<"/"<<nentries<<endl;
		max = GetMax(x, i_max);
		Q = GetQ(t, x, max, i_max, T0, edge_t, cut_t, cut_ratio);
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
void PlotGain(string MeasName){ // e.g.: data/m0
	// TODO: filtrovat pouze slozky s cislem, z nazvu slozky vzit napeti, 
	string GName = MeasName + "/G";
	string f_ImgName = MeasName + "/gain.png";
	std::size_t idx;
	string f_Vname;
	vector<array<double,3>> v_Q;

	int i = 0;
	double V, Q, sig;
	json j_fit;
	ifstream fi_fit;
	for (auto &dir_entry : fs::directory_iterator(GName)){
		f_Vname = dir_entry.path();
		fi_fit.open(f_Vname + "/fit_" + f_Vname.substr(f_Vname.rfind("/")+1) + ".json");
		j_fit = json::parse(fi_fit);
		
		V = stod(f_Vname.substr(f_Vname.rfind("/")+1));
		Q = j_fit["Q1"]["val"];
		sig = j_fit["sig1"]["val"];
		v_Q.push_back({V, Q, sig});
		fi_fit.close();

		// cout<<"V: "<<V<<endl<<"Q: "<<Q<<endl<<"Qe: "<<Qe<<endl<<endl;
		

		// v_Q[i++]
	}
	double a_V[v_Q.size()], a_Ve[v_Q.size()], a_Q[v_Q.size()], a_Qe[v_Q.size()];
	for (int j = 0; j < v_Q.size(); ++j){
		a_V[j] = v_Q[j][0];
		cout<<v_Q[j][2]<<endl;
		a_Ve[j] = 0;
		a_Q[j] = v_Q[j][1];
		a_Qe[j] = v_Q[j][2];
	}
	TCanvas * c_G = new TCanvas("c_G","Gain plot");
	TGraphErrors *g_Q = new TGraphErrors(v_Q.size(), a_V, a_Q, a_Ve, a_Qe);
	g_Q->SetLineColor(1);
	g_Q->SetMarkerColor(1);
	g_Q->SetMarkerSize(1);
	g_Q->SetMarkerStyle(21);
	g_Q->SetTitle("Gain; Bias Voltage (V); Q1 (a.u.)");
	
	g_Q->Draw("AP");
	cout<<f_ImgName<<endl;
	c_G->Print(f_ImgName.data());


}

