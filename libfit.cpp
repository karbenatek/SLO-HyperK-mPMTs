#include <TFitResultPtr.h>
#include <TTree.h>
#include <TObject.h>
#include <TGraph.h>
#include <TVirtualPad.h>
#include <cstdlib>
#include <vector>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TFile.h>

const int Nfpars = 9;
const int Npars = 12;


const array<string,Nfpars> fparnames = {"Q0","sig0","Q1","sig1","Qshift","w","alpha","mu","A"};
const array<string,Npars> parnames = {"Q0", "Q1", "Qv", "Qp", "Qend", "Qshift", "sig0", "sig1", "w", "alpha", "mu", "A"};

void Drawfun(TF1 *f, TLegend * leg, int n = 10){
	double Q0 = f->GetParameter("Q0");
	double sig0 = f->GetParameter("sig0");
	double Q1 = f->GetParameter("Q1");
	double sig1 = f->GetParameter("sig1");
	double w = f->GetParameter("w");
	double alpha = f->GetParameter("alpha");
	double mu = f->GetParameter("mu");
	double A = f->GetParameter("A");
	double Qshift = f->GetParameter("Qshift");
	f->SetLineWidth(2);
	f->SetLineColor(2);
	f->Draw("same");
	TF1 * f_Q0 = new TF1("f_Q0",FunoiseQ0, f->GetXmin(), f->GetXmax(),4);
	f_Q0->SetParameters(Q0, sig0, w, A);
	f_Q0->SetLineColor(3);
	f_Q0->SetLineWidth(1);
	f_Q0->Draw("same");
	leg->AddEntry(f_Q0,"Pedestal Gauss","l");
	TF1 * f_Exp = new TF1("f_Exp",FunoiseExp, Q0, f->GetXmax(),4);
	f_Exp->SetParameters(Q0, w, alpha, A);
	f_Exp->SetLineColor(4);
	f_Exp->SetLineWidth(1);
	f_Exp->Draw("same");
	leg->AddEntry(f_Exp,"Exponential noise","l");
	TF1 *f_s; 
	vector <TF1*> v_f;
	string f_name;
	for (int i = 1; i <= n; ++i){
		f_name = "f_s" + to_string(i);
		v_f.push_back(new TF1( f_name.data(), Funs, f->GetXmin(), f->GetXmax(), 6));
		v_f.back()->SetParameters(Qshift, Q1, sig1, mu, A, i);
		v_f.back()->SetLineColor(6);	
		v_f.back()->SetLineWidth(1);	
		v_f.back()->Draw("same");
	}	
	leg->AddEntry(v_f.back(),"Gauss bunch","l");
}
void  ParEst1(TH1D *h_Q, TH1D *h_Qs, int i_Q0, int i_end, double *pars[12], vector<int> v[4])
{
	double Q0, Q1, Qv, Qp, Qend, Qshift, sig0, sig1, w, alpha, mu, A;


    vector<int> *v_Qp = &v[0], *v_Qv = &v[1], *v_dQp = &v[2], *v_dQv = &v[3];

	h_Q->Draw();
	h_Qs->Draw("same");
	
	Q0 = *pars[0];
	Qv = *pars[2];
	Qp = *pars[3];
    int i_Qv = h_Q->FindBin(Qv), i_Qp = h_Q->FindBin(Qp);
	
	Qend = *pars[4];
	sig0 = (Qv-Q0)/10;
	alpha = 1/(Qv-Q0);
	A = h_Q->GetEntries();

	TF1 *f_Nexp = new TF1("f_N","[0]*[1]*exp(-([1]*(x-[2])))",Q0,Qend);
	f_Nexp->SetParNames("A","alpha","Q0");
	f_Nexp->SetParameters(A,alpha,Q0);
	A = h_Qs->GetBinContent(i_Qv)/f_Nexp->Eval(Qv);
	f_Nexp->SetParameter("A",A);
	f_Nexp->FixParameter(2,Q0);
	f_Nexp->SetLineColor(3);
	f_Nexp->SetLineWidth(1);
	h_Qs->Fit(f_Nexp,"","",Q0,Qv);

	A = f_Nexp->GetParameter("A");
	alpha = f_Nexp->GetParameter("alpha");
	f_Nexp->Draw("same");
	Q1 = Qp -Q0 - w/alpha;
	sig1 = (Qp - Qv)/1.4;


	TH1D * h_Qsig = (TH1D*) h_Qs->Clone("h_Qsig");
	double dQ;
	Qshift = -(Q0 + w/alpha);
	for (int i = 0; i<= h_Q->GetNbinsX();++i)
		if (i >= i_Qv && i <= i_end) {
			dQ = h_Qs->GetBinContent(i) - f_Nexp->Eval(h_Qs->GetBinCenter(i));
			//~ if (dQ < 0) dQ = 0;
			h_Qsig->SetBinContent(i, dQ);
		}
		else h_Qsig->SetBinContent(i,0);
	h_Qsig->ResetStats();
	h_Qsig->Draw("same hist");
	double A1 = h_Qsig->GetMaximum();
	double QA1 = h_Qsig->GetBinCenter(h_Qsig->GetMaximumBin());

	mu = A1/A*sig1*2.5;
	// mu = 2;
	// mu = - utl::LambertW<0>(-A1*sig1*1.8119051);
	
	TF1 *f_sig = new TF1("f_sig",Funsig,Qv,Qend,5);
	f_sig->SetLineColor(2);
	f_sig->SetParNames("Qshift","Q1","sig1","mu","A");
	f_sig->SetParameters(Qshift,Q1,sig1,mu, A);
	A = A*A1/f_sig->Eval(QA1);
	f_sig->SetParameter("A",A);
	f_sig->Draw("same");
	
	
	Qshift = 0;
	w = 1;
	
    

	int i = 0;
    for (double par : {Q0, Q1, Qv, Qp, Qend, Qshift, sig0, sig1, w, alpha, mu, A}){ 
		*pars[i++] = par;
	}
}
void ParEst2(TH1D *h_Q, TH1D *h_Qs, int i_Q0, int i_end, double *pars[12], vector<int> v[4])
{
	double Q0, Q1, Qv, Qp, Qend, Qshift, sig0, sig1, w, alpha, mu, A;

    vector<int> *v_Qp = &v[0], *v_Qv = &v[1], *v_dQp = &v[2], *v_dQv = &v[3];

	h_Q->Draw();
	h_Qs->Draw("same");

	Q0 = *pars[0];
	Qv = *pars[2];
	Qp = *pars[3];
    int i_Qv = h_Q->FindBin(Qv), i_Qp = h_Q->FindBin(Qp);
	
	Qend = *pars[4];
	sig0 = (Qv-Q0)/10;
	w = 0.5;
	// A = h_Q->GetEntries();
	A = h_Q->GetMaximum();

	cout<<"\n\nblah = "<<A<<endl;
	alpha = 1/(Qv-Q0);

	//noise parameters
	TF1 *f_N = new TF1("f_N",Funoise,Q0,Qend,5);
	f_N->SetParNames("Q0","sig0","w","alpha","A");
	//~ f_N->SetLineWidth(1);
	f_N->SetLineStyle(2);
	f_N->SetLineColor(6);

	f_N->SetParameters(Q0,sig0, w,alpha, A);
	cout<<"Qv = "<<Qv<<endl;
	cout<<"i_Qv = "<<i_Qv<<endl;
	// A = A*h_Qs->GetBinContent(i_Qv)/f_N->Eval(Qv);
	f_N->SetParameter("A",A);
	f_N->SetParLimits(0,0,Qv);
	f_N->SetParLimits(1,0,(Qv-Q0)/2);
	f_N->SetParLimits(2,0, 1);
	f_N->SetParLimits(4,0, A);
	cout<<"\nEstimated noise parameters: \n";
	for(int i = 0; i<5;++i)cout<<f_N->GetParName(i)<<" = "<<f_N->GetParameter(i)<<endl;

	cout<<endl<<"NOISE PARAMETER ESTIMATION FIT\n_______________________________\n";
	// h_Q0->Fit(f_N,"","",Q0,Qv);
	f_N->Draw("same");
	
    // subtraction of noise fit result from original spectre
	TH1D * h_Qsig = (TH1D*) h_Qs->Clone("h_Qsig");
	double dQ;
	for (int i = 0; i<= h_Qs->GetNbinsX();++i)
		if (i >= i_Qv && i <= i_end) {
			dQ = h_Qs->GetBinContent(i) - f_N->Eval(h_Qs->GetBinCenter(i));
			//~ if (dQ < 0) dQ = 0;
			h_Qsig->SetBinContent(i, dQ);
		}
		else h_Qsig->SetBinContent(i,0);
	h_Qsig->ResetStats();
	h_Qsig->Draw("same hist");
	
	Q0 = f_N->GetParameter("Q0");
	sig0 = f_N->GetParameter("sig0");
	
	// A = f_N->GetParameter("A");
	alpha = f_N->GetParameter("alpha");
	w = f_N->GetParameter("w");
	f_N->SetParameter("alpha",alpha);
	Qshift = Q0 + w/alpha;	
	
	//gauss bunch parameters -- Q1, sig1, mu
	TF1 *f_sig = new TF1("f_sig",Funsig,Qv,Qend,5);
	f_sig->SetLineColor(3);
	f_sig->SetParNames("Qshift","Q1","sig1","mu","A");
	Q1 = h_Qsig->GetBinCenter(h_Qsig->GetMaximumBin())-Qshift;
	sig1 = (Qp - Qv);
	
	double mu1, mu2;
	TGraph *g_Qsig = new TGraph();
	g_Qsig->SetMarkerStyle(8);
	g_Qsig->SetMarkerColor(3);
	g_Qsig->SetMarkerSize(1);
	
	vector<int> v_Qsig;
	
	for(int i = i_Qv; i < i_end; ++i)
		if (h_Qsig->GetBinContent(i-1) < h_Qsig->GetBinContent(i) && h_Qsig->GetBinContent(i) >= h_Qsig->GetBinContent(i+1)){
			v_Qsig.push_back(i);
			g_Qsig->AddPoint(h_Qsig->GetBinCenter(i), h_Qsig->GetBinContent(i));
		}
	int n = v_Qsig.size();
	if (n > 2) {
		n = 2;
		Q1 = g_Qsig->GetPointX(1) - g_Qsig->GetPointX(0);
	}
	else Q1 = g_Qsig->GetPointX(0) - Qshift;
	Qshift = g_Qsig->GetPointX(0) - Q1;
	A = h_Qsig->GetMaximum();
	mu = TMath::Power(  g_Qsig->GetPointY(n-1)/A*TMath::Factorial(n)*TMath::Sqrt(TMath::TwoPi()*n) ,1./n);
	g_Qsig->Draw("P same");
	
	mu1 = h_Qsig->GetMaximum()*sig1*2.5/A;
	
	f_sig->SetParameters(Qshift,Q1,sig1,mu,A);
	A = A*g_Qsig->GetPointY(0)/f_sig->Eval(g_Qsig->GetPointX(0));
	f_sig->SetParameter("A",A);
	cout<<endl<<"SIGNAL PARAMETER ESTIMATION FIT\n_______________________________\n";
	h_Qsig->Fit(f_sig,"","",Qv,Qend);
	// Q0 = f_sig->GetParameter("Qshift");
	Qshift = f_sig->GetParameter("Qshift");

	Q1 = f_sig->GetParameter("Q1");
	sig1 = f_sig->GetParameter("sig1");
	mu = f_sig->GetParameter("mu");
	A = f_sig->GetParameter("A");
	

	cout<<"h_s parameters: \n";
	for(int i = 0; i<5;++i)cout<<f_sig->GetParName(i)<<" = "<<f_sig->GetParameter(i)<<endl;
	
	//~ h_Qsig->Fit(f_sig,"R");
	//~ TCanvas * cs = new TCanvas();
	//~ h_Qsig->Fit(f_sig,"R");
	f_sig->Draw("same");
    int i = 0;

    for (double par : {Q0, Q1, Qv, Qp, Qend, Qshift, sig0, sig1, w, alpha, mu, A}){ 
		*pars[i++] = par;
	}
	

}
void FindExtrInfl(TH1D * h_Qs, TH1D * h_dQ, TGraph *g[4], vector<int> v[4], int &i_Q0, int &i_end, double *pars[12]){ 
    TGraph *g_Qp = g[0], *g_Qv = g[1], *g_dQp = g[2], *g_dQv = g[3];
    vector<int> *v_Qp = &v[0], *v_Qv = &v[1], *v_dQp = &v[2], *v_dQv = &v[3];
    
	h_Qs->Smooth(5);
	h_Qs->ResetStats();

	//go to Q0 -> pedestal gauss peak
	// while ( h_Qs->GetBinContent(i) < h_Qs->GetBinContent(i+1)) ++i;
	i_Q0 = 1;
	cout<<"i_Q0 = "<<i_Q0<<endl;
    int i = 2;

	//finding local extremes in spectra
	double thr = h_Qs->GetMaximum()*1e-2;
	while (h_Qs->GetBinContent(++i) > thr) {
		if (h_Qs->GetBinContent(i-1) > h_Qs->GetBinContent(i) && h_Qs->GetBinContent(i) <= h_Qs->GetBinContent(i+1)) v_Qv->push_back(i);
		if (h_Qs->GetBinContent(i-1) < h_Qs->GetBinContent(i) && h_Qs->GetBinContent(i) >= h_Qs->GetBinContent(i+1)) v_Qp->push_back(i);
	}
	i_end = i;
	cout<<"i_end = "<<i<<endl;
	//putting local extremes to vec and to TGraph objects
	if (v_Qv->empty()) v_Qv->push_back(i_Q0+1);
	if (v_Qv->begin() > v_Qp->begin()) v_Qv->at(0) = i_Q0+1;
	for(int p, j = 0; j < v_Qp->size(); ++j){
		p = v_Qp->at(j);
		g_Qp->AddPoint(h_Qs->GetBinCenter(p), h_Qs->GetBinContent(p));
	}
    for(int v, j = 0; j < v_Qv->size(); ++j){
		v = v_Qv->at(j);
		g_Qv->AddPoint(h_Qs->GetBinCenter(v), h_Qs->GetBinContent(v));
	}

    //making derivative of spectra
	double dQ;
	for (int j = 1; j < i_end; j++ ) {
		dQ = h_dQ->GetBinContent(j) - h_dQ->GetBinContent(j+1);
		h_dQ->SetBinContent(j,dQ);
	}
	h_dQ->Smooth(10);
	h_dQ->ResetStats();
	
    //finding inflection points
	i = 0;
	while ( h_Qs->GetBinContent(i) < h_Qs->GetBinContent(i+1)) ++i;

	while (++i < i_end) {
		if (h_dQ->GetBinContent(i-1) >= h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) < h_dQ->GetBinContent(i+1)) v_dQv->push_back(i);
		if (h_dQ->GetBinContent(i-1) <= h_dQ->GetBinContent(i) && h_dQ->GetBinContent(i) > h_dQ->GetBinContent(i+1)) v_dQp->push_back(i);
	}
	for(int v,p,j = 0; j < min(v_dQv->size(), v_dQp->size());++j){
		v = v_dQv->at(j);
		p = v_dQp->at(j);
		g_dQv->AddPoint(h_dQ->GetBinCenter(v), h_dQ->GetBinContent(v));
		g_dQp->AddPoint(h_dQ->GetBinCenter(p), h_dQ->GetBinContent(p));
	}
	
	int i_Qv, i_Qp;
	if (v_Qv->size() > 0 && v_Qp->size() > 0){
		i_Qv = v_Qv->at(0);
		i_Qp = v_Qp->at(0);
	}else  //take it from inflection points
	{
		i_Qv = v_dQv->at(0);
		i_Qp = v_dQp->at(0);
	}

	// parameters needed for Fit function Init and for it's parameter estimation
	*pars[0] = h_Qs->GetBinCenter(i_Q0);
	*pars[2] = h_Qs->GetBinCenter(i_Qv);
	*pars[3] = h_Qs->GetBinCenter(i_Qp);
	*pars[4] = h_Qs->GetBinCenter(i_end);

    //plot config
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
}
void DrawFrame(TH1D * h_Q0, int i_Q0, int i_end){
	TGraph *g_W = new TGraph();
	g_W->SetLineColor(1);
	g_W->SetLineWidth(1);
	g_W->SetLineStyle(2);
	
	g_W->AddPoint(h_Q0->GetBinCenter(i_Q0),h_Q0->GetBinContent(i_end));
	g_W->AddPoint(h_Q0->GetBinCenter(i_Q0),h_Q0->GetBinContent(i_Q0));
	g_W->AddPoint(h_Q0->GetBinCenter(i_end),h_Q0->GetBinContent(i_Q0));
	g_W->AddPoint(h_Q0->GetBinCenter(i_end),h_Q0->GetBinContent(i_end));
	g_W->AddPoint(h_Q0->GetBinCenter(i_Q0),h_Q0->GetBinContent(i_end));
	g_W->Draw("same");
}
TF1 *InitFun(double *pars[12]){
	double Q0 = *pars[0];
	// double Q1 = *pars[1];
	double Qv = *pars[2];
	double Qp = *pars[3];
	double Qend = *pars[4];
	// double Qshift = *pars[5];
	// double sig0 = *pars[6];
	// double sig1 = *pars[7];
	// double w = *pars[8];
	// double alpha = *pars[9];
	// double mu = *pars[10];
	// double A = *pars[11];

	TF1 *f = new TF1("f_Q",Fun, Q0, Qend,9);
	f->SetLineColor(1);
	f->SetLineWidth(2);
	f->SetParNames("Q0","sig0","Q1","sig1","Qshift","w","alpha","mu","A");
	f->FixParameter(0, 0);
	f->FixParameter(1,0);
	// "Q0","sig0","Q1","sig1","Qshift","w","alpha","mu","A"

	return f;
}
array<double,9> pars2FitPars(double *pars[12]){
	array<double,9> pars_fit;

	pars_fit[0] = *pars[0];
	pars_fit[1] = *pars[6];
	pars_fit[2] = *pars[1];
	pars_fit[3] = *pars[7];
	pars_fit[4] = *pars[5];
	pars_fit[5] = *pars[8];
	pars_fit[6] = *pars[9];
	pars_fit[7] = *pars[10];
	pars_fit[8] = *pars[11];
	
	return pars_fit;
}
array<double,9> Pars2FitPars(array<double,12> Pars){
	array<double,9> FitPars;

	FitPars[0] = Pars[0];
	FitPars[1] = Pars[6];
	FitPars[2] = Pars[1];
	FitPars[3] = Pars[7];
	FitPars[4] = Pars[5];
	FitPars[5] = Pars[8];
	FitPars[6] = Pars[9];
	FitPars[7] = Pars[10];
	FitPars[8] = Pars[11];
	
	return FitPars;
}
void SetPars(TF1 *f, array<double,12> Pars, bool printout = false){
	array<double,9> FitPars ;
	FitPars = Pars2FitPars(Pars);
	f->SetParameters(FitPars[0], FitPars[1], FitPars[2], FitPars[3], FitPars[4], FitPars[5], FitPars[6], FitPars[7], FitPars[8]);
	cout<<"\n\tSetting parameters\n\n";
	
	if (printout) for (int i = 0; i < Nfpars; ++i) cout<<fparnames[i]<<" = "<<f->GetParameter(i)<<endl;
}
void SetParLimFixs(TF1 *f, array<double,12> Pars, array<bool,9> fix = {1,1,0,0,0,0,0,0,0}, array<bool,9> leave = {0,0,0,0,0,0,0,0,0}){
	double Q0 = Pars[0];
	double Q1 = Pars[1];
	double Qv = Pars[2];
	double Qp = Pars[3];
	double Qend = Pars[4];
	double Qshift = Pars[5];
	double sig0 = Pars[6];
	double sig1 = Pars[7];
	double w = Pars[8];
	double alpha = Pars[9];
	double mu = Pars[10];
	double A = Pars[11];	
	array<array<double,2>,9> l;

	int i = 0;
	l[i++] = {0, abs(Q0)*2}; //Q0
	l[i++] = {sig0/4, sig0*4}; //sig0
	l[i++] = {Q1/4, Q1*4}; //Q1
	l[i++] = {(Qp-Qv)*0.1, (Qp-Qv)*10}; //sig0
	l[i++] = {-Qp, Qp}; //Qshift
	l[i++] = {0, 1}; //w
	l[i++] = {0, Qend*2}; //aplha
	l[i++] = {mu/100, mu*100}; //mu
	l[i++] = {A/100, A*100}; //A

	for (int i = 0; i < 9; ++i){
		if (leave[i]) f->SetParLimits(i, 0, 0);
		else{
			if (!fix[i]) f->SetParLimits(i, l[i][0], l[i][1]);
			else f->FixParameter(i,0);
		}
	}
}



void AddPars(vector<array<double,12>> &v_pars, double *pars[12]){
	array<double,12> ar_pars;

	for (int i = 0; i<12; ++i){

		ar_pars[i] = *pars[i];
	}
	v_pars.push_back(ar_pars);
}

TFitResult * FitCompare(TH1D * h_Q, TF1 * f_Q, vector<array<double,12>> v_Pars){
	int i_best = 0, i = 0;
	double chi = 1e100, Qv = v_Pars[0][2], Qend = v_Pars[0][4];
	TFitResultPtr r;
	TFitResult * r_best = new TFitResult();
	for( auto Pars: v_Pars ){
		SetParLimFixs(f_Q, Pars);
		SetPars(f_Q, Pars, true);
		r = h_Q->Fit(f_Q, "S N", "", Qv, Qend);
		
		
		if (chi > r->Chi2()) {
			chi = r->Chi2();
			*r_best = *r;
			i_best = i;
		}
	
		cout<<"Fit #"<<i<<"\n\tChi^2 = "<<f_Q->GetChisquare()<<endl;
		++i;
	}
	for(int j = 0; j < Nfpars; ++j){
		f_Q->SetParameter(j,r_best->Parameter(j));
	}
	cout<<"blafitres chi2 =" <<r_best->Chi2()<<endl;

	return r_best;
}
void FitHist( TH1D * h_Q, string DataName = "") {
	string Img1Name = "img/fits/" + DataName.substr(DataName.rfind("/")+1) + "extr.png";
	string Img2Name = "img/fits/" + DataName.substr(DataName.rfind("/")+1) + "est.png";
	string Img3Name = "img/fits/" + DataName.substr(DataName.rfind("/")+1) + "fit.png";
	string Img4Name = DataName + "/fit_" + DataName.substr(DataName.rfind("/")+1) + ".png";
	string JsonName = DataName + "/fit_" + DataName.substr(DataName.rfind("/")+1) + ".json";
	double Q0, Q1, Qv, Qp, Qend, Qshift, sig0, sig1, w, alpha, mu, A;
    double *pars[12] = {&Q0, &Q1, &Qv, &Qp, &Qend, &Qshift, &sig0, &sig1, &w, &alpha, &mu, &A};
	int i_Q0, i_Q1, i_Qv, i_Qp, i_end, i_dQ0;

    vector<int> v_all[4];
	vector<int> *v_Qp = &v_all[0], *v_Qv = &v_all[1], *v_dQp = &v_all[2], *v_dQv = &v_all[3];
    TGraph * g_all[4] = {new TGraph(), new TGraph(), new TGraph(), new TGraph()};
	TGraph *g_Qp = g_all[0];
	TGraph *g_Qv = g_all[1];
	TGraph *g_dQp = g_all[2];
	TGraph *g_dQv = g_all[3];
	
	TH1D * h_Q0 = (TH1D*) h_Q->Clone("h_Q0");
	// h_Q0->Scale(1/h_Q0->GetEntries());
	TH1D * h_Qs = (TH1D*) h_Q->Clone("h_Qs");
    TH1D * h_dQ = (TH1D*) h_Q->Clone("h_dQ");

    FindExtrInfl(h_Qs, h_dQ, g_all, v_all, i_Q0, i_end, pars);	

	TCanvas * c_ei = new TCanvas("c","Extremes and inflection points",1000,1400);
	
	c_ei->Divide(1,2);
		
	c_ei->cd(1);
	h_Q0->Draw("HIST");
	h_Q0->GetXaxis()->SetRange(0, i_end*1.2);

	DrawFrame(h_Q0, i_Q0, i_end);
	
	
	
	h_Qs->SetLineColor(2);
	h_Qs->Draw("HIST same");
	if (g_all[0]->GetN() > 0 && g_all[1]->GetN()>0){
		g_all[0]->Draw("P same");
		g_all[1]->Draw("P same");
	}
	c_ei->cd(2);
	h_dQ->Draw("HIST");
	
	g_all[2]->Draw("P same");
	g_all[3]->Draw("P same");
	

	// Fit part
	TCanvas * c_est = new TCanvas("c_est", "Parameter estiamation results");
	gStyle->SetOptStat(0);
	c_est->Divide(1,3);

	TF1 *f_Q = InitFun(pars); // fit function
	vector<array<double,9>> v_FitPars; // vector for estimated parameters from different methods
	vector<array<double,12>> v_Pars;
	c_est->cd(1);
	ParEst1(h_Q0, h_Qs, i_Q0, i_end, pars, v_all);
	AddPars(v_Pars, pars);

	c_est->cd(2);
	ParEst2(h_Q0, h_Qs, i_Q0, i_end, pars, v_all);
	AddPars(v_Pars, pars);	
	
	int i = 0;
	for (auto Pars: v_Pars){
		cout<<"\n\tPars #"<<++i<<":\n\n";
		for (int j = 0; j < Npars; ++j) cout<<parnames[j]<<" = "<<Pars[j]<<endl;

	}
	

	TFitResult r = *FitCompare(h_Q0, f_Q, v_Pars);
	TCanvas * c_fit = new TCanvas("c_fit", "Fit results");
	h_Q0->SetMarkerStyle(7);
	h_Q0->SetMarkerColor(1);
	h_Q0->SetMarkerSize(2);
	h_Q0->Draw("H P");
	f_Q->Draw("same");

	cout<<"\tChi^2/NDF = "<<r.Chi2()/r.Ndf()<<endl;

	// outputs to files
	TLegend * leg_fit = new TLegend(0.6, 0.5, 0.9, 0.9);
	Drawfun(f_Q, leg_fit);
   	leg_fit->AddEntry(h_Q0,"PMT Spectre","L P");
   	leg_fit->AddEntry(f_Q,"FIT","L");
	leg_fit->AddEntry("n",((string) "Chi^{2}/NDF = " + to_string(r.Chi2()/r.Ndf())).data(),"n");
   	leg_fit->Draw();

	cout<<"\nimgname: "<<Img1Name<<endl;
	c_ei->Print(Img1Name.data());
	c_est->Print(Img2Name.data());
	c_fit->Print(Img3Name.data());
	c_fit->Print(Img4Name.data());

	json j_fit;
	for (int i = 0; i < r.NPar(); ++i){
		j_fit[r.GetParameterName(i)]["val"] = r.Parameter(i);
		j_fit[r.GetParameterName(i)]["fix"] = r.IsParameterFixed(i);
		if (!(r.IsParameterFixed(i))) j_fit[r.GetParameterName(i)]["err"] = r.ParError(i);
	}
	cout<<j_fit<<endl;
	ofstream of_fit(JsonName.data());
	of_fit << setw(4)<<j_fit<<endl;
	of_fit.close();

}
void FitHistTB(string DataName){
	string hDataName = DataName + "/OscLecHists.root";
	// TFile * f = TFile::Open("data/dr_smalwfms/TSA.root");
	// TTree * t = (TTree*) f->Get("TSA");
	// int ientry = 20;
	// TH1D * h_Q = new TH1D();
	// t->SetBranchAddress("hQ",&h_Q);
	// t->GetEntry(ientry);
	
	
	TFile * f = TFile::Open(hDataName.data());
	TH1D * h_Q = (TH1D*) f->Get("Q");
	

	FitHist(h_Q, DataName);
	
	
		
	//~ FitHist(h);	
}
void MuPlot(){
	TCanvas *c1 = new TCanvas();
	c1->Divide(1,2);
	
	double Q1 = 1;
	double sig1 = 0.3;


	c1->cd(1);
	TF2 *f_mu = new TF2("Funmu",Fun_mu,0,10,0,4,2);
	f_mu->SetParameters(Q1,sig1);
	f_mu->Draw("surf1");
	
	c1->cd(2);
	TF1 *f = new TF1("f_mu2",Funsig,0,10,5);
	f->SetParameters(0,Q1,sig1,2,1);
	f->Draw();

}
