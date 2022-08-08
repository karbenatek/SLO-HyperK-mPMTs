#ifndef ascii2root_h
#define ascii2root_h

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <filesystem>
#include <TFile.h>
#include <TTree.h>
#include "ascii2root.h"
int n_digit(int number) {
   int count = 0;
   while(number != 0) {
      number = number / 10;
      count++;
   }
   return count;
}
void GetFileInfo(ifstream & file, int & fpos_Wfm, int & fpos_h, int & n_Wfm, int & n_Seg ){
	char line[64];
	file.getline(line,64); 
	string l = line;
	
	while ( l.compare("Time,Ampl")< 0 && !file.eof()){
		if (l.find("Segments,") <= l.length()) {
			n_Wfm = stoi(l.substr(9));
			//cout<<l.substr(9 + n_digit(n_Wfm) + 13)<<endl;
			n_Seg = stoi(l.substr(9 + n_digit(n_Wfm) + 13)) + 1000;
		}
		if (fpos_h <= 0) fpos_h = -file.tellg();
		file.getline(line,64);
		l = line;
		fpos_Wfm = file.tellg();
		if (fpos_h <= 0){
			 if (line[0] == '#') fpos_h = -fpos_h;	
		}
	}
	file.seekg(0,file.beg);
	
	
}
void ParseWfm(ifstream & file, TTree * fChain, int & fpos_Wfm, int & fpos_h, int & i_Wfm, int & i_Seg, string & Date,	double & TimeSinceSeg1, vector<double> & x, vector<float> & t){

	//some parsed values used in process
	double Ampl;//->x
	float Time = -1e10;//->t
	float p_Time;// previous
	string l;
	string::size_type sz;
	char line[64];
	//head process
	file.seekg(fpos_h,file.beg);
	file.getline(line,64);
	l = line;
	l = l.substr(1); //skips "#" character
	i_Seg = std::stoi(l,&sz);
	++i_Wfm;
	l = l.substr(sz+1); //date, time since segment1... date format = 01-Aug-2022 11:55:02 ... 20 chars
	Date = l.substr(0,20);
	l = l.substr(21);
	TimeSinceSeg1 = stod(l);
	fpos_h = file.tellg();
	//wfm process
	file.seekg(fpos_Wfm,file.beg);
	t.clear();
	x.clear();
	while(1) {
		p_Time = Time;
		fpos_Wfm = file.tellg();
		file.getline(line,64);
		if(file.eof()) return;
		l = line;
		Time = stof(l, &sz);
		l = l.substr(sz+1);
		Ampl = stod(l);		
		if (Time >= p_Time){
			t.push_back(Time);
			x.push_back(Ampl);
		}
		else{
			fChain->Fill();
			break;
		}	
	}
}


void ProcessFile(const string fPath, TFile * fOut, TTree * fChain, int & i_Wfm, int & i_Seg, string & Date, double & TimeSinceSeg1, vector<double> & x, vector<float> & t){
	//cout<<Fpath<<endl;

	ifstream file;
	file.open(fPath, std::fstream::in);
	//char line[64];
	//file.getline(line,64);
	

	
	
	
	int fpos_Wfm = 0;
	int fpos_h = 0;
	int n_Wfm;//# of Wfms
	int n_Seg;//# of samples in Wfm aka segment

	GetFileInfo(file, fpos_Wfm, fpos_h, n_Wfm, n_Seg);
	cout<<"Running - "<<fPath<<endl<<"#wfms = "<<n_Wfm<<endl;

	while (!file.eof()) {
		ParseWfm( file, fChain, fpos_Wfm, fpos_h, i_Wfm, i_Seg, Date, TimeSinceSeg1, x, t);
		if (i_Seg % 1000 == 0) cout<<i_Seg<<"/"<<n_Wfm<<endl;
	}

	file.close();
	cout<<"Done"<<endl;
}

void ProcessDirData(string directory,string suffix = ".txt"){
	TFile * fOut = new TFile("OscLecData.root","recreate");
	TTree * fChain = new TTree("OscLecData","Waveforms from LECROYWM806Zi-B");

	//parsed variables	
	int i_Wfm = 0; 
	int i_Seg; 
	string Date; 
	double TimeSinceSeg1; 
	vector<double> x; 
	vector<float> t;
	x.reserve(4000);
	t.reserve(4000);

	//def branches
	auto b_i_Wfm = fChain->Branch("i_Wfm", &i_Wfm);
	auto b_i_Seg = fChain->Branch("i_Seg", &i_Seg);
	auto b_Date = fChain->Branch("Date", &Date,21,0);
	auto b_TimeSinceSeg1 = fChain->Branch("TimeSinceSeg1",&TimeSinceSeg1);
	auto b_x = fChain->Branch("x",&x,4000,0);
	auto b_t = fChain->Branch("t",&t,4000,0);
	
	string fPath;
	const std::filesystem::path dirPath{directory};
	for (auto const& dir_entry : std::filesystem::directory_iterator{dirPath}){
		fPath = dir_entry.path();
		if( fPath.find(suffix,0) <= fPath.length() ) {
			ProcessFile(fPath, fOut, fChain, i_Wfm, i_Seg, Date, TimeSinceSeg1, x, t);
		}
	}
	fOut->Write();
	fOut->Close();
	return;	
}

#endif
