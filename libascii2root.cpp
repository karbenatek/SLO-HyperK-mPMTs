//functions to parse data from osciloscope
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;
#include <TFile.h>
#include <TTree.h>


int n_digit(int number) {
   int count = 0;
   while(number != 0) {
      number = number / 10;
      count++;
   }
   return count;
}
void GetFileInfo(ifstream & file, int & fpos_Wfm, int & fpos_h, int & n_Wfm, int & n_Seg){
	char line[64];
	file.getline(line,64); 
	string l = line;
	
	while ( l.compare("Time,Ampl") < 0 && !file.eof() ){
		if (l.find("Segments,") <= l.length()) {
			n_Wfm = stoi(l.substr(9));
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
void ParseWfm(ifstream & file, int & fpos_Wfm, int & fpos_h, int & i_Wfm, int & i_Seg, char Date[21],	double & TimeSinceSeg1, vector<double> & x, vector<float> & t){
	//some parsed values used in process
	double Ampl;//->x
	float Time = -1e10;//->t
	float p_Time;// previous
	string l;
	string::size_type sz;
	char line[64];
	// wfm head process
	file.seekg(fpos_h,file.beg);
	file.getline(line,64);
	l = line;
	l = l.substr(1); //skips "#" character
	i_Seg = std::stoi(l,&sz);
	++i_Wfm;
	l = l.substr(sz+1); //date, time since segment1... date format = 01-Aug-2022 11:55:02 ... 20 chars
	string sDate = l.substr(0,20);
	for (int i = 0; i < sDate.size(); ++i) Date[i] = sDate[i];
	l = l.substr(21);
	TimeSinceSeg1 = stod(l);
	fpos_h = file.tellg();
	//wfm data process
	file.seekg(fpos_Wfm,file.beg);
	t.clear();
	x.clear();
	int i = 0;
	while(1) {		
		p_Time = Time;
		fpos_Wfm = file.tellg();
		file.getline(line,64);
		l = line;
		if(file.eof() || l.empty()) return;
		Time = stof(l, &sz);
		l = l.substr(sz+1);
		Ampl = stod(l);
		
		if (Time >= p_Time){
			t.push_back(Time);
			x.push_back(Ampl);
			++i;
		}
		else break;
	}
}
void ProcessFile(const string fPath, TFile * fOut, TTree * fChain, int & i_Wfm, int & i_Seg, char Date[21], double & TimeSinceSeg1, vector<double> & x, vector<float> & t){
	ifstream file;
	file.open(fPath, std::fstream::in);

	int fpos_h = 0; //position for reading wfm data
	int fpos_Wfm = 0; //position for reading wfm data
	int n_Wfm;//# of Wfms
	int n_Seg;//# of samples in Wfm ~ segment

	GetFileInfo(file, fpos_Wfm, fpos_h, n_Wfm, n_Seg);
	cout<<"Processing - "<<fPath<<endl;
	while (true) {
		ParseWfm( file, fpos_Wfm, fpos_h, i_Wfm, i_Seg, Date, TimeSinceSeg1, x, t);
		fChain->Fill();
		if (i_Seg % 1000 == 0) cout<<i_Seg<<"/"<<n_Wfm<<endl;
		if (file.eof())break;
	}

	file.close();
}
void ProcessDirData(string DataDir, const char InitLetter = 'C', string Suffix = ".txt", string inCfg = "cfg/cfg.json"){
	// set output filename & config filename to be copied
	string fDataName = DataDir + "/OscLecData.root";
	string outCfg = DataDir + "/cfg.json";
	
	// create output file and copy config file from .
	TFile * fOut = new TFile(fDataName.data(),"recreate");
	TTree * fChain = new TTree("OscLecData","Waveforms from LECROYWM806Zi-B");
	fs::copy(inCfg, outCfg, fs::copy_options::skip_existing);

	// parsed variables
	int i_Wfm = 0; 
	int i_Seg; 
	char Date[22]; 
	double TimeSinceSeg1; 
	vector<double> x; 
	vector<float> t;

	// def branches
	fChain->Branch("i_Wfm", &i_Wfm);
	fChain->Branch("i_Seg", &i_Seg);
	fChain->Branch("Date", Date, "Date/C", 22);
	fChain->Branch("TimeSinceSeg1",&TimeSinceSeg1);
	fChain->Branch("x",&x);
	fChain->Branch("t",&t);
	
	// get files/folders in DataDir
	string fPath;
	set<fs::path> sorted_paths;
	for (auto &dir_entry : fs::directory_iterator(DataDir)){
		sorted_paths.insert(dir_entry.path());
	}
	
	// iterate throught the files/folders in DataDir
	for (auto &filename : sorted_paths ){
		fPath = filename.c_str();
		// process only those with correct suffix & first letter
			// example of osciloscope file output: C1Trace00000.txt
			// we want to process files starting with "C" and ".txt" suffix
		if( ( fPath.find(Suffix, 0) <= fPath.length() ) && ( fPath.substr( fPath.rfind('/')+1)[0] == InitLetter) ) {
			ProcessFile(fPath, fOut, fChain, i_Wfm, i_Seg, Date, TimeSinceSeg1, x, t);
		}
	}
	// write into output file and close it
	fOut->Write();
	fOut->Close();
	cout<<"Data were seccessfully writen to file "<<fDataName<<endl;	
}