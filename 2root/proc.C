#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include "ascii2root.h"

void proc(){
	//dr_data is a folder path with stored ascii files from Oscliloscope LECROYWM806Zi-B
	ProcessDirData("dr_data");
	//will produce file OscLecData.root with wfms from all files in specified folder
	
/* 
	use:
	* change input of ProcessDirData() func to desired directory
	* 
	* run
	* $ root proc.C
	* or
	* $ root proc.C+ ... will be compiled so it might run faster
	*/
}

