/************************************************************************/
/*
* ply2crestlineinput
*
*	Convert the ply format model into the format for crest line code
*
*	Author: Chao Wang
*	Contact: chao.wang3@utdallas.edu
*	Date: 1.7.2013
*/
/************************************************************************/

#include <iostream>
#include "PLYReader.h"

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 4)
	{
		cout<<"Usage: plyconverter --options input.ply outputfilename"<<endl;
		cout<<"options:"<<endl
			<<"--txt/-t: convert a ply file into a txt file (input of Detect-Crest-Line \"setCurvature\" code)"<<endl
			<<"--ply2/-p: convert a ply file into a ply2 file (input of Crest-Line-Viewer \"CrestViewer\" code)"<<endl
			<<"--gts/-g: convert a ply file into a gts (GNU Triangulated Surface) file (input of Smoothing \"smoother\" code)"<<endl;
		return EXIT_FAILURE;
	}

	char* inname = argv[2];
	char* outname = argv[3];
	PLYReader *plyrd = NULL;
	if ( strstr( inname, ".ply" ) != NULL )
	{
		plyrd = new PLYReader(inname);
		plyrd->printInfo();
	}
	else
	{
		cout<<"Unknown file format. Quiting"<<endl;
		return EXIT_FAILURE;
	}
	if ((strcmp(argv[1], "--txt") == 0) || ((strcmp(argv[1], "-t") == 0)) )
	{
		cout<<"Convert ply to txt..."<<endl;
		plyrd->writeCrestLineInputTxt(outname);
		cout<<"Done."<<endl;
	}
	else if ((strcmp(argv[1], "--ply2") == 0) || ((strcmp(argv[1], "-p") == 0)) )
	{
		cout<<"Saving ply to ply2..."<<endl;
		plyrd->writePLY2(outname);
		cout<<"Done."<<endl;
	}
	else if ((strcmp(argv[1], "--gts") == 0) || ((strcmp(argv[1], "-gts") == 0)) )
	{
		cout<<"Converting ply to gts..."<<endl;
		plyrd->writeGTS(outname);
		cout<<"Done."<<endl;
	}

	delete plyrd;
	return EXIT_SUCCESS;
}