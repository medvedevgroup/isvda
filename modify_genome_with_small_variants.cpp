/*
 * Author: Chen Sun(cxs1031@cse.psu.edu)
 *
 *
 * Modify: 9/15/2014 
 *  	   Using a fast algorithm that search all sv, not all site of genome, the running time is worth, haven't debug yet, the usable code is high_confident, which use old algorithm
 *			10/15/2014
 *			Already debuged, this is the up to date version
 * */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <map>
#include <ctime>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
using namespace std;

#define OLDLEN 400000000
#define NEWLEN 400000000
#define TESTLEN  2000

void dsptime()
{
 time_t nowtime;
 nowtime = time(NULL); //get int time number
 struct tm * ptm=localtime(&nowtime);  //get system time
 cout << ptm->tm_mon+1 << "/" << ptm->tm_mday << "/"<< ptm->tm_year+1900 << "," ;
 cout << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec <<" " << endl;
}

/*split function*/
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (! item.empty()){
			elems.push_back(item);
		}
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

string readRef(char* filename, char* & firstLine, int & refLength, int & lineLength)
{
	ifstream refFile;
	refFile.open(filename);
	//char * refSequence = (char*) calloc(OLDLEN , sizeof(char));
	string sequenceString="";
    if (!refFile.good())
    {
        cout << "Error in open reference file!" << endl;
    }
	//int lineNum  = 0;
	while(refFile != NULL){
        string string_line;
        getline(refFile, string_line, '\n');
		//lineNum ++;

        char* line=const_cast<char*>(string_line.c_str());
		if(line[0] == '>')
		{
		    cout << "firstLine" << line << endl;
		    strcpy(firstLine,line);
		    cout<<firstLine << endl;
			continue;
		}
		if (lineLength == 0) lineLength = strlen(line);
		//refSequence = strcat(refSequence, line);
		sequenceString += string_line;
		//cout << lineNum << endl;
	}
    refFile.close();
	//char * refSequence = const_cast<char*>(sequenceString.c_str());
    //refLength = strlen(refSequence);
	refLength = sequenceString.size();
	//cout << lineLength << endl;
	//cout << refSequence << endl;

	return sequenceString;
}

void writeRef(char* filename, char* firstLine, string refSequence, int refLength, int lineLength)
{
    ofstream outFile;
    outFile.open(filename);
    outFile << firstLine << endl;
    int markLength = 0;
    //char * outSequence = (char*) calloc(lineLength+1, sizeof(char));
    string outSequence = "";
	while(markLength + lineLength < refLength)
    {
        //strncpy(outSequence, refSequence+markLength, lineLength);
        outSequence = refSequence.substr(markLength, lineLength);
		outFile << outSequence << endl;
        markLength += lineLength;
    }
    if(markLength < refLength)
    {
        //char* tempSequence = (char*) calloc(refLength - markLength+1, sizeof(char));
        //strncpy(tempSequence, refSequence+markLength, refLength-markLength);
        string tempSequence = refSequence.substr(markLength, refLength-markLength);
		outFile << tempSequence << endl ;
    }
    outFile.close();
}

bool checkInfo(string info){
	//check from info
	vector<string> infoVector = split(info, ';');
	for(int i = 0; i < infoVector.size(); i++){
		string section = infoVector[i];
		vector<string> sectionSplit = split(section, '=');
		if(sectionSplit.size() == 2){
			if(sectionSplit[0] == "AB"){
				if(sectionSplit[1] != "0") return false;
				//else cout << "AB= " << sectionSplit[1] << endl;
			}
			if(sectionSplit[0] == "AF"){
				if(sectionSplit[1] != "1") return false;
				//else cout << "AF= " << sectionSplit[1] << endl;
			}
		}else{
			cout << "Warning: Error when check vcf info";
			for(int i = 0; i < sectionSplit.size(); i++){
				cout << sectionSplit[i] << endl;
			}
		}
	}
	return true;
}

void readSnp(char* filename, map<int, string> & position_ref, map<int, string> & position_alt, char* highConfidentFilename, char* lowConfidentFilename)
{
	ofstream highConfidentFile;
	highConfidentFile.open(highConfidentFilename);
	ofstream lowConfidentFile;
	lowConfidentFile.open(lowConfidentFilename);
	//outFile << "position" << "\t" << "ref" << "\t" << "alt" << "\t" << "abScore" << endl;
	ifstream snpFile;
	snpFile.open(filename);
	if (!snpFile.good())
	{
		cout << "Error: Fail to open snp file " << filename << endl;
	}
	char* columns[31];
	while (snpFile != NULL)
	{
		int section_index = 0;
		string string_line;
		getline(snpFile, string_line, '\n');
		char* line = const_cast<char*>(string_line.c_str());
		string holdLine = line;
		if (line[0] == '#')
		{
			highConfidentFile << line << endl;
			lowConfidentFile << line << endl;
			continue;
		}
		//cout << "breakpoint 1.1" << endl;
		char* section = strtok(line, "\t");
		while (section != NULL)
		{
			columns[section_index] = section;
			section = strtok(NULL, "\t");
			section_index++;
		}
		//cout << "breakpoint 1.2 " << section_index <<endl;
		//cout << line << endl;
		if (section_index == 0){ continue; }
		if (section_index < 5)
		{
			cout << "Error: vcf format error, please check your vcf file" << endl;
			return;
		}
		int position = atoi(columns[1]);
		position--; // vcf is 1 based
		//position *= -1;
		string ref = columns[3];
		string alt = columns[4];
		int qual = atoi(columns[5]);
		string info = columns[7];
		if(qual <= 20){
			lowConfidentFile << holdLine << endl;
			continue;
		}
		//split alt
		//cout << alt << "\t" << info << endl;
		char* alt_charStar = const_cast<char*>(alt.c_str());
		char* alleles[20];
		int allele_index = 0;
		char* alt_section = strtok(alt_charStar, ",");
		while (alt_section != NULL)
		{
			alleles[allele_index] = alt_section;
			alt_section = strtok(NULL, "\t");
			allele_index++;
		}
		if (allele_index > 1) {
			lowConfidentFile << holdLine << endl;
			continue;
		}
		
		if(!checkInfo(info)) {
			//cout << "Error: detecting low confident snp";
			lowConfidentFile << holdLine << endl;
			continue;
		}

		//split info
//		char* info_charStar = const_cast<char*>(info.c_str());
//		char* infoSet[200];
//		int info_index = 0;
//		// in reality, we should traverse all info, but here we know the first is AB
//		char* info_section = strtok(info_charStar, ";");
//		while (info_section != NULL)
//		{
//			infoSet[info_index] = info_section;
//			info_section = strtok(NULL, ";");
//			info_index++;
//		}
//
//		if (info_index < 4)
//		{
//			cout << "Error: not enough info sections." << endl;
//			continue;
//		}
//
//		char* singleInfo = infoSet[0];
//		int abScore = -1;
//		//cout << singleInfo << endl;
//		if (singleInfo[0] == 'A' && singleInfo[1] == 'B' && singleInfo[2] == '='){
//			char* abSign = strtok(singleInfo, "=");
//			char* abScore_charStar = strtok(NULL, "=");
//			abScore = atoi(abScore_charStar);
//			if (abScore != 0) continue;
//		}
//
//		char* afInfo = infoSet[3];
//		int afScore = -1;
//		if (afInfo[0] == 'A' && afInfo[1] == 'F' && afInfo[2] == '='){
//			char* afSign = strtok(afInfo, "=");
//			char* afScore_charStar = strtok(NULL, "=");
//			afScore = atoi(afScore_charStar);
//			if (afScore != 1) continue;
//		}
		//outFile << position << "\t" << ref << "\t" << alt << endl;
		highConfidentFile << holdLine << endl;
		position_ref.insert(pair<int, string>(position, ref));
		position_alt.insert(pair<int, string>(position, alleles[0])); //alt_section will be NULL finally
	}
	snpFile.close();
}

// this function no longer used
void readSnp_old_version(char* filename, map<int, string> & position_ref, map<int, string> & position_alt)
{
	ofstream outFile;
	outFile.open("qualifiedSNP.txt");
	outFile << "position" << "\t" << "ref" << "\t" << "alt" << "\t" << "abScore" << endl;
	ifstream snpFile;
	snpFile.open(filename);
	if (!snpFile.good())
	{
		cout << "Error: Fail to open snp file " << filename << endl;
	}
	char* columns[31];
	while (snpFile != NULL)
	{
		int section_index = 0;
		string string_line;
		getline(snpFile, string_line, '\n');
		char* line = const_cast<char*>(string_line.c_str());
		if (line[0] == '#')
		{
			continue;
		}
		//cout << "breakpoint 1.1" << endl;
		char* section = strtok(line, "\t");
		while (section != NULL)
		{
			columns[section_index] = section;
			section = strtok(NULL, "\t");
			section_index++;
		}
		//cout << "breakpoint 1.2 " << section_index <<endl;
		//cout << line << endl;
		if (section_index == 0){ continue; }
		if (section_index < 5)
		{
			cout << "Error: vcf format error, please check your vcf file" << endl;
			return;
		}
		int position = atoi(columns[1]);
		position--; // vcf is 1 based
		//position *= -1;
		string ref = columns[3];
		string alt = columns[4];
		string info = columns[7];
		//split alt
		//cout << alt << "\t" << info << endl;
		char* alt_charStar = const_cast<char*>(alt.c_str());
		char* alleles[20];
		int allele_index = 0;
		char* alt_section = strtok(alt_charStar, ",");
		while (alt_section != NULL)
		{
			alleles[allele_index] = alt_section;
			alt_section = strtok(NULL, "\t");
			allele_index++;
		}
		if (allele_index > 1) continue;

		//split info
		char* info_charStar = const_cast<char*>(info.c_str());
		char* infoSet[200];
		int info_index = 0;
		// in reality, we should traverse all info, but here we know the first is AB
		char* info_section = strtok(info_charStar, ";");
		while (info_section != NULL)
		{
			infoSet[info_index] = info_section;
			info_section = strtok(NULL, ";");
			info_index++;
		}

		if (info_index < 4)
		{
			cout << "Error: not enough info sections." << endl;
			continue;
		}

		char* singleInfo = infoSet[0];
		int abScore = -1;
		//cout << singleInfo << endl;
		if (singleInfo[0] == 'A' && singleInfo[1] == 'B' && singleInfo[2] == '='){
			char* abSign = strtok(singleInfo, "=");
			char* abScore_charStar = strtok(NULL, "=");
			abScore = atoi(abScore_charStar);
			if (abScore != 0) continue;
		}

		char* afInfo = infoSet[3];
		int afScore = -1;
		if (afInfo[0] == 'A' && afInfo[1] == 'F' && afInfo[2] == '='){
			char* afSign = strtok(afInfo, "=");
			char* afScore_charStar = strtok(NULL, "=");
			afScore = atoi(afScore_charStar);
			if (afScore != 1) continue;
		}
		outFile << position << "\t" << ref << "\t" << alt << "\t" << afScore << endl;
		position_ref.insert(pair<int, string>(position, ref));
		position_alt.insert(pair<int, string>(position, alleles[0])); //alt_section will be NULL finally
	}
	snpFile.close();
}

map<int, int> read_del_vcf(char* filename){
	map<int, int> delMap_start_end;

	ifstream vcfFile;
	vcfFile.open(filename);
	if (!vcfFile.good())
	{
		cout << "Error: Fail to open snp file " << filename << endl;
	}
	//char* columns[31];
	while (vcfFile != NULL)
	{
		int section_index = 0;
		string string_line;
		getline(vcfFile, string_line, '\n');
		//cout << string_line << endl;
		//const string conservedLine = string_line;
		char* line = const_cast<char*>(string_line.c_str());

		if (line[0] == '#')
		{
			//headInfo.push_back(string_line);
			continue;
		}

		vector<string> columns = split(string_line, '\t');
		if (columns.size() == 0){
			cout << "Error: no section in a line." << endl;
			 continue;
		 }
		if (columns.size() < 7)
		{
			cout << "Error: vcf format error, please check your vcf file" << endl;
			continue;
		}
		//cout << conservedLine << endl;
		int startPosition = atoi(columns[1].c_str());
		startPosition--; // vcf is 1 based
		//position *= -1;
		string quality = columns[6];
		if (quality.compare("PASS") != 0){
			//cout << quality << endl;
			continue;
		}
		string info = columns[7];

		vector<string> infoSet = split(info, ';');

		if (infoSet.size() < 4)
		{
			cout << "Error: not enough info sections." << endl;
			continue;
		}

		string singleInfo = infoSet[6];
		int endPosition = -1;
		int length = atoi(infoSet[7].c_str());
		if (singleInfo.substr(0,3) == "END"){
			vector<string> endInfoSet = split(singleInfo, '=');
			endPosition = atoi(endInfoSet[1].c_str());
			endPosition -= 2;
		}
		else{
			cout << "Error in reading end position of SV." << singleInfo.substr(0,3) << endl;
		}
		delMap_start_end.insert(pair<int, int>(startPosition, endPosition));
		//delMap_start_info.insert(pair<int, string>(startPosition, string_line)); //alt_section will be NULL finally
	}
	vcfFile.close();
	cout << "After Read: " << delMap_start_end.size() << endl;
	return delMap_start_end;
}


char* repairGenome(string refSequence, int refLength, map<int, string> position_ref, map<int, string> position_alt)
{
	char* donorSequence = (char*)calloc(NEWLEN, sizeof(char));
	int j = 0;
	for(int i = 0; i < refLength; i++,j++)
	{
		//cout << i << endl;
		if(j >= NEWLEN)
		{
			cout << "Error: donor genome index out of range" << endl;
		}
		if (position_ref.count(i)){
			//cout << i << " " << position_ref[i] << " " << position_alt[i] << endl;
			int refNum = strlen(position_ref[i].c_str());
			int altNum = strlen(position_alt[i].c_str());
			if (refNum == 1 && altNum == 1){
				donorSequence[j] = position_alt[i][0];
				continue;
			}
			donorSequence = strcat(donorSequence, position_alt[i].c_str());
			if (refNum > 1)
			{
				i += (refNum - 1);
			}
			if (altNum > 1)
			{
				j += (altNum - 1);
			}
		}else
		{
			donorSequence[j] = refSequence[i];
		}
	}
	cout << "about to return donor genome" << endl;
	return donorSequence;
}


char* quickRepairing(string refSequence, int refLength, map<int, string> position_ref, map<int, string>position_alt){
	char* donorSequence = (char*)calloc(NEWLEN, sizeof(char));
	int startPosition = 0;
	int endPosition = refSequence.size() - 1;
	//cout << "quick repairing, INDEL num: " << position_ref.size() << endl;
	//int index = 0;
	for(map<int, string>::iterator it = position_ref.begin(); it != position_ref.end(); it++){
		//index ++ ;
		//cout << index << endl;
		int position = it->first;
		int refNum = strlen(position_ref[position].c_str());
		int altNum = strlen(position_alt[position].c_str());
		//cout << position << "\t" << position_ref[position] << "\t" << position_alt[position] << endl;
		string subString = refSequence.substr(startPosition, position - startPosition);
		//cout << startPosition << "\t" << position-startPosition << "\t" << subString << endl;
		char* subCharStar = const_cast<char*>(subString.c_str());
		strcat(donorSequence, subCharStar);
		//if (refNum == 1 && altNum == 1)
		//{
		//	donorSequence[j] = position_alt[position][0];
		//}else{
		strcat(donorSequence, position_alt[position].c_str());
		//}
		startPosition = position + refNum;
	}
	if(startPosition < endPosition)
	{
		string subString = refSequence.substr(startPosition, endPosition - startPosition+1);
		//cout << startPosition << "\t" << endPosition-startPosition+1 << "\t" << subString << endl;
		strcat(donorSequence, subString.c_str());
	}
	return donorSequence;

}

// del info here is 0-based
char* chopDelSection(string refSequence, map<int, int> delMap_start_end){
	char* donorSequence = (char*)calloc(NEWLEN, sizeof(char));
	int startPosition = 0;
	int endPosition = refSequence.size() - 1;
	for(map<int, int>::iterator it = delMap_start_end.begin(); it != delMap_start_end.end(); it++){
		int delStart = it->first;
		int delEnd = it->second;
		if (delStart != 0){
			string subString = refSequence.substr(startPosition, delStart - startPosition);
			char* subCharStar = const_cast<char*>(subString.c_str());
			strcat(donorSequence, subCharStar);
		}
		startPosition = delEnd+1;
	}
	if(startPosition < endPosition){
		string subString = refSequence.substr(startPosition, endPosition - startPosition+1);
		strcat(donorSequence, subString.c_str());
	}
	return donorSequence;
}

//in the future, we may write a generous function to deal with 
/*
string getDonor(string refSequence, int refLength, int & donorLength, map<int, string> position_ref, map<int, string> position_alt)
{
	string donorSequence = "";
	//cout << refSequence.size();
	for(map<int, string>::iterator it = position_ref.begin(); it != position_ref.end(); it++)
	{
		int position = it->first;
		int real_position = position * -1;
		string ref = position_ref[position];
		string alt = position_alt[position];
		int refNum = ref.size();
		int altNum = alt.size();
		//cout << real_position << " " << ref << " " << alt << endl; 
		donorSequence = refSequence.substr(0, real_position) + alt;
		donorSequence += refSequence.substr(real_position+refNum, refLength - real_position - refNum);
	}
	donorLength = donorSequence.size();
	//cout << donorLength << endl;
	return donorSequence;
}
*/
int usage(char* command){
	cout << "\n";
	cout << "\tPlease cite out paper:" << endl;
	cout << "Usage:" << endl;
	cout << "\t" << command << " -r reference -d donor {-s snp_file, -l deletion_file} -h high_confident_vcf  -w low_confident_vcf" << endl;
	cout << "\t reference and donor are mandatory. snp file and deletion file are alternative, you can offer both or just one of them." << endl;
	cout << "\t snp is in vcf format." << endl;
	cout << "\t deletion is now in vcf file format, using Delly's standard." << endl;
}

int test(){
	string info = "AB=0.24;ABP=17.6895;AC=1;AF=0.5;AN=2;AO=6;CIGAR=1X;DP=25;DPB=25;DPRA=0;EPP=16.0391;EPPR=4.03889;GTI=0;LEN=1;MEANALT=1;MQM=23;MQMR=33.7368;NS=1";
	cout << checkInfo(info) << endl;
}

int main(int argc, char* argv[])
{
	//test();
	//return 0;

	if(argc < 2){
		usage(argv[0]);
		return 0;
	}

	char refFilename[3000]={' '};
	char donorFilename[3000]={' '};
	char snpFilename[3000]={' '};
	bool snpModification = false;
	char delFilename[3000]={' '};
	bool delModification = false;
	char highConfidentFilename[3000]={' '};
	char lowConfidentFilename[3000]={' '};
	char blankString[10]={' '};
	
	for(int i = 1; i < argc-1; i++)
	{
		if (!strcmp(argv[i],"-r")){
			strcpy(refFilename, argv[++i]);
		}else if(! strcmp(argv[i], "-d")){
			strcpy(donorFilename, argv[++i]);
		}else if(! strcmp(argv[i], "-s")){
			strcpy(snpFilename, argv[++i]);
			snpModification = true;
		}else if(! strcmp(argv[i], "-l")){
			strcpy(delFilename, argv[++i]);
			delModification = true;
		}else if(! strcmp(argv[i], "-h")){
			strcpy(highConfidentFilename, argv[++i]);
		}else if(! strcmp(argv[i], "-w")){
			strcpy(lowConfidentFilename, argv[++i]);
		}
		
	}

	char * firstLine = (char*)calloc(50, sizeof(char));
    int refLength = 0;
    int lineLength = 0;

	bool donorFileExist = true;
	bool highConfidentExist = true;
	bool lowConfidentExist = true;

	if (!strcmp(donorFilename, blankString)) {
		donorFileExist = false;
		cout << "Error: must have a donor file name"<<endl;
	}
	if (!strcmp(highConfidentFilename, blankString)){
		highConfidentExist = false;
		cout << "Error: must have a high confident vcf file name"<<endl;
	}
	if (!strcmp(lowConfidentFilename, blankString)){ 
		lowConfidentExist = false;
		cout << "Error: must have a low confident vcf file name"<<endl;
	}

	if(!donorFileExist or ! highConfidentExist or !lowConfidentExist){
		usage(argv[0]);
		return 0;
	}

	cout << "Start reading reference file..." ;
	dsptime();
	//char * refSequence = (char*)calloc()
	string refSequence = readRef(refFilename, firstLine, refLength, lineLength);
	cout << "End reading reference file.";
	dsptime();
	int donorLength = 0;
	string donorSequence = "";
	if (snpModification){
    	map<int, string> position_ref;
		map<int, string> position_alt;

    	cout << "Start reading snp file...";
    	dsptime();
		readSnp(snpFilename, position_ref, position_alt, highConfidentFilename, lowConfidentFilename);
    	cout << "End reading snp file.";
    	dsptime();
		//int donorLength = 0;

    	cout << "Start creating donor genome...";
    	dsptime();
		//donorSequence = quickRepairing(refSequence, refLength, position_ref, position_alt);
		donorSequence = repairGenome(refSequence, refLength, position_ref, position_alt);
		cout << "End creating donor genome.";
		dsptime();
	}
	//cout << refSequence << endl;
	//cout << refSequence.size() << endl;

	
	if (delModification){
		cout << "Start reading del vcf file...";
		dsptime();
		map<int, int> delMap_start_end = read_del_vcf(delFilename);
		//cout << delMap_start_end.size() << endl;
		cout << "End reading del vcf file.";
		dsptime();
		cout << "Start chop del section...";
		dsptime();
		if (snpModification){
			donorSequence = chopDelSection(donorSequence, delMap_start_end);
		}else {
			donorSequence = chopDelSection(refSequence, delMap_start_end);
		}
		cout << "End chop del section.";
		dsptime();
	}
	//cout << donorSequence << endl;
	cout << "Start write donor genome...";
	dsptime();
	donorLength = donorSequence.length();
	//writeRef(donorFilename, firstLine, refSequence, refLength, lineLength);
	writeRef(donorFilename, firstLine, donorSequence, donorLength, lineLength);
    cout << "End write donor genome.";
    dsptime();
}
