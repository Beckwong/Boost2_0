#include "GetDataSize.h"
#include "stdafx.h"
#include<vector>
#include<string>
#include<time.h>
#include<fstream>
#include<iostream>
using namespace std;
void count_sample(ifstream& fp_i,vector<vector<int> >&DataSize,int& index,vector<int>&nCase_Ctrl);

void GetDataSize(string filename,vector<vector<int>>& DataSize,vector<int>& nCase_Ctrl){
	ifstream fp,fp_i;
//	int c, ndataset;
	time_t st,ed;
//	int n_sample,nSNP,i,flag,index;
	int index=0;
	fp.open(filename);
	st = time(NULL);
	
	while(!fp.eof()){
		string filename_i;
		getline(fp,filename_i);
		fp_i.open(filename_i);
		count_sample(fp_i,DataSize,index,nCase_Ctrl);
		fp_i.close();
	}
	fp.close();

	ed = time(NULL);
	cout<<"Used time: "<<ed-st<<endl;

}

void count_sample(ifstream& fp_i,vector<vector<int> >&DataSize,int& index,vector<int>&nCase_Ctrl){
	vector<int>data;
	int n=0,p=0;
	int ncase=0,nctrl=0;
	string line;
	while(!fp_i.eof()){
		getline(fp_i,line);
		n++;
		if(index==0){
			if(line.size()==0)continue;
			if(line[0]=='1')ncase++;
			else if(line[0]=='0')nctrl++;
		}

		if(n==1){
			for(int i=0;i<line.size();i++){
				if(isdigit(line[i]))p++;
			}
		}
	}
	
	if(index==0){
	nCase_Ctrl.push_back(ncase);
	nCase_Ctrl.push_back(nctrl);
	}
	index++;
	data.push_back(n);
	data.push_back(p-1);
	DataSize.push_back(data);
	return;
}