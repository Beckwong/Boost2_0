#include "stdafx.h"

#include "LoadData.h"
#include "main.h"
#include<vector>
#include<string>
#include<time.h>
#include<fstream>
#include<iostream>
using namespace std;

void LoadData(string filename,vector<vector<genotype*>>&pgeno,vector<vector<int> >&DataSize){
	
	ifstream filep(filename);
	int col=0,row=0,file_index=0;
	int visited_snp=0;
	unsigned long long mask1 = 0x0000000000000001;
	while(!filep.eof())
	{
		int icase, ictrl;
		string tmp_file;

		file_index++;
		getline(filep,tmp_file);

		ifstream filep_i(tmp_file);
		
		col=0;row=0;
		icase = -1;
		ictrl = -1;
		cout<<"Loading data in file "<<file_index<<" : "<<tmp_file;

		while(!filep_i.eof())
		{
			string s;
			getline(filep_i,s);
			if(s.size()==0)continue;

			if(s[0]=='1')
			{
				icase++;
				int cnt_snp=0;
				for(int i=1;i<s.size();i++)
				{
					if(s[i]==' ')continue;
					else if(s[i]>='0' && s[i]<='2')
					{
						pgeno[visited_snp+cnt_snp][s[i]-'0']->genocase[icase/LengthLongType] |= (mask1 << (icase%LengthLongType));
						cnt_snp++;
					}

				}
			}
			else if(s[0]=='0')
			{
				ictrl++;
				int cnt_snp=0;
				for(int i=1;i<s.size();i++)
				{
					if(s[i]==' ')continue;
					else if(s[i]<='2' &&  s[i]>='0')
					{
						pgeno[visited_snp+cnt_snp][s[i]-'0']->genoctrl[ictrl/LengthLongType] |= (mask1 << (ictrl%LengthLongType));
						cnt_snp++;
					}
				}
			}

		}
		filep_i.close();
		visited_snp+=DataSize[file_index-1][1];
	}
	filep.close();


}