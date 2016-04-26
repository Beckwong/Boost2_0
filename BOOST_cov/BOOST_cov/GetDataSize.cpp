#include"GetDataSize.h"
using namespace std;

void GetDataSize(string filename,vector<int>&DataSize,vector<vector<int>>&nCase_Gender)
{
	ifstream fp,fp_i;
	int index=0;
	fp.open(filename);

	while(!fp.eof())
	{
		string filename_i;
		getline(fp,filename_i);
		fp_i.open(filename_i);
		count_sample(fp_i,DataSize,index,nCase_Gender);
		fp_i.close();
	}
	fp.close();

}

void count_sample(ifstream& fp_i,vector<int>&DataSize,int index,vector<vector<int>>&nCase_Gender)
{
	string line;
	int nsamples=0,nsnps=0;
	while(!fp_i.eof())
	{
		getline(fp_i,line);
		if(line.size()==0)continue;
		
		nsamples++;
		if(index == 0)
		{
			nCase_Gender[line[0]-48][line[2]-48]++;
		}
		if(nsamples==1)
		{
			for(int i=0;i<line.size();i++)
			{
				if(isdigit(line[i]))nsnps++;
			}
		}
	}
	DataSize.push_back(nsnps-2);
	index++;
	return;
}