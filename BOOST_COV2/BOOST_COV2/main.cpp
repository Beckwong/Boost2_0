// BOOST_COV2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "main.h"
#include"CountDigit.h"
#include"GetDataSize.h"

using namespace std;
const int LengthLongType = 64;

int _tmain(int argc, _TCHAR* argv[])
{
	//vector<vector<genotype*>>pgeno;

	vector<vector<int>>nCase_Gender(2,vector<int>(2,0));
	vector<vector<int>>nlongint(2,vector<int>(2,0));
	vector<int>DataSize;
	int nsnps=0,nsamples=0;

	vector<vector<genotype*>>pgeno;
	ifstream fp;
	string filename = "filenamelist.txt";

	CountDigit();
	GetDataSize(filename,DataSize,nCase_Gender);

	nsamples = nCase_Gender[0][0]+nCase_Gender[0][1]+nCase_Gender[1][0]+nCase_Gender[1][1];
	for(int i=0;i<DataSize.size();i++)
	{
		nsnps+=DataSize[i];
	}
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			nlongint[i][j]=ceil((double)nCase_Gender[i][j]/LengthLongType);
		}
	}

	for(int i =0;i<nsnps;i++)
	{
		vector<genotype*>tmp;
		for(int j=0;j<3;j++)
		{
			genotype *tmp_genotype=new genotype;
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					tmp_genotype->geno[m*2+n].resize(nCase_Gender[m][n]);
				}
			}
			tmp.push_back(tmp_genotype);
		}
		pgeno.push_back(tmp);
	}

	/*for(int i=0;i<nsnps;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int m=0;m<2;m++)
			{
				for(int n=0;n<2;n++)
				{
					pgeno[i][j]
				}
			}
		}
	}*/
}

