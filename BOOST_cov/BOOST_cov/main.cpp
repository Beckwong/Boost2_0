#include "main.h"
//#include"CountDigit.h"
#include"GetDataSize.h"
#include"LoadData.h"
#include"CalculateMargnialEntropy.h"
#include"CalculateMarginalDistr.h"
#include"CalculateGenoJointDistr.h"
#include"CalProb.h"
#include"PostCorrection.h"
using namespace std;
//from this line
static unsigned int wordbit[65536];

int popcount(uint64 i)
{
	return (wordbit[i&0xFFFF]+wordbit[(i>>16)&0xFFFF]+wordbit[(i>>32)&0xFFFF]+wordbit[i>>48]);
}

double Abs(double a){
	return (a<0)?-a:a;
}

int bitCount(uint64 i){
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
	i = i + (i >> 8);
	i = i + (i >> 16);
	i = i + (i >> 32);
	return (int)i & 0x7f;
}
//to this line
int main()
{
	//vector<vector<genotype*>>pgeno;

	vector<vector<int>>nCase_Gender(2,vector<int>(2,0));
	vector<vector<int>>nlongint(2,vector<int>(2,0));
	vector<int>DataSize;
	int nsnps=0,nsamples=0;

	vector<vector<genotype*>>pgeno;
	//ifstream fp;
	ofstream fout,fout2;
	string filename = "filenamelist.txt";

	vector<double>MarginalEntropySNP;
	vector<double>MarginalEntropySNP_Y;
	vector<double>MarginalAssociation;
	vector<MarginalDistr*>pMarginalDistr;
	vector<vector<int>>GenoJointDistr(4,vector<int>(9,0));
	
	double maxInteraction = INT_MIN;
	double minInteraction = INT_MAX;
	double thresholdRecord =30;

	vector<pair<int,int>>InteractionPairs;
	vector<double>InteractionMeasurePairs;
	//algorithm start here
	for(uint64 i=0;i<65536;i++)
	{
		wordbit[i] = bitCount(i);
	}

	GetDataSize(filename,DataSize,nCase_Gender);

	

	nsamples = nCase_Gender[0][0]+nCase_Gender[0][1]+nCase_Gender[1][0]+nCase_Gender[1][1];
	for(int i=0;i<DataSize.size();i++)
	{
		nsnps+=DataSize[i];
	}
	//print 1
	cout<<"Number of snps: "<<nsnps<<endl;

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			// print 2
			cout<<nCase_Gender[i][j]<<"  ";
			nlongint[i][j]=ceil((double)nCase_Gender[i][j]/LengthLongType);
		}
		cout<<endl;
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
					tmp_genotype->geno[m*2+n].resize(nlongint[m][n]);
				}
			}
			tmp.push_back(tmp_genotype);
		}
		pgeno.push_back(tmp);
	}

	//load data into bit representation
	LoadData(filename,pgeno,DataSize);
	
	//calculate marginal entropy
	double ptmp1=(double)(nCase_Gender[0][0]+nCase_Gender[0][1])/nsamples;
	double MarginalEntropyY = -ptmp1*log(ptmp1)-(1-ptmp1)*log(1-ptmp1);

	MarginalEntropySNP.resize(nsnps);
	MarginalEntropySNP_Y.resize(nsnps);
	MarginalAssociation.resize(nsnps);
	//print 3
	cout<<"Calculating marginal entropy"<<endl;
	CalculateMarginalEntropy(pgeno,DataSize,nCase_Gender,MarginalEntropySNP,MarginalEntropySNP_Y);

	for(int i=0;i<nsnps;i++)
	{
		MarginalAssociation[i]=(-MarginalEntropySNP_Y[i]+MarginalEntropySNP[i]+MarginalEntropyY)*nsamples*2;
	}

	for(int i=0;i<nsnps;i++)
	{
		MarginalDistr* tmp=new MarginalDistr();
		pMarginalDistr.push_back(tmp);
	}
	//print 4
	cout<<"Calculating marginal distribution"<<endl;
	CalculateMarginalDistr(pgeno,nsnps,pMarginalDistr);

	/*vector<vector<double>> Pk_ij(9,vector<double>(2,0));
	vector<vector<double>> Pj_it(6,vector<double>(3,0));
	vector<vector<double>> Pt_ik(6,vector<double>(2,0));
	vector<vector<double>> Pt_jk(6,vector<double>(2,0));*/

	vector<vector<double>> Pjt(3,vector<double>(2,0));
	vector<vector<double>> Pkt(2,vector<double>(2,0));
	vector<double>Pi(3,0);
	vector<double>Pj(3,0);
	vector<double>Pk(2,0);
	vector<double>Pt(2,0);
	//calculate interaction print 5

	//new ptmp2
	vector<double>Pij(9,0);
	vector<double>Pik(6,0);
	vector<double>Pjk(6,0);
	vector<double>Pit(6,0);
	vector<double>Pijt(18,0);
	cout<<"Cauculating interaction"<<endl;
	ofstream tab;
	tab.open("Contigency table.txt",'w');
	for(int snp1=0;snp1<nsnps-1;snp1++)    //nsnps-1
	{
		for(int snp2=snp1+1;snp2<nsnps;snp2++) //nsnps
		{
			bool flag_bit=0;
			CalculateGenoJointDistr(pgeno,GenoJointDistr,snp1,snp2,pMarginalDistr);
			/*for(int i=0;i<GenoJointDistr.size();i++)
			{
				for(int j=0;j<GenoJointDistr[i].size();j++)
					if(GenoJointDistr[i][j]<=5)
						flag_bit=1;
			}*/
			
			//tab<<"Number: "<<snp1<<" and "<<snp2<<endl;
			//for(int i=0;i<GenoJointDistr.size();i++)
			//{
			//	
			//	for(int j=0;j<GenoJointDistr[i].size();j++)
			//		tab<<GenoJointDistr[i][j]<<",";
			//}
			////tab<<endl;
			//tab<<endl;
			/*if(snp1==0&&snp2==1&&flag_bit)
				cout<<"1st and 2nd SNP has almost empty cells"<<endl;
			else if(snp1==0&&snp2==1&&!flag_bit)
				cout<<"Good for 1st and 2nd SNP pair"<<endl;
			if(flag_bit)
				continue;*/
			
			// Pijk/Pij = P(k|i,j)    i for snp1, j for snp2,k for sample    Pk_ij(9,vector<double>(2,0))
			CalProb(GenoJointDistr,nsamples,Pk_ij,Pj_it,Pt_ik,Pt_jk,Pjt,Pkt,Pi,Pj,Pk,Pt);
			
			//calculate Pij,Pik,Pjk
			for(int i=0;i<3;i++)
			{
				for(int t=0;t<2;t++)
				{
					int sum_it=0;
					for(int j=0;j<3;j++)
					{
						for(int k=0;k<2;k++)
						{
							sum_it+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pit[i*2+t]=(double)sum_it/nsamples;
				}
			}
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					int sum_ij=0;
					for(int k=0;k<2;k++)
					{
						for(int t=0;t<2;t++)
						{
							sum_ij+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pij[i*3+j] = (double)sum_ij/nsamples;
				}
			}

			for(int i=0;i<3;i++)
			{
				for(int k=0;k<2;k++)
				{
					int sum_ik=0;
					for(int j=0;j<3;j++)
					{
						for(int t=0;t<2;t++)
						{
							sum_ik+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pik[i*2+k] = (double)sum_ik/nsamples;
				}
			}

			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					int sum_jk=0;
					for(int i=0;i<3;i++)
					{
						for(int t=0;t<2;t++)
						{
							sum_jk+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pjk[j*2+k] = (double)sum_jk/nsamples;
				}
			}

			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					for(int t=0;t<2;t++)
					{	int sum_ijt=0;
						for(int k=0;k<2;k++)
						{
							sum_ijt+=GenoJointDistr[k*2+t][i*3+j];
						}
						Pijt[i*6+j*2+t] = (double)sum_ijt/nsamples;
					}
				}
			}

			double tao=0.0;
			double InteractionMeasure=0.0;
			double ptmp1=0,ptmp2=0;
			for(int i=0;i<3;i++)
			{
				for(int j=0;j<3;j++)
				{
					for(int k=0;k<2;k++)
					{
						for(int t=0;t<2;t++)
						{
							ptmp1 = (double)(GenoJointDistr[k*2+t][i*3+j])/nsamples;
							if(ptmp1>0)
							{
								InteractionMeasure += ptmp1 * log(ptmp1);
							}
							
							if(Pi[i]>0&&Pj[j]>0&&Pk[k]>0&&Pt[t]>0)
								ptmp2 = (double)(Pik[i*2+k]*Pjk[j*2+k]*Pkt[k][t]*Pijt[i*6+j*2+t])/(Pi[i]*Pj[j]*Pk[k]*Pt[t]);
							else
								ptmp2 = 0;
							/*if(Pjt[j][t]>0&&Pkt[k][t]>0&&Pi[i]&&Pj[j]&&Pk[k]&&Pt[t])
								ptmp2 = (double)(Pk_ij[i*3+j][k]*Pj_it[i*2+t][j]*Pt_ik[i*2+k][t]*Pt_jk[j*2+k][t])/(Pjt[j][t]*Pkt[k][t])*(Pi[i]*Pj[j]*Pk[k]*Pt[t]);
							else
								ptmp2=0;*/
							/*if(Pjt[j][t]>0 && Pkt[k][t]>0 && Pij[i*3+j]>0 && Pik[i*2+k]>0 && Pjk[j*2+k]>0)
							{
								ptmp2 = (double)(Pj_it[i*2+t][j])/(Pjt[j][t]*Pkt[k][t]*Pij[i*3+j]*Pik[i*2+k]*Pjk[j*2+k])*(Pi[i]*Pj[j]*Pk[k]*Pt[t]);
								
							}*/
							if(ptmp2>0)
							{
								InteractionMeasure += -ptmp1*log(ptmp2);
								tao += ptmp2;
							}

						}
					}
				}
			}

			InteractionMeasure += (InteractionMeasure + log(tao))*nsamples*2;
			
			maxInteraction = max(maxInteraction, InteractionMeasure);
			minInteraction = min(minInteraction,InteractionMeasure);

			if(InteractionMeasure > thresholdRecord)
			{
				InteractionPairs.push_back(make_pair(snp1,snp2));
				InteractionMeasurePairs.push_back(InteractionMeasure);
			}

		}
	}
	 tab.close();
	//print 6
	cout<<"Number of interaction pairs: "<<InteractionPairs.size()<<endl;
	cout<<"Max: "<<maxInteraction<<endl;
	cout<<"Min: "<<minInteraction<<endl;
	for(int i=0;i<10 && i<InteractionPairs.size();i++)
		cout<<InteractionMeasurePairs[i]<<endl;
	
	fout2.open("tempRecord.txt",'w');
	int cnt2=0;
	for(int i=0;i<InteractionPairs.size();i++)
	{
		if(InteractionMeasurePairs[i]>thresholdRecord)
		{
			cnt2++;
			fout2<<"No."<<cnt2<<"Pairs:  "<<InteractionPairs[i].first<<" and "<<InteractionPairs[i].second<<" Interaction: "<<InteractionMeasurePairs[i]<<endl;
		}
	}
	
	for(int i=0;i<InteractionPairs.size();i++)
	{
		int snp1 = InteractionPairs[i].first;
		int snp2 = InteractionPairs[i].second;
		CalculateGenoJointDistr(pgeno,GenoJointDistr,snp1,snp2,pMarginalDistr);
		bool flag=0;
		InteractionMeasurePairs[i] = 2*(PostCorrection(GenoJointDistr,nsamples,1)-PostCorrection(GenoJointDistr,nsamples,0));
		if(i%1000==0)
			cout<<"Number "<<i<<" correction";
	}

	fout.open("InteractionRecord.txt",'w');
	int cnt=0;
	for(int i=0;i<InteractionPairs.size();i++)
	{
		if(InteractionMeasurePairs[i]>thresholdRecord)
		{
			cnt++;
			fout<<"No."<<cnt<<"Pairs:  "<<InteractionPairs[i].first<<" and "<<InteractionPairs[i].second<<" Interaction: "<<InteractionMeasurePairs[i]<<endl;
		}
	}
	//while(1);
	return 0;
}
