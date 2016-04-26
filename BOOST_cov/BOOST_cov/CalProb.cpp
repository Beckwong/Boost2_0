#include"CalProb.h"
#include "fstream"
using namespace std;

void CalProb(vector<vector<int>>&GenoJointDistr,int nsamples,vector<vector<double>>&Pk_ij,vector<vector<double>>&Pj_it,vector<vector<double>>&Pt_ik,vector<vector<double>>&Pt_jk,vector<vector<double>>&Pjt,vector<vector<double>>&Pkt,vector<double>&Pi,vector<double>&Pj,vector<double>&Pk,vector<double>&Pt)
{
			for(int i=0;i<3;i++)
			{
				int sum_i=0;
				for(int j=0;j<3;j++)
				{
					int sum_ij=GenoJointDistr[0][i*3+j]+GenoJointDistr[1][i*3+j]+GenoJointDistr[2][i*3+j]+GenoJointDistr[3][i*3+j];
					Pk_ij[i*3+j][0] = (double)(GenoJointDistr[0][i*3+j]+GenoJointDistr[1][i*3+j])/sum_ij;
					Pk_ij[i*3+j][1] = 1-Pk_ij[i*3+j][0];
					sum_i+=sum_ij;
				}
				Pi[i]=(double)sum_i/nsamples;

			}

			//Pijt/Pit = P(j|i,t)
			for(int t=0;t<2;t++)
			{
				int sum_t=0;
				for(int i=0;i<3;i++)
				{
					int sum_it=GenoJointDistr[t][i*3]+GenoJointDistr[t][i*3+1]+GenoJointDistr[t][i*3+2]+GenoJointDistr[t+2][i*3]+GenoJointDistr[t+2][i*3+1]+GenoJointDistr[t+2][i*3+2];
					Pj_it[i*2+t][0] = (double)(GenoJointDistr[t][i*3]+GenoJointDistr[t+2][i*3])/sum_it;
					Pj_it[i*2+t][1] = (double)(GenoJointDistr[t][i*3+1]+GenoJointDistr[t+2][i*3+1])/sum_it;
					Pj_it[i*2+t][2] = 1-Pj_it[i*2+t][0]-Pj_it[i*2+t][1];
					sum_t+=sum_it;
				}
				Pt[t]=(double)sum_t/nsamples;
			}

			//Pikt/Pik = P(t|i,k)

			for(int i=0;i<3;i++)
			{
				for(int k=0;k<2;k++)
				{
					int sum_ik = GenoJointDistr[k][i*3]+GenoJointDistr[k][i*3+1]+GenoJointDistr[k][i*3+2]+GenoJointDistr[k+1][i*3]+GenoJointDistr[k+1][i*3+1]+GenoJointDistr[k+1][i*3+2];
					Pt_ik[i*2+k][0] = (double)(GenoJointDistr[k][i*3]+GenoJointDistr[k][i*3+1]+GenoJointDistr[k][i*3+2])/sum_ik;
					Pt_ik[i*2+k][1] = 1.0-Pt_ik[i*2+k][0];
				}
			}

			// Pjkt/Pjk = P()
			for(int j=0;j<3;j++)
			{
				for(int k=0;k<2;k++)
				{
					int sum_jk = GenoJointDistr[k][j]+ GenoJointDistr[k][1*3+j]+ GenoJointDistr[k][2*3+j]+ GenoJointDistr[k+1][j]+GenoJointDistr[k+1][1*3+j]+ GenoJointDistr[k+1][2*3+j];
					Pt_jk[j*2+k][0] = (double)(GenoJointDistr[k][j]+ GenoJointDistr[k][1*3+j]+ GenoJointDistr[k][2*3+j])/sum_jk;
					Pt_jk[j*2+k][1] = 1-Pt_jk[j*2+k][0];
				}
			}

			// 1/Pjt
			for(int j=0;j<3;j++)
			{
				int sum_t=0;
				for(int t=0;t<2;t++)
				{
					int sum_jt=0;
					for(int i=0;i<3;i++)
					{
						for(int k=0;k<2;k++)
						{
							sum_jt+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pjt[j][t]=(double)sum_jt/nsamples;
					sum_t+=sum_jt;
				}
				Pj[j]=(double)sum_t/nsamples;
			}

			// 1/Pkt
			for(int k=0;k<2;k++)
			{
				int sum_k=0;
				for(int t=0;t<2;t++)
				{
					int sum_kt=0;
					for(int i=0;i<3;i++)
					{
						for(int j=0;j<3;j++)
						{
							sum_kt+=GenoJointDistr[k*2+t][i*3+j];
						}
					}
					Pkt[k][t]=(double)sum_kt/nsamples;
					sum_k+=sum_kt;
				}
				Pk[k]=(double)sum_k/nsamples;
			}







			//////////////////////////////////////
			//output for debug////////////////////
			//////////////////////////////////////
			/*
			vector<vector<double>> Pk_ij(9,vector<double>(2,0));
			vector<vector<double>> Pj_it(6,vector<double>(3,0));
			vector<vector<double>> Pt_ik(6,vector<double>(2,0));
			vector<vector<double>> Pt_jk(6,vector<double>(2,0));
			vector<vector<double>> Pjt(3,vector<double>(2,0));
			vector<vector<double>> Pkt(2,vector<double>(2,0));
			vector<double>Pi(3,0);
			vector<double>Pj(3,0);
			vector<double>Pk(2,0);
			vector<double>Pt(2,0);
			*/
			/*ofstream fcout;
			fcout.open("calculation_probability.txt",'w');
			fcout<<"Output GenoJointDistr:"<<endl;
			for(int i=0;i<GenoJointDistr.size();i++)
			{
				for(int j=0;j<GenoJointDistr[i].size();j++)
				{
					fcout<<GenoJointDistr[i][j]<<"  ";
				}
				fcout<<endl;
			}
			fcout<<"Pk_ij:"<<endl;
			for(int i=0;i<Pk_ij.size();i++)
			{
				for(int j=0;j<Pk_ij[i].size();j++)
				{
					fcout<<Pk_ij[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pj_it"<<endl;
			for(int i=0;i<Pj_it.size();i++)
			{
				for(int j=0;j<Pj_it[i].size();j++)
				{
					fcout<<Pj_it[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pt_ik"<<endl;
			for(int i=0;i<Pt_ik.size();i++)
			{
				for(int j=0;j<Pt_ik[i].size();j++)
				{
					fcout<<Pt_ik[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pt_jk"<<endl;
			for(int i=0;i<Pt_jk.size();i++)
			{
				for(int j=0;j<Pt_jk[i].size();j++)
				{
					fcout<<Pt_jk[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pjt"<<endl;
			for(int i=0;i<Pjt.size();i++)
			{
				for(int j=0;j<Pjt[i].size();j++)
				{
					fcout<<Pjt[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pkt"<<endl;
			for(int i=0;i<Pkt.size();i++)
			{
				for(int j=0;j<Pkt[i].size();j++)
				{
					fcout<<Pkt[i][j]<<" ";
				}
				fcout<<endl;
			}
			fcout<<"Pi:"<<endl;
			for(int i=0;i<Pi.size();i++)
				fcout<<Pi[i]<<" ";

			fcout<<"Pj: "<<endl;
			for(int i=0;i<Pj.size();i++)
				fcout<<Pj[i]<<" ";	

			fcout<<"Pk: "<<endl;
			for(int i=0;i<Pk.size();i++)
				fcout<<Pk[i]<<" ";
			fcout<<"Pt: "<<endl;
			for(int i=0;i<Pt.size();i++)
				fcout<<Pt[i]<<" ";	

			fcout.close();*/

}