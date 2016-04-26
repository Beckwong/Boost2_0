#include"LoadData.h"

using namespace std;

void LoadData(string filename, vector<vector<genotype*>>&pgeno,vector<int>&DataSize)
{
	ifstream fp(filename);
	int col=0, row=0,file_index=0;
	int visited_snps=0;
	unsigned long long mask1 = 0x0000000000000001;

	while(!fp.eof())
	{
		vector<vector<int>>icase_gender(2,vector<int>(2,-1));

		string tmp_file;
		file_index++;
		getline(fp,tmp_file);

		ifstream fp_i(tmp_file);

		col=0;row=0;
		
		while(!fp_i.eof())
		{
			string s;
			getline(fp_i,s);

			if(s.size()==0)continue;

			int case_ctrl=s[0]-48,gender=s[2]-48;
			int cnt_snp=0;
			icase_gender[case_ctrl][gender]++;
			for(int i=3;i<s.size();i++)
			{
				if(s[i]==' ')continue;
				
				pgeno[cnt_snp+visited_snps][s[i]-'0']->geno[case_ctrl*2+gender][icase_gender[case_ctrl][gender]/LengthLongType]|= (mask1 << (icase_gender[case_ctrl][gender]%LengthLongType));
				cnt_snp++;
			}

		}
		fp_i.close();
		visited_snps+=DataSize[file_index-1];
		
	}
	fp.close();
}