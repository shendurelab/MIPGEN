#include <string>
#include "Featurev5.h"
using namespace std;
#define MER_NUM 44
Featurev5::Featurev5 (string chromosome, int start, int stop, int f, string l)
{
	chr = chromosome;
	label = l;
	start_position = start;
	stop_position = stop;
	flank_size = f;
	mip_count = 0;
	start_position_flanked = start - f;
	stop_position_flanked = stop + f;
};
Featurev5::Featurev5()
{}
void Featurev5::get_long_range_content(string extended_sequence, string feature_mers[])
{
	double forward_count;
	double reverse_count;
	string rev_comp;
	for(int i = 0; i < MER_NUM; i++)
	{
		forward_count = 0;
		for (size_t offset = extended_sequence.find(feature_mers[i]); offset != string::npos; offset = extended_sequence.find(feature_mers[i], offset + 1))
		{
			forward_count++;
		}
		string for_seq = feature_mers[i];
		rev_comp = "";
		for(string::reverse_iterator it = for_seq.rbegin(); it != for_seq.rend();it++)
		{
			switch(*it)
			{
				case 'G':rev_comp += "C"; break;
				case 'C':rev_comp += "G"; break;
				case 'A':rev_comp += "T"; break;
				case 'T':rev_comp += "A"; break;
			}
		}
		if(rev_comp.compare(feature_mers[i])!=0)
		{
			reverse_count = 0;
			for (size_t offset = extended_sequence.find(rev_comp); offset != string::npos; offset = extended_sequence.find(rev_comp, offset + 1))
			{
				reverse_count++;
			}
			long_range_content[i] = (forward_count + reverse_count) / (this->chromosomal_sequence_stop_position - this->chromosomal_sequence_start_position + 2001);
		}
		else
		{
			long_range_content[i] = forward_count / (this->chromosomal_sequence_stop_position - this->chromosomal_sequence_start_position + 2001);
		}
	}
}
bool Featurev5::operator < (Featurev5 & b)
{
	if(this->chr == b.chr) return this->start_position < b.start_position;
	else return this->chr < b.chr;
}

