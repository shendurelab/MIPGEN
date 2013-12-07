#ifndef FEATUREV5_H
#define FEATUREV5_H
#include <string>
#define MER_NUM 44

using namespace std;
class Featurev5
{
	public:
	string chr;
	string label;
	int start_position;
	int stop_position;
	int flank_size;
	int start_position_flanked;
	int stop_position_flanked;
	int mip_count;
	int current_scan_start_position;
	string chromosomal_sequence;
	string masked_chromosomal_sequence;
	int chromosomal_sequence_start_position;
	int chromosomal_sequence_stop_position;
	double long_range_content[MER_NUM];
	Featurev5 (string chromosome, int start, int stop, int f, string l);
	Featurev5();
	void get_long_range_content(string extended_sequence, string feature_mers[]);
	bool operator < (Featurev5 & b);
};
#endif
