#include <string>
#include <vector>
#include <map>
using namespace std;

class SVMipv4 
{
	public:
	static map<string, double> junction_scores;	
	int ext_probe_start;
	int ext_probe_stop;
	int lig_probe_start;
	int lig_probe_stop;
	string strand;
	string chr;
	char translocation_failed;
	char snp_failed;
	char mapping_failed;
	char masking_failed;

	int extension_arm_length;
	int ligation_arm_length;
	int scan_size;
	
	string ext_probe_sequence;
	string ext_masked_sequence;
	string lig_probe_sequence;
	string lig_masked_sequence;
	string scan_target_sequence;

	int ext_probe_copy;
	int lig_probe_copy;

	double arm_fraction_masked;

	int scan_start_position;
	int scan_stop_position;

	int snp_count;
	vector<int> snp_positions;

	string mip_seq;
	string ligation_junction;
	double score;
	bool has_snp_mip;
	string snp_ext_sequence;
	string snp_lig_sequence;
	string snp_mip_sequence;

	SVMipv4 (string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length);
	virtual ~SVMipv4(){};
	virtual void set_ext_probe_seq(string seq) =0;
	virtual void set_lig_probe_seq(string seq) =0;
	virtual int get_mip_start() =0;
	double count_ext_mer(string sub);
	double count_lig_mer(string sub);
	double count_insert_mer(string sub);

	void get_parameters(vector<double> & parameters, double long_range_content[]); 
	double get_score (); 
	static void set_junction_scores();
};
