#include <string>

using namespace std;

class MinusSVMipv4: public SVMipv4
{
	public:
	MinusSVMipv4(string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length);
	void set_ext_probe_seq(std::string);
	void set_lig_probe_seq(std::string);
	int get_mip_start();
	void set_scan_target_seq (string seq);
};
