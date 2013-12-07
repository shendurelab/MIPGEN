#include <string>
using namespace std;

class PlusSVMipv4:public SVMipv4
{
	public:
	PlusSVMipv4(string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length);
	void set_scan_target_seq (string seq);
	void set_ext_probe_seq (string seq);
	void set_lig_probe_seq (string seq);
	int get_mip_start();

};
