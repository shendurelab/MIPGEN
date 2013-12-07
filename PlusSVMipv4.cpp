#include <string>
#include "SVMipv4.h"
#include "PlusSVMipv4.h"

using namespace std;

PlusSVMipv4::PlusSVMipv4(string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length):SVMipv4(chromosome, scan_start, scan_stop, ext_length, lig_length)
{
	ext_probe_start = scan_start - ext_length;
	ext_probe_stop = scan_start - 1;
	lig_probe_start = scan_stop + 1;
	lig_probe_stop = scan_stop + lig_length;
	strand = "+";
}

void PlusSVMipv4::set_ext_probe_seq (string seq) {
	this->ext_probe_sequence = seq;
}
void PlusSVMipv4::set_lig_probe_seq (string seq) {
	this->lig_probe_sequence = seq;
	this->ligation_junction = seq.substr(0,2);
}
void PlusSVMipv4::set_scan_target_seq (string seq) {
	this->scan_target_sequence = seq;
}
int PlusSVMipv4::get_mip_start () {
	return this->ext_probe_start;
}
