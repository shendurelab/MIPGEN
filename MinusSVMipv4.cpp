#include <string>
#include "SVMipv4.h"
#include "MinusSVMipv4.h"
using namespace std;

string reverse_comp (string seq) {
   	string rev_comp;
	for (int i = seq.length()-1; i >= 0; i--)
	{
		switch(seq[i])
		{
			case 'G':
				rev_comp.append("C");
				break;
			case 'C':
				rev_comp.append("G");
				break;
			case 'A':
				rev_comp.append("T");
				break;
			case 'T':
				rev_comp.append("A");
				break;
			default:
				rev_comp.append(string(1, seq[i]));
	  	 	 }
		}
       	return rev_comp;
}
MinusSVMipv4::MinusSVMipv4(string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length):SVMipv4(chromosome, scan_start, scan_stop, ext_length, lig_length)
{
	ext_probe_start = scan_stop + 1;
	ext_probe_stop = scan_stop + ext_length;
	lig_probe_start = scan_start - lig_length;
	lig_probe_stop = scan_start - 1;
	strand = "-";
}

void MinusSVMipv4::set_ext_probe_seq (string seq) {
	this->ext_probe_sequence = reverse_comp(seq);
}
void MinusSVMipv4::set_lig_probe_seq (string seq) {
	this->lig_probe_sequence = reverse_comp(seq);
	this->ligation_junction = this->lig_probe_sequence.substr(0,2);
}
void MinusSVMipv4::set_scan_target_seq (string seq) {
	this->scan_target_sequence = reverse_comp(seq);
}
int MinusSVMipv4::get_mip_start () {
	return this->lig_probe_start;
}
