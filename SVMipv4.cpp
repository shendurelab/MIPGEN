#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include "SVMipv4.h"

using namespace std;

#define ARM_DIMER_COUNT 20
#define FEATURE_TRIMER_COUNT 44
#define INSERT_TRIMER_COUNT 84
#define TOTAL_FEATURES 192
#define BINS 100.
SVMipv4::SVMipv4 (string chromosome, int scan_start, int scan_stop, int ext_length, int lig_length)
{
	translocation_failed = '0';
	snp_failed = '0';
	mapping_failed = '0';

	chr = chromosome;
	scan_start_position = scan_start;
	scan_stop_position = scan_stop;
	extension_arm_length = ext_length;
	ligation_arm_length = lig_length;
	snp_count = 0;
	scan_size = scan_stop - scan_start + 1;
	has_snp_mip = false;
}
double SVMipv4::count_ext_mer(string sub)
{
	double count = 0;
	for (size_t offset = this->ext_probe_sequence.find(sub); offset != string::npos; offset = this->ext_probe_sequence.find(sub, offset + 1))
	{
		count++;
	}
	return count;
}
double SVMipv4::count_lig_mer(string sub)
{
	double count = 0;
	for (size_t offset = this->lig_probe_sequence.find(sub); offset != string::npos; offset = this->lig_probe_sequence.find(sub, offset + 1))
	{
		count++;
	}
	return count;
}
double SVMipv4::count_insert_mer(string sub)
{
	double count = 0;
	for (size_t offset = this->scan_target_sequence.find(sub); offset != string::npos; offset = this->scan_target_sequence.find(sub, offset + 1))
	{
		count++;
	}
	return count;
}


void SVMipv4::get_parameters(vector<double> & parameters, double long_range_content[]) 
{
	parameters.clear();
	if (this->ext_probe_sequence.find("N") < (unsigned int) this->extension_arm_length || this->lig_probe_sequence.find("N") < (unsigned int) this->ligation_arm_length || this->mip_seq.find("-") < this->mip_seq.length())
	{
		for(int i = 0; i < TOTAL_FEATURES; i++)
			parameters.push_back(0);
		return;
	}
	string arm_mers [] = {"A","AA","AC","AG","AT","C","CA","CC","CG","CT","G","GA","GC","GG","GT","T","TA","TC","TG","TT"};
	string insert_mers [] = {"A", "AA", "AAA", "AAC", "AAG", "AAT", "AC", "ACA", "ACC", "ACG", "ACT", "AG", "AGA", "AGC", "AGG", "AGT", "AT", "ATA", "ATC", "ATG", "ATT", "C", "CA", "CAA", "CAC", "CAG", "CAT", "CC", "CCA", "CCC", "CCG", "CCT", "CG", "CGA", "CGC", "CGG", "CGT", "CT", "CTA", "CTC", "CTG", "CTT", "G", "GA", "GAA", "GAC", "GAG", "GAT", "GC", "GCA", "GCC", "GCG", "GCT", "GG", "GGA", "GGC", "GGG", "GGT", "GT", "GTA", "GTC", "GTG", "GTT", "T", "TA", "TAA", "TAC", "TAG", "TAT", "TC", "TCA", "TCC", "TCG", "TCT", "TG", "TGA", "TGC", "TGG", "TGT", "TT", "TTA", "TTC", "TTG", "TTT"};
	string junctions [] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
	for(int i = 0; i < ARM_DIMER_COUNT; i++) // extension
	{
		if(arm_mers[i].compare("T") == 0)
		{
			parameters.push_back((this->count_ext_mer("G") + this->count_ext_mer("C")) / (this->extension_arm_length - arm_mers[i].length() + 1));
		}
		parameters.push_back(this->count_ext_mer(arm_mers[i]) / (this->extension_arm_length - arm_mers[i].length() + 1.));
	}
	parameters.push_back(this->extension_arm_length);
	for(int i = 0; i < FEATURE_TRIMER_COUNT; i++) // feature
	{
		parameters.push_back(long_range_content[i]);
	}
	for(int i = 0; i < INSERT_TRIMER_COUNT; i++) // insert
	{
		if(insert_mers[i].compare("T") == 0)
		{
			 parameters.push_back((this->count_insert_mer("G") + this->count_insert_mer("C"))/ (this->scan_size - insert_mers[i].length() + 1.));
		}
		parameters.push_back(this->count_insert_mer(insert_mers[i])/ (this->scan_size - insert_mers[i].length() + 1.));
	}
	parameters.push_back(this->scan_size);
	for(int i = 0; i < ARM_DIMER_COUNT; i++) // ligation
	{
		if(arm_mers[i].compare("T") == 0)
		{
			parameters.push_back((this->count_lig_mer("G") + this->count_lig_mer("C"))/ (this->ligation_arm_length - arm_mers[i].length() + 1.));
		}
		parameters.push_back(this->count_lig_mer(arm_mers[i]) / (this->ligation_arm_length - arm_mers[i].length() + 1.));
	}
	parameters.push_back(this->ligation_arm_length);
	string lj = this->lig_probe_sequence.substr(0,2);
	for(int i = 0; i < 16; i++) // ligation junctions
	{
		parameters.push_back(lj.compare(junctions[i])==0 ? 1 : 0);
	}		
	
	double log_ext_copy = this->ext_probe_copy > 100 ? 2 : log10(this->ext_probe_copy); // log of ext copy
	double log_lig_copy = this->lig_probe_copy > 100 ? 2 : log10(this->lig_probe_copy); // log of lig copy
	parameters.push_back(log_ext_copy);
	parameters.push_back(log_lig_copy);
}
double SVMipv4::get_score () 
{
	if (this->ext_probe_sequence.find("N") < (unsigned int) this->extension_arm_length || this->lig_probe_sequence.find("N") < (unsigned int) this->ligation_arm_length || this->mip_seq.find("-") < this->mip_seq.length()) return -1000.0;

	string last_base = this->scan_target_sequence.substr(0, 1);
	double run_count = 0;
   	for (int i = 1; i< this->scan_size; i++)
	{
		string current_base = this->scan_target_sequence.substr(i, 1);
		if (current_base == "G" || current_base == "C")
		{	
			if (last_base == "G" || last_base == "C") {}
			else
			{
				run_count++;
				last_base = current_base;
			}
		}
		else
		{
	       		if ("A" == last_base || "T" == last_base) {}
	       		else {
				run_count++;
				last_base = current_base;
			}
		}
	}
	run_count++;
	double bases_per_switch = this->scan_size / run_count;
	double ext_g_count = (double) count(this->ext_probe_sequence.begin(), this->ext_probe_sequence.end(), 'G');
	double lig_g_count = (double) count(this->lig_probe_sequence.begin(), this->lig_probe_sequence.end(), 'G');
	double target_g_count = (double) count(this->scan_target_sequence.begin(), this->scan_target_sequence.end(), 'G');

	double ext_gc_count = (double) count(this->ext_probe_sequence.begin(), this->ext_probe_sequence.end(), 'C') + ext_g_count;
	double lig_gc_count = (double) count(this->lig_probe_sequence.begin(), this->lig_probe_sequence.end(), 'C') + lig_g_count;
	double target_gc_count = (double) count(this->scan_target_sequence.begin(), this->scan_target_sequence.end(), 'C') + target_g_count;

	double ext_a_count = (double) count(this->ext_probe_sequence.begin(), this->ext_probe_sequence.end(), 'A');
	double lig_a_count = (double) count(this->lig_probe_sequence.begin(), this->lig_probe_sequence.end(), 'A');
	double target_a_count = (double) count(this->scan_target_sequence.begin(), this->scan_target_sequence.end(), 'A');

	double ext_length = this->extension_arm_length;
	double lig_length = this->ligation_arm_length;
	double target_length = this->scan_size > 250 ? 250 : this->scan_size;

	double ext_gc_content = ext_gc_count/ext_length;
	double lig_gc_content = lig_gc_count/lig_length;
	double target_gc_content = target_gc_count/this->scan_size;

	double ext_g_content = ext_g_count/ext_length;
	double lig_g_content = lig_g_count/lig_length;
	double target_g_content = target_g_count/this->scan_size;

	double ext_a_content = ext_a_count/ext_length;
	double lig_a_content = lig_a_count/lig_length;
	double target_a_content = target_a_count/this->scan_size;

	double junction_score = junction_scores[ligation_junction];
	
	double log_ext_copy = this->ext_probe_copy > 100 ? 2 : log10(this->ext_probe_copy);
	double log_lig_copy = this->lig_probe_copy > 100 ? 2 : log10(this->lig_probe_copy);

	double exponent =
		-35.0464 - 4 +
		-1.974282*bases_per_switch +
		2.63667*bases_per_switch*ext_gc_content +
		2.540741*bases_per_switch*lig_gc_content +
		-0.006488*bases_per_switch*target_length +
		-0.018137*ext_a_content +
		7.795877*pow(ext_a_content, 2) +
		-5.576753*ext_a_content*ext_g_content +
		-0.274062*ext_a_content*ext_length +
		8.695568*ext_a_content*lig_g_content +
		-2.014126*ext_a_content*log_ext_copy +
		3.163087*ext_a_content*log_lig_copy +
		-19.900678*pow(ext_gc_content, 2) +
		35.747084*ext_gc_content +
		-2.082136*ext_g_content +
		7.204324*pow(ext_g_content, 2) +
		11.73888*ext_g_content*lig_g_content +
		-2.173235*ext_g_content*log_lig_copy +
		-7.123214*ext_g_content*target_gc_content +
		1.068617*ext_length +
		-0.008666*pow(ext_length, 2) +
		-0.555599*ext_length*ext_gc_content +
		-0.289857*ext_length*lig_g_content +
		-0.009621*ext_length*lig_length +
		0.119863*ext_length*log_ext_copy +
		2.405833*junction_score +
		-1.764289*junction_score*lig_gc_content +
		2.112564*junction_score*lig_g_content +
		0.656183*junction_score*log_ext_copy +
		-3.099451*junction_score*target_a_content +
		-2.097335*junction_score*target_gc_content +
		6.542827*pow(lig_a_content, 2) +
		-8.885956*lig_a_content +
		5.297562*lig_a_content*ext_g_content +
		13.042436*lig_a_content*lig_gc_content +
		-12.333361*lig_a_content*lig_g_content +
		19.866468*lig_gc_content +
		-13.698276*pow(lig_gc_content, 2) +
		-0.517301*lig_gc_content*lig_length +
		-2.846777*lig_gc_content*log_lig_copy +
		13.64081*lig_gc_content*target_a_content +
		13.614709*lig_g_content +
		-4.759165*lig_g_content*ext_gc_content +
		-7.48883*lig_g_content*lig_gc_content +
		-2.308594*lig_g_content*log_ext_copy +
		-12.640154*lig_g_content*target_a_content +
		1.164626*lig_length +
		-0.010326*pow(lig_length, 2) +
		-0.354197*lig_length*ext_gc_content +
		0.111448*lig_length*log_ext_copy +
		0.238095*lig_length*target_a_content +
		-4.161632*log_ext_copy +
		-2.728864*log_ext_copy*ext_gc_content +
		0.641717*log_ext_copy*log_lig_copy +
		3.738798*log_ext_copy*target_gc_content +
		-1.98457*log_lig_copy +
		2.362253*log_lig_copy*ext_gc_content +
		3.467229*log_lig_copy*target_gc_content +
		-18.443242*target_a_content +
		20.89245*pow(target_a_content, 2) +
		-0.048679*target_a_content*target_length +
		-50.249451*pow(target_gc_content, 2) +
		27.132716*target_gc_content +
		-0.050633*target_gc_content*target_length +
		-20.772366*target_g_content +
		-60.796481*pow(target_g_content, 2) +
		26.630245*target_g_content*target_a_content +
		87.162648*target_g_content*target_gc_content +
		0.030256*target_g_content*target_length +
		0.032811*target_length;	
	return pow(2.71828, exponent)/(1 + pow(2.71828, exponent));
};
void SVMipv4::set_junction_scores()
{
		SVMipv4::junction_scores["AA"] = 0.0;
		SVMipv4::junction_scores["AC"] = 0.35;
		SVMipv4::junction_scores["AG"] = 0.046;
		SVMipv4::junction_scores["AT"] = 0.079;
		SVMipv4::junction_scores["CA"] = 0.34;
		SVMipv4::junction_scores["CC"] = 0.22;
		SVMipv4::junction_scores["CG"] = 0.55;
		SVMipv4::junction_scores["CT"] = -0.071;
		SVMipv4::junction_scores["GA"] = 0.35;
		SVMipv4::junction_scores["GC"] = 0.92;
		SVMipv4::junction_scores["GG"] = 0.24;
		SVMipv4::junction_scores["GT"] = 0.48;
		SVMipv4::junction_scores["TA"] = -0.46;
		SVMipv4::junction_scores["TC"] = -0.35;
		SVMipv4::junction_scores["TG"] = -0.25;
		SVMipv4::junction_scores["TT"] = -0.98;
}
