/*
Written by Evan Boyle, University of Washington
boylee [at] u.washington.edu
Copyright 2014, all rights reserved
*/
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <sstream>
#include <ctype.h>
#include <errno.h>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "Featurev5.h"
#include "SVMipv4.h"
#include "PlusSVMipv4.h"
#include "MinusSVMipv4.h"
#include "svm.h"

using namespace std;

//LIBSVM identifiers
#define TOTAL_FEATURES 192
static string feature_mers [] = {"A","AA","AAA","AAC","AAG","AAT","AC","ACA","ACC","ACG","AG","AGA","AGC","AGG","AGT","AT","ATA","ATC","ATG","CAG","CG","CGG","G","GAC","GAG","GC","GCG","GG","GGC","GGG","GTG","TA","TAA","TAC","TAG","TC","TCC","TCG","TG","TGA","TGC","TGG","TTC","TTG"};
map<string,double> SVMipv4::junction_scores;
static int max_nr_attr = 64;

// comparison function for sorting BED records
bool compare_regions_to_scan (string a, string b)
{
	if (a[0] == '>' || a[0] == '#') // header
		return true;
	else
	{
		string a_chr;
		string b_chr;
		if (a.substr(0,3) == "chr") a_chr = a.substr(3,a.find_first_of(" \t") - 3);
		else a_chr = a.substr(0,a.find_first_of(" \t"));
		if (b.substr(0,3) == "chr") b_chr = b.substr(3,b.find_first_of(" \t") - 3);
		else b_chr = b.substr(0,b.find_first_of(" \t"));
		if(a_chr != b_chr) return a_chr < b_chr;
		else
		{
			int a_start_begin_index = a.find_first_of(" \t", 0);
			a_start_begin_index = a.find_first_not_of(" \t", a_start_begin_index);
			int a_start_size = a.find_first_of(" \t", a_start_begin_index) - a_start_begin_index;
			if(a_start_size < 1)
			{
				a_start_size = a.length() - a_start_begin_index;
			}
			int a_start_position = boost::lexical_cast<int>(a.substr(a_start_begin_index, a_start_size));
			int b_start_begin_index = b.find_first_of(" \t", 0);
			b_start_begin_index = b.find_first_not_of(" \t", b_start_begin_index);
			int b_start_size = b.find_first_of(" \t", b_start_begin_index) - b_start_begin_index;
			int b_start_position = boost::lexical_cast<int>(b.substr(b_start_begin_index, b_start_size));
			return a_start_position < b_start_position;
		}
	}
}
// mipgen object stores all relevant parameters, containers and functions for the purpose of MIP selection
class mipgen{
public:
//parameter and container setup
bool has_arm_length_option;
bool has_arm_length_sum_option;
bool double_tile_strands_separately;
bool double_tile;
map<string,string> args;
set<int> oligo_sizes;
map<int, list<int> > arm_lengths_by_sum;
list<Featurev5> features_to_scan;
int snp_load_count;
map<string, map<int, string> > chr_snp_positions;
map<int, map<string, set<int> > > unmappable_positions;
map<string, map<int, map<int, int> > > copy_chr_start_stop;
ofstream ALLMIPS;
ofstream COLLAPSEDMIPS;
ofstream PICKEDMIPS;
ofstream SNPMIPS;
ofstream PROGRESS;
ofstream GAPS;
ofstream DOUBLEGAPS;
ofstream MINUSGAPS;
ofstream DOUBLEMINUSGAPS;
set<int> arm_length_sum_set;
vector<char> lig_probe_bases;
vector<char> ext_probe_bases;
int bad_design_count;
map<string, map<string, set<int> > > chr_strand_pos_used_arm_bases;
map<int, map<string, boost::shared_ptr<SVMipv4> > > pos_strand_best_mip;
map<int, map<string, boost::shared_ptr<SVMipv4> > > scan_strand_best_mip;
map<int, map<string, list<boost::shared_ptr<SVMipv4> > > > scan_strand_mip_list;
boost::shared_ptr<SVMipv4> null_pointer;
int ext_tag_length;
int lig_tag_length;
string universal_middle_mip_seq;
string file_dir;
//variable and function declarations
int all_mip_counter;
int collapsed_mip_counter;
int picked_mip_counter;
int feature_counter;
int feature_flank;
int max_arm_copy;
int target_arm_copy;
double lower_score_limit;
double upper_score_limit;
string regions_to_scan;
string project_name;
string bwa_genome_index;
int max_mip_overlap;
int starting_mip_overlap;
int max_capture_size;
int min_capture_size;
int capture_increment;

// instantiates a mipgen object with default values for parameters not explicitly set
mipgen(int argc, char * argv[]) {
	has_arm_length_option = false;
	has_arm_length_sum_option = false;
	double_tile_strands_separately = false;
	double_tile = false;
	snp_load_count = 0;
	bad_design_count = 0;
	all_mip_counter = 0;
	collapsed_mip_counter = 0;
	picked_mip_counter = 0;
	feature_counter = 0;
	size_t dir_end = string(argv[0]).find_last_of("/");
	file_dir = dir_end == (size_t) -1 ? "" : string(argv[0]).substr(0, dir_end) + "/";
	set_default_args();
	string parse_status = parse_command_line(argc,argv);
	if (parse_status.length() != 0)
	{
		cerr << parse_status << endl;
		throw 1;
	}
	int a = system(args["-bwa"].c_str());
	if (a != 256)
	{
		cerr << "load bwa" << endl;
		throw 2;
	}
	if(args["-trf"] != "off")
	{
		a = system(args["-trf"].c_str());
		if (a != 65280)
		{
			cerr << "TRF directory invalid" << endl;
			throw 3;
		}
	}
	parse_arg_values();	
}
// sets default values for mipgen object
void set_default_args() {
	args["-bwa"] = "bwa";
	args["-feature_flank"] = "0";
	args["-check_copy_number"] = "on";
	args["-starting_mip_overlap"] = "0";
	args["-double_tile_strand_unaware"] = "off";
	args["-double_tile_strands_separately"] = "off";
	args["-silent_mode"] = "off";
	args["-tag_sizes"] = "5,0";
	args["-trf"] = "off";
	args["-tabix"] = "tabix";
	args["-seal_both_strands"] = "off";
	args["-half_seal_both_strands"] = "off";
	args["-masked_arm_threshold"] = "0.5";
	args["-snp_file"] = "<none>";
	args["-capture_increment"] = "5";
	args["-max_mip_overlap"] = "30";
	args["-max_mip_overlap"] = "30";
	args["-score_method"] = "logistic";
	args["-lig_min_length"] = "18";
	args["-ext_min_length"] = "16";
	args["-max_arm_copy_product"] = "75";
	args["-target_arm_copy"] = "20";
    args["-bwa_threads"] = "1";
}
// turns parameters into their proper types for easy access by field
void parse_arg_values() {
	feature_flank = boost::lexical_cast<int>(args["-feature_flank"]);
	size_t tag_boundary = args["-tag_sizes"].find(",");
	string ext_tag_string = args["-tag_sizes"].substr(0, tag_boundary);
	string lig_tag_string = args["-tag_sizes"].substr(tag_boundary + 1);
	ext_tag_length = boost::lexical_cast<int>(ext_tag_string);
	lig_tag_length = boost::lexical_cast<int>(lig_tag_string);
	max_arm_copy = boost::lexical_cast<int>(args["-max_arm_copy_product"]);
	target_arm_copy = boost::lexical_cast<int>(args["-target_arm_copy"]);
	string universal_constant_mip_seq = "CTTCAGCTTCCCGATATCCGACGGTAGTGT";
	universal_middle_mip_seq = string(lig_tag_length, 'N') + universal_constant_mip_seq + string(ext_tag_length, 'N');
	if (args["-double_tile_strand_unaware"] == "on")
	{
		double_tile = true;
	}
	if (args["-double_tile_strands_separately"] == "on")
	{
		double_tile_strands_separately = true;
	}
	if(args.find("-svr_optimal_score") == args.end())
		args["-svr_optimal_score"] = "2.2";
	if(args.find("-svr_priority_score") == args.end())
		args["-svr_priority_score"] = "1.5";
	if(args.find("-logistic_optimal_score") == args.end())
		args["-logistic_optimal_score"] = "0.98";
	if(args.find("-logistic_priority_score") == args.end())
		args["-logistic_priority_score"] = "0.9";
	if(args["-score_method"] != "logistic" && args["-score_method"] != "svr" && args["-score_method"] != "mixed")
	{
		cerr << "invalid scoring method given" << endl;
		throw 4;
	}
	if (has_arm_length_option)
	{
		string input (args["-arm_lengths"]);
		vector<string> input_parse;
		boost::split(input_parse, input, boost::is_any_of(","));
		for (size_t i = 0; i < input_parse.size(); i++)
		{
			size_t boundary = input_parse.at(i).find(":");
			int extension_arm_length = boost::lexical_cast<int>(input_parse.at(i).substr(0, boundary));
			int ligation_arm_length = boost::lexical_cast<int>(input_parse.at(i).substr(boundary + 1));
			oligo_sizes.insert(extension_arm_length);
			oligo_sizes.insert(ligation_arm_length);
			arm_lengths_by_sum[extension_arm_length + ligation_arm_length].push_back(extension_arm_length);
		}
	}
	if (has_arm_length_sum_option || !(has_arm_length_sum_option || has_arm_length_option))
	{
		string input;
		if(has_arm_length_sum_option)
			input = args["-arm_length_sums"];
		else
			input = "40,41,42,43,44,45";
		vector<string> input_parse;
		boost::split(input_parse, input, boost::is_any_of(","));
		for (size_t i = 0; i<input_parse.size(); i++)
		{
			int arm_sum = boost::lexical_cast<int>(input_parse.at(i));
			for (int ligation_arm_length = boost::lexical_cast<int>(args["-lig_min_length"]); ligation_arm_length <= arm_sum - boost::lexical_cast<int>(args["-ext_min_length"]) && ligation_arm_length <= 30; ligation_arm_length++)
			{
				int extension_arm_length = arm_sum - ligation_arm_length;
				if (extension_arm_length <= 30)
				{
					arm_lengths_by_sum[arm_sum].push_back(extension_arm_length);
					oligo_sizes.insert(extension_arm_length);
					oligo_sizes.insert(ligation_arm_length);
				}
			}
			arm_lengths_by_sum[arm_sum].sort();
		}
	}
	
	
	lower_score_limit = args["-score_method"] == "svr" ? boost::lexical_cast<double>(args["-svr_priority_score"]) : boost::lexical_cast<double>(args["-logistic_priority_score"]);
	upper_score_limit = args["-score_method"] == "svr" ? boost::lexical_cast<double>(args["-svr_optimal_score"]) : boost::lexical_cast<double>(args["-logistic_optimal_score"]);
	regions_to_scan = args["-regions_to_scan"];
	project_name = args["-project_name"];
	bwa_genome_index = args["-bwa_genome_index"];
	max_mip_overlap = boost::lexical_cast<int>(args["-max_mip_overlap"]);
	starting_mip_overlap = boost::lexical_cast<int>(args["-starting_mip_overlap"]);
	max_capture_size = boost::lexical_cast<int>(args["-max_capture_size"]);
	min_capture_size = boost::lexical_cast<int>(args["-min_capture_size"]);
	capture_increment = boost::lexical_cast<int>(args["-capture_increment"]);
	if (capture_increment == 0) capture_increment = 1; // Prevent any strange behavior associated with 0 value

}
// prints version and parameters for mipgen run
void print_header() {
	PROGRESS.open((project_name + ".progress.txt").c_str());
	if(!PROGRESS.is_open())
	{
		cerr << "progress file could not be opened" << endl;
		throw 5;
	}
	PROGRESS << __FILE__ << endl << "last numbered version: 1.1" << endl;
	PROGRESS << "contact: Evan Boyle\nemail: boylee@uw.edu\n";
	for (map<string,string>::iterator it = args.begin(); it!= args.end(); it++) 
	{
		PROGRESS << it->first << " " << it->second << endl;
	}
}
// collects sequence from genome reference, checks for SNPs, copy number and tandem repeats and prepares for selection
void query_sequences(){
	string feature_status = get_features_to_scan(); // creates Feature objects for targeted regions 
	if (feature_status.length() == 0)
	{
		cerr << "[mipgen] region file could not be opened" << endl;
		throw 6;
	}
	PROGRESS << "successfully loaded features for mip design; retrieving chromosomal sequence\n";
	cerr << "[mipgen] features loaded; retrieving chromosomal sequence\n";
	
	if(args.find("-genome_dir") != args.end())
	{
		if (!(get_chr_fasta_sequence_from_genome_dir(args["-genome_dir"])))
		{
			cerr << "[mipgen] chromosome fasta not acquired" << endl;
			throw 7;
		}
	}
	else
	{
		int samtools_result = system("samtools");
		if(samtools_result != 256) {
			cerr << "[mipgen] load samtools or provide genome directory" << endl;
			throw 8;
		}
		if (!(get_chr_fasta_sequence_using_samtools()))
		{
			cerr << "[mipgen] chromosome fasta not acquired" << endl;
			throw 9;
		}
	}
	if (!(get_masked_features_to_scan()))
	{
		cerr << "[mipgen] masked chromosome fasta not acquired; no repetitive bases?" << endl;
	}
	PROGRESS << "successfully acquired chromosomal data for mip design; accessing snp file ...\n";
	cerr << "[mipgen] regions ready; accessing snp file\n";

	map<string, map<int, string> > chr_snp_positions;

    if(!(args["-snp_file"] == "<none>"))
        load_snps(); // loads snp positions and alleles into memory using tabix

	PROGRESS << "all " << snp_load_count << " snps loaded; generating files for bwa\n";
	cerr << "[mipgen] all " << snp_load_count << " snps loaded; generating files for bwa\n";

	string copy_status = check_copy_numbers(); // loads bad mip start positions into memory
	if (copy_status.length() == 0)
	{
		cerr << "error with copy number analysis" << endl;
		throw 11;
	}
	PROGRESS << copy_status;
	find_copy();
	PROGRESS << "bwa copy number analysis finished\n";
	cerr << "[mipgen] bwa copy number analysis finished\n";
	ALLMIPS.open((project_name + ".all_mips.txt").c_str());

	if (ALLMIPS.is_open())
		cerr << "[mipgen] file of all mips ready for write: " << project_name << ".all_mips.txt\n";
	else
	{
		cerr << "[mipgen] all mips file could not be opened" << endl;
		throw 12;
 	}
	ALLMIPS << ">mip_key\t";
	if(args["-score_method"] == "svr")
		ALLMIPS << "svr"; 
	else
		ALLMIPS << "logistic";
	ALLMIPS << "_score\tchr\text_probe_start\text_probe_stop\text_probe_copy\text_probe_sequence\tlig_probe_start\tlig_probe_stop\tlig_probe_copy\tlig_probe_sequence\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tmip_sequence\tfeature_start_position\tfeature_stop_position\tprobe_strand\tfailure_flags\tmip_name\n";
	PROGRESS << "file of all mips ready for write: " << project_name << ".all_mips\n";

	COLLAPSEDMIPS.open((project_name +".collapsed_mips.txt").c_str());
	if (COLLAPSEDMIPS.is_open())
		cerr << "[mipgen] file of collapsed mips ready for write: " << project_name << ".collapsed_mips.txt\n";
	else
	{
		cerr << "[mipgen] file of collapsed mips could not be opened" << endl;
		throw 13;
	}
	COLLAPSEDMIPS << ">mip_key\t";
	if(args["-score_method"] == "svr")
		COLLAPSEDMIPS << "svr";
	else
		COLLAPSEDMIPS << "logistic";
	COLLAPSEDMIPS << "_score\tchr\text_probe_start\text_probe_stop\text_probe_copy\text_probe_sequence\tlig_probe_start\tlig_probe_stop\tlig_probe_copy\tlig_probe_sequence\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tmip_sequence\tfeature_start_position\tfeature_stop_position\tprobe_strand\tfailure_flags\tmip_name\n";
	PROGRESS << "file of collapsed mips ready for write: " << project_name << ".collapsed_mips.txt\n";

	PICKEDMIPS.open((project_name + ".picked_mips.txt").c_str());
	if (!(PICKEDMIPS.is_open()))
		throw 14;
	PICKEDMIPS << ">mip_key\t";
	if(args["-score_method"] == "logistic")
		PICKEDMIPS << "logistic";
	else
		PICKEDMIPS << "svr";
	PICKEDMIPS << "_score\tchr\text_probe_start\text_probe_stop\text_probe_copy\text_probe_sequence\tlig_probe_start\tlig_probe_stop\tlig_probe_copy\tlig_probe_sequence\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tmip_sequence\tfeature_start_position\tfeature_stop_position\tprobe_strand\tfailure_flags\tmip_name\n";
	SNPMIPS.open((project_name + ".snp_mips.txt").c_str());
	if (!(SNPMIPS.is_open()))
		throw 15;
	SNPMIPS << ">mip_key\t";
	if(args["-score_method"] == "logistic")
		SNPMIPS << "logistic";
	else
		SNPMIPS << "svr";
	SNPMIPS << "_score\tchr\text_probe_start\text_probe_stop\text_probe_copy\text_probe_sequence\tlig_probe_start\tlig_probe_stop\tlig_probe_copy\tlig_probe_sequence\tmip_scan_start_position\tmip_scan_stop_position\tscan_target_sequence\tmip_sequence\tfeature_start_position\tfeature_stop_position\tprobe_strand\tfailure_flags\tmip_name\n";
}
// prints MIP data to files and selects a final MIP tiling 
// this is the main process behind MIP design
void tile_regions()
{
	for(map<int, list<int> >::iterator it = arm_lengths_by_sum.begin(); it != arm_lengths_by_sum.end(); it++)
		arm_length_sum_set.insert(it->first);
	
	Featurev5 * feature;
	svm_model * model = svm_load_model((file_dir + "mipgen_svr.model").c_str());
	vector<double> scoring_parameters;
	
	for(list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		feature = &*it;
		feature_counter++;
		PROGRESS << "designing all mips for feature #" << feature_counter << endl;
		cerr << "[mipgen] feature #" << feature_counter << endl;
		string chr = feature->chr;
		//cout << feature->start_position_flanked << "\t" << feature->stop_position_flanked << endl; //comment
			
		feature->current_scan_start_position = feature->start_position_flanked - max_capture_size + *arm_length_sum_set.rbegin(); // build mips starting from when the start of the feature is as far as possible from the extension arm
		if(feature->current_scan_start_position < 0) feature->current_scan_start_position = 0;
		while (feature->current_scan_start_position < feature->stop_position_flanked) // continue until the feature stop position is the first base following the extension arm
		{
			feature->current_scan_start_position++;
			double previous_best_score = 0;
			for (int capture_size = max_capture_size; capture_size >= min_capture_size; capture_size -= capture_increment)
			{
				if (capture_size > feature->stop_position_flanked - feature->start_position_flanked + max_mip_overlap && capture_size - capture_increment >= min_capture_size) continue;
				if (previous_best_score > upper_score_limit) continue;
				for (set<int>::reverse_iterator rit = arm_length_sum_set.rbegin(); rit != arm_length_sum_set.rend(); rit++)
				{
					int arm_length_sum = *rit;
					if (previous_best_score > upper_score_limit && arm_length_sum != *arm_length_sum_set.begin()) continue;
					int previous_minus_score = 0;
					int previous_plus_score = 0;
					bool skip_ahead = false;
					for(list<int>::iterator it2 = arm_lengths_by_sum[arm_length_sum].begin(); it2 != arm_lengths_by_sum[arm_length_sum].end(); it2++)
					{
						if (skip_ahead) continue;
						int extension_arm_length = *it2;
						int ligation_arm_length = arm_length_sum - extension_arm_length;
						if(feature->current_scan_start_position - extension_arm_length <= 0 || feature->current_scan_start_position - ligation_arm_length <= 0) continue;
						if(feature->current_scan_start_position + capture_size - extension_arm_length - 1 > feature->chromosomal_sequence_stop_position || feature->current_scan_start_position + capture_size - ligation_arm_length - 1 > feature->chromosomal_sequence_stop_position) continue;
						//cout << extension_arm_length << ":" << capture_size << ":" << feature->current_scan_start_position << endl; // comment
						boost::shared_ptr<PlusSVMipv4> current_plus_mip (new PlusSVMipv4 (
							feature->chr,
							feature->current_scan_start_position,
							feature->current_scan_start_position + capture_size - arm_length_sum - 1,
							extension_arm_length,
							ligation_arm_length
						));
						boost::shared_ptr<MinusSVMipv4> current_minus_mip (new MinusSVMipv4 (
								feature->chr,
													feature->current_scan_start_position,
							feature->current_scan_start_position + capture_size - arm_length_sum - 1,
							extension_arm_length,
							ligation_arm_length
						));
	
						current_plus_mip->set_scan_target_seq(feature->chromosomal_sequence.substr(current_plus_mip->scan_start_position - feature->chromosomal_sequence_start_position, current_plus_mip->scan_size));		
						current_minus_mip->set_scan_target_seq(feature->chromosomal_sequence.substr(current_minus_mip->scan_start_position - feature->chromosomal_sequence_start_position,current_minus_mip->scan_size));
						// cout << "ready to design\n"; // comment						
						
						boost::shared_ptr<SVMipv4> designed_plus_mip; 
						designed_plus_mip = design_mip(feature, current_plus_mip);
						if(args["-score_method"] == "logistic" || args["-score_method"] == "mixed")
							designed_plus_mip->score = designed_plus_mip->get_score();
						else if (args["-score_method"] == "svr")
						{
							designed_plus_mip->get_parameters(scoring_parameters, feature->long_range_content);
							designed_plus_mip->score = predict_value(scoring_parameters, model);
						}
						all_mip_counter++;
						scan_strand_mip_list[designed_plus_mip->scan_start_position]["+"].push_front(designed_plus_mip); // hashes with keys derived from strand and scan start position and values being arrays of mips at that position
						if(args["-silent_mode"] != "on")
							ALLMIPS << print_details(feature, designed_plus_mip, all_mip_counter, false); 
						
						boost::shared_ptr<SVMipv4> designed_minus_mip;
						designed_minus_mip = design_mip(feature, current_minus_mip);
						if(args["-score_method"] == "logistic" || args["-score_method"] == "mixed")
							designed_minus_mip->score = designed_minus_mip->get_score();
						else if (args["-score_method"] == "svr")
						{
							designed_minus_mip->get_parameters(scoring_parameters, feature->long_range_content);
							designed_minus_mip->score = predict_value(scoring_parameters, model);
						}
						all_mip_counter++;
						scan_strand_mip_list[designed_minus_mip->scan_start_position]["-"].push_front(designed_minus_mip);
						if(args["-silent_mode"] != "on")
							ALLMIPS << print_details(feature, designed_minus_mip, all_mip_counter, false);
						
						//cout << "design done!\n";	
						if (args["-score_method"] == "logistic" && args["-logistic_heuristic"] != "off" && designed_plus_mip->score < previous_plus_score && designed_minus_mip->score < previous_minus_score) skip_ahead = true;
						previous_best_score = (designed_minus_mip->score > designed_plus_mip->score) ? designed_minus_mip->score : designed_plus_mip->score;
						previous_minus_score = designed_minus_mip->score;
						previous_plus_score = designed_plus_mip->score;
					}
				}
			}
		}
		//cout << "condensing mips!\n";
		condense_mips(feature);
		//cout << "collapsing mips!\n";
		collapse_mips(); // collapses mip space by choosing only one (best scoring) hopefully nonfailed mip for every position
		if(args["-silent_mode"] != "on")
			output_collapsed_mips(feature);
		PROGRESS << "mips collapsed! picking mips...\n";
		//cout << "mips collapsed! picking mips...\n";
		if(args["-score_method"] == "mixed")
		{
			lower_score_limit = boost::lexical_cast<double>(args["-svr_priority_score"]);
			upper_score_limit = boost::lexical_cast<double>(args["-svr_optimal_score"]);
		}
		pick_mips(feature, scoring_parameters, model);
		if(args["-score_method"] == "mixed")
		{
			lower_score_limit = boost::lexical_cast<double>(args["-logistic_priority_score"]);
			upper_score_limit = boost::lexical_cast<double>(args["-logistic_optimal_score"]);
		}
		//cout << "mips picked!\n";
		scan_strand_mip_list.clear();
		scan_strand_best_mip.clear();
		pos_strand_best_mip.clear();
	}
	
	
	ALLMIPS.close();
	COLLAPSEDMIPS.close();
	PICKEDMIPS.close();
	SNPMIPS.close();
	PROGRESS << "mip picking complete:\n" << project_name << ".picked_mips.txt\n and \n" << project_name <<".snps_mips.txt\n";
	cerr << "[mipgen] mip picking complete:\n" << project_name << ".picked_mips.txt\n and \n" << project_name << ".snp_mips.txt\n";
	if (bad_design_count > 0)
	{
		PROGRESS << "WARNING: There are " << bad_design_count << " gaps in covering supplied regions\n";
		cerr << "[mipgen] WARNING: There are " << bad_design_count << " gaps in covering supplied regions\n";
	}
	if(GAPS.is_open())
	{
		GAPS.close();
	}
	if(DOUBLEGAPS.is_open())
	{
		DOUBLEGAPS.close();
	}
	if(MINUSGAPS.is_open())
	{
		MINUSGAPS.close();
	}
	if(DOUBLEMINUSGAPS.is_open())
	{
		DOUBLEMINUSGAPS.close();
	}
	PROGRESS.close();
}
// uses BWA to find probe and target copy numbers
void find_copy ()
{
	system((args["-bwa"] + " aln -t " + args["-bwa_threads"] + " " + bwa_genome_index + " " + project_name + ".oligo_copy_count.fq > " + project_name + ".oligo_copy_count.sai").c_str());
	system((args["-bwa"] + " samse " + bwa_genome_index + " " + project_name + ".oligo_copy_count.sai " + project_name + ".oligo_copy_count.fq > " + project_name + ".oligo_copy_count.sam").c_str());
	ifstream OLIGOSAM ((project_name + ".oligo_copy_count.sam").c_str());
	if (OLIGOSAM.is_open()) cerr << "[mipgen] checking oligo copy" << endl;
	else return;
	string line;
	while (OLIGOSAM.good())
	{
		getline(OLIGOSAM, line);
		if (line[0] != '@' && line.length() > 1)
		{
			size_t chr_start = line.find("chr", 0) + 3;
			size_t chr_stop = line.find(":", chr_start);
			size_t start_start = chr_stop + 1;
			size_t start_stop = line.find("-", start_start);
			size_t stop_start = start_stop + 1;
			size_t stop_stop = line.find("\t", stop_start);

			string chr = line.substr(chr_start, chr_stop - chr_start);
			int start = boost::lexical_cast<int>(line.substr(start_start, start_stop - start_start));
			int stop = boost::lexical_cast<int>(line.substr(stop_start, stop_stop - stop_start));
			size_t copy_start = line.find("X0:i:",0);
			if (copy_start < line.length())
			{
				copy_start += 5;
				size_t copy_stop = line.find("\t", copy_start);
				int copy_count = boost::lexical_cast<int>((line.substr(copy_start, copy_stop - copy_start)));
				copy_chr_start_stop[chr][start][stop] = copy_count;
			}
			else
			{
				copy_chr_start_stop[chr][start][stop] = 100;
			}
		}
	}
	OLIGOSAM.close();
}

// gets the basic information about potential mips, used for scoring purposes
boost::shared_ptr<SVMipv4> design_mip (Featurev5 * current_feature, boost::shared_ptr<SVMipv4> current_mip)
{
	string chr = current_feature->chr;
	current_mip->set_ext_probe_seq(current_feature->chromosomal_sequence.substr(current_mip->ext_probe_start - current_feature->chromosomal_sequence_start_position, current_mip->extension_arm_length));
	current_mip->set_lig_probe_seq(current_feature->chromosomal_sequence.substr(current_mip->lig_probe_start - current_feature->chromosomal_sequence_start_position, current_mip->ligation_arm_length));
	
	current_mip->mip_seq = current_mip->lig_probe_sequence + universal_middle_mip_seq + current_mip->ext_probe_sequence;
	string masked_ext_seq = current_feature->masked_chromosomal_sequence.substr(current_mip->ext_probe_start - current_feature->chromosomal_sequence_start_position, current_mip->extension_arm_length);
	string masked_lig_seq = current_feature->masked_chromosomal_sequence.substr(current_mip->lig_probe_start - current_feature->chromosomal_sequence_start_position, current_mip->ligation_arm_length);
	double ext_N_count = count(masked_ext_seq.begin(), masked_ext_seq.end(), 'N');
	double lig_N_count = count(masked_lig_seq.begin(), masked_lig_seq.end(), 'N');
	current_mip->arm_fraction_masked = (ext_N_count + lig_N_count) / (current_mip->ligation_arm_length + current_mip->extension_arm_length);

	current_mip->ext_probe_copy = copy_chr_start_stop[current_mip->chr][current_mip->ext_probe_start][current_mip->ext_probe_stop];
	current_mip->lig_probe_copy = copy_chr_start_stop[current_mip->chr][current_mip->lig_probe_start][current_mip->lig_probe_stop];

	int capture_size = current_mip->scan_size + current_mip->extension_arm_length + current_mip->ligation_arm_length;
	if (	!(unmappable_positions.empty()) 
		&& !(unmappable_positions[capture_size].empty()) 
		&& unmappable_positions[capture_size][chr].find(current_mip->get_mip_start()) != unmappable_positions[capture_size][chr].end()
		&& args["-check_copy_number"] != "off"
	)
	{
		current_mip->mapping_failed = '1';
		// cout << "=(" << endl;
		return current_mip;
	}
	if (current_mip->arm_fraction_masked > boost::lexical_cast<double>(args["-masked_arm_threshold"]))
	{
		current_mip->masking_failed = '1';
	}
	else
	{
		current_mip->masking_failed = '0';
	}
	for (int i = current_mip->ext_probe_start; i <= current_mip->ext_probe_stop; i++)
	{
		if (!(chr_snp_positions.empty()) && !(chr_snp_positions[chr].empty()) && chr_snp_positions[chr].find(i) != chr_snp_positions[chr].end())
		{
			bool flag = false; // used to track whether snp mip creation succeeded
			string alleles = chr_snp_positions[chr][i];
			current_mip->snp_count++;
			string arm_sequence (current_mip->ext_probe_sequence);
			ext_probe_bases.clear();
			ext_probe_bases = vector<char> (arm_sequence.begin(), arm_sequence.end());
			if (alleles.length() == 2 && alleles[0] != 'N' && alleles[1] != 'N' && alleles[0] != '-' && alleles[1] != '-') // only single nucleotide transitions and transversions are supported
			{
				int relative_ext_position;
				if (current_mip->strand == "+")
				{
					relative_ext_position = i - current_mip->ext_probe_start;
				}
				else
				{
					relative_ext_position = current_mip->ext_probe_stop - i;
				}
				if (ext_probe_bases[relative_ext_position] == alleles[0]) 
				{
					ext_probe_bases[relative_ext_position] = alleles[1];
					flag = true;
				}
				else // snp file might have allele on minus strand
				{
					char allele_1;
					switch(alleles[0]) {
						case 'A': allele_1 = 'T'; break;
						case 'T': allele_1 = 'A'; break;
						case 'G': allele_1 = 'C'; break;
						case 'C': allele_1 = 'G'; break;
					}
					char allele_2;
					switch(alleles[1]) {
						case 'A': allele_2 = 'T'; break;
						case 'T': allele_2 = 'A'; break;
						case 'G': allele_2 = 'C'; break;
						case 'C': allele_2 = 'G'; break;
					}
					if (ext_probe_bases[relative_ext_position] == allele_1)
					{
						ext_probe_bases[relative_ext_position] = allele_2;
						flag = true;
					}
				}
			}
			if (flag)
			{	
				current_mip->has_snp_mip = true;
				current_mip->snp_ext_sequence = string(ext_probe_bases.begin(), ext_probe_bases.end());
				current_mip->snp_lig_sequence = current_mip->lig_probe_sequence;
				current_mip->snp_mip_sequence = current_mip->snp_lig_sequence + universal_middle_mip_seq + current_mip->snp_ext_sequence;
			}
			else // this will flag the necesary mip as having a snp with no alternate mip provided
			{
				current_mip->snp_failed = '1';
			}
		}
	}
	for (int i = current_mip->lig_probe_start; i <= current_mip->lig_probe_stop; i++)
	{
		if (!(chr_snp_positions[chr].empty()) && chr_snp_positions[chr].find(i) != chr_snp_positions[chr].end())
		{
			bool flag = false; // used to track whether snp mip creation succeeded
			string alleles = chr_snp_positions[chr][i];
			current_mip->snp_count++;
			string arm_sequence (current_mip->lig_probe_sequence);
			lig_probe_bases.clear();
			lig_probe_bases = vector<char> (arm_sequence.begin(), arm_sequence.end());
			if (alleles.length() == 2 && alleles[0] != 'N' && alleles[1] != 'N' && alleles[0] != '-' && alleles[1] != '-') // only single nucleotide transitions and transversions are supported
			{
				int relative_lig_position;
				if (current_mip->strand == "+")
				{
					relative_lig_position = i - current_mip->lig_probe_start;
				}
				else
				{
					relative_lig_position = current_mip->lig_probe_stop - i;
				}
				if (lig_probe_bases[relative_lig_position] == alleles[0])
				{
					lig_probe_bases[relative_lig_position] = alleles[1];
					flag = true;
				}
				else // snp file might have allele on minus strand
				{
					char allele_1;
					switch(alleles[0]) {
						case 'A': allele_1 = 'T'; break;
						case 'T': allele_1 = 'A'; break;
						case 'G': allele_1 = 'C'; break;
						case 'C': allele_1 = 'G'; break;
					}
					char allele_2;
					switch(alleles[1]) {
						case 'A': allele_2 = 'T'; break;
						case 'T': allele_2 = 'A'; break;
						case 'G': allele_2 = 'C'; break;
						case 'C': allele_2 = 'G'; break;
					}
					if (lig_probe_bases[relative_lig_position] == allele_1)
					{
						lig_probe_bases[relative_lig_position] = allele_2;
						flag = true;
					}
				}
			}
			if (flag)
			{
				current_mip->has_snp_mip = true;
				current_mip->snp_ext_sequence = current_mip->ext_probe_sequence;
				current_mip->snp_lig_sequence = string(lig_probe_bases.begin(), lig_probe_bases.end());
				current_mip->snp_mip_sequence = current_mip->snp_lig_sequence + universal_middle_mip_seq + current_mip->snp_ext_sequence;
			}

			else // this will flag the necesary mip as having a snp with no alternate mip provided
			{
				current_mip->snp_failed = '1';
			}
		}
	}
	if (current_mip->snp_count > 1)
		current_mip->snp_failed = '1';
	return current_mip;
}

// prints standard MIP design file record
string print_details (Featurev5 * feature, boost::shared_ptr<SVMipv4> mip, int mip_index, bool minor)// used throughout the program to standardize output of MIP features
{
	stringstream ss;
	ss << mip->chr << ":";
	ss << (mip->strand == "+" ? mip->ext_probe_start : mip->lig_probe_start) << "-";
	ss << (mip->strand == "+" ? mip->lig_probe_stop : mip->ext_probe_stop) << "/";
	ss << mip->extension_arm_length << ",";
	ss << mip->ligation_arm_length << "/";
	ss << mip->strand << "\t";
	ss << mip->score << "\t";
	ss << mip->chr << "\t";
	ss << mip->ext_probe_start << "\t";
	ss << mip->ext_probe_stop << "\t";
	ss << mip->ext_probe_copy << "\t";
	ss << mip->ext_probe_sequence << "\t";
	ss << mip->lig_probe_start << "\t";
	ss << mip->lig_probe_stop << "\t";
	ss << mip->lig_probe_copy << "\t";
	ss << mip->lig_probe_sequence << "\t";
	ss << mip->scan_start_position << "\t";
	ss << mip->scan_stop_position << "\t";
	ss << mip->scan_target_sequence << "\t";
	ss << mip->mip_seq << "\t";
	ss << feature->start_position - 1 << "\t"; 
	ss << feature->stop_position << "\t";
	ss << mip->strand << "\t";
	ss << mip->mapping_failed << mip->snp_failed << mip->masking_failed << "\t";
	ss << feature->label << "_" << setw(4) << setfill('0') << mip_index << (mip->snp_count == 1 ? ("_SNP_" + string(minor ? "b" : "a")) : "") << "\n";
	return ss.str();
}
// uses bwa to check that targets map uniquely
string check_copy_numbers()
{
	ofstream BWAFQ ((project_name + ".all_sequences.fq").c_str());
	ofstream ARMSFQ ((project_name + ".oligo_copy_count.fq").c_str());
	if (!(BWAFQ.is_open()) || !(ARMSFQ.is_open())) 
	{
		cerr << "[mipgen] copy_check files could not be opened" << endl;
		return "";
	}
	Featurev5 * feature;
	for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		feature = &*it;
		string chr = feature->chr;
		for (int capture_size = max_capture_size; capture_size >= min_capture_size; capture_size -= capture_increment)
		{
			int current_mip_start = feature->start_position_flanked - capture_size;
			while (current_mip_start < feature->stop_position_flanked)
			{
				if(current_mip_start > 0 && current_mip_start + capture_size - 1 <= feature->chromosomal_sequence_stop_position)
				{
					BWAFQ << "@" << capture_size << "_" << feature->chr << "_" << current_mip_start << endl;
					BWAFQ << feature->chromosomal_sequence.substr(current_mip_start - feature->chromosomal_sequence_start_position, capture_size) << endl;
					BWAFQ << "+\n";
					for (int i = 0; i < capture_size; i++) BWAFQ << "#";
					BWAFQ << endl;
				}
				current_mip_start++;
			}
		}
		for (set<int>::iterator it = oligo_sizes.begin(); it != oligo_sizes.end(); it++)
		{
			int oligo_size = *it;
			for(unsigned int relative_start_position = 0; relative_start_position < feature->chromosomal_sequence.length() - oligo_size; relative_start_position++)
			{
				ARMSFQ << "@chr" << chr << ":" << (feature->chromosomal_sequence_start_position + relative_start_position) << "-" << (feature->chromosomal_sequence_start_position + relative_start_position + oligo_size - 1) << endl;
				ARMSFQ << feature->chromosomal_sequence.substr(relative_start_position, oligo_size) << endl;
				ARMSFQ << "+\n";
							for (int i=0; i < oligo_size; i++) ARMSFQ << "#";
							ARMSFQ << endl;
			}
		}
	}
	BWAFQ.close();
	ARMSFQ.close();
	system((args["-bwa"] + " aln -t " + args["-bwa_threads"] + " " + bwa_genome_index + " " + project_name + ".all_sequences.fq > " + project_name + ".all_sequences.sai").c_str());
	system((args["-bwa"] + " samse " + bwa_genome_index + " " + project_name + ".all_sequences.sai " + project_name +".all_sequences.fq > " + project_name + ".all_sequences.sam").c_str());

	int bad_site_counter = 0;
	ifstream CAPTURESAM ((project_name + ".all_sequences.sam").c_str());
	if (CAPTURESAM.is_open()) cerr << "[mipgen] sam file opened" << endl;
	else return "";
	string line;
	bool header = true;
	while (CAPTURESAM.good())
	{
		getline(CAPTURESAM, line);
		if(header && line.find("X0:i:", 0) >= line.length() - 1) continue; 	
		else if (header) header = false;
		unsigned int line_len = line.length();
		if(line_len < 2) continue;
		if(line.find("X0:i:1",0) >= line_len - 1 || line.find("X1:i:0",0) >= line_len - 1) // first regex: ignore header lines. second regex: reject non-unique mappings. third regex: reject mappings that are non-unique with 1 bp edits 
		{
			string size_chr_and_pos = line.substr(0, line.find("\t"));
			size_t sep1 = size_chr_and_pos.find("_", 0);
			size_t sep2 = size_chr_and_pos.rfind("_");
			int size = boost::lexical_cast<int>(size_chr_and_pos.substr(0, sep1));
			string chr = size_chr_and_pos.substr(sep1 + 1, sep2 - sep1 - 1);
			int pos = boost::lexical_cast<int>(size_chr_and_pos.substr(sep2 + 1));
			unmappable_positions[size][chr].insert(pos);
			bad_site_counter++;
		}
	}
	CAPTURESAM.close();
	stringstream i;
	i << bad_site_counter;
	return (i.str() + " ambiguously mapping start positions must be avoided\n");
}
// acquires, prints and load snp data into memory
void load_snps()
{
	string query = "";
	int first_position = features_to_scan.front().start_position_flanked - 500;
	int last_position = features_to_scan.front().stop_position_flanked + 500;
	string current_chr = features_to_scan.front().chr;
	string previous_chr;
	string start_str;
	string stop_str;
	Featurev5 * feature;
	for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		feature = &*it;
		previous_chr = current_chr; 
		current_chr = feature->chr;
		int access_start = feature->start_position_flanked - 500;
		int access_stop = feature->stop_position_flanked + 500;
		if (previous_chr == current_chr && access_stop >= last_position && access_start <= last_position)
		{
			last_position = access_stop;
		}
		else
		{
			stringstream start_ss;
			stringstream stop_ss;
			start_ss << first_position;
			start_str = start_ss.str();
			stop_ss << last_position;
			stop_str = stop_ss.str();
			query.append(previous_chr + ":" + start_str + "-" + stop_str + " ");
			//cout << query << endl;
			first_position = access_start;
			last_position = access_stop;
		}
	}
	
	stringstream start_ss;
	stringstream stop_ss;
	start_ss << first_position;
	stop_ss << last_position;
	start_str = start_ss.str();
	stop_str = stop_ss.str();
	query.append( current_chr + ":" + start_str + "-" + stop_str + " ");
	//PROGRESS << "tabix query: " << endl << query << endl;
	int tabix_test = system(args["-tabix"].c_str());
	if(tabix_test != 256)
	{
		cerr << "[mipgen] tabix not loaded" << endl;
		throw 16;
	}
    if(args["-snp_file"] != "<none>")
    {
        system((args["-tabix"] + " " + args["-snp_file"] + " " + query + " > " + project_name + ".local_snp_data.vcf").c_str()) ;
        parse_vcf(project_name + ".local_snp_data.vcf");
    }
}
    
// reads a VCF and stores SNP positions in memory for evasion purposes
void parse_vcf(string vcf_file)
{
    string line;
    ifstream SNPS (vcf_file.c_str());
    if(SNPS.is_open())
    {
        while (SNPS.good())
        {
            getline(SNPS, line);
            if(line.length() < 2 or line[0] == '#') continue;
            int start_position = 0;
            int stop_position = line.find_first_of(" \t", 0);
            string chr = line.substr(start_position, stop_position);
            start_position = stop_position + 1;
            stop_position = line.find_first_of(" \t", start_position);
            //cout << line.substr(start_position, stop_position - start_position) << endl;
            int position = boost::lexical_cast<int>(line.substr(start_position, stop_position - start_position));
            start_position = stop_position + 1;
            stop_position = line.find_first_of(" \t", start_position);
            start_position = stop_position + 1;
            stop_position = line.find_first_of(" \t", start_position);
            string ref_allele = line.substr(start_position, stop_position - start_position);
            start_position = stop_position + 1;
            stop_position = line.find_first_of(" \t", start_position);
            string alt_allele = line.substr(start_position, stop_position - start_position);
            string alleles_content = ref_allele + alt_allele;
            if (ref_allele.length() > 1) // indel
            {
                for (unsigned int i = 1; i < ref_allele.length(); i++)
                {
                    chr_snp_positions[chr][position + i] = alleles_content;
                }
            }
            else
            {
                chr_snp_positions[chr][position] = alleles_content;
            }
            snp_load_count++;
        }
        SNPS.close();
    }
    else
    {
        cerr << "[mipgen] VCF file could not be opened" << endl;
    }
}

// loads feature details into memory
string get_features_to_scan()
{
	string details;
	list<string> regions_to_scan_contents;
	ifstream FH (regions_to_scan.c_str());
	if(FH.is_open()) cerr << "[mipgen] success on opening region file" << endl;
	else return "";
	while(FH.good())
	{
		getline(FH,details);
		boost::trim(details);
	//	cout << details << endl;
		if (details.length() > 1 && details[0] != '#')
		{
			regions_to_scan_contents.push_back(details);
		}
	}
	FH.close();
	regions_to_scan_contents.sort(compare_regions_to_scan);
	vector<string> bed_fields;
	string default_label = project_name;
	if(default_label.rfind("/") != (unsigned int) -1)
	{
		default_label = default_label.substr(default_label.rfind("/") + 1);
	}
	for (list<string>::iterator it = regions_to_scan_contents.begin(); it != regions_to_scan_contents.end(); it++)
	{
		string line (*it);
		string label;
		if (line[0] == '>' || line[0] == '#' || line.find_first_not_of(" \t\n") > line.length()) // header
			continue;
		else
		{
			boost::split(bed_fields, line, boost::is_any_of(" \t"), boost::token_compress_on);
			label = bed_fields.size() > 3 ? label = bed_fields.at(3) : default_label;
			string chr;
			if (bed_fields.at(0).substr(0,3) == "chr") chr = bed_fields.at(0).substr(3);
			else chr = bed_fields.at(0);
			if(!(features_to_scan.empty()) && features_to_scan.back().chr == chr && boost::lexical_cast<int>(bed_fields.at(1)) - features_to_scan.back().stop_position - 2 * feature_flank < min_capture_size / 2)
			{
				int second_stop = boost::lexical_cast<int>(bed_fields.at(2));
				int later_stop = second_stop > features_to_scan.back().stop_position ? second_stop : features_to_scan.back().stop_position;
				Featurev5 replacement_feature (chr, features_to_scan.back().start_position, later_stop, feature_flank, label);
				features_to_scan.pop_back();
				features_to_scan.push_back(replacement_feature);
			}
			else
			{
				Featurev5 current_feature (chr, boost::lexical_cast<int>(bed_fields.at(1)) + 1,boost::lexical_cast<int>(bed_fields.at(2)), feature_flank, label);
				features_to_scan.push_back(current_feature);
			}
		}
	}
	

	/*for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		cout << it->chr << "\t" << it->start_position << "\t";
		cout << it->stop_position << endl;
	}
	*/
	return "successfully loaded features for mip design; retrieving chromosomal sequence\n";
}
// uses Tandem Repeats Finder to mask tandem repeats
bool get_masked_features_to_scan()
{
	if(!(args["-trf"] == "off"))
	{
		system((args["-trf"] + " " + project_name + ".feature_sequences.fa 2 7 7 80 10 14 100 -m -h").c_str());
	}
	size_t prefix_start = project_name.rfind("/");
	string prefix = project_name.substr(prefix_start + 1);
	
	ifstream MASKEDFEATURES ((prefix + ".feature_sequences.fa.2.7.7.80.10.14.100.mask").c_str());
	if(!(MASKEDFEATURES.is_open()) || args["-trf"]=="off")
	{
		cerr << "[mipgen] no masked feature file found" << endl;
		for(list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
		{
			it->masked_chromosomal_sequence = (*it).chromosomal_sequence;
		}
		return args["-trf"]=="off";
	}
	string line;
	list<Featurev5>::iterator it = features_to_scan.begin();
	Featurev5 * feature;
	while(MASKEDFEATURES.good())
	{
		getline(MASKEDFEATURES, line);
		if(line.length() ==0)
		{
			continue;
		}
		if(line.at(0) == '>')
		{
			feature = &*it;
			it++;
			continue;
		}
		feature->masked_chromosomal_sequence = feature->masked_chromosomal_sequence + line;
	}
	MASKEDFEATURES.close();
	return true;
}
		
// pulls out sequences with faidx		
bool get_chr_fasta_sequence_using_samtools ()
{
	string chr ("0");
	string samtools_query;
	Featurev5 * feature;
	system(("rm -f " + project_name + ".feature_sequences.fa").c_str());
	system(("rm -f " + project_name + ".flanking_sequences.fa").c_str());
	string header;
	string prefix;
	try{
		ifstream GENOMEREF (bwa_genome_index.c_str());
		getline(GENOMEREF, header);
		GENOMEREF.close();
		if(header.length() <= 3)
		{
			prefix = "";
		}
		else
		{
			prefix = header.substr(1,3) == "chr" ? "chr" : ""; //adds "chr" to samtools queries if present in reference being used
		}
		cerr << "[mipgen] first line of reference fasta reads: " << header << endl;
	}
	catch (...) {
		cerr << "genome fasta could not be opened; check file path and permissions?" << endl;
                return false;
	}
	for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		feature = &*it;
		chr = feature->chr;
		feature->chromosomal_sequence_start_position = feature->start_position_flanked - max_capture_size;
		feature->chromosomal_sequence_stop_position = feature->stop_position_flanked + max_capture_size + 14;
		samtools_query = "samtools faidx " + bwa_genome_index + " " + prefix + chr + ":" \
			+ boost::lexical_cast<string>(feature->start_position_flanked - max_capture_size) + "-" \
			+ boost::lexical_cast<string>(feature->stop_position_flanked + max_capture_size + 14) + " " \
			+ ">> " + project_name + ".feature_sequences.fa";
		system(samtools_query.c_str());
		samtools_query = "samtools faidx " + bwa_genome_index + " " + prefix + chr + ":" \
			+ boost::lexical_cast<string>(feature->start_position_flanked - max_capture_size - 1000) + "-" \
			+ boost::lexical_cast<string>(feature->stop_position_flanked + max_capture_size + 14 + 1000) + " " \
			+ ">> " + project_name + ".flanking_sequences.fa";
		system(samtools_query.c_str());
	}
	ifstream FEATURESEQS ((project_name + ".feature_sequences.fa").c_str());
	if (!(FEATURESEQS.is_open())) 
	{
		cerr << "[mipgen] fasta file could not be opened" << endl;
		return false;
	}
	ifstream FLANKINGSEQS ((project_name + ".flanking_sequences.fa").c_str());
	if (!(FLANKINGSEQS.is_open())) 
	{
		cerr << "[mipgen] fasta file could not be opened" << endl;
		return false;
	}
	string line;
	string long_range;
	getline(FEATURESEQS, line);
	getline(FLANKINGSEQS, line);
	for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		long_range="";
		feature = &*it;
		while(FEATURESEQS.good())
		{
			getline(FEATURESEQS, line);
			if(line[0] == '>')
			{
				break;
			}
			for (unsigned int i = 0; i < line.length(); i++) line[i] = toupper(line[i]);
			feature->chromosomal_sequence += line;
		}
		while(FLANKINGSEQS.good())
		{ 
			getline(FLANKINGSEQS, line);
			if(line[0] == '>')
			{
				break;
			}
			for (unsigned int i = 0; i < line.length(); i++) line[i] = toupper(line[i]);
			long_range += line;
		}
		if(args["-score_method"] != "logistic")
			feature->get_long_range_content(long_range, feature_mers);
	}
	FEATURESEQS.close();
	FLANKINGSEQS.close();
	return true;
}		
	
// loads chromosomes into memory one by one
bool get_chr_fasta_sequence_from_genome_dir (string genome_dir)
{
	string chr ("0");
	string chr_fasta;
	string chr_sequence;
	Featurev5 * feature;
	ofstream FEATURESEQS ((project_name + ".feature_sequences.fa").c_str());
	for (list<Featurev5>::iterator it = features_to_scan.begin(); it != features_to_scan.end(); it++)
	{
		feature = &*it;
		if (!(chr == feature->chr))
		{
			chr_sequence = "";
			chr = feature->chr;
			chr_fasta = genome_dir + "/chr" + chr + ".fa";
			//cout << chr_fasta << endl;
			ifstream FH (chr_fasta.c_str());
			if (!(FH.is_open())) 
			{
				cerr << "[mipgen] fasta file could not be opened" << endl;
				return false;
			}
			string line;
			while(FH.good())
			{
				getline(FH, line);
				if (line[0]=='>') continue;
				else {
					for (unsigned int i = 0; i < line.length(); i++) line[i] = toupper(line[i]);
					chr_sequence.append(line);
				}
			}
			FH.close();
		}
		int chr_start_coordinate = feature->start_position_flanked - max_capture_size < 1 ? 1 : feature->start_position_flanked - max_capture_size;
		int chr_stop_coordinate = feature->stop_position_flanked + max_capture_size + 15 > (signed int) chr_sequence.length() ? chr_sequence.length() : feature->stop_position_flanked + max_capture_size + 15; // 15 is the maximum difference between arm lengths, which extends the bases needed slightly
		int feature_length = chr_stop_coordinate - chr_start_coordinate + 1; 

		feature->chromosomal_sequence = chr_sequence.substr(chr_start_coordinate - 1, feature_length); // subtract one to convert from chr position to string position
		feature->chromosomal_sequence_start_position = chr_start_coordinate;
		feature->chromosomal_sequence_stop_position = chr_stop_coordinate;
		
		FEATURESEQS << ">" << feature->chr << ':' << feature->chromosomal_sequence_start_position << "-" << feature->chromosomal_sequence_stop_position << endl;
		FEATURESEQS << feature->chromosomal_sequence << endl;
		if(args["-score_method"] != "logistic")
			feature->get_long_range_content(chr_sequence.substr(feature->start_position_flanked - max_capture_size - 1 - 1000, feature_length + 2000), feature_mers);
	}
	FEATURESEQS.close();
	return true;
}
// prints out regions that have not been tiled in BED format
void print_gaps (ofstream & gap_bed, string file_ext, string progress_note, Featurev5 * feature, set<int> & positions)
{
	if (!(positions.empty()))
	{
		if(!(gap_bed.is_open()))
		{	
			gap_bed.open((args["-project_name"] + file_ext).c_str());
		}
		PROGRESS << progress_note << feature->chr << ":\n";
		int start = *positions.begin();
		int stop = start - 1;
		for (set<int>::iterator it = positions.begin(); it != positions.end(); it++)
		{
			if (*it == stop + 1)
			{
				stop++;
			}
			else
			{
                bad_design_count++;
				gap_bed << feature->chr << "\t" << start - 1 << "\t" << stop << endl;
				start = *it;
				stop = *it;
			}
		}
        bad_design_count++;
		gap_bed << feature->chr << "\t" << start - 1 << "\t" << stop << endl;
	}
}
// prints out regions that have not been tiled in BED format
void create_gap (ofstream & gap_bed, string file_ext, string progress_note, Featurev5 * feature, set<int> & positions)
{
    bad_design_count++;
	if(!(gap_bed.is_open()))
	{	
		gap_bed.open((args["-project_name"] + file_ext).c_str());
	}
	PROGRESS << progress_note << feature->chr << ":\n";
	int start = *positions.begin();
	int stop = start + max_capture_size / 2;
    PROGRESS << feature->chr << "\t" << start - 1 << "\t" << stop << endl;
	for(int i = start; i <= stop; i++)
	{
		positions.erase(i);
	}
	gap_bed << feature->chr << "\t" << start - 1 << "\t" << stop << endl;

}
// intializes all the relevant variables
string parse_command_line(int argc, char * argv[]) 
{
	string all_options = "-regions_to_scan -feature_flank -feature_flank -genome_dir -project_name -bwa -bwa_genome_index \
		-trf -tabix -check_copy_number -common_snps -arm_lengths -arm_length_sums -min_capture_size -download_tabix_index \
		-max_capture_size -capture_increment -max_mip_overlap -starting_mip_overlap -stop_optimizing_scores_above -silent_mode \
		-masked_arm_threshold -seal_both_strands -half_seal_both_strands -tag_sizes -ext_min_length -lig_min_length -bwa_threads \
		-snp_file -double_tile_strand_unaware -double_tile_strands_separately -score_method -logistic_heuristic -file_of_parameters \
		-logistic_priority_score -svr_priority_score -logistic_optimal_score -svr_optimal_score -max_arm_copy_product -target_arm_copy";
	vector<string> option_vector;
	boost::split(option_vector, all_options, boost::is_any_of(" \t"), boost::token_compress_on);
	string splash = "\n\
usage: mipgen (<-parameter_name> <parameter_value>)*\n\
Created by Evan Boyle (boylee@uw.edu)\n\
mipgen -doc for optional parameter descriptions\n\
Required parameters:\n\
\n\
-project_name                   prefix for output\n\
                                ex: /my/output/directory/batch_1\n\
-bwa_genome_index               samtools- and bwa-indexed genome reference file path\n\
-regions_to_scan                BED file of positions GUARANTEED to be captured\n\
                                regions will be merged as needed\n\
-min_capture_size               integer value for length of targeting arms plus insert region to be captured\n\
                                tested down to 120\n\
-max_capture_size               integer value for length of targeting arms plus insert region to be captured\n\
                                tested up to 250\n\
Highly recommended parameter:\n\
\n\
-snp_file			VCF file, either for the whole genome or just select regions\n\
				You must index with tabix in advance so that relevant sequences can be retrieved\n\
				NCBI (at time of writing) has a file of common SNPs for human here:\n\
				ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/common_all.vcf.gz\n";
	string doc = "\n\
usage: mipgen (<-parameter_name> <parameter_value>)*\n\
Created by Evan Boyle (boylee@uw.edu)\n\
Required parameters:\n\
\n\
-project_name                   prefix for output\n\
                                ex: /my/output/directory/batch_1\n\
-bwa_genome_index               samtools- and bwa-indexed genome reference file path\n\
-regions_to_scan                BED file of positions GUARANTEED to be captured\n\
                                regions will be merged as needed\n\
-min_capture_size               integer value for length of targeting arms plus insert region to be captured\n\
                                tested down to 120\n\
-max_capture_size               integer value for length of targeting arms plus insert region to be captured\n\
                                tested up to 250\n\
Oligo control:\n\
\n\
-arm_lengths                    manual setting of individual possible ligation and extension arm lengths\n\
                                format: <extension arm length 1>:<ligation arm length 1>,<extension arm length 2>:<ligation arm length 2>,...\n\
                                ex: 16:24,16:25,16:26\n\
-arm_length_sums                setting of arm lengths to include all pairs of arm lengths adding to a certain length\n\
                                default is 40,41,42,43,44,45\n\
                                format: <sum of extension and ligation arm lengths 1>,<sum of extension and ligation arm lengths 2>,...\n\
                                ex: 40,45\n\
-ext_min_length                 minimum length of extension arm\n\
                                default is 16\n\
-lig_min_length                 minimum length of ligation arm\n\
                                default is 18\n\
-tag_sizes                      specifies degenerate tag (Ns) to place in extension arm and ligation arm\n\
                                more than 8 bases severely reduces on-target rate of capture\n\
                                default is 5,0\n\
                                format: <extension arm tag size>,<ligation arm tag size>\n\
                                ex: 4,4\n\
-masked_arm_threshold           fraction of targeting arms having masked bases (floating point number from 0 to 1)\n\
                                to tolerate before optimizing with respect to masking (requires TRF executable)\n\
                                default is 0.5\n\
-target_arm_copy                threshold over which minimizing copy number to the reference takes priority\n\
                                default is 20\n\
-max_arm_copy_product           maximum permissible product of targeting arm copy number to the reference\n\
                                default is 75\n\
Tool dependencies:\n\
\n\
-tabix                          tabix build\n\
                                default is \"tabix\"\n\
-bwa                            bwa build\n\
                                default is \"bwa\"\n\
-trf                            tandem repeats finder executable\n\
                                default is \"off\"\n\
                                providing a TRF executable enables filtering of simple repeats in probe arms\n\
Input options:\n\
-genome_dir                     genome reference directory helpful for large jobs\n\
                                ex: /my/genome/split/by/chromosomes/\n\
-snp_file                       path to vcf file of snps to avoid during design\n\
-common_snps                    providing \"off\" will disable loading of common SNPs to avoid from NCBI using Tabix\n\
				DEPRECATED:ONLY SNP FILE OPTION NOW SUPPORTED\n\
-file_of_parameters             file containing any of the above parameters in the following format:\n\
                                        -first_parameter_name first_parameter_value\n\
                                        -second_parameter_name second_parameter_value\n\
                                        <etc.>\n\
                                        lines that do not start with '-' are ignored\n\
Tiling control:\n\
\n\
-feature_flank                  integer value for extending BED coordinates\n\
                                default is 0\n\
-capture_increment              integer step size for examining capture sizes\n\
                                default is 5\n\
                                starts at maximum and descends until minimum is reached or\n\
                                acceptable score found, whichever is first\n\
-logistic_heuristic             providing \"off\" will remove assumptions to reduce search space\n\
-max_mip_overlap                integer specifying the maximum number of nucleotides of overlap\n\
                                to test before selecting a MIP\n\
                                default is 30\n\
-starting_mip_overlap           integer specifying the starting number of nucleotides of overlap\n\
                                to use when linearly tiling with MIPs\n\
                                default is 0\n\
-check_copy_number              providing \"off\" enables targeting of non-unique sites in the reference\n\
-seal_both_strands              providing \"on\" disallows targeting arms occupying the same base\n\
                                positions on opposite strands\n\
-half_seal_both_strands         providing \"on\" disallows a second targeting arm from occupying more than half\n\
                                of the positions of a targeting arm placed on the other strand\n\
-double_tile_strand_unaware     providing \"on\" will ensure that all positions in targeted regions will be tiled twice\n\
-double_tile_strands_separately providing \"on\" will ensure that all positions in targeted regions will be tiled twice,\n\
                                at least once on each strand\n\
Scoring parameters:\n\
\n\
-score_method                   \"logistic\" will use the logistic model for scoring (default)\n\
                                \"svr\" will switch behavior to score with the svr model\n\
                                \"mixed\" will first filter with logistic scoring and finalize with the svr model\n\
-logistic_optimal_score         value of score metric at which to stop optimizing MIP score and\n\
                                accept current MIP, lower scores lead to fewer outliers\n\
                                default is 0.98\n\
-svr_optimal_score              analogous to logistic_optimal_score\n\
                                default is 2.2\n\
-logistic_priority_score        value of score metric below which MIPs are placed nonlinearly to\n\
                                optimize score and SNPs are tolerated\n\
                                default is 0.9,lower values enable more efficient linear tiling\n\
-svr_priotity_score             analogous to logistic_priority_score\n\
                                default is 1.5 for SVR\n\
Miscellaneous:\n\
\n\
-silent_mode                    providing \"on\" will reduce volume of text output\n\
-download_tabix_index           providing \"on\" will force redownload of common snp tbi file\n\
                                DEPRECATED\n\
-bwa_threads                    make use of BWA's multithreading option (-t _)\n\
                                default is 1\n";

	if(argc == 1) return splash;
	else if(string(argv[1])=="-doc") return doc;
	bool check_file = false;
	bool option_found;
	for ( int i = 0; i < argc; i++)
	{
		if (argv[i][0] == '-')
		{
			string parameter (argv[i]);
			size_t j = 0;
			option_found = false;
			while(!(option_found))
			{
				if(j == option_vector.size()) {cerr << "not found" << endl; break;}
				else if(option_vector.at(j) == parameter) option_found = true;
				j++;
			}
			if(!(option_found))
			{
				return parameter + " not recognized as valid option\n";
			}
				
			string value (argv[i + 1]);
			args[parameter] = value;
			if (parameter == "-file_of_parameters") check_file = true;
			else if (parameter == "-arm_length_sums") has_arm_length_sum_option = true; 
			else if (parameter == "-arm_lengths") has_arm_length_option = true;
		}
	}
	if (check_file)
	{
		ifstream PARAMETERS (args["-file_of_parameters"].c_str());
		if (PARAMETERS.is_open()==false) 
		{
			cerr << splash << endl;
			return "file of parameters could not be opened\n";
		}
		string line;
		while(PARAMETERS.good())
		{
			getline(PARAMETERS,line);
			if (line[0] == '-')
			{
				int boundary = line.find_first_of(" ", 0);
				size_t j = 0;
				string parameter = line.substr(0,boundary);
				option_found = false;
				while(!(option_found))
				{
					if(j == option_vector.size()) break;
					if(option_vector.at(j) == parameter) option_found = true;
					j++;
				}
				if(!(option_found))
				{
					return parameter + " not recognized as valid option\n";
				}
				string value = line.substr(boundary + 1,line.length()); 
				args[parameter] = value;
				if (parameter == "-arm_length_sums") has_arm_length_sum_option = true;
				else if (parameter == "-arm_lengths") has_arm_length_option = true;
			}
		}
		PARAMETERS.close();
	}
	//for (map<string,string>::iterator it = args.begin(); it!=args.end(); it++)
	//{ cout << it->first << "\t" << it->second << endl; }
	string necessary_parameters [] = 
	{	
		"-regions_to_scan",
		"-project_name",
		"-max_capture_size",
		"-min_capture_size",
		"-bwa_genome_index"
	};
	for (int i = 0; i < 5; i++)
	{
		if (args.find(necessary_parameters[i]) == args.end())
		{
			cerr << splash << endl;
			return "required parameter " + necessary_parameters[i] + " not found";
		}
	}
	//successful parse
	return "";
};



// carries out the actual selection process of final MIPs
void pick_mips (Featurev5 * feature, vector<double> & scoring_parameters, svm_model * model)
{
	set<int> positions_to_scan;
	set<int> positions_to_scan_again;
	set<int> positions_to_scan_minus;
	set<int> positions_to_scan_minus_again;
	int strand_to_use = double_tile_strands_separately ? 0 : -1;
	bool extended_region;

	for(int position = feature->start_position_flanked; position <= feature->stop_position_flanked; position++)
	{
		positions_to_scan.insert(position);
		if (double_tile) positions_to_scan_again.insert(position);
		if (double_tile_strands_separately) positions_to_scan_minus.insert(position);
		if (double_tile && double_tile_strands_separately) positions_to_scan_minus_again.insert(position);
	}
	boost::shared_ptr<SVMipv4> picked_mip = optimize_worst_in_region(feature, positions_to_scan, strand_to_use);
	if(args["-score_method"] == "mixed" && picked_mip != 0)
	{
		picked_mip->get_parameters(scoring_parameters, feature->long_range_content);
		picked_mip->score = predict_value(scoring_parameters, model);
	}
	//if (picked_mip == 0) { cout << "collapsed not defined" << endl; }
	while (!(positions_to_scan.empty()) && picked_mip != 0 && picked_mip->score < lower_score_limit)
	{
		manage_picked_mip(feature, picked_mip, positions_to_scan);
		picked_mip = optimize_worst_in_region(feature, positions_to_scan, strand_to_use);
		if(args["-score_method"] == "mixed" && picked_mip != 0)
		{
			picked_mip->get_parameters(scoring_parameters, feature->long_range_content);
			picked_mip->score = predict_value(scoring_parameters, model);
		}
	}
	if(double_tile_strands_separately) 
	{
		picked_mip = optimize_worst_in_region(feature, positions_to_scan_minus, 1);
		while (!(positions_to_scan_minus.empty()) && picked_mip != 0 && picked_mip->score < lower_score_limit)
        {
            manage_picked_mip(feature, picked_mip, positions_to_scan_minus);
            picked_mip = optimize_worst_in_region(feature, positions_to_scan_minus, 1);
            if(args["-score_method"] == "mixed" && picked_mip != 0)
			{
				picked_mip->get_parameters(scoring_parameters, feature->long_range_content);
				picked_mip->score = predict_value(scoring_parameters, model);
			}	
        }
	}
	// cout << "optimizing done" << endl;
	if (!(positions_to_scan.empty()))
	{
		do 
		{
            // cout << "picking MIP at " << *positions_to_scan.begin() << " to " << *positions_to_scan.rbegin() << endl;
			picked_mip = translocate_down_region(feature, positions_to_scan, scoring_parameters, model, strand_to_use);
			extended_region = *positions_to_scan.rbegin() - *positions_to_scan.begin() > max_capture_size;
            // cout << "extended region? " << extended_region << endl;
			if(picked_mip == null_pointer && extended_region)
			{
				create_gap(GAPS, ".coverage_failed.bed", "GAP INTRODUCED ON CHROMOSOME ", feature, positions_to_scan);
			}
			if (picked_mip != 0) manage_picked_mip(feature, picked_mip, positions_to_scan);
        } while (!(positions_to_scan.empty()) && (picked_mip != 0 || extended_region));
    }
	if (!(positions_to_scan_minus.empty()))
	{
		do
		{
			picked_mip = translocate_down_region(feature, positions_to_scan_minus, scoring_parameters, model, 1);
			extended_region = *positions_to_scan_minus.rbegin() - *positions_to_scan_minus.begin() > max_capture_size;
			if(picked_mip == null_pointer && extended_region)
			{
				create_gap(MINUSGAPS, ".minus_strand_failed.bed", "GAP INTRODUCED ON MINUS STRAND OF CHROMOSOME ", feature, positions_to_scan_minus);
			}
			if (picked_mip != 0) manage_picked_mip(feature, picked_mip, positions_to_scan_minus);
		} while (!(positions_to_scan_minus.empty()) && (picked_mip != 0 || extended_region));
	}
	if (double_tile)
	{
		do
		{
			picked_mip = translocate_down_region(feature,positions_to_scan_again, scoring_parameters, model, strand_to_use);
			extended_region = *positions_to_scan_again.rbegin() - *positions_to_scan_again.begin() > max_capture_size;
			if(picked_mip == null_pointer && extended_region)
			{
				create_gap(DOUBLEGAPS, ".double_tile_failed.bed", "GAP INTRODUCED ON DOUBLE TILING OF CHROMOSOME ", feature, positions_to_scan_again);
			}
			//cout << "iteration, picked mip = " << picked_mip << endl;
			if (picked_mip != 0) manage_picked_mip(feature, picked_mip, positions_to_scan_again);
		} while (!(positions_to_scan_again.empty()) && (picked_mip != 0 || extended_region));
		if (double_tile_strands_separately)
		{
			do
			{
				picked_mip = translocate_down_region(feature,positions_to_scan_again, scoring_parameters, model, 1);
				extended_region = *positions_to_scan_again.rbegin() - *positions_to_scan_again.begin() > max_capture_size;
				if(picked_mip == null_pointer && extended_region)
				{
					create_gap(DOUBLEGAPS, ".minus_strand_double_tile_failed.bed", "GAP INTRODUCED ON MINUS STRAND OF DOUBLE TILING OF CHROMOSOME ", feature, positions_to_scan_minus_again);
				}
				//cout << "iteration, picked mip = " << picked_mip << endl;
				if (picked_mip != 0) manage_picked_mip(feature, picked_mip, positions_to_scan_minus_again);
			} while (!(positions_to_scan_minus_again.empty()) && (picked_mip != 0 || extended_region));
		}
	}
	print_gaps(GAPS, ".coverage_failed.bed", "BASES NOT COVERED ON CHROMOSOME ", feature, positions_to_scan);
	print_gaps(DOUBLEGAPS, ".double_tile_failed.bed", "BASES NOT DOUBLE TILED ON CHROMOSOME ", feature, positions_to_scan_again);
	print_gaps(MINUSGAPS, ".minus_strand_failed.bed", "BASES NOT COVERED ON MINUS STRAND OF CHROMOSOME ", feature, positions_to_scan_minus);
	print_gaps(DOUBLEMINUSGAPS, ".minus_strand_double_tile_failed.bed", "BASES NOT DOUBLE TILED ON MINUS STRAND OF CHROMOSOME ", feature, positions_to_scan_minus_again);
}

// records each position's highest scoring mip possible 
void collapse_mips () // throws out all restriction and mapping failures. retains snp failures if there is no alternative 
{
	PROGRESS << "collapsing feature #" << feature_counter << endl;
	for( map<int, map<string, boost::shared_ptr<SVMipv4> > >::iterator it = scan_strand_best_mip.begin(); it != scan_strand_best_mip.end(); it++)
	{
		int scan = it->first;
		for (map<string, boost::shared_ptr<SVMipv4> >::iterator it2 = scan_strand_best_mip[scan].begin(); it2 != scan_strand_best_mip[scan].end(); it2++)
		{
			boost::shared_ptr<SVMipv4> current_mip;
			string current_strand = it2->first;
			current_mip = it2->second;
			if (current_mip->ext_probe_copy * current_mip->lig_probe_copy > max_arm_copy || current_mip->ext_probe_copy > target_arm_copy || current_mip->lig_probe_copy > target_arm_copy) continue;
			if (current_mip->arm_fraction_masked > boost::lexical_cast<double>(args["-masked_arm_threshold"])) continue;
			for (int interrogated_position = current_mip->scan_start_position; interrogated_position <= current_mip->scan_stop_position; interrogated_position++)
			{
				string strand = it2->first;
				current_mip = it2->second;
				if (pos_strand_best_mip.empty() || pos_strand_best_mip[interrogated_position].find(strand) == pos_strand_best_mip[interrogated_position].end())
				{
					pos_strand_best_mip[interrogated_position][strand] = current_mip;
				}
				else if (current_mip->snp_count < pos_strand_best_mip[interrogated_position][strand]->snp_count)
				{	
					pos_strand_best_mip[interrogated_position][strand] = current_mip;
				}
				else if (current_mip->score > pos_strand_best_mip[interrogated_position][strand]->score && current_mip->snp_count == pos_strand_best_mip[interrogated_position][strand]->snp_count)
				{
					pos_strand_best_mip[interrogated_position][strand] = current_mip;
				}
			}
		}
	}
};
// prints out the collapsed mips for reference
void output_collapsed_mips (Featurev5 * feature)
{
	list<int> positions_collapsed;
	for (map<int, map<string, boost::shared_ptr<SVMipv4> > >::iterator it = pos_strand_best_mip.begin(); it != pos_strand_best_mip.end(); it++)
	{
		positions_collapsed.push_back(it->first);
	}
	positions_collapsed.sort();
	for (list<int>::iterator it = positions_collapsed.begin(); it != positions_collapsed.end(); it++)
	{
		int position = *it;
		for (map<string, boost::shared_ptr<SVMipv4> >::iterator it2 = pos_strand_best_mip[position].begin(); it2 != pos_strand_best_mip[position].end(); it2++)
		{
			collapsed_mip_counter++;
			COLLAPSEDMIPS << print_details(feature, it2->second, collapsed_mip_counter, false);
		}
	}
}
// picks the best mip per scan start				
void condense_mips(Featurev5 * feature)
{
	PROGRESS << "condensing feature #" << feature_counter << endl;
	string strands [] = {"+","-"};
	for(map<int, map<string, list<boost::shared_ptr<SVMipv4> > > >::iterator it = scan_strand_mip_list.begin(); it != scan_strand_mip_list.end(); it++)  
	{
		int position = it->first;
		int chosen_copy_count;
		int current_arm_copy_count;
		double chosen_masked_arm_proportion;
		double current_masked_arm_proportion;
		for (int i = 0; i < 2; i++)
		{
			bool skip_ahead = false;
			string strand = strands[i];
			for (list<boost::shared_ptr<SVMipv4> >::iterator it2 = scan_strand_mip_list[position][strand].begin(); it2 !=  scan_strand_mip_list[position][strand].end(); it2++)
			{
				if (skip_ahead) continue;
				boost::shared_ptr<SVMipv4> current_mip = *it2;
				if (current_mip->ext_probe_copy * current_mip->lig_probe_copy > max_arm_copy) continue;
				if (current_mip->mapping_failed == '0')
				{
					current_arm_copy_count = (current_mip->ext_probe_copy > current_mip->lig_probe_copy) ? current_mip->ext_probe_copy : current_mip->lig_probe_copy;
					current_masked_arm_proportion = current_mip->arm_fraction_masked;
					
					if (scan_strand_best_mip.empty() || scan_strand_best_mip[position].find(strand) == scan_strand_best_mip[position].end())
					{
						scan_strand_best_mip[position][strand] = current_mip;
						chosen_masked_arm_proportion = current_masked_arm_proportion;
						chosen_copy_count = current_arm_copy_count;
					}
					else if (current_masked_arm_proportion > boost::lexical_cast<double>(args["-masked_arm_threshold"]) && current_masked_arm_proportion < chosen_masked_arm_proportion)
					{
						scan_strand_best_mip[position][strand] = current_mip;
						chosen_masked_arm_proportion = current_masked_arm_proportion;
						chosen_copy_count = current_arm_copy_count; 
					}
					else 
					{
						if (current_arm_copy_count > target_arm_copy && current_arm_copy_count < chosen_copy_count)
						{
							scan_strand_best_mip[position][strand] = current_mip;
							chosen_masked_arm_proportion = current_masked_arm_proportion;
							chosen_copy_count = current_arm_copy_count;
						}
						else if (current_arm_copy_count <= target_arm_copy)
						{
							if (current_mip->score < lower_score_limit && current_mip->score > scan_strand_best_mip[position][strand]->score)
							{
								scan_strand_best_mip[position][strand] = current_mip;
								chosen_masked_arm_proportion = current_masked_arm_proportion;
								chosen_copy_count = current_arm_copy_count;
							}
							else if (current_mip->score > lower_score_limit)
							{
								if (current_mip->snp_count < scan_strand_best_mip[position][strand]->snp_count)
								{
									scan_strand_best_mip[position][strand] = current_mip;
									chosen_masked_arm_proportion = current_masked_arm_proportion;
									chosen_copy_count = current_arm_copy_count;
								}
								else if (current_mip->snp_count == scan_strand_best_mip[position][strand]->snp_count)
								{
									if (current_mip->score > scan_strand_best_mip[position][strand]->score)
									{
										scan_strand_best_mip[position][strand] = current_mip;
										if(current_mip->score > upper_score_limit) skip_ahead = true;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}				
// returns the mip with the lowest optimized score for a position
boost::shared_ptr<SVMipv4> optimize_worst_in_region (Featurev5 * feature, set<int> & positions_to_scan, int strand_to_use)
{
	string chr = feature->chr;
	boost::shared_ptr<SVMipv4> mip_for_worst_position;
	bool worst_mip_set = false;
	
	for(set<int>::iterator it = positions_to_scan.begin(); it != positions_to_scan.end(); it++)
	{
		int interrogated_position = *it;
		if (pos_strand_best_mip.find(interrogated_position) == pos_strand_best_mip.end()) continue;
		map<string, boost::shared_ptr<SVMipv4> > mips_at_interrogated_position = pos_strand_best_mip[interrogated_position];	
		boost::shared_ptr<SVMipv4> current_stranded_mip;
		boost::shared_ptr<SVMipv4> current_plus_mip;
		boost::shared_ptr<SVMipv4> current_minus_mip;
		bool minus_mip_set = false;
		bool plus_mip_set = false;

		if (pos_strand_best_mip[interrogated_position].find("+") != pos_strand_best_mip[interrogated_position].end() && strand_to_use != 1)
		{
			plus_mip_set = true;
			current_plus_mip = pos_strand_best_mip[interrogated_position]["+"];
			if (chr_strand_pos_used_arm_bases.find(chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[chr].find("+") != chr_strand_pos_used_arm_bases[chr].end())
			{
				for (int position_to_occupy = current_plus_mip->ext_probe_start; position_to_occupy <= current_plus_mip->ext_probe_stop; position_to_occupy++)
				{
				   	if (chr_strand_pos_used_arm_bases[chr]["+"].count(position_to_occupy) == 1) plus_mip_set = false;
				}
			}
			if (chr_strand_pos_used_arm_bases.find(chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[chr].find("+") != chr_strand_pos_used_arm_bases[chr].end())
			{
				for (int position_to_occupy = current_plus_mip->lig_probe_start; position_to_occupy <= current_plus_mip->lig_probe_stop; position_to_occupy++)
				{
					if (chr_strand_pos_used_arm_bases[chr]["+"].count(position_to_occupy) == 1) plus_mip_set = false;
				}
			}
			if (plus_mip_set) current_stranded_mip = current_plus_mip;
		}
		if (pos_strand_best_mip[interrogated_position].find("-") != pos_strand_best_mip[interrogated_position].end() && strand_to_use != 0)
		{
			current_minus_mip = pos_strand_best_mip[interrogated_position]["-"];
			minus_mip_set = true;
			if (chr_strand_pos_used_arm_bases.find(chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[chr].find("-") != chr_strand_pos_used_arm_bases[chr].end())

			{
				for (int position_to_occupy = current_minus_mip->ext_probe_start; position_to_occupy <= current_minus_mip->ext_probe_stop; position_to_occupy++)
				{
					if (chr_strand_pos_used_arm_bases[chr]["-"].count(position_to_occupy) == 1) minus_mip_set = false;
				}
			}
			if (chr_strand_pos_used_arm_bases.find(chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[chr].find("-") != chr_strand_pos_used_arm_bases[chr].end())
			{
				for (int position_to_occupy = current_minus_mip->lig_probe_start; position_to_occupy <= current_minus_mip->lig_probe_stop; position_to_occupy++)
				{
					if (chr_strand_pos_used_arm_bases[chr]["-"].count(position_to_occupy) == 1) minus_mip_set = false;
				}
			}
			if (minus_mip_set)
			{
				if(plus_mip_set)
					current_stranded_mip = (current_plus_mip->score > current_minus_mip->score) ? current_plus_mip : current_minus_mip;
				else
					current_stranded_mip = current_minus_mip;
			}
		}
		if ((plus_mip_set || minus_mip_set) && (!(worst_mip_set) || current_stranded_mip->score < mip_for_worst_position->score))
		{
			mip_for_worst_position = current_stranded_mip;
			worst_mip_set = true;
		}
	}
	if(worst_mip_set) return mip_for_worst_position;
	else return null_pointer;
}
// linearly tiles down an exon and returns the next mip
boost::shared_ptr<SVMipv4> translocate_down_region (Featurev5 * feature, set<int> & positions_to_scan, vector<double> & scoring_parameters, svm_model * model, int strand_to_use)
{
	int latest_scan_position = *positions_to_scan.begin();
	int positions_to_end = feature->stop_position_flanked - latest_scan_position;
	int min_scan_size = min_capture_size - *arm_length_sum_set.rbegin();
	
	int earliest_scan_position;
	int preliminary_scan_position;
	int direction;
	if (positions_to_end < min_scan_size - 8) // Questionable: how small does the window have to be to switch to this new behavior? Currently within 8 bp of fitting the whole region in one mip
	{
		earliest_scan_position = feature->stop_position_flanked - min_scan_size + 1;
		preliminary_scan_position = earliest_scan_position;
		direction = 1;
	while (scan_strand_best_mip.find(preliminary_scan_position) == scan_strand_best_mip.end() && preliminary_scan_position <= latest_scan_position) preliminary_scan_position++;
	}
	else
	{
		earliest_scan_position = latest_scan_position - max_mip_overlap;
		preliminary_scan_position = latest_scan_position - starting_mip_overlap;
		direction = -1;
		while (scan_strand_best_mip.find(preliminary_scan_position) == scan_strand_best_mip.end() && preliminary_scan_position >= earliest_scan_position) preliminary_scan_position--;
	}
	if (scan_strand_best_mip.find(preliminary_scan_position) == scan_strand_best_mip.end()) 
	{
		return null_pointer;
	}
	string strands[] = {"+","-"};
	boost::shared_ptr<SVMipv4> next_mip;
	int previous_scan_extent = latest_scan_position + 1; // merely for overriding max_mip_overlap when a gap will result
	for (int chosen_scan_position = preliminary_scan_position; (next_mip == 0 && previous_scan_extent > latest_scan_position && chosen_scan_position < feature->stop_position_flanked && chosen_scan_position > feature->start_position_flanked - min_capture_size) || (next_mip != 0 && ((next_mip->score < upper_score_limit || next_mip->snp_count > 0) && chosen_scan_position >= earliest_scan_position && chosen_scan_position <= latest_scan_position && (direction == -1 || chosen_scan_position + next_mip->scan_size > *positions_to_scan.rbegin()))); chosen_scan_position += direction)
	{
		int strand_index; 
		int strand_iterations;
		if(strand_to_use != -1) 
		{
			strand_index = 1 - strand_to_use;
			strand_iterations = 1;
		}
		else 
		{
			strand_index = rand() % 2; // otherwise the MIPs favor the plus strand
			strand_iterations = 2;
		}
		for (int i = 0; i < strand_iterations; i++)
		{
			strand_index = 1 - strand_index;
			if (scan_strand_best_mip[chosen_scan_position].find(strands[strand_index]) == scan_strand_best_mip[chosen_scan_position].end()) continue; // it is possible that condensing eliminated one of the strands
			bool skip_ahead = false;
			string current_strand = strands[strand_index];
			boost::shared_ptr<SVMipv4> test_mip = scan_strand_best_mip[chosen_scan_position][current_strand];
			if(args["-score_method"] == "mixed")
			{
				test_mip->get_parameters(scoring_parameters, feature->long_range_content);
				test_mip->score = predict_value(scoring_parameters, model);
			}
			previous_scan_extent = test_mip->scan_stop_position;
			if (next_mip == 0 || test_mip->score > next_mip->score || test_mip->snp_count < next_mip->snp_count)
			{
				int test_copy_count = test_mip->ext_probe_copy > test_mip->lig_probe_copy ? test_mip->ext_probe_copy : test_mip->lig_probe_copy;
				double test_masked_proportion = test_mip->arm_fraction_masked;
				if (test_copy_count > target_arm_copy && next_mip != 0)
				{
					int next_copy_count = next_mip->ext_probe_copy > next_mip->lig_probe_copy ? next_mip->ext_probe_copy : next_mip->lig_probe_copy;
					if (test_copy_count > next_copy_count) continue;
				}
				if (test_masked_proportion > boost::lexical_cast<double>(args["-masked_arm_threshold"]) && next_mip != 0)
				{
					double next_masked_proportion = next_mip->arm_fraction_masked;
					if (test_masked_proportion > next_masked_proportion) continue;
				}
				for (int position_to_occupy = test_mip->ext_probe_start; position_to_occupy <= test_mip->ext_probe_stop; position_to_occupy++)
				{
					if (chr_strand_pos_used_arm_bases.find(feature->chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[feature->chr].find(current_strand) != chr_strand_pos_used_arm_bases[feature->chr].end() && chr_strand_pos_used_arm_bases[feature->chr][current_strand].count(position_to_occupy) == 1) skip_ahead = true;
				}
				if (skip_ahead) continue;
				for (int position_to_occupy = test_mip->lig_probe_start; position_to_occupy <= test_mip->lig_probe_stop; position_to_occupy++)
				{
					if (chr_strand_pos_used_arm_bases.find(feature->chr) != chr_strand_pos_used_arm_bases.end() && chr_strand_pos_used_arm_bases[feature->chr].find(current_strand) != chr_strand_pos_used_arm_bases[feature->chr].end() && chr_strand_pos_used_arm_bases[feature->chr][current_strand].count(position_to_occupy) == 1) skip_ahead = true;
				}
			if (!(skip_ahead)) next_mip = test_mip;
			}
		}
	}
	if(next_mip == 0) return null_pointer;
	else return next_mip;
}
//processes the new positions covered by a chosen mip and prints out the mip details
void manage_picked_mip (Featurev5 * feature, boost::shared_ptr<SVMipv4> picked_mip, set<int> & positions_to_scan)
{
	picked_mip_counter++;
	PICKEDMIPS << print_details(feature, picked_mip, picked_mip_counter, false);
	if (picked_mip->snp_count == 1 && picked_mip->snp_failed == '0')
	{
		picked_mip->ext_probe_sequence = picked_mip->snp_ext_sequence;
		picked_mip->lig_probe_sequence = picked_mip->snp_lig_sequence;
		picked_mip->mip_seq = picked_mip->snp_mip_sequence;
		SNPMIPS << print_details(feature, picked_mip, picked_mip_counter, true);
	}
	else if(picked_mip->snp_failed == '1')
	{
		SNPMIPS << ">Alternate MIP(s) could not be generated for SNP in arms of MIP #" << picked_mip_counter << endl;
	}
	string other_strand = picked_mip->strand == "+" ? "-" : "+";
	if(args["-seal_both_strands"] == "on")
	{
		for (int occupied_position = picked_mip->ext_probe_start; occupied_position <= picked_mip->ext_probe_stop; occupied_position++)	chr_strand_pos_used_arm_bases[feature->chr][other_strand].insert(occupied_position);
		for (int occupied_position = picked_mip->lig_probe_start; occupied_position <= picked_mip->lig_probe_stop; occupied_position++) chr_strand_pos_used_arm_bases[feature->chr][other_strand].insert(occupied_position);
	}
	else if(args["-half_seal_both_strands"] == "on")
	{
		int occupied_position = (picked_mip->ext_probe_stop + picked_mip->ext_probe_start) / 2; chr_strand_pos_used_arm_bases[feature->chr][other_strand].insert(occupied_position);
		occupied_position = (picked_mip->lig_probe_start + picked_mip->lig_probe_stop) / 2; chr_strand_pos_used_arm_bases[feature->chr][other_strand].insert(occupied_position);
	}
	for (int occupied_position = picked_mip->ext_probe_start; occupied_position <= picked_mip->ext_probe_stop; occupied_position++)	chr_strand_pos_used_arm_bases[feature->chr][picked_mip->strand].insert(occupied_position);
	for (int occupied_position = picked_mip->lig_probe_start; occupied_position <= picked_mip->lig_probe_stop; occupied_position++) chr_strand_pos_used_arm_bases[feature->chr][picked_mip->strand].insert(occupied_position);
	for (int position_scanned = picked_mip->scan_start_position; position_scanned <= picked_mip->scan_stop_position; position_scanned++) positions_to_scan.erase(position_scanned);
}

//LIBSVM method
void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}
//LIBSVM method
double predict_value(vector<double> & parameters, svm_model * model)
{	
//	cout << "predicting" << endl;
	struct svm_node * x = (struct svm_node *) malloc(max_nr_attr * sizeof(struct svm_node));
	int i = 0;
	string line_str("0");
	for(size_t i = 0; i < parameters.size(); i++) //feature count
	{
/*		if((int) i>=max_nr_attr-1)    // need one more for index = -1
		{
			max_nr_attr *= 2;
			x = (struct svm_node *) realloc(x, max_nr_attr * sizeof(struct svm_node));
		}
		x[i].index = (int) i;
		x[i].value = parameters.at(i);*/
//		cout << i + 1<< ":" << parameters.at(i) << " ";
		line_str += " " + boost::lexical_cast<string>(i + 1) + ":" + boost::lexical_cast<string>(parameters.at(i));
	}
//	cout << endl;
//	cout << line_str << endl;
//	x[parameters.size()].index = -1;
	i = 0;
	char * line = new char[line_str.length() + 1];
	line = strcpy(line, line_str.c_str());
	// LIBSVM CODE START 

    int total = 0;

    char *idx, *val, *label, *endptr;
    int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

    label = strtok(line," \t\n");
    if(label == NULL) // empty line
            exit_input_error(total+1);

    strtod(label,&endptr);
    if(endptr == label || *endptr != '\0')
            exit_input_error(total+1);

    while(1)
    {
            if(i>=max_nr_attr-1)    // need one more for index = -1
            {
                    max_nr_attr *= 2;
                    x = (struct svm_node *) realloc(x,max_nr_attr*sizeof(struct svm_node));
            }

            idx = strtok(NULL,":");
            val = strtok(NULL," \t");

            if(val == NULL)
                    break;
            errno = 0;
            x[i].index = (int) strtol(idx,&endptr,10);
            if(endptr == idx || errno != 0 || *endptr != '\0' || x[i].index <= inst_max_index)
                    exit_input_error(total+1);
            else
                    inst_max_index = x[i].index;

            errno = 0;
            x[i].value = strtod(val,&endptr);
            if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
                    exit_input_error(total+1);

            ++i;
    }
    x[i].index = -1;
    //END LIBSVM CODE
	double score =  svm_predict(model, x);
	free(x);
	return score;
}
};
int main(int argc, char * argv[]) {
 
	SVMipv4::set_junction_scores();
	try {
		mipgen * mg = new mipgen(argc, argv);
		mg->print_header();
		mg->query_sequences();
		mg->tile_regions();
	} catch (int e) {
		if(e == 1) { }
		else { cerr << "unable to tile sequences due to circumstance " << e << endl; }
		return 1;
	} catch (exception & e) {
		cerr << "unable to tile sequences" << endl << e.what() << endl;
	}
	
}
