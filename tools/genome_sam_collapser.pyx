# Written by Evan Boyle
# boylee [at] uw.edu
import sys
import re
import numpy as np
cimport numpy as np
cimport cython
from scipy import optimize as optimize
from random import choice
from optparse import OptionParser
from string import maketrans

# container for sam fields and mip information
cdef class sam_read:
  cdef public int read_number, flag, start_coordinate, paired_start_coordinate, tlen, ref_length
  cdef public str sam_line, mip_key, cluster, barcode, mtag, chromosome, cigar, mate_chromosome, seq, qual, read_group
  def __cinit__(self, str sam_line, mip_key_first_coordinate_lookup={}, mip_key_second_coordinate_lookup={}):
    self.sam_line = sam_line.rstrip()
    sam_fields = sam_line.split("\t")
    self.cluster = sam_fields[0].rstrip("/1234")
    barcode_index = sam_fields[0].rfind("#") + 1
    barcode_stop = sam_fields[0].rfind("-")
    if barcode_stop < barcode_index:
      barcode_stop = None
    found_barcode = sam_fields[0][barcode_index:barcode_stop]
    if len(found_barcode) == 0:
      self.barcode = "N"
    else:
      self.barcode = found_barcode
    mtag_index = sam_fields[0].rfind("-") + 1
    if mtag_index < barcode_index:
      self.mtag = ""
    else:
      mtag_length = next((i for i, v in enumerate(sam_fields[0][mtag_index:]) if v not in "ATCGN"), len(sam_fields[0]) - 1)
      self.mtag = sam_fields[0][mtag_index:mtag_index + mtag_length]
    for optional_field in sam_fields[11:]:
      if optional_field[:2] == "RG":
        self.read_group = optional_field[5:]
        break
      self.read_group = None
    self.flag = int(sam_fields[1])
    self.read_number = 1 if self.flag & 128 == 128 else 0
    self.chromosome = sam_fields[2]
    self.start_coordinate = int(sam_fields[3])
    self.cigar = sam_fields[5]
    self.mate_chromosome = sam_fields[6]
    self.paired_start_coordinate = int(sam_fields[7])
    self.tlen = int(sam_fields[8])
    self.ref_length = self.cigar_length()
    self.seq = sam_fields[9]
    self.qual = sam_fields[10]
    self.mip_key = self.lookup_mip(mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup)

  def trim_scan_sequences(self, distance_to_first_scan, distance_to_final_scan, parsed_cigar):
    cdef long d_t_f_s = distance_to_final_scan 
    cdef int read_length = len(self.seq)
    read_scan_start_index = find_first_arm(read_length, distance_to_first_scan, & d_t_f_s, parsed_cigar)
    if read_scan_start_index == -1:
      self.cigar = "0M"
      self.seq = ""
      self.qual = ""
      return 
    read_scan_stop_index = find_second_arm(read_length, & d_t_f_s, read_scan_start_index, parsed_cigar)
    self.adjust_fields(parsed_cigar, read_scan_start_index, read_scan_stop_index, d_t_f_s)

  def add_or_replace_readgroup(self, id):
    if self.read_group != None:
      self.sam_line.replace("RG:Z:" + self.read_group, "RG:Z:" + id)
    else:
      self.sam_line += "\tRG:Z:" + id

  def count_indel(self, wt_mip_counter, mutant_mip_counter):
    if 'N' in self.mtag or 'g' in self.mip_key:
      return
    if(self.cigar != "*"):
      parsed_cigar = []
      trim_and_edit(self)
      if len(self.cigar) <= 1:
        return
      if "I" in self.cigar or "D" in self.cigar:
        mutant_mip_counter[self.mip_key] += 1
      else:
        wt_mip_counter[self.mip_key] += 1
    return

  def verify_cigar(self):
    return self.cigar_seq_length() == len(self.seq)

  def text(self):
    new_fields = self.sam_line.split()
    new_fields[3] = self.start_coordinate
    new_fields[5] = self.cigar
    new_fields[7] = self.paired_start_coordinate
    new_fields[8] = self.tlen
    new_fields[9] = self.seq
    new_fields[10] = self.qual
    return '\t'.join([str(x) for x in new_fields])

  def cigar_seq_length(self):
    if self.cigar == '*':
      return 0
    parsed_cigar = []
    cigar_pattern = re.compile("([\d]+)([A-Z])")
    for pair in cigar_pattern.findall(self.cigar):
      parsed_cigar.append([int(pair[0]), pair[1]])
    seq_length = 0
    for pair in parsed_cigar:
      if pair[1] in "MIS=X":
        seq_length += pair[0]
    return seq_length

  def cigar_length(self):
    if self.cigar == '*':
      return 0
    parsed_cigar = []
    cigar_pattern = re.compile("([\d]+)([A-Z])")
    for pair in cigar_pattern.findall(self.cigar):
      parsed_cigar.append([int(pair[0]), pair[1]])
    read_length = 0
    for pair in parsed_cigar:
      if pair[1] in "MD":
        read_length += pair[0]
    return read_length

  def adjust_fields(self, parsed_cigar, int read_scan_start_index, int read_scan_stop_index, int distance_to_final_scan):
    adjusted_cigar = ""
    self.ref_length = 0
    if(distance_to_final_scan > 0):
      self.seq = self.seq[read_scan_start_index:]
      self.qual = self.qual[read_scan_start_index:]
    else:
      self.seq = self.seq[read_scan_start_index:read_scan_stop_index + 1]
      self.qual = self.qual[read_scan_start_index:read_scan_stop_index + 1]
    # print sam_fields[9]
    for pair in parsed_cigar:
      adjusted_cigar += str(pair[0]) + pair[1]
      self.ref_length += pair[0]
    self.cigar = adjusted_cigar

  cdef str lookup_mip(self, mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup):
    if self.flag & 1 == 1:
      if(self.flag & 163 == 163 and self.chromosome + ':' + str(self.start_coordinate) + ":+" in mip_key_first_coordinate_lookup): # plus paired
        coordinates = self.chromosome + ':' + str(self.start_coordinate) + ":+"
        current_mip_key = mip_key_first_coordinate_lookup[coordinates]
      elif self.flag & 83 == 83 and self.chromosome + ':' + str(self.start_coordinate + self.ref_length - 1) + ":+" in mip_key_second_coordinate_lookup: 
        coordinates = self.chromosome + ':' + str(self.start_coordinate + self.ref_length - 1) + ":+"
        current_mip_key = mip_key_second_coordinate_lookup[coordinates]
      elif self.flag & 99 == 99 and self.chromosome + ':' + str(self.start_coordinate) + ":-" in mip_key_first_coordinate_lookup: # minus paired
        coordinates = self.chromosome + ':' + str(self.start_coordinate) + ":-"
        current_mip_key = mip_key_first_coordinate_lookup[coordinates]
      elif self.flag & 147 == 147 and self.chromosome + ':' + str(self.start_coordinate + self.ref_length - 1) + ":-" in mip_key_second_coordinate_lookup:
        coordinates = self.chromosome + ':' + str(self.start_coordinate + self.ref_length - 1) + ":-"
        current_mip_key = mip_key_second_coordinate_lookup[coordinates]
        #print current_mip_key, "-2"
      else:
        min_position = self.start_coordinate if self.start_coordinate < self.paired_start_coordinate else self.paired_start_coordinate
        max_position = self.tlen if self.tlen > 0 else -1 * self.tlen
        current_mip_key = self.chromosome + ':' + str(min_position) + '-' + str(max_position) + "/0,0/g"
    elif self.flag == 16 and self.chromosome + ':' + str(self.start_coordinate) + ":+" in mip_key_first_coordinate_lookup: # plus single
      coordinates = self.chromosome + ':' + str(self.start_coordinate) + ":+" # SWITCH 0 and 16 FOR PROTON
      current_mip_key = mip_key_first_coordinate_lookup[coordinates]
    elif self.flag == 0 and self.chromosome + ':' + str(self.start_coordinate) + ":-" in mip_key_first_coordinate_lookup: # minus single
      coordinates = self.chromosome + ':' + str(self.start_coordinate) + ":-"
      current_mip_key = mip_key_first_coordinate_lookup[coordinates]
    else:
      current_mip_key = self.chromosome + ':' + str(self.start_coordinate) + '-' + str(self.start_coordinate + self.ref_length) + "/0,0/g"
    return current_mip_key

  cdef int enter_read(self, consensus_details, observed_sites):
    start_and_cigar = (self.start_coordinate, self.cigar)
    if self.mip_key not in observed_sites.keys():
      observed_sites[self.mip_key] = {}
      consensus_details[self.mip_key] = {}
    if self.barcode not in observed_sites[self.mip_key].keys():
      observed_sites[self.mip_key][self.barcode] = {}
      consensus_details[self.mip_key][self.barcode] = {}
    if self.mtag not in observed_sites[self.mip_key][self.barcode].keys():
      observed_sites[self.mip_key][self.barcode][self.mtag] = [{}, {}]
      #print observed_sites[self.mip_key][self.barcode].keys()
      consensus_details[self.mip_key][self.barcode][self.mtag] = [None, None]
    if consensus_details[self.mip_key][self.barcode][self.mtag][self.read_number] == None:
      consensus_details[self.mip_key][self.barcode][self.mtag][self.read_number] = self
    if start_and_cigar not in observed_sites[self.mip_key][self.barcode][self.mtag][self.read_number].keys():
      observed_sites[self.mip_key][self.barcode][self.mtag][self.read_number][start_and_cigar] = [[self.cluster, self.seq, self.qual]]
    else:
      observed_sites[self.mip_key][self.barcode][self.mtag][self.read_number][start_and_cigar].append([self.cluster, self.seq, self.qual])
    return 0
    
cdef int read_barcodes(barcode_labels, barcode_lookup, complexity_by_sample, options):  
  cdef str barcode_label, barcode_seq, rev_seq, for_seq, bases 
  barcode_fh = open(options.barcode_file, 'r')
  for line in barcode_fh:
    if(options.dual_indexed):
      barcode_label, rev_seq, for_seq = line.rstrip().split()[0:3]
      barcode_seq = rev_seq + for_seq
    else:
      barcode_label, barcode_seq = line.rstrip().split()[0:2]
    barcode_labels[barcode_seq] = barcode_label
    barcode_lookup[barcode_seq] = barcode_seq
    complexity_by_sample[barcode_seq] = [0, 0]
    if(options.allow_ambiguous_barcodes):
      bases = "ATCGN"
      for i in range(len(barcode_seq)):
        native_seq = list(barcode_seq)
        for j in range(5):
          mutated_seq = native_seq
          mutated_seq[i] = bases[j]
          barcode_lookup["".join(mutated_seq)] = barcode_seq
  barcode_fh.close()
  return 0

cdef int read_mips(arm_sequences, uniformity_keys, mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup, complexity_by_position, mip_summary, options):
  complement = maketrans("GCATgcat", "CGTAcgta")
  mip_fh = open(options.mip_file, 'r')
  indices = {}
  indices["chr"] = 2
  indices["ext_probe_start"] = 3
  indices["ext_probe_stop"] = 4
  indices["ext_probe_sequence"] = 6
  indices["lig_probe_start"] = 7
  indices["lig_probe_stop"] = 8
  indices["lig_probe_sequence"] = 10
  indices["probe_strand"] = 17
  cdef str ext_arm, lig_arm
  for line in mip_fh:
    if(line.startswith('#')):
      continue
    if(line.startswith('>')):
      header_fields = line.rstrip().split()
      for i, column in enumerate(header_fields):
        if column in indices.keys():
          indices[column] = i
      continue
    mip_fields = line.rstrip().split()
    ext_arm = mip_fields[indices["ext_probe_sequence"]]
    lig_arm = mip_fields[indices["lig_probe_sequence"]]
    
    if mip_fields[indices["probe_strand"]] == '+':
#chr:capture_start-capture_end
      mip_key = mip_fields[indices["chr"]] + ':' + mip_fields[indices["ext_probe_start"]] + '-' + mip_fields[indices["lig_probe_stop"]] + '/' + str(len(ext_arm)) + ',' + str(len(lig_arm)) + "/+"
      arm_sequences[mip_key] = (ext_arm, lig_arm)
      for start in range(int(mip_fields[indices["ext_probe_start"]]) - options.flex_space, int(mip_fields[indices["ext_probe_start"]]) + options.flex_space + 1):
        mip_key_first_coordinate_lookup[mip_fields[indices["chr"]] + ':' + str(start) + ":+"] = mip_key
      for stop in range(int(mip_fields[indices["lig_probe_stop"]]) - options.flex_space, int(mip_fields[indices["lig_probe_stop"]]) + options.flex_space + 1):
        mip_key_second_coordinate_lookup[mip_fields[indices["chr"]] + ':' + str(stop) + ":+"] = mip_key
    else:
      mip_key = mip_fields[indices["chr"]] + ':' + mip_fields[indices["lig_probe_start"]] + '-' + mip_fields[indices["ext_probe_stop"]] + '/' + str(len(ext_arm)) + ',' + str(len(lig_arm)) + "/-"
      #print mip_key
      arm_sequences[mip_key] = (ext_arm[::-1].translate(complement), lig_arm[::-1].translate(complement))      
      for start in range(int(mip_fields[indices["lig_probe_start"]]) - options.flex_space, int(mip_fields[indices["lig_probe_start"]]) + options.flex_space + 1):
        mip_key_first_coordinate_lookup[mip_fields[indices["chr"]] + ':' + str(start) + ":-"] = mip_key
      for stop in range(int(mip_fields[indices["ext_probe_stop"]]) - options.flex_space, int(mip_fields[indices["ext_probe_stop"]]) + options.flex_space + 1):
        mip_key_second_coordinate_lookup[mip_fields[indices["chr"]] + ':' + str(stop) + ":-"] = mip_key
    uniformity_keys.add(mip_key)
    complexity_by_position[mip_key] = [0, 0]
    mip_summary[mip_key] = line.rstrip()
  mip_fh.close()
  #for k,v in mip_key_first_coordinate_lookup.iteritems():
    #print k,v
  return 0

def find_dominant_start_and_cigar(dictionary):
  dominant_start_and_cigar = max(dictionary)
  cdef int current_size = len(dictionary[dominant_start_and_cigar])
  cdef int current_qual_sum = sum([sum([ord(qual) - 33 for qual in clus_bases_quals[2]]) for clus_bases_quals in dictionary[dominant_start_and_cigar]])
  cdef int alt_qual_sum
  for start_and_cigar in dictionary:
    if(len(dictionary[start_and_cigar]) > current_size):
      dominant_start_and_cigar = start_and_cigar
      current_size = len(dictionary[start_and_cigar])
    elif(len(dictionary[start_and_cigar]) == current_size):
      alt_qual_sum = sum([sum([ord(qual) - 33 for qual in clus_bases_quals[2]]) for clus_bases_quals in dictionary[start_and_cigar]])
      if alt_qual_sum > current_qual_sum:
        dominant_start_and_cigar = start_and_cigar
        current_qual_sum = alt_qual_sum
  return dominant_start_and_cigar

cdef find_consensus(read_list): # alternate implementation possible: do not track best base
  cdef long read_length, quality_sum, position
  read_length = min([len(read[1]) for read in read_list])
  consensus = [{} for i in range(read_length)]
  for position in range(read_length):
    for base in "ATGCN":
      consensus[position][base] = 0
  best_bases = ["N" for i in range(read_length)]
  for read in read_list:
    read[1] = read[1][0:read_length]
    for position in range(read_length):
      if(read[1][position] == "N"):
        continue
      consensus[position][read[1][position]] += ord(read[2][position]) - 33
      if(consensus[position][read[1][position]] > consensus[position][best_bases[position]]):
        best_bases[position] = read[1][position]
      elif(consensus[position][read[1][position]] == consensus[position][best_bases[position]] and read[1][position] != best_bases[position]):
        best_bases[position] = "N"
  final_read = [[],[]]
  for position in range(read_length):
    final_read[0].append(best_bases[position])
    quality_sum = 0
    for base in "ATGC":
      if(base == best_bases[position]):
        quality_sum += consensus[position][base]
      else:
        quality_sum -= consensus[position][base]
    if(quality_sum > 40):
      final_read[1].append('I')
    elif(quality_sum < 2):
      final_read[1].append('#')
    else:
      final_read[1].append(chr(quality_sum + 33))
  final_read[0] = "".join(final_read[0])
  final_read[1] = "".join(final_read[1])
  #print final_read, "is final read"
  return final_read

cdef int manage_sites(observed_sites, consensus_details, recent_mip_key, uniformity_keys, barcode_labels, complexity_by_position, complexity_by_sample, mtag_population, probability_list, long * all_total_reads, long * all_unique_reads, file_handles, options):
  #print "managing"
  cdef int barcode_site_total_reads, barcode_site_unique_reads, num_collected
  for active_mip_key in observed_sites.keys():
    active_chr, active_range = active_mip_key.split('/')[0].split(':')
    active_start, active_stop = [int(x) for x in active_range.split('-')]
    recent_chr, recent_range = recent_mip_key.split('/')[0].split(':')
    recent_start, recent_stop = [int(x) for x in recent_range.split('-')]
    if active_chr == recent_chr and active_stop > recent_start - options.flex_space: # not done with site yet
      continue
    for observed_barcode in observed_sites[active_mip_key].keys():
      barcode_site_clusters = set()
      barcode_site_unique_reads = len(observed_sites[active_mip_key][observed_barcode].keys())
      if(barcode_site_unique_reads > len(probability_list)):
        for untaken in range(int(mtag_population) - len(probability_list), int(mtag_population) - barcode_site_unique_reads, -1):
          probability_list.append(probability_list[-1] * untaken / mtag_population) # builds probabilities so that they don't have to be calculated multiple times
      duplicate_free_probability = probability_list[barcode_site_unique_reads - 1]
      for mtag in observed_sites[active_mip_key][observed_barcode].keys():
        for read_index in range(2):
          for s_and_c, read_list in observed_sites[active_mip_key][observed_barcode][mtag][read_index].iteritems():
            for entry in read_list:
              barcode_site_clusters.add(entry[0])
        if len(observed_sites[active_mip_key][observed_barcode][mtag][0]) == 0: # pair not present
          file_handles["unpaired_output"].write(consensus_details[active_mip_key][observed_barcode][mtag][1].text() + "\n")
          continue
        elif not options.single_end and len(observed_sites[active_mip_key][observed_barcode][mtag][1]) == 0:
          file_handles["unpaired_output"].write(consensus_details[active_mip_key][observed_barcode][mtag][0].text() + "\n")
          continue
        if not options.single_end and not consensus_details[active_mip_key][observed_barcode][mtag][1].verify_cigar():
          consensus_details[active_mip_key][observed_barcode][mtag][1].cigar = '*'
          for read_index in range(2):
            file_handles["strange_alignments"].write(consensus_details[active_mip_key][observed_barcode][mtag][read_index].text() + "\n")
          continue
        elif not consensus_details[active_mip_key][observed_barcode][mtag][0].verify_cigar():
          consensus_details[active_mip_key][observed_barcode][mtag][0].cigar = '*'
          file_handles["strange_alignments"].write(consensus_details[active_mip_key][observed_barcode][mtag][0].text() + "\n")
          if not options.single_end:
            file_handles["strange_alignments"].write(consensus_details[active_mip_key][observed_barcode][mtag][1].text() + "\n")
          continue
        consensus_number = 1 if options.single_end else 2
        for read_index in range(consensus_number):
          consensus_read = consensus_details[active_mip_key][observed_barcode][mtag][read_index]
          if(duplicate_free_probability > options.confidence_level): # tags are not saturating
            dominant_start, dominant_cigar = find_dominant_start_and_cigar(observed_sites[active_mip_key][observed_barcode][mtag][read_index])
            consensus_read.start_coordinate = dominant_start
            consensus_read.cigar = dominant_cigar
            consensus_seq = find_consensus(observed_sites[active_mip_key][observed_barcode][mtag][read_index][(dominant_start, dominant_cigar)])
            consensus_read.seq = consensus_seq[0]
            consensus_read.qual = consensus_seq[1]
          else: # too many tags: pick at random to avoid chimeric consensuses
            start_cigar_and_details = []
            for (s, c), d in observed_sites[active_mip_key][observed_barcode][mtag][read_index].iteritems():
              for seq in d:
                start_cigar_and_details.append((s, c, seq))
            random_start, random_cigar, random_details = choice(start_cigar_and_details)
            consensus_read.start_coordinate = random_start
            consensus_read.cigar = random_cigar
            consensus_read.seq = random_details[1]
            consensus_read.qual = random_details[2]

        consensus_read_one = consensus_details[active_mip_key][observed_barcode][mtag][0]
        if not options.single_end:
          consensus_read_two = consensus_details[active_mip_key][observed_barcode][mtag][1]
        if not options.single_end and not (consensus_read_one.verify_cigar() and consensus_read_two.verify_cigar()):
          for read in [consensus_read_one, consensus_read_two]:
            if not read.verify_cigar():
              read.cigar = '*'
          for read in [consensus_read_one, consensus_read_two]:
            file_handles["strange_alignments"].write(read.text() + "\n")
          continue
        elif options.single_end and not consensus_read_one.verify_cigar():
          consensus_read_one.cigar = '*'
          file_handles["strange_alignments"].write(consensus_read_one.text() + "\n")
        
        if not options.single_end:
          first_position = min([consensus_read_two.start_coordinate, consensus_read_one.start_coordinate])
          last_position = max([consensus_read_two.start_coordinate + consensus_read_two.ref_length, consensus_read_one.start_coordinate + consensus_read_one.ref_length])
          consensus_read_two.paired_start_coordinate = consensus_read_one.start_coordinate
          consensus_read_one.paired_start_coordinate = consensus_read_two.start_coordinate
          if consensus_read_two.flag == 99 or consensus_read_one.flag == 163:
            consensus_read_two.tlen = last_position - first_position
            consensus_read_one.tlen = first_position - last_position
          elif consensus_read_one.flag == 99 or consensus_read_two.flag == 163:
            consensus_read_one.tlen = last_position - first_position
            consensus_read_two.tlen = first_position - last_position
        if(options.barcode_file != None and not options.merge_samples):
          file_handles[observed_barcode].write(consensus_read_one.text() + "\n")
          if not options.single_end:
            file_handles[observed_barcode].write(consensus_read_two.text() + "\n")
        else:
          file_handles["merged_output"].write(consensus_read_one.text() + "\n")
          if not options.single_end:
            file_handles["merged_output"].write(consensus_read_two.text() + "\n")
      barcode_site_total_reads = len(barcode_site_clusters)
      if(active_mip_key in uniformity_keys):
        if(options.barcode_file != None):
          complexity_by_sample[observed_barcode][0] += barcode_site_total_reads
          complexity_by_sample[observed_barcode][1] += barcode_site_unique_reads
        if(options.mip_file != None):
          complexity_by_position[active_mip_key][0] += barcode_site_total_reads
          complexity_by_position[active_mip_key][1] += barcode_site_unique_reads
        all_total_reads[0] += barcode_site_total_reads
        all_unique_reads[0] += barcode_site_unique_reads
        if options.barcode_file == None:
          label = observed_barcode
        else:
          label = barcode_labels[observed_barcode]
        file_handles["complexity_summary"].write(label + '\t' + active_mip_key + '\t' + str(barcode_site_total_reads) + '\t' + str(barcode_site_unique_reads) + '\n')
    del observed_sites[active_mip_key]
    del consensus_details[active_mip_key]
  return 0

cdef int find_first_arm(int read_length, int distance_to_first_scan, long * distance_to_final_scan, parsed_cigar):
  cdef int read_scan_start_index = 0
  while(distance_to_first_scan > 0 and read_scan_start_index < read_length):
    if(len(parsed_cigar) == 0):
      distance_to_final_scan[0] = 0
      return -1
    if(parsed_cigar[0][1] == 'D'):
      if(parsed_cigar[0][0] <= distance_to_first_scan):
        distance_to_first_scan -= parsed_cigar[0][0]
        distance_to_final_scan[0] -= parsed_cigar[0][0]
        parsed_cigar.pop(0)
      else:
        parsed_cigar[0][0] -= distance_to_first_scan
        distance_to_final_scan[0] -= distance_to_first_scan
        distance_to_first_scan = 0
    elif(parsed_cigar[0][1] == 'I' or parsed_cigar[0][1] == 'S'):
      read_scan_start_index += parsed_cigar[0][0]
      parsed_cigar.pop(0)
    elif(parsed_cigar[0][1] in "M=X"):
      if(parsed_cigar[0][0] <= distance_to_first_scan):
        distance_to_first_scan -= parsed_cigar[0][0]
        distance_to_final_scan[0] -= parsed_cigar[0][0]
        read_scan_start_index += parsed_cigar[0][0]
        parsed_cigar.pop(0)
      else:
        parsed_cigar[0][0] = parsed_cigar[0][0] - distance_to_first_scan
        read_scan_start_index += distance_to_first_scan

        distance_to_final_scan[0] -= distance_to_first_scan
        distance_to_first_scan = 0
    elif(parsed_cigar[0][1] == 'H'):
      parsed_cigar.pop(0)
    else:
       sys.stderr.write("MORE FUNCTIONALITY NEEDED IN CIGAR PROCESSING, SORRY\n")
       sys.exit()
  #print "returning", read_scan_start_index
  return read_scan_start_index
cdef int find_second_arm(read_length, long * distance_to_final_scan, read_scan_start_index, parsed_cigar):
  cdef int cigar_index, read_scan_stop_index
  cigar_index = 0
  read_scan_stop_index = read_scan_start_index
  while(distance_to_final_scan[0] >= 0 and read_scan_stop_index < read_length):
    if(cigar_index == len(parsed_cigar)): # arm not in read or arm found
      break
    if(parsed_cigar[cigar_index][1] == 'D'):
      if(parsed_cigar[cigar_index][0] <= distance_to_final_scan[0]):
        distance_to_final_scan[0] -= parsed_cigar[cigar_index][0]
      else:
        parsed_cigar[cigar_index][0] = distance_to_final_scan[0] + 1
        read_scan_stop_index -= 1 # must go backward 
        distance_to_final_scan[0] = -1
    elif(parsed_cigar[cigar_index][1] == 'I' or parsed_cigar[cigar_index][1] == 'S'):
      read_scan_stop_index += parsed_cigar[cigar_index][0]
    elif(parsed_cigar[cigar_index][1] in "M=X"):
      if(parsed_cigar[cigar_index][0] <= distance_to_final_scan[0]):
        distance_to_final_scan[0] -= parsed_cigar[cigar_index][0]
        read_scan_stop_index += parsed_cigar[cigar_index][0]
      else:
        parsed_cigar[cigar_index][0] = distance_to_final_scan[0] + 1
        read_scan_stop_index += distance_to_final_scan[0]
        distance_to_final_scan[0] = -1
    elif(parsed_cigar[cigar_index][1] == 'H'):
      pass
    else:
      sys.stderr.write("MORE FUNCTIONALITY NEEDED IN CIGAR PROCESSING, SORRY\n")
      sys.exit()
    cigar_index += 1
  excess_cigar_count = len(parsed_cigar) - cigar_index
  for i in range(excess_cigar_count):
    parsed_cigar.pop()
  #print "final cigar", parsed_cigar
  return read_scan_stop_index


cdef int trim_and_edit(current_read, arm_sequences={}):
  cigar_pattern = re.compile("([\d]+)([A-Z])")
  cdef long read_scan_start_index, distance_to_first_scan, distance_to_final_scan, read_scan_stop_index
  mip_details = current_read.mip_key.split('/')
  mip_start, mip_stop = [int(x) for x in mip_details[0].split(':')[1].split('-')]
  ext_length, lig_length = [int(x) for x in mip_details[1].split(',')]
  read_length = len(current_read.seq)
  current_read_start = int(current_read.start_coordinate)
  if current_read.flag & 1 == 0: 
    if current_read.flag & 16 == 0: # minus mip
      distance_to_first_scan = int(mip_details[1].split(',')[1]) - current_read_start + mip_start
      distance_to_final_scan = mip_stop - current_read_start - int(mip_details[1].split(',')[0])
      if current_read.mip_key in arm_sequences.keys():
        second_arm_reference, first_arm_reference = arm_sequences[current_read.mip_key]
    else: # plus mip
      distance_to_first_scan = int(mip_details[1].split(',')[0]) - current_read_start + mip_start
      distance_to_final_scan = mip_stop - current_read_start - int(mip_details[1].split(',')[1])
      if current_read.mip_key in arm_sequences.keys():
        first_arm_reference, second_arm_reference = arm_sequences[current_read.mip_key]
  elif current_read.flag == 83 or current_read.flag == 163: #plus mips, read 1 and 2 
    distance_to_first_scan = -current_read_start + int(mip_details[1].split(',')[0]) + mip_start
    distance_to_final_scan = mip_stop - current_read_start - int(mip_details[1].split(',')[1])
    if current_read.mip_key in arm_sequences.keys():
      first_arm_reference, second_arm_reference = arm_sequences[current_read.mip_key]
  elif current_read.flag == 99 or current_read.flag == 147: #minus mips, read 1 and 2
    distance_to_first_scan = -current_read_start + int(mip_details[1].split(',')[1]) + mip_start
    distance_to_final_scan = mip_stop - current_read_start - int(mip_details[1].split(',')[0])
    if current_read.mip_key in arm_sequences.keys():
      second_arm_reference, first_arm_reference = arm_sequences[current_read.mip_key]
  else:
    distance_to_final_scan = 0
    return 1
  #print "distance", distance_to_first_scan
  if distance_to_first_scan > 0:
    current_read.start_coordinate = current_read.start_coordinate + distance_to_first_scan
  if(current_read.cigar != "*"):
    parsed_cigar = []
    for pair in cigar_pattern.findall(current_read.cigar):
      parsed_cigar.append([int(pair[0]), pair[1]])
    adjusted_cigar = ""
    read_scan_start_index = find_first_arm(read_length, distance_to_first_scan, & distance_to_final_scan, parsed_cigar)
    if read_scan_start_index == -1:
      current_read.cigar = "0M"
      current_read.seq = ""
      current_read.qual = ""
      return 0
    elif len(arm_sequences) > 0:
      if first_arm_reference[:read_scan_start_index] != current_read.seq[:read_scan_start_index]:
        return 2 
    read_scan_stop_index = find_second_arm(read_length, & distance_to_final_scan, read_scan_start_index, parsed_cigar)    
    if len(arm_sequences) > 0:
      partial_arm_length = min(len(second_arm_reference), len(current_read.seq) - read_scan_stop_index)
      if second_arm_reference[:partial_arm_length] != current_read.seq[read_scan_stop_index + 1: read_scan_stop_index + 1 + partial_arm_length]:
        return 2
    current_read.adjust_fields(parsed_cigar, read_scan_start_index, read_scan_stop_index, distance_to_final_scan)
  return 0

def find_saturation(size_of_mip_pool, mips_sampled, unique_mips_observed):
  return size_of_mip_pool * (1. - pow(1 - 1. / size_of_mip_pool, float(mips_sampled))) - float(unique_mips_observed)

cpdef int output_mip_complexity(complexity_by_position, fh):
  for mip_key in complexity_by_position:
    position_total_reads = complexity_by_position[mip_key][0]
    position_unique_reads = complexity_by_position[mip_key][1]
    if(position_unique_reads == 0):
      position_saturation = "NA"
    elif(position_total_reads > position_unique_reads):
      estimated_sample_size = optimize.brentq(find_saturation, position_unique_reads, position_total_reads * position_unique_reads / float(position_total_reads - position_unique_reads + .1), args = (position_total_reads, position_unique_reads))
      estimated_saturation = int(position_unique_reads / estimated_sample_size * 100)
      if(estimated_saturation < 8):
        position_saturation = "<8%"
      else:
        position_saturation = str(estimated_saturation) + "%"
    else:
      position_saturation = "~0%"
    fh.write('ALL_SAMPLES\t' + mip_key + '\t' + str(complexity_by_position[mip_key][0]) + '\t' + str(complexity_by_position[mip_key][1]) + '\t' + position_saturation + '\n')
  return 0

cpdef int output_sample_complexity(complexity_by_sample, barcode_labels, fh):
  for barcode in complexity_by_sample:
    sample_total_reads = complexity_by_sample[barcode][0]
    sample_unique_reads = complexity_by_sample[barcode][1]
    if(sample_unique_reads == 0):
      sample_saturation = "NA"
    elif(sample_total_reads > sample_unique_reads):
      estimated_sample_size = optimize.brentq(find_saturation, sample_unique_reads, sample_total_reads * sample_unique_reads / float(sample_total_reads - sample_unique_reads + .1), args = (sample_total_reads, sample_unique_reads))
      estimated_saturation = int(sample_unique_reads / estimated_sample_size * 100)
      if(estimated_saturation < 8):
        sample_saturation = "<8%"
      else:
        sample_saturation = str(estimated_saturation) + "%"
    else:
      sample_saturation = "~0%"
    fh.write(barcode_labels[barcode] + '\t' + 'ALL_POSITIONS\t' + str(sample_total_reads) + '\t' + str(sample_unique_reads) + '\t' + sample_saturation + '\n')
  return 0

cdef int output_overall_complexity(long * all_unique_reads, long * all_total_reads, fh):
  if(all_unique_reads[0] == 0):
    all_saturation = "NA"
  elif(all_total_reads[0] == all_unique_reads[0]):
    all_saturation = "~0%"
  elif (find_saturation(all_unique_reads[0], all_total_reads[0], all_unique_reads[0]) < 0 and find_saturation(all_unique_reads[0] * all_total_reads[0] / (all_total_reads[0] - all_unique_reads[0]), all_total_reads[0], all_unique_reads[0]) < 0):
    all_saturation = ">99%"
  elif(all_total_reads[0] > all_unique_reads[0]):
    estimated_size = optimize.brentq(find_saturation, all_unique_reads[0], all_unique_reads[0] * all_total_reads[0] / (all_total_reads[0] - all_unique_reads[0]), args = (all_total_reads[0], all_unique_reads[0]))
    if(estimated_size < 1 or (float(all_unique_reads[0]) / estimated_size) < 0.08):
      all_saturation = "<8%"
    else:
      all_saturation = str(int(float(all_unique_reads[0]) / estimated_size * 100)) + "%"
  else:
    all_saturation = "~0%"
  fh.write('ALL_SAMPLES\tALL_POSITIONS\t' + str(all_total_reads[0]) + '\t' + str(all_unique_reads[0]) + '\t' + all_saturation + '\n')
  return 0

def initialize_and_iterate(options):
  complexity_by_sample = {}
  complexity_by_position = {}
  targeted_sites = {}
  barcode_lookup = {}
  barcode_labels = {}
  mip_summary = {}
  if(options.barcode_file != None):
    read_barcodes(barcode_labels, barcode_lookup, complexity_by_sample, options)
  else:
    barcode_labels['N'] = "anonymous"
  arm_sequences = {}
  uniformity_keys = set()
  mip_key_first_coordinate_lookup = {}
  mip_key_second_coordinate_lookup = {}

  if(options.mip_file != None):
    read_mips(arm_sequences, uniformity_keys, mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup, complexity_by_position, mip_summary, options)
  file_handles = {}
  observed_sites = {}
  consensus_details = {}
  all_total_reads = 0
  all_unique_reads = 0

  mtag_length = int(sys.argv[1])
  probability_list = [1.]
  mtag_population = float(pow(4, mtag_length))
  file_handles["notes"] = open(sys.argv[2] + ".notes.txt", 'w')
  file_handles["unpaired_output"] = open(sys.argv[2] + ".unpaired_reads.sam", 'w')
  file_handles["improper_pairs"] = open(sys.argv[2] + ".improper_pairs.sam", 'w')
  file_handles["strange_alignments"] = open(sys.argv[2] + ".strange_alignments.sam", 'w')
  file_handles["off_target_output"] = open(sys.argv[2] + ".off_target_reads.sam", 'w')
  if options.filter_softclips:
    file_handles["softclipped_output"] = open(sys.argv[2] + ".softclipped.sam", 'w')
  if(options.barcode_file == None or options.merge_samples):
    if options.collapse_free:
      file_handles["merged_output"] =  open(sys.argv[2] + ".all_reads.uncollapsed.sam", 'w')
    else:
      file_handles["merged_output"] =  open(sys.argv[2] + ".all_reads.unique.sam", 'w')
  else:
    for barcode,label in barcode_labels.iteritems():
      if options.collapse_free:
        file_handles[barcode] = open(sys.argv[2] + "." + label + ".uncollapsed.sam",'w')
      else:
        file_handles[barcode] = open(sys.argv[2] + "." + label + ".unique.sam", 'w')
  if options.exact_arms:
    file_handles["imperfect_arms"] = open(sys.argv[2] + ".imperfect_arms.sam", 'w')
  file_handles["complexity_summary"] = open(sys.argv[2] + ".complexity.txt", 'w')
  file_handles["complexity_summary"].write("sample\tmip\ttotal\tunique\n")
  file_handles["mipwise_summary"] = open(sys.argv[2] + ".mipwise_summary.txt", 'w')
  file_handles["mipwise_summary"].write("sample\tmip\ttotal\tunique\tsaturation\n")
  file_handles["samplewise_summary"] = open(sys.argv[2] + ".samplewise_summary.txt", 'w')
  file_handles["samplewise_summary"].write("sample\tmip\ttotal\tunique\tsaturation\n")
  good_flags = frozenset([0,16,147,99,83,163])
  reads_skipped = 0
  off_target_reads = 0
  improper_pairs = 0
  softclippings = 0
  reads_unmapped = 0
  readgroups_printed = not options.add_or_replace_readgroups
  for sam_line in sys.stdin:
    if not readgroups_printed:
      if not sam_line.startswith("@") or sam_line.startswith("@PG") or sam_line.startswith("@CO"):
        for pooled_output in ["imperfect_arms", "improper_pairs", "strange_alignments", "off_target_output", "merged_output", "softclipped_output"]:
          if pooled_output in file_handles.keys():
            for barcode in barcode_labels.keys():
              file_handles[pooled_output].write("@RG\tID:" + barcode_labels[barcode] + "_" + barcode + "\tPL:illumina\tLB:MIPs\tSM:" + barcode_labels[barcode] + "\n")
        if not options.merge_samples:
          for barcode in barcode_labels.keys():
            file_handles[barcode].write("@RG\tID:" + barcode_labels[barcode] + "_" + barcode + "\tPL:illumina\tLB:MIPs\tSM:" + barcode_labels[barcode] + "\n")
        readgroups_printed = True
    if sam_line.startswith("@"):
      for sam_name in ["imperfect_arms", "improper_pairs", "strange_alignments", "off_target_output", "merged_output", "softclipped_output"]:
        if sam_name in file_handles.keys():
          file_handles[sam_name].write(sam_line)
      if options.barcode_file != None:
        for sam_name in barcode_labels.keys():
          if sam_name in file_handles.keys():
            file_handles[sam_name].write(sam_line)
      continue
    current_read = sam_read(sam_line, mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup)
    #print "read parsed"
    if options.barcode_file != None:
      if current_read.barcode not in barcode_lookup:
        reads_skipped += 1
        continue
      else:
        current_read.barcode = barcode_lookup[current_read.barcode]
    if(options.filter_molecular_tags and re.search("[^ATCG]", current_read.mtag)):
      reads_skipped += 1
      continue
    if options.add_or_replace_readgroups:
      current_read.add_or_replace_readgroup(barcode_labels[current_read.barcode] + "_" + current_read.barcode)
    if (not options.single_end and (current_read.mate_chromosome != '=' or (current_read.paired_start_coordinate != 0 and abs(current_read.paired_start_coordinate - current_read.start_coordinate) > 1000))):
      if "discordant_arms" not in file_handles:
        file_handles["discordant_arms"] = open(sys.argv[2] + ".discordant_arms.sam",'w')  
      file_handles["discordant_arms"].write(sam_line)
      continue
    if current_read.flag not in good_flags:
      file_handles["improper_pairs"].write(sam_line)
      improper_pairs += 1
      continue
    if 'S' in current_read.cigar and options.filter_softclips:
      if options.single_end or (re.search("^\d+S", current_read.cigar) and (current_read.flag == 163 or current_read.flag == 99)) or (re.search("S$", current_read.cigar) and (current_read.flag == 83 or current_read.flag == 147)):
        softclippings += 1
        file_handles["softclipped_output"].write(sam_line)
        continue
    if '*' == current_read.chromosome:
      reads_unmapped += 1
      continue
    if options.mip_reference:
      current_read.mip_key = current_read.chromosome
    if current_read.mip_key[-2:] == "/g" and options.mip_file != None:
      file_handles["off_target_output"].write(sam_line)
      off_target_reads += 1
      continue
    mip_start, mip_stop = [int(pos) for pos in current_read.mip_key.split("/")[0].split(':')[1].split('-')]
    if not options.single_end and abs(current_read.tlen) != current_read.ref_length and abs(abs(current_read.tlen) - (mip_stop - mip_start)) > options.flex_space:
      file_handles["off_target_output"].write(sam_line)
      off_target_reads += 1
      continue
    if(options.mip_file and not options.no_trimming):
      if options.exact_arms:
        status = trim_and_edit(current_read, arm_sequences)
      else:
        status = trim_and_edit(current_read)
      if status == 1:
        file_handles["improper_pairs"].write(sam_line)
        improper_pairs += 1
        continue
      if status == 2:
        file_handles["imperfect_arms"].write(sam_line)
        continue
    if options.collapse_free:
      if options.barcode_file == None or options.merge_samples:
        file_handles["merged_output"].write(current_read.text() + "\n")
      else:
        file_handles[current_read.barcode].write(current_read.text() + "\n")
      if options.barcode_file != None:
        complexity_by_sample[current_read.barcode][0] += 1
      if options.mip_file != None and current_read.mip_key in uniformity_keys:
        complexity_by_position[current_read.mip_key][0] += 1
      continue
    if current_read.mip_key not in observed_sites.keys():
      manage_sites(observed_sites, consensus_details, current_read.mip_key, uniformity_keys, barcode_labels, complexity_by_position, complexity_by_sample, mtag_population, probability_list, & all_total_reads, & all_unique_reads, file_handles, options) # try dumping data
    #print "current key", current_read.mip_key
    current_read.enter_read(consensus_details, observed_sites)
  if not options.collapse_free:
    manage_sites(observed_sites, consensus_details, "done:0-0/0,0/0", uniformity_keys, barcode_labels, complexity_by_position, complexity_by_sample, mtag_population, probability_list, & all_total_reads, & all_unique_reads, file_handles, options)
  if(options.mip_file != None):
    output_mip_complexity(complexity_by_position, file_handles["mipwise_summary"])
  if(options.barcode_file != None):
    output_sample_complexity(complexity_by_sample, barcode_labels, file_handles["samplewise_summary"])
  if not options.collapse_free:
    output_overall_complexity(& all_unique_reads, & all_total_reads, file_handles["samplewise_summary"])
  file_handles["notes"].write("%i reads skipped due to barcode or tag filters\n" % reads_skipped)
  file_handles["notes"].write("%i reads with softclipping\n" % softclippings)
  file_handles["notes"].write("%i reads unmapped\n" % reads_unmapped)
  file_handles["notes"].write("%i reads off target\n" % off_target_reads)
  file_handles["notes"].write("%i reads have rejected SAM flags" % improper_pairs)
  for fh in file_handles.values():
    fh.close()

if __name__ == "__main__":

  sys.stderr.write("use other file")
  sys.exit()
