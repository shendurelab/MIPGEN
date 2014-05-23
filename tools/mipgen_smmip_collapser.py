# Written by Evan Boyle
# boylee [at] uw.edu
import sys
import re
import numpy as np
from scipy import optimize as optimize
from random import choice
from optparse import OptionParser
from string import maketrans
from genome_sam_collapser import *

if __name__ == "__main__":

  parser = OptionParser("%prog (STDIN = coordinate_sorted_file.sam) tag_size output_prefix [options]")
  parser.add_option("-b", "--split_barcodes", dest="barcode_file", type="str", help="pulls out only reads with exactly matching barcodes in provided file of label<tab>sequence")
  parser.add_option("-c", "--merge_samples", action="store_true", default=False, help="selects barcodes but does not split into separate files")
  parser.add_option("-d", "--dual_indexed", action="store_true", dest="dual_indexed",  default=False, help="reads barcode file as dual indexed, i.e., with two columns of barcodes")
  parser.add_option("-p", "--picky", action="store_true", dest="filter_molecular_tags", default=False, help="discards reads with non ATCGs in the molecular tag")
  parser.add_option("-t", "--tolerant", action="store_true", dest="allow_ambiguous_barcodes", default=False, help="allows barcodes to be matched with 1bp edit distance")
  parser.add_option("-m", "--mip_design_file", dest="mip_file", type="str", help="only pulls out sequences that are within 2bp of mip sites as determined by mip design file")
  parser.add_option("-n", "--mip_reference", action="store_true", default=False, help="uses chromosome SAM field as MIP key")
  parser.add_option("-C", "--confidence_level", dest="confidence_level", type="float", default=0.9, help="controls consensus calling: confidence refers to the chance of a tag truly representing one distinct haplotype -- high confidence leads to more random sampling to reduce the chances of chimeric consensus and low confidence leads to indiscriminate consensus calling, number refers to probability of ALL site-, barcode-, and tag-stratified reads representing unique captures for that site and barcode sequence (default is 0.9)")
  parser.add_option("-T", "--no_trimming", action="store_true", dest="no_trimming", default=False, help="do not remove number of bases corresponding to mip arm sequences even if mip file is provided")
  parser.add_option("-r", "--add_or_replace_readgroups", action="store_true", default=False, help="use the barcode file (if given) or barcode sequence to generate read groups")
  parser.add_option("-f", "--flex_space", dest="flex_space", type="int", default=0, help="searches given number of bases on either side of read start when looking to assign a read to a known MIP target")
  parser.add_option("-s", "--single_end", action="store_true", default=False, help="single end run")
  parser.add_option("-S", "--no_softclip_filtering", action="store_false", dest="filter_softclips", default=True, help="retains reads with softclipping at the beginning of the read")
  parser.add_option("-w", "--collapse_free", action="store_true", default=False, help="do not run collapsing -- only trim and partition reads")
  parser.add_option("-x", "--exact_arms", action="store_true", default=False, help="only accept MIP reads with exact arm matches, default accepts any read at correct position")

  options, args = parser.parse_args()
 
  if options.merge_samples and not options.barcode_file:
    sys.stderr.write("option 'c' requires option 'b'")
    sys.exit()
  if options.add_or_replace_readgroups and not options.barcode_file:
    sys.stderr.write("option 'r' requires option 'b'")
    sys.exit()
  if options.exact_arms and not options.mip_file:
    sys.stderr.write("option 'x' requires option 'm'")
    sys.exit()
  initialize_and_iterate(options)
  sys.stderr.write("collapsing has terminated\n")
