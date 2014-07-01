# Written by Evan Boyle
# boylee [at] uw.edu
import sys
import re
import gzip
from optparse import OptionParser

parser = OptionParser("%prog read.fq [options]")
parser.add_option("-o", "--output_prefix", type="str", help="directs output to given path")
parser.add_option("-i", "--index_read", dest="index_file", type="str", help="index read file")
parser.add_option("-j", "--index_length", type="int", help="truncate index read to given length")
parser.add_option("-b", "--barcode_file", type="str", help="select barcodes in file of label<tab>sequence")
parser.add_option("-t", "--tolerant", action="store_true", default=False, help="allow 1bp edits from provided barcodes in file")
parser.add_option("-m", "--molecular_tag", type="str", help="molecular tag, template: \"-m 3,2\" for first 3 bases, last 2 bases")
parser.add_option("-l", "--truncated_read_length", type="int", help="truncate reads to given length")
parser.add_option("-L", "--skip_length", type="int", help="skip the number of bases provided (post truncation)")
parser.add_option("-d", "--discard", action="store_true",  default=False, help="discard read pairs for which one of the reads is all '#' quality")
parser.add_option("-e", "--partial_discard", action="store_true", default=False, help="enables truncation of reads to eliminate tailing bases of '#' quality")
parser.add_option("-z", "--gzip_input", action="store_true", default=False, help="input fastq files are gzipped; done by default if file ends in \".gz\"")

(options, args) = parser.parse_args()

if(options.partial_discard and not options.discard):
 print "-e option requires -d option"
if(options.index_file == None and options.index_length != None):
 print "-j option requires -i option"

used_barcodes = set()
bases = "ATCGN"

if options.barcode_file != None:
 with open(options.barcode_file, 'r') as b_in:
  for barcode_line in b_in:
   (barcode_label, barcode_seq) = barcode_line.rstrip().split()
   used_barcodes.add(barcode_seq)
   if(options.tolerant):
    for i in range(len(barcode_seq)):
     native_seq = list(barcode_seq)
     for j in range(5):
      mutated_seq = native_seq
      mutated_seq[i] = bases[j]
      used_barcodes.add("".join(mutated_seq))

if options.molecular_tag != None :
 molecular_tag_specs = options.molecular_tag.split(",")
 molecular_tag_specs = [int(entry) for entry in molecular_tag_specs]
if options.gzip_input or re.search(".gz$", sys.argv[1]):
 fq_in = gzip.open(sys.argv[1], 'rb')
else:
 fq_in = open(sys.argv[1], 'r')
if(options.output_prefix != None):
 outfq = options.output_prefix + ".indexed.fq"
else:
 outfq = sys.argv[1] + ".indexed.fq"
if options.index_file != None:
 if options.gzip_input or re.search(".gz$", options.index_file):
  i_in = gzip.open(options.index_file, 'rb')
 else:
  i_in = open(options.index_file, 'r')
else:
  i_in = None
with open(outfq, 'w') as out:
 while(1):   
  index_block = []
  block = []
  try:
   for i in range(4):
    line = fq_in.next()
    block.append(line)
  except StopIteration:
   break
  block[0] = re.sub("/\d$","",block[0])
  barcode_in_header = re.search("#([ATGCN]+)(-[ATGCN]+)?$", block[0])
  if options.index_file != None:
   for i in range(4):
    index_line = i_in.next()
    index_block.append(index_line)
   barcode = index_block[1].rstrip()
  elif barcode_in_header:
   barcode = barcode_in_header.group(1)
  else:
   barcode = "N"
  if options.index_file != None and options.index_length != None:
   barcode = barcode[:options.index_length]
  if options.barcode_file != None and barcode not in used_barcodes:
   continue
  if options.index_file == None and barcode_in_header:
   pass
  else:
   block[0] = block[0].replace("\n", "#" + barcode + "\n")
   block[0] = block[0].replace(" ", "_")
  tag_in_header = re.search("#[ATGCN]+-[ATGCN]+$", block[0])
  if(options.molecular_tag != None and not tag_in_header):
   tag = block[1][:molecular_tag_specs[0]] + block[1][len(block[1]) - molecular_tag_specs[1] - 1 : -1]
   block[0] = block[0].replace("\n", "-" + tag + "\n")
   block[1] = block[1].rstrip()[molecular_tag_specs[0]:len(block[1]) - 1 - molecular_tag_specs[1]] + "\n"
   block[3] = block[3].rstrip()[molecular_tag_specs[0]:len(block[3]) - 1 - molecular_tag_specs[1]] + "\n"
  if(options.truncated_read_length != None): 
   block[1] = block[1][:options.truncated_read_length].rstrip() + "\n"
   block[3] = block[3][:options.truncated_read_length].rstrip() + "\n"
  if(options.partial_discard):
   match = re.search("#+$", block[3])
   if(match != None):
    block[1] = block[1][0:match.start()] + "\n"
    block[3] = block[3][0:match.start()] + "\n"
  if(options.skip_length != None):
   block[1] = block[1][options.skip_length:].rstrip() + "\n"
   block[3] = block[3][options.skip_length:].rstrip() + "\n"
  for line in block:
   out.write(line)
fq_in.close()

if options.index_file != None:
  i_in.close()

if options.output_prefix != None:
 samfile = options.output_prefix + ".indexed.sam"
else:
 suffix_match = re.search("fq$", sys.argv[1])
 if suffix_match == None:
  suffix_match = re.search("fastq$", sys.argv[1])
 if suffix_match != None:
  samfile = sys.argv[1][:suffix_match.start()] + "indexed.sam"
 else:
  samfile = sys.argv[1] + ".indexed.sam"

sys.stderr.write( "#fq cutting finished\n" )
sys.stderr.write( "#suggested output file for alignment: " + samfile + "\n")
