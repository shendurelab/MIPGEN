import sys
import re
import gzip
from optparse import OptionParser

parser = OptionParser("%prog read.fq [options]")
parser.add_option("-o", "--output_prefix", type="str", help="directs output to given path")
parser.add_option("-i", "--index_read", dest="index_file", type="str", help="truncate index read to given length")
parser.add_option("-j", "--index_length", type="int", help="truncate index read to given length")
parser.add_option("-m", "--molecular_tag", type="str", help="molecular tag, template: \"-m 3,2\" for first 3 bases, last 2 bases")
parser.add_option("-l", "--truncated_read_length", type="int", help="truncate reads to given length")
parser.add_option("-L", "--skip_length", type="int", help="skip the number of bases provided (post truncation)")
parser.add_option("-d", "--discard", action="store_true",  default=False, help="discard read pairs for which one of the reads is all '#' quality")
parser.add_option("-e", "--partial_discard", action="store_true", default=False, help="enables truncation of reads to eliminate tailing bases of '#' quality")
parser.add_option("-z", "--gzip_input", action="store_true", default=False, help="input fastq files are gzipped")

(options, args) = parser.parse_args()

if(options.partial_discard and not options.discard):
 print "-e option requires -d option"
if(options.index_file == None and options.index_length != None):
 print "-j option requires -i option"

if(options.molecular_tag != None):
 molecular_tag_specs = options.molecular_tag.split(",")
 molecular_tag_specs = [int(entry) for entry in molecular_tag_specs]
if(options.gzip_input):
 fq_in = gzip.open(sys.argv[1], 'rb')
else:
 fq_in = open(sys.argv[1], 'r')
if(options.output_prefix != None):
 outfq = options.output_prefix + ".indexed.fq"
else:
 outfq = sys.argv[1] + ".indexed.fq"
if(options.index_file != None):
 if(options.gzip_input):
  i_in = gzip.open(options.index_file, 'rb')
 else:
  i_in = open(options.index_file, 'rb')
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
  if options.index_file != None:
   for i in range(4):
    index_line = i_in.next()
    index_block.append(index_line)
   barcode = index_block[1].rstrip()
  else:
   barcode = "N"
  if(options.index_file != None and options.index_length != None):
   barcode = barcode[:options.index_length]
  block[0] = block[0].replace("\n", "#" + barcode + "\n")
  block[0] = block[0].replace(" ", "_")
  if(options.molecular_tag != None):
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

if(options.index_file != None):
  i_in.close()

if(options.output_prefix != None):
 samfile = options.output_prefix + ".indexed.sam"
else:
 suffix_match = re.search("fq$", sys.argv[1])
 if suffix_match == None:
  suffix_match = re.search("fastq$", sys.argv[1])
 if suffix_match != None:
  samfile = sys.argv[1][:suffix_match.start()] + "indexed.sam"
 else:
  samfile = sys.argv[1] + ".indexed.sam"

print "suggested output file for alignment: " + samfile
