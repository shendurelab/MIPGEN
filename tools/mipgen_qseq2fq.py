import sys
import re
import gzip
from optparse import OptionParser

parser = OptionParser("%prog readgroup_name s_[lane]_[read]_[etc].qseq (x n) [options]\nfiles are concatenated in fastq output\ncontact: boylee@uw.edu")
parser.add_option("-o", "--output_prefix", dest="output_prefix", type="str", help="directs output to given path")
parser.add_option("-z", "--gzip_input", action="store_true", dest="gzip_input", default=False, help="input qseq files are gzipped, not necessary if file ends in '.gz'")
parser.add_option("-Z", "--do_not_gzip_output",action="store_false", dest="gzip_output", default=True, help="do not print fastqs in gzipped format")
(options, args) = parser.parse_args()

if options.output_prefix != None:
  outfqname = options.output_prefix + ".fq"
else:
  if sys.argv[1].find("_1101.qseq.txt") == -1:
    outfqname = sys.argv[1] + ".fq"
  else:
    outfqname = sys.argv[1].replace("_1101.qseq.txt" + ".fq")

def print_with_handle(fq_out):
  for filename in sys.argv[2:]:
    if filename.startswith('-'):
      break
    if(options.gzip_input or re.search(".gz$", filename)):
      qseq_in = gzip.open(filename, 'rb')
    else:
      qseq_in = open(filename, 'r')
    for qseq_line in qseq_in:
      qseq_fields = qseq_line.rstrip().split('\t')
      qseq_fields[8] = qseq_fields[8].replace(".", "N")
      fq_read_cluster = "@" + sys.argv[1] + ":" + qseq_fields[2] + ":" + qseq_fields[3] + ":" + qseq_fields[4] + ":" + qseq_fields[5] + ":" + qseq_fields[10]
      first_line = fq_read_cluster + "\n"
      second_line = qseq_fields[8] + "\n"
      third_line = "+\n"
      fourth_line = ""
      for quality in qseq_fields[9].rstrip():
        fourth_line += chr(ord(quality) - 31)
      fourth_line += "\n"   
      fq_out.write(first_line + second_line + third_line + fourth_line)
    qseq_in.close()

'''
subfields
:sequencing run
:lane number
:tile number
:x coord
:y coord
:pass filter
'''
if options.gzip_output:
  fq_out = gzip.open(outfqname + ".gz", 'wb')
else:
  fq_out = open(outfqname, 'w')
print_with_handle(fq_out)
fq_out.close()
sys.stderr.write( "#conversion to fastqs completed\n" )
