import sys
import re
import gzip
from optparse import OptionParser

parser = OptionParser("%prog read_1.fq read_2.fq [options]")
parser.add_option("-o", "--output_prefix", dest="output_prefix", type="str", help="directs output to given path")
parser.add_option("-i", "--index_read", dest="index_file", type="str", help="index read file to insert barcode into header")
parser.add_option("-j", "--index_length", dest="index_length", type="int", help="truncate index read to given length")
parser.add_option("-b", "--barcode_file", type="str", help="select barcodes in file of label<tab>sequence")
parser.add_option("-t", "--tolerant", action="store_true", default=False, help="allow 1bp edits from provided barcodes in file")
parser.add_option("-m", "--molecular_tag", dest="molecular_tag", type="str", help="molecular tag, template: \"-m 10,0\" for read 1, 10bp")
parser.add_option("-l", "--truncate_read_length", dest="truncated_read_length", type="int", help="truncate reads to given length")
parser.add_option("-L", "--skip_length", dest="skip_length", type="int", help="skip the number of bases provided (post truncation)")
parser.add_option("-d", "--discard", action="store_true", dest="discard", default=False, help="discard read pairs for which one of the reads is all '#' quality")
parser.add_option("-e", "--partial_discard", action="store_true", dest="partial_discard", default=False, help="enables truncation of reads to eliminate tailing bases of '#' quality")
parser.add_option("-z", "--gzip_input", action="store_true", dest="gzip_input", default=False, help="input fastq files are gzipped; done by default if file ends in \".gz\"")

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
 
if options.molecular_tag != None:
 molecular_tag_specs = options.molecular_tag.split(",")
 molecular_tag_specs = [int(entry) for entry in molecular_tag_specs]
if options.gzip_input or re.search(".gz$", sys.argv[1]):
 first_in = gzip.open(sys.argv[1], 'rb')
else:
 first_in = open(sys.argv[1], 'r')
if options.output_prefix != None:
 first_outfq = options.output_prefix + ".r1.indexed.fq"
else:
 first_outfq = sys.argv[1] + ".r1.indexed.fq"
if options.output_prefix != None:
  second_outfq = options.output_prefix + ".r2.indexed.fq"
else:
  second_outfq = sys.argv[2] + ".r2.indexed.fq"
if options.gzip_input or re.search(".gz$", sys.argv[2]):
 second_in = gzip.open(sys.argv[2], 'rb')
else:
 second_in = open(sys.argv[2], 'r')
if(options.index_file != None):
 if options.gzip_input or re.search(".gz$", options.index_file):
  i_in = gzip.open(options.index_file, 'rb')
 else:
  i_in = open(options.index_file, 'rb')
else:
  i_in = None
with open(first_outfq, 'w') as first_out:
 with open(second_outfq, 'w') as second_out:
  while(1):   
   index_block = []
   first_block = []
   second_block = []
   try:
    for i in range(4):
     first_line = first_in.next()
     first_block.append(first_line)
     second_line = second_in.next()
     second_block.append(second_line)
   except StopIteration:
    break
   if(options.discard and (re.match("^#+$", first_block[3].rstrip()) or re.match("^#+$", second_block[3].rstrip()))):
    continue
   if options.index_file:
    for i in range(4):
     index_line = i_in.next()
     index_block.append(index_line)
    barcode = index_block[1].rstrip()
   else:
    barcode = "N"
   if(options.index_file != None and options.index_length != None):
    barcode = barcode[:options.index_length]
   if(options.barcode_file != None and barcode not in used_barcodes):
     continue
   if(options.index_file == None and re.search("#[ATGCN]+$", first_block[0]) and re.search("#[ATGC]+$", second_block[0])):
    pass
   else:
    first_block[0] = first_block[0].replace("\n", "#" + barcode + "\n")
    first_block[0] = first_block[0].replace(" ", "_")
    second_block[0] = second_block[0].replace("\n", "#" + barcode + "\n")
    second_block[0] = second_block[0].replace(" ", "_")
   if(options.molecular_tag != None):
    tag = first_block[1][:molecular_tag_specs[0]] + second_block[1][:molecular_tag_specs[1]]
    first_block[0] = first_block[0].replace("\n", "-" + tag + "\n")
    first_block[1] = first_block[1][molecular_tag_specs[0]:]
    first_block[3] = first_block[3][molecular_tag_specs[0]:]
    second_block[0] = second_block[0].replace("\n", "-" + tag + "\n")
    second_block[1] = second_block[1][molecular_tag_specs[1]:]
    second_block[3] = second_block[3][molecular_tag_specs[1]:]
   if(options.truncated_read_length != None): 
    first_block[1] = first_block[1][:options.truncated_read_length].rstrip() + "\n"
    first_block[3] = first_block[3][:options.truncated_read_length].rstrip() + "\n"
    second_block[1] = second_block[1][:options.truncated_read_length].rstrip() + "\n"
    second_block[3] = second_block[3][:options.truncated_read_length].rstrip() + "\n"
   if(options.partial_discard):
    match = re.search("#+$", first_block[3])
    if(match != None):
     first_block[1] = first_block[1][0:match.start()] + "\n"
     first_block[3] = first_block[3][0:match.start()] + "\n"
    match = re.search("#+$", second_block[3])
    if(match != None):
     second_block[1] = second_block[1][0:match.start()] + "\n"
     second_block[3] = second_block[3][0:match.start()] + "\n"
   if(options.skip_length != None):
    first_block[1] = first_block[1][options.skip_length:].rstrip() + "\n"
    first_block[3] = first_block[3][options.skip_length:].rstrip() + "\n"
    second_block[1] = second_block[1][options.skip_length:].rstrip() + "\n"
    second_block[3] = second_block[3][options.skip_length:].rstrip() + "\n"
   for line in first_block:
    first_out.write(line)
   for line in second_block:
    second_out.write(line)

first_in.close()
second_in.close()

if(options.index_file != None):
  i_in.close()

if(options.output_prefix != None):
 samfile = options.output_prefix + ".indexed.sam"
else:
 f_shared_index = 0
 if(sys.argv[1][-1] == sys.argv[2][-1]):
  r_shared_index = -1
 else:
  r_shared_index = None
 while(sys.argv[1][f_shared_index] == sys.argv[2][f_shared_index]):
  f_shared_index += 1
 while(sys.argv[1][r_shared_index - 1] == sys.argv[2][r_shared_index - 1]):
  r_shared_index -= 1
 root_component = sys.argv[1][:f_shared_index]
 derived_component = sys.argv[1][f_shared_index:r_shared_index] + sys.argv[2][f_shared_index:r_shared_index]
 samfile = root_component + derived_component + ".indexed.sam"

print "suggested output file for alignment: " + samfile
