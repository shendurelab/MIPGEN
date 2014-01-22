import sys
import os
import datetime
import re

a = "" # answer from standard in
print "you will need a mip design file, directories to sequence data, and likely a barcode file of\n\
<label><tab><index sequence> to use this script: use the expected index read sequence rather\n\
than the nucleotides of the oligo (reverse complement as needed)"
merge = False
output_lines = []
memory_lines = []
while(not a.startswith('y') and not a.startswith('n')):
  print "are you starting from qseqs? [y/n]"
  a = sys.stdin.readline()
  memory_lines.append(a)
if a.startswith('y'):
  qdir = ""
  while not os.path.exists(qdir):
    print "where is your directory of qseqs? (validity is checked)"
    qdir = sys.stdin.readline()
    qdir = qdir.rstrip('\n/') + '/'
  memory_lines.append(qdir)
  print "what pattern do your read 1 files match? [blank to accept default]"
  print "default: s_1_1*qseq.txt"
  q1patterndefault = "s_1_1*qseq.txt"
  q1pattern = sys.stdin.readline().rstrip()
  memory_lines.append(q1pattern + "\n")
  q1pattern = q1pattern if len(q1pattern) > 0 else q1patterndefault
  print "what pattern do your index read files match? [blank to accept default]"
  print "default: s_1_2*qseq.txt"
  q2patterndefault = "s_1_2*qseq.txt"
  q2pattern = sys.stdin.readline().rstrip()
  memory_lines.append(q2pattern + "\n")
  q2pattern = q2pattern if len(q2pattern) > 0 else q2patterndefault
  print "what pattern do your read 2 files match? [type se if single end, blank to accept default]"
  print "default: s_1_3*qseq.txt"
  q3patterndefault = "s_1_3*qseq.txt"
  q3pattern = sys.stdin.readline().rstrip()
  memory_lines.append(q3pattern + "\n")
  q3pattern = q3pattern if len(q3pattern) > 0 else q3patterndefault
  index = 0
  se = q3pattern == "se"
  if not se:
    a = ""
    while(not a.startswith('y') and not a.startswith('n')):
      print "do you want to merge reads with PEAR? [y/n]"
      a = sys.stdin.readline()
    memory_lines.append(a)
    merge = a.startswith('y')
  for x,y,z in zip(q1pattern, q2pattern, q3pattern):
    if not (x == y and (x == z or se)):
      break
    index += 1
  fqdefault = qdir + q1pattern[:index] if index > 0 else qdir + "converted"
  fqdefault = re.sub("[\*\?]+.*", "", fqdefault).rstrip("_")
  print "do you want to specify output? [blank to accept default]"
  print "default: " + fqdefault
  fqprefix = sys.stdin.readline().rstrip()
  memory_lines.append(fqprefix + "\n")
  fqprefix = fqprefix if len(fqprefix) > 0 else fqdefault
  print "what would you like to insert into the fastq header? [blank to accept default]"
  fqheaddefault = fqdefault.replace(qdir, "", 1) + "_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
  print "default: " + fqheaddefault
  fqhead = sys.stdin.readline().rstrip()
  memory_lines.append(fqhead + "\n")
  fqhead = fqhead if len(fqhead) > 0 else fqheaddefault
  output_lines.append("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q1pattern + " -o " + fqprefix + ".r1 &\n")
  read1 = fqprefix + ".r1.fq.gz"
  output_lines.append("pid1=$!\n")
  output_lines.append("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q2pattern + " -o " + fqprefix + ".i &\n")
  indexread = fqprefix + ".i.fq.gz"
  indexed = True
  output_lines.append("pid2=$!\n")
  output_lines.append("echo \"waiting for qseq to fastq generation to finish\"\n")
  if not se:
    output_lines.append("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q3pattern + " -o " + fqprefix + ".r2 &\n")
    output_lines.append("pid3=$!\n")
    output_lines.append("wait $pid3\n")
    read2 = fqprefix + ".r2.fq.gz"
  output_lines.append("wait $pid1\n")
  output_lines.append("wait $pid2\n")
  output_lines.append("echo \"done waiting\"\n")
else:
  print "what is the path to read 1?"
  read1 = sys.stdin.readline().rstrip()
  memory_lines.append(read1 + "\n")
  print "what is the path to read 2? [blank if single end]"
  read2 = sys.stdin.readline().rstrip()
  memory_lines.append(read2 + "\n")
  se = len(read2) == 0
  print "what is the path to the index read? [blank if demultiplexed]"
  indexread = sys.stdin.readline().rstrip()
  memory_lines.append(indexread + "\n")
  indexed = len(indexread) != 0
  if not se:
    a = ""
    while not a.startswith('y') and not a.startswith('n'):
      print "do you want to merge reads with PEAR? (files cannot be gzipped) [y/n]"
      a = sys.stdin.readline()
    memory_lines.append(a)
    merge = a.startswith('y')
  outdir = ""
  while not os.path.exists(outdir):
    print "in what directory do you want to place output? (validity is checked)"
    outdir = sys.stdin.readline().rstrip('\n/') + '/'
  memory_lines.append(outdir + "\n")
  print "what file prefix would you like to use for this analyis?"
  fileprefix = sys.stdin.readline().rstrip()
  memory_lines.append(fileprefix + "\n")
  fqprefix = outdir + fileprefix
a = ""
barcodes = "\\"
while len(barcodes) != 0 and not os.path.exists(barcodes):
  print "what is the path to the file of index sequences you would like to select? [blank not to select]"
  barcodes = sys.stdin.readline().rstrip()
memory_lines.append(barcodes + "\n")
a = ""
while not a.startswith('y') and not a.startswith('n'):
  print "do you want to generate read groups automatically?"
  a = sys.stdin.readline()
memory_lines.append(a)
autoreadgroup = a.startswith('y')
if not autoreadgroup:
  readgroup = ""
  while not readgroup.startswith("@RG\tID:"):
    print "what do you want to use as the read group header?\n\
    format = \"@RG<tab>ID:<etc>\""
    readgroup = sys.stdin.readline().rstrip()
  memory_lines.append(readgroup + "\n")
if merge:
  if indexed:
    fqprefix = fqprefix + ".barcoded"
    output_lines.append("python mipgen_fq_cutter_pe.py " + read1 + " " + read2 + " " + \
    "-i " + indexread + " " + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    "-o " + fqprefix + "\n")
    read1 = fqprefix + ".r1.indexed.fq"
    read2 = fqprefix + ".r2.indexed.fq"
  output_lines.append("pear -f " + read1 + " -r " + read2 + " -o " + fqprefix + " &&\n")
  se = True
  read1 = fqprefix + ".assembled.fastq"
mips = ""
while not os.path.exists(mips):
  print "what is the path to the mip design file? (validity is checked)"
  mips = sys.stdin.readline().rstrip()
memory_lines.append(mips + "\n")
mips_fh = open(mips)
header_fields = mips_fh.readline().split()
testmip_fields = mips_fh.readline().split()
mips_fh.close()
seq_index = header_fields.index("mip_sequence")
m = re.search("(N*)CTTCAGCTTCCCGATATCCGACGGTAGTGT(N*)",testmip_fields[seq_index])
ext_tag_size = 0 if m == None else str(len(m.group(2)))
lig_tag_size = 0 if m == None else str(len(m.group(1)))
print "what is your extension arm smmip tag length? [blank to use detected length]"
print "detected: " + str(ext_tag_size)
ext_tag_input = sys.stdin.readline().rstrip()
memory_lines.append(ext_tag_input + "\n")
if len(ext_tag_input) > 0:
  ext_tag_size = ext_tag_input
print "what is your ligation arm smmip tag length? [blank to use detected length]"
print "detected: " + str(lig_tag_size)
lig_tag_input = sys.stdin.readline().rstrip()
memory_lines.append(lig_tag_input + "\n")
if len(lig_tag_input) > 0:
  lig_tag_size = lig_tag_input
if se:
  output_lines.append("python mipgen_fq_cutter_se.py " + read1 + " " + \
    ("-i " + indexread + " " if indexed and not merge else "") + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    ("-m " + lig_tag_size + "," + ext_tag_size + " " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o " + fqprefix + " &&\n")
else:
  output_lines.append("python mipgen_fq_cutter_pe.py " + read1 + " " + read2 + " " + \
    ("-i " + indexread + " " if indexed else "") + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    ("-m " + lig_tag_size + "," + ext_tag_size + " " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o " + fqprefix + " &&\n")
gref = ""
while not os.path.exists(gref):
  print "where is your genome reference?"
  gref = sys.stdin.readline().rstrip()
memory_lines.append(gref + "\n")
if se:
  output_lines.append("bwa mem " + ("-R \"" + readgroup + "\" " if not autoreadgroup else "") + gref + " " + fqprefix + ".indexed.fq > " + fqprefix + ".indexed.sam &&\n")
else:
  output_lines.append("bwa mem " + ("-R \"" + readgroup + "\" " if not autoreadgroup else "") + gref + " " + fqprefix + ".r1.indexed.fq " + fqprefix + ".r2.indexed.fq > " + fqprefix + ".indexed.sam &&\n")
bamprefix = fqprefix + ".indexed.sort"
output_lines.append("samtools view -bS " + fqprefix + ".indexed.sam | samtools sort - " + bamprefix + " &&\n")
a = ""
while(not a.startswith('y') and not a.startswith('n')):
  print "do sample indices need to be split into separate files? [y/n]"
  a = sys.stdin.readline()
memory_lines.append(a)
split_by_barcode = a.startswith('y')
output_lines.append("samtools view -h " + bamprefix + ".bam | python mipgen_smmip_collapser.py " + str(int(ext_tag_size) + int(lig_tag_size)) + " " + bamprefix + ".collapse " + \
  "-m " + mips + " " +\
  "-f 1 " + \
  ("-c " if not split_by_barcode else "") + \
  ("-r " if autoreadgroup else "") + \
  ("-b " + barcodes + " " if len(barcodes) > 0 else "") + \
  ("-s" if se else "") + "\n")
output_lines.append("echo \"analysis commands have terminated\"\n")

output = open("mipgen_analysis_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".sh", 'w')
memory = open("mipgen_memory_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".txt", 'w')
for line in output_lines:
  output.write(line)
output.close()
for line in memory_lines:
  memory.write(line)
memory.close()
print "done building"
print "check that PEAR, BWA and SAMtools can be called"
