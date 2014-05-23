# Written by Evan Boyle
# boylee [at] uw.edu
import sys
import os
import datetime
import re
from collections import OrderedDict

a = "" # answer from standard in
print "you will need a mip design file, directories to sequence data, and likely a barcode file of\n\
<label><tab><index sequence> to use this script: use the expected index read sequence rather\n\
than the nucleotides of the oligo (reverse complement as needed)"
merge = False
output_lines = []
memory_lines = []
variables = OrderedDict({"script_dir" : "<insert directory here>"})

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
  memory_lines.append(qdir + "\n")
  variables["qdir"] = qdir
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
  variables["fqhead"] = fqhead
  variables["qseq_fqprefix"] = fqprefix
  output_lines.append("python ${script_dir}mipgen_qseq2fq.py $fqhead ${qdir}" + q1pattern + " -o ${qseq_fqprefix}.r1 &\n")
  read1 = fqprefix + ".r1.fq.gz"
  output_lines.append("pid1=$!\n")
  output_lines.append("python ${script_dir}mipgen_qseq2fq.py $fqhead ${qdir}" + q2pattern + " -o ${qseq_fqprefix}.i &\n")
  indexed = True
  output_lines.append("pid2=$!\n")
  indexread = fqprefix + ".i.fq.gz"
  variables["indexread"] = indexread
  output_lines.append("echo \"waiting for qseq to fastq generation to finish\"\n")
  if not se:
    output_lines.append("python ${script_dir}mipgen_qseq2fq.py $fqhead ${qdir}" + q3pattern + " -o ${qseq_fqprefix}.r2 &\n")
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
  if indexed:
    variables["indexread"] = indexread
  if not se:
    a = ""
    while not a.startswith('y') and not a.startswith('n'):
      print "do you want to merge reads with PEAR? [y/n]"
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
print "how many threads would you like to run? [blank not to use multithreading]"
threads = sys.stdin.readline().strip()
memory_lines.append(threads + "\n")
if len(threads) > 0:
  try:
    int(threads)
  except:
    sys.stderr.write("multithreading requires integer, quitting")
    sys.exit()
  variables["threads"] = threads
barcodes = "\\"
while len(barcodes) != 0 and not os.path.exists(barcodes):
  print "what is the path to the file of index sequences you would like to select? [blank not to select]"
  barcodes = sys.stdin.readline().rstrip()
memory_lines.append(barcodes + "\n")
if len(barcodes) != 0:
  variables["barcodes"] = barcodes
  bfh = open(barcodes, 'r')
  sample_bline = bfh.readline()
  bfh.close()
  barcode_length = str(len(sample_bline.strip().split()[1]))
a = ""
while not a.startswith('y') and not a.startswith('n'):
  print "do you want to generate read groups automatically? [y/n]"
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
variables["fqprefix"] = fqprefix
if merge:
  if indexed:
    fqprefix = fqprefix + ".barcoded"
    variables["fqprefix"] = fqprefix
    variables["pe_premerge_read1"] = read1
    variables["pe_premerge_read2"] = read2
    output_lines.append("python ${script_dir}mipgen_fq_cutter_pe.py $pe_premerge_read1 $pe_premerge_read2 " + \
    "-i $indexread " + \
    ("-tb $barcodes -j " + barcode_length + " " if len(barcodes) > 0 else "") + \
    "-o $fqprefix &&\n")
    read1 = fqprefix + ".r1.indexed.fq"
    read2 = fqprefix + ".r2.indexed.fq"
  if re.search(".gz$", read1):
    output_lines.append("gunzip $read1\n")
    read1 = re.sub(".gz$", "", read1)
  if re.search(".gz$", read2):
    output_lines.append("gunzip $read2\n")
    read2 = re.sub(".gz$", "", read2)
  variables["pear_read1"] = read1
  variables["pear_read2"] = read2
  output_lines.append("pear " + ("-j $threads " if len(threads) > 0 else "") + "-f $pear_read1 -r $pear_read2 -o $fqprefix &&\n")
  se = True
  read1 = fqprefix + ".assembled.fastq"
mips = ""
while not os.path.exists(mips):
  print "what is the path to the mip design file? (validity is checked)"
  mips = sys.stdin.readline().rstrip()
memory_lines.append(mips + "\n")
variables["mips"] = mips
mips_fh = open(mips)
header_fields = mips_fh.readline().split()
testmip_fields = mips_fh.readline().split()
mips_fh.close()
seq_index = header_fields.index("mip_sequence") if "mip_sequence" in header_fields else 14
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
if int(ext_tag_size) != 0 or int(lig_tag_size) != 0:
  variables["ext_tag_size"] = ext_tag_size
  variables["lig_tag_size"] = lig_tag_size
variables["cut_read1"] = read1
if se:
  output_lines.append("python ${script_dir}mipgen_fq_cutter_se.py $cut_read1 " + \
    ("-i $indexread " if indexed and not merge else "") + \
    ("-tb $barcodes " if len(barcodes) > 0 else "") + \
    ("-j " + barcode_length + " " if len(barcodes) > 0 and indexed and not merge else "") + \
    ("-m ${lig_tag_size},${ext_tag_size} " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o $fqprefix &&\n")
else:
  variables["cut_read2"] = read2
  output_lines.append("python ${script_dir}mipgen_fq_cutter_pe.py $cut_read1 $cut_read2 " + \
    ("-i $indexread " if indexed else "") + \
    ("-tb $barcodes " if len(barcodes) > 0 else "") + \
    ("-j " + barcode_length + " " if len(barcodes) > 0 and indexed else "") + \
    ("-m ${lig_tag_size},${ext_tag_size} " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o $fqprefix &&\n")
gref = ""
while not os.path.exists(gref):
  print "where is your genome reference?"
  gref = sys.stdin.readline().rstrip()
memory_lines.append(gref + "\n")
variables["gref"] = gref
if se:
  output_lines.append("bwa mem " + \
    ("-R \"" + readgroup + "\" " if not autoreadgroup else "") + \
    ("-t $threads " if len(threads) > 0 else  "") + \
    "$gref ${fqprefix}.indexed.fq > ${fqprefix}.indexed.sam &&\n")
else:
  output_lines.append("bwa mem " + \
    ("-R \"" + readgroup + "\" " if not autoreadgroup else "") + \
    ("-t $threads " if len(threads) > 0 else  "") + \
    "$gref ${fqprefix}.r1.indexed.fq ${fqprefix}.r2.indexed.fq > ${fqprefix}.indexed.sam &&\n")
bamprefix = fqprefix + ".indexed.sort"
variables["bamprefix"] = bamprefix
output_lines.append("samtools view -bS ${fqprefix}.indexed.sam | samtools sort - $bamprefix &&\n")
a = ""
while(not a.startswith('y') and not a.startswith('n')):
  print "do sample indices need to be split into separate files? [y/n]"
  a = sys.stdin.readline()
memory_lines.append(a)
split_by_barcode = a.startswith('y')
a = ""
while (not a.startswith('y') and not a.startswith('n')):
  print "do you want to skip collapsing smmip tags? [y/n]"
  a = sys.stdin.readline()
memory_lines.append(a)
uncollapse_tags = a.startswith('y')
output_lines.append("samtools view -h $bamprefix.bam | python ${script_dir}mipgen_smmip_collapser.py " + str(int(ext_tag_size) + int(lig_tag_size)) + " $bamprefix.collapse " + \
  "-m $mips " +\
  "-f 1 " + \
  ("-c " if not split_by_barcode else "") + \
  ("-r " if autoreadgroup else "") + \
  ("-b $barcodes " if len(barcodes) > 0 else "") + \
  ("-w " if uncollapse_tags else "") + \
  ("-s" if se else "") + " &&\n")
output_lines.append("echo \"analysis commands have terminated (successfully or otherwise)\"\n")
output = open("mipgen_analysis_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".sh", 'w')
memory = open("mipgen_memory_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".txt", 'w')
while len(variables) != 0:
  var, val = variables.popitem(last=False)
  output.write(var + "=" + val + "\n")
for line in output_lines:
  output.write(line)
output.close()
for line in memory_lines:
  memory.write(line)
memory.close()
print "done building"
print "check that PEAR, BWA and SAMtools can be called"

