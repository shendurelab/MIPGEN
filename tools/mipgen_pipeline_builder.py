import sys
import os
import datetime
import re

a = "" # answer from standard in
print "you will need a mip design file, directories to sequence data, and possibly a barcode file of\
<label><tab><index sequence> to use this script"
merge = False
output = open("mipgen_analysis_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".sh", 'w')
memory = open("mipgen_memory_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + ".txt", 'w')
while(not a.startswith('y') and not a.startswith('n')):
  print "are you starting from qseqs? [y/n]"
  a = sys.stdin.readline()
  memory.write(a)
if a.startswith('y'):
  qdir = ""
  while not os.path.exists(qdir):
    print "where is your directory of qseqs? (validity is checked)"
    qdir = sys.stdin.readline()
    qdir = qdir.rstrip('\n/') + '/'
  memory.write(qdir)
  print "what pattern do your read 1 files match? [blank to accept default]"
  print "default: s_1_1*qseq.txt"
  q1patterndefault = "s_1_1*qseq.txt"
  q1pattern = sys.stdin.readline().rstrip()
  memory.write(q1pattern + "\n")
  q1pattern = q1pattern if len(q1pattern) > 0 else q1patterndefault
  print "what pattern do your index read files match? [blank to accept default]"
  print "default: s_1_2*qseq.txt"
  q2patterndefault = "s_1_2*qseq.txt"
  q2pattern = sys.stdin.readline().rstrip()
  memory.write(q2pattern + "\n")
  q2pattern = q2pattern if len(q2pattern) > 0 else q2patterndefault
  print "what pattern do your read 2 files match? [type se if single end, blank to accept default]"
  print "default: s_1_3*qseq.txt"
  q3patterndefault = "s_1_3*qseq.txt"
  q3pattern = sys.stdin.readline().rstrip()
  memory.write(q3pattern + "\n")
  q3pattern = q3pattern if len(q3pattern) > 0 else q3patterndefault
  index = 0
  se = q3pattern == "se"
  if not se:
    a = ""
    while(not a.startswith('y') and not a.startswith('n')):
      print "do you want to merge reads with PEAR? [y/n]"
      a = sys.stdin.readline()
    memory.write(a)
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
  memory.write(fqprefix + "\n")
  fqprefix = fqprefix if len(fqprefix) > 0 else fqdefault
  print "what would you like to insert into the fastq header? [blank to accept default]"
  fqheaddefault = fqdefault.replace(qdir, "", 1) + "_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
  print "default: " + fqheaddefault
  fqhead = sys.stdin.readline().rstrip()
  memory.write(fqhead + "\n")
  fqhead = fqhead if len(fqhead) > 0 else fqheaddefault
  output.write("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q1pattern + " -o " + fqprefix + ".r1 &\n")
  read1 = fqprefix + ".r1.fq.gz"
  output.write("pid1=$!\n")
  output.write("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q2pattern + " -o " + fqprefix + ".i &\n")
  indexread = fqprefix + ".i.fq.gz"
  indexed = True
  output.write("pid2=$!\n")
  if not se:
    output.write("python mipgen_qseq2fq.py " + fqhead + " " + qdir + q3pattern + " -o " + fqprefix + ".r2 &\n")
    output.write("pid3=$!\n")
    output.write("wait $pid3\n")
    read2 = fqprefix + ".r2.fq.gz"
  output.write("wait $pid1\n")
  output.write("wait $pid2\n")
else:
  print "what is the path to read 1?"
  read1 = sys.stdin.readline().rstrip()
  memory.write(read1 + "\n")
  print "what is the path to read 2? [blank if single end]"
  read2 = sys.stdin.readline().rstrip()
  memory.write(read2 + "\n")
  se = len(read2) == 0
  print "what is the path to the index read? [blank if demultiplexed]"
  indexread = sys.stdin.readline().rstrip()
  memory.write(indexread + "\n")
  indexed = len(indexread) != 0
  if not se:
    a = ""
    while not a.startswith('y') and not a.startswith('n'):
      print "do you want to merge reads with PEAR? (files cannot be gzipped) [y/n]"
      a = sys.stdin.readline()
    memory.write(a)
    merge = a.startswith('y')
  outdir = ""
  while not os.path.exists(outdir):
    print "in what directory do you want to place output? (validity is checked)"
    outdir = sys.stdin.readline().rstrip('\n/') + '/'
  memory.write(outdir + "\n")
  print "what file prefix would you like to use for this analyis?"
  fileprefix = sys.stdin.readline().rstrip()
  memory.write(fileprefix + "\n")
  fqprefix = outdir + fileprefix
a = ""
barcodes = "\\"
while len(barcodes) != 0 and not os.path.exists(barcodes):
  print "what is the path to the file of index sequences you would like to select? [blank not to select]"
  barcodes = sys.stdin.readline().rstrip()
memory.write(barcodes + "\n")
if merge:
  if indexed:
    fqprefix = fqprefix + ".barcoded"
    output.write("python mipgen_fq_cutter_pe.py " + read1 + " " + read2 + " " + \
    "-i " + indexread + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    "-o " + fqprefix + "\n")
    read1 = fqprefix + ".r1.indexed.fq"
    read2 = fqprefix + ".r2.indexed.fq"
  output.write("pear -f " + read1 + " -r " + read2 + " -o " + fqprefix + " &&\n")
  se = True
  read1 = fqprefix + ".assembled.fastq"
mips = ""
while not os.path.exists(mips):
  print "what is the path to the mip design file? (validity is checked)"
  mips = sys.stdin.readline().rstrip()
memory.write(mips + "\n")
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
memory.write(ext_tag_input + "\n")
if len(ext_tag_input) > 0:
  ext_tag_size = ext_tag_input
print fqprefix + "is the prefix"
print "what is your ligation arm smmip tag length? [blank to use detected length]"
print "detected: " + str(lig_tag_size)
lig_tag_input = sys.stdin.readline().rstrip()
memory.write(lig_tag_input + "\n")
if len(lig_tag_input) > 0:
  lig_tag_size = lig_tag_input
if se:
  output.write("python mipgen_fq_cutter_se.py " + read1 + " " + \
    ("-i " + indexread + " " if indexed and not merge else "") + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    ("-m " + lig_tag_size + "," + ext_tag_size + " " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o " + fqprefix + " &&\n")
else:
  output.write("python mipgen_fq_cutter_pe.py " + read1 + " " + read2 + " " + \
    ("-i " + indexread + " " if indexed else "") + \
    ("-tb " + barcodes + " " if len(barcodes) > 0 else "") + \
    ("-m " + lig_tag_size + "," + ext_tag_size + " " if int(lig_tag_size) != 0 or int(ext_tag_size) !=0 else "") + \
    "-o " + fqprefix + " &&\n")
gref = ""
while not os.path.exists(gref):
  print "where is your genome reference?"
  gref = sys.stdin.readline().rstrip()
memory.write(gref + "\n")
if se:
  output.write("bwa mem " + gref + " " + fqprefix + ".indexed.fq > " + fqprefix + ".indexed.sam &&\n")
else:
  output.write("bwa mem " + gref + " " + fqprefix + ".r1.indexed.fq " + fqprefix + ".r2.indexed.fq > " + fqprefix + ".indexed.sam &&\n")
bamprefix = fqprefix + ".indexed.sort"
output.write("samtools view -bS " + fqprefix + ".indexed.sam | samtools sort - " + bamprefix + " &&\n")
a = ""
while(not a.startswith('y') and not a.startswith('n')):
  print "do sample indices need to be split into separate files? [y/n]"
  a = sys.stdin.readline()
memory.write(a)
split_by_barcode = a.startswith('y')
output.write("samtools view -h " + bamprefix + ".bam | python mipgen_smmip_collapser.py " + str(int(ext_tag_size) + int(lig_tag_size)) + " " + bamprefix + ".collapse " + \
  "-m " + mips + " " +\
  "-f 1 " + \
  ("-c " if not split_by_barcode else "") + \
  ("-b " + barcodes + " " if len(barcodes) > 0 else "") + \
  ("-s" if se else "") + "\n")
output.write("echo \"analysis commands have terminated\"")
output.close()
  
print "done building"
print "check that PEAR, BWA and SAMtools can be called"
