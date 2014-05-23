# Written by Evan Boyle
# boylee [at] uw.edu

import sys
from string import maketrans
complement = maketrans("GCATgcat", "CGTAcgta")
if len(sys.argv) != 2:
  print "usage: PROG <mip design file>"
  sys.exit()
mip_design_fh = open(sys.argv[1], 'r')
for line in mip_design_fh:
  if line.startswith('>') or line.startswith('#'):
    continue
  fields = line.rstrip().split()
  if fields[17] == '+':
    print '>' + fields[2] + ':' + fields[3] + '-' + fields[8] + '/' + str(len(fields[6])) + ',' + str(len(fields[10])) + "/+"
    seq = fields[6] + fields[13] + fields[10]
    while len(seq) > 0:
      print seq[:70]
      seq = seq[70:]
  if fields[17] == '-':
    print '>' + fields[2] + ':' + fields[7] + '-' + fields[4] + '/' + str(len(fields[6])) + ',' + str(len(fields[10])) + "/-"
    seq = fields[6] + fields[13] + fields[10]
    reverse_seq = seq[::-1]
    reverse_complement = reverse_seq.translate(complement)
    while len(reverse_complement) > 0:
      print reverse_complement[:70]
      reverse_complement = reverse_complement[70:]
mip_design_fh.close()

