#!/usr/bin/python

# JMG 3/2017

# find reads that are better aligned in one file than another
# version 2: using pysam, counting alignments
# version 3: editing the BAM
# version 4: removing sec/supp/unmapped alignments, min. MAPQ

import sys
import pysam
if int(pysam.__version__.split('.')[0]) == 0 \
    and int(pysam.__version__.split('.')[1]) < 10:
  sys.stderr.write('Error! pysam version %s (need >= 0.10)\n' \
    % pysam.__version__)
  sys.exit(-1)


def main():
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  <BAM1>  <BAM2>  <out>\n' % sys.argv[0])
    sys.exit(-1)

  minMapq = 10   # min. MAPQ to keep an alignment

  # load AS's from 1st BAM
  d = dict()
  f = pysam.AlignmentFile(args[0], 'rb')
  count = 0
  for line in f:
    if line.is_secondary or line.is_supplementary:
      continue
    count += 1
    if line.is_unmapped:
      continue
    d[line.query_name] = line.get_tag('AS')
  f.close()
  print 'Total reads in %s: %d' % (args[0], count)
  print '  Aligned: %d' % len(d)

  # compare to AS's from 2nd BAM, write output
  count = worse = better = equal = novel = dual = mapq = 0
  f = pysam.AlignmentFile(args[1], 'rb')
  fOut = pysam.AlignmentFile(args[2], 'wb', f)
  for line in f:
    if line.is_secondary or line.is_supplementary:
      #fOut.write(line)
      continue
    count += 1
    if line.is_unmapped:
      #fOut.write(line)
      continue
    if line.query_name in d:
      dual += 1
      aScore = line.get_tag('AS')
      if aScore < d[line.query_name]:
        #line.is_unmapped = True
        #line.reference_id = -1
        #line.reference_start = -1
        #fOut.write(line)
        better += 1
        continue
      elif aScore > d[line.query_name]:
        #fOut.write(line)
        worse += 1
      else:
        #fOut.write(line)
        equal += 1
    else:
      #fOut.write(line)
      novel += 1

    if line.mapping_quality < minMapq:
      mapq += 1
    else:
      fOut.write(line)

  f.close()

  print 'Total reads in %s: %d' % (args[1], count), '(should match)'
  print '  Aligned: %d' % (novel + dual)
  print 'Reads aligned in both     :%10d' % dual
  print '  Better in 2nd (elegans) :%10d' % worse
  print '  Better in 1st (briggsae):%10d' % better
  print '  Same AS                 :%10d' % equal
  print 'Reads with MAPQ < %d: %d' % (minMapq, mapq)

if __name__ == '__main__':
  main()
