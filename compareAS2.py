#!/usr/bin/python

# JMG 3/2017

# find reads that are better aligned in one file than another
# version 2: using pysam, counting alignments

import sys
import pysam
if int(pysam.__version__.split('.')[0]) == 0 \
    and int(pysam.__version__.split('.')[1]) < 10:
  sys.stderr.write('Error! pysam version %s (need >= 0.10)\n' \
    % pysam.__version__)
  sys.exit(-1)

def main():
  # load AS's from CL arg
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  <BAM1>  <BAM2>  <out>\n' % sys.argv[0])
    sys.exit(-1)

  # load AS's from SAM
  d = dict()
  f = pysam.AlignmentFile(args[0], 'rb')
  count = 0
  for line in f:
    if line.is_secondary or line.is_supplementary:
      continue
    count += 1
    if line.is_unmapped:
      continue
    d[line.query_name] = (line, line.get_tag('AS'))

  f.close()
  print 'Total reads in %s: %d' % (args[0], count)
  print '  Aligned: %d' % len(d)

  # compare to AS's from piped-in SAM
  count = worse = better = equal = novel = dual = 0
  f = pysam.AlignmentFile(args[1], 'rb')
  for line in f:
    if line.is_secondary or line.is_supplementary:
      continue
    count += 1
    if line.is_unmapped:
      continue
    if line.query_name in d:
      dual += 1
      aScore = line.get_tag('AS')
      if aScore < d[line.query_name][1]:
        better += 1
      elif aScore > d[line.query_name][1]:
        worse += 1
      else:
        equal += 1
    else:
      novel += 1
  f.close()

  #for read in d:
  #  fOut.write(read + '\tnovel\n')
  #  #print 'Novel:', read, d[read]
  #  novel += 1
  #fOut.close()

  print 'Total reads in %s: %d' % (args[1], count), '(should match)'
  print '  Aligned: %d' % (novel + dual)
  print 'Reads aligned in both     :%10d' % dual
  print '  Better in 2nd (elegans) :%10d' % worse
  print '  Better in 1st (briggsae):%10d' % better
  print '  Same AS                 :%10d' % equal

if __name__ == '__main__':
  main()
