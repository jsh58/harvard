#!/usr/bin/python

# JMG 3/2017

# find reads that are better aligned in one file than another
# version 2: using pysam, counting alignments
# version 3: editing the BAM
# version 4: removing sec/supp/unmapped alignments, min. MAPQ
# version 5: for PE alignments

import sys
import pysam
if int(pysam.__version__.split('.')[0]) == 0 \
    and int(pysam.__version__.split('.')[1]) < 10:
  sys.stderr.write('Error! pysam version %s (need >= 0.10)\n' \
    % pysam.__version__)
  sys.exit(-1)

def loadAS(f):
  '''
  Load alignment scores from given file.
  '''
  d = dict()
  count = mapped = 0
  for line in f:
    if line.is_secondary or line.is_supplementary:
      continue
    count += 1
    if line.is_unmapped:
      d[line.query_name] = d.get(line.query_name, 0) - 200
      #print line
      #print d[line.query_name]
      #raw_input()
      continue
    mapped += 1
    d[line.query_name] = d.get(line.query_name, 0) + line.get_tag('AS')
    #print line
    #print d[line.query_name]
    #raw_input()
  return d, count, mapped


def main():
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  <BAM1>  <BAM2>  <out>\n' % sys.argv[0])
    sys.exit(-1)

  minMapq = 10   # min. MAPQ to keep an alignment

  # load AS's from 1st BAM
  f = pysam.AlignmentFile(args[0], 'rb')
  d, count, mapped = loadAS(f)
  f.close()
  print 'Total reads in %s: %d' % (args[0], count)
  print '  Aligned: %d' % mapped

  # compare to AS's from 2nd BAM, write output
  count = mapped = dual = worse = better \
    = unmapped = equal = mapqCount = printed \
    = paired = unpaired = proper = 0
  f2 = pysam.AlignmentFile(args[1], 'rb')
  f2.fetch(until_eof=True)
  fOut = pysam.AlignmentFile(args[2], 'wb', f2)
  try:
    line = f2.next()
  except StopIteration:
    sys.stderr.write('Error! No reads in %s\n' % args[1])
    sys.exit(-1)
  flag = True
  while flag:
    head = line.query_name
    mapq = line.mapping_quality  # min. MAPQ for this read

    # sum alignment scores
    aScore = 0
    alns = list()
    while head == line.query_name:
      if line.is_secondary or line.is_supplementary:
        try:
          line = f2.next()
        except StopIteration:
          flag = False
          break
        continue
      count += 1
      if line.is_unmapped:
        aScore -= 200
      else:
        mapped += 1
        alns.append(line)
        if line.mapping_quality < mapq:
          mapq = line.mapping_quality
        aScore += line.get_tag('AS')
      #print line
      #print aScore
      #raw_input()
      try:
        line = f2.next()
      except StopIteration:
        flag = False
        break

    # compare alignments
    if head in d:
      dual += 1

      # no alignments
      if len(alns) == 0:
        if aScore < d[head]:
          better += 1  # better in 1st BAM
        else:
          unmapped += 1
        continue

      # compare AS's
      if aScore < d[head]:
        #line.is_unmapped = True
        #line.reference_id = -1
        #line.reference_start = -1
        #fOut.write(line)
        better += 1  # better in 1st BAM
        continue
      elif aScore > d[head]:
        #fOut.write(line)
        worse += 1   # better in 2nd BAM
      else:
        #fOut.write(line)
        equal += 1
    else:
      pass
      #fOut.write(line)

    if mapq < minMapq:
      mapqCount += 1
    else:
      if len(alns) == 1:
        alns[0].is_paired = False
        unpaired += 1
      elif len(alns) == 2:
        paired += 1
        if alns[0].is_proper_pair:
          proper += 1
      else:
        sys.stderr.write('Error! More than 2 alignments for %s\n' \
          % alns[0].query_name)

      for aln in alns:
        fOut.write(aln)
        printed += 1

  f2.close()
  fOut.close()

  print 'Total reads in %s: %d' % (args[1], count), '(should match)'
  print '  Aligned: %d' % mapped
  print 'PE reads compared         :%10d' % dual
  print '  Better in 2nd (elegans) :%10d' % worse
  print '  Better in 1st (briggsae):%10d' % better
  print '  Same AS                 :%10d' % equal
  print '  Both unmapped           :%10d' % unmapped
  print 'C. elegans alignments     :%10d' % (mapqCount + printed)
  print '  MAPQ < %-2d               :%10d' % (minMapq, mapqCount)
  print '  Printed                 :%10d' % printed
  print 'Paired                    :%10d' % paired
  print '  Properly paired         :%10d' % proper
  print 'Unpaired                  :%10d' % unpaired

if __name__ == '__main__':
  main()
