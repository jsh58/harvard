#!/usr/bin/python

# JMG 4/13/16

# Produce a BED file from a SAM.
# Combine paired-end alignments.

import sys
import re

def usage():
  print "Usage: python SAMtoBED.py  [options]  <input>  <output>  \n\
    <input>     SAM file, or '-' for stdin                        \n\
    <output>    BED file                                          \n\
  Options for unpaired alignments:                                \n\
    -n        Do not print unpaired alignments (default)          \n\
    -y        Print unpaired alignments                           \n\
    -a <int>  Print unpaired alignments, with read length         \n\
                increased to specified value                       "
  sys.exit(-1)

def loadChrLen(line, chr):
  '''
  Load chromosome lengths from the SAM file header.
  '''
  mat = re.search(r'@SQ\s+SN:(\S+)\s+LN:(\d+)', line)
  if mat:
    chr[mat.group(1)] = int(mat.group(2))

def parseCigar(cigar):
  '''
  Determine indel offset of alignment from CIGAR.
  '''
  offset = 0
  ins = re.findall(r'(\d+)I', cigar)
  for i in ins:
    offset -= int(i)
  de = re.findall(r'(\d+)D', cigar)
  for d in de:
    offset += int(d)
  return offset

def writeOut(fOut, ref, start, end, read, chr):
  '''
  Write BED output. Adjust any read that extends beyond
    chromosome ends.
  '''
  if start < 0:
    start = 0
    sys.stderr.write('Warning! Read %s prevented from ' % read + \
      'extending below 0 on %s\n' % ref)
  if ref in chr and end > chr[ref]:
    end = chr[ref]
    sys.stderr.write('Warning! Read %s prevented from ' % read + \
      'extending past %d on %s\n' % (chr[ref], ref))
  fOut.write('%s\t%d\t%d\t%s\n' % \
    (ref, start, end, read))

def checkPaired(pos):
  '''
  Check if paired alignments weren't processed.
  '''
  unpaired = 0
  for r in pos:
    if pos[r] >= 0:
      sys.stderr.write('Error! Read %s missing its pair\n' % r)
      unpaired += 1
  return unpaired

def processPaired(spl, flag, start, offset, pos, chr, fOut):
  '''
  Process a properly paired SAM record.
  '''
  # 2nd of PE reads
  if spl[0] in pos:
    if pos[spl[0]] < 0:
      sys.stderr.write('Error! Read %s already analyzed\n' % spl[0])
      sys.exit(-1)

    # save end position
    if flag & 0x10:
      start += offset + len(spl[9])

    writeOut(fOut, spl[2], min(start, pos[spl[0]]), \
      max(start, pos[spl[0]]), spl[0], chr)

    pos[spl[0]] = -1  # records that read was processed

  # 1st of PE reads: save end position
  else:
    if flag & 0x10:
      pos[spl[0]] = start + offset + len(spl[9])
    else:
      pos[spl[0]] = start

def processUnpaired(spl, flag, start, offset, pos, chr, fOut, add):
  '''
  Process an unpaired SAM record.
  '''
  # adjust ends (using parameter 'add')
  end = start + offset + len(spl[9])
  if add != 0:
    if flag & 0x10:
      start = min(end - add, start)
    else:
      end = max(start + add, end)

  writeOut(fOut, spl[2], start, end, spl[0], chr)
  pos[spl[0]] = -2  # records that read was processed

def parseSAM(fIn, fOut, add):
  '''
  Parse the input file, and produce the output file.
  '''
  chr = {}  # chromosome lengths
  pos = {}  # position of alignment (for paired reads)
  line = fIn.readline().rstrip()
  while line:

    # skip header
    if line[0] == '@':
      loadChrLen(line, chr)  # load chromosome length
      line = fIn.readline().rstrip()
      continue

    # save flag and start position
    spl = line.split('\t')
    if len(spl) < 11:
      print 'Error! Poorly formatted SAM record:\n', line
      sys.exit(-1)
    try:
      flag = int(spl[1])
      start = int(spl[3]) - 1
    except ValueError:
      print 'Error parsing SAM record:\n', line
      sys.exit(-1)

    # skip unmapped, secondary, and supplementary
    if flag & 0x904:
      line = fIn.readline().rstrip()
      continue

    # process alignment
    offset = parseCigar(spl[5])
    if flag & 0x2:
      processPaired(spl, flag, start, offset, pos, chr, fOut)
    elif add != -1:
      processUnpaired(spl, flag, start, offset, pos, chr, fOut, add)

    line = fIn.readline().rstrip()

  # check paired alignments that weren't processed
  unpaired = checkPaired(pos)
  if unpaired:
    sys.stderr.write('Reads lacking pairs: %d\n' % unpaired)

def main():
  '''
  Runs the program.
  '''
  # parse command-line args
  args = sys.argv[1:]
  if not args: usage()
  add = -1  # number of bp to add to unpaired reads
  i = 0
  if args[i] == '-y':
    add = 0
    i += 1
  elif args[i] == '-a':
    try:
      add = int(args[i+1])
    except ValueError:
      print 'Error! Cannot convert %s to int' % args[i+1]
      sys.exit(-1)
    i += 2
  elif args[i] == '-n':
    i += 1

  # open files
  if len(args[i:]) < 2: usage()
  try:
    fIn = open(args[i], 'rU')
  except IOError:
    if args[i] == '-':
      fIn = sys.stdin
    else:
      print 'Error! Cannot open', args[i]
      usage()
  fOut = open(args[i+1], 'w')

  parseSAM(fIn, fOut, add)
  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()
