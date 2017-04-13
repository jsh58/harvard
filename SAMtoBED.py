#!/usr/bin/python

# JMG 4/13/16

# Produce a BED file from a SAM.
# Combine paired-end alignments.

import sys
import re

def usage():
  sys.stderr.write('''Usage: python SAMtoBED.py  [options]  -i <input>  -o <output>
    -i <input>    SAM alignment file (can be in any sort order,
                    or unsorted; use '-' for stdin)
    -o <output>   Output BED file
  Options for unpaired alignments:
    -n            Do not print unpaired alignments (default)
    -y            Print unpaired alignments
    -a <int>      Print unpaired alignments, with read length
                    increased to specified value
''')
  sys.exit(-1)

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot read input file %s\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot write to output file %s\n' % filename)
    sys.exit(-1)
  return f

def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to int\n' % arg)
    sys.exit(-1)
  return val

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

def processUnpaired(spl, flag, start, offset, pos, chr, fOut, addBP):
  '''
  Process an unpaired SAM record.
  '''
  # adjust ends (using parameter 'addBP')
  end = start + offset + len(spl[9])
  if addBP != 0:
    if flag & 0x10:
      start = min(end - addBP, start)
    else:
      end = max(start + addBP, end)

  writeOut(fOut, spl[2], start, end, spl[0], chr)
  pos[spl[0]] = -2  # records that read was processed

def parseSAM(fIn, fOut, singleOpt, addBP):
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
      sys.stderr.write('Error! Poorly formatted SAM record:\n' + line)
      sys.exit(-1)
    flag = getInt(spl[1])
    start = getInt(spl[3]) - 1

    # skip unmapped, secondary, and supplementary
    if flag & 0x904:
      line = fIn.readline().rstrip()
      continue

    # process alignment
    offset = parseCigar(spl[5])
    if flag & 0x2:
      processPaired(spl, flag, start, offset, pos, chr, fOut)
    elif singleOpt:
      processUnpaired(spl, flag, start, offset, pos, chr, fOut, addBP)

    line = fIn.readline().rstrip()

  # check paired alignments that weren't processed
  unpaired = checkPaired(pos)
  if unpaired:
    sys.stderr.write('Reads lacking pairs: %d\n' % unpaired)

def main():
  '''
  Main.
  '''
  # Default parameters
  infile = None      # input file
  outfile = None     # output file
  singleOpt = False  # option to print unpaired alignments
  addBP = 0          # number of bp to add to unpaired reads
  verbose = False    # verbose option

  # get command-line args
  args = sys.argv[1:]
  i = 0
  while i < len(args):
    if args[i] == '-h' or args[i] == '--help':
      usage()
    elif args[i] == '--version':
      pass  #printVersion()
    elif args[i] == '-v':
      verbose = True
    elif args[i] == '-n':
      singleOpt = False
    elif args[i] == '-y':
      singleOpt = True
    elif i < len(args) - 1:
      if args[i] == '-i':
        infile = openRead(args[i+1])
      elif args[i] == '-o':
        outfile = openWrite(args[i+1])
      elif args[i] == '-a':
        singleOpt = True
        addBP = max(getInt(args[i+1]), 0)
      else:
        sys.stderr.write('Error! Unknown parameter: %s\n' % args[i])
        usage()
      i += 1
    else:
      sys.stderr.write('Error! Unknown parameter with no arg: ' \
        + '%s\n' % args[i])
      usage()
    i += 1

  # check for I/O errors
  if infile == None or outfile == None:
    sys.stderr.write('Error! Must specify input and output files\n')
    usage()

  parseSAM(infile, outfile, singleOpt, addBP)
  if infile != sys.stdin:
    infile.close()
  if outfile != sys.stdout:
    outfile.close()

if __name__ == '__main__':
  main()
