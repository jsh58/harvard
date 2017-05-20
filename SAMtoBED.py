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
    -a <int>      Print unpaired alignments, with fragment length
                    increased to specified value
    -x            Print unpaired alignments, with fragment length
                    increased to average value calculated from
                    paired alignments
  Other options:
    -v            Run in verbose mode
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

def writeOut(fOut, ref, start, end, read, chr, verbose):
  '''
  Write BED output. Adjust any read that extends beyond
    chromosome ends.
  '''
  if start < 0:
    start = 0
    if verbose:
      sys.stderr.write('Warning! Read %s prevented ' % read \
        + 'from extending below 0 on %s\n' % ref)
  if ref in chr and end > chr[ref]:
    end = chr[ref]
    if verbose:
      sys.stderr.write('Warning! Read %s prevented ' % read \
        + 'from extending past %d on %s\n' % (chr[ref], ref))
  fOut.write('%s\t%d\t%d\t%s\n' % (ref, start, end, read))

def checkPaired(pos, verbose):
  '''
  Check if any paired alignments weren't processed.
  '''
  unpaired = 0
  for r in pos:
    if pos[r] >= 0:
      if verbose:
        sys.stderr.write('Warning! Read %s missing its pair\n' % r)
      unpaired += 1
  return unpaired

def processPaired(header, chrom, rc, start, offset, pos,
    chr, fOut, extendOpt, length, verbose):
  '''
  Process a properly paired SAM record. If first, save end
    position to pos dict; if second, write complete record.
  '''
  # 2nd of PE reads
  if header in pos:
    if pos[header] < 0:
      sys.stderr.write('Error! Read %s already analyzed\n' % header)
      sys.exit(-1)

    # save end position
    if rc:
      start += offset

    writeOut(fOut, chrom, min(start, pos[header]), \
      max(start, pos[header]), header, chr, verbose)

    # keep track of fragment lengths
    length[0] += 1
    length[1] += abs(start - pos[header])

    pos[header] = -1  # records that read was processed

  # 1st of PE reads: save end position
  else:
    if rc:
      pos[header] = start + offset
    else:
      pos[header] = start

def processUnpaired(header, chrom, rc, start, offset, pos,
    chr, fOut, addBP, verbose):
  '''
  Process an unpaired SAM record.
  '''
  # extend 3' end of read (to total length 'addBP')
  end = start + offset
  if addBP != 0:
    if rc:
      start = min(end - addBP, start)
    else:
      end = max(start + addBP, end)

  writeOut(fOut, chrom, start, end, header, chr, verbose)
  pos[header] = -2  # records that read was processed

def processSingle(single, pos, chr, fOut, length, verbose):
  '''
  Process saved singletons (unpaired alignments)
    using calculated extension size.
  '''
  if length[0] == 0:
    sys.stderr.write('Error! Cannot calculate fragment ' \
      + 'lengths: no paired alignments\n')
    sys.exit(-1)

  # calculate average paired alignment length
  avg = int(round(length[1] / length[0]))

  # process reads
  count = 0
  for header in single:
    for idx in range(len(single[header])):
      chrom, rc, start, offset = single[header][idx]
      processUnpaired(header, chrom, rc, start, offset,
        pos, chr, fOut, avg, verbose)
      count += 1
  return count

def parseSAM(fIn, fOut, singleOpt, addBP, extendOpt, verbose):
  '''
  Parse the input file, and produce the output file.
  '''
  chr = {}  # chromosome lengths
  pos = {}  # position of first alignment (for paired reads)
  single = {}  # to save unpaired alignments (for calc.-extension option)
  count = 0    # count of unpaired alignments
  length = [0, 0.0]  # for calculating fragment lengths
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
      sys.stderr.write('Error! Poorly formatted SAM record:\n' \
        + line)
      sys.exit(-1)
    flag = getInt(spl[1])
    start = getInt(spl[3]) - 1

    # skip unmapped, secondary, and supplementary
    if flag & 0x904:
      line = fIn.readline().rstrip()
      continue

    # process alignment
    offset = parseCigar(spl[5]) + len(spl[9])
    if flag & 0x2:
      # properly paired alignment
      processPaired(spl[0], spl[2], flag & 0x10, start, offset,
        pos, chr, fOut, extendOpt, length, verbose)

    elif extendOpt:
      # with calculated-extension option, save unpaired alignments
      #   until after extension length is calculated
      if spl[0] in single:
        single[spl[0]].append((spl[2], flag & 0x10, start, offset))
      else:
        single[spl[0]] = [(spl[2], flag & 0x10, start, offset)]
      pos[spl[0]] = -2  # records that read was processed

    elif singleOpt:
      # process singletons directly (w/o extendOpt)
      processUnpaired(spl[0], spl[2], flag & 0x10, start,
        offset, pos, chr, fOut, addBP, verbose)
      count += 1

    line = fIn.readline().rstrip()

  # check for paired alignments that weren't processed
  unpaired = checkPaired(pos, verbose)

  # for calculated-extension option, process saved unpaired alns
  if extendOpt:
    count = processSingle(single, pos, chr, fOut, length, verbose)

  # log counts
  if verbose:
    sys.stderr.write('Paired alignments (fragments): ' \
      + '%d (%d)\n' % (length[0]*2, length[0]))
    if length[0]:
      sys.stderr.write('  Average fragment length: %.1fbp\n' \
        % ( length[1] / length[0] ))
    if unpaired:
      sys.stderr.write('"Paired" reads missing mates: %d\n' % unpaired)
    if singleOpt or extendOpt:
      sys.stderr.write('Unpaired alignments: %d\n' % count)
      if extendOpt:
        addBP = int(round(length[1] / length[0]))
      if addBP:
        sys.stderr.write('  (extended to length %dbp)\n' % addBP)

def main():
  '''
  Main.
  '''
  # Default parameters
  infile = None      # input file
  outfile = None     # output file
  singleOpt = False  # option to print unpaired alignments
  addBP = 0          # number of bp to add to unpaired reads
  extendOpt = False  # option to calculate extension size
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
    elif args[i] == '-x':
      extendOpt = True
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

  # process files
  parseSAM(infile, outfile, singleOpt, addBP, extendOpt, verbose)
  if infile != sys.stdin:
    infile.close()
  if outfile != sys.stdout:
    outfile.close()

if __name__ == '__main__':
  main()
