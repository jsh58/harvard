#!/usr/bin/python

# JMG 1/2018

# Producing an html summary of the top 20 taxa from
#   centrifuge's kraken-style report.

import sys
import gzip

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
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
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
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

class Node:
  '''
  Node: contains taxon name, taxon number,
    list of child nodes, parent node,
    score (percent of the sample), and
    count (number of reads).
  '''
  def __init__(self, parent, name, taxon, score, count):
    self.child = []
    self.parent = parent
    self.name = name
    self.taxon = taxon
    self.score = float(score)
    self.count = int(count)

def printFooter(f, num, version, date):
  '''
  Print centrifuge/nt info, etc. for the output.
  '''
  f.write('<p>Top %d taxa, identified by ' % num \
    + '<a href="http://www.ccb.jhu.edu/software/centrifuge/manual.shtml">' \
    + 'Centrifuge</a>')
  if version:
    f.write(' (version %s)' % version)
  f.write(', querying the NCBI ' \
    + '<a href="https://www.ncbi.nlm.nih.gov/nucleotide">nt</a> database')
  if date:
    f.write(' (downloaded %s)' % date)
  f.write('.</p>\n')
  f.write('''<p>The value for a given taxon is the percent of all the
  sequence reads assigned to that taxon <strong>or</strong> any lower
  node in its tree.</p>
<h4>Caveats:</h4>
<ul>
  <li>Some sequences in the nt database are contaminated. For example,
    reads with untrimmed Illumina adapters may be assigned to sundry taxa
    (e.g. <strong><i>Cyprinus carpio</i></strong> [in class Actinopteri],
    <strong><i>Eimeria mitis</i></strong> [in phylum Apicomplexa],
    <strong><i>Ralstonia solanacearum</i></strong> [in class
    Betaproteobacteria]) because of the contaminating adapter
    sequences.</li>
  <p>
  <li>Reads derived from one organism may align equally well to other
    related organisms, especially those well-represented in the nt database.
    For example, reads from a sequencing run of a <i>Homo sapiens</i> sample
    may align to <i>Pan</i> or even <i>Mus</i>. Reads that map to multiple
    taxa are counted as a fractional portion to each taxon, rather than to
    the lowest common ancestor (LCA).</li>
  <p>
  <li>The "unclassified" category includes both reads that did not match
    anything in the nt database, plus those that matched a sequence with
    an unspecified or unknown taxonomy.</li>
  <p>
  <li>The depiction of the results above is based on the major levels
    (DKPCOFGS) in the NCBI's
    <a href="https://www.ncbi.nlm.nih.gov/taxonomy">taxonomy</a>
    tree. Some branches in that tree skip a major level; hence, the
    columns in the above table should not be interpreted as
    corresponding to a specific taxonomic level.</li>
</ul>
<p>Questions/concerns/comments/suggestions?
<a href="mailto:jgaspar@fas.harvard.edu">Please let us know.</a></p>
''')

def printLevel(f, n, level, cutoff):
  '''
  Print results for a node (if its count meets cutoff).
    Continue printing for children nodes (recursively).
  '''
  if n.count >= cutoff:  # or level == 0: # to include all children of root
    f.write('  <tr>\n' \
      + '    <td align="right">%.2f&emsp;</td>\n' % n.score \
      + '    <td>%s%s</td>\n' % (level * 2 * '&emsp;', n.name) \
      + '  </tr>\n')
  for m in n.child:
    printLevel(f, m, level + 1, cutoff)

def printOutput(f, unclass, root, num, cutoff, version, date):
  '''
  Begin printing results (header and unclassified).
    Start recursive tree printing.
  '''
  f.write('''<h2>Taxonomy Analysis</h2>
<strong><font color="red" size="4">Warning:</font></strong>
<font size="4"> experimental software; not suitable for publication</font>
<p>
<table style="width:100%;border:1px solid;">
  <tr>
    <th align="right" width=10%>Percent&emsp;</th>
    <th align="left">Taxon</th>
  </tr>
''')
  f.write('  <tr>\n' \
    + '    <td align="right">%.2f&emsp;</td>\n' % unclass \
    + '    <td>unclassified</td>\n' \
    + '  </tr>\n')

  # print tree
  for n in root.child:
    printLevel(f, n, 0, cutoff)
  f.write('</table>\n')

  printFooter(f, num, version, date)

def findCutoff(score, x):
  '''
  Determine threshold for top x scores.
  '''
  if x >= len(score):
    return 0
  res = sorted(score, reverse=True)
  return res[x-1]

def checkTree(node, taxon):
  '''
  Search node and its children for given taxon.
  '''
  if node.taxon == taxon:
    return node
  for n in node.child:
    parent = checkTree(n, taxon)
    if parent:
      return parent
  return None

def loadScores(f, d):
  '''
  Create taxonomic tree (including scores) from a
    Centrifuge report.
  '''
  rank = 'DKPCOFGS'
  unclass = 0.0  # 'unclassified' score
  root = Node(None, 'root', '1', -1, -1) # root of tree
  temp = root    # pointer to previous node
  score = []     # list of scores (read counts)

  for line in f:
    spl = line.split('\t')
    if len(spl) < 6:
      sys.stderr.write('Error! Improperly formatted ' \
        + 'centrifuge-kreport file\n')
      sys.exit(-1)

    # skip non-canonical levels
    if spl[3] == '-':
      if spl[4] not in ['12908', '28384']:
        # make exception for these '-' taxa
        #   (see loadTax(), below)
        continue
    elif spl[3] == 'U':
      # save unclassified value automatically
      unclass = float(spl[0])
      continue
    elif spl[3] not in rank:
      sys.stderr.write('Warning! Unknown taxonomic rank:' \
        + ' %s\n' % spl[3] + '  ' + line)
      continue

    # find parent node in hierarchy
    if spl[4] not in d:
      sys.stderr.write('Warning! Unknown taxon: %s\n' % spl[4])
      continue
    parent = None
    # check current branch first
    while temp != None:
      if temp.taxon == d[spl[4]]:
        parent = temp
        break
      temp = temp.parent
    # if not found, check whole tree from root
    if parent == None:
      parent = checkTree(root, d[spl[4]])
    if parent == None:
      sys.stderr.write('Warning! Cannot find parent for ' \
        + 'taxon %s\n' % spl[4])
      continue

    # save to tree
    name = spl[5].strip()
    if spl[3] in 'GS':
      name = '<i>' + name + '</i>'  # italicize genus/species
    n = Node(parent, name, spl[4], spl[0], spl[1])
    parent.child.append(n)
    temp = n

    # save score
    score.append(int(spl[1]))

  return unclass, root, score

def findParent(d, taxon):
  '''
  Find parent taxon (recursively).
  '''
  if taxon not in d or d[taxon][0] not in d:
    return None
  if d[ d[taxon][0] ][1]:
    return d[taxon][0]
  return findParent(d, d[taxon][0])

def loadTax(f):
  '''
  Load parents of each taxon, keeping only
    canonical (DKPCOFGS) ones.
    (include 1 -> 'root'
         12908 -> 'unclassified sequences'
         28384 -> 'other sequences')
  '''
  # load immediate parent taxa to dict
  temp = {}
  for line in f:
    spl = line.split('|')
    if len(spl) < 3:
      sys.stderr.write('Error! Improperly formatted tree file\n')
      sys.exit(-1)

    level = False
    if spl[2].strip() in ['superkingdom', 'kingdom', 'phylum', \
        'class', 'order', 'family', 'genus', 'species'] \
        or spl[0].strip() in ['1', '12908', '28384']:
      level = True
    temp[spl[0].strip()] = (spl[1].strip(), level)

  # save canonical parent taxa to dict
  d = {}
  for taxon in temp:
    if temp[taxon][1]:
      parent = findParent(temp, taxon)
      if parent:
        d[taxon] = parent
      else:
        sys.stderr.write('Warning! Cannot find parent ' \
          + 'of taxon %s\n' % taxon)

  return d

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: python %s  ' % sys.argv[0] \
      + '<kreport>  <taxTree>  <out> \ \n' \
      + '    [<num>]  [<version>  <date>]\n' \
      + '  <num>      Number of taxa to print (def. 20)\n' \
      + '  <version>  Version of centrifuge\n' \
      + '  <date>     Date of nt download\n')
    sys.exit(-1)

  # load tax tree
  fTax = openRead(args[1])
  d = loadTax(fTax)
  if fTax != sys.stdin:
    fTax.close()

  # load scores and create taxonomic tree
  fIn = openRead(args[0])
  unclass, root, score = loadScores(fIn, d)
  if fIn != sys.stdin:
    fIn.close()

  # find cutoff score for top n taxa
  num = 20
  if len(args) > 3:
    num = int(args[3])
  cutoff = findCutoff(score, num)

  # load centrifuge version, date of nt download
  version = date = ''
  if len(args) > 5:
    version = args[4]
    date = args[5]

  # print output
  fOut = openWrite(args[2])
  printOutput(fOut, unclass, root, num, cutoff, version, date)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
