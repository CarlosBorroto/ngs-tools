import sys
import fileinput
import warnings
import collections

# import python-levenshtein
try:
    import Levenshtein
except ImportError:
    sys.exit("""Error: Python Levenshtein is required by '454-tools'
https://pypi.python.org/pypi/python-Levenshtein""")

# import docopt
try:
    from docopt import docopt
except ImportError:
    sys.exit("""Error: docopt is required by '454-tools'.
https://pypi.python.org/pypi/docopt""")

#import SeqIO from Bio (biopython)
try:
    from Bio import SeqIO
except ImportError:
    sys.exit("""Error: Biopython is required by '454-tools'.
https://pypi.python.org/pypi/biopython""")

def merge(options):
    """
usage: 454-tools merge [options] [--] <fna_file> <qual_file>

options:
    -h --help                       Show this screen.

    -o FILENAME --output=FILENAME   Write FASTQ output to FILENAME.
    """

    options_merge = docopt(merge.__doc__, argv=options)

    try:
        with _open_output_handle(options_merge['--output']) as o:
            try:
                with open(options_merge['<fna_file>']) as f:
                    with open(options_merge['<qual_file>']) as q: 
                        records = SeqIO.QualityIO.PairedFastaQualIterator(f, q)
                        count = SeqIO.write(records, o, "fastq")
            except IOError as e:
                print "Cannot open input file '{0}'. Error: {1}.".format(e.filename, e.strerror)
                print "Please provide valid <fna_file> and <qual_file>."
    except ValueError as e:
        print e

def mid_split(options):
    """
usage: 454-tools mid-split [--barcode BARCODE]... [options] [--] <mid_file> [<input_file>...]

options:
    -h --help                           Show this screen.

    -b BARCODE --barcode=BARCODE        Barcode to use from the <mid_file>. Allow multiple appearances.
                                        By default all barcodes in the <mid_file> are used.
    -q --fastq                          Input file or stdin is in FASTQ format. By default FASTA format
                                        is expected.
    -k --keep-barcode                   Do not trim the barcode.
    -d INTEGER --max-distance=INTEGER   Max Levenshtein's distance when looking for mutated barcodes
                                        [default: 3].
    -s INTEGER --barcode-size=INTEGER   Barcode size [default: 11].
    """

    options_mid_split = docopt(mid_split.__doc__, argv=options)
    m = options_mid_split['<mid_file>']
    i = options_mid_split['<input_file>']
    b = options_mid_split['--barcode']
    f = 'fastq' if options_mid_split['--fastq'] else 'fasta'
    k = options_mid_split['--keep-barcode']
    d = int(options_mid_split['--max-distance'])
    s = int(options_mid_split['--barcode-size'])

    try:
        _split_reads(i, m, b, f, d, s, k)
    except ValueError as e:
        print e

def _open_output_handle(output):
    try:
        if output:
            handle = open(output, "w")
        else:
            handle = sys.stdout
    except IOError as e:
        print "Cannot open output file '{0}'. Error: {1}.".format(e.filename, e.strerror)
        raise ValueError("Please make sure you can write to output file: '{0}'.".format(e.filename))

    return handle

def _split_reads(input_files, barcode_file, barcode_list, format="fasta", max_distance=3, barcode_size=11, keep_barcode=False):
    """
    Given a fasta/fastq set of reads files and a file with the barcode index. Create one 
    fasta/fastq file for each barcode_name based on the closest matching barcode.

    Adapted from:
    https://gist.github.com/dgrtwo/3725741
    """
    try:
        index = BarcodeIndex(barcode_file, barcode_list, max_distance, barcode_size)
    except ValueError:
        raise
 
    counts = collections.defaultdict(int)
 
    try:
        outfs = dict([(g, open("{0}.ld.{1}.{2}".format(g, max_distance, format), "w")) for g in index.barcode_names + ["Unassigned"]])
        try:
            f = fileinput.input(input_files)
            try:
                for r in SeqIO.parse(f, format):
                    barcode_name = index.find_barcode(str(r.seq[:barcode_size]))
                    if barcode_name != None:
                        barcode_name = barcode_name
                        if not keep_barcode:
                            r = r[barcode_size:] 
                    else:
                        barcode_name = "Unassigned"
                    SeqIO.write([r], outfs[barcode_name], format)
                    counts[barcode_name] += 1

                for k, v in sorted(counts.items(), key=lambda t: t[0] == "Unassigned"):
                    print k, v
            finally:
                f.close()
        except IOError as e:
            print "Cannot open input file '{0}'. Error: {1}.".format(e.filename, e.strerror)
            raise ValueError("Please supply valid fasta or fastq <input_file>.")
        finally:            
            for o in outfs.values():
                o.close()
    except IOError as e:
        print "Cannot open output file '{0}'. Error: {1}.".format(e.filename, e.strerror)
        raise ValueError("Please make sure you can write to output file: '{0}'.".format(e.filename))

class BarcodeIndex:
    """
    Represents a set of indices, with the ability to find the closest one
    
    Adapted from:
    https://gist.github.com/dgrtwo/3725741
    """
    def __init__(self, barcode_file, barcode_list, max_distance=3, barcode_size=11):
        self.cache = {}
        self.barcode_size = barcode_size
        try:
            with open(barcode_file) as inf:
                for l in inf:
                    l = l.rstrip('\n')
                    if l.startswith("#") or l == "":
                        continue
                    barcode_name, barcode = l.split("\t")
                    if len(barcode) != self.barcode_size:
                        warnings.warn("Ignoring barcode '{0}' as is not of correct length.".format(barcode))
                        continue
                    if barcode_list and barcode_name not in barcode_list:
                        continue
                    self.cache[barcode] = barcode_name
            self.barcodes = self.cache.keys()
            self.barcode_names = self.cache.values()
        except IOError as e:
            print "Cannot open file '{0}'. Error: {1}".format(e.filename, e.strerror)
            raise ValueError("Please supply a valid <mid_file>.")
        self.max_distance = max_distance
 
    def find_barcode(self, barcode):
        """
        match barcode and return the barcode_name it's supposed to be in. If
        there is none within the Levenshtein distance given, or if there
        is a tie, return None
        """
        # if there's an exact match, return that
        exact = self.cache.get(barcode)
        if exact is not None:
            return exact
 
        # find the Levenshtein distance to each
        distances = [Levenshtein.distance(barcode, b) for b in self.barcodes]
        best = min(distances)
        # check if there's a tie or the distance is too great:
        if best > self.max_distance or distances.count(best) > 1:
            return None
        # otherwise, return the best one, after caching it for future use
        ret = self.barcode_names[distances.index(best)]
        self.cache[barcode] = ret
        return ret

if __name__ == '__main__':
    exit()

