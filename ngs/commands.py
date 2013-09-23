import sys, os
import fileinput
import warnings
import collections

# import python-levenshtein
try:
    import Levenshtein
except ImportError:
    sys.exit("""Error: Python Levenshtein is required by 'ngs-tools'
https://pypi.python.org/pypi/python-Levenshtein""")

# import docopt
try:
    from docopt import docopt
except ImportError:
    sys.exit("""Error: docopt is required by 'ngs-tools'.
https://pypi.python.org/pypi/docopt""")

#import SeqIO from Bio (biopython)
try:
    from Bio import SeqIO
except ImportError:
    sys.exit("""Error: Biopython is required by 'ngs-tools'.
https://pypi.python.org/pypi/biopython""")

def merge_fna_qual(options):
    """
usage: ngs-tools merge-fna-qual [options] [--] <fna_file> <qual_file>

options:
    -h --help                       Show this screen.

    -o FILENAME --output=FILENAME   Write FASTQ output to FILENAME.
    """

    options_merge_fna_qual = docopt(merge_fna_qual.__doc__, argv=options)

    try:
        with _open_output_handle(options_merge_fna_qual['--output']) as o:
            try:
                with open(options_merge_fna_qual['<fna_file>']) as f:
                    with open(options_merge_fna_qual['<qual_file>']) as q: 
                        records = SeqIO.QualityIO.PairedFastaQualIterator(f, q)
                        count = SeqIO.write(records, o, "fastq")
            except IOError as e:
                print "Cannot open input file '{0}'. Error: {1}.".format(e.filename, e.strerror)
                raise ValueError("Please provide valid <fna_file> and <qual_file>.")
    except:
        raise

def split_by_barcode(options):
    """
usage: ngs-tools split-by-barcode [--prefix PREFIX | --galaxy GALAXY_ID]
                    [options] [--] <barcode_file> [<input_file>...]

options:
  -h --help                               Show this screen.

  -b BARCODES --barcodes=BARCODES         Comma separate list of barcodes to use
                                          from the  <barcode_file>. By default
                                          all barcodes in the <barcode_file> are
                                          used.
  -q --fastq                              Input file or stdin is in FASTQ
                                          format. By default FASTA format is
                                          expected.
  -k --keep-barcode                       Do not trim the barcode.
  -d INTEGER --max-distance=INTEGER       Max Levenshtein's distance when
                                          looking for mutated barcodes
                                          [default: 3].
  -s INTEGER --barcode-size=INTEGER       Barcode size [default: 11].
  -p PREFIX --prefix=PREFIX               Prefix for splitted output files.
  -g GALAXY_ID --galaxy=GALAXY_ID         Name splitted output files in away
                                          compatible with galaxy multiple
                                          outputs.
  -r REPORT_OUTPUT --report=REPORT_OUTPUT Write report output to REPORT_OUTPUT.
  -o OUTPUT --output=OUTPUT               Write splitted output files to OUTPUT
                                          directory [default: .]
    """

    options_split_by_barcode = docopt(split_by_barcode.__doc__, argv=options)
    m = options_split_by_barcode['<barcode_file>']
    i = options_split_by_barcode['<input_file>']
    b = list(set(options_split_by_barcode['--barcodes'].split(','))) if options_split_by_barcode['--barcodes'] else []
    f = 'fastq' if options_split_by_barcode['--fastq'] else 'fasta'
    k = options_split_by_barcode['--keep-barcode']
    d = int(options_split_by_barcode['--max-distance'])
    s = int(options_split_by_barcode['--barcode-size'])
    p = options_split_by_barcode['--prefix'] if options_split_by_barcode['--prefix'] else ''
    g = options_split_by_barcode['--galaxy'] if options_split_by_barcode['--galaxy'] else ''
    r = options_split_by_barcode['--report']
    o = options_split_by_barcode['--output'].rstrip('/')

    try:
        report = _split_reads(i, m, b, p, g, o, f, d, s, k)
    except:
        raise

    try:
        with _open_output_handle(r) as fh:
            for k, v in report:
                fh.write("{0}\t{1}\n".format(k, v))
    except:
        raise

def seq_convert(options):
    """
usage: ngs-tools seq-convert [options] [--] [<input_file>...]

options:
    -h --help                       Show this screen.

    -o FILENAME --output=FILENAME       Write output to FILENAME.
    -i FORMAT, --input-format FORMAT    Input format [default: fasta]
    -f FORMAT, --output-format FORMAT   Output format [default: fasta]
    """

    options_seq_convert = docopt(seq_convert.__doc__, argv=options)
    inputs = options_seq_convert['<input_file>']
    output = options_seq_convert['--output']
    in_format = options_seq_convert['--input-format']
    out_format = options_seq_convert['--output-format']

    for i in inputs:
        mode = "rb" if in_format in ['abi'] else "r" 
        try:
            with _open_input_handle(i, mode=mode) as in_fh:
                with _open_output_handle(output) as out_fh:
                    seqs = SeqIO.parse(in_fh, in_format)
                    SeqIO.write(seqs, out_fh, out_format)
        except:
            raise

def _open_input_handle(i, mode="r"):
    try:
        if i:
            handle = open(i, mode)
        else:
            handle = sys.stdin
    except IOError as e:
        raise ValueError("Cannot open input file '{0}'. Error: {1}.".format(e.filename, e.strerror))

    return handle

def _open_output_handle(output):
    try:
        if output:
            handle = open(output, "w")
        else:
            handle = sys.stdout
    except IOError as e:
        raise ValueError("Cannot open output file '{0}'. Error: {1}.".format(e.filename, e.strerror))

    return handle

def _split_reads(input_files, barcode_file, barcode_list, prefix, galaxy_id, output, format="fasta", max_distance=3, barcode_size=11, keep_barcode=False):
    """
    Given a fasta/fastq set of reads files and a file with the barcode index. Create one 
    fasta/fastq file for each barcode_name based on the closest matching barcode.

    Adapted from:
    https://gist.github.com/dgrtwo/3725741
    """
    try:
        index = BarcodeIndex(barcode_file, barcode_list, max_distance, barcode_size)
    except:
        raise
 
    counts = collections.defaultdict(int)
 
    try:
        if galaxy_id:
            outfs = dict([(g, open("{0}/primary_{1}_{2}_visible_{3}".format(output, galaxy_id, g, format), "w")) for g in index.barcode_names + ["Unassigned"]])
        else:
            outfs = dict([(g, open("{0}/{1}{2}.ld.{3}.{4}".format(output, prefix, g, max_distance, format), "w")) for g in index.barcode_names + ["Unassigned"]])
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

                report = sorted(counts.items(), key=lambda t: t[0] == "Unassigned")
            finally:
                f.close()
        except IOError as e:
            raise ValueError("Cannot open input file '{0}'. Error: {1}.".format(e.filename, e.strerror))
        finally:            
            for o in outfs.values():
                o.close()
    except IOError as e:
        raise ValueError("Cannot open output file '{0}'. Error: {1}.".format(e.filename, e.strerror))

    return report

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
                    try:
                        barcode_name, barcode = l.split("\t")
                    except ValueError:
                        raise ValueError("Malformed line on barcode index file. Line content: '{0}'".format(l))
                    if len(barcode) != self.barcode_size:
                        warnings.warn("Ignoring barcode '{0}' as is not of correct length.".format(barcode))
                        continue
                    if barcode_list and barcode_name not in barcode_list:
                        continue
                    self.cache[barcode] = barcode_name
            self.barcodes = self.cache.keys()
            self.barcode_names = self.cache.values()
            barcode_missing = list(set(barcode_list) - set(self.barcode_names))
            if barcode_missing:
                print "Some barcodes could not be found in the <barcode_file> file. Barcodes not present: {0}".format(barcode_missing)
                raise ValueError("Please check if this is the correct <barcode_file> or the barcodes are spelled correctly.")
        except IOError as e:
            print "Cannot open file '{0}'. Error: {1}".format(e.filename, e.strerror)
            raise ValueError("Please supply a valid <barcode_file>.")
        except ValueError:
            raise
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

