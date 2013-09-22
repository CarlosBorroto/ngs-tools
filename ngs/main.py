"""
usage: ngs-tools <command> [<args>...]
          [--version] [--help]

options:
    -h --help   Show this screen.
    --version   Show version.

Available commands are:
    merge-fna-qual      Merge 454's Fasta and Quality files into a Fastq file.
    split-by-barcode    Split Fasta or Fastq files using barcodes.
    seq-convert         Convert one sequence format to another.
"""

import sys
try:
    from docopt import docopt
except ImportError:
    sys.exit("""Error: docopt is required by ngs-tools.
https://pypi.python.org/pypi/docopt""")
from commands import merge_fna_qual, split_by_barcode, seq_convert

def main():
    try:
        options = docopt(__doc__, version='0.1.6', options_first=True)

        if options['<command>'] == 'merge-fna-qual':
            merge_fna_qual([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'split-by-barcode':
            split_by_barcode([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'seq-convert':
            seq_convert([options['<command>']] + options['<args>'])
        else:
            exit("%r is not a ngs-tools command. See 'ngs-tools --help'." % options['<command>'])
    except KeyboardInterrupt:
        pass
#    except ValueError as e:
#        print e

if __name__ == '__main__':
        main()

