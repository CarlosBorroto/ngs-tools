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
    sample              Randomly subsample records.
    coverage            Coverage report.
"""

import sys

try:
    from docopt import docopt
except ImportError:
    sys.stderr.write("ERROR: 'docopt' is required by ngs-tools. https://pypi.python.org/pypi/docopt\n")
    sys.exit(1)
from commands import merge_fna_qual, split_by_barcode, seq_convert, sample, coverage

version = __import__('ngstools.version').get_version()


def main():
    try:
        options = docopt(__doc__, version=version, options_first=True)

        if options['<command>'] == 'merge-fna-qual':
            merge_fna_qual([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'split-by-barcode':
            split_by_barcode([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'seq-convert':
            seq_convert([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'sample':
            sample([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'coverage':
            coverage([options['<command>']] + options['<args>'])
        else:
            exit("%r is not a ngs-tools command. See 'ngs-tools --help'." % options['<command>'])
        return 0
    except KeyboardInterrupt:
        return 1
#    except Exception, e:
#        sys.stderr.write("ERROR: {0}.\n".format(str(e)))
        return 1


if __name__ == '__main__':
    sys.exit(main())
