"""
usage: 454-tools <command> [<args>...]
          [--version] [--help]

options:
    -h --help   Show this screen.
    --version   Show version.

Available commands are:
    merge       Merge 454's Fasta and Quality files into a Fastq file.
    mid-split   Split Fasta or Fastq files using MID barcodes.
"""

import sys
try:
    from docopt import docopt
except ImportError:
    sys.exit("""Error: docopt is required by 454-tools.
https://pypi.python.org/pypi/docopt""")
from tools import merge, mid_split

def main():
    try:
        options = docopt(__doc__, version='0.1.2', options_first=True)

        if options['<command>'] == 'merge':
            merge([options['<command>']] + options['<args>'])
        elif options['<command>'] == 'mid-split':
            mid_split([options['<command>']] + options['<args>'])
        else:
            exit("%r is not a 454-tools command. See '454-tools --help'." % options['<command>'])
    except KeyboardInterrupt:
        pass

if __name__ == '__main__':
        main()

