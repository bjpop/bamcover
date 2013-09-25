--------------------------------------------------------------------------------
bamcover - Plot graphs of read coverage for BAM files.
--------------------------------------------------------------------------------

Version: 0.1.0

Authors: Bernard J Pope (1,2) bjpope@unimelb.edu.au

         (1) Victorian Life Sciences Computation Initiative (VLSCI).
         (2) Department of Computing and Information Systems,
             The University of Melbourne.
         

Web:     https://github.com/bjpop/bamcover

License: BSD

Requirements: Python 2.7, PySam, bx Interval Tree, Matplotlib

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
Command line usage:
--------------------------------------------------------------------------------

usage: bamcover [-h] [--version] --coords COORDS [--log FILE] bams [bams ...]

Generate coverage information for BAM files

positional arguments:
  bams             bam files containing mapped reads

optional arguments:
  -h, --help       show this help message and exit
  --version        show program's version number and exit
  --coords COORDS  DNA coordinates to generate coverage for
  --log FILE       Log progress in FILENAME, defaults to stdout.

Explanation of the arguments:


--------------------------------------------------------------------------------
Example usage (should be all on one line)
--------------------------------------------------------------------------------

bamcover --coords coords.tsv --log logfile *.bam


--------------------------------------------------------------------------------
