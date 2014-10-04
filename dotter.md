
 Dotter - Sequence dotplots with image enhancement tools.

 Usage: dotter [options] <horizontal_sequence> <vertical_sequence>  [X options]

   Where <horizontal_sequence> and <vertical_sequence> are file names for FASTA files
   containing the two sequences.

   Allowed sequence types:   Protein  -   Protein
			     DNA      -   DNA
			     DNA      -   Protein

 Options:
  -h, --help
    Show this usage information

  -b <file>, --batch-save=<file>
    Batch mode; save dot matrix to <file>

  -e <file>, --batch-export=<file>
    Batch mode; export plot to PDF file <file>

  -l <file>, --load
    Load dot matrix from <file>

  -m <float>, --memory-limit=<float>
    Memory usage limit in Mb (default 0.5)

  -z <int>, --zoom
    Set zoom (compression) factor

  -p <int>, --pixel-factor
    Set pixel factor manually (ratio pixelvalue/score)

  -W <int>, --window-size
    Set sliding window size. (K => Karlin/Altschul estimate)

  -M <file>, --matrix-file=<file>
    Read in score matrix from <file> (Blast format; Default: Blosum62).

  -F <file>, --sequence-file
    Read in sequences and data from <file> (replaces sequencefiles).

  -f <file>, --feature-file
    Read feature segments from <file>

  -H, --hsp-mode
    Do not calculate dotplot at startup.

  -R, --reverse-greyramp
    Reversed Greyramp tool at start.

  --greyramp-white=<int>
    Set the white point for the greyramp tool (default 40)
    
  --greyramp-black=<int>
    Set the black point for the greyramp tool (default 100)
    If greyramp-black < greyramp-white then colours are reversed

  -r, --reverse-horizontal
    Reverse and complement horizontal_sequence (if DNA)

  -v, --reverse-vertical
    Reverse and complement vertical_sequence (if DNA)

  --suppress-scale
    Do not print horizontal and vertical scales

  -D, --disable-mirror
    Don't display mirror image in self comparisons.  Now
    equivalent to --triangle=u

  --triangle=u|l
    Display dotplot as upper 'u' or lower 'l' triangle.  

  -w, --watson-only
    For DNA: horizontal_sequence top strand only (Watson)

  -c, --crick-only
    For DNA: horizontal_sequence bottom strand only (Crick)

  -q <int>, --horizontal-offset=<int>
    Horizontal_sequence offset

  -s <int>, --vertical-offset=<int>
    Vertical_sequence offset

  --horizontal-type=p|d
    Horizontal_sequence type ('p' for peptide or 'd' for DNA)

  --vertical-type=p|d
    Vertical_sequence type ('p' for peptide or 'd' for DNA)

  --abbrev-title-on
    Abbreviate window title prefixes

  --abbrev-title-off
    Do not abbreviate window title prefixes

  --labels-off
    Do not add labels to the horizontal and vertical axes

  --labels-size=<int>
    Size labels on horizontal and vertical axes to <int>

  --breakline-colour=<colour_str>
    Set the colour used for breaklines between multiple sequences

  --session_colour=<colour_str>
    Set the background colour of the dotter window

  --compiled
    Show package compile date

  --version
    Show package version

-----
 Written by Gemma Barson <gb10@sanger.ac.uk>
 
 Based on original code by Erik Sonnhammer <Erik.Sonnhammer@sbc.su.se>.

 A few additional options added by Douglas Scofield <douglas.scofield@ebc.uu.se>

 Reference: Sonnhammer ELL & Durbin R (1995). A dot-matrix program
 	    with dynamic threshold control suited for genomic DNA and protein
 	    sequence analysis. Gene 167(2):GC1-10.

 See http://www.sanger.ac.uk/resources/software/seqtools/ for more info.

 Copyright (c) 2010-2011: Genome Research Ltd.
 Dotter is distributed under the GNU Public License; see http://www.gnu.org/copyleft/gpl.txt

 Version 4.28  13:33:30 Sep 18 2014

