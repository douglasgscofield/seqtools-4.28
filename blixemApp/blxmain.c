/*  File: blxmain.c
 *  Author: Erik Sonnhammer, 1999-08-26
 *  Copyright (c) 2009 - 2012 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: main() function for Blixem
 *              BLast matches In an X-windows Embedded Multiple alignment
 *----------------------------------------------------------------------------
 */

#include <blixemApp/blixem_.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxparser.h>
#include <seqtoolsUtils/blxGff3Parser.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>


/* Some globals.... */


/* global debug flag for blixem, set TRUE to get debugging output.           */
gboolean blixem_debug_G = FALSE ;


/* Usage text. This is a duplicate of the text that is in 
 * doc/User_doc/blixem_usage.txt, so ideally we would get rid of this and use
 * the text from the file instead; for now, we must update both. */

#define USAGE_TEXT "\n\
 Blixem - display multiple alignments against a reference sequence.\n\
\n\
 Usage: blixem [options] [<sequencefile>] <datafile> [X options] \n\
\n\
   <sequencefile> contains the reference sequence in FASTA format.\n\
   <datafile> is a GFF v3 file containing alignments and other features.\n\
   If <sequencefile> is ommitted, <datafile> should contain the reference\n\
   sequence in FASTA format, below a comment line that reads ##FASTA.\n\
\n\
   Both <sequencefile> and <datafile> can be substituted by \"-\"\n\
   for reading from stdin (pipe).  If <sequencefile> is piped, the first\n\
   line should contain the sequence name and the second the sequence itself.\n\
\n\
\n\
 Options:\n\
  -t <type>, --display-type=<type>  MANDATORY\n\
    Whether to display sequences in nucleotide or protein mode. Must be one of:\n\
      N = nucleotide\n\
      P = protein\n\
\n\
  -a <names>, --alignment-names=<names>\n\
    Specify a string giving the names of the alignments, e.g. \"EST_mouse EST_human\" etc.\n\
\n\
  -c <file>, --config-file=<file>\n\
    Read configuration options from 'file'.\n\
\n\
  --abbrev-title-on\n\
    Abbreviate window title prefixes\n\
\n\
  --abbrev-title-off\n\
    Do not abbreviate window title prefixes\n\
\n\
  --compiled\n\
    Show package compile date.\n\
\n\
  --dataset\n\
    Optional string to indicate a data-set that the alignments are from.\n\
\n\
  --dotter-first-match\n\
    Call Dotter on the first match to the right of the default start coord.\n\
\n\
  --fetch-server <nodeid:port>\n\
    Causes Blixem to get sequences from a fetch server at machine 'nodeid' on the given\n\
    port (default 22100).\n\
\n\
  -h, --help\n\
    More detailed usage information.\n\
\n\
  --hide-big-picture\n\
    Hide the big picture section on start-up.\n\
\n\
  --hide-inactive-strand\n\
    Hide the inactive strand (i.e. the reverse strand, or the forward strand if the -R option\n\
    is used).\n\
\n\
  --highlight-diffs\n\
    Enable 'highlight differences' mode, where mismatches (rather than matches) are highlighted.\n\
\n\
  --invert-sort\n\
    Invert sorting order\n\
\n\
  -m <from[:to]>, --map-coords=<from[:to]>\n\
    Map the coordinate system so that the given 'from' coordinate maps to the given\n\
    'to' coordinate (or to '1' if 'to' is not given).\n\
\n\
  -n, --negate-coords\n\
    When showing the reverse strand, negate the display coordinates.\n\
\n\
  -o <n>, --offset=<n>\n\
    Offset the reference sequence coordinate system by n.\n\
\n\
  --optional-data\n\
    Parse additional data such as organism and tissue-type on start-up.\n\
\n\
  --remove-input-files\n\
    Delete the input files after they have been parsed.\n\
\n\
  -r, --reverse-strand\n\
    Indicates that the given reference sequence is the reverse strand.\n\
\n\
  --save-temp-files\n\
    Save any temporary files created by Blixem.\n\
\n\
  --show-coverage\n\
    Display the coverage section on start-up.\n\
\n\
  --sort-mode=<mode>\n\
    Default sort mode. Use --help option to see details.\n\
\n\
  --squash-matches\n\
    Compress the alignment lists on start-up.\n\
\n\
  -s <n>, --start-coord=<n>\n\
    Start with the display centred on coordinate n.\n\
\n\
  --start-next-match\n\
    Start with the display centred on the first match to the right of the default start coord.\n\
\n\
  -y <file>, --styles-file=<file>\n\
    Read color options from a key-value file. Use --help option to see details.\n\
\n\
  --version\n\
    Show package version number.\n\
\n\
  -z <start:end>, --zoom-range=<start:end>\n\
    Specify the initial range of coordinates to zoom the big picture in to.\n\
\n\
  --zoom-whole\n\
    Start with the big picture zoomed out to view the full reference sequence range.\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\n"


/* Text to show the authors, version and compile date */
#define FOOTER_TEXT "\
-----\n\
"AUTHOR_TEXT_FULL" \n\
\n\
 Reference: Sonnhammer ELL & Durbin R (1994). A workbench for Large Scale\n\
            Sequence Homology Analysis. Comput. Applic. Biosci. 10:301-307.\n\
\n\
 See http://www.sanger.ac.uk/resources/software/seqtools/ for more info.\n\
\n\
 "BLIXEM_COPYRIGHT_STRING"\n\
 "BLIXEM_LICENSE_STRING"\n\
\n\
 Version "BLIXEM_VERSION_COMPILE"\n\
\n\
"

/* Text to show the version */
#define VERSION_TEXT BLIXEM_PACKAGE_VERSION"\n"

#define HELP_TEXT "\n\
FEATURES\n\
  The prefered file format for <datafile> is GFF v3. (However, Blixem is still compatible with\n\
  older file formats such as exblx and seqbl, as used by MSPcrunch).\n\
\n\
  Blixem is mainly aimed at displaying alignments, but can also show other features such as\n\
  transcripts, variations and polyA tails. The supported SO terms are:\n%s\
\n\
SORT MODE\n\
  The sort mode is specified with the --sort-mode=<mode> argument, where <mode> is\n\
  one of the following:\n\
    s = by Score\n\
    i = by Identity\n\
    n = by Name\n\
    p = by Position\n\
\n\
  If optional data is loaded on start-up using the --optional-data argument, then the following\n\
  sort modes are also valid:\n\
    t = by Tissue type\n\
    m = by Strain\n\
    g = by Gene name\n\
    o = by Organism\n\
\n\
COLOR KEY FILE\n\
  The color key file is specified with the -y <file> or --styles-file=<file> argument. This is a .ini-\n\
  like file that specifies attributes such as fill and line colors for features from particular \n\
  sources (say EST_Human or polya_signal). The file should contain one or more source stanzas followed\n\
  by one or more key=value pairs, i.e. \n\
\n\
    [<source>]\n\
      <key>=<value>\n\
      ...\n\
\n\
  <key> can be one of:\n\
      colours:                 default colours\n\
      transcript-cds-colours:  used to specify a different colour \n\
                               for CDS sections of transcripts\n\
\n\
  <value> is a semi-colon separated list of fill and line colours of the format:\n\
      <normal|selected> <fill|border> <colour>\n\
\n\
  <colour> can be in any of the forms accepted by XParseColor; these include name \n\
           for a colour from rgb.txt, such as DarkSlateGray, or a hex specification \n\
           such as #305050.\n\
\n\
  Example:\n\
    colours=normal border #0000af ; selected border #0000af ; normal fill white ; \\\n\
            selected fill #ffddcc ; \n\
    transcript-cds-colours=normal border #0000af ; selected border #0000af ; \\\n\
            normal fill white ; selected fill #ffddcc ; \n\
\n\
  Note that selection colors will be calculated automatically if they are not\n\
  specified (a darker shade of the default color will be used when the feature is selected).\n\
\n\
MSPcrunch\n\
  To make the datafile from blast output, run MSPcrunch with option -q.\n\
\n\
 o To pipe MSPcrunch output directly to Blixem, use \"-\"\n\
   as the second parameter ([datafile]).  Example:\n\
\n\
   MSPcrunch -q <my.blast_output> | blixem <my.seq> -\n\
\n\
 o The BLAST program (blastp, blastn, blastx, tblastn, tblastx)\n\
   is automatically detected from the Blast output header by MSPcrunch\n\
   and is passed on to Blixem in the seqbl format (-q).\n\n" 


/* set default values for command lines options */
static void initCommandLineOptions(CommandLineOptions *options, char *refSeqName)
{
  options->refSeq = NULL;
  options->refSeqName = refSeqName;
  options->refSeqRange.min = UNSET_INT;
  options->refSeqRange.max = UNSET_INT;
  options->refSeqOffset = 0;
  options->startCoord = 1;
  options->startCoordSet = FALSE;
  options->mspList = NULL;
  options->columnList = NULL;
  options->geneticCode = stdcode1;
  options->activeStrand = BLXSTRAND_FORWARD;
  options->bigPictZoom = 10;          
  options->bigPictRange.min = UNSET_INT;
  options->bigPictRange.max = UNSET_INT;
  
  options->zoomWhole = FALSE;
  options->bigPictON = TRUE;          
  options->hideInactive = FALSE;         
  options->initSortColumn = BLXCOL_ID;
  options->sortInverted = FALSE;	
  options->highlightDiffs = FALSE;   
  options->dotterFirst = FALSE;	
  options->startNextMatch = FALSE;
  options->squashMatches = FALSE;
  options->optionalColumns = FALSE;
  options->saveTempFiles = FALSE;
  options->coverageOn = FALSE;
  options->abbrevTitle = FALSE;
  
  options->blastMode = BLXMODE_UNSET;
  options->seqType = BLXSEQ_NONE;
  options->numFrames = 1;
  options->mapCoords = FALSE;
  options->mapCoordsFrom = UNSET_INT;
  options->mapCoordsTo = 1; /* default to 1-based coordinate system if mapping coords but no 'to' value is specified */
  options->fetchMethods = g_hash_table_new(g_direct_hash, g_direct_equal);
  options->bulkFetchDefault = NULL;
  options->userFetchDefault = NULL;
  options->optionalFetchDefault = NULL;
  options->dataset = NULL;

  options->msgData.titlePrefix = g_strdup(BLIXEM_PREFIX);
  options->msgData.parent = NULL;
  options->msgData.statusBar = NULL;
}


static void validateOptions(CommandLineOptions *options)
{
  options->msgData.titlePrefix = options->abbrevTitle ? g_strdup(BLIXEM_PREFIX_ABBREV) : g_strdup(BLIXEM_PREFIX);
}


/* Determine the sequence type (nucleotide or peptide) from the given char */
static BlxSeqType getSeqTypeFromChar(char seqChar)
{
  BlxSeqType result = BLXSEQ_NONE;
  
  if (seqChar == 'n' || seqChar == 'N')
    result = BLXSEQ_DNA;
  else if (seqChar == 'p' || seqChar == 'P')
    result = BLXSEQ_PEPTIDE;
  else
    g_error("Bad display mode '%c'\n", seqChar);
  
  return result;
}


/* Get the sort mode from a char representing that mode */
static BlxColumnId getSortModeFromChar(char sortChar)
{
  BlxColumnId result = BLXCOL_NONE;
  
  switch (sortChar)
  {
    case 's':
      result = BLXCOL_SCORE;
      break;
    case 'i':
      result = BLXCOL_ID;
      break;
    case 'n':
      result = BLXCOL_SEQNAME;
      break;
    case 'p':
      result = BLXCOL_START;
      break;
    case 't':
      result = BLXCOL_TISSUE_TYPE;
      break;
    case 'm':
      result = BLXCOL_STRAIN;
      break;
    case 'g':
      result = BLXCOL_GENE_NAME;
      break;
    case 'o':
      result = BLXCOL_ORGANISM;
      break;

    default:
      g_error("Bad sort mode: %c\n", sortChar); 
  }
  
  return result;
}


/* Get a list of all the supported GFF types, as a string. The result should
 * be free'd with g_free. */
static char* getSupportedTypesAsString(GSList *supportedTypes)
{
  GString *resultStr = g_string_new(NULL);
  GSList *item = supportedTypes;
  
  for ( ; item; item = item->next)
    {
      BlxGffType *gffType = (BlxGffType*)(item->data);
      g_string_append_printf(resultStr, "    %s\n", gffType->name);
    }
  
  char *result = g_string_free(resultStr, FALSE);
  
  return result;
}


/* Prints usage info to stderr */
static void showUsageText(const int exitCode)
{
  /* Pring usage info followed by authors. */
  /* Send to stderr if shown due to error, otherwise to stdout */
  if (exitCode == EXIT_FAILURE)
    g_message_info("%s%s", USAGE_TEXT, FOOTER_TEXT);
  else
    g_message("%s%s", USAGE_TEXT, FOOTER_TEXT);
}


/* Prints extended usage info to stderr */
static void showHelpText(GSList *supportedTypes, const int exitCode)
{
  /* Print the standard usage text, followed by the additional help text and authors */
  GString *resultStr = g_string_new(USAGE_TEXT);

  char *supported_types_string = getSupportedTypesAsString(supportedTypes);
  
  g_string_append_printf(resultStr, HELP_TEXT, supported_types_string);
  g_string_append(resultStr, FOOTER_TEXT);

  /* Send to stderr if shown due to error, otherwise to stdout */
  if (exitCode == EXIT_FAILURE)
    g_message_info("%s", resultStr->str);
  else
    g_message("%s", resultStr->str);
  
  g_free(supported_types_string);
  g_string_free(resultStr, TRUE);
}


/* Prints version info to stderr */
static void showVersionInfo()
{
  g_message(VERSION_TEXT);  
}

/* Prints compiled date (must go to stdout for our build scripts to work) */
static void showCompiledInfo()
{
  g_message("%s\n", UT_MAKE_COMPILE_DATE());  
}


/* Entry point for blixem standalone program, you should be aware when
 * altering this that blxview.c is also compiled into acedb and that
 * blxview() is called directly by acedb code.
 */
int main(int argc, char **argv)
{
  /* Install error handlers */
  signal(SIGSEGV, errorHandler);
  signal(SIGFPE, errorHandler);

  FILE *seqfile, *FSfile;
  char seqfilename[1000] = {'\0'};
  char FSfilename[1000] = {'\0'};
  
  int install = 1;
  static gboolean showVersion = FALSE;	    /* gets set to true if blixem was called with --version option */
  static gboolean showCompiled = FALSE;	    /* gets set to true if blixem was called with --compiled option */
  
  static gboolean rm_input_files = FALSE ; /* whether to remove input files once we're done with them */
  PfetchParams *pfetch = NULL ;
  gboolean xtra_data = FALSE ;      /* whether we have an extra data file to parse */
  FILE *xtra_file = NULL ;          /* the extra data file */
  char xtra_filename[1000] = {'\0'} ;
  char *align_types = NULL ;        /* string containing alignment types, to display in the title */
  char *config_file = NULL ;        /* optional blixem config file (usually "blixemrc") */
  char *key_file = NULL ;           /* optional keyword file for passing style information */
  GError *error = NULL ;
 
  char refSeqName[FULLNAMESIZE+1] = "";

  static CommandLineOptions options;
  initCommandLineOptions(&options, refSeqName);
 
  /* Set up the GLib message handlers
   * 
   * There are two handlers: the default one for all non-critical messages, which will just log
   * output to the console, and one for critical messages and errors, which will display a 
   * pop-up message (the idea being that we don't bother the user unless it's something serious).
   * So, to get a pop-up message use g_critical, and to log a message or warning use g_message, 
   * g_warning, g_debug etc. Note that g_error is always fatal.
   */
  g_log_set_default_handler(defaultMessageHandler, &options.msgData);
  g_log_set_handler(NULL, (GLogLevelFlags)(G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION), 
                    popupMessageHandler, &options.msgData);

  /* Get the list of supported GFF types, in case we need to print them out in the usage text */
  GSList* supportedTypes = blxCreateSupportedGffTypeList(BLXSEQ_NONE);

  gtk_parse_args(&argc, &argv);

  /* Get the input args. We allow long args, so we need to create a long_options array */
  static struct option long_options[] =
    {
      {"abbrev-title-off",	no_argument,        &options.abbrevTitle, 0},
      {"abbrev-title-on",	no_argument,        &options.abbrevTitle, 1},
      {"compiled",		no_argument,        &showCompiled, 1},
      {"dataset",	        required_argument,  NULL, 0},
      {"dotter-first-match",    no_argument,        &options.dotterFirst, 1},
      {"fetch-server",          required_argument,  NULL, 0},
      {"hide-big-picture",      no_argument,        &options.bigPictON, 0},
      {"hide-inactive-strand",  no_argument,        &options.hideInactive, 1},
      {"highlight-diffs",       no_argument,        &options.highlightDiffs, 1},
      {"invert-sort",		no_argument,        &options.sortInverted, 1},
      {"optional-data",         no_argument,        &options.optionalColumns, 1},
      {"remove-input-files",    no_argument,        &rm_input_files, 1},
      {"save-temp-files",       no_argument,        &options.saveTempFiles, 1},
      {"show-coverage",         no_argument,        &options.coverageOn, 1},
      {"sort-mode",             required_argument,  NULL, 0},
      {"squash-matches",        no_argument,        &options.squashMatches, 1},
      {"start-next-match",      no_argument,        &options.startNextMatch, 1},
      {"version",		no_argument,        &showVersion, 1},
      {"zoom-whole",            no_argument,        &options.zoomWhole, 1},

      {"alignment-names",       required_argument,  0, 'a'},
      {"config-file",           required_argument,  0, 'c'},
      {"help",                  no_argument,        0, 'h'},
      {"disable-install",       no_argument,        0, 'i'}, /* "secret" option (hide from user) */
      {"map-coords",            required_argument,  0, 'm'}, 
      {"negate-coords",         no_argument,        0, 'n'}, 
      {"offset",                required_argument,  0, 'o'},
      {"reverse-strand",        no_argument,        0, 'r'},
      {"start-coord",           required_argument,  0, 's'},
      {"display-type",          required_argument,  0, 't'},
      {"extra-file",            required_argument,  0, 'x'}, /* obsolete? */
      {"styles-file",           required_argument,  0, 'y'},
      {"zoom-range",            required_argument,  0, 'z'},
      {0, 0, 0, 0}
   };

  const char  *optstring="a:c:him:no:rs:t:wx:y:z:";
  extern int   optind;
  extern char *optarg;
  int          optionIndex; /* getopt_long stores the index into the option struct here */
  int          optc;        /* the current option gets stored here */
  gboolean wait = FALSE ;


  while ((optc = getopt_long(argc, argv, optstring, long_options, &optionIndex)) != EOF)
    {
      switch (optc)
	{
        case 0:
            if (long_options[optionIndex].flag != 0)
              {
                /* we get here if getopt_long set a flag; nothing else to do */
              }
            else if (stringsEqual(long_options[optionIndex].name, "fetch-server", TRUE))
              {
                pfetch = (PfetchParams*)g_malloc(sizeof(PfetchParams)) ;
                pfetch->net_id = strtok(optarg, ":") ;
                pfetch->port = atoi(strtok(NULL, ":")) ;
              }                
            else if (stringsEqual(long_options[optionIndex].name, "sort-mode", TRUE))
              {
                options.initSortColumn = getSortModeFromChar(*optarg);
              }
            else if (stringsEqual(long_options[optionIndex].name, "dataset", TRUE))
              {
		options.dataset = g_strdup(optarg);
              }
          break; 
          
        case '?':
          break; /* getopt_long already printed an error message */
          
	case 'a':
	  align_types = g_strdup_printf("%s", optarg) ;
	  break;
	case 'c': 
	  config_file = g_strdup(optarg) ;
	  break;
	case 'h': 
          {
	    showHelpText(supportedTypes, EXIT_SUCCESS);
            exit(EXIT_SUCCESS) ;
            break;
          }
	case 'i':
	  install = 0;
	  break;
        case 'm':
          {
            options.mapCoords = TRUE;
            options.mapCoordsFrom = atoi(optarg); /* will ignore anything after ':', if it exists */
              
            /* Optionally there may be a second number after a ':' character */
            const char *cp = strchr(optarg, ':');
            if (cp)
              options.mapCoordsTo = atoi(cp + 1);
              
            break;
          }
        case 'n':
          options.negateCoords = TRUE;
          break;
	case 'o':
          options.refSeqOffset = convertStringToInt(optarg);
	  break;
        case 'r':
          options.activeStrand = BLXSTRAND_REVERSE;
          break ;
        case 's': 
	  options.startCoord = atoi(optarg);
          options.startCoordSet = TRUE;
	  break;
        case 't':
          options.seqType = getSeqTypeFromChar(*optarg);
          break;
        case 'w':
          wait = TRUE;
          break ;
        case 'x': 
	  xtra_data = TRUE ;
	  strcpy(xtra_filename, optarg);
	  break;
        case 'y':
          key_file = g_strdup(optarg) ;
          break;
        case 'z': 
          {
            int coord1 = atoi(optarg); /* will ignore anything after ':' */
            const char *cp = strchr(optarg, ':');
              
            if (cp)
              {
                int coord2 = atoi(cp + 1);
                intrangeSetValues(&options.bigPictRange, coord1, coord2);
                
                /* If the start coord hasn't already been specified on the
                 * command line, base the default start on the centre of the
                 * big picture range (can still be overridden if start coord
                 * arg is found later) */
                if (!options.startCoordSet)
                  options.startCoord = getRangeCentre(&options.bigPictRange);
              }
            else
              {
                g_warning("Invalid parameters for --zoom-range argument; expected <start:end> but got '%s'. Zoom range will be ignored.\n", optarg);
              }
              
            break;
          }
            
	default : g_error("Illegal option\n");
	}
    }

  if (wait)
    sleep(20);

  if (showVersion)
    {
      /* Just show the version info */
      showVersionInfo();
      exit(EXIT_SUCCESS);
    }

  if (showCompiled)
    {
      /* Just show the version info */
      showCompiledInfo();
      exit(EXIT_SUCCESS);
    }

  validateOptions(&options);

  /* Update the list of supported GFF types, filtering matches by the sequence type.
   * (It's quick and dirty to destroy recreate this but it's a small list and only done once.) */
  blxDestroyGffTypeList(&supportedTypes);
  supportedTypes = blxCreateSupportedGffTypeList(options.seqType);
  
  /* We expect one or two input files */
  const int numFiles = argc - optind;
  if (!(numFiles == 1 || numFiles == 2))
    {
      showUsageText(EXIT_FAILURE);
      exit(EXIT_FAILURE);
    }


  /* Add -install for private colormaps */
  if (install)
    argvAdd(&argc, &argv, "-install");

  gtk_init(&argc, &argv);

  /* mapCoords essentially does the same thing as offset, so we shouldn't be
   * given both.  Get the offset from mapCoords, if given. */
  if (options.mapCoords)
    {
      if (options.refSeqOffset)
        g_error("Error: 'map-coords' and 'offset' arguments are incompatible; please only specify one or the other.\n");
      else
        options.refSeqOffset = options.mapCoordsTo - options.mapCoordsFrom;
    }
  
  /* Set up program configuration. */
  blxInitConfig(config_file, &options, &error);
  reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);

  /* Read in the key file, if there is one */
  GSList *styles = blxReadStylesFile(key_file, &error);
  reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);

  /* Get the file names */
  if (numFiles == 1)
    {
      /* We have a single file containing both the aligments and the ref seq */
      strcpy(FSfilename, argv[optind]);
    }
  else if (numFiles == 2)
    {
      /* The ref seq is in a separate file (the first arg) */
      strcpy(seqfilename, argv[optind]);
      strcpy(FSfilename, argv[optind+1]);
    }
  else
    {
      showUsageText(EXIT_FAILURE);
    }

  /* Parse the data file containing the homol descriptions.                */
  if (!strcmp(FSfilename, "-"))
    {
      FSfile = stdin;
    }
  else if(!(FSfile = fopen(FSfilename, "r")))
    {
      g_error("Cannot open file %s\n", FSfilename);
    }
  
  /* Parser compiles lists of MSPs per type into the following array. Initialise each GList in the array to NULL */
  GArray* featureLists[BLXMSP_NUM_TYPES];
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    featureLists[typeId] = g_array_new(TRUE, FALSE, sizeof(MSP*));
  
  GList *seqList = NULL; /* parser compiles a list of BlxSequences into this list */

  char *dummyseq = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname[FULLNAMESIZE+1] = "";

  /* Create the columns */
  options.columnList = blxCreateColumns(options.optionalColumns, (options.seqType == BLXSEQ_PEPTIDE));

  /* Pass the config file to parseFS */
  GKeyFile *inputConfigFile = blxGetConfig();
  
  /* Create a temporary lookup table for BlxSequences so we can link them on GFF ID */
  GHashTable *lookupTable = g_hash_table_new(g_direct_hash, g_direct_equal);

  /* Set the blast mode from the sequence type, if given. If not, blast mode might
   * get set by parseFS (nasty for backwards compatibility; ideally we'll get rid
   * of blastmode at some point). */
  if (options.blastMode == BLXMODE_UNSET && options.seqType != BLXSEQ_NONE)
    options.blastMode = (options.seqType == BLXSEQ_PEPTIDE ? BLXMODE_BLASTX : BLXMODE_BLASTN);
  
  parseFS(&options.mspList, FSfile, &options.blastMode, featureLists, &seqList, options.columnList, supportedTypes, styles,
          &options.refSeq, options.refSeqName, &options.refSeqRange, &dummyseq, dummyseqname, inputConfigFile, lookupTable, &error) ;
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);  

  /* Now see if blast mode was set and set seqtype from it if not already set... */
  if (options.seqType == BLXSEQ_NONE && options.blastMode != BLXMODE_UNSET)
    options.seqType = (options.blastMode == BLXMODE_BLASTN ? BLXSEQ_DNA : BLXSEQ_PEPTIDE);

  /* Parse the reference sequence, if we have a separate sequence file (and it was
   * not already specified in the features file) */
  if (!options.refSeq && numFiles == 2)
    {
      /* Open the file (or stdin) */
      if (!strcmp(seqfilename, "-"))
        seqfile = stdin;
      else if(!(seqfile = fopen(seqfilename, "r")))
        g_error("Cannot open %s\n", seqfilename);
      
      /* Read in the reference sequence */
      int startCoord = UNSET_INT;
      int endCoord = UNSET_INT;
      options.refSeq = readFastaSeq(seqfile, options.refSeqName, &startCoord, &endCoord, options.seqType);
      
      if (startCoord != UNSET_INT && endCoord != UNSET_INT)
        intrangeSetValues(&options.refSeqRange, startCoord, endCoord);
      
      if (seqfile != stdin)
        fclose(seqfile);
    }
  
  if (!options.refSeq)
    g_error("No reference sequence supplied.");
  
  /* If the ref seq range still has not been set, use 1-based coords */
  if (options.refSeqRange.min == UNSET_INT && options.refSeqRange.max == UNSET_INT)
    {
      options.refSeqRange.min = 1;
      options.refSeqRange.max = strlen(options.refSeq);
    }

  if (FSfile != stdin)
    {
      fclose(FSfile) ;
    }

  /* There may an additional file containing homol data in an alternative format. */
  if (xtra_data)
    {
      if(!(xtra_file = fopen(xtra_filename, "r")))
	{
	  g_error("Cannot open %s\n", xtra_filename) ;
	}
      
      parseFS(&options.mspList, xtra_file, &options.blastMode, featureLists, &seqList, options.columnList, supportedTypes, styles,
              &options.refSeq, options.refSeqName, NULL, &dummyseq, dummyseqname, blxGetConfig(), lookupTable, &error) ;

      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      fclose(xtra_file) ;
    }

  /* Remove the input files if requested.                                  */
  if (rm_input_files)
    {
      if(seqfilename[0] != '\0' && unlink(seqfilename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of sequence input file \"%s\" failed: %s\n", seqfilename, msg) ;
          g_free(msg);
        }
        
      if(FSfilename[0] != '\0' && unlink(FSfilename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of MSP input file \"%s\" failed: %s\n", FSfilename, msg) ;
          g_free(msg);
        }
        
      if (xtra_filename[0] != '\0' && unlink(xtra_filename) != 0)
        {
          char *msg = getSystemErrorText();
          g_warning("Unlink of extra MSP sequence input file \"%s\" failed: %s\n", xtra_filename, msg) ;
          g_free(msg);
        }
    }

  
  /* Now display the alignments. (Note that TRUE signals blxview() that it is being called from
   * this standalone blixem program instead of as part of acedb. */
  if (blxview(&options, featureLists, seqList, supportedTypes, pfetch, align_types, TRUE, styles, lookupTable))
    {
      gtk_main();
    }
 
  g_free(key_file);
  g_free(config_file);
  g_hash_table_unref(lookupTable);

  g_message("Exiting Blixem\n");
  return (EXIT_SUCCESS) ;
}



