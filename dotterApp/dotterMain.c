/*  File: dotterMain.c
 *  Author: esr, 1999-08-26
 *  Copyright (c) 2010 - 2012 Genome Research Ltd
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
 * Description: main() function for Dotter application
 *----------------------------------------------------------------------------
 */

#include <dotterApp/dotter_.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxGff3Parser.h>
#include <seqtoolsUtils/blxparser.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>

#define UNSET_INT  -1

/* Usage text. This is a duplicate of the text that is in 
 * doc/User_doc/dotter_usage.txt, so ideally we would get rid of this and use
 * the text from the file instead; for now, we must update both. */

#define USAGE_TEXT "\
\n\
 Dotter - Sequence dotplots with image enhancement tools.\n\
\n\
 Usage: dotter [options] <horizontal_sequence> <vertical_sequence>  [X options]\n\
\n\
   Where <horizontal_sequence> and <vertical_sequence> are file names for FASTA files\n\
   containing the two sequences.\n\
\n\
   Allowed sequence types:   Protein  -   Protein\n\
			     DNA      -   DNA\n\
			     DNA      -   Protein\n\
\n\
 Options:\n\
  -h, --help\n\
    Show this usage information\n\
\n\
  -b <file>, --batch-save=<file>\n\
    Batch mode; save dot matrix to <file>\n\
\n\
  -e <file>, --batch-export=<file>\n\
    Batch mode; export plot to PDF file <file>\n\
\n\
  -l <file>, --load\n\
    Load dot matrix from <file>\n\
\n\
  -m <float>, --memory-limit=<float>\n\
    Memory usage limit in Mb (default 0.5)\n\
\n\
  -z <int>, --zoom\n\
    Set zoom (compression) factor\n\
\n\
  -p <int>, --pixel-factor\n\
    Set pixel factor manually (ratio pixelvalue/score)\n\
\n\
  -W <int>, --window-size\n\
    Set sliding window size. (K => Karlin/Altschul estimate)\n\
\n\
  -M <file>, --matrix-file=<file>\n\
    Read in score matrix from <file> (Blast format; Default: Blosum62).\n\
\n\
  -F <file>, --sequence-file\n\
    Read in sequences and data from <file> (replaces sequencefiles).\n\
\n\
  -f <file>, --feature-file\n\
    Read feature segments from <file>\n\
\n\
  -H, --hsp-mode\n\
    Do not calculate dotplot at startup.\n\
\n\
  -R, --reverse-greyramp\n\
    Reversed Greyramp tool at start.\n\
\n\
  -r, --reverse-horizontal\n\
    Reverse and complement horizontal_sequence (if DNA)\n\
\n\
  -v, --reverse-vertical\n\
    Reverse and complement vertical_sequence (if DNA)\n\
\n\
  -D, --disable-mirror\n\
    Don't display mirror image in self comparisons\n\
\n\
  -w, --watson-only\n\
    For DNA: horizontal_sequence top strand only (Watson)\n\
\n\
  -c, --crick-only\n\
    For DNA: horizontal_sequence bottom strand only (Crick)\n\
\n\
  -q <int>, --horizontal-offset=<int>\n\
    Horizontal_sequence offset\n\
\n\
  -s <int>, --vertical-offset=<int>\n\
    Vertical_sequence offset\n\
\n\
  --horizontal-type=p|d\n\
    Horizontal_sequence type ('p' for peptide or 'd' for DNA)\n\
\n\
  --vertical-type=p|d\n\
    Vertical_sequence type ('p' for peptide or 'd' for DNA)\n\
\n\
  --abbrev-title-on\n\
    Abbreviate window title prefixes\n\
\n\
  --abbrev-title-off\n\
    Do not abbreviate window title prefixes\n\
\n\
  --breakline_colour=<colour_str>\n\
    Set the colour used for breaklines between multiple sequences\n\
\n\
  --session_colour=<colour_str>\n\
    Set the background colour of the dotter window\n\
\n\
  --compiled\n\
    Show package compile date\n\
\n\
  --version\n\
    Show package version\n\
\n"

/* Text to show the version */
#define VERSION_TEXT DOTTER_PACKAGE_VERSION"\n"


/* Text to show the authors, version and compile date */
#define FOOTER_TEXT "\
-----\n\
"AUTHOR_TEXT_FULL" \n\
\n\
 Reference: Sonnhammer ELL & Durbin R (1995). A dot-matrix program\n\
 	    with dynamic threshold control suited for genomic DNA and protein\n\
 	    sequence analysis. Gene 167(2):GC1-10.\n\
\n\
 See http://www.sanger.ac.uk/resources/software/seqtools/ for more info.\n\
\n\
 "DOTTER_COPYRIGHT_STRING"\n\
 "DOTTER_LICENSE_STRING"\n\
\n\
 Version "DOTTER_VERSION_COMPILE"\n\
\n\
"


static void setDefaultOptions(DotterOptions *options)
{
  options->qoffset = 0;
  options->soffset = 0;
  options->selfcall = FALSE;
  options->qlen = UNSET_INT;
  options->slen = UNSET_INT;
  options->dotterZoom = 0;
  options->install = 1;
  options->pixelFacset = 0;
  options->seqInSFS = 0;
  
  options->memoryLimit = 0.0;
  
  options->savefile = NULL;
  options->exportfile = NULL;
  options->loadfile = NULL;
  options->FSfilename = NULL;
  options->mtxfile = NULL;

  options->winsize = NULL;
  options->qname = NULL;
  options->qseq = NULL;
  options->sname = NULL;
  options->sseq = NULL;
  
  options->mirrorImage = TRUE;
  options->watsonOnly = FALSE;
  options->crickOnly = FALSE;
  options->hspsOnly = FALSE;
  options->swapGreyramp = FALSE;
  options->breaklinesOn = FALSE;
  options->hozScaleRev = FALSE;
  options->vertScaleRev = FALSE;
  options->negateCoords = FALSE;
  options->abbrevTitle = FALSE;
  
  options->msgData.titlePrefix = g_strdup(DOTTER_PREFIX);
  options->msgData.parent = NULL;
  options->msgData.statusBar = NULL;

  options->breaklineColor = NULL;
  options->windowColor = NULL;
}


/* free the memory used by options (doesn't free the options struct itself) */
static void freeOptions(DotterOptions *options)
{
  g_free(options->msgData.titlePrefix);
}



static void strNamecpy(char *dest, char *src)
{
  char *cp;

  while (*src && *src == ' ')
    src++;

  strcpy(dest, src);

  if ((cp = (char *)strchr(dest, ' ')))
    *cp = 0;
  if ((cp = (char *)strchr(dest, '\n')))
    *cp = 0;

  return ;
}


static void addBreakline (MSP **MSPlist, char *name, char *desc, int pos, const int sFrame)
{
  MSP *msp = NULL;
  MSP *lastMsp = NULL;
  char *cp = NULL;

  if (!*MSPlist) 
    {
      /* Create the first msp in the list */
      msp = createEmptyMsp(&lastMsp, MSPlist);
    }
   else
    {
      /* Append a new msp to the end of the list */
      MSP *lastMsp = *MSPlist;
      while(lastMsp->next) 
	lastMsp = lastMsp->next;

      msp = createEmptyMsp(&lastMsp, MSPlist);
    }

  msp->qname = (char*)g_malloc(strlen(name)+1);
  strcpy(msp->qname, name);

  msp->desc = g_strdup(desc);
  if ((cp = (char *)strchr(msp->desc, ' ')))
    *cp = 0;
  if ((cp = (char *)strchr(msp->desc, '\n')))
    *cp = 0;

  msp->qRange.min = msp->qRange.max = pos;
  msp->fsColor = 0;
  msp->type = BLXMSP_FS_SEG;
  msp->score = 100.0;
//   insertFS(msp, "chain_separator");
}		      


/* Print the usage text to stderr */
static void showUsageText(const int exitCode)
{
  /* Send to stderr if shown due to error, otherwise to stdout */
  if (exitCode == EXIT_FAILURE)
    g_message_info("%s%s", USAGE_TEXT, FOOTER_TEXT);
  else
    g_message("%s%s", USAGE_TEXT, FOOTER_TEXT);
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


static void validateOptions(DotterOptions *options)
{
  options->msgData.titlePrefix = options->abbrevTitle ? g_strdup(DOTTER_PREFIX_ABBREV) : g_strdup(DOTTER_PREFIX);
}


int main(int argc, char **argv)
{
  DEBUG_OUT("dotter main\n");
  static DotterOptions options;
  setDefaultOptions(&options);
  
  char   
      line[MAXLINE+1],
      *curChar, *cc, *cq,  
      *qfilename, *sfilename,
      text[MAXLINE+1];

  FILE *qfile, *sfile;
  MSP *MSPlist = NULL;
  GList *seqList = NULL;
  
  /* MSPlist above is obsolete and should be replaced by featureLists, which contains all the MSPs
   * but in GLists in an array indexed by type. Initialise each GList to NULL. */
  GArray* featureLists[BLXMSP_NUM_TYPES];
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      featureLists[typeId] = g_array_new(TRUE, FALSE, sizeof(MSP*));
    }
  
  static char *dotterBinary = NULL;
  static gboolean showHelp = FALSE;
  static gboolean showVersion = FALSE;
  static gboolean showCompiled = FALSE;
  static gboolean hozScaleRev = FALSE;
  static gboolean vertScaleRev = FALSE;
  static BlxSeqType qSeqType = BLXSEQ_NONE;
  static BlxSeqType sSeqType = BLXSEQ_NONE;
  
  /* The strand stuff is a bit hacky, because dotter was originally never designed to deal with
   * reverse match seq strands, so match and ref seq strands work in different ways. If the ref seq
   * strand is reversed then the horizontal scale is reversed as well. (Because the 'active' (top) strand
   * in Blixem is always the reverse strand if the display is reversed. Note that for DNA matches both 
   * strands are always shown anyway, so the -r option essentially just swaps which strand is shown at 
   * the top of the alignment tool, and of course reverses the direction of the display.)
   * The vertical scale was never originally designed to be reversed in dotter. I've 
   * added the vertScaleRev flag in case we might want to do this in the future, but this is currently
   * always set to false, even if we have the reverse match seq strand (which is indicated with the -v option). */
  BlxStrand qStrand = BLXSTRAND_FORWARD;
  BlxStrand sStrand = BLXSTRAND_FORWARD;

  gtk_parse_args(&argc, &argv);
  
  /* Get the input args. We allow long args, so we need to create a long_options array */
  static struct option long_options[] =
    {
      {"abbrev-title-off",	no_argument,        &options.abbrevTitle, 0},
      {"abbrev-title-on",	no_argument,        &options.abbrevTitle, 1},
      {"version",		no_argument,        &showVersion, 1},
      {"compiled",		no_argument,        &showCompiled, 1},
      {"reverse-h-display",	no_argument,        &hozScaleRev, 1},
      {"reverse-v-display",	no_argument,        &vertScaleRev, 1},

      {"help",                  no_argument,        0, 'h'},
      {"batch-save",            required_argument,  0, 'b'},
      {"batch-export",          required_argument,  0, 'e'},
      {"load",                  required_argument,  0, 'l'},
      {"memory-limit",          required_argument,  0, 'm'},
      {"zoom",                  required_argument,  0, 'z'},
      {"pixel-factor",          required_argument,  0, 'p'},
      {"window-size",           required_argument,  0, 'W'},
      {"matrix-file",           required_argument,  0, 'M'},
      {"sequence-file",         required_argument,  0, 'F'},
      {"feature-file",          required_argument,  0, 'f'},
      {"hsp-mode",              no_argument,        0, 'H'},
      {"reverse-greyramp",      no_argument,        0, 'R'},
      {"reverse-horizontal",    no_argument,        0, 'r'},
      {"reverse-vertical",      no_argument,        0, 'v'},
      {"disable-mirror",        no_argument,        0, 'D'},
      {"watson-only",           no_argument,        0, 'w'},
      {"crick-only",            no_argument,        0, 'c'},
      {"horizontal-offset",     required_argument,  0, 'q'},
      {"vertical-offset",       required_argument,  0, 's'},
      {"horizontal-type",       required_argument,  0, 0},
      {"vertical-type",         required_argument,  0, 0},
      {"negate-coords",         no_argument,        0, 'N'},
      {"breakline_colour",      required_argument,  0, 0},
      {"session_colour",        required_argument,  0, 0},
      {0, 0, 0, 0}
    };

  const char  *optstring="b:cDe:f:F:hHil:M:m:Np:q:Rrs:SvW:wz:";
  extern int   optind;
  extern char *optarg;
  int          optionIndex; /* getopt_long stores the index into the option struct here */
  int          optc;        /* the current option gets stored here */
  
  while ((optc = getopt_long(argc, argv, optstring, long_options, &optionIndex)) != EOF)
    {
      switch (optc) 
        {
	  case 0:
            if (long_options[optionIndex].flag != 0)
              {
                /* we get here if getopt_long set a flag; nothing else to do */
              }
            else if (stringsEqual(long_options[optionIndex].name, "horizontal-type", TRUE))
              {
                if (*optarg == 'p')
                  qSeqType = BLXSEQ_PEPTIDE;
                else if (*optarg == 'd')
                  qSeqType = BLXSEQ_DNA;
                else
                  g_critical("Invalid value for horizontal-type argument: expected 'p' or 'd'\n");
              }                
            else if (stringsEqual(long_options[optionIndex].name, "vertical-type", TRUE))
              {
                if (*optarg == 'p')
                  sSeqType = BLXSEQ_PEPTIDE;
                else if (*optarg == 'd')
                  sSeqType = BLXSEQ_DNA;
                else
                  g_critical("Invalid value for vertical-type argument: expected 'p' or 'd'\n");
              }                
            else if (stringsEqual(long_options[optionIndex].name, "breakline_colour", TRUE))
              {
                options.breaklineColor = g_strdup(optarg);
              }
            else if (stringsEqual(long_options[optionIndex].name, "session_colour", TRUE))
              {
                options.windowColor = g_strdup(optarg);
              }
            break;
          
	  case '?':
            break; /* getopt_long already printed an error message */
	  
          case 'b': options.savefile = g_strdup(optarg);   break;
          case 'c': options.crickOnly = TRUE;              break;
          case 'D': options.mirrorImage = FALSE;           break;
          case 'e': options.exportfile = g_strdup(optarg); break;
          case 'f': options.FSfilename = g_strdup(optarg); break;
          case 'F': 
            options.seqInSFS = 1;        
            options.FSfilename = (char*)g_malloc(strlen(optarg)+1);
            strcpy(options.FSfilename, optarg);            break;
	  case 'h': 
            showHelp = TRUE;                               break;
          case 'H': options.hspsOnly = TRUE;               break;
          case 'i': options.install = 0;                   break;
          case 'l': 
            options.loadfile = (char*)g_malloc(strlen(optarg)+1);
            strcpy(options.loadfile, optarg);              break;
          case 'M': 
            options.mtxfile = (char*)g_malloc(strlen(optarg)+1);
            strcpy(options.mtxfile, optarg);               break;
          case 'm': options.memoryLimit = atof(optarg);    break;
	  case 'N': options.negateCoords = TRUE;           break;
          case 'p': options.pixelFacset = atoi(optarg);    break;
          case 'q': options.qoffset = atoi(optarg);        break;
          case 'R': options.swapGreyramp = TRUE;           break;
          case 'r': qStrand = BLXSTRAND_REVERSE;	   break;
          case 's': options.soffset = atoi(optarg);        break;
          case 'S': 
            options.selfcall = TRUE;                       break;
          case 'v': sStrand = BLXSTRAND_REVERSE;	   break;
          case 'W': 
            options.winsize = (char*)g_malloc(strlen(optarg)+1);
            strcpy(options.winsize, optarg);               break;
          case 'w': options.watsonOnly = TRUE;             break;
          case 'z': options.dotterZoom = atoi(optarg);     break;
          default : g_error("Illegal option\n");
        }
    }

  /* We're in batch mode if we've specified a save file or an export file */
  const gboolean batchMode = options.savefile || options.exportfile;

  /* We create the window in batch mode only if exporting (because this needs
   * to print the window); otherwise, don't create the window in batch mode.
   * Obviously we don't need to create the window if just showing usage info either. */
  const gboolean createWindow = (options.exportfile || !batchMode) && !showHelp && !showVersion && !showCompiled;

  /* Initialise gtk */
  if (createWindow)
    gtk_init(&argc, &argv);

  /* Set the message handlers to use our custom handlers. We normally assign a popup
   * message handler for critical messages, but don't do this in batch mode because
   * we can't have user interaction, so use the default handler for all in that case. */
  g_log_set_default_handler(defaultMessageHandler, &options.msgData);

  if (batchMode && !createWindow)
    g_log_set_handler(NULL, (GLogLevelFlags)(G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION), defaultMessageHandler, &options.msgData);
  else
    g_log_set_handler(NULL, (GLogLevelFlags)(G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION), popupMessageHandler, &options.msgData);

  if (showHelp)
    {
      showUsageText(EXIT_SUCCESS);
      exit(EXIT_SUCCESS);
    }
  
  if (showVersion)
    {
      showVersionInfo();
      exit (EXIT_SUCCESS);
    }
  
  if (showCompiled)
    {
      showCompiledInfo();
      exit (EXIT_SUCCESS);
    }
  
  validateOptions(&options);

  options.hozScaleRev = hozScaleRev;
  options.vertScaleRev = vertScaleRev;
  
  /* There's a bug if the rev-scale options are not used with the rev strand and vv, so
   * for now if one option is set, force the other. */
  if (options.hozScaleRev || qStrand == BLXSTRAND_REVERSE)
    {
      options.hozScaleRev = TRUE;
      qStrand = BLXSTRAND_REVERSE;
    }
  
  if (options.vertScaleRev || sStrand == BLXSTRAND_REVERSE)
    {
      options.vertScaleRev = TRUE;
      sStrand = BLXSTRAND_REVERSE;
    }
  
  
  if (options.selfcall) /* Blixem/Dotter calling dotter */
    {
      DEBUG_OUT("Dotter was called internally.\n");
      
      /* The input arguments (following the options) are: qname, qlen, sname, slen, dotterBinary. */
      if (argc - optind < 5 || argc - optind > 5)
        {
          g_error("Incorrect number of arguments passed to dotter from internal program call\n"); 
        }
      
      options.qname = g_strdup(argv[optind]);
      options.qlen = atoi(argv[optind + 1]);
      options.sname = g_strdup(argv[optind + 2]);
      options.slen = atoi(argv[optind + 3]);
      dotterBinary = g_strdup(argv[optind + 4]);
      
      /* Allocate memory for the sequences, now we know their lengths */
      options.qseq = (char*)g_malloc(sizeof(char) * (options.qlen+1));
      options.sseq = (char*)g_malloc(sizeof(char) * (options.slen+1));

      /* Read in the sequences from the piped input */
      DEBUG_OUT("Reading sequences from pipe...\n");

      int l = fread(options.qseq, 1, options.qlen, stdin); 
      if (l != options.qlen) 
        {
          g_error("Only read %d chars to qseq, expected %d\n", l, options.qlen);
        }
      options.qseq[options.qlen] = 0;

      l = fread(options.sseq, 1, options.slen, stdin);
      if (l != options.slen) 
        {
          g_error("Only read %d chars to sseq, expected %d\n", l, options.slen);
        }
      options.sseq[options.slen] = 0;
      DEBUG_OUT("...done.\n");

      /* Read in the features from the piped input */
      DEBUG_OUT("Reading features from pipe...\n");
      MSP *lastMsp = NULL;
      
      while (!feof (stdin))
        { 
          /* read in the blxsequences. they are separated by newlines */
          if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF)
            break;

          int numMsps = 0;
          BlxSequence *blxSeq = readBlxSequenceFromText(text, &numMsps);
          seqList = g_list_append(seqList, blxSeq);
          
          /* read in the msps for this blx sequence and add them to the msplist */
          int i = 0;
          for ( ; i < numMsps; ++i)
            {
              /* msps are also separated by newlines */
              if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF)
                break;

              MSP *msp = createEmptyMsp(&lastMsp, &MSPlist);
              readMspFromText(msp, text);
              
              /* add it to the relevant feature list. */
              featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);
              
              /* add the msp to the blxsequence */
              blxSeq->mspList = g_list_append(blxSeq->mspList, msp);
              msp->sSequence = blxSeq;
              
              /* really horrible hack */
              if (msp->type == BLXMSP_FS_SEG)
                {
                    //insertFS(msp, "chain_separator");
                  options.breaklinesOn = TRUE;
                }
            }
        }
        
      fclose(stdin);
      DEBUG_OUT("...done.\n");
    }
  else if (options.seqInSFS)
    {
      /* The -F option has been used, which replaces the input sequence files. We should therefore
       * only have 0 input arguments*/
      if (argc - optind > 0)
        {
          showUsageText(EXIT_FAILURE);
          exit(EXIT_FAILURE);
        }
    }
  else
    {

      /* The input arguments (following the options) are: qfile, sfile, so we should have 2 arguments */
      if (argc - optind < 2 || argc - optind > 2) 
        {
          showUsageText(EXIT_FAILURE);
          exit(EXIT_FAILURE);
        }

      if(!(qfile = fopen(argv[optind], "r"))) 
        {
          g_error("Cannot open %s\n", argv[optind]); 
        }
    
      qfilename = argv[optind];
      fseek(qfile, 0, SEEK_END);
      options.qlen = ftell(qfile);
      fseek(qfile, 0, SEEK_SET);
    
      if ((curChar = (char *)strrchr(argv[optind], '/')))
        {
          options.qname = g_strdup(curChar+1);
        }
      else
        {
          options.qname = g_strdup(argv[optind]);
        }

      if (!(sfile = fopen(argv[optind+1], "r"))) 
        {
          g_error("Cannot open %s\n", argv[optind+1]); 
        }
    
      sfilename = argv[optind+1];
      fseek(sfile, 0, SEEK_END);
      options.slen = ftell(sfile);
      fseek(sfile, 0, SEEK_SET);
    
      if ((curChar = (char *)strrchr(argv[optind]+1, '/')))
        {
          options.sname = g_strdup(curChar+1);
        }
      else
        {
          options.sname = g_strdup(argv[optind+1]);
        }

      /* Allocate memory for the sequences, now we know their lengths */
      options.qseq = (char*)g_malloc(sizeof(char) * (options.qlen+1));
      options.sseq = (char*)g_malloc(sizeof(char) * (options.slen+1));

      /* Read in the sequences */
      int l = 0, count = 0;
      cc = options.qseq;
      char *firstdesc = NULL;
    
      while (!feof(qfile))
        {
          if (!fgets(line, MAXLINE, qfile))
            {
              break;
            }
      
          if ((cq = (char *)strchr(line, '\n')))
            {
              *cq = 0;
            }

          /* Name headers */
          if ((cq = (char *)strchr(line, '>'))) 
            {

              cq++;
              if (++l == 1) 
                {

                  options.qname = (char*)g_malloc(strlen(cq)+1); strNamecpy(options.qname, cq);
                  firstdesc = g_strdup(cq);
                }
              else
                {

                /* Multiple sequences - add break lines */
                  if (l == 2) 
                    {

                      options.breaklinesOn = TRUE;

                      /* Second sequence - add break line to mark first sequence */
                      addBreakline (&MSPlist, qfilename, firstdesc, options.qoffset, 1);
                      
                      /* change sequence name to filename */
                      options.qname = (char*)g_malloc(strlen(qfilename)+1); strcpy(options.qname, qfilename);
                    }
                  
                  addBreakline (&MSPlist, qfilename, cq, count + options.qoffset, 1);
                }
            }
          else 
            {
              /* Read in sequence data */
              for (cq = line; *cq; cq++) 
                {
                  /* Don't know yet what type of sequence it is, so accept chars for both types */
                  if (isValidIupacChar(*cq, BLXSEQ_DNA) || isValidIupacChar(*cq, BLXSEQ_PEPTIDE))
                    {
                      *cc++ = *cq;
                      count++;
                    }
                }
            }
        }

      *cc = 0;

      if (firstdesc)
        {
          g_free(firstdesc);
          firstdesc = NULL;
        }

      l = 0, count = 0;
      cc = options.sseq;
    
      while (!feof(sfile))
        {
          if (!fgets(line, MAXLINE, sfile))
            break;
          if ((cq = (char *)strchr(line, '\n')))
            *cq = 0;

          /* Name headers */
          if ((cq = (char *)strchr(line, '>'))) {
      	cq++;
      	if (++l == 1) {
          options.sname = (char*)g_malloc(strlen(cq)+1); strNamecpy(options.sname, cq);
      	    firstdesc = g_strdup(cq);
      	}
      	else {
      	  /* Multiple sequences - add break lines */

      	    if (l == 2) {
      		options.breaklinesOn = TRUE;

      	        /* Second sequence - add break line to mark first sequence */
      	        addBreakline (&MSPlist, sfilename, firstdesc, options.soffset, 2);
      		
      		/* change sequence name to filename */
      		options.sname = (char*)g_malloc(strlen(sfilename)+1); strcpy(options.sname, sfilename);
      	    }
      	    addBreakline (&MSPlist, sfilename, cq, count + options.soffset, 2);
      	}
          }
          else 
            {
              for (cq = line; *cq; cq++) 
                {
                  /* Don't know yet what type of sequence it is, so accept chars for both types */
                  if (isValidIupacChar(*cq, BLXSEQ_DNA) || isValidIupacChar(*cq, BLXSEQ_PEPTIDE))
                    {
                      *cc++ = *cq;
                      count++;
                    }
                }
            }
        }
      
      *cc = 0;
      
      if (firstdesc)
        {
          g_free(firstdesc);
          firstdesc = NULL;
        }
    }

  BlxBlastMode blastMode = BLXMODE_UNSET;

  if (options.FSfilename) 
    {
      FILE *file;
      
      if (!strcmp(options.FSfilename, "-")) 
        {
          file = stdin;
        }
      else if(!(file = fopen(options.FSfilename, "r")))
        {
          g_error("Cannot open %s\n", options.FSfilename);
        }
      
      GSList *supportedTypes = blxCreateSupportedGffTypeList(BLXSEQ_NONE);
      GList *columnList = dotterCreateColumns();
      GError *error = NULL;
  
      /* Create a temporary lookup table for BlxSequences so we can link them on GFF ID */
      GHashTable *lookupTable = g_hash_table_new(g_direct_hash, g_direct_equal);

      parseFS(&MSPlist, file, &blastMode, featureLists, &seqList, columnList, supportedTypes, NULL, &options.qseq, options.qname, NULL, &options.sseq, options.sname, NULL, lookupTable, &error);

      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      
      finaliseBlxSequences(featureLists, &MSPlist, &seqList, columnList, 0, BLXSEQ_NONE, -1, NULL, FALSE, lookupTable);
      
      blxDestroyGffTypeList(&supportedTypes);
    }

  /* Determine sequence types */
  if (qSeqType == BLXSEQ_NONE)
    {
      GError *error = NULL;
      qSeqType = determineSeqType(options.qseq, &error);
      prefixError(error, "Error starting dotter; could not determine the sequence type for '%s'.\n", options.qname);
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  if (sSeqType == BLXSEQ_NONE)
    {
      GError *error = NULL;
      sSeqType = determineSeqType(options.sseq, &error);
      prefixError(error, "Error starting dotter; could not determine the sequence type for '%s'.\n", options.sname);
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  if (qSeqType == BLXSEQ_PEPTIDE && sSeqType == BLXSEQ_PEPTIDE) 
    {
      g_message("\nDetected sequence types: Protein vs. Protein\n");
      blastMode = BLXMODE_BLASTP;
    }
  else if (qSeqType == BLXSEQ_DNA && sSeqType == BLXSEQ_DNA) 
    {
      g_message("\nDetected sequence types: DNA vs. DNA\n");
      blastMode = BLXMODE_BLASTN;
    }
  else if (qSeqType == BLXSEQ_DNA && sSeqType == BLXSEQ_PEPTIDE) 
    {
      g_message("\nDetected sequence types: DNA vs. Protein\n");
      blastMode = BLXMODE_BLASTX;;
    }
  else
    {
      g_error("Illegal sequence types: Protein vs. DNA - turn arguments around!\n");
    }
    
  /* Add -install for private colormaps */
  if (options.install) 
    {
      argvAdd(&argc, &argv, "-install");
    }

  /* Create the dot-plot */
  dotter(blastMode, &options, qStrand, sStrand, 0, 0, MSPlist, seqList, 0);

  /* Start the gtk loop to await user interaction */
  if (createWindow)
    {
      gtk_main();
    }
  
  /* destroy the feature lists. note that the stored msps are owned
   * by the msplist, not by the feature lists */
  typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    g_array_free(featureLists[typeId], FALSE);

  freeOptions(&options);

  g_message("Exiting Dotter\n");
  return (EXIT_SUCCESS) ;
}
