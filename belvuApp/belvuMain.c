/*  File: belvuMain.c
 *  Author: Gemma Barson, 2011-04-06
 *  Copyright (c) 2011 - 2012 Genome Research Ltd
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
 * Description: main() function for Belvu application
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvu_.h"
#include "belvuApp/belvuWindow.h"
#include "belvuApp/belvuTree.h"
#include "seqtoolsUtils/utilities.h"
#include "seqtoolsUtils/blxmsp.h"
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>


/* Usage text. This is a duplicate of the text that is in 
 * doc/User_doc/dotter_usage.txt, so ideally we would get rid of this and use
 * the text from the file instead; for now, we must update both. */

#define USAGE_TEXT "\
\n\
 Belvu - View multiple alignments in good-looking colours.\n\
\n\
 Usage: belvu [options] <multiple_alignment>|-\n\
\n\
   <multiple_alignment>|- = alignment file or pipe.\n\
\n\
 Options:\n\
  -c          Print Conservation table.\n\
  -l <file>   Load residue color code file.\n\
  -L <file>   Load markup and organism color code file.\n\
  -m <file>   Read file with matching sequence segments.\n\
  -r          Read alignment in 'raw' format (Name sequence).\n\
  -R          Do not parse coordinates when reading alignment.\n\
  -o <format> Write alignment or tree to stdout in this format and exit.\n\
                Valid formats: MSF, Mul(Stockholm), Selex, \n\
                FastaAlign, Fasta, tree.\n\
  -X <cutoff> Print UPGMA-based subfamilies at cutoff <cutoff>.\n\
  -n <cutoff> Make non-redundant to <cutoff> %identity at startup.\n\
  -Q <cutoff> Remove columns more gappy than <cutoff>.\n\
  -q <cutoff> Remove sequences more gappy than <cutoff>.\n\
  -G          Penalize gaps in pairwise comparisons.\n\
  -i          Ignore gaps in conservation calculation.\n\
  -P          Remove partial sequences at startup.\n\
  -C          Don't write coordinates to saved file.\n\
  -z <char>   Separator char between name and coordinates in saved file.\n\
  -a          Show alignment annotations on screen (Stockholm format only).\n\
  -p          Output random model probabilites for HMMER.\n\
              (Based on all residues.)\n\
  -S <order>  Sort sequences in this order.\n\
                a -> alphabetically\n\
                o -> by Swissprot organism, alphabetically\n\
                s -> by score\n\
                n -> by Neighbor-joining tree\n\
                u -> by UPGMA tree\n\
                S -> by similarity to first sequence\n\
                i -> by identity to first sequence\n\
  -s <file>   Read in file of scores.\n\
  -T <method> Tree options:\n\
                i -> Start up showing tree\n\
                I -> Start up showing only tree\n\
                d -> Show distances in tree\n\
                n -> Neighbor-joining\n\
                u -> UPGMA\n\
                c -> Don't color tree by organism\n\
                o -> Don't display sequence coordinates in tree\n\
                b -> Use Scoredist distance correction (default)\n\
                j -> Use Jukes-Cantor distance correction\n\
                k -> Use Kimura distance correction\n\
                s -> Use Storm & Sonnhammer distance correction\n\
                r -> Use uncorrected distances\n\
                p -> Print distance matrix and exit\n\
                R -> Read distance matrix instead of alignment\n\
              (only in combination with Tree routines)\n\
  -b <n>      Apply boostrap analysis with <n> bootstrap samples\n\
  -B          Print out bootstrap trees and exit\n\
              (Negative value -> display bootstrap trees on screen)\n\
  -O <label>  Read organism info after this label (default OS)\n\
  -t <title>  Set window title.\n\
  -u          Start up with uncoloured alignment (faster).\n\
  -h, --help  Show more detailed usage information\n\
  --compiled  Show package compile date\n\
  --version   Show package version number\n\
  --abbrev-title-on   Abbreviate window title prefixes\n\
  --abbrev-title-off  Do not abbreviate window title prefixes\n\
\n\
"


#define HELP_TEXT "\
 Belvu - View multiple alignments in pretty colours.\n\
\n\
 Usage: belvu [options]  <multiple_alignment>|-\n\
\n\
 <multiple_alignment>|- = file or pipe in Stockholm/Selex/MSF/Fasta format (see below).\n\
\n\
\n\
 Options:\n\
\n\
  -c          Print Conservation table.\n\
  -l <file>   Load residue color code file.\n\
\n\
              Format: <symbol> <color>\n\
              (Lines starting with # are ignored (comment lines))\n\
\n\
              Example of color code file:\n\
\n\
                  # Aroma\n\
                  F YELLOW\n\
                  Y YELLOW\n\
                  W YELLOW\n\
\n\
                  # Yuck\n\
                  D RED \n\
                  N GREEN\n\
                  X BLUE\n\
\n\
              Available colors:\n\
\n\
                    WHITE \n\
                    BLACK \n\
                    LIGHTGRAY\n\
                    DARKGRAY\n\
                    RED \n\
                    GREEN\n\
                    BLUE\n\
                    YELLOW\n\
                    CYAN \n\
                    MAGENTA\n\
                    LIGHTRED \n\
                    LIGHTGREEN \n\
                    LIGHTBLUE\n\
                    DARKRED \n\
                    DARKGREEN \n\
                    DARKBLUE\n\
                    PALERED \n\
                    PALEGREEN \n\
                    PALEBLUE\n\
                    PALEYELLOW \n\
                    PALECYAN \n\
                    PALEMAGENTA\n\
                    BROWN \n\
                    ORANGE \n\
                    PALEORANGE\n\
                    PURPLE \n\
                    VIOLET \n\
                    PALEVIOLET\n\
                    GRAY \n\
                    PALEGRAY\n\
                    CERISE \n\
                    MIDBLUE\n\
\n\
\n\
  -L <file>  Load markup and organism color code file.\n\
             Colour the markup text by residue or colour organism in tree.\n\
\n\
	     Example to set color of letters A and B:\n\
	     A GREEN\n\
	     B YELLOW\n\
\n\
	     Example to set color of organism human:\n\
	     #=OS BLUE human\n\
\n\
\n\
  -m <file>   Read file with matching sequences segments. This is used to\n\
              display a match of  a query sequence to a family.  The format\n\
              of the match is :\n\
\n\
              Line 1: Name/start-end score\n\
              Line 2: Query sequence in matching segment, no pads!\n\
              Line 3: Sequence of matching segments (qstart1 qend1 fstart1\n\
              fend2 qstart2 qend2 fstart2 fend2  etc...).\n\
\n\
              Example:\n\
\n\
              ZK673.9/238-260 21.58\n\
              CPENWVQFTGNGTQYGVCLRGFT\n\
              1 2 1 2  4 7 8 11  \n\
\n\
              NOTE: A sometimes easier way of doing this is to concatenate\n\
              the match to the end of the alignment, after a line with\n\
              exactly this string within the quotes: # matchFooter\n\
\n\
\n\
  -r          Read alignment in 'raw' format (Name sequence).\n\
\n\
              Example of raw alignment file: \n\
\n\
                  seq1_name MFILKTP\n\
                  seq1_name MYI.RTP\n\
\n\
  -R          Do not parse coordinates when reading alignment.\n\
  -o <format> Write alignment or tree to stdout in this format and exit.\n\
                Valid formats: Mul(Stockholm), Selex, MSF, \n\
                FastaAlign, Fasta, tree.\n\
  -X <cutoff> Print UPGMA-based subfamilies at cutoff <cutoff>.\n\
  -n <cutoff> Make non-redundant to <cutoff> %identity at startup.\n\
  -Q <cutoff> Remove columns more gappy than <cutoff>.\n\
  -q <cutoff> Remove sequences more gappy than <cutoff>.\n\
  -G          Penalize gaps in pairwise comparisons.\n\
  -i          Ignore gaps in conservation calculation.\n\
  -P          Remove partial sequences at startup.\n\
  -C          Don't write coordinates to saved file.\n\
  -z <char>   Separator char between name and coordinates in saved file.\n\
  -a          Show alignment annotations on screen (Stockholm format only).\n\
  -p          Output random model probabilites for HMMER.\n\
              (Based on all residues.)\n\
  -S <order>  Sort sequences in this order.\n\
                a -> alphabetically\n\
                o -> by Swissprot organism, alphabetically\n\
                s -> by score\n\
                n -> by Neighbor-joining tree\n\
                u -> by UPGMA tree\n\
                S -> by similarity to first sequence\n\
                i -> by identity to first sequence\n\
  -s <file>   Read in file of scores.A column with scores will\n\
              automatically appear after the coordinates.\n\
\n\
              Format: <score> <sequence_id>\n\
\n\
              Example of score file:\n\
\n\
                  2.78 seq_1/180-206\n\
                  2.78 seq_2/180-206\n\
                  3.79 seq_3/42-94\n\
\n\
\n\
  -T <method> Tree options:\n\
                i -> Start up showing tree\n\
                I -> Start up showing only tree\n\
                d -> Show distances in tree\n\
                n -> Neighbor-joining\n\
                u -> UPGMA\n\
                c -> Don't color tree by organism\n\
                o -> Don't display sequence coordinates in tree\n\
                b -> Use Scoredist distance correction (default)\n\
                j -> Use Jukes-Cantor distance correction\n\
                k -> Use Kimura distance correction\n\
                s -> Use Storm & Sonnhammer distance correction\n\
                r -> Use uncorrected distances\n\
                p -> Print distance matrix and exit\n\
                R -> Read distance matrix instead of alignment\n\
                     (only in combination with Tree routines)\n\
  -b <n>      Apply boostrap analysis with <n> bootstrap samples\n\
  -B          Print out bootstrap trees and exit\n\
              (Negative value -> display bootstrap trees on screen)\n\
  -O <label>  Read organism info after this label (default OS)\n\
  -t <title>  Set window title.\n\
  -u          Start up with uncoloured alignment (faster).\n\
  -h, --help  Show this help information\n\
  --compiled  Show package compile date\n\
  --version   Show package version number\n\
\
\
\
\
\
\
"


/* Common text to show at the bottom of help/usage text */
#define FOOTER_TEXT "\
 setenv BELVU_FETCH to desired sequence fetching program.\n\
 setenv BELVU_FONT_SIZE to specify window font size.\n\
 setenv BELVU_STATUSBAR_SIZE to specify statusbar font size (0 => hide statusbar).\n\
\n\
-----\n\
"AUTHOR_TEXT_FULL" \n\
\n\
 Reference: Scoredist: A simple and robust protein sequence distance estimator.\n\
            Erik LL Sonnhammer and Volker Hollich.\n\
            BMC Bioinformatics 6:108 (2005)\n\
\n\
 See http://www.sanger.ac.uk/resources/software/seqtools/ for more info.\n\
\n\
 "BELVU_COPYRIGHT_STRING"\n\
 "BELVU_LICENSE_STRING"\n\
\n\
 Version "BELVU_VERSION_COMPILE"\n\
\n\
"

/* Text to show the version */
#define VERSION_TEXT BELVU_PACKAGE_VERSION"\n"



/* Prints usage info to stderr */
static void showUsageText(const int exitCode)
{
  /* Pring usage info followed by authors */
  /* Send to stderr if shown due to error, otherwise to stdout */
  if (exitCode == EXIT_FAILURE)
    g_message_info("%s%s", USAGE_TEXT, FOOTER_TEXT);
  else
    g_message("%s%s", USAGE_TEXT, FOOTER_TEXT);
}

/* Prints more detailed usage/help info to stderr */
static void showHelpText(const int exitCode)
{
  /* Pring usage info followed by authors */
  /* Send to stderr if shown due to error, otherwise to stdout */
  if (exitCode == EXIT_FAILURE)
    g_message_info("%s%s", HELP_TEXT, FOOTER_TEXT);
  else
    g_message("%s%s", HELP_TEXT, FOOTER_TEXT);
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


/* Convert swissprot name suffixes to organisms */
static void suffix2organism(GArray *alignArr, GArray *organismArr) 
{
  int i = 0;
  char *cp = NULL;
  
  for (i = 0; i < (int)alignArr->len; ++i) 
    {
      ALN *alnp = g_array_index(alignArr, ALN*, i);
      
      if (!alnp->markup && (cp = strchr(alnp->name, '_'))) 
        {
          char *suffix = (char*)g_malloc(strlen(cp) + 1);
          strcpy(suffix, cp + 1);
          
          /* Add organism to table of organisms.  This is necessary to make all 
           sequences of the same organism point to the same place and to make a 
           non-redundant list of present organisms */
          alnp->organism = suffix;
          
          /* Only insert a new organism if it is not already in the array */
          int ip = 0;
          if (!alnArrayFind(organismArr, alnp, &ip, organism_order))
            {
	      ALN *organism = createEmptyAln();
	      alncpy(organism, alnp);
	      organism->organism = alnp->organism;
	    
              g_array_append_val(organismArr, organism);
              g_array_sort(organismArr, organism_order);
            }
          else
            {
              /* Store pointer to existing organism in ALN struct */
              ALN *alnTmp = g_array_index(organismArr, ALN*, ip);
	      g_free(alnp->organism);
              alnp->organism = alnTmp->organism;
            }
        }
    }
}


/* This is to read in sequence names and count sequences */
static void treeReadDistancesNames(BelvuContext *bc)
{
  char   *cp = NULL;
  
  char line[MAXLENGTH + 1];
  
  if (!fgets (line, MAXLENGTH, bc->treeReadDistancesPipe))
    g_error("Error reading distance matrix\n");
  
  if ((cp = strchr(line, '\n')))
    *cp = 0;
  
  int nseq = 0;
  
  while ((cp = strtok(nseq ? 0 : line, " \t")))
    {
      ALN *aln = createEmptyAln();
      
      strncpy(aln->name, cp, MAXNAMESIZE);
      
      aln->name[MAXNAMESIZE] = 0;
      aln->nr = nseq;
      g_array_insert_val(bc->alignArr, nseq, aln);
      
      nseq++;
      
      /* printf("%d  %s\n", aln.nr, aln.name); */
    }
  
  g_array_sort(bc->alignArr, nrorder);
}


static void readScores(char *filename, BelvuContext *bc)
{
  char line[MAXLENGTH+1], linecp[MAXLENGTH+1], *cp;
  FILE *file;
  int scoreLen;
  
  gboolean found = FALSE;
  gboolean warnings = FALSE;

  ALN aln;
  initAln(&aln);

  if (!(file = fopen(filename, "r")))
    g_error("Cannot open file %s\n", filename);
  
  while (!feof (file))
    { 
      if (!fgets (line, MAXLENGTH, file)) break;
      strcpy(linecp, line);
      
      initAln(&aln);
      
      if (!(cp = strtok(line, " "))) 
	g_error("Error parsing score file %s.\nLine: %s\n", filename, linecp);
      
      if (!(sscanf(cp, "%f", &aln.score)))
        g_error("Error parsing score file %s - bad score.\nLine: %s\n", filename, linecp);
      
      if (!(cp = strtok(0, "/"))) 
	g_error("Error parsing score file %s.\nLine: %s\n", filename, linecp);
      
      strncpy(aln.name, cp, MAXNAMESIZE);
      aln.name[MAXNAMESIZE] = 0;
      
      if (!(cp = strtok(0, "-"))) 
	g_error("Error parsing score file %s.\nLine: %s\n", filename, linecp);
      
      if (!(aln.start = atoi(cp)))
        g_error("Error parsing score file %s - no start coordinate.\nLine: %s\n", filename, linecp);
      
      if (!(cp = strtok(0, "\n"))) 
	g_error("Error parsing score file %s.\nLine: %s\n", filename, linecp);
      
      if (!(aln.end = atoi(cp)))
        g_error("Error parsing score file %s - no end coordinate.\nLine: %s\n", filename, linecp);
      
      int idx = 0;
      if (!alignFind(bc->alignArr, &aln, &idx)) 
        {
	  warnings = TRUE;
          /* printf("Warning: %s/%d-%d (score %.1f) not found in alignment\n", 
           aln.name, aln.start, aln.end, aln.score);*/
        }
      else
        {
	  found = TRUE;
	
          g_array_index(bc->alignArr, ALN*, idx)->score = aln.score;

          char *scoreStr = g_strdup_printf("%.1f", aln.score);
          scoreLen = strlen(scoreStr);
          g_free(scoreStr);

          if (scoreLen > bc->maxScoreLen) 
            bc->maxScoreLen = scoreLen;
        }
    }

  fclose(file);
  
  if (found)
    {
      bc->displayScores = TRUE;
    
      if (warnings)
        g_warning("Some sequences in the scores file were not found in the alignment.\n");
    }
  else
    {
      g_critical("Error reading scores file: no sequences in the scores file were found in the alignment.\n");
    }
}


int main(int argc, char **argv)
{
  FILE    
  *file, *pipe;
  
  char    
  *scoreFile = 0,
  *readMatchFile = 0,
  *colorCodesFile = 0,
  *markupColorCodesFile = 0,
  *output_format = 0,
  *optargc;
  
  int     
  i,
  output_probs = 0,
  show_ann = 0;
  
  double   
  makeNRinit = 0.0,
  init_rmGappyColumns = 0.0,
  init_rmGappySeqs = 0.0;
  
  
  
  /* Initialise GTK - do this early so we can show any errors in a pop-up dialog */
  gtk_init(&argc, &argv);
  
  
  /* Set up the GLib message handlers
   * 
   * There are two handlers: the default one for all non-critical messages, which will just log
   * output to the console, and one for critical messages and errors, which will display a 
   * pop-up message (the idea being that we don't bother the user unless it's something serious).
   * 
   * All errors and warnings will be sent to stderr, as will info messages (g_message_info).
   * Program output destined for stdout should use g_message.
   * g_debug will direct to stdout as well.
   * 
   * In summary:
   *   g_error: pop-up error message (always fatal)
   *   g_critical: pop-up error message
   *   g_warning: error message logged to stderr
   *   g_message_info: program info message sent to stderr
   *   g_message: program output message set to stdout
   *   (g_debug: not sure of usage scenarios but directs to stdout)
   *   
   */
  BlxMessageData msgData;
  msgData.titlePrefix = g_strdup(BELVU_PREFIX);
  msgData.parent = NULL;
  msgData.statusBar = NULL;

  g_log_set_default_handler(defaultMessageHandler, &msgData);
  g_log_set_handler(NULL, (GLogLevelFlags)(G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION), 
                    popupMessageHandler, &msgData);


  /* Initialise up the context with default values */
  BelvuContext *bc = createBelvuContext();


  /* Set up tree defaults */
  char treePickString[50];
  strcpy(treePickString, SWAPstr);
  
  setTreeScaleCorr(bc, bc->treeMethod);
  
  
  char *OrganismLabel = g_strdup("OS");
    
  gboolean verbose = FALSE;
  gboolean init_rmPartial = FALSE;
  
  static gboolean showHelp = FALSE;
  static gboolean showCompiled = FALSE;
  static gboolean showVersion = FALSE;
  static gboolean abbrevTitle = FALSE;

  gtk_parse_args(&argc, &argv);
  
  /* Get the input args. We allow long args, so we need to create a long_options array */
  static struct option long_options[] =
    {
      {"abbrev-title-off",	no_argument,        &abbrevTitle, 0},
      {"abbrev-title-on",	no_argument,        &abbrevTitle, 1},
      {"compiled",		no_argument,        &showCompiled, 1},
      {"version",	        no_argument,        &showVersion, 1},

      {"help",                  no_argument,        0, 'h'},
      {0, 0, 0, 0}
    };
  
  const char  *optstring="aBb:CcGhil:L:m:n:O:o:PpQ:q:RrS:s:T:t:uX:z:";
  extern int   optind;
  extern char *optarg;
  int          optionIndex; /* getopt_long stores the index into the option struct here */
  int          optc;        /* the current option gets stored here */
  
  while ((optc = getopt_long(argc, argv, optstring, long_options, &optionIndex)) != EOF)
    {
      switch (optc) 
        {
          case 0:
            /* we get here if getopt_long set a flag; nothing else to do */
            break; 
            
          case 'a': show_ann = 1;                                       break;
          case 'B': bc->outputBootstrapTrees = TRUE;                    break;
          case 'b': bc->treebootstraps = atoi(optarg);                  break;
          case 'C': bc->saveCoordsOn = FALSE;                           break;
          case 'c': verbose = TRUE;                                     break;
          case 'G': bc->penalize_gaps = TRUE;                           break;
          case 'l': colorCodesFile = g_strdup(optarg);                  break;
          case 'h': showHelp = TRUE;                                    break;
          case 'i': bc->ignoreGapsOn = TRUE;                            break;
          case 'L': markupColorCodesFile = g_strdup(optarg);            break;
          case 'm': readMatchFile = g_strdup(optarg);                   break;
          case 'n':  makeNRinit = atof(optarg);                         break;
          case 'O': strncpy(OrganismLabel, optarg, 2);                  break;
          case 'o': output_format = g_strdup(optarg);                   break;
          case 'P': init_rmPartial = TRUE;                              break;
          case 'Q': init_rmGappyColumns = atof(optarg);                 break;
          case 'q': init_rmGappySeqs = atof(optarg);                    break;
          case 'p': output_probs = 1;                                   break;
          case 'R': bc->stripCoordTokensOn = bc->saveCoordsOn = FALSE;  break;
          case 'r': bc->IN_FORMAT = RAW;                                break;
            
          case 'S': 
            switch (*optarg)
            {
              case 'a': bc->sortType = BELVU_SORT_ALPHA;                break;
              case 'o': bc->sortType = BELVU_SORT_ORGANISM;             break;
              case 's': bc->sortType = BELVU_SORT_SCORE;                break;
              case 'S': bc->sortType = BELVU_SORT_SIM;                  break;
              case 'i': bc->sortType = BELVU_SORT_ID;                   break;
              case 'n': 
                bc->treeMethod = NJ;
                bc->sortType = BELVU_SORT_TREE;                         break;
              case 'u': 
                bc->treeMethod = UPGMA; 
                bc->sortType = BELVU_SORT_TREE;                         break;
              default : g_error("Illegal sorting order: %s\n", optarg);   break;
            };
            break;
            
          case 's': scoreFile = g_strdup(optarg);                       break;
            
          case 'T': 
          for (optargc = optarg; *optargc; optargc++) 
            {
              switch (*optargc)
                {
                  case 'n': 
                    strcpy(bc->treeMethodString, NJstr);
                    bc->treeMethod = NJ;                          break;
                  case 'u': 
                    strcpy(bc->treeMethodString, UPGMAstr);
                    bc->treeMethod = UPGMA;                       break;
                  case 'c':
                    bc->treeColorsOn = FALSE;                     break;
                  case 'd':
                    bc->treeShowBranchlen = TRUE;                 break;
                  case 'I':
                    bc->onlyTree = TRUE; /* fall through */
                  case 'i':
                    bc->initTree = TRUE;                          break;
                  case 'j':
                    strcpy(bc->treeDistString, JUKESCANTORstr);
                    bc->treeDistCorr = JUKESCANTOR;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'k':
                    strcpy(bc->treeDistString, KIMURAstr);
                    bc->treeDistCorr = KIMURA;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'o':
                    bc->treeCoordsOn = FALSE;                     break;
                  case 'p':
                    bc->treePrintDistances = TRUE;  
                    bc->initTree = TRUE;                          break;
                  case 'R':
                    bc->treeReadDistancesOn = TRUE;               break;
                  case 's':
                    strcpy(bc->treeDistString, STORMSONNstr);
                    bc->treeDistCorr = STORMSONN;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'b':
                    strcpy(bc->treeDistString, SCOREDISTstr);
                    bc->treeDistCorr = SCOREDIST;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'r':
                    strcpy(bc->treeDistString, UNCORRstr);
                    bc->treeDistCorr = UNCORR;
                    setTreeScale(bc, 1.0);          break;
                  default : g_error("Illegal sorting order: %s\n", optarg);
                }
            } 
            break;
            
          case 't': strncpy(bc->Title, optarg, 255);		  break;
          case 'u': bc->displayColors = FALSE;                    break;
          case 'X': bc->mksubfamilies_cutoff = atof(optarg);      break;
          case 'z': bc->saveSeparator = *optarg;		  break;
          default : g_error("Illegal option\n");                    break;
        }
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
  
  if (showHelp)
    { 
      showHelpText(EXIT_SUCCESS);
      exit(EXIT_SUCCESS);
    }
  
  if (argc-optind < 1) 
    { 
      showUsageText(EXIT_FAILURE);
      exit(EXIT_FAILURE);
    }
  
  if (!strcmp(argv[optind], "-")) 
    {
      pipe = stdin;
      if (!*bc->Title) strcpy(bc->Title, "stdin");
    }
  else 
    {
      if (!(pipe = fopen(argv[optind], "r")))
	g_error("Cannot open file %s\n", argv[optind]);
      if (!*bc->Title) 
	strncpy(bc->Title, argv[optind], 255);
      
      bc->fileName = g_path_get_basename(argv[optind]);
      bc->dirName = g_path_get_dirname(argv[optind]);
    }
  
  bc->abbrevTitle = abbrevTitle;
  msgData.titlePrefix = abbrevTitle ? g_strdup(BELVU_PREFIX_ABBREV) : g_strdup(BELVU_PREFIX);
    
  if (bc->treeReadDistancesOn) 
    {
      /* Should this really be either or?  Problem: cannot read organism info when reading tree */
      bc->treeReadDistancesPipe = pipe;
      treeReadDistancesNames(bc);
      
      bc->initTree = TRUE;
      bc->onlyTree = TRUE;
      bc->treeCoordsOn = FALSE;
    }
  else
    {
      readFile(bc, pipe);
    }
  
  if (bc->organismArr->len == 0)
    suffix2organism(bc->alignArr, bc->organismArr);
  
  setOrganismColors(bc->organismArr);
  
  if (scoreFile) 
    readScores(scoreFile, bc);
  
  
  /* Sort by the initial sort order. If sorting by similarity or ID, then
   * we must have a selected sequence - select the first one that is not
   * a markup line. */
  if (bc->sortType == BELVU_SORT_SIM || bc->sortType == BELVU_SORT_ID)
    {
      int i = 0;
      bc->selectedAln = g_array_index(bc->alignArr, ALN*, i);

      while (i < (int)bc->alignArr->len && bc->selectedAln && bc->selectedAln->markup)
        {
          ++i;
          bc->selectedAln = g_array_index(bc->alignArr, ALN*, i);
        }
    }

  doSort(bc, bc->sortType, FALSE);
  
  if (!bc->matchFooter && readMatchFile) 
    {
      if (!(file = fopen(readMatchFile, "r"))) 
        g_error("Cannot open file %s\n", readMatchFile);
      
      readMatch(bc, file);
      fclose(file);
    }
  else if (bc->matchFooter) 
    {	 
      readMatch(bc, pipe);
      fclose(pipe);
    }
  
  if (!bc->treeReadDistancesOn) 
    {
      checkAlignment(bc);
      setConsSchemeColors(bc);
    }
  
  if (verbose)
    {
      /* Print conservation statistics */
      int i, j, max, consensus = 0 ;
      double totcons = 0.0;
    
      g_message("\nColumn Consensus        Identity       Conservation\n");
      g_message  ("------ ---------  -------------------  ------------\n");
    
      for (i = 0; i < bc->maxLen; ++i)
        {
          max = 0;
          
          for (j = 1; j < 21; j++)
            {
              if (bc->conservCount[j][i] > max)
                {
                  max = bc->conservCount[j][i];
                  consensus = j;
                }
            }
          
          g_message("%4d       %c      %4d/%-4d = %5.1f %%  %4.1f\n", 
                    i+1, b2aIndex(consensus), max, bc->alignArr->len, (double)max/bc->alignArr->len*100, bc->conservation[i]);
          totcons += bc->conservation[i];
        }
      
      g_message ("\nAverage conservation = %.1f\n", totcons/(bc->maxLen*1.0));
      
      exit(0);
    }
  
  initMarkupColors();
  initCustomColors();

  if (colorCodesFile) 
    {
      if (!(file = fopen(colorCodesFile, "r"))) 
        g_error("Cannot open file %s\n", colorCodesFile);
      
      readResidueColorScheme(bc, file, getColorArray(), TRUE);

      bc->residueScheme = BELVU_SCHEME_CUSTOM;
      bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
      bc->colorByResIdOn = FALSE;
    }
  
  if (markupColorCodesFile) 
    {
      if (!(file = fopen(markupColorCodesFile, "r"))) 
        g_error("Cannot open file %s\n", markupColorCodesFile);

      readResidueColorScheme(bc, file, getMarkupColorArray(), FALSE);

      bc->residueScheme = BELVU_SCHEME_CUSTOM;
      bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
    }

  setResidueSchemeColors(bc);

  if (makeNRinit)
    mkNonRedundant(bc, makeNRinit);
  
  if (init_rmPartial)
    rmPartialSeqs(bc);
  
  if (init_rmGappyColumns)
    rmEmptyColumns(bc, init_rmGappyColumns/100.0);
  
  if (init_rmGappySeqs) {
    rmGappySeqs(bc, init_rmGappySeqs);
    rmFinaliseGapRemoval(bc);
  }
  
  if (output_format) 
    {
      if (!strcasecmp(output_format, "Stockholm") ||
          !strcasecmp(output_format, "Mul") ||
          !strcasecmp(output_format, "Selex"))
        {
          writeMul(bc, stdout);
        }
      else if (!strcasecmp(output_format, "MSF"))
        {
          writeMSF(bc, stdout);
        }
      else if (!strcasecmp(output_format, "FastaAlign")) 
        {
	  bc->saveFormat = BELVU_FILE_ALIGNED_FASTA;
          writeFasta(bc, stdout);
        }
      else if (!strcasecmp(output_format, "Fasta")) 
        {
	  bc->saveFormat = BELVU_FILE_UNALIGNED_FASTA;
          writeFasta(bc, stdout);
        }
      else if (!strcasecmp(output_format, "tree")) 
        {
          separateMarkupLines(bc);
          
          Tree *tree = treeMake(bc, TRUE, TRUE);
          saveTreeNH(tree, tree->head, stdout);
          destroyTree(&tree);
          
          g_message(";\n");
          reInsertMarkupLines(bc);
        }
      else
        {
          g_error("Illegal output format: %s\n", output_format);
        }

      exit(0);
    }
  
  if (bc->outputBootstrapTrees && bc->treebootstraps > 0) 
    {
      treeBootstrap(bc);
      exit(0);
    } 
  
  if (output_probs) 
    {
      outputProbs(bc, stdout);
      exit(0);
    }
  
  if (bc->mksubfamilies_cutoff) 
    {
      mksubfamilies(bc, bc->mksubfamilies_cutoff);
      exit(0);
    }

  g_message_info("\n%d sequences, max len = %d\n", bc->alignArr->len, bc->maxLen);
  
  /* Try to get 8x13 font for menus, if not set on command line */
  for ( i=0; i < argc; i++)
    {
      if (!strcmp(argv[i], "-font")) 
        break;
    }
  
  if (i == argc) 
    {
      argvAdd(&argc, &argv, "-font");
      argvAdd(&argc, &argv, "8x13");
    }
  
  if (show_ann) 
    showAnnotationWindow(bc);

  if (bc->outputBootstrapTrees && bc->treebootstraps < 0)
    {	
      /* Display [treebootstraps] bootstrap trees */
      bc->treebootstrapsDisplay = TRUE;
      
      bc->treebootstraps = -bc->treebootstraps;
      
      treeBootstrap(bc);
    }

  if (!colorCodesFile) 
    {
      bc->schemeType = BELVU_SCHEME_TYPE_CONS;
      bc->consScheme = BELVU_SCHEME_BLOSUM;
      setConsSchemeColors(bc);
    }
  
  /* Create the main belvu graph display of aligned sequences. */
  if (createBelvuWindow(bc, &msgData))
    {
      gtk_main();
    }

  g_free(OrganismLabel) ;
  
  return(EXIT_SUCCESS) ;
}
