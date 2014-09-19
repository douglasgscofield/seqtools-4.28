/*  File: blxparser.c
 *  Author: Erik Sonnhammer, 1993-05-17
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
 * Description: See blxparser.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/blxGff3Parser.h>
#include <seqtoolsUtils/blxparser.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>


/* The following defines give keyword strings for tags in the old file formats */
#define BLX_GAPS_TAG               "Gaps"
#define BLX_DESCRIPTION_TAG        "Description"
#define BLX_SEQUENCE_TAG           "Sequence"


/* Error codes and domain */
#define BLX_PARSER_ERROR g_quark_from_string("Parser")

typedef enum {
  BLX_PARSER_ERROR_NULL_MSP,	      /* MSP is null */
  BLX_PARSER_ERROR_INVALID_TYPE,      /* MSP does not have a valid type */
  BLX_PARSER_ERROR_NO_NAME,           /* MSP does not have a valid name */
  BLX_PARSER_ERROR_INVALID_COORDS,    /* MSP does not have valid ref seq coords */
  BLX_PARSER_ERROR_NO_STRAND,         /* MSP does not have a valid ref seq strand */
  BLX_PARSER_ERROR_INVALID_S_COORDS,  /* MSP does not have valid match coords */
  BLX_PARSER_ERROR_SEQS_DIFFER,       /* parsed data for the same sequence is different on different data lines */
  BLX_PARSER_ERROR_COMPLEMENT_FAILED, /* failed to complement the sequence */
  BLX_PARSER_ERROR_INVALID_FILE,      /* unrecognised file format */
} BlxDotterError;


typedef struct 
{
  GString *seqName;
  GString *seq;
} SeqStruct;



static char *	    nextLineOfFile(FILE *file, GString *line_string);
static char *	    nextLineOfBuffer(const char **buffer_inout, GString *line_string);


static gboolean	    parseHeaderLine(char *line, BlxBlastMode *blastMode, MSP *msp, IntRange *seq1Range, BlxParserState *parserState);

static void	    parseBody(char *line, const int lineNum, BlxBlastMode blastMode, const int resFactor, MSP **msp, GString *line_string, char **seq1, 
                              char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
                              BlxParserState *parserState, GArray* featureLists[], MSP **mspList, GList **seqList, GList *columnList, 
                              GSList *supportedTypes, GSList *styles, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, 
                              GKeyFile *keyFile, GHashTable *lookupTable);

static void	    parseEXBLXSEQBL(GArray* featureLists[], MSP **lastMsp, MSP **mspList, const BlxParserState parserState, BlxBlastMode blastMode, GString *line_string, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseEXBLXSEQBLExtended(GArray* featureLists[], MSP **lastMsp, MSP **mspList, const BlxParserState parserState, BlxBlastMode blastMode, GString *line_string, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseFsHsp(char *line, BlxBlastMode blastMode, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseFsGspHeader(char *line, BlxBlastMode blastMode, GArray* featureLists[], MSP **lastMsp, MSP **mspList, BlxParserState *parserState, GList **seqList, GList *columnList);
static void	    parseFsGspData(char *line, MSP *msp);
static void	    parseFsGff(char *line, BlxBlastMode blastMode, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseFsSeg(char *line, BlxBlastMode blastMode, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseFsXyHeader(char *line, BlxBlastMode blastMode, GArray* featureLists[], MSP **lastMsp, MSP **mspList, char **seq1, char *seq1name, char **seq2, char *seq2name, BlxParserState *parserState, GList **seqList, GList *columnList, GHashTable *lookupTable);
static void	    parseFsXyData(char *line, MSP *msp);
static void         parseFsSeqHeader(char *line, char **seq1, char *seq1name, char **seq2, char *seq2name, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, BlxParserState *parserState);
static void         parseSeqData(char *line, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, const BlxSeqType seqType);

static gboolean	    parseDescription(char **text, MSP *msp_unused) ;
static char*	    parseSequence(char **text, MSP *msp, const BlxBlastMode blastMode) ;
static void	    parseLook(MSP *msp, char *s) ;

static BlxMspType   getMspTypeFromScore(const int score);
static void	    getDesc(MSP *msp, const char *s1, const char *s2) ;
static char*	    prepSeq(const int sStart, char *inputSeq, BlxBlastMode blastMode) ;

static void         getStrandAndFrameFromString(char *text, BlxStrand *strand, int *frame);
static void         checkReversedSubjectAllowed(const MSP *msp, BlxBlastMode blastMode);

/*
 *  Globals
 */

/* These colors match those declared in systags, they must appear in the same order */	
/* If you use more than 256 colours, code WILL break (see for instance the   */
/* 'priority/colour' packing code in griddisp.c). This should not be a       */
/* problem because that's a lot of colours and these colours are NOT used    */
/* for images.                                                               */
/*                                                                           */
enum Colour    {WHITE, BLACK, LIGHTGRAY, DARKGRAY,
                RED, GREEN, BLUE,
                YELLOW, CYAN, MAGENTA,
		LIGHTRED, LIGHTGREEN, LIGHTBLUE,
		DARKRED, DARKGREEN, DARKBLUE,
		PALERED, PALEGREEN, PALEBLUE,
		PALEYELLOW, PALECYAN, PALEMAGENTA,
		BROWN, ORANGE, PALEORANGE,
		PURPLE, VIOLET, PALEVIOLET,
		GRAY, PALEGRAY,
		CERISE, MIDBLUE,
		NUM_TRUECOLORS,
                TRANSPARENT,	/* pseudocolour only for background */
		FORECOLOR,	/* pseudocolor to force box->fcol after graphColor() */
		BACKCOLOR	/* pseudocolor to force box->bcol after graphColor() */
               } ;

static const char *colorNames[NUM_TRUECOLORS] =
  {
    "WHITE", 
    "BLACK", 
    "LIGHTGRAY", 
    "DARKGRAY",
    "RED", 
    "GREEN", 
    "BLUE",
    "YELLOW", 
    "CYAN", 
    "MAGENTA",
    "LIGHTRED", 
    "LIGHTGREEN", 
    "LIGHTBLUE",
    "DARKRED", 
    "DARKGREEN", 
    "DARKBLUE",
    "PALERED", 
    "PALEGREEN", 
    "PALEBLUE",
    "PALEYELLOW", 
    "PALECYAN", 
    "PALEMAGENTA",
    "BROWN", 
    "ORANGE", 
    "PALEORANGE",
    "PURPLE", 
    "VIOLET", 
    "PALEVIOLET",
    "GRAY", 
    "PALEGRAY",
    "CERISE", 
    "MIDBLUE"
};


/* Get the factor to multiply match coords by to get display coords */
static int getResFactorFromMode(const BlxBlastMode blastMode)
{
  return (blastMode == BLXMODE_BLASTX || blastMode == BLXMODE_TBLASTX ? 3 : 1);
}


static void parseLine(char *line, const int lineNum, BlxBlastMode *blastMode, const int resFactor, MSP **msp, GString *line_string,
                      char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
                      BlxParserState *parserState, GArray *featureLists[], MSP **mspList, GList **seqList, GList *columnList, GSList *supportedTypes,
                      GSList *styles, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, GKeyFile *keyFile, GHashTable *lookupTable, GError **error)
{
  if (!line)
    return;

  const int lineLen = strlen(line);
  if (lineLen == 0 || line[0] == '\n')
    {
      return; /* empty file??? */
    }
  
  /* get rid of any trailing '\n', there may not be one if the last line of the file
   * didn't have one. */
  char *charPtr = strchr(line, '\n');
  if (charPtr)
    {
      *charPtr = 0;
    }
  
  /* Check for header info first */
  if (parseHeaderLine(line, blastMode, *msp, seq1Range, parserState))
    {
      return; 
    }
  
  if (*parserState == PARSER_START)
    {
      /* If first line was not a valid header, it's an error */
      *parserState = PARSER_ERROR;
      g_set_error(error, BLX_PARSER_ERROR, BLX_PARSER_ERROR_INVALID_FILE, "Unrecognised file header '%s'.\n", line);
    }
  else
    {
      parseBody(line, lineNum, *blastMode, resFactor, msp, line_string, 
                seq1, seq1name, seq1Range, seq2, seq2name, parserState, featureLists, mspList, seqList, columnList, supportedTypes,
                styles, readSeq, readSeqLen, readSeqMaxLen, keyFile, lookupTable);
    }
}


/* Utility to determine if we're at the end of a file stream (if given) or the 
 * end of the buffer (if given) */
static gboolean endOfFileOrBuffer(FILE *file, const char *buffer)
{
  gboolean result = TRUE;

  if (file)
    result = feof(file);
  else if (buffer)
    result = (*buffer == 0 && *buffer == '\0');
  
  return result;
}

/* Parse input from a FILE stream or a string buffer (don't call this directly - use
 * parseFS or parseBuffer) */
static void parseFileOrBuffer(MSP **MSPlist, FILE *file, const char *buffer_in, BlxBlastMode *blastMode,
                              GArray* featureLists[], GList **seqList, GList *columnList, GSList *supportedTypes, GSList *styles,
                              char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
                              GKeyFile *keyFile, GHashTable *lookupTable, GError **error)
{
  g_return_if_fail(file || buffer_in);

  const int resFactor = getResFactorFromMode(*blastMode);
  const char *buffer = buffer_in;
  
  /* Find the last MSP in the list - new MSPs will be tagged on to the end of this. We also
   * want to keep a record of the last MSP that was added in case additional information needs
   * be be added to the same MSP (e.g. for XY plot information, where we have multiple lines that
   * contain data for a single MSP). */
  MSP *msp = NULL;
  if (*MSPlist)
    {
      msp = *MSPlist;  
      while(msp->next)
	msp = msp->next;
    }

  /* Allocate reusable/extendable string as our buffer..*/
  GString *line_string = g_string_sized_new(MAXLINE + 1); 
  int lineNum = 0;

  BlxParserState parserState = PARSER_START;
  
  /* Variables for reading in SEQ data */
  char **readSeq = NULL;            /* buffer to hold the sequence we're reading in */
  int readSeqMaxLen = UNSET_INT;    /* current max length of the buffer */
  int readSeqLen = UNSET_INT;       /* current end pos of the data in the buffer */

  while (!endOfFileOrBuffer(file, buffer) && parserState != PARSER_ERROR)
    { 
      ++lineNum;
      
      line_string = g_string_truncate(line_string, 0) ;	    /* Reset buffer pointer. */

      char *line = NULL;
      
      if (file)
        line = nextLineOfFile(file, line_string);
      else
        line = nextLineOfBuffer(&buffer, line_string);
  
      parseLine(line, lineNum, blastMode, resFactor, &msp, line_string, 
                seq1, seq1name, seq1Range, seq2, seq2name, &parserState, featureLists, MSPlist, seqList, columnList, supportedTypes,
                styles, &readSeq, &readSeqLen, &readSeqMaxLen, keyFile, lookupTable, error);
    }

  g_string_free(line_string, TRUE) ;			    /* free everything, buffer and all. */

  if (seq1Range)
    {
      if (seq1Range->min == UNSET_INT && seq1Range->max == UNSET_INT && *seq1)
	{
	  /* The seq1 range was not parsed from the file; set the default range to be 1 -> strlen */
	  seq1Range->min = 1;
	  seq1Range->max = strlen(*seq1);
	}
      else if (*seq1)
	{
	  /* Check that the range is the same length as the sequence */
	  int len = strlen(*seq1);
	  if (getRangeLength(seq1Range) > len)
	    {
	      g_warning("Sequence range in features file was %d -> %d (len=%d) but parsed sequence length is %d. Limiting end of sequence range to %d.\n", seq1Range->min, seq1Range->max, getRangeLength(seq1Range), len, seq1Range->min + len - 1);
	      seq1Range->max = seq1Range->min + len - 1;
	    }
	  else if (getRangeLength(seq1Range) < len)
	    {
	      g_warning("Sequence range in features file was %d -> %d (len=%d) but parsed sequence length is %d.\n", seq1Range->min, seq1Range->max, getRangeLength(seq1Range), len);
	    }
	}
    }
  
  return ;
}


/* Function: parse a feature file (GFF or obsolete SFS format) from a FILE stream
 *
 * Assumptions:
 *
 *  *seq1 and *seq2 may or may not be allocated.
 *
 *  *seq1name and *seq2name are allocated at least 255 bytes
 *
 */
void parseFS(MSP **MSPlist, FILE *file, BlxBlastMode *blastMode,
             GArray* featureLists[], GList **seqList, GList *columnList, GSList *supportedTypes, GSList *styles,
	     char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
             GKeyFile *keyFile, GHashTable *lookupTable, GError **error)
{
  parseFileOrBuffer(MSPlist, file, NULL, blastMode, featureLists, seqList, columnList, supportedTypes,
                    styles, seq1, seq1name, seq1Range, seq2, seq2name, keyFile, lookupTable, error);
}


/* Parse file contents from a string buffer */
void parseBuffer(MSP **MSPlist, const char *buffer, BlxBlastMode *blastMode,
                 GArray* featureLists[], GList **seqList, GList *columnList, GSList *supportedTypes, GSList *styles,
                 char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
                 GKeyFile *keyFile, GHashTable *lookupTable, GError **error)
{
  parseFileOrBuffer(MSPlist, NULL, buffer, blastMode, featureLists, seqList, columnList, supportedTypes,
                    styles, seq1, seq1name, seq1Range, seq2, seq2name, keyFile, lookupTable, error);
}


/* Check that the given character is a valid DNA/RNA nucleotide. If not, give
 * a warning.  If it is not even a valid UTF8 character then GTK cannot display
 * it, so replace it with a padding character. */
static void validateIupacChar(char *inputChar, const BlxSeqType seqType)
{
  if (!isValidIupacChar(*inputChar, seqType))
    {
      if (g_utf8_validate(inputChar, 1, NULL))
        {
          g_critical("FASTA input contains invalid %s '%c'\n", (seqType == BLXSEQ_PEPTIDE ? "peptide" : "nucleotide"), *inputChar);
        }
      else
        {
          g_critical("FASTA input contains bad (non-UTF8) character - it will be replaced by '%c'\n", SEQUENCE_CHAR_INVALID);
          *inputChar = SEQUENCE_CHAR_INVALID;
        }
    }
}


/* Read in a FASTA sequence from a FASTA file or stdin */
static char *readFastaSeqFromStdin(FILE *seqfile, char *seqName, int *startCoord, int *endCoord, const BlxSeqType seqType)
{
  char *resultSeq = NULL;
  char line[MAXLINE+1];
  
  if (!fgets(line, MAXLINE, seqfile))
    {
      g_error("Error reading seqFile.\n") ;
    }
  
  sscanf(line, "%s", seqName);
  
  /* Loop through the input text and copy into an auto-expandable string */
  GString *resultStr = g_string_sized_new(5000);
  char currentChar = fgetc(seqfile);
  
  while (currentChar != '\n') 
    {
      /* Ignore whitespace/newlines */
      if (!isWhitespaceChar(currentChar) && !isNewlineChar(currentChar)) 
        {
          validateIupacChar(&currentChar, seqType);
          g_string_append_c(resultStr, currentChar);
        }      
      
      currentChar = fgetc(seqfile);
    }
  
  /* Set the result string and free the GString (but don't free its data) */
  resultSeq = g_string_free(resultStr, FALSE);
  
  return resultSeq;
}


/* Read a fasta sequence from a file */
static GArray *readFastaSeqsFromFile(FILE *seqfile, char *seqName, int *startCoord, int *endCoord, const BlxSeqType seqType)
{
  GArray *resultArr = g_array_new(FALSE, FALSE, sizeof(SeqStruct));
  char line[MAXLINE+1];

  fseek(seqfile, 0, SEEK_END);
  fseek(seqfile, 0, SEEK_SET);

  SeqStruct *currentSeqPtr = NULL;
  int curIdx = -1;

  while (!feof(seqfile))
    {
      /* Get the next line */
      if (!fgets(line, MAXLINE, seqfile)) 
        {
          break;
        }
      
      if (*line == '>')
        {
          /* This line contains the sequence name (and possibly coords), so start a new sequence */
          SeqStruct currentSeq = {g_string_new(NULL), g_string_new(NULL)};
          g_array_append_val(resultArr, currentSeq);

          ++curIdx;
          currentSeqPtr = &g_array_index(resultArr, SeqStruct, curIdx);
          
          char *linePos = line + 1;
          
          for ( ; *linePos && *linePos != ' ' && *linePos != '\n'; ++linePos)
            g_string_append_c(currentSeqPtr->seqName, *linePos);
	
	  if (startCoord && endCoord)
	    {
	      /* See if there coords specified on the line (format is ">seq_name start end") */
              ++linePos;
              sscanf(linePos, "%d %d", startCoord, endCoord);
	    }
        }
      else if (currentSeqPtr)
        {
          /* Ok, this must be a line containing sequence data. Loop through each char in the line. */
          char *linePos = line;
          
          for ( ; *linePos; linePos++)
            {
              /* Ignore whitespace/newlines */
              if (!isWhitespaceChar(*linePos) && !isNewlineChar(*linePos))
                {
                  validateIupacChar(linePos, seqType);
                  g_string_append_c(currentSeqPtr->seq, *linePos);
                }
            }
        }      
    }

  return resultArr;
}


/* Read in multiple FASTA sequences from a file and concatenate them into
 * a single sequence, using the name of the first sequence. */
static char *concatenateFastaSeqs(FILE *seqfile, char *seqName, int *startCoord, int *endCoord, const BlxSeqType seqType)
{
  GString *resultStr = g_string_new(NULL);
  
  GArray *resultArr = readFastaSeqsFromFile(seqfile, seqName, startCoord, endCoord, seqType);
  
  int i = 0;
  for ( ; i < (int)resultArr->len; ++i)
    {
      SeqStruct *curSeq = &g_array_index(resultArr, SeqStruct, i);
      
      if (curSeq && curSeq->seq)
        g_string_append(resultStr, curSeq->seq->str);
      
      if (seqName[0] == 0 && curSeq && curSeq->seqName)
        strcpy(seqName, curSeq->seqName->str);

      g_string_free(curSeq->seqName, TRUE);
      g_string_free(curSeq->seq, TRUE);
    }

  g_array_unref(resultArr);
  char *result = g_string_free(resultStr, FALSE);
  
  return result;
}


/* Read in a FASTA sequence from an input, which can be a file or stdin.
 * 
 * Format: 
 * >seq_name [start end] more words
 *
 *
 * Example with coords (e.g. for blixem):
 * 
 * >chr4-04 43747 43996 more words
 * aaatatgccattttagtattccacaattatgccaccatttggaaagaaag
 * gattgttgacagcaaaagatacccctactaaatatgggtcacagatattt
 * taccctaccctaggctatccacctaccacagaaagctgcagttcattgct
 * ggggctaccagaaatgccaagatgagattacctagagaaacaggagggtc
 * ggccaagcagcaacggccaccaccctttccggggaactccttatggccct
 *
 *
 * Example without coords (e.g. for belvu where coords are included in the name):
 * 
 * >5H1A_HUMAN/53-400 more words
 * GNACVVAAIAL...........ERSLQ.....NVANYLIG..S.LAVTDL
 * MVSVLV..LPMAAL.........YQVL..NKWTL......GQVT.CDL..
 * >5H1A_RAT more words
 * GNACVVAAIAL...........ERSLQ.....NVANYLIG..S.LAVTDL
 * MVSVLV..LPMAAL.........YQVL..NKWTL......GQVT.CDL..
 * ...
 */
char *readFastaSeq(FILE *seqfile, char *seqName, int *startCoord, int *endCoord, const BlxSeqType seqType)
{
  if (seqfile == stdin) 
    return readFastaSeqFromStdin(seqfile, seqName, startCoord, endCoord, seqType);
  else
    return concatenateFastaSeqs(seqfile, seqName, startCoord, endCoord, seqType);
}



/*************************************************
 *               Internal functions.
 *************************************************/

/* Parse a string that contains a color name from our internally-defined list of accepted
 * colors. */
static int parseColor(char *s) 
{
  int i;

  for (i = 0; i < NUM_TRUECOLORS; i++) 
    {
      if (!strcasecmp(colorNames[i], s)) 
        break;
    }
    
  return i;
}

/* Parse a line that contains shape information about a curve. */
static BlxCurveShape parseShape(char *s) 
{
    if (!strcasecmp(s, "interpolate")) 
      return BLXCURVE_INTERPOLATE;
    else if (!strcasecmp(s, "partial")) 
      return BLXCURVE_PARTIAL;
    else 
      return BLXCURVE_BADSHAPE;
}

/* Parse a string that contains information about the way the given MSP
 * should look, i.e. color and shape-interpolation. */
static void parseLook(MSP *msp, char *s)
{
    char *cp, *s2;

    /* Make copy of string to mess up */
    s2 = g_strdup(s);

    cp = strtok(s2, "," );
    while (cp) {
	
	if (parseColor(cp) != NUM_TRUECOLORS) {
	    msp->fsColor = parseColor(cp);
	}
	else if (parseShape(cp) != BLXCURVE_BADSHAPE) {
	    msp->fsShape = parseShape(cp);
	}
	else 
	    g_critical("Unrecognised Look: %s", cp);
	
	cp = strtok(0, "," );
    }
    g_free(s2);
}


/* Copy 'remainder of s1 after word s2' into msp->desc  */
static void getDesc(MSP *msp, const char *s1, const char *s2)
{
    char *cp;
    
    if (!(cp = strstr(s1, s2))) {
	g_critical("Can't find back %s in %s", s2, s1);
	return;
    }

    cp += strlen(s2)+1;

    msp->desc = g_strdup(cp);
}



/* Called when parsing HSP data. For certain modes, prepends dashes onto the start of the match
 * sequence if it does not start at 1. However, this old comment suggests that this isn't necessary:
 * "I think this is a mistake in the code, we get the whole sequence so shouldn't be prepending the
 * '-'s". */
static char* prepSeq(const int sStart, char *inputSeq, BlxBlastMode blastMode)
{
  char *result = NULL;
  
  if (blastMode == BLXMODE_TBLASTN || blastMode == BLXMODE_TBLASTX) 
    {
      result = g_strdup(inputSeq) ;
    }
  else 
    {
      result = (char*)g_malloc(sStart + strlen(inputSeq)+1);
      memset(result, SEQUENCE_CHAR_PAD,sStart); /* Fill up with dashes */
      strcpy(result + sStart - 1, inputSeq);
    }

  return result;
}



/* This routine parses MSP files that have either the exblx or seqbl format.
 * 
 * Format for both exblx and seqbl files is seven tab separated fields:
 * 
 * score reference_strand_and_frame reference_start reference_end match_start match_end match_name
 * 
 * For exblx this is optionally followed by:
 * 
 *             [match_description]
 * 
 * and optionally for seqbl by:
 * 
 *             [match_sequence]
 * 
 * e.g. for exblx
 * 
 *    11 (+2) 49052 49783 102 328 SW:YCF2_MARPO some description or other
 * 
 */
static void parseEXBLXSEQBL(GArray* featureLists[],
                            MSP **lastMsp,
                            MSP **mspList, 
                            BlxParserState parserState, 
                            BlxBlastMode blastMode, 
                            GString *line_string, 
                            GList **seqList,
                            GList *columnList,
                            GHashTable *lookupTable)
{
  int   qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL ;
  
  char sName[MAXLINE+1];
  char qframe[8] = "(+1)";
  int score = UNSET_INT;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sStart = UNSET_INT;
  int sEnd = UNSET_INT;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%d%d%s", 
	     &score, qframe, &qStart, &qEnd, 
	     &sStart, &sEnd, sName) != 7)
    {
      g_error("Incomplete MSP data in input file.\n");
    }

  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (blastMode == BLXMODE_TBLASTN)
    {
      strcpy(qframe, "(+1)");
    }
  
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);
  
  const BlxMspType mspType = getMspTypeFromScore(score);

  /* Create the new MSP */
  GError *error = NULL;
  
  /* Hack for backwards compatibility: remove the 'i' or 'x' postfix from
   * intron and exon names. */
  int len = strlen(sName);
  if (len && mspType == BLXMSP_EXON && toupper(sName[len - 1]) == 'X')
    sName[len - 1] = '\0';
  else if (len && mspType == BLXMSP_INTRON && toupper(sName[len - 1]) == 'I')
    sName[len - 1] = '\0';

  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, mspType, NULL, NULL,
                          score, UNSET_INT, 0,
                          NULL, NULL, qStart, qEnd, qStrand, qFrame,
                          sName, NULL, sStart, sEnd, BLXSTRAND_FORWARD, NULL,
                          0, lookupTable, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  /* Convert subject names to fetchable ones if from NCBI server 
         
  Rule 1: If there is a gi, use that.
  Rule 2: If no gi, use the first and last non-blank field as db:id.
  */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_strdup(msp->sname);
  	
      p = strtok(src, "|");
  	
      if (!strcasecmp(p, "GI"))
	{
	  /* Use only GI number */
	  p = strtok(0, "|");
	  strcpy(msp->sname, "gi");
	  strcat(msp->sname, ":");
	  strcat(msp->sname, p);
	}
      else
	{
	  /* Try to make a proper db:name.  Use last non-blank field */
	  char *db=p, *last;
  	    
	  p = strtok(0, "|");
	  while (p) {
	    if (*p && *p != ' ') last = p;
	    p = strtok(0, "|");
	  }
	  strcpy(msp->sname, db);
	  strcat(msp->sname, ":");
	  strcat(msp->sname, last);
	}
  	
      g_free(src);
    }
  
  if (!(cp = strstr(line, sName)))
    {
      g_error("Line does not include %s\n", sName);
    }
      
  seq_pos = cp + strlen(sName) ;
            
  qlen = mspGetQRangeLen(msp);
  slen = mspGetSRangeLen(msp);

  /* Parse the sequence, if applicable */
  char *sequence = NULL;
  
  if (parserState == EXBLX_BODY)
    {
      /* skip over description */
      while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	seq_pos++;
      if (*seq_pos && !isdigit(*seq_pos))
	{
	  while (*seq_pos && *seq_pos != ' ' && *seq_pos != '\t')
	    seq_pos++;
	  while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	    seq_pos++;
	}
    }
  else if (parserState == SEQBL_BODY)
    {
      if (blastMode == BLXMODE_TBLASTX)
	{
	  slen = (abs(sEnd - sStart) + 1)/3;
	}
      
      
      /* Line contains chars other than sequence so get the starter data...not sure this test
       * is necessary any more now we have a better mechanism of getting a whole line 
       * from a file. All this is a horrible mixture of strtok and sscanf but what else
       * to do.... */
      if (strcspn(line, "acgt") != 0)
	{
	  sequence = (char*)g_malloc(line_string->len + 1) ;
	  
	  if (sscanf(seq_pos, "%s", sequence) != 1)
	    {
	      g_error("Error parsing %s\n", line);
	    }
          
          while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
            seq_pos++;
          while (*seq_pos && *seq_pos != ' ' && *seq_pos != '\t')
            seq_pos++;
	  while (*seq_pos && (*seq_pos == ' ' || *seq_pos == '\t')) 
	    seq_pos++;
	}
    }
  
  if (sequence)
    {
      addBlxSequenceData(msp->sSequence, sequence, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  /* Parse gaps data */
  if (seq_pos)
    {
      if (!blxParseGaps(&seq_pos, msp, FALSE))
        {
          g_error("Incomplete MSP gap data for MSP '%s' [%d - %d]\n", msp->sname, msp->sRange.min, msp->sRange.max) ;
        }
    }
  
  if (parserState == SEQBL_BODY)
    { 
      checkReversedSubjectAllowed(msp, blastMode);
    }
  
  return;
}


/* This routine parses MSP files that have either the _extended_ exblx or
 * seqbl formats which borrow from the GFF format in having a set number
 * of fixed fields followed by a variable number of attribute fields.
 * These formats are distinguished by having exblx_x or seqbl_x in the
 * file comment header:
 *
 * "# seqbl_x"  or  "# exblx_x"
 *
 * Format for both exblx_x and seqbl_x files is eight tab separated fields:
 * 
 * score reference_strand_frame reference_start reference_end match_start match_end match_strand match_name
 * 
 * For exblx_x this is optionally followed by:
 * 
 *             [gaps data] [match_description]
 * 
 * and for seqbl_x by:
 * 
 *             [gaps data] [match_sequence]
 * 
 * The format of these extra fields is "tag value(s) ;", i.e.
 * 
 * "Gaps [ref_start ref_end match_start match_end]+ ;"
 * 
 * "Description the sequence description ;"
 * 
 * "Sequence aaagggtttttcccccc ;"
 * 
 * Currently acedb & ZMap export this format but other programs could also use it.
 * 
 */
static void parseEXBLXSEQBLExtended(GArray* featureLists[],
                                    MSP **lastMsp, 
                                    MSP **mspList, 
                                    BlxParserState parserState, 
                                    BlxBlastMode blastMode, 
                                    GString *line_string, 
                                    GList **seqList,
                                    GList *columnList,
                                    GHashTable *lookupTable)
{
  DEBUG_ENTER("parseEXBLXSEQBLExtended");
  
  gboolean result = FALSE ;
  int  qlen, slen;
  char *cp;
  char *line ;
  char *seq_pos = NULL;
  char *first_pos ;
  
  char sName[MAXLINE+1];
  char qframe[8] ;
  char sframe[8] ;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sStart = UNSET_INT;
  int sEnd = UNSET_INT;
  int score = UNSET_INT;

  line = line_string->str ;

  /* NOTE that sscanf will fail if the sequence name as spaces in it. The name
   * shouldn't have spaces but some do. If it does this function will probably fail
   * in trying to parse the MSP gap data. */
  if (sscanf(line, "%d%s%d%d%s%d%d%s", 
	     &score,
	     qframe, &qStart, &qEnd, 
	     sframe, &sStart, &sEnd, sName) != 8)
    {
      g_error("Incomplete MSP data in input file.\n");
    }

  /* MSPcrunch gives sframe for tblastn - restore qframe */
  if (blastMode == BLXMODE_TBLASTN)
    strcpy(qframe, "(+1)");
  
  const BlxMspType mspType = getMspTypeFromScore(score) ;
  
  /* Extract frame and strand from qframe/sframe text */
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);

  BlxStrand sStrand = BLXSTRAND_NONE;
  getStrandAndFrameFromString(sframe, &sStrand, NULL);

  /* Create the new MSP */
  GError *error = NULL;

  /* Hack for backwards compatibility: remove the 'i' or 'x' postfix from
   * intron and exon names. */
  int len = strlen(sName);
  if (len && typeIsExon(mspType) && toupper(sName[len - 1]) == 'X')
    sName[len - 1] = '\0';
  else if (len && typeIsIntron(mspType) && toupper(sName[len - 1]) == 'I')
    sName[len - 1] = '\0';
  
  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, mspType, NULL, NULL,
                          score, UNSET_INT, 0,
                          NULL, NULL, qStart, qEnd, qStrand, qFrame, 
                          sName, NULL, sStart, sEnd, sStrand, NULL,
                          0, lookupTable, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  
  
  /* Convert subject names to fetchable ones if from NCBI server 
   *  Rule 1: If there is a gi, use that.
   *  Rule 2: If no gi, use the first and last non-blank field as db:id.
   */
  if (strchr(msp->sname, '|'))
    {
      char *p, *src;
  	
      src = g_strdup(msp->sname);
  	
      p = strtok(src, "|");
  	
      if (!strcasecmp(p, "GI"))
	{
	  /* Use only GI number */
	  p = strtok(0, "|");
	  strcpy(msp->sname, "gi");
	  strcat(msp->sname, ":");
	  strcat(msp->sname, p);
	}
      else
	{
	  /* Try to make a proper db:name.  Use last non-blank field */
	  char *db=p, *last;
  	    
	  p = strtok(0, "|");
	  while (p) {
	    if (*p && *p != ' ') last = p;
	    p = strtok(0, "|");
	  }
	  strcpy(msp->sname, db);
	  strcat(msp->sname, ":");
	  strcat(msp->sname, last);
	}
  	
      g_free(src);
    }
  
  if (!(cp = strstr(line, sName)))
    {
      g_error("Line does not include %s\n", sName);
    }

  seq_pos = cp + strlen(sName) ;
            
  qlen = mspGetQRangeLen(msp);
  slen = mspGetSRangeLen(msp);

  /* Now read attributes. */
  char *sequence = NULL;
  
  if (seq_pos && *seq_pos)
    {
      first_pos = seq_pos ;
      result = TRUE ;
      while (result)
	{
	  if (!(seq_pos = strtok(first_pos, " ")))
	    {
	      result = FALSE ;
	      break ;
	    }
	  else
	    {
	      first_pos = NULL ;
              
	      /* Get first word and then parse.... */
	      if ((strstr(seq_pos, BLX_GAPS_TAG)))
		{
		  if (!(result = blxParseGaps(&seq_pos, msp, TRUE)))
		    g_error("Incomplete MSP gap data\n") ;
		}
	      else if (parserState == EXBLX_BODY && (strstr(seq_pos, BLX_DESCRIPTION_TAG)))
		{
		  result = parseDescription(&seq_pos, msp);
                  
		  if (!result)
                    {
                      g_error("Bad description data\n");
                    }
		}
	      else if ((parserState == SEQBL_X_BODY || mspIsVariation(msp)) && (strstr(seq_pos, BLX_SEQUENCE_TAG)))
		{
		  if (blastMode == BLXMODE_TBLASTX)
		    {
		      slen = mspGetSRangeLen(msp) / 3;
		    }
                  
                  sequence = parseSequence(&seq_pos, msp, blastMode);
                  result = (sequence != NULL);
                  
		  if (!result)
                    {
                      g_error("Bad sequence data\n");
                    }
		}
              
	    }
	}
    }
  
  if (sequence)
    {
      addBlxSequenceData(msp->sSequence, sequence, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  if (parserState == SEQBL_X_BODY)
    {
      if (blastMode == BLXMODE_TBLASTX)
	{
	  slen = mspGetSRangeLen(msp) / 3;
	}
	

    }

  if (parserState == SEQBL_X_BODY)
    { 
      checkReversedSubjectAllowed(msp, blastMode);
    }
  
  DEBUG_EXIT("parseEXBLXSEQBLExtended");
  return ;
}


/* Get the next line in a buffer */
static char *nextLineOfBuffer(const char **buffer_inout, GString *line_string)
{
  char *result = NULL;
  g_return_val_if_fail(buffer_inout && line_string, result);

  const char *buffer = *buffer_inout;
  char *cp = strchr(buffer, '\n');
  
  if (cp)
    {
      /* Append the text up to the newline */
      int len = cp - buffer;
      g_string_append_len(line_string, buffer, len);
      
      /* Chop off the bit before the newline, so buffer is ready for the next iteration */
      *buffer_inout = cp + 1;
    }
  else
    {
      /* Just a single line in the buffer, so append the whole thing */
      g_string_append(line_string, buffer);
      *buffer_inout = 0;
    }

  result = line_string->str ;
  return result;
}


/* Read a line from a file, gets the whole line no matter how big...until you run out
 * of memory.
 * Returns NULL if there was nothing to read from the file, otherwise returns the
 * line read which is actually the string held within the GString passed in.
 * Crashes if there is a problem with the file. */
static char *nextLineOfFile(FILE *file, GString *line_string)
{
  enum {BLX_BUF_SIZE = 4096} ;				    /* Vague guess at initial length. */
  char *result = NULL ;
  gboolean line_finished ;
  char buffer[BLX_BUF_SIZE] ;

  g_assert(file) ;

  line_finished = FALSE ;
  while (!line_finished)
    {

      if (!fgets(buffer, BLX_BUF_SIZE, file))
	{
	  if (feof(file))
	    line_finished = TRUE ;
	  else
	    g_error("NULL value returned on reading input file\n") ;
	}
      else
	{
	  if (buffer[0] == '\0')
	    {
	      line_finished = TRUE ;
	      result = line_string->str ;
	    }
	  else
	    {
	      int line_end ;

	      g_string_sprintfa(line_string, "%s", &buffer[0]) ;

	      line_end = strlen(&buffer[0]) - 1 ;

	      if (buffer[line_end] == '\n' || buffer[line_end] == '\r')
		{
		  line_finished = TRUE ;
		  result = line_string->str ;
		}
	    }
	}
    }

  return result ;
}


/* Expects a string in the format "Gaps [ref_start ref_end match_start match_end]+ ; more text....."
 * and parses out the coords. Only spaces allowed in the string, not tabs. Moves text
 * to first char after ';'. 
 * The hasGapsTag is a bit of a hack to also allow this function to work with the old style
 * exblx file format that does not have a "Gaps" tag on the front of the gaps data. (In this case
 * we need to do something slightly different with strtok to start reading in the correct place
 * in the string.) */
gboolean blxParseGaps(char **text, MSP *msp, const gboolean hasGapsTag)
{
  gboolean result = TRUE;

  CoordRange *currentGap = NULL;
  gboolean finished = FALSE;

  char *currentGapStr = hasGapsTag ? strtok(NULL, " ") : strtok(*text, " ") ;

  while (result && !finished && currentGapStr)
    {
      int i  = 0;

      /* Loop through the four coords in this range */
      for (i = 0 ;  i < 4 ; i++)
	{
	  if (!currentGapStr)
	    {
	      result = FALSE;
	      break;
	    }
	  else if (*currentGapStr == ';')
	    {
	      finished = TRUE;
	      break;
	    }
	     
	  switch (i)
	  {
	    case 0:
	    {
	      /* First value is start of subject sequence range. Create the range struct */
              currentGap = (CoordRange*)g_malloc0(sizeof(CoordRange));
              msp->gaps = g_slist_append(msp->gaps, currentGap);
	      currentGap->sStart = convertStringToInt(currentGapStr);
	      break;
	    }
	      
	    case 1:
	    {
	      /* Second value is end of subject sequence range. Order values so that
	       * sStart is less than sEnd if we have the forward strand or v.v. if the reverse. */
	      currentGap->sEnd = convertStringToInt(currentGapStr);
	      sortValues(&currentGap->sStart, &currentGap->sEnd, mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
	      break;
	    }
	      
	    case 2:
	    {
	      /* Third value is start of reference sequence range */
	      currentGap->qStart = convertStringToInt(currentGapStr);
	      break;
	    }
	      
	    case 3:
	    {
	      /* Fourth value is end of reference sequence range. Order values so that
	       * qStart is less than qEnd if ref sequence is forward strand or v.v. if the reverse. */
	      currentGap->qEnd = convertStringToInt(currentGapStr);
	      sortValues(&currentGap->qStart, &currentGap->qEnd, mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
	    }
	  }

	  currentGapStr = strtok(NULL, "\t ,") ; 
	}  
    }

  return result ;
}


/* Description, should just be plain text, format is:
 * 
 * "Description text ; more text....."
 * 
 * text following Description must not contain ';'. Moves text to first
 * char after ';'. */
static gboolean parseDescription(char **text, MSP *msp)
{
  gboolean result = TRUE ;
  char *cp ;

  cp = strtok(NULL, ";") ;

  if (cp && *cp)
    {
      msp->desc = g_strdup(cp) ;
    }

  return result ;
}


/* Find out the sequence type to display (nucleotide or peptide) based on the blast mode */
static BlxSeqType getSeqTypeFromBlastMode(const BlxBlastMode blastMode, GError **error)
{
  BlxSeqType result = BLXSEQ_NONE;
  
  switch (blastMode)
    {
      case BLXMODE_BLASTP: /* fall through */
      case BLXMODE_BLASTX: /* fall through */
      case BLXMODE_TBLASTX:
        result = BLXSEQ_PEPTIDE;
        break;
        
      case BLXMODE_BLASTN: /* fall through */
      case BLXMODE_TBLASTN:
        result = BLXSEQ_DNA;
        break;
        
      default:
        g_set_error(error, BLX_ERROR, 1, "Unknown blast mode '%d'.\n", blastMode);
        break;
    };
  
  return result;
}


/* Description, should just be plain text, format is:
 * 
 * "Sequence text ; "
 * 
 * text following Sequence must not contain ';'. Moves text to first
 * char after ';'. */
static char* parseSequence(char **text, MSP *msp, const BlxBlastMode blastMode)
{
  char *result = NULL;

  char *startPtr = *text ;
  startPtr = strtok(NULL, "\t ") ;				    /* skip "Sequence" */

  const int origLen = (int)strlen(startPtr);
  
  int validLen = 0;
  char *cp = startPtr;
  
  GError *tmpError = NULL;
  BlxSeqType seqType = getSeqTypeFromBlastMode(blastMode, &tmpError);
  reportAndClearIfError(&tmpError, G_LOG_LEVEL_ERROR);
  
  while (isValidIupacChar(*cp, seqType))
    {
      ++cp;
      ++validLen;
    }
  
  if (validLen < 1 || validLen < (int)strlen(startPtr))
    {
      g_error("Error parsing sequence data for MSP '%s'; sequence is only valid up to %d (out of length %d).\n",
		mspGetSName(msp), origLen, validLen) ;
    }
  else
    {
      result = (char*)g_malloc(validLen + 1) ;
	  
      if (sscanf(startPtr, "%s", result) != 1)
	{
	  g_error("Error parsing sequence data '%s' for MSP '%s'\n", startPtr, mspGetSName(msp)) ;
	}
    }

  return result ;
}


/* Utility called by parseFS to parse the header info of a line from a file. Returns true if the
 * line was completely processed, false if further processing on the same line is still required */
static gboolean parseHeaderLine(char *line, BlxBlastMode *blastMode, MSP *msp, IntRange *seq1Range, BlxParserState *parserState)
{
  //DEBUG_ENTER("parseHeaderLine(parserState=%d)", *parserState);
  
  gboolean processed = FALSE;
  
  if (!strncasecmp(line, "##gff-version", 13))
    {
      /* Check it's GFF version 3. Loop past any whitespace first. */
      char *cp = line + 13;
      while (cp && (*cp == ' ' || *cp == '\t'))
	{
	  ++cp;
	} 
    
      if (*cp != '3')
	{
	  *parserState = PARSER_ERROR;
	  g_critical("Error parsing GFF file: GFF version '%s' is not supported (only version 3 is supported).", cp);
	}
      else
	{
	  *parserState = GFF_3_HEADER ;
	  processed = TRUE;
	}
    }
  else if (!strncasecmp(line, "##FASTA", 7))
    {
      *parserState = FASTA_SEQ_HEADER ;
      processed = TRUE;
    }
  else if (!strncasecmp(line, "# seqbl_x", 9))
    {
      *parserState = SEQBL_X_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx_x", 9))
    {
      *parserState = EXBLX_X_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# seqbl", 7))
    {
      /* Only for backwards compatibility */
      *parserState = SEQBL_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# exblx", 7))
    {
      /* Only for backwards compatibility */
      *parserState = EXBLX_BODY ;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=HSP", 13) || 
	   !strncasecmp(line, "# SFS type=HSP", 14))
    {
      *parserState = FS_HSP_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GSP", 13) || 
	   !strncasecmp(line, "# SFS type=GSP", 14))
    {
      *parserState = FS_GSP_HEADER;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# dotter feature format 2", 25) ||
	   !strncasecmp(line, "# FS type=SEG", 13) ||
	   !strncasecmp(line, "# SFS type=SEG", 14))
    {
      *parserState = FS_SEG_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=GFF", 13) || 
	   !strncasecmp(line, "# SFS type=GFF", 14))
    {
      *parserState = FS_GFF_BODY;
      processed = TRUE ;
    }
  else if (!strncasecmp(line, "# FS type=XY", 12) ||
	   !strncasecmp(line, "# SFS type=XY", 13))
    {
      /* More info exists in the header line for this data type, so indicate we need to parse
       * header info and return false to indicate we have not finished with this line */
      *parserState = FS_XY_HEADER;
      processed = FALSE;
    }
  else if (!strncasecmp(line, "# FS type=SEQ", 13) ||
	   !strncasecmp(line, "# SFS type=SEQ", 14))
    {
      *parserState = FS_SEQ_HEADER;
      processed = FALSE; /* more info exists in the header line for this data type */
    }
  else if (!strncasecmp(line, "# FS type=", 10) ||
	   !strncasecmp(line, "# SFS type=", 11))
    {
      g_error("Unrecognised SFS type: %s\n", line);
    }
  else if (*line == '#' && *parserState == GFF_3_HEADER)
    {
      /* If we're in a GFF header we want to take a look at additional comment lines. */
      processed = FALSE ;
    }
  else if (*parserState == GFF_3_HEADER)
    {
      /* We were processing GFF3 header lines (which all start with '#'), but this line does not
       * start with '#', so it must be the start of the GFF3 body. */
      *parserState = GFF_3_BODY;
      processed = FALSE;
    }
  else if (*line == '#')
    {
      /* Very ugly; only for backwards compatibility */
      /* Changed to soft parsing (unknown labels ignored) so that
       any comment can be used */
      if (!strncasecmp(line, "# blastp" , 8))
        *blastMode = BLXMODE_BLASTP;
      else if (!strncasecmp(line, "# tblastn", 9))
        *blastMode = BLXMODE_TBLASTN;
      else if (!strncasecmp(line, "# tblastx", 9))
        *blastMode = BLXMODE_TBLASTX;
      else if (!strncasecmp(line, "# blastn" , 8))
        *blastMode = BLXMODE_BLASTN;
      else if (!strncasecmp(line, "# blastx" , 8))
        *blastMode = BLXMODE_BLASTX;
      else if (!strncasecmp(line, "# hspgaps", 9))
	{
	  /* unused */
	}
      else if (!strncasecmp(line, "# DESC ", 7) &&
	       (*parserState == FS_HSP_BODY || *parserState == FS_GSP_HEADER || 
		*parserState == SEQBL_BODY))
	{
	  if (msp)
	    getDesc(msp, line, mspGetSName(msp));
	}
      else if (!strncasecmp(line, "# sequence-region", 17) &&
               (*parserState == EXBLX_BODY || *parserState == EXBLX_X_BODY || *parserState == SEQBL_BODY || *parserState == SEQBL_X_BODY))
        {
          /* read in the reference sequence name and range */
          char qName[MAXLINE + 1];
          int qStart = UNSET_INT;
          int qEnd = UNSET_INT;
          
          if (sscanf(line, "# sequence-region%s%d%d", qName, &qStart, &qEnd) < 3)
            {
              g_warning("Error parsing sequence-region line in input file. Sequence range was not set. \"%s\"\n", line);
            }
          else
            {
              intrangeSetValues(seq1Range, qStart, qEnd);
            }
        }
	
      processed = TRUE;
    }
  else if (*line == '>' && 
           (*parserState == FASTA_SEQ_BODY || *parserState == FASTA_SEQ_IGNORE))
    {
      /* We can have multiple fasta sequences after the ##FASTA comment line.
       * They don't have another ##FASTA comment so we need to look for the 
       * '>' which indicates the header line of the fasta data. Only do this if
       * we've already seen the ##FASTA comment, though. */
      *parserState = FASTA_SEQ_HEADER;
      processed = FALSE; /* more info exists on the header line that we need to parse out */      
    }

  //DEBUG_EXIT("parseHeaderLine returning processed = %d, (parserState = %d)", processed, *parserState);
  return processed ;
}


/* Parse data from the given line into an MSP of type FS_HSP */
static void parseFsHsp(char *line, 
                       BlxBlastMode blastMode,
                       GArray* featureLists[], 
                       MSP **lastMsp,
                       MSP **mspList,
                       GList **seqList,
                       GList *columnList,
                       GHashTable *lookupTable)
{
  char qName[MAXLINE+1];
  char sName[MAXLINE+1];
  char seq[MAXLINE+1];
  char qframe[8] = "(+1)";
  char sframe[8] = "(+1)";
  int qStart = UNSET_INT;
  int sStart = UNSET_INT;
  int qEnd = UNSET_INT;
  int sEnd = UNSET_INT;
  int score = UNSET_INT;
  
  /* <score> <qname> <qframe> <qstart> <qend> <sname> <sframe> <sstart> <ssend> <sequence> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s%s%d%d%s", 
	     &score, 
	     qName, qframe+1, &qStart, &qEnd, 
	     sName, sframe+1, &sStart, &sEnd,
	     seq) != 10)
    {
      g_error("Error parsing data, type FS_HSP: \"%s\"\n", line);
    }
  
  /* Get the strand and frame from the frame strings */
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);

  BlxStrand sStrand = BLXSTRAND_NONE;
  int sFrame = UNSET_INT;
  getStrandAndFrameFromString(sframe, &sStrand, &sFrame);

  /* Get the sequence data */
  char *sSeq = prepSeq(sStart, seq, blastMode);

  /* Create the new MSP */
  GError *error = NULL;

  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_HSP, NULL,  NULL,
                          score, UNSET_INT, 0, 
                          NULL, qName, qStart, qEnd, qStrand, qFrame, 
                          sName, NULL, sStart, sEnd, sStrand, sSeq, 0, lookupTable, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
  checkReversedSubjectAllowed(msp, blastMode);
}


/* Parse data from the given line, which contains feature-series GSP header info */
static void parseFsGspHeader(char *line, 
                             BlxBlastMode blastMode, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             BlxParserState *parserState, 
                             GList **seqList,
                             GList *columnList)
{
  /* Data goes into a new MSP */
//  msp->type = BLXMSP_GSP;

  /* Will write this as soon as MSPcrunch generates it */
  g_warning("Parser for GSP data not implemented\n");
  
  /* Start parsing the GSP data next */
  *parserState = FS_GSP_BODY ;
}


/* Parse a line that contains feature-series GSP data */
static void parseFsGspData(char *line, MSP *msp)
{
  g_warning("Parser for GSP data not implemented\n");
}


/* Parse data from the given line into an MSP of type FS_SEG (Feature Series Segment) */
static void parseFsSeg(char *line, 
                       BlxBlastMode blastMode, 
                       GArray* featureLists[], 
                       MSP **lastMsp, 
                       MSP **mspList, 
                       GList **seqList,
                       GList *columnList,
                       GHashTable *lookupTable)
{
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  int score = UNSET_INT;
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  
  /* <score> <sequencename> <seriesname> <start> <end> <look> [annotation] */
  if (sscanf(line, "%d%s%s%d%d%s",
	     &score, qName, series, &qStart, &qEnd, look) != 6)
    {
      g_error("Error parsing data, type FS_SEG: \"%s\"\n", line);
    }
  
  
  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_FS_SEG, NULL, NULL,
                          UNSET_INT, UNSET_INT, 0, 
                          NULL, qName, qStart, qEnd, BLXSTRAND_NONE, 1, 
                          series, NULL, qStart, qEnd, BLXSTRAND_NONE, NULL, 0, lookupTable, &error);

  /* Parse in additional feature-series info */
  parseLook(msp, look);
  getDesc(msp, line, look);
//  insertFS(msp, series);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
}


/* Parse data from a line of text that contains feature-series data in GFF format */
static void parseFsGff(char *line, 
                       BlxBlastMode blastMode, 
                       GArray* featureLists[], 
                       MSP **lastMsp, 
                       MSP **mspList, 
                       GList **seqList,
                       GList *columnList,
                       GHashTable *lookupTable)
{
  char scorestring[256];
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  char qframe[8] = "(+1)";
  int qStart = UNSET_INT;
  int qEnd = UNSET_INT;
  
  /* <sequencename> <seriesname> <look> <start> <end> <score> <strand> <transframe> [annotation] */
  if (sscanf(line, "%s%s%s%d%d%s%s%s",
	     qName, series, look, &qStart, &qEnd, scorestring, 
	     qframe+1, qframe+2) != 8)
    {
      g_error("Error parsing data, type FS_GFF: \"%s\"\n", line);
    }
  
  int score = UNSET_INT;
  if (!strcmp(scorestring, ".")) 
    {
      score = 100;
    }
  else 
    {
      score = 50.0*atof(scorestring);
    }
  
  BlxStrand qStrand = BLXSTRAND_NONE;
  int qFrame = UNSET_INT;
  getStrandAndFrameFromString(qframe, &qStrand, &qFrame);
  
  /* Create the new MSP */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_FS_SEG, NULL, NULL,
                          score, UNSET_INT, 0, 
                          NULL, qName, qStart, qEnd, qStrand, qFrame, 
                          series, NULL, qStart, qEnd, BLXSTRAND_FORWARD, NULL, 0, lookupTable, &error);
  

  /* Parse additional feature-series information */
  msp->desc = g_strdup(series);
  parseLook(msp, look);
//  insertFS(msp, series);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
}


/* Parse data from the given line, which contains header info for Feature-Series XY-plot data */
static void parseFsXyHeader(char *line, 
                            BlxBlastMode blastMode, 
                            GArray* featureLists[],
                            MSP **lastMsp, 
                            MSP **mspList, 
                            char **seq1, 
                            char *seq1name, 
                            char **seq2,
                            char *seq2name, 
                            BlxParserState *parserState,
                            GList **seqList,
                            GList *columnList,
                            GHashTable *lookupTable)
{
  int i, seqlen;
  
  char series[MAXLINE+1];
  char qName[MAXLINE+1];
  char look[MAXLINE+1];
  
  /* # FS type=XY <sequencename> <seriesname> <look> [annotation] */
  if (sscanf(line+13, "%s%s%s", qName, series, look) != 3)
    {
      g_error("Error parsing data, type XY: \"%s\"\n", line);
    }
  
  if (!seq1name || !seq2name)
    {
      g_error("Sequencenames not provided\n");
    }
  
  if (!strcasecmp(qName, seq1name) || !strcmp(qName, "@1"))
    {
      if (!seq1 || !*seq1)
        {
          g_error("Sequence for %s not provided\n", qName);
        }
        
      seqlen = strlen(*seq1);
    }
  else if (!strcasecmp(qName, seq2name) || !strcmp(qName, "@2"))
    {
      if (!seq2 || !*seq2)
        {
          g_error("Sequence for %s not provided\n", qName);
        }

      seqlen = strlen(*seq2);
    }
  else
    {
      g_error("Invalid sequence name: %s\n", qName);
    }
  
  if (!seqlen)
    { 
      g_error("Sequence for %s is empty\n", qName);
    }

  /* Create an MSP to put the data in */
  GError *error = NULL;
  
  MSP *msp = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_XY_PLOT, NULL, NULL,
                          UNSET_INT, UNSET_INT, 0, 
                          NULL, qName, UNSET_INT, UNSET_INT, BLXSTRAND_FORWARD, 1,
                          series, NULL, UNSET_INT, UNSET_INT, BLXSTRAND_FORWARD, NULL, 0, lookupTable, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  /* Add additional MSP information */
  msp->xy = g_array_sized_new(TRUE, FALSE, sizeof(int), seqlen);
  
  for (i = 0; i < seqlen; i++)
    {
      int *iPtr = (int*)g_malloc(sizeof(int));
      *iPtr = XY_NOT_FILLED;
      g_array_append_val(msp->xy, *iPtr);
    }
  
  msp->fsShape = BLXCURVE_INTERPOLATE; /* default */
  parseLook(msp, look);
  getDesc(msp, line, look);
//  insertFS(msp, series); 
  
  /* Start parsing the actual XY data next */
  *parserState = FS_XY_BODY;
}


/* Parse a line that contains XY-plot data for a feature-series. Data is populated into
 * the given MSP. */
static void parseFsXyData(char *line, MSP *msp)
{
  /* XY coordinates go in to the XY list in the current MSP. */
  int x, y;
  if (sscanf(line, "%d%d", &x, &y) != 2) 
    {
      g_error("Error parsing data file, type FS_XY_BODY: \"%s\"\n", line);
    }
    
  int *xyVal = &g_array_index(msp->xy, int, x - 1);
  *xyVal = y;
}


/* Parse a line that contains feature-series sequence header info */
static void parseFsSeqHeader(char *line, 
                             char **seq1, char *seq1name, char **seq2, char *seq2name,
                             char ***readSeq, int *readSeqLen, int *readSeqMaxLen, BlxParserState *parserState)
{
  /* We have a header for sequence data. Read in the sequence name from the header line. */
  char series[MAXLINE+1];
  char qname[MAXLINE+1];
  
  if (sscanf(line+14, "%s%s", qname, series) != 2) 
    {
      g_error("Error parsing data file, type FS_SEQ_HEADER: \"%s\"\n", line);
    }
  
  /* We need to allocate a string to read the seqeunce data in to. We want to allocate the releavant
   * return pointer (seq1 or seq2) depending on which string was indicated in the header line. */
  if (!strcmp(qname, "@1")) 
    {
      *readSeq = seq1;
      strcpy(seq1name, series);
    }
  else if (!strcmp(qname, "@2")) 
    {
      *readSeq = seq2;
      strcpy(seq2name, series);
    }
  
  *readSeqMaxLen = MAXLINE;
  **readSeq = (char*)g_malloc(*readSeqMaxLen + 1);
  *readSeqLen = 0;
  
  /* Update the parser state so that we proceed to parse the sequence data next */
  *parserState = FS_SEQ_BODY;
}


/* Parse a line that contains a chunk of sequence data. Concatenates the contents of the
 * line onto readSeq, extending readSeq if necessary. readSeqLen is the length of the current
 * contents of readSeq, and readSeqMaxLen is the currently allocated space for readSeq. */
static void parseSeqData(char *line, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, const BlxSeqType seqType)
{
  if (*readSeqLen == UNSET_INT)
    {
      /* Not initialised for reading in sequence; do nothing */
      return;
    }
  
  /* Read in the actual sequence data. It may span several lines, so concatenate each
   * one on to the end of our result. */
  
  /* First, validate that we have valid input */
  char *cp = line;
  for ( ; cp && *cp; ++cp)
    {
      validateIupacChar(cp, seqType);
    }
  
  /* Also reealloc memory if necessary */
  if (*readSeqLen + (int)strlen(line) > *readSeqMaxLen) 
    {
      char *tmp;
      *readSeqMaxLen += MAXLINE + strlen(line);
      tmp = (char*)g_malloc(*readSeqMaxLen + 1);
      strcpy(tmp, **readSeq);
      g_free(**readSeq);
      **readSeq = tmp;
    }
  
  /* Copy in the contents of the line */
  strcpy(**readSeq + *readSeqLen, line);
  
  /* Update the length counter */
  *readSeqLen += strlen(line);
}


/* Utility to call the relevant parser function to parse the current data type */
static void parseBody(char *line, const int lineNum, BlxBlastMode blastMode, const int resFactor, MSP **lastMsp, GString *line_string,
                      char **seq1, char *seq1name, IntRange *seq1Range, char **seq2, char *seq2name, 
                      BlxParserState *parserState, GArray *featureLists[], MSP **mspList, GList **seqList, GList *columnList, GSList *supportedTypes,
                      GSList *styles, char ***readSeq, int *readSeqLen, int *readSeqMaxLen, GKeyFile *keyFile, GHashTable *lookupTable)
{
  //DEBUG_ENTER("parseBody(parserState=%d, line=%d)", *parserState, lineNum);

  GError *tmpError = NULL;
  
  /* Call the relevant function for the current type of data being parsed */
  switch (*parserState)
  {
    case GFF_3_HEADER:
      parseGff3Header(lineNum, lastMsp, mspList, parserState, line_string, seqList, seq1name, seq1Range, &tmpError);
      break;
      
    case GFF_3_BODY:
      parseGff3Body(lineNum, featureLists, lastMsp, mspList, parserState, line_string, 
                    seqList, columnList, supportedTypes, styles, resFactor, keyFile, seq1Range, lookupTable);
      break;

    case SEQBL_BODY: /* fall through */
    case EXBLX_BODY:
      parseEXBLXSEQBL(featureLists, lastMsp, mspList, *parserState, blastMode, line_string, seqList, columnList, lookupTable) ;
      break;
      
    case SEQBL_X_BODY: /* fall through */
    case EXBLX_X_BODY:
      parseEXBLXSEQBLExtended(featureLists, lastMsp, mspList, *parserState, blastMode, line_string, seqList, columnList, lookupTable) ;
      break;

    case FS_HSP_BODY:
      parseFsHsp(line, blastMode, featureLists, lastMsp, mspList, seqList, columnList, lookupTable);
      break;

    case FS_GSP_HEADER:
      parseFsGspHeader(line, blastMode, featureLists, lastMsp, mspList, parserState, seqList, columnList);
      break;

    case FS_GSP_BODY:
      parseFsGspData(line, *lastMsp);
      break;
      
    case FS_GFF_BODY:
      parseFsGff(line, blastMode, featureLists, lastMsp, mspList, seqList, columnList, lookupTable);
      break;

    case FS_SEG_BODY:
      parseFsSeg(line, blastMode, featureLists, lastMsp, mspList, seqList, columnList, lookupTable);
      break;
      
    case FS_XY_HEADER:
      parseFsXyHeader(line, blastMode, featureLists, lastMsp, mspList, seq1, seq1name, seq2, seq2name, parserState, seqList, columnList, lookupTable);
      break;
      
    case FS_XY_BODY:
      parseFsXyData(line, *lastMsp);
      break;

    case FS_SEQ_HEADER:
      parseFsSeqHeader(line, seq1, seq1name, seq2, seq2name, readSeq, readSeqLen, readSeqMaxLen, parserState);
      break;

    case FASTA_SEQ_HEADER:
      parseFastaSeqHeader(line, lineNum, seq1, seq1name, seq1Range, readSeq, readSeqLen, readSeqMaxLen, parserState);
      break;

    case FS_SEQ_BODY: /* fall through */
    case FASTA_SEQ_BODY:
      parseSeqData(line, readSeq, readSeqLen, readSeqMaxLen, BLXSEQ_DNA); /* ref seq is always dna, not peptide */
      break;

    case FASTA_SEQ_IGNORE:
      /* ignore this line of fasta sequence */
      break;
      
    case PARSER_ERROR:
      break;

    default:
      g_warning("Unexpected state '%d' while parsing input file.\n", (int)(*parserState));
      *parserState = PARSER_ERROR;
      break;
  };
  
  if (tmpError)
    {
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
      *parserState = PARSER_ERROR;
    }

  //DEBUG_EXIT("parseBody");
}


/* For old file types, the hacky way of specifying the type of the MSP was by using a negative score
 * for anything that is not a match. This function determines the type from the given score. */
static BlxMspType getMspTypeFromScore(const int score)
{
  BlxMspType result = BLXMSP_INVALID;
  
  if (score >= 0)
    {
      result = BLXMSP_MATCH;
    }
  else if (score == -1)
    {
      result = BLXMSP_CDS;
    }
  else if (score == -2)
    {
      result = BLXMSP_INTRON;
    }
  else if (score == -3)
    {
      result = BLXMSP_VARIATION;
    }

  return result;  
}


/* Extract the strand and frame from the given string of the format "(+1)" and populate
 * the return values, if non-null. */
static void getStrandAndFrameFromString(char *text, BlxStrand *strand, int *frame)
{
  if (frame)
    {
      *frame = convertStringToInt(&text[2]);
    }

  if (strand)
    {
      if (text[1] == '+')
        {
          *strand = BLXSTRAND_FORWARD;
        }
      else if (text[1] == '-')
        {
          *strand = BLXSTRAND_REVERSE;
        }
      else
        {
          *strand = BLXSTRAND_NONE;
        }
    }
}


/* Check if we have a reversed subject and, if so, if this is allowed. Throws an error if not. */
static void checkReversedSubjectAllowed(const MSP *msp, BlxBlastMode blastMode)
{
  if (mspGetMatchStrand(msp) == BLXSTRAND_REVERSE && (blastMode == BLXMODE_BLASTP || blastMode == BLXMODE_BLASTX))
    {
      g_error("Reversed subjects are not allowed in modes blastp or blastx.\n");
    }
}




