/*  File: blxFetch.c
 *  Author: Ed Griffiths, 2008-06-17
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
 * Description: Blixem functions for control of sequence fetching and
 *              display.
 *              
 *              Compiling with -DPFETCH_HTML includes code to issue
 *              pfetch requests to a proxy pfetch server using html.
 *              This requires libs libpfetch and libcurlobj
 *              which in turn require libcurl.
 *----------------------------------------------------------------------------
 */


#include <sys/socket.h> /* for socket(), connect(), send(), and recv() */
#include <netinet/in.h>
#include <arpa/inet.h>  /* for sockaddr_in and inet_addr() */
#include <netdb.h>                                          /* for gethostbyname() */
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>

#include <seqtoolsUtils/utilities.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blixem_.h>

#ifdef PFETCH_HTML 
#include <libpfetch/libpfetch.h>
#endif



/* States for parsing EMBL files */
typedef enum
  {
    PARSING_NEWLINE,             /* parser is at the start of a new line */
    PARSING_ID,                  /* parsing the 2-letter id at the start of a line */
    PARSING_SEQUENCE_HEADER,     /* first line of the SQ section (does not contain sequence data) */
    PARSING_SEQUENCE,            /* main body of the SQ section (contains sequence data on multiple lines) */
    PARSING_IGNORE,              /* we're parsing an area we can ignore */
    PARSING_FINISHED_SEQ,        /* finished parsing the current sequence */
    PARSING_FINISHED,            /* finished parsing all records */
    PARSING_CANCELLED,           /* cancelled by user */

    PARSING_DATA,                /* parsing data for a particular column */

    PARSING_TAG_SEARCH,          /* parsing a known section, looking for a particular tag name */
    PARSING_TAG_NAME,            /* parsing a tag name */
    PARSING_TAG_IGNORE,          /* parsing the tag data but ignore everything until we get to a quoted section */
    PARSING_DATA_QUOTED          /* parsing the data for a tag; same as PARSING_DATA but we're in a quoted section */
  } BlxSeqParserState;



/* Used to draw/update a progress meter, consider as private. */
typedef struct
{
  GtkWidget *top_level ;
  GtkWidget *progress ;
  GtkWidget *label ;
  GdkColor blue_bar_fg, red_bar_fg ;

  gulong widget_destroy_handler_id;

  gboolean cancelled ;
  int seq_total ;
} ProgressBarStruct, *ProgressBar ;


enum {RCVBUFSIZE = 256} ;               /* size of receive buffer for socket fetch */

/* This struct holds general info about a fetch that is in progress */
typedef struct _GeneralFetchData
{
  const BlxFetchMethod* fetchMethod;      /* details about the fetch method */
  GList *columnList;                      /* list of BlxColumnInfo structs */ 
  
  char *buffer;                           /* receive buffer */
  int lenReceived;                        /* number of chars received into the buffer */
  BlxSequence *currentSeq;                /* the current sequence that we're processing */
  GList *currentSeqItem;                  /* the current sequence that we're processing (list item pointer) */
  GString *currentResult;                 /* this stores the current results that have been extracted for the section we're parsing */
  BlxColumnInfo *currentColumn;           /* the column that the current section relates to */
  ProgressBar bar;                        /* shows the progress of the fetch */
  int numRequested;                       /* the total number of sequences requested */
  int numFetched;                         /* the number of sequences fetched so far */
  int numSucceeded;                       /* the number of sequences fetched successfully so far */
  GString *curLine;                       /* this stores the part of the current line that we've parsed so far */
  char sectionId[3];                      /* EMBL lines start with a two-letter identifier, which will be parsed into this string */
  GString *tagName;                       /* In EMBL FT (feature type) sections, we look for tags of the format /tagname="value".
                                           * The current tag name is parsed into this string. */
  gboolean foundEndQuote;                 /* When we're in a quoted section, this is used to flag that we've found a second quote that
                                           * we think is the end of the quoted section. However, a double quote means an escaped quote, so
                                           * if the next char is also a quote, we know we need to include it in the text and carry on parsing. */
  BlxSeqType seqType;
  BlxSeqParserState parserState;
  gboolean status;                        /* gets set to false if there is a problem */
} GeneralFetchData;



#ifdef PFETCH_HTML 
#define PFETCH_READ_SIZE 80     /* about a line */
#define PFETCH_FAILED_PREFIX "PFetch failed:"


typedef struct
{
  GtkWidget *blxWindow;
  GtkWidget *dialog;
  GtkTextBuffer *text_buffer;
  char *title;
  const BlxSequence *blxSeq;

  gulong widget_destroy_handler_id;
  PFetchHandle pfetch;
  gboolean got_response;
  int attempt;
  const BlxFetchMethod *fetchMethod;
} PFetchDataStruct, *PFetchData;


/* this holds info about an http fetch that is in progress */
typedef struct
{
  gboolean connection_closed;                               /* Gets set to true when pfetch connection is closed */
  char *err_txt ;                                           /* if !status then err message here. */
  gboolean stats ;                                          /* TRUE means record and output stats. */
  int min_bytes, max_bytes, total_bytes, total_reads ;      /* Stats. */
  PFetchHandle pfetch ;
  GList *seqList ;                                          /* List of sequences to fetch */

  GeneralFetchData fetchData ;                              /* general info about a fetch in progress*/
} PFetchSequenceStruct, *PFetchSequence ;

#endif


/* Local function declarations */
#ifdef PFETCH_HTML 
static PFetchStatus pfetch_reader_func(PFetchHandle *handle,
                                       char         *text,
                                       guint        *actual_read,
                                       GError       *error,
                                       gpointer      user_data) ;
static void handle_dialog_close(GtkWidget *dialog, gpointer user_data);
static PFetchStatus pfetch_closed_func(gpointer user_data) ;


static PFetchStatus                sequence_pfetch_reader(PFetchHandle *handle, char *text, guint *actual_read, GError *error, gpointer user_data) ;
static PFetchStatus                sequence_pfetch_closed(PFetchHandle *handle, gpointer user_data) ;
static void                        sequence_dialog_closed(GtkWidget *dialog, gpointer user_data) ;
static gboolean                    parsePfetchHtmlBuffer(const BlxFetchMethod* const fetchMethod, char *read_text, int length, PFetchSequence fetch_data) ;

static void                        httpFetchSequence(const BlxSequence *blxSeq, const BlxFetchMethod* const fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer **text_buffer);
#endif

static int                         socketConstruct(const char *ipAddress, int port, gboolean External, GError **error) ;
static void                        socketSend(int sock, const char *text, GError **error) ;

static ProgressBar                 makeProgressBar(int seq_total) ;
static void                        updateProgressBar(ProgressBar bar, const char *sequence, int numFetched, gboolean fetch_ok) ;
static gboolean                    isCancelledProgressBar(ProgressBar bar) ;
static void                        destroyProgressBar(ProgressBar bar) ;
static void                        destroyProgressCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */
static void                        cancelCB(GtkWidget *widget, gpointer cb_data) ; /* internal to progress bar. */

static void                        readConfigFile(GKeyFile *key_file, CommandLineOptions *options, GError **error) ;

static void                        socketFetchSequence(const BlxSequence *blxSeq, const BlxFetchMethod* const fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer **text_buffer);
static void                        commandFetchSequence(const BlxSequence *blxSeq, const BlxFetchMethod* const fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer **text_buffer);
static void                        internalFetchSequence(const BlxSequence *blxSeq, const BlxFetchMethod* const fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer **text_buffer);
static void                        wwwFetchSequence(const BlxSequence *blxSeq, const BlxFetchMethod* const fetchMethod, const gboolean displayResults, const int attempt, GtkWidget *blxWindow);

static void                        checkFetchMethodExecutable(const BlxFetchMethod* const fetchMethod, GError **error);

/* Pfetch local functions */
static void                        appendCharToString(const char curChar, GString *result);
static void                        appendCharToQuotedString(const char curChar, gboolean *foundEndQuote, GString *result);

static void                        socketFetchInit(const BlxFetchMethod* const fetchMethod, GList *seqsToFetch, gboolean External, int *sock, GError **error);
static void                        checkProgressBar(ProgressBar bar, BlxSeqParserState *parserState, gboolean *status);

static int                         socketFetchReceiveBuffer(GeneralFetchData *fetchData, const int bufferSize, const int sock);

static void                        pfetchGetNextSequence(GeneralFetchData *fetchData, const gboolean pfetch_ok);

static void                        parseRawSequenceBuffer(GeneralFetchData *fetchData, GError **error);

static void                        pfetchGetParserStateFromId(GeneralFetchData *fetchData);

static void                        parseEmblBuffer(GeneralFetchData *fetchData, GError **error);

static gboolean                    pfetchFinishSequence(GeneralFetchData *fetchData);

static void                        pfetchGetParserStateFromTagName(GeneralFetchData *fetchData);



/* global configuration object for blixem. */
static GKeyFile *blx_config_G = NULL ;


/* Utility to convert a fetch mode enum into a string (used in the config file) */
const char *fetchModeStr(const BlxFetchMode fetchMode)
{
  /* Values must be in the same order as BlxFetchMode */
  static const gchar* fetchModeNames[] = 
    {
#ifdef PFETCH_HTML 
      "http",
      "pipe",
#endif
      "socket",
      "www",
      "sqlite",
      "command",
      "internal",
      "none",
      NULL
    }; 
  
  g_assert(g_strv_length((gchar**)fetchModeNames) == BLXFETCH_NUM_MODES);

  const char *result = NULL;
  
  if (fetchMode < BLXFETCH_NUM_MODES)
    result = fetchModeNames[fetchMode];
  
  return result;
}


/* Utility to convert a fetch-method output type into a string (used in the config file) */
const char *outputTypeStr(const BlxFetchOutputType outputType)
{
  /* Values must be in the same order as BlxFetchMode */
  static const gchar* outputNames[] = 
    {
      "<invalid>",
      "raw",
      "fasta",
      "embl",
      "list",
      "gff",
      NULL
    }; 
  
  g_assert(g_strv_length((gchar**)outputNames) == BLXFETCH_NUM_OUTPUT_TYPES);

  const char *result = NULL;
  
  if (outputType < BLXFETCH_NUM_OUTPUT_TYPES)
    result = outputNames[outputType];
  
  return result;
}


/* Get the fetch-method struct containing the details for the 
 * given fetch method */
BlxFetchMethod* getFetchMethodDetails(GQuark fetchMethodQuark, GHashTable *fetchMethods)
{
  BlxFetchMethod *result = (BlxFetchMethod*)g_hash_table_lookup(fetchMethods, GINT_TO_POINTER(fetchMethodQuark));
  return result;
}



/* Fetch the given sequence and optionally display the results. 
 * dialog and text_buffer are only used when recursing via httpFetchSequence;
 * they should be passed as NULL in all other cases. */
void fetchSequence(const BlxSequence *blxSeq, 
                   const gboolean displayResults,
                   const int attempt,
                   GtkWidget *blxWindow,
                   GtkWidget *dialog, 
                   GtkTextBuffer **text_buffer)
{
  g_assert(blxSeq);
  
  /* Look up the fetch method for this sequence */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GQuark fetchMethodQuark = blxSequenceGetFetchMethod(blxSeq, FALSE, FALSE, attempt, bc->userFetchDefault);
  const BlxFetchMethod* const fetchMethod = getFetchMethodDetails(fetchMethodQuark, bc->fetchMethods);

  if (!fetchMethod)
    {
      /* If this is the first attempt then we should have a fetch method; 
       * therefore give a warning if no fetch method was found */
      if (attempt == 0 && !fetchMethodQuark)
        g_warning("No fetch method specified for sequence '%s'\n", blxSequenceGetName(blxSeq));
      else if (fetchMethodQuark)
        g_warning("Error fetching sequence '%s'; could not find details for fetch method '%s'\n", blxSequenceGetName(blxSeq), g_quark_to_string(fetchMethodQuark));

      return;
    }

  if (fetchMethod->mode == BLXFETCH_MODE_NONE && attempt == 0)
    {
      g_message("Fetch method for '%s' is '%s'\n", blxSequenceGetName(blxSeq), fetchModeStr(BLXFETCH_MODE_NONE));
      return;
    }
  
  g_message("Fetching '%s' using method '%s' (attempt %d)\n", blxSequenceGetName(blxSeq), g_quark_to_string(fetchMethodQuark), attempt + 1);

  
  if (fetchMethod->mode == BLXFETCH_MODE_SOCKET)
    {
      socketFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, dialog, text_buffer);
    }
#ifdef PFETCH_HTML 
  else if (fetchMethod->mode == BLXFETCH_MODE_HTTP || fetchMethod->mode == BLXFETCH_MODE_PIPE)
    {
      httpFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, dialog, text_buffer);
    }
#endif
  else if (fetchMethod->mode == BLXFETCH_MODE_COMMAND)
    {
      commandFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, dialog, text_buffer);
    }
  else if (fetchMethod->mode == BLXFETCH_MODE_WWW)
    {
      wwwFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow);
    }
  else if (fetchMethod->mode == BLXFETCH_MODE_SQLITE)
    {
      sqliteFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow);
    }
  else if (fetchMethod->mode == BLXFETCH_MODE_INTERNAL)
    {
      internalFetchSequence(blxSeq, fetchMethod, displayResults, attempt, blxWindow, dialog, text_buffer);
    }
  else
    {
      /* Invalid fetch method. Try again with the next fetch method, if one is specified */
      g_warning("Unknown fetch method: %s\n", g_quark_to_string(fetchMethod->name));
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, dialog, text_buffer);
    }
}


/* Process a single substitution character. Returns false if it's not
 * a known substitution char. */
static gboolean doFetchStringSubstitutionChar(const char substitution_char, 
                                          MatchSequenceData *match_data,
                                          GString *result, 
                                          GString *errorMsg)
{
  gboolean ok = TRUE;
  
  switch (substitution_char)
    {
    case 'p': 
      g_string_append(result, g_get_prgname());
      break;
    case 'h': 
      g_string_append(result, g_get_host_name());
      break;
    case 'u': 
      g_string_append(result, g_get_user_name());
      break;
    case 'm': 
      if (match_data->match_name)
        g_string_append(result, match_data->match_name);
      break;
    case 'r': 
      if (match_data->ref_name)
        g_string_append(result, match_data->ref_name);
      break;
    case 's': 
      g_string_append_printf(result, "%d", match_data->match_start);
      break;
    case 'e': 
      g_string_append_printf(result, "%d", match_data->match_end);
      break;
    case 'd':
      if (match_data->dataset)
        g_string_append(result, match_data->dataset);
      break;
    case 'S':
      if (match_data->source)
        g_string_append(result, match_data->source);
      break;
    case 'f':
      if (match_data->filename)
        g_string_append(result, match_data->filename);
      break;
    case 'g':
      g_string_append_printf(result, "%d", SEQTOOLS_GFF_VERSION);
      break;
    case '%':
      break;
    default:
      ok = FALSE;
      
      if (!errorMsg)
        errorMsg = g_string_new("");
      
      g_string_append_printf(errorMsg, "  Unknown substitution character '%%%c'\n", substitution_char);
      break;
    };

  return ok;
}


/* Process a single substitution keyword. Its format should be:
 *    %(<keyword>) 
 * where <keyword> is a column name e.g. Source or a key value to 
 * look up in the source stanza in the config file. */
static int doFetchStringSubstitutionKeyword(const char* input_string, 
                                            MatchSequenceData *match_data,
                                            GString *result, 
                                            GString *errorMsg)
{
  int len = 0;
  gboolean ok = FALSE;

  /* The input should start with '('. Search for the closing ')' and process
   * the text in between */
  if (input_string && *input_string == '(')
    {
      char *cp = strchr(input_string, ')');
      
      if (cp)
        {
          len = cp - input_string - 1;
          
          if (len > 0)
            {
              char *key = g_strdup_printf("%s", input_string + 1);
              key[len] = '\0';
              
              /*! \todo Ideally we'd allow the keyword to be a known column name in which case
               * we'd need to add a check here. Needs a bit of code reorganisation to 
               * create the columnList before we parse & fetch sequences - currently the columnList
               * is only created when the window is created, i.e. after parsing & fetching is complete. */
              
              /* See if it's a field in this match's source stanza */
              if (match_data->source)
                {
                  GKeyFile *key_file = blxGetConfig();
                  
                  if (key_file)
                    {
                      char *value = g_key_file_get_string(key_file, match_data->source, key, NULL);
                      
                      if (value)
                        {
                          g_string_append(result, value);
                          ok = TRUE;
                          
                          g_free(value);
                        }
                    }
                }
              
              g_free(key);
            }
        }        
    }

  if (!ok)
    {
      if (!errorMsg)
        errorMsg = g_string_new("");
     
      g_string_append_printf(errorMsg, "  Failed to process substitution string '%s'\n", input_string);
    }
  
  return len;
}


/* Process the fetch string, substituting the special substitution characters 
 * for the actual data */
static GString* doFetchStringSubstitutions(const char *command,
                                           MatchSequenceData *match_data,
                                           GError **error)
{
  GString *result = g_string_new("");
  
  /* Loop through the command and substitute any special chars with the relevant info  */
  GString *errorMsg = NULL;
  const char *c = command;
      
  while (c && *c)
    {
      /* If it's preceded by the special char, substitute it for the real value */
      if (*c == '%')
        {
          /* Move to the next char, which should tell us what type of substitution to make */
          ++c;
          
          if (c && *c)
            {
              switch (*c) {
              case '%':
                {
                  g_string_append_c(result, *c);
                  ++c;
                  break;
                }
              case '(':
                {
                  int len = doFetchStringSubstitutionKeyword(c, match_data, result, errorMsg);
                  if (len > 0)
                    c += len + 2; /* progress past the len of the keyword plus the two brackets */
                  else 
                    ++c; /* failed; just increment past the current char */
                  break;
                }
              default:
                {
                  doFetchStringSubstitutionChar(*c, match_data, result, errorMsg);
                  ++c; /* Progress to next char even if failed */
                  break;
                }
              };
            }
        }
      else
        {
          /* Normal char; just append to the result */
          g_string_append_c(result, *c);
          ++c;
        }
    }

  if (errorMsg)
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_SUBSTITUTION, "%s", errorMsg->str);
      prefixError(*error, "Error constructing fetch command:\n");
      g_string_free(errorMsg, TRUE);
    }

  return result;
}


static GString* doGetFetchArgs(const BlxFetchMethod* const fetchMethod,
                               MatchSequenceData *match_data,
                               GError **error)
{
  GString *result = doFetchStringSubstitutions(fetchMethod->args, match_data, error);

  return result;
}


/* Return true if the given fetch method uses http */
static gboolean fetchMethodUsesHttp(const BlxFetchMethod* const fetchMethod)
{
  return (fetchMethod->mode == BLXFETCH_MODE_WWW
#ifdef PFETCH_HTML
          || fetchMethod->mode == BLXFETCH_MODE_HTTP
#endif
          );
}


/* Compile a command path and arguments into a single string */
GString* doGetFetchCommand(const BlxFetchMethod* const fetchMethod,
                           MatchSequenceData *match_data,
                           GError **error)
{
  GString *result = NULL;
  GError *tmpError = NULL;

  if (fetchMethod)
    {
      /* If an executable is required, check that it exists */
      checkFetchMethodExecutable(fetchMethod, &tmpError);

      if (!tmpError)
        {
          /* Compile the command and args into a single string */
          char *command = NULL;
      
          /* For http methods, append the args (i.e. the request) after a '?'.
           * For other methods, append the args after a space. */
          if (fetchMethod->location && fetchMethod->args && fetchMethodUsesHttp(fetchMethod))
            command = g_strdup_printf("%s?%s", fetchMethod->location, fetchMethod->args);
          else if (fetchMethod->location && fetchMethod->args)
            command = g_strdup_printf("%s %s", fetchMethod->location, fetchMethod->args);
          else if (fetchMethod->location)
            command = g_strdup(fetchMethod->location);
          else if (fetchMethod->args)
            command = g_strdup(fetchMethod->args);
          else
            return result;
          
          result = doFetchStringSubstitutions(command, match_data, error);
          
          /* Clean up */
          g_free(command);
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
  
  return result;
}


GString* getFetchArgsMultiple(const BlxFetchMethod* const fetchMethod, GList *seqsToFetch, GError **error)
{
  /* Build up a string containing all the sequence names. */
  GString *seq_string = g_string_sized_new(1000) ;
  GList *seqItem = seqsToFetch;

  const char *separator = fetchMethod->separator ? fetchMethod->separator : " ";
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      g_string_append_printf(seq_string, "%s%s", blxSequenceGetName(blxSeq), separator);
    }

  MatchSequenceData match_data = {seq_string->str, NULL, 0, 0, NULL, NULL, NULL};
  GString *result = doGetFetchArgs(fetchMethod, &match_data, error);

  g_string_free(seq_string, TRUE);
  
  return result;
}


/* Get the command to call for the given fetch method. Parses
 * the command and arguments and substitutes the following
 * characters with values from the given sequence/msp. Note that 
 * either blxSeq or msp must be given, but not both. '%%' is used
 * to represent a normal '%' character.
 * 
 *   %p:      program name
 *   %h:      host name
 *   %u:      user name
 *   %m:      match sequence name
 *   %r:      reference sequence name
 *   %s:      start coord
 *   %e:      end coord
 *   %d:      dataset
 *   %S:      feature source
 *   %f:      file name
 *   %g:      supported GFF version
 * 
 * Returns the command and args compiled into a single string.
 * The caller must free the result with g_string_free.
 * Returns an empty string if the command/args are empty.
*/
GString* getFetchCommand(const BlxFetchMethod* const fetchMethod,
                         const BlxSequence *blxSeq_in,
                         const MSP* const msp,
                         const char *refSeqName,
                         const int refSeqOffset,
                         const IntRange* const refSeqRange,
                         const char *dataset,
                         GError **error)
{
  const BlxSequence *blxSeq = blxSeq_in;

  if (!blxSeq && msp)
    blxSeq = msp->sSequence;

  /* Extract info about the sequence / feature */
  const char *name = blxSequenceGetName(blxSeq);
  if (!name) name = mspGetSName(msp);
  if (!name) name = "";
  
  const char *source = blxSeq ? blxSequenceGetSource(blxSeq) : NULL;
  const char *filename = msp && msp->filename ? g_quark_to_string(msp->filename) : NULL;
  int startCoord = msp ? mspGetQStart(msp) : blxSequenceGetStart(blxSeq, blxSeq->strand);
  int endCoord = msp ? mspGetQEnd(msp) : blxSequenceGetEnd(blxSeq, blxSeq->strand);
  startCoord += refSeqOffset;
  endCoord += refSeqOffset;
  boundsLimitValue(&startCoord, refSeqRange);
  boundsLimitValue(&endCoord, refSeqRange);

  /* Do the substitutions */
  MatchSequenceData match_data = {name, refSeqName, startCoord, endCoord, dataset, source, filename};
  GString *result = doGetFetchCommand(fetchMethod, &match_data, error);
  
  return result;
}


/* Get the arguments string for the given fetch method (replacing
 * any substitution characters with the correct values for the given
 * sequence) */
GString* getFetchArgs(const BlxFetchMethod* const fetchMethod,
                      const BlxSequence *blxSeq,
                      const MSP* const msp,
                      const char *refSeqName,
                      const int refSeqOffset,
                      const IntRange* const refSeqRange,
                      const char *dataset,
                      GError **error)
{
  g_assert(blxSeq || msp);

  /* Extract info about the sequence / feature */

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  const char *name = blxSequenceGetName(blxSeq);
  if (!name) name = mspGetSName(msp);
  if (!name) name = "";
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
  const char *name ;

  if (!(name = mspGetSNameOrig(msp)))
    {
      name = blxSequenceGetName(blxSeq);
      if (!name) name = mspGetSName(msp);
      if (!name) name = "";
    }
  
  const char *source = blxSeq ? blxSequenceGetSource(blxSeq) : NULL;
  const char *filename = msp && msp->filename ? g_quark_to_string(msp->filename) : NULL;
  int startCoord = blxSeq ? blxSequenceGetStart(blxSeq, blxSeq->strand) : mspGetQStart(msp);
  int endCoord = blxSeq ? blxSequenceGetEnd(blxSeq, blxSeq->strand) : mspGetQEnd(msp);
  startCoord += refSeqOffset;
  endCoord += refSeqOffset;
  boundsLimitValue(&startCoord, refSeqRange);
  boundsLimitValue(&endCoord, refSeqRange);

  /* Do the substitutions */
  MatchSequenceData match_data = {name, refSeqName, startCoord, endCoord, dataset, source, filename};
  GString *result = doGetFetchArgs(fetchMethod, &match_data, error);
  
  return result;
}


/* Return true if the given string matches any in the given array
 * of quarks (not case sensitive). */
static gboolean stringInArray(const char *str, GArray *array)
{
  gboolean found = FALSE;

  if (array)
    {
      const int len1 = strlen(str);
      int i = 0;
      
      for ( ; !found && i < (int)array->len; ++i)
        {
          GQuark curQuark = g_array_index(array, GQuark, i);
          const char *curStr = g_quark_to_string(curQuark);
          
          int len = min((int)strlen(curStr), len1);
          
          if (strncasecmp(curStr, str, len) == 0)
            found = TRUE;
        }
    }
  
  return found;
}


/* Check that the given fetch method is not null; set the error if not */
static void checkFetchMethodNonNull(const BlxFetchMethod* const fetchMethod, GError **error)
{
  if (!fetchMethod)
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_FETCH_METHOD, 
                  "Program error: fetch method is null\n");
    }
}


/* Check that the given fetch method's executable is in the path; set 
 * the error if not */
static void checkFetchMethodExecutable(const BlxFetchMethod* const fetchMethod, GError **error)
{
  /* Only command- and socket-fetch require an executable */
  if (fetchMethod &&
      (fetchMethod->mode == BLXFETCH_MODE_COMMAND &&
       fetchMethod->mode == BLXFETCH_MODE_SOCKET))
    {
      if (!fetchMethod->location || !g_find_program_in_path(fetchMethod->location))
        {
          g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_NO_EXE,
                      "[%s]: Executable '%s' not found in path: %s\n", 
                      g_quark_to_string(fetchMethod->name), fetchMethod->location, getenv("PATH"));
        }
    }
}


/* Use the www-fetch method to fetch an entry and optionally display
 * the results in a dialog.
 * Opens a browser to display the results. Does nothing if 
 * not displaying results! */
static void wwwFetchSequence(const BlxSequence *blxSeq,
                             const BlxFetchMethod* const fetchMethod, 
                             const gboolean displayResults, 
                             const int attempt,
                             GtkWidget *blxWindow)
{
  if (displayResults)
    {
      BlxViewContext *bc = blxWindowGetContext(blxWindow);

      GError *error = NULL;
      
      GString *url = getFetchCommand(fetchMethod, 
                                     blxSeq, 
                                     NULL, 
                                     bc->refSeqName, 
                                     bc->refSeqOffset,
                                     &bc->refSeqRange,
                                     bc->dataset,
                                     &error);


      if (!error)
        {
          seqtoolsLaunchWebBrowser(url->str, &error);
        }
      
      if (url)
        {
          g_string_free(url, TRUE);
        }

      /* If failed, re-try with the next-preferred fetch method, if there is one */
      if (error)
        {
          fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, NULL, NULL);
          g_error_free(error);
        }
    }
}


/* Use the command-fetch method to fetch an entry and optionally display
 * the results in a dialog. */
static void commandFetchSequence(const BlxSequence *blxSeq,
                                 const BlxFetchMethod* fetchMethod, 
                                 const gboolean displayResults, 
                                 const int attempt,
                                 GtkWidget *blxWindow,
                                 GtkWidget *dialog,
                                 GtkTextBuffer **text_buffer)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GError *error = NULL;
  GString *command = NULL;
  GString *resultText = NULL;
  
  if (!error)
    checkFetchMethodNonNull(fetchMethod, &error);

  if (!error)
    checkFetchMethodExecutable(fetchMethod, &error);

  if (!error)
    command = getFetchCommand(fetchMethod, blxSeq, NULL, bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, bc->dataset, &error);

  if (!error && command)
    resultText = getExternalCommandOutput(command->str, &error);

  reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);

  if (resultText && resultText->str)
    {
      if (displayResults && !error)
        {
          char *title = g_strdup_printf("%s%s", blxGetTitlePrefix(bc), command->str);
          displayFetchResults(title, resultText->str, blxWindow, dialog, text_buffer);
          g_free(title);
        }

      g_string_free(resultText, TRUE);
    }
  else
    {
      /* Try again with the next-preferred fetch method, if there is one */
      if (resultText)
        g_string_free(resultText, TRUE);
      
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, NULL, NULL);
    }
  
  if (command)
    g_string_free(command, TRUE);
}


/* This "fetch" method doesn't really fetch the sequence: it just
 * returns the internally-stored sequence */
static void internalFetchSequence(const BlxSequence *blxSeq,
                                  const BlxFetchMethod* const fetchMethod, 
                                  const gboolean displayResults, 
                                  const int attempt,
                                  GtkWidget *blxWindow,
                                  GtkWidget *dialog,
                                  GtkTextBuffer **text_buffer)
{
  const char *seq = blxSequenceGetSequence(blxSeq);
  const char *seqName = blxSequenceGetName(blxSeq);

  if (seq)
    {
      char *result = g_strdup_printf(">%s\n%s", seqName ? seqName : "", seq);

      if (displayResults)
        {
          BlxViewContext *bc = blxWindowGetContext(blxWindow);
          char *title = g_strdup_printf("%s%s", blxGetTitlePrefix(bc), seqName ? seqName : "");
          displayFetchResults(title, result, blxWindow, dialog, text_buffer);
          g_free(title);
        }

      g_free(result);
    }
  else
    {
      g_warning("No sequence data found for '%s'\n", seqName ? seqName : "");
      
      /* Try again with the next-preferred fetch method, if there is one */
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, dialog, text_buffer);
    }
}


/* Use the given socket-fetch method to fetch an entry and optionally display the results. */
static void socketFetchSequence(const BlxSequence *blxSeq, 
                                const BlxFetchMethod* const fetchMethod, 
                                const gboolean displayResults, 
                                const int attempt,
                                GtkWidget *blxWindow,
                                GtkWidget *dialog,
                                GtkTextBuffer **text_buffer)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GError *error = NULL;
  GString *resultText = NULL;
  GString *command = NULL;

  if (!error)
    checkFetchMethodNonNull(fetchMethod, &error);

  if (!error)
    checkFetchMethodExecutable(fetchMethod, &error);

  if (!error)
    command = getFetchCommand(fetchMethod, blxSeq, NULL, bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, bc->dataset, &error);  

  if (!error && command)
    resultText = getExternalCommandOutput(command->str, &error);
  
  reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);

  if (resultText && resultText->len && !stringInArray(resultText->str, fetchMethod->errors))  /* Success */
    {
      if (displayResults)
        {
          char *title = g_strdup_printf("%s%s", blxGetTitlePrefix(bc), command->str);
          displayFetchResults(title, resultText->str, blxWindow, dialog, text_buffer);
          g_free(title);
        }

      g_string_free(resultText, TRUE);
    }
  else   /* Failed */
    {
      if (resultText)
        g_string_free(resultText, TRUE);

      /* Try again with the next fetch method, if there is one set */
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, NULL, NULL);
    }

  if (command)
    g_string_free(command, TRUE);
}


#ifdef PFETCH_HTML 
/* Fetch the given list of sequences from an http proxy server. This enables
 * blixem to be run and get sequences from anywhere that can see the http 
 * proxy server. */
static gboolean httpFetchList(GList *seqsToFetch, 
                              const BlxFetchMethod* const fetchMethod,
                              GList *seqList, 
                              GList *columnList,
                              const BlxSeqType seqType,
                              GError **error)
{
  gboolean status = FALSE ;
  gboolean debug_pfetch = FALSE ;

  GType pfetch_type = PFETCH_TYPE_HTTP_HANDLE ;

  if (fetchMethod->mode == BLXFETCH_MODE_PIPE)
    pfetch_type = PFETCH_TYPE_PIPE_HANDLE ;

  PFetchSequenceStruct fetch_data = {FALSE} ;

  fetch_data.connection_closed = FALSE;
  fetch_data.err_txt = NULL;
  fetch_data.stats = FALSE ;
  fetch_data.min_bytes = INT_MAX ;
  fetch_data.max_bytes = 0 ;
  fetch_data.total_bytes = 0 ;
  fetch_data.total_reads = 0 ;
  fetch_data.seqList = seqsToFetch;
  fetch_data.pfetch = PFetchHandleNew(pfetch_type);

  fetch_data.fetchData.fetchMethod = fetchMethod;
  fetch_data.fetchData.columnList = columnList;
  fetch_data.fetchData.currentColumn = NULL;
  fetch_data.fetchData.currentSeq = seqsToFetch ? (BlxSequence*)(seqsToFetch->data) : NULL;
  fetch_data.fetchData.currentSeqItem = seqsToFetch;
  fetch_data.fetchData.currentResult = g_string_new("");
  fetch_data.fetchData.numRequested = g_list_length(seqsToFetch);
  fetch_data.fetchData.numFetched = 0;
  fetch_data.fetchData.numSucceeded = 0;
  fetch_data.fetchData.bar = makeProgressBar(fetch_data.fetchData.numRequested) ;
  fetch_data.fetchData.curLine = g_string_new("");
  fetch_data.fetchData.sectionId[0] = ' ';
  fetch_data.fetchData.sectionId[1] = ' ';
  fetch_data.fetchData.sectionId[2] = '\0';
  fetch_data.fetchData.tagName = g_string_new("");
  fetch_data.fetchData.foundEndQuote = FALSE;
  fetch_data.fetchData.seqType = seqType ;
  fetch_data.fetchData.parserState = PARSING_NEWLINE;
  fetch_data.fetchData.status = TRUE ;
  
  g_signal_connect(G_OBJECT(fetch_data.fetchData.bar->top_level), "destroy",
               G_CALLBACK(sequence_dialog_closed), &fetch_data) ;
  
  if (PFETCH_IS_HTTP_HANDLE(fetch_data.pfetch))
    {
      PFetchHandleSettings(fetch_data.pfetch, 
                           "port",       fetchMethod->port,
                           "debug",      debug_pfetch,
                           "pfetch",     fetchMethod->location,
                           "cookie-jar", fetchMethod->cookie_jar,
                           NULL);
    }
  else
    {
      PFetchHandleSettings(fetch_data.pfetch, 
                           "pfetch",     fetchMethod->location,
                           NULL);
    }
  
  g_signal_connect(G_OBJECT(fetch_data.pfetch), "reader", G_CALLBACK(sequence_pfetch_reader), &fetch_data) ;

  g_signal_connect(G_OBJECT(fetch_data.pfetch), "closed", G_CALLBACK(sequence_pfetch_closed), &fetch_data) ;

  GString *request = getFetchArgsMultiple(fetchMethod, seqsToFetch, error);
  
  /* Set up pfetch/curl connection routines, this is non-blocking so if connection
   * is successful we block using our own flag. */
  if (PFetchHandleFetch(fetch_data.pfetch, request->str) == PFETCH_STATUS_OK)
    {
      status = TRUE ;

      while (!(fetch_data.connection_closed) && status && fetch_data.fetchData.parserState != PARSING_CANCELLED)
        {
          checkProgressBar(fetch_data.fetchData.bar, &fetch_data.fetchData.parserState, &status);
          gtk_main_iteration() ;
        }

      status = fetch_data.fetchData.status ;
      
      if (!status)
        {
          if (fetch_data.fetchData.parserState != PARSING_CANCELLED)
            {
              g_critical("Sequence fetch from http server failed: %s\n", fetch_data.err_txt) ;
              
              if (fetch_data.err_txt)
                { 
                  g_free(fetch_data.err_txt) ;
                }
            }
        }
    }
  else
    {
      status = FALSE ;
    }
  
  destroyProgressBar(fetch_data.fetchData.bar) ;
  g_string_free(fetch_data.fetchData.currentResult, TRUE);
  g_string_free(fetch_data.fetchData.tagName, TRUE);
  g_string_free(fetch_data.fetchData.curLine, TRUE);

  if (request)
    g_string_free(request, FALSE);
  
  return status ;
}


/* Use the http proxy to pfetch an entry */
/* Note that this uses a callback to update the display 
 * window; it cannot return the result immediately and the 
 * code is not currently structured to allow the callback to
 * do anything other than update the display window, so if
 * we're just requesting the sequence, this currently just 
 * returns null. */
static void httpFetchSequence(const BlxSequence *blxSeq,
                              const BlxFetchMethod* const fetchMethod,
                              const gboolean displayResults,
                              const int attempt,
                              GtkWidget *blxWindow, 
                              GtkWidget *dialog, 
                              GtkTextBuffer **text_buffer)
{
  if (!displayResults)
    {
      g_warning("Program error: http-fetch expected to display results but displayResults is false.\n");
      return;
    }

  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  gboolean debug_pfetch = FALSE ;
  PFetchData pfetch_data ;
  GError *tmpError = NULL;
  GString *command = NULL;
  GString *request = NULL;

  if (fetchMethod->location == NULL)
    g_set_error(&tmpError, BLX_ERROR, 1, "%s", "Failed to obtain preferences specifying how to pfetch.\n");
  
  if (!tmpError)
    {
      GType pfetch_type = PFETCH_TYPE_HTTP_HANDLE ;

      if (fetchMethod->mode == BLXFETCH_MODE_PIPE)
        pfetch_type = PFETCH_TYPE_PIPE_HANDLE ;
        
      pfetch_data = g_new0(PFetchDataStruct, 1);

      pfetch_data->pfetch = PFetchHandleNew(pfetch_type);
      
      pfetch_data->blxWindow = blxWindow;
      pfetch_data->blxSeq = blxSeq;
      pfetch_data->attempt = attempt;
      pfetch_data->fetchMethod = fetchMethod;
      
      command = getFetchCommand(fetchMethod, blxSeq, NULL, bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, bc->dataset, &tmpError);
    }
  
  if (!tmpError)
    {
      request = getFetchArgs(fetchMethod, blxSeq, NULL, 
                             bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange, 
                             bc->dataset, &tmpError);
    }
  
  if (tmpError)
    {
      /* Couldn't initiate the fetch; try again with a different fetch method */
      g_free(pfetch_data->pfetch);
      g_free(pfetch_data);
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
      fetchSequence(blxSeq, TRUE, attempt + 1, blxWindow, dialog, text_buffer);
    }
  else
    {
      pfetch_data->title = g_strdup_printf("%s%s", blxGetTitlePrefix(bc), command->str);
      
      if (pfetch_data->title)
        {
          if (dialog && text_buffer && *text_buffer)
            {
              gtk_window_set_title(GTK_WINDOW(dialog), pfetch_data->title);
              pfetch_data->dialog = dialog;
              pfetch_data->text_buffer = *text_buffer;
            }
          else
            {
              pfetch_data->dialog = displayFetchResults(pfetch_data->title, "pfetching...\n", blxWindow, dialog, &pfetch_data->text_buffer);
            }

          pfetch_data->widget_destroy_handler_id = 
            g_signal_connect(G_OBJECT(pfetch_data->dialog), "destroy", 
                             G_CALLBACK(handle_dialog_close), pfetch_data); 
        }
      
      if (PFETCH_IS_HTTP_HANDLE(pfetch_data->pfetch))
        {
            PFetchHandleSettings(pfetch_data->pfetch, 
                                 "port",       fetchMethod->port,
                                 "debug",      debug_pfetch,
                                 "pfetch",     fetchMethod->location,
                                 "cookie-jar", fetchMethod->cookie_jar,
                                 NULL);
        }
      else
        {
            PFetchHandleSettings(pfetch_data->pfetch, 
                                 "pfetch",     fetchMethod->location,
                                 NULL);
        }
      
      g_signal_connect(G_OBJECT(pfetch_data->pfetch), "reader", G_CALLBACK(pfetch_reader_func), pfetch_data);
      g_signal_connect(G_OBJECT(pfetch_data->pfetch), "closed", G_CALLBACK(pfetch_closed_func), pfetch_data);

      PFetchHandleFetch(pfetch_data->pfetch, request->str) ;
      
    }

  if (request)
    {
      g_string_free(request, FALSE);
    }
}


#endif


/* Fetch a list of sequences using sockets.
 * 
 * adapted from Tony Cox's code pfetch.c
 *
 *  - this version incorporates a progress monitor as a window,
 *    much easier for user to control + has a cancel button.
 *  - can be called after the fasta sequence data is already populated:
 *    in that case it will ignore the sequence data and just populate
 *    the additional data.
 *  - sequence data will also be ignored for sequences that do not 
 *    require sequence data
 */
gboolean socketFetchList(GList *seqsToFetch, 
                         const BlxFetchMethod* const fetchMethod,
                         GList *seqList, 
                         GList *columnList,
                         gboolean External, 
                         const BlxSeqType seqType, 
                         GError **error)
{
  /* Initialise and send the requests */
  int sock;
  GError *tmpError = NULL;

  socketFetchInit(fetchMethod, seqsToFetch, External, &sock, &tmpError);

  /* Create a struct to pass around data releated to the fetch */
  GeneralFetchData fetchData;
  fetchData.fetchMethod = fetchMethod;
  fetchData.columnList = columnList;
  fetchData.currentColumn = NULL;
  fetchData.seqType = seqType,
  fetchData.status = (tmpError == NULL);

  /* Get the sequences back. They will be returned in the same order that we asked for them, i.e. 
   * in the order they are in our list. */
  fetchData.currentSeqItem = seqsToFetch;
  fetchData.currentSeq = (BlxSequence*)(fetchData.currentSeqItem->data);

  if (!tmpError && fetchData.currentSeq)
    {
      char buffer[RCVBUFSIZE + 1];
      fetchData.buffer = buffer;
      fetchData.numRequested = g_list_length(seqsToFetch); /* total number of sequences requested */
      fetchData.bar = makeProgressBar(fetchData.numRequested);
      fetchData.numFetched = 0;
      fetchData.numSucceeded = 0;
      fetchData.parserState = PARSING_NEWLINE;
      fetchData.curLine = g_string_new("");
      fetchData.sectionId[0] = ' ';
      fetchData.sectionId[1] = ' ';
      fetchData.sectionId[2] = '\0';
      fetchData.tagName = g_string_new("");
      fetchData.currentResult = g_string_new("");
      fetchData.foundEndQuote = FALSE;      
     
      while (fetchData.status &&
             !tmpError &&
             fetchData.parserState != PARSING_CANCELLED && 
             fetchData.parserState != PARSING_FINISHED)
        {
          /* Receive and parse the next buffer */
          checkProgressBar(fetchData.bar, &fetchData.parserState, &fetchData.status);
          fetchData.lenReceived = socketFetchReceiveBuffer(&fetchData, RCVBUFSIZE, sock);
          
          if (fetchMethod->outputType == BLXFETCH_OUTPUT_EMBL)
            {              
              parseEmblBuffer(&fetchData, &tmpError);
            }
          else if (fetchMethod->outputType == BLXFETCH_OUTPUT_RAW)
            {
              parseRawSequenceBuffer(&fetchData, &tmpError);
            }
          else
            {
              g_set_error(error, BLX_ERROR, 1, "Invalid output format for fetch method %s (expected '%s' or '%s')\n", 
                          g_quark_to_string(fetchMethod->name),
                          outputTypeStr(BLXFETCH_OUTPUT_RAW),
                          outputTypeStr(BLXFETCH_OUTPUT_EMBL));
            }
        }
      
      /* Finish up */
      shutdown(sock, SHUT_RDWR);
      destroyProgressBar(fetchData.bar);
      fetchData.bar = NULL ;
      
      if (fetchData.tagName)
        g_string_free(fetchData.tagName, TRUE);

      if (fetchData.currentResult)
        g_string_free(fetchData.currentResult, TRUE);

      if (fetchData.status && !tmpError && fetchData.numSucceeded != fetchData.numRequested)
        {
          double proportionOk = (float)fetchData.numSucceeded / (float)fetchData.numRequested;

          /* We don't display a critical error message when fetching the full EMBL file because we're
           * going to re-try fetching just the fasta data anyway. Display a warning, or just an info
           * message if a small proportion failed */
          if (proportionOk < 0.5)
            {
              g_warning("pfetch sent back %d when %d requested\n", fetchData.numSucceeded, fetchData.numRequested) ;
            }
          else
            {
              g_message("pfetch sent back %d when %d requested\n", fetchData.numSucceeded, fetchData.numRequested) ;
            }
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);

  return fetchData.status ;
}


/* Get the contents of the given config file. Returns the contents in
 * a string which must be freed by the caller. Optionally return the 
 * length of the string. */
static gchar* getConfigFileContent(const char *config_file, gsize *len)
{
  gchar *content = NULL;

  if (g_file_test(config_file, G_FILE_TEST_EXISTS))
    {
      GKeyFile *keyFile = g_key_file_new();
      
      if (g_key_file_load_from_file(keyFile, config_file, G_KEY_FILE_NONE, NULL))
        content = g_key_file_to_data(keyFile, len, NULL);
      
      g_key_file_free(keyFile);
    }

  return content;
}



/* Set/Get global config, necessary because we don't have some blixem context pointer....
 * To do: we do have a context now, so this should be moved to there. 
 * Sets the error if there were any problems. Note that the error is not set if the
 * config file does not exist or is empty  */
void blxInitConfig(const char *config_file, CommandLineOptions *options, GError **error)
{
  g_assert(!blx_config_G) ;

  /* Create the result and set it in the global. We always return something
   * even if it's empty. */
  GKeyFile *key_file = g_key_file_new();
  blx_config_G = key_file;

  /* Get the content of the supplied config file, if any */
  gsize len1 = 0;
  gchar *content1 = getConfigFileContent(config_file, &len1);

  /* Get the content of the local settings file, if any */
  char *settings_file = g_strdup_printf("%s/%s", g_get_home_dir(), BLIXEM_SETTINGS_FILE);
  gsize len2 = 0;
  gchar *content2 = getConfigFileContent(settings_file, &len2);

  /* Merge both file contents into a single string, if applicable,
   * or just use whichever content is non-null, if any */
  gchar *content = NULL;
  gsize len = 0;

  if (content1 && content2)
    {
      /* Merge both file contents. Note that the second file supplied here
       * will take priority if there are duplicate fields across the two files */
      len = len1 + len2 + 2; /* need 2 extra chars for the sprintf below (newline and terminating nul) */
      content = (gchar*)g_malloc(len);
      sprintf(content, "%s\n%s", content1, content2);
      
      /* Free the original strings */
      g_free(content1);
      g_free(content2);
    }
  else if (content1)
    {
      len = len1;
      content = content1;
    }
  else if (content2)
    {
      len = len2;
      content = content2;
    }

  /* Now load the content into the key file, if we have any */  
  if (content)
    {
      if (g_key_file_load_from_data(key_file, content, len, G_KEY_FILE_NONE, NULL))
        {
          /* Parse the config file to read in the blixem options */
          readConfigFile(key_file, options, error);

          if (error && *error && settings_file && config_file)
            prefixError(*error, "Errors found while reading config files '%s' and '%s':\n", config_file, settings_file);
          else if (error && *error && settings_file)
            prefixError(*error, "Errors found while reading config file '%s':\n", settings_file);
          else if (error && *error && config_file)
            prefixError(*error, "Errors found while reading config file '%s':\n", config_file);
        }

      g_free(content);
    }

  g_free(settings_file);
}



GKeyFile *blxGetConfig(void)
{
  return blx_config_G ;
}


/* 
 *                        Internal functions.
 */


static int socketConstruct(const char *ipAddress, int port, gboolean External, GError **error)
{
  int sock ;                       /* socket descriptor */
  struct sockaddr_in *servAddr ;   /* echo server address */
  struct hostent *hp ;


  /* Create a reliable, stream socket using TCP */
  if ((sock = socket(PF_INET, SOCK_STREAM, 0)) < 0)
    {
      g_set_error(error, BLX_FETCH_ERROR, BLX_FETCH_ERROR_SOCKET,
                  "Error creating socket\n") ;
      return -1 ;
    }


  /* Construct the server address structure */
  servAddr = (struct sockaddr_in *) g_malloc (sizeof (struct sockaddr_in)) ;
  hp = gethostbyname(ipAddress) ;
  
  if (!hp)
    {
      g_set_error(error, BLX_FETCH_ERROR, BLX_FETCH_ERROR_HOST,
                  "Unknown host \"%s\"\n", ipAddress);
      return -1;
    }
  
  servAddr->sin_family = AF_INET ;                          /* Internet address family */
  bcopy((char*)hp->h_addr, (char*) &(servAddr->sin_addr.s_addr), hp->h_length) ;
                                                            /* Server IP address */
  servAddr->sin_port = htons ((unsigned short) port) ;      /* Server port */


  /* Establish the connection to the server */
  if (connect(sock, (struct sockaddr *) servAddr, sizeof(struct sockaddr_in)) < 0)
    {
      g_set_error(error, BLX_FETCH_ERROR, BLX_FETCH_ERROR_CONNECT, 
                  "Error connecting socket to host '%s'\n", ipAddress) ;
      sock = -1 ;
    }

  g_free (servAddr) ;

  return sock ;
}


static void socketSend (int sock, const char *text, GError **error)
{
  int len, bytes_to_send, bytes_written ;
  char *tmp ;
  struct sigaction oursigpipe, oldsigpipe ;

  /* The adding of 0x20 to the end looks wierd but I think it's because the  */
  /* server may not hold strings in the way C does (i.e. terminating '\0'),  */
  /* so 0x20 is added to mark the end of the string and the C string term-   */
  /* inator is moved up one and is not actually sent.                        */
  len = strlen(text) ;
  tmp = (char*)g_malloc(len + 2) ;
  strcpy(tmp, text) ;
  tmp[len] = 0x20 ;                                         /* Add new string terminator. */
  tmp[len + 1] = 0 ;
  bytes_to_send = len + 1 ;

  /* send() can deliver a SIGPIPE if the socket has been disconnected, by    */
  /* ignoring it we will receive -1 and can look for EPIPE as the errno.     */
  oursigpipe.sa_handler = SIG_IGN ;
  sigemptyset(&oursigpipe.sa_mask) ;
  oursigpipe.sa_flags = 0 ;
  if (sigaction(SIGPIPE, &oursigpipe, &oldsigpipe) < 0)
    g_error("Cannot set SIG_IGN for SIGPIPE for socket write operations.\n") ;

  bytes_written = send(sock, tmp, bytes_to_send, 0) ;
  if (bytes_written == -1)
    {
      if (errno == EPIPE || errno == ECONNRESET || errno == ENOTCONN)
        {
          char *msg = getSystemErrorText();
          g_set_error(error, BLX_FETCH_ERROR, BLX_FETCH_ERROR_CONNECT,
                      "Socket connection to server has failed, error was: %s\n", msg) ;
          g_free(msg);
        }
      else
      {
        char *msg = getSystemErrorText();
        g_error("Fatal error on socket connection to pfetch server, error was: %s\n", msg) ;
        g_free(msg);
      }
    }
  else if (bytes_written != bytes_to_send)
    {
      g_error("send() call should have written %d bytes, but actually wrote %d.\n", bytes_to_send, bytes_written) ;
    }

  /* Reset the old signal handler.                                           */
  if (sigaction(SIGPIPE, &oldsigpipe, NULL) < 0)
    {
      g_error("Cannot reset previous signal handler for signal SIGPIPE for socket write operations.\n") ;
    }

  g_free(tmp) ;
}


#ifdef PFETCH_HTML 

static PFetchStatus pfetch_reader_func(PFetchHandle *handle,
                                       char         *text,
                                       guint        *actual_read,
                                       GError       *error,
                                       gpointer      user_data)
{
  PFetchData pfetch_data = (PFetchData)user_data;
  PFetchStatus status    = PFETCH_STATUS_OK;

  if (actual_read && *actual_read > 0 && pfetch_data)
    {
      GtkTextBuffer *text_buffer = pfetch_data->text_buffer;

      /* clear the buffer the first time... */
      if(pfetch_data->got_response == FALSE)
        {
          gtk_text_buffer_set_text(text_buffer, "", 0);
        }

      gtk_text_buffer_insert_at_cursor(text_buffer, text, *actual_read);
      pfetch_data->got_response = TRUE;

      /* If we tried fetching the full entry and failed, try again
       * with the next fetch method, if there is one */
      if (stringInArray(text, pfetch_data->fetchMethod->errors))
        {
          fetchSequence(pfetch_data->blxSeq, TRUE, pfetch_data->attempt + 1, pfetch_data->blxWindow, pfetch_data->dialog, &pfetch_data->text_buffer);
        }
    }

  return status;
}



static PFetchStatus pfetch_closed_func(gpointer user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK;

  DEBUG_OUT("pfetch closed\n");

  return status;
}


/* GtkObject destroy signal handler */
static void handle_dialog_close(GtkWidget *dialog, gpointer user_data)
{
  PFetchData pfetch_data   = (PFetchData)user_data;
  pfetch_data->text_buffer = NULL;
  pfetch_data->widget_destroy_handler_id = 0; /* can we get this more than once? */

  if(pfetch_data->pfetch)
    pfetch_data->pfetch = PFetchHandleDestroy(pfetch_data->pfetch);

  return ;
}



/* 
 *             Multiple sequence fetch callbacks.
 */

static PFetchStatus sequence_pfetch_reader(PFetchHandle *handle,
                                           char         *text,
                                           guint        *actual_read,
                                           GError       *error,
                                           gpointer      user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK ;
  PFetchSequence fetch_data = (PFetchSequence)user_data ;
  ProgressBar bar = fetch_data->fetchData.bar ;

  if (fetch_data->fetchData.parserState != PARSING_FINISHED && 
      fetch_data->fetchData.parserState != PARSING_CANCELLED)
    {
      if (!(*text) || *actual_read <= 0)
        {
          fetch_data->fetchData.parserState = PARSING_FINISHED ;
          fetch_data->fetchData.status = FALSE ;
          fetch_data->err_txt = g_strdup("No data returned by http proxy server.") ;
        }
      else if (*actual_read > 0)
        {
          if (isCancelledProgressBar(bar))
            {
              fetch_data->fetchData.status = FALSE ;
              fetch_data->fetchData.parserState = PARSING_CANCELLED ;
      
              status = PFETCH_STATUS_FAILED ;
            }
          else if (!parsePfetchHtmlBuffer(fetch_data->fetchData.fetchMethod, text, *actual_read, fetch_data))
            {
              status = PFETCH_STATUS_FAILED ;
            }
        }
    }

  return status ;
}

static PFetchStatus sequence_pfetch_closed(PFetchHandle *handle, gpointer user_data)
{
  PFetchStatus status = PFETCH_STATUS_OK;
  PFetchSequence fetch_data = (PFetchSequence)user_data ;


  DEBUG_OUT("pfetch closed\n");


  fetch_data->fetchData.parserState = PARSING_FINISHED ;
  fetch_data->connection_closed = TRUE;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
  destroyProgressBar(fetch_data->bar) ;
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */

  if (fetch_data->stats)
    {
      g_message("Stats for %d reads:\tmin = %d\tmax = %d\tmean = %d\n",
                fetch_data->total_reads,
                fetch_data->min_bytes, fetch_data->max_bytes,
                (fetch_data->total_bytes / fetch_data->total_reads)) ;
    }

  return status ;
}


/* Progress meter has been destroyed so free all pfetch stuff. */
static void sequence_dialog_closed(GtkWidget *dialog, gpointer user_data)
{
  PFetchSequence fetch_data = (PFetchSequence)user_data ;

  fetch_data->pfetch = PFetchHandleDestroy(fetch_data->pfetch) ;

  return ;
}


/* Utility to record the stats for the current read buffer, if stats are enabled */
static void pfetchHtmlRecordStats(const char *read_text, const int length, PFetchSequence fetch_data)
{
  if (fetch_data->stats)
    {
      if (length < fetch_data->min_bytes)
        fetch_data->min_bytes = length ;

      if (length > fetch_data->max_bytes)
        fetch_data->max_bytes = length ;

      fetch_data->total_bytes += length ;

      fetch_data->total_reads++ ;
  
      if (length < 100)
        {
          char *str ;

          str = g_strndup(read_text, length) ;
          g_message("Read %d: \"%s\"\n", fetch_data->total_reads, str) ;
          g_free(str) ;
        }
    }
}


/* Parse the buffer sent back by proxy server. The pfetch server sends back data separated
 * by newlines. The data is either a valid IUPAC dna sequence or a valid IUPAC peptide
 * sequence or an error message. The problem here is that the data is not returned to
 * this function in complete lines so we have to reconstruct the lines as best we can. It's
 * even possible for very long sequences that they may span several buffers. Note also
 * that the buffer is _not_ null terminated, we have to use length to know when to stop reading.
 */
static gboolean parsePfetchHtmlBuffer(const BlxFetchMethod* const fetchMethod,
                                      char *read_text,
                                      int length,
                                      PFetchSequence fetch_data)
{
  gboolean status = TRUE ;

  /* Validate input and record stats, if requested */
  g_assert(fetch_data && fetch_data->fetchData.currentSeqItem && fetch_data->fetchData.currentSeqItem->data);
  pfetchHtmlRecordStats(read_text, length, fetch_data);
  
  GError *error = NULL;
  
  fetch_data->fetchData.buffer = read_text;
  fetch_data->fetchData.lenReceived = length;

  if (fetch_data->fetchData.fetchMethod->outputType == BLXFETCH_OUTPUT_EMBL)
    {
      /* We're fetching the full EMBL entries */
      parseEmblBuffer(&fetch_data->fetchData, &error);
    }
  else if (fetch_data->fetchData.fetchMethod->outputType == BLXFETCH_OUTPUT_RAW)
    {
      /* The fetched entries just contain the FASTA sequence */
      parseRawSequenceBuffer(&fetch_data->fetchData, &error);
    }
  else 
    {
      g_set_error(&error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT, 
                  "Invalid output format specified for fetch method '%s'; expected '%s' or '%s'\n",
                  g_quark_to_string(fetch_data->fetchData.fetchMethod->name), 
                  outputTypeStr(BLXFETCH_OUTPUT_EMBL), 
                  outputTypeStr(BLXFETCH_OUTPUT_RAW));
    }
  
  if (error)
    {
      fetch_data->err_txt = g_strdup(error->message);
      g_error_free(error);
      error = NULL;
    }
  
  status = fetch_data->fetchData.status ;
  
  return status;
}



#endif /* PFETCH_HTML */



/* Functions to display, update, cancel and remove a progress meter. */

static ProgressBar makeProgressBar(int seq_total)
{
  ProgressBar bar = g_new0(ProgressBarStruct, 1) ;

  bar->seq_total = seq_total ;
  bar->cancelled = FALSE;
  bar->widget_destroy_handler_id = 0;
  
  gdk_color_parse("blue", &(bar->blue_bar_fg)) ;
  gdk_color_parse("red", &(bar->red_bar_fg)) ;

  bar->top_level = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  char *title = g_strdup_printf("Blixem - pfetching %d sequences...", seq_total) ;
  gtk_window_set_title(GTK_WINDOW(bar->top_level), title) ;
  g_free(title) ;
  g_signal_connect(G_OBJECT(bar->top_level), "destroy", G_CALLBACK(destroyProgressCB), bar) ;


  gtk_window_set_default_size(GTK_WINDOW(bar->top_level), 350, -1) ;

  GtkWidget *frame = gtk_frame_new(NULL) ;
  gtk_container_add(GTK_CONTAINER(bar->top_level), frame) ;

  GtkWidget *vbox = gtk_vbox_new(FALSE, 5) ;
  gtk_container_add(GTK_CONTAINER(frame), vbox) ;

  bar->progress = gtk_progress_bar_new() ; 

  gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->blue_bar_fg)) ;

  gtk_box_pack_start(GTK_BOX(vbox), bar->progress, TRUE, TRUE, 0);

  bar->label = gtk_label_new("") ;
  gtk_box_pack_start(GTK_BOX(vbox), bar->label, TRUE, TRUE, 0);

  GtkWidget *hbox = gtk_hbox_new(FALSE, 0) ;
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
  gtk_container_border_width(GTK_CONTAINER(hbox), 5);

  GtkWidget *cancel_button = gtk_button_new_with_label("Cancel") ;
  gtk_box_pack_start(GTK_BOX(hbox), cancel_button, TRUE, TRUE, 0) ;

  gtk_signal_connect(GTK_OBJECT(cancel_button), "clicked", GTK_SIGNAL_FUNC(cancelCB), (gpointer)bar) ;

  gtk_widget_show_all(bar->top_level) ;
  
  while (gtk_events_pending())
    {
      gtk_main_iteration() ;
    }

  return bar ;
}

static void updateProgressBar(ProgressBar bar, const char *sequence, int numFetched, gboolean fetch_ok)
{
  char *label_text ;

  if (fetch_ok)
    gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->blue_bar_fg)) ;
  else
    gtk_widget_modify_fg(bar->progress, GTK_STATE_NORMAL, &(bar->red_bar_fg)) ;

  label_text = g_strdup_printf("%s - %s", sequence, fetch_ok ? "Fetched" : "Not found") ;
  gtk_label_set_text(GTK_LABEL(bar->label), label_text) ;
  g_free(label_text) ;
  
  gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(bar->progress),
                                (double)((double)numFetched / (double)(bar->seq_total))) ;

  while (gtk_events_pending())
    gtk_main_iteration();

  return ;
}


static gboolean isCancelledProgressBar(ProgressBar bar)
{
  return bar->cancelled ;
}


static void destroyProgressBar(ProgressBar bar)
{
  gtk_widget_destroy(bar->top_level) ;
  
  return ;
}


/* called from progress meter when being destroyed, do not call directly. */
static void destroyProgressCB(GtkWidget *widget, gpointer cb_data)
{
  ProgressBar bar = (ProgressBar)cb_data ;

  g_free(bar) ;

  return ;
}


/* Called from progress meter. */
static void cancelCB(GtkWidget *widget, gpointer cb_data)
{
  ProgressBar bar = (ProgressBar)cb_data ;

  bar->cancelled = TRUE ;

  return ;
}





/* 
 *            Configuration file reading.
 */

/* Get default blixem options from the "blixem" stanza */
static void readBlixemStanza(GKeyFile *key_file,
                             const char *group, 
                             CommandLineOptions *options,
                             GError **error)
{
  /* Note that all values are optional, so don't set the error if they
   * do not exist */

  /* Get the window background color if set */
  options->windowColor = g_key_file_get_string(key_file, group, SEQTOOLS_WINDOW_COLOR, NULL);

  /* Get the comma-separated list of possible fetch methods */
  options->bulkFetchDefault = keyFileGetCsv(key_file, group, SEQTOOLS_BULK_FETCH, NULL);
  options->userFetchDefault = keyFileGetCsv(key_file, group, SEQTOOLS_USER_FETCH, NULL);
  options->optionalFetchDefault = keyFileGetCsv(key_file, group, SEQTOOLS_OPTIONAL_FETCH, NULL);

  /* If the bulk-fetch key wasn't found, try the old
   * default-fetch-mode key, for backwards compatibility */
  if (!options->bulkFetchDefault)
    options->bulkFetchDefault = keyFileGetCsv(key_file, group, BLIXEM_OLD_BULK_FETCH, NULL);

  /* Get the default values for the MSP flags, if they're specified in the blixem stanza. */
  int flag = MSPFLAG_MIN + 1;
  for ( ; flag < MSPFLAG_NUM_FLAGS; ++flag)
    {
      GError *tmpError = NULL;
      const char *key = mspFlagGetConfigKey((MspFlag)flag);
      gboolean value = g_key_file_get_boolean(key_file, group, key, &tmpError);
      
      /* If found, set it as the default */
      if (!tmpError)
        mspFlagSetDefault((MspFlag)flag, value);
    }
}


/* Create a fetch method struct */
static BlxFetchMethod* createBlxFetchMethod(const char *fetchName, 
                                            const GQuark fetchMode)
{
  BlxFetchMethod *result = (BlxFetchMethod*)g_malloc(sizeof *result);

  result->name = g_quark_from_string(fetchName);
  result->mode = BLXFETCH_MODE_NONE;

  result->location = NULL;
  result->node = NULL;
  result->port = 0;
  result->cookie_jar = NULL;
  result->args = NULL;
  result->columns = NULL;

  result->separator = NULL;
  result->errors = NULL;
  result->outputType = BLXFETCH_OUTPUT_INVALID;

  return result;
}


/* Get the output type from the given fetch stanza. Returns BLXFETCH_OUTPUT_INVALID if
 * not found. The input error value can be non-null, in which case any new error message
 * will be appended. */
static BlxFetchOutputType readFetchOutputType(GKeyFile *key_file, const char *group, GError **error)
{
  BlxFetchOutputType result = BLXFETCH_OUTPUT_INVALID;

  GError *tmpError = NULL;
  char *outputTypeName = removeDelimiters(g_key_file_get_string(key_file, group, FETCH_OUTPUT, &tmpError));

  if (!tmpError)
    {
      /* Loop through all output types looking for one with this name */
      int outputType = 0;
      
      for ( ; outputType < BLXFETCH_NUM_OUTPUT_TYPES; ++outputType)
        {
          if (stringsEqual(outputTypeName, outputTypeStr((BlxFetchOutputType)outputType), FALSE))
            {
              result = (BlxFetchOutputType)outputType;
              break;
            }
        }
      
      if (result == BLXFETCH_OUTPUT_INVALID)
        {
          g_set_error(&tmpError, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT,
                      "Invalid output type '%s' for fetch method '%s'", outputTypeName, group);
        }
    }
  
  if (tmpError && error && *error)
    postfixError(*error, "; %s", tmpError->message);
  else if (tmpError)
    g_propagate_error(error, tmpError);
 
  g_free(outputTypeName);

  return result;
}


/* Utilities to get a value from a key-file and remove delimiters and whitespace (where applicable).
 * The input error may be non-null, in which case any new error message will be appended to it. */
static char* configGetString(GKeyFile *key_file, const char *group, const char *key, GError **error)
{
  GError *tmpError = NULL;
  
  char *result = g_key_file_get_string(key_file, group, key, &tmpError);

  if (result)
    {
      /* remove any leading/trailing whitespace */
      result = g_strchug(g_strchomp(result));

      /* remove any delimiters */
      result = removeDelimiters(result);
    }
  
  if (tmpError && error && *error)
    postfixError(*error, "; %s", tmpError->message);
  else if (tmpError)
    g_propagate_error(error, tmpError);

  return result;
}

static int configGetInteger(GKeyFile *key_file, const char *group, const char *key, GError **error)
{
  GError *tmpError = NULL;
  int result = g_key_file_get_integer(key_file, group, SOCKET_FETCH_PORT, &tmpError);

  if (tmpError && error && *error)
    postfixError(*error, "; %s", tmpError->message);
  else if (tmpError)
    g_propagate_error(error, tmpError);

  return result;
}


/* Get details about the given fetch method stanza and add it to 
 * the list of fetch methods in 'options' */
static void readFetchMethodStanza(GKeyFile *key_file, 
                                  const char *group, 
                                  CommandLineOptions *options, 
                                  GError **error)
{
  BlxFetchMethod *result = NULL;

  GError *tmpError = NULL;
  char *fetchMode = configGetString(key_file, group, FETCH_MODE_KEY, &tmpError);

  if (!tmpError)
    result = createBlxFetchMethod(group, g_quark_from_string(fetchMode));

  /* Set the relevant properties for this type of fetch method. (Append a
   * warning to tmpError for mandatory arguments if they are not found.) */
  if (!tmpError)
    {
      if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_SOCKET), FALSE))
        {
          result->mode = BLXFETCH_MODE_SOCKET;
          result->location = configGetString(key_file, group, SOCKET_FETCH_LOCATION, &tmpError);
          result->node = configGetString(key_file, group, SOCKET_FETCH_NODE, &tmpError);
          result->port = configGetInteger(key_file, group, SOCKET_FETCH_PORT, &tmpError);
          result->args = configGetString(key_file, group, SOCKET_FETCH_ARGS, NULL);
          result->errors = keyFileGetCsv(key_file, group, FETCH_ERRORS, NULL);
          result->separator = configGetString(key_file, group, FETCH_SEPARATOR, NULL);
          result->outputType = readFetchOutputType(key_file, group, &tmpError);
        }
    #ifdef PFETCH_HTML
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_HTTP), FALSE))
        {
          result->mode = BLXFETCH_MODE_HTTP;
          result->location = configGetString(key_file, group, HTTP_FETCH_LOCATION, &tmpError);
          result->port = configGetInteger(key_file, group, HTTP_FETCH_PORT, &tmpError);
          result->cookie_jar = configGetString(key_file, group, HTTP_FETCH_COOKIE_JAR, &tmpError);
          result->args = configGetString(key_file, group, HTTP_FETCH_ARGS, NULL);
          result->errors = keyFileGetCsv(key_file, group, FETCH_ERRORS, NULL);
          result->separator = configGetString(key_file, group, FETCH_SEPARATOR, NULL);
          result->outputType = readFetchOutputType(key_file, group, &tmpError);
        }
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_PIPE), FALSE))
        {
          result->mode = BLXFETCH_MODE_PIPE;
          result->location = configGetString(key_file, group, PIPE_FETCH_LOCATION, &tmpError);
          result->args = configGetString(key_file, group, PIPE_FETCH_ARGS, NULL);
          result->errors = keyFileGetCsv(key_file, group, FETCH_ERRORS, NULL);
          result->outputType = readFetchOutputType(key_file, group, &tmpError);
        }
    #endif
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_WWW), FALSE))
        {
          result->mode = BLXFETCH_MODE_WWW;
          result->location = configGetString(key_file, group, WWW_FETCH_LOCATION, &tmpError);
          result->args = configGetString(key_file, group, WWW_FETCH_ARGS, NULL);
        }
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_INTERNAL), FALSE))
        {
          result->mode = BLXFETCH_MODE_INTERNAL;
        }
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_SQLITE), FALSE))
        {
          result->mode = BLXFETCH_MODE_SQLITE;
          result->location = configGetString(key_file, group, DB_FETCH_LOCATION, &tmpError);
          result->args = configGetString(key_file, group, DB_FETCH_QUERY, &tmpError);
          result->errors = keyFileGetCsv(key_file, group, FETCH_ERRORS, NULL);
          result->outputType = readFetchOutputType(key_file, group, &tmpError);
          
          /* Hard code the separator */
          result->separator = g_strdup("','");
        }
      else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_COMMAND), FALSE))
        {
          result->mode = BLXFETCH_MODE_COMMAND;
          result->location = configGetString(key_file, group, COMMAND_FETCH_SCRIPT, &tmpError);
          result->args = configGetString(key_file, group, COMMAND_FETCH_ARGS, NULL);
          result->errors = keyFileGetCsv(key_file, group, FETCH_ERRORS, NULL);
          result->outputType = readFetchOutputType(key_file, group, &tmpError);
        }
        else if (stringsEqual(fetchMode, fetchModeStr(BLXFETCH_MODE_NONE), FALSE))
        {
          result->mode = BLXFETCH_MODE_NONE;
        }
      else
        {
          g_set_error(&tmpError, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_INVALID_FETCH_MODE, "Unrecognised fetch mode '%s'", fetchMode);
          result->mode = BLXFETCH_MODE_NONE;
        }
    }
  
  /* Add result to list */
  if (!tmpError && result)
    {
      GQuark fetchName = g_quark_from_string(group);
      g_hash_table_insert(options->fetchMethods, GINT_TO_POINTER(fetchName), result);
    }
  else if (result)
    {
      /* Fetch method details are incomplete, so delete it. */
      g_free(result);
    }
  
  /* Clean up */
  g_free(fetchMode);

  /* Error handling */
  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Tries to load values for this group if it is a standard, recognised group or 
 * group type. If any required values are missing, sets the error. */
static void readConfigGroup(GKeyFile *key_file, const char *group, CommandLineOptions *options, GError **error)
{
  /* Check for known group names first */
  if (stringsEqual(group, BLIXEM_GROUP, FALSE))
    {
      readBlixemStanza(key_file, group, options, error);
    }
  else
    {
      /* Check if this is a fetch method; if it is it the group will have a fetch-mode key */
      char *fetchMode = configGetString(key_file, group, FETCH_MODE_KEY, NULL);
      
      if (fetchMode)
        readFetchMethodStanza(key_file, group, options, error);

      g_free(fetchMode);
    }
}


/* Read standard stanzas from the given config file. Sets the error if any problems */
static void readConfigFile(GKeyFile *key_file, CommandLineOptions *options, GError **error)
{
  if (!key_file)
    return;

  /* Loop through each stanza in the config file */
  gsize num_groups ;
  char **groups = g_key_file_get_groups(key_file, (gsize*)(&num_groups));
  char **group = groups;
  int i = 0;
  
  for ( ; i < (int)num_groups ; i++, group++)
    {
      /* Read in the data in this stanza */
      GError *tmpError = NULL;
      readConfigGroup(key_file, *group, options, &tmpError) ;
      
      if (tmpError && error)
        {
          /* Compile all errors into one message */
          if (*error == NULL)
            g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_WARNINGS, "  [%s]: %s\n", *group, tmpError->message);
          else
            postfixError(*error, "  [%s]: %s\n", *group, tmpError->message);
        }
    }
  
  /* For backwards compatibility with config files that do not set the optional-fetch
   * field: If the default optional-fetch method is not set then find any fetch methods that
   * have an output type of 'embl'. gb10: we may wish to remove this once we no longer have these
   * config files. */
  if (!options->optionalFetchDefault)
    {
      GHashTableIter iter;
      gpointer key, value;

      g_hash_table_iter_init (&iter, options->fetchMethods);
      while (g_hash_table_iter_next (&iter, &key, &value))
        {
          BlxFetchMethod* fetchMethod = (BlxFetchMethod*)value;

          if (fetchMethod->outputType == BLXFETCH_OUTPUT_EMBL)
            {
              if (!options->optionalFetchDefault)
                options->optionalFetchDefault = g_array_sized_new(FALSE, TRUE, sizeof(GQuark), 1);

              g_array_append_val(options->optionalFetchDefault, fetchMethod->name);
            }
        }
    }

  g_strfreev(groups);
}


/*******************************************************************
 *                      Pfetch utility functions                   *
 *******************************************************************/


/* Initialise a pfetch connection with the given options and calls pfetch on each of the 
 * sequences in the given list. */
static void socketFetchInit(const BlxFetchMethod* const fetchMethod,
                            GList *seqsToFetch, 
                            gboolean External, 
                            int *sock,
                            GError **error)
{
  GError *tmpError = NULL;

  /* open socket connection */
  if (!tmpError)
    {
      *sock = socketConstruct(fetchMethod->node, fetchMethod->port, External, &tmpError) ;
    }
  
  /* send the command/names to the server */
  if (!tmpError)
    {
      GString *command = getFetchArgsMultiple(fetchMethod, seqsToFetch, &tmpError);
      socketSend(*sock, command->str, &tmpError);
      
      if (command)
        g_string_free(command, TRUE);
    }
  
//  if (!tmpError)
//    {
//      /* For each sequence, send a command to fetch that sequence, in the order that they are in our list */
//      GList *seqItem = seqsToFetch;
//      
//      for ( ; seqItem && !tmpError ; seqItem = seqItem->next)
//        {
//          BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
//          socketSend(*sock, blxSequenceGetName(blxSeq), &tmpError);
//        }
//    }
  
  if (!tmpError)
    {
      /* send a final newline to flush the socket */
      if (send(*sock, "\n", 1, 0) != 1)
        {
          g_set_error(&tmpError, BLX_FETCH_ERROR, BLX_FETCH_ERROR_SEND,
                      "Failed to send final newline to socket\n") ;
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Check the progress bar to see if the user cancelled. If so, sets the state to cancelled
 * and the status to false (meaning don't continue). */
static void checkProgressBar(ProgressBar bar, BlxSeqParserState *parserState, gboolean *status)
{
  if (isCancelledProgressBar(bar))
    {
      *status = FALSE;
      *parserState = PARSING_CANCELLED;
    }
}


/* Receive the next buffer back from a socket server. Only does anything if the status is ok
 * (i.e. true). Checks the length etc and sets the state to finished if there is no more data
 * to receive. Sets 'status' to false if there was an error. Returns the number of chars received. */
static int socketFetchReceiveBuffer(GeneralFetchData *fetchData, const int bufferSize, const int sock)
{
  int lenReceived = 0;
  
  if (fetchData->status == FALSE)
    {
      return lenReceived;
    }
  
  /* Ask for the next chunk of data to be put into our buffer */
  lenReceived = recv(sock, fetchData->buffer, bufferSize, 0);
  
  if (lenReceived < 0)
    {
      /* Problem with this one - skip to the next */
      fetchData->status = FALSE;
      char *msg = getSystemErrorText();
      g_critical("Could not retrieve sequence data from pfetch server, error was: %s\n", msg) ;
      g_free(msg);
    }
  else if (lenReceived == 0)
    {
      /* No more data, so quit out of the loop */
      fetchData->parserState = PARSING_FINISHED;
    }
  
  return lenReceived;
}


/* Update the given sequence and sequence item to the next item in the list (unless the state is
 * 'finished' in which case there is nothing to do). Sets status to false if there was an
 * error (i.e. if we ran out of sequences). */
static void pfetchGetNextSequence(GeneralFetchData *fetchData, const gboolean pfetch_ok)
{
  if (fetchData->parserState != PARSING_FINISHED)
    {
      /* Check that another BlxSequence item exists in the list */
      if (fetchData->currentSeqItem->next)
        {
          updateProgressBar(fetchData->bar, blxSequenceGetName(fetchData->currentSeq), fetchData->numFetched, pfetch_ok) ;
          
          /* Move to the next BlxSequence */
          fetchData->currentSeqItem = fetchData->currentSeqItem->next;
          fetchData->currentSeq = (BlxSequence*)(fetchData->currentSeqItem->data);
        }
      else
        {
          fetchData->status = FALSE ;
          g_critical("Unexpected data from pfetch server: received too many lines. (%d sequences were requested.)\n", fetchData->numRequested);
        }
    }
}


/* Parse the given buffer that contains an arbitrary section
 * of data from a fasta file (or concatenation of multiple 
 * fasta files). The parserState indicates on entry what 
 * state we are in and is updated on exit with the new state,
 * if it has changed. */
static void parseRawSequenceBuffer(GeneralFetchData *fetchData, GError **error)
{
  if (fetchData->status == FALSE || fetchData->parserState == PARSING_FINISHED || fetchData->parserState == PARSING_CANCELLED)
    {
      return;
    }
  
  /* Loop through each character in the buffer */
  int i = 0;
  
  for ( ; i < fetchData->lenReceived && fetchData->status; ++i)
    {
      /* Check for user cancellation again */
      checkProgressBar(fetchData->bar, &fetchData->parserState, &fetchData->status);
      
      if (fetchData->parserState == PARSING_CANCELLED)
        {
          break;
        }

      const char curChar = fetchData->buffer[i];
      
      if (curChar == '\n')
        {
          /* finish up this sequence and move to the next one */
          gboolean pfetch_ok = pfetchFinishSequence(fetchData);
          pfetchGetNextSequence(fetchData, pfetch_ok);
        }
      else
        {
          /* Append this character to the current result string */
          g_string_append_c(fetchData->currentResult, curChar);
        }
    }
}


/* Returns true if the given string is the terminating string of
 * an embl file, i.e. '//' */
static gboolean isEmblTerminator(GQuark value)
{
  static GQuark terminator = 0;
  
  if (!terminator)
    terminator = g_quark_from_string("//");

  return (value == terminator);
}


/* Get the parser state from the given two-letter at at the start of an EMBL file line. */
static void pfetchGetParserStateFromId(GeneralFetchData *fetchData)
{
  GQuark sectionId = g_quark_from_string(fetchData->sectionId);

  /* First, check if the section id is the terminator string, in 
   * which case finish this sequence */
  if (isEmblTerminator(sectionId))
    {
      fetchData->parserState = PARSING_FINISHED_SEQ;
      return;
    }

  /* Loop through all the columns and check if there's a column with
   * an embl ID that matches the current section ID. If not, then there's 
   * nothing to do for this tag (so default to parsing_ignore if not found) */
  GList *item = fetchData->columnList;
  fetchData->parserState = PARSING_IGNORE;
  
  for ( ; item; item = item->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
      
      if (columnInfo->emblId && columnInfo->emblId == sectionId)
        {
          /* First, check if we need to bother getting the info for this column */
          if (blxSequenceRequiresColumnData(fetchData->currentSeq, columnInfo->columnId) &&
              !blxSequenceGetValueAsString(fetchData->currentSeq, columnInfo->columnId))
            {
              if (columnInfo->emblTag)
                {
                  /* The column has an embl tag; continue parsing to find the
                   * tag. (Note that multiple columns may have this section ID so
                   * the next step may end up finding a different column depending
                   * on the embl tag; at this point, we're just finding out if any 
                   * column has an embl tag that we need to continue parsing for,
                   * so we don't set the currentColumn yet). */
                  fetchData->parserState = PARSING_TAG_SEARCH;
                }
              else
                {
                  /* No tag to look for so start parsing the data for this column immediately */
                  fetchData->currentColumn = columnInfo;

                  /* If this is the sequence column, we need to parse the sequence header
                   * first, so set a special parser state. Other columns can be parsed directly. */
                  if (columnInfo->columnId == BLXCOL_SEQUENCE)
                    fetchData->parserState = PARSING_SEQUENCE_HEADER;
                  else
                    fetchData->parserState = PARSING_DATA;
                }

              break; /* we found a column, so exit the loop */
            }
        }
    }
}


/* Process the current character in the current buffer of EMBL data. The char might be read into
 * various data locations depending on the current parser state, or might cause the parser state
 * to change. */
static void pfetchProcessEmblBufferChar(GeneralFetchData *fetchData, const char curChar)
{
  switch (fetchData->parserState)
    {
      case PARSING_NEWLINE:
      {
        /* First char should be the start of the two-letter ID */
        fetchData->sectionId[0] = curChar;
        fetchData->parserState = PARSING_ID;
        break;
      }
        
      case PARSING_ID:
      {
        /* Read the second (and last) character from the ID. */
        fetchData->sectionId[1] = curChar;
        pfetchGetParserStateFromId(fetchData);
        break;
      }
        
      case PARSING_SEQUENCE:
      {
        /* Look out for the terminating "//" characters. We ignore newlines in the sequence
         * body so won't catch these by the normal method. We can assume that if we see a single
         * slash that it's the terminator, because we should never see one here otherwise */
        if (curChar == '/')
          {
            fetchData->sectionId[0] = curChar;
            fetchData->parserState = PARSING_ID; /* continue parsing to read the other '/' */
          }
        else
          {
            /* Add the character to the end of the sequence, but only if it is a valid IUPAC character
             * (because we want to ignore whitespace, newlines and digits). */
            if (isValidIupacChar(curChar, fetchData->seqType))
              {
                appendCharToString(curChar, fetchData->currentResult);
              }
          }
        
        break;
      }
        
      case PARSING_DATA:
      {
        appendCharToString(curChar, fetchData->currentResult);
        break;
      }
        
      case PARSING_TAG_SEARCH:
      {
        /* If we find a forward slash, it's the start of a tag, so we'll parse the tag name next */
        if (curChar == '/')
          fetchData->parserState = PARSING_TAG_NAME;
        
        break;
      }
        
      case PARSING_TAG_NAME:
      {
        if (curChar == '=') /* signals the end of the tag */
          pfetchGetParserStateFromTagName(fetchData);
        else
          g_string_append_c(fetchData->tagName, curChar); /* read in tag name */
        
        break;
      }
        
      case PARSING_TAG_IGNORE:
      {
        /* Look for the start quote before we start parsing the tag properly */
        if (curChar == '"')
          pfetchGetParserStateFromTagName(fetchData);
        
        break;
      }
        
      case PARSING_DATA_QUOTED:
      {
        /* Like PARSING_DATA but we need to watch out for end quotes */
        appendCharToQuotedString(curChar, &fetchData->foundEndQuote, fetchData->currentResult);
        break;
      }
        
      default:
        break;
        
    };
}


static void checkParserState(GeneralFetchData *fetchData, BlxSeqParserState origState)
{
  /* If we were parsing data for a column but the parser state has now changed to
   * something else, then we've finished parsing the string, so save it in the 
   * sequence */
  if (fetchData->parserState != origState && fetchData->currentResult->str)
    {
      if (origState == PARSING_DATA || origState == PARSING_DATA_QUOTED)
        {
          blxSequenceSetValueFromString(fetchData->currentSeq, 
                                        fetchData->currentColumn->columnId, 
                                        fetchData->currentResult->str); 

          /* Reset the result string and current column */
          g_string_truncate(fetchData->currentResult, 0);
          fetchData->currentColumn = NULL;
        }
    }

  /* If parsing for this sequence has finished, tidy it up and move to the next sequence */
  if (fetchData->parserState == PARSING_FINISHED_SEQ)
    {
      gboolean pfetch_ok = pfetchFinishSequence(fetchData);
      pfetchGetNextSequence(fetchData, pfetch_ok);
      fetchData->parserState = PARSING_NEWLINE;
    }
}


/* Parse the given buffer, which contains an arbitrary section of
 * data from an EMBL entry (or combination of multiple embl entries).
 * The parserState indicates what state we are in on entry, and gets
 * updated with the new state on exit. */
static void parseEmblBuffer(GeneralFetchData *fetchData, GError **error)
{
  if (fetchData->status == FALSE || 
      fetchData->parserState == PARSING_FINISHED || 
      fetchData->parserState == PARSING_CANCELLED)
    {
      return;
    }

  /* Loop through each character in the buffer */
  int i = 0;
  
  for ( ; i < fetchData->lenReceived && fetchData->status; ++i)
    {
      checkProgressBar(fetchData->bar, &fetchData->parserState, &fetchData->status);

      if (fetchData->parserState == PARSING_CANCELLED)
        {
          break;
        }

      /* Remember the state prior to processing this character */
      BlxSeqParserState origState = fetchData->parserState;

      const char curChar = fetchData->buffer[i];
      g_string_append_c(fetchData->curLine, curChar);

      /* Special treatment if we've previously found an end quote: if this char is NOT also
       * a quote, then it means it genuinely was an end quote, so we finish reading the tag data that
       * we were reading (by setting the tag name to null). We return to parsing the section
       * that contained that tag in case there are any more tags in it that we're interested in. */
      if (fetchData->foundEndQuote && curChar != '"')
        {
          g_string_truncate(fetchData->tagName, 0);
          fetchData->foundEndQuote = FALSE;
          fetchData->parserState = PARSING_TAG_SEARCH;
        }
      

      /* Newline character means we need to re-read the 2-letter ID at the start of the
       * line to find out what section we're in. Ignore newlines in the sequence body, though,
       * because it spans several lines and does not have an ID on each line. */
      if (curChar == '\n' && fetchData->parserState != PARSING_SEQUENCE)
        {
          if (fetchData->parserState == PARSING_SEQUENCE_HEADER)
            {
              /* We were previously in the sequence header, so the next line is the sequence body */
              fetchData->parserState = PARSING_SEQUENCE;
            }
          else if (stringInArray(fetchData->curLine->str, fetchData->fetchMethod->errors))
            {
              /* The line is an error message. Finish this sequence and move to next. */
              g_warning("[%s] Error fetching sequence '%s': %s", g_quark_to_string(fetchData->fetchMethod->name), blxSequenceGetName(fetchData->currentSeq), fetchData->curLine->str);
              fetchData->parserState = PARSING_FINISHED_SEQ;
            }
          else
            {
              fetchData->parserState = PARSING_NEWLINE;
            }

          g_string_truncate(fetchData->curLine, 0);
        }
      else
        {
          pfetchProcessEmblBufferChar(fetchData, curChar);
        }

      checkParserState(fetchData, origState);
    }
}


/* Check the given tag name to see if it contains data for a column we're interested
 * in, and update the parser state accordingly. Sets the parser state to PARSING_IGNORE 
 * if this is not a tag we care about. */
static void pfetchGetParserStateFromTagName(GeneralFetchData *fetchData)
{
  if (!fetchData->tagName || !fetchData->tagName->str)
    return;

  /* Check if there's a column with this section ID and tag. (Assumes
   * that only one column at most will match.) */
  GQuark foundSection = g_quark_from_string(fetchData->sectionId);
  GQuark foundTag = g_quark_from_string(fetchData->tagName->str);
  GList *item = fetchData->columnList;
  fetchData->currentColumn = NULL;

  for ( ; item && !fetchData->currentColumn; item = item->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
      GQuark currentSection = columnInfo->emblId;
      GQuark currentTag = columnInfo->emblTag;

      if (currentSection == foundSection && currentTag == foundTag)
        {
          fetchData->currentColumn = columnInfo;
          
          if (fetchData->parserState == PARSING_TAG_NAME)
            {
              /* Start parsing the data within the tag but initially ignore
               * everything until we find a quoted section */
              fetchData->parserState = PARSING_TAG_IGNORE;
            }
          else
            {
              /* Should only get here if we found a quoted section. Start 
               * parsing the actual data for the current column. */
              fetchData->parserState = PARSING_DATA_QUOTED;
            }
        }
    }

  if (!fetchData->currentColumn)
    {
      /* No column was found, so ignore this tag. Reset the tag name to
       * zero-length so we know to start again */
      fetchData->parserState = PARSING_IGNORE;
      g_string_truncate(fetchData->tagName, 0);
    }
}


/* Append the given char to the given GString, but only if it is not an end quote. Assumes the
 * start quote has already been seen. Treats double quotes ("") as a single escaped quote and
 * includes it in the text - foundEndQuote is set for the first quote found and if the next char
 * this is called with is also a quote the quote is included and foundEndQuote reset.
 * To do: this should really identify newlines and ignore the FT identifier and whitespace at the
 * start of the line. Since none of our tags currently have multi-line data that's not worth doing just yet... */
static void appendCharToQuotedString(const char curChar, gboolean *foundEndQuote, GString *result)
{
  if (curChar == '"')
    {
      if (*foundEndQuote)
        {
          /* The previous char was also a quote, which means this is an 
           * escaped quote, so we do include the char in the text. */
          appendCharToString(curChar, result);
          *foundEndQuote = FALSE;
        }
      else
        {
          /* This looks like an end quote but we won't know until the next loop, so just flag it. */
          *foundEndQuote = TRUE;
        }
    }
  else
    {
      appendCharToString(curChar, result);
    }
}


/* Append the given char to the given GString. Doesn't
 * add whitespace chars to the start of the string. Also 
 * filters out multiple whitespace chars. */
static void appendCharToString(const char curChar, GString *result)
{
  if (result)
    {
      /* Don't add multiple consecutive whitespace characters, or
       * white space at the start of the string */
      const int len = result->len;
      const char *str = result->str;

      if (!isWhitespaceChar(curChar) || (len > 0 && !isWhitespaceChar(str[len - 1])))
        {
          g_string_append_c(result, curChar);
        }
    }
}


/* This is called when we've finished parsing a given sequence. It checks that the sequence
 * data is valid (i.e. not an error message) and complements it if necessary. It updates the parser
 * state to finished if we've got all the sequences we requested. Returns true if the pfetch
 * was successful, false if not */
static gboolean pfetchFinishSequence(GeneralFetchData *fetchData)
{
  /* We issue warnings if a sequence fetch failed but set a limit on the number so we don't end
   * up with pages and pages of terminal output */
  static int numWarnings = 0;
  const int maxWarnings = 50;

  fetchData->numFetched += 1;
  
  if (fetchData->numFetched >= fetchData->numRequested)
    {
      /* We've fetched all of the sequences we requested. */
      fetchData->parserState = PARSING_FINISHED;
    }
  
  /* The pfetch failed if our sequence is null or equal to an error string. */
  gboolean pfetch_ok = FALSE;
  
  /* If the sequence isn't required, or is already set, then we don't need 
   * to check if it was fetched */
  if (!blxSequenceRequiresColumnData(fetchData->currentSeq, BLXCOL_SEQUENCE) ||
      blxSequenceGetValueAsString(fetchData->currentSeq, BLXCOL_SEQUENCE))
    {
      pfetch_ok = TRUE;
    }
  else 
    {
      if (!fetchData->currentResult || fetchData->currentResult->len == 0)
        {
          if (numWarnings < maxWarnings)
            {
              ++numWarnings;
              g_warning("No sequence data fetched for '%s'\n", blxSequenceGetName(fetchData->currentSeq));
            }
        }
      else if (stringInArray(fetchData->currentResult->str, fetchData->fetchMethod->errors))
        {
          if (numWarnings < maxWarnings)
            {
              ++numWarnings;
              g_warning("Sequence fetch failed for '%s': %s\n", blxSequenceGetName(fetchData->currentSeq), fetchData->currentResult->str);
            }
        }
      else if (!isValidIupacChar(fetchData->currentResult->str[0], fetchData->seqType))
        {
          if (numWarnings < maxWarnings)
            {
              ++numWarnings;
              g_warning("Sequence fetch failed for '%s': invalid character '%c' found in fetched text: %s\n", blxSequenceGetName(fetchData->currentSeq), fetchData->currentResult->str[0], fetchData->currentResult->str);
            }
        }
      else
        {
          /* Succeeded: save the result in the blxsequence */
          pfetch_ok = TRUE;
          blxSequenceSetValueFromString(fetchData->currentSeq, BLXCOL_SEQUENCE, fetchData->currentResult->str);
        }
    }
  
  if (pfetch_ok)
    {
      fetchData->numSucceeded += 1;
    }

  /* Reset the result string */
  g_string_truncate(fetchData->currentResult, 0);
  
  return pfetch_ok;
}


/*******************************************************************
 *                          Bulk fetch                             *
 *******************************************************************/

/* Returns true if the given fetch method retrieves sequence 
 * data. Note that fetch methods that return gff files for 
 * re-parsing will cause this function to return true. */
static gboolean fetchMethodReturnsSequence(const BlxFetchMethod* const fetchMethod)
{
  gboolean result = FALSE;
  
  if (fetchMethod)
    {
      result = 
        fetchMethod->outputType == BLXFETCH_OUTPUT_RAW ||
        fetchMethod->outputType == BLXFETCH_OUTPUT_FASTA ||
        fetchMethod->outputType == BLXFETCH_OUTPUT_EMBL || 
        fetchMethod->outputType == BLXFETCH_OUTPUT_GFF;
    }
  
  return result;
}


/* Returns true if the given fetch method retrieves data for
 * optional columns. Note that fetch methods that return gff files for 
 * re-parsing will cause this function to return true. */
static gboolean fetchMethodReturnsOptionalColumns(const BlxFetchMethod* const fetchMethod)
{
  gboolean result = FALSE;
  
  if (fetchMethod)
    {
      /* - EMBL files can be parsed for optional columns like tissue-type and strain.
       * - DB fetch methods that return a list of columns can also return optional columns.
       * - GFF files will be re-parsed and the results will be checked again, so we return
       * true for those too. */
      result = 
        fetchMethod->outputType == BLXFETCH_OUTPUT_EMBL || 
        fetchMethod->outputType == BLXFETCH_OUTPUT_LIST || 
        fetchMethod->outputType == BLXFETCH_OUTPUT_GFF;
    }
  
  return result;
}


static const BlxFetchMethod* findFetchMethod(const BlxFetchMode mode, 
                                             const BlxFetchOutputType outputType,
                                             GHashTable *fetchMethods)
{
  const BlxFetchMethod* result = NULL;

  GHashTableIter iter;
  gpointer key, value;

  g_hash_table_iter_init (&iter, fetchMethods);
  while (g_hash_table_iter_next (&iter, &key, &value) && !result)
    {
      const BlxFetchMethod* fetchMethod = (const BlxFetchMethod*)value;
      
      if (fetchMethod->mode == mode && fetchMethod->outputType == outputType)
        result = fetchMethod;
    }

  return result;
}


/* Checks the list of sequences for blixem to display to see which ones
 * need fetching.
 * Returns lists of all the sequences to be fetched, categorised by the fetch 
 * mode that should be used to fetch them. The return value is a map of
 * a GQuark (representing the fetch-mode string) to the GList of sequences
 * to be fetched.
 * This function can be called multiple times on the same sequences to
 * re-try fetching sequences with different fetch methods if the original
 * fetch method fails. Pass 'attempt' as 0 for the first try, 1 for 
 * the second etc. 
 * If optionalColumns is true, then this 'forces' optional data to be 
 * loaded even if the bulk fetch method for a sequence does not return
 * embl data. We do this by looking for another fetch method of the same
 * mode that does return embl data. This is perhaps a bit hacky but avoids 
 * us having to have yet another set of fetch methods in the config for the
 * user-triggered load-optional-data action. */
static GHashTable* getSeqsToPopulate(GList *inputList, 
                                     const GArray *defaultFetchMethods,
                                     const int attempt,
                                     GHashTable *fetchMethods,
                                     const gboolean optionalColumns)
{
  GHashTable *resultTable = g_hash_table_new(g_direct_hash, g_direct_equal);
  
  /* Loop through the input list */
  GList *inputItem = inputList;
  
  for ( ; inputItem; inputItem = inputItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(inputItem->data);
      GQuark fetchMethodQuark = blxSequenceGetFetchMethod(blxSeq, TRUE, optionalColumns, attempt, defaultFetchMethods);

      if (fetchMethodQuark)
        {
          const BlxFetchMethod* fetchMethod = getFetchMethodDetails(fetchMethodQuark, fetchMethods);

          if (fetchMethod)
            {
              /* Check if sequence data is required and is not already set.
               * Also only attempt to fetch the sequence if this fetch method
               * can return it! */
              gboolean getSeq = (blxSequenceRequiresSeqData(blxSeq) && 
                                 fetchMethodReturnsSequence(fetchMethod) &&
                                 !blxSequenceGetSequence(blxSeq));
              
              /* Check if full embl data data is required and is not already set. */
              gboolean getEmbl = (blxSequenceRequiresOptionalData(blxSeq) &&
                                  (!blxSequenceGetOrganism(blxSeq) ||
                                   !blxSequenceGetGeneName(blxSeq) ||
                                   !blxSequenceGetTissueType(blxSeq) ||
                                   !blxSequenceGetStrain(blxSeq)));
              
              /* Only attempt to fetch the embl data if this fetch method can
               * return it, or, if optionalColumns is true, then search for a 
               * fetch method that does return embl data */ 
              if (getEmbl)
                {
                  if (!fetchMethodReturnsOptionalColumns(fetchMethod) && optionalColumns)
                    {
                      fetchMethod = findFetchMethod(fetchMethod->mode, BLXFETCH_OUTPUT_EMBL, fetchMethods);

                      if (fetchMethod)
                        fetchMethodQuark = fetchMethod->name;
                      else
                        fetchMethodQuark = 0;
                    }

                  getEmbl = fetchMethodReturnsOptionalColumns(fetchMethod);
                }
              
              getSeq |= getEmbl;

              if (getSeq)
                {
                  /* Get the result list for this fetch method. It's ok if it is 
                   * null because the list will be created by g_list_prepend. */
                  GList *resultList = (GList*)g_hash_table_lookup(resultTable, GINT_TO_POINTER(fetchMethodQuark));
                  resultList = g_list_prepend(resultList, blxSeq);
                  
                  /* Update the existing (or insert the new) list */
                  g_hash_table_insert(resultTable, GINT_TO_POINTER(fetchMethodQuark), resultList);
                }
            }
          else
            {
              /* Warn user that fetch method details are missing (only warn once per fetch method) */
              static GHashTable *seenMethods = NULL;

              if (!seenMethods)
                seenMethods = g_hash_table_new(g_direct_hash, g_direct_equal);

              if (!g_hash_table_lookup(seenMethods, GINT_TO_POINTER(fetchMethodQuark)))
                {
                  g_warning("Error fetching sequence '%s'; could not find details for fetch method '%s'\n", blxSequenceGetName(blxSeq), g_quark_to_string(fetchMethodQuark));
                  g_hash_table_insert(seenMethods, GINT_TO_POINTER(fetchMethodQuark), GINT_TO_POINTER(fetchMethodQuark));
                }
            }
        }
    }

  return resultTable;
}


void sendFetchOutputToFile(GString *command,
                           GKeyFile *keyFile,
                           BlxBlastMode *blastMode,
                           GArray* featureLists[],
                           GSList *supportedTypes, 
                           GSList *styles,
                           GList **seqList,
                           MSP **mspListIn,
                           const char *fetchName,
                           const gboolean saveTempFiles,
                           MSP **newMsps,
                           GList **newSeqs,
                           GList *columnList,
                           GHashTable *lookupTable,
                           const int refSeqOffset,
                           const IntRange* const refSeqRange,
                           GError **error)
{
  /* Create a temp file for the results */
  const char *tmpDir = getSystemTempDir();
  char *fileName = g_strdup_printf("%s/%s_%s", tmpDir, MKSTEMP_CONST_CHARS_GFF, MKSTEMP_REPLACEMENT_CHARS);
  int fileDesc = mkstemp(fileName);
  GError *tmpError = NULL;
  
  if (!fileName || fileDesc == -1)
    {
      g_set_error(&tmpError, BLX_ERROR, 1, "  %s: Error creating temp file for fetch results (filename=%s)\n", fetchName, fileName);
    }

  if (fileDesc != -1)
    {
      close(fileDesc);
    }

  if (!tmpError)
    {
      /* Send the output to the temp file */
      g_string_append_printf(command, " > %s", fileName);
      
      FILE *outputFile = fopen(fileName, "w");
      
      g_debug("Fetch command:\n%s\n", command->str);
      
      GtkWidget *dialog = gtk_message_dialog_new(NULL, GTK_DIALOG_MODAL, GTK_MESSAGE_INFO, GTK_BUTTONS_NONE, "Fetching features...");
      gtk_widget_show_all(dialog);

      g_message_info("Executing fetch...\n");
      const gboolean success = (system(command->str) == 0);

      gtk_widget_destroy(dialog);
      
      fclose(outputFile);

      if (success)
        {
          /* Parse the results */
          g_message_info("... ok.\n");
          g_message_info("Parsing fetch results...");

          loadNativeFile(fileName, NULL, keyFile, blastMode, featureLists, 
                         supportedTypes, styles, newMsps, newSeqs, columnList, 
                         lookupTable, refSeqOffset, refSeqRange, &tmpError);

          if (!tmpError)
            {
              g_message_info("... ok.\n");
            }
          else
            {
              g_message_info("... failed.\n");
            }
        }
      else
        {
          g_message_info("... failed.\n");
          g_set_error(&tmpError, BLX_ERROR, 1, "  %s: Command failed.\n", fetchName);
        }
    }

  /* Delete the temp file (unless the 'save temp files' option is on) */
  if (!saveTempFiles)
    {
      if (unlink(fileName) != 0)
        g_warning("Error removing temp file '%s'.\n", fileName);
    }

  g_free(fileName);

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Called by regionFetchList. Fetches the sequences for a specific region.
 *   tmpDir is the directory in which to place the temporary files
 *   script is the script to call to do the fetch
 *   dataset will be passed as the -dataset argument to the script if it is not null */
static void regionFetchFeature(const MSP* const msp, 
                               const BlxSequence* const blxSeq,
                               const BlxFetchMethod* const fetchMethod,
                               const char *script,
                               const char *dataset,
                               const char *tmpDir,
                               const int refSeqOffset,
                               BlxBlastMode *blastMode,
                               GList **seqList,
                               GList *columnList,
                               MSP **mspListIn,
                               GArray* featureLists[],
                               GSList *supportedTypes, 
                               GSList *styles,
                               const gboolean saveTempFiles,
                               const IntRange* const refSeqRange,
                               GHashTable *lookupTable,
                               GError **error)
{
  GKeyFile *keyFile = blxGetConfig();
  const char *fetchName = g_quark_to_string(fetchMethod->name);
  GError *tmpError = NULL;

  /* Get the command string, including the args */
  GString *command = getFetchCommand(fetchMethod, NULL, 
                                     msp, mspGetRefName(msp), 
                                     refSeqOffset, refSeqRange, 
                                     dataset, &tmpError);

  if (tmpError)
    prefixError(tmpError, "  %s: Error constructing fetch command:\n", fetchName);
  
  if (!tmpError)
    {
      MSP *newMsps  = NULL;
      GList *newSeqs = NULL;
      
      sendFetchOutputToFile(command, keyFile, blastMode, 
                            featureLists, supportedTypes, styles, 
                            seqList, mspListIn, fetchName, saveTempFiles, 
                            &newMsps, &newSeqs, columnList, lookupTable, refSeqOffset, refSeqRange, &tmpError);

      blxMergeFeatures(newMsps, newSeqs, mspListIn, seqList);
    }
  
  if (command)
    g_string_free(command, TRUE);

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Fetch sequences for a given region. This uses the config file to find the
 * script and arguments to call to fetch the sequences.
 * The input GList contains a list of BlxSequences that are parent objects for
 * MSPs that identify regions. For each region, the script is called to fetch
 * all sequences that lie within that region, and the results are placed in 
 * a temporary GFF file, which is then parsed to get the results. The GFF file
 * is deleted when finished, unless the saveTempFiles argument is true. */
static void regionFetchList(GList *regionsToFetch, 
                            GList **seqList, 
                            GList *columnList,
                            const BlxFetchMethod* const fetchMethod, 
                            MSP **mspListIn,
                            BlxBlastMode *blastMode,
                            GArray* featureLists[],
                            GSList *supportedTypes, 
                            GSList *styles,
                            gboolean External,
                            const gboolean saveTempFiles,
                            const BlxSeqType seqType,
                            const int refSeqOffset,
                            const char *dataset,
                            const IntRange* const refSeqRange,
                            GHashTable *lookupTable,
                            GError **error)
{
  /* Get the command to run */
  const gchar *script = fetchMethod->location;

  if (!script)
    {
      g_set_error(error, BLX_ERROR, 1, "Error fetching sequences; no command given for fetch-method [%s].\n", g_quark_to_string(fetchMethod->name));
      return;
    }

  const char *tmpDir = getSystemTempDir();

  /* Loop through each region, creating a GFF file with the results for each region */
  GList *regionItem = regionsToFetch;
  GError *tmpError = NULL;

  for ( ; regionItem && !tmpError; regionItem = regionItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(regionItem->data);
      GList *mspItem = blxSeq->mspList;
    
      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);

          /* Only fetch regions that are at least partly inside our display range */
          if (!rangesOverlap(&msp->qRange, refSeqRange))
            continue;
          
          regionFetchFeature(msp, blxSeq, fetchMethod, script, dataset, tmpDir, refSeqOffset,
                             blastMode, seqList, columnList, mspListIn, featureLists, supportedTypes,
                             styles, saveTempFiles, refSeqRange, lookupTable, &tmpError);
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Fetch sequences using a given command-line script */
static void commandFetchList(GList *regionsToFetch, 
                             GList **seqList, 
                             GList *columnList,
                             const BlxFetchMethod* const fetchMethod, 
                             MSP **mspListIn,
                             BlxBlastMode *blastMode,
                             GArray* featureLists[],
                             GSList *supportedTypes, 
                             GSList *styles,
                             gboolean External,
                             const gboolean saveTempFiles,
                             const BlxSeqType seqType,
                             const int refSeqOffset,
                             const char *dataset,
                             const IntRange* const refSeqRange,
                             GHashTable *lookupTable,
                             GError **error)
{
  /* Currently we only support an output type of gff */
  if (fetchMethod->outputType == BLXFETCH_OUTPUT_GFF)
    {
      regionFetchList(regionsToFetch, seqList, columnList, fetchMethod, mspListIn,
                      blastMode, featureLists, supportedTypes, styles, 
                      External, saveTempFiles, seqType, refSeqOffset,
                      dataset, refSeqRange, lookupTable, error);
    }
  else
    {
      g_warning("Fetch mode 'command' currently only supports an output type of 'gff'\n");
    }
}


/* This function determines which function to call to fetch the 
 * given list of features, based on the given fetch method.
 * The first argument is the list  of features to fetch and 
 * the second is the list of all features. */
static gboolean fetchList(GList *seqsToFetch, 
                          GList **seqList,
                          GList *columnList,
                          const BlxFetchMethod* const fetchMethod,
                          const BlxSeqType seqType,
                          const gboolean saveTempFiles,
                          const gboolean External,
                          MSP **mspList,
                          BlxBlastMode *blastMode,
                          GArray* featureLists[],
                          GSList *supportedTypes, 
                          GSList *styles,
                          const int refSeqOffset,
                          const IntRange* const refSeqRange,
                          const char *dataset,
                          GHashTable *lookupTable,
                          GError **error)
{
  gboolean success = TRUE;
  
  if (g_list_length(seqsToFetch) > 0)
    {
      if (fetchMethod)
        {
          if (fetchMethod->mode == BLXFETCH_MODE_SOCKET)
            {  
              success = socketFetchList(seqsToFetch,
                                        fetchMethod,
                                        *seqList,
                                        columnList,
                                        External,
                                        seqType,
                                        error);
            }
#ifdef PFETCH_HTML 
          else if (fetchMethod->mode == BLXFETCH_MODE_HTTP || fetchMethod->mode == BLXFETCH_MODE_PIPE)
            {
              success = httpFetchList(seqsToFetch,
                                      fetchMethod,
                                      *seqList,
                                      columnList,
                                      seqType,
                                      error);
            }
#endif
          else if (fetchMethod->mode == BLXFETCH_MODE_SQLITE)
            {
              sqliteFetchSequences(seqsToFetch, fetchMethod, columnList, error);
            }
          else if (fetchMethod->mode == BLXFETCH_MODE_COMMAND)
            {
              commandFetchList(seqsToFetch, 
                               seqList, 
                               columnList,
                               fetchMethod, 
                               mspList,
                               blastMode,
                               featureLists,
                               supportedTypes,
                               styles,
                               External,
                               saveTempFiles,
                               seqType,
                               refSeqOffset,
                               dataset,
                               refSeqRange,
                               lookupTable,
                               error);
            }
          else if (fetchMethod->mode == BLXFETCH_MODE_NONE)
            {
              /* do nothing */
            }
          else
            {
              g_set_error(error, BLX_ERROR, 1, "Bulk fetch is not implemented yet in %s mode.\n", g_quark_to_string(fetchMethod->name));
            }
        }
      else
        {
          g_set_error(error, BLX_ERROR, 1, "Fetch mode not specified.\n");
        }
    }
  
  return success;
}


/* Find out if we need to fetch any sequences (they may all be 
 * contained in the input files so there might not be anything to 
 * fetch). If we do need to, then fetch them by the preferred method.
 * If the preferred fetch method fails, recusively try any other
 * fetch methods set up for each sequence until we have either fetched 
 * everything or run out of fetch methods to try.
 * 'attempt' should be passed as 0 for the first call. */
gboolean bulkFetchSequences(const int attempt,
                            gboolean External, 
                            const gboolean saveTempFiles,
                            const BlxSeqType seqType,
                            GList **seqList, /* list of BlxSequence structs for all required sequences */
                            GList *columnList,
                            const GArray *defaultFetchMethods,
                            GHashTable *fetchMethods,
                            MSP **mspList,
                            BlxBlastMode *blastMode,
                            GArray* featureLists[],
                            GSList *supportedTypes, 
                            GSList *styles,
                            const int refSeqOffset,
                            const IntRange* const refSeqRange,
                            const char *dataset,
                            const gboolean optionalColumns,
                            GHashTable *lookupTable)
{
  gboolean success = FALSE; /* will get set to true if any of the fetch methods succeed */
  
  /* Fetch any sequences that do not have their sequence data
   * already populated. If this is a re-try attempt, then use
   * a secondary fetch method, if one is given; otherwise, exclude
   * from the list (i.e. when we run out of fetch methods or everything
   * has been successfully fetched, then this table will be empty). */
  GHashTable *seqsTable = getSeqsToPopulate(*seqList, defaultFetchMethods, attempt, fetchMethods, optionalColumns);

  if (g_hash_table_size(seqsTable) < 1)
    {
      /* Nothing to fetch */
      success = TRUE;
    }
  else
    {
      /* Loop through each fetch mode */
      GHashTableIter iter;
      g_hash_table_iter_init(&iter, seqsTable);
      gpointer key, value;
      GError *error = NULL;
      
      while (g_hash_table_iter_next(&iter, &key, &value))
        {
          GList *seqsToFetch = (GList*)value;          
          GQuark fetchMethodQuark = GPOINTER_TO_INT(key);

          if (!fetchMethodQuark)
            continue;
          
          g_message_info("Fetching %d items using method '%s' (attempt %d)\n", g_list_length(seqsToFetch), g_quark_to_string(fetchMethodQuark), attempt + 1);

          const BlxFetchMethod* const fetchMethod = getFetchMethodDetails(fetchMethodQuark, fetchMethods);

          if (!fetchMethod)
            {
              g_warning("Fetch method '%s' not found\n", g_quark_to_string(fetchMethodQuark));
              continue;
            }
          
          GError *tmpError = NULL;
          
          if (fetchList(seqsToFetch, seqList, columnList, fetchMethod, seqType, saveTempFiles, External, mspList, blastMode, featureLists, supportedTypes, styles, refSeqOffset, refSeqRange, dataset, lookupTable, &tmpError))
            {
              success = TRUE;
              
              /* Compile all errors into a single error */
              if (error && tmpError)
                {
                  prefixError(error, tmpError->message);
                  g_error_free(tmpError);
                  tmpError = NULL;
                }
              else if (tmpError)
                {
                  error = tmpError;
                }
            }
          
          /* We're done with this list now, so free the memory. Don't delete it from
           * the table yet, though, because that will invalidate the iterators. */
          g_list_free(seqsToFetch);
        }
      
      if (success && error)
        {
          /* Some fetches succeeded, so just issue a warning */
          prefixError(error, "Error fetching sequences:\n");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
      else if (error)
        {
          /* All failed, so issue a critical warning */
          prefixError(error, "Error fetching sequences:\n");
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }

      /* Recurse to re-try any that failed. This returns straight
       * away everything was fetched successfully, or if there are
       * no more fetch methods to try. */
      success = bulkFetchSequences(attempt + 1, External, saveTempFiles, seqType, seqList, columnList,
                                   defaultFetchMethods, fetchMethods, mspList,
                                   blastMode, featureLists, supportedTypes, 
                                   styles, refSeqOffset, refSeqRange, dataset, optionalColumns, lookupTable);
    }

  /* Clean up */
  g_hash_table_unref(seqsTable);

  return success;
}

