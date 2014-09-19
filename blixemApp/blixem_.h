/*  File: blixem_.h
 *  Author: Ed Griffiths, 2001-11-29
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
 * Description: Internal header for Blixem package-wide code
 *----------------------------------------------------------------------------
 */

#ifndef DEF_BLIXEM_P_H
#define DEF_BLIXEM_P_H

#include <gtk/gtk.h>
#include <config.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/version.h>


/*            blixem program version and information.                        */
#define BLIXEM_TITLE   "Blixem"
#define BLIXEM_PREFIX  BLIXEM_TITLE" - "  /* The prefix to use for window titles etc. */
#define BLIXEM_PREFIX_ABBREV  "B: "       /* An abbreviated version of the prefix to save space e.g. good for displaying in the taskbar */
#define BLIXEM_DESC    "Multiple alignment visualisation tool."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define BLIXEM_VERSION_STRING      SEQTOOLS_VERSION
#define BLIXEM_PACKAGE_VERSION     UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, SEQTOOLS_VERSION)
#define BLIXEM_TITLE_STRING        UT_MAKE_TITLE_STRING(BLIXEM_TITLE, BLIXEM_VERSION_STRING)
#define BLIXEM_VERSION_COMPILE     BLIXEM_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define BLIXEM_COPYRIGHT_STRING    UT_MAKE_COPYRIGHT_STRING("2009-2011")
#define BLIXEM_WEBSITE_STRING      "http://www.sanger.ac.uk/resources/software/seqtools/"
#define BLIXEM_LICENSE_STRING      UT_MAKE_LICENCE_STRING(BLIXEM_TITLE)

#define BLIXEM_COMMENTS_STRING()                                \
"("BLIXEM_TITLE_STRING", "                                      \
UT_COMPILE_PHRASE " " UT_MAKE_COMPILE_DATE() ")\n\n"            \
AUTHOR_TEXT "\n"




#define COLUMN_WIDTHS_GROUP             "column-widths"   /* group name in the config file */
#define COLUMN_SUMMARY_GROUP            "summary-columns" /* group name in the config file for
                                                             columns to include in the summary info */

/* Define the columns' default widths and titles. */
#define BLXCOL_DEFAULT_WIDTH            50    /* default width for generic columns */
#define BLXCOL_SEQNAME_WIDTH            120   /* default width for the name column */
#define BLXCOL_SCORE_WIDTH              40    /* default width for the score column */
#define BLXCOL_ID_WIDTH                 45    /* default width for the ID column */
#define BLXCOL_SOURCE_WIDTH             85    /* default width for source column  */
#define BLXCOL_GROUP_WIDTH              58    /* default width for group column  */
#define BLXCOL_START_WIDTH              50    /* default width for the start coord column */
#define BLXCOL_END_WIDTH                80    /* default width for end coord column (bigger because it also spans the scrollbar) */
#define BLXCOL_SEQUENCE_WIDTH           40    /* default width for sequence column */
#define BLXCOL_ORGANISM_WIDTH           28    /* default width for organism column */
#define BLXCOL_GENE_NAME_WIDTH          58    /* default width for gene-name column  */
#define BLXCOL_STRAIN_WIDTH             85    /* default width for strain column  */
#define BLXCOL_TISSUE_TYPE_WIDTH        100   /* default width for tissue-type column  */


/* 
 * config file groups/keywords, these should not be changed willy nilly as they
 * are used external programs and users when constructing config files.
 * If you change things or add additional groups or values, then you should
 * update the getConfig function.
 */

/* For overall blixem settings. */
#define BLIXEM_OLD_BULK_FETCH      "default-fetch-mode"      /* for compatibility with old config files (new config files use SEQTOOLS_BULK_FETCH) */

/* Fetch settings */
#define FETCH_MODE_KEY             "fetch-mode"  /* any group with this key is a fetch method, and this specifies what type of fetch to do */

/* These are the supported fetch modes. ***If you add anything here, also add it in fetchModeStr*** */
typedef enum
  {
#ifdef PFETCH_HTML 
    BLXFETCH_MODE_HTTP,
    BLXFETCH_MODE_PIPE,
#endif
    BLXFETCH_MODE_SOCKET,
    BLXFETCH_MODE_WWW,
    BLXFETCH_MODE_SQLITE,
    BLXFETCH_MODE_COMMAND,
    BLXFETCH_MODE_INTERNAL,
    BLXFETCH_MODE_NONE,

    BLXFETCH_NUM_MODES /* must be last in list */
  } BlxFetchMode;


/* Required keys for http-and pipe-fetch groups */
#ifdef PFETCH_HTML
#define HTTP_FETCH_LOCATION       "url"
#define HTTP_FETCH_PORT           "port"
#define HTTP_FETCH_ARGS           "request"
#define HTTP_FETCH_COOKIE_JAR     "cookie-jar"

#define PIPE_FETCH_LOCATION       "command"
#define PIPE_FETCH_ARGS           "args"
#endif 

/* Required keys for socket-fetch groups */
#define SOCKET_FETCH_LOCATION     "command"
#define SOCKET_FETCH_NODE         "node"
#define SOCKET_FETCH_PORT         "port"
#define SOCKET_FETCH_ARGS         "args"

/* Required keys for www-fetch groups */
#define WWW_FETCH_LOCATION       "url"
#define WWW_FETCH_ARGS           "request"

/* Required keys for db-fetch groups */
#define DB_FETCH_LOCATION        "location"
#define DB_FETCH_QUERY           "query"
#define DB_FETCH_COLUMNS         "columns"

/* Required keys for command-fetch groups */
#define COMMAND_FETCH_SCRIPT        "command"
#define COMMAND_FETCH_ARGS          "args"


/* Generic fetch keys (applicable to more than one method) */
#define FETCH_OUTPUT       "output"      /* output format */
#define FETCH_SEPARATOR    "separator"   /* separator, when combining multiple sequences */
#define FETCH_ERRORS       "errors"      /* list of messages that indicate errors */


/* For settings */
#define BLIXEM_SETTINGS_FILE             ".blixemrc"  /* default file name for saving blixem settings to */
#define SETTINGS_GROUP                   "user-settings"
#define STYLES_FILE_KEY                  "stylesfile" /* styles-file key in the [blixem] group */  

#define SETTING_NAME_INVERT_SORT "invert-sort"
#define SETTING_NAME_HIGHLIGHT_DIFFS "highlight-diffs"
#define SETTING_NAME_HIGHLIGHT_VARIATIONS "highlight-variations"
#define SETTING_NAME_SHOW_VARIATION_TRACK "show-variations-track"
#define SETTING_NAME_SHOW_UNALIGNED "show-unaligned"
#define SETTING_NAME_SHOW_UNALIGNED_SELECTED "show-unaligned-selected-seq"
#define SETTING_NAME_LIMIT_UNALIGNED_BASES "limit-unaligned"
#define SETTING_NAME_SHOW_POLYA_SITE "show-polya-site"
#define SETTING_NAME_SHOW_POLYA_SITE_SELECTED "show-poly-site-selected-seq"
#define SETTING_NAME_SHOW_POLYA_SIG "show-poly-sig"
#define SETTING_NAME_SHOW_POLYA_SIG_SELECTED "show-polya-sig-selected-seq"
#define SETTING_NAME_SHOW_SPLICE_SITES "show-splice-sites"
#define SETTING_NAME_SQUASH_MATCHES "squash-matches"
#define SETTING_NAME_SHOW_COLINEARITY "show-colinearity"
#define SETTING_NAME_SHOW_COLINEARITY_SELECTED "show-colinearity-selected"


/* would be good to get rid of this.... */
#define FULLNAMESIZE               255


#define MKSTEMP_CONST_CHARS_GFF               "BLIXEM_gff"        /* the prefix to use when creating a temp file name */
#define MKSTEMP_REPLACEMENT_CHARS             "XXXXXX"            /* the required string that will be replaced by unique chars when creating a temp file name */




/* Blixem config/fetch error domain */
#define BLX_CONFIG_ERROR g_quark_from_string("Blixem config")
#define BLX_FETCH_ERROR g_quark_from_string("Blixem config")

/* Error codes */
typedef enum
  {
    BLX_CONFIG_ERROR_NO_GROUPS,             /* no groups in config file */
    BLX_CONFIG_ERROR_INVALID_KEY_TYPE,      /* invalid key type given in config file */
    BLX_CONFIG_ERROR_MISSING_KEY,           /* a required key is missing */
    BLX_CONFIG_ERROR_INVALID_FETCH_MODE,    /* invalid fetch mode specified */
    BLX_CONFIG_ERROR_INVALID_OUTPUT_FORMAT, /* invalid output format specified for fetch mode */
    BLX_CONFIG_ERROR_WARNINGS,              /* warnings found while reading config file */
    BLX_CONFIG_ERROR_SUBSTITUTION,          /* error with substitution string */
    BLX_CONFIG_ERROR_INVALID_FETCH_METHOD,  /* null fetch method */
    BLX_CONFIG_ERROR_NO_EXE,                /* fetch method executable does not exist */
    BLX_CONFIG_ERROR_NULL_FETCH,            /* fetch method is null */
    BLX_CONFIG_ERROR_NO_ARGS               /* mandatory args weren't specified */
  } BlxConfigError;


/* Error codes */
typedef enum
  {
    BLX_FETCH_ERROR_SOCKET,                /* error creating socket */
    BLX_FETCH_ERROR_HOST,                  /* unknown host */
    BLX_FETCH_ERROR_CONNECT,               /* error connecting to host */
    BLX_FETCH_ERROR_SEND,                  /* error sending to socket */
    BLX_FETCH_ERROR_PATH                   /* executable not found in path */
  } BlxFetchError;



/* Function pointer to a function that performs a fetch of a particular sequence */
typedef void(*FetchFunc)(const char *seqName, gpointer fetchMethod, const gboolean bulk, GtkWidget *blxWindow);


/* output types for fetch modes. *** if you add anything here, also add it in outputTypeStr *** */
typedef enum
{
  BLXFETCH_OUTPUT_INVALID,
  BLXFETCH_OUTPUT_RAW,      /* raw sequence data, separated by newlines */
  BLXFETCH_OUTPUT_FASTA,    /* sequence data in FASTA format */
  BLXFETCH_OUTPUT_EMBL,     /* the sequence's EMBL entry */
  BLXFETCH_OUTPUT_LIST,     /* a list of named columns is returned */
  BLXFETCH_OUTPUT_GFF,      /* a new gff for re-parsing is returned */

  BLXFETCH_NUM_OUTPUT_TYPES
} BlxFetchOutputType;


/* struct to hold info about a fetch method */
typedef struct _BlxFetchMethod
{
  GQuark name;                      /* fetch method name */
  BlxFetchMode mode;                /* the type of fetch method */
  
  char *location;                   /* e.g. url, script, command, db location etc. */
  char *node;                       /* for socket fetch mode */
  int port;                         /* for socket and http/pipe fetch modes */
  char *cookie_jar;                 /* for http/pipe fetch mode */
  char *args;                       /* arguments/query/request */
  GArray *columns;                  /* for db-fetch, the list of columns the query will populate */

  char *separator;                  /* separator when combining multiple sequence names into a list */
  GArray *errors;                   /* array of messages (as GQuarks) that indicate that an error occurred, e.g. "no match" */
  BlxFetchOutputType outputType;    /* the output format to expect from the fetch command */
} BlxFetchMethod;


/* The following are used to define default colors for certain types of features in Blixem.
 * One of several different actual colors from the BlxColor struct may be used depending 
 * on state, e.g. we use a different color if "print colors" (i.e. black and 
 * white mode) is on. */

typedef enum 
  {
    BLXCOLOR_MIN,           /* dummy value so that we don't get a zero ID */
  
    BLXCOLOR_BACKGROUND,    /* background color of the widgets */
    BLXCOLOR_REF_SEQ,       /* default background color for the reference sequence */  
    BLXCOLOR_MATCH,         /* background color for an exact match */
    BLXCOLOR_CONS,          /* background color for a conserved match */
    BLXCOLOR_MISMATCH,      /* background color for a mismatch */
    BLXCOLOR_INSERTION,     /* color for an insertion marker */
    BLXCOLOR_EXON_START,    /* color for the start boundary line of an exon */
    BLXCOLOR_EXON_END,      /* color for the end boundary line of an exon */
    BLXCOLOR_CODON,         /* color in which to highlight the nucleotides for the currently-selected codon */
    BLXCOLOR_MET,           /* background color for MET codons in the three frame translation */
    BLXCOLOR_STOP,          /* background color for STOP codons in the three frame translation */
    BLXCOLOR_GRID_LINE,     /* color of the gridlines in the big picture grids */
    BLXCOLOR_GRID_TEXT,     /* color of the text in the big picture grids */
    BLXCOLOR_HIGHLIGHT_BOX, /* color of the highlight box in the big picture */
    BLXCOLOR_PREVIEW_BOX,   /* color of the preview box in the big picture */
    BLXCOLOR_MSP_LINE,      /* color of the MSP lines in the big picture */
    BLXCOLOR_SNP,           /* background color for SNPs */
    BLXCOLOR_GROUP,         /* default highlight color for generic groups */
    BLXCOLOR_MATCH_SET,     /* default highlight color for the special match-set group */
    BLXCOLOR_EXON_FILL,     /* fill color for an exon in the big picture */
    BLXCOLOR_EXON_LINE,     /* line color for an exon in the big picture */
    BLXCOLOR_CDS_FILL,      /* fill color for a CDS in the big picture (coding region) */
    BLXCOLOR_CDS_LINE,      /* line color for a CDS in the big picture (coding region) */
    BLXCOLOR_UTR_FILL,      /* fill color for an UTR in the big picture (non-coding/untranslated region) */
    BLXCOLOR_UTR_LINE,      /* line color for an UTR in the big picture (non-coding/untranslated region) */
    BLXCOLOR_UNALIGNED_SEQ, /* color in which to show additional sequence in the match that is not part of the alignment */
    BLXCOLOR_CANONICAL,     /* background highlight color for canonical intron bases */
    BLXCOLOR_NON_CANONICAL, /* background highlight color for non-canonical intron bases */
    BLXCOLOR_POLYA_TAIL,    /* background color for polyA tails in the detail view */
    BLXCOLOR_POLYA_SIGNAL,  /* background color for non-annotated polyA signals in the detail view */
    BLXCOLOR_POLYA_SIGNAL_ANN,/* background color for annotated polyA signals in the detail view */
    BLXCOLOR_POLYA_SITE_ANN,/* background color for annotated polyA sites in the detail view */
    BLXCOLOR_TREE_GRID_LINES,/* color of the tree grid lines (i.e. column separator lines) */
    BLXCOLOR_CLIP_MARKER,   /* color of the marker line used to indicate a match has been clipped */
    BLXCOLOR_COVERAGE_PLOT, /* color of the coverage plot */
    BLXCOLOR_ASSEMBLY_GAP,  /* highlight color for assembly gaps */
    BLXCOLOR_SELECTION,     /* highlight color for selections */
    BLXCOLOR_PARTIAL_EXON_CROSSHATCH, /* color of cross-hatch lines for partial exons */
    BLXCOLOR_COLINEAR_PERFECT, /* color of lines joining alignment blocks that are perfectly colinear */
    BLXCOLOR_COLINEAR_IMPERFECT, /* color of lines joining alignment blocks that are imperfectly colinear */
    BLXCOLOR_COLINEAR_NOT,  /* color of lines joining alignment blocks that are not colinear */

    BLXCOL_NUM_COLORS
  } BlxColorId;
  

/* Utility struct to pass around a set of colors from a styles file */
typedef struct _BlxStyleColors
{
  char *fillColor;
  char *lineColor; 
  char *fillColorSelected;
  char *lineColorSelected;
  char *fillColorCds;
  char *lineColorCds;
  char *fillColorCdsSelected;
  char *lineColorCdsSelected;
  gboolean normalFound;
  gboolean cdsFound;
} BlxStyleColors;

  
/* This enum contains a list of all the boolean options that the user can toggle on/off */
typedef enum
  {
    BLXFLAG_MIN,                    /* Start index for looping through flags */
  
    BLXFLAG_INVERT_SORT,            /* Inverts the default sort order */
    BLXFLAG_HIGHLIGHT_DIFFS,        /* Hides matching bases and highlights mis-matching ones */
    BLXFLAG_HIGHLIGHT_VARIATIONS,   /* Whether to highlight bases that have variations in the reference sequence */
    BLXFLAG_SHOW_VARIATION_TRACK,   /* Shows the Variations track */
    BLXFLAG_SHOW_UNALIGNED,         /* Shows additional bits of the match sequence that are not part of the aligned section */
    BLXFLAG_SHOW_UNALIGNED_SELECTED,/* Only show unaligned bits of sequence for the currently-selected sequence(s) */
    BLXFLAG_LIMIT_UNALIGNED_BASES,  /* If the show-unaligned-sequence option is on, limits how many bases from the unaligned sequence are shown */
    BLXFLAG_SHOW_POLYA_SITE,        /* Show polyA tails */
    BLXFLAG_SHOW_POLYA_SITE_SELECTED,/* Only show polyA tails for the currently-selected sequence(s) */
    BLXFLAG_SHOW_POLYA_SIG,         /* Show polyA signals in the reference sequence */
    BLXFLAG_SHOW_POLYA_SIG_SELECTED,/* Only show polyA signals for the currently-selected sequence(s) */
    BLXFLAG_SHOW_SPLICE_SITES,      /* Highlights splice sites in the reference sequence for the currently-selected MSPs */
    BLXFLAG_OPTIONAL_COLUMNS,       /* Gets set to true if the optional columns have been loaded */
    BLXFLAG_SHOW_CDS,               /* True if CDS/UTR regions should be shown; false if plain exons should be shown */
    BLXFLAG_NEGATE_COORDS,          /* True if coords should be negated when display is reversed (so coords appear to increase left-to-right when really they decrease) */
    BLXFLAG_HIDE_UNGROUPED_SEQS,    /* Hide all sequences that are not in a group (unless their group is also hidden) */
    BLXFLAG_HIDE_UNGROUPED_FEATURES,/* Hide all features that are not in a group (unless their group is also hidden) */
    BLXFLAG_SAVE_TEMP_FILES,        /* save any temporary files that blixem creates, e.g. the GFF file created by the region-fetch fetch mode */
    BLXFLAG_ABBREV_TITLE,           /* whether to abbreviate the window titles to save space */
    BLXFLAG_LINK_FEATURES,          /* whether features with the same name should be linked */
    BLXFLAG_SHOW_COLINEARITY,       /* whether to show colinearity lines between alignment blocks */
    BLXFLAG_SHOW_COLINEARITY_SELECTED, /* whether to show colinearity lines for selected sequence only */
    
    BLXFLAG_NUM_FLAGS               /* Total number of flags e.g. for creating arrays and loops etc */
  } BlxFlag;


/* Structure representing a group of sequences. This groups several BlxSequences together and 
 * sets various flags so that we can hide/highlight/etc all of the sequences in a group. */
typedef struct _SequenceGroup
  {
    char *groupName;               /* user-friendly name for the group (should be unique to save confusion) */
    int groupId;                   /* unique ID number for the group */
    int order;                     /* field for sorting - lower numbers will be listed first */
    GList *seqList;                /* list of BlxSequences */
    gboolean ownsSeqNames;         /* If true, the group will free the sequence names when it is destroyed */
    gboolean hidden;               /* true if the group should be hidden from the detail view */
    gboolean highlighted;          /* true if the group should be highlighted */
    GdkColor highlightColor;       /* the color to highlight the group's sequences in (in both the big picture and the detail view) */
  } SequenceGroup;


/* This enum gives a more meaningful way of indexing the "opts" string */
typedef enum
  {
    BLXOPT_MODE = 0,              /* Blast mode: L = tblastx, N = blastn, P = blastp, T = tblastn, X = blastx */
    BLXOPT_START_NEXT_MATCH = 1,  /* 'M' means start at next match. Blank to ignore. */
    BLXOPT_SHOW_BIG_PICT = 2,     /* 'B' to show big picture, 'b' to hide it */
    BLXOPT_REV_BIG_PICT = 3,      /* 'R' to show the reverse strand in the big picture; 'r' just show the forward strand (blxparser.c) */
    BLXOPT_FULL_EMBL_INFO = 4,    /* 'F' to parse full embl info on startup; blank to just parse the sequence data (quicker) */
    BLXOPT_FULL_ZOOM = 5,         /* 'Z' to use full zoom by default; blank otherwise (blxparser.c) */
    BLXOPT_INVERT_SORT = 6,       /* 'I' to invert the default sort order; blank otherwise */
    BLXOPT_HSP_GAPS = 7,          /* 'G' for HSP gaps mode; blank otherwise (blxparser.c) */
    BLXOPT_SORT_MODE = 8,         /* Initial sort mode: i=ID, n=name, p=position, s=score, t=tissue type, m=strain, g=gene name, o=organism;
                                   * OR: set to 'd' to automatically dotter the first sequence; blank otherwise (blxselect.c) */
    
    BLXOPT_NUM_OPTS               /* Total number of options (must always be last in list) */
  } BlxOptsIdx ;


/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the BlxViewContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    BLXDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */
    
    BLXDIALOG_HELP,                 /* The Help dialog */
    BLXDIALOG_SETTINGS,             /* The Settings dialog */
    BLXDIALOG_SORT,                 /* The Sort dialog */
    BLXDIALOG_FIND,                 /* The Find dialog */
    BLXDIALOG_GROUPS,               /* The Groups dialog */
    BLXDIALOG_VIEW,                 /* The View dialog */
    BLXDIALOG_DOTTER,               /* The Dotter dialog */
    
    BLXDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } BlxDialogId;



/* Struct to hold all the settings that come from the command line options */
typedef struct _CommandLineOptions
{
  char *refSeq;                   /* the section of reference sequence we're viewing */
  char *refSeqName;               /* the name of the reference sequence */
  IntRange refSeqRange;           /* the range of the reference sequence (before any offset is applied) */
  int refSeqOffset;               /* if non-zero, all parsed coords will be offset by this amount */
  int startCoord;                 /* which coord to start the initial display range at */
  gboolean startCoordSet;         /* true if the start coord has been specified on the command line */
  MSP *mspList;                   /* the list of alignments */
  GList *columnList;              /* the list of display columns */
  char **geneticCode;             /* the genetic code */
  
  BlxStrand activeStrand;         /* which strand will initially be the active one */
  gboolean negateCoords;          /* if this option is true, the display will show coords as negative when the reverse strand is active */
  gboolean zoomWhole;             /* whether to zoom out to view the entire big picture range on startup */
  int bigPictZoom;                /* initial zoom level for the big picture (as a multiple of the initial detail view range) */
  IntRange bigPictRange;          /* initial range for the big picture (will override bigPictZoom if both are set) */
  gboolean bigPictON;             /* whether to show the big picture by default */
  gboolean hideInactive;          /* whether to show the inactive strand in the big picture/detail view */
  BlxColumnId initSortColumn;     /* initial column to sort by */
  gboolean sortInverted;          /* whether initial sort order should be inverted */
  gboolean highlightDiffs;        /* whether the initial display should highlight mismatches rather than matches */
  gboolean dotterFirst;           /* open dotter when blixem starts */
  gboolean startNextMatch;        /* start at the coord of the next match from the default start coord */
  gboolean squashMatches;         /* start with the 'squash matches' option on */
  gboolean optionalColumns;       /* load data for optional columns on startup to populate additional info like tissue-type and strain */
  gboolean saveTempFiles;         /* save any temporary files that blixem creates */
  gboolean coverageOn;            /* show the coverage view on start-up */
  gboolean abbrevTitle;           /* if true, use a abbreviated window titles to save space */
  
  gboolean mspFlagDefaults[MSPFLAG_NUM_FLAGS]; /* default values for MSP flags */     
  
  BlxBlastMode blastMode;         /* the blast match mode */
  BlxSeqType seqType;             /* whether the display shows sequences as peptides or nucleotides */
  int numFrames;                  /* the number of reading frames */
  GHashTable *fetchMethods;       /* table of fetch methods (keyed on name as a GQuark) */
  GArray *bulkFetchDefault;       /* the default method(s) for bulk fetching sequences (can be overridden by an MSPs data-type properties) */
  GArray *userFetchDefault;       /* the default method(s) for fetching individual sequences interactively */
  GArray *optionalFetchDefault;   /* the default method(s) for bulk fetching optional sequence data */
  char *dataset;                  /* the name of a dataset, e.g. 'human' */
  BlxMessageData msgData;         /* data to be passed to the message handlers */
  gboolean mapCoords;             /* whether the map-coords command-line argument was specified */
  int mapCoordsFrom;              /* the coord to map from */
  int mapCoordsTo;                /* the coord to map to */
  char *windowColor;              /* if not null, set the main window background color to this */
} CommandLineOptions;



/* A Blixem View context, containing all status information required to draw the blixem view */
typedef struct _BlxViewContext
  {
    GtkWidget *statusBar;                   /* The Blixem window's status bar */
    
    char *refSeq;                           /* The reference sequence (always forward strand, always DNA sequence) */
    const char *refSeqName;                 /* The name of the reference sequence */
    IntRange refSeqRange;                   /* The range of the reference sequence */
    IntRange fullDisplayRange;              /* The range of the displayed sequence */
    int refSeqOffset;                       /* how much the coordinate system has been offset from the input coords */

    BlxBlastMode blastMode;                 /* The type of blast matching that was used */
    BlxSeqType seqType;                     /* The type of the match sequences, e.g. DNA or peptide */
    char **geneticCode;                     /* The genetic code used to translate DNA <-> peptide */
    int numFrames;                          /* The number of reading frames */

    GArray *bulkFetchDefault;               /* The default method(s) of bulk fetching sequences (can be overridden by an MSPs data-type properties) */
    GArray *userFetchDefault;               /* The default method(s) for interactively fetching individual sequences */
    GArray *optionalFetchDefault;           /* The default method(s) for bulk fetching optional sequence data */
    GHashTable *fetchMethods;               /* List of fetch methods, keyed on name as a GQuark */

    char* dataset;                          /* the name of a dataset, e.g. 'human' */
    gboolean optionalColumns;               /* load data for optional columns on startup to populate additional info like tissue-type */
    gboolean saveTempFiles;                 /* whether to save temporary files created to store results of file conversions */

    MSP *mspList;                           /* List of all MSPs. Obsolete - use featureLists array instead */
    GArray* featureLists[BLXMSP_NUM_TYPES]; /* Array indexed by the BlxMspType enum. Each array entry contains a zero-terminated array of all the MSPs of that type. */
    
    GList *matchSeqs;                       /* List of all match sequences (as BlxSequences). */
    GSList *supportedTypes;                 /* List of supported GFF types */
    const char *paddingSeq;                 /* A sequence of padding characters, used if the real sequence could not be found. All padded MSPs
                                             * use this same padding sequence - it is constructed to be long enough for the longest required seq. */
    
    gboolean displayRev;                    /* True if the display is reversed (i.e. coords decrease as you read from left to right, rather than increase). */
    gboolean external;                      /* True if Blixem was run externally or false if it was run internally from another program */
    
    GList *selectedSeqs;                    /* A list of sequences that are selected (as BlxSequences) */
    GList *sequenceGroups;                  /* A list of SequenceGroups */
    SequenceGroup *matchSetGroup;           /* A special group that can be created/deleted quickly from the 'toggle match set' shortcuts */
    
    gboolean autoDotter;                    /* Whether to use automatic dotter params */
    gboolean dotterSelf;                    /* Whether the dotter "call on self" option is on by default */
    gboolean dotterHsps;                    /* Whether the dotter "HSPs only" option is on by default */
    int dotterStart;                        /* Start coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterEnd;                          /* End coord to call dotter on, or UNSET_INT to calculate automatically */
    int dotterZoom;                         /* Zoom param to call dotter with, if using manual params */
    
    GArray *defaultColors;                  /* Default colors used by Blixem */
    gboolean usePrintColors;                /* Whether to use print colors (i.e. black and white) */
    char *windowColor;                      /* If not null, background color for the window */

    GList *columnList;                      /* A list of details about all the columns in the detail view (might have been better to use an array here but it's a short list so not important) */
    GSList *styles;
    
    gboolean flags[BLXFLAG_NUM_FLAGS];              /* Array of all the flags the user can toggle. Indexed by the BlxFlags enum. */
    GtkWidget *dialogList[BLXDIALOG_NUM_DIALOGS];   /* Array of all the persistent dialogs in the application */
    GSList *spawnedProcesses;                       /* List of processes spawned by Blixem */
    BlxModelId modelId;                             /* which tree model to use (i.e. normal or squashed) */
  
    int *depthArray;                        /* this array holds the depth (num alignments) at each coord of the ref seq */
    int minDepth;                           /* minimum value in the depthArray */
    int maxDepth;                           /* maximum value in the depthArray */
} BlxViewContext;



/* blixem can use either efetch (default) or a pfetch server to get
 * sequences, to use pfetch node/port information must be specified. */
typedef struct
  {
    char *net_id ;
    int port ;
  } PfetchParams ;


/* Used to pass around info about a match sequence when getting fetch command arguments */
typedef struct
{
  const char *match_name;
  const char *ref_name;
  int match_start;
  int match_end;
  const char *dataset;
  const char *source;
  const char *filename;
} MatchSequenceData ;


/* blxview.c */
/* Function to show blixem window, can be called from any application. */
gboolean                            blxview(CommandLineOptions *options,
                                            GArray* featureLists[],
                                            GList *seqList, 
                                            GSList *supportedTypes,
                                            PfetchParams *pfetch, 
                                            char *align_types, 
                                            gboolean External,
                                            GSList *styles,
                                            GHashTable *lookupTable) ;

BlxColumnInfo*                     getColumnInfo(GList *columnList, const BlxColumnId columnId);
int                                getColumnWidth(GList *columnList, const BlxColumnId columnId);
const char*                        getColumnTitle(GList *columnList, const BlxColumnId columnId);
void                               getColumnXCoords(GList *columnList, const BlxColumnId columnId, IntRange *xRange);
void                               saveColumnWidths(GList *columnList, GKeyFile *key_file);
void                               saveSummaryColumns(GList *columnList, GKeyFile *key_file);
gboolean                           showColumn(BlxColumnInfo *columnInfo);
void                               resetColumnWidths(GList *columnList);

void                               blviewRedraw(void);
GtkWidget*                         getBlixemWindow(void);
const IntRange*                    mspGetFullSRange(const MSP* const msp, const gboolean seqSelected, const BlxViewContext* const bc);
const IntRange*                    mspGetDisplayRange(const MSP* const msp);
const IntRange*                    mspGetFullDisplayRange(const MSP* const msp, const gboolean seqSelected, const BlxViewContext* const bc);
void                               mspCalculateFullExtents(MSP *msp, const BlxViewContext* const bc, const int numUnalignedBases);
void                               cacheMspDisplayRanges(const BlxViewContext* const bc, const int numUnalignedBases);

int                                mspGetMatchCoord(const MSP *msp, 
                                                    const int qIdx, 
                                                    const gboolean seqSelected,
                                                    const int numUnalignedBases,
                                                    BlxViewContext *bc);


void                               drawAssemblyGaps(GtkWidget *widget,
                                                    GdkDrawable *drawable,
                                                    GdkColor *color,
                                                    const gboolean displayRev,
                                                    GdkRectangle *rect, 
                                                    const IntRange* const dnaRange,
                                                    const GArray *mspArray);

GList*                             blxCreateColumns(const gboolean optionalColumns, const gboolean customSeqHeader);

GSList*                            blxReadStylesFile(const char *keyFileName_in, GError **error);

const char*                        blxGetAppName();
const char*                        blxGetTitlePrefix(const BlxViewContext* const bc);
const char*                        blxGetCopyrightString();
const char*                        blxGetWebSiteString();
const char*                        blxGetCommentsString();
const char*                        blxGetLicenseString();
const char*                        blxGetVersionString();        

/* dotter.c */
//void                               selectFeatures(void);

/* blxparser.c */
gboolean                           mspHasFs(const MSP *msp);  
char*                              readFastaSeq(FILE *seqfile, char *qname, int *startCoord, int *endCoord, const BlxSeqType seqType);

/* blxFetch.c */
BlxFetchMethod*                    getFetchMethodDetails(GQuark fetchMethodQuark, GHashTable *fetchMethods);
GString*                           getFetchCommand(const BlxFetchMethod* const fetchMethod, const BlxSequence *blxSeq, const MSP* const msp, const char *refSeqName, const int refSeqOffset, const IntRange* const refSeqRange, const char *dataset, GError **error);
GString*                           doGetFetchCommand(const BlxFetchMethod* const fetchMethod,MatchSequenceData *match_data, GError **error);
GString*                           getFetchArgs(const BlxFetchMethod* const fetchMethod, const BlxSequence *blxSeq,const MSP* const msp,const char *refSeqName,const int refSeqOffset,const IntRange* const refSeqRange,const char *dataset,GError **error);
GString*                           getFetchArgsMultiple(const BlxFetchMethod* const fetchMethod, GList *seqsToFetch, GError **error);
void                               fetchSequence(const BlxSequence *blxSeq, const gboolean displayResults, const int attempt, GtkWidget *blxWindow, GtkWidget *dialog, GtkTextBuffer **text_buffer) ;
void                               finaliseFetch(GList *seqList, GList *columnList);
void                               sendFetchOutputToFile(GString *command, GKeyFile *keyFile, BlxBlastMode *blastMode,GArray* featureLists[],GSList *supportedTypes, GSList *styles, GList **seqList, MSP **mspListIn,const char *fetchName, const gboolean saveTempFiles, MSP **newMsps, GList **newSeqs, GList *columnList, GHashTable *lookupTable, const int refSeqOffset, const IntRange* const refSeqRange, GError **error);
const char*                        outputTypeStr(const BlxFetchOutputType outputType);

void                               fetchSeqsIndividually(GList *seqsToFetch, GtkWidget *blxWindow);
gboolean                           populateSequenceDataHtml(GList *seqsToFetch, const BlxSeqType seqType, const BlxFetchMethod* const fetchMethod) ;
gboolean                           populateFastaDataPfetch(GList *seqsToFetch, BlxFetchMethod *fetchMethod, gboolean External, const BlxSeqType seqType, GError **error) ;

gboolean                           populateFullDataPfetch(GList *seqsToFetch, BlxFetchMethod *fetchMethod, gboolean External, const BlxSeqType seqType, GError **error);
void                               blxInitConfig(const char *config_file, CommandLineOptions *options, GError **error) ;
GKeyFile*                          blxGetConfig(void) ;

void                               loadNativeFile(const char *filename, const char *buffer, GKeyFile *keyFile, BlxBlastMode *blastMode, GArray* featureLists[], GSList *supportedTypes, GSList *styles, MSP **newMsps, GList **newSeqs, GList *columnList, GHashTable *lookupTable, const int refSeqOffset, const IntRange* const refSeqRange, GError **error);
void                               blxMergeFeatures(MSP *newMsps, GList *newSeqs, MSP **mspList, GList **seqList);

/* Create/destroy sequences and MSPs */
void                               blviewResetGlobals();

BlxStyle*                          createBlxStyle(const char *styleName, const char *fillColor, const char *fillColorSelected, const char *lineColor, const char *lineColorSelected, const char *fillColorUtr, const char *fillColorUtrSelected, const char *lineColorUtr, const char *lineColorUtrSelected, GError **error);
void                               destroyBlxStyle(BlxStyle *style);

void                               createPfetchDropDownBox(GtkBox *box, GtkWidget *blxWindow);

gboolean                           bulkFetchSequences(const int attempt, 
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
                                                      GHashTable *lookupTable);


/* Dotter/Blixem Package-wide variables...........MORE GLOBALS...... */
extern char      *stdcode1[];      /* 1-letter amino acid translation code */
extern int       aa_atob[];
extern int       PAM120[23][23];
extern int       dotterGraph;
extern float     fsPlotHeight;
extern GtkWidget *blixemWindow;


/* blxFetchDb.c */
void sqliteFetchSequences(GList *seqsToFetch, const BlxFetchMethod* const fetchMethod, GList *columnList, GError **error);
void sqliteFetchSequence(const BlxSequence* const blxSeq, const BlxFetchMethod* const fetchMethod,const gboolean displayResults, const int attempt,GtkWidget *blxWindow);


#endif /*  !defined DEF_BLIXEM_P_H */
