/*  File: blxmsp.h
 *  Author: Gemma Barson, 2010-09-02
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
 *      Ed  fGriffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: Defines the MSP data struct and related functions. I think 
 *              MSP stood for Matching Segment Pair and originally represented
 *              an alignment. However, it is now used for other feature
 *              types too, so should be renamed to "BlxFeature" or something.
 *              Ideally also it would be separated out into a base feature and
 *              derived features for the different features types.
 *----------------------------------------------------------------------------
 */ 

#ifndef _blxmsp_included_
#define _blxmsp_included_

#include <seqtoolsUtils/utilities.h>
#include <stdlib.h>
#include <gtk/gtk.h>


//extern GArray    *fsArr;                 /* in dotter.c */

/* These are used by both blixem and dotter */
#define selectFeaturesStr          "Feature series selection tool"
#define XY_NOT_FILLED -1000        /* Magic value meaning "value not provided" */


/* Key names for values in data-type stanzas in the config file */
#define BLIXEM_GROUP                 "blixem"                   /* [blixem] stanza */
#define SEQTOOLS_BULK_FETCH          "bulk-fetch"               /* methods(s) for batch-fetching sequences on start */
#define SEQTOOLS_USER_FETCH          "user-fetch"               /* method(s) for interactively fetching sequences by the user */
#define SEQTOOLS_OPTIONAL_FETCH      "optional-fetch"           /* method(s) for batch-fetching additional data on user request */
#define SEQTOOLS_GFF_FILENAME_KEY    "file"
#define SEQTOOLS_WINDOW_COLOR        "session-colour"           /* color for the main window background */

/* Main Blixem error domain */
#define BLX_ERROR g_quark_from_string("Blixem")

/* Error codes */
typedef enum
{
  BLX_ERROR_SEQ_SEGMENT,              /* error finding sequence segment */
  BLX_ERROR_EMPTY_STRING,           /* error code for when user entered a zero-length string */
  BLX_ERROR_STRING_NOT_FOUND,       /* error code for when a search string is not found */
  BLX_ERROR_SEQ_NAME_NOT_FOUND,     /* the sequence name(s) being searched for were not found */
  BLX_ERROR_SEQ_DATA_MISMATCH,      /* same sequence was parsed more than once and data does not match */
  BLX_ERROR_INVALID_COLUMN          /* error when an invalid column is requested */
} BlxError;



/* Supported types of MSP */
typedef enum 
{
  BLXMSP_INVALID,                /* No valid type was set */
  
  BLXMSP_MATCH,                  /* A match (i.e. alignment) */
  BLXMSP_MATCH_SET,              /* The parent of a set of matches. Can be used to specify generic properties such as color. */
  BLXMSP_CDS,                    /* CDS (coding) region of an exon */
  BLXMSP_UTR,                    /* UTR (untranslated) region of an exon */
  BLXMSP_INTRON,                 /* Intron */
  BLXMSP_EXON,                   /* Exon (should appear AFTER CDS and UTR for sorting, as required by constructTranscriptData) */
  BLXMSP_POLYA_SITE,             /* polyA tail site */
  BLXMSP_POLYA_SIGNAL,           /* polyA signal */
  
  BLXMSP_VARIATION,              /* SNP, substitution, deletion, insertion */
  
  BLXMSP_HSP,                    /* obsolete? */
  BLXMSP_GSP,                    /* obsolete? */
  
  BLXMSP_FS_SEG,                 /* Feature Series Segment - obsolete? */
  BLXMSP_XY_PLOT,                /* x/y coordinates - for plotting feature-series curves - obsolete? */

  BLXMSP_REGION,                 /* Region */
  BLXMSP_GAP,                    /* Gap, e.g. assembly gap */

  
  
  BLXMSP_NUM_TYPES,               /* the number of valid MSP types - any types following this may be used
                                   * e.g. for parsing, but no real MSP will be created from them */
  
  BLXMSP_TRANSCRIPT              /* Transcript */
} BlxMspType;


/* Type definition for BlxSequences */
typedef enum
  {
    BLXSEQUENCE_UNSET,
    BLXSEQUENCE_TRANSCRIPT,         /* transcript (i.e. collection of exons and introns) */
    BLXSEQUENCE_MATCH,              /* match sequence (i.e. collection of matches) */
    BLXSEQUENCE_VARIATION,          /* variation (i.e. insertion, deletion or substitution) */
    BLXSEQUENCE_REGION              /* region */
  } BlxSequenceType;


/* This enum provides an ID for each type of data model that we use for the
 * detail-view trees. */
typedef enum
  {
    BLXMODEL_NORMAL,                /* the normal model, where each row contains one feature */
    BLXMODEL_SQUASHED,              /* the "squashed" model, where all MSPs from the same sequence appear on the same row */
    
    BLXMODEL_NUM_MODELS             /* the number of model IDs. MUST BE LAST IN LIST */
  } BlxModelId;

/* This enum contains a list of all the boolean flags in the BlxDataType */
/* YOU MUST UPDATE g_MspFlagConfigKeys AFTER CHANGING THIS ENUM */
typedef enum
  {
    MSPFLAG_MIN,                        /* Start index for looping through flags */
  
    MSPFLAG_LINK_FEATURES_BY_NAME,      /* whether features with the same name are part of the same parent */
    MSPFLAG_SQUASH_LINKED_FEATURES,     /* whether features with the same parent should be compressed onto the same line when you do 'squash matches' */
    MSPFLAG_SQUASH_IDENTICAL_FEATURES,  /* whether alignments that are identical should be compressed onto the same line when you do 'squash matches' */
    MSPFLAG_STRAND_SPECIFIC,            /* if false, show all features on the forward strand; else only show forward features on forward strand */
    MSPFLAG_SHOW_REVERSE_STRAND,        /* if true, and strand-specific, show rev strand features in rev strand area of display */

    /* Add new items above here. */
    /* YOU MUST UPDATE g_MspFlagConfigKeys AFTER CHANGING THIS ENUM */

    MSPFLAG_NUM_FLAGS                   /* Total number of flags e.g. for creating arrays and loops etc */
  } MspFlag;


/* Defines a data type for sequences. The data type contains properties applicable
 * to multiple sequences, e.g. which fetch method to use. */
typedef struct _BlxDataType
  {
    GQuark name;           /* the name of the data-type */
    GArray *bulkFetch;     /* list of fetch methods (by name as a GQuark) to use when bulk fetching sequences, in order of priority */
    GArray *userFetch;     /* list of fetch methods (by name as a GQuark) to use when user fetches a sequence, in order of priority */
    GArray *optionalFetch; /* list of fetch methods (by name as a GQuark) to use when user requests optional data to be loaded */
    gboolean flags[MSPFLAG_NUM_FLAGS];  /* boolean flags */
  } BlxDataType;


/* COLUMNS: To add a new column you must do the following:
 *    - create the column in blxCreateColumns or dotterCreateColumns
 * and optionally:
 *    - add a BlxColumnId enum, if you need a convenient way of referring specifically to the column in the code;
 *    - add a custom data function in createTreeColumn;
 *    - add a custom header widget and/or header refresh function in createTreeColHeader;
 *    - specify specific sort behaviour in sortColumnCompareFunc. */

/* This enum declares identifiers for known data columns. It is
 * NOT A COMPREHENSIVE LIST of columns because further columns
 * can be added dynamically to the column list in BlxCreateColumns
 * or dotterCreateColumns. This enum is just used as a convenient
 * way of referring to the main set of known columns. */
typedef enum
  {
    BLXCOL_NONE=-1,             /* Used for sorting to indicate that no sorting is required; negative value => not a valid column ID in the trees */

    BLXCOL_SEQNAME=0,           /* The match sequence's name */
    BLXCOL_SOURCE,              /* The match sequence's source */

    BLXCOL_GROUP,               /* The group that this alignment belongs to */
    BLXCOL_SCORE,               /* The alignment's score */
    BLXCOL_ID,                  /* The alignment's %ID */
    BLXCOL_START,               /* The start coord of the alignment on the match sequence */
    BLXCOL_SEQUENCE,            /* This column will display the part of the alignment currently in the display range. */
    BLXCOL_END,                 /* The end coord of the alignment on the match sequence */

    /* The following columns are optional */
    BLXCOL_ORGANISM,
    BLXCOL_GENE_NAME,
    BLXCOL_TISSUE_TYPE,
    BLXCOL_STRAIN,

    BLXCOL_NUM_COLUMNS          /* Number of main columns. Further additional columns may be
                                 * added dynamically with a columnId greater than or equal to this. */
  } BlxColumnId;


/* This struct describes a column in the detail view. Multiple widgets (i.e. headers
 * and tree columns) in the detail view must all have columns that share the same
 * properties (namely the column width). */
typedef struct _BlxColumnInfo
  {
    BlxColumnId columnId;       /* the column identifier */
    int columnIdx;              /* 0-based index of columns in display order */
    GType type;                 /* the type of data, e.g. G_TYPE_STRING */
    GtkWidget *headerWidget;    /* the header widget for this column (in the detail-view header) */
    GtkCallback refreshFunc;    /* the function that will be called on the header widget when columns are refreshed */
    const char *title;          /* the default column title */
    const char *propertyName;   /* the property name (used to set the data for the SequenceCellRenderer) */
    const char *sortName;       /* the name to display in the sort-by drop-down box (NULL if the view is not sortable on this column) */

    GQuark emblId;              /* 2-char embl line ID, e.g. 'SQ' for sequence or 'OS' for organism */
    GQuark emblTag;             /* tag name within an embl line, e.g. the 'tissue_type' tag within the 'FT' section. 
                                 * Only supports tags of the following format (i.e. like those in the 'FT' section):
                                 *     /tissue_type="testis"    */

    int width;                  /* the column width */
    gboolean dataLoaded;        /* whether the data for this column has been loaded from the EMBL file (or tried to be loaded, if it doesn't exist) */
    gboolean showColumn;        /* whether the column should be shown in the detail view */
    gboolean showSummary;       /* whether the column should be shown in the summary info (i.e. the mouse-over feedback bar) */
    gboolean canShowSummary;    /* whether it's possible to show summary info for this column */
    gboolean searchable;        /* whether searching sequences by data in this column is supported */
  } BlxColumnInfo;


/* Structure that contains information about a sequence */
typedef struct _BlxSequence
{
  BlxSequenceType type;            /* What type of collection of MSPs this is */
  BlxDataType *dataType;           /* Optional data type that specifies additional properties for this type of sequence data */
  char *idTag;                     /* Unique identifier e.g. from ID tag in GFF files */

  GArray* values;                  /* Array of values (as GValue) for the columns */
 
  BlxStrand strand;                /* which strand of the sequence this is */
  gboolean sequenceReqd;           /* whether the sequence data is required (e.g. it is not needed for exons/introns etc.) */

  IntRange qRangeFwd;              /* the extent of alignments from this sequence on the ref sequence forward strand */ 
  IntRange qRangeRev;              /* the extent of alignments from this sequence on the ref sequence reverse strand */ 
  char *organismAbbrev;            /* internally-calculated abbreviation for the Organism column value, if set */
  
  GList *mspList;                  /* list of MSPs from this sequence */
} BlxSequence;


typedef struct _FeatureSeries
{
  char        *name;             /* Name of the Feature Series */
  int         on;                /* Flag to show/hide the Feature Series */
  int         order;             /* Order number used for sorting */
  float       x;                   /* Series offset on x axis, to bump series on the screen */
  float       y;                   /* Series offset on y axis */
  int         xy;                  /* Flag for XY plot series */
} FeatureSeries;


/* Shapes of XY curves */
typedef enum
{ 
  BLXCURVE_PARTIAL, 
  BLXCURVE_INTERPOLATE, 
  BLXCURVE_BADSHAPE
} BlxCurveShape;


/* Structure holding information about a feature (see note at the top of this
 * file about the naming of this struct). */
typedef struct _MSP
{
  gchar* treePaths[BLXMODEL_NUM_MODELS]; /* identifies the row in the tree data model this msp is in (for the model given by the modelId) */

  struct _MSP       *next;
  GList             *childMsps;    /* Child MSPs of this MSP if it has them, e.g. an exon has CDS and UTR children (part_of relationship). */
  
  BlxMspType        type;          /* The type of the MSP, e.g. match, exon, SNP etc. */
  gdouble           score;         /* Score as a percentage. Technically this should be a weighted score taking into account gaps, length of the match etc., but for unknown reasons the ID has always been passed instead of score and the ID gets stored in here */
  gdouble           id;            /* Identity as a percentage. A simple comparison of bases within the match, ignoring gaps etc. Currently this is calculated internally by blixem. */
  int               phase;         /* phase: q start coord is offset by this amount to give the first base in the first complete codon (only relevant to CDSs) */
  GQuark            filename;      /* optional filename, e.g. for features used to fetch data from a bam file */
  
  char              *qname;        /* For Dotter, the MSP can belong to either sequence */
  IntRange          qRange;        /* the range of coords on the ref sequence where the alignment lies */
  BlxStrand         qStrand;       /* which strand on the reference sequence the match is on */
  int               qFrame;        /* which frame on the reference sequence the match is on */
  
  BlxSequence       *sSequence;    /* pointer to a struct holding info about the sequence/strand this match is from */
  char              *sname;        /* sequence name (could be different to the sequence name in
                                      the blxSequence e.g. exons have a postfixed 'x') */
  char              *sname_orig;   /* sequence name, original case version of sname. */
  IntRange          sRange;        /* the range of coords on the match sequence where the alignment lies */
  
  /* The following ranges are all calculated from the above but are
   * cached in the MSP because they are used a lot. Note that these
   * these are not current used by dotter so dotter does not bother to
   * initialise them. */
  IntRange          displayRange;  /* the same range as qRange but in display coords */
  IntRange          fullRange;     /* the full range of display coords to show this match against (includes any unaligned portions of sequence that we're showing) */
  IntRange          fullSRange;    /* the full range of coords on the match sequence that we're showing (including any unaligned portions of sequence) */
  
  char              *desc;         /* Optional description text for the MSP */
  GSList            *gaps;         /* Array of "gaps" in this homolgy (this is a bit of a misnomer because the array
                                    * gives the ranges of the bits that align rather than the ranges of the gaps in between */
  
  BlxStyle          *style;        /* Specifies drawing style for this MSP, e.g. fill color and line color */
  
  /* obsolete? */
  FeatureSeries     *fs;           /* Feature series that this MSP belongs to */
  int               fsColor;       /* Color to draw this MSP in the feature series */
  BlxCurveShape     fsShape;       /* Shape data for drawing feature series curves, i.e. XY type PARTIAL or INTERPOLATE shapes */
  GArray            *xy;            /* For XY plot feature series */
} MSP ;


/* MSP functions */
gboolean              typeIsExon(const BlxMspType mspType);
gboolean              typeIsIntron(const BlxMspType mspType);
gboolean              typeIsTranscript(const BlxMspType mspType);
gboolean              typeIsMatch(const BlxMspType mspType);
gboolean              typeIsVariation(const BlxMspType mspType);
gboolean              typeIsRegion(const BlxMspType mspType);
gboolean              typeShownInDetailView(const BlxMspType mspType);
gboolean              blxSequenceShownInDetailView(const BlxSequence *blxSeq);
gboolean              blxSequenceShownInGrid(const BlxSequence *blxSeq);

const char*           mspGetRefName(const MSP* const msp);
int                   mspGetRefFrame(const MSP* const msp, const BlxSeqType seqType);
BlxStrand             mspGetRefStrand(const MSP* const msp);
BlxStrand             mspGetMatchStrand(const MSP* const msp);
const char*           mspGetMatchSeq(const MSP* const msp);
const char*           mspGetSource(const MSP* const msp);
const char*           mspGetSName(const MSP *msp);
const char*           mspGetSNameOrig(const MSP *msp);
const IntRange*       mspGetRefCoords(const MSP* const msp);
const IntRange*       mspGetMatchCoords(const MSP* const msp);
int                   mspGetQStart(const MSP* const msp);
int                   mspGetQEnd(const MSP* const msp);
int                   mspGetSStart(const MSP* const msp);
int                   mspGetSEnd(const MSP* const msp);
int                   mspGetQRangeLen(const MSP* const msp);
int                   mspGetSRangeLen(const MSP* const msp);
int                   mspGetMatchSeqLen(const MSP* const msp);

const GdkColor*       mspGetColor(const MSP* const msp, 
                                  GArray *defaultColors,
                                  const int defaultColorId,
                                  const BlxSequence *blxSeq, 
                                  const gboolean selected, 
                                  const gboolean usePrintColors, 
                                  const gboolean fill,
                                  const int exonFillColorId,
                                  const int exonLineColorId,
                                  const int cdsFillColorId,
                                  const int cdsLineColorId,
                                  const int utrFillColorId,
                                  const int utrLineColorId);

const char*           mspGetColumn(const MSP* const msp, const BlxColumnId columnId);
const char*           mspGetOrganism(const MSP* const msp);
const char*           mspGetOrganismAbbrev(const MSP* const msp);
const char*           mspGetGeneName(const MSP* const msp);
const char*           mspGetTissueType(const MSP* const msp);
const char*           mspGetStrain(const MSP* const msp);
char*                 mspGetCoordsAsString(const MSP* const msp);
gchar*                mspGetTreePath(const MSP* const msp, BlxModelId modelId);

MSP*                  mspArrayIdx(const GArray* const array, const int idx);
gint                  compareFuncMspPos(gconstpointer a, gconstpointer b);
gint                  compareFuncMspArray(gconstpointer a, gconstpointer b);

gboolean              mspLayerIsVisible(const MSP* const msp);
gboolean              mspIsExon(const MSP* const msp);
gboolean              mspIsIntron(const MSP* const msp);
gboolean              mspIsSnp(const MSP* const msp);
gboolean              mspIsBlastMatch(const MSP* const msp);
gboolean              mspIsPolyASite(const MSP* const msp);
gboolean              mspIsVariation(const MSP* const msp);
gboolean              mspIsZeroLenVariation(const MSP* const msp);

gboolean              mspHasSName(const MSP* const msp);
gboolean              mspHasSSeq(const MSP * const msp);
gboolean              mspHasSCoords(const MSP* const msp);
gboolean              mspHasSStrand(const MSP* const msp);
gboolean              mspHasPolyATail(const MSP* const msp);
gboolean              mspCoordInPolyATail(const int coord, const MSP* const msp);
gboolean              mspGetFlag(const MSP* const msp, const MspFlag flag);
const char*           mspFlagGetConfigKey(const MspFlag flag);
gboolean              mspFlagGetDefault(const MspFlag flag);
void                  mspFlagSetDefault(const MspFlag flag, const gboolean value);

ColinearityType       mspIsColinear(const MSP* const msp1, const MSP* const msp2);

int                   getMaxMspLen();
void                  setMaxMspLen(const int len);

void                  writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *range1, IntRange *range2);
BlxSequence*          readBlxSequenceFromText(char *text, int *numMsps);
void                  writeMspToOutput(FILE *pipe, const MSP* const msp);
void                  readMspFromText(MSP *msp, char *text);

void                  destroyMspList(MSP **mspList);
void                  destroyBlxSequenceList(GList **seqList);
void                  destroyMspData(MSP *msp);
MSP*                  createEmptyMsp(MSP **lastMsp, MSP **mspList);
MSP*                  createNewMsp(GArray* featureLists[],
                                   MSP **lastMsp, MSP **mspList, GList **seqList, GList *columnList,
                                   const BlxMspType mspType, BlxDataType *dataType, const char *source,
                                   const gdouble score, const gdouble percentId, const int phase, const char *idTag,
                                   const char *qName, const int qStart, const int qEnd,
                                   const BlxStrand qStrand, const int qFrame,
                                   const char *sName, const char *const sName_orig, int sStart, const int sEnd, 
                                   const BlxStrand sStrand, char *sequence,
                                   const GQuark filename, GHashTable *lookupTable, GError **error);  
MSP*                  copyMsp(const MSP* const src, GArray* featureLists[], MSP **lastMsp, MSP **mspList, const gboolean addToParent);

//void                  insertFS(MSP *msp, char *series);

void                  finaliseBlxSequences(GArray* featureLists[], MSP **mspList, GList **seqList, GList *columnList, const int offset, const BlxSeqType seqType, const int numFrames, const IntRange* const refSeqRange, const gboolean calcFrame, GHashTable *lookupTable);
int                   findMspListSExtent(GList *mspList, const gboolean findMin);
int                   findMspListQExtent(GList *mspList, const gboolean findMin, const BlxStrand strand);

/* Feature series */
gint                  fsSortByNameCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);
gint                  fsSortByOrderCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);

/* Columns */
gint                  columnIdxCompareFunc(gconstpointer a, gconstpointer b);

/* BlxSequence */
char*                 blxSequenceGetSummaryInfo(const BlxSequence* const blxSeq, GList *columnList);
BlxDataType*          createBlxDataType();
void                  destroyBlxDataType(BlxDataType **blxDataType);
const char*           getDataTypeName(BlxDataType *blxDataType);
BlxSequence*          createEmptyBlxSequence();
void                  addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error);
BlxSequence*          addBlxSequence(const char *name, const char *name_orig, const char *idTag,
                                     BlxStrand strand, BlxDataType *dataType, const char *source,
                                     GList **seqList, GList *columnList, char *sequence, MSP *msp,
                                     GHashTable *lookupTable, GError **error);
GList*                blxSequenceConstructCdsList(BlxSequence *seq);
void                  blxSequenceSetValue(const BlxSequence *seq, const int columnId, GValue *value);
void                  blxSequenceSetValueFromString(const BlxSequence *seq, const int columnId, const char *inputStr);
void                  blxSequenceSetColumn(BlxSequence *seq, const char *colName, const char *value, GList *columnList);
const char*           blxSequenceGetName(const BlxSequence *seq);
GQuark                blxSequenceGetFetchMethod(const BlxSequence *seq, const gboolean bulk, const gboolean optionalColumns, const int index, const GArray *defaultMethods);
int                   blxSequenceGetLength(const BlxSequence *seq);
const char*           blxSequenceGetSequence(const BlxSequence *seq);
gboolean              blxSequenceRequiresSeqData(const BlxSequence *seq);
gboolean              blxSequenceRequiresOptionalData(const BlxSequence *seq);
gboolean              blxSequenceRequiresColumnData(const BlxSequence *seq, const BlxColumnId columnId);
BlxSequence*          blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs);
char*                 blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, GList *columnList);
int                   blxSequenceGetStart(const BlxSequence *seq, const BlxStrand strand);
int                   blxSequenceGetEnd(const BlxSequence *seq, const BlxStrand strand);
const char*           blxSequenceGetSource(const BlxSequence *seq);
gboolean              blxSequenceGetLinkFeatures(const BlxSequence *seq, const gboolean defaultLinkFeatures);

GValue*               blxSequenceGetValue(const BlxSequence *seq, const int columnId);

const char*           blxSequenceGetValueAsString(const BlxSequence *seq, const int columnId);
const char*           blxSequenceGetColumn(const BlxSequence* const blxSeq, const BlxColumnId columnId);
const char*           blxSequenceGetOrganism(const BlxSequence *seq);
const char*           blxSequenceGetOrganismAbbrev(const BlxSequence *seq);
const char*           blxSequenceGetGeneName(const BlxSequence *seq);
const char*           blxSequenceGetTissueType(const BlxSequence *seq);
const char*           blxSequenceGetStrain(const BlxSequence *seq);
char*                 blxSequenceGetFasta(const BlxSequence *seq);
gboolean              blxSequenceGetFlag(const BlxSequence* const blxSeq, const MspFlag flag);

void                  destroyBlxSequence(BlxSequence *seq);

void                  blxColumnCreate(BlxColumnId columnId, const gboolean createHeader, const char *title,
                                      GType type, const char *propertyName, const int defaultWidth,
                                      const gboolean dataLoaded, const gboolean showColumn,
                                      const gboolean showSummary, const gboolean canShowSummary, const gboolean searchable,
                                      const char *sortName, const char *emblId, const char *emblTag,
                                      GList **columnList);

/* BlxDataType */
gboolean              dataTypeGetFlag(const BlxDataType* const dataType, const MspFlag flag);

#endif /* _blxmsp_included_ */








