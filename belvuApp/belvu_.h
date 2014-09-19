/*  File: belvu_.h
 *  Author: Gemma Barson, 2011-03-06
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
 * Description: Internal header for Belvu package-wide code
 *----------------------------------------------------------------------------
 */

#ifndef DEF_BELVU_P_H
#define DEF_BELVU_P_H

#include "seqtoolsUtils/utilities.h"
#include "seqtoolsUtils/blxmsp.h"
#include "seqtoolsUtils/version.h"
#include <sys/param.h>
#include <gtk/gtk.h>
#include <config.h>


/*            blixem program version and information.                        */
#define BELVU_TITLE   "Belvu"
#define BELVU_PREFIX  "Belvu - " /* window title prefix */
#define BELVU_PREFIX_ABBREV  "Bv: " /* abbreviated window title prefix */
#define BELVU_DESC    "Multiple alignment visualisation tool."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define BELVU_VERSION_STRING	   SEQTOOLS_VERSION
#define BELVU_PACKAGE_VERSION	   UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, SEQTOOLS_VERSION)
#define BELVU_TITLE_STRING	   UT_MAKE_TITLE_STRING(BELVU_TITLE, BELVU_VERSION_STRING)
#define BELVU_VERSION_COMPILE	   BELVU_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define BELVU_COPYRIGHT_STRING	   UT_MAKE_COPYRIGHT_STRING("2011")
#define BELVU_WEBSITE_STRING	   "http://www.sanger.ac.uk/resources/software/seqtools/"
#define BELVU_LICENSE_STRING	   UT_MAKE_LICENCE_STRING(BELVU_TITLE)

#define BELVU_COMMENTS_STRING()                                \
"("BELVU_TITLE_STRING", "					\
UT_COMPILE_PHRASE " " UT_MAKE_COMPILE_DATE() ")\n\n"            \
AUTHOR_TEXT "\n"


#define DEFAULT_BELVU_WINDOW_WIDTH_FRACTION     0.95   /* default width of belvu window (as fraction of screen width) */
#define DEFAULT_BELVU_WINDOW_HEIGHT_FRACTION    0.45   /* default height of belvu window (as fraction of screen height) */
#define DIALOG_XPAD                             12      /* default x padding around dialog widgets */
#define DIALOG_YPAD                             8       /* default y padding around dialog widgets */
#define MAX_PIXMAP_WIDTH                        15000 /* max width in pixels of a pixmap */
#define MAX_PIXMAP_HEIGHT                       15000 /* max width in pixels of a pixmap */
#define FONT_SIZE_ENV_VAR                       "BELVU_FONT_SIZE"  /* optional environment variable to specify the default font size in points */
#define STATUSBAR_SIZE_ENV_VAR                  "BELVU_STATUSBAR_SIZE"  /* optional environment variable to specify the font size for the statusbar in points */



/* _MAX_PATH is 260 in WIN32 but each path component can be max. 256 in size */
#define FIL_BUFFER_SIZE 256
#define DIR_BUFFER_SIZE MAXPATHLEN

#define MAXNAMESIZE  256
#define MAXLENGTH 100000


#define UPGMAstr "UPGMA"
#define NJstr "NJ"

#define SWAPstr "Swap"
#define ROOTstr "Reroot"

#define UNCORRstr "Uncorrected differences"
#define KIMURAstr "Kimura correction"
#define JUKESCANTORstr "Jukes-Cantor correction"
#define STORMSONNstr "Storm & Sonnhammer correction"
#define SCOREDISTstr "Scoredist correction"

#define NOCOLOR -1
#define NN NOCOLOR
#define BG WHITE
#define BOXCOLOR  WHITE /* LIGHTGRAY */


/* List of color ids. Must be kept in sync with colorNames in belvu.c */
enum Colour    
{
  WHITE, BLACK, LIGHTGRAY, DARKGRAY,
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
};


enum { MUL, RAW };		/* Input formats for IN_FORMAT */
enum { dummy, GF, GC, GS, GR };	/* Markup types in Stockholm format */


/* Define the different file formats that belvu supports. This list
 * must be kept in sync with fileFormatNames in belvu.c */
typedef enum _BelvuFileFormat
{
  BELVU_FILE_MUL,		/* Stockholm or selex file format */
  BELVU_FILE_MSF,		/* MSF file format */
  BELVU_FILE_ALIGNED_FASTA,	/* Aligned FASTA file format */
  BELVU_FILE_UNALIGNED_FASTA,	/* Unaligned FASTA file format */
  
  BELVU_NUM_FILE_FORMATS	/* must be last in list */
} BelvuFileFormat;  


/* Tree picking methods */
typedef enum _BelvuPickMode
{
  NODESWAP,
  NODEROOT
} BelvuPickMode;


/* Tree building methods */
typedef enum _BelvuBuildMethod
{
  UPGMA,
  NJ
} BelvuBuildMethod;


/* Distance correction methods for tree building */
typedef enum _BelvuDistCorr
{
  UNCORR,
  KIMURA,
  STORMSONN,
  SCOREDIST,
  JUKESCANTOR
} BelvuDistCorr; 


/* The following are used to define default colors for certain types of features in Belvu.
 * One of several different actual colors from the BlxColor struct may be used depending 
 * on state, e.g. we use a different color if "print colors" (i.e. black and 
 * white mode) is on. */
typedef enum _BelvuColorId
  {
    BELCOLOR_MIN,                             /* dummy value so that we don't get a zero ID */
    
    BELCOLOR_BACKGROUND,                      /* default background color for general widgets */
    BELCOLOR_ALIGN_TEXT,                      /* text color for alignments */
    BELCOLOR_COLUMN_HIGHLIGHT,                /* highlight colour for the currently-selected column */
    
    BELCOLOR_TREE_BACKGROUND,                 /* background color for trees */
    BELCOLOR_TREE_LINE,                       /* default line color for the tree */
    BELCOLOR_TREE_TEXT,                       /* default text color for the tree */
    BELCOLOR_TREE_BOOTSTRAP,
    BELCOLOR_CONS_PLOT,                       /* line color of the plot on the conservation profile */
    BELCOLOR_CONS_PLOT_AVG,                   /* line color of the average line on the conservation plot */
    BELCOLOR_CONS_PLOT_SCALE,                 /* line color of the scale for the conservation plot */

    BELCOLOR_NUM_COLORS
  } BelvuColorId;


/* This defines the possible sort orders */
typedef enum _BelvuSortType
  {
    BELVU_UNSORTED,	       /* Not sorted */
  
    BELVU_SORT_SCORE,          /* Sort by score */
    BELVU_SORT_ALPHA,          /* Sort alphabetically */
    BELVU_SORT_ORGANISM,       /* Sort by organism */
    BELVU_SORT_TREE,           /* Sort by tree order */
    BELVU_SORT_SIM,            /* Sort by similarity to selected sequence */
    BELVU_SORT_ID,             /* Sort by identity to selected sequence */
    BELVU_SORT_CONS            /* Sort by conservation order */
  } BelvuSortType;


/* This defines the types of color scheme that we deal with */
typedef enum _BelvuSchemeType
  {
    BELVU_SCHEME_TYPE_RESIDUE,          /* Color by residue */
    BELVU_SCHEME_TYPE_CONS              /* Color by conservation */
  } BelvuSchemeType;


/* This defines the color schemes when coloring by residue */
typedef enum _BelvuColorSchemes
  {
    BELVU_SCHEME_NONE,              /* Clean slate (no coloring) */
    BELVU_SCHEME_ERIK,              /* Erik's original scheme */
    BELVU_SCHEME_GIBSON,            /* Toby's */
    BELVU_SCHEME_CGP,               /* Cys/Gly/Pro */
    BELVU_SCHEME_CGPH,              /* Cys/Gly/Pro/His */
    BELVU_SCHEME_CUSTOM,            /* Custom color scheme (this is activated after colors have been edited) */
    
    NUM_RESIDUE_SCHEMES,            /* this allows us to identify whether a scheme is a color-by-residue or -by-conservation mode scheme */
    
    BELVU_SCHEME_BLOSUM,            /* Average similarity by BLOSUM62 */
    BELVU_SCHEME_ID,                /* Percent identity */
    BELVU_SCHEME_ID_BLOSUM          /* Percent ID and BLOSUM62 */
  } BelvuColorScheme;


/* This defines the different levels of conservation that we colour by
 * in color-by-conservation mode */
typedef enum _BelvuConsLevel
  {
    CONS_LEVEL_MAX,                /* maximum conservation */
    CONS_LEVEL_MID,                /* medium conservation */
    CONS_LEVEL_LOW                 /* lowest conservation */
  } BelvuConsLevel;


/* Structs */
typedef struct alnStruct {
  char name[MAXNAMESIZE + 1];
  int  start;
  int  end;
  GString *sequenceStr;         /* The sequence data */
  int  nr;			/* Ordinal number in array - must be sorted on this */
  char fetch[MAXNAMESIZE + 11];
  float score;
  int  color;			/* Background color of name */
  int  markup;                  /* Markup line */
  gboolean  hide;               /* Hide this line */
  gboolean nocolor;		/* Exclude from coloring */
  char *organism;
  int startColIdx;              /* 0-based index indicating which column the sequence data starts in */
} ALN;


typedef struct _TreeNode 
{
  double dist;			/* Absolute distance position */
  double branchlen;		/* Length of branch to higher node */
  double boot;                  /* Bootstrap value */
  struct _TreeNode *left;
  struct _TreeNode *right;
  struct _TreeNode *parent;
  char *name;
  char *organism;
  ALN *aln;
  int box;
  int color;
} TreeNode;


typedef struct _Tree
{
  TreeNode *head;     /* Root node of the tree */
} Tree;


/* Struct to store bootstrap group */
typedef struct BootstrapGroupStruct
{   
  TreeNode *node;    /* Points to node in original tree (for incrementing) */
  char *s;           /* Sorted concatenation of all sequences in node, to be inserted in list */
} BootstrapGroup;


typedef struct SegStruct
{
  int  start;
  int  end;
  int  qstart;
  int  qend;
  struct SegStruct *next;
} SEG;


/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the BelvuContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    BELDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */
    
    BELDIALOG_MAKE_TREE,            /* The make-tree dialog */
    BELDIALOG_EDIT_RESIDUE_COLORS,  /* The edit-residue-colors dialog */
    BELDIALOG_EDIT_CONS_COLORS,     /* The edit-conservation-colors dialog */
    BELDIALOG_FIND,                 /* The find dialog */

    BELDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } BelvuDialogId;


typedef struct BelvuContextStruct
{
  GtkWidget *belvuWindow;          /* Pointer to the main belvu window, or NULL if it has not been created yet */
  GSList *spawnedWindows;          /* List of all top-level windows spawned from the main window */
  GtkWidget *belvuTree;            /* The tree window */
  GtkWidget *belvuAlignment;       /* The widget that draws the alignments for the main window */
  GtkWidget *consPlot;             /* The conservation-plot window */
  GtkWidget *orgsWindow;           /* The organisms window */

  GdkCursor *defaultCursor;        /* default cursor */
  GdkCursor *busyCursor;           /* cursor to use when busy */
  GdkCursor *removeSeqsCursor;     /* cursor to use when removing sequences */

  GArray *defaultColors;           /* Default colors used by Belvu */
  
  GArray *alignArr;
  GArray *organismArr;
  GArray *markupAlignArr;
  GArray *bootstrapGroups;

  ALN *selectedAln;                /* The currently-selected alignment */
  GSList *highlightedAlns;         /* List of all currently-highlighted alignments
                                    * (generally, all ALNs with the same name as
                                    * the selectedAln are highlighted) */

  Tree *mainTree;                  /* global current tree */
  TreeNode *treeBestBalancedNode;

  FILE *treeReadDistancesPipe;

  int IN_FORMAT;
  int maxScoreLen;
  int alignYStart;
  int treebootstraps;              /* Number of bootstrap trees to be made */
  int maxLen;                      /* number of columns in alignment */
  int maxTreeWidth;
  int maxNameLen;                  /* Max string length of any sequence name */
  int maxStartLen;                 /* Max string length of any sequence start */
  int maxEndLen;                   /* Max string length of any sequence end */
  int maxScoreen;                  /* Max string length of any sequence score */
  int selectedCol;                 /* The currently-selected column index (0 for unset) */
  int highlightedCol;              /* The currently-highlighted column index (0 for unset) */

  int maxfgColor;                  /* Foreground color for max conservation */
  int midfgColor;                  /* Foreground color for mid conservation */
  int lowfgColor;                  /* Foreground color for lowmax conservation */
  int maxbgColor;                  /* Background color for max conservation */
  int midbgColor;                  /* Background color for mid conservation */
  int lowbgColor;                  /* Background color for low conservation */
  int maxfgPrintColor;             /* As above but for when the 'print colors' option is enabled */
  int midfgPrintColor;
  int lowfgPrintColor;
  int maxbgPrintColor;
  int midbgPrintColor;
  int lowbgPrintColor;
  
  BelvuSchemeType schemeType;	      /* Current colour scheme mode (color-by-residue or -by-conservation) */
  BelvuColorScheme residueScheme;     /* Which color-by-residue color scheme is selected */
  BelvuColorScheme consScheme;	      /* Which color-by-conservation color scheme is selected */
  
  BelvuBuildMethod treeMethod;        /* Default building method for trees */
  BelvuDistCorr treeDistCorr;         /* Default distance correction method for trees */
  BelvuPickMode treePickMode;         /* Default action when picking a node in a tree */
  
  BelvuSortType sortType;             /* What data to sort the alignments by */
  BelvuFileFormat saveFormat;	      /* Which file format to use for saving alignments */
  
  double treeBestBalance;
  double treeBestBalance_subtrees;
  double tree_y;
  double lowIdCutoff;                 /* %id cutoff for lowest colour */
  double midIdCutoff;                 /* %id cutoff for medium colour */
  double maxIdCutoff;                 /* %id cutoff for maximum colour */
  double lowSimCutoff;                /* %id cutoff for lowest colour */
  double midSimCutoff;                /* %id cutoff for medium colour */
  double maxSimCutoff;                /* %id cutoff for maximum colour */
  double colorByResIdCutoff;          /* Cutoff when only coloring residues above a given %ID */
  double mksubfamilies_cutoff; 
  double treeScale;                   /* Default scale to use for drawing the tree */
  double treeLineWidth;               /* Default line width of the branch lines in trees */
  
  char gapChar;
  char saveSeparator;
  char treeDistString[50];
  char treeMethodString[50];
  char Title[256];
  char *fileName;		   /* Default file name for file browser */
  char *dirName;		   /* Default directory for file browser */
  char organismLabel[3];
  
  int **conservCount;              /* Matrix of conservation values (1st index is amino acid code; 2nd index is column index; value is the number of that residue in that column) - 21 x maxLen */
  int **colorMap;                  /* Matrix of conservation colours - 21 x maxLen */
  int *conservResidues;            /* Array of number of residues present in each column */
  double *conservation;            /* The max conservation in each column [0..maxLen] */

  GSList *annotationList;	   /* List of annotation lines from the input file */
  
  gboolean treeCoordsOn;
  gboolean treeReadDistancesOn;
  gboolean treePrintDistances;
  gboolean penalize_gaps;
  gboolean stripCoordTokensOn;
  gboolean saveCoordsOn;
  gboolean displayScores;
  gboolean outputBootstrapTrees;   /* Output the individual bootstrap trees */
  gboolean treebootstrapsDisplay;  /* Display bootstrap trees on screen */
  gboolean treeColorsOn;           
  gboolean treeShowOrganism;       /* whether to display the organism name in the tree */
  gboolean treeShowBranchlen;      /* whether to display the branch length in the tree */
  gboolean matchFooter;
  gboolean saved;
//  gboolean color_by_similarity;    /* FALSE => by id */
//  gboolean color_by_conserv;       /* TRUE => by id or similarity; FALSE => by residue  */
  gboolean ignoreGapsOn;
  gboolean colorByResIdOn;         /* colour by residue type above identity cutoff */
  gboolean id_blosum;              /* Paint similar residues too */
  gboolean rmEmptyColumnsOn;       /* if true then remove empty columns after deleting sequences */
  gboolean lowercaseOn;		   /* Set to true to highlight lowercase characters */
  gboolean removingSeqs;	   /* Set to true if in the 'removing sequences' mode */
  gboolean displayColors;	   /* Whether to display colors (faster without) */
  gboolean haveCustomColors;       /* Whether the custom colors have been set */
  gboolean printColorsOn;          /* Whether to use greyscale colors for printing */
  gboolean highlightOrthologs;     /* Whether to highlight orthologs or not in the tree */
  gboolean useWWWFetch;            /* Whether to fetch sequences via a web URL rather than a local program */
  gboolean initTree;               /* Start up showing the tree */
  gboolean onlyTree;               /* Start up showing only the tree */
  gboolean abbrevTitle;            /* Abbreviate window title prefixes */

  GtkWidget *dialogList[BELDIALOG_NUM_DIALOGS];   /* Array of all the persistent dialogs in the application */
  
} BelvuContext;


/* Functions */
const char*                               belvuGetAppName(void);
const char*                               belvuGetTitlePrefix(BelvuContext *bc);
const char*                               belvuGetCopyrightString(void);
const char*                               belvuGetWebSiteString(void);
const char*                               belvuGetCommentsString(void);
const char*                               belvuGetLicenseString(void);
const char*                               belvuGetVersionString(void);       

BelvuContext*                             createBelvuContext();
void                                      destroyBelvuContext(BelvuContext **bc);

void                                      greyOutInvalidActions(BelvuContext *bc);
void                                      greyOutInvalidActionsForGroup(BelvuContext *bc, GtkActionGroup *action_group);

void					  setTreeScaleCorr(BelvuContext *bc, const int treeMethod);
void					  setTreeScale(BelvuContext *bc, const double newScale) ;

gint                                      alphaorder(gconstpointer xIn, gconstpointer yIn);
gint                                      organism_order(gconstpointer xIn, gconstpointer yIn);
gint                                      organismorder(gconstpointer xIn, gconstpointer yIn);
gint                                      scoreorder(gconstpointer xIn, gconstpointer yIn);
gint                                      nrorder(gconstpointer xIn, gconstpointer yIn);

void                                      highlightScoreSort(char mode, BelvuContext *bc);
void                                      treeSortBatch(BelvuContext *bc);
void					  doSort(BelvuContext *bc, const BelvuSortType sortType, const gboolean showTree);


void                                      arrayOrder(GArray *alignArr);
gboolean                                  alignFind(GArray *alignArr, ALN *obj, int *idx);
void                                      initAln(ALN *alnp);
ALN*                                      createEmptyAln();

void                                      setOrganismColors(GArray *organismArr);

void                                      parseMulLine(BelvuContext *bc, char *line, ALN *aln);

void                                      readMatch(BelvuContext *bc, FILE *fil);                
void                                      checkAlignment(BelvuContext *bc);
void                                      setConsSchemeColors(BelvuContext *bc);
void					  updateSchemeColors(BelvuContext *bc);
void                                      saveCustomColors(BelvuContext *bc);
void                                      initResidueColors(BelvuContext *bc);
void                                      initMarkupColors(void);              
void                                      initCustomColors(void);              
gboolean				  colorByConservation(BelvuContext *bc);
gboolean				  colorByResidue(BelvuContext *bc);
gboolean				  colorBySimilarity(BelvuContext *bc);
gboolean                                  colorByResId(BelvuContext *bc);
void                                      setResidueSchemeColors(BelvuContext *bc);
const char*                               getColorNumName(const int colorNum);
const char*				  getFileFormatString(const int formatId);
int*                                      getConsColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground);
void                                      setExcludeFromConsCalc(BelvuContext *bc, const gboolean exclude);

void                                      readLabels(BelvuContext *bc, FILE *fil);

void                                      mkNonRedundant(BelvuContext *bc, double cutoff);
void                                      rmPartialSeqs(BelvuContext *bc);         
void                                      rmEmptyColumns(BelvuContext *bc, double cutoff);
void                                      rmGappySeqs(BelvuContext *bc, double cutoff);
void                                      rmFinaliseGapRemoval(BelvuContext *bc);
void					  rmOutliers(BelvuContext *bc, const double cutoff);
void					  rmScore(BelvuContext *bc, const double cutoff);
void                                      rmColumn(BelvuContext *bc, const int from, const int to);
void                                      rmColumnCutoff(BelvuContext *bc, const double from, const double to);
void                                      rmFinaliseColumnRemoval(BelvuContext *bc);

void                                      readFile(BelvuContext *bc, FILE *pipe);
void                                      writeMul(BelvuContext *bc, FILE *fil);
void                                      writeFasta(BelvuContext *bc, FILE *pipe);
void                                      writeMSF(BelvuContext *bc, FILE *pipe);

void                                      separateMarkupLines(BelvuContext *bc);
void                                      reInsertMarkupLines(BelvuContext *bc);
Tree*                                     treeMake(BelvuContext *bc, const gboolean doBootstrap, const gboolean displayFeedback);

void                                      outputProbs(BelvuContext *bc, FILE *fil);
void                                      mksubfamilies(BelvuContext *bc, double cutoff);        

void                                      treeDisplay(BelvuContext *bc);       

void                                      colorSim(BelvuContext *bc);          

char                                      b2aIndex(const int idx);
int                                       getMarkupColor(const char inputChar);
int                                       getConservColor(BelvuContext *bc, const char inputChar, const int idx);
int                                       getColor(const char inputChar);
void                                      setColor(const char inputChar, const int colorNum);
int*                                      getColorArray();
int*                                      getMarkupColorArray();
void                                      saveResidueColorScheme(BelvuContext *bc, FILE *fil);
void                                      readResidueColorScheme(BelvuContext *bc, FILE *fil, int *colorarr, const gboolean storeCustomColors);

char*                                     alnGetSeq(ALN *aln);
int                                       alnGetSeqLen(ALN *aln);
gboolean                                  isGap(char c);
int                                       strcmp_(gconstpointer xIn, gconstpointer yIn);
gboolean                                  alnArrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer));
gboolean                                  bsArrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer));
GArray*                                   copyAlignArray(GArray *inputArr);
void                                      columnCopy(GArray *alignArrDest, int destIdx, GArray *alignArrSrc, int srcIdx);
double                                    identity(char *s1, char *s2, const gboolean penalize_gaps);

void                                      convertColorNumToGdkColor(const int colorNum, const gboolean isSelected, GdkColor *result);
void                                      drawText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const char *text, int *textWidth, int *textHeight);
void                                      drawIntAsText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const int value);
void                                      drawDoubleAsText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const double value);

void                                      treeTraverse(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *treeNode));
void                                      treeSort(BelvuContext *bc, const gboolean showTree);
void                                      saveTreeNH(Tree *tree, TreeNode *node, FILE *file);

void                                      listIdentity(BelvuContext *bc);
void                                      belvuContextSetTree(BelvuContext *bc, Tree **tree);

void                                      fetchAln(BelvuContext *bc, ALN *alnp);
gboolean                                  alignmentHighlighted(BelvuContext *bc, ALN *alnp);
void					  str2aln(BelvuContext *bc, char *src, ALN *alnp) ;
void					  alncpy(ALN *dest, ALN *src);

void                                      setBusyCursor(BelvuContext *bc, const gboolean busy);

#endif /* DEF_BELVU_P_H */
