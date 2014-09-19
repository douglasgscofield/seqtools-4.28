/*  File: utilities.h
 *  Author: Gemma Barson, 2010-01-05
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
 * Description: Utility functions for the SeqTools package.
 *----------------------------------------------------------------------------
 */
 
#ifndef _utilities_h_included_
#define _utilities_h_included_

#include <gtk/gtk.h>

#define UNSET_INT                     -1   /* this value indicates an unset integer */
#define DEFAULT_LABEL_X_PAD           0    /* default x padding to use for header labels */
#define DEFAULT_PRINT_RESOLUTION      300  /* default resolution for printing in DPI */
#define TABLE_XPAD                    12      /* default x padding around table elements */
#define TABLE_YPAD                    2       /* default y padding around table elements */
#define MIN_FONT_SIZE                 2
#define MAX_FONT_SIZE                 20
#define MACHINE_RES                   0.000000001
#define SEQTOOLS_GFF_VERSION          3    /* the default GFF version */


/* Really the buffers that use this should be dynamic but I'm not going to do that, this
 * code is so poor that it doesn't warrant the effort.... */
#define NAMESIZE      12
#define MAXLINE       10000


#define max(a,b)        (((a) > (b)) ? (a) : (b))
#define min(a,b)        (((a) < (b)) ? (a) : (b))


/* Debug logging macros. #define DEBUG to enable debug output. */
#ifdef DEBUG
#define DEBUG_OUT(format, args...) debugLogLevel(0); printf(format, ##args); fflush(stdout);
#else
#define DEBUG_OUT(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_ENTER(format, args...) debugLogLevel(1); printf("-->"); printf(format, ##args); printf("\n"); fflush(stdout);
#else
#define DEBUG_ENTER(format, args...)
#endif

#ifdef DEBUG
#define DEBUG_EXIT(format, args...) debugLogLevel(-1); printf("<--"); printf(format, ##args); printf("\n"); fflush(stdout);
#else
#define DEBUG_EXIT(format, args...)
#endif


/* Generic SeqTools error domain */
#define SEQTOOLS_ERROR g_quark_from_string("SeqTools")

/* Error codes */
typedef enum
{
  SEQTOOLS_ERROR_PARSING_COLOR,       /* error parsing color string */
  SEQTOOLS_ERROR_SEQ_SEGMENT,         /* error getting the requested segment of a sequence */
  SEQTOOLS_ERROR_NO_STYLE,            /* style does not exist */
  SEQTOOLS_ERROR_EXECUTING_CMD,       /* error executing command */
  SEQTOOLS_ERROR_SEQ_TYPE,            /* error determining sequence type */
  SEQTOOLS_ERROR_NO_EXE       /* executable does not exist */
} SeqToolsError;


/* This enum defines the columns for a combo box with the common combination
 * of two columns: an int relating to an enum, and a text description. */
typedef enum _BlxComboColumns
  {
    COMBO_ENUM_COL,
    COMBO_TEXT_COL,
    N_COMBO_COLUMNS
  } BlxComboColumns;


/* Printing scale type - we can scale to fit the page width, or height, or both */
typedef enum
{
  PRINT_FIT_WIDTH,
  PRINT_FIT_HEIGHT,
  PRINT_FIT_BOTH
} PrintScaleType;



#if defined(MACINTOSH)
#define      SYSERR_FORMAT             "system error %d"
#else
#define      SYSERR_FORMAT             "system error %d - %s"
#endif


/* Special characters for displaying in sequences */
#define SEQUENCE_CHAR_DELETION '.'   /* represents a deletion */
#define SEQUENCE_CHAR_GAP      '-'   /* represents a gap. Note that gaps in the reference sequence will be highlighted as assembly gaps */
#define SEQUENCE_CHAR_PAD      '-'   /* used for padding when the sequence is unavailable */
#define SEQUENCE_CHAR_BLANK    '-'   /* used to display a blank when we're not interested in what the actual base is */
#define SEQUENCE_CHAR_INVALID  ' '   /* used to replace non-utf8 characters (which GTK can't display) */
#define SEQUENCE_CHAR_STOP     '*'   /* STOP codon */
#define SEQUENCE_CHAR_MET      'M'   /* MET codon */


/* Color strings that can be passed to create a GdkColor */
#define BLX_BLACK             "#000000"
#define BLX_WHITE             "#ffffff"

#define BLX_GREY              "#bebebe"
#define BLX_LIGHT_GREY        "#cacaca"
#define BLX_VERY_LIGHT_GREY   "#eeeeee"
#define BLX_DARK_GREY         "#929292"
#define BLX_VERY_DARK_GREY    "#5f5f5f"

#define BLX_YELLOW            "#ffff00" 
#define BLX_DARK_YELLOW       "#d0d000"
#define BLX_PALE_YELLOW       "#ffffcc"

#define BLX_BLUE              "#0000ff"
#define BLX_DARK_BLUE         "#000080"
#define BLX_ROYAL_BLUE        "#4169e1"
#define BLX_SKY_BLUE          "#87cefa"
#define BLX_PALE_BLUE         "#c0d8f0"
#define BLX_CYAN              "#15ced2"
#define BLX_LIGHT_CYAN        "#6defe9"
#define BLX_VIOLET            "#78b4f0"
#define BLX_DARK_VIOLET       "#5c98d5"
#define BLX_MID_BLUE          "#6495ED"
#define BLX_SLATE_BLUE        "#7C7FFF"

#define BLX_RED               "#ff0000"
#define BLX_LIGHT_RED         "#ff7373"
#define BLX_SALMON_RED        "#ffa07a"
#define BLX_DARK_RED          "#800000"
#define BLX_VERY_DARK_RED     "#400000"
#define BLX_ORANGE_RED        "#ff4500"
#define BLX_ORANGE            "#ffa500"
#define BLX_PALE_ORANGE       "#ff9c00"
#define BLX_MAGENTA           "#ff00ff"
#define BLX_PALE_MAGENTA      "#ffccff"
#define BLX_PURPLE            "#A020F0"
#define BLX_PALE_VIOLET       "#EE82EE"
#define BLX_CERISE            "#DE3163"

#define BLX_GREEN             "#00ff00"
#define BLX_LIGHT_GREEN       "#C1FFC1" 
#define BLX_LAWN_GREEN        "#7cfc00"
#define BLX_DARK_GREEN        "#00bb00"
#define BLX_VERY_DARK_GREEN   "#015800"

#define BLX_BROWN             "#A52A2A"

typedef struct _BlxColor
  {
    char *name;                   /* meaningful name for the color e.g. "Match" */
    char *desc;                   /* meaningful description for what the color is used for e.g. "background color for exact matches" */
    gboolean transparent;         /* if this is true, the colors below are not specified and the background color should be used instead */
    
    GdkColor normal;              /* the color in normal operation */
    GdkColor selected;            /* the color in a selected state */
    GdkColor print;               /* the color used for printing */
    GdkColor printSelected;       /* the selected-state color used for printing */
  } BlxColor;
  
  
/* This handle holds a list of pointers to all memory allocated via this handle. 
 * Use handleDestroy to free the handle and all its allocated memory. */
typedef struct _BlxHandle
  {
    GSList *memoryList;
  } BlxHandleStruct, *BlxHandle;
  
  
/* This struct is used to pass user data to the message handlers */
typedef struct _BlxMessageData
{
  char *titlePrefix;
  GtkWindow *parent;
  GtkStatusbar *statusBar;
} BlxMessageData;
  
  
/* Define a drawing style for an MSP */
typedef struct _BlxStyle
  {
    char *styleName;
    BlxColor fillColor;
    BlxColor lineColor;
    BlxColor fillColorUtr; /* used only for transcript features; fillColor and lineColor are for CDS portions and these are for UTR portions */
    BlxColor lineColorUtr; /* used only for transcript features; fillColor and lineColor are for CDS portions and these are for UTR portions */
  } BlxStyle;

/* Utility structs to hold a range of integers/doubles */
typedef struct _IntRange
  {
    int min;
    int max;
  } IntRange ;

typedef struct _DoubleRange
  {
    gdouble min;
    gdouble max;
  } DoubleRange ;
  

/* Fundamental strand direction. */
typedef enum
  {
    BLXSTRAND_NONE, 
    BLXSTRAND_FORWARD, 
    BLXSTRAND_REVERSE
  } BlxStrand ;
  
/* Fundamental type of sequence (DNA really means nucleotide, because it could be RNA as well). */
typedef enum
  {
    BLXSEQ_NONE, 
    BLXSEQ_DNA, 
    BLXSEQ_PEPTIDE
  } BlxSeqType ;

/* Fundamental blast mode used */
typedef enum
  {
    BLXMODE_UNSET, 
    BLXMODE_BLASTX, 
    BLXMODE_TBLASTX, 
    BLXMODE_BLASTN, 
    BLXMODE_TBLASTN, 
    BLXMODE_BLASTP
  } BlxBlastMode ;
  

/* Function pointer for callback functions used by widgets on dialog boxes. */
typedef gboolean (*BlxResponseCallback)(GtkWidget *widget, const gint responseId, gpointer data);

/* This struct holds generic callback information. It can be stored as a widget
 * property and called on the widget on request (e.g. by a dialog when it applies changes). */
typedef struct _CallbackData
  {
    BlxResponseCallback func;     /* Callback function to be called */
    gpointer data;                /* User data to pass to the callback function */
  } CallbackData;

/* Custom dialog response types, which can be used in addition
 * to the default types specified by GtkResponseType */
typedef enum
  {
    BLX_RESPONSE_FORWARD, 
    BLX_RESPONSE_BACK,
    BLX_RESPONSE_RESET,
    BLX_RESPONSE_CLEAR
  } BlxResponseType;


/* Struct to hold a pair of coordinate ranges, one holding ref seq coords, the other match coords */
typedef struct _CoordRange
  {
    int qStart;     /* start coord on the reference (Query) seq */
    int qEnd;       /* end coord on the reference (Query) seq */
    int sStart;     /* start coord on the match (Subject) seq */
    int sEnd;       /* end coord on the match (Subject) seq */
  } CoordRange;
  

/* Used to specify the degree of colinearity between two alignment blocks. */
typedef enum
  {
    COLINEAR_INVALID,
    COLINEAR_NOT,                                          /* blocks not colinear. */
    COLINEAR_IMPERFECT,                                    /* blocks colinear but not contiguous. */
    COLINEAR_PERFECT,                                      /* blocks colinear and contiguous. */
    COLINEARITY_N_TYPE
  } ColinearityType ;


/* Supported drag-and-drop target types */
typedef enum {TARGET_STRING, TARGET_URL} DragDropTargetType;


GdkDrawable*          widgetGetDrawable(GtkWidget *widget);
void                  widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable);
gboolean              widgetGetHidden(GtkWidget *widget);
void                  widgetSetHidden(GtkWidget *widget, const gboolean hidden);
void                  hideUserHiddenWidget(GtkWidget *widget, gpointer data);
void                  widgetClearCachedDrawable(GtkWidget *widget, gpointer data);
void                  callFuncOnAllChildWidgets(GtkWidget *widget, gpointer data);

void                  blxSetWidgetColor(GtkWidget* widget, char *colorName);
gboolean              onExposePrintable(GtkWidget *widget, GdkEventExpose *event, gpointer data);
GtkWidget*            createLabel(const char *text, const gdouble xalign, const gdouble yalign, const gboolean ellipsize, const gboolean enableCopyPaste, const gboolean showWhenPrinting);
GtkWidget*            getLabelWidget(GtkWidget *widget);
void                  labelSetFont(GtkWidget *widget, PangoFontDescription *fontDesc);
GdkDrawable*          createBlankPixmap(GtkWidget *widget);
GdkDrawable*          createBlankSizedPixmap(GtkWidget *widget, GdkDrawable *window, const int width, const int height);

gboolean              clickedInRect(GdkEventButton *event, GdkRectangle *rect, const int minWidth);

BlxSeqType            determineSeqType(char *seq, GError **error);
void                  argvAdd(int *argc, char ***argv, const char *s);
char*                 getSystemErrorText();
gpointer              handleAlloc(BlxHandle *handle, size_t numBytes);
BlxHandle             handleCreate();
void                  handleDestroy(BlxHandle *handle);
BlxStyle*             getBlxStyle(const char *styleName, GSList *styles, GError **error);

void                  sortValues(int *val1, int *val2, gboolean forwards);
int                   numDigitsInInt(int val);
gboolean              getColorFromString(const char *colorStr, GdkColor *color, GError **error);
void                  getSelectionColor(const GdkColor* const origColor, GdkColor *result);
void                  getDropShadowColor(const GdkColor* const origColor, GdkColor *result);
void                  convertToGrayscale(const GdkColor* const origColor, GdkColor *result);
void                  adjustColorBrightness(const GdkColor* const origColor, const double factor, GdkColor *result);

void                  getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax);

gboolean              isDelimiter(const char c);
char*                 removeDelimiters(char *text);
GArray*               keyFileGetCsv(GKeyFile *keyFile, const char *group, const char *key, GError **error);

int                   getRangeLength(const IntRange* const range);
int                   getRangeCentre(const IntRange* const range);
void                  centreRangeOnCoord(IntRange *range, const int coord, const int length);
gboolean              valueWithinRange(const int value, const IntRange* const range);
gboolean              rangesOverlap(const IntRange* const range1, const IntRange* const range2);
gboolean              rangesAdjacent(const IntRange* const range1, const IntRange* const range2);
gboolean              rangesEqual(const IntRange* const range1, const IntRange* const range2);
void                  boundsLimitValue(int *value, const IntRange* const range);
void                  boundsLimitRange(IntRange *range, const IntRange* const limit, const gboolean maintainLen);
gboolean              pointInRect(const int x, const int y, const GdkRectangle* const rect);
char                  convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType);

void                  convertDisplayRangeToDnaRange(const IntRange* const displayRange, 
                                                    const BlxSeqType displaySeqType,
                                                    const int numFrames,
                                                    const gboolean displayRev,
                                                    const IntRange* const refSeqRange,
                                                    IntRange *result);

int                   convertDisplayIdxToDnaIdx(const int inputIdx, 
                                                const BlxSeqType inputIdxType,
                                                const int frame, 
                                                const int baseNum, 
                                                const int numFrames,
                                                const gboolean displayRev,
                                                const IntRange* const dnaIdxRange);

int                   convertDnaIdxToDisplayIdx(const int dnaIdx, 
                                                const BlxSeqType displaySeqType,
                                                const int frame,
                                                const int numFrames, 
                                                const gboolean displayRev,
                                                const IntRange* const dnaIdxRange,
                                                int *baseNum);

char                  getStrandAsChar(const BlxStrand strand);

int                   roundNearest(const double val);
int                   roundToValue(const int inputVal, const int roundTo);
int                   roundToValueFromList(const int inputVal, GSList *roundValues, int *roundedTo);
int                   roundUpToValueFromList(const int inputVal, GSList *roundValues, int *roundedTo);

char                  getSequenceIndex(char *seq, 
                                       const int qIdx, 
                                       const gboolean complement, 
                                       const IntRange* const seqRange,
                                       const BlxSeqType seqType);

int                   getStartDnaCoord(const IntRange* const displayRange, 
                                       const int frame,
                                       const BlxSeqType displaySeqType, 
                                       const gboolean displayRev, 
                                       const int numFrames,
                                       const IntRange* const refSeqRange);

int                   getEndDnaCoord(const IntRange* const displayRange, 
                                     const int frame,
                                     const BlxSeqType displaySeqType, 
                                     const gboolean displayRev, 
                                     const int numFrames,
                                     const IntRange* const refSeqRange);

int                   wildcardSearch(const char *textToSearch, const char *searchStr);

char*                 convertIntToString(const int value);
char*                 convertDoubleToString(const gdouble value, const int numDp);
int                   convertStringToInt(const char *inputStr);
gboolean              isWhitespaceChar(const char curChar);
gboolean              isNewlineChar(const char curChar);
char*                 abbreviateText(const char *inputStr, const int maxLen);
gboolean              stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive);
gboolean              isValidIupacChar(const char inputChar, const BlxSeqType seqType);
void                  stringProtect(FILE *file, const char *string);
char*                 stringUnprotect(char **textp, char *target);

int                   invertCoord(const int coord, const IntRange* const range, const gboolean invert);

void                  popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void                  defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data);
void                  g_message_info(const char *formatStr, ...);

GtkWidget*            showMessageDialog(const char *title,  
                                        const char *messageText,
                                        GtkWidget *parent,
                                        const int initWidth,
                                        const int maxHeight,
                                        const gboolean wrapText,
                                        const gboolean useMarkup,
                                        PangoFontDescription *fontDesc,
                                        GtkTextView **textView);

void                    destroyMessageList();

GtkWidget*              createScrollableTextView(const char *messageText,
                                                 const gboolean wrapText,
                                                 PangoFontDescription *fontDesc,
                                                 const gboolean useMarkup,
                                                 int *width,
                                                 int *height,
                                                 GtkTextView **textViewOut);
                                    
void                  widgetSetCallbackData(GtkWidget *widget, BlxResponseCallback callbackFunc, gpointer callbackData);
gboolean              widgetCallAllCallbacks(GtkWidget *widget, gpointer data);
void                  onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data);
void                  onCloseDialog(GtkDialog *dialog, gpointer data);
void                  dialogClearContentArea(GtkDialog *dialog);
GtkWidget*            getPersistentDialog(GtkWidget* dialogList[], const int dialogId);
void                  addPersistentDialog(GtkWidget* dialogList[], const int dialogId, GtkWidget *widget);

void                  setPrimaryClipboardText(const char *text);
void                  setDefaultClipboardText(const char *text);
void                  requestPrimaryClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);
void                  requestDefaultClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data);

int                   parseMatchLine(const char *inputText,
                                     char **matchNameOut,
                                     int *matchStartOut, 
                                     int *matchEndOut, 
                                     int *matchLenOut);

GList*                parseMatchList(const char *inputText);

GtkWidget*            getNamedChildWidget(GtkWidget *widget, const gchar *searchName);

gint                  runConfirmationBox(GtkWidget *blxWindow, const char *title, const char *messageText);

void                  prefixError(GError *error, const char *prefixStr, ...);
void                  postfixError(GError *error, const char *formatStr, ...);
void                  reportAndClearIfError(GError **error, GLogLevelFlags log_level);

void                  intrangeSetValues(IntRange *range, const int val1, const int val2);

#ifdef DEBUG
void                  debugLogLevel(const int increaseAmt);
#endif

void                  drawRect(GdkDrawable *drawable, GdkColor *color, const int x, const int y, const int width, const int height, const double alpha, cairo_operator_t op);
void                  drawHighlightBox(GdkDrawable *drawable, const GdkRectangle* const rect, const gint minWidth, GdkColor *color);

void                  blxPrintWidget(GtkWidget *widget, GdkDrawable *drawable, GtkWindow *window, GtkPrintSettings **printSettings, GtkPageSetup **pageSetup, const char *filename, const gboolean printCachedOnly, const PrintScaleType scaleType);
void                  setStatusBarShadowStyle(GtkWidget *statusBar, const char *shadowStyle);

BlxColor*             getBlxColor(GArray *defaultColors, const int colorId);
GdkColor*             getGdkColor(int colorId, GArray *defaultColors, const gboolean selected, const gboolean usePrintColors);
const GdkColor*       blxColorGetColor(const BlxColor *blxColor, const gboolean selected, const gboolean usePrintColors);
char*                 convertColorToString(GdkColor *color);
void                  destroyBlxColor(BlxColor *blxColor);
gboolean              colorsEqual(GdkColor *color1, GdkColor *color2);

void                  createBlxColor(GArray *defaultColors,
                                     int colorId,
                                     const char *name, 
                                     const char *desc, 
                                     const char *normalCol, 
                                     const char *printCol,
                                     const char *normalColSelected,
                                     const char *printColSelected);

gdouble               pixelsPerBase(const gint displayWidth, 
                                    const IntRange* const displayRange);

gdouble               convertBaseIdxToRectPos(const gint dnaIdx, 
                                              const GdkRectangle* const rect, 
                                              const IntRange* const dnaDispRange,
                                              const gboolean horizontal,
                                              const gboolean displayRev,
                                              const gboolean clip);

gchar*                    getSequenceSegment(const char* const dnaSequence,
                                             IntRange *qRange,
                                             const BlxStrand strand,
                                             const BlxSeqType srcSeqType,
                                             const BlxSeqType destSeqType,
                                             const int frame,
                                             const int numFrames,
                                             const IntRange* const refSeqRange,
                                             const BlxBlastMode blastMode,
                                             char **geneticCode,
                                             const gboolean displayRev,
                                             const gboolean reverseResult,
                                             const gboolean allowComplement,
                                             GError **error);

const char*                        findFixedWidthFont(GtkWidget *widget);
const char*                        findFixedWidthFontFamily(GtkWidget *widget, GList *pref_families);
void                               getFontCharSize(GtkWidget *widget, PangoFontDescription *fontDesc, gdouble *width, gdouble *height);

GtkWidget*                         createToolbarHandle();
GtkToolItem*                       addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget, const int position);

void                               onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data);
void                               collatePixmaps(GtkWidget *widget, gpointer data);
void                               onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data);
void                               setWidgetBackgroundColor(GtkWidget *widget, gpointer data);

void                               forceResize(GtkWidget *widget);

gboolean                           onComboChanged(GtkWidget *combo, const gint responseId, gpointer data);
GtkComboBox*                       createComboBox();
void                               addComboItem(GtkComboBox *combo, GtkTreeIter *parent, const int val, const char *text, const int initVal);

const char*                        getSaveFileName(GtkWidget *widget, const char *currentName, const char *defaultPath, const char *defaultExtension, const char *title);
char*                              getLoadFileName(GtkWidget *widget, const char *defaultPath, const char *title);

void                               enableMenuAction(GtkActionGroup *action_group, const char *actionName, const gboolean enable);
void                               setToggleMenuStatus(GtkActionGroup *action_group, const char *actionName, const gboolean active);
gboolean                           getToggleMenuStatus(GtkActionGroup *action_group, const char *actionName);
void                               setRadioMenuStatus(GtkActionGroup *action_group, const char *actionName, const gint value);

GtkWidget*                         externalCommand(const char *command, const char *progName, GtkWidget *widget, GError **error);
GString*                           getExternalCommandOutput(const char *command, GError **error);
GtkWidget*                         displayFetchResults(const char *title, const char *displayText, GtkWidget *widget, GtkWidget *dialog, GtkTextBuffer **textBuffer);

int                                scrollBarWidth();
void                               getTextSize(GtkWidget *widget, const char *text, int *width, int *height);
int                                getTextWidth(GtkWidget *widget, const char *text);
int                                getTextHeight(GtkWidget *widget, const char *text);

void                               getScreenSizeFraction(GtkWidget *widget, const double widthFraction, const double heightFraction, int *widthOut, int *heightOut);

GtkWidget*                         createTextEntryFromInt(GtkWidget *widget,
                                                          GtkTable *table, 
                                                          const int row,
                                                          const int col,
                                                          const int xpad,
                                                          const int ypad,
                                                          const char *mnemonic,
                                                          const int value,
                                                          BlxResponseCallback callback);

GtkWidget*                         createTextEntryFromDouble(GtkWidget *widget,
                                                             GtkTable *table, 
                                                             const int row,
                                                             const int col,
                                                             const int xpad,
                                                             const int ypad,
                                                             const char *mnemonic,
                                                             const double value,
                                                             BlxResponseCallback callback);

void                               drawHScale(GtkWidget *widget, 
                                              GdkDrawable *drawable,
                                              const IntRange* const range,
                                              const GdkRectangle* const rect,
                                              const int widthPerVal,       
                                              const int majorTickInterval, 
                                              const int labelInterval,     
                                              const int minorTickHeight,
                                              const int majorTickHeight);

const char*                        getStringFromTextEntry(GtkEntry *entry);
void                               widgetSetFontSize(GtkWidget *widget, gpointer data);
void                               widgetSetFontSizeAndCheck(GtkWidget *belvuAlignment, const int newSize);

GtkRadioButton*                    createRadioButton(GtkTable *table,
                                                     const int col,
                                                     const int row,
                                                      GtkRadioButton *existingButton,
                                                     const char *mnemonic,
                                                     const gboolean isActive,
                                                     const gboolean createTextEntry,
                                                     const gboolean multiline,
                                                     BlxResponseCallback callbackFunc,
                                                     GtkWidget *blxWindow,
                                                     GSList **entryList);

const char*                        getSystemTempDir();     

void                               errorHandler(const int sig); 

/* seqtoolsWebBrowser.c */
gboolean                           seqtoolsLaunchWebBrowser(const char *link, GError **error);


/* translate.c */
char*                              blxTranslate(const char *seq, char **code);
void                               blxComplement(char *seq) ;    
char*                              revComplement(char *comp, char *seq) ;
char                               complementChar(const char inputChar, GError **error);


void    gtk_text_buffer_insert_markup             (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup);

void    gtk_text_buffer_insert_markup_with_tag    (GtkTextBuffer *buffer,
                                                   GtkTextIter   *iter,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup_with_tag       (GtkTextBuffer *buffer,
                                                   const gchar   *markup,
                                                   GtkTextTag    *tag);

void    gtk_text_buffer_set_markup                (GtkTextBuffer *buffer,
                                                   const gchar   *markup);


#endif /* _utilities_h_included_ */
