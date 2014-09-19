/*  File: detailview.c
 *  Author: Gemma Barson, 2009-11-23
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
 * Description: The "detail view" section shows the sequence data in detail,
 *              down to the actual nucleotides and peptides. There is a 
 *              separate section to show the alignments for each strand and
 *              (in protein mode) reading-frame. Each section is a list of
 *              alignments (called a "tree" in the code - see detailviewtree.h)
 *
 *              Since there are 6 combinations of strand and reading frame for
 *              proteins, only one strand is shown at the time and the user can
 *              "toggle" which strand is currently shown (active). In DNA mode,
 *              we show both strands together, and toggling the active strand
 *              controls which strand is shown at the top (and which direction
 *              sequences are displayed).
 *----------------------------------------------------------------------------
 */

#ifndef _detail_view_included_
#define _detail_view_included_

#include <gtk/gtk.h>
#include <seqtoolsUtils/utilities.h>
#include <blixemApp/blxwindow.h>


#define HEADER_CONTAINER_NAME           "header container"
#define SNP_TRACK_HEADER_NAME           "SNP track header"
#define DNA_TRACK_HEADER_NAME           "DNA track header"
#define DETAIL_VIEW_STATUSBAR_CONTEXT   "statusBarCtx"


/* Column name to use when multiple, duplicate reads with different names are
 * shown on the same row. Note that this is a printf format that must take the
 * number of reads as an integer argument. */
#define DUPLICATE_READS_COLUMN_NAME               "(%d) matches" 
#define DUPLICATE_READS_COLUMN_NAME_SGL           "(%d) match" /* as above but for when there is just one read */


/* This struct contains info about canonical splice sites */
typedef struct _BlxSpliceSite
  {
    char donorSite[3];                  /* The bases expected at the donor splice site */
    char acceptorSite[3];               /* The bases expected at the acceptor splice site */

    char donorSiteRev[3];               /* Same as donorSite but reversed */
    char acceptorSiteRev[3];            /* Same as acceptorSite but reversed */
    
    gboolean bothReqd;                  /* Whether both donor and acceptor sites must match in order to be considered canonical */
  } BlxSpliceSite;


/* Essential info required by the the detail view */
typedef struct _DetailViewProperties
  {
    GtkWidget *blxWindow;         /* The main blixem window that this view belongs to */
    GtkCellRenderer *renderer;    /* The cell renderer that renders the sequences */
    GtkAdjustment *adjustment;    /* The scroll adjustment control for the detail view */

    GtkWidget *header;            /* Contains all the widgets in the detail view header */
    GtkWidget *feedbackBox;       /* A text box that feeds back info to the user about the currently selected items */
    GtkWidget *statusBar;         /* A status bar that feeds back info to the user about the currently moused-over items */
    GList *columnList;            /* A list of details about all the columns in the detail view */    
    BlxColumnId* sortColumns;     /* Array of columns to sort by, in order of priority. The length of this array will be set to the same length as columnList */
    
    GList *fwdStrandTrees;        /* A list of all the trees that show the forward strand of the ref seq */
    GList *revStrandTrees;        /* A list of all the trees that show the reverse strand of the ref seq */
    
    int cellXPadding;             /* The x padding between the tree cell background area and their drawing area */
    int cellYPadding;             /* The y padding between the tree cell background area and their drawing area */
        
    IntRange displayRange;        /* The currently-displayed range of bases in the reference sequence */
    
    gboolean selectedBaseSet;     /* indicates whether the selected base index below is set or not */
    int selectedBaseIdx;          /* The currently-selected index in the display range */
    int selectedFrame;            /* The reading frame to display selected bases for */
    int selectedBaseNum;          /* The currently-selected base within the selected reading frame */
    int selectedDnaBaseIdx;       /* The currently-selected index in terms of the DNA sequence */

    int clickedBaseIdx;           /* Stores the index the user right clicked on (used when copying a range of ref seq) */

    BlxStrand selectedStrand;     /* BlxStrand of the tree that the last-selected  */
    PangoFontDescription *fontDesc; /* The fixed-width font that will be used to display the alignments */

    int snpConnectorHeight;       /* The height of the connector between the SNP track and the DNA base track */
    int numUnalignedBases;        /* If the display-unaligned-sequence option is on, this specifies how many additional bases to show at each end of the alignment */

    /* Cached font sizes, needed often for calculations. */
    gdouble charHeight;
    gdouble charWidth;
        
    int exonBoundaryLineWidth;                 /* line width for exon boundaries */
    GdkLineStyle exonBoundaryLineStyle;        /* line style for exon boundaries */
    GdkLineStyle exonBoundaryLineStylePartial; /* line style for exon boundaries (where the boundary is part-way through a codon) */
    
    GSList *spliceSites;           /* List of splice sites that can be found and highlighted by Blixem */
  } DetailViewProperties;


typedef struct _DrawBaseData
{
  int dnaIdx;
  char baseChar;
  BlxStrand strand;
  int frame;
  BlxSeqType seqType;
  gboolean topToBottom;         /* true if we're displaying bases top-to-bottom instead of left-to-right */
  gboolean displayIdxSelected;  /* whether to use default background color or leave blank */
  gboolean dnaIdxSelected;      /* true if this base is the currently-selected dna index */
  gboolean showBackground;      /* true if we should draw the background colour */
  gboolean highlightSnps;       /* whether to highlight variations */
  gboolean showCodons;          /* whether to highlight DNA bases within the selected codon, for protein matches */
  BlxColorId defaultBgColor;    /* the default background color for the header */
  GdkColor *fillColor;
  GdkColor *outlineColor;
  gboolean drawStart;
  gboolean drawEnd;
  gboolean drawJoiningLines;
  gboolean shadeBackground;
} DrawBaseData;



/* Public function declarations */
int                     detailViewGetNumFrames(GtkWidget *detailView);
IntRange*               detailViewGetDisplayRange(GtkWidget *detailView);
int                     detailViewGetClickedBaseIdx(GtkWidget *detailView);
int                     detailViewGetSelectedBaseIdx(GtkWidget *detailView);
gboolean                detailViewGetSelectedBaseSet(GtkWidget *detailView);
int                     detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView);
int                     detailViewGetOldSelectedBaseIdx(GtkWidget *detailView);
GtkAdjustment*          detailViewGetAdjustment(GtkWidget *detailView);
GtkWidget*              detailViewGetTree(GtkWidget *detailView, const BlxStrand strand, const int frame);
GtkWidget*              detailViewGetTreeContainer(GtkWidget *detailView, const BlxStrand strand, const int frame);
gboolean                detailViewGetDisplayRev(GtkWidget *detailView);
PangoFontDescription*   detailViewGetFontDesc(GtkWidget *detailView);
GtkCellRenderer*        detailViewGetRenderer(GtkWidget *detailView);
int                     detailViewGetCellXPadding(GtkWidget *detailView);
int                     detailViewGetCellYPadding(GtkWidget *detailView);
BlxSeqType              detailViewGetSeqType(GtkWidget *detailView);
IntRange*               detailViewGetRefSeqRange(GtkWidget *detailView);
GtkWidget*              detailViewGetBlxWindow(GtkWidget *detailView);
gdouble                 detailViewGetCharWidth(GtkWidget *detailView);
gdouble                 detailViewGetCharHeight(GtkWidget *detailView);
int                     detailViewGetNumUnalignedBases(GtkWidget *detailView);
BlxColumnId*            detailViewGetSortColumns(GtkWidget *detailView);
GList*                  detailViewGetColumnList(GtkWidget *detailView);
GType*                  columnListGetTypes(GList *columnList);
int                     detailViewGetActiveFrame(GtkWidget *detailView);
BlxStrand               detailViewGetSelectedStrand(GtkWidget *detailView);
void                    detailViewSetSelectedStrand(GtkWidget *detailView, BlxStrand strand);

DetailViewProperties*   detailViewGetProperties(GtkWidget *widget);

int                     detailViewGetColumnWidth(GtkWidget *detailView, const BlxColumnId columnId);
const char*             detailViewGetColumnTitle(GtkWidget *detailView, const BlxColumnId columnId);
void                    detailViewGetColumnXCoords(DetailViewProperties *properties, const BlxColumnId columnId, IntRange *xRange);
gboolean                detailViewShowColumn(BlxColumnInfo *columnInfo);
void                    detailViewSaveProperties(GtkWidget *detailView, GKeyFile *key_file);
void                    detailViewResetColumnWidths(GtkWidget *detailView);

int                     getBaseIndexAtColCoords(const int x, const int y, const gdouble charWidth, const IntRange* const displayRange);

MSP*                    prevMatch(GtkWidget *detailView, GList *seqList);
MSP*                    nextMatch(GtkWidget *detailView, GList *seqList);
MSP*                    firstMatch(GtkWidget *detailView, GList *seqList);
MSP*                    lastMatch(GtkWidget *detailView, GList *seqList);
void                    goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType);
void                    setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void                    setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType);
void                    scrollDetailViewLeft1(GtkWidget *detailView);
void                    scrollDetailViewRight1(GtkWidget *detailView);
void                    scrollDetailViewLeftStep(GtkWidget *detailView);
void                    scrollDetailViewRightStep(GtkWidget *detailView);
void                    scrollDetailViewLeftPage(GtkWidget *detailView);
void                    scrollDetailViewRightPage(GtkWidget *detailView);

void                    detailViewSetSortColumn(GtkWidget *detailView, const BlxColumnId sortColumn);

void                    zoomDetailView(GtkWidget *detailView, const gboolean zoomIn);

void                    updateDynamicColumnWidths(GtkWidget *detailView);
void                    refilterDetailView(GtkWidget *detailView, const IntRange* const oldRange);

void                    detailViewSetSelectedBaseIdx(GtkWidget *detailView, 
                                                     const int selectedBaseIdx, 
                                                     const int frame, 
                                                     const int baseNum, 
                                                     const gboolean allowScroll, 
                                                     const gboolean scrollMinimum);

void                    detailViewSetSelectedDnaBaseIdx(GtkWidget *detailView,
                                                        const int selectedDnaBaseIdx,
                                                        const int frame,
                                                        const gboolean allowScroll,
                                                        const gboolean scrollMinimum);

void                    detailViewScrollToKeepInRange(GtkWidget *detailView, const IntRange* const range);

void                    detailViewUnsetSelectedBaseIdx(GtkWidget *detailView);
void                    detailViewSetActiveFrame(GtkWidget *detailView, const int frame);
void                    detailViewResortTrees(GtkWidget *detailView);

void                    updateFeedbackBox(GtkWidget *detailView);
void                    clearFeedbackArea(GtkWidget *detailView);
void                    updateFeedbackAreaNucleotide(GtkWidget *detailView, const int dnaIdx, const BlxStrand strand);
void                    toggleStrand(GtkWidget *detailView);

void                    detailViewAddMspData(GtkWidget *detailView, MSP *mspList, GList *seqList);

void                    updateDetailViewFontDesc(GtkWidget *detailView);
void                    updateDetailViewRange(GtkWidget *detailView);
void                    resizeDetailViewHeaders(GtkWidget *detailView);
void                    detailViewRedrawAll(GtkWidget *detailView);
void                    detailViewRefreshAllHeaders(GtkWidget *detailView);

void                    detailViewUpdateSquashMatches(GtkWidget *detailView, const gboolean squash);
void                    detailViewUpdateSortInverted(GtkWidget *detailView, const gboolean invert);
void                    detailViewUpdateShowSnpTrack(GtkWidget *detailView, const gboolean showSnpTrack);
void                    detailViewUpdateMspLengths(GtkWidget *detailView, const int numUnalignedBases);

void                    detailViewSetNumUnalignedBases(GtkWidget *detailView, const int numBases);
void                    detailViewToggleSnpTrack(GtkWidget *detailView);

GtkWidget*              createSnpTrackHeader(GtkBox *parent, GtkWidget *detailView, const BlxStrand strand);
void                    refreshTextHeader(GtkWidget *widget, gpointer data);
gboolean                onExposeDnaTrack(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);

void                    selectClickedSnp(GtkWidget *snpTrack,
                                         GtkWidget *colHeader,
                                         GtkWidget *detailView, 
                                         const int xIn, 
                                         const int yIn, 
                                         const gboolean expandSnps,
                                         const int clickedBase);

void                    seqColHeaderSetRow(GtkWidget *header, const int frame);
int                     seqColHeaderGetRow(GtkWidget *header);
int                     seqColHeaderGetBase(GtkWidget *header, const int frame, const int numFrames);

GHashTable*             getRefSeqBasesToHighlight(GtkWidget *detailView, const IntRange* const displayRange, const BlxSeqType seqType, const BlxStrand strand);

void                    drawColumnSeparatorLine(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const BlxViewContext *bc);
gboolean                onExposeGenericHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);

gint                    sortByColumnCompareFunc(GList *mspGList1,
                                                GList *mspGList2,
                                                GtkWidget *detailView, 
                                                const BlxColumnId sortColumn);

void                    drawHeaderChar(BlxViewContext *bc,
                                       DetailViewProperties *properties,
                                       GdkDrawable *drawable,
                                       GdkGC *gc,
                                       const int x,
                                       const int y,
                                       GHashTable *intronBases,
                                       DrawBaseData *baseData);

GtkWidget*              createDetailView(GtkWidget *blxWindow,
                                         GtkContainer *parent,
                                         GtkWidget *toolbar,
                                         GtkAdjustment *adjustment, 
                                         GtkWidget *fwdStrandGrid, 
                                         GtkWidget *revStrandGrid,
                                         MSP *mspList,
                                         GList *columnList,
                                         BlxBlastMode mode,
                                         BlxSeqType seqType,
                                         int numFrames,
                                         const char* const refSeqName,
                                         const int startCoord,
                                         const gboolean sortInverted,
                                         const BlxColumnId sortColumn,
                                         const gboolean optionalDataLoaded,
                                         char *windowColor);

GtkWidget*              createDetailViewScrollBar(GtkAdjustment *adjustment, 
                                                  GtkWidget *blxWindow);


#endif /* _detail_view_included_ */
