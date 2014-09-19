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
 * Description: See detailview.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/detailview.h>
#include <blixemApp/detailviewtree.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/exonview.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <gtk/gtk.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define DETAIL_VIEW_TOOLBAR_NAME        "DetailViewToolbarName"
#define DETAIL_VIEW_WIDGET_NAME         "DetailViewWidget"
#define SORT_BY_NAME_STRING             "Name"
#define SORT_BY_SCORE_STRING            "Score"
#define SORT_BY_ID_STRING               "Identity"
#define SORT_BY_POS_STRING              "Position"
#define SORT_BY_GROUP_ORDER_STRING      "Group"
#define FONT_INCREMENT_SIZE             1
#define NO_SUBJECT_SELECTED_TEXT        "<no subject selected>"
#define MULTIPLE_SUBJECTS_SELECTED_TEXT "<multiple subjects selected>"
#define MULTIPLE_VARIATIONS_TEXT        "<multiple variations>"
#define MULTIPLE_POLYA_SIGNALS_TEXT     "<multiple polyA signals>"
#define MULTIPLE_POLYA_SITES_TEXT       "<multiple polyA sites>"
#define DEFAULT_SNP_CONNECTOR_HEIGHT    0
#define DEFAULT_NUM_UNALIGNED_BASES     5     /* the default number of additional bases to show if displaying unaligned parts of the match sequence */
#define POLYA_SIG_BASES_UPSTREAM        50    /* the number of bases upstream from a polyA tail to search for polyA signals */
#define POLYA_SIGNAL                    "aataaa"

#define SETTING_NAME_NUM_UNALIGNED_BASES "num-unaligned-bases"


typedef struct 
  {
    const int startDnaIdx;              /* the DNA coord to start searching from */
    const gboolean searchRight;         /* search towards the right or left */
    const int searchDirection;          /* multiplier to add/subtract depending on whether searching right/left */
    const gboolean displayRev;          /* true if the display is reversed */
    const BlxSeqType seqType;           /* whether viewing DNA or peptide seqs */
    const int numFrames;                /* number of reading frames */
    const IntRange* const refSeqRange;  /* the full range of the reference sequence */
    GList *seqList;                     /* only search matches in these sequences */

    int offset;                         /* the offset of the found MSP */
    int foundFrame;                     /* which ref seq frame the MSP we chose is in */
    MSP *foundMsp;                      /* the found msp */
  } MatchSearchData;


/* Utility struct to pass data to a recursive function */
typedef struct _RecursiveFuncData
{
  const char *widgetName;   /* the name of widgets that the function will be called on */
  GtkCallback callbackFunc; /* the function to call */
  gpointer callbackData;    /* user data to pass to the function */
} RecursiveFuncData;


/* Local function declarations */
static BlxViewContext*        detailViewGetContext(GtkWidget *detailView);
static GtkWidget*             detailViewGetBigPicture(GtkWidget *detailView);
static GtkWidget*             detailViewGetHeader(GtkWidget *detailView);
static GtkWidget*             detailViewGetFeedbackBox(GtkWidget *detailView);
static int                    detailViewGetSnpConnectorHeight(GtkWidget *detailView);
static void                   refreshDetailViewHeaders(GtkWidget *detailView);

static void                   snpTrackSetStrand(GtkWidget *snpTrack, const BlxStrand strand);
static BlxStrand              snpTrackGetStrand(GtkWidget *snpTrack, GtkWidget *detailView);
static void                   getVariationDisplayRange(const MSP *msp, const gboolean expand, const BlxSeqType seqType, const int numFrames, const gboolean displayRev, const int activeFrame, const IntRange* const refSeqRange, IntRange *displayRange, IntRange *expandedRange);
static void                   recalculateSnpTrackBorders(GtkWidget *snpTrack, gpointer data);

static void                   detailViewCacheFontSize(GtkWidget *detailView, gdouble charWidth, gdouble charHeight);
static gboolean               widgetIsTree(GtkWidget *widget);
static gboolean               widgetIsTreeContainer(GtkWidget *widget);
static void                   updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc);
static void                   setDetailViewScrollPos(GtkAdjustment *adjustment, int value);
static const char*            spliceSiteGetBases(const BlxSpliceSite *spliceSite, const gboolean donor, const gboolean reverse);
static int                    getNumSnpTrackRows(const BlxViewContext *bc, DetailViewProperties *properties, const BlxStrand strand, const int frame);
static int                    getVariationRowNumber(const IntRange* const rangeIn, const int numRows, GSList **rows);
static void                   freeRowsList(GSList *rows);

static gboolean               coordAffectedByVariation(const int dnaIdx,
                                                       const BlxStrand strand, 
                                                       BlxViewContext *bc,
                                                       const MSP **msp,
                                                       gboolean *drawStartBoundary, 
                                                       gboolean *drawEndBoundary, 
                                                       gboolean *drawJoiningLines, 
                                                       gboolean *drawBackground,
                                                       gboolean *multipleVariations);

static gboolean               coordAffectedByPolyASignal(const int dnaIdx,
                                                         const BlxStrand strand, 
                                                         BlxViewContext *bc,
                                                         const MSP **mspOut,
                                                         gboolean *multipleVariations);

static gboolean                coordAffectedByPolyASite(const int dnaIdx,
                                                        const BlxStrand strand, 
                                                        BlxViewContext *bc,
                                                        const MSP **mspOut,
                                                        gboolean *multipleVariations);


/***********************************************************
 *                     Utility functions                   *
 ***********************************************************/

/* Recursively call a given function on a given widget and its
 * children if it/they have the given name. The callback function,
 * callback data and widget name are passed in a RecursiveFuncData
 * struct */
void callFuncOnChildren(GtkWidget *widget, gpointer data)
{
  RecursiveFuncData *funcData = (RecursiveFuncData*)data;

  if (stringsEqual(gtk_widget_get_name(widget), funcData->widgetName, TRUE))
    funcData->callbackFunc(widget, funcData->callbackData);
  
  if (GTK_IS_CONTAINER(widget))
    gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnChildren, data);
}


/* Refresh headers for the detail view and its children */
void detailViewRefreshAllHeaders(GtkWidget *detailView)
{
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  
  RecursiveFuncData data = {SNP_TRACK_HEADER_NAME, recalculateSnpTrackBorders, detailView};    
  callFuncOnChildren(detailView, &data);

  gtk_widget_queue_draw(detailView);
}


void detailViewRedrawAll(GtkWidget *detailView)
{
  /* Redraw all the trees */
  callFuncOnAllDetailViewTrees(detailView, widgetClearCachedDrawable, NULL);

  /* Recalculate the size of the snp track headers */
  RecursiveFuncData data = {SNP_TRACK_HEADER_NAME, recalculateSnpTrackBorders, detailView};  
  callFuncOnChildren(detailView, &data);

  gtk_widget_queue_draw(detailView);
}


/* Save any user-settings that are stored in the detail-view properties */
void detailViewSaveProperties(GtkWidget *detailView, GKeyFile *key_file)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  g_key_file_set_integer(key_file, SETTINGS_GROUP, SETTING_NAME_NUM_UNALIGNED_BASES, properties->numUnalignedBases);
  saveColumnWidths(detailViewGetColumnList(detailView), key_file);
  saveSummaryColumns(detailViewGetColumnList(detailView), key_file);
}


/* Scroll the detail view so that the given coord is at the start of the display
 * range (within bounds). the coord should be in terms of display coords */
void setDetailViewStartIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, coord);
}


/* Scroll the detail view so that the given coord is at the end of the display
 * range (within bounds). the coord should be in terms of display coords */
void setDetailViewEndIdx(GtkWidget *detailView, int coord, const BlxSeqType coordSeqType)
{
  /* Get the new start coord */
  const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
  const int displayLen = getRangeLength(displayRange);
  int newStart = coord - displayLen + 1;

  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, newStart);
}


/* Update the scroll position of the adjustment to the given value. Does bounds checking. */
static void setDetailViewScrollPos(GtkAdjustment *adjustment, int value)
{  
  /* bounds checking */
  int maxValue = adjustment->upper - adjustment->page_size + 1;
  
  if (value > maxValue)
    {
      value = maxValue;
    }
  
  if (value < adjustment->lower)
    {
      value = adjustment->lower;
    }
  
  adjustment->value = value;
  
  /* Emit notification that the scroll pos has changed */
   gtk_adjustment_value_changed(adjustment);
}


/* Scroll left/right by one base. */
void scrollDetailViewLeft1(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value - 1);
}

void scrollDetailViewRight1(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + 1);
}

/* Scroll by one step increment */
void scrollDetailViewLeftStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value - adjustment->step_increment);
}

void scrollDetailViewRightStep(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + adjustment->step_increment);
}

/* Scroll by one page size */
void scrollDetailViewLeftPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value - adjustment->page_increment);
}

void scrollDetailViewRightPage(GtkWidget *detailView)
{
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);
  setDetailViewScrollPos(adjustment, adjustment->value + adjustment->page_increment);
}


/* Calculate the number of bases that can be displayed in the sequence column */
static int calcNumBasesInSequenceColumn(DetailViewProperties *properties)
{
  int numChars = 0;
  
  /* Find the width of the sequence column */
  int colWidth = UNSET_INT;
  GList *listItem = blxWindowGetColumnList(properties->blxWindow);
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);
      if (columnInfo && columnInfo->columnId == BLXCOL_SEQUENCE)
        {
          colWidth = columnInfo->width;
          break;
        }
    }
  
  if (colWidth >= 0)
    {
      /* Don't include the cell padding area */
      GtkCellRenderer *renderer = properties->renderer;
      colWidth -= (2 * renderer->xpad) + (2 * renderer->xalign);
      
      /* Return the number of whole characters that fit in the column. */
      numChars = (int)((gdouble)colWidth / properties->charWidth);
    }
  
  return numChars;
}


/* This should be called when the width of the sequence column has changed (or the
 * size of the font has changed). This function will adjust the scroll range of our
 * custom scroll adjustment so that it represents the range that can be displayed 
 * in the new column width. */
void updateDetailViewRange(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  if (properties->adjustment)
    {
      int newPageSize = calcNumBasesInSequenceColumn(properties);
      
      /* Only trigger updates if things have actually changed */
      if (newPageSize != properties->adjustment->page_size)
        {
          properties->adjustment->page_size = newPageSize;
          properties->adjustment->page_increment = newPageSize;
          
          ///* Reset the display range so that it is between the scrollbar min and max. Try to keep
//         * it centred on the same base. The base index is in terms of the display range coords, 
//         * so the sequence type of the coord is whatever the display sequence type is. */
//        const IntRange const *displayRange = &properties->displayRange;
//
//        /* First time through, both coords are set to the initial start coord. So, if
//         * they are the same, just use that coord as the start coord */
//        int newStart = displayRange->min;
//        
//        if (displayRange->min != displayRange->max)
//          {
//            int centre = getRangeCentre(displayRange);
//            int offset = roundNearest((double)properties->adjustment->page_size / 2.0);
//            newStart = centre - offset;
//          }
//            
//        const BlxSeqType seqType = blxWindowGetSeqType(properties->blxWindow);
//        setDetailViewStartIdx(detailView, newStart, seqType);
          
          gtk_adjustment_changed(properties->adjustment); /* signal that the scroll range has changed */
        }
    }
}


/* Find a child widget of this widget that is a paned widget. Returns null if it does not
 * have one. If there are two, it returns the first one it finds. */
static GtkWidget *findNestedPanedWidget(GtkContainer *parentWidget)
{
  GtkWidget *result = NULL;
  
  GList *children = gtk_container_get_children(parentWidget);
  GList *child = children;
  
  while (child)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);
      if (GTK_IS_PANED(childWidget))
        {
          result = childWidget;
          break;
        }
      
      child = child->next;
    }
  
  g_list_free(children);

  return result;
}


/* Add all trees in the given list to the detail view. Note that the treeList may not contain
 * actual trees, it may contain their containers (e.g. if they live in scrolled windows). 'first'
 * indicates that this is the first (set of) tree(s) added. This information is used to determine
 * which pane to put the tree in and/or whether the first tree in the list should have headers. */
static void addTreesToDetailViewPane(GtkPaned *panedWin, 
                                     GList *treeList, 
                                     const gboolean first)
{
  int numTrees = g_list_length(treeList);
  
  if (numTrees == 1)
    {
      /* Two panes, one tree. Use 'first' flags to decide which pane to put it in */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      
      if (first)
        {
          gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);
        }
      else
        {
          gtk_paned_pack2(panedWin, tree1, TRUE, TRUE);
        }
    }
  else if (numTrees == 2)
    {
      /* Two panes, two trees. Easy - put one in each. */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      GtkWidget *tree2 = GTK_WIDGET(treeList->next->data);
      
      gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);
      gtk_paned_pack2(panedWin, tree2, TRUE, TRUE);
    }
  else if (numTrees > 2)
    {
      /* Two panes, three trees. Put one tree in pane 1. There should be a
       * nested widget in pane 2 that we can put the remaining trees in. */
      GtkWidget *tree1 = GTK_WIDGET(treeList->data);
      gtk_paned_pack1(panedWin, tree1, TRUE, TRUE);

      GtkWidget *nestedPanedWidget = findNestedPanedWidget(GTK_CONTAINER(panedWin));
      
      if (nestedPanedWidget)
        {
          /* Create a new list containing the remaining trees and call this 
           * function again on the nested paned widget. */
          GList *remainingTrees = NULL;
          GList *nextTree = treeList->next;
          while (nextTree)
            {
              remainingTrees = g_list_append(remainingTrees, nextTree->data);
              nextTree = nextTree->next;
            }
          
          /* Pass 'first' as false to make sure we don't add any more headers */
          addTreesToDetailViewPane(GTK_PANED(nestedPanedWidget), remainingTrees, FALSE);
        }
    }
}


/* Returns true if this widget is a detail-view-tree */
static gboolean widgetIsTree(GtkWidget *widget)
{
  const char *name = gtk_widget_get_name(widget);
  return (strcmp(name, DETAIL_VIEW_TREE_NAME) == 0);
}


/* Returns true if this widget is a detail-view-tree-container (that is, 
 * a container that has only one child that is a detail-view-tree.) */
static gboolean widgetIsTreeContainer(GtkWidget *widget)
{
  const char *name = gtk_widget_get_name(widget);
  return (strcmp(name, DETAIL_VIEW_TREE_CONTAINER_NAME) == 0);
}


/* Given a widget contains a single child that is a detail-view-tree, return
 * that tree. */
static GtkWidget *treeContainerGetTree(GtkContainer *container)
{
  GtkWidget *result = NULL;
  
  GList *children = gtk_container_get_children(container);
  GList *child = children;
  
  for ( ; child && !result; child = child->next)
    {
      if (GTK_IS_WIDGET(child->data))
        {
          GtkWidget *childWidget = GTK_WIDGET(child->data);
          
          if (widgetIsTree(childWidget))
            {
              result = childWidget;
            }
          else if (GTK_IS_CONTAINER(childWidget))
           {
             result = treeContainerGetTree(GTK_CONTAINER(childWidget));
           }
        }
    }

  g_list_free(children);
  
  return result;
}


/* This function removes all detail-view-trees from the given container widget. */
static void removeAllTreesFromContainer(GtkWidget *widget, gpointer data)
{
  /* See if this widget is a tree (or a container containing only a single tree - the
   * detail view contains the outermost single-tree container, e.g. typically the trees
   * will live in a vbox or scrolled window, and it is this container that the detail view
   * knows about) */
  if (widgetIsTree(widget) || widgetIsTreeContainer(widget))
    {
      GtkContainer *parent = GTK_CONTAINER(data);
      gtk_container_remove(parent, widget);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      /* It is a nested container, containing multiple children. Recurse. */
      gtk_container_foreach(GTK_CONTAINER(widget), removeAllTreesFromContainer, widget);
    }
}


/* This function removes the trees from the detail view and re-adds them in the
 * correct order according to the displayRev flag. It should be called every
 * time the strands are toggled. It assumes the trees are already in the 
 * detailView container, and that the properties have been set for all 3 widgets. */
static void refreshTreeOrder(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  BlxViewContext *bc = detailViewGetContext(detailView);
  gboolean toggled = detailViewGetDisplayRev(detailView);

  /* Get the direct parent of the trees */
  GtkWidget *detailViewContainer = getNamedChildWidget(detailView, DETAIL_VIEW_WIDGET_NAME);

  gtk_container_foreach(GTK_CONTAINER(detailView), removeAllTreesFromContainer, detailViewContainer);

  /* Extract the named detail-view widget, which is the direct parent of the trees. This 
   * should be a paned window. */
  GtkWidget *panedWin = getNamedChildWidget(detailView, DETAIL_VIEW_WIDGET_NAME);

  if (!GTK_IS_PANED(panedWin))
    {
      g_error("Unexpected detail view type: expected a paned widget. Could not add child trees.\n");
    }

  if (bc->seqType == BLXSEQ_DNA)
    {
      /* Add both trees - the order they're added will depend on whether the display is toggled or not. */
      addTreesToDetailViewPane(GTK_PANED(panedWin), properties->fwdStrandTrees, !toggled);
      addTreesToDetailViewPane(GTK_PANED(panedWin), properties->revStrandTrees, toggled);
    }
  else if (bc->seqType == BLXSEQ_PEPTIDE)
    {
      /* Only add one set of trees - the reverse strand if toggled, the forward strand if not. */
      GList *treeList = toggled ? properties->revStrandTrees : properties->fwdStrandTrees;
      addTreesToDetailViewPane(GTK_PANED(panedWin), treeList, TRUE);
    }
  
  /* Must show all child widgets because some of them may not have been in this parent before.
   * (Just calling gtk_widget_show on the individual trees doesn't seem to work.)
   * However, we then need to re-hide any that may have been previously hidden by the user. */
  gtk_widget_show_all(detailView);
  gtk_container_foreach(GTK_CONTAINER(detailView), hideUserHiddenWidget, NULL);
}


/* Refresh a header that contains fixed-width text that needs to be resized whenever a zoom
 * happens; updates the height of the widget and the font description, and clears its cached
 * drawable, if it has one. This just sets things up ready for the next expose event, where 
 * the actual drawing will take place. */
void refreshTextHeader(GtkWidget *header, gpointer data)
{
  const char *widgetName = gtk_widget_get_name(header);
  
  if (GTK_IS_LABEL(header))
    {
      /* For straightforward labels, just update the font size. */
      GtkWidget *detailView = GTK_WIDGET(data);
      
      PangoFontDescription *dvFont = detailViewGetFontDesc(detailView);
      int newSize = pango_font_description_get_size(dvFont);

      pango_font_description_set_size(header->style->font_desc, newSize);
      gtk_widget_modify_font(header, header->style->font_desc);
    }
  else if (!strcmp(widgetName, SNP_TRACK_HEADER_NAME) || !strcmp(widgetName, DNA_TRACK_HEADER_NAME))
    {
      /* Clear the cached drawable so that it gets recreated on the next expose */
      widgetClearCachedDrawable(header, NULL);

      /* Update the font and the widget height, in case the font-size has changed. */
      GtkWidget *detailView = GTK_WIDGET(data);
      gtk_widget_modify_font(header, detailViewGetFontDesc(detailView));

      const gdouble charHeight = detailViewGetCharHeight(detailView);
      
      if (!strcmp(widgetName, SNP_TRACK_HEADER_NAME))
        {
          /* SNP track */
          BlxViewContext *bc = detailViewGetContext(detailView);
        
          if (!bc->flags[BLXFLAG_SHOW_VARIATION_TRACK] || !bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS])
            {
              /* SNP track is hidden, so set the height to 0 */
              gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, 0);
              gtk_widget_set_size_request(header, -1, 0);
            }
          else
            {
              /* Set the height to the character height plus the connector line height */
              const int height = roundNearest(charHeight + (gdouble)detailViewGetSnpConnectorHeight(detailView));
              gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, height);
              gtk_widget_set_size_request(header, -1, height);
            }
        }
      else if (GTK_IS_LAYOUT(header))
        {
          /* Normal text header. Set the height to the character height */
          gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, roundNearest(charHeight));
          gtk_widget_set_size_request(header, -1, roundNearest(charHeight));
        }

      gtk_widget_queue_draw(header);
    }
  
  /* If this is a container of header widgets, recurse over each child */
  if (!strcmp(widgetName, HEADER_CONTAINER_NAME) && GTK_IS_CONTAINER(header))
    {
      gtk_container_foreach(GTK_CONTAINER(header), refreshTextHeader, data);
    }
}


/* Refresh the headers for the detail view. NOte: does not refresh tree view
 * headers; use detailViewRefreshAllHeaders if you want to do that */
static void refreshDetailViewHeaders(GtkWidget *detailView)
{
  /* Loop through all widgets in the header and call refreshTextHeader. This
   * updates the font etc. if it is a type of widget that requires that. */
  GtkWidget *header = detailViewGetHeader(detailView);
  gtk_container_foreach(GTK_CONTAINER(header), refreshTextHeader, detailView);
  
  /* Loop through all columns and call the individual refresh callbacks for
   * each of their headers. This updates the specific information within the header. */
  GList *columnList = detailViewGetColumnList(detailView);
  GList *column = columnList;
  
  for ( ; column; column = column->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)column->data;
      if (columnInfo)
        {
          if (columnInfo->headerWidget && columnInfo->refreshFunc)
            {
              columnInfo->refreshFunc(columnInfo->headerWidget, detailView);
            }
        }
      else
        {
          g_critical("Invalid column data for detail view header; header may not be refreshed correctly.\n");
        }
    }
}


/* Resize the detail view header widgets. Should be called whenever a column is resized. */
void resizeDetailViewHeaders(GtkWidget *detailView)
{
  GList *listItem = detailViewGetColumnList(detailView);
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)listItem->data;

      /* For the sequence column, don't set the size request, or we won't be able
       * to shrink the window. (The sequence col header will be resized dynamically.) */
      if (columnInfo->columnId != BLXCOL_SEQUENCE)
        {
          /* For other columns, we can set the size request: they're small enough
           * that we can live without the need to shrink below their sum widths. */
          if (showColumn(columnInfo))
            {
              widgetSetHidden(columnInfo->headerWidget, FALSE);
              gtk_widget_set_size_request(columnInfo->headerWidget, columnInfo->width, -1);
            }
          else
            {
              /* Zero width: we must hide the widget because it always seems to have
               * some width even if we set the size request to 0 */
              widgetSetHidden(columnInfo->headerWidget, TRUE);
            }
        }
    }
}


/* Update the font description for all relevant components of the detail view */
void updateDetailViewFontDesc(GtkWidget *detailView)
{
  PangoFontDescription *fontDesc = detailViewGetFontDesc(detailView);
  
  gtk_widget_modify_font(detailView, fontDesc);
  callFuncOnAllDetailViewTrees(detailView, (GtkCallback)gtk_widget_modify_font, fontDesc);
  
  updateCellRendererFont(detailView, fontDesc);
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  callFuncOnAllDetailViewTrees(detailView, treeUpdateFontSize, NULL);
  detailViewRedrawAll(detailView);
}


static void incrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) + FONT_INCREMENT_SIZE;
  
  if (newSize <= MAX_FONT_SIZE)
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView);
    }
}


static void decrementFontSize(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int newSize = (pango_font_description_get_size(properties->fontDesc) / PANGO_SCALE) - FONT_INCREMENT_SIZE;
  
  /* Note that the font needs to be big enough to cover the vertical separators and padding around the
   * cell, otherwise we will end up with visible gaps between rows. */
  if (newSize >= MIN_FONT_SIZE && newSize >= (detailViewGetCellYPadding(detailView) * 2))
    {
      pango_font_description_set_size(properties->fontDesc, newSize * PANGO_SCALE);
      updateDetailViewFontDesc(detailView);
    }
}


/* Zoom the detail view in/out */
void zoomDetailView(GtkWidget *detailView, const gboolean zoomIn)
{
  if (zoomIn)
    {
      incrementFontSize(detailView);
    }
  else
    {
      decrementFontSize(detailView);
    }
  
  updateDetailViewRange(detailView);
}


/* Get the text displayed in the user feedback box based on the given MSPs sequence name
 * (if an MSP is given), and also the currently-selected base index (if there is one). 
 * The string returned by this function must be free'd with g_free. */
static char* getFeedbackText(GtkWidget *detailView, const BlxSequence *seq, const int numSeqsSelected)
{
  /* The info we need to find... */
  int qIdx = UNSET_INT; /* index into the ref sequence. Ref seq is always a DNA seq */
  int sIdx = UNSET_INT; /* index into the match sequence. Will be coords into the peptide sequence if showing peptide matches */
  int sLen = UNSET_INT; /* the length of the match sequence */

  DetailViewProperties *properties = detailViewGetProperties(detailView);
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  /* Get the selected base index. */
  if (properties->selectedBaseSet)
    qIdx = properties->selectedDnaBaseIdx;
  
  /* Find the sequence name text (or some default text to indicate that a sequence is not selected) */
  const char *noSeqText = numSeqsSelected > 0 ? MULTIPLE_SUBJECTS_SELECTED_TEXT : NO_SUBJECT_SELECTED_TEXT;
  
  if (seq)
    {
      if (g_list_length(seq->mspList) > 0)
        {
          MSP *firstMsp = (MSP*)(seq->mspList->data);
          
          if (mspGetMatchSeq(firstMsp))
            {
              sLen = strlen(mspGetMatchSeq(firstMsp));
            }

          /* If a q index is selected, see if there is a valid base at that index 
           * for any of the MSPs for the selected sequence. */
          if (properties->selectedBaseSet)
            {
              GList *mspListItem = seq->mspList;
              const int numUnalignedBases = detailViewGetNumUnalignedBases(detailView);
              
              for ( ; mspListItem; mspListItem = mspListItem->next)
                {
                  MSP *msp = (MSP*)(mspListItem->data);
                  
                  sIdx = mspGetMatchCoord(msp, qIdx, TRUE, numUnalignedBases, bc);

                  if (sIdx != UNSET_INT)
                    {
                      break;
                    }
                }
            }
        }
    }

  /* Add all the bits into a text string */
  GString *resultString = g_string_sized_new(200); /* will be extended if we need more space */
  
  if (properties->selectedBaseSet)
    {
      /* Negate the coord for the display, if necessary */
      int coord = (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS] ? -1 * qIdx : qIdx);
      g_string_printf(resultString, "%d   ", coord);
    }
  
  if (seq)
    {
      const char *seqName = blxSequenceGetName(seq);
      const char *sequence = blxSequenceGetSequence(seq);

      if (seqName)
        g_string_append_printf(resultString, "%s", seqName);
        
      if (seq->type == BLXSEQUENCE_VARIATION && sequence)
        g_string_append_printf(resultString, " : %s", sequence);
    }
  else if (qIdx != UNSET_INT)
    {
      g_string_append_printf(resultString, "%s", noSeqText); 
    }
    
  if (sLen != UNSET_INT && (!seq || seq->type != BLXSEQUENCE_VARIATION))
    {
      g_string_append_printf(resultString, "(%d)", sLen);
    }

  if (sIdx != UNSET_INT && (!seq || seq->type != BLXSEQUENCE_VARIATION))
    {
      g_string_append_printf(resultString, " : %d", sIdx);
    }
  
  char *messageText = g_string_free(resultString, FALSE);
  
  return messageText;
}


/* Updates the feedback box with info about any currently-selected MSPs. This 
 * currently assumes single-MSP selections, but could be extended in the future 
 * to display, say, summary information about multiple MSPs. */
void updateFeedbackBox(GtkWidget *detailView)
{
  char *messageText = NULL;

  BlxViewContext *bc = detailViewGetContext(detailView);
  const int numSeqsSelected = g_list_length(bc->selectedSeqs);
  
  if (numSeqsSelected == 1) /* currently we only properly handle single sequence selection */
    {
      const BlxSequence *seq = (const BlxSequence*)(bc->selectedSeqs->data);
      messageText = getFeedbackText(detailView, seq, numSeqsSelected);
    }
  else
    {
      /* 0 or multiple MSPs selected. Just see if a base index is selected. */
      messageText = getFeedbackText(detailView, NULL, numSeqsSelected);
    }
  
  GtkWidget *feedbackBox = detailViewGetFeedbackBox(detailView);
  gtk_entry_set_text(GTK_ENTRY(feedbackBox), messageText);
  gtk_widget_queue_draw(feedbackBox);
  
  g_free(messageText);
}


/* Clear the feedback area */
void clearFeedbackArea(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  if (properties)
    {
      GtkStatusbar *statusBar = GTK_STATUSBAR(properties->statusBar);
      
      if (statusBar)
        {
          guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
          gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, "");
        }
    }
}


/* This updates the detail-view feedback area with info about a given nucleotide in the reference sequence */
void updateFeedbackAreaNucleotide(GtkWidget *detailView, const int dnaIdx, const BlxStrand strand)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);

  /* First clear the existing message if there is one */
  GtkStatusbar *statusBar = GTK_STATUSBAR(properties->statusBar);
  
  if (statusBar)
    {
      guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
      gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
      
      /* See if there's a variation at the given coord */
      const MSP *msp = NULL;
      
      gboolean multiple = FALSE;
      
      if (coordAffectedByVariation(dnaIdx, strand, bc, &msp, NULL, NULL, NULL, NULL, &multiple))
        {
          if (msp && mspGetSName(msp))
            {
              char *displayText = NULL;
              
              /* If we're displaying coords negated, negate it now */
              int coord = (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS] ? -1 * dnaIdx : dnaIdx);

              /* If there are multiple variations on this coord, display some summary text.
               * Otherwise, check we've got the sequence info to display. We should have, but 
               * if not just display the name. */
              if (multiple)
                displayText = g_strdup_printf("%d  %s", coord, MULTIPLE_VARIATIONS_TEXT);
              else if (mspGetMatchSeq(msp))
                displayText = g_strdup_printf("%d  %s : %s", coord, mspGetSName(msp), mspGetMatchSeq(msp));
              else
                displayText = g_strdup_printf("%d  %s", coord, mspGetSName(msp));
              
              /* Send the message to the status bar */
              gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
              g_free(displayText);
            }
        }
      else if (coordAffectedByPolyASignal(dnaIdx, strand, bc, &msp, &multiple) && msp)
        {
          char *displayText = NULL;

          /* If we're displaying coords negated, negate it now */
          const int negate = (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS] ? -1 : 1);
          int coord = negate * dnaIdx;

          /* If there are multiple signals on this coord, display some summary text.
           * Otherwise, check we've got the sequence info to display. If not just display the name. */
          if (multiple)
            {
              displayText = g_strdup_printf("%d  %s", coord, MULTIPLE_POLYA_SIGNALS_TEXT);
            }
          else
            {
              if (bc->displayRev) /* swap coords and negage if applicable */
                displayText = g_strdup_printf("%d %s: %d,%d", coord, "polyA signal", negate * msp->qRange.max, negate * msp->qRange.min);
              else
                displayText = g_strdup_printf("%d %s: %d,%d", coord, "polyA signal", msp->qRange.min, msp->qRange.max);
            }
              
          /* Send the message to the status bar */
          gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
          g_free(displayText);
        }
      else if (coordAffectedByPolyASite(dnaIdx, strand, bc, &msp, &multiple) && msp)
        {
          char *displayText = NULL;

          /* If we're displaying coords negated, negate it now */
          const int negate = (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS] ? -1 : 1);
          int coord = negate * dnaIdx;

          /* If there are multiple sites on this coord, display some summary text.
           * Otherwise, check we've got the sequence info to display. If not just display the name. */
          if (multiple)
            {
              displayText = g_strdup_printf("%d  %s", coord, MULTIPLE_POLYA_SITES_TEXT);
            }
          else
            {
              if (bc->displayRev) /* swap coords and negate if applicable */
                displayText = g_strdup_printf("%d %s: %d,%d", coord, "polyA site", negate * msp->qRange.max, negate * msp->qRange.min);
              else
                displayText = g_strdup_printf("%d %s: %d,%d", coord, "polyA site", msp->qRange.min, msp->qRange.max);
            }
              
          /* Send the message to the status bar */
          gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
          g_free(displayText);
        }
    }
}


/* Scroll the detail-view if necessary to keep it within the given range */
void detailViewScrollToKeepInRange(GtkWidget *detailView, const IntRange* const range)
{
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const BlxSeqType seqType = detailViewGetSeqType(detailView);
      
  if (displayRange->min < range->min)
    {
      setDetailViewStartIdx(detailView, range->min, seqType);
    }
  else if (displayRange->max > range->max)
    {
      setDetailViewEndIdx(detailView, range->max, seqType);
    }
}


/* If the selected base index is outside the current display range, scroll to
 * keep it in range. We scroll by the minimum number of bases possible if 
 * scrollMinimum is true; otherwise we re-centre on the selection. */
static void scrollToKeepSelectionInRange(GtkWidget *detailView, const gboolean scrollMinimum)
{
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);

  if (selectedBaseIdx != UNSET_INT && !valueWithinRange(selectedBaseIdx, displayRange))
    {
      const BlxSeqType seqType = detailViewGetSeqType(detailView);
      
      if (scrollMinimum)
        {
          if (selectedBaseIdx < displayRange->min)
            {
              setDetailViewStartIdx(detailView, selectedBaseIdx, seqType);
            }
          else if (selectedBaseIdx > displayRange->max)
            {
              setDetailViewEndIdx(detailView, selectedBaseIdx, seqType);
            }
        }
      else
        {
          const int newStart = selectedBaseIdx - (getRangeLength(displayRange) / 2);
          setDetailViewStartIdx(detailView, newStart, seqType);
        }
    }
}


/* This is called when the Squash matches option has been toggled. It switches to the condensed 
 * view of the trees (where multiple MSPs in the same sequence appear on the same row) if squash
 * is true, or reverts to the expanded version (where each MSP has its own row) if squash is false. */
void detailViewUpdateSquashMatches(GtkWidget *detailView, const gboolean squash)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  if (squash && bc->modelId != BLXMODEL_SQUASHED)
    {
      /* Set the "squashed" model to be active */
      bc->modelId = BLXMODEL_SQUASHED;
      callFuncOnAllDetailViewTrees(detailView, treeUpdateSquashMatches, NULL);
    }
  else if (!squash && bc->modelId != BLXMODEL_NORMAL)
    {
      /* Set the normal model to be active */
      bc->modelId = BLXMODEL_NORMAL;
      callFuncOnAllDetailViewTrees(detailView, treeUpdateSquashMatches, NULL);
    }

  detailViewRedrawAll(detailView);
}



/* Sort comparison function for sorting by group */
static gint sortByGroupCompareFunc(const MSP *msp1, const MSP *msp2, GtkWidget *blxWindow)
{
  gint result = 0;
  
  /* Get the order number out of the group and sort on that. If the sequence
   * is not in a group, its order number is UNSET_INT, and it gets sorted after
   * any sequences that are in groups. */
  const int msp1Order = sequenceGetGroupOrder(blxWindow, msp1->sSequence);
  const int msp2Order = sequenceGetGroupOrder(blxWindow, msp2->sSequence);
  
  if (msp1Order == UNSET_INT && msp2Order != UNSET_INT)
    {
      result = 1;
    }
  else if (msp1Order != UNSET_INT && msp2Order == UNSET_INT)
    {
      result = -1;
    }
  else
    {
      result = msp1Order - msp2Order;
    }

  return result;
}

/* Sort comparison function for sorting by the start position on the 
 * reference sequence. (Does a secondary sort by the alignment length) */
static gint sortByStartCompareFunc(const MSP *msp1, const MSP *msp2, const gboolean displayRev)
{
  gint result = 0;

  if (displayRev)
    {
      /* Display is reversed (i.e. numbers shown descending) so use compare on the max coord
       * and look for the max */
      result = msp2->qRange.max - msp1->qRange.max;
    }
  else 
    {
      result = msp1->qRange.min - msp2->qRange.min;
    }
  
  if (result == 0)
    {
      /* If the MSPs have the same start coord, do a secondary sort 
       * by alignment length */
      result = getRangeLength(&msp1->qRange) - getRangeLength(&msp2->qRange);
    }

  return result;
}

/* Sort comparison function for sorting by doubles */
static gint sortByDoubleCompareFunc(const MSP *msp1, const MSP *msp2)
{
  gint result = 0;
  
  gdouble dResult = msp1->id - msp2->id;
  
  if (dResult == 0)
    result = 0;
  else if (dResult > 0)
    result = 1;
  else 
    result = -1;
  
  return result;
}

/* Sort comparison function for sorting by the start position on the reference sequence
 * when we have multiple MSPs to compare (does a secondary sort by alignment length) */
static gint sortByStartCompareFuncMultiple(GList *mspList1, GList *mspList2, const gboolean msp1Fwd, const gboolean msp2Fwd, const gboolean displayRev)
{
  gint result = 0;

  const int coord1 = findMspListQExtent(mspList1, !displayRev, BLXSTRAND_NONE); /* find min coord unless display rev */
  const int coord2 = findMspListQExtent(mspList2, !displayRev, BLXSTRAND_NONE);

  result = coord1 - coord2;

  if (displayRev)
    {
      /* Display is reversed (i.e. numbers shown descending) so reverse the result */
      result *= -1;
    }

  return result;
}

/* Sort comparison function for sorting by string values. Allows NULL values and
 * sorts them AFTER non-null values. Comparison is case-insensitive. */
static int sortByStringCompareFunc(const char *str1, const char *str2)
{
  int result = 0;

  if (!str1 && !str2)
    {
      result = 0;
    }
  else if (!str1)
    {
      result = 1;
    }
  else if (!str2)
    {
      result = -1;
    }
  else
    {
      const int len = min(strlen(str1), strlen(str2));
      result = g_ascii_strncasecmp(str1, str2, len);
    }
  
  return result;
}


/* Determine the sort order for this column. Columns may have a different default
 * sort order (e.g. name is sorted ascending, but score is sorted descending.)
 * However, the opposite sort order is returned if the detail view's invert-sort-
 * order flag is set. */
static GtkSortType getColumnSortOrder(BlxViewContext *bc, const BlxColumnId columnId)
{
  GtkSortType result = GTK_SORT_ASCENDING;

  switch (columnId)
    {
      case BLXCOL_SCORE: /* fall through */
      case BLXCOL_ID:
        result = GTK_SORT_DESCENDING;
        break;
        
      default: /* all others ascending */
        break;
    };

  if (bc->flags[BLXFLAG_INVERT_SORT])
    {
      result = (result == GTK_SORT_ASCENDING) ? GTK_SORT_DESCENDING : GTK_SORT_ASCENDING;
    }

  return result;
}


  /* We're only interested in sorting exons and matches */
static gboolean mspIsSortable(const MSP* const msp)
{
  return (msp && (typeIsMatch(msp->type) || mspIsExon(msp)));
}


/* Sort comparison function for sorting by a particular column of the tree view. */
gint sortByColumnCompareFunc(GList *mspGList1,
                             GList *mspGList2,
                             GtkWidget *detailView, 
                             const BlxColumnId sortColumn)
{
  gint result = 0;
  
  /* Get the first MSP in each list. */
  MSP *msp1 = mspGList1 ? (MSP*)(mspGList1->data) : NULL;
  MSP *msp2 = mspGList2 ? (MSP*)(mspGList2->data) : NULL;

  /* If an msp is of a type that we don't bother sorting, place it before any that we do sort */
  if (!mspIsSortable(msp1) && !mspIsSortable(msp2))
    return 0;
  else if (!mspIsSortable(msp1))
    return -1;
  else if (!mspIsSortable(msp2))
    return 1;

  /* Check whether either row has more than one MSP. If so, it means some options
   * aren't applicable (\to do: unless they're short reads, which should be identical if
   * they're in the same row, so we can treat those as singular). */
  const gboolean multipleMsps = 
    (g_list_length(mspGList1) > 1 || g_list_length(mspGList2) > 1);
  
  BlxViewContext *bc = detailViewGetContext(detailView);
  gboolean displayRev = bc->displayRev;

  switch (sortColumn)
  {
    case BLXCOL_NONE:
      result = 0;
      break;
      
    case BLXCOL_SCORE:
      result = multipleMsps ? 0 : (int)(msp1->score - msp2->score);
      break;
      
    case BLXCOL_ID:
      result = multipleMsps ? 0 : sortByDoubleCompareFunc(msp1, msp2);
      break;
      
    case BLXCOL_START:
      if (multipleMsps)
        result = sortByStartCompareFuncMultiple(mspGList1, mspGList2, msp1->qStrand == BLXSTRAND_FORWARD, msp2->qStrand == BLXSTRAND_FORWARD, displayRev);
      else
        result = sortByStartCompareFunc(msp1, msp2, displayRev);
      break;
      
    case BLXCOL_GROUP:
      result = sortByGroupCompareFunc(msp1, msp2, detailViewGetBlxWindow(detailView));
      break;
      
    default:
      /* Generic string column */
      result = sortByStringCompareFunc(mspGetColumn(msp1, sortColumn), mspGetColumn(msp2, sortColumn));
      break;
  };

  /* Invert the result if we need to sort descending instead of ascending. */
  if (getColumnSortOrder(bc, sortColumn) == GTK_SORT_DESCENDING)
    {
      result *= -1;
    }
  
  return result;
}


/* This is the main sort comparison function for comparing two BlxSequences.
 * The sort criteria are specified in the detailView properties;
 * we may sort by multiple columns.
 * 
 * Returns a negative value if the first row appears before the second,
 * positive if the second appears before the first, or 0 if they are 
 * equivalent. */
static gint detailViewSortByColumns(gconstpointer a, gconstpointer b)
{
  gint result = UNSET_INT;

  const BlxSequence *seq1 = (const BlxSequence*)a;
  const BlxSequence *seq2 = (const BlxSequence*)b;

  GtkWidget *blxWindow = getBlixemWindow();
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GList *columnList = detailViewGetColumnList(detailView);
  const int numColumns = g_list_length(columnList);

  /* Sort by each requested column in order of priority */
  BlxColumnId *sortColumns = detailViewGetSortColumns(detailView);
  
  if (sortColumns)
    {
      int priority = 0;

      for ( ; priority < numColumns; ++priority)
        {
          BlxColumnId sortColumn = sortColumns[priority];
        
          /* NONE indicates an unused entry in the priority array; if we reach
           * an unset value, there should be no more values after it */
          if (sortColumn == BLXCOL_NONE)
            break;

          /* Extract the MSPs for the two rows that we're comparing */
          GList *mspGList1 = seq1->mspList;
          GList *mspGList2 = seq2->mspList;
  
          /* Do the comparison on this column */
          result = sortByColumnCompareFunc(mspGList1, mspGList2, detailView, sortColumn);
          
          /* If rows are equal, continue to sort; otherwise we're done */
          if (result != 0)
            break;
        }
    }
  
  return result;
}



/* Re-sort all trees */
void detailViewResortTrees(GtkWidget *detailView)
{
  /* Sort the data for each tree */
  callFuncOnAllDetailViewTrees(detailView, resortTree, NULL);

  /* Sort the list of BlxSequences (uesd by the exon view) */
  BlxViewContext *bc = detailViewGetContext(detailView);
  bc->matchSeqs = g_list_sort(bc->matchSeqs, detailViewSortByColumns);

  bigPictureRedrawAll(detailViewGetBigPicture(detailView));
}


/* Set the value of the 'invert sort order' flag */
void detailViewUpdateSortInverted(GtkWidget *detailView, const gboolean invert)
{
  detailViewResortTrees(detailView);
}


/* Set the value of the 'Show SNP track' flag */
void detailViewUpdateShowSnpTrack(GtkWidget *detailView, const gboolean showSnpTrack)
{
  refreshDetailViewHeaders(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  detailViewRedrawAll(detailView);
  
  refreshDialog(BLXDIALOG_SETTINGS, detailViewGetBlxWindow(detailView));
}


/* This performs required updates required after editing anything that affects
 * the display length of MSPs in the detail-view, e.g. the 'show unaligned sequence'
 * or 'show polyA tails' options. */
void detailViewUpdateMspLengths(GtkWidget *detailView, const int numUnalignedBases)
{
  /* Re-calculate the full extent of all MSPs, and the max msp length */
  setMaxMspLen(0);
  
  BlxViewContext *bc = detailViewGetContext(detailView);
  MSP *msp = bc->mspList;
  
  for ( ; msp; msp = msp->next)
    {
      mspCalculateFullExtents(msp, bc, numUnalignedBases);
    }
  
  /* Do a full re-sort and re-filter because the lengths of the displayed match
   * sequences may have changed (and we need to make sure they're sorted by start pos) */
  detailViewResortTrees(detailView);
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  detailViewRedrawAll(detailView);
}


/* Set the number of additional bases to show for the show-unaligned-sequence option */
void detailViewSetNumUnalignedBases(GtkWidget *detailView, const int numBases)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->numUnalignedBases = numBases;
  detailViewUpdateMspLengths(detailView, properties->numUnalignedBases);
}


int getBaseIndexAtColCoords(const int x, const int y, const gdouble charWidth, const IntRange* const displayRange)
{
  int result = UNSET_INT;
  
  const int leftEdge = 0; /* to do: allow for padding? check upper bound? */
  
  if (x > leftEdge)
    {
      result = (int)(((gdouble)x - leftEdge) / charWidth);
    }
  
  result += displayRange->min;
  
  return result;
}


/* Get the nucleotide at the given x/y position in the given sequence-column header (protein
 * matches only). */
static void getSeqColHeaderClickedNucleotide(GtkWidget *header, GtkWidget *detailView, const int x, const int y, int *coordOut, int *frameOut, int *baseOut)
{
  int coord = getBaseIndexAtColCoords(x, y, detailViewGetCharWidth(detailView), detailViewGetDisplayRange(detailView));
  
  if (coord != UNSET_INT)
    {   
      /* Get the base number of the clicked base within the active frame. The
       * base number is determined by which row in the header the mouse pointer 
       * is over: for frame 1, row 1 will give base 1; for frame 2, row 2 will
       * give base 1, etc. Start by getting the frame number for the clicked row: */
      int frame = detailViewGetActiveFrame(detailView);
      int row = seqColHeaderGetRow(header);
      const int numFrames = detailViewGetNumFrames(detailView);
      
      /* The header widget passed to this function is the originally-clicked widget.
       * If the pointer has dragged onto another row in the header, we can work out
       * which one based on the y coord. If it is outside altogether, take the 
       * topmost/bottommost row. (NB this assumes all rows in the header have the 
       * same height as this one - should be true because they're just simple labels). */
      if (y < 0 || y > header->allocation.height)
        {
          const int offsetRows = floor((double)y / (double)header->allocation.height); /* can be negative */
          row += offsetRows;
          
          if (row < 1)
            {
              row = 1;
            }
          else if (row > numFrames)
            {
              row = numFrames;
            }
        }
      
      /* Calculate the base number based on the row and the currently-active frame */
      int baseNum = row - frame + 1;
      
      if (baseNum < 1)
        {
          /* Cycle round if gone below the min base number */
          baseNum += numFrames;
          --coord;
        }
      
      if (coordOut)
        *coordOut = coord;
      
      if (frameOut)
        *frameOut = frame;
      
      if (baseOut)
        *baseOut = baseNum;
    }
}


/* Select the coord at the given x,y position in the given sequence column header */
static void selectClickedNucleotide(GtkWidget *header, GtkWidget *detailView, const int x, const int y)
{
  int coord, frame, baseNum;
  getSeqColHeaderClickedNucleotide(header, detailView, x, y, &coord, &frame, &baseNum);
  
  detailViewSetSelectedBaseIdx(detailView, coord, frame, baseNum, FALSE, TRUE);
}


/* If the user clicked on a SNP in the SNP header, select it. This updates the 
 * feedback box with the SNP name and coord, selects the SNP coord and recentres 
 * on it (if allowScroll is true), clearing any previous selections. If expandSnps
 * is true it means that the SNPs are drawn full-width (i.e. showing all the sequence
 * in the MSP) and we should take into account when clicking; otherwise, the SNP
 * is just shown as a marker at its actual coord, and we only need to check the
 * position of that coord to see if the user clicked there. snpTrack is the widget
 * that was clicked on and the optional 'colHeader' argument can be passed if coordinates
 * should be converted to the coordinate system of that column header (i.e. if the SNP track
 * widget is not the same width of the sequence column, the coords will need converting). */
void selectClickedSnp(GtkWidget *snpTrack,
                      GtkWidget *colHeader,
                      GtkWidget *detailView, 
                      const int xIn, 
                      const int yIn, 
                      const gboolean expandSnps,
                      const int clickedBase)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  
  /* Convert x coord to sequence-column coords */
  int x = xIn, y = yIn;
  
  if (colHeader)
    {
      gtk_widget_translate_coordinates(snpTrack, colHeader, xIn, 0, &x, NULL);
    }

  int clickedDisplayIdx = getBaseIndexAtColCoords(x, y, properties->charWidth, &properties->displayRange);

  /* Get the clicked row index (0-based) */
  const int rowHeight = ceil(properties->charHeight);
  const int clickedRow = y / rowHeight;
  const int numRows = getNumSnpTrackRows(bc, properties, snpTrackGetStrand(snpTrack, detailView), activeFrame);

  if (clickedDisplayIdx != UNSET_INT)
    {
      /* Loop through all variations and see if there are any at this displayIdx and row */
      GList *snpList = NULL;
      GSList *rows = NULL;
      int i = 0;
      const MSP *msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], i);
      
      for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], ++i))
        {
          /* Get the variation coords in terms of display coords, and also the 'expanded' range
           * of coords where the variation is displayed */
          IntRange mspExpandedRange;
          getVariationDisplayRange(msp, expandSnps, bc->seqType, bc->numFrames, bc->displayRev, activeFrame, &bc->refSeqRange, NULL, &mspExpandedRange);
          
          const int clickedDnaIdx = convertDisplayIdxToDnaIdx(clickedDisplayIdx, bc->seqType, 1, clickedBase, bc->numFrames, bc->displayRev, &bc->refSeqRange);
          gboolean found = FALSE;
          int dnaIdxToSelect = UNSET_INT;
          
          /* Get the row this variation is in */
          int mspRow = 1;
          if (expandSnps)
            mspRow = getVariationRowNumber(&mspExpandedRange, numRows, &rows) - 1;

          if (expandSnps && valueWithinRange(clickedDisplayIdx, &mspExpandedRange) && mspRow == clickedRow)
            {
              /* We clicked inside this MSP on the variation track. Select the first coord in
               * the MSP. */
              found = TRUE;
              dnaIdxToSelect = msp->qRange.min;
            }
          else if (!expandSnps && valueWithinRange(clickedDnaIdx, &msp->qRange))
            {
              /* We clicked on a coord in the ref seq header that is affected by this variation.
               * Select the clicked coord. */
              found = TRUE;
              dnaIdxToSelect = clickedDnaIdx;
            }
          
          if (found)
            {
              snpList = g_list_prepend(snpList, msp->sSequence);
              detailViewSetSelectedDnaBaseIdx(detailView, dnaIdxToSelect, activeFrame, TRUE, FALSE);
            }
        }
      
      freeRowsList(rows);

      /* Clear any existing selections and select the new SNP(s) */
      blxWindowSetSelectedSeqList(blxWindow, snpList);
    }
}


/***********************************************************
 *                          Drawing                        *
 ***********************************************************/

/* Draw a vertical separator line at the edges of the given widget */
void drawColumnSeparatorLine(GtkWidget *widget, 
                             GdkDrawable *drawable, 
                             GdkGC *gc, 
                             const BlxViewContext *bc)
{
  const int lineWidth = 1; /* width of the separator lines */
  
  if (GTK_WIDGET_VISIBLE(widget) && widget->allocation.width > lineWidth)
    {
      GdkColor *color = getGdkColor(BLXCOLOR_TREE_GRID_LINES, bc->defaultColors, FALSE, bc->usePrintColors);
      gdk_gc_set_foreground(gc, color);
      
      /* We want dashed lines */
      gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      /* Make the dashes very short and closely packed (i.e. dash length of 1 pixel and gaps of 1 pixel) */
      int listLen = 1;
      gint8 dashList[listLen];
      dashList[0] = 1;
      gdk_gc_set_dashes(gc, 0, dashList, listLen);
      
      /* Draw vertical lines. Draw the rightmost edge of each widget (because we don't really want a
       * line at the leftmost edge of the first cell of the tree view.) */
      const int x = widget->allocation.x + widget->allocation.width - lineWidth;
      const int y1 = widget->allocation.y;
      const int y2 = y1 + widget->allocation.height;
      
      gdk_draw_line(drawable, gc, x, y1, x, y2);
    }
}


/* Populates the 'bases' string with the two bases from the ref seq adjacent to the start/end
 * of the given MSP (if it is an exon/match). 'Start' here means at the lower value coord end
 * of the MSP. 'revStrand' is true if this MSP is on the reverse strand of the ref seq */
static void mspGetAdjacentBases(const MSP* const msp, char *bases, const gboolean start, const gboolean revStrand, const BlxViewContext *bc)
{
  if (start)
    {
      bases[0] = getSequenceIndex(bc->refSeq, msp->qRange.min - 2, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[1] = getSequenceIndex(bc->refSeq, msp->qRange.min - 1, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[2] = '\0';
    }  
  else
    {
      bases[0] = getSequenceIndex(bc->refSeq, msp->qRange.max + 1, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[1] = getSequenceIndex(bc->refSeq, msp->qRange.max + 2, revStrand, &bc->refSeqRange, BLXSEQ_DNA);
      bases[2] = '\0';
    }
}


/* This function returns the prev/next MSP in the given sequence in order of start position on
 * the ref seq. Returns the next MSP (with the next-highest start coord) if 'getNext' is true; 
 * otherwise returns the previous one. The return value will be NULL if there is no next/prev MSP.
 * The error will be set if there was any problem (e.g. if the given MSP does not exist the the sequence).
 * Only considers exons and matches. */
static const MSP* sequenceGetNextMsp(const MSP* const msp, 
                                     const BlxSequence *blxSeq, 
                                     const gboolean getNext, 
                                     GError **error)
{
  const MSP *result = NULL;

  GList *mspItem = blxSeq->mspList;
  const MSP *prevMsp = NULL;
  gboolean found = FALSE;
  
  /* Loop through all MSPs and look for the given one. */
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP* const curMsp = (const MSP*)(mspItem->data);
      
      if (curMsp == msp)
        {
          found = TRUE;
          
          /* If looking for the previous MSP, we already know it. */
          if (!getNext)
            {
              result = prevMsp;
              break;
            }
          else
            {
              /* Continue looping to find the next MSP. */
              continue;
            }
        }
      
      if (found && getNext)
        {
          /* We've found the input MSP and we're on the next one. See if it is an exon/match,
           * otherwise continue looping. */
          if (mspIsExon(curMsp) || mspIsBlastMatch(msp))
            {
              result = curMsp;
              break;
            }
        }
      
      /* Set the previous MSP (but only if it's an exon/match) */
      if (mspIsExon(curMsp) || mspIsBlastMatch(curMsp))
        {
          prevMsp = curMsp;
        }
    }
  
  if (!found && error)
    {
      g_set_error(error, BLX_ERROR, 1, "The given MSP '%s' was not found in the given sequence '%s'.\n", mspGetSName(msp), blxSequenceGetName(blxSeq));
    }
  
  return result;
}


/* Determine whether the bases at the other end of the intron for the start/end of 
 * the given MSP are canonical. 'canonicalStart' and 'canonicalEnd' contain the two bases to look
 * for at the start/end of the next/previous MSP to determine if the result is canonical. */
static gboolean mspIsSpliceSiteCanonicalOtherEnd(const MSP* const msp,
                                                 const BlxSequence *blxSeq,
                                                 const gboolean isMinCoord, 
                                                 const gboolean revStrand,
                                                 const gboolean donor,
                                                 const BlxViewContext *bc,
                                                 const BlxSpliceSite *spliceSite)
{
  gboolean result = FALSE;
 
  /* Get the previous MSP if we're looking at the min coord of the current MSP, or 
   * the next MSP if we're at the end */
  const MSP *nextMsp = sequenceGetNextMsp(msp, blxSeq, !isMinCoord, NULL);
  
  if (nextMsp)
    {
      /* Get the two bases at the end of the previous MSP / start of the next MSP */
      char bases[3];
      mspGetAdjacentBases(nextMsp, bases, !isMinCoord, revStrand, bc);
      
      /* If the original site was a donor site, look for an acceptor site (or vice versa) */
      const char *canonicalBases = spliceSiteGetBases(spliceSite, !donor, revStrand);

      if (stringsEqual(bases, canonicalBases, FALSE))
        {
          result = TRUE;
        }
    }
  
  return result;
}


/* Create a BlxSpliceSite and add it to the given list */
static void addBlxSpliceSite(GSList **spliceSites, const char *donorSite, const char *acceptorSite, const gboolean bothReqd)
{
  BlxSpliceSite *spliceSite = (BlxSpliceSite*)g_malloc(sizeof(BlxSpliceSite));

  if (strlen(donorSite) < 2 || strlen(acceptorSite) < 2)
    {
      g_critical("Error adding splice site info ['%s', '%s'].\n", donorSite, acceptorSite);
      return;
    }
  
  spliceSite->donorSite[0] = donorSite[0];
  spliceSite->donorSite[1] = donorSite[1];
  spliceSite->donorSite[2] = '\0';

  spliceSite->donorSiteRev[0] = donorSite[1];
  spliceSite->donorSiteRev[1] = donorSite[0];
  spliceSite->donorSiteRev[2] = '\0';

  spliceSite->acceptorSite[0] = acceptorSite[0];
  spliceSite->acceptorSite[1] = acceptorSite[1];
  spliceSite->acceptorSite[2] = '\0';
  
  spliceSite->acceptorSiteRev[0] = acceptorSite[1];
  spliceSite->acceptorSiteRev[1] = acceptorSite[0];
  spliceSite->acceptorSiteRev[2] = '\0';
  
  spliceSite->bothReqd = bothReqd;
  
  /* Add it to the list */
  *spliceSites = g_slist_append(*spliceSites, spliceSite);
}


/* Frees the memory used by a BlxSpliceSite */
static void destroyBlxSpliceSite(gpointer listItemData, gpointer data)
{
  BlxSpliceSite *spliceSite = (BlxSpliceSite*)listItemData;
  
  g_free(spliceSite);
}


/* Return the canonical bases at the donor/acceptor end of the given BlxSpliceSite. Returns them in
 * reverse order if 'reverse' is true, e.g. for a GC-AG intron, the dono is GC and acceptor is AG.
 * If reverse is true the donor returns CG and acceptor returns GA. The result is a pointer to the
 * string in the splice site, which is owned by the Detail View properties. */
static const char* spliceSiteGetBases(const BlxSpliceSite *spliceSite, const gboolean donor, const gboolean reverse)
{
  const char *result = NULL;
  
  if (donor)
    {
      result = reverse ? spliceSite->donorSiteRev : spliceSite->donorSite;
    }
  else
    {
      result = reverse ? spliceSite->acceptorSiteRev : spliceSite->acceptorSite;
    }
  
  return result;
}


/* Determines whether the intron splice sites for the given MSP are canonical or not.
 * We look for GC-AG or AT-AC introns. The latter must have both ends of the intron matching
 * to be canonical. GC is the donor site and AG is the acceptor site. */
static gboolean isMspSpliceSiteCanonical(const MSP* const msp, 
                                         const BlxSequence *blxSeq, 
                                         const gboolean isMinCoord, 
                                         const gboolean revStrand,
                                         const BlxViewContext *bc,
                                         GSList *spliceSites)
{
  gboolean result = FALSE;

  GSList *item = spliceSites;
  for ( ; item; item = item->next)
    {
      const BlxSpliceSite *spliceSite = (const BlxSpliceSite*)(item->data);

      /* Get the two adjacent bases at the start/end of the given MSP */
      char bases[3];
      mspGetAdjacentBases(msp, bases, isMinCoord, revStrand, bc);
      
      /* The min coord end of the MSP is the acceptor site of the intron and the max coord the 
       * donor site (or vice versa if the strand is reversed). */
      const gboolean donor = (isMinCoord == revStrand);
      const char *canonicalBases = spliceSiteGetBases(spliceSite, donor, revStrand);

      if (stringsEqual(bases, canonicalBases, FALSE))
        {
          if (spliceSite->bothReqd)
            {
              /* Must also check if the other end of the intron matches */
              result = mspIsSpliceSiteCanonicalOtherEnd(msp, blxSeq, isMinCoord, revStrand, donor, bc, spliceSite);
            }
          else
            {
              result = TRUE;
            }
        }
    }

  return result;
}


/* If the given MSP is an exon/match, get the 2 nucleotides from the reference sequence at the
 * splice sites of the adjacent introns. Inserts the results into the given hash table with
 * an enum indicating whether the nucleotides are canonical/non-canonical. Only considers
 * MSPs that are within the given ref seq range. */
static void mspGetSpliceSiteCoords(const MSP* const msp, 
                                   const BlxSequence *blxSeq, 
                                   const IntRange* const qRange, 
                                   const BlxViewContext *bc, 
                                   GSList *spliceSites,
                                   GHashTable *result)
{
  const gboolean revStrand = (mspGetRefStrand(msp) != BLXSTRAND_FORWARD);
  
  /* Ignore the termini (i.e. 5' end of first exon and 3' end of last exon) */
  const gboolean getMin = (msp->sSequence && msp->qRange.min != blxSequenceGetStart(msp->sSequence, msp->qStrand));
  const gboolean getMax = (msp->sSequence && msp->qRange.max != blxSequenceGetEnd(msp->sSequence, msp->qStrand));
  
  /* See if the min coord is within the given range */
  if (getMin && valueWithinRange(msp->qRange.min, qRange) && msp->qRange.min >= bc->refSeqRange.min + 2)
    {
      /* Find out if the adjacent bases are canonical/non-canonical. */
      const gboolean canonical = isMspSpliceSiteCanonical(msp, blxSeq, TRUE, revStrand, bc, spliceSites);
      BlxColorId colorId = canonical ? BLXCOLOR_CANONICAL : BLXCOLOR_NON_CANONICAL;

      /* Insert the two coords into the hash table */
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.min - 2), GINT_TO_POINTER(colorId));
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.min - 1), GINT_TO_POINTER(colorId));
    }
  
  /* See if the max coord is within the given range */
  if (getMax && valueWithinRange(msp->qRange.max, qRange) && msp->qRange.max <= bc->refSeqRange.max - 2)
    {
      /* Find out if they are canonical/non-canonical */
      const gboolean canonical = isMspSpliceSiteCanonical(msp, blxSeq, FALSE, revStrand, bc, spliceSites);
      BlxColorId colorId = canonical ? BLXCOLOR_CANONICAL : BLXCOLOR_NON_CANONICAL;
      
      /* Insert the two coords into the hash table */
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.max + 1), GINT_TO_POINTER(colorId));
      g_hash_table_insert(result, GINT_TO_POINTER(msp->qRange.max + 2), GINT_TO_POINTER(colorId));
    }
}


/* This looks through all of the polyA signals and sees if any of them are in the given display
 * range (in nucleotide coords) and whether we want to display them. If the 'display polyA signals for
 * selected sequences only' option is enabled, then any selected sequence MSPs that have polyA tails
 * should be passed in the given GSList. Any coords within a relevant polyA signal range are added to
 * the given hash table with the BlxColorId that they should be drawn with. */
static void getAnnotatedPolyASignalBasesToHighlight(const BlxViewContext *bc,
                                                    GSList *polyATailMsps,
                                                    const BlxStrand qStrand,
                                                    const IntRange* const qRange,
                                                    GHashTable *result)
{
  /* We've more work to do if showing polyA signals (unless we're only showing them for selected
   * sequences and there were no relevant selected sequences) */
  if (bc->flags[BLXFLAG_SHOW_POLYA_SIG] && (!bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED] || g_slist_length(polyATailMsps) > 0))
    {
      BlxColorId colorId = BLXCOLOR_POLYA_SIGNAL_ANN;

      /* Loop through all polyA signals */
      int i = 0;
      const MSP *sigMsp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SIGNAL], i);
      
      for ( ; sigMsp; sigMsp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SIGNAL], ++i))
        {
          /* Only interested the polyA signal has the correct strand and is within the display range. */
          if (sigMsp->qStrand == qStrand && rangesOverlap(&sigMsp->qRange, qRange))
            {
              gboolean addSignal = FALSE;
              
              if (bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED])
                {
                  /* We're only interested in polyA signals that are 50 bases upstream of one of our polyA-tail MSPs */
                  GSList *item = polyATailMsps;
                  
                  for ( ; item && !addSignal; item = item->next)
                    {
                      const MSP* const tailMsp = (const MSP*)(item->data);
                      const int qEnd = tailMsp->qRange.max;
                      const int qStart = qEnd - POLYA_SIG_BASES_UPSTREAM;
                      
                      IntRange upstreamRange;
                      intrangeSetValues(&upstreamRange, qStart, qEnd); /* sorts out which is min and which is max */
                                        
                      addSignal = rangesOverlap(&sigMsp->qRange, &upstreamRange);
                    }
                }
              else
                {
                  /* Add all signals that are in the display range */
                  addSignal = TRUE;
                }
            
              if (addSignal)
                {
                  /* Add each base in the polyA signal range to the hash table. This may overwrite
                   * splice site bases that we previously found, which is fine (something has to take priority). */
                  int i = sigMsp->qRange.min;
                  for ( ; i <= sigMsp->qRange.max; ++i)
                    {
                      g_hash_table_insert(result, GINT_TO_POINTER(i), GINT_TO_POINTER(colorId));
                    }
                }
            }
        }
    }
}


/* Scan upstream from any polyA tails to see if there are any polyA signals in the reference
 * sequence. Actually, because the visible range in the detail-view is generally quite small and
 * we don't want to re-scan lots of times, lets just scan through the whole visible range once
 * and show all signals - but only if there is at least one visible polyA tail in the range. */
static void getPolyASignalBasesToHighlight(GtkWidget *detailView,
                                           const BlxViewContext *bc,
                                           GSList *polyATailMsps,
                                           const IntRange* const qRange,
                                           GHashTable *result)
{
  /* If only showing polyAs for selected sequences, check that there's a (selected) msp with a
     polyA tail in range (i.e. in the input list). Otherwise show all polyA signals. */
  if (bc->flags[BLXFLAG_SHOW_POLYA_SIG] && (!bc->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED] || g_slist_length(polyATailMsps) > 0))
    {
      BlxColorId colorId = BLXCOLOR_POLYA_SIGNAL;
                                        
      /* Loop through each base in the visible range */
      const char *seq = bc->refSeq;
      int idx = qRange->min - bc->refSeqRange.min; /* convert to 0-based */
      const int max_idx = qRange->max - bc->refSeqRange.min; /* convert to 0-based */
      const char *cp = seq + idx;
      const char *comparison = POLYA_SIGNAL;
      const int comparison_len = strlen(comparison);

      for ( ; idx < max_idx && cp && *cp; ++idx, ++cp)
        {
          if (!strncasecmp(cp, comparison, comparison_len))
            {
              /* Add each base in the polyA signal range to the hash table. This may overwrite
               * splice site bases that we previously found, which is fine (something has to take priority). */
              int i = idx + bc->refSeqRange.min; /* convert back to coords */
              const int max = i + comparison_len;

              for ( ; i < max; ++i)
                {
                  g_hash_table_insert(result, GINT_TO_POINTER(i), GINT_TO_POINTER(colorId));
                }
            }
        }
    }
}


/* This looks through all of the polyA sites and sees if any of them are in the given display
 * range (in nucleotide coords) and whether we want to display them. */
static void getAnnotatedPolyASiteBasesToHighlight(const BlxViewContext *bc,
                                                  const BlxStrand qStrand,
                                                  const IntRange* const qRange,
                                                  GHashTable *result)
{
  /* We've more work to do if showing polyA signals (unless we're only showing them for selected
   * sequences and there were no relevant selected sequences) */
  if (bc->flags[BLXFLAG_SHOW_POLYA_SIG])
    {
      BlxColorId colorId = BLXCOLOR_POLYA_SITE_ANN;

      /* Loop through all polyA signals */
      int i = 0;
      const MSP *siteMsp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SITE], i);
      
      for ( ; siteMsp; siteMsp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SITE], ++i))
        {
          /* Only interested the polyA site has the correct strand and is within the display range. */
          if (siteMsp->qStrand == qStrand && rangesOverlap(&siteMsp->qRange, qRange))
            {
              /* Add just the first base in the polyA signal range to the hash table (because the
               * site is between the two bases in the range and we only want to draw it once). */
              g_hash_table_insert(result, GINT_TO_POINTER(siteMsp->qRange.min), GINT_TO_POINTER(colorId));
            }
        }
    }
}


/* This function looks for special bases in the reference sequence header to highlight and stores their
 * coords in the returned hash table with the BlxColorId (converted to a gpointer with GINT_TO_POINTER)
 * of the fill color they should be drawn with. These special coords include splice sites and polyA signals. */
GHashTable* getRefSeqBasesToHighlight(GtkWidget *detailView, 
                                      const IntRange* const qRange, 
                                      const BlxSeqType seqType,
                                      const BlxStrand qStrand)
{
  GHashTable *result = g_hash_table_new(g_direct_hash, g_direct_equal);

  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  const BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* We only highlight nucleotides, so there is nothing to do if we're showing a peptide sequence,
   * or if the show-splice-sites or show-polyA-signals options are disabled. */
  if (seqType == BLXSEQ_PEPTIDE || 
      (!bc->flags[BLXFLAG_SHOW_SPLICE_SITES] && !bc->flags[BLXFLAG_SHOW_POLYA_SIG]))
    {
      return result;
    }
  
  /* Loop through the selected sequences. If showing polyA signals for selected sequences, we'll 
   * compile a list of MSPs we're interested in seeing polyA signals for (i.e. those with polyA tails) */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GList *seqItem = blxWindowGetSelectedSeqs(blxWindow);
  GSList *polyATailMsps = NULL;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      const BlxSequence *blxSeq = (const BlxSequence*)(seqItem->data);
      GList *mspItem = blxSeq->mspList;
     
      /* Loop through all MSPs for this sequence */ 
      for ( ; mspItem; mspItem = mspItem->next)
        {
          MSP *msp = (MSP*)(mspItem->data);
          
          /* Only look at matches/exons on the correct strand */
          if ((mspIsBlastMatch(msp) || msp->type == BLXMSP_EXON) && mspGetRefStrand(msp) == qStrand)
            {
              if (bc->flags[BLXFLAG_SHOW_SPLICE_SITES])
                mspGetSpliceSiteCoords(msp, blxSeq, qRange, bc, properties->spliceSites, result);
              
              if (bc->flags[BLXFLAG_SHOW_POLYA_SIG] && mspHasPolyATail(msp))
                polyATailMsps = g_slist_append(polyATailMsps, msp);
            }
        }
    }
  
  /* Now check the polyA signals and sites and see if any of them are in range. */
  getPolyASignalBasesToHighlight(detailView, bc, polyATailMsps, qRange, result);
  getAnnotatedPolyASignalBasesToHighlight(bc, polyATailMsps, qStrand, qRange, result);
  getAnnotatedPolyASiteBasesToHighlight(bc, qStrand, qRange, result);

  return result;
}


static void drawDnaTrack(GtkWidget *dnaTrack, GtkWidget *detailView, const BlxStrand strand, const int frame)
{
  GdkDrawable *drawable = createBlankPixmap(dnaTrack);
  
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Find the segment of the ref sequence to display (complemented if this tree is
   * displaying the reverse strand, and reversed if the display is toggled). Ref seq
   * coords are in nucleotides, so convert display range to nucleotide coords. */
  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  
  const int qIdx1 = convertDisplayIdxToDnaIdx(displayRange->min, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);         /* 1st base in frame */
  const int qIdx2 = convertDisplayIdxToDnaIdx(displayRange->max, bc->seqType, frame, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange); /* last base in frame */
  IntRange qRange = {min(qIdx1, qIdx2), max(qIdx1, qIdx2)};
  
  GError *error = NULL;
  
  gchar *segmentToDisplay = getSequenceSegment(bc->refSeq,
                                               &qRange,
                                               strand, 
                                               BLXSEQ_DNA,      /* ref seq is always in nucleotide coords */
                                               BLXSEQ_DNA,      /* required segment is in nucleotide coords */
                                               frame, 
                                               bc->numFrames,
                                               &bc->refSeqRange,
                                               bc->blastMode,
                                               bc->geneticCode,
                                               bc->displayRev,
                                               bc->displayRev,
                                               bc->displayRev,
                                               &error);
  
  if (!segmentToDisplay)
    {
      g_assert(error);
      prefixError(error, "Could not draw DNA header. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
    }
    
    
  GdkGC *gc = gdk_gc_new(drawable);
  
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  const BlxStrand activeStrand = blxWindowGetActiveStrand(blxWindow);
  const gboolean highlightSnps = bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS];

  gtk_layout_set_size(GTK_LAYOUT(dnaTrack), dnaTrack->allocation.width, roundNearest(properties->charHeight));
  
  /* Find out if there are any bases in the introns that need highlighting. */
  GHashTable *basesToHighlight = getRefSeqBasesToHighlight(detailView, &qRange, BLXSEQ_DNA, activeStrand);
  
  /* We'll loop forward/backward through the display range depending on which strand we're viewing */
  int incrementValue = bc->displayRev ? -bc->numFrames : bc->numFrames;
  int displayLen = qRange.max - qRange.min + 1;
  char displayText[displayLen + 1];
  int displayTextPos = 0;
  
  int qIdx = bc->displayRev ? qRange.max : qRange.min;
  int displayIdx = convertDnaIdxToDisplayIdx(qIdx, bc->seqType, activeFrame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int y = 0;
  DrawBaseData baseData = {qIdx, 0, strand, UNSET_INT, BLXSEQ_DNA, TRUE, FALSE, FALSE, FALSE, highlightSnps, TRUE, BLXCOLOR_BACKGROUND, NULL, NULL, FALSE, FALSE, FALSE, FALSE};
  
  while (qIdx >= qRange.min && qIdx <= qRange.max)
    {
      /* Get the character to display at this index, and its position */
      displayText[displayTextPos] = getSequenceIndex(bc->refSeq, qIdx, bc->displayRev, &bc->refSeqRange, BLXSEQ_DNA);
      const int x = (int)((gdouble)displayTextPos * properties->charWidth);
      
      baseData.dnaIdx = qIdx;
      baseData.baseChar = displayText[displayTextPos];
      baseData.displayIdxSelected = (displayIdx == properties->selectedBaseIdx);
      baseData.dnaIdxSelected = (qIdx == properties->selectedDnaBaseIdx);

      /* Color the base depending on whether it is selected or affected by a SNP */
      drawHeaderChar(bc, properties, drawable, gc, x, y, basesToHighlight, &baseData);
      
      /* Increment indices */
      ++displayTextPos;
      ++displayIdx;
      qIdx += incrementValue;
    }
  
  /* Make sure the string is terminated properly */
  displayText[displayTextPos] = '\0';
  
  /* Draw the text */
  PangoLayout *layout = gtk_widget_create_pango_layout(detailView, displayText);
  pango_layout_set_font_description(layout, detailViewGetFontDesc(detailView));
  
  if (layout)
    {
      gtk_paint_layout(dnaTrack->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, 0, 0, layout);
      g_object_unref(layout);
    }
  
  drawColumnSeparatorLine(dnaTrack, drawable, gc, bc);
  
  g_hash_table_unref(basesToHighlight);
  g_free(segmentToDisplay);
  g_object_unref(gc);
}


/* Utility to find the display coord(s) that a variation lies on. If expand is true, 
 * it means that if the variation contains more than one alternative, and these will be 
 * displayed horizontally across the display, taking up more width than where 
 * the actual variation coords lie. expandedRange will take this into account. displayRange 
 * always returns the actual msp start/end in display coords. */
static void getVariationDisplayRange(const MSP *msp, 
                                     const gboolean expand,
                                     const BlxSeqType seqType, 
                                     const int numFrames,
                                     const gboolean displayRev, 
                                     const int activeFrame,
                                     const IntRange* const refSeqRange,
                                     IntRange *displayRange,
                                     IntRange *expandedRange)
{
  /* Convert the MSP coords to display coords. We want to display the variation
   * in the same position regardless of reading frame, so always use frame 1. */
  const IntRange* const mspRange = mspGetDisplayRange(msp);

  if (displayRange)
    intrangeSetValues(displayRange, mspRange->min, mspRange->max);
  
  if (!expandedRange)
    return;

  /* Work out the expanded range (it will be the same as the MSP range if unexpanded) */
  intrangeSetValues(expandedRange, mspRange->min, mspRange->max);
  
  if (expand && mspGetMatchSeq(msp))
    {
      /* Expand the variation range so that we display its entire sequence. We'll 
       * position the variation so that the middle of its sequence lies at the centre
       * coord of its ref seq range */
      const int numChars = strlen(mspGetMatchSeq(msp));

      int offset = (int)((double)numChars / 2.0);
  
      if (numChars % 2 == 0)
        {
          /* Shift to the right by one if we've got an even number of chars */
          --offset;
        }
  
      expandedRange->min = getRangeCentre(expandedRange) - offset;
      expandedRange->max = expandedRange->min + numChars - 1;
    }
}


/* Utility that returns true if the given bool is either true or not 
 * requested (i.e. the pointer is null) */
static gboolean trueOrNotRequested(const gboolean* const ptr)
{
  return (ptr == NULL || *ptr == TRUE);
}


/* Determine whether the given coord in the given frame/strand is affected by
 * a variation. Sets whether to draw the start/end/top/bottom boundaries for the coord
 * when outlining it by the extent of the found variation. Also sets drawBackground to true
 * if the background of the base should be highlighted in the variation color (i.e. if the base
 * is part of a selected variation). Sets multipleVariations to true if there is more than
 * one variation at this coord */
static gboolean coordAffectedByVariation(const int dnaIdx,
                                         const BlxStrand strand, 
                                         BlxViewContext *bc,
                                         const MSP **mspOut, /* the variation we found */
                                         gboolean *drawStartBoundary, 
                                         gboolean *drawEndBoundary, 
                                         gboolean *drawJoiningLines, 
                                         gboolean *drawBackground,
                                         gboolean *multipleVariations)
{
  gboolean result = FALSE;
  
  /* Loop through all variations */
  int i = 0;
  const MSP *msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], i);
  
  for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], ++i))
    {
      if (mspGetRefStrand(msp) == strand && valueWithinRange(dnaIdx, &msp->qRange))
        {
          /* If result has already been set, then there are multiple variations on this coord */
          if (result && multipleVariations)
            *multipleVariations = TRUE;

          result = TRUE;
          
          if (mspOut)
            *mspOut = msp;
          
          /* We draw the start boundary of this base if it's the first base in the variation's range 
           * and the end boundary if it's the last, unless it's a zero-length feature, i.e. an insertion
           * site; in this case we draw a line to the right of the coord in the direction of the
           * reference sequence (i.e. to give just a vertical line at the site, not an outline around the bases) */
          const gboolean isFirst = (dnaIdx == msp->qRange.min);
          const gboolean isLast = (dnaIdx == msp->qRange.max);
          const gboolean isZeroLen = mspIsZeroLenVariation(msp);

          /* Note that if any of the boundaries are already set to be drawn (by a previous snp
           * found in this loop) then we keep the current value */
          if (drawStartBoundary)
            *drawStartBoundary |= (isZeroLen != bc->displayRev ? isLast : isFirst);
            
          if (drawEndBoundary)
            *drawEndBoundary |= (isZeroLen != bc->displayRev ? isFirst : isLast);

          /* Don't draw joining lines between the start and end boundaries if it's zero-length */
          if (drawJoiningLines)
            *drawJoiningLines |= !isZeroLen;
            
          /* Highlight the background if the variation is selected (unless it's a zero-length variation) */
          if (drawBackground)
            *drawBackground |= !isZeroLen && blxContextIsSeqSelected(bc, msp->sSequence);
          
          /* If we only want to know if this coord is affected by a snp, we can return now.
           *
           * If we also need the boundary-lines info for this coord, then we need to check 
           * all snps that are at this coord, because some may have different boundaries to
           * others. We can break if all boundar-lines have already been set to true, though.
           *
           * If we need to return whether there are multiple snps at this coord, we can return
           * after we've found the second on (i.e. after multipleVariations is set to TRUE). */
          if (trueOrNotRequested(drawStartBoundary) &&
              trueOrNotRequested(drawEndBoundary) && 
              trueOrNotRequested(drawJoiningLines) && 
              trueOrNotRequested(drawBackground) &&
              trueOrNotRequested(multipleVariations))
            break;
        }
    }
  
  return result;
}


/* Determine whether the given coord in the given frame/strand is affected by
 * a polyA signal */
static gboolean coordAffectedByPolyASignal(const int dnaIdx,
                                           const BlxStrand strand, 
                                           BlxViewContext *bc,
                                           const MSP **mspOut, /* the variation we found */
                                           gboolean *multiple)
{
  gboolean result = FALSE;
  
  /* Loop through all variations */
  int i = 0;
  const MSP *msp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SIGNAL], i);
  
  for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SIGNAL], ++i))
    {
      if (mspGetRefStrand(msp) == strand && valueWithinRange(dnaIdx, &msp->qRange))
        {
          /* If result has already been set, then there are multiple signals on this coord */
          if (result && multiple)
            *multiple = TRUE;

          result = TRUE;
          
          if (mspOut)
            *mspOut = msp;
          
          /* If we only want to know if this coord is affected by a polyA signal, we can return now.
           *
           * If we need to return whether there are multiple sites at this coord, we can return
           * after we've found the second on (i.e. after multiple is set to TRUE). */
          if (trueOrNotRequested(multiple))
            break;
        }
    }
  
  return result;
}


/* Determine whether the given coord in the given frame/strand is affected by
 * a polyA site */
static gboolean coordAffectedByPolyASite(const int dnaIdx,
                                         const BlxStrand strand, 
                                         BlxViewContext *bc,
                                         const MSP **mspOut, /* the variation we found */
                                         gboolean *multiple)
{
  gboolean result = FALSE;
  
  /* Loop through all variations */
  int i = 0;
  const MSP *msp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SITE], i);
  
  for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_POLYA_SITE], ++i))
    {
      if (mspGetRefStrand(msp) == strand && valueWithinRange(dnaIdx, &msp->qRange))
        {
          /* If result has already been set, then there are multiple sites on this coord */
          if (result && multiple)
            *multiple = TRUE;

          result = TRUE;
          
          if (mspOut)
            *mspOut = msp;
          
          /* If we only want to know if this coord is affected by a polyA signal, we can return now.
           *
           * If we need to return whether there are multiple sites at this coord, we can return
           * after we've found the second on (i.e. after multiple is set to TRUE). */
          if (trueOrNotRequested(multiple))
            break;
        }
    }
  
  return result;
}


/* Draw a rectangle with an outline and a fill color */
static void drawRectangle(GdkDrawable *drawable, 
                          GdkGC *gc, 
                          GdkColor *fillColor,
                          GdkColor *outlineColor,
                          const int x, 
                          const int y,
                          const int width, 
                          const int height,
                          const gboolean drawLeft,
                          const gboolean drawRight,
                          const gboolean drawTop,
                          const gboolean drawBottom)
{
  const int lineWidth = 1;
  
  if (fillColor)
    {
      gdk_gc_set_foreground(gc, fillColor);
      gdk_draw_rectangle(drawable, gc, TRUE, x, y, width, height);
    }

  if (outlineColor)
    {
      gdk_gc_set_foreground(gc, outlineColor);
      gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      int w = width - lineWidth;
      int h = height - lineWidth;
      
      if (drawLeft)
        gdk_draw_line(drawable, gc, x, y, x, y + h);
        
      if (drawRight)
        gdk_draw_line(drawable, gc, x + w, y, x + w, y + h);

      if (drawTop)
        gdk_draw_line(drawable, gc, x, y, x + w, y);

      if (drawBottom)
        gdk_draw_line(drawable, gc, x, y + h, x + w, y + h);

      //gdk_draw_rectangle(drawable, gc, FALSE, x, y, width - lineWidth, height - lineWidth);
    }
}


/* Returns true if the given coord (which must be in nucleotide coords) is within the
 * range of any MSP that is selected. Only considers MSPs that are on the given ref seq
 * frame and strand (if given; otherwise considers all MSPs). */
static gboolean isCoordInSelectedMspRange(const BlxViewContext *bc,
                                          const int dnaIdx, 
                                          const BlxStrand refSeqStrand, 
                                          const int refSeqFrame,
                                          const BlxSeqType seqType)
{
  gboolean inSelectedMspRange = FALSE;
  
  /* Loop through all the selected sequences */
  GList *seqItem = bc->selectedSeqs;
  
  for ( ; seqItem && !inSelectedMspRange; seqItem = seqItem->next)
    {
      /* Loop through all the MSPs in this sequence */
      const BlxSequence* const blxSeq = (const BlxSequence*)(seqItem->data);
      GList *mspItem = blxSeq->mspList;
      
      for ( ; mspItem && !inSelectedMspRange; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);
          const int mspFrame = mspGetRefFrame(msp, seqType);
          
          if ((mspIsBlastMatch(msp) || mspIsExon(msp)) && 
              (refSeqStrand == BLXSTRAND_NONE || mspGetRefStrand(msp) == refSeqStrand) &&
              (refSeqFrame == UNSET_INT || mspFrame == refSeqFrame))
            {
              /* Passed-in coord is the nucleotide for base 1 in frame 1 of the current display index.
               * We need to compare against base1 for the frame for the individual MSP, though, so 
               * we need to convert it: just add 1 base if we're in frame 2 or 2 bases if we're in frame3. */
              int idx = dnaIdx;
              
              if (bc->displayRev)
                {
                  idx -= (mspFrame - 1);
                }
              else
                {
                  idx += (mspFrame - 1);
                }
              
              inSelectedMspRange = valueWithinRange(idx, &msp->qRange);
            }
        }
    }

  return inSelectedMspRange;
}


/* Determine whether any highlighting needs to be done for snps in the ref seq */
static void getSnpHighlighting(DrawBaseData *data,
                               BlxViewContext *bc)
{     
  gboolean drawBackground = FALSE;
  
  if (coordAffectedByVariation(data->dnaIdx, data->strand, bc, NULL,
                               &data->drawStart, &data->drawEnd, &data->drawJoiningLines, &drawBackground, NULL))
    {
      /* The coord is affected by a SNP. Outline it in the "selected" SNP color
       * (which is darker than the normal color) */
      data->outlineColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, TRUE, bc->usePrintColors);
      
      /* If the SNP is selected, also fill it with the SNP color (using the
       * "unselected" SNP color, which is lighter than the outline). */
      if (drawBackground)
        {
          data->fillColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, data->shadeBackground, bc->usePrintColors);
        }
    }
}


/* Draw a given nucleotide or peptide. Determines the color depending on various
 * parameters */
void drawHeaderChar(BlxViewContext *bc,
                    DetailViewProperties *properties,
                    GdkDrawable *drawable,
                    GdkGC *gc,
                    const int x,
                    const int y,
                    GHashTable *basesToHighlight,
                    DrawBaseData *data)
{
  /* Shade the background if the base is selected XOR if the base is within the range of a 
   * selected sequence. (If both conditions are true we don't shade, to give the effect of an
   * inverted selection color.) */
  gboolean inSelectedMspRange = isCoordInSelectedMspRange(bc, data->dnaIdx, data->strand, data->frame, data->seqType);
  data->shadeBackground = (data->displayIdxSelected != inSelectedMspRange);

  /* Reset background color and outlines to defaults */
  data->outlineColor = NULL;
  data->fillColor = NULL;
  data->drawStart = FALSE;
  data->drawEnd = FALSE;
  data->drawJoiningLines = FALSE;
 
  /* Check if this coord already has a special color stored for it */
  gpointer hashValue = g_hash_table_lookup(basesToHighlight, GINT_TO_POINTER(data->dnaIdx));

  if (hashValue)
    {
      BlxColorId colorId = (BlxColorId)GPOINTER_TO_INT(hashValue);
      GdkColor *color = getGdkColor(colorId, bc->defaultColors, data->shadeBackground, bc->usePrintColors);

      if (colorId == BLXCOLOR_POLYA_SITE_ANN)
        {
          /* For polyA sites we don't shade the background; instead we want to draw a bar after the
           * coord (or before if the display is reversed), so set the outline. 
           * (NB Bit of a hack to use the colorId to identify polyA sites here but it works fine.) */ 
          data->outlineColor = color;
          
          if (bc->displayRev)
            data->drawStart = TRUE;
          else
            data->drawEnd = TRUE;
        }
      else
        {
          /* Normal highlighting: just fill in the background with the specified color */
          data->fillColor = color;
        }
    }

  /* Check if this base is in the currently-selected codon and needs highlighting */
  if (data->seqType == BLXSEQ_DNA && data->showCodons && (data->dnaIdxSelected || data->displayIdxSelected || data->shadeBackground))
    {
      if (data->dnaIdxSelected || data->displayIdxSelected)
        {
          /* The coord is a nucleotide in the currently-selected codon. The color depends
           * on whether the actual nucleotide itself is selected, or just the codon that it 
           * belongs to. */
          data->fillColor = getGdkColor(BLXCOLOR_CODON, bc->defaultColors, data->dnaIdxSelected, bc->usePrintColors);
          
        }
      else if (!data->fillColor)
        {
          /* The coord is not selected but this coord is within the range of a selected MSP, so 
           * shade the background. */
          data->fillColor = getGdkColor(data->defaultBgColor, bc->defaultColors, data->shadeBackground, bc->usePrintColors);
        }
    }

  /* Check if this base is a SNP (or other variation) and needs highlighting */
  if (data->highlightSnps)
    {
      getSnpHighlighting(data, bc);
    }
  
  /* If the base is not already assigned some special highlighting, then
   * check whether it should be highlighted as a stop or MET. Otherwise,
   * give it the default background colour. */
  if (!data->fillColor)
    {
      if (data->seqType == BLXSEQ_PEPTIDE && data->baseChar == SEQUENCE_CHAR_MET)
        {
          /* The coord is a MET codon */
          data->fillColor = getGdkColor(BLXCOLOR_MET, bc->defaultColors, data->shadeBackground, bc->usePrintColors);
        }
      else if (data->seqType == BLXSEQ_PEPTIDE && data->baseChar == SEQUENCE_CHAR_STOP)
        {
          /* The coord is a STOP codon */
          data->fillColor = getGdkColor(BLXCOLOR_STOP, bc->defaultColors, data->shadeBackground, bc->usePrintColors);
        }
      else if (data->showBackground)
        {
          /* Use the default background color for the reference sequence */
          data->fillColor = getGdkColor(data->defaultBgColor, bc->defaultColors, data->shadeBackground, bc->usePrintColors);
        }
    }
  
  /* Ok, now we know all the colours and outlines, draw the background for the base */
  if (data->topToBottom)
    {
      /* We're drawing nucleotides from top-to-bottom instead of left-to-right, so the start border is
       * the top border and the bottom border is the end border. */
      drawRectangle(drawable, gc, data->fillColor, data->outlineColor, x, y, ceil(properties->charWidth), roundNearest(properties->charHeight),
                    data->drawJoiningLines, data->drawJoiningLines, data->drawStart, data->drawEnd);
    }
  else
    {
      drawRectangle(drawable, gc, data->fillColor, data->outlineColor, x, y, ceil(properties->charWidth), roundNearest(properties->charHeight),
                    data->drawStart, data->drawEnd, data->drawJoiningLines, data->drawJoiningLines);
    }
}


/* This works out the row number (1-based) to draw a variation with the given range
 * in. Variations are drawn in row 1 if possible, but if the given range overlaps
 * another variation in row 1, then we will try row 2, then row 3 etc.
 * If numRows is UNSET_INT, then the straightforward row calculation is returned.
 * If numRows is given, then the result is inverted (i.e. if there are three rows, rows
 * 1 is translated to 3, 2->2 and 3->1. This is so that we can draw the first row at the
 * bottom). */
static int getVariationRowNumber(const IntRange* const rangeIn, const int numRows, GSList **rows)
{
  /* Loop through each row and add it to the first row where it will not overlap any other */
  GSList *row = *rows;
  int rowNum = 1;

  for ( ; row; row = row->next, ++rowNum)
    {
      /* See if it overlaps an existing range */
      GSList *ranges = (GSList*)(row->data);
      GSList *rangeItem = ranges;
      
      gboolean overlaps = FALSE;
  
      for ( ; rangeItem && !overlaps; rangeItem = rangeItem->next)
        {
          IntRange *range = (IntRange*)(rangeItem->data);
          overlaps = rangesOverlap(range, rangeIn);
        }

      if (!overlaps)
        {
          /* Doesn't overlap anything in this row; add it to the row and exit */
          IntRange *range = (IntRange*)g_malloc(sizeof *range);
          range->min = rangeIn->min;
          range->max = rangeIn->max;
          ranges = g_slist_append(ranges, range);
          
          row->data = ranges;

          if (numRows != UNSET_INT)
            rowNum = numRows - rowNum + 1; /* invert row order */
          
          return rowNum;
        }
    }

  /* If we get here, we didn't find a row that we don't overlap, so add a new row */
  GSList *ranges = NULL;
  IntRange *range = (IntRange*)g_malloc(sizeof *range);
  range->min = rangeIn->min;
  range->max = rangeIn->max;
  ranges = g_slist_append(ranges, range);
  *rows = g_slist_append(*rows, ranges);

  if (numRows != UNSET_INT)
    rowNum = numRows - rowNum + 1; /* invert row order */

  return rowNum;
}


/* This deletes the memory allocated to the 'rows' array created
 * when drawing the snp track */
static void freeRowsList(GSList *rows)
{
  GSList *rowItem = rows;
  for ( ; rowItem; rowItem = rowItem->next)
    {
      GSList *ranges = (GSList*)(rowItem->data);
      GSList *rangeItem = ranges;
      
      for ( ; rangeItem; rangeItem = rangeItem->next)
        g_free(rangeItem->data);

      g_slist_free(ranges);
    }

  g_slist_free(rows);
}


/* Utility to calculate the number of rows we need in the SNP track to avoid 
 * any overlapping variations. */
static int getNumSnpTrackRows(const BlxViewContext *bc, 
                              DetailViewProperties *properties, 
                              const BlxStrand strand, 
                              const int frame)
{
  int numRows = 1; /* Always show one row, even if it is empty */

  int i = 0;
  const MSP* msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], i);
  GSList *rows = NULL;

  for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], ++i))
    {
      if (mspGetRefStrand(msp) == strand && mspGetMatchSeq(msp))
        {
          /* Get the range of display coords where this SNP will appear */
          IntRange mspDisplayRange, mspExpandedRange;
          getVariationDisplayRange(msp, TRUE, bc->seqType, bc->numFrames, bc->displayRev, frame, &bc->refSeqRange, &mspDisplayRange, &mspExpandedRange);

          /* See if the variation is in the current display range */
          if (rangesOverlap(&mspExpandedRange, &properties->displayRange))
            {
              /* Get the row number (1-based) to draw this variation in. */
              int row = getVariationRowNumber(&mspExpandedRange, UNSET_INT, &rows);
              
              /* Keep track of how many rows we've seen so we can size the widget accordingly */
              if (row > numRows)
                numRows = row;
            }
        }
    }

  freeRowsList(rows);
  
  return numRows;
}


/* Recalculate the size of the variations track. This needs to be done whenever
 * we turn the track on/off, and also whenever we we scroll, because more/less
 * variations can appear and if they overlap then we add extra rows in which to
 * draw them. */
static void recalculateSnpTrackBorders(GtkWidget *snpTrack, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *blxWindow = properties->blxWindow;
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  if (bc->flags[BLXFLAG_SHOW_VARIATION_TRACK])
    {
      const int activeFrame = detailViewGetActiveFrame(detailView);
      BlxStrand strand = snpTrackGetStrand(snpTrack, detailView);

      const int rowHeight = ceil(properties->charHeight);
      const int numRows = getNumSnpTrackRows(bc, properties, strand, activeFrame);

      gtk_widget_set_size_request(snpTrack, -1, numRows * rowHeight);
    }
  else
    {
      gtk_widget_set_size_request(snpTrack, -1, 0);
    }
}


/* Function that does the drawing for the variations track */
static void drawVariationsTrack(GtkWidget *snpTrack, GtkWidget *detailView)
{
  /* Create the drawable for the widget (whether we're actually going to do any drawing or not) */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GdkDrawable *drawable = createBlankSizedPixmap(snpTrack, snpTrack->window, snpTrack->allocation.width, ceil(properties->charHeight) * 10);

  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (!bc->flags[BLXFLAG_SHOW_VARIATION_TRACK] || !bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS])
    {
      return;
    }
  
  const int activeFrame = detailViewGetActiveFrame(detailView);
  
  BlxStrand strand = snpTrackGetStrand(snpTrack, detailView);
  
  GdkGC *gc = gdk_gc_new(drawable);
  
  /* Find the left margin. It will be at the same x coord as the left edge of
   * the sequence column header. */
  int leftMargin = UNSET_INT;
  BlxColumnInfo *seqColInfo = getColumnInfo(detailViewGetColumnList(detailView), BLXCOL_SEQUENCE);
  gtk_widget_translate_coordinates(seqColInfo->headerWidget, snpTrack, 0, 0, &leftMargin, NULL);
  
  /* Maintain lists for each row where the variations are drawn; remember their display
   * ranges and don't allow any to overlap. */
  const int numRows = getNumSnpTrackRows(bc, properties, strand, activeFrame);

  /* Loop through all variations and see if any are in the current display range */
  const int y = 0;
  const int rowHeight = ceil(properties->charHeight);
  int i = 0;
  const MSP *msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], i);
  GSList *rows = NULL;

  for ( ; msp; msp = mspArrayIdx(bc->featureLists[BLXMSP_VARIATION], ++i))
    {
      if (mspGetRefStrand(msp) == strand && mspGetMatchSeq(msp))
        {
          /* Get the range of display coords where this SNP will appear */
          IntRange mspDisplayRange, mspExpandedRange;
          getVariationDisplayRange(msp, TRUE, bc->seqType, bc->numFrames, bc->displayRev, activeFrame, &bc->refSeqRange, &mspDisplayRange, &mspExpandedRange);

          /* See if the variation is in the current display range */
          if (rangesOverlap(&mspExpandedRange, &properties->displayRange))
            {
              /* Get the row number to draw this variation in. */
              const int rowNum = getVariationRowNumber(&mspExpandedRange, numRows, &rows);
              const int rowIdx = rowNum - 1;
              
              int x = leftMargin + (int)((gdouble)(mspExpandedRange.min - properties->displayRange.min) * properties->charWidth);
              const int width = ceil((gdouble)strlen(mspGetMatchSeq(msp)) * properties->charWidth);
              const gboolean isSelected = blxWindowIsSeqSelected(blxWindow, msp->sSequence);
              
              /* Draw the outline in the default SNP color. If the SNP is selected, also
               * fill in the rectangle in the SNP color (use the selected color for the
               * outline and the unselected color for the fill, so that the outline is darker). */
              GdkColor *outlineColor = getGdkColor(BLXCOLOR_SNP, bc->defaultColors, TRUE, bc->usePrintColors);
              GdkColor *fillColor = isSelected ? getGdkColor(BLXCOLOR_SNP, bc->defaultColors, FALSE, bc->usePrintColors) : NULL;
              
              /* Draw the background rectangle for the char */
              drawRectangle(drawable, gc, fillColor, outlineColor, x, y + (rowHeight * rowIdx), width, rowHeight, TRUE, TRUE, TRUE, TRUE);
              
              /* Draw the text */
              PangoLayout *layout = gtk_widget_create_pango_layout(detailView, mspGetMatchSeq(msp));
              pango_layout_set_font_description(layout, properties->fontDesc);

              if (layout)
                {
                  gtk_paint_layout(snpTrack->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, x, y + (rowHeight * rowIdx), layout);
                  g_object_unref(layout);
                }             
            }
        }
    }

  freeRowsList(rows);
  
  drawColumnSeparatorLine(snpTrack, drawable, gc, bc);
  
  g_object_unref(gc);
}


/* This function centres the detail view display on the currently-selected base index. Does 
 * nothing if there isn't a selected base. */
static void detailViewCentreOnSelection(GtkWidget *detailView)
{
  const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
  
  if (selectedBaseIdx != UNSET_INT)
    {
      /* The coord is in terms of the display coords, i.e. whatever the displayed seq type is. */
      const BlxSeqType seqType = detailViewGetSeqType(detailView);
      const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
      
      int newStart = selectedBaseIdx - (getRangeLength(displayRange) / 2);
      setDetailViewStartIdx(detailView, newStart, seqType);
    }
}


/* In the detail view, get the base index at the given coords, if those coords lie within the
 * sequence column; returns UNSET_INT otherwise. */
static int getBaseIndexAtDetailViewCoords(GtkWidget *detailView, const int x, const int y)
{
  int baseIdx = UNSET_INT;
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GList *columnList = detailViewGetColumnList(detailView);

  /* Get the x coords at the start/end of the sequence column */
  IntRange xRange;
  getColumnXCoords(columnList, BLXCOL_SEQUENCE, &xRange);
  
  /* See if our x coord lies inside the sequence column */
  if (x >= xRange.min && x <= xRange.max)
    {
      /* Get the 0-based char index at x */
      gdouble charWidth = properties->charWidth;
      int charIdx = (int)(((gdouble)x - xRange.min) / charWidth);

      /* Add the start of the scroll range to convert this to the display index */
      GtkAdjustment *adjustment = properties->adjustment;
      baseIdx = charIdx + adjustment->value;
    }
  
  return baseIdx;
}


/***********************************************************
 *              Detail view scrollbar events               *
 ***********************************************************/

/* Callback called when the detail-view scrollbar's range has changed */
static void onScrollRangeChangedDetailView(GtkObject *object, gpointer data)
{
  DEBUG_ENTER("onScrollRangeChangedDetailView");

  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  IntRange *displayRange = detailViewGetDisplayRange(detailView);

  /* First time round, set the adjusment range to be centred on the centre of
   * the detail-view range*/
  if (adjustment->value == UNSET_INT)
    {
      adjustment->value = getRangeCentre(displayRange) - (adjustment->page_size / 2);

      if (adjustment->value < adjustment->lower)
        adjustment->value = adjustment->lower;
      else if (adjustment->value > adjustment->upper - adjustment->page_size + 1)
        adjustment->value = adjustment->upper - adjustment->page_size + 1;
    }

  int newStart = adjustment->value;
  int newEnd = newStart + adjustment->page_size - 1;
  
  /* Only update if something has changed */
  if (displayRange->min != newStart || displayRange->max != newEnd)
    {
      IntRange oldRange = {displayRange->min, displayRange->max};

      displayRange->min = newStart;
      displayRange->max = newEnd;
      
      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      refilterDetailView(detailView, &oldRange);

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
      detailViewRedrawAll(detailView);

      /* Update the big picture because the highlight box has moved (and changed size) */
      GtkWidget *bigPicture = detailViewGetBigPicture(detailView);
      refreshBigPictureDisplayRange(bigPicture, FALSE);
      calculateBigPictureCellSize(bigPicture, bigPictureGetProperties(bigPicture));
    }
  
  DEBUG_EXIT("onScrollRangeChangedDetailView returning");
}


/* Callback called when the detail-view scrollbar's position has changed */
static void onScrollPosChangedDetailView(GtkObject *object, gpointer data)
{
  DEBUG_ENTER("onScrollPosChangedDetailView");

  GtkAdjustment *adjustment = GTK_ADJUSTMENT(object);
  GtkWidget *detailView = GTK_WIDGET(data);
  
  /* Set the display range so that it starts at the new scroll pos */
  int newStart = adjustment->value;
  int newEnd = adjustment->value + adjustment->page_size - 1;

  /* Only update if something has changed */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  IntRange *displayRange = &properties->displayRange;
  
  if (displayRange->min != newStart || displayRange->max != newEnd)
    {
      IntRange oldRange = {displayRange->min, displayRange->max};

      displayRange->min = newStart;
      displayRange->max = newEnd;

      /* Refilter the data for all trees in the detail view because rows may have scrolled in/out of view */
      refilterDetailView(detailView, &oldRange);

      /* Refresh the detail view header (which may contain the DNA sequence), and 
       * the headers for all the trees (which contains the reference sequence) */
      refreshDetailViewHeaders(detailView);
      callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
      detailViewRedrawAll(detailView);

      /* Update the big picture because the highlight box has moved */
      GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
      refreshBigPictureDisplayRange(bigPicture, FALSE);
    }
  
  DEBUG_EXIT("onScrollPosChangedDetailView returning");
}


/***********************************************************
 *            Detail View toolbar events                   *
 ***********************************************************/

/* Set the font for the detail view cell renderer. Should be called after
 * the font size is changed by zooming in/out of the detail view. */
static void updateCellRendererFont(GtkWidget *detailView, PangoFontDescription *fontDesc)
{
  /* Calculate the row height from the font size */
  gdouble charWidth, charHeight;
  getFontCharSize(detailView, fontDesc, &charWidth, &charHeight);
  
  /* Cache these results, because we use them often for calculations */
  detailViewCacheFontSize(detailView, charWidth, charHeight);
  
  /* Set the row height. Subtract the padding between the cell's actual area and
   * its background area. We will render at the background area's height, so that
   * we draw over the "gaps" between the cells, giving the impression of no gaps. */
  gint rowHeight = roundNearest(charHeight) - (detailViewGetCellYPadding(detailView) * 2);

  GtkCellRenderer *renderer = detailViewGetRenderer(detailView);
  gtk_cell_renderer_set_fixed_size(renderer, 0, rowHeight);
}


/* Refilter the tree row that this msp is in */
static void refilterMspRow(MSP *msp, GtkWidget *detailView, BlxViewContext *bc)
{
  DEBUG_ENTER("refilterMspRow(msp=[%d,%d])", msp->fullRange.min, msp->fullRange.max);

  /* Find the tree row that this MSP is in and force that row to update
   * its visibility status. */
  gchar *pathStr = mspGetTreePath(msp, bc->modelId);
  
  if (pathStr)
    {
      GtkTreePath *path = gtk_tree_path_new_from_string(pathStr);
      
      if (path)
        {
          GtkWidget *tree = detailViewGetTree(detailView, mspGetRefStrand(msp), mspGetRefFrame(msp, bc->seqType));
          GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
          GtkTreeIter iter;
          
          if (gtk_tree_model_get_iter(model, &iter, path))
            {
              gtk_tree_model_row_changed(model, path, &iter);
            }
        }

      gtk_tree_path_free(path);
    }
  
  DEBUG_EXIT("refilterMspRow returning ");
}


/* Utility to see if the start value of the first range is within the second
 * range. If 'rev' is true, the 'start' of range1 is the max value rather than
 * the min value of this range. */
static gboolean startValueWithinRange(const IntRange* const range1,
                                      const IntRange* const range2,
                                      const gboolean rev)
{
  gboolean result = FALSE;
  
  if (rev)
    result = valueWithinRange(range1->max, range2);
  else 
    result = valueWithinRange(range1->min, range2);
  
  return result;
}


/* Quick search to find any MSP in the given array that whose start coord
 * lies within the given range. Sets the found array index in 'idx' or returns 
 * FALSE if no MSP was found in this range. */
static gboolean getAnyMspInRange(GArray *mspArray, 
                                 const IntRange* const range, 
                                 const gboolean displayRev, 
                                 int *idx)
{
  DEBUG_ENTER("getAnyMspInRange");

  gboolean result = FALSE;
  
  /* Do a binary search based on start coord until we find a start coord 
   * that lies within the given range. */
  int iMax = mspArray->len - 1;
  int iMin = 0;

  while (1)
    {
      const int i = iMin + (iMax - iMin) / 2;  
      const MSP *msp = mspArrayIdx(mspArray, i);
      
      if (!msp)
        break;
    
      if (startValueWithinRange(&msp->displayRange, range, displayRev))
        {
          result = TRUE;
          *idx = i;
          break;
        }
    
      if (iMax - iMin < 1)
        break;
    
      if ((msp->displayRange.min < range->min) != displayRev)
        {
          iMin = i + 1;
        }
      else
        {
          iMax = i - 1;
        }
    }
  
  
  DEBUG_EXIT("getAnyMspInRange returning %d (idx=%d)", result, *idx);
  return result;
}


/* Refilter the tree rows for all MSPs in the given array whose start coords
 * lie within the given range. 'startIdx' should be an index in the array that
 * gives the position of an MSP that is known to lie within this range (i.e. 
 * it gives us a rough starting point so that we don't have to search through
 * the entire array; it does not necessarily have to be the FIRST MSP that lies
 * within range). */
static void refilterMspList(const int startIdx, 
                            GArray *array,
                            const IntRange* const range, 
                            GtkWidget *detailView, 
                            BlxViewContext *bc)
{
  DEBUG_ENTER("refilterMspList(startIdx=%d, range=[%d,%d])", startIdx, range->min, range->max);

  /* Loop forwards through the array from the start index, refiltering the rows
   * for each MSP until we find an MSP that is out of range. */
  int i = startIdx;
  MSP *msp = mspArrayIdx(array, i);
  
  for ( ; msp; msp = mspArrayIdx(array, ++i))
    {
      if (!startValueWithinRange(&msp->displayRange, range, bc->displayRev))
        break;
      
      refilterMspRow(msp, detailView, bc);
    }
  
  /* Also loop backwards from the start index, because startIdx may not have 
   * been the very first MSP in the array that was in range. */
  i = startIdx - 1;
  msp = mspArrayIdx(array, i);
  
  for ( ; msp; msp = mspArrayIdx(array, --i))
    {
      if (!startValueWithinRange(&msp->displayRange, range, bc->displayRev))
        break;
      
      refilterMspRow(msp, detailView, bc);
    }
  
  DEBUG_EXIT("refilterMspList returning ");
}


/* Refilter rows in the detail-view trees for MSPs of the given type. Only filter
 * the rows of MSPs that lie within the given display range. */
static void refilterMspArrayByRange(BlxMspType mspType, 
                                   BlxViewContext *bc, 
                                   const IntRange* const displayRange, 
                                   GtkWidget *detailView)
{
  DEBUG_ENTER("refilterMspArrayByRange(mspType=%d)", mspType);

  if (!displayRange)
    return;
  
  /* We only want to update MSPs that are within the given display range.
   * 
   * We want to avoid searching through the entire MSP array, because it may 
   * contain many thousands of MSPs. However, it is difficult to do a quick
   * search to find an item from the list based on two values (i.e. the start and
   * end of its range). 
   * What we do is extend the start of the display range by the maximum MSP
   * length so that we can do a binary search just on the start coord of the MSP 
   * to check that it lies within the extended range. For any MSP that lies within
   * the display range, it must be true that its start coord lies within this 
   * extended range.
   * Because the extended range is based on the maximum MSP length, it is likely
   * to include some MSPs that do not actually lie within the display range.
   * However, there is no harm in filtering a few extra rows, and this should 
   * still be more efficient that searching through the entire array of MSPs. */

  const int maxLen = getMaxMspLen();

  /* If the display is reversed, the 'start' is the maximum coord end 
   * rather than the minimum coord */
  IntRange extendedRange;
  
  if (bc->displayRev)
    {
      extendedRange.min = displayRange->min;
      extendedRange.max = displayRange->max + maxLen;
    }
  else
    {
      extendedRange.min = displayRange->min - maxLen;
      extendedRange.max = displayRange->max;
    }
  
  int startIdx = 0;
  
  if (getAnyMspInRange(bc->featureLists[mspType], &extendedRange, bc->displayRev, &startIdx))
    {
      refilterMspList(startIdx, bc->featureLists[mspType], &extendedRange, detailView, bc);
    }
  
  DEBUG_EXIT("refilterMspArrayByRange returning ");
}


/* Refilter the rows in the detail-view trees. Only updates the rows in the
 * current display range and (if given) the old display range. */
void refilterDetailView(GtkWidget *detailView, const IntRange* const oldRange)
{
  DEBUG_ENTER("refilterDetailView(oldRange=[%d,%d])", oldRange ? oldRange->min : 0, oldRange ? oldRange->max : 0);
  
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  if (bc->flags[BLXFLAG_SHOW_POLYA_SITE] || bc->flags[BLXFLAG_SHOW_UNALIGNED])
    {
      /* If these options are enabled, then determining which rows need to be 
       * refiltered becomes much more complex because visible alignments lengths
       * may differ from the actual alignment coords, and may differ depending
       * on whether the sequence is selected or not. For now, play safe and
       * refilter everything. */
      callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
    }
  else
    {
      /* We only want to update MSPs that are in the old detail-view range
       * and the new detail-view range. (Because updating every row in all
       * trees can be very slow if there are a lot of MSPs.) */
      const IntRange* const newRange = detailViewGetDisplayRange(detailView);

      /* Only consider MSPs that are shwon in the detail-view */
      int mspType = 0;
      for ( ; mspType < BLXMSP_NUM_TYPES; ++mspType)
        {
          if (typeShownInDetailView((BlxMspType)mspType))
            {
              refilterMspArrayByRange((BlxMspType)mspType, bc, oldRange, detailView);
              refilterMspArrayByRange((BlxMspType)mspType, bc, newRange, detailView);
            }
        }
      
      detailViewRedrawAll(detailView);
    }

  DEBUG_EXIT("refilterDetailView returning ");
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static void assertDetailView(GtkWidget *detailView)
{
  /* Check it's a valid detail-view tree type */
  if (!detailView)
    g_error("Detail-view widget is null\n");
  
  if (!GTK_IS_WIDGET(detailView))
    g_error("Detail-view is not a valid widget [%p]\n", detailView);
  
  if (!GTK_IS_CONTAINER(detailView))
    g_error("Detail-view is not a valid container [%p]\n", detailView);
  
  if (!detailViewGetProperties(detailView))
    g_error("Tree properties not set [widget=%p]\n", detailView);
}

GtkWidget* detailViewGetBlxWindow(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->blxWindow : NULL;
}

GtkAdjustment* detailViewGetAdjustment(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->adjustment : NULL;
}

GtkCellRenderer* detailViewGetRenderer(GtkWidget *detailView)
{
  assertDetailView(detailView);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->renderer : NULL;
}

/* Get the list of columns */
GList* detailViewGetColumnList(GtkWidget *detailView)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  return bc ? bc->columnList : NULL;
}

/* For the given list of columns, extract their types into an
 * array of 'GType's */
GType* columnListGetTypes(GList *columnList)
{
  GType *result = NULL;
  const int len = g_list_length(columnList);

  if (len > 0)
    {
      result = (GType*)g_malloc(len * sizeof(GType));

      GList *item = columnList;
      int i = 0;
      
      for ( ; item; item = item->next, ++i)
        {
          BlxColumnInfo* columnInfo = (BlxColumnInfo*)(item->data);
          result[i] = columnInfo->type;
        }
    }
  
  return result;
}

gboolean detailViewGetDisplayRev(GtkWidget *detailView)
{
  return blxWindowGetDisplayRev(detailViewGetBlxWindow(detailView));
}

PangoFontDescription *detailViewGetFontDesc(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->fontDesc : NULL;
}

int detailViewGetCellXPadding(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->cellXPadding : 0;
}

int detailViewGetCellYPadding(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->cellYPadding : 0;
}

/* Get the character width. */
gdouble detailViewGetCharWidth(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->charWidth : 0.0;
}

gdouble detailViewGetCharHeight(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->charHeight : 0.0;
}

static void detailViewCacheFontSize(GtkWidget *detailView, gdouble charWidth, gdouble charHeight)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->charWidth = charWidth;
  properties->charHeight = charHeight;
}


/* Get the outermost container for the GtkTreeView for the given frame on the given strand */
GtkWidget* detailViewGetTreeContainer(GtkWidget *detailView, const BlxStrand activeStrand, const int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GtkWidget *result = NULL;
  
  /* Get the list of trees for the relevant strand */
  GList *list = (activeStrand == BLXSTRAND_FORWARD) ? properties->fwdStrandTrees : properties->revStrandTrees;

  /* The list should be sorted in order of frame number, and we should not be requesting a frame number
   * greater than the number of items in the list */
  if (frame <= (int)g_list_length(list))
    {
      GList *listItem = list;
      int count = 1;
      
      for ( ; listItem && count < frame; ++count)
        {
          listItem = listItem->next;
        }

      if (count == frame && 
          GTK_IS_WIDGET(listItem->data) && 
          (widgetIsTree(GTK_WIDGET(listItem->data)) || widgetIsTreeContainer(GTK_WIDGET(listItem->data))))
        {
          result = GTK_WIDGET(listItem->data);
        }
    }
  
  return result;
}


/* Extract the actual GtkTreeView for the given frame on the given strand */
GtkWidget* detailViewGetTree(GtkWidget *detailView, const BlxStrand activeStrand, const int frame)
{
  GtkWidget *result = NULL;
  GtkWidget *treeContainer = detailViewGetTreeContainer(detailView, activeStrand, frame);
  
  if (treeContainer && widgetIsTree(treeContainer))
    {
      result = treeContainer;
    }
  else if (treeContainer)
    {
      result = treeContainerGetTree(GTK_CONTAINER(treeContainer));
    }
  
  if (!result)
    {
      g_warning("Tree not found for '%s' strand, frame '%d'. Returning NULL.\n", ((activeStrand == BLXSTRAND_FORWARD) ? "forward" : "reverse"), frame);
    }
  
  return result;
}

int detailViewGetNumFrames(GtkWidget *detailView)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  return bc->numFrames;
}

static GtkWidget* detailViewGetHeader(GtkWidget *detailView)
{
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  return detailViewProperties->header;
}

static GtkWidget* detailViewGetFeedbackBox(GtkWidget *detailView)
{
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  return detailViewProperties->feedbackBox;
}

IntRange* detailViewGetDisplayRange(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? &properties->displayRange : NULL;
}

IntRange* detailViewGetRefSeqRange(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetRefSeqRange(blxWindow);
}

BlxSeqType detailViewGetSeqType(GtkWidget *detailView)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  return bc->seqType;
}

int detailViewGetSelectedBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedBaseIdx : UNSET_INT;
}

gboolean detailViewGetSelectedBaseSet(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedBaseSet : FALSE;
}

int detailViewGetSelectedDnaBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedDnaBaseIdx : UNSET_INT;
}

/* Get the active frame. Returns the last-selected frame, or 1 if no frame is selected. */
int detailViewGetActiveFrame(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  int result = 1;
  
  if (properties && properties->selectedFrame != UNSET_INT)
    {
      result = properties->selectedFrame;
}

  return result;
}

/* Get the strand of the tree that was last selected (which defaults to the active strand if none is selected). */
BlxStrand detailViewGetSelectedStrand(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->selectedStrand : BLXSTRAND_FORWARD;
}

/* Set the strand of the tree that was last selected. */
void detailViewSetSelectedStrand(GtkWidget *detailView, BlxStrand strand)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  properties->selectedStrand = strand;
}


int detailViewGetNumUnalignedBases(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->numUnalignedBases : FALSE;
}


static int detailViewGetSnpConnectorHeight(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->snpConnectorHeight : UNSET_INT;
}


BlxColumnId* detailViewGetSortColumns(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties ? properties->sortColumns : NULL;
}


/* Perform required updates after the selected base has changed. */
static void updateFollowingBaseSelection(GtkWidget *detailView,
                                         const gboolean allowScroll,
                                         const gboolean scrollMinimum)
{
  if (allowScroll)
    {
      scrollToKeepSelectionInRange(detailView, scrollMinimum);
    }
  
  /* Update the feedback box */
  updateFeedbackBox(detailView);
  
  /* Update the headers so that the newly-selected index is highlighted */
  detailViewRefreshAllHeaders(detailView);
}


/* Cancel the current base index selection */
void detailViewUnsetSelectedBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  properties->selectedBaseSet = FALSE;
  properties->selectedDnaBaseIdx = UNSET_INT;
  properties->selectedBaseIdx = UNSET_INT;
  properties->selectedBaseNum = UNSET_INT;
  
  updateFollowingBaseSelection(detailView, FALSE, FALSE);
}


/* Set the selected base index to a specific DNA index. Performs any required 
 * refreshes. Scrolls the view to keep the selected base in view if allowScroll 
 * is true. (Such scrolling is by the minimum number of bases necessary if 
 * scrollMinimum is true.) */
void detailViewSetSelectedDnaBaseIdx(GtkWidget *detailView,
                                     const int selectedDnaBaseIdx,
                                     const int frame,
                                     const gboolean allowScroll,
                                     const gboolean scrollMinimum)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  properties->selectedBaseSet = TRUE;
  properties->selectedDnaBaseIdx = selectedDnaBaseIdx;
  properties->selectedFrame = frame;
  
  /* For protein matches, calculate the display index and base number of this dna idx */
  BlxViewContext *bc = detailViewGetContext(detailView);
  properties->selectedBaseIdx = convertDnaIdxToDisplayIdx(selectedDnaBaseIdx, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &properties->selectedBaseNum);
  
  updateFollowingBaseSelection(detailView, allowScroll, scrollMinimum);
}


/* Set which frame is the active (currently-selected) reading frame */
void detailViewSetActiveFrame(GtkWidget *detailView, const int frame)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  properties->selectedFrame = frame;
  
  /* Keep the selected DNA coord the same and update the base number for the
   * new reading frame */
  if (properties->selectedBaseSet)
    {
      properties->selectedBaseIdx = convertDnaIdxToDisplayIdx(properties->selectedDnaBaseIdx, bc->seqType, frame, bc->numFrames, 
                                                              bc->displayRev, &bc->refSeqRange, &properties->selectedBaseNum);
    }
  
  updateFollowingBaseSelection(detailView, FALSE, FALSE);
}


/* Set the selected base index to the given display index and base/frame number.
 *  Performs any required refreshes. Scrolls the view to keep the selected base 
 * in view if allowScroll is true. (Such scrolling is by the minimum
 * number of bases necessary if scrollMinimum is true.) */
void detailViewSetSelectedBaseIdx(GtkWidget *detailView, 
                                  const int selectedBaseIdx, 
                                  const int frame, 
                                  const int baseNum, 
                                  const gboolean allowScroll,
                                  const gboolean scrollMinimum)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  properties->selectedBaseSet = TRUE;
  properties->selectedBaseIdx = selectedBaseIdx;
  properties->selectedFrame = frame;
  properties->selectedBaseNum = baseNum;
  
  /* For protein matches, calculate the base index in terms of the DNA sequence and cache it */
  BlxViewContext *bc = detailViewGetContext(detailView);
  properties->selectedDnaBaseIdx = convertDisplayIdxToDnaIdx(selectedBaseIdx, bc->seqType, frame, baseNum, bc->numFrames, bc->displayRev, &bc->refSeqRange);

  updateFollowingBaseSelection(detailView, allowScroll, scrollMinimum);
}


int detailViewGetClickedBaseIdx(GtkWidget *detailView)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  return properties->clickedBaseIdx;
}


/* Set the selected base index to the given display index and base/frame number.
 *  Performs any required refreshes. Scrolls the view to keep the selected base 
 * in view if allowScroll is true. (Such scrolling is by the minimum
 * number of bases necessary if scrollMinimum is true.) */
void detailViewSetClickedBaseIdx(GtkWidget *detailView, 
                                  const int clickedBaseIdx, 
                                  const int frame, 
                                  const int baseNum, 
                                  const gboolean allowScroll,
                                  const gboolean scrollMinimum)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  /* For protein matches, calculate the base index in terms of the DNA sequence and cache it */
  BlxViewContext *bc = detailViewGetContext(detailView);
  properties->clickedBaseIdx = convertDisplayIdxToDnaIdx(clickedBaseIdx, bc->seqType, frame, baseNum, bc->numFrames, bc->displayRev, &bc->refSeqRange);
}


static GtkWidget *detailViewGetBigPicture(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetBigPicture(blxWindow);
}


DetailViewProperties* detailViewGetProperties(GtkWidget *widget)
{
  /* optimisation: cache result, because we know there is only ever one detail view */
  static DetailViewProperties *properties = NULL;

  if (!properties && widget)
    properties = (DetailViewProperties*)(g_object_get_data(G_OBJECT(widget), "DetailViewProperties"));
  
  return properties;
}

static BlxViewContext* detailViewGetContext(GtkWidget *detailView)
{
  GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
  return blxWindowGetContext(blxWindow);
}

static void onDestroyDetailView(GtkWidget *widget)
{
  DetailViewProperties *properties = detailViewGetProperties(widget);

  /* N.B. Don't free the cell renderer, or it causes memory corruption. I'm not
   * sure what owns it - the columns it is added to? */
    
  if (properties->fwdStrandTrees)
    {
      g_list_free(properties->fwdStrandTrees);
      properties->fwdStrandTrees = NULL;
    }
  
  if (properties->revStrandTrees)
    {
      g_list_free(properties->revStrandTrees);
      properties->revStrandTrees = NULL;
    }

  if (properties->fontDesc)
    {
      pango_font_description_free(properties->fontDesc);
      properties->fontDesc = NULL;
    }
  
  if (properties->spliceSites)
    {
      g_slist_foreach(properties->spliceSites, destroyBlxSpliceSite, NULL);
      g_slist_free(properties->spliceSites);
      properties->spliceSites = NULL;
    }
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "DetailViewProperties", NULL);
    }
}


static void detailViewCreateProperties(GtkWidget *detailView,
                                       GtkWidget *blxWindow,
                                       GtkCellRenderer *renderer,
                                       GList *fwdStrandTrees,
                                       GList *revStrandTrees,
                                       GtkWidget *header,
                                       GtkWidget *feedbackBox,
                                       GtkWidget *statusBar,
                                       GList *columnList,
                                       GtkAdjustment *adjustment, 
                                       const int startCoord,
                                       const BlxColumnId sortColumn)
{
  if (detailView)
    { 
      DetailViewProperties *properties = (DetailViewProperties*)g_malloc(sizeof *properties);

      /* Find a fixed-width font */
      const char *fontFamily = findFixedWidthFont(detailView);
      PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(fontDesc, pango_font_description_get_size(detailView->style->font_desc));
      
      properties->blxWindow = blxWindow;
      properties->renderer = renderer;
      properties->fwdStrandTrees = fwdStrandTrees;
      properties->revStrandTrees = revStrandTrees;
      properties->header = header;
      properties->feedbackBox = feedbackBox;
      properties->statusBar = statusBar;
      properties->adjustment = adjustment;
      properties->selectedBaseSet = FALSE;
      properties->selectedBaseIdx = UNSET_INT;
      properties->selectedBaseNum = UNSET_INT;
      properties->selectedFrame = UNSET_INT;
      properties->selectedDnaBaseIdx = UNSET_INT;
      properties->fontDesc = fontDesc;
      properties->charWidth = 0.0;
      properties->charHeight = 0.0;
      properties->snpConnectorHeight = DEFAULT_SNP_CONNECTOR_HEIGHT;
      properties->numUnalignedBases = DEFAULT_NUM_UNALIGNED_BASES;
      
      /* The numunalignedbases may be set in the config file; if so, override the default */
      GKeyFile *key_file = blxGetConfig();
      if (key_file)
        {
          GError *error = NULL;
          int numUnaligned = g_key_file_get_integer(key_file, SETTINGS_GROUP, SETTING_NAME_NUM_UNALIGNED_BASES, &error);
          
          if (!error) /* we don't care if it wasn't found */
            properties->numUnalignedBases = numUnaligned;
        }

      /* Add the splice sites that we want Blixem to identify as canonical */
      properties->spliceSites = NULL;
      addBlxSpliceSite(&properties->spliceSites, "GT", "AG", FALSE);
      addBlxSpliceSite(&properties->spliceSites, "GC", "AG", FALSE);
      addBlxSpliceSite(&properties->spliceSites, "AT", "AC", TRUE);

      /* We don't know the display range yet, so set an arbitrary range centred
       * on the start coord. Set the adjustment value to be unset so that we know
       * we need to calculated it first time round. */
      properties->adjustment->value = UNSET_INT;
      properties->displayRange.min = startCoord - 1;
      properties->displayRange.max = startCoord + 1;
      
      /* Find the padding between the background area of the tree cells and the actual
       * drawing area. This is used to render the full height of the background area, so
       * that we don't have gaps between rows. */ 
      if (fwdStrandTrees)
        {
          GtkWidget *widget = GTK_WIDGET(fwdStrandTrees->data);

          GtkWidget *tree = widgetIsTree(widget) ? widget : treeContainerGetTree(GTK_CONTAINER(widget));
          gtk_widget_realize(tree); /* must realize tree to pick up any overriden style properties */
            
          /* I can't get this to work properly using gtk_tree_view_get_cell_area etc. The cell
           * area and background area come out wrong - perhaps because I'm not using a real
           * row path. After a bit of experimentation it seems that the y padding is always
           * related to the vertical-separator as follows: ypad = trunc(vseparator / 2) + 1 */
          gint vertical_separator, horizontal_separator;
          gtk_widget_style_get (tree, "vertical-separator", &vertical_separator, NULL);
          gtk_widget_style_get (tree, "horizontal-separator", &horizontal_separator, NULL);
          
          properties->cellXPadding = (horizontal_separator / 2) + 1;
          properties->cellYPadding = (vertical_separator / 2) + 1;
        }
      properties->exonBoundaryLineWidth   = 1;
      properties->exonBoundaryLineStyle = GDK_LINE_SOLID;
      properties->exonBoundaryLineStylePartial = GDK_LINE_ON_OFF_DASH;
      
      /* Allocate sortColumns array to be same length as columnList */
      const int numColumns = g_list_length(columnList);
      properties->sortColumns = (BlxColumnId*)g_malloc(numColumns * sizeof(BlxColumnId));
      
      int i = 0;
      for ( ; i < numColumns; ++i)
        properties->sortColumns[i] = BLXCOL_NONE;

      /* Sort by the default/input column, and then by name and then position */
      properties->sortColumns[0] = sortColumn;
      properties->sortColumns[1] = BLXCOL_SEQNAME;
      properties->sortColumns[2] = BLXCOL_START;
      
      g_object_set_data(G_OBJECT(detailView), "DetailViewProperties", properties);
      g_signal_connect(G_OBJECT(detailView), "destroy", G_CALLBACK(onDestroyDetailView), NULL); 
    }
}


/* Get/set functions for the sequence column header's row number. (There is only
 * one property to set, so it's not worth having a struct for the properties.) */
void seqColHeaderSetRow(GtkWidget *header, const int frame)
{
  g_object_set_data(G_OBJECT(header), "seqColHeaderFrameNumber", GINT_TO_POINTER(frame));
}

int seqColHeaderGetRow(GtkWidget *header)
{
  return header ? GPOINTER_TO_INT(g_object_get_data(G_OBJECT(header), "seqColHeaderFrameNumber")) : UNSET_INT;
}

/* Gets the base number within the given reading frame that this header row corresponds to*/
int seqColHeaderGetBase(GtkWidget *header, const int frame, const int numFrames)
{
  const int row = seqColHeaderGetRow(header);
  
  int baseNum = UNSET_INT;
  
  if (row != UNSET_INT)
    {
      baseNum = row - (frame - 1);

      /* Cyclic decrement */
      if (baseNum < 1)
        {
          baseNum += numFrames;
        }
    }
  
  return baseNum;
}


/* Get/set functions for the SNP track header strand. (There is only one property
 * to set, so it's not worth having a struct for the properties.) */
static void snpTrackSetStrand(GtkWidget *snpTrack, const BlxStrand strand)
{
  g_object_set_data(G_OBJECT(snpTrack), "snpTrackStrand", GINT_TO_POINTER(strand));
}

static BlxStrand snpTrackGetStrand(GtkWidget *snpTrack, GtkWidget *detailView)
{
  BlxStrand strand = BLXSTRAND_NONE;
  
  if (snpTrack)
    {
      strand = (BlxStrand)GPOINTER_TO_INT(g_object_get_data(G_OBJECT(snpTrack), "snpTrackStrand"));
    }
  
  if (strand == BLXSTRAND_NONE)
    {
      /* use the active strand instead */
      strand = blxWindowGetActiveStrand(detailViewGetBlxWindow(detailView));
    }
  
  return strand;
}


/***********************************************************
 *                     Events                              *
 ***********************************************************/

/* Expose function for generic text headers. Draws vertical separator lines and then lets the
 * default handler continue.  */
gboolean onExposeGenericHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkGC *gc = gdk_gc_new(headerWidget->window);
  
  GtkWidget *detailView = GTK_WIDGET(data);
  const BlxViewContext *bc = detailViewGetContext(detailView);
  
  drawColumnSeparatorLine(headerWidget, headerWidget->window, gc, bc);
  
  g_object_unref(gc);
  
  return FALSE;
}


/* expose function to push a cached bitmap to screen */
gboolean onExposeDnaTrack(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkWindow *window = GTK_IS_LAYOUT(headerWidget) ? GTK_LAYOUT(headerWidget)->bin_window : headerWidget->window;
  
  if (window)
    {
      /* See if there's a cached drawable and, if not, create it */
      GdkDrawable *bitmap = widgetGetDrawable(headerWidget);

      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          GtkWidget *detailView = GTK_WIDGET(data);
          const BlxStrand activeStrand = detailViewGetDisplayRev(detailView) ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
          
          drawDnaTrack(headerWidget, detailView, activeStrand, seqColHeaderGetRow(headerWidget));
          bitmap = widgetGetDrawable(headerWidget);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
    }
  
  return FALSE;
}


/* Expose handler for the variations-track header */
static gboolean onExposeVariationsTrack(GtkWidget *snpTrack, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(snpTrack)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(snpTrack);
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          GtkWidget *detailView = GTK_WIDGET(data);
          drawVariationsTrack(snpTrack, detailView);
          
          bitmap = widgetGetDrawable(snpTrack);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
      else
        {
          g_critical("Failed to draw SNP track [%p] - could not create bitmap.\n", snpTrack);
        }
    }
  
  return TRUE;
}


/* Handler for when the mouse button is pressed in the SNP track header */
static gboolean onButtonPressSnpTrack(GtkWidget *snpTrack, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  switch (event->button)
  {
    case 1:
      {
        GtkWidget *detailView = GTK_WIDGET(data);
        GtkWidget *blxWindow = detailViewGetBlxWindow(detailView);
        GList *columnList = detailViewGetColumnList(detailView);
        
        if (event->type == GDK_BUTTON_PRESS) /* first click */
          {
            /* Select the variation that was clicked on.  */
            blxWindowDeselectAllSeqs(blxWindow);

            /* The SNP track is not the same width as the sequence column, so pass the
             * sequence column header so that we can convert to the correct coords */
            BlxColumnInfo *seqColInfo = getColumnInfo(columnList, BLXCOL_SEQUENCE);

            selectClickedSnp(snpTrack, seqColInfo->headerWidget, detailView, event->x, event->y, TRUE, UNSET_INT); /* SNPs are always expanded in the SNP track */
            
            detailViewRefreshAllHeaders(detailView);
          }      
        else if (event->type == GDK_2BUTTON_PRESS) /* double-click */
          {
            /* If a variation was double-clicked, open its URL in a browser. 
             * If multiple are selected, use the last-selected one. */
            GList *seqItem = g_list_last(blxWindowGetSelectedSeqs(blxWindow));
            
            if (seqItem)
              {
                BlxSequence *seq = (BlxSequence*)(seqItem->data);
                fetchSequence(seq, TRUE, 0, blxWindow, NULL, NULL);
              }
          }

        handled = TRUE;
        break;
      }
  }
  
  return handled;
}


/* Parent handler for button presses anywhere on the detail view. */
static gboolean onButtonPressDetailView(GtkWidget *detailView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
  {
    case 2:
    {
      /* Middle button: select the base index at the clicked coords */
      
      /* First make sure the detail view trees are fully drawn (so that the
       * cached drawable is up to date - see notes in onExposeDetailViewTree)
       * and then set the mouse-drag flag. */
      detailViewRedrawAll(detailView);
      gdk_window_process_all_updates();
      setMouseDragMode(TRUE);
      
      int baseIdx = getBaseIndexAtDetailViewCoords(detailView, event->x, event->y);
      
      if (baseIdx != UNSET_INT)
        {
          /* For protein matches, select the first base in the triplet */
          const int baseNum = 1;
          const int frame = detailViewGetActiveFrame(detailView);
          detailViewSetSelectedBaseIdx(detailView, baseIdx, frame, baseNum, FALSE, TRUE);
        }
      
      handled = TRUE;
      break;
    }
      
    case 3:
    {
      /* Right button: store the clicked coord (required for copying ref seq),
       * and then show the context menu. */
      int baseIdx = getBaseIndexAtDetailViewCoords(detailView, event->x, event->y);
      if (baseIdx != UNSET_INT)
        {
          const int baseNum = 1;
          const int frame = detailViewGetActiveFrame(detailView);
          detailViewSetClickedBaseIdx(detailView, baseIdx, frame, baseNum, FALSE, TRUE);
        }
      

      GtkWidget *mainmenu = blxWindowGetMainMenu(detailViewGetBlxWindow(detailView));
      gtk_menu_popup (GTK_MENU(mainmenu), NULL, NULL, NULL, NULL, event->button, event->time);
      handled = TRUE;
      break;
    }
      
    default:
      break;
  };
  
  return handled;
}


static gboolean onMouseMoveDetailView(GtkWidget *detailView, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base
       * (but don't re-centre on the selected base until the mouse button is released). */
      const int baseIdx = getBaseIndexAtDetailViewCoords(detailView, event->x, event->y);
      
      if (baseIdx != UNSET_INT)
        {
          /* For protein matches, get the 1st base in the peptide */
          const int baseNum = 1;
          const int frame = detailViewGetActiveFrame(detailView);
          detailViewSetSelectedBaseIdx(detailView, baseIdx, frame, baseNum, FALSE, TRUE);
        }
      
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseDetailView(GtkWidget *detailView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Right button: show context menu (handled in button-press event) */
  if (event->button == 3)
    {
      handled = TRUE;
    }
  
  /* Middle button: scroll the selected base index to the centre (unless CTRL is pressed) */
  if (event->button == 2)
    {
      /* Cancel middle-drag mode */
      setMouseDragMode(FALSE);
      
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
      
      if (!ctrlModifier)
        {
          /* Move the scrollbar so that the currently-selected base index is at the centre */
          const int selectedBaseIdx = detailViewGetSelectedBaseIdx(detailView);
          
          if (selectedBaseIdx != UNSET_INT)
            {
              /* The coord is in terms of the display coords, i.e. whatever the displayed seq type is. */
              const BlxSeqType seqType = detailViewGetSeqType(detailView);
              const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
              
              int newStart = selectedBaseIdx - (getRangeLength(displayRange) / 2);
              setDetailViewStartIdx(detailView, newStart, seqType);
            }
        }

      detailViewRedrawAll(detailView);
      
      handled = TRUE;
    }
  
  return handled;
}


/* Handler for when a mouse button is pressed on a sequence column header widget.
 * This signal only gets connected for protein matches, where there are 3 header
 * widgets, each showing the DNA nucleotides for a specific reading frame. */
static gboolean onButtonPressSeqColHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
  {
    case 1:
      {
        GtkWidget *detailView = GTK_WIDGET(data);

        if (event->type == GDK_BUTTON_PRESS)
          {
            /* Select the SNP that was clicked on.  */
            blxWindowDeselectAllSeqs(detailViewGetBlxWindow(detailView));
            const int clickedBase = seqColHeaderGetBase(header, 1, detailViewGetNumFrames(detailView));

            selectClickedSnp(header, NULL, detailView, event->x, event->y, FALSE, clickedBase); /* SNPs are always un-expanded in the DNA track */
            
            detailViewRefreshAllHeaders(detailView);
          }
        else if (event->type == GDK_2BUTTON_PRESS)
          {
            /* Double-click: toggle variations-track visibility */
            BlxViewContext *bc = detailViewGetContext(detailView);
            const gboolean showTrack = !bc->flags[BLXFLAG_SHOW_VARIATION_TRACK];
            
            /* If we're enabling the variations track, also make sure that variations are highlighted in the header */
            if (showTrack)
              bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS] = TRUE;
            
            bc->flags[BLXFLAG_SHOW_VARIATION_TRACK] = showTrack;
            detailViewUpdateShowSnpTrack(detailView, showTrack);
          }

        handled = TRUE;
        break;
      }
      
    case 2:
      {
        /* Middle button: select the nucleotide that was clicked on. */
        GtkWidget *detailView = GTK_WIDGET(data);
        selectClickedNucleotide(header, detailView, event->x, event->y);
        handled = TRUE;
        break;
      }
      
    default:
      break;
  }
  
  return handled;
}


static gboolean onButtonReleaseSeqColHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* Middle button: scroll the selected base index to the centre (unless CTRL is pressed) */
  if (event->button == 2)
    {
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
      
      if (!ctrlModifier)
        {
          /* Move the scrollbar so that the currently-selected base index is at the centre */
          GtkWidget *detailView = GTK_WIDGET(data);
          detailViewCentreOnSelection(detailView);
        }
      
      handled = TRUE;
    }

  return handled;
}


static gboolean onMouseMoveSeqColHeader(GtkWidget *header, GdkEventMotion *event, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);

  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Moving mouse with middle mouse button down. Update the currently-selected base */
      selectClickedNucleotide(header, detailView, event->x, event->y);
    }
  else
    {
      /* Get the nucelotide we're hovering over and feedback info about that nucelotide in the
       * feedback area */
      int baseIdx, frame, baseNum;
      getSeqColHeaderClickedNucleotide(header, detailView, event->x, event->y, &baseIdx, &frame, &baseNum);
      
      if (baseIdx != UNSET_INT && frame != UNSET_INT && baseNum != UNSET_INT)
        {
          BlxViewContext *bc = detailViewGetContext(detailView);
          const int dnaIdx = convertDisplayIdxToDnaIdx(baseIdx, bc->seqType, frame, baseNum, bc->numFrames, bc->displayRev, &bc->refSeqRange);
          const BlxStrand strand = blxWindowGetActiveStrand(detailViewGetBlxWindow(detailView));
          
          updateFeedbackAreaNucleotide(detailView, dnaIdx, strand);
        }
    }
  
  return TRUE;
}


/* Updates the widths of any dynamic-sized columns (really they are fixed-width
 * columns so we update their widths manually here). */
void updateDynamicColumnWidths(GtkWidget *detailView)
{
  /* Currently, only the sequence column has dynamic width, so just sum all of the
   * other (visible) column widths and subtract from the allocation width. */
  int width = detailView->allocation.width;
  BlxColumnInfo *seqColInfo = NULL;
  
  GList *listItem = detailViewGetColumnList(detailView);
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);
      
      if (columnInfo->columnId == BLXCOL_SEQUENCE)
        seqColInfo = columnInfo;
      else if (columnInfo->showColumn && columnInfo->dataLoaded)
        width -= columnInfo->width;
    }
      
  if (seqColInfo && seqColInfo->width != width)
    {
      seqColInfo->width = width;
  
      callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns, NULL);
      resizeDetailViewHeaders(detailView);

      updateDetailViewRange(detailView);
    }
}


/***********************************************************
 *                     Callbacks                           *
 ***********************************************************/

/* If there are two trees and one is visible and the other hidden, toggle their hidden states */
static void swapTreeVisibility(GtkWidget *detailView)
{
  /* Only do anything for DNA matches, where we have 1 frame and both strands are visible. */
  if (detailViewGetSeqType(detailView) == BLXSEQ_DNA)
    {
      GtkWidget *tree1 = detailViewGetTreeContainer(detailView, BLXSTRAND_FORWARD, 1);
      GtkWidget *tree2 = detailViewGetTreeContainer(detailView, BLXSTRAND_REVERSE, 1);
      
      if (widgetGetHidden(tree1) != widgetGetHidden(tree2))
        {
          widgetSetHidden(tree1, !widgetGetHidden(tree1));
          widgetSetHidden(tree2, !widgetGetHidden(tree2));
        }
    }
}


/* If one grid is visible and the other hidden, toggle their hidden states */
static void swapGridVisibility(GtkWidget *bigPicture)
{
  GtkWidget *grid1 = bigPictureGetFwdGrid(bigPicture);
  GtkWidget *grid2 = bigPictureGetRevGrid(bigPicture);
  
  if (widgetGetHidden(grid1) != widgetGetHidden(grid2))
    {
      widgetSetHidden(grid1, !widgetGetHidden(grid1));
      widgetSetHidden(grid2, !widgetGetHidden(grid2));
    }
}


/* If one exon view is visible and the other hidden, toggle their hidden states */
static void swapExonViewVisibility(GtkWidget *bigPicture)
{
  GtkWidget *exonView1 = bigPictureGetFwdExonView(bigPicture);
  GtkWidget *exonView2 = bigPictureGetRevExonView(bigPicture);
  
  if (widgetGetHidden(exonView1) != widgetGetHidden(exonView2))
    {
      widgetSetHidden(exonView1, !widgetGetHidden(exonView1));
      widgetSetHidden(exonView2, !widgetGetHidden(exonView2));
    }
}


void toggleStrand(GtkWidget *detailView)
{
  DEBUG_ENTER("toggleStrand");

  BlxViewContext *blxContext = detailViewGetContext(detailView);
  GtkWidget *bigPicture = detailViewGetBigPicture(detailView);

  /* Update the flag */
  blxContext->displayRev = !blxContext->displayRev;
  
  /* Invert the display range for both detail view and big picture */
  IntRange *bpRange = bigPictureGetDisplayRange(bigPicture);
  const IntRange* const fullRange = &blxContext->fullDisplayRange;
  const int bpStart = fullRange->max - bpRange->min + fullRange->min;
  const int bpEnd = fullRange->max - bpRange->max + fullRange->min;
  intrangeSetValues(bpRange, bpStart, bpEnd);

  IntRange *displayRange = detailViewGetDisplayRange(detailView);
  const int newStart = fullRange->max - displayRange->max + fullRange->min;
  setDetailViewStartIdx(detailView, newStart, blxContext->seqType);

  /* Re-select the currently-selected index, if there is one, because the display coords
   * have changed. */
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  if (properties->selectedBaseSet)
    {
      const int activeFrame = detailViewGetActiveFrame(detailView); 
      detailViewSetSelectedDnaBaseIdx(detailView, properties->selectedDnaBaseIdx, activeFrame, FALSE, TRUE);
    }
  
  /* Re-calculate the cached display ranges for the MSPs */
  cacheMspDisplayRanges(blxContext, properties->numUnalignedBases);
  
  /* If one grid/tree is hidden and the other visible, toggle which is hidden */
  swapTreeVisibility(detailView);
  swapGridVisibility(bigPicture);
  swapExonViewVisibility(bigPicture);
  
  /* Toggle the order of the trees and grids. */
  refreshTreeOrder(detailView);
  refreshGridOrder(bigPicture);

  /* Redraw the tree and detail-view headers. Also need to resort because if the display is
   * reversed it affects sorting by position. */
  refreshDetailViewHeaders(detailView);
  detailViewResortTrees(detailView);
  callFuncOnAllDetailViewTrees(detailView, refreshTreeHeaders, NULL);
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  detailViewRedrawAll(detailView);
  
  /* Redraw the grids and grid headers */
  refreshBigPictureDisplayRange(bigPicture, FALSE);
  
  DEBUG_EXIT("toggleStrand returning");
}


/* Go to a user-specified coord on the reference sequence (in terms of the DNA
 * or peptide sequence as specified by coordSeqType). Selects the base index
 * and centres on it. */
void goToDetailViewCoord(GtkWidget *detailView, const BlxSeqType coordSeqType)
{
  static gchar defaultInput[32] = "";
  static gboolean toplevelCoords = FALSE;

  BlxViewContext *bc = detailViewGetContext(detailView);
  
  /* Pop up a dialog to request a coord from the user */
  char *title = g_strdup_printf("%sGo to position", blxGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title,
                                                  NULL, 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_NO_SEPARATOR),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;
  const int padding = 4;

  /* Create a text-entry widget for the user to enter the coord */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(contentArea), hbox, FALSE, FALSE, padding);

  GtkWidget *label = gtk_label_new("Enter coord: ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, padding);

  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, padding);

  /* Create a tick-box to toggle to top-level coords (default is local coords) */
  GtkWidget *checkBox = gtk_check_button_new_with_mnemonic("Use _top-level coords");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkBox), toplevelCoords);
  gtk_box_pack_start(GTK_BOX(contentArea), checkBox, FALSE, FALSE, padding);

  /* Show the dialog and wait for a response */
  gtk_widget_show_all(dialog);

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const gchar *inputText = gtk_entry_get_text(GTK_ENTRY(entry));
      toplevelCoords = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkBox));

      /* Convert the input string to an int */
      int requestedCoord = atoi(inputText);
      
      if (requestedCoord)
        {
          /* Remember this input for next time */
          sprintf(defaultInput, "%d", requestedCoord);
          
          int coord = requestedCoord;
          
          /* If display coords are negated, assume the user has entered a 
           * negative coord too, and un-negate it. */
          const gboolean negate = bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS];
          
          if (negate)
            {
              coord *= -1;
            }

          /* If the user entered a top-level coord, apply the offset to convert to local coords */
          if (toplevelCoords)
            {
              coord += bc->refSeqOffset;
            }

          /* Check if it's in range. If not, try negating it in case the user
           * entered the wrong sign */
          if (!valueWithinRange(coord, &bc->refSeqRange))
            {
              g_debug("Coord '%d' does not lie within the reference sequence range [%d %d]; trying '%d'\n", requestedCoord, negate ? bc->refSeqRange.max * -1 : bc->refSeqRange.min, negate ? bc->refSeqRange.min * -1 : bc->refSeqRange.max, requestedCoord * -1);
              coord *= -1;
            }
          
          if (valueWithinRange(coord, &bc->refSeqRange))
            {
              const int activeFrame = detailViewGetActiveFrame(detailView);
              detailViewSetSelectedDnaBaseIdx(detailView, coord, activeFrame, TRUE, FALSE);
              detailViewRedrawAll(detailView);
            }
          else
            {
              g_debug("Coord '%d' does not lie within the reference sequence range [%d %d]\n", negate ? coord * -1 : coord, negate? bc->refSeqRange.max * -1 : bc->refSeqRange.min, negate ? bc->refSeqRange.min * -1 : bc->refSeqRange.max);
              g_critical("Coord was not in the reference sequence range [%d %d]\n", negate ? bc->refSeqRange.max * -1 : bc->refSeqRange.min, negate ? bc->refSeqRange.min * -1 : bc->refSeqRange.max);
            }
        }
    }

  gtk_widget_destroy(dialog);
}


/* This sets the main sort column to be the given column (it does not change
 * any secondary sort columns). */
void detailViewSetSortColumn(GtkWidget *detailView, const BlxColumnId sortColumn)
{
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(properties->blxWindow);
  const int numColumns = g_list_length(columnList);

  if (numColumns > 0)
    {      
      if (properties->sortColumns[0] != sortColumn)
        {
          properties->sortColumns[0] = sortColumn;
          detailViewResortTrees(detailView);
        }
    }
}


/* Find the next MSP (out the MSPs in this tree) whose start/end is the next closest to the
 * current start position of the display, searching only in the direction specified by the
 * search criteria passed in the user data. Updates the searchData with the offset found. */
static gboolean findNextMatchInTree(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  MatchSearchData *searchData = (MatchSearchData*)data;

  /* Loop through all MSPs in this tree row. */
  GList* mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      /* Check its in the sequence list (if given), is a valid match or exon, and is in a visible layer */
      if (mspLayerIsVisible(msp) &&
          (mspIsExon(msp) || mspIsBlastMatch(msp)) &&
          (!searchData->seqList || g_list_find(searchData->seqList, msp->sSequence)))
        {
          /* Get the offset of the msp coords from the given start coord and find the smallest,
           * ignorning zero and negative offsets (negative means its the wrong direction) */
          const int offset1 = (msp->qRange.min - searchData->startDnaIdx) * searchData->searchDirection;
          const int offset2 = (msp->qRange.max - searchData->startDnaIdx) * searchData->searchDirection;
          
          int currentFrame = mspGetRefFrame(msp, searchData->seqType);
          int currentBest = UNSET_INT;
          
          if (offset1 > 0 && (offset2 <= 0 || offset2 >= offset1))
            {
              currentBest = offset1;
            }
          else
            {
              currentBest = offset2;
            }
          if (currentBest > 0)
            {
              gboolean useNew = FALSE;

              /* Check if this is the first one we've looked at */
              useNew |= (searchData->offset == UNSET_INT);
              
              /* Check if this offset is smaller than the previous best */
              useNew |= (currentBest < searchData->offset);
              
              /* If it's the same as the previous best, use the one with the smallest 
               * frame number, i.e. give preference to the topmost frame (in the hope
               * that this will be the least confusing!) */
              useNew |= (currentBest == searchData->offset && currentFrame < searchData->foundFrame);
              
              if (useNew)
                {
                  searchData->offset = currentBest;
                  searchData->foundFrame = currentFrame;
                  searchData->foundMsp = msp;
                }
            }
        }
    }
  
  return FALSE;
}


/* Find and go to the next match (either left or right depending on the search flag).
 * If a list of sequence names is given, only look at matches in those sequences.
 * startDnaIdx determines where to start searching from. Sets the found match idx as
 * the currently-selected base index */
static MSP* goToNextMatch(GtkWidget *detailView, const int startDnaIdx, const gboolean searchRight, GList *seqList)
{
  BlxViewContext *bc = detailViewGetContext(detailView);
  
  const int searchDirection = (searchRight != bc->displayRev) ? 1 : -1;

  MatchSearchData searchData = {startDnaIdx, 
                                searchRight, 
                                searchDirection, 
                                bc->displayRev, 
                                bc->seqType,
                                bc->numFrames,
                                &bc->refSeqRange,
                                seqList,
                                UNSET_INT,
                                UNSET_INT,
                                NULL};

  /* Loop through the MSPs in all visible trees */
  int frame = 1;

  for (  ; frame <= searchData.numFrames; ++frame)
    {
      GtkWidget *treeContainer = detailViewGetTreeContainer(detailView, BLXSTRAND_FORWARD, frame);
      GtkWidget *tree = treeContainerGetTree(GTK_CONTAINER(treeContainer));
      
      if (GTK_WIDGET_VISIBLE(tree) && gtk_widget_get_parent(treeContainer)) /* ignore if not currently included in view */
        {
          GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
          gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
        }

      treeContainer = detailViewGetTreeContainer(detailView, BLXSTRAND_REVERSE, frame);
      tree = treeContainerGetTree(GTK_CONTAINER(treeContainer));
      
      if (GTK_WIDGET_VISIBLE(tree) && gtk_widget_get_parent(treeContainer)) /* ignore if not currently included in view */
        {
          GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
          gtk_tree_model_foreach(model, findNextMatchInTree, &searchData);
        }
    }
  
  if (searchData.offset != UNSET_INT)
    {
      /* Offset the start coord by the found amount. */
      int newDnaIdx = searchData.startDnaIdx + (searchData.offset * searchDirection);
      detailViewSetSelectedDnaBaseIdx(detailView, newDnaIdx, searchData.foundFrame, TRUE, FALSE);
      callFuncOnAllDetailViewTrees(detailView, treeScrollSelectionIntoView, NULL);
      detailViewRedrawAll(detailView);
    }
  
  return searchData.foundMsp;
}


/* Go to the previous match (optionally limited to matches in the given list)  */
MSP* prevMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startDnaIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
  
  if (!valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
      
      /* Use base 1 within the currently selected frame for this display coord */
      int frame = detailViewGetActiveFrame(detailView);
      BlxViewContext *bc = detailViewGetContext(detailView);
      startDnaIdx = convertDisplayIdxToDnaIdx(startCoord, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  return goToNextMatch(detailView, startDnaIdx, FALSE, seqList);
}


/* Go to the next match (optionally limited to matches in the given list)  */
MSP* nextMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the currently selected base index, if there is
   * one and if it is currently visible. Otherwise use the current display centre. */
  int startDnaIdx = detailViewGetSelectedDnaBaseIdx(detailView);
  int startCoord = detailViewGetSelectedBaseIdx(detailView);
  const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
  
  if (!valueWithinRange(startCoord, displayRange))
    {
      startCoord = getRangeCentre(displayRange);
      
      /* Use base 1 within the currently selected frame for this display coord */
      int frame = detailViewGetActiveFrame(detailView);
      BlxViewContext *bc = detailViewGetContext(detailView);
      startDnaIdx = convertDisplayIdxToDnaIdx(startCoord, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  return goToNextMatch(detailView, startDnaIdx, TRUE, seqList);
}


/* Go to the first match (optionally limited to matches in the given list)  */
MSP* firstMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the start of the ref seq */
  const IntRange* const refSeqRange = detailViewGetRefSeqRange(detailView);
  const int startIdx = detailViewGetDisplayRev(detailView) ? refSeqRange->max : refSeqRange->min;
  
  return goToNextMatch(detailView, startIdx, TRUE, seqList);
}


/* Go to the last match (optionally limited to matches in the given list)  */
MSP* lastMatch(GtkWidget *detailView, GList *seqList)
{
  /* Jump to the nearest match to the end of the reference sequence */
  const IntRange* const refSeqRange = detailViewGetRefSeqRange(detailView);
  const int startIdx = detailViewGetDisplayRev(detailView) ? refSeqRange->min : refSeqRange->max;

  return goToNextMatch(detailView, startIdx, FALSE, seqList);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* This loops through all the columns and adds each column's header widget to the header bar. */
static void addColumnsToHeaderBar(GtkBox *headerBar, GList *columnList)
{
  /* Loop through each column and add its header widget to the header bar */
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(columnItem->data);

      /* The sequence column is a special one that wants to fill any additional space, so
       * set the 'expand' property to true for that column only. */
      const gboolean expand = (columnInfo->columnId == BLXCOL_SEQUENCE );
      gtk_box_pack_start(headerBar, columnInfo->headerWidget, expand, TRUE, 0);
    }
}


/* Create the header bar for the detail view. This contains the labels for the
 * detail-view trees (since we only want one label at the top, rather than 
 * individual labels for each tree). For protein sequence matches, the header
 * for the sequence column will also show the DNA sequence (separated into reading
 * frames). 
 * Column data is compiled into the detailViewColumns return argument. */
static GtkWidget* createDetailViewHeader(GtkWidget *detailView, 
                                         const BlxSeqType seqType, 
                                         const int numFrames,
                                         GList *columnList,
                                         const gboolean includeSnpTrack)
{
  GtkBox *header = NULL; /* outermost container for the header */
  GtkBox *headerBar = GTK_BOX(gtk_hbox_new(FALSE, 0)); /* hbox for the column headers */
  
  /* Create a SNP track, if requested */
  if (includeSnpTrack)
    {
      header = GTK_BOX(gtk_vbox_new(FALSE, 0));
      createSnpTrackHeader(GTK_BOX(header), detailView, BLXSTRAND_NONE);
      gtk_box_pack_start(header, GTK_WIDGET(headerBar), FALSE, FALSE, 0);
      gtk_widget_set_name(GTK_WIDGET(headerBar), HEADER_CONTAINER_NAME);
    }
  else
    {
      /* Only need the column header bar */
      header = headerBar;
    }

  /* Add all the column headers to the header bar */
  addColumnsToHeaderBar(headerBar, columnList);

  return GTK_WIDGET(header);
}


/* Create the variations track header widget. If the strand is BLXSTRAND_NONE
 * then the active strand will be used. */
GtkWidget* createSnpTrackHeader(GtkBox *parent, GtkWidget *detailView, const BlxStrand strand)
{
  GtkWidget *snpTrack = gtk_layout_new(NULL, NULL);

  gtk_widget_set_name(snpTrack, SNP_TRACK_HEADER_NAME);  
  gtk_box_pack_start(parent, snpTrack, FALSE, TRUE, 0);
  snpTrackSetStrand(snpTrack, strand);

  gtk_widget_add_events(snpTrack, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(snpTrack), "expose-event", G_CALLBACK(onExposeVariationsTrack), detailView);
  g_signal_connect(G_OBJECT(snpTrack), "button-press-event", G_CALLBACK(onButtonPressSnpTrack), detailView);
  
  return snpTrack;
}


/* Create a custom header widget for the sequence column for protein matches (this is where
 * we will display the triplets that make up the codons.) Returns NULL for DNA matches. */
static void createSeqColHeader(GtkWidget *detailView,
                               const BlxSeqType seqType,
                               const int numFrames,
                               GList *columnList)
{
  if (seqType == BLXSEQ_PEPTIDE)
    {
      GtkWidget *header = gtk_vbox_new(FALSE, 0);
      gtk_widget_set_name(header, HEADER_CONTAINER_NAME);
      
      int frame = 0;
      for ( ; frame < numFrames; ++frame)
        {
          GtkWidget *line = gtk_layout_new(NULL, NULL);
          gtk_box_pack_start(GTK_BOX(header), line, FALSE, TRUE, 0);

          gtk_widget_set_name(line, DNA_TRACK_HEADER_NAME);
          seqColHeaderSetRow(line, frame + 1);

          /* Disable double buffering to try to speed up drawing */
          GTK_WIDGET_UNSET_FLAGS (line, GTK_DOUBLE_BUFFERED);

          gtk_widget_add_events(line, GDK_BUTTON_PRESS_MASK);
          gtk_widget_add_events(line, GDK_BUTTON_RELEASE_MASK);
          gtk_widget_add_events(line, GDK_POINTER_MOTION_MASK);
          g_signal_connect(G_OBJECT(line), "expose-event", G_CALLBACK(onExposeDnaTrack), detailView);
          g_signal_connect(G_OBJECT(line), "button-press-event", G_CALLBACK(onButtonPressSeqColHeader), detailView);
          g_signal_connect(G_OBJECT(line), "button-release-event", G_CALLBACK(onButtonReleaseSeqColHeader), detailView);
          g_signal_connect(G_OBJECT(line), "motion-notify-event", G_CALLBACK(onMouseMoveSeqColHeader), detailView);
        }

      /* Set the header widget (and its refresh function) in the column data */
      BlxColumnInfo* columnInfo = getColumnInfo(columnList, BLXCOL_SEQUENCE);
      columnInfo->headerWidget = header;
      columnInfo->refreshFunc = refreshTextHeader;
    }
}


/* Create the horizontal scrollbar for custom scrolling. This takes ownership of
 * the GtkAdjustment. */
 GtkWidget* createDetailViewScrollBar(GtkAdjustment *adjustment, GtkWidget *detailView)
{
  GtkWidget *scrollBar = gtk_hscrollbar_new(adjustment);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollRangeChangedDetailView), detailView);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollPosChangedDetailView), detailView);
  
  return scrollBar;
}


/* Create the feedback box. (This feeds back info to the user about the currently-
 * selected base/sequence.) */
static GtkWidget* createFeedbackBox(GtkToolbar *toolbar, char *windowColor)
{
  GtkWidget *feedbackBox = gtk_entry_new() ;

  blxSetWidgetColor(feedbackBox, windowColor);

  /* User can copy text out but not edit contents */
  gtk_editable_set_editable(GTK_EDITABLE(feedbackBox), FALSE);

  /* want fixed width because feedback area needs as much space as possible - 
   * could do with a way to make sure this box is always wide enough though */
  const int numChars = 36; /* guesstimate of max number of chars we'll need */
  const int charWidth = 8; /* guesstimate of char width for default font */
  
  gtk_widget_set_size_request(feedbackBox, numChars * charWidth, -1) ;
  //GtkToolItem *item = addToolbarWidget(toolbar, feedbackBox, 0) ;
  //gtk_tool_item_set_expand(item, FALSE); 
  
  /* We want the box to be printed, so connect the expose function that will 
   * draw to a pixmap for printing */
  g_signal_connect(G_OBJECT(feedbackBox), "expose-event", G_CALLBACK(onExposePrintable), NULL);
  
  return feedbackBox;
}


/* Create the status bar for the detail-view toolbar. (This feeds back info to the user 
 * about the currently-moused-over sequence.) */
static GtkWidget* createStatusBar(GtkToolbar *toolbar, char *windowColor)
{
  GtkWidget *statusBar = gtk_statusbar_new() ;
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), FALSE);
  blxSetWidgetColor(statusBar, windowColor);

  /* Make it expandable so we use all available space. Set minimum size to be 0
   * because it's better to show it small than not at all. */
  gtk_widget_set_size_request(statusBar, 0, -1) ;
  //GtkToolItem *item = addToolbarWidget(toolbar, statusBar, 0) ;
  //gtk_tool_item_set_expand(item, TRUE); /* make as big as possible */

  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");

  return statusBar;
}


/* Create the detail view toolbar */
static GtkWidget* createDetailViewButtonBar(GtkWidget *detailView, 
                                            GtkWidget *toolbarIn,
                                            BlxBlastMode mode,
                                            const BlxColumnId sortColumn,
                                            GList *columnList,
                                            char *windowColor,
                                            GtkWidget **feedbackBox,
                                            GtkWidget **statusBar)
{
  GtkToolbar *toolbar = GTK_TOOLBAR(toolbarIn);

  blxSetWidgetColor(GTK_WIDGET(toolbar), windowColor);
  
  gtk_toolbar_set_style(toolbar, GTK_TOOLBAR_ICONS);
  gtk_toolbar_set_icon_size(toolbar, GTK_ICON_SIZE_SMALL_TOOLBAR);

  *feedbackBox = createFeedbackBox(toolbar, windowColor);
  *statusBar = createStatusBar(toolbar, windowColor);

  /* Create a parent hbox which will hold the main toolbar
   * and the detail-view tools. (We do this rather than putting
   * the detail-view tools directly into the toolbar so that
   * we have more control over spacing. In particular, the 
   * toolbar tools go into a nice drop-down box when there isn't
   * enough space, but the detail-view tools won't display if 
   * that happens, so we keep them separate.) */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, toolbarIn, TRUE, TRUE, 0);
  gtk_box_pack_start(hbox, *feedbackBox, FALSE, FALSE, 0);
  gtk_box_pack_start(hbox, *statusBar, TRUE, TRUE, 0);

  /* Put the toolbar in a handle box so that it can be torn off */
  GtkWidget *toolbarContainer = createToolbarHandle();
  gtk_container_add(GTK_CONTAINER(toolbarContainer), GTK_WIDGET(hbox));
  blxSetWidgetColor(toolbarContainer, windowColor);

  return toolbarContainer;
}


/* Create two detail-view trees and place them in the 2 panes of the container (if one is
 * given). The first tree gets associated with grid1 and appended to list1, and the second
 * with grid2 and list2. If hasHeaders is true, the first tree will have visible headers. */
static void createTwoPanedTrees(GtkWidget *detailView,
                                GtkPaned *panedWin, 
                                GtkCellRenderer *renderer,
                                GtkWidget *grid1, 
                                GtkWidget *grid2,
                                GList **list1,
                                GList **list2,
                                BlxSeqType seqType,
                                GList *columnList,
                                const char* const refSeqName,
                                const int frame1,
                                const int frame2,
                                const gboolean includeSnpTrack)
{
  GtkWidget *tree1 = createDetailViewTree(grid1, detailView, renderer, list1, columnList, seqType, refSeqName, frame1, includeSnpTrack);
  GtkWidget *tree2 = createDetailViewTree(grid2, detailView, renderer, list2, columnList, seqType, refSeqName, frame2, includeSnpTrack);
  
  if (panedWin)
    {
      gtk_paned_pack1(GTK_PANED(panedWin), tree1, TRUE, TRUE);
      gtk_paned_pack2(GTK_PANED(panedWin), tree2, TRUE, TRUE);
    }
}


/* Create three detail-view trees and place them in the paned widget. We only have
 * two panes in the widget, so two of the trees will be placed in a second, nested
 * paned widget. */
static void createThreePanedTrees(GtkWidget *detailView,
                                  GtkPaned *panedWin,
                                  GtkCellRenderer *renderer,
                                  GtkWidget *grid,
                                  GList **list,
                                  BlxSeqType seqType,
                                  const gboolean addToDetailView,
                                  GList *columnList,
                                  const char* const refSeqName,
                                  const gboolean includeSnpTrack)
{
  const int frame1 = 1, frame2 = 2, frame3 = 3;
  
  /* Create a tree for pane1 (but only add it to the detailView if instructed to).
   * The first tree has headers. */
  GtkWidget *tree1 = createDetailViewTree(grid, detailView, renderer, list, columnList, seqType, refSeqName, frame1, includeSnpTrack);

  GtkPaned *nestedPanedWidget = NULL;
  if (addToDetailView)
    {
      gtk_paned_pack1(GTK_PANED(panedWin), tree1, TRUE, TRUE);
      
      /* Create another paned widget and place it in pane 2 */
      nestedPanedWidget = GTK_PANED(gtk_vpaned_new());
      gtk_paned_pack2(panedWin, GTK_WIDGET(nestedPanedWidget), TRUE, TRUE);
    }
  
  /* Create two more trees (and place them in the nested paned widget, if it is not null). 
   * Neither of these trees should have headers. */
  createTwoPanedTrees(detailView, nestedPanedWidget, renderer, grid, grid, list, list, seqType, columnList, refSeqName, frame2, frame3, includeSnpTrack);
}


/* Create the trees for the detail view, creating sub-panes if necessary depending
 * on how many trees we need */
static void createDetailViewPanes(GtkWidget *detailView, 
                                  GtkPaned *panedWin,
                                  GtkCellRenderer *renderer,
                                  GtkWidget *fwdStrandGrid, 
                                  GtkWidget *revStrandGrid, 
                                  const int numFrames,
                                  GList **fwdStrandTrees,
                                  GList **revStrandTrees,
                                  BlxSeqType seqType,
                                  GList *columnList,
                                  const char* const refSeqName,
                                  const gboolean includeSnpTrack)
{
  if (numFrames == 1)
    {
      /* DNA matches: we need 2 trees, one for the forward strand and one for the reverse. */
      createTwoPanedTrees(detailView, 
                          panedWin, 
                          renderer,
                          fwdStrandGrid, 
                          revStrandGrid, 
                          fwdStrandTrees, 
                          revStrandTrees, 
                          seqType,
                          columnList,
                          refSeqName,
                          1,
                          1,
                          includeSnpTrack);
    }
  else if (numFrames == 3)
    {
      /* Protein matches: we need 3 trees for the 3 reading frames for EACH strand (although only
       * one set of trees will be displayed at a time). */
      createThreePanedTrees(detailView, panedWin, renderer, fwdStrandGrid, fwdStrandTrees, seqType, TRUE, columnList, refSeqName, includeSnpTrack);
      createThreePanedTrees(detailView, panedWin, renderer, revStrandGrid, revStrandTrees, seqType, FALSE, columnList, refSeqName, includeSnpTrack);
    }
}



/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* Add the MSPs to the detail-view trees. Calculates and returns the lowest ID 
 * out of all the blast matches. modelId specifies which tree data model should
 * be active at the start. If create is true, the tree data stores are created;
 * otherwise we assume they already exist. */
void detailViewAddMspData(GtkWidget *detailView, MSP *mspList, GList *seqList)
{
  BlxViewContext *bc = detailViewGetContext(detailView);

  /* First, make sure there is a data model for each tree */
  callFuncOnAllDetailViewTrees(detailView, treeCreateBaseDataModel, NULL);

  /* Loop through each MSP, and add it to the correct tree based on its strand and
   * reading frame. Also find the lowest ID value out of all the matches. */
  MSP *msp = mspList;
  
  for ( ; msp; msp = msp->next)
    {
      /* Only add matches/exons to trees. For exons, only add the parent exon;
       * it's child UTR/CDSs will be added automatically */
      if (mspIsBlastMatch(msp) || msp->type == BLXMSP_EXON)
        {
          /* Find the tree that this MSP should belong to based on its reading frame and strand */
          BlxStrand strand = mspGetRefStrand(msp);
          const int frame = mspGetRefFrame(msp, bc->seqType);
          GtkWidget *tree = detailViewGetTree(detailView, strand, frame);

          if (tree)
            {
              GtkListStore *store = GTK_LIST_STORE(treeGetBaseDataModel(GTK_TREE_VIEW(tree)));
              addMspToTree(msp, tree, store);
            }
          else
            {
              g_warning("Could not determine alignment list. Sequence may not be shown. (sequence = '%s', q range [%d-%d], s range [%d-%d])\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
            }
        }
    }
    
  /* Finally, create a custom-filtered version of the data store for each tree. We do 
   * this AFTER adding the data so that it doesn't try to re-filter every time we add a row.
   * To do: Note that if we're calling this function a second time to add additional rows, the filter
   * will already exist so this may be slow and we may need to look into improving this. */
  callFuncOnAllDetailViewTrees(detailView, treeCreateFilteredDataModel, NULL);
  
  /* Also create a second data store that will store one sequence per row (as opposed to one
   * MSP per row). This data store will be switched in when the user selects 'squash matches'. */
  callFuncOnAllDetailViewTrees(detailView, addSequencesToTree, seqList);
}


GtkWidget* createDetailView(GtkWidget *blxWindow,
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
                            char *windowColor)
{
  /* We'll group the trees in their own container so that we can pass them all around
   * together (so that operations like zooming and scrolling can act on the group). The
   * trees might be a direct child of this or a grandchild/grand-grandchild, so we will need
   * to look at all children recursively and check if they're the correct type. (We'll give 
   * all of our detail-view trees the same name so that we can identify them.) */
  GtkWidget *detailView = gtk_vbox_new(FALSE, 0);
  gtk_container_add(parent, detailView);
  
  GtkWidget *panedWin = gtk_vpaned_new();
  gtk_widget_set_name(panedWin, DETAIL_VIEW_WIDGET_NAME);

  /* Update the sequence column with a custom header type. (Must be done
   * before calling createDetailViewHeader) */
  createSeqColHeader(detailView, seqType, numFrames, columnList);

  /* Create the header bar. If viewing protein matches include one SNP track in the detail 
   * view header; otherwise create SNP tracks in each tree header. */
  const gboolean singleSnpTrack = (seqType == BLXSEQ_PEPTIDE);
  GtkWidget *header = createDetailViewHeader(detailView, seqType, numFrames, columnList, singleSnpTrack);

  /* Create the toolbar. We need to remember the feedback box and status bar so we can set them in the properties. */
  GtkWidget *feedbackBox = NULL;
  GtkWidget *statusBar = NULL;
  GtkWidget *buttonBar = createDetailViewButtonBar(detailView, toolbar, mode, sortColumn, columnList, windowColor, &feedbackBox, &statusBar);

  /* Create a custom cell renderer to render the sequences in the detail view */
  GtkCellRenderer *renderer = sequence_cell_renderer_new();

  /* Create the trees. */
  GList *fwdStrandTrees = NULL, *revStrandTrees = NULL;
  createDetailViewPanes(detailView, 
                        GTK_PANED(panedWin),
                        renderer,
                        fwdStrandGrid, 
                        revStrandGrid, 
                        numFrames, 
                        &fwdStrandTrees,
                        &revStrandTrees,
                        seqType,
                        columnList,
                        refSeqName,
                        !singleSnpTrack);
    
  gtk_box_pack_start(GTK_BOX(detailView), buttonBar, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(detailView), header, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(detailView), panedWin, TRUE, TRUE, 0);

  /* Connect signals */
  gtk_widget_add_events(detailView, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(detailView), "button-press-event", G_CALLBACK(onButtonPressDetailView), NULL);
  g_signal_connect(G_OBJECT(detailView), "button-release-event", G_CALLBACK(onButtonReleaseDetailView), NULL);
  g_signal_connect(G_OBJECT(detailView), "motion-notify-event", G_CALLBACK(onMouseMoveDetailView), NULL);

  detailViewCreateProperties(detailView, 
                             blxWindow, 
                             renderer,
                             fwdStrandTrees,
                             revStrandTrees,
                             header,
                             feedbackBox, 
                             statusBar,
                             columnList,
                             adjustment, 
                             startCoord,
                             sortColumn);

  return detailView;
}

