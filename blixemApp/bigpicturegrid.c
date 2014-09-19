/*  File: bigpicturegrid.c
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
 * Description: See bigpicturegrid.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/bigpicturegrid.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/detailview.h>
#include <blixemApp/detailviewtree.h>
#include <blixemApp/blxwindow.h>
#include <seqtoolsUtils/utilities.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#define BIG_PICTURE_MSP_LINE_NAME	"BigPictureMspLine"
#define DEFAULT_MSP_LINE_HEIGHT		3	  /* the height of the MSP lines in the grid */
#define DEFAULT_GRID_Y_PADDING		5	  /* this provides space between the grid and the edge of the widget */
#define MIN_MSP_LINE_WIDTH		1	  /* used to make sure that MSP lines never shrink to nothing */

typedef struct _DrawGridData
{
  GtkWidget *grid;
  GdkDrawable *drawable;
  GdkGC *gc;
  GdkColor *color;
  GdkColor *shadowColor;
  gboolean drawColinearityLines;
} DrawGridData;

  

/* Local function declarations */
static BlxViewContext*	    gridGetContext(GtkWidget *grid);
static IntRange*	    gridGetDisplayRange(GtkWidget *grid);
static GdkColor*	    gridGetMspLineColor(GtkWidget *grid, const gboolean selected);
static GtkWidget*	    gridGetBlxWindow(GtkWidget *grid);
static void                 drawBigPictureGrid(GtkWidget *grid, GdkDrawable *drawable);

/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

/* Calls the given function (passed as the data pointer) on the given widget 
 * if it is a grid in the big picture view, or, if it is a container, 
 * calls the function on all children/grandchildren/etc that are grids */
void callFuncOnAllBigPictureGrids(GtkWidget *widget, gpointer data)
{
  const gchar *name = gtk_widget_get_name(widget);
  if (strcmp(name, BIG_PICTURE_GRID_NAME) == 0)
    {
      GtkCallback func = (GtkCallback)data;
      func(widget, NULL);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllBigPictureGrids, data);
    }
}


static int gridGetNumVCells(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetNumVCells(bigPicture);
}


/* Convert an ID% value to the y coord on the given grid */
gint convertValueToGridPos(GtkWidget *grid, const gdouble value)
{
  /* The top line of the grid is drawn one cell height down from the top of the grid border */
  GridProperties *properties = gridGetProperties(grid);
  
  DoubleRange *valRange = bigPictureGetPercentIdRange(properties->bigPicture);
  gdouble percent = (valRange->max - value) / (valRange->max - valRange->min);
  
  /* Make sure we do the multiplication on doubles before rounding to int */
  gint result = properties->gridRect.y + roundNearest((gdouble)properties->gridRect.height * percent);
  return result;
}


/* Calculates the size and position of an MSP line in the given grid. Return
 * args can be null if not required. */
static void calculateMspLineDimensions(GtkWidget *grid, 
                                       const MSP* const msp, 
                                       gboolean clip,
                                       int *x, 
                                       int *y, 
                                       int *width, 
                                       int *height)
{
  BlxViewContext *bc = gridGetContext(grid);
  GridProperties *gridProperties = gridGetProperties(grid);

  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(gridGetDisplayRange(grid), bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dnaDispRange);

  /* The grid pos for coords gives the left edge of the coord, so draw to max + 1 to be inclusive */
  const int x1 = convertBaseIdxToRectPos(msp->qRange.min, &gridProperties->gridRect, &dnaDispRange, TRUE, bc->displayRev, TRUE);
  const int x2 = convertBaseIdxToRectPos(msp->qRange.max + 1, &gridProperties->gridRect, &dnaDispRange, TRUE, bc->displayRev, TRUE);
  
  const int xMin = min(x1, x2);
  const int xMax = max(x1, x2);
  
  if (x)
    {
      *x = xMin;
    }
    
  if (width)
    {
      *width = max((xMax - xMin), MIN_MSP_LINE_WIDTH);
    }
  
  /* Find where in the y axis we should draw the line, based on the %ID value */
  if (y)
    {
      *y = convertValueToGridPos(grid, msp->id);
    }
    
  if (height)
    {
      *height = gridProperties->mspLineHeight;
    }
}


/* Returns true if sequences in the given group should be shown. */
static gboolean isGroupVisible(const SequenceGroup* const group)
{
  gboolean result = TRUE;
  
  result = (!group || !group->hidden);
  
  return result;
}


/* Returns true if the given msp is displayed in the given grid, i.e. is the correct
 * strand, is not an intron or exon, and (if checkRange is true) is not out of range */
static gboolean mspShownInGrid(const MSP* const msp, GtkWidget *grid, gboolean checkRange)
{
  gboolean result = FALSE ; 

  SequenceGroup *group = NULL ;
  BlxViewContext *bc = NULL ;

  /* First check strand and type */
  result = (grid && msp && mspIsBlastMatch(msp) && mspGetRefStrand(msp) == gridGetStrand(grid)) ;

  /* If in a group, check that the group is visible */
  if (result)
    {
      bc = gridGetContext(grid) ;
      group = blxContextGetSequenceGroup(bc, msp->sSequence) ;
      result = isGroupVisible(group) ;
    }

  if (result)
    {
      /* If hiding ungrouped sequences and this is a sequence (i.e. match) without a group, hide it */
      if (!group && bc->flags[BLXFLAG_HIDE_UNGROUPED_SEQS] && msp->sSequence->type == BLXSEQUENCE_MATCH)
        result = FALSE ;
    }
  
  if (result)
    {
      if (checkRange)
        {
          /* See if the msp lies within the grid's display range */
          const IntRange* const displayRange = gridGetDisplayRange(grid);
          const IntRange* const mspRange = mspGetDisplayRange(msp);
          result = rangesOverlap(mspRange, displayRange);
        }
      else
        {
          /* We're just checking the correct type and strand, so return true */
          result = TRUE;
        }
    }
  
  return result;
}


/* Draw a line on the given grid to represent the given match */
static void drawMspLine(const MSP* const msp, DrawGridData *drawData)
{
  /* Ignore introns. */
  if (mspShownInGrid(msp, drawData->grid, TRUE))
    {
      gdk_gc_set_subwindow(drawData->gc, GDK_INCLUDE_INFERIORS);
      gdk_gc_set_foreground(drawData->gc, drawData->color);
      
      /* Calculate where it should go */
      int x, y, width, height;
      calculateMspLineDimensions(drawData->grid, msp, TRUE, &x, &y, &width, &height);
      
      /* Draw a block rectangle */
      gdk_draw_rectangle(drawData->drawable, drawData->gc, TRUE, x, y, width, height);
      
      /* Draw a drop-shadow, to make sure it is always visible even in the paler highlight colors */
      const int shadowHt = 1;
      gdk_gc_set_foreground(drawData->gc, drawData->shadowColor);
      gdk_gc_set_line_attributes(drawData->gc, shadowHt, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
      
      int yBottom = y + height - shadowHt;
      gdk_draw_line(drawData->drawable, drawData->gc, x, yBottom, x + width, yBottom);
      
      /* only draw the right-hand-side if the width is great enough to fit it and still be visible */
      if (width > shadowHt + 1)
      {
	int xRight = x + width - shadowHt;
	gdk_draw_line(drawData->drawable, drawData->gc, xRight, y, xRight, y + height);
      }
    }
}


/* Draw colinearity lines between the given MSP list item and the next one (if there is one) */
static void drawColinearityLines(GList *mspListItem, DrawGridData *drawData)
{
  g_return_if_fail(drawData);

  BlxViewContext *bc = gridGetContext(drawData->grid);
  
  if (drawData->drawColinearityLines && mspListItem && mspListItem->next && bc->flags[BLXFLAG_SHOW_COLINEARITY])
    {
      const IntRange* const displayRange = gridGetDisplayRange(drawData->grid);

      const MSP* msp1 = (const MSP*)mspListItem->data;
      const MSP* msp2 = (const MSP*)mspListItem->next->data;

      /* If the display is reversed we need to swap the order of the msps */
      if (bc->displayRev)
        {
          msp1 = (const MSP*)mspListItem->next->data;
          msp2 = (const MSP*)mspListItem->data;
        }

      /* Check that the line will lie within the visible range */
      if ((msp1->displayRange.max < displayRange->max && msp2->displayRange.min > displayRange->min) ||
          (msp1->displayRange.max > displayRange->max && msp2->displayRange.min < displayRange->min))
        {
          ColinearityType colinearityType = COLINEAR_INVALID;

          if ((bc->displayRev && msp1->qStrand == BLXSTRAND_FORWARD) || (!bc->displayRev && msp1->qStrand != BLXSTRAND_FORWARD))
            colinearityType = mspIsColinear(msp2, msp1);
          else
            colinearityType = mspIsColinear(msp1, msp2);

          if (colinearityType != COLINEAR_INVALID)
            {
              /* get the line color */
              GdkColor *color = NULL;

              if (colinearityType == COLINEAR_PERFECT)
                color = getGdkColor(BLXCOLOR_COLINEAR_PERFECT, bc->defaultColors, FALSE, bc->usePrintColors);
              else if (colinearityType == COLINEAR_IMPERFECT)
                color = getGdkColor(BLXCOLOR_COLINEAR_IMPERFECT, bc->defaultColors, FALSE, bc->usePrintColors);
              else
                color = getGdkColor(BLXCOLOR_COLINEAR_NOT, bc->defaultColors, FALSE, bc->usePrintColors);

              gdk_gc_set_foreground(drawData->gc, color);

              /* get coords of the msp lines */
              int x1, y1, width1, height1;
              int x2, y2, width2, height2;
              calculateMspLineDimensions(drawData->grid, msp1, FALSE, &x1, &y1, &width1, &height1);
              calculateMspLineDimensions(drawData->grid, msp2, FALSE, &x2, &y2, &width2, &height2);

              /* get coords of line to draw */
              IntRange lineXRange = {x1 + width1, x2};
              IntRange lineYRange = {y1 + height1 / 2, y2 + height2 / 2};

              gdk_draw_line(drawData->drawable, drawData->gc, 
                            lineXRange.min, lineYRange.min, lineXRange.max, lineYRange.max);
            }
        }
    }
}


/* Draw the MSPs for the given sequence in the given color. */
static void drawSequenceMspLines(gpointer listItemData, gpointer data)
{
  const BlxSequence *seq = (const BlxSequence*)listItemData;
  
  /* Only applicable to alignments and transcripts */
  if (!blxSequenceShownInGrid(seq))
    return;
  
  DrawGridData *drawData = (DrawGridData*)data;  
  GList *mspListItem = seq->mspList;

  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);

      /* Draw the msp (if it's on this grid and is in the visible range) */
      if (mspShownInGrid(msp, drawData->grid, TRUE))
        drawMspLine(msp, drawData);

      /* Draw colinearity lines between this and the next msp. Note that we need to call this
       * even if this msp is out of range, in case the next msp is in range. */
      if (mspShownInGrid(msp, drawData->grid, FALSE))
        drawColinearityLines(mspListItem, drawData);
    }
}


/* Draw MSP lines that are in groups. Use the highlight color of the group if one is set */
static void drawGroupedMspLines(gpointer listItemData, gpointer data)
{
  SequenceGroup *seqGroup = (SequenceGroup*)listItemData;
  DrawGridData *drawData = (DrawGridData*)data;
  GdkColor *origColor = drawData->color;
  GdkColor *origShadowColor = drawData->shadowColor;

  /* Use the group's highlight color, if this sequence should be highlighted */
  if (seqGroup->highlighted)
    {
      /* Use the group's highlight color */
      drawData->color = &seqGroup->highlightColor;
      GdkColor shadowColor;
      getDropShadowColor(drawData->color, &shadowColor);
      drawData->shadowColor = &shadowColor;
      
      g_list_foreach(seqGroup->seqList, drawSequenceMspLines, drawData);

      /* Set the draw data back to its original values */
      drawData->color = origColor;
      drawData->shadowColor = origShadowColor;
    }
  else
    {
      g_list_foreach(seqGroup->seqList, drawSequenceMspLines, drawData);
    }
  
}


/* Draw a line for each MSP in the given grid */
static void drawMspLines(GtkWidget *grid, GdkDrawable *drawable)
{
  BlxViewContext *bc = gridGetContext(grid);
  GdkGC *gc = gdk_gc_new(drawable);
  
  DrawGridData drawData = {
    grid, 
    drawable, 
    gc, 
    gridGetMspLineColor(grid, FALSE), 
    gridGetMspLineColor(grid, FALSE),
    FALSE
  };
  
  /* Draw all MSPs for this grid */ 
  g_list_foreach(bc->matchSeqs, drawSequenceMspLines, &drawData);  

  /* Now draw MSPs that are in groups (to do: it would be good to do this in reverse
   * Sort Order, so that those ordered first get drawn last and therefore appear on top) */
  g_list_foreach(bc->sequenceGroups, drawGroupedMspLines, &drawData);
  
  /* Finally, draw selected sequences. These will appear on top of everything else and will have
   * colinearity lines displayed (we only ever display colinearity lines for selected sequences
   * in the grid because it would be too cluttered otherwise). */
  drawData.color = gridGetMspLineColor(grid, TRUE);
  GdkColor shadowColor;
  getDropShadowColor(drawData.color, &shadowColor);
  drawData.shadowColor = &shadowColor;
  drawData.drawColinearityLines = TRUE;
  g_list_foreach(bc->selectedSeqs, drawSequenceMspLines, &drawData);
  
  g_object_unref(gc);
}


/* Main function that does the drawing for the grid. Drawing is done onto a pixmap
 * which is then stored in the grid properties */
static void drawBigPictureGrid(GtkWidget *grid, GdkDrawable *drawable)
{
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);

  /* Calculate some factors for scaling */
  const gdouble percentPerCell = bigPictureGetIdPerCell(properties->bigPicture);
  const gint numVCells = gridGetNumVCells(grid);

  /* Highlight any regions that are assembly gaps */
  /* Get the display range in dna coords */
  IntRange bpRange;
  convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &bpRange);
  
  GdkColor *color = getGdkColor(BLXCOLOR_ASSEMBLY_GAP, bc->defaultColors, FALSE, bc->usePrintColors);

  drawAssemblyGaps(grid, drawable, color, bc->displayRev,
                   &properties->gridRect, &bpRange, 
                   bc->featureLists[BLXMSP_GAP]);
  
  /* Draw the grid lines */
  drawVerticalGridLines(&properties->gridRect, &properties->highlightRect, properties->gridYPadding, bc, bpProperties, drawable);
  drawHorizontalGridLines(grid, properties->bigPicture, &properties->gridRect, bc, bpProperties, drawable, numVCells, percentPerCell, bpProperties->percentIdRange.max, FALSE, "%");
  
  /* Draw lines corresponding to the MSPs */
  drawMspLines(grid, drawable);
}


/* Prepare the grid for printing (draws the transient hightlight box
 * onto the cached drawable). */
void gridPrepareForPrinting(GtkWidget *grid)
{
  GdkDrawable *drawable = widgetGetDrawable(grid);
  
  if (drawable)
    {
      GridProperties *properties = gridGetProperties(grid);
      BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
      BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
      
      GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
      drawHighlightBox(drawable, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);
    }
}


void calculateGridHighlightBoxBorders(GtkWidget *grid)
{
  DEBUG_ENTER("calculateGridHighlightBoxBorders(grid)");

  GridProperties *properties = gridGetProperties(grid);
  calculateHighlightBoxBorders(&properties->gridRect, &properties->highlightRect, properties->bigPicture, properties->mspLineHeight);
    
  DEBUG_EXIT("calculateGridHighlightBoxBorders returning");
}


/* Calculate the borders for the big picture view */
void calculateGridBorders(GtkWidget *grid)
{
  DEBUG_ENTER("calculateGridBorders");

  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  int numVCells = gridGetNumVCells(grid);
  
  /* Get some info about the size of the layout and the font used */
  guint layoutWidth, layoutHeight;
  gtk_layout_get_size(GTK_LAYOUT(grid), &layoutWidth, &layoutHeight);
  
  properties->displayRect.x = 0;
  properties->displayRect.y = 0;
  properties->displayRect.width = grid->allocation.width;
  
  /* Get the boundaries of the grid */
  properties->gridRect.x = roundNearest(bigPictureProperties->charWidth * (gdouble)bigPictureProperties->leftBorderChars);
  properties->gridRect.y = bigPictureProperties->highlightBoxYPad + properties->gridYPadding;
  properties->gridRect.width = properties->displayRect.width - properties->gridRect.x;
  properties->gridRect.height = bigPictureGetCellHeight(properties->bigPicture) * numVCells;
  
  /* Get the boundaries of the highlight box */
  calculateGridHighlightBoxBorders(grid);
  
  /* Get the total display height required. Set the layout size to fit. */
  const int newHeight = properties->highlightRect.height;
  
  if (newHeight != properties->displayRect.height)
    {
      DEBUG_OUT("Setting new grid height = %d\n", newHeight);
      properties->displayRect.height = newHeight;
      gtk_layout_set_size(GTK_LAYOUT(grid), properties->displayRect.width, properties->displayRect.height);
  
      /* Set the size request to our desired height. We want a fixed heigh but don't set the
       * width, because we want the user to be able to resize that. */
      gtk_widget_set_size_request(grid, 0, properties->displayRect.height);
    }
  
  DEBUG_EXIT("calculateGridBorders returning");
}


/***********************************************************
 *			    Events			   *
 ***********************************************************/

/* Expose handler for the grid. If the grid has a bitmap cached from the previous
 * call, this function just pushes that back to screen and then (optionally) draws
 * the preview box directly to the window over the top of it. If the cached bitmap
 * is null, this function first calls drawBigPictureGrid to create and cache
 * the bitmap first. */
static gboolean onExposeGrid(GtkWidget *grid, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(grid)->bin_window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(grid);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(grid);
          drawBigPictureGrid(grid, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);

          /* Draw the highlight box on top of it */
          GridProperties *properties = gridGetProperties(grid);
          BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
          BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
          
          GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
          drawHighlightBox(window, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);
          
          /* Draw the preview box too, if set */
          drawPreviewBox(properties->bigPicture, window, &properties->gridRect, &properties->highlightRect);
        }
      else
	{
	  g_warning("Failed to draw grid [%p] - could not create bitmap.\n", grid);
	}
    }
  
  return TRUE;
}


/* Re-calculate size info about the given grid following a resize */
static void onSizeAllocateBigPictureGrid(GtkWidget *grid, GtkAllocation *allocation, gpointer data)
{
  DEBUG_ENTER("onSizeAllocateBigPictureGrid");

  /* Recalculate the borders for the grids and the header */
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  calculateGridHeaderBorders(bigPictureProperties->header);
  calculateGridBorders(grid);
  
  DEBUG_EXIT("onSizeAllocateBigPictureGrid returning");
}


/* Mark the given MSP's sequence as selected if this MSP line contains the given coords.
 * Returns true if it was selected. */
static gboolean selectMspIfContainsCoords(GtkWidget *grid, 
					  const MSP *msp,
					  const int x,
					  const int y,
					  gboolean deselectOthers)
{
  gboolean wasSelected = FALSE;
  
  if (mspShownInGrid(msp, grid, TRUE))
    {
      int mspX, mspY, mspWidth, mspHeight;
      calculateMspLineDimensions(grid, msp, TRUE, &mspX, &mspY, &mspWidth, &mspHeight);
      
      if (x >= mspX && x <= mspX + mspWidth && y >= mspY && y <= mspY + mspHeight)
	{
	  /* It's a hit. Select this sequence. */
	  GtkWidget *blxWindow = gridGetBlxWindow(grid);
	  
	  if (deselectOthers)
	    {
	      blxWindowDeselectAllSeqs(blxWindow);
	    }
	  
	  blxWindowSelectSeq(blxWindow, msp->sSequence);
	  
	  /* Update the selected strand */
	  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
	  detailViewSetSelectedStrand(detailView, gridGetStrand(grid));

	  /* Scroll the detail view trees to bring the new selection into view */
	  callFuncOnAllDetailViewTrees(detailView, treeScrollSelectionIntoView, NULL);

	  wasSelected = TRUE;
	}
    }
  
  return wasSelected;
}


/* Loop through all the msp lines for this grid and mark them as selected
 * if they contain the coords of the mouse press */
static gboolean selectClickedMspLines(GtkWidget *grid, GdkEventButton *event, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  /* Loop through all the MSPs until we find one under the click coords */
  GtkWidget *blxWindow = gridGetBlxWindow(grid);
  const MSP *msp = blxWindowGetMspList(blxWindow);
  const gboolean deselectOthers = !ctrlModifier && !shiftModifier; /* whether to deselect all others first */
  gboolean found = FALSE;

  for ( ; msp && !found; msp = msp->next)
    {
      found = selectMspIfContainsCoords(grid, msp, event->x, event->y, deselectOthers);
    }

  return found;
}


static gboolean onButtonPressGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  /* left button */
  if (event->button == 1)
    {
      /* If we clicked on top of an msp line, select that msp */
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
      const gboolean shiftModifier = ((event->state & modifiers) == GDK_SHIFT_MASK);
      
      handled = selectClickedMspLines(grid, event, ctrlModifier, shiftModifier);
    }
  
  /* Middle button: always show the preview box; left button: show
   * preview box if we double-clicked, or if we clicked in the highlight
   * box (i.e. left button selects and drags the highlight box; middle 
   * button or double-click makes the highlight box jump) */
  GridProperties *properties = gridGetProperties(grid);
  BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
  
  if (event->button == 2 || 
      (event->button == 1 && !handled && 
       (event->type == GDK_2BUTTON_PRESS || 
        clickedInRect(event, &properties->highlightRect, bpProperties->highlightBoxMinWidth))))
    {
      /* If dragging the highlight box (left button), then centre the preview 
       * box on the existing highlight box centre; otherwise, centre it on the click pos */
      int x = event->x;
      
      if (event->button == 1 && event->type == GDK_BUTTON_PRESS)
        x = properties->highlightRect.x + properties->highlightRect.width / 2;

      showPreviewBox(gridGetBigPicture(grid), event->x, TRUE, x - event->x);
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseGrid(GtkWidget *grid, GdkEventButton *event, gpointer data)
{
  if (event->button == 1 || event->button == 2) /* left or middle button */
    {
      GridProperties *properties = gridGetProperties(grid);
      acceptAndClearPreviewBox(gridGetBigPicture(grid), event->x, &properties->gridRect, &properties->highlightRect);
    }
  
  return TRUE;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the grid.
 * This scrolls the position of the highlight box, i.e. it scrolls the display
 * range in the detail view. */
static gboolean onScrollGrid(GtkWidget *grid, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
    {
      case GDK_SCROLL_LEFT:
	{
          scrollBigPictureLeftStep(gridGetBigPicture(grid));
	  handled = TRUE;
	  break;
	}
	
      case GDK_SCROLL_RIGHT:
	{
          scrollBigPictureRightStep(gridGetBigPicture(grid));
	  handled = TRUE;
	  break;
	}

      default:
	{
	  handled = FALSE;
	  break;
	}
    };
  
  return handled;
}


static gboolean onMouseMoveGrid(GtkWidget *grid, GdkEventMotion *event, gpointer data)
{
  if ((event->state & GDK_BUTTON1_MASK) || /* left button */
      (event->state & GDK_BUTTON2_MASK))   /* middle button */
    {
      /* Draw a preview box at the mouse pointer location */
      showPreviewBox(gridGetBigPicture(grid), event->x, FALSE, 0);
    }
  
  return TRUE;
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

GridProperties* gridGetProperties(GtkWidget *widget)
{
  return widget ? (GridProperties*)(g_object_get_data(G_OBJECT(widget), "GridProperties")) : NULL;
}

static BlxViewContext* gridGetContext(GtkWidget *grid)
{
  GtkWidget *blxWindow = gridGetBlxWindow(grid);
  return blxWindowGetContext(blxWindow);
}

static void onDestroyGrid(GtkWidget *widget)
{
  GridProperties *properties = gridGetProperties(widget);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "GridProperties", NULL);
    }
}


static void gridCreateProperties(GtkWidget *widget, 
				 GtkWidget *bigPicture,
				 gulong exposeHandlerId,
				 BlxStrand strand)
{
  if (widget)
    { 
      /* Work out how many characters we need to display in the left border for
       * the vertical grid axis label.  This is the number of digits in the 
       * largest label value plus one (for a '%' character)
       * Add a fudge factor to give more space to allow for the fact that 
       * the calculated char width is approximate and may not give enough space */
      GridProperties *properties = (GridProperties*)g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->strand = strand;
      
      properties->gridYPadding = DEFAULT_GRID_Y_PADDING;
      
      properties->mspLineHeight = DEFAULT_MSP_LINE_HEIGHT;
      properties->exposeHandlerId = exposeHandlerId;
      properties->ignoreSizeAlloc = FALSE;
      
      g_object_set_data(G_OBJECT(widget), "GridProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyGrid), NULL); 
    }
}


static GtkWidget* gridGetBlxWindow(GtkWidget *grid)
{
  GridProperties *properties = grid ? gridGetProperties(grid) : NULL;
  return properties ? bigPictureGetBlxWindow(properties->bigPicture) : NULL;
}


BlxStrand gridGetStrand(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties->strand;
}

GtkWidget* gridGetBigPicture(GtkWidget *grid)
{
  GridProperties *properties = gridGetProperties(grid);
  return properties ? properties->bigPicture : NULL;
}

static IntRange* gridGetDisplayRange(GtkWidget *grid)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  return bigPictureGetDisplayRange(bigPicture);
}

static GdkColor *gridGetMspLineColor(GtkWidget *grid, const gboolean selected)
{
  GtkWidget *bigPicture = gridGetBigPicture(grid);
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  return getGdkColor(BLXCOLOR_MSP_LINE, bc->defaultColors, selected, bc->usePrintColors);
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

GtkWidget* createBigPictureGrid(GtkWidget *bigPicture, BlxStrand strand)
{
  /* Create a layout area for the grid */
  GtkWidget *grid = gtk_layout_new(NULL, NULL);
  
  /* Grid style properties */
  gtk_widget_set_redraw_on_allocate(grid, FALSE);
  gtk_widget_set_name(grid, BIG_PICTURE_GRID_NAME);
  
  /* Connect signals */
  gtk_widget_add_events(grid, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(grid, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(grid, GDK_POINTER_MOTION_MASK);
  
  gulong exposeHandlerId = g_signal_connect(G_OBJECT(grid), "expose-event", G_CALLBACK(onExposeGrid), NULL);  
  g_signal_connect(G_OBJECT(grid), "size-allocate",	    G_CALLBACK(onSizeAllocateBigPictureGrid), NULL);
  g_signal_connect(G_OBJECT(grid), "button-press-event",    G_CALLBACK(onButtonPressGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "button-release-event",  G_CALLBACK(onButtonReleaseGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "motion-notify-event",   G_CALLBACK(onMouseMoveGrid),	      NULL);
  g_signal_connect(G_OBJECT(grid), "scroll-event",	    G_CALLBACK(onScrollGrid),		      NULL);
  
  /* Set required data in the grid. We can't set the tree yet because it hasn't been created yet. */
  gridCreateProperties(grid, bigPicture, exposeHandlerId, strand);
  
  return grid;
}


