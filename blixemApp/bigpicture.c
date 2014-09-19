/*  File: bigpicture.c
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
 * Description: See bigpicture.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/bigpicture.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/exonview.h>
#include <blixemApp/coverageview.h>
#include <seqtoolsUtils/utilities.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEFAULT_PREVIEW_BOX_LINE_WIDTH  1
#define DEFAULT_GRID_NUM_HEADER_LINES   1	  /* the default number of lines of text in the grid header */
#define DEFAULT_GRID_HEADER_Y_PAD	0	  /* the default y padding around the grid header */
#define DEFAULT_LABEL_X_PADDING		5	  /* padding around the grid labels */
#define DEFAULT_LABEL_Y_PADDING		-2	  /* padding around the grid labels */
#define DEFAULT_GRID_CELL_WIDTH		100	  /* the default cell width of the grids */
#define DEFAULT_GRID_NUM_HOZ_CELLS	5	  /* the default number of cells to show horizontally in the grids */
#define DEFAULT_PERCENT_ID_PER_CELL	20	  /* the default %ID per vertical cell to show in the grids */
#define DEFAULT_GRID_PERCENT_ID_MAX	100	  /* default maximum %ID to show on the scale */
#define MIN_NUM_V_CELLS			1	  /* minimum number of vertical cells to show in the grid */
#define DEFAULT_HIGHLIGHT_BOX_Y_PAD	2	  /* this provides space between highlight box and the top/bottom of the grid */
#define MIN_HIGHLIGHT_BOX_WIDTH         5         /* minimum width of the highlight box */
#define GRID_SCALE_MIN_ID_PER_CELL      0.1       /* minimum %ID per grid cell */
#define GRID_SCALE_MIN                  0         /* minimum possible value for grid scale */
#define GRID_SCALE_MAX                  100       /* maximum possible value for grid scale */
#define BIG_PICTURE_GRID_HEADER_NAME	"BigPictureGridHeader"
#define BIG_PICTURE_WIDGET_NAME		"BigPictureWidget" /* name of the direct parent of the grids etc. */
#define MAX_BIG_PICTURE_HEIGHT_RATIO	0.4	  /* max height of the big picture wrt the whole window size */
#define BLX_SCROLL_INCREMENT_RATIO      20        /* determines speed of scrolling of big picture
                                                   * (we use 1/nth of the display range as the scroll step) */

/* Local function declarations */
static GridHeaderProperties*	    gridHeaderGetProperties(GtkWidget *gridHeader);
static int			    bigPictureGetInitialZoom(GtkWidget *bigPicture);

static void                         drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc);

/***********************************************************
 *                     Utility functions	           *
 ***********************************************************/

/* This causes a complete redraw of all the big picture components (including
 * the coverage view, which depends on some big picture properties). */
void bigPictureRedrawAll(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  widgetClearCachedDrawable(properties->header, NULL);
  widgetClearCachedDrawable(properties->fwdStrandGrid, NULL);
  widgetClearCachedDrawable(properties->revStrandGrid, NULL);
  widgetClearCachedDrawable(properties->fwdExonView, NULL);
  widgetClearCachedDrawable(properties->revExonView, NULL);
  coverageViewRedraw(properties->coverageView);
  
  gtk_widget_queue_draw(bigPicture);
}


/* This just causes an expose on all the big picture child widgets, and also
 * on the coverage view */
static void bigPictureRefreshAll(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);

  gtk_widget_queue_draw(bigPicture);
  gtk_widget_queue_draw(properties->coverageView);
}


/* This function calculates the cell size and number of cells for the big picture grids.
 * It should be called whenever the big picture is resized or its display range changes. */
void calculateBigPictureCellSize(GtkWidget *bigPicture, BigPictureProperties *properties)
{
  DEBUG_ENTER("calculateBigPictureCellSize");

  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *header = properties->header;
  GridHeaderProperties *headerProperties = gridHeaderGetProperties(header);
  
  /* Calculate the number of bases per cell and round it to a "nice" value from our stored list */
  const int displayStart = convertDisplayIdxToDnaIdx(properties->displayRange.min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int displayEnd = convertDisplayIdxToDnaIdx(properties->displayRange.max, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int displayWidth = abs(displayEnd - displayStart); /* use abs because can be negative if display reversed */
  
  const int defaultBasesPerCell = ceil((double)(displayWidth) / DEFAULT_GRID_NUM_HOZ_CELLS);
  properties->basesPerCell = roundToValueFromList(defaultBasesPerCell, properties->roundValues, &properties->roundTo);

  /* Adjust the cell width (and hence number of cells) proportionally wrt the rounded number of bases per cell */
  const double defaultCellWidth = (double)(headerProperties->headerRect.width) / DEFAULT_GRID_NUM_HOZ_CELLS;
  const double actualCellWidth = (defaultCellWidth * properties->basesPerCell) / defaultBasesPerCell;
  properties->numHCells = ceil((double)headerProperties->headerRect.width / actualCellWidth);
  
  DEBUG_EXIT("calculateBigPictureCellSize returning");
}


/* Utility to convert the given decimal number (percent) to text 
 * returned in 'text'). showDecimal indicates whether it should be
 * shown as a decimal or not. If abbrev is true, 1000 is abbreviated 
 * as 1k etc. 'unit' is a string displayed after the text and can be
 * an empty string. */
static void drawNumericLabel(char *text, 
                             const gdouble percent,
                             const gboolean showDecimal, 
                             const gboolean abbrev, 
                             const char *unit)
{
  if (!text)
    {
      return;
    }
  
  if (showDecimal)
    {
      sprintf(text, "%1.1f%s", percent, unit);
    }
  else 
    {
      sprintf(text, "%d%s", (int)percent, unit);
    
      if (abbrev)
        {
          /* Abbreviate the number so 1000 becomes 1k, 1000000 becomes 1M etc. */
          const int len = strlen(text);
          int i = 3;

          for ( ; i < len; i += 3)
            {
              char suffix = 0;

              switch (i)
                {
                  case 3:  suffix = 'k'; break;
                  case 6:  suffix = 'M'; break;
                  case 9:  suffix = 'G'; break;
                  case 12: suffix = 'T'; break;
                  case 15: suffix = 'P'; break;
                  case 18: suffix = 'E'; break;
                  case 21: suffix = 'Z'; break;
                  case 24: suffix = 'Y'; break;
                  default: break;
                };
  
              if (suffix && text[len-i]=='0' && text[len-i+1]=='0' && text[len-i+2]=='0')
                {
                  text[len-i] = suffix;
                  text[len-i+1] = '\0';
                }
            }
        }
    }
}


static void drawVerticalGridLineHeaders(GtkWidget *header, 
					GtkWidget *bigPicture, 
                                        GdkDrawable *drawable,
					GdkGC *gc, 
					const GdkColor* const textColor, 
					const GdkColor* const lineColor,
                                        const gboolean abbrev)
{
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  GridHeaderProperties *headerProperties = gridHeaderGetProperties(header);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
  
  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dnaDispRange);

  const int direction = bc->displayRev ? -1 : 1; /* to subtract instead of add when display reversed */
  
  /* Get the first base index and round it to a nice round number. We'll offset all of the gridlines 
   * by the distance between this and the real start coord. */
  const int firstBaseIdx = roundToValue(bc->displayRev ? dnaDispRange.max : dnaDispRange.min, bpProperties->roundTo);
  
  /* Calculate the top and bottom heights for the lines. */
  const gint bottomBorder = headerProperties->headerRect.y + headerProperties->headerRect.height;
  const gint topBorder = bottomBorder - headerProperties->markerHeight;

  const int minX = headerProperties->headerRect.x;
  const int maxX = headerProperties->headerRect.x + headerProperties->headerRect.width;

  /* Loop through each grid cell */
  gint hCell = 0;
  for ( ; hCell <= bpProperties->numHCells; ++hCell)
    {
      /* Draw the label, showing which base index is at this x coord */
      int numBasesFromLeft = bpProperties->basesPerCell * hCell;
      int baseIdx = firstBaseIdx + (numBasesFromLeft * direction);

      const int x = convertBaseIdxToRectPos(baseIdx, &headerProperties->headerRect, &dnaDispRange, TRUE, bc->displayRev, TRUE);
      
      if (x > minX && x < maxX)
	{
          /* If we're displaying negative coords, negate the base index */
          if (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS])
            baseIdx *= -1;
            
	  gdk_gc_set_foreground(gc, textColor);
	  gchar text[numDigitsInInt(baseIdx) + 1];
          drawNumericLabel(text, baseIdx, FALSE, abbrev, "");

	  PangoLayout *layout = gtk_widget_create_pango_layout(header, text);
	  gdk_draw_layout(drawable, gc, x, 0, layout);
	  g_object_unref(layout);
	  
	  /* Draw a marker line (small cosmetic touch to better join the label to the grid below) */
	  gdk_gc_set_foreground(gc, lineColor);
	  gdk_draw_line (drawable, gc, x, topBorder, x, bottomBorder);
	}
    }
}


/* Draw the vertical gridlines for a big picture component. This is a generic routine used
 * by the grids and the coverage view, so that they all have their gridlines spaced the same */
void drawVerticalGridLines(GdkRectangle *drawingRect,
			   GdkRectangle *highlightRect,
			   const int yPadding,
			   BlxViewContext *bc,
			   BigPictureProperties *bpProperties,
			   GdkDrawable *drawable)
{
  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dnaDispRange);
  
  const int direction = bc->displayRev ? -1 : 1; /* to subtract instead of add when display reversed */
  
  /* Get the first base index (in terms of the nucleotide coords) and round it to a nice round
   * number. We'll offset all of the gridlines by the distance between this and the real start coord. */
  const int realFirstBaseIdx = convertDisplayIdxToDnaIdx(bpProperties->displayRange.min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int firstBaseIdx = roundToValue(realFirstBaseIdx, bpProperties->roundTo);
  
  /* Calculate the top and bottom heights for the lines. */
  const gint topBorder = highlightRect->y - yPadding;
  const gint bottomBorder = drawingRect->y + drawingRect->height;
  
  const int minX = drawingRect->x;
  const int maxX = drawingRect->x + drawingRect->width;
  
  gint hCell = 0;
  for ( ; hCell <= bpProperties->numHCells; ++hCell)
    {
      /* Get the base index for this grid line and calc its x coord */
      int numBasesFromLeft = bpProperties->basesPerCell * hCell;
      int baseIdx = firstBaseIdx + (numBasesFromLeft * direction);
      
      const int x = convertBaseIdxToRectPos(baseIdx, drawingRect, &dnaDispRange, TRUE, bc->displayRev, TRUE);
      
      if (x > minX && x < maxX)
	{
	  GdkGC *gc = gdk_gc_new(drawable);

	  GdkColor *lineColor = getGdkColor(BLXCOLOR_GRID_LINE, bc->defaultColors, FALSE, bc->usePrintColors);
	  gdk_gc_set_foreground(gc, lineColor);
	  gdk_draw_line (drawable, gc, x, topBorder, x, bottomBorder);
          
          g_object_unref(gc);
	}
    }
}


/* Get the standard height for grid cells */
gint bigPictureGetCellHeight(GtkWidget *bigPicture)
{
  /* Base the cell height on the font height */
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return roundNearest(properties->charHeight + (gdouble)(2 * DEFAULT_LABEL_Y_PADDING));
}


/* Draw the horizontal grid lines for the big picture view */
void drawHorizontalGridLines(GtkWidget *widget,
			     GtkWidget *bigPicture,
			     GdkRectangle *drawingRect,
			     BlxViewContext *bc,
			     BigPictureProperties *bpProperties,
			     GdkDrawable *drawable,
			     const gint numCells, 
			     const gdouble rangePerCell, 
			     const gdouble maxVal,
                             const gboolean abbrev,
			     const char *unit)
{
  const gint rightBorder = drawingRect->x + drawingRect->width;
  
  GdkColor *textColor = getGdkColor( BLXCOLOR_GRID_TEXT, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkGC *textGc = gdk_gc_new(drawable);
  gdk_gc_set_foreground(textGc, textColor);
  
  GdkColor *lineColor = getGdkColor(BLXCOLOR_GRID_LINE, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkGC *lineGc = gdk_gc_new(drawable);
  gdk_gc_set_foreground(lineGc, lineColor);
  
  /* Show decimal places if the range per cell is a fraction of a percent */
  const gboolean showDecimal = (rangePerCell < 1.0);
  const int cellHeight = bigPictureGetCellHeight(bigPicture);
  
  gint vCell = 0;
  for ( ; vCell <= numCells; ++vCell)
    {
      gint y = drawingRect->y + (gint)((gdouble)vCell * cellHeight);
      gint x = drawingRect->x - DEFAULT_LABEL_X_PADDING;
      
      /* Label this gridline with the %ID */
      gdouble percent = maxVal - (rangePerCell * vCell);
      char text[bpProperties->leftBorderChars + 3]; /* +3 to include decimal point, 1dp, and terminating nul */

      drawNumericLabel(text, percent, showDecimal, abbrev, unit);
      
      PangoLayout *layout = gtk_widget_create_pango_layout(widget, text);
      
      int width = UNSET_INT;
      pango_layout_get_pixel_size(layout, &width, NULL);

      gdk_draw_layout(drawable, textGc, x - width, y - cellHeight/2, layout);
      g_object_unref(layout);
      
      /* Draw the gridline */
      gdk_draw_line (drawable, lineGc, drawingRect->x, y, rightBorder, y);
    }
  
  g_object_unref(lineGc);
  g_object_unref(textGc);
}



/* Refresh the header - clears and redraws its bitmap */
static void redrawBigPictureGridHeader(GtkWidget *header)
{
  /* Check that the header is shown on screen */
  if (GTK_LAYOUT(header)->bin_window)
    {
      /* Create a new bitmap to draw on to and set it in the widget. (This automatically
       * deletes the old one, if there is one.) */
      GdkDrawable *bitmap = gdk_pixmap_new(GTK_LAYOUT(header)->bin_window, header->allocation.width, header->allocation.height, -1);
      gdk_drawable_set_colormap(bitmap, gdk_colormap_get_system());
      widgetSetDrawable(header, bitmap);

      /* Clear the bitmap to the background color */
      GdkGC *gc = gdk_gc_new(bitmap);
      
      GtkStyle *style = gtk_widget_get_style(header);
      GdkColor *bgColor = &style->bg[GTK_STATE_NORMAL];
      gdk_gc_set_foreground(gc, bgColor);
      gdk_draw_rectangle(bitmap, gc, TRUE, 0, 0, header->allocation.width, header->allocation.height);

      /* Draw the header */
      drawBigPictureGridHeader(header, bitmap, gc);
      
      g_object_unref(gc);
    }
}



/* Draw the big picture header onto the given drawable */
static void drawBigPictureGridHeader(GtkWidget *header, GdkDrawable *drawable, GdkGC *gc)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  
  /* Set the drawing properties */
  gdk_gc_set_subwindow(gc, GDK_INCLUDE_INFERIORS);
  
  /* First, highlight any assembly gaps */
  /* Get the display range in dna coords */
  const IntRange* const displayRange = bigPictureGetDisplayRange(properties->bigPicture);
  IntRange bpRange;
  convertDisplayRangeToDnaRange(displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &bpRange);
  
  GdkColor *gapColor = getGdkColor(BLXCOLOR_ASSEMBLY_GAP, bc->defaultColors, FALSE, bc->usePrintColors);
  drawAssemblyGaps(header, drawable, gapColor, bc->displayRev, &properties->headerRect, &bpRange, bc->featureLists[BLXMSP_GAP]);
  
  /* Draw the grid headers */
  drawVerticalGridLineHeaders(header, 
			      properties->bigPicture, 
                              drawable,
			      gc,
			      getGdkColor(BLXCOLOR_GRID_TEXT, bc->defaultColors, FALSE, bc->usePrintColors), 
			      getGdkColor(BLXCOLOR_GRID_LINE, bc->defaultColors, FALSE, bc->usePrintColors),
                              FALSE);
}


/* Recalculate the borders for the grid header. Should be called after a resize */
void calculateGridHeaderBorders(GtkWidget *header)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(header);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Calculate the size of the grid header (zero height if it does not have one) */
  properties->headerRect.x = roundNearest(bigPictureProperties->charWidth * (gdouble)bigPictureProperties->leftBorderChars);
  properties->headerRect.y = 0;
  properties->headerRect.width = header->allocation.width - properties->headerRect.x;
  properties->headerRect.height = properties->refButton->allocation.height + (properties->headerYPad * 2);
  
  properties->markerHeight = properties->headerRect.height - roundNearest(bigPictureProperties->charHeight * (gdouble)properties->numHeaderLines);
  if (properties->markerHeight < 0)
    properties->markerHeight = 0;
  
  gtk_layout_set_size(GTK_LAYOUT(header), header->allocation.width, properties->headerRect.height);
  gtk_widget_set_size_request(header, 0, properties->headerRect.height);
}


void calculateHighlightBoxBorders(GdkRectangle *drawingRect, 
				  GdkRectangle *highlightRect,
				  GtkWidget *bigPicture,
				  const int yPadding)
{
  DEBUG_ENTER("calculateGridHighlightBoxBorders(grid)");
  
  /* Calculate how many pixels from the left edge of the widget to the first base in the range. Truncating
   * the double to an int after the multiplication means we can be up to 1 pixel out, but this should be fine. */
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (adjustment)
    {
      BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
      BlxViewContext *bc = bigPictureGetContext(bigPicture);
      
      /* Get the big picture display range in dna coords */
      IntRange bpRange;
      convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &bpRange);
      
      /* Get the detail view display range in dna coords */
      IntRange dvRange;
      convertDisplayRangeToDnaRange(detailViewGetDisplayRange(detailView), bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dvRange);
      
      /* Get the x coords for the start and end of the detail view display range */
      const int x1 = convertBaseIdxToRectPos(dvRange.min, drawingRect, &bpRange, TRUE, bc->displayRev, TRUE);
      const int x2 = convertBaseIdxToRectPos(dvRange.max + 1, drawingRect, &bpRange, TRUE, bc->displayRev, TRUE);
      
      highlightRect->x = min(x1, x2);
      highlightRect->y = 0;
      
      highlightRect->width = abs(x1 - x2);
      highlightRect->height = drawingRect->height + roundNearest(bpProperties->charHeight / 2.0) + yPadding + (2 * bpProperties->highlightBoxYPad);
    }
  
  DEBUG_EXIT("calculateGridHighlightBoxBorders returning");
}


/* Add the given child to the big picture. This is separated out into its
 * own function in case we want to change the container type in the future. */
static void addChildToBigPicture(GtkWidget *container, GtkWidget *child, gboolean expand)
{
  if (GTK_IS_BOX(container))
    {
      gtk_box_pack_start(GTK_BOX(container), child, expand, FALSE, 0);
    }
}


/* This function removes the grids from the big picture and re-adds them in the
 * correct order according to the displayRev flag. It should be called every
 * time the strands are toggled. It assumes the two grids are both already in the 
 * bigPicture container, and that the properties have been set for all 3 widgets. */
void refreshGridOrder(GtkWidget *bigPicture)
{
  GtkWidget *fwdStrandGrid = bigPictureGetFwdGrid(bigPicture);
  GtkWidget *revStrandGrid = bigPictureGetRevGrid(bigPicture);
  GtkWidget *fwdExonView = bigPictureGetFwdExonView(bigPicture);
  GtkWidget *revExonView = bigPictureGetRevExonView(bigPicture);
  
  /* Increase the reference count to make sure the widgets aren't destroyed when we remove them. */
  g_object_ref(fwdStrandGrid);
  g_object_ref(revStrandGrid);
  g_object_ref(fwdExonView);
  g_object_ref(revExonView);
  
  /* Find the direct container of the grids etc and remove them */
  GtkWidget *bpContainer = getNamedChildWidget(bigPicture, BIG_PICTURE_WIDGET_NAME);
  gtk_container_remove(GTK_CONTAINER(bpContainer), fwdStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bpContainer), revStrandGrid);
  gtk_container_remove(GTK_CONTAINER(bpContainer), fwdExonView);
  gtk_container_remove(GTK_CONTAINER(bpContainer), revExonView);
  
  /* Add them back, with the forward-strand grid at the top and the reverse-strand grid
   * at the bottom, or vice versa if the strands are toggled. */
  if (bigPictureGetDisplayRev(bigPicture))
    {
      addChildToBigPicture(bpContainer, revStrandGrid, FALSE);
      addChildToBigPicture(bpContainer, revExonView, FALSE);
      addChildToBigPicture(bpContainer, fwdExonView, FALSE);
      addChildToBigPicture(bpContainer, fwdStrandGrid, FALSE);
    }
  else
    {
      addChildToBigPicture(bpContainer, fwdStrandGrid, FALSE);
      addChildToBigPicture(bpContainer, fwdExonView, FALSE);
      addChildToBigPicture(bpContainer, revExonView, FALSE);
      addChildToBigPicture(bpContainer, revStrandGrid, FALSE);
    }
  
  /* Decrease the ref count again */
  g_object_unref(fwdStrandGrid);
  g_object_unref(revStrandGrid);
  g_object_unref(fwdExonView);
  g_object_unref(revExonView);
  
  /* Must show all child widgets because some of them may not have been in this parent before.
   * (Just calling gtk_widget_show on the individual trees doesn't seem to work.)
   * However, we then need to re-hide any that may have been previously hidden by the user. */
  gtk_widget_show_all(bigPicture);
  gtk_container_foreach(GTK_CONTAINER(bigPicture), hideUserHiddenWidget, NULL);  
}


/* This recalculates the highlight box positions for all relevant child
 * widgets of the big picture */
static void updateHighlightBox(GtkWidget *bigPicture, BigPictureProperties *properties)
{
  callFuncOnAllBigPictureGrids(bigPicture, (gpointer)calculateGridHighlightBoxBorders);
  callFuncOnAllBigPictureExonViews(bigPicture, (gpointer)calculateExonViewHighlightBoxBorders);
  calculateCoverageViewHighlightBoxBorders(properties->coverageView);
}


/* This should be called to do the required updates after the big picture range has changed */
static void onBigPictureRangeChanged(GtkWidget *bigPicture, BigPictureProperties *properties)
{
  /* Recalculate the exon view height, because it may have changed with more/less
   * exons being scrolled into view */
  calculateExonViewHeight(properties->fwdExonView);
  calculateExonViewHeight(properties->revExonView);

  /* We must force a resize, because the size-allocate signal does not
   * get emitted if the exon views have shrunk, only if they have expanded. */
  forceResize(bigPicture);
      
  /* Do a complete redraw */
  bigPictureRedrawAll(bigPicture);

  /* Refresh the dotter dialog, if it's open, because it may be tracking the big picture range */
  refreshDialog(BLXDIALOG_DOTTER, properties->blxWindow);
}


/* Set the display range for the big picture, based on the given width (i.e. number of
 * bases wide). Keeps the display centred on the same range that is shown in the detail view.
 * If recalcHighlightBox is true, the highlight box borders are recalculated. */
static void setBigPictureDisplayRange(GtkWidget *bigPicture, 
                                      BigPictureProperties *properties, 
                                      int width, 
                                      const gboolean keepCentered)
{
  DEBUG_ENTER("setBigPictureDisplayRange");

  GtkWidget *detailView = blxWindowGetDetailView(properties->blxWindow);
  IntRange *detailViewRange = detailViewGetDisplayRange(detailView);

  IntRange *displayRange = &properties->displayRange;
  IntRange *fullRange = blxWindowGetFullRange(properties->blxWindow);
  
  int detailViewWidth = detailViewRange->max - detailViewRange->min;
  int maxWidth = fullRange->max - fullRange->min;
  
  gboolean changedRange = FALSE;

  if (width < detailViewWidth)
    {
      /* Don't display less than the detail view range */
      displayRange->min = detailViewRange->min;
      displayRange->max = detailViewRange->max;
      width = displayRange->max - displayRange->min;
      changedRange = TRUE;
    }
  else if (width >= maxWidth)
    {
      /* Don't display more than the full range of the reference sequence */
      changedRange = displayRange->min != fullRange->min || displayRange->max != fullRange->max;
      
      displayRange->min = fullRange->min;
      displayRange->max = fullRange->max;
    }
  else if (keepCentered)
    {
      /* Keep the view centred on what's visible in the detail view */
      const int newcentre = getRangeCentre(detailViewRange);

      /* Try to display an equal amount either side of the centre */
      const int offset = roundNearest((double)width / 2.0);

      displayRange->min = newcentre - offset;
      displayRange->max = displayRange->min + width;
      
      changedRange = TRUE;
    }
  
  if (!changedRange) /* i.e. not already done */
    {
      /* See if we need to scroll the big picture range so that it completely 
       * contains the detail-view range.
       * NB We use a shorter big picture range so that we never get the highlight
       * box bumped right up against the edge of the big picture; UNLESS the
       * detail view range is larger than that limited range, in which case use
       * the original range. NOTE: disabled the border because it messes up the
       * initial big picture range when it is specified using the zoom-range arg */
      const gdouble borderFraction = 0.0;
      int border = getRangeLength(displayRange) * borderFraction;

      if (getRangeLength(detailViewRange) > width - (2 * border))
        border = 0;
      
      int offset = (displayRange->min + border) - detailViewRange->min;
      
      if (offset > 0)
        {
          displayRange->min -= offset;
          displayRange->max = displayRange->min + width;
          changedRange = TRUE;
        }
      
      offset = detailViewRange->max - (displayRange->max - border);
      
      if (offset > 0)
        {
          displayRange->max += offset;
          displayRange->min = displayRange->max - width;
          changedRange = TRUE;
        }
    }

  if (changedRange)
    {
      boundsLimitRange(displayRange, fullRange, TRUE);
      onBigPictureRangeChanged(bigPicture, properties);
    }
  
  updateHighlightBox(bigPicture, properties);
  bigPictureRefreshAll(bigPicture);

  DEBUG_EXIT("setBigPictureDisplayRange returning");
}


/* This function should be called every time the detail view display range has 
 * changed. It updates the position of the highlight box on the big picture and,
 * if necessary, scrolls to keep the highlight box in range. */
void refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean keepCentered)
{
  DEBUG_ENTER("refreshBigPictureDisplayRange");

  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  IntRange *bpRange = &properties->displayRange;
  
  GtkWidget *detailView = blxWindowGetDetailView(properties->blxWindow);
  IntRange *dvRange = detailViewGetDisplayRange(detailView);
  
  if (bpRange->min == UNSET_INT && bpRange->max == UNSET_INT) /* check both in case UNSET_INT is a real coord for one end! */
    {
      /* This is the first time we've refreshed the detail view range. Set the
       * initial big picture range width to be a ratio of the detail view range width. */
      const int width = (dvRange->max - dvRange->min) * bigPictureGetInitialZoom(bigPicture);
      setBigPictureDisplayRange(bigPicture, properties, width, TRUE);
    }
  else
    {
      /* Call set-width (but just use the existing width) to force the required updates. */
      const int width = bpRange->max - bpRange->min;
      setBigPictureDisplayRange(bigPicture, properties, width, keepCentered);
    }
  
  DEBUG_EXIT("refreshBigPictureDisplayRange returning");
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the 
 * left edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
int getLeftCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int leftCoord = centreCoord - roundNearest((double)width / 2.0);
  
  if (outerRect)
    {
      if (leftCoord < outerRect->x) 
	leftCoord = outerRect->x;
      else
	{
	  int leftCoordMax = outerRect->x + outerRect->width - width;
	  
	  if (leftCoord > leftCoordMax) 
	    {
	      leftCoord = leftCoordMax;
	    }
	}
    }
  
  return leftCoord;
}


/* Given the centre x coord of a rectangle and its width, find the x coord of the 
 * right edge. If an outer rectangle is given, limit the coord so that the
 * rectangle lies entirely within the outer rect. */
int getRightCoordFromCentre(const int centreCoord, const int width, const GdkRectangle *outerRect)
{
  int rightCoord = centreCoord + roundNearest((double)width / 2.0);
  
  if (outerRect)
    {
      if (rightCoord > outerRect->x + outerRect->width)
	rightCoord = outerRect->x + outerRect->width;
      else
	{
	  int rightCoordMin = outerRect->x + width;
	  
	  if (rightCoord < rightCoordMin) 
	    {
	      rightCoord = rightCoordMin;
	    }
	}
    }
  
  return rightCoord;
}


/* Zoom the big picture in/out. Just halves or doubles the current display range. */
void zoomBigPicture(GtkWidget *bigPicture, const gboolean zoomIn)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  IntRange *displayRange = &properties->displayRange;
  int newWidth = UNSET_INT;
  
  if (zoomIn)
    {
      newWidth = roundNearest(((double)(displayRange->max - displayRange->min)) / 2.0);
    }
  else
    {
      newWidth = (displayRange->max - displayRange->min) * 2;
    }

  setBigPictureDisplayRange(bigPicture, properties, newWidth, TRUE);
  calculateBigPictureCellSize(bigPicture, properties);
}


/* Zoom the big picture out to view the whole reference sequence */
void zoomWholeBigPicture(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  IntRange *displayRange = &properties->displayRange;
  IntRange *fullRange = blxWindowGetFullRange(properties->blxWindow);
  
  /* Check we're not already showing the whole range */
  if (displayRange->min != fullRange->min || displayRange->max != fullRange->max)
    {
      setBigPictureDisplayRange(bigPicture, properties, fullRange->max - fullRange->min, TRUE);
      calculateBigPictureCellSize(bigPicture, properties);
    }
}


/* Calculate the number of cells to show vertically in the grids, and update the
 * grid's percent-ID range to be a full number of cells. */
void calculateNumVCells(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (properties->idPerCell == 0.0)
    {
      g_warning("%%ID per cell setting is 0. Cannot calculate number of cells.\n");
      return;
    }
  
  const double idRangeLen = properties->percentIdRange.max - properties->percentIdRange.min;
  properties->numVCells = ceil(idRangeLen / properties->idPerCell); 

  if (properties->numVCells < MIN_NUM_V_CELLS)
    {
      properties->numVCells = MIN_NUM_V_CELLS;
    }
  
  properties->percentIdRange.min = properties->percentIdRange.max - (properties->numVCells * properties->idPerCell);
}


/* update function to be called when the percent ID range or percent-ID-per-cell
 * values have been changed */
static void updateOnPercentIdChanged(GtkWidget *bigPicture)
{
  if (bigPicture)
    {
      BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
      if (properties)
        {
          calculateNumVCells(bigPicture);
  
          callFuncOnAllBigPictureGrids(bigPicture, (gpointer)calculateGridBorders);
          callFuncOnAllBigPictureGrids(bigPicture, (gpointer)calculateGridHighlightBoxBorders);
          calculateCoverageViewBorders(properties->coverageView);
  
          bigPictureRedrawAll(bigPicture);
        }
    }
}


/* Prepare the big picture for printing - this draws normally-transient 
 * components onto the cached drawable so that they get included in the print.
 * bigPictureRedrawAll should be called afterwards to remove the transient 
 * components. */
void bigPicturePrepareForPrinting(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  callFuncOnAllBigPictureGrids(bigPicture, (gpointer)gridPrepareForPrinting);
  callFuncOnAllBigPictureExonViews(bigPicture, (gpointer)exonViewPrepareForPrinting);
  coverageViewPrepareForPrinting(properties->coverageView);
}


/***********************************************************
 *			    Events			   *
 ***********************************************************/

static void onZoomInBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomBigPicture(bigPicture, TRUE);
}


static void onZoomOutBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomBigPicture(bigPicture, FALSE);
}


static void onZoomWholeBigPicture(GtkButton *button, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);
  zoomWholeBigPicture(bigPicture);
}


static gboolean onExposeGridHeader(GtkWidget *header, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(header)->bin_window;
  
  if (window)
    {
      /* Just push the stored bitmap onto the screen */
      GdkDrawable *bitmap = widgetGetDrawable(header);

      if (!bitmap)
        {
          /* If the cache is empty, create it now. */
	  redrawBigPictureGridHeader(header);
	  bitmap = widgetGetDrawable(header);
        }
      
      if (bitmap)
        {
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
    }
  
  return TRUE;
}


/* Scroll the big picture left by one increment */
void scrollBigPictureLeftStep(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  
  IntRange *displayRange = &properties->displayRange;
  
  if (displayRange->min > bc->fullDisplayRange.min)
    {
      /* Check we can scroll the full increment amount. If not, scroll to the end of the full range */
      int diff = getRangeLength(displayRange) / BLX_SCROLL_INCREMENT_RATIO;

      if (displayRange->min - diff < bc->fullDisplayRange.min)
        diff = displayRange->min - bc->fullDisplayRange.min;

      /* Update the range */
      displayRange->min -= diff;
      displayRange->max -= diff;

      /* Update */
      onBigPictureRangeChanged(bigPicture, properties);
      updateHighlightBox(bigPicture, properties);
      bigPictureRefreshAll(bigPicture);

      /* Scroll the detail view too if necessary to keep it visible */
      detailViewScrollToKeepInRange(bigPictureGetDetailView(bigPicture), displayRange);
    }
}


/* Scroll the big picture right by one increment */
void scrollBigPictureRightStep(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  
  IntRange *displayRange = &properties->displayRange;
  
  /* Check we're not already at the max of the full range. */
  if (displayRange->max < bc->fullDisplayRange.max)
    {
      /* Check we can scroll the full increment amount. If not, scroll to the end of the full range */
      int diff = getRangeLength(displayRange) / BLX_SCROLL_INCREMENT_RATIO;

      if (displayRange->max + diff > bc->fullDisplayRange.max)
        diff = bc->fullDisplayRange.max - displayRange->max;

      /* Adjust the range */
      displayRange->min += diff;
      displayRange->max += diff;

      /* Update */
      onBigPictureRangeChanged(bigPicture, properties);
      updateHighlightBox(bigPicture, properties);
      bigPictureRefreshAll(bigPicture);

      /* Scroll the detail view too if necessary to keep it visible */
      detailViewScrollToKeepInRange(bigPictureGetDetailView(bigPicture), displayRange);
    }
}


/* Convert an x coord in the given rectangle to a base index (in nucleotide coords) */
static gint convertRectPosToBaseIdx(const gint x, 
                                    const GdkRectangle* const displayRect,  
                                    const IntRange* const dnaDispRange,
                                    const gboolean displayRev)
{
  gint result = UNSET_INT;
  
  gdouble distFromEdge = (gdouble)(x - displayRect->x);
  int basesFromEdge = (int)(distFromEdge / pixelsPerBase(displayRect->width, dnaDispRange));

  if (displayRev)
    {
      result = dnaDispRange->max - basesFromEdge;
    }
  else
    {
      result = dnaDispRange->min + basesFromEdge;
    }
  
  return result;
}


/* Draw the preview box on the given drawable within the boundaries of the given displayRect.
 * The boundaries of the preview box are given by highlightRect.
 * Only does anything if the preview box centre is set. */
void drawPreviewBox(GtkWidget *bigPicture, 
                    GdkDrawable *drawable, 
                    GdkRectangle *displayRect, 
                    GdkRectangle *highlightRect)
{
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  
  if (!bpProperties->displayPreviewBox)
    {
      return;
    }

  /* Get the display range in dna coords */
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dnaDispRange);
  
  /* Find the x coord for the left edge of the preview box (or the right edge, if
   * the display is right-to-left). */
  int x = getLeftCoordFromCentre(bpProperties->previewBoxCentre, highlightRect->width, displayRect);
  
  /* Convert it to the base index and back again so that we get it rounded to the position of
   * the nearest base. */
  int baseIdx = convertRectPosToBaseIdx(x, displayRect, &dnaDispRange, bc->displayRev);
  int xRounded = convertBaseIdxToRectPos(baseIdx, displayRect, &dnaDispRange, TRUE, bc->displayRev, TRUE);
  
  /* The other dimensions of the preview box are the same as the current highlight box. */
  GdkRectangle previewRect = {xRounded, highlightRect->y, highlightRect->width, highlightRect->height};

  GdkColor *previewBoxColor = getGdkColor(BLXCOLOR_PREVIEW_BOX, bc->defaultColors, FALSE, bc->usePrintColors);

  drawHighlightBox(drawable, &previewRect, bpProperties->previewBoxLineWidth, previewBoxColor);
}


/* Show a preview box centred on the given x coord. If init is true, we're initialising
 * a drag; otherwise, we're already dragging */
void showPreviewBox(GtkWidget *bigPicture, const int x, const gboolean init, const int offset)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (init)
    {
      properties->displayPreviewBox = TRUE;
      properties->previewBoxOffset = offset;

      static GdkCursor *cursor = NULL;
      cursor = gdk_cursor_new(GDK_FLEUR);

      GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
      gdk_window_set_cursor(blxWindow->window, cursor);
    }

  /* We might get called by a drag operation where the preview box drag has not been
   * initialised; just ignore it */
  if (!properties->displayPreviewBox)
    return;

  /* Whether dragging or initialising, we need to update the position */
  properties->previewBoxCentre = x + properties->previewBoxOffset;

  /* Refresh all child widgets, and also the coverage view (which may not
   * be a child of the big picture) */
  gtk_widget_queue_draw(bigPicture);
  gtk_widget_queue_draw(properties->coverageView);
}


/* Scroll the big picture so that it is centred on the current preview box position, and clear
 * the preview box. (Do nothing if preview box is not currently shown.)  */
void acceptAndClearPreviewBox(GtkWidget *bigPicture, const int xCentreIn, GdkRectangle *displayRect, GdkRectangle *highlightRect)
{
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);

  gdk_window_set_cursor(blxWindow->window, NULL);

  if (!bpProperties->displayPreviewBox)
    return;

  /* Apply any offset required to get the real centre coord */
  const int xCentre = xCentreIn + bpProperties->previewBoxOffset;
  
  /* Get the display range in dna coords */
  BlxViewContext *bc = bigPictureGetContext(bigPicture);
  IntRange dnaDispRange;
  convertDisplayRangeToDnaRange(&bpProperties->displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dnaDispRange);
  
  /* Find the base index where the new scroll range will start. This is the leftmost
   * edge of the preview box if numbers increase in the normal left-to-right direction, 
   * or the rightmost edge if the display is reversed. */
  const int x = getLeftCoordFromCentre(xCentre, highlightRect->width, displayRect);
  int baseIdx = convertRectPosToBaseIdx(x, displayRect, &dnaDispRange, bc->displayRev);
  
  /* Subtract 1 if the display is reversed to give the base to the right of x, rather than the base to the left of x */
  if (bc->displayRev)
    --baseIdx;
  
  /* Reset the preview box status. We don't need to un-draw it because we
   * will redraw the whole big picture below. */
  bpProperties->displayPreviewBox = FALSE;
  bpProperties->previewBoxOffset = 0;
  
  /* Update the detail view's scroll pos to start at the new base. The base index is in terms of
   * the nucleotide coords so we need to convert to display coords */
  const int displayIdx = convertDnaIdxToDisplayIdx(baseIdx, bc->seqType, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  
  GtkWidget *detailView = bigPictureGetDetailView(bigPicture);
  setDetailViewStartIdx(detailView, displayIdx, bc->seqType);
  
  /* Re-centre the big picture */
  refreshBigPictureDisplayRange(bigPicture, TRUE);
  
  gtk_widget_queue_draw(bigPicture);
}


/* Recursively loop through the children of the given widget and sum the heights of any
 * that are grid or exon-view widgets. */
static int getBigPictureChildrenHeights(GtkWidget *widget, const int heightIn)
{
  GList *children = gtk_container_get_children(GTK_CONTAINER(widget));
  GList *child = children;
  int height = heightIn;
  
  for ( ; child; child = child->next)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);

      if (GTK_WIDGET_VISIBLE(childWidget))
	{
	  const gchar *name = gtk_widget_get_name(childWidget);
      
	  if (stringsEqual(name, BIG_PICTURE_GRID_NAME, TRUE) || 
	      stringsEqual(name, BIG_PICTURE_EXON_VIEW_NAME, TRUE) ||
	      stringsEqual(name, BIG_PICTURE_GRID_HEADER_NAME, TRUE))
	    {
	      height += childWidget->allocation.height;
	    }
	  else if (GTK_IS_CONTAINER(childWidget))
	    {
	      height += getBigPictureChildrenHeights(childWidget, height);
	    }
	}
    }

  g_list_free(children);
    
  return height;
}


/* Calculate the size of the big picture's children and set its size request to fit them in
 * but don't go greater than the maximum. This is called when the exon view or grids change size
 * etc. to make sure we only use as much space as necessary. Note that because the big picture
 * is in a paned window this does not have an effect if the user has changed (i.e. specifically
 * set) the pane size, which makes sense but it would be nice to have some way to be able to 
 * revert to the original behaviour */
static void bigPictureRecalculateSize(GtkWidget *bigPicture)
{
  DEBUG_ENTER("bigPictureRecalculateSize");

  int height = getBigPictureChildrenHeights(bigPicture, 0);
  height += 4; /* not sure where this extra space comes from (padding or something?) */

  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  int maxHeight = blxWindow->allocation.height * MAX_BIG_PICTURE_HEIGHT_RATIO;
  
  height = min(height, maxHeight);
  gtk_widget_set_size_request(bigPicture, -1, height);
  
  DEBUG_EXIT("bigPictureRecalculateSize returning");
}


/* Callback called when the big picture widget's size has changed */
static void onSizeAllocateBigPicture(GtkWidget *bigPicture, GtkAllocation *allocation, gpointer data)
{
  DEBUG_ENTER("onSizeAllocateBigPicture");

  /* Recalculate the widget size based on its child widget sizes. Note that we
   * don't need to do anything else here because onScrollPosChanged will be 
   * called, which does all the updating of the big picture range. */
  bigPictureRecalculateSize(bigPicture);
  
  bigPictureRedrawAll(bigPicture);
  
  DEBUG_EXIT("onSizeAllocateBigPicture returning");
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

BigPictureProperties* bigPictureGetProperties(GtkWidget *bigPicture)
{
  /* optimisation: cache result, because we know there is only ever one big picture */
  static BigPictureProperties *properties = NULL;
  
  if (!properties && bigPicture)
    properties = (BigPictureProperties*)(g_object_get_data(G_OBJECT(bigPicture), "BigPictureProperties"));
  
  return properties;
}

BlxViewContext* bigPictureGetContext(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetContext(blxWindow);
}

static GridHeaderProperties* gridHeaderGetProperties(GtkWidget *gridHeader)
{
  /* optimisation: cache result, because we know there is only ever one grid header */
  static GridHeaderProperties *properties = NULL;
  
  if (!properties && gridHeader)
    properties = (GridHeaderProperties*)(g_object_get_data(G_OBJECT(gridHeader), "GridHeaderProperties"));
  
  return properties;
}


static void onDestroyBigPicture(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (properties)
    {
      if (properties->roundValues)
	{
	  g_slist_free(properties->roundValues);
	  properties->roundValues = NULL;
	}
      
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", NULL);
    }
}


static void bigPictureCreateProperties(GtkWidget *bigPicture, 
				       GtkWidget *blxWindow, 
                                       GtkWidget *coverageView, 
				       GtkWidget *header, 
				       GtkWidget *fwdStrandGrid,
				       GtkWidget *revStrandGrid,
				       GtkWidget *fwdExonView,
				       GtkWidget *revExonView,
				       int previewBoxCentre,
                                       const IntRange* const initRange,
                                       const IntRange* const fullRange,
				       const int initialZoom,
				       const gdouble lowestId)
{
  if (bigPicture)
    { 
      BigPictureProperties *properties = (BigPictureProperties*)g_malloc(sizeof *properties);
      
      properties->blxWindow = blxWindow;
      properties->header = header;
      properties->coverageView = coverageView;
      properties->fwdStrandGrid = fwdStrandGrid;
      properties->revStrandGrid = revStrandGrid;
      properties->fwdExonView = fwdExonView;
      properties->revExonView = revExonView;
      
      properties->numHCells = UNSET_INT;
      properties->basesPerCell = UNSET_INT;
      properties->roundTo = 25;
      
      properties->numVCells = UNSET_INT;
      properties->idPerCell = DEFAULT_PERCENT_ID_PER_CELL;
      properties->percentIdRange.min = lowestId;
      properties->percentIdRange.max = (gdouble)DEFAULT_GRID_PERCENT_ID_MAX;

      properties->displayPreviewBox = FALSE;
      properties->previewBoxOffset = 0;
      properties->previewBoxCentre = previewBoxCentre;
      properties->leftBorderChars = numDigitsInInt(DEFAULT_GRID_PERCENT_ID_MAX) + 2; /* Extra fudge factor because char width is approx */
      properties->highlightBoxMinWidth = MIN_HIGHLIGHT_BOX_WIDTH;
      properties->previewBoxLineWidth = DEFAULT_PREVIEW_BOX_LINE_WIDTH;
      properties->highlightBoxYPad = DEFAULT_HIGHLIGHT_BOX_Y_PAD;
      properties->initialZoom = initialZoom;
      
      /* Set the initial big picture range from that requested. Note that the
       * intput range may be UNSE_INTs, in which case these values will be
       * initialized when the detail view size is first allocated, based on
       * the "initial zoom level" which is a ratio of the detail view range. */
      properties->displayRange.min = initRange->min;
      properties->displayRange.max = initRange->max;

      /* Calculate the font size */
      getFontCharSize(bigPicture, bigPicture->style->font_desc, &properties->charWidth, &properties->charHeight);
      
      /* Create the list of "nice round values" to round the grid header values to.
       * Create the list in reverse order (i.e. highest values first). */
      properties->roundValues = NULL;
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(25));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(50));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(100));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(125));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(150));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(250));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(500));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(1000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(2500));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(5000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(10000));
      properties->roundValues = g_slist_prepend(properties->roundValues, GINT_TO_POINTER(25000));
      
      g_object_set_data(G_OBJECT(bigPicture), "BigPictureProperties", properties);
      g_signal_connect(G_OBJECT(bigPicture), "destroy", G_CALLBACK(onDestroyBigPicture), NULL); 
    }
}


static void onDestroyGridHeader(GtkWidget *bigPicture)
{
  GridHeaderProperties *properties = gridHeaderGetProperties(bigPicture);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(bigPicture), "GridHeaderProperties", NULL);
    }
}

static void gridHeaderCreateProperties(GtkWidget *gridHeader, GtkWidget *bigPicture, GtkWidget *refButton)
{
  if (gridHeader)
    {
      GridHeaderProperties *properties =(GridHeaderProperties*)g_malloc(sizeof *properties);
      
      properties->bigPicture = bigPicture;
      properties->refButton = refButton;
      properties->numHeaderLines = DEFAULT_GRID_NUM_HEADER_LINES;
      properties->headerYPad = DEFAULT_GRID_HEADER_Y_PAD;
      
      g_object_set_data(G_OBJECT(gridHeader), "GridHeaderProperties", properties);
      g_signal_connect(G_OBJECT(gridHeader), "destroy", G_CALLBACK(onDestroyGridHeader), NULL);
    }
}




GtkWidget* bigPictureGetBlxWindow(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->blxWindow : NULL;
}

GtkWidget* bigPictureGetDetailView(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetDetailView(blxWindow);
}


GtkWidget* bigPictureGetFwdGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->fwdStrandGrid : NULL;
}

GtkWidget* bigPictureGetRevGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->revStrandGrid : NULL;
}

/* Get the active grid (forward strand grid by default, reverse strand grid if display toggled) */
GtkWidget* bigPictureGetActiveGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->revStrandGrid : properties->fwdStrandGrid;
}

/* Get the in-active grid (reverse strand grid by default, forward strand grid if display toggled) */
GtkWidget* bigPictureGetInactiveGrid(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->fwdStrandGrid : properties->revStrandGrid;
}

GtkWidget* bigPictureGetFwdExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->fwdExonView : NULL;
}

GtkWidget* bigPictureGetRevExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->revExonView : NULL;
}

GtkWidget* bigPictureGetCoverageView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties ? properties->coverageView : NULL;
}

GtkWidget* bigPictureGetActiveExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->revExonView : properties->fwdExonView;
}

GtkWidget* bigPictureGetInactiveExonView(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return bigPictureGetDisplayRev(bigPicture) ? properties->fwdExonView : properties->revExonView;
}

gboolean bigPictureGetDisplayRev(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindow ? blxWindowGetDisplayRev(blxWindow) : FALSE;
}

IntRange* bigPictureGetDisplayRange(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->displayRange;
}

static int bigPictureGetInitialZoom(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->initialZoom;
}

GtkWidget* bigPictureGetGridHeader(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->header;
}

BlxSeqType bigPictureGetSeqType(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetSeqType(blxWindow);
}

int bigPictureGetNumFrames(GtkWidget *bigPicture)
{
  GtkWidget *blxWindow = bigPictureGetBlxWindow(bigPicture);
  return blxWindowGetNumFrames(blxWindow);
}

gdouble bigPictureGetIdPerCell(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->idPerCell;
}

gboolean bigPictureSetIdPerCell(GtkWidget *bigPicture, const gdouble idPerCell)
{
  gboolean result = FALSE;
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (idPerCell < GRID_SCALE_MIN_ID_PER_CELL)
    {
      g_critical("Cannot set ID per cell less than %1.1f.\n", GRID_SCALE_MIN_ID_PER_CELL);
    }
  else
    {
      properties->idPerCell = idPerCell;
      updateOnPercentIdChanged(bigPicture);
      result = TRUE;
    }
  
  return result;
}

int bigPictureGetNumVCells(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return properties->numVCells;
}

DoubleRange* bigPictureGetPercentIdRange(GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  return &properties->percentIdRange;
}

gboolean bigPictureSetMaxPercentId(GtkWidget *bigPicture, const gdouble newValue)
{
  gboolean result = FALSE;
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  if (newValue < GRID_SCALE_MIN)
    {
      g_critical("Cannot set grid scale less than %d.\n", GRID_SCALE_MIN);
    }
  else if (newValue > GRID_SCALE_MAX)
    {
      g_critical("Cannot set grid scale greater than %d.\n", GRID_SCALE_MAX);
    }
  else if (newValue < properties->percentIdRange.min)
    {
      g_critical("Cannot set grid maximum to %1.1f because this is less than the grid minimum of %1.1f.\n", (double)newValue, (double)properties->percentIdRange.min);
    }
  else
    {
      properties->percentIdRange.max = newValue;
      updateOnPercentIdChanged(bigPicture);
      result = TRUE;
    }
  
  return result;
}

gboolean bigPictureSetMinPercentId(GtkWidget *bigPicture, const gdouble newValue)
{
  gboolean result = FALSE;
  
  if (bigPicture && newValue != -1) /* -1 means unset */
    {
      BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
      
      if (properties)
        {
          if (newValue < GRID_SCALE_MIN)
            {
              g_critical("Cannot set grid scale less than %d.\n", GRID_SCALE_MIN);
            }
          else if (newValue > GRID_SCALE_MAX)
            {
              g_critical("Cannot set grid scale greater than %d.\n", GRID_SCALE_MAX);
            }
          else if (newValue > properties->percentIdRange.max)
            {
              g_critical("Cannot set grid minimum to %1.1f because this is greater than the grid maximum of %1.1f.\n", (double)newValue, (double)properties->percentIdRange.max);
            }
          else
            {    
              properties->percentIdRange.min = newValue;
              updateOnPercentIdChanged(bigPicture);
              result = TRUE;
            }
        }
    }

  return result;
}

/***********************************************************
 *                     Initialization                      *
 ***********************************************************/

/* Create a tool button for the big picture header */
static GtkWidget* createButton(GtkWidget *container, const char *label, const char *tooltip, GtkSignalFunc callback_func, gpointer data)
{
  GtkWidget *eventBox = gtk_event_box_new();
  gtk_box_pack_start(GTK_BOX(container), eventBox, FALSE, FALSE, 0);
  
  GtkWidget *button = gtk_button_new_with_label(label);
  gtk_container_add(GTK_CONTAINER(eventBox), button);
  
  /* gtk_widget_set_tooltip_text(button, tooltip); */ /* not in pre-GTK-2.12 */
  
  g_signal_connect(GTK_OBJECT(button), "clicked", callback_func, data);
  
  return button;
}


/* Create the header for the big picture section. This is a custom header
 * that contains some tool buttons and the header labels for the grid. We
 * squash them all into the same horizontal section to save space, and also
 * so that we can hide any of the grids but still keep the header visible. */
static GtkWidget *createBigPictureGridHeader(GtkWidget *bigPicture)
{
  GtkWidget *header = gtk_layout_new(NULL, NULL);
  gtk_widget_set_name(header, BIG_PICTURE_GRID_HEADER_NAME);
  
  /* Create the header buttons. Put them in an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_layout_put(GTK_LAYOUT(header), hbox, 0, 0);
  
  /* Store a ref to the first button so we can find out the height of the
   * buttons when we come to calculate the grid borders */
  GtkWidget *refButton = createButton(hbox, "Zoom in", "Zoom in (Ctrl-=)", (GtkSignalFunc)onZoomInBigPicture, bigPicture);
  createButton(hbox, "Zoom out", "Zoom out (Ctrl--)", (GtkSignalFunc)onZoomOutBigPicture, bigPicture);
  createButton(hbox, "Whole", "Zoom out to whole width (Shift-Ctrl--_", (GtkSignalFunc)onZoomWholeBigPicture, bigPicture);
  
  /* Create the header properties */
  gridHeaderCreateProperties(header, bigPicture, refButton);
  
  /* Conect signals */
  g_signal_connect(G_OBJECT(header), "expose-event", G_CALLBACK(onExposeGridHeader), NULL);
  
  return header;
}


GtkWidget* createBigPicture(GtkWidget *blxWindow, 
			    GtkContainer *parent,
                            GtkWidget *coverageView,
			    GtkWidget **fwdStrandGrid, 
			    GtkWidget **revStrandGrid,
                            const IntRange* const initRange,
                            const IntRange* const fullRange,
			    const int initialZoom,
			    const gdouble lowestId)
{
  /* Create the main big picture widget, which will contain all of the 
   * individual big-picture grids, plus a header, in a scrollable vbox. */
  GtkWidget *bigPicture = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(bigPicture), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_container_add(parent, bigPicture);
   
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(bigPicture), vbox);
  gtk_widget_set_name(vbox, BIG_PICTURE_WIDGET_NAME);
  
  /* Our big picture needs to have a header plus 3 panes:
   * 1. the top grid (fwd strand);
   * 2. the exon view; and
   * 3. bottom grid (reverse strand) */
  
  GtkWidget *header = createBigPictureGridHeader(bigPicture);
  addChildToBigPicture(vbox, header, FALSE);

  *fwdStrandGrid = createBigPictureGrid(bigPicture, BLXSTRAND_FORWARD);
  *revStrandGrid = createBigPictureGrid(bigPicture, BLXSTRAND_REVERSE);

  GtkWidget *fwdExonView = createExonView(bigPicture, BLXSTRAND_FORWARD);
  GtkWidget *revExonView = createExonView(bigPicture, BLXSTRAND_REVERSE);
  
  /* By default, make the forward strand the top grid */
  addChildToBigPicture(vbox, *fwdStrandGrid, FALSE);
  addChildToBigPicture(vbox, fwdExonView, FALSE);
  addChildToBigPicture(vbox, revExonView, FALSE);
  addChildToBigPicture(vbox, *revStrandGrid, FALSE);

  g_signal_connect(G_OBJECT(bigPicture), "size-allocate", G_CALLBACK(onSizeAllocateBigPicture), NULL);

  /* Set the big picture properties */
  bigPictureCreateProperties(bigPicture, 
			     blxWindow,
                             coverageView,
			     header, 
			     *fwdStrandGrid,
			     *revStrandGrid,
			     fwdExonView,
			     revExonView,
			     UNSET_INT,
                             initRange,
                             fullRange,
			     initialZoom,
			     lowestId);
  
  return bigPicture;
}



