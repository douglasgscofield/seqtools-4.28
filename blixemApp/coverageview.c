/*  File: coverageview.c
 *  Author: Gemma Barson, 2011-03-21
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
 * Description: This view shows a plot of the number of alignments at each
 *              coordinate in the reference sequence.
 *----------------------------------------------------------------------------
 */

#include "blixemApp/coverageview.h"
#include "blixemApp/bigpicture.h"
#include "blixemApp/detailview.h"
#include "blixemApp/blixem_.h"
#include <gtk/gtk.h>
#include <math.h>


#define DEFAULT_COVERAGE_VIEW_Y_PADDING		3	  /* this provides space between the drawing area and the edge of the widget */
#define DEFAULT_NUM_V_CELLS			4	  /* number of vertical cells to show on the grid */
#define MIN_LINE_WIDTH				0.5	  /* this provides space between the drawing area and the edge of the widget */
#define COVERAGE_VIEW_NAME                      "CoverageView"


typedef struct _CoverageViewProperties
  {
    GtkWidget *blxWindow;   /* The main blixem window */

    int viewYPadding;	     /* The y padding around the view rect */
    double numVCells;	     /* The number of cells to show vertically */
    gdouble rangePerCell;    /* The range of depth values shown per grid cell on the plot */
    
    GdkRectangle viewRect;   /* The rectangle we draw in */
    GdkRectangle displayRect; /* The total display area */
    GdkRectangle highlightRect; /* The area that is highlighted (which indicates the detail-view range) */
  } CoverageViewProperties;




/***********************************************************
 *                       Properties                        *
 ***********************************************************/

CoverageViewProperties* coverageViewGetProperties(GtkWidget *widget)
{
  /* optimisation: cache result, because we know there is only ever one coverage view */
  static CoverageViewProperties *properties = NULL;
  
  if (!properties && widget)
    properties = (CoverageViewProperties*)(g_object_get_data(G_OBJECT(widget), "CoverageViewProperties"));
  
  return properties;
}

static void onDestroyCoverageView(GtkWidget *widget)
{
  CoverageViewProperties *properties = coverageViewGetProperties(widget);
  
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", NULL);
    }
}

static void coverageViewCreateProperties(GtkWidget *widget, 
                                         GtkWidget *blxWindow,
					 BlxViewContext *bc)
{
  if (widget)
    { 
      CoverageViewProperties *properties = (CoverageViewProperties*)g_malloc(sizeof *properties);
      
      properties->blxWindow = blxWindow;
      properties->viewYPadding = DEFAULT_COVERAGE_VIEW_Y_PADDING;
      properties->numVCells = DEFAULT_NUM_V_CELLS;
      properties->rangePerCell = 0;
      
      g_object_set_data(G_OBJECT(widget), "CoverageViewProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCoverageView), NULL); 
    }
}


/* This function should be called whenever the coverage depth data has changed */
void updateCoverageDepth(GtkWidget *coverageView, BlxViewContext *bc)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);

  /* Set up a list of 'nice' values to round to for displaying labels */
  static GSList *roundValues = NULL;
  
  if (!roundValues)
    {
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(200));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2500));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(250000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(1000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(2500000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(5000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(10000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(25000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(50000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(100000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(250000000));
      roundValues = g_slist_prepend(roundValues, GINT_TO_POINTER(500000000));

      /* First time round, calculate the range per cell, aiming for 
       * around 5 cells. (If we enter this function again, it's because the
       * user has manually entered the range per cell so we just need to calculate
       * the relevant number of cells) */
      properties->numVCells = 5;
      properties->rangePerCell = ceil((gdouble)bc->maxDepth / (gdouble)properties->numVCells);

      /* Round the result and recalculate the number of cells */
      properties->rangePerCell = roundUpToValueFromList(properties->rangePerCell, roundValues, NULL);
      
      if (properties->rangePerCell < 1)
        properties->rangePerCell = 1;
    }
  
  properties->numVCells = (gdouble)bc->maxDepth / properties->rangePerCell;
  
  coverageViewRecalculate(coverageView);
}


static GtkWidget *coverageViewGetBigPicture(GtkWidget *coverageView)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  return (properties ? blxWindowGetBigPicture(properties->blxWindow) : NULL);
}

double coverageViewGetDepthPerCell(GtkWidget *coverageView)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  return properties ? properties->rangePerCell : 0.0;
}

gboolean coverageViewSetDepthPerCell(GtkWidget *coverageView, const double depthPerCell)
{
  if (depthPerCell <= 0.0)
    return FALSE;
  
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  
  if (properties)
    {
      properties->rangePerCell = depthPerCell;
      BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
      updateCoverageDepth(coverageView, bc);
    }
  
  return TRUE;
}

/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Clear the cached drawable and re-draw the coverage view */
void coverageViewRedraw(GtkWidget *coverageView)
{
  widgetClearCachedDrawable(coverageView, NULL);
  gtk_widget_queue_draw(coverageView);
}


/* Recalculate the size of the coverage view widget and redraw */
void coverageViewRecalculate(GtkWidget *coverageView)
{
  calculateCoverageViewBorders(coverageView);
  coverageViewRedraw(coverageView);
}


/* Draw a single bar of the coverage bar graph */
static void drawCoverageBar(const double x1,
			    const double x2,
			    double y1,
			    const double y2,
			    cairo_t *cr)
{
  double height = y2 - y1;
  double width = x2 - x1;

  if (height <= MIN_LINE_WIDTH)
    {
      height = MIN_LINE_WIDTH;
      y1 = y2 - height;
    }
  
  if (width <= 0)
    return;
  
  cairo_rectangle(cr, x1 - MIN_LINE_WIDTH, y1, width + (2 * MIN_LINE_WIDTH), height);
  cairo_fill(cr);
}


/* Utility to get the max depth that the coverage view shows, based on 
 * the number of cells and the range per cell. Note that this returns the
 * max lable value, i.e. the value of the top gridline; the real max
 * depth may be slightly greater than this, and may extend above the top
 * gridline (the height of the widget is made big enough to accommodate this). */
static int coverageViewGetMaxLabeledDepth(CoverageViewProperties *properties)
{
  /* to do: ideally we would round numcells to the nearest int rather than
   * truncating, but there is a bug with that where the horizontal grid lines
   * are sometimes not drawn with the correct labels */
  int numCells = (int)(properties->numVCells);
  const int result = properties->rangePerCell * numCells;
  return result;
}


/* Draw the actual coverage data as a bar chart */
static void drawCoveragePlot(GtkWidget *coverageView, GdkDrawable *drawable)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  
  if (!bc || bc->maxDepth <= 0)
    return;

  cairo_t *cr = gdk_cairo_create(drawable);

  const GdkColor *color = getGdkColor(BLXCOLOR_COVERAGE_PLOT, bc->defaultColors, FALSE, bc->usePrintColors);
  gdk_cairo_set_source_color(cr, color);
  
  const int maxDepth = coverageViewGetMaxLabeledDepth(properties);
  const double pixelsPerVal = (double)properties->viewRect.height / (double)maxDepth;
  const int bottomBorder = properties->viewRect.y + properties->viewRect.height;
  
  /* Loop through each coord in the display range */
  const IntRange* const displayRange = bigPictureGetDisplayRange(bigPicture);
  
  double startX = -1.0;
  double prevX = -1.0;
  double prevY = -1.0;
  int coord = displayRange->min;
  
  for ( ; coord <= displayRange->max; ++coord)
    {
      /* Get the x position for this coord (always pass displayRev as false because
       * display coords are already inverted if the display is reversed). */
      const double x = convertBaseIdxToRectPos(coord, &properties->viewRect, displayRange, TRUE, FALSE, TRUE);
      
      /* Convert the display coord to a zero-based coord in the full ref seq
       * display range, for indexing the depth array. Note that we need to
       * un-invert the display coord if the display is reversed. */
      const int idx = invertCoord(coord, &bc->fullDisplayRange, bc->displayRev);
      const int depth = bc->depthArray[idx - bc->fullDisplayRange.min];

      /* Calculate the y position based on the depth */
      const double height = (pixelsPerVal * (double)depth);
      const double y = (double)bottomBorder - height;

      /* First time round, don't draw the line - just get the starting x and y.
       * Also, don't draw the line if y is the same as prevY; drawing is quite
       * expensive, so if we have multiple y values that are the same we remember
       * the first x coord at this y value and draw from there all in one go. */
      if (prevX == -1.0 && prevY == -1.0)
        {
          startX = x;
        }
      else if (y != prevY || coord == displayRange->max)
        {
          /* If we had multiple positions where y was the same, draw a horizontal
           * line at that y position. */
          if (prevX != startX)
	    drawCoverageBar(startX, prevX, prevY, bottomBorder, cr);
	    
          /* Now draw the sloped line from the previous y to the new y. */
	  drawCoverageBar(prevX, x, y, bottomBorder, cr);

          /* Reset the starting point */
          startX = x;
        }

      prevX = x;
      prevY = y;
    }
  
  cairo_destroy(cr);
}


/* Main function for drawing the coverage view */
static void drawCoverageView(GtkWidget *coverageView, GdkDrawable *drawable)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);

  drawVerticalGridLines(&properties->viewRect, &properties->highlightRect, 
			properties->viewYPadding, bc, bpProperties, drawable);
  
  const int maxDepth = coverageViewGetMaxLabeledDepth(properties);
  
  drawHorizontalGridLines(coverageView, bigPicture, &properties->viewRect, bc, bpProperties, drawable,
			  (int)(properties->numVCells), properties->rangePerCell, (gdouble)maxDepth, TRUE, "");
  
  drawCoveragePlot(coverageView, drawable);
}


/* Calculate the borders of the highlight box (the shaded area that indicates the
 * detail-view range). (This is just a convenience way to call calculateHighlightBoxBorders
 * from an external function.) */
void calculateCoverageViewHighlightBoxBorders(GtkWidget *coverageView)
{  
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);

  calculateHighlightBoxBorders(&properties->displayRect, &properties->highlightRect, bigPicture, 0);
}


/* Calculate the borders of the view */
void calculateCoverageViewBorders(GtkWidget *coverageView)
{
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
  BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
  
  /* Calculate the height based on the number of cells */
  const int height = ceil(properties->numVCells * (double)bigPictureGetCellHeight(bigPicture));
  const int gridHeight = (int)properties->numVCells * bigPictureGetCellHeight(bigPicture);
  
  properties->displayRect.x = roundNearest(bpProperties->charWidth * (gdouble)bpProperties->leftBorderChars);
  properties->displayRect.y = height - gridHeight;
  
  properties->viewRect.x = properties->displayRect.x;
  properties->viewRect.y = properties->displayRect.y + bpProperties->highlightBoxYPad + DEFAULT_COVERAGE_VIEW_Y_PADDING;

  properties->displayRect.width = coverageView->allocation.width - properties->viewRect.x;
  properties->displayRect.height = height + 2 * (bpProperties->highlightBoxYPad + DEFAULT_COVERAGE_VIEW_Y_PADDING);

  properties->viewRect.width = properties->displayRect.width;
  properties->viewRect.height = gridHeight;
  
  /* Get the boundaries of the highlight box */
  calculateHighlightBoxBorders(&properties->displayRect, &properties->highlightRect, bigPicture, 0);
  
  /* Set the size request to our desired height. We want a fixed heigh but don't set the
   * width, because we want the user to be able to resize that. */
  gtk_widget_set_size_request(coverageView, 0, properties->displayRect.height);
}


/* Prepare the coverage view for printing (draws the transient hightlight box
 * onto the cached drawable). */
void coverageViewPrepareForPrinting(GtkWidget *coverageView)
{
  GdkDrawable *drawable = widgetGetDrawable(coverageView);
  
  if (drawable)
    {
      CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
      GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
      BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
      BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
      
      GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
      drawHighlightBox(drawable, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);
    }
}


/***********************************************************
 *                         Events                          *
 ***********************************************************/

/* Expose handler. */
static gboolean onExposeCoverageView(GtkWidget *coverageView, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(coverageView)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(coverageView);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
	  bitmap = createBlankPixmap(coverageView);
          drawCoverageView(coverageView, bitmap);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
          
          /* Draw the highlight box on top of it */
          CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
          GtkWidget *bigPicture = blxWindowGetBigPicture(properties->blxWindow);
          BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
          BigPictureProperties *bpProperties = bigPictureGetProperties(bigPicture);
          
          GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
          drawHighlightBox(window, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);
          
          /* Draw the preview box too, if set */
          drawPreviewBox(bigPicture, window, &properties->viewRect, &properties->highlightRect);
        }
      else
	{
	  g_warning("Failed to draw coverageView [%p] - could not create bitmap.\n", coverageView);
	}
    }
  
  return TRUE;
}


static void onSizeAllocateCoverageView(GtkWidget *coverageView, GtkAllocation *allocation, gpointer data)\
{
  calculateCoverageViewBorders(coverageView);
}


static gboolean onButtonPressCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
  BigPictureProperties *bpProperties = bigPictureGetProperties(coverageViewGetBigPicture(coverageView));
  
  if (event->button == 2 ||
      (event->button == 1 && !handled && 
       (event->type == GDK_2BUTTON_PRESS || 
        clickedInRect(event, &properties->highlightRect, bpProperties->highlightBoxMinWidth))))
    {
      /* Draw the preview box (draw it on the other big picture components as well) */
      int x = event->x;
      
      if (event->button == 1 && event->type == GDK_BUTTON_PRESS)
        x = properties->highlightRect.x + properties->highlightRect.width / 2;
      
      showPreviewBox(coverageViewGetBigPicture(coverageView), event->x, TRUE, x - event->x);
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseCoverageView(GtkWidget *coverageView, GdkEventButton *event, gpointer data)
{
  if (event->button == 1 || event->button == 2) /* left or middle button */
    {
      CoverageViewProperties *properties = coverageViewGetProperties(coverageView);
      acceptAndClearPreviewBox(coverageViewGetBigPicture(coverageView), event->x, &properties->viewRect, &properties->highlightRect);
    }
  
  return TRUE;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the coverageView.
 * This scrolls the position of the highlight box, i.e. it scrolls the display
 * range in the detail view. */
static gboolean onScrollCoverageView(GtkWidget *coverageView, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
  {
    case GDK_SCROLL_LEFT:
    {
      scrollBigPictureLeftStep(coverageViewGetBigPicture(coverageView));
      handled = TRUE;
      break;
    }
      
    case GDK_SCROLL_RIGHT:
    {
      scrollBigPictureRightStep(coverageViewGetBigPicture(coverageView));
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


static gboolean onMouseMoveCoverageView(GtkWidget *coverageView, GdkEventMotion *event, gpointer data)
{
  if ((event->state & GDK_BUTTON1_MASK) || /* left or middle button */
      (event->state & GDK_BUTTON2_MASK))
    {
      /* Draw a preview box at the mouse pointer location */
      showPreviewBox(coverageViewGetBigPicture(coverageView), event->x, FALSE, 0);
    }
  
  return TRUE;
}


/***********************************************************
 *                     Initialisation                      *
 ***********************************************************/

GtkWidget* createCoverageView(GtkWidget *blxWindow, BlxViewContext *bc)
{
  GtkWidget *coverageView = gtk_layout_new(NULL, NULL);

  /* Style properties */
  gtk_widget_set_redraw_on_allocate(coverageView, FALSE);
  gtk_widget_set_name(coverageView, COVERAGE_VIEW_NAME);

  /* Connect signals */
  gtk_widget_add_events(coverageView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(coverageView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(coverageView, GDK_POINTER_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(coverageView), "expose-event",          G_CALLBACK(onExposeCoverageView),                 NULL);  
  g_signal_connect(G_OBJECT(coverageView), "size-allocate",	    G_CALLBACK(onSizeAllocateCoverageView),           NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-press-event",    G_CALLBACK(onButtonPressCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "button-release-event",  G_CALLBACK(onButtonReleaseCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "motion-notify-event",   G_CALLBACK(onMouseMoveCoverageView),	      NULL);
  g_signal_connect(G_OBJECT(coverageView), "scroll-event",	    G_CALLBACK(onScrollCoverageView),                 NULL);

  /* Set required data in the coverageView. */
  coverageViewCreateProperties(coverageView, blxWindow, bc);
  
  return coverageView;
}

