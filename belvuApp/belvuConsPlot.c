/*  File: belvuConsPlot.h
 *  Author: Gemma Barson, 2011-07-01
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
 * Description: See belvuConsPlot.h
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvuConsPlot.h"
#include "belvuApp/belvuWindow.h"
#include "seqtoolsUtils/utilities.h"
#include <gtk/gtk.h>


#define DEFAULT_CONS_PLOT_SCALE_HEIGHT      200   /* default height of the conservation plot scale */
#define CONS_PLOT_XPAD                      10    /* Default x padding for the conservation plot */
#define CONS_PLOT_YPAD                      10    /* Default y padding for the conservation plot */
#define CONS_PLOT_LABEL_PAD                 4     /* Padding around labels in the conservation plot */
#define TICKMARK_INTERVAL                   10    /* how many values between tickmark labels */
#define MAJOR_TICKMARK_HEIGHT               6     /* height of major tick marks in the sequence area header */
#define MINOR_TICKMARK_HEIGHT               3     /* height of minor tick marks in the sequence area header */
#define AVG_CONSERVATION_LABEL              "Average conservation"  /* label to show on the plot's "average conservation" line */
#define CONS_PLOT_WINDOW_YPAD               30    /* extra y padding allocated to the plot's window on initialisation */


/* Local function declarations */
static void                         calculateConservation(GtkWidget *consPlot);
static void                         calculateConsPlotBorders(GtkWidget *consPlot);



/***********************************************************
 *                         Properties                      *
 ***********************************************************/

/* Properties specific to the belvu window */
typedef struct _ConsPlotProperties
  {
    BelvuContext *bc;                   /* The belvu context */
    GtkActionGroup *actionGroup;
    
    GtkWidget *drawingArea;             /* The drawing widget */
    GdkRectangle headerRect;            /* Space for drawing the header line */
    GdkRectangle plotRect;              /* Space for drawing the main plot */
    GdkRectangle yScaleRect;            /* Space for drawing the y scale */
    GdkRectangle xScaleRect;            /* Space for drawing the y scale */
    
    int windowSize;                     /* Size of the sliding window used for smothing the profile */
    int lineWidth;                      /* Line width of the plot */
    double xScale;                      /* Used for scaling the x axis */
    double yScale;                      /* Used for scaling the y axis */
    
    double maxcons;                     /* maximum conservation */
    double mincons;                     /* minimum conservation */
    double avgcons;                     /* average conservation */
    double *smooth;                     /* Used for calculating the smoothed profile */
  } ConsPlotProperties;


/* Properties for generic windows */
static ConsPlotProperties* consPlotGetProperties(GtkWidget *consPlot)
{
  return consPlot ? (ConsPlotProperties*)(g_object_get_data(G_OBJECT(consPlot), "BelvuConsPlotProperties")) : NULL;
}


/* Does the job of destroying the window */
static void onDestroyConsPlot(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  if (properties)
    {
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(consPlot), "BelvuConsPlotProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void consPlotCreateProperties(GtkWidget *consPlot, 
                                     BelvuContext *bc, 
                                     GtkActionGroup *actionGroup,
                                     GtkWidget *drawingArea)
{
  if (consPlot)
    {
      ConsPlotProperties *properties = (ConsPlotProperties*)g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->actionGroup = actionGroup;

      properties->drawingArea = drawingArea;
      
      properties->windowSize = 1;
      properties->lineWidth = 1;
      properties->xScale = 10.0;
      properties->yScale = 20.0;
      
      properties->smooth = NULL;
      
      g_object_set_data(G_OBJECT(consPlot), "BelvuConsPlotProperties", properties);
      g_signal_connect(G_OBJECT(consPlot), "destroy", G_CALLBACK (onDestroyConsPlot), NULL);
    }
}


BelvuContext* consPlotGetContext(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  return (properties ? properties->bc : NULL);
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Clear cached drawables and redraw everything */
static void belvuConsPlotRedrawAll(GtkWidget *consPlot)
{
  if (!consPlot)
    return;
  
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  if (properties && properties->drawingArea)
    {
      widgetClearCachedDrawable(properties->drawingArea, NULL);
      gtk_widget_queue_draw(properties->drawingArea);
    }
}


/* Recalculate the conservation array and redraw everything */
void belvuConsPlotRecalcAll(GtkWidget *consPlot)
{
  if (consPlot)
    {
      calculateConservation(consPlot);
      calculateConsPlotBorders(consPlot);
    }
}


static void drawConsPlot(GtkWidget *widget, GdkDrawable *drawable, ConsPlotProperties *properties)
{
  if (!properties->drawingArea)
    return;
  
  BelvuContext *bc = properties->bc;
  GdkGC *gc = gdk_gc_new(drawable);
  gdk_gc_set_line_attributes(gc, properties->lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  
  GtkAdjustment *hAdjustment = gtk_layout_get_hadjustment(GTK_LAYOUT(properties->drawingArea));
  const double xMin = hAdjustment->value;
  const double xMax = min(hAdjustment->value + hAdjustment->page_size, hAdjustment->upper);
  int iMin = (xMin - (double)properties->plotRect.x) / properties->xScale;
  int iMax = (xMax - (double)properties->plotRect.x) / properties->xScale;
  
  if (iMin < 0)
    iMin = 0;
  
  if (iMax > bc->maxLen - 1)
    iMax = bc->maxLen - 1;
  
  /* Draw the header line */
  char *tmpStr = g_strdup_printf("Window = %d", properties->windowSize);
  drawText(widget, drawable, gc, properties->headerRect.x, properties->headerRect.y, tmpStr, NULL, NULL);
  g_free(tmpStr);
  
  /* Draw x scale */
  GdkColor *scaleColor = getGdkColor(BELCOLOR_CONS_PLOT_SCALE, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, scaleColor);

  const int majorXTickInterval = 10;
  const int minorXTickInterval = 5;
  
  /* Get the start and end x coords of the plot. Note that we offset by the
   * minimum x coord so that the leftmost edge is 0 rather than xMin */
  double xStart = max(xMin, properties->xScaleRect.x) - xMin;
  const double xEnd = min(xMax, properties->xScaleRect.x + properties->xScaleRect.width) - xMin;
  
  /* Draw the base line for the x scale */
  double y = properties->xScaleRect.y;
  gdk_draw_line(drawable, gc, xStart, y, xEnd, y);

  /* Loop through each visible column that will have a major tickmark. For
   * each major tickmark we also draw the minor tickmark previous to it, so 
   * loop until (i - minorXTickInterval) is outside the visible range. */
  int i = roundToValue(iMin, majorXTickInterval);
  double x = 0;
  
  for ( ; i - minorXTickInterval <= iMax; i += majorXTickInterval) 
    {
      /* Get the x position of this column */
      x = properties->xScaleRect.x + (i * properties->xScale) - xMin;
      
      /* Draw the major tickmark (if in range; the very last one might not be) */
      if (i <= iMax)
        gdk_draw_line(drawable, gc, x, y, x, y + MAJOR_TICKMARK_HEIGHT);

      /* Decrement the column number to give the previous minor tickmark
       * position and draw the minor tickmark there (if it is still in the
       * visible area). */
      if (i + minorXTickInterval <= iMax)
        {
          double x2 = properties->xScaleRect.x + ((i + minorXTickInterval) * properties->xScale) - xMin;
          gdk_draw_line(drawable, gc, x2, y, x2, y + MINOR_TICKMARK_HEIGHT);
        }
      
      tmpStr = g_strdup_printf("%d", i);
      drawText(widget, drawable, gc, x - minorXTickInterval, y + MAJOR_TICKMARK_HEIGHT + CONS_PLOT_LABEL_PAD, tmpStr, NULL, NULL);
      g_free(tmpStr);
    }
  
  /* Draw the y scale. Note that the scale is inversed; that is, the lowest number
   * is at the bottom, i.e. with the highest y position. */
  const double majorYTickInterval = 1.0;
  const double minorYTickInterval = 0.5;

  x = max(xMin, properties->yScaleRect.x + properties->yScaleRect.width) - xMin;;
  const double yStart = properties->yScaleRect.y + properties->yScaleRect.height;
  const double yEnd = properties->yScaleRect.y;
  
  gdk_draw_line(drawable, gc, x, yEnd, x, yStart);

  double f = properties->mincons;
  y = 0;

  for ( ; f < properties->maxcons; f += majorYTickInterval) 
    {
      y = yStart - (f * properties->yScale);
      
      gdk_draw_line(drawable, gc, x, y, x - MAJOR_TICKMARK_HEIGHT, y);
      
      if (f + minorYTickInterval < properties->maxcons) 
        {
          double y2 = yStart - ((f + minorYTickInterval) * properties->yScale);
          gdk_draw_line(drawable, gc, x, y2, x - MINOR_TICKMARK_HEIGHT, y2);
        }
      
      tmpStr = g_strdup_printf("%.0f", f);
      PangoLayout *layout = gtk_widget_create_pango_layout(widget, tmpStr);
      g_free(tmpStr);

      int textWidth = 0, textHeight = 0;
      pango_layout_get_size(layout, &textWidth, &textHeight);
      textWidth /= PANGO_SCALE;
      textHeight /= PANGO_SCALE;
      
      gdk_draw_layout(drawable, gc, 
                      x - textWidth - MAJOR_TICKMARK_HEIGHT - CONS_PLOT_LABEL_PAD, 
                      y - textHeight / 2, layout);
      
      g_object_unref(layout);
    }
  
  /* Draw average line */
  GdkColor *avgLineColor = getGdkColor(BELCOLOR_CONS_PLOT_AVG, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, avgLineColor);
  const double yAvg = yStart - (properties->avgcons * properties->yScale);
  
  gdk_draw_line(drawable, gc, xStart, yAvg, xEnd, yAvg);
  drawText(widget, drawable, gc, xEnd + CONS_PLOT_XPAD, yAvg, AVG_CONSERVATION_LABEL, NULL, NULL);
  
  /* Plot the conservation value for each column */
  GdkColor *plotColor = getGdkColor(BELCOLOR_CONS_PLOT, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, plotColor);
  
  for (i = iMin + 1; i <= iMax; ++i) 
    {
      x = properties->xScaleRect.x + (i * properties->xScale) - xMin;
      
      const double y1 = yStart - (properties->smooth[i - 1] * properties->yScale);
      const double y2 = yStart - (properties->smooth[i] * properties->yScale);
      
      gdk_draw_line(drawable, gc, x, y1, x + properties->xScale, y2);
    }

  /* Clean up */
  g_object_unref(gc);
}


/***********************************************************
 *                   Settings dialog                       *
 ***********************************************************/

static void showSettingsDialog(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  BelvuContext *bc = properties->bc;
  
  char *title = g_strdup_printf("%sPlot Settings", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(consPlot), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  const int numRows = 3;
  const int numCols = 4;
  const int xpad = TABLE_XPAD;
  const int ypad = TABLE_YPAD;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), GTK_WIDGET(table), FALSE, FALSE, 12);

  GtkWidget *winEntry = createTextEntryFromInt(dialog, table, 0, 1, xpad, ypad, "_Window size: ", properties->windowSize, NULL);
  GtkWidget *xEntry = createTextEntryFromInt(dialog, table, 1, 1, xpad, ypad, "_X scale: ", properties->xScale, NULL);
  GtkWidget *yEntry = createTextEntryFromInt(dialog, table, 1, 3, 0, ypad, "_Y scale: ", properties->yScale, NULL);
  GtkWidget *lineEntry = createTextEntryFromInt(dialog, table, 2, 1, xpad, ypad, "_Line width: ", properties->lineWidth, NULL);

  gtk_widget_show_all(dialog);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const int newWin = convertStringToInt(gtk_entry_get_text(GTK_ENTRY(winEntry)));
      const int newX = convertStringToInt(gtk_entry_get_text(GTK_ENTRY(xEntry)));
      const int newY = convertStringToInt(gtk_entry_get_text(GTK_ENTRY(yEntry)));
      const int newLine = convertStringToInt(gtk_entry_get_text(GTK_ENTRY(lineEntry)));

      gboolean changed = FALSE;
      
      /* Recalculate the conservation array if the window size has changed */
      if (newWin != properties->windowSize)
        {
          changed = TRUE;
          properties->windowSize = newWin;
          calculateConservation(consPlot);
        }
      
      /* Recalculate the window size if either scale has changed */
      if (newX != properties->xScale || newY != properties->yScale)
        {
          changed = TRUE;
          properties->xScale = newX;
          properties->yScale = newY;
          calculateConsPlotBorders(consPlot);
        }

      /* Redraw if anything has changed */
      if (changed || newLine != properties->lineWidth)
        {
          changed = TRUE;
          properties->lineWidth = newLine;
          belvuConsPlotRedrawAll(consPlot);
        }
    }
  
  gtk_widget_destroy(dialog);  
}


/***********************************************************
 *                     Calculations                        *
 ***********************************************************/

/* Calculate the max, min and average conservation. */
static void calculateConservation(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  BelvuContext *bc = properties->bc;
  
  /* Smooth the conservation profile by applying a window */
  if (properties->smooth)
    g_free(properties->smooth);
  
  properties->smooth = (double*)g_malloc(bc->maxLen * sizeof(double));
  
  double sum = 0.0;
  
  int i = 0;
  for (i = 0; i < properties->windowSize; ++i) 
    sum += bc->conservation[i];
  
  properties->smooth[properties->windowSize / 2] = sum / properties->windowSize;
  
  for ( ; i < bc->maxLen; ++i) 
    {
      sum -= bc->conservation[i-properties->windowSize];
      sum += bc->conservation[i];
      properties->smooth[i - properties->windowSize / 2] = sum/properties->windowSize;
    }
  
  /* Find max and min and avg conservation */
  properties->maxcons = -1;
  properties->mincons = 10000;
  properties->avgcons = 0;
  
  for (i = 0; i < bc->maxLen; ++i) 
    {
      if (properties->smooth[i] > properties->maxcons) 
        properties->maxcons = properties->smooth[i];
      
      if (properties->smooth[i] < properties->mincons) 
        properties->mincons = properties->smooth[i];
      
      properties->avgcons += bc->conservation[i];
    }
  
  properties->avgcons /= bc->maxLen * 1.0;
}


/* Calculate the size of the drawing areas */
static void calculateConsPlotBorders(GtkWidget *consPlot)
{  
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  /* Calculate plot size */
  properties->plotRect.width = properties->xScale * properties->bc->maxLen;
  properties->plotRect.height = (properties->maxcons - properties->mincons) * properties->yScale;

  /* y scale size (based on width required to display max conservation) */
  char *tmpStr = g_strdup_printf("%.0f", properties->maxcons);
  int textWidth, textHeight;
  getTextSize(properties->drawingArea, tmpStr, &textWidth, &textHeight);
  g_free(tmpStr);
  
  properties->yScaleRect.width = textWidth + MAJOR_TICKMARK_HEIGHT + CONS_PLOT_XPAD;
  properties->yScaleRect.height = properties->plotRect.height;
  
  /* header line size (one character high) */
  properties->headerRect.width = properties->plotRect.width;
  properties->headerRect.height = textHeight;

  /* x scale size */
  properties->xScaleRect.width = properties->plotRect.width;
  properties->xScaleRect.height = textHeight + MAJOR_TICKMARK_HEIGHT + CONS_PLOT_YPAD;

  /* Calculate x and y coords */
  properties->headerRect.x = CONS_PLOT_XPAD;
  properties->headerRect.y = CONS_PLOT_YPAD;
  
  properties->yScaleRect.x = CONS_PLOT_XPAD;
  properties->yScaleRect.y = properties->headerRect.y + properties->headerRect.height + CONS_PLOT_YPAD;
  
  properties->plotRect.x = properties->yScaleRect.x + properties->yScaleRect.width + CONS_PLOT_XPAD;
  properties->plotRect.y = properties->yScaleRect.y;
  
  properties->xScaleRect.x = properties->plotRect.x;
  properties->xScaleRect.y = properties->plotRect.y + properties->plotRect.height + CONS_PLOT_YPAD;

  /* Set the size of the layout. This must include the rightmost extent of the plot
   * rectangle, and also some extra space on the right for the 'average conservation'
   * label. */
  const int width = properties->plotRect.x + properties->plotRect.width + 
                    getTextWidth(properties->drawingArea, AVG_CONSERVATION_LABEL) +
                    (CONS_PLOT_XPAD * 2);
  
  gtk_layout_set_size(GTK_LAYOUT(properties->drawingArea), 
                      width, 
                      properties->xScaleRect.y + properties->xScaleRect.height);
  
  belvuConsPlotRedrawAll(consPlot);
}

/***********************************************************
 *                        Events                           *
 ***********************************************************/

static gboolean onExposeConsPlot(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *consPlot = GTK_WIDGET(data);
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, 
                                          widget->allocation.width, 
                                          properties->xScaleRect.y + properties->xScaleRect.height);

          drawConsPlot(widget, bitmap, properties);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          GtkAdjustment *hAdjustment = gtk_layout_get_hadjustment(GTK_LAYOUT(widget));
          gdk_draw_drawable(window, gc, bitmap, 0, 0, hAdjustment->value, 0, -1, -1);
          g_object_unref(gc);
        }
      else
        {
          g_warning("Failed to draw conservation plot [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}


static void onSizeAllocateConsPlot(GtkWidget *widget, GtkAllocation *allocation, gpointer data)
{
  GtkWidget *consPlot = GTK_WIDGET(data);
  belvuConsPlotRecalcAll(consPlot);
}


void onPlotOptsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *consPlot = GTK_WIDGET(data);
  showSettingsDialog(consPlot);
}


static gboolean onScrollPosChangedConsPlot(GtkObject *object, gpointer data)
{
  GtkWidget *consPlot = GTK_WIDGET(data);
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  if (properties->drawingArea)
    {
      belvuConsPlotRedrawAll(consPlot);
    }
  
  return FALSE;
}

/***********************************************************
 *                     Initialisation                      *
 ***********************************************************/

static void setConsPlotStyleProperties(GtkWidget *window, BelvuContext *bc, const int height)
{
  gtk_widget_set_name(window, BELVU_CONS_PLOT_WINDOW_NAME);
  
  /* Just hide the widget when it is closed rather than destroy it */
  g_signal_connect(window, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
  
  /* Set default size based on scale height and alignment window width */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_BELVU_WINDOW_WIDTH_FRACTION;
  
  gtk_window_set_default_size(GTK_WINDOW(window), width, height + CONS_PLOT_WINDOW_YPAD);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


void createConsPlot(BelvuContext *bc)
{
  /* Create the window */
  bc->consPlot = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  
  char *title = g_strdup_printf("%sConservation Profile", belvuGetTitlePrefix(bc));
  gtk_window_set_title(GTK_WINDOW(bc->consPlot), title);
  g_free(title);
  
  /* We must add all toplevel windows to the list of spawned windows */
  bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, bc->consPlot);
  
  /* Create the context menu and set a callback to show it */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(bc->consPlot, bc, &actionGroup);
  GtkWidget *contextmenu = createBelvuMenu(bc->consPlot, "/PlotContextMenu", uiManager);
  
  gtk_widget_add_events(bc->consPlot, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(bc->consPlot), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(bc->consPlot), vbox);
  
  /* We'll place the drawing area in a scrolled window */
  GtkWidget *drawing = gtk_layout_new(NULL, NULL);
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrollWin), drawing);
  
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_ALWAYS, GTK_POLICY_ALWAYS);
  gtk_box_pack_start(GTK_BOX(vbox), scrollWin, TRUE, TRUE, 0);
  
  /* Set default background color */
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);
  gtk_widget_modify_bg(drawing, GTK_STATE_NORMAL, bgColor);
  
  /* Set properties */
  consPlotCreateProperties(bc->consPlot, bc, actionGroup, drawing);
  
  /* Calculate the plot size so we can set the initial window height correctly.
   * Must be done after setting properties, but before showing the widgets. */
  calculateConservation(bc->consPlot);
  calculateConsPlotBorders(bc->consPlot);
  ConsPlotProperties *properties = consPlotGetProperties(bc->consPlot);
  setConsPlotStyleProperties(bc->consPlot, bc, properties->xScaleRect.y + properties->xScaleRect.height + scrollBarWidth());

  /* Connect signals */
  GtkAdjustment *hAdjustment = gtk_layout_get_hadjustment(GTK_LAYOUT(drawing));
  g_signal_connect(G_OBJECT(drawing), "size-allocate", G_CALLBACK(onSizeAllocateConsPlot), bc->consPlot);
  g_signal_connect(G_OBJECT(drawing), "expose-event", G_CALLBACK(onExposeConsPlot), bc->consPlot);
  g_signal_connect(G_OBJECT(hAdjustment), "value-changed",  G_CALLBACK(onScrollPosChangedConsPlot), bc->consPlot);

  gtk_widget_show_all(bc->consPlot);
  gtk_window_present(GTK_WINDOW(bc->consPlot));
}
