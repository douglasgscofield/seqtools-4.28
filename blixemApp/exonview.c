/*  File: exonview.c
 *  Author: Gemma Barson, 2009-12-24
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
 * Description: See exonview.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/exonview.h>
#include <blixemApp/detailviewtree.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/bigpicturegrid.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blxwindow.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <stdlib.h>

#define DEFAULT_EXON_HEIGHT			 10
#define DEFAULT_EXON_HEIGHT_BUMPED		 7
#define	DEFAULT_EXON_YPAD			 4
#define	DEFAULT_EXON_YPAD_BUMPED		 4

typedef struct _ExonViewProperties
  {
    GtkWidget *bigPicture;	      /* The big picture that this view belongs to */
    BlxStrand currentStrand;	      /* Which strand of the ref seq this view displays exons for */
    
    gboolean expanded;		      /* Whether the exon view is expanded or compressed */
    
    int yPad;			      /* y padding */
    
    GdkRectangle exonViewRect;	      /* The drawing area for the exon view */
    GdkRectangle highlightRect;       /* The area that the highlight box will cover (indicating the current detail-view display range) */
    
    int exonHeight;                   /* the height of an individual exon */ 
  } ExonViewProperties;


/* Utility struct to pass around data required for drawing exons */
typedef struct _DrawData
  {
    GdkDrawable *drawable;
    GdkGC *gc;
    GdkRectangle *exonViewRect;
    GtkWidget *blxWindow;
    BlxViewContext *bc;
    const BlxStrand strand;
    const IntRange* const displayRange;
    const IntRange* const refSeqRange;
    const gboolean displayRev;
    const int numFrames;
    const BlxSeqType seqType;
    gboolean expanded;		      /* whether exon view is expaned or compressed */
    gboolean normalOnly;	      /* Only draw "normal" MSPs, i.e. thse that are not selected and are not in a group */
    int yPad;			      /* y padding */
    int y;			      /* y position to draw this exon at (constant if view compressed; increases if view bumped) */
    int height;			      /* height of exon box */
  } DrawData;



/* Local function declarations */
static GtkWidget*		exonViewGetBigPicture(GtkWidget *exonView);
static GtkWidget*		exonViewGetBlxWindow(GtkWidget *exonView);
static ExonViewProperties*	exonViewGetProperties(GtkWidget *exonView);


/***********************************************************
 *                       Utility functions                 *
 ***********************************************************/

/* Calls the given function (passed as the data pointer) on the given widget 
 * if it is an exon view in the big picture view, or, if it is a container, 
 * calls the function on all children/grandchildren/etc that are exon views */
void callFuncOnAllBigPictureExonViews(GtkWidget *widget, gpointer data)
{
  const gchar *name = gtk_widget_get_name(widget);
  if (strcmp(name, BIG_PICTURE_EXON_VIEW_NAME) == 0)
    {
      GtkCallback func = (GtkCallback)data;
      func(widget, NULL);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllBigPictureExonViews, data);
    }
}


static gboolean calculateExonIntronDimensions(const MSP* const msp, 
                                              BlxViewContext *bc,
                                              const IntRange* const displayRange,
                                              GdkRectangle *exonViewRect,
                                              int *x, 
                                              int *width)
{
  gboolean drawn = FALSE;
  
  /* Find the coordinates of the start and end base in this msp, converting to display coords. Note
   * that display coords always increase from left-to-right, even if the actual coords are inverted. */
  const IntRange* const mspDisplayRange = mspGetDisplayRange(msp);

  if (rangesOverlap(mspDisplayRange, displayRange))
    {
      drawn = TRUE;

      /* Get the display range in dna coords */
      IntRange dnaDispRange;
      convertDisplayRangeToDnaRange(displayRange, bc->seqType,
                                    bc->numFrames, bc->displayRev,
                                    &bc->refSeqRange, &dnaDispRange);
      
      /* The grid pos for coords gives the left edge of the coord, so draw to max + 1 to be inclusive */
      const int qStart = msp->qRange.min;
      const int qEnd = msp->qRange.max + 1;
      
      const gint x1 = convertBaseIdxToRectPos(qStart, exonViewRect, &dnaDispRange,
                                              TRUE, bc->displayRev, FALSE);
      const gint x2 = convertBaseIdxToRectPos(qEnd, exonViewRect, &dnaDispRange, 
                                              TRUE, bc->displayRev, FALSE);
      
      *x = min(x1, x2);
      *width = abs(x1 - x2);
    }
 
  return drawn;
}


/* Returns true if sequences in the given group should be shown. */
static gboolean isGroupVisible(const SequenceGroup* const group)
{
  gboolean result = TRUE;
  
  result = (!group || !group->hidden);
  
  return result;
}


/* Returns true if the given MSP should be shown in this exon view */
static gboolean showMspInExonView(const MSP *msp, const BlxStrand strand, BlxViewContext *bc)
{
  /* Check it's an exon or intron */
  gboolean showMsp = mspIsExon(msp) || mspIsIntron(msp);
  
  /* Check it's in a visible layer */
  showMsp &= mspLayerIsVisible(msp);
  
  /* Check it's the correct strand */
  showMsp &= (mspGetRefStrand(msp) == strand);

  /* Check if its in a group that's hidden */
  if (showMsp)
    {
      SequenceGroup *group = blxContextGetSequenceGroup(bc, msp->sSequence) ;
      showMsp &= isGroupVisible(group) ;

      /* If hiding ungrouped features and this is a feature (i.e. not a match) without a group, hide it */
      if (!group && bc->flags[BLXFLAG_HIDE_UNGROUPED_FEATURES] && msp->sSequence->type != BLXSEQUENCE_MATCH)
        showMsp = FALSE ;
    }

  return showMsp;
}


/* Mark the given exon's transcript as selected if this exon contains the given coords.
 * Returns true if it was selected. */
static gboolean selectExonIfContainsCoords(GtkWidget *exonView, 
                                           ExonViewProperties *properties,
                                           BlxViewContext *bc,
                                           const IntRange* const displayRange,
                                           const MSP *msp,
                                           const int x,
                                           const int y,
                                           gboolean deselectOthers,
                                           const int mspY,
                                           gboolean *drawn)
{
  gboolean wasSelected = FALSE;
  
  if (showMspInExonView(msp, properties->currentStrand, bc) && 
      rangesOverlap(displayRange, mspGetDisplayRange(msp)))
    {
      *drawn = TRUE;
      
      int mspX = UNSET_INT, mspWidth = UNSET_INT;
      calculateExonIntronDimensions(msp, bc, displayRange, &properties->exonViewRect, &mspX, &mspWidth);

      if (x >= mspX && x <= mspX + mspWidth && y >= mspY && y <= mspY + properties->exonHeight)
	{
	  /* It's a hit. Select this sequence. */
	  GtkWidget *blxWindow = exonViewGetBlxWindow(exonView);
	  
	  if (deselectOthers)
	    {
	      blxWindowDeselectAllSeqs(blxWindow);
	    }
	  
	  blxWindowSelectSeq(blxWindow, msp->sSequence);
	  
	  /* Update the selected strand */
	  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
          ExonViewProperties *properties = exonViewGetProperties(exonView);
	  detailViewSetSelectedStrand(detailView, properties->currentStrand);

	  /* Scroll the detail view trees to bring the new selection into view */
	  callFuncOnAllDetailViewTrees(detailView, treeScrollSelectionIntoView, NULL);

	  wasSelected = TRUE;
	}
    }
  
  return wasSelected;
}


/* Loop through all the msp lines for this grid and mark them as selected
 * if they contain the coords of the mouse press */
static gboolean selectClickedExon(GtkWidget *exonView, 
                                  GdkEventButton *event, 
                                  const gboolean ctrlModifier, 
                                  const gboolean shiftModifier)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  const gboolean deselectOthers = !ctrlModifier && !shiftModifier; /* whether to deselect all others first */

  /* Loop through all the MSPs until we find one under the click coords */
  GtkWidget *blxWindow = exonViewGetBlxWindow(exonView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GList *seqItem = blxWindowGetAllMatchSeqs(blxWindow);

  gboolean found = FALSE;
  int y = properties->yPad;
  
  for ( ; seqItem && !found; seqItem = seqItem->next)
    {
      const BlxSequence* const seq = (const BlxSequence*)(seqItem->data);
      GList *mspItem = seq->mspList;
      gboolean drawn = FALSE;
            
      for ( ; mspItem && !found; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);
          found = selectExonIfContainsCoords(exonView, properties, bc, 
                                             bigPictureGetDisplayRange(properties->bigPicture),
                                             msp, event->x, event->y, deselectOthers, y, &drawn);
        }

      if (drawn && properties->expanded)
        y += properties->exonHeight + properties->yPad;
    }

  return found;
}


void calculateExonViewHeight(GtkWidget *exonView)
{
  DEBUG_ENTER("calculateExonViewHeight");

  ExonViewProperties *properties = exonViewGetProperties(exonView);

  BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
  const IntRange* const displayRange = &bpProperties->displayRange;
  
  BlxViewContext *bc = blxWindowGetContext(bpProperties->blxWindow);

  /* Calculate the height based on how many exon lines will actually be drawn */
  int numExons = 0;
  int maxExons = properties->expanded ? UNSET_INT : 1; /* unset means no limit */
  
  /* Loop through all sequences */
  GList *seqItem = bc->matchSeqs;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      /* Loop through all msps */
      const BlxSequence *seq = (BlxSequence*)(seqItem->data);
      GList *mspItem = seq->mspList;
      
      for ( ; mspItem; mspItem = mspItem->next)
	{
	  const MSP *msp = (const MSP*)(mspItem->data);
          
          if (showMspInExonView(msp, properties->currentStrand, bc))
	    {
	      const IntRange* const mspDisplayRange = mspGetDisplayRange(msp);
              
              if (rangesOverlap(mspDisplayRange, displayRange))
		{
		  ++numExons;
		  break; /* break inner loop and move to next sequence */
		}
	    }
	}
      
      /* Break after we've found the maximum number of lines, if a max is specified */
      if (maxExons != UNSET_INT && numExons >= maxExons)
	{
	  break;
	}
    }
  
  const int newHeight = (numExons * (properties->exonHeight + properties->yPad)) + (2 * properties->yPad);
  
  if (newHeight != properties->exonViewRect.height)
    {
      DEBUG_OUT("Setting new height = %d\n", newHeight);
      properties->exonViewRect.height = newHeight;
      gtk_widget_set_size_request(exonView, -1, properties->exonViewRect.height);
    }
  
  DEBUG_EXIT("calculateExonViewHeight returning");
}


void calculateExonViewHighlightBoxBorders(GtkWidget *exonView)
{
  DEBUG_ENTER("calculateExonViewHighlightBoxBorders");

  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
  
  /* Get the big picture display range in dna coords */
  IntRange bpRange;
  convertDisplayRangeToDnaRange(bigPictureGetDisplayRange(properties->bigPicture), bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &bpRange);

  /* Get the detail view display range in dna coords */
  IntRange dvRange;
  GtkWidget *detailView = bigPictureGetDetailView(properties->bigPicture);
  convertDisplayRangeToDnaRange(detailViewGetDisplayRange(detailView), bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &dvRange);
  
  /* Calculate how many pixels from the left edge of the widget to the first base in the range. */
  const int x1 = convertBaseIdxToRectPos(dvRange.min, &properties->exonViewRect, &bpRange, TRUE, bc->displayRev, TRUE);
  const int x2 = convertBaseIdxToRectPos(dvRange.max + 1, &properties->exonViewRect, &bpRange, TRUE, bc->displayRev, TRUE);
  
  properties->highlightRect.x = min(x1, x2);
  properties->highlightRect.y = 0;
  
  properties->highlightRect.width = abs(x1 - x2);
  properties->highlightRect.height = exonView->allocation.height;
  
  DEBUG_EXIT("calculateExonViewHighlightBoxBorders returning");
}


static void calculateExonViewBorders(GtkWidget *exonView)
{
  DEBUG_ENTER("calculateExonViewBorders");

  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BigPictureProperties *bigPictureProperties = bigPictureGetProperties(properties->bigPicture);
  
  /* Calculate the size of the exon view */
  properties->exonViewRect.x = roundNearest(bigPictureProperties->charWidth * (gdouble)bigPictureProperties->leftBorderChars);
  properties->exonViewRect.width = exonView->allocation.width - properties->exonViewRect.x;
  
  /* Calculate the size of the highlight box */
  calculateExonViewHighlightBoxBorders(exonView);
  
  DEBUG_EXIT("calculateExonViewBorders returning");
}


/***********************************************************
 *                       Drawing functions                 *
 ***********************************************************/

/* Draw an exon */
static void drawExon(const MSP* const msp, 
                     DrawData *data, 
                     const BlxSequence *blxSeq, 
                     const gboolean isSelected, 
                     const gint x, 
                     const gint y,
                     const gint widthIn,
                     const gint height)
{
  /* Clip to the drawing area (or just beyond, so we don't get end lines where the exon doesn't really end) */
  gint xStart = x;
  gint xEnd = x + widthIn;
  gint width = widthIn;
  
  const gint xMin = data->exonViewRect->x;
  const gint xMax = data->exonViewRect->x + data->exonViewRect->width;
  
  if (xStart <= xMax && xEnd >= xMin)
    {
      if (xStart < xMin)
        {
          xStart = xMin - 1;
          width = xEnd - xStart;
        }

      if (xEnd > xMax)
        {
          xEnd = xMax + 1;
          width = xEnd - xStart;
        }

      /* Draw the fill rectangle */
      const GdkColor *fillColor = mspGetColor(msp, data->bc->defaultColors, BLXCOLOR_BACKGROUND, blxSeq, isSelected, data->bc->usePrintColors, TRUE, BLXCOLOR_EXON_FILL, BLXCOLOR_EXON_LINE, BLXCOLOR_CDS_FILL, BLXCOLOR_CDS_LINE, BLXCOLOR_UTR_FILL, BLXCOLOR_UTR_LINE);
      gdk_gc_set_foreground(data->gc, fillColor);
      gdk_draw_rectangle(data->drawable, data->gc, TRUE, xStart, y, width, height);
      
      /* Draw outline (exon box outline always the same (unselected) color; only intron lines change when selected) */
      const GdkColor *lineColor = mspGetColor(msp, data->bc->defaultColors, BLXCOLOR_BACKGROUND, blxSeq, isSelected, data->bc->usePrintColors, FALSE, BLXCOLOR_EXON_FILL, BLXCOLOR_EXON_LINE, BLXCOLOR_CDS_FILL, BLXCOLOR_CDS_LINE, BLXCOLOR_UTR_FILL, BLXCOLOR_UTR_LINE);
      gdk_gc_set_foreground(data->gc, lineColor);
      gdk_draw_rectangle(data->drawable, data->gc, FALSE, xStart, y, width, height);
    }
}
  

/* Utility to actually draw the line for an intron. Clips it if necessary, maintaining the same
 * angle for the line */
static void drawIntronLine(DrawData *data, const gint x1, const gint y1, const gint x2, const gint y2, GdkRectangle *clipRect)
{
  /* Only draw anything if at least part of the line is within range. We are only ever called with
   * y values that are in range so don't bother checking them. */
  const gint xMax = clipRect->x + clipRect->width;
  
  if (x1 <= xMax && x2 >= clipRect->x)
    {
      int xStart = x1;
      int xEnd = x2;
      int yStart = y1;
      int yEnd = y2;
      
      /* Clip the start/end x values if out of range */
      if (xStart < clipRect->x)
        {
          const int origWidth = abs(xEnd - xStart);
        
          xStart = clipRect->x;

          const int newWidth = abs(xEnd - xStart);
          const int newHeight = roundNearest((double)(yEnd - yStart) * (double)newWidth / (double)origWidth); /* negative if yend < ystart */
          
          yStart = yEnd - newHeight;
        }

      if (xEnd > xMax)
        {
          const int origWidth = abs(xEnd - xStart);

          xEnd = xMax;

          const int newWidth = abs(xEnd - xStart);
          const int newHeight = roundNearest((double)(yEnd - yStart) * (double)newWidth / (double)origWidth);
          
          yEnd = yStart + newHeight;
        }
        
      gdk_draw_line(data->drawable, data->gc, xStart, yStart, xEnd, yEnd);
    }
}


/* Draw an intron */
static void drawIntron(const MSP* const msp, 
                       DrawData *data, 
                       const BlxSequence *blxSeq, 
                       const gboolean isSelected, 
                       const gint x, 
                       const gint y, 
                       const gint width, 
                       const gint height)
{
  const GdkColor *lineColor = mspGetColor(msp, data->bc->defaultColors, BLXCOLOR_BACKGROUND, blxSeq, isSelected, data->bc->usePrintColors, FALSE, BLXCOLOR_EXON_FILL, BLXCOLOR_EXON_LINE, BLXCOLOR_CDS_FILL, BLXCOLOR_CDS_LINE, BLXCOLOR_UTR_FILL, BLXCOLOR_UTR_LINE);
  gdk_gc_set_foreground(data->gc, lineColor);

  int yTop = y;
  int yBottom = y + roundNearest((double)height / 2.0);
  
  /* Draw the first section, from the given x to the mid point, sloping up */
  int xStart = x;
  int xEnd = x + roundNearest((double)width / 2.0);
  drawIntronLine(data, xStart, yBottom, xEnd, yTop, data->exonViewRect);

  /* Draw the second section, from the mid point to the end, sloping down */
  xStart  = xEnd;
  xEnd = x + width;
  drawIntronLine(data, xStart, yTop, xEnd, yBottom, data->exonViewRect);
}


/* Draw the given exon/intron, if it is in range. Returns true if it was drawn */
static gboolean drawExonIntron(const MSP *msp, 
                               DrawData *data,
                               const gboolean isSelected,
                               const BlxSequence *blxSeq)
{
  int x = UNSET_INT, width = UNSET_INT;
  
  gboolean drawn = calculateExonIntronDimensions(msp, data->bc, data->displayRange, 
                                                 data->exonViewRect, &x, &width);

  if (drawn)
    {
      if (mspIsExon(msp))
	{
	  drawExon(msp, data, blxSeq, isSelected, x, data->y, width, data->height);
	}
      else if (mspIsIntron(msp))
	{
	  drawIntron(msp, data, blxSeq, isSelected, x, data->y, width, data->height);
	}
    }
  
  return drawn;
}


/* Draw the msps in the given sequence, if they are exons/introns. Use the color
 * specified in the user data */
static void drawExonIntronItem(gpointer listItemData, gpointer data)
{
  const BlxSequence *seq = (const BlxSequence*)listItemData;
  DrawData *drawData = (DrawData*)data;

  const gboolean isSelected = blxWindowIsSeqSelected(drawData->blxWindow, seq);
  SequenceGroup *group = blxContextGetSequenceGroup(drawData->bc, seq);
  gboolean seqDrawn = FALSE;
  
  if (!drawData->normalOnly || (!isSelected && !group))
    {
      /* Loop through all msps in this sequence */
      GList *mspListItem = seq->mspList;
  
      for ( ; mspListItem; mspListItem = mspListItem->next)
	{
	  MSP *msp = (MSP*)(mspListItem->data);
      
          if (showMspInExonView(msp, drawData->strand, drawData->bc))
	    {
	      seqDrawn |= drawExonIntron(msp, drawData, isSelected, seq);
	    }
	}
    }
  
  /* If the view is expanded, increase the y-coord for the next sequence */
  if (seqDrawn && drawData->expanded)
    {
      drawData->y += drawData->height + drawData->yPad;
    }
}


/* Draw the exon view */
static void drawExonView(GtkWidget *exonView, GdkDrawable *drawable)
{
  GtkWidget *blxWindow = exonViewGetBlxWindow(exonView);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  const IntRange* const displayRange = bigPictureGetDisplayRange(properties->bigPicture);

  /* First, highlight any assembly gaps */
  /* Get the display range in dna coords */
  IntRange bpRange;
  convertDisplayRangeToDnaRange(displayRange, bc->seqType, bc->numFrames, bc->displayRev, &bc->refSeqRange, &bpRange);
  
  GdkColor *gapColor = getGdkColor(BLXCOLOR_ASSEMBLY_GAP, bc->defaultColors, FALSE, bc->usePrintColors);
  drawAssemblyGaps(exonView, drawable, gapColor, bc->displayRev, &properties->exonViewRect, &bpRange, bc->featureLists[BLXMSP_GAP]);
  
  /* Set a clip rectangle for drawing the exons and introns (because they are drawn "over the
   * edges" to make sure intron lines have the correct slope etc.) */
  GdkGC *gc = gdk_gc_new(drawable);
  
  gdk_gc_set_clip_origin(gc, 0, 0);
  gdk_gc_set_clip_rectangle(gc, &properties->exonViewRect);
  
  /* Draw the exons and introns. Since we could have a lot of them in the loop, extract all the
   * info we need now and pass it around so we don't have to look for this stuff each time. */
  
  DrawData drawData = {
    drawable,
    gc,
    &properties->exonViewRect,
    blxWindow,
    bc,
    properties->currentStrand,
    displayRange,
    &bc->refSeqRange,
    bc->displayRev,
    bc->numFrames,
    bc->seqType,
    properties->expanded,
    FALSE,
    properties->yPad,
    properties->exonViewRect.y,
    properties->exonHeight
  };
  
  /* If the view is compressed (i.e. exons will overlap each other), then
   * only draw "normal" MSPs the first time round, and draw grouped/selected
   * MSPs afterwards, so that they appear on top. If the view is expanded, 
   * we can draw them all in a single loop, because they will not overlap. */
  drawData.normalOnly = !properties->expanded;
  
  /* Loop through all sequences, drawing all msps that are exons/introns */
  GList *seqList = blxWindowGetAllMatchSeqs(blxWindow);
  g_list_foreach(seqList, drawExonIntronItem, &drawData);

  if (!properties->expanded)
    {
      drawData.normalOnly = FALSE;
  
      /* Draw all selected msps */
      g_list_foreach(bc->selectedSeqs, drawExonIntronItem, &drawData);
      
      /* Increment the y value when finished, because we calculate the view height based on this */
      drawData.y += drawData.height + drawData.yPad;
    }

  /* Set the height based on the height of the exons that were actually drawn */
  const int newHeight = drawData.y - properties->exonViewRect.y + drawData.yPad;
  gtk_layout_set_size(GTK_LAYOUT(exonView), exonView->allocation.width, newHeight);
  
  g_object_unref(gc);
}


/* Prepare the exon view for printing (draws the transient hightlight box
 * onto the cached drawable). */
void exonViewPrepareForPrinting(GtkWidget *exonView)
{
  GdkDrawable *drawable = widgetGetDrawable(exonView);
  
  if (drawable)
    {
      ExonViewProperties *properties = exonViewGetProperties(exonView);
      BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
      BlxViewContext *bc = bigPictureGetContext(properties->bigPicture);
      
      GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
      drawHighlightBox(drawable, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);
    }
}



/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static ExonViewProperties* exonViewGetProperties(GtkWidget *exonView)
{
  return exonView ? (ExonViewProperties*)(g_object_get_data(G_OBJECT(exonView), "ExonViewProperties")) : NULL;
}

static void onDestroyExonView(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  if (properties)
    {
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(exonView), "ExonViewProperties", NULL);
    }
}

static void exonViewCreateProperties(GtkWidget *exonView, 
				     GtkWidget *bigPicture, 
				     const BlxStrand currentStrand)
{
  if (exonView)
    {
      ExonViewProperties *properties = (ExonViewProperties*)g_malloc(sizeof *properties);
      
      properties->bigPicture	      = bigPicture;
      properties->currentStrand	      = currentStrand;
      
      properties->expanded	      = FALSE;
      properties->yPad		      =	DEFAULT_EXON_YPAD;
      
      properties->exonViewRect.x      = 0;
      properties->exonViewRect.y      = DEFAULT_EXON_YPAD;
      properties->exonViewRect.width  = 0;
      properties->exonViewRect.height = DEFAULT_EXON_HEIGHT;

      properties->exonHeight          = DEFAULT_EXON_HEIGHT;
      
      gtk_widget_set_size_request(exonView, 0, DEFAULT_EXON_HEIGHT + (2 * DEFAULT_EXON_YPAD));

      g_object_set_data(G_OBJECT(exonView), "ExonViewProperties", properties);
      g_signal_connect(G_OBJECT(exonView), "destroy", G_CALLBACK(onDestroyExonView), NULL);
    }
}

static GtkWidget* exonViewGetBigPicture(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  return properties->bigPicture;
}

static GtkWidget* exonViewGetBlxWindow(GtkWidget *exonView)
{
  GtkWidget *bigPicture = exonViewGetBigPicture(exonView);
  return bigPictureGetBlxWindow(bigPicture);
}

gboolean exonViewGetExpanded(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  return properties->expanded;
}

/* Set whether the view is expanded or not */
void exonViewSetExpanded(GtkWidget *exonView, const gboolean expanded)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  properties->expanded = expanded;

  if (expanded)
    {
      properties->yPad = DEFAULT_EXON_YPAD_BUMPED;
      properties->exonHeight = DEFAULT_EXON_HEIGHT_BUMPED;
    }
  else 
    {
      properties->yPad = DEFAULT_EXON_YPAD;
      properties->exonHeight = DEFAULT_EXON_HEIGHT;
    }
  
  calculateExonViewHeight(exonView);
  forceResize(properties->bigPicture);
  bigPictureRedrawAll(properties->bigPicture);
}


/* Toggle whether the view is expanded or not */
void exonViewToggleExpanded(GtkWidget *exonView)
{
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  exonViewSetExpanded(exonView, !properties->expanded);
}


/***********************************************************
 *                       Events                            *
 ***********************************************************/

static gboolean onExposeExonView(GtkWidget *exonView, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *drawable = widgetGetDrawable(exonView);
  
  if (!drawable)
    {
      /* Create a pixmap and draw the exon view onto it */
      drawable = createBlankPixmap(exonView);
      drawExonView(exonView, drawable);
    }
  
  if (drawable)
    {  
      /* Push the pixmap onto the screen */
      GdkDrawable *window = GTK_LAYOUT(exonView)->bin_window;
      
      GdkGC *gc = gdk_gc_new(window);
      gdk_draw_drawable(window, gc, drawable, 0, 0, 0, 0, -1, -1);
      g_object_unref(gc);
      
      /* Draw the highlight box on top of it */
      ExonViewProperties *properties = exonViewGetProperties(exonView);
      BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);
      BlxViewContext *bc = blxWindowGetContext(bpProperties->blxWindow);
      
      GdkColor *highlightBoxColor = getGdkColor(BLXCOLOR_HIGHLIGHT_BOX, bc->defaultColors, FALSE, bc->usePrintColors);
      drawHighlightBox(window, &properties->highlightRect, bpProperties->highlightBoxMinWidth, highlightBoxColor);

      /* Draw the preview box too, if it is set */
      drawPreviewBox(properties->bigPicture, window, &properties->exonViewRect, &properties->highlightRect);
    }
  
  return TRUE;
}

static void onSizeAllocateExonView(GtkWidget *exonView, GtkAllocation *allocation, gpointer data)
{
  DEBUG_ENTER("onSizeAllocateExonView");

  calculateExonViewBorders(exonView);
  
  DEBUG_EXIT("onSizeAllocateExonView returning");
}


static gboolean onButtonPressExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->button == 1) /* left button */
    {
      /* If we clicked on top of an exon, select that msp */
      guint modifiers = gtk_accelerator_get_default_mod_mask();
      const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
      const gboolean shiftModifier = ((event->state & modifiers) == GDK_SHIFT_MASK);

      handled = selectClickedExon(exonView, event, ctrlModifier, shiftModifier);
    }
  
  ExonViewProperties *properties = exonViewGetProperties(exonView);
  BigPictureProperties *bpProperties = bigPictureGetProperties(properties->bigPicture);

  if (event->button == 2 ||
      (event->button == 1 && !handled && 
       (event->type == GDK_2BUTTON_PRESS || 
        clickedInRect(event, &properties->highlightRect, bpProperties->highlightBoxMinWidth))))
    {
      /* Draw the preview box (draw it on the other big picture components as well) */
      int x = event->x;
      
      if (event->button == 1 && event->type == GDK_BUTTON_PRESS)
        x = properties->highlightRect.x + properties->highlightRect.width / 2;
      
      showPreviewBox(exonViewGetBigPicture(exonView), event->x, TRUE, x - event->x);
      handled = TRUE;
    }
  
  return handled;
}


static gboolean onButtonReleaseExonView(GtkWidget *exonView, GdkEventButton *event, gpointer data)
{
  if (event->button == 1 || event->button == 2) /* left or middle button */
    {
      ExonViewProperties *properties = exonViewGetProperties(exonView);
      acceptAndClearPreviewBox(exonViewGetBigPicture(exonView), event->x, &properties->exonViewRect, &properties->highlightRect);
    }
  
  return TRUE;
}


static gboolean onMouseMoveExonView(GtkWidget *exonView, GdkEventMotion *event, gpointer data)
{
  if ((event->state & GDK_BUTTON1_MASK) || (event->state & GDK_BUTTON2_MASK)) /* left or middle button */
    {
      /* Draw a preview box at the mouse pointer location */
      showPreviewBox(exonViewGetBigPicture(exonView), event->x, FALSE, 0);
    }
  
  return TRUE;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the grid.
 * This scrolls the position of the highlight box, i.e. it scrolls the display
 * range in the detail view. */
static gboolean onScrollExonView(GtkWidget *exonView, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
    {
      case GDK_SCROLL_LEFT:
	{
          scrollBigPictureLeftStep(exonViewGetBigPicture(exonView));
	  handled = TRUE;
	  break;
	}
	
      case GDK_SCROLL_RIGHT:
	{
          scrollBigPictureRightStep(exonViewGetBigPicture(exonView));
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



/***********************************************************
 *                       Initialisation                    *
 ***********************************************************/

/* Create the part of the view that will show the exons */
GtkWidget *createExonView(GtkWidget *bigPicture, const BlxStrand currentStrand)
{
  GtkWidget *exonView = gtk_layout_new(NULL, NULL);
  gtk_widget_set_name(exonView, BIG_PICTURE_EXON_VIEW_NAME);

  /* Connect signals */
  gtk_widget_add_events(exonView, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(exonView, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(exonView, GDK_POINTER_MOTION_MASK);
  
  g_signal_connect(G_OBJECT(exonView),	"expose-event",		G_CALLBACK(onExposeExonView),	      NULL);
  g_signal_connect(G_OBJECT(exonView),	"size-allocate",	G_CALLBACK(onSizeAllocateExonView),   NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-press-event",   G_CALLBACK(onButtonPressExonView),    NULL);
  g_signal_connect(G_OBJECT(exonView),	"button-release-event", G_CALLBACK(onButtonReleaseExonView),  NULL);
  g_signal_connect(G_OBJECT(exonView),	"motion-notify-event",  G_CALLBACK(onMouseMoveExonView),      NULL);
  g_signal_connect(G_OBJECT(exonView),  "scroll-event",	        G_CALLBACK(onScrollExonView),         NULL);

  exonViewCreateProperties(exonView, bigPicture, currentStrand);

  return exonView;
}
