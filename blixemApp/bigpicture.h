/*  File: bigpicture.h
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
 * Description: The "big picture" section displays an overview of the entire
 *              portion of reference sequence that was passed to Blixem. The
 *              reference sequence coordinates are shown along the top.
 *
 *              Alignments are shown as lines and are placed in a scaled grid
 *              (one for each strand) against their %ID (see bigpicturegrid.h).
 *
 *              Exons and introns are drawn in standard notation in the 
 *              transcript section (see exonview.h).
 *----------------------------------------------------------------------------
 */

#ifndef _big_picture_included_
#define _big_picture_included_

#include <blixemApp/bigpicturegrid.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/blixem_.h>


typedef struct _BigPictureProperties
  {
    GtkWidget *blxWindow;	/* The main blixem window that this grid belongs to */
    GtkWidget *header;		/* The grid header */
    GtkWidget *coverageView;    /* The depth-coverage view */
    GtkWidget *fwdStrandGrid;	/* The grid that displays the forward ref seq strand */
    GtkWidget *revStrandGrid;	/* The grid that displays the reverse ref seq strand */
    GtkWidget *fwdExonView;	/* The section showing the exons for the forward ref seq strand */
    GtkWidget *revExonView;	/* The section showing the exons for the reverse ref seq strand */
    
    GSList *roundValues;	/* List of "nice" values to round to, for the display values in the grid header */
    int initialZoom;		/* Multiple to multiply the detail view display range by to get the initial big picture display range */
    
    IntRange displayRange;	/* The currently-displayed range in the big picture */
    
    int numHCells;		/* The number of cells in the grid horizontally */
    int basesPerCell;		/* The number of bases show per horizontal cell */
    int roundTo;		/* The number of bases to round grid lines to the nearest multiple of */

    int numVCells;		/* The number of cells in the grid vertically */
    gdouble idPerCell;		/* The percent ID to show per vertical cell */
    DoubleRange percentIdRange;	/* The max and min %ID values displayed */

    gboolean displayPreviewBox; /* Whether to display the preview box */
    int previewBoxCentre;	/* The base that the preview box is centered on */
    int previewBoxOffset;       /* Can be used to offset the actual preview box position from previewBoxCentre */

    int leftBorderChars;	/* The number of characters in the left border of the big picture grids */
    gdouble charWidth;		/* The width of the characters in the big picture grids */
    gdouble charHeight;		/* The height of the characters in the big picture grids */
    
    int highlightBoxMinWidth;
    int previewBoxLineWidth;
    int highlightBoxYPad;       /* Vertical padding between the highlight box and the grid */
  } BigPictureProperties;

typedef struct _GridHeaderProperties
  {
    GtkWidget *bigPicture;  /* The big picture view that this header belongs to */
    GtkWidget *refButton;   /* A reference button, so we can query properties like its height */
    
    GdkRectangle headerRect; /* The actual drawing area where we'll draw the labels */
    int numHeaderLines;	    /* The number of lines of text in the header labels */
    int markerHeight;	    /* The height of the marker lines between the grid and the header labels */
    int headerYPad;	    /* Y padding around the header buttons */
  } GridHeaderProperties;


/* Public function declarations */
void			      setGdkColorBlack(GdkColor *color);
void			      setGdkColorBlue(GdkColor *color);
void			      setGdkColorYellow(GdkColor *color);
void			      setGdkColorCyan(GdkColor *color);

BigPictureProperties*	      bigPictureGetProperties(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetBlxWindow(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveGrid(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetFwdExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetRevExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetCoverageView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetActiveExonView(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetInactiveExonView(GtkWidget *bigPicture);
gboolean		      bigPictureGetDisplayRev(GtkWidget *bigPicture);
IntRange*		      bigPictureGetDisplayRange(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetGridHeader(GtkWidget *bigPicture);
GtkWidget*		      bigPictureGetDetailView(GtkWidget *bigPicture);
BlxSeqType		      bigPictureGetSeqType(GtkWidget *bigPicture);
int			      bigPictureGetNumFrames(GtkWidget *bigPicture);
gdouble			      bigPictureGetIdPerCell(GtkWidget *bigPicture);
DoubleRange*		      bigPictureGetPercentIdRange(GtkWidget *bigPicture);
int			      bigPictureGetNumVCells(GtkWidget *bigPicture);
BlxViewContext*		      bigPictureGetContext(GtkWidget *bigPicture);

gboolean                      bigPictureSetIdPerCell(GtkWidget *bigPicture, const gdouble idPerCell);
gboolean		      bigPictureSetMaxPercentId(GtkWidget *bigPicture, const gdouble newValue);
gboolean		      bigPictureSetMinPercentId(GtkWidget *bigPicture, const gdouble newValue);

void			      calculateGridHeaderBorders(GtkWidget *header);
void			      refreshBigPictureDisplayRange(GtkWidget *bigPicture, const gboolean keepCentered);
void                          calculateBigPictureCellSize(GtkWidget *bigPicture, BigPictureProperties *properties);
void			      calculateNumVCells(GtkWidget *bigPicture);
void			      bigPictureRedrawAll(GtkWidget *bigPicture);
void                          bigPicturePrepareForPrinting(GtkWidget *bigPicture);

void                          scrollBigPictureLeftStep(GtkWidget *bigPicture);
void                          scrollBigPictureRightStep(GtkWidget *bigPicture);
void                          drawPreviewBox(GtkWidget *bigPicture, GdkDrawable *drawable, GdkRectangle *displayRect, GdkRectangle *highlightRect);
void                          showPreviewBox(GtkWidget *bigPicture, const int x, const gboolean bOffset, const int offset);
void                          acceptAndClearPreviewBox(GtkWidget *bigPicture, const int xCentre, GdkRectangle *displayRect, GdkRectangle *highlightRect);

gint			      bigPictureGetCellHeight(GtkWidget *bigPicture);

void			      drawVerticalGridLines(GdkRectangle *drawingRect, GdkRectangle *highlightRect,
						    const int yPadding, BlxViewContext *bc, 
						    BigPictureProperties *bpProperties, GdkDrawable *drawable);

void			      drawHorizontalGridLines(GtkWidget *widget, GtkWidget *bigPicture,
						      GdkRectangle *drawingRect, BlxViewContext *bc,
						      BigPictureProperties *bpProperties, GdkDrawable *drawable,
						      const gint numCells, const gdouble rangePerCell, 
						      const gdouble maxVal, const gboolean abbrev, const char *unit);

void			      calculateHighlightBoxBorders(GdkRectangle *drawingRect, GdkRectangle *highlightRect,
							   GtkWidget *bigPicture, const int yPadding);

void			      zoomBigPicture(GtkWidget *bigPicture, const gboolean zoomIn);
void			      zoomWholeBigPicture(GtkWidget *bigPicture);

int			      getLeftCoordFromCentre(const int centreCoord, 
						     const int width, 
						     const GdkRectangle *outerRect);

int			      getRightCoordFromCentre(const int centreCoord, 
						      const int width, 
						      const GdkRectangle *outerRect);

void			      refreshGridOrder(GtkWidget *bigPicture);

GtkWidget*		      createBigPicture(GtkWidget *blxWindow,
					       GtkContainer *parent,
					       GtkWidget *coverageView,
					       GtkWidget **fwdStrandGrid, 
					       GtkWidget **revStrandGrid,
                                               const IntRange* const initRange,
                                               const IntRange* const fullRange,
					       const int bigPictZoom,
					       const gdouble lowestId);


#endif /* _big_picture_included_ */
