/*  File: bigpicturegrid.h
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
 * Description: Creates a grid on the Big Picture section. A grid shows the
 *              alignments on a particular strand of the reference sequence.
 *              The x coord of the grid is the reference sequence coord and the
 *              y coord is the %ID of the alignment. An alignment is shown as
 *              a horizontal line between the start and end coords of the 
 *              match, at the y coord corresponding to the relevant %ID.
 *----------------------------------------------------------------------------
 */

#ifndef _big_picture_grid_included_
#define _big_picture_grid_included_


#include <gtk/gtk.h>
#include <blixemApp/blixem_.h>
#include <seqtoolsUtils/utilities.h>


#define BIG_PICTURE_GRID_NAME		"BigPictureGrid"


typedef struct _GridProperties
  {
    GtkWidget *bigPicture;   /* The big picture that this grid belongs to */
    
    BlxStrand strand;	     /* Whether this grid shows the fwd or rev strand of the ref sequence. */
    
    int mspLineHeight;	     /* The height of the msp lines */
    
    int gridYPadding;	     /* The y padding around the grid */
    
    gulong exposeHandlerId;  /* The handler ID for the expose event */
    gboolean ignoreSizeAlloc; /* Flag to indicate that we should ignore size allocation events */
    
    GdkRectangle gridRect;
    GdkRectangle displayRect;
    GdkRectangle highlightRect;
  } GridProperties;


/* Public function declarations */
GridProperties*	    gridGetProperties(GtkWidget *widget);
BlxStrand	    gridGetStrand(GtkWidget *grid);
GtkWidget*	    gridGetBigPicture(GtkWidget *grid);

void		    calculateGridBorders(GtkWidget *grid);
void		    calculateGridHighlightBoxBorders(GtkWidget *grid);

void		    callFuncOnAllBigPictureGrids(GtkWidget *widget, 
						 gpointer data);

gint		    convertValueToGridPos(GtkWidget *grid, 
					  const gdouble value);

void                gridPrepareForPrinting(GtkWidget *grid);

GtkWidget*	    createBigPictureGrid(GtkWidget *bigPicture, 
					 BlxStrand strand);

#endif /* _big_picture_grid_included_ */
