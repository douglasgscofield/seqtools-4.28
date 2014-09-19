/*  File: sequencecellrenderer.h
 *  Author: Gemma Barson, 2009-10-15
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
 * Description: A custom renderer to render cells in the detail-view trees.
 *              Each row in the tree contains a match sequence, and the 
 *              renderer draws the portion of that sequence that is within the
 *              current display range.
 *
 *              Each base is colored to indicate whether it is a match, 
 *              mismatch or conserved (similar) match.
 *              Markers are drawn to indicate deletions or insertions.
 *
 *              This renderer was originally designed for renderer the cell
 *              that draws the sequence data, hence the name, but it has since
 *              been extended to draw the plain text cells too because it
 *              colors the background and sets the correct font.
 *----------------------------------------------------------------------------
 */


#ifndef _sequence_cell_renderer_included_
#define _sequence_cell_renderer_included_

#include <blixemApp/blixem_.h>
#include <gtk/gtk.h>
#include <gtk/gtkcellrenderertext.h>

/* Some boilerplate GObject type check and type cast macros.
 *  'klass' is used here instead of 'class', because 'class'
 *  is a c++ keyword */

#define SEQUENCE_CELL_RENDERER_TYPE             (sequence_cell_renderer_get_type())
#define SEQUENCE_CELL_RENDERER(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRenderer))
#define SEQUENCE_CELL_RENDERER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRendererClass))
#define IS_SEQUENCE_CELL_RENDERER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), SEQUENCE_CELL_RENDERER_TYPE))
#define IS_SEQUENCE_CELL_RENDERER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass),  SEQUENCE_CELL_RENDERER_TYPE))
#define SEQUENCE_CELL_RENDERER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRendererClass))


/* Define the property names that we will use to set the data in the renderer */
#define RENDERER_TEXT_PROPERTY          "text"
#define RENDERER_SEQUENCE_PROPERTY      "sequence"
#define RENDERER_DATA_PROPERTY          "data"



/* SequenceCellRenderer: Our custom cell renderer
 *   structure. Extend according to need */

typedef struct _SequenceCellRenderer
{
  GtkCellRenderer   parent;
  
  /* The cell renderer can be used to render a match sequence or plain text */
  char *text;       /* generic text property */
  GList *mspGList;  /* property for the sequence column. Contains the MSP(s) to be displayed in this row */
  GList *data;      /* property for data that is set for every column */
  
} SequenceCellRenderer;


typedef struct _SequenceCellRendererClass
{
  GtkCellRendererTextClass  parent_class;
} SequenceCellRendererClass;


GType                sequence_cell_renderer_get_type (void);
GtkCellRenderer     *sequence_cell_renderer_new (void);
int		     rendererGetCellBackgroundPadding(GtkCellRenderer *cell);

#endif /* _sequence_cell_renderer_included_ */
