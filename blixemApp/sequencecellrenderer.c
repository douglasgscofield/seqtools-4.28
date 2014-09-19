/*  File: sequencecellrenderer.c
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
 * Description: See sequencecellrenderer.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/sequencecellrenderer.h>
#include <blixemApp/detailview.h>
#include <blixemApp/detailviewtree.h>
#include <blixemApp/blxwindow.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/utilities.h>
#include <gtk/gtkcellrenderertext.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>


#define SEQUENCE_CELL_RENDERER_NAME	"SequenceCellRenderer"
#define GAP_WIDTH_AS_FRACTION		0.375	/* multiplier used to get the width of the "gap" marker based on a fraction of char width */
#define MIN_GAP_WIDTH			2

typedef struct _RenderData
  {
    BlxViewContext *bc;
    GdkRectangle *cell_area;
    GdkWindow *window;
    GtkStateType state;
    GdkGC *gc;
    const BlxStrand qStrand;
    const int qFrame;
    const int selectedBaseIdx;
    gboolean seqSelected;
    const int cellXPadding;
    const int cellYPadding;
    const gdouble charWidth;
    const gdouble charHeight;
    const IntRange* const displayRange;
    const gboolean highlightDiffs;
    GdkDrawable *drawable;
    GtkWidget *blxWindow;
    GdkColor *exonColor;
    GdkColor *exonColorSelected;
    GdkColor *cdsColor;
    GdkColor *cdsColorSelected;
    GdkColor *utrColor;
    GdkColor *utrColorSelected;
    GdkColor *crosshatchColor;
    GdkColor *crosshatchColorSelected;
    GdkColor *insertionColor;
    GdkColor *insertionColorSelected;
    GdkColor *matchColor;
    GdkColor *matchColorSelected;
    GdkColor *consColor;
    GdkColor *consColorSelected;
    GdkColor *mismatchColor;
    GdkColor *mismatchColorSelected;
    GdkColor *unalignedSeqColor;
    GdkColor *unalignedSeqColorSelected;
    GdkColor *exonBoundaryColorStart;
    GdkColor *exonBoundaryColorEnd;
    GdkColor *polyAColor;
    GdkColor *polyAColorSelected;
    GdkColor *clipMarkerColor;
    int exonBoundaryWidth;
    GdkLineStyle exonBoundaryStyle;
    GdkLineStyle exonBoundaryStylePartial;
    gboolean limitUnalignedBases;
    int numUnalignedBases;
  } RenderData;



/* Properties */
enum 
{
  PROP_0,
  
  PROP_TEXT,
  PROP_MSP,
  PROP_DATA
};


/* Some boring function declarations: GObject type system stuff */
static void     sequence_cell_renderer_init       (SequenceCellRenderer      *cellrenderersequence);
static void     sequence_cell_renderer_class_init (SequenceCellRendererClass *klass);
static void     sequence_cell_renderer_finalize (GObject *gobject);

static gboolean mspGetVisibleRange(MSP *msp, RenderData *data, IntRange *result);
void		drawAllVisibleExonBoundaries(GtkWidget *tree, RenderData *data);

static void mspDrawSequenceText(GtkWidget *tree,
                                gchar *displayText, 
                                const IntRange* const segmentRange,
                                RenderData *data);


static void segmentGetCoordsForBaseIdx(const int segmentIdx, 
                                       const IntRange* const segmentRange,
                                       RenderData *data,
                                       int *x, 
                                       int* y);


/* These functions are the heart of our custom cell renderer: */
static void     sequence_cell_renderer_get_size   (GtkCellRenderer            *cell,
                                                          GtkWidget                  *widget,
                                                          GdkRectangle               *cell_area,
                                                          gint                       *x_offset,
                                                          gint                       *y_offset,
                                                          gint                       *width,
                                                          gint                       *height);
static void     sequence_cell_renderer_render     (GtkCellRenderer            *cell,
                                                          GdkDrawable                *window,
                                                          GtkWidget                  *tree,
                                                          GdkRectangle               *background_area,
                                                          GdkRectangle               *cell_area,
                                                          GdkRectangle               *expose_area,
                                                          GtkCellRendererState       flags);
static void sequence_cell_renderer_get_property (GObject      *object,
					 guint         param_id,
					 GValue	      *value,
					 GParamSpec   *pspec);
static void sequence_cell_renderer_set_property (GObject      *object,
					 guint         param_id,
					 const GValue *value,
					 GParamSpec   *pspec);

static   gpointer parent_class;



/***************************************************************************/

/* Custom version of gtk_paint_layout that draws the same thing to two GdkDrawables */
static void paintLayout2(GtkStyle *style,
			 GdkDrawable *drawable1,
			 GdkDrawable *drawable2,
			 GtkStateType state_type,
			 gboolean use_text,
			 GdkRectangle *area,
			 GtkWidget *widget,
			 gchar *detail,
			 gint x,
			 gint y,
			 PangoLayout *layout)
{
  if (drawable1)
    gtk_paint_layout (style, drawable1, state_type, use_text, area, widget, detail, x, y, layout);
  
  if (drawable2)
    gtk_paint_layout (style, drawable2, state_type, use_text, area, widget, detail, x, y, layout);
}

/* Custom version of gdk_draw_line that draws the same thing to two GdkDrawables */
void drawLine2(GdkDrawable *drawable1,
               GdkDrawable *drawable2,
               GdkGC *gc,
               gint x1,
               gint y1,
               gint x2,
               gint y2)
{
  if (drawable1)
    {
//      cairo_t *cr = gdk_cairo_create(drawable1);
//      cairo_set_line_width(cr, 0.5);
//      cairo_move_to(cr, x1 + 0.5, y1);
//      cairo_line_to(cr, x2 + 0.5, y2);
//      cairo_stroke(cr); 
//      cairo_destroy(cr);
      gdk_draw_line(drawable1, gc, x1, y1, x2, y2);
    }
  
  if (drawable2)
    {
//      cairo_t *cr = gdk_cairo_create(drawable2);
//      cairo_set_line_width(cr, 0.5);
//      cairo_move_to(cr, x1 + 0.5, y1);
//      cairo_line_to(cr, x2 + 0.5, y2);
//      cairo_stroke(cr); 
//      cairo_destroy(cr);
      
      gdk_draw_line(drawable2, gc, x1, y1, x2, y2);
    }
}

/* Custom version of gdk_draw_rectangle that draws the same thing to two GdkDrawables */
void drawRectangle2(GdkDrawable *drawable1,
		    GdkDrawable *drawable2,
		    GdkGC *gc,
		    gboolean filled,
		    gint x,
		    gint y,
		    gint width,
		    gint height)
{
//  cairo_t *cr = gdk_cairo_create(drawable1);
//  cairo_rectangle(cr, x, y, width, height);
//  cairo_fill(cr);
//  cairo_destroy(cr);
  
  if (drawable1)
    {
      gdk_draw_rectangle(drawable1, gc, filled, x, y, width, height);
    }
  
  if (drawable2)
    {
      gdk_draw_rectangle(drawable2, gc, filled, x, y, width, height);
    }
}


/* Fill a rectangle with crosshatch lines */
void drawCrosshatchRectangle(GdkDrawable *drawable1,
                             GdkDrawable *drawable2,
                             GdkGC *gc,
                             gint x,
                             gint y,
                             gint width,
                             gint height)
{
  /* Draw diagonal lines from top-right to bottom-left. */
  drawLine2(drawable1, drawable2, gc, x + width, y, x, y + height);
  drawLine2(drawable1, drawable2, gc, x + width / 2, y, x, y + height / 2);
  drawLine2(drawable1, drawable2, gc, x + width, y + height / 2, x + width / 2, y + height);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_get_type: here we register our type with
 *                                          the GObject type system if we
 *                                          haven't done so yet. Everything
 *                                          else is done in the callbacks.
 *
 ***************************************************************************/

GType
sequence_cell_renderer_get_type (void)
{
  static GType cell_sequence_type = 0;
  
  if (cell_sequence_type == 0)
    {
      static const GTypeInfo cell_sequence_info =
      {
	sizeof (SequenceCellRendererClass),
	NULL,                                                     /* base_init */
	NULL,                                                     /* base_finalize */
	(GClassInitFunc) sequence_cell_renderer_class_init,
	NULL,                                                     /* class_finalize */
	NULL,                                                     /* class_data */
	sizeof (SequenceCellRenderer),
	0,                                                        /* n_preallocs */
	(GInstanceInitFunc) sequence_cell_renderer_init,
      };
      
      /* Derive from GtkCellRenderer */
      cell_sequence_type = g_type_register_static (GTK_TYPE_CELL_RENDERER,
						   "SequenceCellRenderer",
						   &cell_sequence_info,
                                                   0);
    }
  
  return cell_sequence_type;
}


/***************************************************************************
 *
 *  sequence_cell_renderer_init: set some default properties of the
 *                                      parent (GtkCellRendererText).
 *
 ***************************************************************************/

static void
sequence_cell_renderer_init (SequenceCellRenderer *cellrenderersequence)
{
  GTK_CELL_RENDERER(cellrenderersequence)->mode = GTK_CELL_RENDERER_MODE_INERT;
  GTK_CELL_RENDERER(cellrenderersequence)->xalign = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->yalign = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->xpad = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->ypad = 0;

  cellrenderersequence->data = NULL;
  cellrenderersequence->mspGList = NULL;
  cellrenderersequence->text = NULL;
}


/***************************************************************************
 *
 *  sequence_cell_renderer_class_init:
 *
 *  If you want cells that can be activated on their own (ie. not
 *  just the whole row selected) or cells that are editable, you
 *  will need to override 'activate' and 'start_editing'.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_class_init (SequenceCellRendererClass *klass)
{
  GtkCellRendererClass *cell_class   = GTK_CELL_RENDERER_CLASS(klass);
  GObjectClass         *object_class = G_OBJECT_CLASS(klass);
  
  parent_class           = g_type_class_peek_parent (klass);
  object_class->finalize = sequence_cell_renderer_finalize;

  /* Override the two crucial functions that are the heart
   *   of a cell renderer in the parent class */
  cell_class->get_size = sequence_cell_renderer_get_size;
  cell_class->render   = sequence_cell_renderer_render;
  
  /* Override functions to get/set properties */
  object_class->get_property = sequence_cell_renderer_get_property;
  object_class->set_property = sequence_cell_renderer_set_property;
  
  g_object_class_install_property (object_class,
                                   PROP_DATA,
                                   g_param_spec_pointer (RENDERER_DATA_PROPERTY,
                                                         "Data",
                                                         "Pointer to the msp",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_MSP,
                                   g_param_spec_pointer (RENDERER_SEQUENCE_PROPERTY,
                                                         "Sequence",
                                                         "Pointer to an msp whose sequence to display",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_TEXT,
                                   g_param_spec_string (RENDERER_TEXT_PROPERTY,
                                                        "Text",
                                                        "Text field",
							NULL,
                                                        G_PARAM_WRITABLE));
}


/***************************************************************************
 *
 *  sequence_cell_renderer_finalize: free any resources here
 *
 ***************************************************************************/

static void
sequence_cell_renderer_finalize (GObject *object)
{
  /*
   SequenceCellRenderer *cellrenderersequence = SEQUENCE_CELL_RENDERER(object);
   */
  
  /* Free any dynamically allocated resources here */
  
  (* G_OBJECT_CLASS (parent_class)->finalize) (object);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_new: return a new cell renderer instance
 *
 ***************************************************************************/

GtkCellRenderer *
sequence_cell_renderer_new (void)
{
  return (GtkCellRenderer*)g_object_new(SEQUENCE_CELL_RENDERER_TYPE, NULL);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_set_property: 
 *
 ***************************************************************************/

static void setAllPropertiesNull(SequenceCellRenderer *renderer)
{
  if (renderer->text)
    g_free(renderer->text);

  renderer->text = NULL;
  renderer->mspGList = NULL;
}

static void
sequence_cell_renderer_set_property (GObject      *object,
				     guint         param_id,
				     const GValue *value,
				     GParamSpec   *pspec)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(object);
  
  switch (param_id)
  {
    case PROP_DATA:
      /* Additional data. DON'T reset other properties. */
      renderer->data = (GList*)g_value_get_pointer(value);
      break;
      
    case PROP_MSP:
      /* We set either the msp data or the text data, not both;
       * so make sure nothing else is set */
      setAllPropertiesNull(renderer);
      renderer->mspGList = (GList*)g_value_get_pointer(value);
      break;

    case PROP_TEXT:
      /* We set either the msp data or the text data, not both;
       * so make sure nothing else is set */
      setAllPropertiesNull(renderer);
      renderer->text = g_strdup(g_value_get_string(value));
      break;
      
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
      break;
  };
}


/***************************************************************************
 *
 *  sequence_cell_renderer_get_property: 
 *
 ***************************************************************************/

static void
sequence_cell_renderer_get_property (GObject      *object,
					 guint         param_id,
					 GValue	      *value,
					 GParamSpec   *pspec)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(object);
  switch (param_id)
  {
    case PROP_DATA:
      g_value_set_pointer(value, renderer->data);
      break;

    case PROP_MSP:
      g_value_set_pointer(value, renderer->mspGList);
      break;
      
    case PROP_TEXT:
      g_value_set_string(value, renderer->text);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
      break;
  };
}


static void
get_size (GtkCellRenderer *cell,
	  GtkWidget       *widget,
	  GdkRectangle    *cell_area,
	  gint            *x_offset,
	  gint            *y_offset,
	  gint            *width,
	  gint            *height)
{
  if (height)
    *height = cell->height;
  
  if (width)
    *width = cell->width;
  
  if (cell_area)
    cell_area->height = cell->height;
}


/* Render function for a cell that contains simple text */
static void rendererDrawSimpleText(SequenceCellRenderer *renderer, 
                                   GtkWidget *tree,
                                   GdkWindow *window, 
                                   GtkStateType state, 
                                   GdkRectangle *cell_area)
{
  gchar *displayText = renderer->text;
  PangoLayout *layout = gtk_widget_create_pango_layout(tree, displayText);

  PangoFontDescription *font_desc = pango_font_description_copy(tree->style->font_desc);
  pango_layout_set_font_description(layout, font_desc);
  pango_font_description_free(font_desc);

  paintLayout2(tree->style,
	       window,
	       widgetGetDrawable(tree),
	       state,
	       TRUE,
	       NULL,
	       tree,
	       NULL, 
	       cell_area->x - treeGetCellXPadding(tree),
	       cell_area->y - treeGetCellYPadding(tree),
	       layout);
  
  g_object_unref(layout);
}


/* The given renderer is an MSP. This function checks if there is a base index
 * selected and, if so, colors the background for that base with the given color. */
static void highlightSelectedBase(const int selectedBaseIdx,
				  GdkColor *highlightColor,
				  RenderData *data)
{
  if (selectedBaseIdx != UNSET_INT && valueWithinRange(selectedBaseIdx, data->displayRange))
    {
      /* Convert the display-range index to a 0-based index for the section of sequence displayed */
      const int segmentIdx = selectedBaseIdx - data->displayRange->min;
      
      int x, y;
      segmentGetCoordsForBaseIdx(segmentIdx, data->displayRange, data, &x, &y);

      gdk_gc_set_foreground(data->gc, highlightColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, roundNearest(data->charWidth), roundNearest(data->charHeight));
    }
}


/* Utility to get the exon color based on whether it is CDS/UTR and whether it is selected */
static GdkColor* exonGetFillColor(const MSP* const msp, const gboolean isSelected, RenderData *data)
{
  GdkColor *result = NULL;

  if (msp->type == BLXMSP_EXON)
    {
      result = (isSelected ? data->exonColorSelected : data->exonColor);
    }
  else if (msp->type == BLXMSP_CDS)
    {
      result = (isSelected ? data->cdsColorSelected : data->cdsColor);
    }
  else
    {
      result = (isSelected ? data->utrColorSelected : data->utrColor);
    }
    
  return result;
}


/* Returns true if the min coord (or max coord, if 'start' is false) is a
 * partial codon, i.e. if it does not start at base 1 (if start) or end at
 * base 3 (if end)... or vice versa if the display is reversed */
static gboolean exonCoordIsPartialCodon(const MSP* const msp, const gboolean start, RenderData *data, int *displayIdxOut)
{
  /* To be a complete codon, if we're at the start the the coord must be base 1 
   * or if we're at the end then coord must be base 3 */
  const int reqdBase = (start != data->bc->displayRev) ? 1 : data->bc->numFrames;
  
  /* Calculate the actual base number of the start/end coord */
  const int dnaIdx = start ? msp->qRange.min : msp->qRange.max;
  
  int baseNum = UNSET_INT;
  const int frame = mspGetRefFrame(msp, data->bc->seqType);
  
  const int displayIdx = convertDnaIdxToDisplayIdx(dnaIdx, data->bc->seqType, frame, 
    data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange, &baseNum);
  
  if (displayIdxOut)
    *displayIdxOut = displayIdx;

  gboolean isPartial = (baseNum != reqdBase);
  return isPartial;
}


/* Highlight the start (min) peptide of the given msp if it is a partial codon, i.e. if it
 * does not start at base 1 (or end at base 3, if the display is reversed). OR highlight the end
 * (max) codon if 'start' is FALSE. The x/y coords give the top left corner of the peptide. */
static void exonHighlightPartialCodons(const MSP* const msp, 
                                       const gboolean start, 
                                       const int x, 
                                       const int y, 
                                       RenderData *data)
{
  int displayIdx = UNSET_INT;
  gboolean isPartial = exonCoordIsPartialCodon(msp, start, data, &displayIdx);
  
  if (isPartial && valueWithinRange(displayIdx, data->displayRange)) 
    {
      /* It's not a complete codon, so highlight this base */
      const gboolean coordSelected = (data->selectedBaseIdx == displayIdx);
      const gboolean isSelected = (data->seqSelected != coordSelected);

      GdkColor *color = (isSelected ? data->crosshatchColorSelected : data->crosshatchColor);
      gdk_gc_set_foreground(data->gc, color);

      drawCrosshatchRectangle(data->window, data->drawable, data->gc, x, y, data->charWidth, data->charHeight);
    }
}



/* The given renderer is an MSP that is an exon. This function draws the exon
 * or part of the exon that is in view, if it is within the current display range. */
static void drawExon(SequenceCellRenderer *renderer,
		     MSP *msp,
		     GtkWidget *tree,
		     RenderData *data)
{
  if (mspLayerIsVisible(msp))
    {
      IntRange segmentRange = {UNSET_INT, UNSET_INT};
      
      if (!mspGetVisibleRange(msp, data, &segmentRange))
        {
          return;
        }

      const int segmentLen = segmentRange.max - segmentRange.min + 1;

      int x, y;
      segmentGetCoordsForBaseIdx(0, &segmentRange, data, &x, &y);
      const int width = ceil((gdouble)segmentLen * data->charWidth);

      /* Just draw one big rectangle the same color for the whole thing. Color depends if row selected. */
      GdkColor *color = exonGetFillColor(msp, data->seqSelected, data);
      gdk_gc_set_foreground(data->gc, color);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, width, roundNearest(data->charHeight));
      
      /* If a base is selected, highlight it. Its color depends on whether it the base is within the exon range or not. */
      if (data->selectedBaseIdx != UNSET_INT && valueWithinRange(data->selectedBaseIdx, &segmentRange))
        {
          /* Negate the color if double-selected (i.e. if the row is selected as well) */
          GdkColor *color = exonGetFillColor(msp, !data->seqSelected, data);
          highlightSelectedBase(data->selectedBaseIdx, color, data);
        }
      
      /* If the start or end index is not a full codon, highlight it in a different color */
      exonHighlightPartialCodons(msp, !data->bc->displayRev, x, y, data);
      exonHighlightPartialCodons(msp, data->bc->displayRev, x + width - data->charWidth, y, data);
      
      //drawAllVisibleExonBoundaries(tree, data);
    }
} 


/* Get the base from the given sequence at the given index. Converts to
 * upper or lower case as appropriate for the sequence type. The given index is
 * assumed to be 1-based. If the given BlxSequence has null sequence data, then
 * this function returns a padding character instead */
static char blxSeqGetMatchSeqBase(BlxSequence *blxSeq, const int sIdx, const BlxSeqType seqType)
{
  char result = SEQUENCE_CHAR_PAD;
  
  const char *sequence = blxSequenceGetSequence(blxSeq);
  
  if (sequence && sIdx <= (int)strlen(sequence))
    {
      result = sequence[sIdx - 1];
      result = convertBaseToCorrectCase(result, seqType);
    }

  return result;
}


/* Color in the background of a particular base in the given match sequence. Returns
 * the calculated index into the match sequence (for effiency, so that we don't have to
 * recalculate it later on). Returns the equivalent index in the subject sequence, or 
 * UNSET_INT if there is none. */
static void mspDrawBaseBg(MSP *msp,
                          const int segmentIdx, 
                          const IntRange* const segmentRange,
                          char *refSeqSegment, 
                          RenderData *data,
                          const int x,
                          const int y,
                          gchar *displayText,
                          int *sIdx,
                          int *qIdx)
{
  char sBase = '\0';
  GdkColor *baseBgColor = NULL;
  
  /* From the segment index, find the display index and the ref seq coord */
  const int displayIdx = segmentRange->min + segmentIdx;
  *qIdx = convertDisplayIdxToDnaIdx(displayIdx, data->bc->seqType, data->qFrame, 1, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange);
  
  /* Find the match-sequence coord at this ref-seq coord */
  *sIdx = mspGetMatchCoord(msp, *qIdx, data->seqSelected, data->numUnalignedBases, data->bc);
  
  /* Highlight the base if its base index is selected, or if its sequence is selected.
   * (If it is selected in both, show it in the normal color) */
  gboolean selected = (displayIdx == data->selectedBaseIdx) != data->seqSelected;

  if (!valueWithinRange(*qIdx, &msp->qRange))
    {
      /* We're outside the alignment range. There might still be a base to display if
       * we're displaying unaligned parts of the match sequence or polyA tails; otherwise, we 
       * show nothing. */
      if (*sIdx != UNSET_INT)
	{
          sBase = blxSeqGetMatchSeqBase(msp->sSequence, *sIdx, data->bc->seqType);

          if (data->bc->flags[BLXFLAG_SHOW_POLYA_SITE] && 
              (!data->bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED] || data->seqSelected) &&
              mspCoordInPolyATail(*qIdx, msp))
            {
              baseBgColor = selected ? data->polyAColorSelected : data->polyAColor;
            }
          else
            {
              baseBgColor = selected ? data->unalignedSeqColorSelected : data->unalignedSeqColor;
            }
	}
    }
  else if (*sIdx == UNSET_INT)
    {
      /* We're inside the alignment range but there is no base to display: we must be in a deletion. */
      sBase = SEQUENCE_CHAR_DELETION;
      baseBgColor = selected ? data->mismatchColorSelected : data->mismatchColor;
    }
  else
    {
      /* There is a base in the match sequence. See if it matches the ref sequence */
      sBase = blxSeqGetMatchSeqBase(msp->sSequence, *sIdx, data->bc->seqType);
      char qBase = refSeqSegment[segmentIdx];

      if (tolower(sBase) == tolower(qBase))
	{
	  /* Match */
	  baseBgColor = selected ? data->matchColorSelected : data->matchColor;
	  
	  /* If we're highlighting differences, don't show this base (show a dash instead) */
	  if (data->highlightDiffs)
	    {
	      sBase = SEQUENCE_CHAR_BLANK;
	    }
	}
      else if (data->bc->blastMode != BLXMODE_BLASTN && PAM120[aa_atob[(unsigned int)qBase]-1 ][aa_atob[(unsigned int)sBase]-1 ] > 0)
	{
	  /* 'Conserved' match (i.e. similar amino acid) */
	  baseBgColor = selected ? data->consColorSelected : data->consColor;
	}
      else
	{
	  /* Mismatch */
	  baseBgColor = selected ? data->mismatchColorSelected : data->mismatchColor;
	}
    }

  /* Draw the background color */
  if (baseBgColor)
    {
      gdk_gc_set_foreground(data->gc, baseBgColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, ceil(data->charWidth), roundNearest(data->charHeight));
    }
  
  if (sBase != '\0')
    {
      /* Add this character into the display text */
      displayText[segmentIdx] = sBase;
    }
  else
    {
      displayText[segmentIdx] = ' ';
    }
}


static PangoLayout* pangoGetLayoutFromText(gchar *displayText, GtkWidget *tree, PangoFontDescription *font_desc)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(tree, displayText);
  pango_layout_set_font_description(layout, font_desc);
  return layout;
}


/* Return the x/y coords for the top-left corner where we want to draw the base
 * with the given index in the segment, where the segment starts at the given 
 * index in the display. */
static void segmentGetCoordsForBaseIdx(const int segmentIdx, 
                                       const IntRange* const segmentRange,
                                       RenderData *data,
                                       int *x, 
                                       int* y)
{
  /* Find the start of the segment with respect to the display range */
  const int startPos = segmentRange->min - data->displayRange->min;
  
  /* Find the position of the character within the segment */
  int charIdx = startPos + segmentIdx;

  /* Calculate the coords */
  *x = data->cell_area->x - data->cellXPadding + (int)((gdouble)charIdx * data->charWidth);
  *y = data->cell_area->y - data->cellYPadding;
}


/* Draw the start/end boundaries for the given MSP if it is an exon whose start/end
 * coords are within the current display range */
static gboolean exonDrawBoundary(const MSP *msp, RenderData *rd)
{
  if (msp && msp->type == BLXMSP_EXON)
    {
      /* Get the msp's start/end in terms of the display coords */
      const IntRange* const mspRange = mspGetDisplayRange(msp);
      
      if (valueWithinRange(mspRange->min, rd->displayRange))
	{
	  /* Draw the lower index. The color and line style depend on whether it's the start or end index. */
	  GdkColor *color = rd->bc->displayRev ? rd->exonBoundaryColorEnd : rd->exonBoundaryColorStart;
	  gdk_gc_set_foreground(rd->gc, color);
	  
          /* Check if it's the boundary of a partial codon - we'll draw a dotted
           * line if it is, to indicate that the boundary is not exactly at this position. */
          const gboolean isPartial = exonCoordIsPartialCodon(msp, !rd->bc->displayRev, rd, NULL);
	  GdkLineStyle lineStyle = isPartial ? rd->exonBoundaryStylePartial : rd->exonBoundaryStyle;
          gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, lineStyle, GDK_CAP_BUTT, GDK_JOIN_MITER);

          const int idx = mspRange->min - rd->displayRange->min;
          
	  int x = UNSET_INT, y = UNSET_INT;
	  segmentGetCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);
          
          drawLine2(rd->window, rd->drawable, rd->gc, x, y, x, y + roundNearest(rd->charHeight));
	}
      
      if (valueWithinRange(mspRange->max, rd->displayRange))
	{
	  /* Draw the upper index. The color and line style depend on whether it's the start or end index. */
	  GdkColor *color = rd->bc->displayRev ? rd->exonBoundaryColorStart : rd->exonBoundaryColorEnd;
	  gdk_gc_set_foreground(rd->gc, color);
	  
          /* Check if it's the boundary of a partial codon - we'll draw a dotted
           * line if it is, to indicate that the boundary is not exactly at this position. */
          const gboolean isPartial = exonCoordIsPartialCodon(msp, rd->bc->displayRev, rd, NULL);
	  GdkLineStyle lineStyle = isPartial ? rd->exonBoundaryStylePartial : rd->exonBoundaryStyle;
          gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, lineStyle, GDK_CAP_BUTT, GDK_JOIN_MITER);
	  
	  const int idx = mspRange->max + 1 - rd->displayRange->min;

	  int x = UNSET_INT, y = UNSET_INT;
	  segmentGetCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);

          drawLine2(rd->window, rd->drawable, rd->gc, x, y, x, y + roundNearest(rd->charHeight));
	}
    }
  
  return FALSE;
}


/* Draw the boundaries of all exons in the given tree that are within the current
 * display range */
void drawAllVisibleExonBoundaries(GtkWidget *tree, RenderData *data)
{
  /* Loop through all MSPs. */
  const MSP *msp = blxWindowGetMspList(data->blxWindow);

  for ( ; msp; msp = msp->next)
    {
      if (mspIsExon(msp) && mspGetRefFrame(msp, data->bc->seqType) == data->qFrame && mspGetRefStrand(msp) == data->qStrand)
	{
	  exonDrawBoundary(msp, data);
	}
    }
}


/* Draw a vertical yellow line to indicate an insertion between two bases in the match sequence */
static void mspDrawInsertionMarker(int sIdx, 
                                   int lastFoundSIdx, 
                                   int qIdx, 
                                   int lastFoundQIdx, 
                                   int x, 
                                   int y, 
                                   RenderData *data)
{
  if (sIdx != UNSET_INT && lastFoundSIdx != UNSET_INT && abs(sIdx - lastFoundSIdx) > 1)
    {
      /* There is a gap between this index and the previous one (in a left-to-right sense),
       * so draw an insertion marker between these two bases (at the position given by x and y). */
      
      /* This is not very sophisticated - just uses a fudge factor to find a suitable width and
       * draws it half over the current base and half over the previous one. */
      int gapWidth = roundNearest(data->charWidth * GAP_WIDTH_AS_FRACTION);

      if (gapWidth < MIN_GAP_WIDTH)
	{
	  gapWidth = MIN_GAP_WIDTH;
        }

      /* No selections to worry about. Draw the whole thing in the normal color */
      gdk_gc_set_foreground(data->gc, data->insertionColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x - gapWidth/2, y, gapWidth, roundNearest(data->charHeight));
    }
}


/* Draw the given sequence text at the given coords */
static void mspDrawSequenceText(GtkWidget *tree,
                                gchar *displayText, 
                                const IntRange* const segmentRange,
                                RenderData *data)
{
  if (g_utf8_validate(displayText, -1, NULL))
    {
      /* Get the coords for the first base. The display text should have been
       * was constructed such that everything else will line up from here. */
      int x, y;
      segmentGetCoordsForBaseIdx(0, segmentRange, data, &x, &y);
      
      PangoFontDescription *font_desc = pango_font_description_copy(tree->style->font_desc);
      PangoLayout *layout = pangoGetLayoutFromText(displayText, tree, font_desc);
      pango_font_description_free(font_desc);

      if (layout)
	{
	  paintLayout2(tree->style, data->window, data->drawable, data->state, TRUE, NULL, tree, NULL, x, y, layout);
	  g_object_unref(layout);
	}
      else
	{
	  g_warning("Error creating layout while trying to display sequence:\n%s\n", displayText);
	}
    }
  else
    {
      g_critical("Invalid string constructed when trying to display sequence.\n");
    }
  
}


/* Get the range of the msp that is inside the currently displayed range. The
 * return value is FALSE if no part of the msp is visible. The result coords 
 * are in display coords. */
static gboolean mspGetVisibleRange(MSP *msp, RenderData *data, IntRange *result)
{
  gboolean found = FALSE;
  
  /* Get the full display range of the MSP (including any portions of unaligned sequence etc.) */
  const IntRange *fullRange = mspGetFullDisplayRange(msp, data->seqSelected, data->bc);
  result->min = fullRange->min;
  result->max = fullRange->max;
  
  if (rangesOverlap(result, data->displayRange))
    {
      /* Limit the returned range to the display range. */
      boundsLimitRange(result, data->displayRange, FALSE);
      found = TRUE;
    }
  else
    {
      /* No portion of the MSP range is in the display range, so return UNSET_INTs */
      result->min = UNSET_INT;
      result->max = UNSET_INT;
    }

  return found;
}


/* If the current MSP has been clipped (i.e. extends outside the current reference
 * sequence range) then draw a marker at the clip point so that the user knows there
 * is more data for the match that isn't shown. */
static void mspDrawClippedMarker(const MSP* const msp,
                                 const int qIdx, 
                                 const int segmentIdx, 
                                 const IntRange* const segmentRange,
                                 const int x, 
                                 const int y, 
                                 RenderData *data)
{
  gboolean clipStart = FALSE;
  gboolean clipEnd = FALSE;
  
  /* Check if we're at the very first/last index in the reference sequence range */
  if (qIdx == data->bc->refSeqRange.min && msp->qRange.min < qIdx)
    {
      if (data->bc->displayRev)
        {
          clipEnd = TRUE;
        }
      else
        {
          clipStart = TRUE;
        }
    }
  else if (qIdx == data->bc->refSeqRange.max && msp->qRange.max > qIdx)
    {
      if (data->bc->displayRev)
        {
          clipStart = TRUE;
        }
      else
        {
          clipEnd = TRUE;
        }
    }

  if (clipStart || clipEnd)
    {
      gdouble gapWidth = data->charWidth * GAP_WIDTH_AS_FRACTION;
      const int xPos = clipStart ? x : roundNearest((gdouble)x + data->charWidth - gapWidth);

      gdk_gc_set_foreground(data->gc, data->clipMarkerColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, xPos, y, gapWidth, roundNearest(data->charHeight));
    }
}


/* Draw the the given MSP. As well as the sequence text this also draws the background
 * color for each base based on how well it matches the reference sequence, and also
 * draws insertion markers and "clipped" markers. 
 * Optionally returns the visible range of the MSP in segmentRange_out. */
static void drawMsp(SequenceCellRenderer *renderer,
                    MSP *msp,
                    GtkWidget *tree,
                    RenderData *data)
{
  /* Extract the section of the reference sequence that we're interested in. */
  IntRange segmentRange = {UNSET_INT, UNSET_INT};
  
  if (!mspGetVisibleRange(msp, data, &segmentRange))
    {
      return;
    }

  /* The ref seq is in nucleotide coords, so convert the segment coords to nucleotide coords */
  const int coord1 = convertDisplayIdxToDnaIdx(segmentRange.min, data->bc->seqType, data->qFrame, 1, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange);
  const int coord2 = convertDisplayIdxToDnaIdx(segmentRange.max, data->bc->seqType, data->qFrame, data->bc->numFrames, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange);

  IntRange qRange;
  intrangeSetValues(&qRange, coord1, coord2);

  GError *error = NULL;
  gchar *refSeqSegment = getSequenceSegment(data->bc->refSeq,
                                            &qRange,
					    data->qStrand, 
                                            BLXSEQ_DNA,               /* ref seq is always in nucleotide coords */
					    data->bc->seqType,        /* required segment is in display coords */
					    data->qFrame, 
					    data->bc->numFrames,
					    &data->bc->refSeqRange,
					    data->bc->blastMode,
					    data->bc->geneticCode,
					    data->bc->displayRev,
					    data->bc->displayRev,
					    TRUE,
					    &error);

  if (!refSeqSegment)
    {
      g_assert(error);
      prefixError(error, "Could not draw alignment for sequence '%s'. ", mspGetSName(msp));
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
    }
    
  /* We'll populate a string with the characters we want to display as we loop through the indices. */
  const int segmentLen = strlen(refSeqSegment);
  gchar displayText[segmentLen + 1];
  displayText[0] = '\0';
  
  int lastFoundSIdx = UNSET_INT;  /* remember the last index where we found a valid base */
  int lastFoundQIdx = UNSET_INT;  /* remember the last index where we found a valid base */

  int segmentIdx = 0;
  for ( ; segmentIdx < segmentLen; ++segmentIdx)
    {
      int x = UNSET_INT, y = UNSET_INT;
      segmentGetCoordsForBaseIdx(segmentIdx, &segmentRange, data, &x, &y);
      
      /* Find the base in the match sequence and draw the background color according to how well it matches */
      int sIdx = UNSET_INT, qIdx = UNSET_INT;
      mspDrawBaseBg(msp, segmentIdx, &segmentRange, refSeqSegment, data, x, y, displayText, &sIdx, &qIdx);
      
      /* If there is an insertion (i.e. extra bases on the match sequence) between this 
       * and the previous coord, draw a marker */
      mspDrawInsertionMarker(sIdx, lastFoundSIdx, qIdx, lastFoundQIdx, x, y, data);

      /* If this match has been clipped, draw a marker to indicate as such */
      mspDrawClippedMarker(msp, qIdx, segmentIdx, &segmentRange, x, y, data);

      if (sIdx != UNSET_INT)
	{
	  lastFoundSIdx = sIdx;
	  lastFoundQIdx = qIdx;
	}
    }

  /* Null-terminate the string */
  displayText[segmentLen] = '\0';
//  insertChar(displayText, &segmentIdx, '\0', msp);

  /* Draw the sequence text */
  mspDrawSequenceText(tree, displayText, &segmentRange, data);

  g_free(refSeqSegment);
}    


/* Draw a colinearity line between the two given MSPs, if applicable */
static void mspDrawColinearityLine(const MSP* msp1, const MSP* msp2, const gboolean selected, RenderData *data)
{
  if (msp1 && msp2 && data && data->bc)
    {
      /* If the display is reversed we need to swap the order of the msps */
      if (data->bc->displayRev)
        {
          const MSP* tmp = msp1;
          msp1 = msp2;
          msp2 = tmp;
        }

      ColinearityType colinearityType = COLINEAR_INVALID;

      if ((data->bc->displayRev && msp1->qStrand == BLXSTRAND_FORWARD) || (!data->bc->displayRev && msp1->qStrand != BLXSTRAND_FORWARD))
        colinearityType = mspIsColinear(msp2, msp1);
      else
        colinearityType = mspIsColinear(msp1, msp2);

      if (colinearityType != COLINEAR_INVALID)
        {
          /* get the line color */
          GdkColor *color = NULL;

          if (colinearityType == COLINEAR_PERFECT)
            color = getGdkColor(BLXCOLOR_COLINEAR_PERFECT, data->bc->defaultColors, FALSE, data->bc->usePrintColors);
          else if (colinearityType == COLINEAR_IMPERFECT)
            color = getGdkColor(BLXCOLOR_COLINEAR_IMPERFECT, data->bc->defaultColors, FALSE, data->bc->usePrintColors);
          else
            color = getGdkColor(BLXCOLOR_COLINEAR_NOT, data->bc->defaultColors, FALSE, data->bc->usePrintColors);

          gdk_gc_set_foreground(data->gc, color);

          /* Get the coords of the line */
          const int y = data->cell_area->y + (data->charHeight / 2);
          int x1 = data->cell_area->x + ((msp1->displayRange.max + 1 - data->displayRange->min) * data->charWidth); /* +1 to get rightmost edge of char */
          int x2 = data->cell_area->x + ((msp2->displayRange.min - data->displayRange->min) * data->charWidth) - 1; /* -1 offset by 1 pixel so we don't overdraw the char */

          if (x1 < data->cell_area->x + data->cell_area->width && x2 > data->cell_area->x)
            {
              if (x1 < data->cell_area->x)
                x1 = data->cell_area->x;
      
              if (x2 > data->cell_area->x + data->cell_area->width)
                x2 = data->cell_area->x + data->cell_area->width;

              drawLine2(data->window, data->drawable, data->gc, x1, y, x2, y);
            }
        }
    }
}


/* Draw a colinearity line between the given msp and the previous/next one in the BlxSequence's list */
static void mspDrawColinearityLineAdjacent(const MSP* msp, const gboolean selected, RenderData *data, gboolean prev)
{
  if (data && msp && msp->sSequence)
    {
      /* Find the current msp in the list, and get the one before it */
      GList *mspItem = msp->sSequence->mspList;

      for ( ; mspItem; mspItem = mspItem->next)
        {
          if (mspItem->data == msp)
            break;
        }
      
      if (mspItem) /* found it */
        {
          GList *adjacentItem = (prev ? mspItem->prev : mspItem->next);
          
          if (adjacentItem)
            {
              const MSP* adjacentMsp = (const MSP*)(adjacentItem->data);

              if (prev)
                mspDrawColinearityLine(adjacentMsp, msp, selected, data);
              else
                mspDrawColinearityLine(msp, adjacentMsp, selected, data);
            }  
        }
    }
}


/* Draw colinearity lines between the two given, adjacent MSPs in the same row, if applicable */
static void mspDrawColinearityLines(const MSP* cur_msp, const MSP* prev_msp, const gboolean selected, RenderData *data)
{
  g_return_if_fail(data);

  if (data && data->bc && data->bc->flags[BLXFLAG_SHOW_COLINEARITY] &&
      (selected || !data->bc->flags[BLXFLAG_SHOW_COLINEARITY_SELECTED]))
    {
      if (data->bc->seqType == BLXSEQ_DNA)
        {
          /* Draw a line between the adjacent MSPs*/
          mspDrawColinearityLine(prev_msp, cur_msp, selected, data);

          /* If it's the first MSP in the row, then check for adjacent MSPs in the sequence that may be
           * in different rows. */
          if (prev_msp == NULL)
            mspDrawColinearityLineAdjacent(cur_msp, selected, data, TRUE);
          
          /* If it's the last MSP in the row, then check for adjacent MSPs in the sequence that may be
           * in different rows. */
          if (cur_msp == NULL)
            mspDrawColinearityLineAdjacent(prev_msp, selected, data, FALSE);
        }
      else
        {
          /* Protein mode. The two given MSPs are in the same row and hence the same frame, but
           * we therefore can't assume they're adjacent because the real next/previous MSP might
           * be in a different frame. Therefore we have to draw the line to the previous AND next
           * MSP for every MSP (i.e. we draw the joining lines twice; if they are in the same frame
           * they will overlap but otherwise they need to be shown in each frame) */
          if (cur_msp)
            mspDrawColinearityLineAdjacent(cur_msp, selected, data, TRUE);
          
          if (prev_msp)
            mspDrawColinearityLineAdjacent(prev_msp, selected, data, FALSE);
        }

    }
}


/* There can be multiple MSPs in the same cell. This function loops through them
 * and draws each one (IF it is in the correct strand/frame for this tree). */
static void rendererDrawMsps(SequenceCellRenderer *renderer,
                             GtkWidget *tree,
                             GdkWindow *window, 
                             GtkStateType state,
                             GdkRectangle *cell_area)
{
  /* Extract all the info from the tree that we'll need repeatedly. */
  TreeProperties *treeProperties = treeGetProperties(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(treeProperties->detailView);
  BlxViewContext *bc = blxWindowGetContext(detailViewProperties->blxWindow);
  
  const gboolean highlightDiffs = bc->flags[BLXFLAG_HIGHLIGHT_DIFFS]; /* swap match/mismatch colors if this is true */
  const MSP *firstMsp = (const MSP*)(renderer->mspGList->data);
  const BlxSequence *seq = firstMsp ? firstMsp->sSequence : NULL;

  GdkColor *matchColor = getGdkColor(BLXCOLOR_MATCH, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkColor *matchColorSelected = getGdkColor(BLXCOLOR_MATCH, bc->defaultColors, TRUE, bc->usePrintColors);
  GdkColor *mismatchColor = getGdkColor(BLXCOLOR_MISMATCH, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkColor *mismatchColorSelected = getGdkColor(BLXCOLOR_MISMATCH, bc->defaultColors, TRUE, bc->usePrintColors);
  GdkColor *backgroundColorSelected = getGdkColor(BLXCOLOR_BACKGROUND, bc->defaultColors, TRUE, bc->usePrintColors);
  
  GdkGC *gc = gdk_gc_new(window);
  
  /* Make the dashes of the partial-boundary lines very short and closely
   * packed (i.e. dash length of 2 pixels and gaps of 1 pixel) */
  int listLen = 2;
  gint8 dashList[listLen];
  dashList[0] = 2;
  dashList[1] = 1;
  gdk_gc_set_dashes(gc, 1, dashList, listLen);
  
  RenderData data = {
    bc,
    cell_area,
    window,
    state,
    gc,
    treeGetStrand(tree),
    treeProperties->readingFrame,
    detailViewProperties->selectedBaseIdx,
    blxWindowIsSeqSelected(detailViewProperties->blxWindow, seq),
    detailViewProperties->cellXPadding,
    detailViewProperties->cellYPadding,
    detailViewProperties->charWidth,
    detailViewProperties->charHeight,
    &detailViewProperties->displayRange,
    highlightDiffs,
    widgetGetDrawable(tree),
    detailViewProperties->blxWindow,
    getGdkColor(BLXCOLOR_EXON_FILL, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_FILL, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_CDS_FILL, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_CDS_FILL, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_UTR_FILL, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_UTR_FILL, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_PARTIAL_EXON_CROSSHATCH, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_PARTIAL_EXON_CROSSHATCH, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_INSERTION, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_INSERTION, bc->defaultColors, TRUE, bc->usePrintColors),
    highlightDiffs ? mismatchColor : matchColor,
    highlightDiffs ? mismatchColorSelected : matchColorSelected,
    getGdkColor(BLXCOLOR_CONS, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_CONS, bc->defaultColors, TRUE, bc->usePrintColors),
    highlightDiffs ? matchColor : mismatchColor,
    highlightDiffs ? matchColorSelected : mismatchColorSelected,
    getGdkColor(BLXCOLOR_UNALIGNED_SEQ, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_UNALIGNED_SEQ, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_START, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_END, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_POLYA_TAIL, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_POLYA_TAIL, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_CLIP_MARKER, bc->defaultColors, FALSE, bc->usePrintColors),
    detailViewProperties->exonBoundaryLineWidth,
    detailViewProperties->exonBoundaryLineStyle,
    detailViewProperties->exonBoundaryLineStylePartial,
    bc->flags[BLXFLAG_LIMIT_UNALIGNED_BASES],
    detailViewProperties->numUnalignedBases
  };  
  
  /* If a base is selected highlight it now in case we don't come to process it (in 
   * which case this will get drawn over). */
  highlightSelectedBase(data.selectedBaseIdx, backgroundColorSelected, &data);

  /* Draw all MSPs in this row */
  GList *mspListItem = renderer->mspGList;
  MSP *savedMsp = NULL;
  MSP *prevMsp = NULL;
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if (mspIsExon(msp))
        {
          drawExon(renderer, msp, tree, &data);
        }
      else if (mspIsBlastMatch(msp))
        {
          gboolean selected = blxWindowIsSeqSelected(detailViewProperties->blxWindow, msp->sSequence);

          if (mspGetFlag(msp, MSPFLAG_SQUASH_IDENTICAL_FEATURES) && !mspGetFlag(msp, MSPFLAG_SQUASH_LINKED_FEATURES))
            {
              /* The first condition here means that identical matches are placed in the
               * same row and the second means that matches linked in any other way aren't. Hence
               * if we get here we know that if we have multiple MSPs in this row then they
               * must be identical (i.e. duplicate) matches and we only need to draw one of them.
               * This is an important optimisation for BAM data, where we can have hundreds of
               * identical reads on the same row.
               * gb10 2014: Ideally we should add optimisation for the case where we have both 
               * linked features and identical matches in the same row, but this isn't used at the
               * moment. */
              
              /* If any of the MSPs is selected, we want to draw the row as selected, so loop through
               * checking if any are selected. If a selected one is found then draw it; otherwise,
               * save the first MSP so we can go back and draw that. */
              if (selected)
                {
                  data.seqSelected = TRUE;
                  drawMsp(renderer, msp, tree, &data);
                  savedMsp = NULL;
                  break;
                }
              else if (!savedMsp)
                {
                  savedMsp = msp;
                }
            }
          else
            {
              /* Ordinary row: draw all MSPs */
              drawMsp(renderer, msp, tree, &data);
              mspDrawColinearityLines(msp, prevMsp, selected, &data);
            }
        }
      
      prevMsp = msp;
    }
  
  /* Draw colinearity lines for the last msp. Passing curMsp as null indicates it's the last msp
   * so we search for adjacent msps in different rows. */
  if (prevMsp)
    mspDrawColinearityLines(NULL, prevMsp, blxWindowIsSeqSelected(detailViewProperties->blxWindow, prevMsp->sSequence), &data);

  if (savedMsp)
    drawMsp(renderer, savedMsp, tree, &data);
  
  drawAllVisibleExonBoundaries(tree, &data);
  
  g_object_unref(gc);
}


static GtkStateType getState(GtkWidget *widget, GtkCellRendererState flags)
{
  GtkStateType state;
  if ((flags & GTK_CELL_RENDERER_SELECTED) == GTK_CELL_RENDERER_SELECTED)
    {
      state = GTK_WIDGET_HAS_FOCUS(widget) ? GTK_STATE_SELECTED : GTK_STATE_ACTIVE;
    }
  else if ((flags & GTK_CELL_RENDERER_PRELIT) == GTK_CELL_RENDERER_PRELIT && GTK_WIDGET_STATE (widget) == GTK_STATE_PRELIGHT)
    {
      state = GTK_STATE_PRELIGHT;
    }
  else if (GTK_WIDGET_STATE (widget) == GTK_STATE_INSENSITIVE)
    {
      state = GTK_STATE_INSENSITIVE;
    }
  else
    {
      state = GTK_STATE_NORMAL;
    }
  
  return state;
}


/* Utility function that returns true if any of the MSPs in the given list
 * is selected. */
static gboolean mspListContainsSelectedMsp(GList *mspList, const BlxViewContext* const bc)
{
  gboolean isSelected = FALSE;
  GList *mspItem = mspList;
  
  for ( ; mspItem && !isSelected; mspItem = mspItem->next)
    {
      const MSP* const msp = (const MSP*)(mspItem->data);
      isSelected = blxContextIsSeqSelected(bc, msp->sSequence);
    }
  
  return isSelected;
}


/* Utility to determine if any MSP in the given list is in a group and, if so
 * to return that group.  Returns the first group found and ignores any 
 * subsequent MSPs in the list that also have groups. Returns null if no group
 * was found. */
static SequenceGroup* mspListContainsGroupedMsp(GList *mspList, const BlxViewContext* const bc)
{
  SequenceGroup *group = NULL;
  GList *mspItem = mspList;
  
  for ( ; mspItem && !group; mspItem = mspItem->next)
    {
      const MSP* const msp = (const MSP*)(mspItem->data);
      group = blxContextGetSequenceGroup(bc, msp->sSequence);
    }
  
  return group;
  
}


/* This function sets the background color for the row. */
static void rendererSetBackgroundColor(GtkCellRenderer *cell, GtkWidget *tree, GdkWindow *window, GdkRectangle *background_area)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  if (renderer->data)
    {
      /* Find out whether the MSP(s) that this cell is displaying are in 
       * a grouped sequence or are selected. */
      const BlxViewContext* const bc = blxWindowGetContext(treeGetBlxWindow(tree));
      GList *mspList = renderer->data;
      
      const gboolean isSelected = mspListContainsSelectedMsp(mspList, bc);
      const SequenceGroup* const group = mspListContainsGroupedMsp(mspList, bc);
      
      GdkGC *gc = gdk_gc_new(window);

      if (isSelected || (group && group->highlighted))
	{
	  if (group && group->highlighted && isSelected)
	    {
	      /* Use the group's highlight color but darken it because it is also selected */
	      GdkColor color;
	      getSelectionColor(&group->highlightColor, &color);
	      gdk_gc_set_foreground(gc, &color);
	    }
	  else if (group && group->highlighted)
	    {
	      /* Use the group's highlight color */
	      gdk_gc_set_foreground(gc, &group->highlightColor);
	    }
	  else
	    {
	      /* Not in a group but is selected. Use the background color but darken it. */
	      GdkColor color;
	      getSelectionColor(&tree->style->bg[GTK_STATE_NORMAL], &color);
	      gdk_gc_set_foreground(gc, &color);
	    }
	}
      else
        {
  	  gdk_gc_set_foreground(gc, &tree->style->bg[GTK_STATE_NORMAL]);
        }

      drawRectangle2(window, widgetGetDrawable(tree), gc, TRUE, background_area->x, background_area->y, background_area->width, background_area->height);
      
      g_object_unref(gc);
    }
}


/* Draw the tree grid lines (i.e. column separators) */
static void rendererDrawGridLines(GtkWidget *tree, GdkWindow *window, GdkRectangle *cell_area)
{
  GdkDrawable *drawable = widgetGetDrawable(tree);
  BlxViewContext *bc = treeGetContext(tree);

  GdkGC *gc = gdk_gc_new(window);

  /* Set the line color */
  GdkColor *color = getGdkColor(BLXCOLOR_TREE_GRID_LINES, bc->defaultColors, FALSE, bc->usePrintColors);
  gdk_gc_set_foreground(gc, color);
  
  /* We want dashed lines */
  const int lineWidth = 1;
  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_MITER);

  /* Make the dashes very short and closely packed (i.e. dash length of 1 pixel and gaps of 1 pixel) */
  int listLen = 1;
  gint8 dashList[listLen];
  dashList[0] = 1;
  gdk_gc_set_dashes(gc, 0, dashList, listLen);
  
  /* Draw vertical lines. Draw the right edge of each cell (because we don't really want a
   * line at the leftmost edge of the first cell.) */
  const int x = cell_area->x + cell_area->width - lineWidth;
  drawLine2(window, drawable, gc, x, cell_area->y, x, cell_area->y + cell_area->height);
  
  g_object_unref(gc);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_get_size: crucial - calculate the size
 *                                          of our cell, taking into account
 *                                          padding and alignment properties
 *                                          of parent.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_get_size (GtkCellRenderer *cell,
                                        GtkWidget       *widget,
                                        GdkRectangle    *cell_area,
                                        gint            *x_offset,
                                        gint            *y_offset,
                                        gint            *width,
                                        gint            *height)
{
  get_size(cell, widget, cell_area, x_offset, y_offset, width, height);  
}


/***************************************************************************
 *
 *  sequence_cell_renderer_render: crucial - do the rendering.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_render (GtkCellRenderer *cell,
			       GdkDrawable     *window,
			       GtkWidget       *tree,
			       GdkRectangle    *background_area,
			       GdkRectangle    *cell_area,
			       GdkRectangle    *expose_area,
			       GtkCellRendererState flags)
{
  rendererSetBackgroundColor(cell, tree, window, background_area);
  
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  MSP *msp = renderer->mspGList ? (MSP*)(renderer->mspGList->data) : NULL;
  
  if (msp)
    {
      rendererDrawMsps(renderer, tree, window, getState(tree, flags), cell_area);
    }
  else
    {
      rendererDrawSimpleText(renderer, tree, window, getState(tree, flags), cell_area);
    }
  
  rendererDrawGridLines(tree, window, background_area);
}



