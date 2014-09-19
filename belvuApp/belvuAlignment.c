/*  File: belvuAlignment.c
 *  Author: Gemma Barson, 2011-04-12
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
 * Description: see belvuAlignment.h
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvuAlignment.h"
#include "belvuApp/belvuWindow.h"
#include <math.h>


#define DEFAULT_XPAD                            2
#define DEFAULT_YPAD                            4
#define DEFAULT_PADDING_CHARS                   1  /* number of char widths to use to pad between columns */
#define SEQ_AREA_PADDING                        2  /* padding around the sequence area, in pixels */
#define DEFAULT_NAME_COLUMN_PADDING_CHARS       2  /* number of char widths to use to pad after the name column */
#define WRAP_DISPLAY_PADDING_CHARS              4
#define TICKMARK_INTERVAL                       10 /* number of coords between each tick marker in the sequence area header */
#define MAJOR_TICKMARK_HEIGHT                   6  /* height of major tick marks in the sequence area header */
#define MINOR_TICKMARK_HEIGHT                   3  /* height of minor tick marks in the sequence area header */

/* Local function declarations */
static void               bg2fgColor(BelvuContext *bc, GdkColor *bgColor, GdkColor *result);


/* Properties specific to the belvu alignment */
typedef struct _BelvuAlignmentProperties
{
  BelvuContext *bc;                 /* The belvu context */
  
  GtkWidget *columnsArea;           /* Drawing widget for the columns (i.e. name and coord columns) */
  GtkWidget *columnsHeader;         /* Drawing widget for the header for the the columns area */
  GtkWidget *seqArea;               /* Drawing widget for the actual sequence */
  GtkWidget *seqHeader;             /* Drawing widget for the header for the sequence area */
  GdkRectangle columnsRect;         /* Drawing area for the columns (i.e. name and coord columns) */
  GdkRectangle seqRect;             /* Drawing area for the actual sequence */
  GdkRectangle seqHeaderRect;       /* Drawing area for the sequence header */
  GtkAdjustment *hAdjustment;       /* Controls horizontal scroll bar */
  GtkAdjustment *vAdjustment;       /* Controls vertical scroll bar */
  
  int columnPadding;                /* Padding in between columns, in pixels */
  int nameColumnPadding;            /* Padding after name column, in pixels */
  char *title;                /* title to display at the top of the alignments */
  int wrapWidth;                    /* Number of characters after which to wrap (or UNSET_INT for no wrapping) */
  
  gdouble charWidth;                /* The width of each character in the display */
  gdouble charHeight;               /* The height of each character in the display */
} BelvuAlignmentProperties;


/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuAlignmentProperties* belvuAlignmentGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuAlignmentProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuAlignmentProperties")) : NULL;
}

static void onDestroyBelvuAlignment(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  if (properties)
    {
      if (properties->title)
        {
          g_free(properties->title);
          properties->title = NULL;
        }
      
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void belvuAlignmentCreateProperties(GtkWidget *belvuAlignment, 
                                           BelvuContext *bc,
                                           GtkWidget *columnsArea,
                                           GtkWidget *columnsHeader,
                                           GtkWidget *seqArea,
                                           GtkWidget *seqHeader,
                                           GtkAdjustment *hAdjustment,
                                           GtkAdjustment *vAdjustment,
                                           const char *title,
                                           const int wrapWidth)
{
  if (belvuAlignment)
    {
      BelvuAlignmentProperties *properties = (BelvuAlignmentProperties*)g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->columnsArea = columnsArea;
      properties->columnsHeader = columnsHeader;
      properties->seqArea = seqArea;
      properties->seqHeader = seqHeader;
      properties->columnPadding = 0; /* calculated in calculate-borders */
      properties->hAdjustment = hAdjustment;
      properties->vAdjustment = vAdjustment;
      properties->title = g_strdup(title);
      properties->wrapWidth = wrapWidth;
      
      properties->seqRect.x = 0;
      properties->seqRect.y = DEFAULT_YPAD;
      properties->seqHeaderRect.x = properties->seqRect.x;
      properties->seqHeaderRect.y = DEFAULT_YPAD;
      properties->columnsRect.x = DEFAULT_XPAD;
      properties->columnsRect.y = DEFAULT_YPAD;

      /* Find a fixed-width font */
      const char *fontFamily = findFixedWidthFont(properties->seqArea);
      PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
      gtk_widget_modify_font(seqArea, fontDesc);

      if (columnsArea)
        gtk_widget_modify_font(columnsArea, fontDesc);

      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", properties);
      g_signal_connect(G_OBJECT(belvuAlignment), "destroy", G_CALLBACK (onDestroyBelvuAlignment), NULL);
    }
}

/***********************************************************
 *                         Misc                         *
 ***********************************************************/

/* This should be called after the font has changed. It caches
 * the new font size in the properties and recalculates all sizes */
void onBelvuAlignmentFontSizeChanged(GtkWidget *belvuAlignment)
{
  if (!belvuAlignment)
    return;

  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (!properties || !properties->seqArea)
    return;
  
  getFontCharSize(properties->seqArea, properties->seqArea->style->font_desc, &properties->charWidth, &properties->charHeight);

  /* Update the column padding, which is based on the character width. We use
   * slightly different padding if the view is wrapped because we don't draw 
   * column separators */
  if (properties->wrapWidth != UNSET_INT) /* wrapped */
    {
      properties->columnPadding = (DEFAULT_PADDING_CHARS * properties->charWidth);
      properties->nameColumnPadding = (DEFAULT_NAME_COLUMN_PADDING_CHARS * properties->charWidth);
    }
  else /* unwrapped */
    {
      properties->columnPadding = (DEFAULT_PADDING_CHARS * properties->charWidth) / 2;
      properties->nameColumnPadding = 0;
    }
  
  
  /* The font size affects the height of the header line */
  properties->seqHeaderRect.height = properties->charHeight + DEFAULT_YPAD;

  if (properties->seqHeader)
    gtk_widget_set_size_request(properties->seqHeader, -1, properties->seqHeaderRect.y + properties->seqHeaderRect.height);

  /* Set the size of the header columns widget */
  updateHeaderColumnsSize(belvuAlignment);

  /* Recalculate the size of the drawing areas */
  calculateDrawingSizes(belvuAlignment);

  belvuAlignmentRedrawAll(belvuAlignment);
}



/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Refresh all */
void belvuAlignmentRefreshAll(GtkWidget *belvuAlignment)
{
  gtk_widget_queue_draw(belvuAlignment);
}


/* Clear cached drawables and redraw all */
void belvuAlignmentRedrawAll(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  if (properties->seqArea)
    widgetClearCachedDrawable(properties->seqArea, NULL);
  
  if (properties->seqHeader)
    widgetClearCachedDrawable(properties->seqHeader, NULL);
  
  if (properties->columnsArea)
    widgetClearCachedDrawable(properties->columnsArea, NULL);

  if (properties->columnsHeader)
    widgetClearCachedDrawable(properties->columnsHeader, NULL);
  
  gtk_widget_queue_draw(belvuAlignment);
}


/* Clear cached drawables for the sequence area only */
static void belvuAlignmentRedrawSequenceArea(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->seqArea)
    widgetClearCachedDrawable(properties->seqArea, NULL);
  
  if (properties->seqHeader)
    widgetClearCachedDrawable(properties->seqHeader, NULL);
  
  gtk_widget_queue_draw(belvuAlignment);
}


/* Find the background color of the given residue index i in the given alignment */
static void findResidueBGcolor(BelvuContext *bc, ALN* alnp, int i, const gboolean isSelected, GdkColor *result) 
{
  int colorNum = 0;
  char *alnpSeq = alnGetSeq(alnp);
  
  if (!alnpSeq || i >= alnGetSeqLen(alnp))
    colorNum = WHITE;
  else if (bc->lowercaseOn && alnpSeq[i] >= 'a' && alnpSeq[i] <= 'z') 
    colorNum = ORANGE;
  else if (alnp->nocolor)
    colorNum = BOXCOLOR;
  else if (alnp->markup)
    colorNum = getMarkupColor(alnpSeq[i]);
  else if (colorByConservation(bc) || colorByResId(bc))
    colorNum = getConservColor(bc, alnpSeq[i], i);
  else
    colorNum = getColor(alnpSeq[i]);

  convertColorNumToGdkColor(colorNum, isSelected, result);
}


/* Draw a single row in the columns area */
static void drawSingleColumn(GtkWidget *widget,
                             GdkDrawable *drawable, 
                             BelvuAlignmentProperties *properties,
                             ALN *alnp, 
                             const int lineNum)
{
  int x = 0;
  const int y = properties->columnsRect.y + (properties->charHeight * lineNum);
  
  const gboolean highlightAln = alignmentHighlighted(properties->bc, alnp);
  GdkGC *gc = gdk_gc_new(drawable);
  
  if (highlightAln || alnp->color != WHITE)
    {
      /* Highlight the background */
      GdkColor bgColor;
      convertColorNumToGdkColor(alnp->color, highlightAln, &bgColor);
      
      gdk_gc_set_foreground(gc, &bgColor);
      gdk_draw_rectangle(drawable, gc, TRUE, properties->columnsRect.x, y, properties->columnsRect.width, properties->charHeight);
    }
  
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, highlightAln, FALSE);
  gdk_gc_set_foreground(gc, textColor);

  /* Draw the name */
  if (alnp->name)
    {
      x += properties->columnPadding;
      drawText(widget, drawable, gc, x, y, alnp->name, NULL, NULL);
      x += (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding;
    }
  
  /* Draw the coords (but only if this is not GC markup; note that we always
   * draw the column separator lines for the coords, though) */
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  x += properties->columnPadding + properties->bc->maxStartLen * properties->charWidth; /* right-hand edge of column for right-aligned text */

  if (alnp->markup != GC)
    drawIntAsText(widget, drawable, gc, x, y, alnp->start);
  
  x += properties->columnPadding;
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  
  x += properties->columnPadding + properties->bc->maxEndLen * properties->charWidth; /* right-aligned */

  if (alnp->markup != GC)
    drawIntAsText(widget, drawable, gc, x, y, alnp->end);

  x += properties->columnPadding;
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  
  /* Draw the score, if displaying scores (and if not a markup row) */
  if (properties->bc->displayScores)
    {
      x += (properties->bc->maxScoreLen * properties->charWidth) + properties->columnPadding;
      
      if (!alnp->markup)
        drawDoubleAsText(widget, drawable, gc, x, y, alnp->score);
      
      x += properties->columnPadding;
      gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
    }
  
  g_object_unref(gc);
}


/* Utility to return true if the given alignment is selected */
static gboolean alignmentSelected(BelvuContext *bc, ALN *alnp)
{
  /* Return true if this alignment has the same name as the selected alignment */
  return (bc->selectedAln && bc->selectedAln == alnp);
}


/* Draw a single character in the given sequence. In standard mode, fgColorsOnly
 * should be false and the function will draw the background color for the 
 * base (and will not draw any text). This is because drawing each character
 * separately is slow, so all text in the default color is drawn at once
 * by the calling function.
 * However, if fgColorsOnly is true then this function checks if the required
 * foreground color is different to the default and, if so, draws the character
 * in that color. (It does not draw background colors in this mode.) */
static void drawSequenceChar(BelvuAlignmentProperties *properties, 
                             ALN *alnp,
                             const int colIdx,
                             const gboolean rowHighlighted,
                             GtkWidget *widget,
                             GdkDrawable *drawable,
                             GdkGC *gc,
                             GdkColor *defaultFgColor,
                             const gboolean fgColorsOnly,
                             const int x, 
                             const int y)
{
  if (!fgColorsOnly)
    {
      /* Draw the background */
      GdkColor bgColor;
      findResidueBGcolor(properties->bc, alnp, colIdx, rowHighlighted, &bgColor);
      gdk_gc_set_foreground(gc, &bgColor);
      
      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
    }
  else if (colIdx < alnGetSeqLen(alnp) && colorByConservation(properties->bc))
    {
      /* Draw the text if it is in a different color to the default. This is only
       * possible in color by conservation mode. We get the text colour from the 
       * background color */
      GdkColor fgColor, bgColor;
      findResidueBGcolor(properties->bc, alnp, colIdx, rowHighlighted, &bgColor);
      bg2fgColor(properties->bc, &bgColor, &fgColor);
      
      if (!colorsEqual(&fgColor, defaultFgColor))
        {
          gdk_gc_set_foreground(gc, &fgColor);
          
          char displayText[2];
          displayText[0] = alnGetSeq(alnp)[colIdx];
          displayText[1] = '\0';
          drawText(widget, drawable, gc, x, y, displayText, NULL, NULL);
        }
    }
  
}


/* Draw a single line in the sequence area.
 * 
 * This is rather convoluted in an attempt to maximise speed.  Drawing of text in
 * GTK is particularly slow, and is best done a whole line at a time, rather than
 * drawing individual characters.  The characters may each have a different background
 * colour, so we draw the background colours first and then draw the text over the top
 * a whole line at a time.  A further problem, though, is that some characters may
 * have a different foreground colour to others on the same line (in colour-by-conservation
 * mode), so some characters do need to be drawn individually.  There are usually
 * relatively few characters that this affects, though, so we continue to draw the
 * whole line in the default colour first, and then we look to see if any individual
 * characters have a different foreground colour and just draw those over the top.
*/
static void drawSingleSequence(GtkWidget *widget,
                                GdkDrawable *drawable, 
                                BelvuAlignmentProperties *properties,
                                ALN *alnp, 
                                const int lineNum)
{
  if (!alnGetSeq(alnp))
    return;
  
  const int startX = properties->seqRect.x;
  int x = startX;
  const int y = properties->seqRect.y + (properties->charHeight * lineNum);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkAdjustment *hAdjustment = properties->hAdjustment;
  
  GdkColor *defaultFgColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, FALSE, FALSE);
  
  /* Loop through each column in the current display range and color
   * the according to the relevant highlight color */
  int colIdx = hAdjustment->value;
  int iMax = min(properties->bc->maxLen, hAdjustment->value + hAdjustment->page_size);
  const gboolean rowHighlighted = alignmentSelected(properties->bc, alnp);
  
  if (properties->bc->displayColors)
    {
      for ( ; colIdx < iMax; ++colIdx)
        {
          drawSequenceChar(properties, alnp, colIdx, rowHighlighted, widget, drawable, gc, defaultFgColor, FALSE, x, y);
          x += properties->charWidth;
        }
    }  
  else if (rowHighlighted)
    {
      /* We're not displaying colors for individual chars, but we still need to
       * highlight the background if the row is selected. */
      GdkColor color;
      convertColorNumToGdkColor(WHITE, TRUE, &color);
      gdk_gc_set_foreground(gc, &color);
      
      gdk_draw_rectangle(drawable, gc, TRUE, startX, y, properties->seqRect.width, properties->charHeight);
    }
  
  
  /* Draw the sequence text (current display range only) */
  if (alnGetSeq(alnp))
    {
      GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, FALSE, FALSE);
      gdk_gc_set_foreground(gc, textColor);

      /* Start at the number of characters into the string where the horizontal 
       * scrollbar indicates we are (making sure that's not out of the end of the
       * string) */
      if ((int)hAdjustment->value < alnGetSeqLen(alnp))
        {
          char *cp = alnGetSeq(alnp) + (int)hAdjustment->value;
          char *displayText = g_strndup(cp, iMax - properties->hAdjustment->value);
      
          drawText(widget, drawable, gc, startX, y, displayText, NULL, NULL);
      
          g_free(displayText);
        }
    }
  
  /* Loop again and draw any characters that are not in the default text color.
   * We don't draw every character individually because the pango layout calls
   * are expensive, but we can't really avoid it when some characters are in
   * difference colors.  (We could put markup in the pango layout instead, but
   * generally the number of characters in a non-default color is small, so 
   * this should be good enough for most cases.) */
  if (properties->bc->displayColors)
    {
      for (x = startX, colIdx = hAdjustment->value; colIdx < iMax; ++colIdx)
        {
          drawSequenceChar(properties, alnp, colIdx, rowHighlighted, widget, drawable, gc, defaultFgColor, TRUE, x, y);
          x += properties->charWidth;
        }
    } 
  
  g_object_unref(gc);
}


/* Return the foreground color selected for a given
   background color in color-by-conservation mode 
 */
static void bg2fgColor(BelvuContext *bc, GdkColor *bgColor, GdkColor *result)
{
  int fgColorNum = BLACK; /* default colour for uncoloured or markup */
  
  /* See if the given bgcolor matches the max-conservation background colour */
  int testColor = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
  convertColorNumToGdkColor(testColor, FALSE, result);
  
  if (colorsEqual(bgColor, result))
    {
      fgColorNum = *getConsColor(bc, CONS_LEVEL_MAX, TRUE);
    }
  else
    {
      /* Try the mid-conservation colour */
      int testColor = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
      convertColorNumToGdkColor(testColor, FALSE, result);
      
      if (colorsEqual(bgColor, result))
        {
          fgColorNum = *getConsColor(bc, CONS_LEVEL_MID, TRUE);
        }
      else
        {
          /* Lastly, try the low-conservation colour */
          int testColor = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
          convertColorNumToGdkColor(testColor, FALSE, result);
          
          if (colorsEqual(bgColor, result))
            fgColorNum = *getConsColor(bc, CONS_LEVEL_LOW, TRUE);
        }
    }

  /* Convert the foreground colour number into the GdkColor result. */
  convertColorNumToGdkColor(fgColorNum, FALSE, result);
}


static void drawWrappedSequences(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  /*  Column makeup:
   space, maxNameLen, 2xspace, maxEndLen, space, maxLen, maxEndLen
   i.e. 4 spaces
   */
  
  int i, oldpos, collapseRes = 0, collapsePos = 0, collapseOn = 0;
  int numSpaces = WRAP_DISPLAY_PADDING_CHARS;
  char collapseStr[10],
  ch[2] = " ";
  static int *pos=0;                    /* Current residue position of sequence j */
  
  GdkColor bgColor;
  GdkColor *pBgColor = &bgColor;
  GdkGC *gcText = gdk_gc_new(drawable);
  GdkGC *gc = gdk_gc_new(drawable);
  
  BelvuContext *bc = properties->bc;
  
  
  /* Initialise the sequence and position array */
  char *seq = (char*)g_malloc(properties->wrapWidth + 1);
  
  if (!pos) 
    pos = (int *)g_malloc(bc->alignArr->len * sizeof(int *));
  
  int j = 0;
  for (j = 0; j < (int)bc->alignArr->len; ++j) 
    pos[j] = g_array_index(bc->alignArr, ALN*, j)->start;
  
  
  int paragraph = 0;
  int totCollapsed = 0;
  int line = 1; /* leave a blank line (line=0) at the top for spacing */
  
  if (properties->title)
    {
      const int x = properties->columnPadding;
      const int y = line * properties->charHeight;
      
      drawText(widget, drawable, gcText, x, y, properties->title, NULL, NULL);
      line += 2; /* increment, and also add another blank line for spacing */
    }
  
  gboolean quit = FALSE;
  while (!quit && paragraph * properties->wrapWidth + totCollapsed < bc->maxLen)
    {
      gboolean emptyPara = TRUE;
      
      for (j = 0; j < (int)bc->alignArr->len; ++j)
        {
          ALN *alnp = g_array_index(bc->alignArr, ALN*, j);
          char *alnpSeq = alnGetSeq(alnp);
          
          if (alnp->hide) 
            continue;
          
          int alnstart = paragraph*properties->wrapWidth +totCollapsed; 
          int alnlen = ( (paragraph+1)*properties->wrapWidth +totCollapsed < bc->maxLen ? properties->wrapWidth : bc->maxLen - alnstart );
          int alnend = min(alnstart + alnlen, alnGetSeqLen(alnp));
          
          gboolean emptyLine = TRUE;
          for (i = alnstart; i < alnend; ++i) 
            {
              if (alnpSeq && !isGap(alnpSeq[i]) && alnpSeq[i] != ' ') 
                {
                  emptyLine = FALSE;
                  emptyPara = FALSE;
                  break;
                }
            }
          
          if (!emptyLine) 
            {
              const int y = line * properties->charHeight;
              
              if (y > properties->seqRect.y + properties->seqRect.height)
                {
                  quit = TRUE;
                  break;
                }

              for (collapsePos = 0, oldpos = pos[j], i = alnstart; i < alnend; i++) 
                {       
                  const int xpos = bc->maxNameLen + bc->maxEndLen + numSpaces + i - alnstart - collapsePos;
                  const int x = xpos * properties->charWidth;
                  
                  if (alnpSeq[i] == '[') 
                    {
                      /* Experimental - collapse block: "[" -> Collapse start, "]" -> Collapse end, e.g.
                       1 S[..........]FKSJFE
                       2 L[RVLKLKYKNS]KDEJHF
                       3 P[RL......DK]FKEJFJ
                                 |
                                 V
                       1 S[  0]FKSJFE
                       2 L[ 10]KDEJHF
                       3 P[  4]FKEJFJ
                       
                       Minimum collapse = 5 chars. The system depends on absolute coherence in format; 
                       very strange results will be generated otherwise. Edges need to be avoided manually.
                       */
                      collapseOn = 1;
                      collapsePos += 1;
                      collapseRes = 0;
                      continue;
                    }
                  
                  if (alnpSeq[i] == ']') 
                    {
                      collapseOn = 0;
                      collapsePos -= 4;
                      alnend -= 3;
                      
                      sprintf(collapseStr, "[%3d]", collapseRes);
                      drawText(widget, drawable, gcText, x, y, collapseStr, NULL, NULL);

                      continue;
                    }
                  
                  if (collapseOn) 
                    {
                      collapsePos++;
                      if (!isGap(alnpSeq[i])) 
                        {
                          collapseRes++;
                          pos[j]++;
                        }
                      alnend++;
                      continue;
                    }
                  
                  
                  if (!isGap(alnpSeq[i]) && alnpSeq[i] != ' ') 
                    {
                      findResidueBGcolor(bc, alnp, i, FALSE, pBgColor);
                      gdk_gc_set_foreground(gc, pBgColor);
                      
                      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
                      
                      pos[j]++;
                    }
                  
                  /* Foreground color */
                  GdkColor fgColor;
                  if (colorByConservation(bc))
                    bg2fgColor(bc, pBgColor, &fgColor);
                  else
                    convertColorNumToGdkColor(BLACK, FALSE, &fgColor);
                  
                  gdk_gc_set_foreground(gc, &fgColor);

                  *ch = alnpSeq[i];
                  drawText(widget, drawable, gc, x, y, ch, NULL, NULL);
                }
          
              drawText(widget, drawable, gcText, properties->charWidth, y, alnp->name, NULL, NULL);
              
              if (!alnp->markup) 
                {
                  char *tmpStr = g_strdup_printf("%*d", bc->maxEndLen, oldpos);
                  drawText(widget, drawable, gcText, (bc->maxNameLen + 3) * properties->charWidth, y, tmpStr, NULL, NULL);
                  g_free(tmpStr);
                  
                  if (alnend == bc->maxLen) 
                    {
                      char *tmpStr = g_strdup_printf("%-d", pos[j] - 1);
                      drawText(widget, drawable, gcText, (bc->maxNameLen + bc->maxEndLen + alnlen + 5) * properties->charWidth, y, tmpStr, NULL, NULL);
                      g_free(tmpStr);
                    }
                }
              
              line++;
            }
        }

      /* Move to next paragraph */
      paragraph++;
      
      /* Add a blank line (unless nothing was drawn for this paragraph) */
      if (!emptyPara)
        line++;
      
      totCollapsed += collapsePos;
    }
  
  g_free(seq);
  g_object_unref(gc);
  g_object_unref(gcText);

  /* Set the layout size now we know how big it is */
  properties->seqRect.height = line * properties->charHeight;
  
  gtk_layout_set_size(GTK_LAYOUT(properties->seqArea),
                      properties->seqRect.x + properties->seqRect.width,
                      properties->seqRect.y + properties->seqRect.height);
  
  if (properties->hAdjustment)
    gtk_adjustment_changed(properties->hAdjustment); /* signal that the scroll range has changed */
  
  if (properties->vAdjustment)
    gtk_adjustment_changed(properties->vAdjustment);
}


/* Draw the header line for the columns section */
static void drawBelvuColumnsHeader(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  GdkGC *gc = gdk_gc_new(drawable);

  /* Draw the summary line, which shows the number of sequences & alignment length */
  int x = properties->columnPadding;
  int y = DEFAULT_YPAD;
  
  char *tmpStr = g_strdup_printf("(%dx%d)", bc->alignArr->len, bc->maxLen);
  
  drawText(widget, drawable, gc, x, y, tmpStr, NULL, NULL);
  y += properties->charHeight + DEFAULT_YPAD - 1;
  
  g_free(tmpStr);
  
  /* Draw a horizontal separator line */
  x = 0;
  gdk_draw_line(drawable, gc, x, y, x + properties->columnsRect.width, y);
  
  g_object_unref(gc);
}


/* Draw the alignment view */
static void drawBelvuColumns(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;

  /* Draw each visible alignment */
  GtkAdjustment *vAdjustment = properties->vAdjustment;
  
  int i = vAdjustment->value;
  int lineNum = 0;
  const int iMax = min(bc->alignArr->len, vAdjustment->value + vAdjustment->page_size);
  
  for ( ; i < iMax; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (alnp && !alnp->hide)
        {
          drawSingleColumn(widget, drawable, properties, alnp, lineNum);
          ++lineNum;
        }
    }
}


/* Draw the header line in the sequence area */
static void drawBelvuSequenceHeader(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  /* The range starts at the current adjustment value (add one to make it a 1-based column number) */
  IntRange displayRange = {properties->hAdjustment->value + 1,
                           min(properties->bc->maxLen, properties->hAdjustment->value + properties->hAdjustment->page_size)};
  
  drawHScale(widget, 
             drawable,
             &displayRange,
             &properties->seqHeaderRect,
             properties->charWidth, /* width between each minor tickmark */
             TICKMARK_INTERVAL / 2, /* major markers half way between labels */
             TICKMARK_INTERVAL,     /* interval between main tickmarks with labels */
             MINOR_TICKMARK_HEIGHT,
             MAJOR_TICKMARK_HEIGHT);
}


/* This draws some shading where the currently-selected column is in order to
 * highlight it.  If there is no column selected, does nothing. */
static void drawSelectedColumn(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  /* Quit if there is no column selected */
  if (bc->highlightedCol < 1)
    return;

  GtkAdjustment *hAdjustment = properties->hAdjustment;
  
  const int iMin = properties->hAdjustment->value + 1;
  const int iMax = min(bc->maxLen, hAdjustment->value + hAdjustment->page_size + 1);

  /* Quit if the selected column is not in the visible range */
  if (bc->highlightedCol < iMin || bc->highlightedCol > iMax)
    return;

  /* Get the column number as a zero-based index in the current range and 
   * convert this to a distance from the left edge, which will be our x coord.
   * Add half a char's width to get the x coord at the centre of the char. */
  const int colIdx = bc->highlightedCol - iMin;
  const int x = properties->seqRect.x + (colIdx * properties->charWidth) + (properties->charWidth / 2);
  

  /* Draw a vertical line at this index */
  GdkGC *gc = gdk_gc_new(drawable);

  GdkColor *color = getGdkColor(BELCOLOR_COLUMN_HIGHLIGHT, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, color);
  gdk_gc_set_function(gc, GDK_INVERT);
  gdk_draw_line(drawable, gc, x, 0, x, properties->seqRect.height);
  
  g_object_unref(gc);
  
//  /* Draw a transparent box at this index */
//  cairo_t *cr = gdk_cairo_create(drawable);
//  gdk_cairo_set_source_color(cr, color);
//  cairo_rectangle(cr, x, 0, properties->charWidth, properties->seqRect.height);
//  cairo_clip(cr);
//  cairo_paint_with_alpha(cr, 0.2);
//  cairo_destroy(cr);
}


/* Draw the alignment view */
static void drawBelvuSequence(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  /* Draw each visible alignment */
  GtkAdjustment *vAdjustment = properties->vAdjustment;
  
  int i = vAdjustment->value;
  int lineNum = 0;
  const int iMax = min(bc->alignArr->len, vAdjustment->value + vAdjustment->page_size);
  
  for ( ; i < iMax; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (alnp && !alnp->hide)
        {
          drawSingleSequence(widget, drawable, properties, alnp, lineNum);
          ++lineNum;
        }
    }
}


/* Expose handler for the alignment section */
static gboolean onExposeBelvuColumns(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->columnsRect.width, properties->columnsRect.height);
          drawBelvuColumns(widget, bitmap, properties);
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
          g_warning("Failed to draw Belvu alignment [%p] - could not create bitmap.\n", widget);
        }
    }

  return TRUE;
}

/* Expose handler for the alignment section */
static gboolean onExposeBelvuSequence(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->seqRect.width, properties->seqRect.height);
          
          if (properties->wrapWidth == UNSET_INT)
            drawBelvuSequence(widget, bitmap, properties);
          else
            drawWrappedSequences(widget, bitmap, properties);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
          
          /* Draw the currently-selected column on top */
          drawSelectedColumn(widget, window, properties);
        }
      else
        {
          g_warning("Failed to draw Belvu alignment [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}


/* Expose handler for the header widget of the alignment section */
static gboolean onExposeBelvuSequenceHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  GdkDrawable *window = widget->window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->seqRect.width, properties->seqRect.height);
          drawBelvuSequenceHeader(widget, bitmap, properties);
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
          g_warning("Failed to draw Belvu alignment header [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}


/* Expose handler for the header widget of the columns section */
static gboolean onExposeBelvuColumnsHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  GdkDrawable *window = widget->window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->seqRect.width, properties->seqRect.height);
          drawBelvuColumnsHeader(widget, bitmap, properties);
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
          g_warning("Failed to draw Belvu columns header [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}

/***********************************************************
 *                         Scrolling                       *
 ***********************************************************/

/* Called when the horizontal scroll position has changed */
static void onHScrollPosChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->wrapWidth == UNSET_INT)
    {
      /* We need to clear and re-draw the sequence area (unless we're displaying
       * wrapped sequences, in which case we've already drawn the full sequence
       * width in the drawing area */
      belvuAlignmentRedrawSequenceArea(belvuAlignment);
    }
}

/* Called when the vertical scroll position has changed */
static void onVScrollPosChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->wrapWidth == UNSET_INT)
    {
      belvuAlignmentRedrawAll(belvuAlignment);
    }
}


/* Called when the horizontal scroll range has changed */
static void onHScrollRangeChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  /* Set the horizontal adjustment page size to be the number of characters
   * that will fit in the width */
  properties->hAdjustment->page_size = (int)(properties->seqRect.width / properties->charWidth);
  properties->hAdjustment->page_increment = properties->hAdjustment->page_size;

  /* Make sure the scroll position still lies within upper - page_size. */
  if (properties->hAdjustment->value > properties->hAdjustment->upper - properties->hAdjustment->page_size)
    properties->hAdjustment->value = properties->hAdjustment->upper - properties->hAdjustment->page_size;

  /* Make sure we're still within the lower limit */
  if (properties->hAdjustment->value < properties->hAdjustment->lower)
    properties->hAdjustment->value = properties->hAdjustment->lower;

  belvuAlignmentRedrawSequenceArea(belvuAlignment);
}

/* Called when the vertical scroll range has changed */
static void onVScrollRangeChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  /* Set the horizontal adjustment page size to be the number of characters
   * that will fit in the width */
  properties->vAdjustment->page_size = (int)(properties->seqArea->allocation.height / properties->charHeight);
  properties->vAdjustment->page_increment = properties->vAdjustment->page_size;
  
  /* Make sure the scroll position still lies within upper - page_size. */
  if (properties->vAdjustment->value > properties->vAdjustment->upper - properties->vAdjustment->page_size)
    properties->vAdjustment->value = properties->vAdjustment->upper - properties->vAdjustment->page_size;
  
  /* Make sure we're still within the lower limit */
  if (properties->vAdjustment->value < properties->vAdjustment->lower)
    properties->vAdjustment->value = properties->vAdjustment->lower;

  belvuAlignmentRedrawAll(belvuAlignment);
}


static void updateOnAlignmentSizeChanged(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->vAdjustment)
    {
      properties->vAdjustment->upper = properties->bc->alignArr->len + 1;
      gtk_adjustment_changed(properties->vAdjustment);
    }
  
  if (properties->hAdjustment)
    {
      properties->hAdjustment->upper = properties->bc->maxLen;
      gtk_adjustment_changed(properties->hAdjustment);
    }

  belvuAlignmentRedrawAll(belvuAlignment);
}


/* Called when the vertical scroll range upper value has changed (i.e. the length
 * of the alignments array has changed). */
void updateOnVScrollSizeChaged(GtkWidget *belvuAlignment)
{
  /* Recalculate the borders, because the max name length etc. may have changed */
  calculateDrawingSizes(belvuAlignment);
  
  /* Update the scroll ranges because the number of rows and columns may have changed */
  updateOnAlignmentSizeChanged(belvuAlignment);
}


/* This should be called when the alignment length has changed, i.e. when columns
 * have been deleted. */
void updateOnAlignmentLenChanged(GtkWidget *belvuAlignment)
{
  updateOnAlignmentSizeChanged(belvuAlignment);
}

/***********************************************************
 *                         Sizing                          *
 ***********************************************************/

int belvuAlignmentGetWidth(GtkWidget *belvuAlignment)
{
  if (!belvuAlignment)
    return 0;

  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  return properties ? properties->seqRect.width : 0;
}


/* Get the width of the drawing area that shows the sequences */
static int getAlignmentDisplayWidth(BelvuAlignmentProperties *properties)
{
  int result = 0;
  
  if (properties->wrapWidth == UNSET_INT)
    {
      /* Drawing can be very slow if we draw the full alignment, so we only
       * draw what's in the current display width. */
      result = properties->seqArea->allocation.width - SEQ_AREA_PADDING;
    }
  else
    {
      /* The sequences are wrapped, so we only have a short width to deal with
       * and can therefore draw the full width. (Although specify a size limit
       * in case the user specifies a very large width because this can cause
       * a badalloc crash if too big.) */
      result =
        properties->columnPadding + (properties->bc->maxNameLen * properties->charWidth) + 
        properties->nameColumnPadding + (properties->bc->maxEndLen * properties->charWidth) +
        properties->columnPadding + (properties->wrapWidth * properties->charWidth) + 
        properties->columnPadding + (properties->bc->maxEndLen * properties->charWidth) + 
        properties->columnPadding + (properties->bc->maxScoreLen * properties->charWidth);
      
      if (result > MAX_PIXMAP_WIDTH)
        {
          result = MAX_PIXMAP_WIDTH;
          g_warning("The alignment window is too large and will be clipped.\n");
        }
    }
  
  return result;
}


static int getAlignmentDisplayHeight(BelvuAlignmentProperties *properties)
{
  int result = 0;
  
  if (properties->wrapWidth == UNSET_INT)
    {
      result = properties->seqArea->allocation.height - DEFAULT_YPAD;
    }  
  else
    {
      /* Divide the full sequence length by the display width to give the number
       * of paragraphs we require */
      const int numParagraphs = ceil((double)properties->bc->maxLen / (double)properties->wrapWidth);
      const int paragraphHeight = (properties->bc->alignArr->len + 1) * properties->charHeight;
    
      result = paragraphHeight * numParagraphs;
      
      /* Add space for the title, if given (one line for the title and one as a spacer) */
      if (properties->title)
        result += properties->charHeight * 2;
      
      if (result > MAX_PIXMAP_HEIGHT)
        {
          result = MAX_PIXMAP_HEIGHT;
          g_warning("The alignment window is too large and will be clipped.\n");
        }
    }
  
  return result;
}


/* Update the size of the header-columns widget. Should be called after any
 * change that affects the length of the text displayed in the columns area. */
void updateHeaderColumnsSize(GtkWidget *belvuAlignment)
{
  if (!belvuAlignment)
    return;
  
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->columnsArea)
    {
      /* There is a separate drawing area for the columns that contain the row
       * headers; calculate its size. */
      properties->columnsRect.width = (properties->bc->maxNameLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxStartLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxEndLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      properties->columnPadding;
      
      if (properties->bc->displayScores)
        properties->columnsRect.width += (properties->bc->maxScoreLen * properties->charWidth) + (2 * properties->columnPadding);

      /* Set the height of the header line to the the same as the sequence header line */
      gtk_widget_set_size_request(properties->columnsHeader, -1, properties->seqHeaderRect.y + properties->seqHeaderRect.height);
      
      /* Set the height and width of the drawing area */
      gtk_layout_set_size(GTK_LAYOUT(properties->columnsArea), properties->columnsRect.width, properties->columnsRect.height);
      
      /* Force the width of the drawing widget to be the width of the drawing area */
      if (properties->columnsRect.width != properties->columnsArea->allocation.width)
        gtk_widget_set_size_request(properties->columnsArea, properties->columnsRect.width, -1);
    }
}


/* This function calculates the size of the drawing area for sequences based
 * on either a) the current allocated size of the sequence area or b) for a 
 * wrapped alignment, the full size required to draw everything.
 * It also updates the height of the columns area, if applicable, so that it is
 * the same height as the sequence area. */
void calculateDrawingSizes(GtkWidget *belvuAlignment)
{
  if (!belvuAlignment)
    return;
  
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  /* Recalculate the size of the drawing area for the sequences (both height
   * and width may have changed on a size-allocate) */
  properties->seqRect.width = getAlignmentDisplayWidth(properties);
  properties->seqRect.height = getAlignmentDisplayHeight(properties);
  
  /* Update the size of the drawing widget and signal that the scroll range
   * has changed (n/a for wrapped view because this gets done by the draw function. */
  if (properties->wrapWidth == UNSET_INT)
    {
      gtk_layout_set_size(GTK_LAYOUT(properties->seqArea), properties->seqRect.x + properties->seqRect.width, properties->seqRect.y + properties->seqRect.height);

      if (properties->hAdjustment)
        gtk_adjustment_changed(properties->hAdjustment); /* signal that the scroll range has changed */
      
      if (properties->vAdjustment)
        gtk_adjustment_changed(properties->vAdjustment);
    }
  
  /* Keep the sequence-header area the same width as the sequence area */
  properties->seqHeaderRect.width = properties->seqRect.width;

  if (properties->columnsArea)
    {
      /* Update the height of the drawing area (the width shouldn't have
       * changed for a size-allocate because the table cell does not expand) */
      properties->columnsRect.height = properties->seqRect.height;
      gtk_layout_set_size(GTK_LAYOUT(properties->columnsArea), properties->columnsRect.width, properties->columnsRect.height);
    }
}


static void onSizeAllocateBelvuAlignment(GtkWidget *widget, GtkAllocation *allocation, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  
  /* Update the size of the drawing areas */
  calculateDrawingSizes(belvuAlignment);
}


/***********************************************************
 *                Removing sequences                       *
 ***********************************************************/

void removeSelectedSequence(BelvuContext *bc, GtkWidget *belvuAlignment)
{ 
  if (!bc->selectedAln) 
    {
      g_critical("Please select a sequence to remove.\n");
      return;
    }
  
  const int idx = bc->selectedAln->nr - 1;
  g_array_remove_index(bc->alignArr, idx);
  arrayOrder(bc->alignArr);
  
  bc->saved = FALSE;
  
  g_message("Removed %s/%d-%d.  %d sequences left.\n\n", bc->selectedAln->name, bc->selectedAln->start, bc->selectedAln->end, bc->alignArr->len);
  
  bc->selectedAln = NULL;
  
  rmFinaliseGapRemoval(bc);
  updateOnVScrollSizeChaged(belvuAlignment);  
  onRowSelectionChanged(bc);
}


/* Get rid of seqs that are too gappy. */
void removeGappySeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmGappySeqs(bc, cutoff);
  rmFinaliseGapRemoval(bc);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of partial seqs. */
void removePartialSeqs(BelvuContext *bc, GtkWidget *belvuAlignment)
{
  rmPartialSeqs(bc);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of redundant seqs (those that are more than the given % identical). */
void removeRedundantSeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  mkNonRedundant(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of outlier seqs (those that are less than the given % identical to any other). */
void removeOutliers(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmOutliers(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of seqs that have a score below the given value */
void removeByScore(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmScore(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);
  centerHighlighted(bc, belvuAlignment);
}


/* Scroll the alignment list if necessary to make sure that the currently
 * highlighted sequence is visible */
void centerHighlighted(BelvuContext *bc, GtkWidget *belvuAlignment)
{
  if (bc->selectedAln) 
    {
      BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
      GtkAdjustment *vAdjustment = properties->vAdjustment;
      
      /* Create the range of y coords that we want to be in range */
      const int newIdx = bc->selectedAln->nr - 1;
      gtk_adjustment_clamp_page(vAdjustment, newIdx, newIdx + 1);
    }
}


/* Utility to set the value of an adjustment and to bounds-limit it to lie
 * between 'lower' and 'upper - page_size' */
static void setAdjustmentValue(GtkAdjustment *adjustment, const int value)
{
  int newValue = value;

  /* _set_value checks the lower limit, but we need to check the upper... */
  if (newValue > adjustment->upper - adjustment->page_size)
    newValue = adjustment->upper - adjustment->page_size;
                     
  gtk_adjustment_set_value(adjustment, newValue);
}


/* Utility to scroll an adjustment to the start/end of its range */
static void adjustmentScrollExtremity(GtkAdjustment *adjustment, const gboolean start)
{
  if (start)
    setAdjustmentValue(adjustment, adjustment->lower);
  else
    setAdjustmentValue(adjustment, adjustment->upper);
}

/* Utility to scroll an adjustment up or down by one page */
static void adjustmentScrollPage(GtkAdjustment *adjustment, const gboolean lower)
{
  if (lower)
    setAdjustmentValue(adjustment, adjustment->value - adjustment->page_size);
  else
    setAdjustmentValue(adjustment, adjustment->value + adjustment->page_size);
}

/* Utility to scroll an adjustment up or down by the given value */
static void adjustmentScrollValue(GtkAdjustment *adjustment, const gboolean lower, const int value)
{
  if (lower)
    setAdjustmentValue(adjustment, adjustment->value - value);
  else
    setAdjustmentValue(adjustment, adjustment->value + value);
}


/* Scroll vertically to the start (top) or end (bottom) of the alignment list */
void vScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollExtremity(properties->vAdjustment, start);
}

/* Scroll horizontally to the start (leftmost) or end (rightmost) coord */
void hScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollExtremity(properties->hAdjustment, start);
}

/* Scroll vertically one page up or down */
void vScrollPageUpDown(GtkWidget *belvuAlignment, const gboolean up)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollPage(properties->vAdjustment, up);
}

/* Scroll horizontally one page left or right */
void hScrollPageLeftRight(GtkWidget *belvuAlignment, const gboolean left)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollPage(properties->hAdjustment, left);
}

/* Scroll vertically the given number of rows up or down */
void vScrollUpDown(GtkWidget *belvuAlignment, const gboolean up, const int numRows)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollValue(properties->vAdjustment, up, numRows);
}

/* Scroll horizontally the given number of characters left or right */
void hScrollLeftRight(GtkWidget *belvuAlignment, const gboolean left, const int numChars)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  adjustmentScrollValue(properties->hAdjustment, left, numChars);
}




/* Select the row at the given y coord */
static void selectRowAtCoord(BelvuAlignmentProperties *properties, const int y)
{
  BelvuContext *bc = properties->bc;
  
  const int rowIdx = properties->vAdjustment->value + ((y - properties->seqRect.y) / properties->charHeight);
  
  /* Set the selected alignment. Note that some alignments are hidden so
   * we need to count how many visible alignments there are up to this row. */
  int i = 0;
  int count = -1;
  bc->selectedAln = NULL;
  
  for ( i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (!alnp->hide)
        {
          ++count;
          if (count == rowIdx)
            {
              bc->selectedAln = alnp;
              break;
            }
        }
    }
}


/* Utility to find the column index (0-based from left edge) at the given
 * x coordinate in the sequence area of the alignment view. */
static int getColumnAtCoord(BelvuAlignmentProperties *properties, const int x)
{
  int result = (x - properties->seqRect.x) / properties->charWidth;
  return result;
}


/* Select the column at the given x coord */
static void selectColumnAtCoord(BelvuAlignmentProperties *properties, const int x, const gboolean highlightCol)
{
  BelvuContext *bc = properties->bc;

  const int colIdx = getColumnAtCoord(properties, x);
  
  if (colIdx >= 0 && colIdx < bc->maxLen)
    {
      bc->selectedCol = colIdx + properties->hAdjustment->value + 1;
      
      if (highlightCol)
        bc->highlightedCol = bc->selectedCol;
      else
        bc->highlightedCol = 0;
    }
}


/* Mouse button handler for the columns area of the alignment window */
static gboolean onButtonPressColumnsArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 1)  /* single click left button */
    {
      /* Select the clicked row */
      selectRowAtCoord(properties, event->y);
      onRowSelectionChanged(properties->bc);
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
        {
          /* Removed the clicked sequence (which will be the selected one) */
          removeSelectedSequence(properties->bc, belvuAlignment);
        }
      else 
        {
          /* Fetch the clicked sequence (i.e. the currently selected one) */
          fetchAln(properties->bc, properties->bc->selectedAln);
        }
    }
  
  return handled;
}


/* Mouse button handler for the sequence area of the alignment window */
static gboolean onButtonPressSeqArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->type == GDK_BUTTON_PRESS &&
      (event->button == 1 || event->button == 2))  /* single click left or middle buttons */
    {
      /* If the middle button was pressed, highlight the selected column */
      const gboolean highlightCol = (event->button == 2);
      
      /* If the left button was pressed, select the clicked row */
      if (event->button == 1)
        {
          selectRowAtCoord(properties, event->y);
          onRowSelectionChanged(properties->bc);
        }
      
      /* Select the clicked column */
      selectColumnAtCoord(properties, event->x, highlightCol);
      onColSelectionChanged(properties->bc);
      
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
        {
          /* Removed the clicked sequence (i.e. the currently selected one) */
          removeSelectedSequence(properties->bc, belvuAlignment);
        }
      else 
        {
          /* Fetch the clicked sequence (i.e. the currently selected one) */
          fetchAln(properties->bc, properties->bc->selectedAln);
        }
    }
  
  return handled;
}


/* Mouse button handler release for the sequence area of the alignment window */
static gboolean onButtonReleaseSeqArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->button == 2)  /* released middle button */
    {
      /* Scroll to centre on the current column */
      int colIdx = getColumnAtCoord(properties, event->x);
      double newValue = colIdx + properties->hAdjustment->value - (properties->hAdjustment->page_size / 2.0);
      
      if (newValue > properties->hAdjustment->upper - properties->hAdjustment->page_size)
        newValue = properties->hAdjustment->upper - properties->hAdjustment->page_size;
        
      gtk_adjustment_set_value(properties->hAdjustment, newValue);
      
      /* Clear the current column highlihting */
      properties->bc->highlightedCol = 0;
      onColSelectionChanged(properties->bc);
      
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
        {
          /* Removed the clicked sequence (which will be the selected one) */
          removeSelectedSequence(properties->bc, belvuAlignment);
        }
    }
  
  return handled;
}


/* Mouse move handler for the sequence area of the alignment window */
static gboolean onMouseMoveSeqArea(GtkWidget *widget, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->state & GDK_BUTTON2_MASK) /* middle button pressed */
    {
      /* Update the selected/highlighted column */
      selectColumnAtCoord(properties, event->x, TRUE);
      onColSelectionChanged(properties->bc);
      
      handled = TRUE;
    }
  
  return handled;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the alignment. */
static gboolean onScrollAlignment(GtkWidget *widget, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  switch (event->direction)
    {
      /* left/right scrolling for sequence area or sequence header only */
      case GDK_SCROLL_LEFT:
        {
          if (widget == properties->seqArea || widget == properties->seqHeader)
            hScrollLeftRight(belvuAlignment, TRUE, properties->hAdjustment->step_increment);
          handled = TRUE;
          break;
        }
        
      case GDK_SCROLL_RIGHT:
        {
          if (widget == properties->seqArea || widget == properties->seqHeader)
            hScrollLeftRight(belvuAlignment, FALSE, properties->hAdjustment->step_increment);
          handled = TRUE;
          break;
        }

      /* up/down scrolling for sequence area or columns area only */
      case GDK_SCROLL_UP:
        {
          if (widget == properties->seqArea || widget == properties->columnsArea)
            vScrollUpDown(belvuAlignment, TRUE, properties->vAdjustment->step_increment);
          handled = TRUE;
          break;
        }
        
      case GDK_SCROLL_DOWN:
        {
          if (widget == properties->seqArea || widget == properties->columnsArea)
            vScrollUpDown(belvuAlignment, FALSE, properties->vAdjustment->step_increment);
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
 *                      Initialisation                     *
 ***********************************************************/

/* Set style properties for the belvu alignment widgets */
static void setBelvuAlignmentStyle(BelvuContext *bc, GtkWidget *seqArea, GtkWidget *columnsArea)
{
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);

  gtk_widget_modify_bg(seqArea, GTK_STATE_NORMAL, bgColor);

  if (columnsArea)
    gtk_widget_modify_bg(columnsArea, GTK_STATE_NORMAL, bgColor);
}


/* Create a widget that will draw the alignments. This type of widget is
 * used for both the standard and wrapped views - pass wrapWidth as UNSET_INT
 * for the standard view, or pass wrapWidth as the number of characters after
 * which to wrap for the wrapped view. */
GtkWidget* createBelvuAlignment(BelvuContext *bc, const char* title, const int wrapWidth)
{
  GtkWidget *belvuAlignment = NULL;      /* the outermost container */
  
  GtkWidget *seqArea = NULL;             /* the "sequence area", which will draw the actual peptide sequences */
  GtkWidget *columnsArea = NULL;         /* optional "headers" column for drawing names, coords etc. */
  GtkWidget *seqHeader = NULL;       /* header for the sequence area */
  GtkWidget *columnsHeader = NULL;   /* header for the columns area */
  
  GtkAdjustment *hAdjustment = NULL;
  GtkAdjustment *vAdjustment = NULL;
  
  if (wrapWidth == UNSET_INT)
    {
      /* The standard (unwrapped) view. Two columns that share the same vertical
       * scrollbar but only the sequence area will have a horizontal scrollbar. */
      
      /* Place everything in a table */
      belvuAlignment = gtk_table_new(3, 3, FALSE);
      const int xpad = 0, ypad = 0;

      /* Create the scrollbars */
      hAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, bc->maxLen, 5, 0, 0));
      vAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, bc->alignArr->len + 1, 1, 0, 0));
      GtkWidget *hScrollbar = gtk_hscrollbar_new(hAdjustment);
      GtkWidget *vScrollbar = gtk_vscrollbar_new(vAdjustment);
      
      /* Create the drawing areas */
      seqArea = gtk_layout_new(NULL, NULL);
      columnsArea = gtk_layout_new(NULL, NULL);
      
      /* Create a header widget for each drawing area */
      columnsHeader = gtk_drawing_area_new();
      seqHeader = gtk_drawing_area_new();
      
      /* Place all the widgets in the table */
      gtk_table_attach(GTK_TABLE(belvuAlignment), columnsHeader,0, 1, 0, 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), columnsArea,  0, 1, 1, 2, GTK_FILL, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), seqHeader,    1, 2, 0, 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), seqArea,      1, 2, 1, 2, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), hScrollbar,   1, 2, 2, 3, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), vScrollbar,   2, 3, 1, 2, GTK_SHRINK, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), xpad, ypad);

      /* Connect signals */
      gtk_widget_add_events(columnsArea, GDK_BUTTON_PRESS_MASK);
      g_signal_connect(G_OBJECT(columnsArea), "expose-event", G_CALLBACK(onExposeBelvuColumns), belvuAlignment);  
      g_signal_connect(G_OBJECT(columnsArea), "button-press-event", G_CALLBACK(onButtonPressColumnsArea), belvuAlignment);
      g_signal_connect(G_OBJECT(columnsArea), "scroll-event",  G_CALLBACK(onScrollAlignment), belvuAlignment);
      
      g_signal_connect(G_OBJECT(columnsHeader), "expose-event", G_CALLBACK(onExposeBelvuColumnsHeader), belvuAlignment);  

      gtk_widget_add_events(seqHeader, GDK_BUTTON_PRESS_MASK);
      g_signal_connect(G_OBJECT(seqHeader), "expose-event", G_CALLBACK(onExposeBelvuSequenceHeader), belvuAlignment);  
      g_signal_connect(G_OBJECT(seqHeader), "scroll-event",  G_CALLBACK(onScrollAlignment), belvuAlignment);

      gtk_widget_add_events(seqArea, GDK_BUTTON_RELEASE_MASK);
      gtk_widget_add_events(seqArea, GDK_POINTER_MOTION_MASK);
      g_signal_connect(G_OBJECT(seqArea), "button-press-event", G_CALLBACK(onButtonPressSeqArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "button-release-event", G_CALLBACK(onButtonReleaseSeqArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "motion-notify-event", G_CALLBACK(onMouseMoveSeqArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "scroll-event",  G_CALLBACK(onScrollAlignment), belvuAlignment);

      g_signal_connect(G_OBJECT(hAdjustment), "changed", G_CALLBACK(onHScrollRangeChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(vAdjustment), "changed", G_CALLBACK(onVScrollRangeChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(hAdjustment), "value-changed", G_CALLBACK(onHScrollPosChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(vAdjustment), "value-changed", G_CALLBACK(onVScrollPosChangedBelvuAlignment), belvuAlignment);
    }
  else
    {
      /* For the wrapped view, we just have a single layout area. We'll need to
       * draw the entire contents at once so that it can be printed. We can
       * therefore just place it in a standard scrolled window. We also don't
       * need to worry about button-press events because this widget is not
       * interactive; it's just for printing. */
      seqArea = gtk_layout_new(NULL, NULL);

      belvuAlignment = gtk_scrolled_window_new(NULL, NULL);
      gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(belvuAlignment), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
      gtk_container_add(GTK_CONTAINER(belvuAlignment), seqArea);
      
      hAdjustment = gtk_layout_get_hadjustment(GTK_LAYOUT(seqArea));
      vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(seqArea));
    }

  gtk_widget_add_events(seqArea, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(seqArea), "expose-event",  G_CALLBACK(onExposeBelvuSequence), belvuAlignment);  
  g_signal_connect(G_OBJECT(seqArea), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), belvuAlignment);
  
  /* Set the style and properties */
  g_assert(belvuAlignment);
  setBelvuAlignmentStyle(bc, seqArea, columnsArea);
  belvuAlignmentCreateProperties(belvuAlignment, bc, columnsArea, columnsHeader, seqArea, seqHeader, hAdjustment, vAdjustment, title, wrapWidth);

  return belvuAlignment;
}





