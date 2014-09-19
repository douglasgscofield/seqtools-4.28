/*  File: belvuAlignment.h
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
 * Description: Draws the alignment section in the main Belvu window
 *----------------------------------------------------------------------------
 */

#ifndef _belvualignment_h_included_
#define _belvualignment_h_included_

#include <belvuApp/belvu_.h>
#include <gtk/gtk.h>

GtkWidget*              createBelvuAlignment(BelvuContext *bc, const char *title, const int wrapWidth);
void                    belvuAlignmentRedrawAll(GtkWidget *belvuAlignment);
void                    belvuAlignmentRefreshAll(GtkWidget *belvuAlignment);
void                    updateOnAlignmentLenChanged(GtkWidget *belvuAlignment);

void                    onBelvuAlignmentFontSizeChanged(GtkWidget *belvuAlignment);

void                    updateOnVScrollSizeChaged(GtkWidget *belvuAlignment);
void                    centerHighlighted(BelvuContext *bc, GtkWidget *belvuAlignment);

void                    removeSelectedSequence(BelvuContext *bc, GtkWidget *belvuAlignment);
void                    removeGappySeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removePartialSeqs(BelvuContext *bc, GtkWidget *belvuAlignment);
void                    removeRedundantSeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removeOutliers(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);
void                    removeByScore(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff);

void                    vScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start);
void                    hScrollStartEnd(GtkWidget *belvuAlignment, const gboolean start);
void                    vScrollPageUpDown(GtkWidget *belvuAlignment, const gboolean up);
void                    hScrollPageLeftRight(GtkWidget *belvuAlignment, const gboolean left);
void                    vScrollUpDown(GtkWidget *belvuAlignment, const gboolean up, const int numRows);
void                    hScrollLeftRight(GtkWidget *belvuAlignment, const gboolean left, const int numChars);

int                     belvuAlignmentGetWidth(GtkWidget *belvuAlignment);

void                    calculateDrawingSizes(GtkWidget *belvuAlignment);
void                    updateHeaderColumnsSize(GtkWidget *belvuAlignment);


#endif /* _belvualignment_h_included_ */
