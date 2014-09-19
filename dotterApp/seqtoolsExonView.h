/*  File: seqtoolsExonView.h
 *  Author: Gemma Barson, 2009-12-24
 *  Copyright (c) 2010 - 2012 Genome Research Ltd
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
 * Description: Widget to draw the exons and introns in Dotter.
 *
 *              seqtoolsExonView.c/.h were copied from exonview.c/.h. They
 *              contain much of the same code so the aim was to make these
 *              files generic enough for both Blixem and Dotter to use and to
 *              put them in the seqtoolsUtils library. Unfortunately due to
 *              time constraints this didn't happen, so these files still 
 *              contain some Dotter- and Blixem- specific stuff: in Dotter
 *              we need to draw the crosshair over the exon view; in Blixem
 *              we need to draw the highlight box to indicate the current
 *              detail-view range.
 *----------------------------------------------------------------------------
 */

#ifndef _seqtools_exon_view_included_
#define _seqtools_exon_view_included_

#include <gtk/gtk.h>
#include <seqtoolsUtils/blxmsp.h>
#include <dotterApp/dotter_.h>

#define SEQTOOLS_EXON_VIEW_NAME		"SeqtoolsExonView"
#define DEFAULT_EXON_HEIGHT			 10
#define DEFAULT_EXON_HEIGHT_BUMPED		 7
#define	DEFAULT_EXON_YPAD			 10
#define	DEFAULT_EXON_YPAD_BUMPED		 4


/* Public function declarations */
GtkWidget *createDotterExonView(GtkWidget *parent, 
				GtkCallback refreshFunc,
                                const gboolean horizontal,
				const BlxStrand strand, 
				DotterWindowContext *dwc,
				const int width,
                                const int height,
				const IntRange* const qRange,
                                const gboolean drawCrosshair,
                                GtkWidget **exonViewOut);

//gboolean	exonViewGetBumped(GtkWidget *exonView);
//void		exonViewSetBumped(GtkWidget *exonView, const gboolean bumped);
void		exonViewToggleBumped(GtkWidget *exonView);
//
//void            callFuncOnAllChildExonViews(GtkWidget *widget, gpointer data);
void		calculateDotterExonViewHeight(GtkWidget *exonView);
void		calculateDotterExonViewBorders(GtkWidget *exonView, const int width, const int height);
void            exonViewSetShowCrosshair(GtkWidget *exonView, const gboolean showCrosshair);
void            exonViewPrepareForPrinting(GtkWidget *exonView);

#endif /* _seqtools_exon_view_included_ */
