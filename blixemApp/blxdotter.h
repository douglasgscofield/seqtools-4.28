/*  File: blxdotter.h
 *  Author: Gemma Barson, 2010-02-03
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
 * Description: Provides functionality to call Dotter as an external 
 *              application.
 *              
 *              Determines which sequences to call Dotter with based on the
 *              current selection in Blixem. Provides a dialog allowing the
 *              user to set certain parameters and optionally calculates a 
 *              sensible coordinate range to call Dotter with.
 *
 *              Generally Dotter will be called with one alignment sequence vs
 *              a relevant portion of the reference sequence. It can also be
 *              called on the reference sequence vs itself. This is useful if
 *              the reference sequence has been constructed from several 
 *              sequences concatenated together and one wishes to identify
 *              alignments between those sequences.
 *----------------------------------------------------------------------------
 */

#ifndef _blx_dotter_h_included_
#define _blx_dotter_h_included_


#include <gtk/gtk.h>


void			showDotterDialog(GtkWidget *blxWindow, const gboolean bringToFront);
gboolean		callDotter(GtkWidget *blxWindow, const gboolean hspsOnly, char *dotterSSeq, GError **error);
char*                   getDotterSSeq(GtkWidget *blxWindow, GError **error);


#endif /* _blx_dotter_h_included_ */
