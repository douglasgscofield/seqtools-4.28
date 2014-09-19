/*  File: blxWindow.h
 *  Author: Gemma Barson, 2009-11-24
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
 * Description: Creates the main Blixem window. Also creates a "context", 
 *              which contains all of the variables associated with a Blixem
 *              session. 
 *
 *              The context could live somewhere else, but was just included
 *              here because there is one context for each Blixem window.
 *              Ideally the context would be created earlier (perhaps in the
 *              main function) but so far that has not been practical.
 *----------------------------------------------------------------------------
 */

#ifndef _blxwindow_included_
#define _blxwindow_included_

#include <gtk/gtk.h>
#include <blixemApp/blixem_.h>
#include <seqtoolsUtils/utilities.h>


/* Public function declarations */
BlxViewContext*		  blxWindowGetContext(GtkWidget *widget);
GList*                    blxWindowGetColumnList(GtkWidget *blxWindow);
gboolean		  blxWindowGetDisplayRev(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetBigPicture(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetDetailView(GtkWidget *blxWindow);
GtkWidget*                blxWindowGetCoverageView(GtkWidget *blxWindow);
GtkWidget*		  blxWindowGetMainMenu(GtkWidget *blxWindow);
BlxBlastMode		  blxWindowGetBlastMode(GtkWidget *blxWindow);
IntRange*		  blxWindowGetFullRange(GtkWidget *blxWindow);
IntRange*		  blxWindowGetRefSeqRange(GtkWidget *blxWindow);
const char*		  blxWindowGetRefSeqName(GtkWidget *blxWindow);
BlxSeqType		  blxWindowGetSeqType(GtkWidget *blxWindow);
char**			  blxWindowGetGeneticCode(GtkWidget *blxWindow);
char*			  blxWindowGetRefSeq(GtkWidget *blxWindow);
int			  blxWindowGetNumFrames(GtkWidget *blxWindow);
int			  blxWindowGetDotterStart(GtkWidget *blxWindow);
int			  blxWindowGetDotterEnd(GtkWidget *blxWindow);
int			  blxWindowGetDotterZoom(GtkWidget *blxWindow);
int			  blxWindowGetAutoDotter(GtkWidget *blxWindow);
MSP*			  blxWindowGetMspList(GtkWidget *blxWindow);
GList*			  blxWindowGetAllMatchSeqs(GtkWidget *blxWindow);
GList*			  blxWindowGetSequenceGroups(GtkWidget *blxWindow);
SequenceGroup*		  blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind);
const char*		  blxWindowGetPaddingSeq(GtkWidget *blxWindow);
int			  blxWindowGetOffset(GtkWidget *blxWindow);
BlxStrand		  blxWindowGetActiveStrand(GtkWidget *blxWindow);
gboolean                  blxWindowGetNegateCoords(GtkWidget *blxWindow);

GList*                    blxWindowGetSelectedSeqs(GtkWidget *blxWindow);
void                      blxWindowSelectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void                      blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList);
void                      blxWindowDeselectSeq(GtkWidget *blxWindow, BlxSequence *seq);
void                      blxWindowDeselectAllSeqs(GtkWidget *blxWindow);
gboolean                  blxWindowIsSeqSelected(GtkWidget *blxWindow, const BlxSequence *seq);
void                      blxWindowSetSeqSelected(GtkWidget *blxWindow, BlxSequence *seq, const gboolean selected);
void                      blxWindowSelectionChanged(GtkWidget *blxWindow);
BlxSequence*              blxWindowGetLastSelectedSeq(GtkWidget *blxWindow);

gboolean                  blxContextIsSeqSelected(const BlxViewContext* const bc, const BlxSequence *seq);
SequenceGroup*            blxContextGetSequenceGroup(const BlxViewContext *bc, const BlxSequence *seqToFind);

int                       sequenceGetGroupOrder(GtkWidget *blxWindow, const BlxSequence *seq);
void                      copySelectionToClipboard(GtkWidget *blxWindow);
void                      findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data);

void                      refreshDialog(const BlxDialogId dialogId, GtkWidget *blxWindow);
void                      showHelpDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showSettingsDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showSortDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showViewPanesDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups, const gboolean bringToFront);
void                      showFindDialog(GtkWidget *blxWindow, const gboolean bringToFront);
void                      showAboutDialog(GtkWidget *blxWindow);
void                      showInfoDialog(GtkWidget *blxWindow);

void                      blxWindowRedrawAll(GtkWidget *blxWindow);
  
GtkWidget*                createBlxWindow(CommandLineOptions *options, 
                                          const char *paddingSeq, 
                                          GArray* featureLists[], 
                                          GList *seqList, 
                                          GSList *supportedTypes,
                                          const gboolean External,
                                          GSList *styles);


#endif /* _blxwindow_included_ */
