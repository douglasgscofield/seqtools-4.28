/*  File: detailviewtree.h
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
 * Description: A detail-view tree shows all of the alignments for a particular 
 *              reading-frame and strand of the reference sequence. It 
 *              constitutes one pane in the detail-view. It has a header
 *              showing the reference sequence and, lined up below this, shows
 *              the sequence data for each match sequence within the current
 *              display range.
 *
 *              One row in the tree represents one match sequence. Only matches
 *              that lie within the current display range are visible in the
 *              tree. The display range is specified by selecting a range from
 *              the big picture or by scrolling etc.
 *
 *              All detail-view trees share the same scrollbar/display-range
 *              so that they all show the same section of reference sequence.
 *              The cell contents of the tree are drawn by a custom cell 
 *              renderer - see sequencecellrenderer.h.
 *----------------------------------------------------------------------------
 */

#ifndef _detail_view_tree_included_
#define _detail_view_tree_included_

#include <gtk/gtk.h>
#include <seqtoolsUtils/utilities.h>
#include <blixemApp/sequencecellrenderer.h>
#include <blixemApp/detailview.h>

#define DETAIL_VIEW_TREE_NAME		  "DetailViewTreeName"
#define DETAIL_VIEW_TREE_CONTAINER_NAME	  "DetailViewTreeContainerName"
#define START_TRUE                        1
#define START_FALSE                       2

/* This struct holds info about a tree header widget. */
typedef struct _TreeColumnHeaderInfo
  {
    GtkWidget *headerWidget;	/* The actual header widget */
    GtkWidget *tree;		/* The tree view that this header belongs to */
    GList *columnIds;		/* A list of columns spanned by this header. (Columns must be adjacent to each other.) */
    GtkCallback refreshFunc;    /* Function to be called when a refresh is requested */
  } TreeColumnHeaderInfo;


/* This struct holds info about a detail-view tree pane */
typedef struct _TreeProperties 
  {
    GtkWidget *grid;		    /* The grid that this tree corresponds to */
    GtkWidget *detailView;	    /* The detail view that this tree belongs to */

    int readingFrame;		    /* Which reading frame this tree displays */
    GtkWidget *treeHeader;	    /* The container that contains all the widgets for the tree header */
    GList *treeColumnHeaderList;    /* List of info about the tree column headers */
    gboolean hasSnpHeader;	    /* Whether a SNP track is shown above this tree */

    GtkTreeModel *treeModels[BLXMODEL_NUM_MODELS];  /* The tree data store(s) */
  } TreeProperties;


/* Public function declarations */
BlxViewContext*	  treeGetContext(GtkWidget *tree);
TreeProperties*	  treeGetProperties(GtkWidget *widget);
BlxStrand	  treeGetStrand(GtkWidget *tree);
GList*		  treeGetMsps(GtkTreeModel *model, GtkTreeIter *iter);
GtkTreeModel*	  treeGetBaseDataModel(GtkTreeView *tree);
GtkWidget*	  treeGetBlxWindow(GtkWidget *tree);
int		  treeGetCellXPadding(GtkWidget *tree);
int		  treeGetCellYPadding(GtkWidget *tree);
int               treeGetFrame(GtkWidget *tree);

void		  callFuncOnAllDetailViewTrees(GtkWidget *widget, GtkCallback func, gpointer data);

void		  refilterTree(GtkWidget *tree, gpointer data);
void		  resortTree(GtkWidget *tree, gpointer data);
void		  refreshTreeHeaders(GtkWidget *tree, gpointer data);
void		  resizeTreeColumns(GtkWidget *tree, gpointer data);
void		  treeUpdateFontSize(GtkWidget *tree, gpointer data);

void		  treeUpdateSquashMatches(GtkWidget *tree, gpointer data);

gboolean	  treeMoveRowSelection(GtkWidget *tree, const gboolean moveUp, const gboolean shiftModifier);
void		  treeScrollSelectionIntoView(GtkWidget *tree, gpointer data);

void              addMspToTree(MSP *msp, GtkWidget *tree, GtkListStore *store);
void		  addSequencesToTree(GtkWidget *tree, gpointer data);

void              treeDrawCachedBitmap(GtkWidget *tree, gpointer data);

void              setMouseDragMode(const gboolean value);

GtkWidget*	  createDetailViewTree(GtkWidget *grid, 
				       GtkWidget *detailView,
				       GtkCellRenderer *renderer,
				       GList **treeList,
				       GList *columnList,
				       BlxSeqType seqType,
				       const char* const refSeqName,
				       const int frame,
				       const gboolean includeSnpTrack);

void		   treeCreateBaseDataModel(GtkWidget *tree, gpointer data);
void		   treeCreateFilteredDataModel(GtkWidget *tree, gpointer data);

#endif /* _detail_view_tree_included_ */
