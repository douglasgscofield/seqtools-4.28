/*  File: detailviewtree.c
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
 * Description: See detailviewtree.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/detailviewtree.h>
#include <blixemApp/detailview.h>
#include <blixemApp/bigpicturegrid.h>
#include <blixemApp/sequencecellrenderer.h>
#include <blixemApp/blxwindow.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>
#include <math.h>


/* Global variables */

/* This enables "mouse-drag" mode: see the expose function for how this is used. */
static gboolean g_mouse_drag_mode = FALSE;  


/* Local function declarations */
static GtkWidget*	treeGetDetailView(GtkWidget *tree);
static gboolean		onSelectionChangedTree(GObject *selection, gpointer data);
static gint		sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data);
static int		calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree);
static gboolean		isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data);
static gboolean		onExposeRefSeqHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data);
static GList*		treeGetSequenceRows(GtkWidget *tree, const BlxSequence *clickedSeq);
static BlxSequence*	treeGetSequence(GtkTreeModel *model, GtkTreeIter *iter);
static void             destroyTreePathList(GList **list);

/***********************************************************
 *                Tree - utility functions                 *
 ***********************************************************/

/* This sets the "mouse drag mode" flag */
void setMouseDragMode(const gboolean value)
{
  g_mouse_drag_mode = value;
}

/* Check the given widget is a detail-view tree and has its properties set */
static void assertTree(GtkWidget *tree)
{
  if (!tree)
    g_error("Tree is null\n");
  
  if (!GTK_IS_WIDGET(tree))
    g_error("Tree is not a valid widget [%p]\n", tree);
    
  if (!GTK_IS_TREE_VIEW(tree))
    g_error("Tree is not a valid tree view [%p]\n", tree);

  if (strcmp(gtk_widget_get_name(tree), DETAIL_VIEW_TREE_NAME) != 0)
    g_error("Tree is not a valid detail-view tree [%p]\n", tree);

  if (!treeGetProperties(tree))
    g_error("Tree properties not set [widget=%p]\n", tree);
}


static GtkWidget *treeGetDetailView(GtkWidget *tree)
{
  assertTree(tree);
  
  /* optimisation: cache result, because we know there is only ever one detail view */
  static GtkWidget *detailView = NULL;

  if (!detailView)
    {
      TreeProperties *properties = treeGetProperties(tree);
      detailView = properties ? properties->detailView : NULL;
    }
  
  return detailView;
}

GtkWidget *treeGetBlxWindow(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetBlxWindow(detailView);
}

BlxStrand treeGetStrand(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return gridGetStrand(properties->grid);
}

static gboolean treeGetDisplayRev(GtkWidget *tree)
{
  assertTree(tree);
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetDisplayRev(detailView);
}

int treeGetFrame(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties->readingFrame;
}

static gboolean treeHasSnpHeader(GtkWidget *tree)
{
  TreeProperties *properties = treeGetProperties(tree);
  return properties->hasSnpHeader;
}


/* Returns the currently-displayed range in the tree view */
static IntRange* treeGetDisplayRange(GtkWidget *tree)
{
  assertTree(tree);
  TreeProperties *properties = treeGetProperties(tree);
  return properties ? detailViewGetDisplayRange(properties->detailView) : NULL;
}

static PangoFontDescription* treeGetFontDesc(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetFontDesc(detailView);
}

int treeGetCellXPadding(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCellXPadding(detailView);
}

int treeGetCellYPadding(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCellYPadding(detailView);
}

static gdouble treeGetCharWidth(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  return detailViewGetCharWidth(detailView);
}

static TreeColumnHeaderInfo* treeColumnGetHeaderInfo(GtkWidget *tree, BlxColumnId columnId)
{
  /* Loop through each headerinfo item until we find the one that contains this column id */
  TreeProperties *properties = treeGetProperties(tree);
  GList *headerInfoItem = properties->treeColumnHeaderList;
  TreeColumnHeaderInfo *result = NULL;
  
  for ( ; headerInfoItem && !result; headerInfoItem = headerInfoItem->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)(headerInfoItem->data);
      
      /* Loop through each column id included under this header */
      GList *columnIdItem = headerInfo->columnIds;
      
      for ( ; columnIdItem && !result; columnIdItem = columnIdItem->next)
	{
	  if (GPOINTER_TO_INT(columnIdItem->data) == columnId)
	    {
	      result = headerInfo;
	    }
	}
    }
  
  return result;
}

/***********************************************************
 *			utility functions		   *
 ***********************************************************/

/* Call the given function on all trees in the detail view */
void callFuncOnAllDetailViewTrees(GtkWidget *detailView, GtkCallback func, gpointer data)
{
  int numFrames = detailViewGetNumFrames(detailView);
  
  /* Call the function on the forward strand tree and reverse strand tree
   * for each frame. */
  int frame = 1;
  for ( ; frame <= numFrames; ++frame)
    {
      GtkWidget *fwdTree = detailViewGetTree(detailView, BLXSTRAND_FORWARD, frame);
      GtkWidget *revTree = detailViewGetTree(detailView, BLXSTRAND_REVERSE, frame);
      
      if (fwdTree)
	{
	  func(fwdTree, data);
	}
      
      if (revTree)
	{
	  func(revTree, data);
	}
    }
}


/* Add the msps in the given BlxSequence to a single row in the given tree store */
static void addSequenceMspsToSingleRow(BlxSequence *blxSeq, GtkWidget *tree, GtkListStore *store)
{
  /* Only add msps that are in the correct strand for this tree (since the same 
   * sequence may have matches against both ref seq strands) */
  GList *mspsToAdd = NULL;
  GList *mspItem = blxSeq->mspList;
  const BlxStrand treeStrand = treeGetStrand(tree);
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      MSP *msp  = (MSP*)(mspItem->data);
      if (typeShownInDetailView(msp->type) && msp->qStrand == treeStrand && msp->qFrame == treeGetFrame(tree))
        {
          mspsToAdd = g_list_append(mspsToAdd, msp);
        }
    }

  /* Now add a row to the tree store, if there is anything to add */
  if (g_list_length(mspsToAdd) > 0)
    {
      GtkWidget *detailView = treeGetDetailView(tree);
      GList *columnList = detailViewGetColumnList(detailView);

      GtkTreeIter iter;
      gtk_list_store_append(store, &iter);

      /* If there is only one msp, then we can add specific info about that MSP */
      MSP *msp = (g_list_length(mspsToAdd) == 1 ? (MSP*)mspsToAdd->data : NULL);

      /* Add the hard-coded column data */
      const double score = msp ? msp->score : 0.0;
      const double id = msp ? msp->id : 0.0;
      const int start = msp ? msp->sRange.min : blxSequenceGetStart(blxSeq, treeStrand);
      const int end = msp ? msp->sRange.max : blxSequenceGetEnd(blxSeq, treeStrand);

      /* Loop through the rest of the columns */
      GList *item = columnList;
      
      for ( ; item; item = item->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);

          if (columnInfo->columnId == BLXCOL_SCORE)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, score, -1);
          else if (columnInfo->columnId == BLXCOL_ID)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, id, -1);
          else if (columnInfo->columnId == BLXCOL_START)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, start, -1);
          else if (columnInfo->columnId == BLXCOL_SEQUENCE)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, mspsToAdd, -1);
          else if (columnInfo->columnId == BLXCOL_END)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, end, -1);
          else
            {
              GValue *val = blxSequenceGetValue(blxSeq, columnInfo->columnId);
              if (val) gtk_list_store_set_value(store, &iter, columnInfo->columnIdx, val);
            }
        }
      
      /* Remember which row these msps are in */
      GtkTreePath *path = gtk_tree_model_get_path(GTK_TREE_MODEL(store), &iter);
      GList *mspItem = mspsToAdd;
      
      for ( ; mspItem; mspItem = mspItem->next)
        {
          MSP *curMsp = (MSP*)(mspItem->data);
          curMsp->treePaths[BLXMODEL_SQUASHED] = gtk_tree_path_to_string(path);
        }

      gtk_tree_path_free(path);
    }
}


/* Add the msps in the given BlxSequence to the given tree store with a separate
 * row for each msp */
static void addSequenceMspsToSeparateRows(BlxSequence *blxSeq, GtkWidget *tree, GtkListStore *store)
{
  const BlxStrand treeStrand = treeGetStrand(tree);
  GList *mspItem = blxSeq->mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      MSP *msp  = (MSP*)(mspItem->data);

      /* Only add msps that are in the correct strand for this tree (since the same 
       * sequence may have matches against both ref seq strands) */
      if (typeShownInDetailView(msp->type) && 
          msp->qStrand == treeStrand && 
          msp->qFrame == treeGetFrame(tree))
        {
          addMspToTree(msp, tree, store);
        }
    }
}


/* Add the given BlxSequence as a row in the given tree store */
static void addSequenceToTree(BlxSequence *blxSeq, GtkWidget *tree, GtkListStore *store)
{
  /* Only add matches and transcripts to the detail-view. Also, 
   * we exclude sequences with squash-identical-features set because
   * these are added separately. */
  if (!blxSequenceShownInDetailView(blxSeq) ||
      blxSequenceGetFlag(blxSeq, MSPFLAG_SQUASH_IDENTICAL_FEATURES))
    {
      return;
    }

  /* If the squash-linked-features property is set, add all msps in this
   * sequence to the same row; otherwise, add them to separate rows*/
  if (blxSequenceGetFlag(blxSeq, MSPFLAG_SQUASH_LINKED_FEATURES))
    addSequenceMspsToSingleRow(blxSeq, tree, store);
  else
    addSequenceMspsToSeparateRows(blxSeq, tree, store);
}


/* Comparison function for sorting/determining if two alignments are identical;
 * they are identical if the DNA sequence, start positions, source, score and ID
 * are all identical */
static int sortByDnaCompareFunc(gconstpointer a, gconstpointer b)
{
  double result = 0;
  
  const MSP* const msp1 = *((const MSP**)a);
  const MSP* const msp2 = *((const MSP**)b);

  const char *sequence1 = mspGetMatchSeq(msp1);
  const char *sequence2 = mspGetMatchSeq(msp2);

  const gboolean msp1HasSeq = (sequence1 != NULL);
  const gboolean msp2HasSeq = (sequence2 != NULL);

  if (msp1HasSeq && msp2HasSeq)
    {
      result = msp1->qRange.min - msp2->qRange.min;

      if (result == 0) result = getRangeLength(&msp1->sRange) - getRangeLength(&msp2->sRange);
      if (result == 0) result = msp1->score - msp2->score;
      if (result == 0) result = msp1->id - msp2->id;
      if (result == 0) result = strncmp(sequence1 + msp1->sRange.min - 1, sequence2 + msp2->sRange.min - 1, getRangeLength(&msp1->sRange));
    }
  else if (!msp1HasSeq && !msp2HasSeq)
    {
      result = 0;
    }
  else if (!msp1HasSeq)
    {
      result = -1;
    }
  else 
    {
      result = 1;
    }

  /* For fractional results, round up (or down if negative), so that we still
   * recognise a fracional difference as a difference (e.g. a difference of 
   * -0.2 in the scores would mean a < b, so we return -1). */
  if (result > 0)
    return ceil(result);
  else if (result < 0)
    return floor(result);
  else 
    return 0;
}


/* For matches that should be squashed if identical, add them to the given tree, 
 * placing duplicate sequences on the same row as each other to create a compact tree. */
static void addFeaturesToCompactTree(GtkWidget *tree, GtkListStore *store, GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxStrand treeStrand = treeGetStrand(tree);

  /* Extract matches that should be squashed if identical into their own array. */
  GArray *matchArray = bc->featureLists[BLXMSP_MATCH];
  GArray *array = g_array_sized_new(FALSE, FALSE, sizeof(MSP*), matchArray->len);

  int i = 0;
  for ( ; i < (int)matchArray->len; ++i)
    {
      MSP *val = g_array_index(matchArray, MSP*, i);

      if (mspGetFlag(val, MSPFLAG_SQUASH_IDENTICAL_FEATURES))
        g_array_prepend_val(array, val);
    }

  /* Sort the array by DNA sequence so that duplicates are adjacent. */ 
  g_array_sort(array, sortByDnaCompareFunc);

  /* Loop through each MSP, compiling a list of MSPs to add to the
   * current row. */
  i = 0;
  MSP *prevMsp = NULL;
  MSP *msp = mspArrayIdx(array, i);
  GList *mspsToAdd = NULL;
  
  for ( ; msp || prevMsp; msp = mspArrayIdx(array, ++i))
    {
      if (msp && mspGetRefStrand(msp) != treeStrand)
        continue;
      
      if (msp && prevMsp && sortByDnaCompareFunc(&prevMsp, &msp) == 0)
        {
          /* Same as previous sequence; add to current row */
          mspsToAdd = g_list_prepend(mspsToAdd, msp);
        }
      else
        {
          if (g_list_length(mspsToAdd) > 0)
            {
              /* Add the current msp list to the tree */
              GtkTreeIter iter;
              gtk_list_store_append(store, &iter);

              /* Loop through the columns */
              GList *item = blxWindowGetColumnList(blxWindow);
      
              for ( ; item; item = item->next)
                {
                  BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
                  
                  if (columnInfo->columnId == BLXCOL_SCORE)
                    gtk_list_store_set(store, &iter, columnInfo->columnIdx, prevMsp->score, -1);
                  else if (columnInfo->columnId == BLXCOL_ID)
                    gtk_list_store_set(store, &iter, columnInfo->columnIdx, prevMsp->id, -1);
                  else if (columnInfo->columnId == BLXCOL_START)
                    gtk_list_store_set(store, &iter, columnInfo->columnIdx, prevMsp->sRange.min, -1);
                  else if (columnInfo->columnId == BLXCOL_SEQUENCE)
                    gtk_list_store_set(store, &iter, columnInfo->columnIdx, mspsToAdd, -1);
                  else if (columnInfo->columnId == BLXCOL_END)
                    gtk_list_store_set(store, &iter, columnInfo->columnIdx, prevMsp->sRange.max, -1);
                  else
                    {
                      GValue *val = blxSequenceGetValue(msp->sSequence, columnInfo->columnId);
                      if (val) gtk_list_store_set_value(store, &iter, columnInfo->columnIdx, val);
                    }
                }
              
              /* Remember which row these msps are in */
              GtkTreePath *path = gtk_tree_model_get_path(GTK_TREE_MODEL(store), &iter);
              GList *mspItem = mspsToAdd;
          
              for ( ; mspItem; mspItem = mspItem->next)
                {
                  MSP *curMsp = (MSP*)(mspItem->data);
                  curMsp->treePaths[BLXMODEL_SQUASHED] = gtk_tree_path_to_string(path);
                }

              gtk_tree_path_free(path);
          
              /* Reset the list pointer ready for the next row. Note that
               * we do not free the list because it is now owned by the tree */
              mspsToAdd = NULL;
            }
          
          /* Add the new msp */
          if (msp)
            mspsToAdd = g_list_prepend(mspsToAdd, msp);
        }
      
      prevMsp = msp;
    }
  
  /* Re-sort by default sort order (required for filtering) if this msp type
   * appears in the detail-view */
  if (typeShownInDetailView(BLXMSP_MATCH))
    g_array_sort(bc->featureLists[BLXMSP_MATCH], compareFuncMspArray);
}


/* Add all of the BlxSequences that are in the given tree's strand to the tree */
void addSequencesToTree(GtkWidget *tree, gpointer data)
{
  GList *seqList = (GList*)data;
  
  /* Create the data store for the tree. We must know how many columns
   * we have, and their data types. */
  GtkWidget *detailView = treeGetDetailView(tree);  
  GList *columnList = detailViewGetColumnList(detailView);
  const int numCols = g_list_length(columnList);
  GType *typeList = columnListGetTypes(columnList);

  GtkListStore *store = gtk_list_store_newv(numCols, typeList);
  
  /* Set the sort function for each column */
  int colNum = 0;
  for ( ; colNum < numCols; ++colNum)
    {
      gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), colNum,	sortColumnCompareFunc, tree, NULL);
    }

  /* Add the rows - one row per sequence. Use the list we've already compiled of all
   * sequences as BlxSequences */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  GList *seqItem = seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      addSequenceToTree(blxSeq, tree, store);
    }

  /* Also add one row for each match that has duplcate DNA to the compact tree */
  addFeaturesToCompactTree(tree, store, blxWindow);
  
  /* Create a filtered version which will only show sequences that are in the display range */
  GtkTreeModel *filter = gtk_tree_model_filter_new(GTK_TREE_MODEL(store), NULL);
  gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(filter), (GtkTreeModelFilterVisibleFunc)isTreeRowVisible, tree, NULL);
  
  /* Lose the local reference to 'store' because the tree filter has a reference to it. */
  g_object_unref(G_OBJECT(store));
  
  /* Remember the base tree store in the properties so we can switch between this and the 'unsquashed' tree model */
  TreeProperties *properties = treeGetProperties(tree);
  properties->treeModels[BLXMODEL_SQUASHED] = GTK_TREE_MODEL(filter);
  
  /* Note that we don't decrement the ref count to 'filter' even though we're losing
   * the local pointer to it, because we've also added a pointer to it from the tree
   * properties. */
}


/* Return the MSP(s) in a given tree row */
GList* treeGetMsps(GtkTreeModel *model, GtkTreeIter *iter)
{
  GList *mspGList = NULL;  

  /* Get the sequence column index (note this may be different to the column id) */
  GtkWidget *blxWindow = getBlixemWindow();
  GList *columnList = blxWindowGetColumnList(blxWindow);
  BlxColumnInfo *columnInfo = getColumnInfo(columnList, BLXCOL_SEQUENCE);

  if (columnInfo)
    gtk_tree_model_get(model, iter, columnInfo->columnIdx, &mspGList, -1);

  return mspGList;
}


/* For the given tree view, return the data model used for the visible view
 * (which may be a filtered tree model). */
static GtkTreeModel* treeGetVisibleDataModel(GtkTreeView *tree)
{
  /* Just return the model stored directly in the tree view. This will be the filter
   * if there is one, otherwise it will be the original data model. */
  assertTree(GTK_WIDGET(tree));
  return gtk_tree_view_get_model(tree);
}

/* Get the status bar that this tree should report details about moused-over rows to.
 * This may be different to the main window status bar. */
static GtkStatusbar* treeGetStatusBar(GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  return GTK_STATUSBAR(dvProperties->statusBar);
  
  //GTK_STATUSBAR(treeGetContext(tree)->statusBar);
}

/* For the given tree view, return the current base data model i.e. all data in the tree
 * without any filtering. Returns the same as treeGetVisibleDataModel if the tree does
 * not have a filter. */
GtkTreeModel* treeGetBaseDataModel(GtkTreeView *tree)
{
  assertTree(GTK_WIDGET(tree));
  
  GtkTreeModel *result = NULL;
  
  if (tree)
    {
      GtkTreeModel *model = gtk_tree_view_get_model(tree);
      if (model && GTK_IS_TREE_MODEL_FILTER(model))
        {
          result = gtk_tree_model_filter_get_model(GTK_TREE_MODEL_FILTER(model));
        }
      else
        {
          result = model;
        }
    }
  
  return result;
}


/* Decrease the font size in the detail view trees (i.e. effectively zoom out) */
void treeUpdateFontSize(GtkWidget *tree, gpointer data)
{
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);  
  gtk_widget_modify_font(tree, fontDesc);
}


/* This function updates the tree following a change in which tree model we're viewing */
void treeUpdateSquashMatches(GtkWidget *tree, gpointer data)
{
  BlxViewContext *bc = treeGetContext(tree);
  TreeProperties *properties = treeGetProperties(tree);
  
  /* Find the new model */
  GtkTreeModel *newModel = properties->treeModels[bc->modelId];
  
  if (newModel)
    gtk_tree_view_set_model(GTK_TREE_VIEW(tree), newModel);
  
  /* Re-sort and re-filter, because the new one might not be up to date */
  /* Note that we sort the base data model, not the filtered one (i.e. sort all rows, 
   * not just visible ones) */
  newModel = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
  
  if (newModel)
    {
      resortTree(tree, NULL);
      refilterTree(tree, NULL);
    }
}


/***********************************************************
 *			      Updates			   *
 ***********************************************************/


/* Refresh the tree header widgets */
void refreshTreeHeaders(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);

  /* Loop through all widgets in the header and call refreshTextHeader. This
   * updates the font etc. if it is a type of widget that requires that. */
  gtk_container_foreach(GTK_CONTAINER(properties->treeHeader), refreshTextHeader, treeGetDetailView(tree));
  
  /* Loop through all column headers and call their individual refresh functions
   * (this sets the specific data for the column headers) */
  GList *headerItem = properties->treeColumnHeaderList;
  PangoFontDescription *fontDesc = treeGetFontDesc(tree);
  
  for ( ; headerItem; headerItem = headerItem->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)headerItem->data;
      
      if (headerInfo && headerInfo->headerWidget)
	{
	  /* Set the background color  */
	  BlxViewContext *bc = treeGetContext(tree);
	  GdkColor *bgColor = getGdkColor(BLXCOLOR_REF_SEQ, bc->defaultColors, FALSE, bc->usePrintColors);
	  gtk_widget_modify_bg(headerInfo->headerWidget, GTK_STATE_NORMAL, bgColor);
	  gtk_widget_modify_bg(headerInfo->headerWidget, GTK_STATE_NORMAL, bgColor);

	  /* Update the font, in case its size has changed */
	  labelSetFont(headerInfo->headerWidget, fontDesc);

	  /* Call the refresh function, if there is one */
	  if (headerInfo->refreshFunc)
	    {
	      headerInfo->refreshFunc(headerInfo->headerWidget, tree);
	    }
	}
    }
}


/* Resize tree header widgets (should be called after the column width has changed) */
static void resizeTreeHeaders(GtkWidget *tree, gpointer data)
{
  TreeProperties *properties = treeGetProperties(tree);
  GList *header = properties->treeColumnHeaderList;
  
  for ( ; header; header = header->next)
    {
      TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)header->data;
      
      if (headerInfo && headerInfo->headerWidget && g_list_length(headerInfo->columnIds) > 0)
	{
	  int firstColId = GPOINTER_TO_INT(headerInfo->columnIds->data);

          /* For the sequence column, don't set the size request to the real size, or we 
           * won't be able to shrink the window. (The sequence col header will be resized
           * dynamically anyway.) */
	  if (firstColId != BLXCOL_SEQUENCE)
	    {
	      /* For other columns, we can set the size request: they're small enough
	       * that we can live without the need to shrink below their sum widths. */
	      const int width = calculateColumnWidth(headerInfo, tree);
	      gtk_widget_set_size_request(headerInfo->headerWidget, width, -1);
	    }
	}
    }
  
  refreshTreeHeaders(tree, NULL);
}


/* Resize the columns for this tree. Should be called after the column width
 * is changed manually by the user. */
void resizeTreeColumns(GtkWidget *tree, gpointer data)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  
  GList *listItem = detailViewGetColumnList(detailView);

  /* Loop through each column in the tree and set the column width and visibility
   * based on what's stored in our column info. */
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);

      /* We don't set the width of the sequence column - this is an autosize column, so it will 
       * be updated dynamically when any of the other columns change. */
//      if (columnInfo->columnId != BLXCOL_SEQUENCE)
	{
	  GtkTreeViewColumn *treeColumn = gtk_tree_view_get_column(GTK_TREE_VIEW(tree), columnInfo->columnIdx);
      
	  int width = columnInfo->width;
	  if (columnInfo->columnId == BLXCOL_END)
	    {
	      width -= scrollBarWidth();
	    }

	  if (width > 0)
	    {
              const gboolean displayColumn = showColumn(columnInfo);
	      gtk_tree_view_column_set_visible(treeColumn, displayColumn);
	      gtk_tree_view_column_set_fixed_width(treeColumn, width);
	    }
	  else
	    {
	      /* Can't have 0 width, so hide the column instead */
	      gtk_tree_view_column_set_visible(treeColumn, FALSE);
	    }
	}
    }

  resizeTreeHeaders(tree, NULL);
}


/***********************************************************
 *			  Selections			   *
 ***********************************************************/


/* Scrolls the currently-selected row(s) into view */
void treeScrollSelectionIntoView(GtkWidget *tree, gpointer data)
{
  /* Get the last selected sequence */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  const BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
  
  /* Get the tree row(s) that contain that sequence. (If multiple, just use 1st one) */
  GList *rows = treeGetSequenceRows(tree, lastSelectedSeq);
  
  if (g_list_length(rows) > 0)
    {
      GtkTreePath *path = (GtkTreePath*)(rows->data);
      gtk_tree_view_scroll_to_cell(GTK_TREE_VIEW(tree), path, NULL, FALSE, 0.0, 0.0);
    }

  destroyTreePathList(&rows);
}


/* We handle row selections ourselves */
static gboolean onSelectionChangedTree(GObject *selection, gpointer data)
{
  return TRUE;
}


/***********************************************************
 *                   Sorting and filtering                 *
 ***********************************************************/

/* Refilter the data for the given tree */
void refilterTree(GtkWidget *tree, gpointer data)
{
  assertTree(tree);
  
  GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(tree));
  gtk_tree_model_filter_refilter(GTK_TREE_MODEL_FILTER(model));
}


/* Update the cached path pointer for each msp in this tree row. Should be called
 * after the paths have changed, e.g. after a sort. */
static gboolean updateMspPaths(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  BlxModelId modelId = (BlxModelId)GPOINTER_TO_INT(data);
  GList *mspItem = treeGetMsps(model, iter);
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      MSP *msp = (MSP*)(mspItem->data);
      
      /* clear any existing path string */
      if (msp->treePaths[modelId])
        g_free(msp->treePaths[modelId]);
      
      msp->treePaths[modelId] = gtk_tree_path_to_string(path);
    }
  
  return FALSE;
}


/* Re-sort the data for the given tree */
void resortTree(GtkWidget *tree, gpointer data)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);
  const int numColumns = g_list_length(columnList);

  if (numColumns < 1)
    return; 

  /* Find the main column to sort by. We set this as the sort column on the tree.
   * It actually doesn't make a lot of difference which column we set except that 
   * we call the correct sort-by function; the sort-by function will actually sort 
   * by multiple columns based on the detail-view properties. */
  int sortColumn = dvProperties->sortColumns[0];
  
  /* The column ID enum includes some values that are not valid tree columns, i.e. those
   * outside the enum for the max number of columns. If we've got one of these, then
   * set the tree to be unsorted. */
  if (sortColumn >= numColumns)
    sortColumn = GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID;
  
  /* Note that the sort function takes care of whether it's asc or desc, so we
   * always set asc. */
  GtkSortType sortOrder = GTK_SORT_ASCENDING;
  
  /* Not sure if there's a better way to do this, but we can force a re-sort by 
   * setting the sort column to something else and then back again. */
  GtkTreeModel *model = treeGetBaseDataModel(GTK_TREE_VIEW(tree));
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), GTK_TREE_SORTABLE_UNSORTED_SORT_COLUMN_ID, sortOrder);
  gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(model), sortColumn, sortOrder);

  
  /* Update the cached path held by each MSP about the row it is in. */
  BlxViewContext *bc = treeGetContext(tree);
  gtk_tree_model_foreach(model, updateMspPaths, GINT_TO_POINTER(bc->modelId));
}


/* Utility that returns true if the given MSP is currently shown in the tree with the given
 * strand/frame */
static gboolean isMspVisible(const MSP* const msp, 
			     const BlxViewContext *bc, 
			     const int frame, 
			     const IntRange* const displayRange,
			     const int numUnalignedBases,
                             const gboolean seqSelected,
                             const SequenceGroup *group)
{
  g_return_val_if_fail(msp && msp->sSequence, FALSE) ;

  gboolean result = TRUE ;

  /* If hiding ungrouped sequences and this is a sequence (i.e. match) without a group, hide it */
  if (!group && bc->flags[BLXFLAG_HIDE_UNGROUPED_SEQS] && msp->sSequence->type == BLXSEQUENCE_MATCH)
    result = FALSE ;

  /* If hiding ungrouped features and this is feature (i.e. anything except a sequence) without
   * a group, hide it */
  if (!group && bc->flags[BLXFLAG_HIDE_UNGROUPED_FEATURES] && msp->sSequence->type != BLXSEQUENCE_MATCH)
    result = FALSE ;

  if (result)
    {
      /* Check the MSP in the current display range. Get the full MSP display range including
       * any portions outside the actual alignment. */
      const IntRange *mspDisplayRange = mspGetFullDisplayRange(msp, seqSelected, bc);
      result = rangesOverlap(mspDisplayRange, displayRange);
    }

  return result;
}


/* Returns true if sequences in the given group should be shown. */
static gboolean isGroupVisible(const SequenceGroup* const group)
{
  gboolean result = TRUE;
  
  result = (!group || !group->hidden);
  
  return result;
}


/* Filter function. Returns true if the given row in the given tree model should be visible. */
static gboolean isTreeRowVisible(GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  gboolean bDisplay = FALSE;
  
  GList *mspList = treeGetMsps(model, iter);
  
  if (g_list_length(mspList) > 0)
    {
      GtkWidget *tree = GTK_WIDGET(data);
      TreeProperties *properties = treeGetProperties(tree);
      DetailViewProperties *dvProperties = detailViewGetProperties(properties->detailView);
      BlxViewContext *bc = blxWindowGetContext(dvProperties->blxWindow);

      /* Check the first msp to see if this sequence is in a group that's hidden.
       * (Note that all MSPs in the same row should be in the same sequence - we 
       * don't check this here because this function is called many times so we
       * avoid any unnecessary checks.) */
      const MSP *firstMsp = (const MSP*)(mspList->data);
      SequenceGroup *group = blxContextGetSequenceGroup(bc, firstMsp->sSequence);
      
      if (isGroupVisible(group))
	{
	  BlxViewContext *bc = treeGetContext(tree);
          GtkWidget *detailView = treeGetDetailView(tree);
          DetailViewProperties *dvProperties = detailViewGetProperties(detailView);

	  const int frame = properties->readingFrame;
	  const IntRange* const displayRange = &dvProperties->displayRange;
          const gboolean seqSelected = blxContextIsSeqSelected(bc, firstMsp->sSequence);

	  /* Show the row if any MSP in the list is an exon or blast match in the correct frame/strand
	   * and within the display range */
	  GList *mspListItem = mspList;
	
	  for ( ; mspListItem; mspListItem = mspListItem->next)
	    {
	      const MSP* msp = (const MSP*)(mspListItem->data);
	      
              if (isMspVisible(msp, bc, frame, displayRange, dvProperties->numUnalignedBases, seqSelected, group))
		{
		  bDisplay = TRUE;
		  break;
		}
	    }	  
	}
    }
  
  return bDisplay;
}


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

TreeProperties* treeGetProperties(GtkWidget *widget)
{
  return widget ? (TreeProperties*)(g_object_get_data(G_OBJECT(widget), "TreeProperties")) : NULL;
}

BlxViewContext* treeGetContext(GtkWidget *tree)
{
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  return blxWindowGetContext(blxWindow);
}

static void onDestroyTree(GtkWidget *widget)
{
  TreeProperties *properties = treeGetProperties(widget);
  
  if (properties)
    {
      if (properties->treeColumnHeaderList)
	{
	  g_list_free(properties->treeColumnHeaderList);
	  properties->treeColumnHeaderList = NULL;
	}
    
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "TreeProperties", NULL);
    }
}

static void treeCreateProperties(GtkWidget *widget, 
				 GtkWidget *grid, 
				 GtkWidget *detailView, 
				 const int frame,
				 GtkWidget *treeHeader,
				 GList *treeColumnHeaderList,
				 const gboolean hasSnpHeader)
{
  if (widget)
    { 
      TreeProperties *properties = (TreeProperties*)g_malloc(sizeof *properties);
      
      properties->grid = grid;
      properties->detailView = detailView;
      properties->readingFrame = frame;
      properties->treeHeader = treeHeader;
      properties->treeColumnHeaderList = treeColumnHeaderList;
      properties->hasSnpHeader = hasSnpHeader;
      
      int i = 0;
      for ( ; i < BLXMODEL_NUM_MODELS; ++i)
        properties->treeModels[i] = NULL;

      g_object_set_data(G_OBJECT(widget), "TreeProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyTree), NULL); 
    }
}


/* returns true if display coords should be negated */
static gboolean treeGetNegateCoords(GtkWidget *tree)
{
  /* We negate coords (literally just stick a '-' on the front) for display purposes if the
   * display is reversed and the negate-coords option is enabled. This gives the effect that coords
   * always increase left-to-right, whereas when the display is reversed they really decrease. */
  return blxWindowGetNegateCoords(treeGetBlxWindow(tree));
}


/***********************************************************
 *                       Tree events                       *
 ***********************************************************/

/* Handler for when vertical scrollbar for the tree has changed (either value or range) */
static void onScrollChangedTree(GtkObject *object, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  /* Remove any message from the detail-view status bar, because this shows info for the
   * currently-moused-over item, which may have changed now we've scrolled. */
   GtkStatusbar *statusBar = treeGetStatusBar(tree);
   
   if (statusBar)
     {
       guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
       gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
     }
}

/* Find the position of the left edge of the source widget wrt the dest widget and offset 
 * the event by this amount, then propagate the event to the dest widget. This allows us
 * to get the correct x coord on the parent widget for a click on the child widget. */
static void propagateEventButton(GtkWidget *srcWidget, GtkWidget *destWidget, GdkEventButton *event)
{
  int xOffset;
  gtk_widget_translate_coordinates(srcWidget, destWidget, 0, 0, &xOffset, NULL);
  
  event->x += xOffset;
  
  gtk_propagate_event(destWidget, (GdkEvent*)event);
}


/* Utility to translate the coords in the given event from coords relative to the source
 * widget to coords relative to the destination widget, and then propagate the event to the
 * dest widget. */
static void propagateEventMotion(GtkWidget *srcWidget, GtkWidget *destWidget, GdkEventMotion *event)
{
  /* Find the position of 0,0 of the source widget wrt the dest widget and offset by this amount. */
  int xOffset, yOffset;
  gtk_widget_translate_coordinates(srcWidget, destWidget, 0, 0, &xOffset, &yOffset);
  
  event->x += xOffset;
  event->y += yOffset;  
  
  gtk_propagate_event(destWidget, (GdkEvent*)event);
}


/* Draw the part of the tree header that shows the reference sequence (either as
 * nucleotide sequence, or a peptide sequence if viewing protein matches) */
static void drawRefSeqHeader(GtkWidget *headerWidget, GtkWidget *tree)
{
  GdkDrawable *drawable = createBlankPixmap(headerWidget);
  
  BlxViewContext *bc = treeGetContext(tree);
  const BlxStrand strand = treeGetStrand(tree);
  const int frame = treeGetFrame(tree);
  
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  
  const gboolean highlightSnps = treeHasSnpHeader(tree) && bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS];
  
  /* Find the segment of the ref seq to display. Ref seq is in nucleotide coords so convert the 
   * display range to nucleotide coords. */
  const int qIdx1 = convertDisplayIdxToDnaIdx(properties->displayRange.min, bc->seqType, frame, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  const int qIdx2 = convertDisplayIdxToDnaIdx(properties->displayRange.max, bc->seqType, frame, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange); 
  IntRange qRange = {min(qIdx1, qIdx2), max(qIdx1, qIdx2)};

  /* The q range may be outside the ref seq range if we are at the start/end and we have included
   * "missing" bases to make up complete codons. Adjust to within the range, but maintain the same
   * base number within the reading frame. */
  int offsetMin = 0;
  int offsetMax = 0;
  
  while (qRange.min < bc->refSeqRange.min)
    {
      qRange.min += bc->numFrames;
      ++offsetMin;
    }

  while (qRange.max > bc->refSeqRange.max)
    {
      qRange.max -= bc->numFrames;
      ++offsetMax;
    }
  
  GError *error = NULL;
  
  gchar *segmentToDisplay = getSequenceSegment(bc->refSeq,
                                               &qRange,
					       strand, 
					       BLXSEQ_DNA,      /* input ref seq is always in nucleotide coords */
                                               bc->seqType,     /* required segment is in display coords */
					       frame, 
					       bc->numFrames,
					       &bc->refSeqRange,
					       bc->blastMode,
					       bc->geneticCode,
					       bc->displayRev,
					       bc->displayRev,	/* show backwards if display reversed */
					       TRUE,		/* always complement reverse strand */
					       &error);
  
  if (!segmentToDisplay)
    {
      g_assert(error);
      prefixError(error, "Could not draw reference sequence header. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
    }
  
  GdkGC *gc = gdk_gc_new(drawable);

  /* Offset the x coord where we'll start drawing if we did not start at the beginning of the display range. */
  const int offset = bc->displayRev ? offsetMax : offsetMin;
  gdouble xStart = (gdouble)offset * properties->charWidth;
  const int yStart = 0;
  
  /* Find out if there are any special bases that need highlighting. */
  GHashTable *basesToHighlight = getRefSeqBasesToHighlight(detailView, &qRange, bc->seqType, strand);

  const int incrementValue = bc->displayRev ? -1 * bc->numFrames : bc->numFrames;
  int displayIdx = properties->displayRange.min;
  DrawBaseData baseData = {qIdx1, 0, strand, frame, bc->seqType, FALSE, FALSE, FALSE, TRUE, highlightSnps, FALSE, BLXCOLOR_REF_SEQ, NULL, NULL, FALSE, FALSE, FALSE, FALSE};

  while (displayIdx >= properties->displayRange.min && displayIdx <= properties->displayRange.max)
    {
      baseData.displayIdxSelected = (displayIdx == properties->selectedBaseIdx - offset);
      baseData.dnaIdxSelected = baseData.displayIdxSelected;
      baseData.baseChar = segmentToDisplay[displayIdx - properties->displayRange.min];
      
      const int x = (int)((gdouble)xStart + (gdouble)(displayIdx - properties->displayRange.min) * properties->charWidth);

      /* Draw the character, seting the background color and outline depending on whether this base is selected or
       * is affected by a SNP or polyA signal etc. */
      drawHeaderChar(bc, properties, drawable, gc, x, yStart, basesToHighlight, &baseData);
     
      baseData.dnaIdx += incrementValue;
      ++displayIdx;
    }
  
  /* Mark up the text to highlight the selected base, if there is one */
  PangoLayout *layout = gtk_widget_create_pango_layout(detailView, segmentToDisplay);
  pango_layout_set_font_description(layout, detailViewGetFontDesc(detailView));
  
  if (layout)
    {
      gtk_paint_layout(headerWidget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, detailView, NULL, xStart, yStart, layout);
      g_object_unref(layout);
    }
  
  drawColumnSeparatorLine(headerWidget, drawable, gc, bc);
 
  g_hash_table_unref(basesToHighlight);
  g_free(segmentToDisplay);
  g_object_unref(gc);
}



/* expose function to push a cached bitmap to screen */
static gboolean onExposeRefSeqHeader(GtkWidget *headerWidget, GdkEventExpose *event, gpointer data)
{
  GdkWindow *window = GTK_IS_LAYOUT(headerWidget) ? GTK_LAYOUT(headerWidget)->bin_window : headerWidget->window;
  
  if (window)
    {
      /* See if there's a cached drawable and, if not, create it */
      GdkDrawable *bitmap = widgetGetDrawable(headerWidget);
      
      if (!bitmap)
	{
	  /* There isn't a bitmap yet. Create it now. */
	  GtkWidget *tree = GTK_WIDGET(data);
	  drawRefSeqHeader(headerWidget, tree);
	  bitmap = widgetGetDrawable(headerWidget);
	}
      
      if (bitmap)
	{
	  /* Push the bitmap onto the window */
	  GdkGC *gc = gdk_gc_new(window);
	  gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
	}
    }
  
  return FALSE;
}


/* The given renderer is an MSP. This function checks if there is a base index
 * selected and, if so, colors the background for that base with the given color. */
static void treeHighlightSelectedBase(GtkWidget *tree, GdkDrawable *drawable)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  GList *columnList = detailViewGetColumnList(detailView);
  
  if (properties->selectedBaseIdx != UNSET_INT && valueWithinRange(properties->selectedBaseIdx, &properties->displayRange))
    {
      /* Convert the display-range index to a 0-based index in the display range */
      const int charIdx = properties->selectedBaseIdx - properties->displayRange.min;
      
      /* Get the x coords for the start and end of the sequence column */
      IntRange xRange;
      getColumnXCoords(columnList, BLXCOL_SEQUENCE, &xRange);
      
      const int x = xRange.min + (charIdx * properties->charWidth);
      const int y = 0;
      
      BlxViewContext *bc = blxWindowGetContext(properties->blxWindow);
      GdkColor *color = getGdkColor(BLXCOLOR_SELECTION, bc->defaultColors, FALSE, bc->usePrintColors);
      
      drawRect(drawable, color, x, y, roundNearest(properties->charWidth), tree->allocation.height, 0.3, CAIRO_OPERATOR_XOR);
    }
}


/* Expose function for a detail-view tree 
 * 
 * There is a bit of hacky code in here to try to speed up the redraw of the
 * tree, which can be very slow on some systems (I'm not sure why but it seems
 * quite a common problem with GtkTreeView, sadly).
 * 
 * The way the drawing works is as follows:
 *  - First time round, there is no cached drawable, so it creates the drawable
 *    and saves it in the widget. It is blank at this point. The default handler
 *    then continues.
 *  - The default handler gets the cell renderer to draw each row. The cell
 *    renderer draws to the cached drawable as well as the window.
 *  - Even though we have a cached drawable, we always let the default handler
 *    do the drawing and overwrite the cached drawable (see notes below) ...
 *  - ... EXCEPT when middle-dragging. The default handler would have to re-draw
 *    everything when middle-dragging, because the highlighted column on
 *    every row changes, and this can be slow. We therefore use the cached 
 *    drawable instead and draw the highlighted column over the top of it.
 *
 * Notes:
 *  - Ideally on subsequent calls we would just push the cached drawable to screen,
 *    (unless it has been cleared to force a re-draw). However, this does not
 *    work well for vertical scrolling, because it slows things down if we clear 
 *    and re-draw everything (I think the renderer normally just draws the
 *    relevant rows, so is quicker). This is why we let the default handler do 
 *    the drawing in most cases.  (We could perhaps implement some cleverer
 *    caching for vertical scrolling. Ideally caching needs to be done by the
 *    renderer, not by the tree.)
 *  - It should be safe to use the cached drawable while middle-dragging because the
 *    user should not be doing any other operations that will change what is
 *    shown in the trees (although we might need to introduce some blocks to avoid 
 *    other inputs happening by accident). We make sure the cached drawable is 
 *    up to date by forcing a re-draw when the user clicks the mouse, and then
 *    set the "mouse drag mode" flag to indicate that we want to use it (and
 *    to also stop it being overwritten).
 *  - The previously-highlighted column will still be highlighted for the duration 
 *    of a drag, but the highlighting is in a slightly different color, so this is
 *    actually quite a useful effect because you can see where you started dragging
 *    from.
 * */
static gboolean onExposeDetailViewTree(GtkWidget *tree, GdkEventExpose *event, gpointer data)
{
  gboolean handled = FALSE;
  GdkDrawable *drawable = widgetGetDrawable(tree);

  if (g_mouse_drag_mode && drawable)
    {
      /* Push the cached drawable to the window */
      GdkGC *gc = gdk_gc_new(drawable);
      GdkWindow *window = gtk_tree_view_get_bin_window(GTK_TREE_VIEW(tree));
      gdk_draw_drawable(window, gc, drawable, 0, 0, 0, 0, -1, -1);
      g_object_unref(gc);
      
      treeHighlightSelectedBase(tree, window);

      handled = TRUE;
    }
  else
    {
      /* Create a new drawable to draw to. Our custom cell renderer will draw 
       * to this as well as the widget's window. The cached drawable is used
       * for printing and for mouse-drag mode. */
      GdkDrawable *drawable = gdk_pixmap_new(tree->window, tree->allocation.width, tree->allocation.height, -1);
      gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
      widgetSetDrawable(tree, drawable);

      /* Draw a blank rectangle of the required widget background color */
      GdkGC *gc = gdk_gc_new(drawable);
      
      GdkColor *bgColor = tree->style->bg;
      gdk_gc_set_foreground(gc, bgColor);
      gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, tree->allocation.width, tree->allocation.height);
      
      g_object_unref(gc);
      
      /* Let the default handler continue to do the actual drawing */
      handled = FALSE;
    }
  
  return handled;
}


/* Select the range of rows in the given tree from the first path to the last inclusive */
static void treeSelectRowRange(GtkWidget *blxWindow, GtkTreeModel *model, GtkTreePath *firstPath, GtkTreePath *lastPath)
{
  /* We currently only ever have a list, so this function assumes the model is a list. If 
   * the alignment lists are ever changed to be trees, this function will need updating. */
  if (GTK_IS_TREE_STORE(model))
    {
      g_warning("Error selecting range of rows: function not implemented (expected alignments to be in a list store, not a tree store).");
      return;
    }
  
  /* Loop through all the rows from the last selected one to the one that was clicked
   * on. The direction depends on which one appears first in the list. */
  gint *firstIdx = gtk_tree_path_get_indices(firstPath);
  gint *lastIdx = gtk_tree_path_get_indices(lastPath);
  
  gint currentIdx = *firstIdx;
  const gint incrementValue = *firstIdx < *lastIdx ? 1 : -1;
  gboolean done = FALSE;
  
  while (!done)
    {
      /* Select the sequence in the current row */
      GtkTreePath *currentPath = gtk_tree_path_new_from_indices(currentIdx, -1);
      
      GtkTreeIter currentIter;
      if (gtk_tree_model_get_iter(model, &currentIter, currentPath))
	{
	  BlxSequence *currentSeq = treeGetSequence(model, &currentIter);
	  blxWindowSelectSeq(blxWindow, currentSeq);
	}
      else
	{
	  g_warning("Invalid iterator found when selecting a range of rows\n");
	  done = TRUE;
	}
      
      /* Get the next path, unless we're already at the last one */
      if (currentIdx == *lastIdx)
	{
	  done = TRUE;
	}
      else
	{
	  currentIdx += incrementValue;
	}

      gtk_tree_path_free(currentPath);
    }
}


static gboolean treeSelectRow(GtkWidget *tree, GdkEventButton *event)
{
  guint modifiers = gtk_accelerator_get_default_mod_mask();
  const gboolean ctrlModifier = ((event->state & modifiers) == GDK_CONTROL_MASK);
  const gboolean shiftModifier = ((event->state & modifiers) == GDK_SHIFT_MASK);
  
  if (!ctrlModifier && !shiftModifier)
    {
      blxWindowDeselectAllSeqs(treeGetBlxWindow(tree));
    }
  
  /* If we've clicked a different tree to the previously selected one, make this tree's 
   * reading frame the active one */
  GtkWidget *detailView = treeGetDetailView(tree);
  detailViewSetActiveFrame(detailView, treeGetFrame(tree));
  
  /* Find which row was clicked */
  GtkTreePath *clickedPath = NULL;
  GtkTreeViewColumn *clickedCol = NULL;
  int cell_x, cell_y;
  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), event->x, event->y, &clickedPath, &clickedCol, &cell_x, &cell_y);

  if (clickedPath)
    {
      GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
      GtkTreeIter clickedIter;
      gtk_tree_model_get_iter(model, &clickedIter, clickedPath);
      
      /* Get the sequence of the MSP(s) in this row */
      BlxSequence *clickedSeq = treeGetSequence(model, &clickedIter);

      if (clickedSeq)
	{
	  GtkWidget *blxWindow = treeGetBlxWindow(tree);
	  
	  if (!ctrlModifier && !shiftModifier)
	    {
	      /* No modifiers: select the sequence */
	      blxWindowSelectSeq(blxWindow, clickedSeq);
	    }
	  else if (ctrlModifier && !shiftModifier)
	    {
	      /* Ctrl pressed: toggle the selection state of the sequence */
	      blxWindowSetSeqSelected(blxWindow, clickedSeq, !blxWindowIsSeqSelected(blxWindow, clickedSeq));
	    }
	  else if (!ctrlModifier && shiftModifier)
	    {
	      /* Shift pressed: select all rows between the last-selected sequence and the clicked row */
	      const BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
	      GList *lastSelectedRows = treeGetSequenceRows(tree, lastSelectedSeq);
	      
	      if (g_list_length(lastSelectedRows) > 0) /* do nothing if no rows were previously selected */
		{
		  /* Get the tree row for the last selected sequence. (If there are multiple rows for this
		   * sequence, the desired behaviour here is a bit ambiguous, so for now just use any of them */
		  GtkTreePath *lastSelectedPath = (GtkTreePath*)(lastSelectedRows->data);
		  treeSelectRowRange(blxWindow, model, lastSelectedPath, clickedPath);
		}

              destroyTreePathList(&lastSelectedRows);
	    }
	}

      GtkWidget *detailView = treeGetDetailView(tree);
      detailViewSetSelectedStrand(detailView, treeGetStrand(tree));
      detailViewRedrawAll(detailView);

      gtk_tree_path_free(clickedPath);
    }
  
  return TRUE; /* handled */
}


static gboolean treePfetchRow(GtkWidget *tree)
{
  /* Get the selected sequence (assumes that only one is selected) */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (selectedSeqs)
    {
      const BlxSequence *clickedSeq = (const BlxSequence*)selectedSeqs->data;
      fetchSequence(clickedSeq, TRUE, 0, blxWindow, NULL, NULL);
    }

  return TRUE;
}


static gboolean onButtonPressTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
    {
    case 1:
      {
	/* Left button: selects a row, or shows the pfetch window if this was a double-click. */
	if (event->type == GDK_BUTTON_PRESS)
	  {
	    handled = treeSelectRow(tree, event);
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
	    handled = treePfetchRow(tree);
	  }

	break;
      }
      
    default:
      break;
    };
  
  if (!handled)
    { 
      propagateEventButton(tree, treeGetDetailView(tree), event);
    }
  
  return handled;
}


static gboolean onButtonPressTreeHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  GtkWidget *tree = GTK_WIDGET(data);

  switch (event->button)
    {
      case 1:
      {
	GtkWidget *detailView = treeGetDetailView(tree);
	
	if (event->type == GDK_BUTTON_PRESS)
	  {
            /* Set the active frame to the tree that was clicked on */
            detailViewSetActiveFrame(detailView, treeGetFrame(tree));
            
	    /* If variations are highlighted, select the variation that was clicked on, if any.  */
	    BlxViewContext *bc = treeGetContext(tree);

            if (bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS])
              {
                blxWindowDeselectAllSeqs(detailViewGetBlxWindow(detailView));
                int clickedBase = 1; /* only get in here for DNA matches; so frame/base number is always one */

                selectClickedSnp(header, NULL, detailView, event->x, event->y, FALSE, clickedBase); /* SNPs are always un-expanded in the DNA track */
	    
                detailViewRefreshAllHeaders(detailView);
              }
	  }
	else if (event->type == GDK_2BUTTON_PRESS)
	  {
            /* Double click. Show/hide the variations track (and highlight variations in the ref
             * sequence, if not already highighted). */
	    BlxViewContext *bc = treeGetContext(tree);
            const gboolean showTrack = !bc->flags[BLXFLAG_SHOW_VARIATION_TRACK];
	    bc->flags[BLXFLAG_SHOW_VARIATION_TRACK] = showTrack;
            
            if (showTrack)
              bc->flags[BLXFLAG_HIGHLIGHT_VARIATIONS] = TRUE;
            
            detailViewUpdateShowSnpTrack(detailView, showTrack);
	  }
	
	handled = TRUE;
	break;
      }
      
    default:
      break;
    };
  
  if (!handled)
    {
      /* Propagate to the detail view, translating event coords to detail view widget coords */
      GtkWidget *detailView = treeGetDetailView(tree);
      propagateEventButton(header, detailView, event);
    }
  
  return handled;
}


static BlxSequence* treeGetSequence(GtkTreeModel *model, GtkTreeIter *iter)
{
  BlxSequence *result = NULL;
  GList *mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      const MSP *firstMsp = (const MSP*)(mspGList->data);
      result = firstMsp->sSequence;
    }
  
  return result;
}


/* utility to destroy a GList of GtkTreePaths: destroys the data and the list
 * and sets the list pointer to null */
static void destroyTreePathList(GList **list)
{
  GList *item = *list;

  for ( ; item; item = item->next)
    {
      GtkTreePath *path = (GtkTreePath*)(item->data);
      gtk_tree_path_free(path);
    }

  g_list_free(*list);
  *list = NULL;
}


/* Get all the tree rows that contain MSPs from the given sequence. Returns
 * a GList of GtkTreePaths; the result must be freed by calling destroyTreePathList */
static GList *treeGetSequenceRows(GtkWidget *tree, const BlxSequence *clickedSeq)
{
  GList *resultList = NULL;
  
  /* Loop through all tree rows looking for any with a matching sequence */
  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
  
  GtkTreeIter iter;
  gboolean validIter = gtk_tree_model_get_iter_first(model, &iter);
  
  while (validIter)
    {
      /* Loop through all msps in this row and see if any belong
       * to the clicked sequence. (For normal matches we could just
       * check one msp because they should all belong to the same
       * BlxSequence, but for short reads we squash matches from 
       * different sequences onto the same row, so we can't assume
       * this) */
      GList *mspItem = treeGetMsps(model, &iter);
      
      for ( ; mspItem; mspItem = mspItem->next)
	{
	  const MSP *firstMsp = (const MSP*)(mspItem->data);
          
	  if (firstMsp->sSequence == clickedSeq)
	    {
	      GtkTreePath *path = gtk_tree_model_get_path(model, &iter);
	      resultList = g_list_append(resultList, path);
              break;
	    }
	}
      
      validIter = gtk_tree_model_iter_next(model, &iter);
    }

  return resultList;
}


/* Move the current row selection up/down by one row */
gboolean treeMoveRowSelection(GtkWidget *tree, const gboolean moveUp, const gboolean shiftModifier)
{
  /* Get the last selected sequence */
  GtkWidget *blxWindow = treeGetBlxWindow(tree);
  BlxSequence *lastSelectedSeq = blxWindowGetLastSelectedSeq(blxWindow);
  
  /* Get the row in this tree that contain MSPs from that sequence (if there are multiple 
   * rows behaviour is a bit ambiguous; for now just act on the first one that was found) */
  GList *rows = treeGetSequenceRows(tree, lastSelectedSeq);

  if (g_list_length(rows) > 0)
    {
      GtkTreePath *path = (GtkTreePath*)(rows->data);
      GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
      
      GtkTreeIter newRowIter;
      gboolean newRowExists = FALSE;
      
      /* We have to use a path to get the previous row and an iter to get the next row (because
       * gtk_tree_path_next is cyclic, so we can't tell if we're at the bottom already). */
      if (moveUp && gtk_tree_path_prev(path))
	{
	  newRowExists = gtk_tree_model_get_iter(model, &newRowIter, path);
	}
      else if (!moveUp)
	{
	  newRowExists = gtk_tree_model_get_iter(model, &newRowIter, path);
	  newRowExists &= gtk_tree_model_iter_next(model, &newRowIter);
	}
      
      if (newRowExists)
	{
	  BlxSequence *newSeq = treeGetSequence(model, &newRowIter);

	  if (shiftModifier && blxWindowIsSeqSelected(blxWindow, newSeq))
	    {
	      /* The new row is already selected; deselect the last-selected row */
	      blxWindowDeselectSeq(blxWindow, lastSelectedSeq);
	    }
	  else if (shiftModifier)
	    {
	      /* Add new row to current selection */
	      blxWindowSelectSeq(blxWindow, newSeq);
	    }
	  else
	    {
	      /* New row becomes sole selection */
	      blxWindowDeselectAllSeqs(treeGetBlxWindow(tree));
	      blxWindowSelectSeq(blxWindow, newSeq);
	    }
	  
	  treeScrollSelectionIntoView(tree, NULL);
	}
    }

  destroyTreePathList(&rows);
  
  return TRUE;
}


static gboolean onButtonReleaseTree(GtkWidget *tree, GdkEventButton *event, gpointer data)
{
  propagateEventButton(tree, treeGetDetailView(tree), event);
  return FALSE;
}


static gboolean onButtonReleaseTreeHeader(GtkWidget *header, GdkEventButton *event, gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  propagateEventButton(tree, treeGetDetailView(tree), event);
  return FALSE;
}


/* Implement custom scrolling for horizontal mouse wheel movements over the tree.
 * This scrolls our custom horizontal scrollbar for the whole Detail View. Leaves
 * vertical scrolling to the default handler. */
static gboolean onScrollTree(GtkWidget *tree, GdkEventScroll *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->direction)
    {
      case GDK_SCROLL_LEFT:
	{
	  scrollDetailViewLeftStep(treeGetDetailView(tree));
	  handled = TRUE;
	  break;
	}
	
      case GDK_SCROLL_RIGHT:
	{
	  scrollDetailViewRightStep(treeGetDetailView(tree));
	  handled = TRUE;
	  break;
	}

      default:
	{
	  /* Default handler can handle vertical scrolling */
	  handled = FALSE;
	  break;
	}
    };
  
  return handled;
}


static gboolean onMouseMoveTree(GtkWidget *tree, GdkEventMotion *event, gpointer data)
{
  if (!(event->state & GDK_BUTTON1_MASK) &&
      !(event->state & GDK_BUTTON2_MASK) &&
      !(event->state & GDK_BUTTON3_MASK) &&
      !(event->state & GDK_BUTTON4_MASK) &&
      !(event->state & GDK_BUTTON5_MASK))
    {
      /* Feed back details about the currently-hovered over row in the detail-view statusbar area */
      GtkStatusbar *statusBar = treeGetStatusBar(tree);
      
      if (statusBar)
        {
          /* Remove any previous message */
          guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
          gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
          
          /* Add details about the current row that is being hovered over */
          int bin_x = event->x, bin_y = event->y;
    //      gtk_tree_view_convert_widget_to_bin_window_coords(GTK_TREE_VIEW(tree), event->x, event->y, &bin_x, &bin_y); //only from GTK v2.12
          
          GtkTreePath *path = NULL;
          GtkTreeViewColumn *column = NULL;
          
          if (gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), bin_x, bin_y, &path, &column, NULL, NULL))
            {
              GtkTreeIter iter;
              GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
              gtk_tree_model_get_iter(model, &iter, path);

              GList *mspList = treeGetMsps(model, &iter);
              
              if (g_list_length(mspList) > 0)
                {
                  const MSP* const msp = (const MSP*)(mspList->data);
                  
                  /* Note: this assumes that all MSPs in this row are in the same BlxSequence.
                   * If there are more, it will just show data for the first BlxSequence found. */
                  if (msp->sSequence)
                    {
                      GtkWidget *detailView = treeGetDetailView(tree);
                      GList *columnList = detailViewGetColumnList(detailView);

                      char *displayText = blxSequenceGetSummaryInfo(msp->sSequence, columnList);
                  
                      if (displayText)
                        {
                          gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
                          g_free(displayText);
                        }
                    }
                }
            }

          gtk_tree_path_free(path);
        }
	
      return TRUE;
    }
  else
    {
      propagateEventMotion(tree, treeGetDetailView(tree), event);
      return FALSE;
    }
}

/* Get the display index at the given position in the tree header */
static int treeHeaderGetCoordAtPos(GtkWidget *header, GtkWidget *tree, const int x, const int y)
{
  int baseIdx = UNSET_INT;
  
  GtkWidget *detailView = treeGetDetailView(tree);
  GtkAdjustment *adjustment = detailViewGetAdjustment(detailView);

  if (x >= 0 && x <= header->allocation.width)
    {
      /* Get the 0-based char index at x */
      gdouble charWidth = detailViewGetCharWidth(detailView);
      int charIdx = (int)((gdouble)x / charWidth);
      
      /* Add the start of the scroll range to convert this to the display index */
      baseIdx = charIdx + adjustment->value;
    }
  else if (x < 0)
    {
      baseIdx = adjustment->value;
    }
  else if (x > header->allocation.width)
    {
      baseIdx = adjustment->value = adjustment->page_size;
    }
  
  return baseIdx;
}


static gboolean onMouseMoveTreeHeader(GtkWidget *header, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  GtkWidget *tree = GTK_WIDGET(data);

  if (event->state & GDK_BUTTON2_MASK)
    {
      /* Propagate the event to the detail view. The start of the header widget is at 
       * the start of the sequence column, so offset the coords before propagating. */
      GtkWidget *detailView = treeGetDetailView(tree);
      IntRange xRange;
      getColumnXCoords(detailViewGetColumnList(detailView), BLXCOL_SEQUENCE, &xRange);
      event->x += xRange.min;

      propagateEventMotion(tree, treeGetDetailView(tree), event);
      handled = FALSE;
    }
  else
    {
      /* If we're hovering over a base that's affected by a variation, then feed back info
       * about the variation to the user. Only applicable if we're showing a nucleotide sequence. */
      BlxViewContext *bc = treeGetContext(tree);
      
      if (bc->seqType == BLXSEQ_DNA)
        {
          /* Get the index we're hovering over */
          GtkWidget *detailView = treeGetDetailView(tree);
          const int displayIdx = treeHeaderGetCoordAtPos(header, tree, event->x, event->y);
          const int dnaIdx = convertDisplayIdxToDnaIdx(displayIdx, bc->seqType, detailViewGetActiveFrame(detailView), treeGetFrame(tree), bc->numFrames, bc->displayRev, &bc->refSeqRange);
          
          updateFeedbackAreaNucleotide(detailView, dnaIdx, treeGetStrand(tree));
        }
    }
      
  return handled;
}


//static void onDragBeginTree(GtkWidget *widget, GdkDragContext *event, gpointer data)
//{
//}
//
//static void onDragEndTree(GtkWidget *widget, GdkDragContext *event, gpointer data)
//{
//}
//
//static gboolean onDragMotionTree(GtkWidget *widget, GdkDragContext *event, gint x, gint y, guint time, gpointer data)
//{
//  return FALSE;
//}

/* Returns false => not in a drop zone */
static gboolean onDragDropTree(GtkWidget *widget,
                               GdkDragContext *drag_context,
                               gint            x,
                               gint            y,
                               guint           time,
                               gpointer        user_data)
{
  /* This call is required to stop a warning from gtk about drag-and-drop not
   * being supported for tree views */
  g_signal_stop_emission_by_name(widget, "drag-drop");
  return FALSE;
}

static gboolean onEnterTree(GtkWidget *tree, GdkEventCrossing *event, gpointer data)
{
  /* Return true to stop the default handler re-drawing when the focus changes */
  return TRUE;
}


/* Clear the detail-view feedback area */
static void clearStatusbar(GtkWidget *tree)
{
  GtkStatusbar *statusBar = treeGetStatusBar(tree);
  
  if (statusBar)
    {
      guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), DETAIL_VIEW_STATUSBAR_CONTEXT);
      gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
    }
}


static gboolean onLeaveTree(GtkWidget *tree, GdkEventCrossing *event, gpointer data)
{
  /* Remove any statusbar message that was added by mousing over the tree rows */
  clearStatusbar(tree);
  
  /* Return true to stop the default handler re-drawing when the focus changes */
  return TRUE;
}

static gboolean onLeaveTreeHeader(GtkWidget *header, GdkEventCrossing *event, gpointer data)
{
  /* Remove any statusbar message that was added by mousing over the tree header */
  GtkWidget *tree = GTK_WIDGET(data);
  clearStatusbar(tree);
  
  return TRUE;
}


/* Add a row to the given tree containing the given MSP */
void addMspToTree(MSP *msp, GtkWidget *tree, GtkListStore *store)
{
  if (tree)
    {
      GtkWidget *detailView = treeGetDetailView(tree);
      GList *columnList = detailViewGetColumnList(detailView);

      GtkTreeIter iter;
      gtk_list_store_append(store, &iter);
      
      /* The SequenceCellRenderer expects a GList of MSPs, so put our MSP in a list. For exons,
       * we want to add the child CDS/UTRs rather than the exon itself, so use the child list.
       * Note that this means we can have multiple MSPs on the same row even when the 'squash
       * matches' option is not on (which makes sense because realistically they are the same
       * object). */
      GList *mspGList = NULL;
      
      if (msp->type == BLXMSP_EXON && g_list_length(msp->childMsps) > 0)
        {
          mspGList = msp->childMsps;
        }
      else
        {
          mspGList = g_list_append(NULL, msp);
        }

      /* Loop through the rest of the columns */
      GList *item = columnList;
      
      for ( ; item; item = item->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);

          if (columnInfo->columnId == BLXCOL_SCORE)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, msp->score, -1);
          else if (columnInfo->columnId == BLXCOL_ID)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, msp->id, -1);
          else if (columnInfo->columnId == BLXCOL_START)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, msp->sRange.min, -1);
          else if (columnInfo->columnId == BLXCOL_SEQUENCE)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, mspGList, -1);
          else if (columnInfo->columnId == BLXCOL_END)
            gtk_list_store_set(store, &iter, columnInfo->columnIdx, msp->sRange.max, -1);
          else
            {
              GValue *val = blxSequenceGetValue(msp->sSequence, columnInfo->columnId);
              if (val) gtk_list_store_set_value(store, &iter, columnInfo->columnIdx, val);
            }
        }
      
      /* Remember the path to this tree row for each MSP */
      GtkTreePath *path = gtk_tree_model_get_path(GTK_TREE_MODEL(store), &iter);
      GList *mspItem = mspGList;
      
      for ( ; mspItem; mspItem = mspItem->next)
        {
          MSP *curMsp = (MSP*)(mspItem->data);
          curMsp->treePaths[BLXMODEL_NORMAL] = gtk_tree_path_to_string(path);
        }

      gtk_tree_path_free(path);
    }
}


/* Cell data function for the "name" column. This displays the sequence name, but
 * abbreviated to fit the cell, and with a symbol appended to indicate the s strand */
static void cellDataFunctionNameCol(GtkTreeViewColumn *column,
				     GtkCellRenderer *renderer, 
				     GtkTreeModel *model, 
				     GtkTreeIter *iter, 
				     gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);
  
  /* Get the MSP(s) for this row. They should all have the same sequence name. */
  GList	*mspGList = treeGetMsps(model, iter);
  const int numMsps = g_list_length(mspGList);

  if (GTK_WIDGET_VISIBLE(tree) && numMsps > 0)
    {
      MSP *msp = (MSP*)(mspGList->data);

      const int colWidth = gtk_tree_view_column_get_width(column);
      const gdouble charWidth = treeGetCharWidth(tree);
      const int maxLen = (int)((gdouble)colWidth / charWidth);

      if (maxLen > 2)
	{
	  /* Get the display name */
          const char *name = mspGetSName(msp);
          char *displayName = NULL;

          if (!name)
            name = "<no name>";
	
	  /* If the display is squashed and identical matches are on 
           * the same linke, we need to create a name that includes the
           * number of duplicate matches. */
          BlxViewContext *bc = treeGetContext(tree);

	  if (bc->modelId == BLXMODEL_SQUASHED && mspGetFlag(msp, MSPFLAG_SQUASH_IDENTICAL_FEATURES))
	    {
              char *name2 = g_strdup_printf(numMsps == 1 ? DUPLICATE_READS_COLUMN_NAME_SGL : DUPLICATE_READS_COLUMN_NAME, numMsps);
              displayName = abbreviateText(name2, maxLen - 2);
              g_free(name2);
            }
          
          if (!displayName)
            displayName = abbreviateText(name, maxLen - 2);
          
          if (displayName)
            {
              char displayText[maxLen + 1];
              sprintf(displayText, "%s", displayName);
	  
              int i = strlen(displayName);
              for ( ; i < maxLen - 1; ++i)
                {
                  displayText[i] = ' ';
                }

              displayText[maxLen - 1] = getStrandAsChar(msp && msp->sSequence ? msp->sSequence->strand : BLXSTRAND_NONE);
              displayText[maxLen] = 0;

              g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
              
              g_free(displayName);
            }
	}
      else
	{
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
	}
    }
}


/* Cell data function for the "start" column. This displays the start coord of the match 
 * sequence in normal left-to-right display, but the end coord if the display is reversed */
static void cellDataFunctionStartCol(GtkTreeViewColumn *column,
				     GtkCellRenderer *renderer, 
				     GtkTreeModel *model, 
				     GtkTreeIter *iter, 
				     gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);

  /* Get the MSPs in this row */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {  
      /* We want to display the min coord if we're in the same direction as the q strand,
       * or the max coord if we're in the opposite direction (unless the display is reversed,
       * in which case it's vice versa). */
      const MSP* const msp = (const MSP*)(mspGList->data);
      
      if (mspIsBlastMatch(msp))
        {
          const gboolean sameDirection = (treeGetStrand(tree) == mspGetMatchStrand(msp));
          const gboolean findMin = (treeGetDisplayRev(tree) != sameDirection);
          
          int coord = findMspListSExtent(mspGList, findMin);

          char displayText[numDigitsInInt(coord) + 1];
          sprintf(displayText, "%d", coord);
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
        }
      else
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
        }
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Cell data function for the "end" column. This displays the end coord of the match 
 * sequence in normal left-to-right display, but the start coord if the display is reversed */
static void cellDataFunctionEndCol(GtkTreeViewColumn *column, 
				   GtkCellRenderer *renderer, 
				   GtkTreeModel *model, 
				   GtkTreeIter *iter, 
				   gpointer data)
{
  GtkWidget *tree = GTK_WIDGET(data);

  /* Get the MSPs in this row */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {  
      /* We want to display the max coord if we're in the same direction as the q strand,
       * or the min coord if we're in the opposite direction (unless the display is reversed,
       * in which case it's vice versa). */
      const MSP* const msp = (const MSP*)(mspGList->data);
      
      if (mspIsBlastMatch(msp))
        {
          const gboolean sameDirection = (treeGetStrand(tree) == mspGetMatchStrand(msp));
          const gboolean findMin = (treeGetDisplayRev(tree) == sameDirection);
          
          int coord = findMspListSExtent(mspGList, findMin);

          char displayText[numDigitsInInt(coord) + 1];
          sprintf(displayText, "%d", coord);
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
        }
      else
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
        }
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Cell data function for the score column. */
static void cellDataFunctionScoreCol(GtkTreeViewColumn *column, 
				    GtkCellRenderer *renderer, 
				    GtkTreeModel *model, 
				    GtkTreeIter *iter, 
				    gpointer data)
{
  /* Get the MSP(s) for this row. Do not display coords if the row contains
   * multiple MSPs (unless they're short-reads, where we can assume that they
   * all have the same score because duplicate short-reads are grouped together). */
  GList	*mspGList = treeGetMsps(model, iter);

  if (g_list_length(mspGList) < 1)
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
  else
    {
      const MSP* const msp = (const MSP*)(mspGList->data);
    
      /* If the score is negative it means do not show */
      if (msp->score >= 0.0 && (g_list_length(mspGList) == 1 || mspGetFlag(msp, MSPFLAG_SQUASH_IDENTICAL_FEATURES)))
	{
	  const gdouble score = msp->score;
	  char displayText[numDigitsInInt((int)score) + 3]; /* +3 to include decimal point, 1 dp, and terminating nul */
	  
    //      sprintf(displayText, "%1.1f", score); 
	  sprintf(displayText, "%d", (int)score); 
	  
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
	}
      else
	{
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
	}
	
    }
}


/* Cell data function for the id column. */
static void cellDataFunctionIdCol(GtkTreeViewColumn *column, 
				  GtkCellRenderer *renderer, 
				  GtkTreeModel *model, 
				  GtkTreeIter *iter, 
				  gpointer data)
{
  /* Get the MSP(s) for this row. Do not display coords if the row contains multiple MSPs */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) < 1)
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
  else
    {
      const MSP* const msp = (const MSP*)(mspGList->data);

      /* If the score is negative it means do not show. Also only display the
       * ID if we only have one msp in the row (unless they are short reads, in 
       * which case the ID should be the same for all of them) */
      if (msp->id >= 0.0 && (g_list_length(mspGList) == 1 || mspGetFlag(msp, MSPFLAG_SQUASH_IDENTICAL_FEATURES)))
	{
	  const gdouble id = msp->id;
	  char displayText[numDigitsInInt((int)id) + 3]; /* +3 to include decimal point, 1 dp, and terminating nul */
	  
          sprintf(displayText, "%1.1f", id); 
//	  sprintf(displayText, "%d", (int)id); 

	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, displayText, NULL);
	}
      else
	{
	  g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
	}
    }
}


/* Cell data function for the Group column. */
static void cellDataFunctionGroupCol(GtkTreeViewColumn *column, 
                                     GtkCellRenderer *renderer, 
                                     GtkTreeModel *model, 
                                     GtkTreeIter *iter, 
                                     gpointer data)
{
  /* Get the MSP(s) for this row and find out they are in a group. All MSPs in a row should
   * be in the same sequence. */
  GList	*mspGList = treeGetMsps(model, iter);
  
  if (g_list_length(mspGList) > 0)
    {
      const MSP* const msp = (const MSP*)(mspGList->data);

      GtkWidget *tree = GTK_WIDGET(data);
      GtkWidget *blxWindow = treeGetBlxWindow(tree);

      SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, msp->sSequence);
      
      if (group)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, group->groupName, NULL);
        }
    }
}


/* Cell data function for the Organism column. */
static void cellDataFunctionOrganismCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  if (g_list_length(mspGList) > 0)
    {
      const MSP* const msp = (const MSP*)(mspGList->data);

      /* Use the abbreviation if available. */
      const char *text = mspGetOrganismAbbrev(msp);

      if (!text)
        {
          text = mspGetOrganism(msp);
        }
      
      if (text)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, text, NULL);
        }
    }
}


/* Cell data function for generic text columns. */
static void cellDataFunctionGenericCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  BlxColumnId columnId = (BlxColumnId)GPOINTER_TO_INT(data);
  GList	*mspGList = treeGetMsps(model, iter);

  if (g_list_length(mspGList) > 0)
    {
      const MSP* const msp = (const MSP*)(mspGList->data);
      const char *text = mspGetColumn(msp, columnId);
      
      if (text)
        {
          g_object_set(renderer, RENDERER_TEXT_PROPERTY, text, NULL);
        }
    }
}


static const char* mspListGetSource(GList *mspList)
{
  const char *result = NULL;
  GList *item = mspList;
  
  for ( ; item; item = item->next)
    {
      const MSP* const msp = (const MSP*)(item->data);
      const char *source = mspGetSource(msp);

      if (!source)
        {
          /* If any source is null, return null */
          result = NULL;
          break;
        }
      else if (!result)
        {
          /* First one found; set result */
          result = source;
        }
      else if (strcmp(source, result) != 0)
        {
          /* All sources are not the same. Return null. */
          result = NULL;
          break;
        }      
    }

  return result;
}


/* Cell data function for the Source column. */
static void cellDataFunctionSourceCol(GtkTreeViewColumn *column, GtkCellRenderer *renderer, GtkTreeModel *model, GtkTreeIter *iter, gpointer data)
{
  GList	*mspGList = treeGetMsps(model, iter);
  const char *source = mspListGetSource(mspGList);
  
  if (source)
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, source, NULL);
    }
  else
    {
      g_object_set(renderer, RENDERER_TEXT_PROPERTY, "", NULL);
    }
}


/* Callback called when the width of the sequence column has changed. */
static void onSeqColWidthChanged(GtkTreeViewColumn *column, GParamSpec *paramSpec, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  updateDynamicColumnWidths(detailView);
}


/* Create a single column in the tree. */
static GtkTreeViewColumn* createTreeColumn(GtkWidget *tree, 
                                           GtkWidget *detailView,
                                           GtkCellRenderer *renderer, 
                                           BlxColumnInfo *columnInfo,
                                           BlxColumnInfo *seqColInfo)
{
  /* Create the column in the tree */
  GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes(
    columnInfo->title, 
    renderer, 
    columnInfo->propertyName, columnInfo->columnIdx,  /* set the relevant property in the tree for the given column */
    RENDERER_DATA_PROPERTY, seqColInfo->columnIdx,  /* set the 'data' property in the tree to be the 'sequence' column (i.e. the MSP data) */
    NULL);

  /* Reduce the width of the end col by the scrollbar width. (This is so it matches the width
   * of the header, which does not have a scrollbar. ) */
  int width = columnInfo->width;
  if (columnInfo->columnId == BLXCOL_END)
    {
      width = width - scrollBarWidth();
    }
  
  if (columnInfo->columnId == BLXCOL_SEQUENCE)
    {
      g_signal_connect(G_OBJECT(column), "notify::width", G_CALLBACK(onSeqColWidthChanged), detailView);
    }
  
  /* Set the column properties and add the column to the tree */
  if (width > 0)
    {
      gboolean displayColumn = showColumn(columnInfo);
      gtk_tree_view_column_set_visible(column, displayColumn);
      gtk_tree_view_column_set_fixed_width(column, width);
    }
  else
    {
      /* Can't have 0 width, so hide the column instead */
      gtk_tree_view_column_set_visible(column, FALSE);
    }
  
  gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_FIXED);
  gtk_tree_view_column_set_resizable(column, TRUE);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tree), column);
  
  /* Special treatment for specific columns */
  switch (columnInfo->columnId)
  {
    case BLXCOL_SEQUENCE:
      gtk_tree_view_column_set_expand(column, TRUE);
      gtk_tree_view_column_set_sizing(column, GTK_TREE_VIEW_COLUMN_AUTOSIZE);
      break;
      
    case BLXCOL_SEQNAME:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionNameCol, tree, NULL);
      break;

    case BLXCOL_START:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionStartCol, tree, NULL);
      break;
      
    case BLXCOL_END:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionEndCol, tree, NULL);
      break;
      
    case BLXCOL_SCORE:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionScoreCol, tree, NULL);
      break;

    case BLXCOL_ID:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionIdCol, tree, NULL);
      break;

    case BLXCOL_GROUP:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionGroupCol, tree, NULL);
      break;
      
    case BLXCOL_SOURCE:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionSourceCol, tree, NULL);
      break;

    case BLXCOL_ORGANISM:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionOrganismCol, tree, NULL);
      break;
    
    default:
      gtk_tree_view_column_set_cell_data_func(column, renderer, cellDataFunctionGenericCol, GINT_TO_POINTER(columnInfo->columnId), NULL);
      break;
  }
  
  return column;
}


/* Refresh the name column header. This displays an abbreviated version of the
 * reference sequence name. It needs refreshing when the columns change size,
 * so we can re-abbreviate with more/less text as required. */
static void refreshNameColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *label = getLabelWidget(headerWidget);
  
  if (GTK_IS_LABEL(label))
    {
      GtkWidget *tree = GTK_WIDGET(data);

      /* Update the font, in case its size has changed */
      gtk_widget_modify_font(label , treeGetFontDesc(tree));

      TreeColumnHeaderInfo *headerInfo = treeColumnGetHeaderInfo(tree, BLXCOL_SEQNAME);
      const int colWidth = calculateColumnWidth(headerInfo, tree);

      /* Abbreviate the name */
      const char *refSeqName = blxWindowGetRefSeqName(treeGetBlxWindow(tree));
      const int maxLen = (int)((gdouble)colWidth / treeGetCharWidth(tree));
      
      char strandChar = (treeGetStrand(tree) == BLXSTRAND_FORWARD ? '+' : '-');
      char *stringToAppend = g_strdup_printf("(%c%d)", strandChar, treeGetFrame(tree));
      const int numCharsToAppend = strlen(stringToAppend);
      
      gchar *displayText = NULL;
      
      if (maxLen > numCharsToAppend)
	{
	  /* Abbreviate the name and then append the strand/frame */
	  gchar *displayName = abbreviateText(refSeqName, maxLen - numCharsToAppend);
	  displayText = g_strconcat(displayName, stringToAppend, NULL);
	  g_free(displayName);
	}
      else
	{
	  /* No space to concatenate the frame and strand. Just include whatever of the name we can */
	  displayText = abbreviateText(refSeqName, maxLen);
	}

      g_free(stringToAppend);

      if (displayText)
	{
	  gtk_label_set_text(GTK_LABEL(label), displayText);
	  g_free(displayText);
	}
    }
  else
    {
      g_warning("Unexpected widget type for Name column header; header may not refresh properly.\n");
    }
}


/* Refresh the start column header. This displays the start index of the current display range */
static void refreshStartColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *label = getLabelWidget(headerWidget);
  
  if (GTK_IS_LABEL(label))
    {
      GtkWidget *tree = GTK_WIDGET(data);
      BlxViewContext *bc = treeGetContext(tree);

      /* Update the font, in case its size has changed */
      gtk_widget_modify_font(label, treeGetFontDesc(tree));

      int displayVal = getStartDnaCoord(treeGetDisplayRange(tree), 
					treeGetFrame(tree),
					bc->seqType, 
					bc->displayRev, 
					bc->numFrames,
					&bc->refSeqRange);
      
      if (treeGetNegateCoords(tree))
        displayVal *= -1;

      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(label), displayText);
    }
  else
    {
      g_warning("Unexpected widget type for Start column header; header may not refresh properly.\n");
    }
}


/* Refresh the start column header. This displays the end index of the current display range */
static void refreshEndColHeader(GtkWidget *headerWidget, gpointer data)
{
  GtkWidget *label = getLabelWidget(headerWidget);

  if (GTK_IS_LABEL(label))
    {
      GtkWidget *tree = GTK_WIDGET(data);
      BlxViewContext *bc = treeGetContext(tree);
    
      int displayVal = getEndDnaCoord(treeGetDisplayRange(tree),
				      treeGetFrame(tree),
				      bc->seqType, 
				      bc->displayRev, 
				      bc->numFrames,
				      &bc->refSeqRange);
      
      if (treeGetNegateCoords(tree))
        displayVal *= -1;
      
      const int displayTextLen = numDigitsInInt(displayVal) + 1;
      
      gchar displayText[displayTextLen];
      sprintf(displayText, "%d", displayVal);
      displayText[displayTextLen - 1] = '\0';
      
      gtk_label_set_text(GTK_LABEL(label), displayText);
    }
  else
    {
      g_warning("Unexpected widget type for End column header; header may not refresh properly.\n");
    }
}


static int calculateColumnWidth(TreeColumnHeaderInfo *headerInfo, GtkWidget *tree)
{
  GtkWidget *detailView = treeGetDetailView(tree);
  GList *columnList = detailViewGetColumnList(detailView);
  
  /* Sum the width of all columns that this header includes */
  GList *listItem = headerInfo->columnIds;
  int width = 0;
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnId columnId = (BlxColumnId)GPOINTER_TO_INT(listItem->data);
      
      BlxColumnInfo *columnInfo = getColumnInfo(columnList, columnId);
      
      if (showColumn(columnInfo))
        {
          width = width + columnInfo->width;
        }
    }
  
  return width;
}


/* Utility function to create a TreeColumnHeaderInfo struct and initialise its values */
static TreeColumnHeaderInfo* createTreeColumnHeaderInfo(GtkWidget *headerWidget, 
							GtkWidget *tree, 
							GList *columnIds, 
							GtkCallback refreshFunc)
{
  TreeColumnHeaderInfo *headerInfo = (TreeColumnHeaderInfo*)g_malloc(sizeof(TreeColumnHeaderInfo));
  
  headerInfo->headerWidget = headerWidget;
  headerInfo->tree = tree;
  headerInfo->columnIds = columnIds;
  headerInfo->refreshFunc = refreshFunc;
  
  return headerInfo;
}


/* Create the column header widget for the given column in the tree header. It gets placed in
 * a TreecolumnHeaderInfo struct along with other header properties. The tree header shows
 * information about the reference sequence. */
static TreeColumnHeaderInfo* createTreeColHeader(GList **columnHeaders, 
                                                 GtkTreeViewColumn *treeColumn,
                                                 BlxColumnInfo *columnInfo,
                                                 TreeColumnHeaderInfo* firstTreeCol,
                                                 GtkWidget *headerBar,
                                                 GtkWidget *tree,
                                                 GtkWidget *detailView,
                                                 const char* const refSeqName,
                                                 const int frame,
                                                 const BlxStrand strand)
{
  /* Create a header widget for this column, if required. Also create a list of other 
   * column IDs we wish to merge under the same header (i.e. the 'Name' header spans over
   * the name column as well as the score and Id columns). */
  GtkWidget *columnHeader = NULL;
  GList *columnIds = NULL;
  GtkCallback refreshFunc = NULL;
  
  switch (columnInfo->columnId)
    {
      case BLXCOL_SEQNAME:
	{
	  /* The header above the name column will display the reference sequence name.
	   * This header will also span the score and id columns (and any optional columns)
           * seeing as we don't need to show any info anout the ref seq in those columns. */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE, TRUE);
	  refreshFunc = refreshNameColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
	  g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}
	
      case BLXCOL_SEQUENCE:
	{
	  /* The sequence column header contains the reference sequence. */
	  columnHeader = gtk_layout_new(NULL, NULL);
	  
	  seqColHeaderSetRow(columnHeader, frame);
	  gtk_widget_set_name(columnHeader, DNA_TRACK_HEADER_NAME);
	  g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeRefSeqHeader), tree);
	  
	  refreshFunc = refreshTextHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));

          gtk_widget_add_events(columnHeader, GDK_BUTTON_PRESS_MASK);
          gtk_widget_add_events(columnHeader, GDK_BUTTON_RELEASE_MASK);
          gtk_widget_add_events(columnHeader, GDK_POINTER_MOTION_MASK);
          gtk_widget_add_events(columnHeader, GDK_LEAVE_NOTIFY);
          g_signal_connect(G_OBJECT(columnHeader), "button-press-event",    G_CALLBACK(onButtonPressTreeHeader), tree);
          g_signal_connect(G_OBJECT(columnHeader), "button-release-event",  G_CALLBACK(onButtonReleaseTreeHeader), tree);
          g_signal_connect(G_OBJECT(columnHeader), "motion-notify-event",   G_CALLBACK(onMouseMoveTreeHeader), tree);
          g_signal_connect(G_OBJECT(columnHeader), "leave-notify-event",    G_CALLBACK(onLeaveTreeHeader), tree);
          
	  break;
	}

      case BLXCOL_START:
	{
	  /* The start column header displays the start index of the current display range */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE, TRUE);
	  refreshFunc = refreshStartColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
          g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}
	
      case BLXCOL_END:
	{
	  /* The end column header displays the start index of the current display range */
	  columnHeader = createLabel("", 0.0, 1.0, TRUE, TRUE, TRUE);
	  refreshFunc = refreshEndColHeader;
	  columnIds = g_list_append(columnIds, GINT_TO_POINTER(columnInfo->columnId));
          g_signal_connect(G_OBJECT(columnHeader), "expose-event", G_CALLBACK(onExposeGenericHeader), detailView);
	  break;
	}

      default:
        /* Any columns not specified above will not have their own header, so add them under the
         * first column's header instead. */
        {
          if (firstTreeCol)
            {
              firstTreeCol->columnIds = g_list_append(firstTreeCol->columnIds, GINT_TO_POINTER(columnInfo->columnId));
            }
          
          break;
        }
    };

  TreeColumnHeaderInfo *headerInfo = NULL;
  
  if (columnHeader)
    {
      /* Create a header info struct for each column, even if its contents are null */
      headerInfo = createTreeColumnHeaderInfo(columnHeader, tree, columnIds, refreshFunc);
      *columnHeaders = g_list_append(*columnHeaders, headerInfo);
    }
  
  return headerInfo;
}


///* Tooltip function for detail-view trees */
//static gboolean onQueryTooltip(GtkWidget *tree, gint x, gint y, gboolean keyboard_mode, GtkTooltip *tooltip, gpointer data)
//{
//  int binx, biny;
//  gtk_tree_view_convert_widget_to_bin_window_coords(GTK_TREE_VIEW(tree), x, y, &binx, &biny);
//  
//  GtkTreeViewColumn *column = NULL;
//  gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(tree), binx, biny, NULL, &column, NULL, NULL);
//  
//  if (column == gtk_tree_view_get_column(GTK_TREE_VIEW(tree), BLXCOL_SEQNAME))
//    {
//        gtk_tooltip_set_text(tooltip, "hello");
//
//    }
//  
//  return TRUE;
//}


/* Create the columns. Returns a list of header info for the column headers */
static GList* createTreeColumns(GtkWidget *tree, 
                                GtkWidget *detailView,
                                GtkCellRenderer *renderer, 
                                const BlxSeqType seqType,
                                GList *columnList,
                                GtkWidget *columnHeaderBar,
                                const char* const refSeqName,
                                const int frame,
                                const BlxStrand strand)
{
  /* We'll create a list of tree column headers */
  GList *treeColumns = NULL;

  /* The columns are defined by the columnList from the detail view. This list contains the
   * info such as the column width and title for each column. Loop through and create the tree
   * column header for each detail-view column. If any columns do not require their own header
   * in the tree, add them under the first column's header instead. */
  TreeColumnHeaderInfo* firstTreeCol = NULL;
  GList *column = columnList;
  BlxColumnInfo *seqColInfo = getColumnInfo(columnList, BLXCOL_SEQUENCE);
  
  for ( ; column; column = column->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)column->data;
      
      if (columnInfo)
	{
	  GtkTreeViewColumn *treeColumn = createTreeColumn(tree, detailView, renderer, columnInfo, seqColInfo);
          
	  TreeColumnHeaderInfo* headerInfo = createTreeColHeader(&treeColumns, 
                                                                 treeColumn, 
                                                                 columnInfo, 
                                                                 firstTreeCol, 
                                                                 columnHeaderBar, 
                                                                 tree, 
                                                                 detailView, 
                                                                 refSeqName, 
                                                                 frame, 
                                                                 strand);
          
          if (!firstTreeCol)
            {
              firstTreeCol = headerInfo;
            }
	}
      else
	{
	  g_warning("Error creating column; invalid column info in detail-view column-list.\n");
	}
    }
  
  /* Set a tooltip to display the sequence name. To do: at the moment this shows
   * the tooltip when hovering anywhere over the tree: really we probably just
   * want to show it when hovering over the name column. */
//  gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(tree), BLXCOL_SEQNAME); /* only in GTK 2.12 and higher */
//  gtk_widget_set_has_tooltip(tree, TRUE);
//  g_signal_connect(G_OBJECT(tree), "query-tooltip", G_CALLBACK(onQueryTooltip), NULL);

  return treeColumns;
}


/* Add the header widget for each column to the tree's header bar */
static void addColumnsToTreeHeader(GtkWidget *headerBar, GList *columnList)
{
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      TreeColumnHeaderInfo *columnInfo = (TreeColumnHeaderInfo*)(columnItem->data);
      
      /* Put the header widget in an event box so that we can color its background. */
      GtkWidget *eventBox = gtk_event_box_new();
      gtk_container_add(GTK_CONTAINER(eventBox), columnInfo->headerWidget);
      
      /* Put the event box into the header bar. If it's the sequence column, set the expand property. */
      gboolean expand = (g_list_find(columnInfo->columnIds, GINT_TO_POINTER(BLXCOL_SEQUENCE)) != NULL);
      
      gtk_box_pack_start(GTK_BOX(headerBar), eventBox, expand, TRUE, 0);
    }
}


/* This is the main sort comparison function for comparing two rows of a 
 * tree view. The sort criteria are specified in the detailView properties;
 * we may sort by multiple columns.
 * 
 * Returns a negative value if the first row appears before the second,
 * positive if the second appears before the first, or 0 if they are 
 * equivalent. */
static gint sortColumnCompareFunc(GtkTreeModel *model, GtkTreeIter *iter1, GtkTreeIter *iter2, gpointer data)
{
  gint result = UNSET_INT;

  GtkWidget *tree = GTK_WIDGET(data);
  GtkWidget *detailView = treeGetDetailView(tree);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);
  const int numColumns = g_list_length(columnList);

  /* Sort by each requested column in order of priority */
  if (dvProperties->sortColumns)
    {
      int priority = 0;

      for ( ; priority < numColumns; ++priority)
        {
          BlxColumnId sortColumn = dvProperties->sortColumns[priority];
          
          /* NONE indicates an unused entry in the priority array; if we reach
           * an unset value, there should be no more values after it */
          if (sortColumn == BLXCOL_NONE)
            break;

          /* Extract the MSPs for the two rows that we're comparing */
          GList *mspGList1 = treeGetMsps(model, iter1);
          GList *mspGList2 = treeGetMsps(model, iter2);
  
          /* Do the comparison on this column */
          result = sortByColumnCompareFunc(mspGList1, mspGList2, detailView, sortColumn);
          
          /* If rows are equal, continue to sort; otherwise we're done */
          if (result != 0)
            break;
        }
    }
  
  return result;
}


/* Create the base data store for a detail view tree */
void treeCreateBaseDataModel(GtkWidget *tree, gpointer data)
{
  /* Create the data store for the tree view (unless it already exists) */
  if (treeGetBaseDataModel(GTK_TREE_VIEW(tree)))
    return;
  
  /* Create the data store for the tree. We must know how many columns
   * we have, and their data types. */
  GtkWidget *detailView = treeGetDetailView(tree);  
  GList *columnList = detailViewGetColumnList(detailView);
  const int numCols = g_list_length(columnList);
  GType *typeList = columnListGetTypes(columnList);

  GtkListStore *store = gtk_list_store_newv(numCols, typeList);

  /* Set the sort function for each column */
  int colNum = 0;
  for ( ; colNum < numCols; ++colNum)
    {
      gtk_tree_sortable_set_sort_func(GTK_TREE_SORTABLE(store), colNum,	sortColumnCompareFunc, tree, NULL);
    }
  
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), GTK_TREE_MODEL(store));
  
  /* gtk_tree_view_set_model increments the reference count to 'store', so we should
   * decrement the reference count when we lose our local reference to it. */
  g_object_unref(G_OBJECT(store));  
}


/* Create a filtered version of the data store to only show rows that are within the display range. */
void treeCreateFilteredDataModel(GtkWidget *tree, gpointer data)
{
  GtkTreeModel *baseModel = treeGetBaseDataModel(GTK_TREE_VIEW(tree));

  /* If there is already a filtered model, there's nothing to do. (If the visible
   * model is different to the base model, we can assume it's the filtered model.)*/
  if (treeGetVisibleDataModel(GTK_TREE_VIEW(tree)) != baseModel)
    return;
  
  GtkTreeModel *filter = gtk_tree_model_filter_new(GTK_TREE_MODEL(baseModel), NULL);

  gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(filter), 
					 (GtkTreeModelFilterVisibleFunc)isTreeRowVisible, 
					 tree, 
					 NULL);
  
  /* Add the filtered store to the tree view */
  gtk_tree_view_set_model(GTK_TREE_VIEW(tree), GTK_TREE_MODEL(filter));
  
  /* Keep a reference to the model in the properties so we can switch between this and the 'squashed' model */
  TreeProperties *properties = treeGetProperties(tree);
  properties->treeModels[BLXMODEL_NORMAL] = gtk_tree_view_get_model(GTK_TREE_VIEW(tree));
  
  /* Note that we would normally decrement the reference count to 'filter' because we 
   * will lose the local pointer to it.  However, we've also added a pointer to it
   * from our properties, so we would need to increment the reference count again */
}


static void setTreeStyle(GtkTreeView *tree)
{
  gtk_widget_set_name(GTK_WIDGET(tree), DETAIL_VIEW_TREE_NAME);
  gtk_widget_set_redraw_on_allocate(GTK_WIDGET(tree), FALSE);
  
  gtk_tree_view_set_grid_lines(tree, GTK_TREE_VIEW_GRID_LINES_NONE);
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(tree), GTK_SELECTION_MULTIPLE);
  gtk_tree_view_set_reorderable(tree, TRUE);
  gtk_tree_view_set_headers_visible(tree, FALSE);
  gtk_tree_view_set_enable_search(tree, FALSE);
  
  /* Set the background color for the rows to be the same as the widget's background color */
  gtk_widget_modify_base(GTK_WIDGET(tree), GTK_STATE_NORMAL, GTK_WIDGET(tree)->style->bg);

  /* The default text color when rows are selected is white. This doesn't work 
   * well against our default background color of cyan, so use the same text color
   * as unselected rows. */
  gtk_widget_modify_text(GTK_WIDGET(tree), GTK_STATE_SELECTED, GTK_WIDGET(tree)->style->text);
  gtk_widget_modify_text(GTK_WIDGET(tree), GTK_STATE_ACTIVE, GTK_WIDGET(tree)->style->text);
  
  /* Set the expander size to 0 so that we can have tiny rows (otherwise the min is 12pt).
   * Also set the vertical separator to 0 so that we can have the option of the smallest
   * fonts possible. (The vertical separator causes gaps between rows. We want rows to be
   * flush, so we render over this gap. However, the bigger the vertical separator, the
   * bigger the row height (and hence font size) we must have to cover the gaps.) */
  char parseString[500];
  sprintf(parseString, "style \"packedTree\"\n"
	  "{\n"
	  "GtkTreeView::expander-size	      = 0\n"
	  "GtkTreeView::vertical-separator    = 0\n"
	  "GtkTreeView::horizontal-separator  = 0\n"
	  "}"
	  "widget \"*%s*\" style \"packedTree\"", DETAIL_VIEW_TREE_NAME);
  gtk_rc_parse_string(parseString);  
}


/* Create the widget that will contain all the header widgets. */
static GtkWidget *createDetailViewTreeHeader(GtkWidget *detailView, 
					     const gboolean includeSnpTrack,
					     const BlxStrand strand,
					     GtkWidget *columnHeaderBar)
{
  GtkWidget *treeHeader = NULL; /* the outermost container for the header widgets */
  
  if (includeSnpTrack)
    {
      /* Pack the snp track and the column header bar into a vbox */
      treeHeader = gtk_vbox_new(FALSE, 0);
      gtk_widget_set_name(treeHeader, HEADER_CONTAINER_NAME);
      
      createSnpTrackHeader(GTK_BOX(treeHeader), detailView, strand);
      gtk_box_pack_start(GTK_BOX(treeHeader), columnHeaderBar, FALSE, FALSE, 0);
    }
  else
    {
      /* Only include the column headers */
      treeHeader = columnHeaderBar;
    }
  
  return treeHeader;
}


GtkWidget* createDetailViewTree(GtkWidget *grid, 
				GtkWidget *detailView, 
				GtkCellRenderer *renderer,
				GList **treeList,
				GList *columnList,
				BlxSeqType seqType,
				const char* const refSeqName,
				const int frame,
				const gboolean includeSnpTrack)
{
  /* Find the strand the tree belongs to from the grid */
  const BlxStrand strand = gridGetStrand(grid);
  
  /* Create a tree view for the list of match sequences */
  GtkWidget *tree = gtk_tree_view_new();
  setTreeStyle(GTK_TREE_VIEW(tree));

  /* Put it in a scrolled window for vertical scrolling only (hoz scrolling will be via our
   * custom adjustment). Always display the scrollbars because we assume the column widths 
   * are the same for all trees and they won't be if one shows a scrollbar and another doesn't. */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_widget_set_name(scrollWin, DETAIL_VIEW_TREE_CONTAINER_NAME);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_NEVER, GTK_POLICY_ALWAYS);
  gtk_container_add(GTK_CONTAINER(scrollWin), tree);

  GtkAdjustment *adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrollWin));
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollChangedTree), tree);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollChangedTree), tree);

  /* Create a header, and put the tree and header in a vbox. This vbox is the outermost
   * container for the tree */
  GtkWidget *columnHeaderBar = gtk_hbox_new(FALSE, 0);
  GtkWidget *treeHeader = createDetailViewTreeHeader(detailView, includeSnpTrack, strand, columnHeaderBar);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_widget_set_name(vbox, DETAIL_VIEW_TREE_CONTAINER_NAME);
  gtk_box_pack_start(GTK_BOX(vbox), treeHeader, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), scrollWin, TRUE, TRUE, 0);
  
  /* Create the tree columns */
  GList *treeColumnHeaderList = createTreeColumns(tree, detailView, renderer, seqType, columnList, columnHeaderBar, refSeqName, frame, strand);
  
  /* Add the columns to the tree header */
  addColumnsToTreeHeader(columnHeaderBar, treeColumnHeaderList);
  
  /* Set the essential tree properties */
  treeCreateProperties(tree, grid, detailView, frame, treeHeader, treeColumnHeaderList, includeSnpTrack);
  
  /* Connect signals */
  gtk_widget_add_events(tree, GDK_FOCUS_CHANGE_MASK);
  g_signal_connect(G_OBJECT(tree), "button-press-event",    G_CALLBACK(onButtonPressTree),	NULL);
  g_signal_connect(G_OBJECT(tree), "button-release-event",  G_CALLBACK(onButtonReleaseTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "motion-notify-event",   G_CALLBACK(onMouseMoveTree),	detailView);
  g_signal_connect(G_OBJECT(tree), "scroll-event",	    G_CALLBACK(onScrollTree),		detailView);
  g_signal_connect(G_OBJECT(tree), "enter-notify-event",    G_CALLBACK(onEnterTree),		NULL);
  g_signal_connect(G_OBJECT(tree), "leave-notify-event",    G_CALLBACK(onLeaveTree),		NULL);
  //g_signal_connect(G_OBJECT(tree), "drag-begin",	    G_CALLBACK(onDragBeginTree),	NULL);
  //g_signal_connect(G_OBJECT(tree), "drag-end",		    G_CALLBACK(onDragEndTree),		NULL);
  //g_signal_connect(G_OBJECT(tree), "drag-motion",	    G_CALLBACK(onDragMotionTree),	NULL);
  g_signal_connect(G_OBJECT(tree), "drag-drop",		    G_CALLBACK(onDragDropTree),		NULL);
  g_signal_connect(G_OBJECT(tree), "expose-event",	    G_CALLBACK(onExposeDetailViewTree), NULL);

  GtkTreeSelection *selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
  g_signal_connect(G_OBJECT(selection), "changed", G_CALLBACK(onSelectionChangedTree), tree);
  
  /* Add the tree's outermost container to the given tree list. Also increase its ref count so that
   * we can add/remove it from its parent (which we do to switch panes when we toggle strands)
   * without worrying about it being destroyed. */
  *treeList = g_list_append(*treeList, vbox);
  g_object_ref(vbox);
  
  return vbox;
}

