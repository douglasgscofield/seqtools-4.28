/*  File: belvuTree.c
 *  Author: Gemma Barson, 2011-05-06
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
 * Description: See belvuTree.h
 *----------------------------------------------------------------------------
 */


#include "belvuApp/belvuTree.h"
#include "belvuApp/belvuWindow.h"
#include "seqtoolsUtils/utilities.h"
#include <gdk/gdkkeysyms.h>
#include <string.h>
#include <ctype.h> /* for isspace etc. */
#include <math.h>


#define DEFAULT_TREE_WINDOW_WIDTH_FRACTION      0.6    /* default width of tree window (as fraction of screen width) */
#define DEFAULT_TREE_WINDOW_HEIGHT_FRACTION     0.35   /* default height of tree window (as fraction of screen height) */
#define DEFAULT_XPAD                            10
#define DEFAULT_YPAD                            10
#define DIALOG_XPAD                             12      /* default x padding around dialog widgets */
#define DIALOG_YPAD                             8       /* default y padding around dialog widgets */
#define TABLE_XPAD                              12      /* default x padding around table elements */
#define TABLE_YPAD                              2       /* default y padding around table elements */
#define TITLE_MAIN_TREE_PREFIX                  ""                  /* prefix for tree window title for the main belvu tree */
#define TITLE_BOOTSTRAP_TREE_PREFIX             "Bootstrap "        /* prefix for tree window title for bootstrap trees */
#define TITLE_NJ_TREE_DESCRIPTION               "Neighbor-joining " /* window title tree description for NJ trees */
#define TITLE_UPGMA_TREE_DESCRIPTION            "UPGMA "            /* window title tree description for NJ trees */



/* Globals; original build methods in the tree dialog */
static BelvuBuildMethod origBuildMethod;
static BelvuDistCorr origDistCorr;


/* Utility struct to hold data for the build-changed signal */
typedef struct _BuildMethodChangedData
{
  BelvuContext *bc;
  GtkWidget **treeScaleEntry;
  GtkWidget *buildMethodCombo;
  GtkWidget *distCorrCombo;
} BuildMethodChangedData;


/*  BLOSUM62 930809
 
 #  Matrix made by matblas from blosum62.iij
 #  * column uses minimum score
 #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
 #  Blocks Database = /data/blocks_5.0/blocks.dat
 #  Cluster Percentage: >= 62
 #  Entropy =   0.6979, Expected =  -0.5209
 
 Note: to use with a2b[], always subtract 1 from the values !!!!
 
 A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
static int BLOSUM62[24][24] = {
{4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1}
};


/* ASCII-to-binary translation table
 Note: to use with BLOSUM62[], always subtract 1 from the values !!!! */
#undef NA
#define NA 23
static int a2b[] =
{
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA, 1,21, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,NA,19,NA,NA,NA,NA,NA,NA,
NA, 1,NA, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
15, 6, 2,16,17,NA,20,18,23,19,22,NA,NA,NA,NA,NA,

NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
};


/* This struct represents an area on the main tree drawing area
 * where a node is drawn. It can be used to find a node that was 
 * clicked on. */
typedef struct _ClickableRect
{
  GdkRectangle rect;                /* The area that this clickable rect covers on the drawing area */
  TreeNode *node;                   /* The node associated with this area */
  gboolean isBranch;                /* True if this is a tree branch line (false if it's the actual sequence name) */
} ClickableRect;


/* Properties specific to the belvu tree */
typedef struct _BelvuTreeProperties
  {
    BelvuContext *bc;	            /* The belvu context */
    GtkActionGroup *actionGroup;
    
    Tree *tree;                     /* The underlying tree struct */
    gboolean isMainTree;            /* True if tree is the main tree's root node. */
    
    GtkWidget *treeArea;            /* Drawing widget for the tree */
    GdkRectangle treeRect;          /* Specifies the actual area in which we'll draw the tree within treeAre */
    
    gdouble charWidth;
    gdouble charHeight;
    
    BelvuBuildMethod buildMethod;   /* The build method used to build the tree */
    BelvuDistCorr distCorr;         /* The distance-correction method used to build the tree */
    
    double treeScale;               /* The tree scale (not used if main tree - uses the values in the context instead) */
    double lineWidth;               /* Line width (not used if main tree) */
    gboolean showBranchLen;         /* Whether to display branch lengths (not used if main tree) */
    gboolean showOrganism;          /* Whether to show organism names (not used if main tree) */
    
    GArray *clickableRects;         /* Array of rectangles that associate clickable areas in treeArea to TreeNodes. */
  } BelvuTreeProperties;



/* Local function declarations */
static Tree*                        createEmptyTree();
static void                         calculateBelvuTreeBorders(GtkWidget *belvuTree);


/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuTreeProperties* belvuTreeGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuTreeProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuTreeProperties")) : NULL;
}

static void onDestroyBelvuTree(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);

  /* We must remove the tree from the list of spawned windows */
  properties->bc->spawnedWindows = g_slist_remove(properties->bc->spawnedWindows, belvuTree);
  
  /* If it is the main tree, also set the pointer in the context to null */
  if (properties->isMainTree)
    properties->bc->belvuTree = NULL;
  
  if (properties)
    {
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuTree), "BelvuTreeProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void belvuTreeCreateProperties(GtkWidget *belvuTree, 
                                      BelvuContext *bc,
                                      Tree *tree,
                                      const gboolean isMainTree,
                                      GtkActionGroup *actionGroup,
                                      GtkWidget *treeArea,
                                      BelvuBuildMethod buildMethod,
                                      BelvuDistCorr distCorr)
{
  if (belvuTree)
    {
      BelvuTreeProperties *properties = (BelvuTreeProperties*)g_malloc(sizeof *properties);
      
      properties->bc = bc;
      
      properties->tree = tree;
      properties->isMainTree = isMainTree;
      
      properties->treeArea = treeArea;
      properties->actionGroup = actionGroup;
      
      properties->buildMethod = buildMethod;
      properties->distCorr = distCorr;
      properties->treeScale = bc->treeScale;
      properties->lineWidth = bc->treeLineWidth;
      properties->showBranchLen = bc->treeShowBranchlen;
      properties->showOrganism = bc->treeShowOrganism;
      
      properties->charHeight = 0;
      properties->charWidth = 0;

      properties->clickableRects = g_array_new(FALSE, FALSE, sizeof(ClickableRect));
      
      g_object_set_data(G_OBJECT(belvuTree), "BelvuTreeProperties", properties);
      g_signal_connect(G_OBJECT(belvuTree), "destroy", G_CALLBACK (onDestroyBelvuTree), NULL);
    }
}


BelvuContext* belvuTreeGetContext(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  return properties ? properties->bc : NULL;
}


GtkActionGroup* belvuTreeGetActionGroup(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  return (properties ? properties->actionGroup : NULL);
}


/* Get the tree scale for this tree. Uses the value from the main context
 * if this is the main belvu tree, or the value stored in the properties 
 * for other trees. */
static double treeGetScale(BelvuTreeProperties *properties)
{
  return (properties->isMainTree ? properties->bc->treeScale : properties->treeScale);
}

static double treeGetLineWidth(BelvuTreeProperties *properties)
{
  return (properties->isMainTree ? properties->bc->treeLineWidth : properties->lineWidth);
}

static gboolean treeGetShowBranchLen(BelvuTreeProperties *properties)
{
  return (properties->isMainTree ? properties->bc->treeShowBranchlen : properties->showBranchLen);
}

static gboolean treeGetShowOrganism(BelvuTreeProperties *properties)
{
  return (properties->isMainTree ? properties->bc->treeShowOrganism : properties->showOrganism);
}


/***********************************************************
 *                   Tree bootstrapping                    *
 ***********************************************************/

static int BSorder(gconstpointer xIn, gconstpointer yIn)
{
  const BootstrapGroup *x = (const BootstrapGroup*)xIn;
  const BootstrapGroup *y = (const BootstrapGroup*)yIn;
  
  int result = 0;
  
  if (x && x->s && y && y->s)
    result = strcmp(x->s, y->s);
  else if (x && x->s)
    result = 1;
  else if (y && y->s)
    result = -1;
  else
    result = 0;
  
  return result;
}


static int BSptrorder(gconstpointer xIn, gconstpointer yIn)
{
  const BootstrapGroup *x = *((const BootstrapGroup**)xIn);
  const BootstrapGroup *y = *((const BootstrapGroup**)yIn);
  
  return BSorder(x, y);
}


static BootstrapGroup* createEmptyBootstrapGroup()
{
  BootstrapGroup *result = (BootstrapGroup*)g_malloc(sizeof *result);
  
  result->node = NULL;
  result->s = NULL;
  
  return result;
}


/* Combines left and right sequence groups and insert to bootstrap group list */
static GArray* fillBootstrapGroups(BelvuContext *bc, Tree *tree, TreeNode *node, const gboolean isMainTree) 
{
  GArray *result = NULL;
  
  BootstrapGroup *BS = createEmptyBootstrapGroup();
  
  if (!node->name)
    {
      /* Internal node */
      /* Combine left node sequences into right node array */
      
      GArray *left = fillBootstrapGroups(bc, tree, node->left, isMainTree);
      GArray *right = fillBootstrapGroups(bc, tree, node->right, isMainTree);
      
      /* Nothing to do for root */
      if (node != tree->head)
        {      
          /* Combine left and right groups */
          int i = 0;
          for (i = 0 ; i < (int)left->len; ++i) 
            {
              char *s = g_array_index(left, char*, i);
              g_array_append_val(right, s);
              g_array_sort(right, strcmp_);
            }
          
          g_array_unref(left);
          
          /* Create string with group members */
          int ssize = 1;
          for (i = 0 ; i < (int)right->len ; ++i) 
            ssize += (strlen(g_array_index(right, char*, i)) + 1);
          
          BS->s = (char*)g_malloc(ssize+1);
          BS->s[0] = 0;
          
          for (i = 0 ; i < (int)right->len ; ++i) 
            {
              char *ch = g_array_index(right, char*, i);
              strcat(BS->s, ch);
              strcat(BS->s, " ");
            }
          
          /* printf("   New bootstrap group: %s\n", BS->s); */
          
          if (isMainTree) 
            {
              /* Associate this node with the group string */
              BS->node =  node;
              
              /* Add group string to array of bootstrap groups */
              g_array_append_val(bc->bootstrapGroups, BS);
            }
          else
            {
              /* Find group string and increment counter if exists */
              BootstrapGroup *BS2;
              
              int ip = 0;
              if (bsArrayFind(bc->bootstrapGroups, BS, &ip, BSptrorder)) 
                {
                  BS2 = g_array_index(bc->bootstrapGroups, BootstrapGroup*, ip);
                  BS2->node->boot++;
                  DEBUG_OUT("Found bootgroup %s (%d)\n", BS->s, ip);
                }
              else
                {
                  DEBUG_OUT("Did not find bootgroup %s\n", BS->s);
                }
              
              if (BS && BS->s)
                g_free(BS->s);
              
              if (BS)
                g_free(BS);
            }
          
          result = right;
        }
    }
  else
    {
      /* Leaf node - return as embryonic array with one entry */
      GArray *leaf = g_array_sized_new(FALSE, FALSE, sizeof(char*), bc->alignArr->len);
      g_array_append_val(leaf, node->name);
      g_array_sort(leaf, strcmp_);
      result = leaf;
    }
  
  return result;
}


static void normaliseBootstraps(BelvuContext *bc, TreeNode *node) 
{
  node->boot = 100.0 * (double)node->boot/(double)bc->treebootstraps;
}


/* Bootstrap the internal nodes in a tree.
 1. Set up a list (array) of all internal nodes with leaf content and pointer to node
 2. In bootstrap tree, for each internal node, check if its contents exists in list. If so, increment node's bootstrap count
 3. Turn increments to percentages
 */
static void treeBootstrapStats(BelvuContext *bc, Tree *tree)
{
  /* Traverse tree, fill array bootstrapGroups */
  bc->bootstrapGroups = g_array_sized_new(FALSE, TRUE, sizeof(BootstrapGroup*), bc->alignArr->len);
  fillBootstrapGroups(bc, tree, tree->head, TRUE);
  g_array_sort(bc->bootstrapGroups, BSptrorder);

#ifdef DEBUG
  DEBUG_OUT("Created bootstrap groups:\n");
  int i = 0;
  for ( ; i < bc->bootstrapGroups->len; ++i)
    {
      BootstrapGroup *bs = g_array_index(bc->bootstrapGroups, BootstrapGroup*, i);
      DEBUG_OUT("%d:  %s\n", i, (bs && bs->s ? bs->s : "null"));
    }
#endif
  
  treeBootstrap(bc);
  
  /* Turn increments to percentages */
  treeTraverse(bc, tree->head, normaliseBootstraps);
}


void treeBootstrap(BelvuContext *bc)
{
  separateMarkupLines(bc);
  
  /* We will change the sequence strings in the alignments, so first
   * we must make a copy so we can revert them. */
  GArray *alignArrTmp = copyAlignArray(bc->alignArr);
  
#if defined(__CYGWIN__) || defined(DARWIN)
  srand(time(0));
#else
  srand48(time(0));
#endif
  
  int iter = 0;
  for (iter = 0; iter < bc->treebootstraps; ++iter) 
    {
      /* Make bootstrapped Alignment */
      int i = 0;
      for (i = 0; i < bc->maxLen; ++i) 
        {
          int src = 0;
          
#if defined(__CYGWIN__) || defined(DARWIN)
          src = (int)((((double)rand())/RAND_MAX) * bc->maxLen);
#else
          src = (int)(drand48() * bc->maxLen);
#endif
          
          columnCopy(bc->alignArr, i, alignArrTmp, src);
        }
      
      /* Create the bootstrap tree */
      Tree* tree = treeMake(bc, FALSE, FALSE);
      
      if (bc->outputBootstrapTrees) 
        {
          if (bc->treebootstrapsDisplay) 
            {
              /* The tree window takes ownership of the tree struct */
              createBelvuTreeWindow(bc, tree, FALSE);
            }
          else
            {
              /* Print the tree, then destroy it */
              g_message("\n");
              saveTreeNH(tree, tree->head, stdout);
              g_message(";\n");
              
              destroyTree(&tree);
            }
        }
      else
        {
          /* Collect the bootstrap statistics then destroy the tree */
          fillBootstrapGroups(bc, tree, tree->head, FALSE);
          destroyTree(&tree);
        }
    }
  
  /* Restore alignment */
  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *srcAln = g_array_index(alignArrTmp, ALN*, i);
      ALN *destAln = g_array_index(bc->alignArr, ALN*, i);
      
      if (destAln->sequenceStr)
        g_string_free(destAln->sequenceStr, TRUE);
      
      if (srcAln->sequenceStr)
        destAln->sequenceStr = srcAln->sequenceStr; /* take ownership of the sequence string */
      else 
        destAln->sequenceStr = NULL;
    }
  
  g_array_free(alignArrTmp, TRUE); /* destroys the array and all its ALN structs */
  
  reInsertMarkupLines(bc);
}


/***********************************************************
 *                   Business logic                        *
 ***********************************************************/

/*
 * Expect format:
 *
 * name1 name2 name3
 * 1-1   1-2   1-3
 * 2-1   2-2   2-3
 * 3-1   3-2   3-3
 */
static void treeReadDistances(BelvuContext *bc, double **pairmtx)
{
  char   *p;
  int    i = 0, j = 0 ;
  double  d;
  char line[MAXLENGTH+1];
  
  while (!feof (bc->treeReadDistancesPipe))
    { 
      if (!fgets (line, MAXLENGTH, bc->treeReadDistancesPipe))
	break;
      
      j = 0;
      
      while ( (p = strtok(j ? 0 : line, " \t")) && !isspace(*p))
	{
	  d = atof(p);
	  DEBUG_OUT("%d  %d  %f\n", i, j, d);
	  pairmtx[i][j] = d;
	  j++;
	}
      
      if (j != (int)bc->alignArr->len)
	g_error("nseq = %d, but read %d distances in row %d\n", bc->alignArr->len, j, i);
      
      i++;
    }
  
  if (j != (int)bc->alignArr->len)
    g_error("nseq = %d, but read %d distance rows\n", bc->alignArr->len, i);
  
  return ;
}


/* Correct an observed distance to an estimated 'true' distance.
 
 Adapted from Lasse Arvestad's lapd program
 
 od = observed distance in percent 
 
 */
static double treeJUKESCANTOR(double od)
{
  double cd;
  
  od /= 100;
  
  od = 20.0/19.0*od;
  if  (od > 0.95) od = 0.95; /* Limit to 300 PAM */
  
  cd = -19.0/20.0 * log(1-od);
  
  return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.
 
 Adapted from Lasse Arvestad's lapd program
 
 od = observed distance in percent 
 
 */
static double treeKimura(double od)
{
  double adjusted, cd;
  
  od /= 100;
  
  adjusted = od + 0.2 * od*od;
  
  if (adjusted > 0.95) adjusted = 0.95;   /* Limit to 300 PAM */
  
  cd = -log(1 - adjusted);
  
  return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.
 
 Based on Christian Storm's 5th order polynomial curve fitting of ROSE simulated data
 
 od = observed distance in percent 
 
 */
static double treeSTORMSONN(double od)
{
  double cd;			/* Corrected distance */
  
  if (od > 91.6) return 1000.0;	/* Otherwise log of negative value below */
  
  od /= 100;
  
  cd= -log(1 
           -0.95844*od 
           -0.69957*od*od
           +2.4955*od*od*od
           -4.6353*od*od*od*od
           +2.8076*od*od*od*od*od);
  
  /* printf(" od=%f  cd=%f\n", od, cd); */
  
  if (cd > 3.0) cd=3.0; /* Limit to 300 PAM */
  
  return cd*100;
}


static double treeSCOREDIST(char *seq1, char *seq2, BelvuContext *bc)
{
  char *s1 = seq1;
  char *s2 = seq2;
  
  
  /* Calc scores */
  int len = 0;
  int sc = 0;
  int s1sc = 0;
  int s2sc = 0;
  
  for (len=0; *s1; s1++, s2++) 
    {
      if (bc->penalize_gaps) 
        {
          if (isGap(*s1) || isGap(*s2)) sc -= 0.6;	
        }
      
      if (!isGap(*s1) && !isGap(*s2)) 
        {
          int val1 = a2b[(unsigned char)(*s1)];
          int val2 = a2b[(unsigned char)(*s2)];
          
          if (val1 > 0 && val2 > 0)
            {
              sc += (double) BLOSUM62[val1 - 1][val2 - 1];
              s1sc += (double) BLOSUM62[val1 - 1][val1 - 1];
              s2sc += (double) BLOSUM62[val2 - 1][val2 - 1];
            }

          len++;
        }
    }
  
  double maxsc = (s1sc + s2sc) / 2.0;
  
  /* Calc expected score */
  double expect =  -0.5209 * len;
  
  double od = ((double)sc - expect) / (maxsc - expect);
  
  if (!len || od < 0.05) od = 0.05;  /* Limit to 300 PAM;  len==0 if no overlap */
  if (od > 1.0) od = 1.0; 
  
  double cd = -log(od);
  
  DEBUG_OUT("SCOREDIST: len=%d  sc=%d  maxsc=%.2f  expect=%.2f  maxsc-expect=%.2f  od=%.3f\n", len, sc, maxsc, expect, maxsc-expect, od);
  
  cd = cd * 100.0 * 1.337 /* Magic scaling factor optimized for Dayhoff data */ ;
  
  if (cd > 300) cd = 300;  /* Limit to 300 PAM */
  
  /*    if (!bc->penalize_gaps) {
   g_free(s1);
   g_free(s2);
   } */
  
  return cd;
}



/* Sum branchlengths, allow going up parents
 
 If going up parents, arg1 = parent  ;  arg2 = child.
 
 If you're not interested in going up parents, simply call me with
 the same node in both arguments.
 
 */
static double treeSize3way(TreeNode *node, TreeNode *fromnode) 
{
  int root = 0;
  TreeNode *left = NULL;
  TreeNode *right = NULL;
  double len = 0.0;
  
  if (!node) 
    return 0.0;
  
  /* Get the correct branch length */
  if (node->left == fromnode || node->right == fromnode) /* Going up the tree */
    len = fromnode->branchlen;
  else
    len = node->branchlen;
  
  if (node->left == fromnode) 
    {
      left = node->parent;
      if (!left) 
        root = 1;
    }
  else
    {
      left = node->left;
    }
  
  if (node->right == fromnode) 
    {
      right = node->parent;
      if (!right) 
        root = 1;
    }
  else 
    {
      right = node->right;
    }
  
  if (root)
    {
      double retval;
      
      /* go across tree root */
      if (left) /* Coming from right */
        retval = treeSize3way(left, left);
      else  /* Coming from left */
        retval = treeSize3way(right, right);
      
      DEBUG_OUT("Returning (%.1f + %.1f) = %.1f\n", fromnode->branchlen, retval, fromnode->branchlen + retval);
      
      return fromnode->branchlen + retval;
    }
  else 
    {
      double 
      l = treeSize3way(left, node),
      r = treeSize3way(right, node),
      retval = (l + r)/2.0 + len;
      
      DEBUG_OUT("Returning (%.1f + %.1f)/2 + %.1f = %.1f\n", l, r, len, retval);
      
      return retval;
    }
}


/* Calculate the difference between left and right trees if the tree
 were to be rerooted at this node.
 
 What is calculated (bal) is the difference between 'left' and
 'right' subtrees minus the length of the branch itself. If this
 difference is negative, a perfectly balanced tree can be made.  For
 imperfectly balanced trees we want to root at the branch that gives
 the best balance.  However, perfectly balanced trees are all
 'perfect' so here we chose the branch with most equal subtrees.
 
 Actually it is not "left" and "right" but "down" and "up" subtrees.  */
static void treeCalcBalance(BelvuContext *bc, TreeNode *node) 
{
  double bal, lweight, rweight;
  
  if (node == bc->treeBestBalancedNode) 
    return;
  
  DEBUG_OUT("Left/Downstream weight\n");
  
  lweight = treeSize3way(node, node);
  
  DEBUG_OUT("Right/Upstream weight\n");
  
  rweight = treeSize3way(node->parent, node);
  
  bal = fabsf(lweight - rweight) - node->branchlen;
  
  DEBUG_OUT("Node=%s (branchlen = %.1f).  Weights = %.1f  %.1f. Bal = %.1f\n", 
            node->name, node->branchlen, lweight, rweight, bal);
  
  if (bal < bc->treeBestBalance) 
    { /* better balance */
      if (bc->treeBestBalance > 0.0 ||
          /* If previous tree was not perfectly balanced, or
       If previous tree was perfectly balanced - choose root with best subtree balance */
          fabsf(lweight - rweight) < bc->treeBestBalance_subtrees)
        {
          DEBUG_OUT("            %s has better balance %.1f < %.1f\n", node->name, bal, bc->treeBestBalance);
          
          bc->treeBestBalancedNode = node;
          bc->treeBestBalance = bal;
          bc->treeBestBalance_subtrees = fabsf(lweight - rweight);
        }
    }
}


static TreeNode *treeParent2leaf(TreeNode *newparent, TreeNode *curr)
{
  if (!curr->parent) 
    { /* i.e. the old root */
      if (curr->left == newparent)
        {
          newparent->branchlen += curr->right->branchlen; /* Add second part of this vertex */
          return curr->right;
        }
      else
        {
          newparent->branchlen += curr->left->branchlen;
          return curr->left;
        }
    } 
  else
    {
      if (curr->left == newparent) 
        {
          /* Change the link to the new parent to the old parent */
          curr->left = treeParent2leaf(curr, curr->parent);
          curr->left->branchlen = curr->branchlen;
        }
      else
        {
          curr->right = treeParent2leaf(curr, curr->parent);
          curr->right->branchlen = curr->branchlen;
        }
    }
  
  return curr;
}


/* Balance the two sides of a tree branch by the weights of the two subtrees.
 The branch has two sides: llen and rlen.
 
 Rationale: The correct center point balances the treesizes on both sides.
 
 Method: Calculate half the unbalance, add it to the appropriate side.
 
 No real theory for this, but it seems to work in easy cases 
 */
static void treeBalanceByWeight(TreeNode *lnode, TreeNode *rnode, double *llen, double *rlen)
{
  double adhocRatio = 0.95;
  
  double 
  halfbal, 
  branchlen = *rlen+*llen;
  double lweight = treeSize3way(lnode, lnode) /*- lnode->branchlen*/;
  double rweight = treeSize3way(rnode, rnode) /*- rnode->branchlen*/;
  
  halfbal = fabsf(lweight-rweight) / 2.0;
  
  if (halfbal < *llen && halfbal < *rlen) {
    
    if (lweight < rweight) {
      *llen += halfbal;
      *rlen -= halfbal;
    }
    else {
      *llen -= halfbal;
      *rlen += halfbal;
    }
  }
  else {
    /* The difference is larger than the final branch -
     give nearly all weight to the shorter one (cosmetic hack) */
    if (lweight < rweight) {
      *llen = branchlen*adhocRatio;
      *rlen = branchlen*(1.0 - adhocRatio);
    }
    else {
      *rlen = branchlen*adhocRatio;
      *llen = branchlen*(1.0 - adhocRatio);
    }
  }
}


static void fillParents(TreeNode *parent, TreeNode *node)
{
  if (!node) 
    return;
  
  node->parent = parent;
  fillParents(node, node->left);
  fillParents(node, node->right);
}


static char *fillOrganism(TreeNode *node)
{
  char 
  *leftorganism,
  *rightorganism;
  
  if (node->name) 
    return node->organism;
  
  leftorganism = fillOrganism(node->left);
  rightorganism = fillOrganism(node->right);
  
  /* printf("\nFill: (left=%s, right=%s):  ", leftorganism, rightorganism); */
  
  node->organism = (leftorganism == rightorganism ? leftorganism : NULL);
  
  return node->organism;
}


/* Allocate memory for a new TreeNode and initialise its contents to empty values */
static TreeNode* createEmptyTreeNode()
{
  TreeNode *result = (TreeNode*)g_malloc(sizeof *result);
  
  result->dist = 0.0;
  result->branchlen = 0.0;
  result->boot = 0.0;
  result->left = NULL;
  result->right = NULL;
  result->parent = NULL;
  result->name = NULL;
  result->organism =NULL;
  result->aln = NULL;
  result->box = 0;
  result->color = 0;
  
  return result;
}


/* Recursively destroy the tree node and all its child nodes */
static void destroyTreeNode(TreeNode **node)
{
  if (node == NULL || *node == NULL)
    return;

  /* Destroy all child nodes first */
  destroyTreeNode(&((*node)->left));
  destroyTreeNode(&((*node)->right));

  /* Free memory owned by the node */
  if ((*node)->name)
    {
      g_free((*node)->name);
      (*node)->name = NULL;
    }

  /* Free the node struct itself, and set the pointer to it to null */
  g_free(*node);
  *node = NULL;
}


/* Swap the left and right branches of the given node */
static void treeSwapNode(TreeNode *node)
{
  void *tmp;
  
  tmp = node->left;
  node->left = node->right;
  node->right = (TreeNode*)tmp;
}


/* Rerooting works roughly this way:
 
 - A new node is created with one child being the node chosen as new root 
 and the other child the previous parent of it.  The sum of left and right
 branch lengths of the new root should equal the branch length of the chosen node.
 - All previous parent nodes are visited and converted so that:
 1. The previous parent node becomes a child.
 2. The new parent is the node that calls.
 3. The branchlength of the previous parent node is assigned to its parent.
 - When the previous root is reached, it is deleted and the other child
 becomes the child of the calling node.
 
 Note that treeReroot destroys the old tree, since it reuses the nodes.  
 Use treecpy if you still need it for something later.
 */
static TreeNode *treeReroot(TreeNode *node)
{
  TreeNode *newroot = createEmptyTreeNode();
  
  newroot->left = node;
  newroot->right = treeParent2leaf(node, node->parent);
  
  fillParents(newroot, newroot->left);
  fillParents(newroot, newroot->right);
  
  newroot->left->branchlen = newroot->right->branchlen = node->branchlen / 2.0;
  treeBalanceByWeight(newroot->left, newroot->right, 
                      &newroot->left->branchlen, &newroot->right->branchlen);
  
  return newroot;
}


static gboolean doublesEqual(const double a, const double b)
{
  return (b < a + MACHINE_RES && b > a - MACHINE_RES);
}


/* Find the node which has most equal balance, return new tree with this as root.
 */
static TreeNode *treeFindBalance(BelvuContext *bc, TreeNode *node) 
{
  double lweight = treeSize3way(node->left, node->left);
  double rweight = treeSize3way(node->right, node->right);
  
  bc->treeBestBalancedNode = node;
  bc->treeBestBalance = fabsf(lweight - rweight);
  
  bc->treeBestBalance_subtrees = 
  fabsf((lweight - node->left->branchlen) - (rweight - node->right->branchlen));
  
  DEBUG_OUT("Initial weights = %.1f  %.1f. Bal = %.1f\n", lweight, rweight, bc->treeBestBalance);
  
  treeTraverse(bc, node, treeCalcBalance);
  
  if (bc->treeBestBalancedNode == node)
    return node;
  else
    return treeReroot(bc->treeBestBalancedNode);
}


/* Calculate the pairwise tree distances */
static void calcPairwiseDistMatrix(BelvuContext *bc, double **pairmtx)
{
  /* Calculate pairwise distance matrix. Note that this only calculates
   * the portion of the array above the diagonal (where j > i); the other
   * values are left uninitialised and should not be used. */
  arrayOrder(bc->alignArr);
  
  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len - 1; ++i)
    {
      ALN *aln_i = g_array_index(bc->alignArr, ALN*, i);
      char *alniSeq = alnGetSeq(aln_i);
      
      int j = i+1;
      for (j = i+1; j < (int)bc->alignArr->len; ++j)
        {
          ALN *aln_j = g_array_index(bc->alignArr, ALN*, j);
          char *alnjSeq = alnGetSeq(aln_j);
          
          pairmtx[i][j] = 100.0 - identity(alniSeq, alnjSeq, bc->penalize_gaps);
          
          if (bc->treeDistCorr == KIMURA) 
            pairmtx[i][j] = treeKimura(pairmtx[i][j]);
          else if (bc->treeDistCorr == JUKESCANTOR) 
            pairmtx[i][j] = treeJUKESCANTOR(pairmtx[i][j]);
          else if (bc->treeDistCorr == STORMSONN) 
            pairmtx[i][j] = treeSTORMSONN(pairmtx[i][j]);
          else if (bc->treeDistCorr == SCOREDIST) 
            pairmtx[i][j] = treeSCOREDIST(alniSeq, alnjSeq, bc);
        }
    }
}


/* Output the tree distances */
static void printTreeDistances(BelvuContext *bc, double **pairmtx)
{
  double dist;

  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; i++) 
    {
      if (!bc->treeCoordsOn) 
        {
          g_message("%s\t", g_array_index(bc->alignArr, ALN*, i)->name);
        }
      else
        {
          g_message("%s/%d-%d\t",
                    g_array_index(bc->alignArr, ALN*, i)->name,
                    g_array_index(bc->alignArr, ALN*, i)->start,
                    g_array_index(bc->alignArr, ALN*, i)->end);
        }
    }
  g_message ("\n");
  
  for (i = 0; i < (int)bc->alignArr->len; ++i) 
    {
      int j = 0;
      for (j = 0; j < (int)bc->alignArr->len; ++j) 
        {
          if (i == j) 
            dist = 0.0;
          else if (i < j)
            dist = pairmtx[i][j];
          else
            dist = pairmtx[j][i];
          
          g_message("%7.3f\t", dist);
        }
      g_message ("\n");
    }
}


#ifdef DEBUG
/* Print the given matrix */
static void printMtx(BelvuContext *bc, double **mtx) 
{
  int i, j;
  
  g_message ("\n");
  for (i = 0; i < bc->alignArr->len; i++) 
    {
      for (j = 0; j < bc->alignArr->len; j++)
        g_message("%6.2f ", mtx[i][j]);
      
      g_message ("\n");
  }
}


/* Print staisitics about the tree construction */
static void printTreeStats(BelvuContext *bc, double **pairmtx, double *avgdist, double **Dmtx, TreeNode **node)
{
  g_message("Node status, Avg distance:\n");
  
  int i = 0;
  for (i = 0; i < bc->alignArr->len; i++)
    g_message("%6d ", (node[i] ? 1 : 0)); 
  
  g_message("\n");
  
  for (i = 0; i < bc->alignArr->len; i++)
    g_message("%6.2f ", avgdist[i]); 
  
  g_message("\n\nPairdistances, corrected pairdistances:");
  printMtx(bc, pairmtx);
  printMtx(bc, Dmtx);
  g_message("\n");
}
#endif


/* Construct the tree using the NJ method. Populates Dmtx */
static void constructTreeNJ(BelvuContext *bc, double **pairmtx, double *avgdist, TreeNode **node, double **Dmtx)
{
  int count = 0;
  
  /* Calculate vector r (avgdist) */
  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      if (!node[i]) 
        continue;
      
      avgdist[i] = 0.0;
      
      int j = 0;
      for (count=0, j = 0; j < (int)bc->alignArr->len; j++)
        {
          if (!node[j]) 
            continue;
          
          if (i < j)
            avgdist[i] += pairmtx[i][j];
          else if (i > j)
            avgdist[i] += pairmtx[j][i];
          
          count++;
        }
      
      if (count == 2)	/* Hack, to accommodate last node */
        avgdist[i] = 1;
      else
        avgdist[i] /= 1.0*(count - 2);
    }
  
  /* Calculate corrected matrix Dij. Note again that only values 
   * above the diagonal (where j > i) are populated; other values
   * are uninitialised and should not be used. */
  if (1 /* gjm */)
    {
      for (i = 0; i < (int)bc->alignArr->len - 1; ++i)
        {
          if (!node[i])
            continue;
          
          int j = i+1;
          for (j = i+1; j < (int)bc->alignArr->len; ++j)
            {
              if (!node[j]) 
                continue;
              
              Dmtx[i][j] = pairmtx[i][j] - (avgdist[i] + avgdist[j]);
            }
        }
    }
  else
    {		/* Felsenstein */
      double Q = 0.0;
      
      for (i = 0; i < (int)bc->alignArr->len - 1; ++i)
        {
          if (!node[i])
            continue;
          
          int j = i+1;
          for (j = i+1; j < (int)bc->alignArr->len; j++)
            {
              if (!node[j]) 
                continue;
              
              Q += pairmtx[i][j];
            }
        }
      
      for (i = 0; i < (int)bc->alignArr->len - 1; ++i)
        {
          if (!node[i])
            continue;
          
          int j = i+1;
          for (j = i+1; j < (int)bc->alignArr->len; j++)
            {
              if (!node[j])
                continue;
              
              Dmtx[i][j] = (pairmtx[i][j] +
                            (2.0*Q)/(count-(count == 2 ? 1 : 2)) - 
                            avgdist[i] - avgdist[j]) / 2.0;
            }
        }
    }
}


/* Find smallest distance pair in pairmtx */
static double treeFindSmallestDist(BelvuContext *bc, double **pairmtx, double **curMtx, double **Dmtx, TreeNode **node, int *maxiOut, int *maxjOut)
{
  int maxi = -1;
  int maxj = -1;
  double maxid = 1000000;
  double pmaxid = 1000000;

  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len - 1; ++i) 
    {
      if (!node[i]) 
        continue;
      
      int j = i+1;
      for (j = i+1; j < (int)bc->alignArr->len; ++j) 
        {
          if (!node[j]) 
            continue;
          
          /* printf("iter %d, i=%d. j=%d, dist= %f\n", iter, i, j, curMtx[i][j]);*/
          if (curMtx[i][j] < maxid) 
            {
              maxid = curMtx[i][j];
              pmaxid = pairmtx[i][j];
              maxi = i;
              maxj = j;
            }
          else if (bc->treeMethod == NJ && doublesEqual(Dmtx[i][j], maxid) && pairmtx[i][j] < pmaxid) 
            {
              /* To resolve ties - important for tree look! */
              maxi = i;
              maxj = j;
              pmaxid = pairmtx[i][j];
            }
        }
    }
  
  maxid = pairmtx[maxi][maxj]; /* Don't want to point to Dmtx in NJ */
  
  if (maxiOut)
    *maxiOut = maxi;
  
  if (maxjOut)
    *maxjOut = maxj;
  
  return maxid;
}


/* This does the work to create all the nodes in a tree. All the memory for
 * the nodes etc. is allocated using a BlxHandle which is stored in the tree.
 * To free the memory used by the tree, the handle should be destroyed. */
Tree* treeMake(BelvuContext *bc, const gboolean doBootstrap, const gboolean displayFeedback)
{
  /* This can take a long time, so let the user know we're doing something.
   * Only display feedback text if asked, though (e.g. we don't want this each
   * time if calculating a lot of bootstrap trees) */
  if (displayFeedback)
    {
      g_message_info("Calculating tree...\n");
    }

  setBusyCursor(bc, TRUE);

  TreeNode *newnode = NULL ;
  int maxi = -1, maxj = -1;
  double maxid = 0.0, **pairmtx, **Dmtx, **curMtx, *src, *trg, 
  *avgdist,		/* vector r in Durbin et al */
  llen = 0, rlen = 0;
  TreeNode **nodes;   /* Array of (primary) nodes.  Value=0 => stale column */

  /* Create the tree struct */
  Tree *tree = createEmptyTree();

  /* Allocate memory */
  BlxHandle localHandle = handleCreate(); /* handles local memory that will be free'd before we return */
  
  nodes = (TreeNode**)handleAlloc(&localHandle, bc->alignArr->len * sizeof(TreeNode *));
  pairmtx = (double**)handleAlloc(&localHandle, bc->alignArr->len * sizeof(double *));
  Dmtx = (double**)handleAlloc(&localHandle, bc->alignArr->len * sizeof(double *));
  avgdist = (double*)handleAlloc(&localHandle, bc->alignArr->len * sizeof(double));

  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *aln_i = g_array_index(bc->alignArr, ALN*, i);

      pairmtx[i] = (double*)handleAlloc(&localHandle, bc->alignArr->len*sizeof(double));
      Dmtx[i] = (double*)handleAlloc(&localHandle, bc->alignArr->len*sizeof(double));

      nodes[i] = createEmptyTreeNode();
      nodes[i]->name =  (char*)g_malloc(strlen(aln_i->name) + 50);
      
      if (!bc->treeCoordsOn) 
        sprintf(nodes[i]->name, "%s", aln_i->name);
      else
        sprintf(nodes[i]->name, "%s/%d-%d", aln_i->name, aln_i->start, aln_i->end);
      
      nodes[i]->aln = aln_i;
      nodes[i]->organism = aln_i->organism;
    }

  /* Get the pairwise tree distances (from file if given, or calculate them) */
  if (bc->treeReadDistancesOn)
    treeReadDistances(bc, pairmtx);
  else
    calcPairwiseDistMatrix(bc, pairmtx);

  /* If requested (or if debug is on), print the distance matrix */
  if (bc->treePrintDistances) 
    {
      printTreeDistances(bc, pairmtx);
      exit(0);
    }  

#ifdef DEBUG
  printTreeDistances(bc, pairmtx);
#endif
  
  /* Construct the tree */
  int iter = 0;
  for (iter = 0; iter < (int)bc->alignArr->len - 1; ++iter)
    {
      if (bc->treeMethod == NJ)
        {
          constructTreeNJ(bc, pairmtx, avgdist, nodes, Dmtx);
          curMtx = Dmtx;
        }
      else 
        {
          curMtx = pairmtx;
        }
      
#ifdef DEBUG
      printTreeStats(bc, pairmtx, avgdist, Dmtx, nodes);
#endif
      
      /* Find smallest distance pair in pairmtx */
      maxid = treeFindSmallestDist(bc, pairmtx, curMtx, Dmtx, nodes, &maxi, &maxj);

      /* Merge rows & columns of maxi and maxj into maxi
       Recalculate distances to other nodes */
      for (i = 0; i < (int)bc->alignArr->len; ++i)
        {
          if (!nodes[i])
            continue;
          
          if (i < maxi) 
            trg = &pairmtx[i][maxi];
          else if (i > maxi) 
            trg = &pairmtx[maxi][i];
          else continue;
          
          if (i < maxj) 
            src = &pairmtx[i][maxj];
          else if (i > maxj) 
            src = &pairmtx[maxj][i];
          else continue;
          
          if (bc->treeMethod == UPGMA) 
            *trg = (*trg + *src) / 2.0;
          else
            *trg = (*trg + *src - maxid) / 2.0;
        }
      
      /* Create node for maxi and maxj */
      newnode = createEmptyTreeNode();
      
      if (bc->treeMethod == UPGMA)
        {
          /* subtract lower branch lengths from absolute distance
           Horribly ugly, only to be able to share code UPGMA and NJ */
          TreeNode *tmpnode = nodes[maxi]->left;
          
          for (llen = maxid; tmpnode; tmpnode = tmpnode->left)
            llen -= tmpnode->branchlen;
          
          tmpnode = nodes[maxj]->right;
          
          for (rlen = maxid; tmpnode; tmpnode = tmpnode->right)
            rlen -= tmpnode->branchlen;
        }
      else
        {
          llen = (maxid + avgdist[maxi] - avgdist[maxj]) / 2.0;
          rlen = maxid - llen;
          
          if (iter == (int)bc->alignArr->len - 2)
            { /* Root node */
              
              /* Not necessary anymore, the tree is re-balanced at the end which calls this too
               treeBalanceByWeight(node[maxi], node[maxj], &llen, &rlen);*/
              
              /* Put entire length of root branch in one leg so the rebalancing 
               will work properly (otherwise it is hard to take this branch into account */
              rlen += llen;
              llen = 0;
            }
        }
      
      DEBUG_OUT("Iter %d: Merging %d and %d, dist= %f\n", iter, maxi+1, maxj+1, curMtx[maxi][maxj]);
      DEBUG_OUT("maxid= %f  llen= %f  rlen= %f\n", maxid, llen, rlen);
      DEBUG_OUT("avgdist[left]= %f  avgdist[right]= %f\n\n", avgdist[maxi], avgdist[maxj]);
      
      newnode->left = nodes[maxi];
      newnode->left->branchlen = llen;
      
      newnode->right = nodes[maxj];
      newnode->right->branchlen = rlen;
      
      newnode->organism = (nodes[maxi]->organism == nodes[maxj]->organism ?
                           nodes[maxi]->organism : NULL);
      
      nodes[maxi] = newnode;
      nodes[maxj] = NULL;
    }

  fillParents(newnode, newnode->left);
  fillParents(newnode, newnode->right);
  
  if (bc->treeMethod == UPGMA)
    newnode->branchlen = 100 - maxid ;
  
  if (bc->treeMethod == NJ)
    newnode = treeFindBalance(bc, newnode) ;
  
  fillOrganism(newnode);

  /* Set the newnode to be the root node of the tree */
  tree->head = newnode;

  /* Clean up locally allocated memory */
  handleDestroy(&localHandle);

  if (doBootstrap && bc->treebootstraps) 
    treeBootstrapStats(bc, tree);
  
  setBusyCursor(bc, FALSE);
  
  if (displayFeedback)
    {
      g_message_info("Finished calculating tree.\n");
    }
  
  return tree;
}


/***********************************************************
 *                   Find Orthologs                        *
 ***********************************************************/

static void treePrintNode(BelvuContext *bc, TreeNode *node) 
{
  if (node->name) 
    g_message("%s ", node->name);
}


static gboolean treePrintOrthologsRecur(BelvuContext *bc, TreeNode *node) 
{
  gboolean found = FALSE;
  
  if (!node || !node->left || !node->right) 
    return found;
  
  DEBUG_OUT("\n 1 (%s, seq=%s):  ", node->left->organism, node->left->name);
  DEBUG_OUT("\n 2 (%s, seq=%s)\n: ", node->right->organism, node->right->name);
  
  if (node->left->organism && node->right->organism &&
      node->left->organism != node->right->organism) 
    {
      found = TRUE;
      
      g_message("\nSpecies 1 (%s):  ", node->left->organism);
      treeTraverse(bc, node->left, treePrintNode);
      g_message("\nSpecies 2 (%s): ", node->right->organism);
      treeTraverse(bc, node->right, treePrintNode);
      g_message("\n");
    }
  else 
    {
      if (treePrintOrthologsRecur(bc, node->left))
        found = TRUE;
      
      if (treePrintOrthologsRecur(bc, node->right))
        found = TRUE;
    }
  
  return found;
}


void treePrintOrthologs(BelvuContext *bc, GtkWidget *treeWindow) 
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(treeWindow);
  
  if (properties->tree == NULL)
    {
      g_critical("Tree has not been calculated.\n");
      return;
    }
  else if (treePrintOrthologsRecur(bc, properties->tree->head))
    {
      g_message_info("Found orthologs\n");
    }
  else
    {
      g_message_info("No orthologs\n");
    }
}


/***********************************************************
 *                        Drawing                          *
 ***********************************************************/

/* This is called to re-make the underlying tree for the given tree
 * widget. It must be called after a change that renders the existing
 * tree invalid, e.g. deleting an alignment or changing the build method */
void belvuTreeRemakeTree(GtkWidget *belvuTree)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  BelvuContext *bc = properties->bc;

  /* Re-make the tree */
  separateMarkupLines(bc);
  
  Tree *newTree = treeMake(bc, TRUE, TRUE);
  
  /* Update the contents of our tree with the new member variables (we don't 
   * change the actual tree pointer if it already exists because the context 
   * may point to the same tree struct). */
  if (properties->tree)
    {
      destroyTreeContents(properties->tree);
      
      if (newTree)
        {
          properties->tree->head = newTree->head;
          g_free(newTree);
          newTree = properties->tree;
        }
    }
  else
    {
      properties->tree = newTree;
    }
  
  reInsertMarkupLines(bc);
  
  /* Make sure our properties are up to date with the data used to create
   * the new tree. */
  properties->buildMethod = bc->treeMethod;
  properties->distCorr = bc->treeDistCorr;
  
  calculateBelvuTreeBorders(belvuTree);
  
  /* If the alignment is sorted by tree and we've changed the main belvu tree's
   * order, we need to refresh the sort order */
  if (properties->isMainTree)
    onTreeOrderChanged(bc);
}


/* This is called when the tree settings have been changed. */
static void belvuTreeUpdateSettings(BelvuContext *bc, GtkWidget *treeWindow)
{
  if (!treeWindow)
    return;
  
  /* The tree window exists, so must be updated */
  BelvuTreeProperties *properties = belvuTreeGetProperties(treeWindow);
  g_assert(properties);
  
  if (properties->tree == NULL || properties->tree->head == NULL)
    {
      /* The underlying tree has been invalidated, so we'll need to re-make
       * the whole tree. */
      belvuTreeRemakeTree(treeWindow);
    }
  else
    {
      /* Just redraw the existing tree */
      calculateBelvuTreeBorders(treeWindow);
      belvuTreeRedrawAll(treeWindow, NULL);
    }
}


/* Clear any cached drawables and redraw everything. It also recalculates
 * the tree if the build method has changed. */
void belvuTreeRedrawAll(gpointer widget, gpointer data)
{
  if (!widget || !stringsEqual(gtk_widget_get_name(GTK_WIDGET(widget)), BELVU_TREE_WINDOW_NAME, TRUE))
    return;
  
  GtkWidget *belvuTree = GTK_WIDGET(widget);
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);

  widgetClearCachedDrawable(properties->treeArea, NULL);
  gtk_widget_queue_draw(belvuTree);
}


static void createClickableRect(BelvuTreeProperties *properties,
                                TreeNode *node,
                                const int x,
                                const int y,
                                const int width,
                                const int height,
                                const gboolean isBranch)
{
  ClickableRect clickRect;

  clickRect.node = node;
  clickRect.isBranch = isBranch;
  
  clickRect.rect.x = x;
  clickRect.rect.y = y;
  clickRect.rect.width = width;
  clickRect.rect.height = height;
  
  g_array_append_val(properties->clickableRects, clickRect);
}


/* The actual tree drawing routine.
 Note: must be in sync with treeDrawNodeBox, which draws clickable
 boxes first.
 */
static double treeDrawNode(BelvuContext *bc, 
                           GtkWidget *widget,
                           GdkDrawable *drawable, 
                           GdkGC *gc, 
                           BelvuTreeProperties *properties,
                           GdkColor *defaultColor, 
                           TreeNode *node, 
                           double x) 
{
  double y, yl, yr;
  
  if (!node) 
    return 0.0;
  
  const int curX = x + roundNearest(node->branchlen * (double)(treeGetScale(properties) * properties->charWidth));

  GdkGC *leftGc = gdk_gc_new(drawable);
  gdk_gc_copy(leftGc, gc);
  yl = treeDrawNode(bc, widget, drawable, leftGc, properties, defaultColor, node->left, curX);
  
  gdk_gc_set_foreground(gc, defaultColor);
  yr = treeDrawNode(bc, widget, drawable, gc, properties, defaultColor, node->right, curX);
  
  GdkGC *gcTmp = gdk_gc_new(drawable);
  gdk_gc_set_foreground(gcTmp, defaultColor);

  if (yl) 
    {
      /* internal node */
      y = (yl + yr) / 2.0;
      
      /* connect children */
      gdk_draw_line(drawable, gc, curX, yr, curX, y);
      gdk_gc_copy(gc, leftGc);
      gdk_draw_line(drawable, gc, curX, yl, curX, y);
      
      if (node->left->organism != node->right->organism)
        {
          /* Reset color */
          gdk_gc_set_foreground(gc, defaultColor);
        }
    }
  else 
    {
      /* Sequence name */
      const gboolean isSelected = (bc->selectedAln && bc->selectedAln == node->aln);
      
      y = bc->tree_y * properties->charHeight;
      bc->tree_y++;
      
      if (bc->treeColorsOn && node->organism) 
        {
          /* Get the color for this organism */
          ALN aln;
          initAln(&aln);
          
          aln.organism = node->organism;
          
          int ip = 0;
          if (alnArrayFind(bc->organismArr, &aln, &ip, organism_order)) 
            {
              GdkColor color;
              int colorNum = g_array_index(bc->organismArr, ALN*, ip)->color;
              convertColorNumToGdkColor(colorNum, FALSE, &color); /* we currently don't change the text color when the node is selected */
              
              gdk_gc_set_foreground(gc, &color);
	    }
	}
      else
        {
          /* Just use the default color */
          gdk_gc_set_foreground(gc, defaultColor);
        }
      
      /* Draw the sequence name */
      int nameWidth = 0, nameHeight = 0;
      const int textX = curX + DEFAULT_XPAD;
      const int textY = y - properties->charHeight / 2;

      drawText(widget, drawable, gcTmp, textX, textY, node->name, &nameWidth, &nameHeight);

      if (isSelected || (node->aln && node->aln->color != WHITE))
        {
          /* Highlight this node */
          GdkColor color;
          convertColorNumToGdkColor(node->aln->color, isSelected, &color);
          
          gdk_gc_set_foreground(gcTmp, &color);
          gdk_draw_rectangle(drawable, gcTmp, TRUE, textX - DEFAULT_XPAD/2, textY - 1, nameWidth + DEFAULT_XPAD, nameHeight + 2);
          
          /* This is a bit hacky because it re-draws text we've already drawn because it was covered by the
           * background. We should rearrange things slightly so that we can get nameWidth and nameHeight another
           * way so we can avoidthis, but it's a very small performance hit so not worth worrying about. */
          GdkColor *fgColor = getGdkColor(BELCOLOR_TREE_TEXT, bc->defaultColors, isSelected, FALSE);
          gdk_gc_set_foreground(gcTmp, fgColor);
          drawText(widget, drawable, gcTmp, curX + DEFAULT_XPAD, y - properties->charHeight / 2, node->name, NULL, NULL);
        }
      
      /* Make a clickable box for the sequence name */
      createClickableRect(properties, node, textX, textY, nameWidth, properties->charHeight, FALSE);

      if (treeGetShowOrganism(properties) && node->organism) 
        {
          drawText(widget, drawable, gc, curX + nameWidth + DEFAULT_XPAD * 2, y - properties->charHeight / 2, node->organism, NULL, NULL);
        }
    }
  
  /* Horizontal branches. If highlighting orthologs, draw a shaded area underneath
   * the branch. */
  if (bc->highlightOrthologs && 
      node->left && node->right &&
      node->left->organism && node->right->organism &&
      node->left->organism != node->right->organism) 
    {
      gdk_gc_set_foreground(gcTmp, defaultColor);
      gdk_draw_rectangle(drawable, gcTmp, TRUE, x, y - properties->charHeight/2, curX - x, properties->charHeight);
    }
  
  gdk_draw_line(drawable, gc, curX, y, x, y);
  createClickableRect(properties, node, x, y - properties->charHeight/2, curX - x, properties->charHeight, TRUE);
  
  if (treeGetShowBranchLen(properties) && node->branchlen) 
    {
      /* Draw the branch label, which shows the branch length */
      char *tmpStr = g_strdup_printf("%.1f", node->branchlen);

      const int textWidth = getTextWidth(widget, tmpStr);
      double pos = min(curX, x) + (abs(curX - x) / 2) - (textWidth * 0.5); /* centre text at middle of branch */

      drawText(widget, drawable, gc, pos, y, tmpStr, NULL, NULL);
      
      g_free(tmpStr);
    }
  
  /* This prints the tree bootstrap statistics */
  if (bc->treebootstraps && !node->name && node != properties->tree->head && !bc->treebootstrapsDisplay) 
    {
      GdkColor *color = getGdkColor(BELCOLOR_TREE_BOOTSTRAP, bc->defaultColors, FALSE, FALSE);
      gdk_gc_set_foreground(gc, color);

      char *tmpStr = g_strdup_printf("%.0f", node->boot);
      const int textWidth = getTextWidth(widget, tmpStr);
      double pos = curX - textWidth - 0.5;
      
      if (pos < 0.0) 
        pos = 0;
      
      g_message("%f  %f   \n", node->boot, pos / properties->charWidth);
      drawText(widget, drawable,  gc, pos, y, tmpStr, NULL, NULL);

      g_free(tmpStr);
      gdk_gc_set_foreground(gc, defaultColor);
    }
  
  g_object_unref(leftGc);

  return y;
}


/* Destroy the contents of the given tree and reset all member variables to
 * empty values. Does not destroy the tree struct itself */
void destroyTreeContents(Tree *tree)
{
  if (tree == NULL)
    return;

  /* Free all memory used by the tree nodes */
  if (tree->head)
    destroyTreeNode(&tree->head);

  tree->head = NULL;
}


/* Destroy the given tree's contents and the tree itself. */
void destroyTree(Tree **tree)
{
  if (tree == NULL)
    return;
  
  destroyTreeContents(*tree);
  
  /* Free the tree struct */
  g_free(*tree);
  *tree = NULL;
}


/* Create a tree struct, initialised with sensible empty values */
static Tree* createEmptyTree()
{
  Tree *result = (Tree*)g_malloc(sizeof *result);
  
  result->head = NULL;
  
  return result;
}


static void drawBelvuTree(GtkWidget *widget, GdkDrawable *drawable, BelvuTreeProperties *properties)
{
  BelvuContext *bc = properties->bc;

  if (!properties->tree || !properties->tree->head)
    return;
  
  /* Clear any previous clickable rects that were created */
  g_array_unref(properties->clickableRects);
  properties->clickableRects = g_array_new(FALSE, FALSE, sizeof(ClickableRect));
  
  bc->tree_y = 1;

  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *defaultColor = getGdkColor(BELCOLOR_TREE_LINE, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_line_attributes(gc, treeGetLineWidth(properties) * properties->charWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);
  
  treeDrawNode(bc, widget, drawable, gc, properties, defaultColor, properties->tree->head, properties->treeRect.x);

  int xscale = treeGetScale(properties) * properties->charWidth;
  int yscale = properties->charHeight;
  
  const int markerHt = 0.5 * yscale;
  
  const int xStart = properties->treeRect.x;
  const int xWidth = 10 * xscale;
  const int xHalfWidth = xWidth / 2;
  int xMax = 0;
  
  /* Draw scale */
  if (bc->treeMethod == UPGMA) 
    {
      xMax = xStart + (100 * xscale);

      bc->tree_y += 2;
      int y = bc->tree_y * yscale;

      gdk_draw_line(drawable, gc, xStart, y, xMax, y);
      
      int x = xStart;
      int i = 1;
      for (; i <= 101; i += 10, x += xWidth) 
        {
          if (i < 101) 
            {
              /* Draw half-way marker lines */
              gdk_draw_line(drawable, gc, x + xHalfWidth, y, x + xHalfWidth, y + markerHt);
            }
          
          /* Draw the full marker line, with a label */
          gdk_draw_line(drawable, gc, x, y, x, y + markerHt);

          char *tmpStr = g_strdup_printf("%d", i - 1);
          drawText(widget, drawable, gc,  x - (0.5 * xscale), y + (2 * markerHt), tmpStr, NULL, NULL);
          g_free(tmpStr);
        }
    }
  else
    {
      xMax = xStart + (10 * xscale);

      bc->tree_y += 3;
      int y = bc->tree_y * yscale;
      int x = xStart;
      
      drawText(widget, drawable, gc, x + xHalfWidth, y - 3 * markerHt, "0.1", NULL, NULL);
      
      gdk_draw_line(drawable, gc, x, y, x + xWidth, y);
      gdk_draw_line(drawable, gc, x, y - markerHt, x, y + markerHt);
      gdk_draw_line(drawable, gc, x + xWidth, y - markerHt, x + xWidth, y + markerHt);
    }
  
  if (bc->treeMethod == NJ) 
    {	
      int y = bc->tree_y * yscale;

      double lweight, rweight;
      
      lweight = treeSize3way(properties->tree->head->left, properties->tree->head);
      rweight = treeSize3way(properties->tree->head->right, properties->tree->head);
      
      char *tmpStr = NULL;
#ifdef DEBUG
      tmpStr = g_strdup_printf("Tree balance = %.1f (%.1f-%.1f)", fabsf(lweight - rweight), lweight, rweight);
#else
      tmpStr = g_strdup_printf("Tree balance = %.1f", fabsf(lweight - rweight));
#endif
      
      drawText(widget, drawable, gc, xMax + 2 * DEFAULT_XPAD, y - markerHt, tmpStr, NULL, NULL);
      g_free(tmpStr);
    }
  
  g_object_unref(gc);
}


/* Expose function for the drawing area. The main belvuTree widget is passed as
 * the user data. */
static gboolean onExposeBelvuTree(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuTree = GTK_WIDGET(data);
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          if (properties->treeRect.width > 0 && properties->treeRect.height > 0)
            {
              bitmap = createBlankSizedPixmap(widget, window, properties->treeRect.x * 2 + properties->treeRect.width, 
                                              properties->treeRect.y * 2 + properties->treeRect.height);
              
              separateMarkupLines(properties->bc);
              drawBelvuTree(widget, bitmap, properties);
              reInsertMarkupLines(properties->bc);
            }
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
	  g_warning("Failed to draw Belvu tree [%p] - could not create bitmap.\n", widget);
	}
    }
  
  return TRUE;
}


/***********************************************************
 *                       Mouse events                      *
 ***********************************************************/

/* Called when the tree has been left-clicked */
static void onLeftClickTree(GtkWidget *belvuTree, const int x, const int y)
{
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  BelvuContext *bc = properties->bc;

  /* We store a list of rectangles in our properties that tell us where 
   * clickable items lie. Loop through them and see if the click point lies
   * inside any of them. */
  ClickableRect *foundRect = NULL;

  int i = 0;
  for ( ; i < (int)properties->clickableRects->len; ++i)
    {
      ClickableRect *clickRect = &g_array_index(properties->clickableRects, ClickableRect, i);
      
      if (pointInRect(x, y, &clickRect->rect))
        {
          foundRect = clickRect;
          break; /* we shouldn't have overlapping items, so exit once we have found one */
        }
    }

  if (foundRect)
    {
      if (foundRect->isBranch)
        {
          BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);

          /* We clicked on a tree branch - swap or re-root */
          if (bc->treePickMode == NODESWAP)
            {
              treeSwapNode(foundRect->node);
            }
          else if (bc->treePickMode == NODEROOT)
            {
              /* Re-routing changes the tree's head node. */
              properties->tree->head = treeReroot(foundRect->node);
              
              /* Re-routing can also affect the drawing area size, so recalculate borders */
              calculateBelvuTreeBorders(belvuTree);
            }
          else
            {
              g_warning("Program error: unrecognised tree selection mode '%d'.\n", bc->treePickMode);
            }
          
          belvuTreeRedrawAll(belvuTree, NULL);

          /* If this is the main tree, then we need to update the main alignment
           * to reflect the new tree order. */
          if (properties->isMainTree)
            onTreeOrderChanged(bc);
        }
      else if (foundRect->node)
        {
          bc->selectedAln = foundRect->node->aln;
          onRowSelectionChanged(bc);
        }
    }
}

static gboolean onButtonPressBelvuTree(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      GtkWidget *belvuTree = GTK_WIDGET(data);
      onLeftClickTree(belvuTree, event->x, event->y);
      handled = TRUE;
    }
  
  return handled;
}


/* This should be called if the font size of the tree window has changed. */
void onBelvuTreeFontSizeChanged(GtkWidget *belvuTree)
{
  if (!belvuTree)
    return;
  
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);
  
  getFontCharSize(properties->treeArea, properties->treeArea->style->font_desc, &properties->charWidth, &properties->charHeight);
  
  calculateBelvuTreeBorders(belvuTree);
  belvuTreeRedrawAll(belvuTree, NULL);
}


static gboolean onZoomBelvuTree(BelvuContext *bc, const gboolean zoomIn)
{
  if (bc->belvuTree)
    {
      int size = pango_font_description_get_size(bc->belvuTree->style->font_desc) / PANGO_SCALE;
      
      if (zoomIn)
        widgetSetFontSizeAndCheck(bc->belvuTree, size + 1);
      else
        widgetSetFontSizeAndCheck(bc->belvuTree, size - 1);
      
      onBelvuTreeFontSizeChanged(bc->belvuTree);
    }  
  
  return TRUE;
}


/* Key press handler */
gboolean onKeyPressBelvuTree(GtkWidget *window, GdkEventKey *event, gpointer data)
{
  gboolean handled = FALSE;
  
  BelvuContext *bc = (BelvuContext*)data;
  
  switch (event->keyval)
  {
    case GDK_plus:      /* fall through */
    case GDK_equal:     handled = onZoomBelvuTree(bc, TRUE);    break;
    case GDK_minus:     /* fall through */
    case GDK_underscore: handled = onZoomBelvuTree(bc, FALSE);  break;
      
    default: break;
  };
  
  return handled;
}


/***********************************************************
 *                    Settings dialog                      *
 ***********************************************************/

/* Callback to set the value of a boolean pointer (passed as the user data)
 * when the given toggle button has changed value. */
static gboolean onBoolChangedCallback(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean *value = (gboolean*)data;
  *value = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  return TRUE;
}

/* Callback to set the value of an integer (passed as the user data) when the
 * given text entry (which must contain a double or int as text) has changed value */
static gboolean onDoubleChangedCallback(GtkWidget *entry, const gint responseId, gpointer data)
{
  double *value = (double*)data;
  
  const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
  *value = g_strtod(text, NULL);
  
  return TRUE;
}

/* Utility to create a check button and place it in the given parent box. The
 * button will control the setting of the given boolean value, which will be 
 * updated when the button's parent dialog gets a response. */
static void createCheckButton(const char *mnemonic, gboolean *value, GtkTable *table, const int row, const int col)
{
  g_assert(value);
  
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_table_attach(table, button, col, col + 2, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), *value);
  widgetSetCallbackData(button, onBoolChangedCallback, value);
}

/* Utility to create a text entry box for updating the given integer value.
 * The text entry is added to the given box, which should belong to a dialog.
 * The value will be updated when the dialog gets a response. */
static GtkWidget* createDoubleTextEntry(const char *labelText, double *value, GtkTable *table, int row, int col)
{
  g_assert(value);

  /* Create the label */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  /* Create the text entry */
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col + 1, col + 2, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  widgetSetCallbackData(entry, onDoubleChangedCallback, value);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

  char *defaultInput = convertDoubleToString(*value, 2);
  gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);

  
  const int defaultLen = min(strlen(defaultInput) * 8, 500);
  gtk_widget_set_size_request(entry, defaultLen, -1);

  g_free(defaultInput);
  
  return entry;
}

static GtkWidget* createTreeDisplayOptsButtons(GtkBox *box, 
                                               double *treeScale,
                                               double *lineWidth,
                                               gboolean *showBranchLen, 
                                               gboolean *showOrganism,
					       GtkWidget **treeScaleEntry)
{
  /* Put the display options in a table in a frame */
  GtkContainer *frame = GTK_CONTAINER(gtk_frame_new("Display options"));
  gtk_box_pack_start(box, GTK_WIDGET(frame), FALSE, FALSE, DIALOG_YPAD);
  
  GtkTable *table = GTK_TABLE(gtk_table_new(2, 2, FALSE));
  gtk_container_add(frame, GTK_WIDGET(table));

  *treeScaleEntry = createDoubleTextEntry("Tree scale:", treeScale, table, 0, 0);
  createDoubleTextEntry("Line width:", lineWidth, table, 1, 0);
  createCheckButton("Display branch lengths", showBranchLen, table, 2, 0);
  createCheckButton("Display organism", showOrganism, table, 3, 0);
  
  return GTK_WIDGET(frame);
}


static void createTreeInteractionButtons(GtkBox *box, BelvuPickMode *pickMode)
{
  GtkWidget *frame = gtk_frame_new("Interactions");
  gtk_box_pack_start(box, GTK_WIDGET(frame), FALSE, FALSE, DIALOG_YPAD);
  
  /* Create the tree pick-method drop-down box */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(hbox));

  GtkWidget *label = gtk_label_new("Action when picking a node:");
  gtk_box_pack_start(hbox, label, FALSE, FALSE, DIALOG_XPAD);
  
  GtkComboBox *combo = createComboBox();
  
  GtkTreeIter *iter = NULL;
  BelvuPickMode initMode = *pickMode;
  
  addComboItem(combo, iter, NODESWAP, SWAPstr, initMode);
  addComboItem(combo, iter, NODEROOT, ROOTstr, initMode);

  widgetSetCallbackData(GTK_WIDGET(combo), onComboChanged, pickMode);
  
  gtk_box_pack_start(hbox, GTK_WIDGET(combo), FALSE, FALSE, DIALOG_XPAD);
}


/* Returns true if the given tree setup is distance-correcting or not */
static gboolean correctingDistances(const BelvuBuildMethod buildMethod, const BelvuDistCorr distCorr)
{
  return (buildMethod != UPGMA && distCorr != UNCORR);
}


/* Callback for when the build method has been changed. */
gboolean onBuildMethodChanged(GtkWidget *combo, const gint responseId, gpointer data)
{
  BelvuContext *bc = (BelvuContext*)data;
  
  const int origVal = bc->treeMethod;
  
  onComboChanged(combo, responseId, &bc->treeMethod);
  
  /* This change invalidates the tree, so set the tree head to null to indicate this */
  if (origVal != bc->treeMethod)
    {
      belvuContextSetTree(bc, NULL);
    }
  
  return TRUE;
}


/* Callback for when the distance-correction method has been changed. */
gboolean onDistCorrChanged(GtkWidget *combo, const gint responseId, gpointer data)
{
  BelvuContext *bc = (BelvuContext*)data;
  
  const int origVal = bc->treeDistCorr;

  onComboChanged(combo, responseId, &bc->treeDistCorr);
  
  /* This change invalidates the tree, so set the tree head to null to indicate this */
  if (origVal != bc->treeDistCorr)
    {
      belvuContextSetTree(bc, NULL);
    }
  
  return TRUE;
}


/* Utility to create a standard 2-column combo box and place it in the given
 * table with a label with the given text. */
static GtkComboBox* createComboWithLabel(GtkTable *table, 
                                         const char *labelText, 
                                         const int col, 
                                         const int row,
                                         BlxResponseCallback callbackFunc,
                                         gpointer data)
{
  /* Create the label, and right-align it */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  GtkComboBox *combo = createComboBox();
  gtk_table_attach(table, GTK_WIDGET(combo), col + 1, col + 2, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  widgetSetCallbackData(GTK_WIDGET(combo), callbackFunc, data);
  
  return combo;
} 

static void updateTreeScaleEntry(GtkComboBox *combo, gpointer data)
{
  BuildMethodChangedData *userData = (BuildMethodChangedData*)data;

  int newBuildMethod, newDistCorr;
  onComboChanged(userData->buildMethodCombo, 0, &newBuildMethod);
  onComboChanged(userData->distCorrCombo, 0, &newDistCorr);
  
  if (correctingDistances(origBuildMethod, origDistCorr) == correctingDistances((BelvuBuildMethod)newBuildMethod, (BelvuDistCorr)newDistCorr))
    return;

  /* If we're now distance-correcting and weren't previously (or vice versa)
   * then reset the tree scale to the appropriate default value */
  if (correctingDistances((BelvuBuildMethod)newBuildMethod, (BelvuDistCorr)newDistCorr))
    {
      char *tmpStr = g_strdup_printf("%.2f", DEFAULT_TREE_SCALE_CORR);
      gtk_entry_set_text(GTK_ENTRY(*userData->treeScaleEntry), tmpStr);
      g_free(tmpStr);
    }
  else
    {
      char *tmpStr = g_strdup_printf("%.2f", DEFAULT_TREE_SCALE_NON_CORR);
      gtk_entry_set_text(GTK_ENTRY(*userData->treeScaleEntry), tmpStr);
      g_free(tmpStr);
    }    
  
  origBuildMethod = (BelvuBuildMethod)newBuildMethod;
  origDistCorr = (BelvuDistCorr)newDistCorr;
}


/* Create the drop-down boxes for the build methods on the settings dialog */
static void createTreeBuildMethodButtons(BelvuContext *bc, GtkBox *box, BelvuBuildMethod *buildMethod, BelvuDistCorr *distCorr, GtkWidget **treeScaleEntry)
{
  /* We'll put everything in a table inside a frame */
  GtkWidget *frame = gtk_frame_new("Build methods");
  gtk_box_pack_start(box, GTK_WIDGET(frame), FALSE, FALSE, DIALOG_YPAD);
  
  GtkTable *table = GTK_TABLE(gtk_table_new(2, 2, FALSE));
  gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(table));
  
  /* Create the build-method drop-down box */
  static BuildMethodChangedData data;
  data.bc = bc;
  data.treeScaleEntry = treeScaleEntry;
  
  GtkComboBox *buildMethodCombo = createComboWithLabel(table, "Tree building method:", 0, 0, onBuildMethodChanged, bc);
  data.buildMethodCombo = GTK_WIDGET(buildMethodCombo);
  
  GtkTreeIter *iter = NULL;
  int initMode = *buildMethod;
  addComboItem(buildMethodCombo, iter, NJ, NJstr, initMode);
  addComboItem(buildMethodCombo, iter, UPGMA, UPGMAstr, initMode);
  
  /* Create the distance-correction method drop-down box */
  GtkComboBox *distCorrCombo = createComboWithLabel(table, "Distance correction method:", 0, 1, onDistCorrChanged, bc);
  data.distCorrCombo = GTK_WIDGET(distCorrCombo);
  
  iter = NULL;
  initMode = *distCorr;
  addComboItem(distCorrCombo, iter, UNCORR, UNCORRstr, initMode);
  addComboItem(distCorrCombo, iter, JUKESCANTOR, JUKESCANTORstr, initMode);
  addComboItem(distCorrCombo, iter, KIMURA, KIMURAstr, initMode);
  addComboItem(distCorrCombo, iter, STORMSONN, STORMSONNstr, initMode);
  addComboItem(distCorrCombo, iter, SCOREDIST, SCOREDISTstr, initMode);
  
  origBuildMethod = bc->treeMethod;
  origDistCorr = bc->treeDistCorr;
  g_signal_connect(G_OBJECT(buildMethodCombo), "changed", G_CALLBACK(updateTreeScaleEntry), &data);
  g_signal_connect(G_OBJECT(distCorrCombo), "changed", G_CALLBACK(updateTreeScaleEntry), &data);

}


/* Utility function to create the content for the tree settings dialog */
GtkWidget* createTreeSettingsDialogContent(BelvuContext *bc, 
                                           GtkWidget *dialog, 
                                           const gboolean isMainTree,
                                           double *treeScale,
                                           double *lineWidth,
                                           gboolean *showBranchLen, 
                                           gboolean *showOrganism,
                                           BelvuPickMode *pickMode,
                                           BelvuBuildMethod *buildMethod, 
                                           BelvuDistCorr *distCorr)
{
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), GTK_WIDGET(vbox), FALSE, FALSE, 0);
  
  static GtkWidget *treeScaleEntry = NULL;
  
  /* Only give the option to change the build method for the main belvu tree */
  if (isMainTree)
    createTreeBuildMethodButtons(bc, vbox, buildMethod, distCorr, &treeScaleEntry);

  /* For all trees, give the option to change display settings */
  GtkWidget *content = createTreeDisplayOptsButtons(vbox, treeScale, lineWidth, showBranchLen, showOrganism, &treeScaleEntry);
  createTreeInteractionButtons(vbox, pickMode);

  /* Set the focus on the main content area, because this has widgets that can
   * activate the default response, thereby allowing the user to create the tree
   * very quickly just by pressing Enter. */
  gtk_container_set_focus_child(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), content);
  
  return GTK_WIDGET(vbox);
}


/* Called when the user makes a response on the tree settings dialog. The tree
 * window is passed as the user data. */
void onResponseTreeSettingsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  BelvuContext *bc = (BelvuContext*)data;
  gboolean destroy = TRUE;

  /* Get the toplevel window that this dialog was created from and if if it's
   * not a tree window then use the main belvu tree window. */
  GtkWidget *window = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  GtkWidget *treeWindow = NULL;
  gboolean isMainTree = FALSE;
  
  if (window && stringsEqual(gtk_widget_get_name(window), BELVU_TREE_WINDOW_NAME, TRUE))
    {
      treeWindow = window;
      BelvuTreeProperties *properties = belvuTreeGetProperties(treeWindow);
      isMainTree = properties->isMainTree;
    }
  else 
    {
      treeWindow = bc->belvuTree;
      isMainTree = TRUE;
    }
  
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Call all of the callbacks for each individual widget to update the 
       * properties. Then refresh the window and show the tree. Destroy dialog
       * if successful. */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      belvuTreeUpdateSettings(bc, treeWindow);
      
      if (!treeWindow)
        createAndShowBelvuTree(bc, isMainTree);
      else
        gtk_window_present(GTK_WINDOW(treeWindow));
        
      break;
      
    case GTK_RESPONSE_APPLY:
      /* Never destroy */
      destroy = FALSE;
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      belvuTreeUpdateSettings(bc, treeWindow);
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/* Dialog to allow the user to edit the settings for a tree */
void showTreeSettingsDialog(GtkWidget *window, BelvuContext *bc)
{
  char *title = g_strdup_printf("%sTree Settings", belvuGetTitlePrefix(bc));
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title,
                                       GTK_WINDOW(window), 
                                       (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_MODAL),
                                       GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                       GTK_STOCK_APPLY, GTK_RESPONSE_APPLY,
                                       GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                       NULL);

  g_free(title);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseTreeSettingsDialog), bc);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkWidget *treeWindow = NULL;
  gboolean isMainTree = FALSE;
  
  /* If this is a tree window then check the properties to see if it is the main
   * window. Otherwise, use the main tree from the context. If it's the main window
   * then updating the tree settings will change them in the context; otherwise
   * we use locally-stored properties so that we can change them and only affect
   * the local tree. (Note that the build method properties are only available
   * on the main tree.) */
  double *treeScale = &bc->treeScale;
  double *lineWidth = &bc->treeLineWidth;
  gboolean *showBranchLen = &bc->treeShowBranchlen;
  gboolean *showOrganisms = &bc->treeShowOrganism;
  
  if (window && stringsEqual(gtk_widget_get_name(window), BELVU_TREE_WINDOW_NAME, TRUE))
    {
      treeWindow = window;
      BelvuTreeProperties *properties = belvuTreeGetProperties(treeWindow);
      isMainTree = properties->isMainTree;

      if (!isMainTree)
        {
          /* Use locally-stored properties */
          treeScale = &properties->treeScale;
          lineWidth = &properties->lineWidth;
          showBranchLen = &properties->showBranchLen;
          showOrganisms = &properties->showOrganism;
        }
    }
  else
    {
      treeWindow = bc->belvuTree;
      isMainTree = TRUE;
    }

  createTreeSettingsDialogContent(bc, dialog, isMainTree,
                                  treeScale, lineWidth,
                                  showBranchLen, showOrganisms,
                                  &bc->treePickMode, &bc->treeMethod, &bc->treeDistCorr);
  
  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
}

/***********************************************************
 *                          Sizing                         *
 ***********************************************************/

static void calculateNodeWidth(BelvuTreeProperties *properties, TreeNode *node, const int x)
{
  if (!node)
    return;
  
  const int curX = x + (node->branchlen * treeGetScale(properties) * properties->charWidth);
  
  /* Recurse left and right */
  calculateNodeWidth(properties, node->left, curX);
  calculateNodeWidth(properties, node->right, curX);
  
  /* Check if the end of this text is outside the maximum x position found so far */
  int textWidth = 0;
  textWidth += getTextWidth(properties->treeArea, node->name) + DEFAULT_XPAD;
  textWidth += getTextWidth(properties->treeArea, node->organism) + DEFAULT_XPAD;
  
  int pos = curX + textWidth;
  
  if (pos > properties->bc->maxTreeWidth) 
    properties->bc->maxTreeWidth = pos;
}

static void calculateBelvuTreeBorders(GtkWidget *belvuTree)
{
  if (!belvuTree)
    return;
  
  BelvuTreeProperties *properties = belvuTreeGetProperties(belvuTree);

  if (!properties || !properties->tree || !properties->tree->head)
    return;
  
  /* This loops through all nodes and calculates the max tree width */
  properties->bc->maxTreeWidth = 0;
  calculateNodeWidth(properties, properties->tree->head, properties->treeRect.x);
  
  int treeHeight = (properties->bc->alignArr->len + 7) * properties->charHeight;
  
  if (treeHeight > MAX_PIXMAP_HEIGHT)
    {
      treeHeight = MAX_PIXMAP_HEIGHT;
      g_warning("The tree window is too large and will be clipped.\n");
    }

  int treeWidth = properties->bc->maxTreeWidth;
  
  if (treeWidth > MAX_PIXMAP_WIDTH)
    {
      treeWidth = MAX_PIXMAP_WIDTH;
      g_warning("The tree window is too large and will be clipped.\n");
    }
  
  properties->treeRect.x = DEFAULT_XPAD;
  properties->treeRect.y = DEFAULT_YPAD;
  properties->treeRect.width = treeWidth;
  properties->treeRect.height = treeHeight;
  
  gtk_layout_set_size(GTK_LAYOUT(properties->treeArea), DEFAULT_XPAD * 2 + properties->treeRect.width, 
                      DEFAULT_YPAD * 2 + properties->treeRect.height);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set style properties for the belvu alignment widgets */
static void setBelvuTreeStyle(BelvuContext *bc, GtkWidget *belvuTree)
{
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);
  
  gtk_widget_modify_bg(belvuTree, GTK_STATE_NORMAL, bgColor);
}


/* Create a widget for drawing a tree. Returns the GtkLayout widget that 
 * does the drawing. Places the widget into the given parent box. */
static GtkWidget* createBelvuTreeWidget(BelvuContext *bc, GtkBox *box)
{
  /* Create the drawing area */
  GtkWidget *treeArea = gtk_layout_new(NULL, NULL);
  
  /* Make it scrollable. The scroll-window will be the outermost container and
   * therefore will be the widget that we treat as the "belvu tree". */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_box_pack_start(box, scrollWin, TRUE, TRUE, 0);
  
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(scrollWin), treeArea);

  /* Set the style */
  setBelvuTreeStyle(bc, treeArea);
  
  return treeArea;
}


static void setTreeWindowStyleProperties(GtkWidget *window)
{
  gtk_widget_set_name(window, BELVU_TREE_WINDOW_NAME);
  
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_TREE_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_TREE_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


/* Display an existing tree. Opens in a new window. */
GtkWidget* createBelvuTreeWindow(BelvuContext *bc, Tree *tree, const gboolean isMainTree)
{
  /* Create the window */
  GtkWidget *belvuTree = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setTreeWindowStyleProperties(belvuTree);

  const char *titlePrefix = (isMainTree ? TITLE_MAIN_TREE_PREFIX : TITLE_BOOTSTRAP_TREE_PREFIX);
  const char *titleDesc = (bc->treeMethod == NJ ? TITLE_NJ_TREE_DESCRIPTION : TITLE_UPGMA_TREE_DESCRIPTION);
  
  char *title = g_strdup_printf("%s%s%stree using %s distances of %s", 
                                belvuGetTitlePrefix(bc),
                                titlePrefix, 
                                titleDesc,
                                bc->treeDistString,
                                bc->Title);
  
  gtk_window_set_title(GTK_WINDOW(belvuTree), title);
  g_free(title);
  
  /* We must add all toplevel windows to the list of spawned windows */
  bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, belvuTree);

  /* If this is the main tree, then set it in the context */
  if (isMainTree)
    bc->belvuTree = belvuTree;

  /* Create the context menu and set a callback to show it */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(belvuTree, bc, &actionGroup);
  GtkWidget *contextmenu = createBelvuMenu(belvuTree, "/TreeContextMenu", uiManager);
  
  gtk_widget_add_events(belvuTree, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(belvuTree), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(belvuTree), vbox);
  
  /* Add the alignment section */
  GtkWidget *treeArea = createBelvuTreeWidget(bc, GTK_BOX(vbox));

  /* Set the properties on the tree window */
  belvuTreeCreateProperties(belvuTree, bc, tree, isMainTree, actionGroup, treeArea, bc->treeMethod, bc->treeDistCorr);

  /* Set the initial size (must be called after properties are set) */
  calculateBelvuTreeBorders(belvuTree);

  /* Connect signals. Pass the window as the data. */
  gtk_widget_add_events(treeArea, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(treeArea), "expose-event", G_CALLBACK(onExposeBelvuTree), belvuTree);  
  g_signal_connect(G_OBJECT(treeArea), "button-press-event", G_CALLBACK(onButtonPressBelvuTree), belvuTree);  
  g_signal_connect(G_OBJECT(belvuTree), "key-press-event", G_CALLBACK(onKeyPressBelvuTree), bc);
  
  /* Set the font size to be the same as the main alignment window */
  if (bc->belvuAlignment)
    widgetSetFontSize(belvuTree, GINT_TO_POINTER(pango_font_description_get_size(bc->belvuAlignment->style->font_desc) / PANGO_SCALE));
  
  onBelvuTreeFontSizeChanged(belvuTree);

  gtk_widget_show_all(belvuTree);
  gtk_window_present(GTK_WINDOW(belvuTree));

  return belvuTree;
}


/* Create a new tree and create a tree-window to display it. If isMainTree is 
 * true then set this as the main belvu tree. */
GtkWidget* createAndShowBelvuTree(BelvuContext *bc, const gboolean isMainTree)
{
  GtkWidget *belvuTree = NULL;
  
  /* For the main tree, use the existing tree nodes if any; otherwise create
   * one. For other trees, we always need to create the tree nodes. */
  if (!isMainTree || !bc->mainTree || !bc->mainTree->head)
    {
      separateMarkupLines(bc);
      Tree *tree = treeMake(bc, TRUE, TRUE);
      
      if (isMainTree)
        belvuContextSetTree(bc, &tree);
          
      belvuTree = createBelvuTreeWindow(bc, tree, isMainTree);
      reInsertMarkupLines(bc);
    }
  else
    {
      belvuTree = createBelvuTreeWindow(bc, bc->mainTree, isMainTree);
    }
  
  return belvuTree;
}



