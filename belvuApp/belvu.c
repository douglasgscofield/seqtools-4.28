/*  File: belvu.c
 *  Author: Erik Sonnhammer
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
 * You should have received a copy of the GNU General Public License2 * along with this program; if not, write to the Free Software
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
 * Description: Logic functions for the belvu applications
 *----------------------------------------------------------------------------
 */



/* 

    Pending:

        Internal bootstrapping.
	Rationale for basic bootstrapping:
	     0. struct Subtree {nodep, seqstring}.  nodep points to node in orig tree, seqnamestring a sorted concatenation of seqence ordinal numbers
	     1. Traverse original tree from root to 2-seq nodes, insert Subtree structs at each node into a list (sorted on seqstring)
	     2. For each Boostrap tree {
	         Traverse tree, for each internal node except root: if seqstring found in Subtreelist, increment node's boostrap counter
	     }
	        Traverse bootstrap and original tree
	Rationale for group bootstrapping:
		   For each sequence, store array of sequences in merging order in bootstrap tree.
		   For each node in original tree, check if group of sequences 
		   corresponds to the first n sequences in bootstraptree array.
	Rationale for exact bootstrapping:
		   Traverse bootstrap tree from leaf to root:
		      oritree nodes = {same, different, 0}
		      if issequence(parent->other_child) {
		         if (parent->other_child == boottree(parent->other_child))
			    parent = same
			 else
			    parent = different
		      if isseen(parent->other_child) {
		         if (parent->other_child == boottree(parent->other_child))
			    parent = same
			 else
			    parent = diferent
		      if same, continue to parent

	Clean up the parsing code with str2aln calls.

	read in other trees with bootstraps
	use confidence cutoff in find orthologs
	
        make alignment collapsing easier to use, see above.

        Keyboard edit of one sequence.  How to find the right place ? !!

        Undo edits.

	Would like to abort parsing if two sequence lines accidentally have same name -
	How? - this is used in selex format...

	Koonin ideas:
	    Consensus pattern lines (Default 0% and 100% lines, maybe 90% too).
	    Add/remove any nr of lines cutoffs 0 - 100%.
	    Tool to define groups of residues and corresponding symbol.  These
	    are shown on pattern lines (priority to smaller groups).

	Low priority:

	   "Add sequences with matching segments" for more than one sequence.
	    Will probably be very deleterious to the alignment.

	    Ability to choose both foreground and background color - makes it twice
	    as slow - see "Black/white for printing.

*/


/*  Description of color_by_similarity algorithm in setConsSchemeColors():

    ( Corresponds to summing all pairwise scores)

    for each residue i {
        for each residue j {
	    if (i == j) 
	        score(i) += (count(i)-1)*count(j)*matrix(i,j)
	    else
	        score(i) += count(i)*count(j)*matrix(i,j)
	    }
	}

	if (ignore gaps)
	    n = nresidues(pos)
	else 
	    n = nsequences
		    
	if (n == 1)
	    id = 0.0
	else
	    id = score/(n*(n-1))
    }
		
*/


#include <belvuApp/belvu_.h>
#include <belvuApp/belvuWindow.h>
#include <belvuApp/belvuTree.h>
#include <belvuApp/belvuConsPlot.h>
#include <belvuApp/belvuAlignment.h>

#include <stdarg.h>
/*#include <stdlib.h> / * Needed for RAND_MAX but clashes with other stuff */
#include <sys/types.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>


#define FETCH_PROG_ENV_VAR        "BELVU_FETCH"       /* environment variable used to specify the fetch program */
#define FETCH_URL_ENV_VAR         "BELVU_FETCH_WWW"   /* environment variable used to specify the WWW-fetch URL */
#define DEFAULT_FETCH_PROG        "efetch"            /* default program for fetching sequences */
#define DEFAULT_FETCH_URL         "http://www.sanger.ac.uk/cgi-bin/otter/52/nph-pfetch?request=%s"
//#define DEFAULT_FETCH_URL         "http://www.sanger.ac.uk/cgi-bin/seq-query?%s"  /* default url for fetching sequences in WWW-fetch mode */



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
 *  This converts an ascii char to a 1-based index that can be used in
 *  the BLOSUM matrix - note that you need to subtract 1 from the values
 *  to get a 0-based index for use in BLOSUM62.
 *  
 *  It specifies a 1-based index for the 20 standard amino acids. For any
 *  character that is not a residue, NA is returned.
 *
 * Note: to use with BLOSUM62[], always subtract 1 from the values !!!!
 */

#undef NA
#define NA 23 /* not a residue - eaual to X */

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


/* ASCII-to-binary translation table
 *  Similar to a2b but returns NR (0) for characters that are
 *  not residues, rather than NA. (This is used in our color
 *  arrays, which are 1-based, with the 0 index being reserved
 *  for unknown characters.
 */
#undef NR
#define NR 0

static int n2b[] =
  {
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR, 1,NR, 5, 4, 7,14, 8, 9,10,NR,12,11,13, 3,NR,
    15, 6, 2,16,17,NR,20,18,NR,19,NR,NR,NR,NR,NR,NR,
    NR, 1,NR, 5, 4, 7,14, 8, 9,10,NR,12,11,13, 3,NR,
    15, 6, 2,16,17,NR,20,18,NR,19,NR,NR,NR,NR,NR,NR,

    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
  };


#ifdef OLD_BELVU_CODE

/* Note: this treecpy does not reallocate any strings, assuming they will 
   always exist and never change.  Serious problems await if this is
   not true
*/
static treeNode *treecpy(treeNode *node)
{
    treeNode *newnode;

    if (!node) return 0;

    newnode = g_malloc(sizeof(treeNode));

    newnode->dist      = node->dist;
    newnode->branchlen = node->branchlen;	
    newnode->boot      = node->boot;
    newnode->name      = node->name;
    newnode->organism  = node->organism;
    newnode->aln       = node->aln;
    newnode->box       = node->box;
    newnode->color     = node->color;

    newnode->left      = treecpy(node->left);
    newnode->right     = treecpy(node->right);

    return newnode;
}

/* Not used - the idea was to find balancing point of a tree by walking to
   neighbors and checking if they are better.  However, the best balanced
   node may be beyond less balanced nodes.  (This is because the subtree weights
   are averaged branchlengths and not the sum of the total branchlengths).
*/
static void treeBalanceTreeRecurse(treeNode *node, double *bestbal, treeNode **bestNode)
{
    double newbal, lweight, rweight;

    /* Probe left subtree */

    if (node->left) {
	lweight = treeSize3way(node, node->left) - node->branchlen;
	rweight = treeSize3way(node->left, node) - node->left->branchlen;
	newbal = fabsf(lweight - rweight);
	
	/*
	printf("Subtree weights = %.1f  %.1f\n", lweight, rweight);
	printf("oldbal = %.1f   newbal = %.1f\n", *bestbal, newbal);
	*/

	if (newbal < *bestbal) { /* better balance */

	    g_message("Better balance: %.1f (oldbest = %.1f)\n", newbal, *bestbal);
	    
	    *bestbal = newbal;
	    *bestNode = node;
	    
	    treeBalanceTreeRecurse(node->right, bestbal, bestNode);
	}
    }

    if (node->right) {

	lweight = treeSize3way(node, node->right) - node->branchlen;
	rweight = treeSize3way(node->right, node) - node->right->branchlen;
	newbal = fabsf(lweight - rweight);
	
	/*
	printf("Subtree weights = %.1f  %.1f\n", lweight, rweight);
	printf("oldbal = %.1f   newbal = %.1f\n", *bestbal, newbal);
	*/
	
	if (newbal < *bestbal) { /* better balance */
	    
	    g_message("Better bal: %f\n", newbal);
	    
	    *bestbal = newbal;
	    *bestNode = node;
	    
	    treeBalanceTreeRecurse(node->right, bestbal, bestNode);
	}
    }
}

#endif /* OLD_BELVU_CODE */



/* These values define the defaults for the thresholds when coloring by
 * conservation; the first three are when coloring by %ID and the last 
 * three when coloring by similarity (i.e. BLOSUM62) */
#define DEFAULT_LOW_ID_CUTOFF           0.4
#define DEFAULT_MID_ID_CUTOFF           0.6
#define DEFAULT_MAX_ID_CUTOFF           0.8
#define DEFAULT_LOW_SIM_CUTOFF          0.5
#define DEFAULT_MID_SIM_CUTOFF          1.5
#define DEFAULT_MAX_SIM_CUTOFF          3.0

/* These values define the default color IDs for the conservation colors */
#define DEFAULT_MAX_FG_COLOR            BLACK
#define DEFAULT_MID_FG_COLOR            BLACK
#define DEFAULT_LOW_FG_COLOR            BLACK
#define DEFAULT_MAX_BG_COLOR            CYAN
#define DEFAULT_MID_BG_COLOR            MIDBLUE
#define DEFAULT_LOW_BG_COLOR            LIGHTGRAY  
#define DEFAULT_MAX_FG_PRINT_COLOR      WHITE
#define DEFAULT_MID_FG_PRINT_COLOR      BLACK
#define DEFAULT_LOW_FG_PRINT_COLOR      BLACK
#define DEFAULT_MAX_BG_PRINT_COLOR      BLACK
#define DEFAULT_MID_BG_PRINT_COLOR      GRAY
#define DEFAULT_LOW_BG_PRINT_COLOR      LIGHTGRAY  

/* Global variables */

/* Color names (must be in same order as Color enum) */
static const char *colorNames[NUM_TRUECOLORS] = {
"WHITE", 
"BLACK", 
"LIGHTGRAY", 
"DARKGRAY",
"RED", 
"GREEN", 
"BLUE",
"YELLOW", 
"CYAN", 
"MAGENTA",
"LIGHTRED", 
"LIGHTGREEN", 
"LIGHTBLUE",
"DARKRED", 
"DARKGREEN", 
"DARKBLUE",
"PALERED", 
"PALEGREEN", 
"PALEBLUE",
"PALEYELLOW", 
"PALECYAN", 
"PALEMAGENTA",
"BROWN", 
"ORANGE", 
"PALEORANGE",
"PURPLE", 
"VIOLET", 
"PALEVIOLET",
"GRAY", 
"PALEGRAY",
"CERISE", 
"MIDBLUE"
};


/* Define the old-style acedb colors as hex values that can be parsed by GDK */
static const char* colorTable[NUM_TRUECOLORS]= {
"#ffffff", 	   /* WHITE           */
"#000000",	   /* BLACK           */
"#c8c8c8",	   /* LIGHTGRAY       */
"#646464", 	   /* DARKGRAY        */
"#ff0000",         /* RED             */
"#00ff00",         /* GREEN           */
"#0000ff",	   /* BLUE            */
"#ffff00",	   /* YELLOW          */
"#00ffff",	   /* CYAN            */
"#ff00ff", 	   /* MAGENTA         */
"#ffa0a0",	   /* LIGHTRED        */
"#a0ffa0",	   /* LIGHTGREEN      */
"#a0c8ff",	   /* LIGHTBLUE       */
"#af0000", 	   /* DARKRED         */
"#00af00",	   /* DARKGREEN       */
"#0000af",         /* DARKBLUE        */
"#ffe6d2",	   /* PALERED         */ 
"#d2ffd2", 	   /* PALEGREEN       */
"#d2ebff",	   /* PALEBLUE        */
"#ffffc8",	   /* PALEYELLOW      */
"#c8ffff",	   /* PALECYAN        */
"#ffc8ff",	   /* PALEMAGENTA     */
"#a05000",	   /* BROWN           */
"#ff8000",	   /* ORANGE          */
"#ffdc6e",         /* PALEORANGE      */
"#c000ff",	   /* PURPLE          */
"#c9aaff",	   /* VIOLET          */
"#ebd7ff",	   /* PALEVIOLET      */
"#969696",	   /* GRAY            */
"#ebebeb",	   /* PALEGRAY        */
"#ff0080",	   /* CERISE          */
"#56b2de"	   /* MIDBLUE         */
} ;


/* File format names (must be in same order as BelvuFileFormat enum)  */
static const char *fileFormatNames[BELVU_NUM_FILE_FORMATS] = {
  "Stockholm (Pfam/HMMER)",
  "MSF",
  "Aligned Fasta",
  "Unaligned Fasta"
};


/* Current residue colors */
static int color[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};


/* Saved custom residue colors */
static int customColor[] = {
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};


/* Residue colors */
static int markupColor[] = {
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,

        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG
};


static const char b2a[] = ".ARNDCQEGHILKMFPSTWYV" ;

static int a2b_sean[] =
  {
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR, 1,NR, 2, 3, 4, 5, 6, 7, 8,NR, 9,10,11,12,NR,
    13,14,15,16,17,NR,18,19,NR,20,NR,NR,NR,NR,NR,NR,
    NR, 1,NR, 2, 3, 4, 5, 6, 7, 8,NR, 9,10,11,12,NR,
    13,14,15,16,17,NR,18,19,NR,20,NR,NR,NR,NR,NR,NR,

    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR
  };


 

/* Local function declarations */
static double		   score(char *s1, char *s2, const gboolean penalize_gaps);
static void		   initConservMtx(BelvuContext *bc);
static int		   countResidueFreqs(BelvuContext *bc);
static int                 stripCoordTokens(char *cp, BelvuContext *bc);
int*                       getConsColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground);
static void                rmFinalise(BelvuContext *bc) ;


/***********************************************************
 *		          Sorting			   *
 ***********************************************************/

/* Sort comparision function to sort ALNs by name */
gint alphaorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

  int retval=0;
  
  if (!(retval = strcmp(x->name, y->name)))
    {
      if (x->start == y->start)
	{
	  if (x->end < y->end)
	    retval = -1;
	  else if (x->end > y->end)
	    retval = 1;
	}
      else if (x->start < y->start)
	retval = -1;
      else if (x->start > y->start)
	retval = 1;
    }
  
  /* printf("Comparing %10s %4d %4d with %10s %4d %4d = %d\n", 
   x->name, x->start, x->end, y->name, y->start, y->end, retval); */
  
  return retval;
}


/* Sort comparision function to sort ALNs by organism */
gint organism_order(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

  if (!x->organism && !y->organism)
    return 0;
  else if (!x->organism)
    return -1;
  else if (!y->organism)
    return 1;
  else
    return strcmp(x->organism, y->organism);
}


/* Sort comparision function to sort ALNs by organism */
gint organismorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

  int retval=0;
  char *p1 = strchr(x->name, '_'), 
       *p2 = strchr(y->name, '_');

  if (!p1 && !p2) return alphaorder(xIn, yIn);
  if (!p1) return 1;
  if (!p2) return -1;
  
  if (!(retval = strcmp(p1, p2))) 
    return alphaorder(xIn, yIn);

  return retval;
}


/* Sort comparison function to sort ALNs by score */
gint scoreorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

  if (x->score < y->score)
    return -1;
  else if (x->score > y->score)
    return 1;
  else return 0;
}

static gint scoreorderRev(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

    if (x->score > y->score)
	return -1;
    else if (x->score < y->score)
	return 1;
    else return 0;
}


/* Sort comparison function to sort ALNs by nr */
gint nrorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = *((const ALN**)xIn);
  const ALN *y = *((const ALN**)yIn);

  if (x->nr < y->nr)
    return -1;
  else if (x->nr > y->nr)
    return 1;
  else return 0;
}


void scoreSort(BelvuContext *bc)
{
  if (!bc->displayScores) 
    { 
      g_critical("No scores available.\n");
      return;
    }

  g_array_sort(bc->alignArr, scoreorder);
    
  arrayOrder(bc->alignArr);
}


static void alphaSort(BelvuContext *bc)
{
  g_array_sort(bc->alignArr, alphaorder);
  arrayOrder(bc->alignArr);
}


static void organismSort(BelvuContext *bc)
{
  g_array_sort(bc->alignArr, organismorder);
  arrayOrder(bc->alignArr);
}


void highlightScoreSort(char mode, BelvuContext *bc)
{
  if (!bc->selectedAln) 
    {
      g_critical("Please highlight a sequence first\n");
      return;
    }
  
  if (bc->selectedAln->markup) 
    {
      g_critical("Please do not highlight a markup line\n");
      return;
    }
  
  separateMarkupLines(bc);
  
  /* if (displayScores) {
   if (!(graphQuery("This will erase the current scores, do you want to continue?")))
   return;
   }*/
  
  bc->displayScores = TRUE;
  
  /* Calculate score relative to highlighted sequence */
  int i = 0;
  
  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *curAln = g_array_index(bc->alignArr, ALN*, i);

      if (mode == 'P')
	{
	  curAln->score = score(alnGetSeq(bc->selectedAln), alnGetSeq(curAln), bc->penalize_gaps);
	}
      else if (mode == 'I')
	{
	  curAln->score = identity(alnGetSeq(bc->selectedAln), alnGetSeq(curAln), bc->penalize_gaps);
	}
      
      char *scoreStr = g_strdup_printf("%.1f", curAln->score);
      int len = strlen(scoreStr);
      
      if (len > bc->maxScoreLen)
	bc->maxScoreLen = len;
      
      g_free(scoreStr);
    }
  
  g_array_sort(bc->alignArr, scoreorderRev);
  arrayOrder(bc->alignArr);
  
  reInsertMarkupLines(bc);
  
  bc->alignYStart = 0;
  
  if (bc->belvuAlignment)
    {
      updateHeaderColumnsSize(bc->belvuAlignment);
      calculateDrawingSizes(bc->belvuAlignment);
      belvuAlignmentRedrawAll(bc->belvuAlignment);
    }
}


/* This sorts the alignment list by the current sort mode specified in the context.
 * If the mode is to sort by tree, then this function will create the tree if it
 * does not already exist, and will also display the tree to the user if showTree
 * is true. */
void doSort(BelvuContext *bc, const BelvuSortType sortType, const gboolean showTree)
{
  g_array_sort(bc->alignArr, nrorder);
  
  switch(sortType) 
  {
    case BELVU_SORT_ALPHA :	      alphaSort(bc);                        break;
    case BELVU_SORT_ORGANISM :        organismSort(bc);                     break;
    case BELVU_SORT_SCORE :	      scoreSort(bc);                        break;
    case BELVU_SORT_TREE  :	      treeSort(bc, showTree);               break;
    case BELVU_SORT_CONS  :	      /* sort by nrorder - already done */  break;
    case BELVU_SORT_SIM :             highlightScoreSort('P', bc);          break;
    case BELVU_SORT_ID :              highlightScoreSort('I', bc);          break;
    case BELVU_UNSORTED : break;
    
    default: 
      g_warning("Initial sort order '%d' not recognised.\n", sortType);
      break;
  }
}

/***********************************************************
 *		          Trees				   *
 ***********************************************************/

void setTreeScaleCorr(BelvuContext *bc, const int treeMethod) 
{
  if (treeMethod == UPGMA)
      bc->treeScale = 1.0;
  else if (treeMethod == NJ)
      bc->treeScale = 0.3;
}


void setTreeScale(BelvuContext *bc, const double newScale) 
{
  bc->treeScale = newScale;
}


static int treeOrder(TreeNode *node, const int treeOrderNrIn) 
{
  int treeOrderNr = treeOrderNrIn;
  
  if (node) 
    {
      treeOrderNr = treeOrder(node->left, treeOrderNr);
      
      if (node->aln)
        node->aln->nr = treeOrderNr++;
      
      treeOrderNr = treeOrder(node->right, treeOrderNr);    
    }

  return treeOrderNr;
}


void treeSortBatch(BelvuContext *bc)
{
  separateMarkupLines(bc);
  
  if (!bc->mainTree || !bc->mainTree->head)
    {
      Tree *tree = treeMake(bc, FALSE, TRUE);
      belvuContextSetTree(bc, &tree);
    }
  
  treeOrder(bc->mainTree->head, 1); /* Set nr field according to tree order */
  
  g_array_sort(bc->alignArr, nrorder);
  
  reInsertMarkupLines(bc);
}


void treeTraverse(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *treeNode)) 
{
  if (!node) 
    return;
  
  treeTraverse(bc, node->left, func);
  func(bc, node);
  treeTraverse(bc, node->right, func);
}


/* General purpose routine to convert a string to ALN struct.
   Note: only fields Name, Start, End are filled!
 */
void str2aln(BelvuContext *bc, char *src, ALN *alnp) 
{
  char *tmp = g_strdup(src);
  stripCoordTokens(tmp, bc);

  if (sscanf(tmp, "%s%d%d", alnp->name, &alnp->start, &alnp->end) != 3)
    {
      g_critical("Name to field conversion failed for %s (%s).\n", src, tmp);
      return;
    }
  
  if (strlen(alnp->name) > MAXNAMESIZE)
    g_error("buffer overrun in %s !!!!!!!!!\n", "str2aln") ;
  
  g_free(tmp);
}


/* Sort the alignments by tree order (creates the tree if it does not exist,
 * and also displays it to the user if showTree is true). */
void treeSort(BelvuContext *bc, const gboolean showTree)
{
  treeSortBatch(bc);
  
  if (showTree)
    {
      /* Show the tree window (create it if necessary) */
      if (bc->belvuTree)
        gtk_window_present(GTK_WINDOW(bc->belvuTree));
      else
        createBelvuTreeWindow(bc, bc->mainTree, TRUE);
    }
}


static void treeTraverseLRfirst(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *node)) 
{
  if (!node) 
    return;
  
  treeTraverseLRfirst(bc, node->left, func);
  treeTraverseLRfirst(bc, node->right, func);
  func(bc, node);
}


static void subfamilyTrav(BelvuContext *bc, TreeNode *node) 
{
  static double dist = 0.0;
  static int 
  newgroup = 1,
  groupnr = 0;
  
  if (!node) 
    return;
  
  if (node->name) 
    {
      dist = node->branchlen;
    
      if (newgroup) 
        {
          g_message("\nGroup nr %d:\n", ++groupnr);
          newgroup = 0;
        }
      
      g_message("%s\n", node->name);
    }
  else
    {
      /* internal node */
      dist += node->branchlen;
    }
  
  if ( bc->mksubfamilies_cutoff > (100.0-dist) )
    { 
      newgroup = 1; 
    }
  
  /* printf("abs=%.1f  branch=%.1f\n", dist, node->branchlen); */
}


void mksubfamilies(BelvuContext *bc, double cutoff)
{
  separateMarkupLines(bc);
  
  strcpy(bc->treeMethodString, UPGMAstr);
  bc->treeMethod = UPGMA;
  
  Tree *tree = treeMake(bc, FALSE, TRUE);
  
  treeTraverseLRfirst(bc, tree->head, subfamilyTrav);
}


/***********************************************************
 *		          Alignments			   *
 ***********************************************************/

/* This finalises the addition of an alignment from a fasta file - the given
 * alignment is appended into the array, and the array takes ownership of it.
 * If there is already a sequence with this name in the array then it will not
 * be added, and it will be deleted; therefore, the caller should be wary of
 * maintaining a pointer to it. */
static void readFastaAlnFinalise(BelvuContext *bc, ALN *aln)
{
  if (bc->maxLen) 
    {
      if (alnGetSeqLen(aln) != bc->maxLen) 
        g_error("Differing sequence lengths: %d %d\n", bc->maxLen, alnGetSeqLen(aln));
    }
  else
    {
      bc->maxLen = alnGetSeqLen(aln);
    }

  int ip = 0;
  if (alnArrayFind(bc->alignArr, &aln, &ip, alphaorder))
    {
      g_error("Sequence name occurs more than once: %s%c%d-%d\n", 
              aln->name, bc->saveSeparator, aln->start, aln->end);
      
      g_string_free(aln->sequenceStr, TRUE);
      g_free(aln);
    }
  else
    {
      aln->nr = bc->alignArr->len + 1;
      g_array_append_val(bc->alignArr, aln);
      g_array_sort(bc->alignArr, alphaorder);
    }
}


/* initialise an ALN with empty values */
void initAln(ALN *alnp)
{
  alnp->name[0] = '\0';
  alnp->start = 0;
  alnp->end = 0;
  alnp->sequenceStr = NULL;
  alnp->nr = 0;
  alnp->fetch[0] = '\0';
  alnp->score = 0.0;
  alnp->color = 0;
  alnp->markup = 0;
  alnp->hide = FALSE;
  alnp->nocolor = FALSE;
  alnp->organism = NULL; 
  alnp->startColIdx = 0;
}


ALN* createEmptyAln()
{
  ALN *aln = (ALN*)g_malloc(sizeof *aln);
  initAln(aln);
  return aln;
}


/* Read in fasta sequences from a file and create a sequence in the given
 * alignments array for each of the fasta sequences. */
static void readFastaAln(BelvuContext *bc, FILE *pipe)
{
  char line[MAXLENGTH+1];

  ALN *currentAln = NULL;

  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe))
	break;

      char *cp = strchr(line, '\n');
      
      /* Cut off the newline char at the end, if it has one */
      if (cp)
	*cp = 0;

      if (*line == '>')
        {
          if (currentAln)
            {
              /* Finish off the previous sequence */
              readFastaAlnFinalise(bc, currentAln);
              currentAln = NULL;
            }

          /* Create a new sequence */
          currentAln = createEmptyAln();
          currentAln->sequenceStr = g_string_new(NULL);
          
          /* Parse the new line. Note that this resets the ALN struct for us. */
          parseMulLine(bc, line + 1, currentAln);
        }
      else if (currentAln)
        {
          /* Part-way through reading a sequnce; append the current line to it */
          g_string_append(currentAln->sequenceStr, line);
        }
    }
  
  if (currentAln)
    {
      readFastaAlnFinalise(bc, currentAln);
    }

  bc->saveFormat = BELVU_FILE_ALIGNED_FASTA;
  
  return ;
}


/* Light copy of an ALN struct - the sequence and organism strings are pointers to the
 * original values; they are not duplicated */
void alncpy(ALN *dest, ALN *src)
{
  strncpy(dest->name, src->name, MAXNAMESIZE);
  dest->start = src->start;
  dest->end = src->end;
  /* dest->sequenceStr = src->sequenceStr ? g_string_new(src->sequenceStr->str) : NULL; */
  dest->sequenceStr = src->sequenceStr;
  dest->nr = src->nr;			
  strncpy(dest->fetch, src->fetch, MAXNAMESIZE+10);
  dest->score = src->score;
  dest->color = src->color;
  dest->markup = src->markup;
  dest->hide = src->hide;
  dest->nocolor = src->nocolor;
  dest->organism = src->organism;
  dest->startColIdx = src->startColIdx;
}


/***********************************************************
 *		          Parsing			   *
 ***********************************************************/

/* Convenience routine for converting "name/start-end" to "name start end".
 Used by parsing routines.
 
 Return 1 if both tokens were found, otherwise 0.
 */
static int stripCoordTokens(char *cp, BelvuContext *bc)
{
  char *end = cp;
  
  while (*end && !isspace(*end)) end++;
  
  if ((cp = strchr(cp, bc->saveSeparator)) && cp < end) {
    *cp = ' ';
    if ((cp = strchr(cp, '-')) && cp < end) {
      *cp = ' ';
      bc->IN_FORMAT = MUL;
      return 1;
    }
  }

  bc->IN_FORMAT = RAW;
  return 0;
}


/*
 Parse name, start and end of a Mul format line
 
 Convenience routine, part of readMul and other parsers
 */
void parseMulLine(BelvuContext *bc, char *line, ALN *aln)
{
  char line2[MAXLENGTH+1], *cp=line2, *cq, GRfeat[MAXNAMESIZE+1];
  GRfeat[0] = 0;
  
  strncpy(cp, line, MAXLENGTH);
  
  if (!strncmp(cp, "#=GC", 4))
    {
      aln->markup = GC;
      cp += 5;
    }

  if (!strncmp(cp, "#=RF", 4)) 
    {
      aln->markup = GC;
    } 

  if (!strncmp(cp, "#=GR", 4)) 
    {
      aln->markup = GR;
      cp += 5;
    }
  
  if (bc->stripCoordTokensOn) 
    stripCoordTokens(cp, bc);
  
  /* Name */
  strncpy(aln->name, cp, MAXNAMESIZE);
  aln->name[MAXNAMESIZE] = 0;

  if ((cq = strchr(aln->name, ' '))) 
    *cq = 0;
  
  /* Add Start and End coords */
  if (bc->stripCoordTokensOn && (bc->IN_FORMAT != RAW) )
    sscanf(strchr(cp, ' ') + 1, "%d%d", &aln->start, &aln->end);
  
  if (aln->markup == GR) 
    {
      /* Add GR markup names to name 
       
         #=GR O83071 192 246 SA 999887756453524252..55152525....36463774777.....948472782969685958
         ->name = O83071_SA
      */
      if (bc->IN_FORMAT == MUL) 
        {
          sscanf(cp + strlen(aln->name), "%d%d%s", &aln->start, &aln->end, GRfeat);
        }
      else
        {
          sscanf(cp + strlen(aln->name), "%s", GRfeat);
        }
      
      if (strlen(aln->name)+strlen(GRfeat)+2 > MAXNAMESIZE)
        g_critical("Too long name or/and feature name\n");

      strcat(aln->name, " ");
      strncat(aln->name, GRfeat, MAXNAMESIZE-strlen(aln->name));
      
      /* printf("%s, %d chars\n", aln->name, strlen(aln->name)); fflush(stdout); */
    }
}


/***********************************************************
 *		          Colours			   *
 ***********************************************************/

/* Set the default colors of organisms to something somewhat intelligent */
void setOrganismColors(GArray *organismArr) 
{
    static int treeColors[16] = {
      RED, 
      BLUE,
      DARKGREEN, 
      ORANGE, 
      MAGENTA,
      BROWN, 
      PURPLE, 
      CYAN, 
      VIOLET, 
      MIDBLUE,
      CERISE, 
      LIGHTBLUE,
      DARKRED, 
      GREEN, 
      DARKBLUE,
      GRAY
  };
  
  int i = 0;
  int color = 0;
  
  for (i = 0; i < (int)organismArr->len; ++i) 
    {
      color = treeColors[i % 16];
      /*if (i > 15) color = treeColors[15];*/
     
      g_array_index(organismArr, ALN*, i)->color = color;
    }
}



/* Return the markup color for the given char */
int getMarkupColor(const char inputChar)
{
  return markupColor[(unsigned char)(inputChar)];
}

/* Return the conservation color for the given char at the given index */
int getConservColor(BelvuContext *bc, const char inputChar, const int idx)
{
  return bc->colorMap[n2b[(unsigned char)(inputChar)]][idx];
}

/* Return the color from the colors array for the given char */
int getColor(const char inputChar)
{
  return color[(unsigned char)(inputChar)];
}

/* Set the color in the colors array for the given char */
void setColor(const char inputChar, const int colorNum)
{
  color[(unsigned char)(toupper(inputChar))] = colorNum;
  color[(unsigned char)(tolower(inputChar))] = colorNum;
}

/* Returns the 'color' array */
int* getColorArray()
{
  return color;
}

/* Returns the 'markupColor' array */
int* getMarkupColorArray()
{
  return markupColor;
}


/* Convert one of the old acedb-style color numbers to a hex string */
static const char* convertColorNumToStr(const int colorNum)
{
  const char *result = colorTable[WHITE];
  
  if (colorNum >= 0 && colorNum < NUM_TRUECOLORS)
    result = colorTable[colorNum];
  
  return result;
}


/* Convert an old-style ACEDB color number to a GdkColor */
void convertColorNumToGdkColor(const int colorNum, const gboolean isSelected, GdkColor *result)
{
  const char *colorStr = convertColorNumToStr(colorNum);
  getColorFromString(colorStr, result, NULL);
  
  /* If an item is selected, we use a slightly different shade of the same color. */
  if (isSelected)
    getSelectionColor(result, result);
}


//static void colorCons(BelvuContext *bc)
//{
//  setConsSchemeColors(bc);
//  
//  //  menuSetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);
//  //  menuUnsetFlags(menuItem(colorMenu, printColorsStr), MENUFLAG_DISABLED);
//  //  menuUnsetFlags(menuItem(colorMenu, ignoreGapsStr), MENUFLAG_DISABLED);
//  
//  bc->colorByResIdOn = FALSE;
//  //  belvuRedraw();
//}


/* Returns true if we're coloring by conservation */
gboolean colorByConservation(BelvuContext *bc)
{
  return (bc->schemeType == BELVU_SCHEME_TYPE_CONS);
}

/* Returns true if we're coloring by residue */
gboolean colorByResidue(BelvuContext *bc)
{
  return (bc->schemeType == BELVU_SCHEME_TYPE_RESIDUE);
}

/* Returns true if we're coloring by similarity */
gboolean colorBySimilarity(BelvuContext *bc)
{
  return (colorByConservation(bc) && bc->consScheme == BELVU_SCHEME_BLOSUM);
}


/* Returns true if we should only color residues above a set threshold */
gboolean colorByResId(BelvuContext *bc)
{
  /* This only applies in color-by-residue mode */
  return (colorByResidue(bc) && bc->colorByResIdOn);
}


/* Utility to clear all residue colors to white */
static void clearResidueColors(BelvuContext *bc)
{
  color['A'] = color['a'] = WHITE;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = WHITE;
  color['D'] = color['d'] = WHITE;
  color['E'] = color['e'] = WHITE;
  color['F'] = color['f'] = WHITE;
  color['G'] = color['g'] = WHITE;
  color['H'] = color['h'] = WHITE;
  color['I'] = color['i'] = WHITE;
  color['J'] = color['j'] = WHITE;
  color['K'] = color['k'] = WHITE;
  color['L'] = color['l'] = WHITE;
  color['M'] = color['m'] = WHITE;
  color['N'] = color['n'] = WHITE;
  color['O'] = color['o'] = WHITE;
  color['P'] = color['p'] = WHITE;
  color['Q'] = color['q'] = WHITE;
  color['R'] = color['r'] = WHITE;
  color['S'] = color['s'] = WHITE;
  color['T'] = color['t'] = WHITE;
  color['V'] = color['v'] = WHITE;
  color['U'] = color['u'] = WHITE;
  color['W'] = color['w'] = WHITE;
  color['X'] = color['x'] = WHITE;
  color['Y'] = color['y'] = WHITE;
  color['Z'] = color['z'] = NOCOLOR;
}


/* Utility to clear all custom colors to white. Should be called during initialisation. */
void initCustomColors()
{
  customColor['A'] = customColor['a'] = WHITE;
  customColor['B'] = customColor['b'] = NOCOLOR;
  customColor['C'] = customColor['c'] = WHITE;
  customColor['D'] = customColor['d'] = WHITE;
  customColor['E'] = customColor['e'] = WHITE;
  customColor['F'] = customColor['f'] = WHITE;
  customColor['G'] = customColor['g'] = WHITE;
  customColor['H'] = customColor['h'] = WHITE;
  customColor['I'] = customColor['i'] = WHITE;
  customColor['J'] = customColor['j'] = WHITE;
  customColor['K'] = customColor['k'] = WHITE;
  customColor['L'] = customColor['l'] = WHITE;
  customColor['M'] = customColor['m'] = WHITE;
  customColor['N'] = customColor['n'] = WHITE;
  customColor['O'] = customColor['o'] = WHITE;
  customColor['P'] = customColor['p'] = WHITE;
  customColor['Q'] = customColor['q'] = WHITE;
  customColor['R'] = customColor['r'] = WHITE;
  customColor['S'] = customColor['s'] = WHITE;
  customColor['T'] = customColor['t'] = WHITE;
  customColor['V'] = customColor['v'] = WHITE;
  customColor['U'] = customColor['u'] = WHITE;
  customColor['W'] = customColor['w'] = WHITE;
  customColor['X'] = customColor['x'] = WHITE;
  customColor['Y'] = customColor['y'] = WHITE;
  customColor['Z'] = customColor['z'] = NOCOLOR;
}


static void colorSchemeCGP(BelvuContext *bc)
{
  clearResidueColors(bc);
  
  color['C'] = color['c'] = CYAN;
  color['G'] = color['g'] = RED;
  color['P'] = color['p'] = GREEN;
}


static void colorSchemeCGPH(BelvuContext *bc)
{
  clearResidueColors(bc);
  
  color['C'] = color['c'] = CYAN;
  color['G'] = color['g'] = RED;
  color['P'] = color['p'] = GREEN;
  color['H'] = color['p'] = YELLOW;
}


static void colorSchemeEmpty(BelvuContext *bc)
{
  clearResidueColors(bc);
}


static void colorSchemeErik(BelvuContext *bc)
{
  /* Erik's favorite colours:
   
   C        - MIDBLUE
   GP       - CYAN
   HKR      - GREEN
   AFILMVWY - YELLOW
   BDENQSTZ - LIGHTRED
   */
  
  color['A'] = color['a'] = YELLOW;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = MIDBLUE;
  color['D'] = color['d'] = LIGHTRED;
  color['E'] = color['e'] = LIGHTRED;
  color['F'] = color['f'] = YELLOW;
  color['G'] = color['g'] = CYAN;
  color['H'] = color['h'] = GREEN;
  color['I'] = color['i'] = YELLOW;
  color['K'] = color['k'] = GREEN;
  color['L'] = color['l'] = YELLOW;
  color['M'] = color['m'] = YELLOW;
  color['N'] = color['n'] = LIGHTRED;
  color['P'] = color['p'] = CYAN;
  color['Q'] = color['q'] = LIGHTRED;
  color['R'] = color['r'] = GREEN;
  color['S'] = color['s'] = LIGHTRED;
  color['T'] = color['t'] = LIGHTRED;
  color['V'] = color['v'] = YELLOW;
  color['W'] = color['w'] = YELLOW;
  color['Y'] = color['y'] = YELLOW;
  color['Z'] = color['z'] = NOCOLOR;
}


static void colorSchemeGibson(BelvuContext *bc)
{
  /* Colour scheme by Gibson et. al (1994) TIBS 19:349-353
   
   Listed in Figure 1:
   
   
   Gibson      AA        Here
   
   orange      G         ORANGE (16-colours: LIGHTRED)
   yellow      P         YELLOW
   blue        ACFILMVW  MIDBLUE
   light blue  Y         LIGHTBLUE (16-colours: CYAN)
   green       NQST      GREEN
   purple      DE        PURPLE (16-colours: MAGENTA)
   red         RK        RED
   pink        H         LIGHTRED (16-colours: DARKRED)
   */
  
  color['A'] = color['a'] = MIDBLUE;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = MIDBLUE;
  color['D'] = color['d'] = PURPLE;
  color['E'] = color['e'] = PURPLE;
  color['F'] = color['f'] = MIDBLUE;
  color['G'] = color['g'] = ORANGE;
  color['H'] = color['h'] = LIGHTRED;
  color['I'] = color['i'] = MIDBLUE;
  color['K'] = color['k'] = RED;
  color['L'] = color['l'] = MIDBLUE;
  color['M'] = color['m'] = MIDBLUE;
  color['N'] = color['n'] = GREEN;
  color['P'] = color['p'] = YELLOW;
  color['Q'] = color['q'] = GREEN;
  color['R'] = color['r'] = RED;
  color['S'] = color['s'] = GREEN;
  color['T'] = color['t'] = GREEN;
  color['V'] = color['v'] = MIDBLUE;
  color['W'] = color['w'] = MIDBLUE;
  color['Y'] = color['y'] = LIGHTBLUE;
  color['Z'] = color['z'] = NOCOLOR;
}


/* Save the current residue colors as the 'custom' color scheme */
void saveCustomColors(BelvuContext *bc)
{
  bc->haveCustomColors = TRUE;
  
  customColor['A'] = customColor['a'] = color['a'];
  customColor['B'] = customColor['b'] = color['b'];
  customColor['C'] = customColor['c'] = color['c'];
  customColor['D'] = customColor['d'] = color['d'];
  customColor['E'] = customColor['e'] = color['e'];
  customColor['F'] = customColor['f'] = color['f'];
  customColor['G'] = customColor['g'] = color['g'];
  customColor['H'] = customColor['h'] = color['h'];
  customColor['I'] = customColor['i'] = color['i'];
  customColor['K'] = customColor['k'] = color['k'];
  customColor['L'] = customColor['l'] = color['l'];
  customColor['M'] = customColor['m'] = color['m'];
  customColor['N'] = customColor['n'] = color['n'];
  customColor['P'] = customColor['p'] = color['p'];
  customColor['Q'] = customColor['q'] = color['q'];
  customColor['R'] = customColor['r'] = color['r'];
  customColor['S'] = customColor['s'] = color['s'];
  customColor['T'] = customColor['t'] = color['t'];
  customColor['V'] = customColor['v'] = color['v'];
  customColor['W'] = customColor['w'] = color['w'];
  customColor['Y'] = customColor['y'] = color['y'];
  customColor['Z'] = customColor['z'] = color['z'];
  
}


/* Set the current colors the the last-saved 'custom' color scheme */
static void colorSchemeCustom(BelvuContext *bc)
{
  color['A'] = color['a'] = customColor['a'];
  color['B'] = color['b'] = customColor['b'];
  color['C'] = color['c'] = customColor['c'];
  color['D'] = color['d'] = customColor['d'];
  color['E'] = color['e'] = customColor['e'];
  color['F'] = color['f'] = customColor['f'];
  color['G'] = color['g'] = customColor['g'];
  color['H'] = color['h'] = customColor['h'];
  color['I'] = color['i'] = customColor['i'];
  color['K'] = color['k'] = customColor['k'];
  color['L'] = color['l'] = customColor['l'];
  color['M'] = color['m'] = customColor['m'];
  color['N'] = color['n'] = customColor['n'];
  color['P'] = color['p'] = customColor['p'];
  color['Q'] = color['q'] = customColor['q'];
  color['R'] = color['r'] = customColor['r'];
  color['S'] = color['s'] = customColor['s'];
  color['T'] = color['t'] = customColor['t'];
  color['V'] = color['v'] = customColor['v'];
  color['W'] = color['w'] = customColor['w'];
  color['Y'] = color['y'] = customColor['y'];
  color['Z'] = color['z'] = customColor['z'];
}


/* This is called when the color scheme has been set to 'by residue'. It updates
 * the colors according to the active color scheme. */
void setResidueSchemeColors(BelvuContext *bc)
{
  switch (bc->residueScheme)
    {
      case BELVU_SCHEME_ERIK:     colorSchemeErik(bc);	    break;
      case BELVU_SCHEME_GIBSON:   colorSchemeGibson(bc);    break;
      case BELVU_SCHEME_CGP:      colorSchemeCGP(bc);	    break;
      case BELVU_SCHEME_CGPH:     colorSchemeCGPH(bc);	    break;
      case BELVU_SCHEME_NONE:     colorSchemeEmpty(bc);	    break;
      case BELVU_SCHEME_CUSTOM:   colorSchemeCustom(bc);    break;
      
      default: 
	g_warning("Program error: unrecognised color scheme '%d'.\n", bc->residueScheme);
	break;
    }
}


/* This should be called when the color scheme has changed. It updates the
 * colors according to the active scheme. */
void updateSchemeColors(BelvuContext *bc)
{
  /* Set the color scheme if coloring by conservation or if applying a 
   * threshold when coloring by residue */
  if (colorByConservation(bc) || colorByResId(bc))
    setConsSchemeColors(bc);
}


/* Return 1 if c1 has priority over c2, 0 otherwise */
static int colorPriority(BelvuContext *bc, int c1, int c2) 
{
  if (c2 == WHITE) 
    return 1;
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_MAX, FALSE))
    return 0;
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_LOW, FALSE)) 
    {
      if (c1 ==*getConsColor(bc, CONS_LEVEL_LOW, FALSE)) 
        return 0;
      else 
        return 1;
    }
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_MID, FALSE)) 
    {
      if (c1 == *getConsColor(bc, CONS_LEVEL_MAX, FALSE)) 
        return 1;
      else
        return 0;
    }
  
  g_critical("Program error: invalid background colour '%s' when calculating color priority.\n", colorNames[c2]);
  
  return 0 ;
}


/* This is called when the color scheme type has been changed to 'by conservation'. It updates
 * the colors according to the active color scheme. */
void setConsSchemeColors(BelvuContext *bc)
{
  int i, j, k, l, colornr, simCount, n;
  double id, maxid;
  
  if (!bc->conservCount) 
    initConservMtx(bc);
  
  int totalNumSeqs = countResidueFreqs(bc);
  
  for (i = 0; i < bc->maxLen; ++i)
    {
      for (k = 1; k < 21; ++k)
	{
	  bc->colorMap[k][i] = WHITE;
	}
    }
  
  for (i = 0; i < bc->maxLen; ++i) 
    {
      maxid = -100.0;
    
      for (k = 1; k < 21; k++) 
	{
	  if (colorBySimilarity(bc)) 
	    {
              /* Convert counts to similarity counts */
              simCount = 0;
              for (j = 1; j < 21; j++) 
                {
                  /* Get the blosum comparison score of the two residues */
                  int score_k_vs_j = BLOSUM62[j-1][k-1];
                  
                  /* This comparison score applies for each occurance of k vs
                   * each occurance of j, e.g. if there are 3 occurances of k 
                   * and 2 occurances of j, we have:
                   *   k1 vs j1 = score_k_vs_j
                   *   k1 vs j2 = score_k_vs_j
                   *   k2 vs j1 = score_k_vs_j
                   *   k2 vs j2 = score_k_vs_j
                   *   k3 vs j1 = score_k_vs_j
                   *   k4 vs j2 = score_k_vs_j
                   * 
                   * i.e. score_k_vs_j occurs (count_k * count_j) times.
                   */
                  int count_k = bc->conservCount[k][i];
                  int count_j = bc->conservCount[j][i];

                  /* Don't compare the same amino acid against itself, i.e. if
                   * there are three occurances of k then we compare:
                   *   k1 vs k2 
                   *   k1 vs k3 
                   * 
                   * but NOT
                   *   k1 vs k1
                   *
                   * so in this case score_k_vs_j occurs (count_k * (count_k - 1)) times.
                   */
                  if (j == k) 
                    --count_k;

                  simCount += count_k * count_j * score_k_vs_j;
                }
	    
              if (bc->ignoreGapsOn) 
                n = bc->conservResidues[i]; /* total number of residues in this column */
              else 
                n = totalNumSeqs;  /* total number of sequences */
              
              if (n < 2)
                {
                  id = 0.0;
                }
              else 
                {
                  /* Divide the similarity count by the total number of comparisons 
                   * made for each column; we made n * (n - 1) comparisons because
                   * we compared each of the n residues in the column to each other
                   * residue in the column except itself. */
                  id = (double)simCount / (n * (n-1));
                }
              
              /* printf("%d, %c:  simCount= %d, id= %.2f\n", i, b2a[k], simCount, id); */
	    
              /* Colour this residue if it is above the %ID threshold */
              if (id > bc->lowSimCutoff) 
                {
                  /* Choose the colour based on the 3 specified levels */
                  if (id > bc->maxSimCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
                  else if (id > bc->midSimCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
                  else
                    colornr = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
                  
                  /* Set the colour for this residue, unless it has already been
                   * given a colour with a higher priority than this one (i.e. it
                   * has already been marked as more conserved) */
                  if (colorPriority(bc, colornr, bc->colorMap[k][i]))
                    bc->colorMap[k][i] = colornr;
	      
                  /* Color all similar residues too; that is, any residue that has
                   * a positive blosum score when compared to the current residue
                   * should be coloured with same level of conservation in this 
                   * column; again, we only set the colour if it doesn't already 
                   * have a higher priority colour set on it. */
                  for (l = 1; l < 21; l++) 
                    {
                      if (BLOSUM62[k-1][l-1] > 0 && colorPriority(bc, colornr, bc->colorMap[l][i])) 
                        {
                          /*printf("%d: %c -> %c\n", i, b2a[k], b2a[l]);*/
                          bc->colorMap[l][i] = colornr;
                        }
                    }
                }
	    }
	  else 
	    {
              /* We are colouring by %ID */
              
              /* First, get the %ID; this is the count of this residue divided
               * by the total number of residues (or the total number of sequences,
               * if we are including gaps). 
               * If ignoring gaps but there is only one residue in this column
               * then the ID takes into account the total number of sequences;
               * I'm not sure why - perhaps because there are no other residues
               * to compare it to; it seems a bit inconsistent, though. */
              if (bc->ignoreGapsOn && bc->conservResidues[i] != 1)
                id = (double)bc->conservCount[k][i]/bc->conservResidues[i];
              else
                id = (double)bc->conservCount[k][i]/totalNumSeqs;
              
              if (colorByResId(bc)) 
                {
                  /* We're colouring by residue type, but only colouring the residues
                   * if their %ID is above the set threshold */
                  if (id * 100.0 > bc->colorByResIdCutoff)
                    bc->colorMap[k][i] = color[(unsigned char)(b2a[k])];
                  else
                    bc->colorMap[k][i] = WHITE;
                }
              else if (id > bc->lowIdCutoff) 
                {
                  /* We're colouring by conservation, using the %ID to determine
                   * the colour according to the three thresholds: */
                  if (id > bc->maxIdCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
                  else if (id > bc->midIdCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
                  else
                    colornr = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
                  
                  /* Set the colour in the array (to do: should this use 
                   * colorPriority to check if it's already been set? At the moment
                   * it overrides any previous (possibly better) colour set 
                   * from a similar residue's result)  */
                  bc->colorMap[k][i] = colornr;
                  
                  if (bc->consScheme == BELVU_SCHEME_ID_BLOSUM) 
                    {
                      /* Colour all similar residues too; that is, any residues
                       * that have a positive blosum score when compared to the 
                       * current residue should be given the same colour in this
                       * column (unless a higher priority colour has already been
                       * set). */
                      for (l = 1; l < 21; l++) 
                        {
                          if (BLOSUM62[k-1][l-1] > 0 && colorPriority(bc, colornr, bc->colorMap[l][i])) 
                            {
                              /*printf("%d: %c -> %c\n", i, b2a[k], b2a[l]);*/
                              bc->colorMap[l][i] = colornr;
                            }
                        }
                    }
                }
	    }
	
	  if (id > maxid) 
	    {
	      maxid = id;
	    }
	}
      
      bc->conservation[i] = maxid;
    }
}


void initMarkupColors(void)
{
  markupColor['0'] = DARKBLUE;
  markupColor['1'] = BLUE;
  markupColor['2'] = MIDBLUE;
  markupColor['3'] = LIGHTBLUE;
  markupColor['4'] = VIOLET;
  markupColor['5'] = PALEBLUE;
  markupColor['6'] = PALECYAN;
  markupColor['7'] = CYAN;
  markupColor['8'] = CYAN;
  markupColor['9'] = CYAN;
  markupColor['A'] = markupColor['a'] = WHITE;
  markupColor['B'] = markupColor['b'] = RED;
  markupColor['C'] = markupColor['c'] = PALEYELLOW;
  markupColor['D'] = markupColor['d'] = WHITE;
  markupColor['E'] = markupColor['e'] = RED;
  markupColor['F'] = markupColor['f'] = WHITE;
  markupColor['G'] = markupColor['g'] = DARKGREEN;
  markupColor['H'] = markupColor['h'] = DARKGREEN;
  markupColor['I'] = markupColor['i'] = DARKGREEN;
  markupColor['J'] = markupColor['j'] = WHITE;
  markupColor['K'] = markupColor['k'] = WHITE;
  markupColor['L'] = markupColor['l'] = WHITE;
  markupColor['M'] = markupColor['m'] = WHITE;
  markupColor['N'] = markupColor['n'] = WHITE;
  markupColor['O'] = markupColor['o'] = WHITE;
  markupColor['P'] = markupColor['p'] = WHITE;
  markupColor['Q'] = markupColor['q'] = WHITE;
  markupColor['R'] = markupColor['r'] = WHITE;
  markupColor['S'] = markupColor['s'] = YELLOW;
  markupColor['T'] = markupColor['t'] = YELLOW;
  markupColor['V'] = markupColor['v'] = WHITE;
  markupColor['U'] = markupColor['u'] = WHITE;
  markupColor['W'] = markupColor['w'] = WHITE;
  markupColor['X'] = markupColor['x'] = WHITE;
  markupColor['Y'] = markupColor['y'] = WHITE;
  markupColor['Z'] = markupColor['z'] = NOCOLOR;
}


/* Save the current color-by-residue color scheme */
void saveResidueColorScheme(BelvuContext *bc, FILE *fil)
{
  if (!fil) 
    return;
  
  int i = 1;
  for (i = 1; i < 21; i++)
    {
      fprintf(fil, "%c %s\n", b2a[i], colorNames[color[(unsigned char)(b2a[i])]]);
    }
  
  fclose(fil);
}


/* This reads in residue colors from the given file and places them into
 * the given array, which is typically one of the active color scheme arrays
 * 'color' or 'markupColor'. 
 * If storeCustomColors is true, then it also saves the colors to the 'customColor'
 * array so that they can be retrieved later (because 'color' gets overwritten
 * with the current colors every time the user toggles between different color
 * schemes). */
void readResidueColorScheme(BelvuContext *bc, FILE *fil, int *colorarr, const gboolean storeCustomColors)
{
  if (!fil) 
    return;

  char *cp=NULL, line[MAXLINE+1], setColor[MAXLINE+1];
  unsigned char c ;
  int i=0, colornr=0;
  
  while (!feof(fil)) 
    {
      if (!fgets (line, MAXLINE, fil)) 
	break;
    
      /* remove newline */
      if ((cp = strchr(line, '\n'))) 
        *cp = 0 ;
    
      /* Parse color of organism in tree 
         Format:  #=OS BLUE D. melanogaster*/
      if (!strncmp(line, "#=OS ", 5)) 
        {
          cp = line+5;
          sscanf(cp, "%s", setColor);
          
          for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
            {
              if (!strcasecmp(colorNames[i], setColor)) 
                colornr = i;
            }
          
          if (colornr == -1) 
            {
              g_warning("Unrecognized color: %s, using black instead.\n", setColor);
              colornr = BLACK;
            }
          
          while(*cp == ' ') cp++;
          while(*cp != ' ') cp++;
          while(*cp == ' ') cp++;
          
          /* Find organism and set its colour */
          ALN aln;
          initAln(&aln);
          aln.organism = cp;
          
          int ip = 0;
          if (!alnArrayFind(bc->organismArr, &aln, &ip, organism_order))
            g_critical("Cannot find organism \"%s\", specified in color code file. Hope that's ok\n", aln.organism);
          else
            g_array_index(bc->organismArr, ALN*, ip)->color = colornr;
        }
      
      /* Ignore comments */
      if (*line == '#') 
        continue;
      
      /* Parse character colours */
      if (sscanf(line, "%c%s", &c, setColor) == 2) 
        {
          c = toupper(c);
          for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
            {
              if (!strcasecmp(colorNames[i], setColor)) 
                colornr = i;
            }
          
          if (colornr == -1) 
            {
              g_warning("Unrecognized color: %s\n", setColor);
              colornr = 0;
            }
          
          colorarr[(unsigned char)(c)] = colornr;
          
          if (c > 64 && c <= 96)
            colorarr[(unsigned char)(c+32)] = colorarr[(unsigned char)(c)];
          else if (c > 96 && c <= 128)
            colorarr[(unsigned char)(c-32)] = colorarr[(unsigned char)(c)];
        }
    }
  
  fclose(fil);
  
  /* Store the custom colors */
  saveCustomColors(bc);
}


/* Utility function to get the color display name for a color number */
const char* getColorNumName(const int colorNum)
{
  g_assert(colorNum < NUM_TRUECOLORS);
  return colorNames[colorNum];
}

/* Utility function to get the file format name for a format id */
const char* getFileFormatString(const int formatId)
{
  g_assert(formatId < BELVU_NUM_FILE_FORMATS);
  return fileFormatNames[formatId];
}


int* getConsPrintColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground)
{
  int *result = &bc->lowbgPrintColor;

  switch (consLevel)
    {
      case CONS_LEVEL_MAX:
        result = foreground ? &bc->maxfgPrintColor : &bc->maxbgPrintColor;
        break;
        
      case CONS_LEVEL_MID:
        result = foreground ? &bc->midfgPrintColor : &bc->midbgPrintColor;
        break;
        
      case CONS_LEVEL_LOW:
        result = foreground ? &bc->lowfgPrintColor : &bc->lowbgPrintColor;
        break;

      default:
        g_critical("Program error: invalid conservation level '%d' when getting conservation colour.\n", consLevel);
        break;
    };

  return result;
}


/* Utility to get foreground/background the colour number for the given
 * conservation level. Returns a pointer to the value in the context. */
int* getConsColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground)
{
  if (bc->printColorsOn)
    return getConsPrintColor(bc, consLevel, foreground);

  int *result = &bc->lowbgColor;

  switch (consLevel)
    {
      case CONS_LEVEL_MAX:
        result = foreground ? &bc->maxfgColor : &bc->maxbgColor;
        break;
        
      case CONS_LEVEL_MID:
        result = foreground ? &bc->midfgColor : &bc->midbgColor;
        break;
        
      case CONS_LEVEL_LOW:
        result = foreground ? &bc->lowfgColor : &bc->lowbgColor;
        break;

      default:
        g_critical("Program error: invalid conservation level '%d' when getting conservation colour.\n", consLevel);
        break;
    };

  return result;
}


/* Set whether the currently-selected sequence should be excluded
 * from the conservation colours calculations */
void setExcludeFromConsCalc(BelvuContext *bc, const gboolean exclude)
{
  if (!bc->selectedAln) 
    {
      g_critical("Please select a sequence first.\n");
      return;
  }
  
  if ((exclude && bc->selectedAln->nocolor) ||
     (!exclude && !bc->selectedAln->nocolor))
    {
      /* Nothing to do */
      return;
    }
  
  /* Store orignal state in the nocolor field:
   1 = normal sequence
   2 = markup line
   
   This is needed to restore markup lines (nocolor lines are always markups)
   */
  
  if (exclude) 
    {
      if (bc->selectedAln->markup)
        bc->selectedAln->nocolor = 2;
      else
        bc->selectedAln->nocolor = 1;
      
      bc->selectedAln->markup = 1;
    }
  else 
    {
      if (bc->selectedAln->nocolor == 1) 
        bc->selectedAln->markup = 0;
      
      bc->selectedAln->nocolor = 0;
    }
}

/***********************************************************
 *		          Arrays			   *
 ***********************************************************/

/* Finds Entry s from an array of ALNs  a
 * sorted in ascending order of order()
 * If found, returns TRUE and sets *ip
 * if not, returns FALSE and sets *ip one step left
 */
gboolean alnArrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer))
{
  int ord;
  int i = 0 , j = a->len, k;

  if (!j || (ord = orderFunc(&s, &g_array_index(a, ALN*, 0))) < 0)
    { 
      if (ip)
	*ip = -1; 
      return FALSE;
    }   /* not found */

  if (ord == 0)
    { 
      if (ip)
	*ip = 0;
      return TRUE;
    }

  if ((ord = orderFunc(&s, &g_array_index(a, ALN*, --j))) > 0)
    {
      if (ip)
	*ip = j; 
      return FALSE;
    }
  
  if (ord == 0)
    { 
      if (ip)
	*ip = j;
      return TRUE;
    }

  while (TRUE)
    { 
      k = i + ((j-i) >> 1) ; /* midpoint */

      if ((ord = orderFunc(&s, &g_array_index(a, ALN*, k))) == 0)
	{ 
          if (ip)
	    *ip = k; 
	  return TRUE;
	}
      
      if (ord > 0) 
	(i = k);
      else
	(j = k) ;
      
      if (i == (j-1) )
        break;
    }
  
  if (ip)
    *ip = i ;
  
  return FALSE;
}


/* As alnArrayFind, but for an array of BootstrapGroup*'s */
gboolean bsArrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer))
{
  int ord;
  int i = 0 , j = a->len, k;
  
  if (!j || (ord = orderFunc(&s, &g_array_index(a, BootstrapGroup*, 0))) < 0)
    { 
      if (ip)
	*ip = -1; 
      return FALSE;
    }   /* not found */
  
  if (ord == 0)
    { 
      if (ip)
	*ip = 0;
      return TRUE;
    }
  
  if ((ord = orderFunc(&s, &g_array_index(a, BootstrapGroup*, --j))) > 0)
    {
      if (ip)
	*ip = j; 
      return FALSE;
    }
  
  if (ord == 0)
    { 
      if (ip)
	*ip = j;
      return TRUE;
    }
  
  while (TRUE)
    { 
      k = i + ((j-i) >> 1) ; /* midpoint */
      
      if ((ord = orderFunc(&s, &g_array_index(a, BootstrapGroup*, k))) == 0)
	{ 
          if (ip)
	    *ip = k; 
	  return TRUE;
	}
      
      if (ord > 0) 
	(i = k);
      else
	(j = k) ;
      
      if (i == (j-1) )
        break;
    }
  
  if (ip)
    *ip = i ;
  
  return FALSE;
}


/* Linear search for an exact matching sequence name and coordinates,
 typically to find back highlighted row after sorting 
 */
gboolean alignFind(GArray *alignArr, ALN *obj, int *idx)
{
  for (*idx = 0; *idx < (int)alignArr->len; (*idx)++) 
    {
      ALN *currentAln = g_array_index(alignArr, ALN*, *idx);
      
      if (alphaorder(&currentAln, &obj) == 0)
        return TRUE;
    }

  return FALSE;
}


/* Set the ->nr field in Align array in ascending order */
void arrayOrder(GArray *alignArr)
{
  int i = 0;
  
  for (i = 0; i < (int)alignArr->len; ++i) 
    {
      ALN *alnp = g_array_index(alignArr, ALN*, i);
      alnp->nr = i + 1;
    }
}


///* Set the ->nr field in Align array in ascending order with increments of 10
//   Necessary if one wants to insert items.
// */
//static void arrayOrder10(GArray *alignArr)
//{
//  int i;
//
//  for (i = 0; i < alignArr->len; ++i) 
//    g_array_index(alignArr, ALN*, i)->nr = (i+1)*10;
//}



/* Separate markuplines to another array before resorting
 */
void separateMarkupLines(BelvuContext *bc)
{
  bc->markupAlignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN*), 100);
  
  arrayOrder(bc->alignArr);
  
  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; ) 
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (alnp->markup) 
        {
          /* printf ("Moving line %d, %s/%d-%d, nseq=%d\n", i, alnp->name, alnp->start, alnp->end, nseq);*/
          g_array_append_val(bc->markupAlignArr, alnp);
          g_array_sort(bc->markupAlignArr, alphaorder);
          g_array_remove_index(bc->alignArr, i);
        }
      else
        {
          ++i;
        }
    }
  
  arrayOrder(bc->alignArr);
}


/* Reinsert markuplines after mother sequence or if orphan at bottom */
void reInsertMarkupLines(BelvuContext *bc)
{
  int i, j;
  char tmpname[MAXNAMESIZE+1], *cp;
  
  g_array_sort(bc->alignArr, nrorder); /* to do: can we move this out of the loop ? */

  for (i = bc->markupAlignArr->len - 1; i >=0 ; --i)
    {
      ALN *alnp = g_array_index(bc->markupAlignArr, ALN*, i);

      /* Insert the markup line after the non-markup entry with the same name
       * minus the postfix (after a space); or at the end of the array if none
       * exists. */
      strcpy(tmpname, alnp->name);
      
      if ((cp = strchr(tmpname, ' ')))
        *cp = 0;
      
      for (j = 0; j < (int)bc->alignArr->len; ++j)
        {
          if (!strcmp(tmpname, g_array_index(bc->alignArr, ALN*, j)->name))
            break;
        }

      if (j >= (int)bc->alignArr->len)
	g_array_append_val(bc->alignArr, alnp);
      else
	g_array_insert_val(bc->alignArr, j + 1, alnp);
    }

  arrayOrder(bc->alignArr);
}


int strcmp_(gconstpointer xIn, gconstpointer yIn)
{
  const char *x = (const char*)xIn;
  const char *y = (const char*)yIn;
  
  int retval = strcmp(x, y);
  return retval;
}


GArray *copyAlignArray(GArray *inputArr)
{
  GArray *result = g_array_sized_new(FALSE, FALSE, sizeof(char*), inputArr->len);
  
  int i = 0;
  for ( ; i < (int)inputArr->len; ++i)
    {
      ALN *inputAln = g_array_index(inputArr, ALN*, i);
      ALN *destAln = createEmptyAln();
      alncpy(destAln, inputAln); /* shallow copy; doesn't copy sequence string */
      
      g_array_append_val(result, destAln);

      /* Duplicate the sequence string, if not null */
      if (alnGetSeq(inputAln))
        destAln->sequenceStr = g_string_new(alnGetSeq(inputAln));
      else
        destAln->sequenceStr = NULL;
    }
  
  return result;
}


void columnCopy(GArray *alignArrDest, int destIdx, GArray *alignArrSrc, int srcIdx)
{
  int i;
  
  for (i = 0; i < (int)alignArrSrc->len; ++i)
    {
      ALN *srcAln = g_array_index(alignArrSrc, ALN*, i);
      ALN *destAln = g_array_index(alignArrDest, ALN*, i);

      char *srcSeq = alnGetSeq(srcAln);
      char *destSeq = alnGetSeq(destAln);
      
      if (srcSeq && destSeq && destIdx < alnGetSeqLen(destAln) && srcIdx < alnGetSeqLen(srcAln))
        destSeq[destIdx] = srcSeq[srcIdx];
    }
}



/***********************************************************
 *		          Context			   *
 ***********************************************************/

/* Create the colors that belvu will use for various specific purposes */
static void createBelvuColors(BelvuContext *bc)
{
  /* Initialise the array with empty BlxColor structs */
  bc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), BELCOLOR_NUM_COLORS);
  int i = BELCOLOR_MIN + 1;
  
  for ( ; i < BELCOLOR_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = (BlxColor*)g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(bc->defaultColors, *blxColor);
    }
  
  createBlxColor(bc->defaultColors, BELCOLOR_BACKGROUND, "Background", "Background color", BLX_WHITE, BLX_WHITE, "#bdbdbd", NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_ALIGN_TEXT, "Text color for alignments", "Text color for alignments", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_COLUMN_HIGHLIGHT, "Highlight color for selected column", "Highlight color for selected column",  "#dddddd", BLX_BLACK, NULL, NULL);
  
  /* Trees */
  createBlxColor(bc->defaultColors, BELCOLOR_TREE_BACKGROUND, "Tree background", "Tree background color", BLX_WHITE, BLX_WHITE, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_TREE_LINE, "Default tree line color", "Default tree line color", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_TREE_TEXT, "Default tree text color", "Default tree text color", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_TREE_BOOTSTRAP, "Tree boostrap line color", "Tree boostrap line color", BLX_BLUE, BLX_BLUE, NULL, NULL);
  
  /* Conservation plot */
  createBlxColor(bc->defaultColors, BELCOLOR_CONS_PLOT, "Line color of the conservation plot", "Line color of the conservation plot", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_CONS_PLOT_AVG, "Average-conservation line color on the conservation profile", "Average-conservation line color on the conservation profile", BLX_RED, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BELCOLOR_CONS_PLOT_SCALE, "Scale color for the conservation profile", "Scale color for the conservation profile", BLX_DARK_GREY, BLX_DARK_GREY, NULL, NULL);
}


/* Create the context, which contains all program-wide variables */
BelvuContext* createBelvuContext()
{
  BelvuContext *bc = (BelvuContext*)g_malloc(sizeof *bc);
  
  bc->belvuWindow = NULL;
  bc->spawnedWindows = NULL;
  bc->belvuTree = NULL;
  bc->belvuAlignment = NULL;
  bc->consPlot = NULL;
  bc->orgsWindow = NULL;
  
  bc->defaultCursor = NULL; /* get from gdkwindow once it is shown */
  bc->removeSeqsCursor = gdk_cursor_new(GDK_PIRATE);
  bc->busyCursor = gdk_cursor_new(GDK_WATCH);

  bc->defaultColors = NULL;
  createBelvuColors(bc);
  
  bc->alignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN*), 100);
  bc->organismArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN*), 100);
  bc->markupAlignArr = NULL;
  bc->bootstrapGroups = NULL;
  
  bc->selectedAln = NULL;
  bc->highlightedAlns = NULL;
  
  bc->mainTree = NULL;
  bc->treeBestBalancedNode = NULL;
  
  bc->treeReadDistancesPipe = NULL;
  
  bc->IN_FORMAT = MUL;
  bc->maxScoreLen = 0;
  bc->alignYStart = 0;
  bc->treebootstraps = 0; 
  bc->maxLen = 0;
  bc->maxTreeWidth = 0;
  bc->maxNameLen = 0;   
  bc->maxStartLen = 0; 
  bc->maxEndLen = 0; 
  bc->maxScoreLen = 0; 
  bc->selectedCol = 0;
  bc->highlightedCol = 0;
  
  bc->maxfgColor = DEFAULT_MAX_FG_COLOR;
  bc->midfgColor = DEFAULT_MID_FG_COLOR,
  bc->lowfgColor = DEFAULT_LOW_FG_COLOR;
  bc->maxbgColor = DEFAULT_MAX_BG_COLOR;
  bc->midbgColor = DEFAULT_MID_BG_COLOR;
  bc->lowbgColor = DEFAULT_LOW_BG_COLOR;
  bc->maxfgPrintColor = DEFAULT_MAX_FG_PRINT_COLOR;
  bc->midfgPrintColor = DEFAULT_MID_FG_PRINT_COLOR,
  bc->lowfgPrintColor = DEFAULT_LOW_FG_PRINT_COLOR;
  bc->maxbgPrintColor = DEFAULT_MAX_BG_PRINT_COLOR;
  bc->midbgPrintColor = DEFAULT_MID_BG_PRINT_COLOR;
  bc->lowbgPrintColor = DEFAULT_LOW_BG_PRINT_COLOR;
  
  bc->schemeType = BELVU_SCHEME_TYPE_CONS;
  bc->residueScheme = BELVU_SCHEME_ERIK;
  bc->consScheme = BELVU_SCHEME_BLOSUM;
  
  bc->treeMethod = NJ;
  bc->treeDistCorr = SCOREDIST;
  bc->treePickMode = NODESWAP;
  
  bc->sortType = BELVU_UNSORTED;
  
  bc->annotationList = NULL;
  
  bc->treeBestBalance = 0.0;
  bc->treeBestBalance_subtrees = 0.0;
  bc->tree_y = 0.3;
  bc->lowIdCutoff = DEFAULT_LOW_ID_CUTOFF;
  bc->midIdCutoff = DEFAULT_MID_ID_CUTOFF;
  bc->maxIdCutoff = DEFAULT_MAX_ID_CUTOFF;
  bc->lowSimCutoff = DEFAULT_LOW_SIM_CUTOFF;
  bc->midSimCutoff = DEFAULT_MID_SIM_CUTOFF;
  bc->maxSimCutoff = DEFAULT_MAX_SIM_CUTOFF;
  bc->colorByResIdCutoff = 20.0;
  bc->mksubfamilies_cutoff = 0.0;
  bc->treeScale = DEFAULT_TREE_SCALE_CORR;
  bc->treeLineWidth = 0.3;
  
  bc->gapChar = '.';
  bc->saveSeparator = '/';
  strcpy(bc->treeDistString, SCOREDISTstr);
  strcpy(bc->treeMethodString, NJstr);
  bc->Title[0] = '\0';
  bc->saveFormat = BELVU_FILE_MUL;
  bc->fileName = NULL;
  bc->dirName = NULL;
  bc->organismLabel[0] = 'O';
  bc->organismLabel[0] = 'S';   
  bc->organismLabel[0] = '\0'; 
  
  bc->conservCount = NULL;
  bc->colorMap = NULL;
  bc->conservResidues = NULL;
  bc->conservation = NULL;
  
  bc->treeCoordsOn = TRUE;
  bc->treeReadDistancesOn = FALSE;
  bc->treePrintDistances = FALSE;
  bc->penalize_gaps = FALSE;
  bc->stripCoordTokensOn = TRUE;
  bc->saveCoordsOn = TRUE;
  bc->displayScores = FALSE;
  bc->outputBootstrapTrees = FALSE;
  bc->treebootstrapsDisplay = FALSE;
  bc->treeColorsOn = TRUE;
  bc->treeShowOrganism = TRUE;
  bc->treeShowBranchlen = FALSE;
  bc->matchFooter = FALSE;
  bc->saved = TRUE;
  bc->ignoreGapsOn = FALSE;
  bc->colorByResIdOn = FALSE;
  bc->rmEmptyColumnsOn = TRUE;
  bc->lowercaseOn = FALSE;
  bc->removingSeqs = FALSE;
  bc->displayColors = TRUE;
  bc->haveCustomColors = FALSE;
  bc->printColorsOn = FALSE;
  bc->highlightOrthologs = FALSE;
  bc->useWWWFetch = FALSE;
  bc->initTree = FALSE;
  bc->onlyTree = FALSE;
  
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < BELDIALOG_NUM_DIALOGS; ++dialogId)
    {
      bc->dialogList[dialogId] = NULL;
    }
  
  return bc;
}


/* Destroy the context */
void destroyBelvuContext(BelvuContext **bc)
{
  if (bc && *bc)
    {
    if ((*bc)->alignArr)
      g_array_unref((*bc)->alignArr);
    
    if ((*bc)->organismArr)
      g_array_unref((*bc)->organismArr);
    
    if ((*bc)->markupAlignArr)
      g_array_unref((*bc)->markupAlignArr);
    
    if ((*bc)->bootstrapGroups)
      g_array_unref((*bc)->bootstrapGroups);
    
    g_free(*bc);
    *bc = NULL;
    }
}


/***********************************************************
 *                           Utilities                     *
 ***********************************************************/

/* Returns a string which is the name of the Blixem application. */
const char *belvuGetAppName(void)
{
  return BELVU_TITLE ;
}

/* Returns a string which is the prefix to use for window titles. */
const char *belvuGetTitlePrefix(BelvuContext *bc)
{
  return bc->abbrevTitle ? BELVU_PREFIX_ABBREV : BELVU_PREFIX ;
}

/* Returns a copyright string for the Blixem application. */
const char *belvuGetCopyrightString(void)
{
  return BELVU_COPYRIGHT_STRING ;
}

/* Returns the Blixem website URL. */
const char *belvuGetWebSiteString(void)
{
  return BELVU_WEBSITE_STRING ;
}

/* Returns a comments string for the Blixem application. */
const char *belvuGetCommentsString(void)
{
  return BELVU_COMMENTS_STRING() ;
}

/* Returns a license string for the belvu application. */
const char *belvuGetLicenseString(void)
{
  return BELVU_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Blixem code. */
const char *belvuGetVersionString(void)
{
  return BELVU_VERSION_STRING ;
}

/* Utility to return the sequence data in the given alignment; returns null if
 * not set. */
char* alnGetSeq(ALN *aln)
{
  return (aln && aln->sequenceStr ? aln->sequenceStr->str : NULL);
}


/* Utility to return the length sequence data in the given alignment; returns
 * 0 if not set. */
int alnGetSeqLen(ALN *aln)
{
  return (aln && aln->sequenceStr ? aln->sequenceStr->len : 0);
}


/* This function just returns the value in the b2a array at the given index */
char b2aIndex(const int idx)
{
  return b2a[idx];
}


void drawText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const char *text, int *textWidth, int *textHeight)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, text);
  
  if (drawable)
    gdk_draw_layout(drawable, gc, x, y, layout);
  
  /* Return the width and height of the layout, if requested */
  pango_layout_get_size(layout, textWidth, textHeight);
  
  if (textWidth)
    *textWidth /= PANGO_SCALE;
  
  if (textHeight)
    *textHeight /= PANGO_SCALE;
  
  g_object_unref(layout);
}


/* Utility to draw the given integer as text. The text is right-aligned, so 
 * the input x coord must be the RIGHTMOST EDGE of where you want the text. */
void drawIntAsText(GtkWidget *widget, 
                   GdkDrawable *drawable, 
                   GdkGC *gc, 
                   const int x, 
                   const int y, 
                   const int value)
{
  char *tmpStr = g_strdup_printf("%d", value);
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, tmpStr);
  g_free(tmpStr);

  int textWidth;
  pango_layout_get_pixel_size(layout, &textWidth, NULL);

  gdk_draw_layout(drawable, gc, x - textWidth, y, layout);
  g_object_unref(layout);
}


/* Utility to draw the given double as text. The text is right-aligned, so 
 * the input x coord must be the RIGHTMOST EDGE of where you want the text. */
void drawDoubleAsText(GtkWidget *widget, 
                      GdkDrawable *drawable, 
                      GdkGC *gc, 
                      const int x, 
                      const int y, 
                      const double value)
{
  char *tmpStr = g_strdup_printf("%.1f", value);
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, tmpStr);
  g_free(tmpStr);
  
  int textWidth;
  pango_layout_get_pixel_size(layout, &textWidth, NULL);
  
  gdk_draw_layout(drawable, gc, x - textWidth, y, layout);
  g_object_unref(layout);
}


/* Calculate percent identity of two strings */
double identity(char *s1, char *s2, const gboolean penalize_gaps)
{
    int n, id;

    for (n = id = 0; *s1 && *s2; s1++, s2++) 
      {
	if (isGap(*s1) && isGap(*s2)) 
          continue;
        
	if (isGap(*s1) || isGap(*s2))
          {
	    if (!penalize_gaps) 
              continue;
          }
        
	n++;
	if (toupper(*s1) == toupper(*s2)) 
          id++;
      }

    if (n)
      return (double)id/n*100;
    else
      return 0.0;
}


/* Calculate score of two sequences */
static double score(char *s1, char *s2, const gboolean penalize_gaps)
{
    double sc=0.0;

    for (;*s1 && *s2; s1++, s2++) 
      {
	if (isGap(*s1) && isGap(*s2)) 
          {
	    continue;
          }
	else if (isGap(*s1) || isGap(*s2)) 
          {
	    if (penalize_gaps) sc -= 0.6;
          }
	else
          {
            int val1 = a2b[(unsigned char)(*s1)];
            int val2 = a2b[(unsigned char)(*s2)];
            
            if (val1 > 0 && val2 > 0)
              sc += (double) BLOSUM62[val1 - 1][val2 - 1];
          }
      }
    
    return sc;
}


gboolean isGap(char c) 
{
  if (c == '.' || c == '-' ||
      c == '[' || c == ']' /* Collapse-control chars */ ) 
    return TRUE;
  else
    return FALSE;
}


static gboolean isAlign(char c)
{
  if (isalpha(c) || isGap(c) || c == '*')
    return TRUE;
  else 
    return FALSE;
}


int GCGchecksum(BelvuContext *bc, char *seq)
{
  int  
  check = 0, 
  count = 0, 
  i=0;
  
  for (i = 0; i < bc->maxLen; ++i) 
    {
    ++count;
    check += count * toupper((int) seq[i]);
    
    if (count == 57) 
      count = 0;
    }
  
  return (check % 10000);
}


int GCGgrandchecksum(BelvuContext *bc)
{
  int 
  i=0,
  grand_checksum=0;
  
  for(i=0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      grand_checksum += GCGchecksum(bc, alnGetSeq(alnp));
    }
  
  return (grand_checksum % 10000);
}


/* Read labels of highlighted sequence and spread them */
void readLabels(BelvuContext *bc, FILE *fil)
{
  char *labelseq = (char*)g_malloc(bc->maxLen + 1); /* The raw sequence of labels */
  char *label = (char*)g_malloc(bc->maxLen + 1);    /* The mapped sequence of labels, 1-maxlen */
  char line[MAXLENGTH+1];
  
  /* read file */
  char *cq = labelseq;
  while (!feof (fil)) 
    {
      if (!fgets (line, MAXLENGTH, fil)) 
        break;
      
      char *cp = line;
      while (*cp) 
        {
          if (isalpha(*cp)) 
            *cq++ = *cp;
          cp++;
        }
    }

  fclose(fil);
  
  /* Warn if seq too long, return if too short */
  int seqlen = bc->selectedAln->end - bc->selectedAln->start + 1;
  if ((int)strlen(labelseq) > seqlen)
    {
      g_critical("The sequence of labels is longer (%d) than the sequence (%d).\nHope that's ok\n",
		 (int)strlen(labelseq), seqlen);
    }
  else if ((int)strlen(labelseq) < seqlen) 
    {
      g_critical("The sequence of labels is shorter (%d) than the sequence (%d).\nAborting\n",
		 (int)strlen(labelseq), seqlen);
      return;
    }
  
  /* map labels to alignment */
  int col = 0;
  int seqpos = 0;
  
  for (col = 0, seqpos = 0; col < bc->maxLen; col++) 
    {
      label[col] = labelseq[seqpos];
      const char *selectedSeq = alnGetSeq(bc->selectedAln);
      
      if (selectedSeq && isalpha(selectedSeq[col]) && labelseq[seqpos+1]) 
        seqpos++;
    }
  
  int row = 0;
  for (row = 0; row < (int)bc->alignArr->len; ++row) 
    {
      ALN *alnrow = g_array_index(bc->alignArr, ALN*, row);
      
      for (col = 0; col < bc->maxLen; col++)
        {
          char *rowSeq = alnGetSeq(alnrow);
          
          if (rowSeq && isalpha(rowSeq[col]))
            rowSeq[col] = label[col];
        }
    }
  
  g_free(labelseq);
}

/***********************************************************
 *		          			   *
 ***********************************************************/

static void makeSegList(BelvuContext *bc, SEG **SegList, char *line)
{
  int n, i;
  char *p, *linecopy;
  SEG *seg, *prevseg=0;

  /* Count coordinates - has to be multiple of 4 */
  linecopy = (char*)g_malloc(strlen(line)+1);
  strcpy(linecopy, line);
  n = 0;
  
  if (atoi(strtok(linecopy, " "))) 
    n++;;
  
  while ( (p = strtok(0, " ")) && atoi(p) ) 
    n++;
  
  g_free(linecopy);
    
  if (!n || n % 4) 
    g_error("Segments not multiple of 4 ints (%s)\n", line);

  for (i = 0; i < n/4; i++) 
    {
      seg = (SEG *)g_malloc(sizeof(SEG));
      if (prevseg) 
        prevseg->next = seg;
      else
        *SegList = seg;
      
      prevseg = seg;

      seg->qstart = atoi(i ? strtok(0, " ") : strtok(line, " "));
      seg->qend = atoi(strtok(0, " "));
      seg->start = atoi(strtok(0, " "));
      seg->end = atoi(strtok(0, " "));
      seg->next = NULL;
    }

  for (seg = *SegList; seg; seg = seg->next) 
    {
      DEBUG_OUT("%d %d %d %d\n", seg->qstart, seg->qend, seg->start, seg->end);

      if (seg == *SegList && seg->qstart != 1)
        g_error("Bad qstart: Must start on 1\n");

      if (seg->qstart < 1 || seg->qstart > bc->maxLen)
        g_error("Bad qstart: %d.  Range: 1-%d\n", seg->qstart, bc->maxLen);
      if (seg->qend < 1 || seg->qend > bc->maxLen)
        g_error("Bad qend: %d.  Range: 1-%d\n", seg->qend, bc->maxLen);
      if (seg->start < 1 || seg->start > bc->maxLen)
        g_error("Bad start: %d.  Range: 1-%d\n", seg->start, bc->maxLen);
      if (seg->end < 1 || seg->end > bc->maxLen)
        g_error("Bad end: %d.  Range: 1-%d\n", seg->end, bc->maxLen);

      if (seg->qstart > seg->qend)
        g_error("qstart > qend  (%d > %d)\n", seg->qstart, seg->qend);
      if (seg->start > seg->end)
        g_error("start > end  (%d > %d)\n", seg->start, seg->end);
    }
}


static int countInserts(SEG *seg)
{
  int   
    Align_gap, Query_gap, gap, insert_counter=0;

  if (!seg)
    return insert_counter ;

  while (seg->next) 
    {
      Align_gap = seg->next->start - seg->end - 1;
      if (Align_gap < 0) 
        g_error("Negative Align_gap: %d (%d-%d)\n", Align_gap, seg->start, seg->next->end );
		
      Query_gap = seg->next->qstart - seg->qend - 1;
      if (Query_gap < 0) 
        g_error("Negative Query_gap: %d (%d-%d)\n", Query_gap, seg->qstart, seg->next->qend );

      gap = Query_gap - Align_gap;
      if (gap > 0) 
        {
          insert_counter += gap;
	}
      
      seg = seg->next;
    }
  
  return insert_counter;
}


/* Inserts n gap columns after position p, in sequence coordinate, 1...maxLen 
 */
static void insertColumns(BelvuContext *bc, int p, int n)
{
  int
    i, j;

  ALN 
    *alni;
  char
    *dest, *src, *seq;

  g_message("Inserting %d columns after column %d\n", n, p);

  bc->maxLen += n;

  for (i = 0; i < (int)bc->alignArr->len; ++i) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
	
      seq = (char*)g_malloc(bc->maxLen + 1);

      dest = seq;
      src = alnGetSeq(alni);

      for (j = 0; j < p;  j++) 
        {
          *dest++ = *src++;
	}
      for (; j < p+n;  j++) 
        {
          *dest++ = '.';
	}
      for (; j < bc->maxLen;  j++) 
        {
          *dest++ = *src++;
	}

      g_string_free(alni->sequenceStr, TRUE);
      alni->sequenceStr = g_string_new(seq);
    }

  bc->saved = FALSE;
}


/* readMatch displays a matching sequences to the alignment so that
   the match can be displayed in relation to the whole family.

   Considerations:

   - Must work both from ACEDB and command line (and from Web server).

   Display:

   - The score should be visible (turn on score column.)

   - Should matches be displayed on separate lines or in the middle of the alignment?
     On top would be nice, but would require a lot of extra programming (status bar etc).
     Instead add to alignment (SEED usually fits on one screen) 
     Draw names with red background.  
     Add on top.


   Steps needed on command line:

   hmmb -W tmp HMM align
   selex2alignmap tmp > alignmap

   hmm?s -F HMM query | ?smapback map=alignmap > query.mapmatch  ( -> acedb)
   belvu -m query.mapmatch align

   The two last lines should be scripted in hmm?sBelvu


*/
void readMatch(BelvuContext *bc, FILE *fil)
{
  /* Format:
       
  Using segments is more difficult than residues, but is
  necessary to display in ACEDB Pepmap.

  seg_match/3-30                   (matching segment name)
  WLPLHTLinsertAACGEFYLVDSLLKH     (only matching part, no pads!)
  1 7 3 9  14 28 10 24

  seq_match2...
  */
  int	orig_maxLen, tmp,
    Align_gap, Query_gap, gap, inserts, seglen,
    done_one=0, len,  was_saved=bc->saved, insert_counter;
  char *cp, *seq, *rawseq = NULL, *seqp;
  SEG	*SegList = NULL, *seg;
  char line[MAXLENGTH+1];
  
  while (!done_one)
    {
      if (!fgets(line, MAXLENGTH, fil))
	break;

      /*printf("%s\n", line); continue;*/

      if (*line != '\n' && *line != '#')
	{
	  if (done_one)
	    {
	      g_critical("Sorry, can only do one match\n");
	      continue;
	    }

          ALN *aln = createEmptyAln();

	  /* Name */
	  if (!strtok(line, "/"))
	    g_error("Bad format: %s\n", line);

	  strncpy(aln->name, line, MAXNAMESIZE);
	  aln->name[MAXNAMESIZE] = 0;
          
	  if ((cp = strchr(aln->name, ' ')))
	    *cp = 0;

	  if (bc->maxNameLen  < (int)strlen(aln->name))
	    bc->maxNameLen = strlen(aln->name);

	  /* Start & End */
	  if (!(cp = strtok(0, "-"))) 
            g_error("Bad start: %s\n", cp);

	  aln->start = atoi(cp);
          
	  if ((len=strlen(cp)) > bc->maxStartLen)
	    bc->maxStartLen = len;
          
          char *tmpStr = g_strdup_printf("%d", aln->start);

	  if (bc->maxStartLen < (tmp = strlen(tmpStr)))
	    bc->maxStartLen = tmp;
                                   
          g_free(tmpStr);
          tmpStr = NULL;

	  if (!(cp = strtok(0, " ")))
	    g_error("Bad end: %s\n", cp);
          
	  aln->end = atoi(cp);
          
	  if ((len=strlen(cp)) > bc->maxEndLen)
	    bc->maxEndLen = len;
          
          tmpStr = g_strdup_printf("%d", aln->end);
          
	  if (bc->maxEndLen < (tmp = strlen(tmpStr)))
	    bc->maxEndLen = tmp;
    
          g_free(tmpStr);
          tmpStr = NULL;
          
          int ip = 0;
	  if (alignFind(bc->alignArr, aln, &ip))
	    {
	      g_critical("Sorry, already have sequence %s/%d-%d\n",
                         aln->name, aln->start, aln->end);
	      return;
	    }

	  /* Score */
	  if (!(cp = strtok(0, "\n")))
	    g_error("Bad score: %s\n", cp);

	  aln->score = atof(cp);
          
          tmpStr = g_strdup_printf("%.1f", aln->score);
          
	  if ((len=strlen(tmpStr)) > bc->maxScoreLen)
	    bc->maxScoreLen = len;
          
          g_free(tmpStr);
          tmpStr = NULL;

	  bc->displayScores = TRUE;

	  /* Sequence */
	  if (!fgets(line, MAXLENGTH, fil))
	    break;

	  if ((cp = strchr(line, '\n')))
	    *cp = 0 ;
          
	  rawseq = (char*)g_malloc(strlen(line)+1);
	  strcpy(rawseq, line);

	  /* Matching segments */
          
	  if (!fgets(line, MAXLENGTH, fil))
	    break;

	  if ((cp = strchr(line, '\n')))
            *cp = 0;

	  makeSegList(bc, &SegList, line);

	  inserts = countInserts(SegList);
	  seq = (char*)g_malloc(bc->maxLen + inserts + 1);
	  memset(seq, '.', bc->maxLen + inserts);
	    
	  orig_maxLen = bc->maxLen;

	  /* Add first segment */
	  seg = SegList;

	  seqp = seq + seg->start-1;
	  seglen = seg->end - seg->start+1;
	  strncpy(seqp, rawseq, seglen);
	  seqp += seglen; 

	  /* Add subsequent segments */
	  insert_counter = 0;
	  while (seg->next)
	    {
	      Align_gap = seg->next->start - seg->end - 1;
	      Query_gap = seg->next->qstart - seg->qend - 1;


	      if (Query_gap > 0)
		{ 
		  strncpy(seqp, rawseq + seg->qend, Query_gap);
		  seqp += Query_gap; 
		}

	      gap = Query_gap - Align_gap;
	      if (gap > 0)
		{
		  insertColumns(bc, seg->end+insert_counter, gap);
		  insert_counter += gap;
		}
	      else if (gap < 0)
		{
		  seqp += -gap;
		  DEBUG_OUT("inserting %d gaps at %d\n", -gap, (int)(seqp-seq));
		}

	      seg = seg->next;

	      /* Add sequence segment */
	      seglen = seg->end - seg->start+1;
	      strncpy(seqp, rawseq + seg->qstart-1, seglen);
	      seqp += seglen;
	    }

	  /* Add final pads */
	  memset(seqp, '.', orig_maxLen - seg->end);

	  aln->sequenceStr = g_string_new(seq);
          g_free(seq);
          
	  aln->color = RED;
	  aln->nr = 0;

          g_array_append_val(bc->alignArr, aln);
          
	  done_one = 1;

	} /* End of this matchSeq */
    }

  g_array_sort(bc->alignArr, nrorder);

  g_free(rawseq);

  arrayOrder(bc->alignArr);
  bc->saved = was_saved;

  return ;
}


void checkAlignment(BelvuContext *bc)
{
  int i, g, cres, nres, tmp;

  bc->maxNameLen = bc->maxStartLen = bc->maxEndLen = 0;

  for (i = 0; i < (int)bc->alignArr->len; ++i) 
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

      if (!alnp->markup) 
        {
          char *alnSeq = alnGetSeq(alnp);
          
          /* Count residues */
          for (cres = g = 0; g < bc->maxLen; g++) 
            {
              if (alnSeq && isAlign(alnSeq[g]) && !isGap(alnSeq[g])) 
                cres++;
	    }
	    
          if (!alnp->start) 
            {
              /* No coords provided - reconstruct them */
              alnp->start = 1;
              alnp->end = cres;
	    }
          else 
            {
              /* Check if provided coordinates are correct */
              nres = abs(alnp->end - alnp->start) + 1;
              
              if ( nres != cres)
                {
                  g_warning("Found wrong number of residues in %s/%d-%d: %d instead of %d\n",
                            alnp->name, alnp->start, alnp->end,
                            cres,
                            nres);
                }
	    }
	}

      /* Find max string length of name, start, and end */
      if (bc->maxNameLen  < (int)strlen(alnp->name))
        bc->maxNameLen = strlen(alnp->name);

      char *tmpStr = g_strdup_printf("%d", alnp->start);
      
      if (bc->maxStartLen < (tmp = strlen(tmpStr))) 
        bc->maxStartLen = tmp;

      g_free(tmpStr);
      tmpStr = NULL;

      tmpStr= g_strdup_printf("%d", alnp->end);

      if (bc->maxEndLen   < (tmp = strlen(tmpStr))) 
        bc->maxEndLen = tmp;
    }
}    


static void initConservMtx(BelvuContext *bc)
{
  int i;

  bc->conservCount = (int **)g_malloc(21*sizeof(int *));
  bc->colorMap = (int **)g_malloc(21*sizeof(int *));
    
  for (i = 0; i < 21; ++i) 
    {
      bc->conservCount[i] = (int *)g_malloc(bc->maxLen*sizeof(int));
      bc->colorMap[i] = (int *)g_malloc(bc->maxLen*sizeof(int));
    }

  bc->conservResidues = (int*)g_malloc(bc->maxLen*sizeof(int));
  bc->conservation = (double*)g_malloc(bc->maxLen*sizeof(double));
}	


/* This populates conservCount (the count of how many of each residue there is
 * in each column) and conservResidues (the count of how many residues in total
 * there are in each column). 
 * Returns: the total number of sequences (excluding markup lines) */
static int countResidueFreqs(BelvuContext *bc)
{
  int   i, j, totalNumSeqs=0;

  if (!bc->conservCount) 
    initConservMtx(bc);

    /* Must reset since this routine may be called many times */
  for (i = 0; i < bc->maxLen; ++i)
    {
      for (j = 0; j < 21; j++)
        {
          bc->conservCount[j][i] = 0;
          bc->colorMap[j][i] = 0;
        }
      
      bc->conservResidues[i] = 0;
    }
    
  for (j = 0; j < (int)bc->alignArr->len; ++j)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, j);
      
      if (alnp->markup)
        continue;

      totalNumSeqs++;

      for (i = 0; i < bc->maxLen; ++i)
        {
          char *alnSeq = alnGetSeq(alnp);
          int val = n2b[(unsigned char)(alnSeq[i])];

          bc->conservCount[val][i]++;

          if (isalpha(alnSeq[i]) || alnSeq[i] == '*')
            bc->conservResidues[i]++;
        }
    }

  return totalNumSeqs;
}


/* Calculate overhang between two aligned sequences 
   Return nr of residues at the ends unique to s1 (overhanging s2).
*/
static int alnOverhang(char *s1, char *s2)
{
  int overhang = 0,
    s1started = 0,
    s2started = 0;
  char *s1save = s1, *s2save = s2;

  for (; *s1 && *s2; s1++, s2++) 
    {
      if (!isGap(*s1)) s1started = 1;
      if (!isGap(*s2)) s2started = 1;
      if (s1started && !s2started) overhang++;
    }	
  
  s1--; s2--;
  s1started = s2started = 0;
  for (; s1>=s1save && s2>=s2save; s1--, s2--) 
    {
      if (!isGap(*s1)) s1started = 1;
      if (!isGap(*s2)) s2started = 1;
      if (s1started && !s2started) overhang++;
    }	

  /* printf ("%s\n%s\nOverhang=%d\n", s1save, s2save, overhang);*/

  return overhang;
}


/***********************************************************
 *		Removing alignments/columns		   *
 ***********************************************************/

void rmFinaliseColumnRemoval(BelvuContext *bc)
{
  rmGappySeqs(bc, 100.0);
  rmFinalise(bc);
}

/* Removes all columns between from and to, inclusively.
 * Note that from and to are sequence coordinates, 1...maxLen !!!!!  */
void rmColumn(BelvuContext *bc, const int from, const int to)
{
  int
    i=0, j=0, len = to - from + 1;
  ALN 
    *alni=NULL;
  
  g_message_info("Removing Columns %d-%d.\n", from, to);

  for (i = 0; i < (int)bc->alignArr->len; i++) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
      char *alnSeq = alnGetSeq(alni);

      /* If N or C terminal trim, change the coordinates */
      if (from == 1)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            
            if (!isGap(alnSeq[j-1]))
              (alni->start < alni->end ? alni->start++ : alni->start--); 
          }

      if (to == bc->maxLen)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            if (!isGap(alnSeq[j-1]))
              (alni->start < alni->end ? alni->end-- : alni->end++); 
          }
	
      /* Remove the columns */
      for (j = 0; j < bc->maxLen-to+1 /* Including terminator 0 */;  j++) 
        {
          alnSeq[from+j-1] = alnSeq[to+j];
	}
      j--;
      
      if (alnSeq[from+j-1] || alnSeq[to+j])
        g_warning("Still a bug in rmColumn(): End=%c, Oldend=%c\n", alnSeq[from+j-1], alnSeq[to+j]);
    }

  bc->maxLen -= len;

  bc->saved = FALSE;
}


/* Remove columns whose conservation is below the given cutoff */
void rmColumnCutoff(BelvuContext *bc, const double from, const double to)
{
  int 
    i, j, max, removed=0, oldmaxLen=bc->maxLen;
  static double cons ;
  
  
  for (i = bc->maxLen-1; i >= 0; i--) 
    {
      if (colorBySimilarity(bc)) 
        {
          cons = bc->conservation[i];
        }
      else 
        {
          max = 0;
          for (j = 1; j < 21; j++)
            {
              if (bc->conservCount[j][i] > max) 
                {
                  max = bc->conservCount[j][i];
                }
            }	
          cons = (double)max / bc->alignArr->len;
        }
      
      if (cons > from && cons <= to) 
        {
          g_message("removing %d, cons= %.2f\n", i+1, cons);
          rmColumn(bc, i+1, i+1);
          if (++removed == oldmaxLen) 
            {
              g_critical("You have removed all columns.  Prepare to exit Belvu\n");
              exit(EXIT_SUCCESS);
            }
        }
    }
  
  bc->saved = FALSE;
  rmFinaliseColumnRemoval(bc);
}


void rmEmptyColumns(BelvuContext *bc, double cutoff)
{
  int
    i=0, j=0, gaps=0, totseq=0, removed=0, oldmaxLen=bc->maxLen;
  ALN 
    *alni=NULL;
  char
    c='\0';

  for (totseq = i = 0; i < (int)bc->alignArr->len; ++i) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
      
      if (!alni->markup) 
        totseq++;
    }

  for (j = bc->maxLen-1; j >= 0; j--) 
    {
      for (gaps = i = 0; i < (int)bc->alignArr->len; i++) 
        {
          alni = g_array_index(bc->alignArr, ALN*, i);
          char *alnSeq = alnGetSeq(alni);
          
          if (!alni->markup) 
            {
              c = alnSeq[j];
              
              if (isGap(c) || c == ' ') 
                gaps++;
	    }
	}
      
      if ((double)gaps/totseq >= cutoff - MACHINE_RES) 
        {
          rmColumn(bc, j+1, j+1);
          if (++removed == oldmaxLen)
            {
              g_critical("You have removed all columns.  Prepare to exit Belvu\n");
              exit(0);
            }
	}
    }
}


/* Certain actions need to be enabled/disabled depending on certain properties */
void greyOutInvalidActionsForGroup(BelvuContext *bc, GtkActionGroup *action_group)
{
  if (!action_group)
    return;
  
  enableMenuAction(action_group, "Output", bc->displayScores);
  enableMenuAction(action_group, "rmScore", bc->displayScores);
  enableMenuAction(action_group, "scoreSort", bc->displayScores);
  
  enableMenuAction(action_group, "toggleColorByResId", !colorByConservation(bc));
  enableMenuAction(action_group, "colorByResId", !colorByConservation(bc));
  enableMenuAction(action_group, "saveColorScheme", !colorByConservation(bc));
  enableMenuAction(action_group, "loadColorScheme", !colorByConservation(bc));
  
  enableMenuAction(action_group, "colorSchemeCustom", bc->haveCustomColors);
  
  /* ignoreGaps used to be greyed out in old belvu when in color-by-residue
   * mode; it does affect the display when in that mode if using a %ID
   * threshold, though, so it should probably be enabled (or alternatively the 
   * logic should be adjusted so that it does not affect the residue colors, if 
   * that was the intent...) */
  /* enableMenuAction(action_group, "ignoreGaps", colorByConservation(bc)); */
  enableMenuAction(action_group, "printColors", colorByConservation(bc));
  
  enableMenuAction(action_group, "excludeHighlighted", bc->selectedAln != NULL);
}


/* Certain actions need to be enabled/disabled depending on certain properties */
void greyOutInvalidActions(BelvuContext *bc)
{
  /* Need to do the action groups for the main window and the tree (wrapped
   * windows currently don't have anything that gets greyed out but if they
   * do in future they will need updating here too) */
  
  GtkActionGroup *actionGroup = belvuWindowGetActionGroup(bc->belvuWindow);
  greyOutInvalidActionsForGroup(bc, actionGroup);

  actionGroup = belvuTreeGetActionGroup(bc->belvuTree);
  greyOutInvalidActionsForGroup(bc, actionGroup);
}


/* Set the main tree in the belvu context (null to reset). Takes ownership
 * of the given tree. Destroys the existing tree contents first, if there is one. */
void belvuContextSetTree(BelvuContext *bc, Tree **tree)
{
  if (bc->mainTree)
    {
      /* We don't destroy the actual tree struct because there may be pointers
       * to it from other places; instead, we just replace its contents. */
      destroyTreeContents(bc->mainTree);
      
      if (tree && *tree)
        {
          bc->mainTree->head = (*tree)->head;
          g_free(*tree);
          
          /* Update input pointer so that it points to the tree that actually
           * got set. */
          *tree = bc->mainTree;
        }
    }
  else if (tree && *tree)
    {
      bc->mainTree = *tree;
    }
  else
    {
      bc->mainTree = NULL;
    }
  
  /* Whether a tree exists affects whether some menu items are greyed out */
  greyOutInvalidActions(bc);
}


/* This function should be called after removing sequences or columns */
static void rmFinalise(BelvuContext *bc) 
{
  /*    ruler[maxLen] = 0;*/
  checkAlignment(bc);
  setConsSchemeColors(bc);
  
  /* Removing seqs/cols invalidates the tree, so set the tree head to NULL. */
  belvuContextSetTree(bc, NULL);
  
  /* Removing seqs/cols invalidates the conservation plot, so recalculate it */
  belvuConsPlotRecalcAll(bc->consPlot);
  
  /* Removing sequences can change the size of the columns in the alignment view,
   * so recalculate them */
  updateHeaderColumnsSize(bc->belvuAlignment);
  calculateDrawingSizes(bc->belvuAlignment);
}


/* Get rid of seqs that start or end with a gap.
 */
void rmPartialSeqs(BelvuContext *bc)
{
  int i=0, n=0;
  ALN *alni = NULL;

  for (i = 0; i < (int)bc->alignArr->len; ) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
      char *alnSeq = alnGetSeq(alni);
      
      if (isGap(alnSeq[0]) || isGap(alnSeq[bc->maxLen-1])) 
        {
          /* Remove entry */
          n++;

          if (bc->selectedAln == alni) 
            bc->selectedAln = 0;
          
          g_array_remove_index(bc->alignArr, i);
          bc->saved = 0;
	}
      else 
        {
          ++i;
        }
    }

  g_message_info("%d partial sequences removed.  %d seqs left.\n\n", n, bc->alignArr->len);

  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Get rid of seqs that are too gappy.
 */
void rmGappySeqs(BelvuContext *bc, const double cutoff)
{
  int i=0, j=0, n=0, gaps;
  ALN *alni = NULL;

  for (i = 0; i < (int)bc->alignArr->len; ) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
      char *alnSeq = alnGetSeq(alni);

      for (gaps = j = 0; j < bc->maxLen; j++)
        if (isGap(alnSeq[j])) 
          gaps++;

      if ((double)gaps/bc->maxLen >= cutoff/100.0) 
        {
          /* Remove entry */
          n++;
          
          if (bc->selectedAln == alni) 
            bc->selectedAln = 0;
          
          g_array_remove_index(bc->alignArr, i);
          bc->saved = 0;
	}
      else 
	{
	  i++;
	}
    }

    g_message_info("%d gappy sequences removed.  %d seqs left.\n\n", n, bc->alignArr->len);

    arrayOrder(bc->alignArr);
}


/* Remove empty (gappy) columns if the 'remove empty columns' option
 * is enabled. */
void rmFinaliseGapRemoval(BelvuContext *bc)
{
  if (bc->rmEmptyColumnsOn) 
    rmEmptyColumns(bc, 1.0);
  
  rmFinalise(bc);
}


/* Get rid of seqs that are more than x% identical with others. 
 * Keep the  first one.
 */
void mkNonRedundant(BelvuContext *bc, const double cutoff)
{
  int i=0,j=0, n=0, overhang;
  ALN *alni=NULL, *alnj=NULL;
  double id = 0.0;

  for (i = 0; i < (int)bc->alignArr->len; i++) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
      char *alniSeq = alnGetSeq(alni);

      for (j = 0; j < (int)bc->alignArr->len; j++) 
        {
          if (i == j) 
            continue;

          alnj = g_array_index(bc->alignArr, ALN*, j);
          char *alnjSeq = alnGetSeq(alnj);
	    
          overhang = alnOverhang(alnjSeq, alniSeq);
          id = identity(alniSeq, alnjSeq, bc->penalize_gaps);

          if (!overhang && id > cutoff)
	    {
              g_message_info("%s/%d-%d and %s/%d-%d are %.1f%% identical. "
                             "The first includes the latter which was removed.\n",
                             alni->name, alni->start, alni->end,
                             alnj->name, alnj->start, alnj->end,
                             id);
		
              /* Remove entry j */
              n++;
		
              if (bc->selectedAln == alnj) 
                bc->selectedAln = NULL;
              
              g_array_remove_index(bc->alignArr, j);
              bc->saved = FALSE;
              
              if (j < i) 
                {
                  i--;
                  alni = g_array_index(bc->alignArr, ALN*, i);
		}
              j--;
	    }
	}
    }

  g_message_info("%d sequences removed at the %.0f%% level.  %d seqs left.\n\n", n, cutoff, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Get rid of seqs that are less than x% identical with any of the others. 
 */
void rmOutliers(BelvuContext *bc, const double cutoff)
{
  int i=0,j=0, n=0;
  ALN *alni=NULL, *alnj=NULL;
  
  for (i = 0; i < (int)bc->alignArr->len - 1; ) 
    {
      alni = g_array_index(bc->alignArr, ALN*, i);
    
      double maxid = 0;
      for (maxid=0, j = 0; j < (int)bc->alignArr->len; ++j) 
	{
	  if (i == j) 
	    continue;
	
	  alnj = g_array_index(bc->alignArr, ALN*, j);
	  double id = identity(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
	
	  if (id > maxid) 
	    maxid = id;
	}
    
      if (maxid < cutoff) 
	{
	  g_message("%s/%d-%d was max %.1f%% identical to any other sequence and was removed.\n",
		    alni->name, alni->start, alni->end, maxid);
	
	  /* Remove entry */
	  n++;
	
	  if (bc->selectedAln == alni) 
	    bc->selectedAln = NULL;
	
	  g_array_remove_index(bc->alignArr, i);
	  bc->saved = FALSE;
	}
      else 
	{
	  i++;
	}
    }
  
  g_message("%d sequences removed at the %.0f%% level.  %d seqs left.\n\n", n, cutoff, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Remove sequences that have a score below the given value */
void rmScore(BelvuContext *bc, const double cutoff)
{
  scoreSort(bc);
  
  /* Save bc->selectedAln */
  ALN aln;
  initAln(&aln);
  
  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  int numRemoved = 0;
  int i = 0;
  
  for (i = 0; i < (int)bc->alignArr->len; ) 
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
    
      if (alnp->score < cutoff) 
	{ 
	  ++numRemoved;
	  
	  g_message_info("Removing %s/%d-%d (score %.1f)\n",
                         alnp->name, alnp->start, alnp->end, alnp->score);
	
	  if (bc->selectedAln == alnp) 
	    bc->selectedAln = NULL;
	  
	  g_array_remove_index(bc->alignArr, i);
	  bc->saved = FALSE;
	}
      else
	{
	  i++;
	}
    }
  
  arrayOrder(bc->alignArr);
  
  g_message_info("%d sequences with score < %.1f removed.  %d seqs left.\n\n", numRemoved, cutoff, bc->alignArr->len);
  
  /* Find bc->selectedAln in new array */
  if (bc->selectedAln) 
    { 
      int ip = 0;
      ALN aln;
      initAln(&aln);
      
      if (!alnArrayFind(bc->alignArr, &aln, &ip, scoreorder)) 
	{
	  bc->selectedAln = NULL;
	}
      else
	{
	  bc->selectedAln = g_array_index(bc->alignArr, ALN*, ip);
	}
  }
  
  bc->alignYStart = 0;
  
  rmFinaliseGapRemoval(bc);
}


/***********************************************************
 *		      Read/write files 			   *
 ***********************************************************/

/* Save the given tree to the given file in New Hampshire format */
void saveTreeNH(Tree *tree, TreeNode *node, FILE *file)
{
  if (!node) 
    return;
  
  if (node->left && node->right) 
    {
      fprintf(file, "(\n");
      saveTreeNH(tree, node->left, file);
      fprintf(file, ",\n");
      saveTreeNH(tree, node->right, file);
      fprintf(file, ")\n");
      
      if (node != tree->head)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
        fprintf(file, "%.0f", node->boot+0.5);
    }
  else
    {
      fprintf(file, "%s", node->name);
    }
  
  if (node != tree->head)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
    fprintf(file, ":%.3f", node->branchlen/100.0);
}


void writeMSF(BelvuContext *bc, FILE *pipe) /* c = separator between name and coordinates */
{
  int i=0, j=0;
  int maxfullnamelen = bc->maxNameLen + 1 + bc->maxStartLen + 1 + bc->maxEndLen;
  int paragraph=0, alnstart=0, alnend=0, alnlen=0, linelen=50, blocklen=10;
  
  if (bc->saveCoordsOn)
    maxfullnamelen = bc->maxNameLen + 1 + bc->maxStartLen + 1 + bc->maxEndLen;
  else
    maxfullnamelen = bc->maxNameLen;
  
  /* Title */
  fprintf(pipe, "PileUp\n\n %s  MSF: %d  Type: X  Check: %d  ..\n\n",
          bc->Title, bc->maxLen, GCGgrandchecksum(bc));

  /* Names and checksums */
  for(i=0; i < (int)bc->alignArr->len; ++i) 
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (alnp->markup) 
        continue;

      char *tmpStr = bc->saveCoordsOn 
        ? g_strdup_printf("%s%c%d-%d", alnp->name, bc->saveSeparator, alnp->start, alnp->end)
        : g_strdup_printf("%s", alnp->name);

      fprintf(pipe, "  Name: %-*s  Len:  %5d  Check:  %5d  Weight: %.4f\n",
              maxfullnamelen,
              tmpStr,
              bc->maxLen,
              GCGchecksum(bc, alnGetSeq(alnp)),
              1.0);

      g_free(tmpStr);
    }
  
  fprintf(pipe, "\n//\n\n");

  /* Alignment */
  while (paragraph * linelen < bc->maxLen)
    {
      for (i = 0; i < (int)bc->alignArr->len; ++i)
	{
          alnstart = paragraph * linelen; 
          alnlen = ( (paragraph+1)*linelen < bc->maxLen ? linelen : bc->maxLen - alnstart );
          alnend = alnstart + alnlen;

          ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
          
          if (alnp->markup) 
            continue;

          char *tmpStr = bc->saveCoordsOn 
            ? g_strdup_printf("%s%c%d-%d", alnp->name, bc->saveSeparator, alnp->start, alnp->end)
            : g_strdup_printf("%s", alnp->name);
          
          fprintf(pipe, "%-*s  ", maxfullnamelen, tmpStr);
          g_free(tmpStr);
          
          for (j = alnstart; j < alnend; ) 
            {
              char *alnpSeq = alnGetSeq(alnp);
              fprintf(pipe, "%c", alnpSeq[j]);
              j++;
              
              if ( !((j-alnstart) % blocklen) ) 
                fprintf(pipe, " ");
	    }
          
          fprintf(pipe, "\n");
	}
      
      fprintf(pipe, "\n");
      paragraph++;
    }

  fclose(pipe);
  fflush(pipe);
  
  bc->saved = TRUE;
}


static void readMSF(BelvuContext *bc, FILE *pipe)
{
  char seq[1001], *cp=NULL, *cq=NULL;
  char line[MAXLENGTH + 1];
  line[0] = 0;
  seq[0] = 0;
  
  g_message_info("\nDetected MSF format\n");

  /* Read sequence names */
  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe)) 
        break;

      if (!strncmp(line, "//", 2)) 
        {
          break;
        }
      else if (strstr(line, "Name:") && (cp = strstr(line, "Len:")) && strstr(line, "Check:")) 
	{
          int len;

          sscanf(cp+4, "%d", &len);

          if (bc->maxLen < len) 
            bc->maxLen = len;

          cp = strstr(line, "Name:")+6;
          
          while (*cp && *cp == ' ') 
            cp++;

          /* Read sequence name and coords */
          ALN *aln = createEmptyAln();
          parseMulLine(bc, cp, aln);

          /* Create a string to contain the sequence data */
          aln->sequenceStr = g_string_new("");
          
          /* Add to the array */
          aln->nr = bc->alignArr->len + 1;

          g_array_append_val(bc->alignArr, aln);
          g_array_sort(bc->alignArr, alphaorder);
	}
    }

  /* Read sequence alignment */
  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe)) 
        break;

      cp = line;
      cq = seq;
      while (*cp && *cp == ' ') cp++;
      while (*cp && *cp != ' ') cp++; /* Spin to sequence */
      
      while (*cp && cq-seq < 1000) 
        {
          if (isAlign(*cp)) *cq++ = *cp;
          cp++;
	}
	
      *cq = 0;

      if (*seq) 
        {
          cp = line;
          while (*cp && *cp == ' ') 
            cp++;

          /* Get the name of the sequence  */
          ALN aln;
          initAln(&aln);
          parseMulLine(bc, cp, &aln);

          /* Get the alignment with this name from the array, and append this
           * bit of sequence to it. */
          int ip = 0;
          if (alignFind(bc->alignArr, &aln, &ip))
            {
              ALN *alnp = g_array_index(bc->alignArr, ALN*, ip);
              g_string_append(alnp->sequenceStr, seq);
            }
          else
            {
              g_warning("Alignment does not appear in file header section (%s %d %d).\n", aln.name, aln.start, aln.end);
            }
	}
    }
  
  bc->saveFormat = BELVU_FILE_MSF;
}


/* Parse annotation lines from #=GS lines in the Mul file format */
static void parseMulAnnotationLine(BelvuContext *bc, const char *seqLine)
{
  const char *cp = seqLine;

  char
    *namep,		/* Seqname */
    name[MAXNAMESIZE+1],/* Seqname */
    *labelp,		/* Label (OS) */
    *valuep;		/* Organism (c. elegans) */

  
  /* Ignore anything that's not a 'GS' line */
  if (strncmp(cp, "#=GS", 4) != 0)
    return;
      
  /* Copy the portion after  the '#=GS' into 'name' */
  strcpy(name, cp+4);
  namep = name;
  
  /* Ignore leading whitespace */
  while (*namep == ' ') 
    namep++;
      
  /* Find the end of the name */
  labelp = namep;
  while (*labelp != ' ') 
    labelp++;
  
  *labelp = 0;		/* Terminate namep */
  ++labelp;		/* Walk across terminator */
  
  while (*labelp == ' ') /* Find the start of the label */
    ++labelp;
      
  valuep = labelp;
  while (*valuep != ' ') valuep++; /* Find the end of the label */
  while (*valuep == ' ') valuep++; /* Find the start of the value */
  
  if (strncasecmp(labelp, "LO", 2) == 0)
    {
      int i, colornr;
      
      /* Check the value to see if it matches one of our colour names */
      for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
        {
          if (strcasecmp(colorNames[i], valuep) == 0) 
            {
              colornr = i;
              break;
            }
        }
      
      if (colornr == -1)
        {
          g_message("Unrecognized color: %s\n", valuep);
          colornr = 0;
        }
          
      /* Create an ALN struct from the name */
      ALN aln;
      initAln(&aln);
      str2aln(bc, namep, &aln);

      /* Find the corresponding sequence */
      if (!alnArrayFind(bc->alignArr, &aln, &i, alphaorder))
        {
          g_critical("Cannot find '%s' [%d %d] in alignment.\n", aln.name, aln.start, aln.end);
          return;
        }
      
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      alnp->color = colornr;
    }
  else if (strncasecmp(labelp, bc->organismLabel, 2) == 0)
    {
      /* Add organism to table of organisms.  This is necessary to make all 
         sequences of the same organism point to the same place and to make a 
         non-redundant list of present organisms */
      
      /* Find string in permanent stack */
      if (!(valuep = strstr(cp, valuep)))
        {
          g_critical("Failed to parse organism properly.\n");
          return;
        }
      
      /* Create an ALN struct from this organism */
      ALN *aln = createEmptyAln();
      
      aln->organism = valuep; /* Find organism string in permanent stack */
      g_array_append_val(bc->organismArr, aln);
      g_array_sort(bc->organismArr, organism_order);
          
      if (strchr(cp, '/') && strchr(cp, '-'))
        {
          str2aln(bc, namep, aln);
              
          /* Find the corresponding sequence */
          int ip = 0;
          
          if (!alnArrayFind(bc->alignArr, aln, &ip, alphaorder)) 
            {
              g_critical("Cannot find '%s' [%d %d] in alignment.\n", aln->name, aln->start, aln->end);
              return;
            }
             
          ALN *alnp = g_array_index(bc->alignArr, ALN*, ip);
              
          /* Store pointer to unique organism in ALN struct */
          aln->organism = valuep;
          ip = 0;
          
          alnArrayFind(bc->organismArr, aln, &ip, organism_order);
          alnp->organism = g_array_index(bc->organismArr, ALN*, ip)->organism;
        }
      else 
        {
          /* Organism specified for entire sequences.
             Find all ALN instances of this sequences.
          */
          int i = 0;
          for (i = 0; i < (int)bc->alignArr->len; ++i) 
            {
              ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

              if (strcmp(alnp->name, namep) == 0) 
                {
                  /* Store pointer to unique organism in ALN struct */
                  ALN aln;
                  initAln(&aln);
                  aln.organism = valuep;
                  int ip = 0;

                  alnArrayFind(bc->organismArr, &aln, &ip, organism_order);
                  alnp->organism = g_array_index(bc->organismArr, ALN*, ip)->organism;
                }
            }
        }
    }
}


/* Add sequence data to an alignment from a line of text of the format:
 * SEQ_NAME     SEQUENCE_DATA,
 * where alnstart gives the position of the start of SEQUENCE_DATA. */
static void appendSequenceDataToAln(BelvuContext *bc, char *line, const int alnstart)
{
  /* Find the alignment name */
  ALN aln;
  initAln(&aln);
  parseMulLine(bc, line, &aln);

  /* See if this alignment is in the alignments array */
  int ip = 0;
  if (alnArrayFind(bc->alignArr, &aln, &ip, alphaorder))
    {
      /* Append this bit of sequence to the existing alignment */
      ALN *alnp = g_array_index(bc->alignArr, ALN*, ip);
      g_string_append(alnp->sequenceStr, line + alnstart);
      
      /* Recalculate the max alignment length */
      if ((int)alnp->sequenceStr->len > bc->maxLen)
        bc->maxLen = alnp->sequenceStr->len;
    }
  else
    {
      /* Create a new alignment and add this bit of sequence */
      ALN *alnp = createEmptyAln();
      alncpy(alnp, &aln);
      
      alnp->sequenceStr = g_string_new(line + alnstart);
      alnp->nr = bc->alignArr->len + 1;
      
      g_array_append_val(bc->alignArr, alnp);
      g_array_sort(bc->alignArr, alphaorder);
      
      /* Recalculate the max alignment length */
      if ((int)alnp->sequenceStr->len > bc->maxLen)
        bc->maxLen = alnp->sequenceStr->len;
    }
}


/* ReadMul 
 * parses an alignment file in mul (stockholm) or selex format and puts it in the Align array
 *
 * Assumes header contains ' ' or '#'
 *
 * Alignment can have one of the following formats:
 *  CSW_DROME  VTHIKIQNNGDFFDLYGGEKFATLP
 *  CSW_DROME  51   75   VTHIKIQNNGDFFDLYGGEKFATLP
 *  CSW_DROME  51   75   VTHIKIQNNGDFFDLYGGEKFATLP P29349
 *  KFES_MOUSE/458-539    .........WYHGAIPW.....AEVAELLT........HTGDFLVRESQG
 *
 */
static void readMul(BelvuContext *bc, FILE *pipe)
{
  char line[MAXLENGTH+1];
  line[0] = 0;

  /* Read raw alignment into stack
   *******************************/
  
  int alnstart = MAXLENGTH;
  GSList *alnList = NULL;
  
  while (!feof (pipe))
    { 
      /* EOF checking to make acedb calling work */
      if (!fgets (line, MAXLENGTH, pipe) || (unsigned char)*line == (unsigned char)EOF ) 
        break;
      
      if (!strncmp(line, "PileUp", 6)) 
        {
          readMSF(bc, pipe);
          return;
        }

      /* Remove any trailing newline */
      char *cp = strchr(line, '\n');
      if (cp) 
        *cp = 0;
      
      int lineLen = strlen(line);

      if (lineLen > 0 && *line != '#' && strcmp(line, "//") != 0)
	{
          /* Sequence line */
          char *cq = strchr(line, ' ');
          
          if (!cq) 
            g_error("Error reading selex file; no spacer between name and sequence in the following line:\n%s", line);
          
          /* Find which column the alignment starts in */
          int i = 0;
          for (i = 0; line[i] && line[i] != ' '; ++i);
          for (; line[i] && !isAlign(line[i]); ++i);
          
          /* Remember the leftmost start position of any alignment. We'll assume
           * all alignments start in the same column. */
          if (i < alnstart) 
            alnstart = i;

          /* Remove optional accession number at end of alignment */
          /* FOR PRODOM STYLE ALIGNMENTS W/ ACCESSION TO THE RIGHT - MAYBE MAKE OPTIONAL
           * This way it's incompatible with alignments with ' ' gapcharacters */
          for (; line[i] && isAlign(line[i]); i++);
          line[i] = 0;

          /* Store the line for processing later (once alnstart has been calculated) */
          alnList = g_slist_prepend(alnList, g_strdup(line));
	}
      else if (!strncmp(line, "#=GF ", 5) || 
               !strncmp(line, "#=GS ", 5)) 
        {
	  /* Store all annotation lines in a list. Prepend the items because that
	   * is more efficient, and then reverse the list at the end */
	  bc->annotationList = g_slist_prepend(bc->annotationList, g_strdup(line));
          parseMulAnnotationLine(bc, line);
        }
      else if (!strncmp(line, "#=GC ", 5) || 
               !strncmp(line, "#=GR ", 5) || 
               !strncmp(line, "#=RF ", 5)) 
        {
          /* These are markup lines that are shown in the alignment list */
          alnList = g_slist_prepend(alnList, g_strdup(line));
        }
      else if (!strncmp(line, "# matchFooter", 13)) 
        {
          /* Match Footer  */
          bc->matchFooter = TRUE;
          break;
        }
    }
  
  /* Reverse the lists, because we prepended items instead of appending them */
  bc->annotationList = g_slist_reverse(bc->annotationList);
  alnList = g_slist_reverse(alnList);
  
  /* Loop through all of the alignment lines and extract the sequence string */
  GSList *alnItem = alnList;
  for ( ; alnItem; alnItem = alnItem->next)
    {
      char *alnStr = (char*)(alnItem->data);
      appendSequenceDataToAln(bc, alnStr, alnstart);
      
      /* Free the line string, now we're done with it */
      g_free(alnStr);
      alnItem->data = NULL;
    }

  g_slist_free(alnList);
  alnList = NULL;
  
/* For debugging * /
   for (i = 0; i < nseq; i++) {
   alnp = arrp(Align, i, ALN);
   printf("\n%-10s %4d %4d %s %d %s\n", 
   alnp->name, alnp->start, alnp->end, alnp->seq, 
   alnp->len, alnp->organism);
   }
   for (i=0; i < arrayMax(organismArr); i++) 
   printf("%s\n", arrp(organismArr, i, ALN)->organism);
   */
  
  if (bc->alignArr->len == 0 || bc->maxLen == 0) 
    g_error("Unable to read sequence data\n");
  
  bc->saveFormat = BELVU_FILE_MUL;
}


/* 
 Format the name/start-end string 
 
 For convenience, used by writeMul. Returns a newly allocated string
 which must be freed with g_free
 */
static char *writeMulName(BelvuContext *bc, ALN *aln) 
{  
  char *name = NULL; /* result */
  
  char *cp, *namep, GRname[MAXNAMESIZE*2+2], GRfeat[MAXNAMESIZE+1];
  
  if (aln->markup == GC) 
    {
      name = g_strdup_printf("#=GC %s", aln->name);
      return name;
    }
  
  /* NOTE: GR lines have the feature concatenated in aln->name - must separate */
  if (aln->markup == GR) 
    {
      namep = GRname;
      strncpy(namep, aln->name, 50);
      cp = strchr(namep, ' ');
      strncpy(GRfeat, cp+1, 50);
      *cp = 0;
    }
  else
    {
      namep = aln->name;
    }
  
  if (!bc->saveCoordsOn) 
    {
      name = g_strdup_printf("%s", namep);
    }
  else 
    {
      name = g_strdup_printf("%s%c%d-%d", namep, bc->saveSeparator, aln->start, aln->end);
    }
  
  if (aln->markup == GR)
    name = g_strdup_printf("#=GR %s %s", name, GRfeat);
  
  return name;
}


void writeMul(BelvuContext *bc, FILE *fil)
{
  /* Write Annotation */
  GSList *annotationItem = bc->annotationList;
  
  for ( ; annotationItem; annotationItem = annotationItem->next)
    {
      const char *line = (const char*)(annotationItem->data);
      fprintf(fil, "%s\n", line);
    }

  /* Find max width of name column */
  int i = 0;
  int maxWidth = 0;
  
  for ( ; i < (int)bc->alignArr->len; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

      char *mulName = writeMulName(bc, alnp);
      int curWidth = strlen(mulName);
      g_free(mulName);
    
      if ( curWidth > maxWidth) 
	maxWidth = curWidth;
    }
    
  /* Write alignment */
  for (i = 0; i < (int)bc->alignArr->len; i++)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      char *mulName = writeMulName(bc, alnp);
      const char *alnpSeq = alnGetSeq(alnp);
    
      fprintf(fil, "%-*s %s\n", maxWidth, mulName, (alnpSeq ? alnpSeq : ""));
    
      g_free(mulName);
    }
  
  fprintf(fil, "//\n");
  
  fclose(fil);
  fflush(fil);

  bc->saved = TRUE;
}


/* readFile 
 * Determines the format of the input file and calls the appropriate 
 * parser. Valid file formats are fasta, MSF, mul (stockholm) or selex.
 */
void readFile(BelvuContext *bc, FILE *pipe)
{
  char   ch = '\0';
  char line[MAXLENGTH+1];
  line[0] = 0;
  
  /* Parse header to check for MSF or Fasta format */
  while (!feof (pipe))
    { 
      ch = fgetc(pipe);
      
      if (!isspace(ch)) 
        {
          if (ch == '>') 
            {
              ungetc(ch, pipe);
              return readFastaAln(bc, pipe);
            }
          else
            {
              break;
            }
        }
      else if (ch == '\n')
        {
          if (!fgets (line, MAXLENGTH, pipe)) 
            break;
        }
      
      if (strstr(line, "MSF:") && strstr(line, "Type:") && strstr(line, "Check:")) 
        {
          return readMSF(bc, pipe);
        }
    }
  
  if (!feof(pipe))
    ungetc(ch, pipe);

  return readMul(bc, pipe);
}


void writeFasta(BelvuContext *bc, FILE *pipe)
{
  int i=0, n=0;
  char *cp = NULL;

  for (i = 0; i < (int)bc->alignArr->len; ++i) 
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      
      if (alnp->markup) 
        continue;
	
      if (bc->saveCoordsOn)
        fprintf(pipe, ">%s%c%d-%d\n", alnp->name, bc->saveSeparator, alnp->start, alnp->end);
      else
        fprintf(pipe, ">%s\n", alnp->name);
	
      for (n=0, cp = alnGetSeq(alnp); *cp; cp++)
        {
	  if (bc->saveFormat == BELVU_FILE_ALIGNED_FASTA)
            {
              fputc(*cp, pipe);
              n++;
	    }
          else
            {
              if (!isGap(*cp)) 
                {
                  fputc(*cp, pipe);
                  n++;
		}
	    }
	    
          if (n && !(n % 80) ) 
            {
              fputc('\n', pipe);
              n = 0;
	    }
	}
	
      if (n)
        fputc('\n', pipe);
    }

  fclose(pipe);
  fflush(pipe);
  
  bc->saved = TRUE;
}


static int getMatchStates(BelvuContext *bc)
{
  int i=0, j=0, n=0, retval=0;

  for (j = 0; j < bc->maxLen; ++j) 
    {
      n = 0;
      for (i = 0; i < (int)bc->alignArr->len; ++i) 
        {
          ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
          char *alnpSeq = alnGetSeq(alnp);
          
          if (isalpha(alnpSeq[j])) 
            n++;
	}
      
      if (n > (int)bc->alignArr->len / 2) 
        retval++;
    }

  return retval;
}


/* Output a.a. probabilities in HMMER format.

   Disabled counting of match-state residues only (script HMM-random does that already)

   Disabled pseudocounts - not sure what the best way is.
*/
void outputProbs(BelvuContext *bc, FILE *fil)
{
  /* From swissprot 33:

     Ala (A) 7.54   Gln (Q) 4.02   Leu (L) 9.31   Ser (S) 7.19
     Arg (R) 5.15   Glu (E) 6.31   Lys (K) 5.94   Thr (T) 5.76
     Asn (N) 4.54   Gly (G) 6.86   Met (M) 2.36   Trp (W) 1.26
     Asp (D) 5.29   His (H) 2.23   Phe (F) 4.06   Tyr (Y) 3.21
     Cys (C) 1.70   Ile (I) 5.72   Pro (P) 4.91   Val (V) 6.52
  */
  
  double f_sean[21] = {
    0.0, 
    .08713, .03347, .04687, .04953, .03977, .08861, .03362, .03689, .08048, .08536, 
    .01475, .04043, .05068, .03826, .04090, .06958, .05854, .06472, .01049, .02992
  };

  double p = 0.0;
  int i = 0, c[21], n=0, nmat = 0;
  char *cp = NULL;

  for (i = 1; i <= 20; i++) 
    c[i] = 0;

  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
      cp = alnGetSeq(alnp);

      for (; *cp; cp++)
        {
          if (a2b_sean[(unsigned char)(*cp)])
            {
              c[a2b_sean[(unsigned char)(*cp)]]++;
              n++;
            }
        }
    }

  if (!n)
    g_error("No residues found\n");

  if (0) 
    {
      nmat = getMatchStates(bc);	/* Approximation, HMM may differ slightly */

      if (!nmat) 
        g_error("No match state columns found\n");
      
      g_message("Amino\n");
      
      for (i = 1; i <= 20; i++) 
        {
          /* One way:  p = ((double)c[i]/nmat + 20*f_sean[i]) / (double)(n/nmat+20);*/
          /* Other way: */
          p = ((double)c[i] + 20*f_sean[i]) / (double)(n+20);

          if (p < 0.000001) 
            p = 0.000001; /* Must not be zero */

          g_message("%f ", p);
          /* printf("   %d %f %d \n", c[i], f[i], n);*/
        }
    }
  else
    {
      g_message("Amino\n");
      for (i = 1; i <= 20; ++i) 
        {
          p = (double)c[i] / (double)n;
          if (p < 0.000001) p = 0.000001; /* Must not be zero */
          g_message("%f ", p);
	}
    }
    
  g_message("\n");
  
  fclose(fil);
  fflush(fil);

  bc->saved = TRUE;
}



/* Calculate the identity of each sequence against each other, and print the
 * results to stdout. */
void listIdentity(BelvuContext *bc)
{
  /* This may take some time, so feed back to the user what is happening. */
  g_message_info("Outputting identities...\n");
  setBusyCursor(bc, TRUE);

  int i=0,j=0,n=0 ;
  double totsc=0, maxsc=0, minsc=1000000,
         totid=0.0, maxid=0.0, minid=100.0;
  
  for (i = n = 0; i < (int)bc->alignArr->len - 1; ++i) 
    {
      /* Update the display periodically, otherwise the busy cursor might not
       * get updated (e.g. if we were outside the main window when we set the 
       * cursor and the user later enters the window). */
      while (gtk_events_pending())
        gtk_main_iteration();

      ALN *alni = g_array_index(bc->alignArr, ALN*, i);
      
      if (alni->markup > 0) /* ignore markup lines */
        continue;

      for (j = i+1; j < (int)bc->alignArr->len; ++j) 
        {
          ALN *alnj = g_array_index(bc->alignArr, ALN*, j);

          if (alnj->markup > 0) /* ignore markup lines */
            continue;

          double id = identity(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
          totid += id;

          if (id > maxid) 
            maxid = id;
          
          if (id < minid) 
            minid = id;
      
          double sc = score(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
          totsc += sc;
          
          if (sc > maxsc) 
            maxsc = sc;
          
          if (sc < minsc) 
            minsc = sc;
          
          g_message("%s/%d-%d and %s/%d-%d are %.1f%% identical, score=%f\n",
                    alni->name, alni->start, alni->end,
                    alnj->name, alnj->start, alnj->end,
                    id, sc);
          
          ++n;
        }
      g_message("\n");
    }

  g_message("Maximum %%id was: %.1f\n", maxid);
  g_message("Minimum %%id was: %.1f\n", minid);
  g_message("Mean    %%id was: %.1f\n", totid/n);
  g_message("Maximum score was: %.1f\n", maxsc);
  g_message("Minimum score was: %.1f\n", minsc);
  g_message("Mean    score was: %.1f\n", (double)totsc/n);

  setBusyCursor(bc, FALSE);
  g_message_info("Finished outputting identities.\n");
}


static void onDestroySpawnedWindow(GtkWidget *window, gpointer data)
{
  /* We must remove the window from the list of spawned windows */
  BelvuContext *bc = (BelvuContext*)data;
  bc->spawnedWindows = g_slist_remove(bc->spawnedWindows, window);
}


void fetchAln(BelvuContext *bc, ALN *alnp)
{
  GError *error = NULL;

  if (bc->useWWWFetch)
    {
      /* Get the URL from the enviroment var (or use the hard-coded default
       * value) */
      char  *env, url[1025]="";
      
      if ((env = getenv(FETCH_URL_ENV_VAR)) )
	strcpy(url, env);
      else
	strcpy(url, DEFAULT_FETCH_URL);

      /* Add the sequence name to the url. The URL should be a format string
       * containing '%s' for the name. */
      char *cp = strchr(url, '%');
      
      if (cp)
        ++cp;
      
      if (!cp || *cp != 's')
        {
          g_critical("Invalid URL string %s.\nThe URL must contain the search string \"%%s\", which will be replaced with the sequence name.\n", url);
        }
      else
        {
          char *link = g_strdup_printf(url, alnp->name);
          
          g_message_info("Opening URL: %s\n", link);
          seqtoolsLaunchWebBrowser(link, &error);
          
          g_free(link);
        }
    }
  else
    {
      char  *env, fetchProg[1025]="";
      
      if ((env = getenv(FETCH_PROG_ENV_VAR)) )
        strcpy(fetchProg, env);
      else
        strcpy(fetchProg, DEFAULT_FETCH_PROG);

      char *path = g_find_program_in_path(fetchProg);
      
      if (!path)
        {
          g_warning("Executable '%s' not found in path: %s\n", fetchProg, getenv("PATH"));
        }
      else
        {
          g_free(path);
          
          char *cmd = g_strdup_printf("%s '%s' &", fetchProg, alnp->name);
      
          GtkWidget *pfetchWin = externalCommand(cmd, BELVU_TITLE, bc->belvuAlignment, &error);
      
          if (pfetchWin)
            {
              const gchar *env = g_getenv(FONT_SIZE_ENV_VAR);
              if (env)
                widgetSetFontSizeAndCheck(pfetchWin, convertStringToInt(env));
              
              /* Add the window to our list of spawned windows */
              bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, pfetchWin);
              g_signal_connect(G_OBJECT(pfetchWin), "destroy", G_CALLBACK(onDestroySpawnedWindow), bc);
            }
      
          g_free(cmd);
        }
    }
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

}



/* Utility to return true if the given alignment is highlighted (i.e. has the
 * same name as the selected alignment) */
gboolean alignmentHighlighted(BelvuContext *bc, ALN *alnp)
{
  /* Return true if this alignment has the same name as the selected alignment */
  return (bc->selectedAln && stringsEqual(alnp->name, bc->selectedAln->name, TRUE));
}


/***********************************************************
 *                        Cursors                          *
 ***********************************************************/

/* Utility to set/unset the busy cursor on all open toplevel windows */
void setBusyCursor(BelvuContext *bc, const gboolean busy)
{
  GdkCursor *cursor = bc->defaultCursor;

  if (busy)
    cursor = bc->busyCursor;
  else if (bc->removingSeqs)
    cursor = bc->removeSeqsCursor;

  if (bc->belvuWindow && bc->belvuWindow->window)
    gdk_window_set_cursor(bc->belvuWindow->window, cursor);

  if (bc->belvuTree && bc->belvuTree->window)
    gdk_window_set_cursor(bc->belvuTree->window, cursor);

  if (bc->consPlot && bc->consPlot->window)
    gdk_window_set_cursor(bc->consPlot->window, cursor);

  if (bc->orgsWindow && bc->orgsWindow->window)
    gdk_window_set_cursor(bc->orgsWindow->window, cursor);

  /* Force cursor to change immediately */
  while (gtk_events_pending())
    gtk_main_iteration();
}

