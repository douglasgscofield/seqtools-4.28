/*  File: blxview.c
 *  Author: Erik Sonnhammer, 1993-02-20
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
 * Description: Setup and miscellaneous functions for Blixem.
 *               
 *              This file is a bit of a mish-mash. It originally contained all
 *              of the graphical stuff for Blixem, but that has been moved to
 *              other files. This file now mostly contains setup functions,
 *              but it also contains other functions that don't really belong
 *              anywhere else.
 *----------------------------------------------------------------------------
 */


/*
MSP score codes (for obsolete exblx file format):
----------------
-1 exon                  -> Big picture + Alignment
-2 intron                -> Big picture + Alignment
-3 Any colored segment   -> Big picture
-4 stringentSEGcolor     -> Big picture
-5 mediumSEGcolor        -> Big picture
-6 nonglobularSEGcolor   -> Big picture
-10 hidden by hand
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <gdk/gdkkeysyms.h>

#include <blixemApp/blixem_.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/blxdotter.h>
#include <blixemApp/sequencecellrenderer.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/blxparser.h>
#include <seqtoolsUtils/utilities.h>


#define MAXALIGNLEN			      10000
#define ORGANISM_PREFIX_SEPARATOR	      "  "
#define GFF3_VERSION_HEADER		      "##gff-version 3"	  /* the header line of a GFF v3 file */
#define GFF3_SEQUENCE_REGION_HEADER	      "##sequence-region" /* the start comment of the sequence-region comment line in a GFF v3 file */
#define MIN_GAP_HIGHLIGHT_WIDTH		      5			  /* minimum width of assembly gaps markers */


static void            blviewCreate(char *align_types, const char *paddingSeq, GArray* featureLists[], GList *seqList, GSList *supportedTypes, CommandLineOptions *options, const gboolean External, GSList *styles) ;
static void            processGeneName(BlxSequence *blxSeq);
static void            processOrganism(BlxSequence *blxSeq);


/* GLOBAL VARIABLES... sigh... */

GtkWidget *blixemWindow = NULL ;
static char *padseq = 0;



/*
 *                Internal routines
 */


//static int possort(MSP *msp1, MSP *msp2) {
//    return ( (msp1->qstart > msp2->qstart) ? 1 : -1 );
//}
//
//static int namesort(MSP *msp1, MSP *msp2) {
//    if (strcmp(msp1->sname, msp2->sname))
//	return strcmp(msp1->sname, msp2->sname);
//    else
//	return possort(msp1, msp2);
//}
//static int scoresort(MSP *msp1, MSP *msp2) {
//    return ( (msp1->score < msp2->score) ? 1 : -1 );
//}
//
//static int idsort(MSP *msp1, MSP *msp2)
//{
//  int result = 0 ;
//
//  result = ( (msp1->id < msp2->id) ? 1 : -1 ) ;
//
//  return result ;
//}


//static void wholePrint(void)
//{
   // int
//	tmp,
//	dispstart_save = dispstart,
//	BigPictON_save = BigPictON;
//
//    static int
//	start=1, end=0;
//    ACEIN zone_in;
//
//    if (!end) end = qlen;
//
//    /* Swap coords if necessary */
//    if ((plusmin < 0 && start < end) || (plusmin > 0 && start > end )) {
//	tmp = start;
//	start = end;
//	end = tmp;
//    }
//
//    /* Apply max limit MAXALIGNLEN */
//    if ((abs(start-end)+1) > MAXALIGNLEN*symbfact) {
//	start = dispstart - plusmin*MAXALIGNLEN*symbfact;
//	if (start < 1) start = 1;
//	if (start > qlen) start = qlen;
//
//	end = start + plusmin*MAXALIGNLEN*symbfact;
//	if (end > qlen) end = qlen;
//	if (end < 1) end = 1;
//    }
//
//    if (!(zone_in = messPrompt("Please state the zone you wish to print",
//			       messprintf("%d %d", start, end),
//			       "iiz", 0)))
//      return;
//
//    aceInInt(zone_in, &start);
//    aceInInt(zone_in, &end);
//    aceInDestroy (zone_in);
//
//    dispstart = start;
//    displen = abs(end-start)+1;
//
//    /* Validation */
//    if (plusmin > 0 && start > end) {
//	g_message("Please give a range where from: is less than to:");
//	return;
//    }
//    else if (plusmin < 0 && start < end) {
//	g_message("Please give a range where from: is more than to:");
//	return;
//    }
//    if (displen/symbfact > MAXALIGNLEN) {
//	g_message("Sorry, can't print more than %d residues.  Anyway, think of the paper!", MAXALIGNLEN*symbfact);
//	return;
//    }
//
//    wholePrintOn = 1;
//    BigPictON = 0;
////    oneGraph = 1;
//
//    blviewRedraw();
//    graphPrint();
//
//    /* Restore */
//    wholePrintOn = 0;
//    dispstart = dispstart_save;
//    displen = dispstart_save;
////    oneGraph = 0;
//    BigPictON = BigPictON_save;
//
//    blviewRedraw();
//}

//static void blxPrint(void)
//{
//  oneGraph = 1;
//  blviewRedraw();
//  graphPrint();
//  oneGraph = 0;
//
//  blviewRedraw();
//}


/* Check we have everything we need */
static void validateInput(CommandLineOptions *options)
{
  if (!options->refSeq)
    {
      g_error("Error: reference sequence not specified.\n");
    }
  
  /* if we were given the blast mode instead of seq type, determine the display sequence type */
  if (options->seqType == BLXSEQ_NONE && options->blastMode != BLXMODE_UNSET)
    {
      if (options->blastMode == BLXMODE_BLASTN)
        options->seqType = BLXSEQ_DNA;
      else
        options->seqType = BLXSEQ_PEPTIDE;
    }
  
  /* Check we have a valid seq type */
  if (options->seqType == BLXSEQ_NONE)
    {
      g_message("\nNo sequence type specified. Detected ");
      
      GError *error = NULL;
      BlxSeqType seqType = determineSeqType(options->refSeq, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
      
      if (seqType == BLXSEQ_PEPTIDE)
	{
	  g_message("protein sequence. Will try to run Blixem in protein mode.\n");
	  options->seqType = BLXSEQ_PEPTIDE;
	}
      else
	{
	  g_message("nucleotide sequence. Will try to run Blixem in nucelotide mode.\n");
	  options->seqType = BLXSEQ_DNA;
	}
    }
  
  /* Ideally we'd get rid of blast mode and just deal with real sequence types (ref seq type, match
   * seq type and display type), but it's too in-built for now so make sure it's set */
  if (options->blastMode == BLXMODE_UNSET)
    options->blastMode = (options->seqType == BLXSEQ_DNA ? BLXMODE_BLASTN : BLXMODE_BLASTX);
  
  /* Set the number of reading frames */
  if (options->seqType == BLXSEQ_PEPTIDE)
    options->numFrames = 3;
  else
    options->numFrames = 1;
  
  /* Now we've got the ref seq, we can work out the zoom, if zooming to view the whole sequence */
  if (options->zoomWhole)
    {
      options->bigPictZoom = strlen(options->refSeq);
    }
  else
    {
      if (options->seqType == BLXSEQ_DNA)
        options->bigPictZoom = 30;
      else
        options->bigPictZoom = 10;
    }  
}




/* Utility to get the parent sequence of the given variant, if it is not already set.
 * returns true if the returned parent is not null. */
static gboolean getParent(BlxSequence *variant, BlxSequence **parent, GList *allSeqs)
{
  if (*parent == NULL)
    {
      *parent = blxSequenceGetVariantParent(variant, allSeqs);
    }
  
  return (*parent != NULL);
}


/* This function checks if the given sequence is missing its optional data (such as organism
 * and gene name) and, if so, looks for the parent sequence and copies the data from there. */
static void populateMissingDataFromParent(BlxSequence *curSeq, GList *seqList, GList *columnList)
{
  BlxSequence *parent = NULL;
  GList *item = columnList;
  
  for ( ; item; item = item->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
      GValue *value = blxSequenceGetValue(curSeq, columnInfo->columnId);
      
      if (!value)
        {
          getParent(curSeq, &parent, seqList);
          GValue *parentValue = blxSequenceGetValue(parent, columnInfo->columnId);
          blxSequenceSetValue(curSeq, columnInfo->columnId, parentValue);
        }
    }
}


/* Merge new features into our existing context: merges the newMsps list into mspList and newSeqs
 * list into seqList. Takes ownership of the contents of both newMsps and newSeqs. */
void blxMergeFeatures(MSP *newMsps, GList *newSeqs, MSP **mspList, GList **seqList)
{
  /* Append new MSPs to MSP list */
  MSP *lastMsp = *mspList;
  
  while (lastMsp->next)
    lastMsp = lastMsp->next;
  
  lastMsp->next = newMsps;
  
  /* Append new sequences to sequence list */
  *seqList = g_list_concat(*seqList, newSeqs);
}


/* This function to loads the contents of a natively-supported features file
 * (e.g. GFF) into blixem. The new features are appended onto the existing
 * sequence and MSP lists.
 * A non-local or non-native file can be passed to check if it is loadable
 * and, if not, this function returns early and sets the error. */
void loadNativeFile(const char *filename,
                    const char *buffer,
                    GKeyFile *keyFile,
                    BlxBlastMode *blastMode,
                    GArray* featureLists[],
                    GSList *supportedTypes, 
                    GSList *styles,
                    MSP **newMsps,
                    GList **newSeqs,
                    GList *columnList,
                    GHashTable *lookupTable,
                    const int refSeqOffset,
                    const IntRange* const refSeqRange,
                    GError **error)
{
  if (!filename && !buffer)
    {
      g_set_error(error, BLX_ERROR, 1, "No file or buffer provided.");
      return;
    }

  char *dummyseq1 = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname1[FULLNAMESIZE+1] = "";
  char *dummyseq2 = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char dummyseqname2[FULLNAMESIZE+1] = "";

  /* The range passed to the parser must be in the original parsed coords, so subtract the
   * offset. (It is used to validate that the range in the GFF file we're reading in is within 
   * blixem's known range.)  */
  IntRange toplevelRange = {UNSET_INT, UNSET_INT};
  
  if (refSeqRange)
    {
      toplevelRange.min = refSeqRange->min - refSeqOffset;
      toplevelRange.max = refSeqRange->max - refSeqOffset;
    }
  
  if (filename)
    {
      /* Open the file for reading */
      FILE *file = fopen(filename, "r");

      if (!file)
        {
          g_set_error(error, BLX_ERROR, 1, "Error opening file '%s' for reading.\n", filename);
        }
      else
        {
          parseFS(newMsps, file, blastMode, featureLists, newSeqs, columnList, supportedTypes, styles,
                  &dummyseq1, dummyseqname1, &toplevelRange, &dummyseq2, dummyseqname2, keyFile, lookupTable, error) ;      
          
          fclose(file);
        }
    }
  else if (buffer)
    {
      parseBuffer(newMsps, buffer, blastMode, featureLists, newSeqs, columnList, supportedTypes, styles,
                  &dummyseq1, dummyseqname1, &toplevelRange, &dummyseq2, dummyseqname2, keyFile, lookupTable, error) ;      
    }
}


/* Once we've fetched all the sequences we need to do some post-processing. Loop 
 * twice: once to modify any fields in our own custom manner, and once more to
 * see if any with missing data can copy it from their parent. (Need to do these in
 * separate loops or we don't know if the data we're copying is processed or not.) */
void finaliseFetch(GList *seqList, GList *columnList)
{
  GList *seqItem = seqList;

  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      processGeneName(blxSeq);
      processOrganism(blxSeq);
    }
  
  for (seqItem = seqList; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      populateMissingDataFromParent(blxSeq, seqList, columnList);
    }
}


/* Find any gaps in the reference sequence */
static void findAssemblyGaps(const char *refSeq, GArray *featureLists[], MSP **mspList, const IntRange* const refSeqRange)
{
  /* Find the last msp in the list so we can append new ones to the end */
  MSP *lastMsp = *mspList;
  while (lastMsp && lastMsp->next)
    lastMsp = lastMsp->next;
  
  /* Scan for the gap character */
  const char *cp = strchr(refSeq, SEQUENCE_CHAR_GAP);

  while (cp && *cp)
    {
      /* Found the start of a gap; remember the start coord */
      const int startCoord = cp - refSeq + refSeqRange->min;
    
      /* Loop until we find a non-gap character (or the end of the string) */
      while (cp && *cp == SEQUENCE_CHAR_GAP)
	++cp;
    
      const int endCoord = cp - refSeq - 1 + refSeqRange->min;
    
      MSP *msp = createEmptyMsp(&lastMsp, mspList);
      msp->type = BLXMSP_GAP;
      intrangeSetValues(&msp->qRange, startCoord, endCoord);
      
      featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);
    
      /* Continue looking for more gaps */
      cp = strchr(cp, SEQUENCE_CHAR_GAP);
    }
}


/* blxview() can be called either from other functions in the Blixem
 * program itself or directly by functions in other programs such as
 * xace.
 *
 * Interface
 *    pfetch:  if non-NULL, then we use pfetch instead of efetch for
 *             _all_ sequence fetching (use the node/port info. in the
 *             pfetch struct to locate the pfetch server).
 *
 */
gboolean blxview(CommandLineOptions *options,
                 GArray* featureLists[],
                 GList *seqList,
                 GSList *supportedTypes,
                 PfetchParams *pfetch, 
                 char *align_types, 
                 gboolean External,
                 GSList *styles,
                 GHashTable *lookupTable)
{
  if (blixemWindow)
    gtk_widget_destroy(blixemWindow) ;

  if (!External)
    {
      /* When being called internally, the config file will not have been
       * initialised yet, so do it now. */
      GError *error = NULL;
      blxInitConfig(NULL, options, &error);
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
    }
  
  validateInput(options);
  
  /* Find any assembly gaps (i.e. gaps in the reference sequence) */
  findAssemblyGaps(options->refSeq, featureLists, &options->mspList, &options->refSeqRange);
  
  gboolean status = bulkFetchSequences(
    0, External, options->saveTempFiles, options->seqType, &seqList, options->columnList,
    options->bulkFetchDefault, options->fetchMethods, &options->mspList, &options->blastMode, 
    featureLists, supportedTypes, NULL, 0, &options->refSeqRange, 
    options->dataset, FALSE, lookupTable); /* offset has not been applied yet, so pass offset=0 */

  if (status)
    {
      finaliseFetch(seqList, options->columnList);

      /* Construct missing data and do any other required processing now we have all the sequence data */
      finaliseBlxSequences(featureLists, &options->mspList, &seqList, options->columnList, 
                           options->refSeqOffset, options->seqType, 
                           options->numFrames, &options->refSeqRange, TRUE, lookupTable);

    }

  /* Note that we create a blxview even if MSPlist is empty.
   * But only if it's an internal call.  If external & anything's wrong, we die. */
  if (status || !External)
    {
      blviewCreate(align_types, padseq, featureLists, seqList, supportedTypes, options, External, styles) ;
    }

  return status;
}


/* Find all of the different alignment types (i.e. source names) and 
 * compile them into a single string, separated by the given separator string */
static char* findAlignTypes(GArray* featureLists[], const char *separatorStr)
{
  GSList *sourceList = NULL;
  GString *resultStr = g_string_new(NULL);
  
  /* Loop through all MSPs of type 'match' */
  const GArray* const mspList = featureLists[BLXMSP_MATCH];
  
  int i = 0;
  const MSP *msp = mspArrayIdx(mspList, i);
  
  for ( ; msp; msp = mspArrayIdx(mspList, ++i))
    {
      const char *source = mspGetSource(msp);
      
      if (source)
        {
          /* See if we've already seen this source. If not, add it to the list,
           * and append it to the string. */
          GQuark quark = g_quark_from_string(source);
          
          if (!g_slist_find(sourceList, GINT_TO_POINTER(quark)))
            {
              sourceList = g_slist_prepend(sourceList, GINT_TO_POINTER(quark));
              
              if (resultStr->len)
                g_string_append(resultStr, separatorStr);
              
              g_string_append(resultStr, source);
            }
        }
    }
  
  char *result = resultStr->str;
  g_string_free(resultStr, FALSE);
  
  return result;
}


/* Initialize the display and the buttons */
static void blviewCreate(char *align_types, 
			 const char *paddingSeq,
                         GArray* featureLists[],
                         GList *seqList,
                         GSList *supportedTypes,
			 CommandLineOptions *options,
                         const gboolean External,
                         GSList *styles)
{
  if (!blixemWindow)
    {
      /* Create the window */
      blixemWindow = createBlxWindow(options, paddingSeq, featureLists, seqList, supportedTypes, External, styles);

      /* Set the window title. Get a description of all the alignment types 
       * (unless already supplied) */
      if (!align_types)
        align_types = findAlignTypes(featureLists, ", ");

      /* If no alignment description was set, create a generic description */
      if (!align_types)
        align_types = g_strdup_printf("%s", options->seqType == BLXSEQ_PEPTIDE ? "peptide alignment" : "nucleotide alignment");
      
      BlxViewContext *bc = blxWindowGetContext(blixemWindow);
      
      char *title = g_strdup_printf("%s(%s) %s %s",
                                    blxGetTitlePrefix(bc),
                                    align_types,
                                    (options->dataset ? options->dataset : ""),
                                    options->refSeqName);
      
      gtk_window_set_title(GTK_WINDOW(blixemWindow), title);
      g_free(title);
    }

  char *nameSeparatorPos = (char *)strrchr(options->refSeqName, '/');
  if (nameSeparatorPos)
    {
      options->refSeqName = nameSeparatorPos + 1;
    }
  
  if (options->dotterFirst && options->mspList && options->mspList->sname &&
      (mspIsBlastMatch(options->mspList) || options->mspList->type == BLXMSP_HSP || options->mspList->type == BLXMSP_GSP))
    {
      /* We must select the sequence before calling callDotter. Get the first
       * sequence to the right of the start coord. */
      MSP *msp = nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
      
      if (msp)
        {
          blxWindowSelectSeq(blixemWindow, msp->sSequence);
          
          GError *error = NULL;
          char *dotterSSeq = getDotterSSeq(blixemWindow, &error);
          
          if (!error)
            {
              callDotter(blixemWindow, FALSE, dotterSSeq, NULL);
            }
            
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
      else 
        {
          g_error("Could not run Dotter on first sequence; no sequences found to the right of the start coord.\n");
        }
    }

  if (options->startNextMatch)
    {
      /* Set the start coord to be the start of the next MSP on from the default start coord */
      nextMatch(blxWindowGetDetailView(blixemWindow), NULL);
    }
}

/***********************************************************
 *               Sequences and MSPs
 ***********************************************************/

/* The parsed gene name generally contains extra info that we're not interested in. This function
 * tries to extract just the part of the info that we're interested in. */
static void processGeneName(BlxSequence *blxSeq)
{
  /* Data is usually in the format: "Name=SMARCA2; Synonyms=BAF190B, BRM, SNF2A, SNF2L2;"
   * Therefore, look for the Name tag and extract everything up to the ; or end of line.
   * If there is no name tag, just include everything */
  const char *geneName = blxSequenceGetGeneName(blxSeq);
  
  if (geneName)
    {
      char *startPtr = strstr(geneName, "Name=");
      
      if (!startPtr)
        {
          startPtr = strstr(geneName, "name=");
        }
        
      if (startPtr)
        {
          startPtr += 5;
          char *endPtr = strchr(startPtr, ';');
          const int numChars = endPtr ? endPtr - startPtr : strlen(startPtr);
          
          char *result = (char*)g_malloc((numChars + 1) * sizeof(char));
          
          g_utf8_strncpy(result, startPtr, numChars);
          result[numChars] = 0;
          
          blxSequenceSetValueFromString(blxSeq, BLXCOL_GENE_NAME, result);

          g_free(result);
        }
    }
}


/* Add a prefix to the given BlxSequence's organism field, e.g. change "Homo sapiens (human)" 
 * to "Hs - Homo sapiens (human)". This is so that we can have a very narrow column that just 
 * displays the prefix part of the field. */
static void processOrganism(BlxSequence *blxSeq)
{
  const char *organism = blxSequenceGetOrganism(blxSeq);
  
  if (organism)
    {
      /* To create the alias, we'll take the first letter of each word until we hit a non-alpha
       * character, e.g. the alias for "Homo sapiens (human)" would be "Hs" */
      int srcIdx = 0;
      int srcLen = strlen(organism);
      
      int destIdx = 0;
      int destLen = 5; /* limit the length of the alias */
      char alias[destLen + 1];
      
      gboolean startWord = TRUE;
      
      for ( ; srcIdx < srcLen && destIdx < destLen; ++srcIdx)
        {
          if (organism[srcIdx] == ' ')
            {
              startWord = TRUE;
            }
          else if (g_ascii_isalpha(organism[srcIdx]))
            {
              if (startWord)
                {
                  alias[destIdx] = organism[srcIdx];
                  ++destIdx;
                  startWord = FALSE;
                }
            }
          else
            {
              /* Finish if we've hit a char that's not an alpha char or space. */
              break;
            }
        }

      alias[destIdx] = '\0';
      
      if (destIdx > 0)
        {
          blxSeq->organismAbbrev = g_strdup(alias);
        }
    }
}


/* Returns true if this msp is part of a feature series */
gboolean mspHasFs(const MSP *msp)
{
  gboolean result = (msp->type == BLXMSP_FS_SEG || msp->type == BLXMSP_XY_PLOT);
  return result;
}


/* Return the (cached) full extent of the match that we're showing in match seq coords */
const IntRange* mspGetFullSRange(const MSP* const msp, const gboolean seqSelected, const BlxViewContext* const bc)
{
  const IntRange *result = NULL;
  
  if (seqSelected || !bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED] || !bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED])
    result = &msp->fullSRange;
  else
    result = &msp->sRange;
  
  return result;
}

/* Return the (cached) full extent of the match on the ref sequence in display coords
 * (including any portions of unaligned sequence that we're showing). Depending on the
 * options, the range may depend on whether the sequence is selected or not. */
const IntRange* mspGetFullDisplayRange(const MSP* const msp, const gboolean seqSelected, const BlxViewContext* const bc)
{
  const IntRange *result = NULL;
  
  /* Check if showing unaligned sequence or polya tails for all sequences, or just the selected
     sequence */
  if ((bc->flags[BLXFLAG_SHOW_UNALIGNED] || bc->flags[BLXFLAG_SHOW_POLYA_SITE]) && 
      (seqSelected || 
       (bc->flags[BLXFLAG_SHOW_UNALIGNED] && !bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED]) || 
       (bc->flags[BLXFLAG_SHOW_POLYA_SITE] && !bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED])))
    {
      result = &msp->fullRange;
    }
  else
    {
      result = &msp->displayRange;
    }
  
  return result;
}

/* Return the (cached) extent of the alignment that we're showing in display coords
 * (excluding unaligned sequence) */
const IntRange* mspGetDisplayRange(const MSP* const msp)
{
  return &msp->displayRange;
}


/* Get the full range of the given MSP that we want to display, in s coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
static void mspCalcFullSRange(const MSP* const msp, 
			      const gboolean *flags,
			      const int numUnalignedBases, 
			      const GArray* const polyASiteList,
			      IntRange *result)
{
  /* Normally we just display the part of the sequence in the alignment */
  result->min = msp->sRange.min;
  result->max = msp->sRange.max;
  
  if (mspIsBlastMatch(msp))
    {
      if (flags[BLXFLAG_SHOW_UNALIGNED] && mspGetMatchSeq(msp))
	{
	  /* We're displaying additional unaligned sequence outside the alignment range. Get 
	   * the full range of the match sequence */
	  result->min = 1;
	  result->max = mspGetMatchSeqLen(msp);
	  
	  if (flags[BLXFLAG_LIMIT_UNALIGNED_BASES])
	    {
	      /* Only include up to 'numUnalignedBases' each side of the MSP range (still limited
	       * to the range we found above, though). */
	      result->min = max(result->min, msp->sRange.min - numUnalignedBases);
	      result->max = min(result->max, msp->sRange.max + numUnalignedBases);
	    }
	}
      
      if (flags[BLXFLAG_SHOW_POLYA_SITE] && mspHasPolyATail(msp))
	{
	  /* We're displaying polyA tails, so override the 3' end coord with the full extent of
	   * the s sequence if there is a polyA site here. The 3' end is the min q coord if the
	   * match is on forward ref seq strand or the max coord if on the reverse. */
	  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
	  
	  if (sameDirection)
	    {
	      result->max = mspGetMatchSeqLen(msp);
	    }
	  else
	    {
	      result->min = 1;
	    }
	}
    }
}


/* Get the full range of the given MSP that we want to display, in q coords. This will generally 
 * be the coords of the alignment but could extend outside this range we are displaying unaligned 
 * portions of the match sequence or polyA tails etc. */
static void mspCalcFullQRange(const MSP* const msp, 
			      const gboolean *flags,
			      const int numUnalignedBases, 
			      const GArray* const polyASiteList,
			      const int numFrames,
			      const IntRange* const fullSRange,
			      IntRange *result)
{
  /* Default to the alignment range so we can exit quickly if there are no special cases */
  result->min = msp->qRange.min;
  result->max = msp->qRange.max;
  
  if (mspIsBlastMatch(msp) && (flags[BLXFLAG_SHOW_UNALIGNED] || flags[BLXFLAG_SHOW_POLYA_SITE]))
    {
      /* Find the offset of the start and end of the full range compared to the alignment range and
       * offset the ref seq range by the same amount. We need to multiply by the number of reading
       * frames because q coords are in nucleotides and s coords are in peptides. */
      const int startOffset = (msp->sRange.min - fullSRange->min) * numFrames;
      const int endOffset = (fullSRange->max - msp->sRange.max) * numFrames;
      
      const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
      result->min -= sameDirection ? startOffset : endOffset;
      result->max += sameDirection ? endOffset : startOffset;
    }
}


/* Calculate the full extent of the match sequence to display, in display coords,
 * and cache the result in the msp. Includes any portions of unaligned sequence that we're
 * displaying */
void mspCalculateFullExtents(MSP *msp, const BlxViewContext* const bc, const int numUnalignedBases)
{
  mspCalcFullSRange(msp, bc->flags, numUnalignedBases, bc->featureLists[BLXMSP_POLYA_SITE], &msp->fullSRange);
  mspCalcFullQRange(msp, bc->flags, numUnalignedBases, bc->featureLists[BLXMSP_POLYA_SITE], bc->numFrames, &msp->fullSRange, &msp->fullRange);
 
  /* convert the Q range to display coords */
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int coord1 = convertDnaIdxToDisplayIdx(msp->fullRange.min, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->fullRange.max, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  intrangeSetValues(&msp->fullRange, coord1, coord2);
  
  /* Remember the max len of all the MSPs in the detail-view */
  if (typeShownInDetailView(msp->type) && getRangeLength(&msp->fullRange) > getMaxMspLen())
    {
      setMaxMspLen(getRangeLength(&msp->fullRange));
    }
}


/* Convert the ref-seq range of the given msp in display coords and cache it in the msp */
static void mspCalculateDisplayRange(MSP *msp, const BlxViewContext* const bc)
{
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  intrangeSetValues(&msp->displayRange, coord1, coord2);  
}


/* This caches the display range (in display coords rather than dna coords,
 * and inverted if the display is inverted) for each MSP */
void cacheMspDisplayRanges(const BlxViewContext* const bc, const int numUnalignedBases)
{
  /* This also calculates the max msp len */
  setMaxMspLen(0);
  
  MSP *msp = bc->mspList;
  for ( ; msp; msp = msp->next)
    {
      mspCalculateDisplayRange(msp, bc);
      mspCalculateFullExtents(msp, bc, numUnalignedBases);
    }
}


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the MSP is a gapped MSP and the ref-seq coord is known to lie within the
 * MSP's alignment range. */
static int mspGetGappedAlignmentCoord(const MSP *msp, const int qIdx, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));

  /* Gapped alignment. Look to see if x lies inside one of the "gaps" ranges. */
  GSList *rangeItem = msp->gaps;
  
  for ( ; rangeItem ; rangeItem = rangeItem->next)
    {
      CoordRange *curRange = (CoordRange*)(rangeItem->data);
      
      int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
      getCoordRangeExtents(curRange, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
      
      /* We've "found" the value if it's in or before this range. Note that the
       * the range values are in decreasing order if the q strand is reversed. */
      gboolean found = qForward ? qIdx <= qRangeMax : qIdx >= qRangeMin;
      
      if (found)
        {
          gboolean inRange = (qForward ? qIdx >= qRangeMin : qIdx <= qRangeMax);
          
          if (inRange)
            {
              /* It's inside this range. Calculate the actual index. */
              int offset = (qIdx - qRangeMin) / bc->numFrames;
              result = sameDirection ? sRangeMin + offset : sRangeMax - offset;
            }
          
          break;
        }
    }
  
  return result;
}


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the MSP is an ungapped MSP and the ref-seq coord is known to lie within the
 * MSP's alignment range. */
static int mspGetUngappedAlignmentCoord(const MSP *msp, const int qIdx, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  /* If strands are in the same direction, find the offset from qRange.min and add it to 
   * sRange.min. If strands are in opposite directions, find the offset from qRange.min and
   * subtract it from sRange.max. Note that the offset could be negative if we're outside
   * the alignment range. */
  int offset = (qIdx - msp->qRange.min) / bc->numFrames ;
  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));

  result = (sameDirection) ? msp->sRange.min + offset : msp->sRange.max - offset ;
  
  if (result < msp->sRange.min || result > msp->sRange.max)
    {
      result = UNSET_INT;
    }
  
  return result;
}  


/* Return the match-sequence coord of an MSP at the given reference-sequence coord,
 * where the ref-seq coord is known to lie outside the MSP's alignment range. The 
 * result will be unset unless the option to display unaligned portions of 
 * sequence is enabled. */
static int mspGetUnalignedCoord(const MSP *msp, const int qIdx, const gboolean seqSelected, const int numUnalignedBases, const BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  /* First convert to display coords */
  const int frame = mspGetRefFrame(msp, bc->seqType);
  const int displayIdx = convertDnaIdxToDisplayIdx(qIdx, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, NULL);
  const IntRange* const mspRange = mspGetDisplayRange(msp);

  /* Note that because we have converted to display coords the ref seq coords are
   * now in the direction of the display, regardless of the ref seq strand. */
  const gboolean sameDirection = (mspGetMatchStrand(msp) == mspGetRefStrand(msp));
  const gboolean sForward = (sameDirection != bc->displayRev);
  
  if (displayIdx < mspRange->min)
    {
      /* Find the offset backwards from the low coord and subtract it from the low end
       * of the s coord range (or add it to the high end, if the directions are opposite).
       * We're working in display coords here. */
      const int offset = mspRange->min - displayIdx;
      result = sForward ? msp->sRange.min - offset : msp->sRange.max + offset;
    }
  else
    {
      /* Find the offset forwards from the high coord and add it to the high end of the
       * s coord range (or subtract it from the low end, if the directions are opposite). */
      const int offset = displayIdx - mspRange->max;
      result = sForward ? msp->sRange.max + offset : msp->sRange.min - offset;
    }
  
  /* Get the full display range of the match sequence. If the result is still out of range
   * then there's nothing to show for this qIdx. */
  const IntRange *fullSRange = mspGetFullSRange(msp, seqSelected, bc);
  
  if (!valueWithinRange(result, fullSRange))
    {
      result = UNSET_INT;
    }
  
  return result;
}



/* Given a base index on the reference sequence, find the corresonding base 
 * in the match sequence. The return value is always UNSET_INT if there is not
 * a corresponding base at this position. */
int mspGetMatchCoord(const MSP *msp, 
                     const int qIdx, 
                     const gboolean seqSelected,
                     const int numUnalignedBases,
                     BlxViewContext *bc)
{
  int result = UNSET_INT;
  
  if (mspIsBlastMatch(msp) || mspIsExon(msp))
    {
      const gboolean inMspRange = valueWithinRange(qIdx, &msp->qRange);
      
      if (msp->gaps && g_slist_length(msp->gaps) >= 1 && inMspRange)
        {
          result = mspGetGappedAlignmentCoord(msp, qIdx, bc);
        }
      else if (!inMspRange && mspIsBlastMatch(msp))
        {
          /* The q index is outside the alignment range but if the option to show
           * unaligned sequence is enabled we may still have a valie result. */
          result = mspGetUnalignedCoord(msp, qIdx, seqSelected, numUnalignedBases, bc);
        }
      else
        {
          result = mspGetUngappedAlignmentCoord(msp, qIdx, bc);
        }
    }
  
  return result;
}



/***********************************************************
 *               General
 ***********************************************************/

/* Redraw the entire blixem window. (Call blxWindowRedrawAll directly where possible,
 * rather than this function, which relies on the global variable 'blixemWindow'). */
void blviewRedraw(void)
{
  if (blixemWindow)
    {
      blxWindowRedrawAll(blixemWindow);
    }
}

GtkWidget* getBlixemWindow()
{
  return blixemWindow;
}


/* Reset any global variables / singleton instances */
void blviewResetGlobals()
{
  blixemWindow = NULL;
  destroyMessageList();
}

/***********************************************************
 *               Styles and colours
 ***********************************************************/

/* Set the given BlxColor elements from the given color string(s). Works out some good defaults
 * for blxColor.selected blxcolor.print etc if these colors are not given */
static void setBlxColorValues(const char *normal, const char *selected, BlxColor *blxColor, GError **error)
{
  GError *tmpError = NULL;
  
  if (normal)
    {
      getColorFromString(normal, &blxColor->normal, &tmpError);

      if (!tmpError)
        {
          getSelectionColor(&blxColor->normal, &blxColor->selected); /* will use this if selection color is not given/has error */
          
          /* Calculate print coords as a greyscale version of normal colors */
          convertToGrayscale(&blxColor->normal, &blxColor->print);            /* calculate print colors, because they are not given */
          getSelectionColor(&blxColor->print, &blxColor->printSelected);
        }

      /* Treat white as transparent. Not ideal but blixem can read zmap styles
       * files where this assumption is made */
      blxColor->transparent = (!strncasecmp(normal, "white", 5) || !strncasecmp(normal, "#ffffff", 7));
    }
  else
    {
      blxColor->transparent = TRUE;
    }

  /* Use the selected-color string, if given */
  if (!tmpError && selected)
    {
      getColorFromString(selected, &blxColor->selected, &tmpError);
    }
  
  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Create a BlxStyle. For transcripts, CDS and UTR features can have different colors, but they
 * are from the same source and hence have the same style object. We therefore use the default
 * fillColor, lineColor etc. variables for CDS features and provide additional options to supply 
 * UTR colors as well. */
BlxStyle* createBlxStyle(const char *styleName, 
			 const char *fillColor, 
			 const char *fillColorSelected, 
			 const char *lineColor, 
			 const char *lineColorSelected, 
			 const char *fillColorUtr, 
			 const char *fillColorUtrSelected, 
			 const char *lineColorUtr, 
			 const char *lineColorUtrSelected, 
			 GError **error)
{
  BlxStyle *style = (BlxStyle*)g_malloc(sizeof(BlxStyle));
  GError *tmpError = NULL;

  if (!styleName)
    {
      g_set_error(error, BLX_ERROR, 1, "Style name is NULL.\n");
    }
    
  if (!tmpError)
    {
      style->styleName = g_strdup(styleName);
    }
  
  if (!tmpError)
    {      
      setBlxColorValues(fillColor ? fillColor : fillColorUtr,
                        fillColorSelected ? fillColorSelected : fillColorUtrSelected,
                        &style->fillColor, &tmpError);
    }
  
  if (!tmpError)
    {
      setBlxColorValues(lineColor ? lineColor : lineColorUtr,
                        lineColorSelected ? lineColorSelected : lineColorUtrSelected,
                        &style->lineColor, &tmpError);
    }
  
  if (!tmpError)
    {
      setBlxColorValues(fillColorUtr ? fillColorUtr : fillColor,
                        fillColorUtrSelected ? fillColorUtrSelected : fillColorSelected,
                        &style->fillColorUtr, &tmpError);
    }

  if (!tmpError)
    {
      setBlxColorValues(lineColorUtr ? lineColorUtr : lineColor,
                        lineColorUtrSelected ? lineColorUtrSelected : lineColorSelected,
                        &style->lineColorUtr, &tmpError);
    }

  if (tmpError)
    {
      g_free(style);
      style = NULL;
      g_propagate_error(error, tmpError);
    }
  
  return style;
}


/* Destroy a BlxStyle */
void destroyBlxStyle(BlxStyle *style)
{
  if (style)
    {
      g_free(style->styleName);
    }
}


/* This function highlights any regions within the given display range that
 * are assembly gaps (which are BLXMSP_GAP type MSPs in the given MSP array).
 * The given rect defines the drawing area on the given drawable that corresponds
 * to the given display range. Gaps that lie within the display range are drawn
 * within the rect at the appropriate positions. */
void drawAssemblyGaps(GtkWidget *widget,
                      GdkDrawable *drawable,
                      GdkColor *color,
                      const gboolean displayRev,
                      GdkRectangle *rect, 
                      const IntRange* const dnaRange,
                      const GArray *mspArray)
{
  /* See if any gaps lie within the display range. */
  int i = 0;
  MSP *gap = mspArrayIdx(mspArray, i);
  
  for ( ; gap; gap = mspArrayIdx(mspArray, ++i))
    {
      if (rangesOverlap(&gap->qRange, dnaRange))
        {
	  /* Draw to max coord plus one (or min coord minus one if reversed) to be inclusive */
	  const int x1 = convertBaseIdxToRectPos(displayRev ? gap->qRange.min - 1 : gap->qRange.min, rect, dnaRange, TRUE, displayRev, TRUE);
	  const int x2 = convertBaseIdxToRectPos(gap->qRange.max + 1, rect, dnaRange, TRUE, displayRev, TRUE);
          
	  const int width = max(MIN_GAP_HIGHLIGHT_WIDTH, abs(x2 - x1));
	
          cairo_t *cr = gdk_cairo_create(drawable);
          gdk_cairo_set_source_color(cr, color);
          cairo_rectangle(cr, min(x1, x2), 0, width, widget->allocation.height);
          cairo_clip(cr);
          cairo_paint_with_alpha(cr, 0.1);
          cairo_destroy(cr);
        }
    }
  
}


/* Get the color strings for the given group from the given 
 * styles file. If the string is not found then recurse
 * through any parent styles */
static void getStylesFileColorsRecursive(GKeyFile *keyFile,
                                         const char *group,
                                         BlxStyleColors *colors,
                                         GError **error)
{
  /* Loop through the normal colors, if we still need to find any */
  if (!colors->fillColor || !colors->lineColor ||
      !colors->fillColorSelected || !colors->lineColorSelected)
    {
      char **normalColors = g_key_file_get_string_list(keyFile, group, "colours", NULL, NULL);
      char **color = normalColors;
      
      if (normalColors)
        colors->normalFound = TRUE;
      
      for ( ; color && *color; ++color)
        {
          /* Ignore leading whitespace */
          char *c = *color;
          while (*c == ' ')
            ++c;
          
          if (c && *c)
            {
              if (!strncasecmp(c, "normal fill ", 12))
                {
                  if (!colors->fillColor) 
                    colors->fillColor = g_strchug(g_strchomp(g_strdup(c + 12)));
                }
              else if (!strncasecmp(c, "normal border ", 14))
                {
                  if (!colors->lineColor)
                    colors->lineColor = g_strchug(g_strchomp(g_strdup(c + 14)));
                }
              else if (!strncasecmp(c, "selected fill ", 14))
                {
                  if (!colors->fillColorSelected)
                    colors->fillColorSelected = g_strchug(g_strchomp(g_strdup(c + 14)));
                }
              else if (!strncasecmp(c, "selected border ", 16))
                {
                  if (!colors->lineColorSelected)
                    colors->lineColorSelected = g_strchug(g_strchomp(g_strdup(c + 16)));
                }
              else if (error && *error)
                {
                  postfixError(*error, "  '%s' is not a valid color field\n", c);
                }
              else
                {
                  g_set_error(error, BLX_ERROR, 1, "  '%s' is not a valid color field\n", c);
                }
            }
        }

      if (normalColors)
        g_strfreev(normalColors);
    }

  /* Loop through the CDs colors, if we still need to find any */
  if (!colors->fillColorCds || !colors->lineColorCds ||
      !colors->fillColorCdsSelected || !colors->lineColorCdsSelected)
    {
      char **cdsColors = g_key_file_get_string_list(keyFile, group, "transcript-cds-colours", NULL, NULL);
      char **color = cdsColors;

      if (cdsColors)
        colors->cdsFound = TRUE;

      for ( ; color && *color; ++color)
        {
          /* Ignore leading whitespace */
          char *c = *color;
          while (*c == ' ')
            ++c;
          
          if (c && *c)
            {
              if (!strncasecmp(c, "normal fill ", 12))
                {
                  if (!colors->fillColorCds) 
                    colors->fillColorCds = g_strchug(g_strchomp(g_strdup(c + 12)));
                }
              else if (!strncasecmp(c, "normal border ", 14))
                {
                  if (!colors->lineColorCds)
                    colors->lineColorCds = g_strchug(g_strchomp(g_strdup(c + 14)));
                }
              else if (!strncasecmp(c, "selected fill ", 14))
                {
                  if (!colors->fillColorCdsSelected)
                    colors->fillColorCdsSelected = g_strchug(g_strchomp(g_strdup(c + 14)));
                }
              else if (!strncasecmp(c, "selected border ", 16))
                {
                  if (!colors->lineColorCdsSelected)
                    colors->lineColorCdsSelected = g_strchug(g_strchomp(g_strdup(c + 16)));
                }
              else if (error && *error)
                {
                  postfixError(*error, "  '%s' is not a valid color field\n", c);
                }
              else
                {
                  g_set_error(error, BLX_ERROR, 1, "  '%s' is not a valid color field\n", c);
                }
            }
        }
      
      if (cdsColors)
        g_strfreev(cdsColors);
    }

  /* If there are still any outstanding, recurse to the next parent */
  if (!colors->fillColor || !colors->lineColor || 
      !colors->fillColorSelected || !colors->lineColorSelected || 
      !colors->fillColorCds || !colors->lineColorCds ||
      !colors->fillColorCdsSelected || !colors->lineColorCdsSelected)
    {
      char *parent = g_key_file_get_string(keyFile, group, "parent-style", NULL);
      
      if (parent)
        {
          getStylesFileColorsRecursive(keyFile, parent, colors, error);
          g_free(parent);
        }
    }
}


/* For backwards compatibility, read the old-style color fields
 * for the given source (group) from the given key file. (The 
 * key names were changed to be consistent with zmap). */
static void readStylesFileColorsOld(GKeyFile *keyFile, 
                                    const char *group, 
                                    GSList **stylesList,
                                    GError **error)
{
  char *fillColor = g_key_file_get_string(keyFile, group, "fill_color", NULL);
  char *lineColor = g_key_file_get_string(keyFile, group, "line_color", NULL);
  char *fillColorSelected = g_key_file_get_string(keyFile, group, "fill_color_selected", NULL);
  char *lineColorSelected = g_key_file_get_string(keyFile, group, "line_color_selected", NULL);
  char *fillColorUtr = g_key_file_get_string(keyFile, group, "fill_color_utr", NULL);
  char *lineColorUtr = g_key_file_get_string(keyFile, group, "line_color_utr", NULL);
  char *fillColorUtrSelected = g_key_file_get_string(keyFile, group, "fill_color_utr_selected", NULL);
  char *lineColorUtrSelected = g_key_file_get_string(keyFile, group, "line_color_utr_selected", NULL);
  
  /* If there was an error, skip this group. Otherwise, go ahead and create the style */
  if (fillColor && lineColor)
    {
      BlxStyle *style = createBlxStyle(group, 
                                       fillColor, fillColorSelected,
                                       lineColor, lineColorSelected,
                                       fillColorUtr, fillColorUtrSelected,
                                       lineColorUtr, lineColorUtrSelected,
                                       error);

      if (style)
        *stylesList = g_slist_append(*stylesList, style);
    }
  else if (!fillColor)
    {
      g_set_error(error, BLX_ERROR, 1, "Style '%s' does not contain required field fill_color", group);
    }
  else if (!lineColor)
    {
      g_set_error(error, BLX_ERROR, 1, "Style '%s' does not contain required field line_color", group);
    }
  
  
  if (fillColor) g_free(fillColor);
  if (lineColor) g_free(lineColor);
  if (fillColorSelected) g_free(fillColorSelected);
  if (lineColorSelected) g_free(lineColorSelected);
  if (fillColorUtr) g_free(fillColorUtr);
  if (lineColorUtr) g_free(lineColorUtr);
  if (fillColorUtrSelected) g_free(fillColorUtrSelected);
  if (lineColorUtrSelected) g_free(lineColorUtrSelected);
}


/* Read the new-style color fields for the given source (group) 
 * from the given key file. Returns true if colors were found. */
static gboolean readStylesFileColors(GKeyFile *keyFile, 
                                     const char *group, 
                                     GSList **stylesList,
                                     GError **error)
{
  gboolean result = FALSE;
  GError *tmpError = NULL;
  g_key_file_set_list_separator(keyFile, ';');

  /* Get the colors from the styles. We may need to recurse through multiple
   * style parents before we find them. */
  BlxStyleColors colors = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, FALSE, FALSE};
  getStylesFileColorsRecursive(keyFile, group, &colors, &tmpError);
  
  /* If colors were found, create the style */
  if (!tmpError && (colors.normalFound || colors.cdsFound))
    {
      /* Unfortunately blixem stores cds/utr colors differently to zmap, which
       * leads to the following messy logic:
       * If the optional cds colors are set then the default 'colors' field 
       * in the styles file gives our utr colors and the 'transcript-cds-colors'
       * field gives our normal colors. Otherwise, the 'colors' field gives our
       * normal colors and we don't set the utr colors (because they default 
       * to the normal colors). */
      BlxStyle *style = NULL;
      
      if (colors.cdsFound)
        {
          style = createBlxStyle(group, 
                                 colors.fillColorCds, colors.fillColorCdsSelected,
                                 colors.lineColorCds, colors.lineColorCdsSelected,
                                 colors.fillColor, colors.fillColorSelected,
                                 colors.lineColor, colors.lineColorSelected,
                                 &tmpError);
        }
      else
        {
          style = createBlxStyle(group, 
                                 colors.fillColor, colors.fillColorSelected,
                                 colors.lineColor, colors.lineColorSelected,
                                 NULL, NULL,
                                 NULL, NULL,
                                 &tmpError);
        }

      if (style)
        *stylesList = g_slist_append(*stylesList, style);

      if (tmpError)
        prefixError(tmpError, "  "); /* indent the error message */
    }
  
  /* Clean up */
  if (colors.fillColor) g_free(colors.fillColor);
  if (colors.lineColor) g_free(colors.lineColor);
  if (colors.fillColorSelected) g_free(colors.fillColorSelected);
  if (colors.lineColorSelected) g_free(colors.lineColorSelected);
  if (colors.fillColorCds) g_free(colors.fillColorCds);
  if (colors.lineColorCds) g_free(colors.lineColorCds);
  if (colors.fillColorCdsSelected) g_free(colors.fillColorCdsSelected);
  if (colors.lineColorCdsSelected) g_free(colors.lineColorCdsSelected);

  /* Error handling */
  if (tmpError)
    g_propagate_error(error, tmpError);

  return result;
}


/* Read in the styles file. Returns a list of style structs for each style
 * found. Uses the given key file if specified, otherwise checks to see
 * if there is a key file specified in the config file. */
GSList* blxReadStylesFile(const char *keyFileName_in, GError **error)
{
  GSList *result = NULL;

  char *keyFileName = NULL;

  if (keyFileName_in)
    {
      keyFileName = g_strdup(keyFileName_in);
    }
  else
    {
      /* See if the styles file is specified in the config */
      GKeyFile *keyFile = blxGetConfig();

      if (keyFile)
        keyFileName = g_key_file_get_string(keyFile, BLIXEM_GROUP, STYLES_FILE_KEY, NULL);
    }
  
  if (!keyFileName)
    {
      return result;
    }
  
  /* Load the key file */
  GKeyFile *keyFile = g_key_file_new();
  GKeyFileFlags flags = G_KEY_FILE_NONE ;

  if (g_key_file_load_from_file(keyFile, keyFileName, flags, error))
    {
      /* Get all the groups (i.e. style names) from the file */
      gsize num_groups;
      char **groups = g_key_file_get_groups(keyFile, &num_groups) ;

      /* Loop through each style */
      char **group;
      int i;
      GError *tmpError = NULL;
      
      for (i = 0, group = groups ; i < (int)num_groups ; i++, group++)
	{
          gboolean found = readStylesFileColors(keyFile, *group, &result, &tmpError);

          /* If it wasn't found, look for old-style colors */
          if (!found && !tmpError)
            readStylesFileColorsOld(keyFile, *group, &result, NULL);

          if (tmpError)
            {
              /* Compile all errors into one */
              prefixError(tmpError, "[%s]\n", *group);
              
              if (error && *error)
                postfixError(*error, "%s", tmpError->message);
              else
                g_set_error(error, BLX_ERROR, 1, "%s", tmpError->message);

              g_error_free(tmpError);
              tmpError = NULL;
            }
        }
      
      if (tmpError)
        {          
          g_propagate_error(error, tmpError);
        }
      
      g_strfreev(groups);
    }
  
  if (error && *error)
    {
      prefixError(*error, "Errors found while reading styles file '%s'\n", keyFileName);
      postfixError(*error, "\n");
    }

  g_free(keyFileName);
  g_key_file_free(keyFile) ;
  keyFile = NULL ;
  
  return result;
}


/***********************************************************
 *                        Utilities
 ***********************************************************/

/* Returns a string which is the name of the Blixem application. */
const char *blxGetAppName()
{
  return BLIXEM_TITLE ;
}

/* Returns a string which is used to prefix window titles (full or abbrev
 * depending on settings). */
const char *blxGetTitlePrefix(const BlxViewContext * const bc)
{
  return bc->flags[BLXFLAG_ABBREV_TITLE] ? BLIXEM_PREFIX_ABBREV : BLIXEM_PREFIX ;
}

/* Returns a copyright string for the Blixem application. */
const char *blxGetCopyrightString()
{
  return BLIXEM_COPYRIGHT_STRING ;
}

/* Returns the Blixem website URL. */
const char *blxGetWebSiteString()
{
  return BLIXEM_WEBSITE_STRING ;
}

/* Returns a comments string for the Blixem application. */
const char *blxGetCommentsString()
{
  return BLIXEM_COMMENTS_STRING() ;
}

/* Returns a license string for the blx application. */
const char *blxGetLicenseString()
{
  return BLIXEM_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Blixem code. */
const char *blxGetVersionString()
{
  return BLIXEM_VERSION_STRING ;
}


/***********************************************************
 *                      Columns
 ***********************************************************/

/* Check if the column has a column-width set in the config and if so override the default value */
static void getColumnWidthsConfig(BlxColumnInfo *columnInfo)
{
  /* Do nothing for the sequence column, because it has dynamic width */
  if (columnInfo->columnId == BLXCOL_SEQUENCE)
    return;
    
  GKeyFile *key_file = blxGetConfig();
  GError *error = NULL;

  if (key_file && g_key_file_has_group(key_file, COLUMN_WIDTHS_GROUP) && g_key_file_has_key(key_file, COLUMN_WIDTHS_GROUP, columnInfo->title, &error))
    {
      if (!error)
        {
          const int newWidth = g_key_file_get_integer(key_file, COLUMN_WIDTHS_GROUP, columnInfo->title, &error);
          
          if (error)
            {
              reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
            }
          else if (newWidth > 0)
            {
              columnInfo->width = newWidth;
              columnInfo->showColumn = TRUE;
            }
          else
            {
              columnInfo->showColumn = FALSE;
            }
        }
    }
}


/* Check if the column has a column-summary flag set in the config and if so 
 * override the default value */
static void getColumnSummaryConfig(BlxColumnInfo *columnInfo)
{
  /* Do nothing if we can't show the summary for this column */
  if (!columnInfo->canShowSummary)
    return;
    
  GKeyFile *key_file = blxGetConfig();

  if (key_file && g_key_file_has_group(key_file, COLUMN_SUMMARY_GROUP))
    {
      /* The value is true if we should include the column in summary info, false if we should
       * not. If the group exists but the column is not listed then that also indicates we should not
       * show it (this call will correctly return false in that case). */
      columnInfo->showSummary = g_key_file_get_boolean(key_file, COLUMN_SUMMARY_GROUP, columnInfo->title, NULL);
    }
}


/* This checks if the column has any properties specified for it in the config
 * file and, if so, overrides the default values in the given column with 
 * the config file values. */
static void getColumnConfig(BlxColumnInfo *columnInfo)
{
  getColumnWidthsConfig(columnInfo);
  getColumnSummaryConfig(columnInfo);
}

//static void destroyColumnList(GList **columnList)
//{
//  GList *column = *columnList;
//  
//  for ( ; column; column = column->next)
//    {
//      g_free(column->data);
//    }
//    
//  g_list_free(*columnList);
//  *columnList = NULL;
//}


/* This creates BlxColumnInfo entries for each column required in the detail view. It
 * returns a list of the columns created. */
GList* blxCreateColumns(const gboolean optionalColumns, const gboolean customSeqHeader)
{
  GList *columnList = NULL;
  
  /* Create the columns' data structs. The columns appear in the order
   * that they are added here. */
  blxColumnCreate(BLXCOL_SEQNAME,     TRUE,             "Name",       G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_SEQNAME_WIDTH,        TRUE,            TRUE,  TRUE,  TRUE,  TRUE,   "Name",        NULL, NULL, &columnList);
  blxColumnCreate(BLXCOL_SOURCE,      TRUE,             "Source",     G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_SOURCE_WIDTH,         TRUE,            TRUE,  TRUE,  TRUE,  TRUE,   "Source",      NULL, NULL, &columnList);
                                                                                                                                                                                                              
  blxColumnCreate(BLXCOL_ORGANISM,    TRUE,             "Organism",   G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_ORGANISM_WIDTH,       optionalColumns, TRUE,  TRUE,  TRUE,  TRUE,   "Organism",    "OS", NULL, &columnList);
  blxColumnCreate(BLXCOL_GENE_NAME,   TRUE,             "Gene Name",  G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_GENE_NAME_WIDTH,      optionalColumns, FALSE, TRUE,  TRUE,  TRUE,   "Gene name",   "GN", NULL, &columnList);
  blxColumnCreate(BLXCOL_TISSUE_TYPE, TRUE,             "Tissue Type",G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_TISSUE_TYPE_WIDTH,    optionalColumns, FALSE, TRUE,  TRUE,  TRUE,   "Tissue type", "FT", "tissue_type", &columnList);
  blxColumnCreate(BLXCOL_STRAIN,      TRUE,             "Strain",     G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_STRAIN_WIDTH,         optionalColumns, FALSE, TRUE,  TRUE,  TRUE,   "Strain",      "FT", "strain", &columnList);
                                                                                                                                                                             
  /* Insert optional columns here, with a dynamically-created IDs that are >= BLXCOL_NUM_COLS */                                                                             
  int columnId = BLXCOL_NUM_COLUMNS;                                                                                                                                         
  blxColumnCreate((BlxColumnId)columnId++, TRUE,        "Description",G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_DEFAULT_WIDTH,        optionalColumns, FALSE, TRUE,  TRUE,  TRUE,   "Description", "DE", NULL, &columnList);
                                                                                                                                                                                                             
  blxColumnCreate(BLXCOL_GROUP,       TRUE,             "Group",      G_TYPE_STRING, RENDERER_TEXT_PROPERTY,     BLXCOL_GROUP_WIDTH,          TRUE,            FALSE, FALSE, FALSE, TRUE,   "Group",       NULL, NULL, &columnList);
  blxColumnCreate(BLXCOL_SCORE,       TRUE,             "Score",      G_TYPE_DOUBLE, RENDERER_TEXT_PROPERTY,     BLXCOL_SCORE_WIDTH,          TRUE,            TRUE,  FALSE, FALSE, FALSE,  "Score",       NULL, NULL, &columnList);
  blxColumnCreate(BLXCOL_ID,          TRUE,             "%Id",        G_TYPE_DOUBLE, RENDERER_TEXT_PROPERTY,     BLXCOL_ID_WIDTH,             TRUE,            TRUE,  FALSE, FALSE, FALSE,  "Identity",    NULL, NULL, &columnList);
  blxColumnCreate(BLXCOL_START,       TRUE,             "Start",      G_TYPE_INT,    RENDERER_TEXT_PROPERTY,     BLXCOL_START_WIDTH,          TRUE,            TRUE,  FALSE, FALSE, FALSE,  "Position",    NULL, NULL, &columnList);
  blxColumnCreate(BLXCOL_SEQUENCE,    !customSeqHeader, "Sequence",   G_TYPE_POINTER,RENDERER_SEQUENCE_PROPERTY, BLXCOL_SEQUENCE_WIDTH,       TRUE,            TRUE,  FALSE, FALSE, FALSE,  NULL,          "SQ", NULL, &columnList);
  blxColumnCreate(BLXCOL_END,         TRUE,             "End",        G_TYPE_INT,    RENDERER_TEXT_PROPERTY,     BLXCOL_END_WIDTH,            TRUE,            TRUE,  FALSE, FALSE, FALSE,  NULL,          NULL, NULL, &columnList);

  /* For each column, check if it is configured in the config file and if so update accordingly */
  GList *columnItem  = columnList;

  for ( ; columnItem; columnItem = columnItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(columnItem->data);
      getColumnConfig(columnInfo);
    }

  return columnList;
}


/* Return the width of the column with the given column id */
int getColumnWidth(GList *columnList, const BlxColumnId columnId)
{
  int result = 0;

  BlxColumnInfo *columnInfo = getColumnInfo(columnList, columnId);

  if (columnInfo)
    {
      result = columnInfo->width;
    }
  
  return result;
}


/* Return the width of the column with the given column id */
const char* getColumnTitle(GList *columnList, const BlxColumnId columnId)
{
  const char *result = NULL;
  
  BlxColumnInfo *columnInfo = getColumnInfo(columnList, columnId);

  if (columnInfo)
    {
      result = columnInfo->title;
    }
  
  return result;
}


/* Gets the x coords at the start/end of the given column and populate them into the range
 * return argument. */
void getColumnXCoords(GList *columnList, const BlxColumnId columnId, IntRange *xRange)
{
  xRange->min = 0;
  xRange->max = 0;
  
  /* Loop through all visible columns up to the given column, summing their widths. */
  GList *columnItem = columnList;
  
  for ( ; columnItem; columnItem = columnItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(columnItem->data);
      
      if (columnInfo->columnId != columnId)
        {
          if (showColumn(columnInfo))
            xRange->min += columnInfo->width;
        }
      else
        {
          /* We've got to the required column. Calculate the x coord at the end
           * of this column and then break. */
          if (showColumn(columnInfo))
            xRange->max = xRange->min + columnInfo->width;
          else
            xRange->max = xRange->min; /* return zero-width if column is not visible */
          
          break;
        }
    }
}


/* Returns true if the given column should be visible */
gboolean showColumn(BlxColumnInfo *columnInfo)
{
  return (columnInfo->showColumn && columnInfo->dataLoaded && columnInfo->width > 0);
}


/* Save the column widths to the given config file. */
void saveColumnWidths(GList *columnList, GKeyFile *key_file)
{
  /* Loop through each column */
  GList *listItem = columnList;
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);

      if (columnInfo && columnInfo->columnId != BLXCOL_SEQUENCE)
        {
          /* If column is not visible, set width to zero to indicate that it should be hidden */
          if (columnInfo->showColumn)
            g_key_file_set_integer(key_file, COLUMN_WIDTHS_GROUP, columnInfo->title, columnInfo->width);
          else
            g_key_file_set_integer(key_file, COLUMN_WIDTHS_GROUP, columnInfo->title, 0);
        }
    }
}


/* Save the summary columns to the given config file. */
void saveSummaryColumns(GList *columnList, GKeyFile *key_file)
{
  /* Loop through each column */
  GList *listItem = columnList;
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);

      if (columnInfo && columnInfo->canShowSummary)
        {
          g_key_file_set_boolean(key_file, COLUMN_SUMMARY_GROUP, columnInfo->title, columnInfo->showSummary);
        }
    }
}


/* Reset column widths to default values */
void resetColumnWidths(GList *columnList)
{
  /* Quick and dirty: just set the width and visibility manually. This should be
   * done differently so we can just loop through and find the correct values somehow;
   * currently we run the risk of getting out of sync with the initial values set in
   * blxCreateColumns */
  getColumnInfo(columnList, BLXCOL_SEQNAME)->width = BLXCOL_SEQNAME_WIDTH;
  getColumnInfo(columnList, BLXCOL_SEQNAME)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_SCORE)->width = BLXCOL_SCORE_WIDTH;
  getColumnInfo(columnList, BLXCOL_SCORE)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_ID)->width = BLXCOL_ID_WIDTH;
  getColumnInfo(columnList, BLXCOL_ID)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_START)->width = BLXCOL_START_WIDTH;
  getColumnInfo(columnList, BLXCOL_START)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_SEQUENCE)->width = BLXCOL_SEQUENCE_WIDTH;
  getColumnInfo(columnList, BLXCOL_SEQUENCE)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_END)->width = BLXCOL_END_WIDTH;
  getColumnInfo(columnList, BLXCOL_END)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_SOURCE)->width = BLXCOL_SOURCE_WIDTH;
  getColumnInfo(columnList, BLXCOL_SOURCE)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_GROUP)->width = BLXCOL_GROUP_WIDTH;
  getColumnInfo(columnList, BLXCOL_GROUP)->showColumn = FALSE;

  getColumnInfo(columnList, BLXCOL_ORGANISM)->width = BLXCOL_ORGANISM_WIDTH;
  getColumnInfo(columnList, BLXCOL_ORGANISM)->showColumn = TRUE;

  getColumnInfo(columnList, BLXCOL_GENE_NAME)->width = BLXCOL_GENE_NAME_WIDTH;
  getColumnInfo(columnList, BLXCOL_GENE_NAME)->showColumn = FALSE;

  getColumnInfo(columnList, BLXCOL_TISSUE_TYPE)->width = BLXCOL_TISSUE_TYPE_WIDTH;
  getColumnInfo(columnList, BLXCOL_TISSUE_TYPE)->showColumn = FALSE;

  getColumnInfo(columnList, BLXCOL_STRAIN)->width = BLXCOL_STRAIN_WIDTH;
  getColumnInfo(columnList, BLXCOL_STRAIN)->showColumn = FALSE;

}



/***************** end of file ***********************/
