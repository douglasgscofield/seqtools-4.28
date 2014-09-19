/*  File: blxmsp.c
 *  Author: Gemma Barson, 2010-09-02
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
 * Description: See blxmsp.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/utilities.h>
#include <string.h>


#define POLYA_TAIL_BASES_TO_CHECK -1 /* number of bases to check when looking for a polyA tail (-1
                                        means check all of the unaligned sequence) */


/* Globals */
static int g_MaxMspLen = 0;                   /* max length in display coords of all MSPs in the detail-view */
static BlxDataType *g_DefaultDataType = NULL; /* data type containing default values; used if sequences do not have a data-type specified */

/* The config value keys for each flag in BlxDataType. 
 * Use NULL if you don't want the value to be configurable via the config file.
 * THIS ARRAY MUST BE UPDATED IF YOU ADD ITEMS TO THE MspFlag ENUM */
static const char* g_MspFlagConfigKeys[] = 
  {
    "dummy", /* dummy value for MSPFLAG_MIN */

    "link-features-by-name",
    "squash-linked-features",
    "squash-identical-features",
    "strand-specific",
    "show-reverse-strand",
    
    "dummy" /* dummy value for MSPFLAG_NUM_FLAGS */
  };


static void addBlxSequences(const char *name, const char *name_orig, const char *idTag, 
                            BlxStrand strand, BlxDataType *dataType, const char *source, 
                            GArray *featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, 
                            GList *columnList, char *sequence, 
                            MSP *msp, GHashTable *lookupTable, GError **error);
static void findSequenceExtents(BlxSequence *blxSeq);
static MSP* createMissingMsp(const BlxMspType newType,
                             const int newStart,
                             const int newEnd,
                             const char *qname,
                             const int newFrame,
                             BlxStyle *newStyle,
                             BlxSequence *blxSeq, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             GList **seqList,
                             GList *columnList,
                             GHashTable *lookupTable,
                             GError **error);

/* Get/set the max MSP length */
int getMaxMspLen()
{
  return g_MaxMspLen;
}

void setMaxMspLen(const int len)
{
  g_MaxMspLen = len;
}


/* Type determination methods */
gboolean typeIsExon(const BlxMspType mspType)
{
  return (mspType == BLXMSP_CDS || mspType == BLXMSP_UTR || mspType == BLXMSP_EXON);
}

gboolean typeIsIntron(const BlxMspType mspType)
{
  return (mspType == BLXMSP_INTRON);
}

gboolean typeIsTranscript(const BlxMspType mspType)
{
  return (mspType == BLXMSP_TRANSCRIPT);
}

gboolean typeIsMatch(const BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH);
}

gboolean typeIsVariation(const BlxMspType mspType)
{
  return (mspType == BLXMSP_VARIATION);
}

gboolean typeIsRegion(const BlxMspType mspType)
{
  return (mspType == BLXMSP_REGION);
}

/* This returns true if the given type should be shown in the detail-view */
gboolean typeShownInDetailView(const BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH || mspType == BLXMSP_CDS || mspType == BLXMSP_UTR);
}

/* This returns true if the given sequence should be shown in the detail-view */
gboolean blxSequenceShownInDetailView(const BlxSequence *blxSeq)
{
  return (blxSeq->type == BLXSEQUENCE_MATCH || blxSeq->type == BLXSEQUENCE_TRANSCRIPT);
}

/* This returns true if the given sequence should be shown in the big picture grids */
gboolean blxSequenceShownInGrid(const BlxSequence *blxSeq)
{
  return (blxSeq->type == BLXSEQUENCE_MATCH);
}

gboolean mspIsExon(const MSP* const msp)
{
  return (msp && typeIsExon(msp->type));
}


/* Determine whether the given msp is in a visible layer */
gboolean mspLayerIsVisible(const MSP* const msp)
{
  gboolean result = TRUE;
  
  /* Currently only applicable to exons. Show plain exons OR their CDS/UTR sections,
   * but not both. The plan is to add some options to toggle layers on and off, but for 
   * now just hard code this. */
  if (msp->type == BLXMSP_EXON /* msp->type == BLXMSP_CDS || msp->type == BLXMSP_UTR */)
    {
    result = FALSE;
    }
  
  return result;
}


/* Determine whether the given MSP is in a coding region or untranslated region. For 
 * exons, this is determined by the exon type. For introns, we have to look at the
 * adjacent exons to determine whether to show them as CDS or UTR - only show it as
 * CDS if there is a CDS exon on both sides of the intron. */
static const GdkColor* mspGetIntronColor(const MSP* const msp, 
                                         GArray *defaultColors,
                                         const int defaultColorId,
                                         const BlxSequence *blxSeq,
                                         const gboolean selected,
                                         const gboolean usePrintColors,
                                         const gboolean fill,
                                         const int exonFillColorId,
                                         const int exonLineColorId,
                                         const int cdsFillColorId,
                                         const int cdsLineColorId,
                                         const int utrFillColorId,
                                         const int utrLineColorId)
{
  const GdkColor *result = NULL;
  
  /* Find the nearest exons before and after this MSP */
  const MSP *prevExon = NULL;
  const MSP *nextExon = NULL;
  
  GList *mspItem = blxSeq->mspList;
  for ( ; mspItem; mspItem = mspItem->next)
    {
    const MSP *curMsp = (const MSP *)(mspItem->data);
    
    if (mspIsExon(curMsp) && mspLayerIsVisible(curMsp))
      {
      const int curOffset = mspGetQStart(curMsp) - mspGetQStart(msp);
      
      if (curOffset < 0 && (!prevExon || curOffset > mspGetQStart(prevExon) - mspGetQStart(msp)))
	{
	/* Current MSP is before our MSP and is the smallest offset so far */
	prevExon = curMsp;
	}
      else if (curOffset > 0 && (!nextExon || curOffset < mspGetQStart(nextExon) - mspGetQStart(msp)))
	{
	/* Current MSP is after our MSP and is the smallest offset so far */
	nextExon = curMsp;
	}
      }
    }
  
  gboolean prevIsUtr = prevExon && prevExon->type == BLXMSP_UTR;
  gboolean nextIsUtr = nextExon && nextExon->type == BLXMSP_UTR;
  
  /* if either exon is UTR, the intron is UTR */
  if (prevIsUtr)
    {
      result = mspGetColor(prevExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextIsUtr)
    {
      result = mspGetColor(nextExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (prevExon)
    {
      /* Both exons (or the sole exon, if only one exists) are CDS, so the intron is CDS */
      result = mspGetColor(prevExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextExon)
    {
      /* This is the only exon and it is CDS, so make the intron CDS */
      result = mspGetColor(nextExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else
    {
    /* No exon exists adjacent to this intron: default to generic exon color for want of anything better to do. */
    if (fill)
      result = getGdkColor(exonFillColorId, defaultColors, selected, usePrintColors);
    else
      result = getGdkColor(exonLineColorId, defaultColors, selected, usePrintColors);
    }
  
  return result;
}

gboolean mspIsIntron(const MSP* const msp)
{
  return (msp && msp->type == BLXMSP_INTRON);
}

gboolean mspIsBlastMatch(const MSP* const msp)
{
  return (msp && (msp->type == BLXMSP_MATCH));
}

gboolean mspIsPolyASite(const MSP* const msp)
{
  return (msp && msp->type == BLXMSP_POLYA_SITE);
}

gboolean mspIsVariation(const MSP* const msp)
{
  return (msp && typeIsVariation(msp->type));
}

gboolean mspIsZeroLenVariation(const MSP* const msp)
{
  gboolean result = mspIsVariation(msp);
  
  if (result)
    {
      const char *seq = mspGetMatchSeq(msp);
      result = (seq && seq[0] == '-');
    }

  return result;
}

/* Whether the msp has a Target name (i.e. Subject sequence name) */
gboolean mspHasSName(const MSP* const msp)
{
  return TRUE;
}

/* Whether the MSP requires subject sequence coords to be set. Only matches 
 * and exons have coords on the subject sequence. (to do: is this optional for exons?) */
gboolean mspHasSCoords(const MSP* const msp)
{
  return mspIsExon(msp) || mspIsBlastMatch(msp);
}

/* Whether the MSP requires subject sequence strand to be set. Only matches 
 * require strand on the subject sequence, although exons may have them set. */
gboolean mspHasSStrand(const MSP* const msp)
{
  return mspIsBlastMatch(msp);
}

/* Whether the MSP requires the actual sequence data for the subject sequence. Only
 * matches require the sequence data. */
gboolean mspHasSSeq(const MSP* const msp)
{
  return mspIsBlastMatch(msp);
}


/***********************************************************
 *		MSP data access functions		   * 
 ***********************************************************/

/* Get the range of coords of the alignment on the reference sequence */
const IntRange* mspGetRefCoords(const MSP* const msp)
{
  return &msp->qRange;
}

/* Get the range of coords of the alignment on the match sequence */
const IntRange* mspGetMatchCoords(const MSP* const msp)
{
  return &msp->sRange;
}

/* Return the length of the range of alignment coords on the ref seq */
int mspGetQRangeLen(const MSP* const msp)
{
  return msp->qRange.max - msp->qRange.min + 1;
}

/* Return the length of the range of alignment coords on the match seq */
int mspGetSRangeLen(const MSP* const msp)
{
  return msp->sRange.max - msp->sRange.min + 1;
}

/* Get the start (5 prime) coord of the alignment on the reference sequence. This is
 * the lowest value coord if the strand is forwards or the highest if it is reverse. */
int mspGetQStart(const MSP* const msp)
{
  return (mspGetRefStrand(msp) == BLXSTRAND_REVERSE ? msp->qRange.max : msp->qRange.min);
}

/* Get the end (3 prime) coord of the alignment on the reference sequence. This is
 * the highest value coord if the strand is forwards or the lowest if it is reverse. */
int mspGetQEnd(const MSP* const msp)
{
  return (mspGetRefStrand(msp) == BLXSTRAND_REVERSE ? msp->qRange.min : msp->qRange.max);
}

/* Get the start coord of the alignment on the match sequence. This is
 * the lowest value coord if the match strand is the same direction as the ref seq strand,
 * or the highest value coord otherwise. */
int mspGetSStart(const MSP* const msp)
{
  return (mspGetMatchStrand(msp) == mspGetRefStrand(msp) ? msp->sRange.min : msp->sRange.max);
}

/* Get the end coord of the alignment on the match sequence. This is
 * the highest value coord if the match strand is in the same direction as the ref seq strand, 
 * or the lowest value coord otherwise. */
int mspGetSEnd(const MSP* const msp)
{
  return (mspGetMatchStrand(msp) == mspGetRefStrand(msp) ? msp->sRange.max : msp->sRange.min);
}


/* Return the match sequence name. (Gets it from the BlxSequence if the MSP itself doesn't have
 * a name) */
const char *mspGetSName(const MSP* const msp)
{
  const char *result = NULL;
  
  if (msp)
    {
      if (msp->sname && msp->sname[0] != 0)
        {
          result = msp->sname;
        }
      else if (msp->sSequence)
        {
          result = blxSequenceGetName(msp->sSequence);
        }
    }
  
  return result;
}

/* Return the match sequence name original name. */
const char *mspGetSNameOrig(const MSP* const msp)
{
  const char *result = NULL;
  
  if (msp)
    {
      if (msp->sname_orig && *(msp->sname_orig))
        {
          result = msp->sname_orig;
        }
    }
  
  return result;
}


/* Return the length of the match sequence that the given MSP lies on */
int mspGetMatchSeqLen(const MSP* const msp)
{
  return blxSequenceGetLength(msp->sSequence);
}

/* Return the reading frame of the ref sequence that the given MSP is a match against */
int mspGetRefFrame(const MSP* const msp, const BlxSeqType seqType)
{
  int result = UNSET_INT;
  
  if (seqType == BLXSEQ_DNA)
    {
    /* Ignore the frame in  the msp. For DNA matches we only have one frame on each strand. */
    result = 1;
    }
  else
    {
    result = msp->qFrame;
    }
  
  return result;
}

/* Return the reference sequence name that this msp aligns against */
const char* mspGetRefName(const MSP* const msp)
{
  return msp->qname;
}

/* Return the strand of the ref sequence that the given MSP is a match against */
BlxStrand mspGetRefStrand(const MSP* const msp)
{
  BlxStrand result = msp->qStrand;
  
  /* If not strand specific, always return the forward strand */
  if (!mspGetFlag(msp, MSPFLAG_STRAND_SPECIFIC))
    result = BLXSTRAND_FORWARD;
  
  return result;
}

/* Return the strand of the match sequence that the given MSP is a match on */
BlxStrand mspGetMatchStrand(const MSP* const msp)
{
  BlxStrand result = (msp->sSequence ? msp->sSequence->strand : BLXSTRAND_NONE);

  /* If not strand specific, always return the forward strand */
  if (!mspGetFlag(msp, MSPFLAG_STRAND_SPECIFIC))
    result = BLXSTRAND_FORWARD;
  
  return result;
}

/* Get the match sequence for the given MSP */
const char* mspGetMatchSeq(const MSP* const msp)
{
  return (msp ? blxSequenceGetSequence(msp->sSequence) : NULL);
}

/* Get the source of the MSP */
const char* mspGetSource(const MSP* const msp)
{
  return blxSequenceGetSource(msp->sSequence);
}


/* Return the color matching the given properties from the given style */
static const GdkColor *styleGetColor(const BlxStyle *style, 
                                     const gboolean selected, 
                                     const gboolean usePrintColors,
                                     const gboolean fill, 
                                     const gboolean utr,
                                     GArray *defaultColors,
                                     const int defaultColorId)
{
  const GdkColor *result = NULL;
  const BlxColor *blxColor = NULL;
  
  if (fill)
    {
      if (utr)
        blxColor = &style->fillColorUtr;
      else
        blxColor = &style->fillColor;
    }
  else
    {
      if (utr)
        blxColor = &style->lineColorUtr;
      else
        blxColor = &style->lineColor;
    }

  /* If it's transparent, use the background color instead, unless 
   * selected is true in which case we need to use the highlight color */
  if (blxColor && blxColor->transparent && !selected)
    result = getGdkColor(defaultColorId, defaultColors, selected, usePrintColors);
  else if (blxColor)
    result = blxColorGetColor(blxColor, selected, usePrintColors);

  return result;
}


/* Get the color for drawing the given MSP (If 'selected' is true, returns
 * the color when the MSP is selected.). Returns the fill color if 'fill' is 
 * true, otherwise the line color */
const GdkColor* mspGetColor(const MSP* const msp, 
			    GArray *defaultColors, 
                            const int defaultColorId,
			    const BlxSequence *blxSeq,
			    const gboolean selected, 
			    const gboolean usePrintColors, 
			    const gboolean fill,
			    const int exonFillColorId,
			    const int exonLineColorId,
			    const int cdsFillColorId,
			    const int cdsLineColorId,
			    const int utrFillColorId,
			    const int utrLineColorId)
{
  const GdkColor *result = NULL;
  
  if (msp->style)
    {
      result = styleGetColor(msp->style, selected, usePrintColors, fill, msp->type == BLXMSP_UTR, defaultColors, defaultColorId);
    }
  
  if (!result)
    {
    /* Use the default color for this MSP's type */
    switch (msp->type)
      {
	case BLXMSP_EXON:
	result = getGdkColor(fill ? exonFillColorId : exonLineColorId, defaultColors, selected, usePrintColors);
	break;
	
	case BLXMSP_CDS:
	result = getGdkColor(fill ? cdsFillColorId : cdsLineColorId, defaultColors, selected, usePrintColors);
	break;
	
	case BLXMSP_UTR:
	result = getGdkColor(fill ? utrFillColorId : utrLineColorId, defaultColors, selected, usePrintColors);
	break;
	
        /* to do: mspGetIntronColor() is non-trivial, because it has to work out the color
         * from the adjacent exons. Since mspGetColor() is called many times on re-draw, it
         * would be better to work out whether an intron is CDS or UTR during initialisation 
         * and use different types (e.g. BLXMSP_INTRON_CDS) that can be queried here to quickly
         * determine what color to use. */
	case BLXMSP_INTRON:
          result = mspGetIntronColor(msp, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
	break;
	
	default:
	break;
      };
    }
  
  return result;
}


/* Get functions from msps for various sequence properties. Result is owned by the sequence
 * and should not be free'd. Returns NULL if the property is not set. */
const char *mspGetColumn(const MSP* const msp, const BlxColumnId columnId)
{
  return (msp ? blxSequenceGetColumn(msp->sSequence, columnId) : NULL);
}

const char *mspGetOrganism(const MSP* const msp)
{
  return (msp ? blxSequenceGetOrganism(msp->sSequence) : NULL);
}

const char *mspGetOrganismAbbrev(const MSP* const msp)
{
  return (msp ? blxSequenceGetOrganismAbbrev(msp->sSequence) : NULL);
}

const char *mspGetGeneName(const MSP* const msp)
{
  return (msp ? blxSequenceGetGeneName(msp->sSequence) : NULL);
}

const char *mspGetTissueType(const MSP* const msp)
{
  return (msp ? blxSequenceGetTissueType(msp->sSequence) : NULL);
}

const char *mspGetStrain(const MSP* const msp)
{
  return (msp ? blxSequenceGetStrain(msp->sSequence) : NULL);
}


/* Return the coords of an MSP as a string. The result should be free'd with g_free */
char *mspGetCoordsAsString(const MSP* const msp)
{
  char *result = NULL;
  
  if (msp)
    {
      GString *resultStr = g_string_new("");

      /* If both s coords are UNSET_INT then they are not relevant, so exclude them */
      if (msp->sRange.min == UNSET_INT && msp->sRange.max == UNSET_INT)
        g_string_append_printf(resultStr, "%d,%d", msp->qRange.min, msp->qRange.max);
      else
        g_string_append_printf(resultStr, "%d,%d[%d,%d]", msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
      
      result = g_string_free(resultStr, FALSE);
    }
  
  return result;
}


/* Return the path in the given tree model that this MSP lies in */
gchar* mspGetTreePath(const MSP* const msp, BlxModelId modelId)
{
  return msp->treePaths[modelId];
}


///* Returns true if a feature-series by the given name exists in the feature-series array and
// * and, if so, sets index_out with its index. */
//static gboolean fsArrayFindByName(GArray *fsArray, FeatureSeries *fs, int *index_out)
//{
//  gboolean result = FALSE;
//  
//  if (fsArray)
//    {
//    int i = 0;
//    for ( ; i < fsArray->len; ++i)
//      {
//      FeatureSeries *compareFs = &g_array_index(fsArray, FeatureSeries, i);
//      
//      if (!fsSortByNameCompareFunc(fs, compareFs))
//	{
//	result = TRUE;
//	*index_out = i;
//	break;
//	}
//      }
//    }
//  
//  return result;
//}


/* Comparison function to sort two Feature Series by the order number stored in the FeatureSeries
 * struct. Returns -1 if the first item is before the second, 1 if the second is first, or 0 if 
 * they are equal.  */
gint fsSortByOrderCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in)
{
  int result = 0;
  
  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;
  
  if (fs1->order < fs2->order)
    result = -1;
  else if (fs1->order > fs2->order)
    result = 1;
  
  return result;
}


/* Comparison function to sort two Features Series by name. */
gint fsSortByNameCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in)
{
  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;
  
  /*printf("%s - %s : %d\n", fs1->name, fs2->name,  strcmp(fs1->name, fs2->name));*/
  
  return strcmp(fs1->name, fs2->name);
}


/* Insert the given MSP into the Feature Series of the given name. If the Feature Series
 * does not exist yet, create it. Also, if the Feature Series array does not exist yet, create
 * that too. */
//void insertFS(MSP *msp, char *series)
//{
//  if (!fsArr) 
//    {
//    fsArr = g_array_sized_new(TRUE, FALSE, sizeof(FeatureSeries), 50);
//    }
//  
//  static int orderNum = 0; /* will increment this each time we add a feature series to the array */
//  
//  FeatureSeries *fs = g_malloc(sizeof(FeatureSeries));
//  fs->on = 1;
//  fs->y = 0.0;
//  fs->xy = (msp->type == BLXMSP_XY_PLOT ? 1 : 0);
//  
//  fs->name = g_strdup(series);
//  
//  int i;
//  if (fsArrayFindByName(fsArr, fs, &i))
//    {
//    msp->fs = &g_array_index(fsArr, FeatureSeries, i);
//    g_free(fs->name);
//    g_free(fs);
//    }
//  else
//    {
//    /* Remember the order we added them so we can sort by it again later. */
//    orderNum++;
//    fs->order = orderNum;
//    
//    g_array_append_val(fsArr, *fs);
//    //      g_array_sort(fsArr, fsSortByNameCompareFunc);
//    
//    msp->fs = fs;
//    }
//}
//


/* Get the MSP at the given index in the given array. Returns null if out of bounds */
MSP* mspArrayIdx(const GArray* const array, const int idx)
{
  MSP *msp = NULL;
  
  if (idx >= 0 && idx < (int)array->len)
    msp = g_array_index(array, MSP*, idx);
  
  return msp;
}


/* Check if the given character is a polyA character (i.e. 'a', or if the strand is reverse 't') */
static gboolean isPolyAChar(const char c, const BlxStrand strand)
{
  gboolean result = FALSE;

  if (strand == BLXSTRAND_FORWARD)
    result = (c == 'a' || c == 'A');
  else
    result = (c == 't' || c == 'T');

  return result;
}


/* Get the number of bases to check when looking for a polyA tail. We may want to make this configurable? */
static int getNumPolyATailBasesToCheck()
{
  static int result = POLYA_TAIL_BASES_TO_CHECK;
  return result;
}


/* Returns true if there is a polyA site at the 3' end of this MSP's alignment range. */
gboolean mspHasPolyATail(const MSP* const msp)
{
  gboolean found = FALSE;
  
  /* Only matches have polyA tails. */
  if (mspIsBlastMatch(msp))
    {
      const char *seq = mspGetMatchSeq(msp);
      
      if (seq)
        {
          const int numRequired = getNumPolyATailBasesToCheck();
          const int minRequired = 3; /* check at least 3 bases */
          const int len = strlen(seq);
          BlxStrand sStrand = mspGetMatchStrand(msp);
          BlxStrand qStrand = mspGetRefStrand(msp);
          int sCoord = mspGetSEnd(msp);

          if (qStrand == sStrand) 
            {
              ++sCoord; /* next coord after alignment block end */
              int sMax = mspGetSStart(msp);

              if (numRequired > 0) /* -1 means check all of the unaligned sequence */
                sMax = sCoord + numRequired;

              if (sMax <= len && sMax - sCoord >= minRequired)
                {
                  found = TRUE;

                  for ( ; sCoord <= sMax; ++sCoord)
                    {
                      if (!isPolyAChar(seq[sCoord - 1], qStrand))
                        {
                          found = FALSE;
                          break;
                        }
                    }
                }
            }
          else
            {
              --sCoord; /* next coord after alignment block end */
              int sMin = 1;
              
              if (numRequired > 0)
                sMin = sCoord - numRequired;

              if (sMin >= 1 && sCoord - sMin >= minRequired)
                {
                  found = TRUE;

                  for ( ; sCoord >= sMin; --sCoord)
                    {
                      if (!isPolyAChar(seq[sCoord - 1], qStrand))
                        {
                          found = FALSE;
                          break;
                        }
                    }
                }
            }
        }
    }
  
  return found;
}


/* Returns true if the given MSP coord (in ref seq nucleotide coords) is inside a polyA tail, if
 * this MSP has one. */
gboolean mspCoordInPolyATail(const int coord, const MSP* const msp)
{
  gboolean result = mspHasPolyATail(msp);
  
  /* See if the coord is outside the 3' end of the alignment range (i.e. is greater than the
   * max coord if we're on the forward strand or less than the min coord if on the reverse). */
  //result &= ((mspGetRefStrand(msp) == BLXSTRAND_FORWARD && coord > msp->displayRange.max) ||
  //           (mspGetRefStrand(msp) == BLXSTRAND_REVERSE && coord < msp->displayRange.min));
  result &= coord > msp->displayRange.max;

  return result;
}


/***********************************************************
 *		      BlxSequence			   * 
 ***********************************************************/

/* Append the contents of the given text to the GString, if the text is non-null,
 * including the given separator before the text is appended. Pass the separator
 * as an empty string if no separation is required. */
static void appendTextIfNonNull(GString *gstr, const char *separator, const char *text)
{
  if (text)
    {
      g_assert(separator);
      g_string_append_printf(gstr, "%s%s", separator, text);
    }
}

/* Return summary info about a given BlxSequence (e.g. for displaying in the status bar). The
 * result should be free'd with g_free. */
char* blxSequenceGetSummaryInfo(const BlxSequence* const blxSeq, GList *columnList)
{
  char *result = NULL;
  
  if (blxSeq)
    {
      GString *resultStr = g_string_new("");
      const char *separator = "";

      /* Loop through all the columns, appending any that should be shown
       * to the result string */
      GList *item = columnList;

      for ( ; item; item = item->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
          
          if (columnInfo->showSummary)
            {
              const char *valueText = blxSequenceGetValueAsString(blxSeq, columnInfo->columnId);

              if (valueText)
                {
                  appendTextIfNonNull(resultStr, separator, valueText);

                  /* The first item has no separator, but subsequence items do,
                   * so set the separator after the first item has been added. */
                  separator = ";  ";
                }
            }
        }
      
      result = g_string_free(resultStr, FALSE);
    }
  
  return result;
}


/* Return the full name of a BlxSequence (including prefix and variant) */
const char *blxSequenceGetName(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq)
    {
      result = blxSequenceGetValueAsString(seq, BLXCOL_SEQNAME);
      
      if (!result && seq->idTag)
        result = seq->idTag;

      if (!result)
        g_warning("Sequence does not have a name specified.\n");
    }
  
  return result;
}

/* Return the Source text of a BlxSequence, if it has one (note that it gets
 * this from the first MSP and does no checking whether other MSPs have the 
 * same source or not). */
const char *blxSequenceGetSource(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq)
    result = blxSequenceGetValueAsString(seq, BLXCOL_SOURCE);
  
  return result;
}

/* Return the fetch method of a BlxSequence. If 'bulk' is true,
 * get the bulk-fetch method, otherwise the user-fetch method.
 * 'index' indicates which method to choose if multiple methods
 * are available; 0 is the first (preferred) method, 1 the second 
 * etc. If no fetch method is set, return the given default method 
 * instead (for index==0 only). Pass 'defaultMethod' as 0 if N/A. */
GQuark blxSequenceGetFetchMethod(const BlxSequence *seq, 
                                 const gboolean bulk,
                                 const gboolean optionalColumns,
                                 const int index,
                                 const GArray *defaultMethods)
{
  GQuark result = 0;

  if (seq && seq->dataType)
    {
      GArray *array = NULL;

      if (bulk && optionalColumns)
        array = seq->dataType->optionalFetch;
      else if (bulk)
        array = seq->dataType->bulkFetch;
      else
        array = seq->dataType->userFetch;
      
      if (array && index >= 0 && index < (int)array->len)
        result = g_array_index(array, GQuark, index);
    }

  if (!result && defaultMethods && index >= 0 && index < (int)defaultMethods->len) 
    {
      result = g_array_index(defaultMethods, GQuark, index);
    }
  
  return result;
}

/* Return the length of the given blxsequence's sequence data */
int blxSequenceGetLength(const BlxSequence *seq)
{
  int result = 0;
  
  if (seq)
    {
      const char *sequence = blxSequenceGetSequence(seq);

      if (sequence)
        result = strlen(sequence);
    }
  
  return result;
}

/* Get the start extent of the alignments from this match sequence on the
 * given reference sequence strand */
int blxSequenceGetStart(const BlxSequence *seq, const BlxStrand strand)
{
  return (strand == BLXSTRAND_REVERSE ? seq->qRangeRev.min : seq->qRangeFwd.min);
}

/* Get the end extend of the sequence on the ref sequence */
int blxSequenceGetEnd(const BlxSequence *seq, const BlxStrand strand)
{
  return (strand == BLXSTRAND_FORWARD ? seq->qRangeRev.max : seq->qRangeFwd.max);
}

/* Get the sequence data for the given blxsequence */
const char *blxSequenceGetSequence(const BlxSequence *seq)
{
  const char *result = blxSequenceGetValueAsString(seq, BLXCOL_SEQUENCE);
  return result;
}

/* This returns true if the given sequence object requires sequence data
 * (i.e. the actual dna or peptide sequence string) */
gboolean blxSequenceRequiresSeqData(const BlxSequence *seq)
{
  return (seq && 
          (seq->type == BLXSEQUENCE_MATCH || 
           seq->type == BLXSEQUENCE_VARIATION || 
           seq->type == BLXSEQUENCE_REGION));
}

/* This returns true if the given sequence object requires optional data such
 * as organism and tissue-type, if this data is available. */
gboolean blxSequenceRequiresOptionalData(const BlxSequence *seq)
{
  return (seq && seq->type == BLXSEQUENCE_MATCH);
}

/* This returns true if the given sequence object requires the data for
 * the given column */
gboolean blxSequenceRequiresColumnData(const BlxSequence *seq, const BlxColumnId columnId)
{
  gboolean result = FALSE;

  if (seq)
    {
      if (columnId == BLXCOL_ORGANISM || columnId == BLXCOL_TISSUE_TYPE ||
          columnId == BLXCOL_GENE_NAME || columnId == BLXCOL_STRAIN ||
          columnId >= BLXCOL_NUM_COLUMNS)
        {
          result = blxSequenceRequiresOptionalData(seq);
        }
      else if (columnId == BLXCOL_SEQUENCE)
        {
          result = blxSequenceRequiresSeqData(seq);
        }
    }
  
  return result;
}

/* Get the value for the given column */
GValue* blxSequenceGetValue(const BlxSequence *seq, const int columnId)
{
  GValue *result = NULL;
  
  if (seq && seq->values && columnId < (int)seq->values->len)
    result = &g_array_index(seq->values, GValue, columnId);

  return result;
}

void blxSequenceSetValue(const BlxSequence *seq, const int columnId, GValue *value_in)
{
  if (seq && seq->values)
    {
      GValue *value = blxSequenceGetValue(seq, columnId);
      g_value_reset(value);
      g_value_copy(value_in, value);
    }
}

/* Clear the idTag stored in the BlxSequence. This gets done after importing a GFF file is
 * complete because we don't want the IDs to clash with further files which may be imported
 * later. If the name is not set, this sets the name to the ID before clearing the ID (really we
 * should probably be trying to construct a unique name from the source, type, strand and coords,
 * but we don't currently store the GFF type). */
static void blxSequenceClearGFFIds(BlxSequence *seq)
{
  if (seq)
    {
      if (!blxSequenceGetName(seq))
        {
          if (seq->idTag)
            blxSequenceSetValueFromString(seq, BLXCOL_SEQNAME, seq->idTag);
          else
            blxSequenceSetValueFromString(seq, BLXCOL_SEQNAME, "<no name>");
        }

      if (seq->idTag)
        {
          g_free(seq->idTag);
          seq->idTag = NULL;
        }
    }
}

/* Set the value for the given column. The given string is converted to
 * the relevant value type */
void blxSequenceSetValueFromString(const BlxSequence *seq, const int columnId, const char *inputStr)
{
  if (inputStr && seq && seq->values)
    {
      GValue *value = blxSequenceGetValue(seq, columnId);

      if (value)
        {
          g_value_reset(value);

          if (G_VALUE_HOLDS_STRING(value))
            {
              g_value_take_string(value, g_strdup(inputStr));
            }
          else if (G_VALUE_HOLDS_INT(value))
            {
              const int tmp = atoi(inputStr);
              g_value_set_int(value, tmp);
            }
          else if (G_VALUE_HOLDS_DOUBLE(value))
            {
              const gdouble tmp = g_strtod(inputStr, NULL);
              g_value_set_double(value, tmp);
            }
          else
            {
              g_warning("Tried to set value of unknown type for column '%d'\n", columnId);
            }
        }
    }
}

/* Get the string value for the given column. Returns null if the
 * value is not a string type. */
const char* blxSequenceGetValueAsString(const BlxSequence *seq, const int columnId)
{
  const char *result = NULL;
  
  GValue *value = blxSequenceGetValue(seq, columnId);

  if (value)
    {
      if (G_VALUE_HOLDS_STRING(value))
        result = g_value_get_string(value);
      else if (G_VALUE_HOLDS_INT(value))
        result = convertIntToString(g_value_get_int(value));
      else if (G_VALUE_HOLDS_DOUBLE(value))
        result = convertDoubleToString(g_value_get_double(value), 2);
    }
  
  /* Return null if it's an empty value (i.e. if it's unset) */
  if (result && *result == 0)
    result = NULL;

  return result;
}

/* Get the relevant data about a sequence for the given column ID (at the moment
 * this only supports string data) */
const char* blxSequenceGetColumn(const BlxSequence* const blxSeq, const BlxColumnId columnId)
{
  const char *result = NULL;

  GValue *value = blxSequenceGetValue(blxSeq, columnId);
    
  if (value && G_VALUE_HOLDS_STRING(value))
    {
      result = g_value_get_string(value);
    }
  
  return result;
}

const char *blxSequenceGetOrganism(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_ORGANISM);
}

const char *blxSequenceGetOrganismAbbrev(const BlxSequence *seq)
{
  return (seq ? seq->organismAbbrev : NULL);
}

const char *blxSequenceGetGeneName(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_GENE_NAME);
}

const char *blxSequenceGetTissueType(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_TISSUE_TYPE);
}

const char *blxSequenceGetStrain(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_STRAIN);
}

/* Return the sequence data as a string in fasta format.
 * Result should be free'd with g_free, unless the name
 * or sequence couldn't be found, in which case the result
 * is null */
char *blxSequenceGetFasta(const BlxSequence *seq)
{
  char *result = NULL;
  
  if (seq)
    {
      const char *name = blxSequenceGetName(seq);
      const char *sequence = blxSequenceGetSequence(seq);
      
      if (name && sequence)
        {
          result = g_strdup_printf(">%s\n%s", name, sequence);
        }
    }
  
  return result;
}


/* Return all the stored info about a blx sequenece (description, organism, tissue type etc.) 
 * in a single string. The result should be free'd by the caller using g_free. If 'allowNewlines'
 * is true the data is separated with newline characters, otherwise with tabs (i.e. pass as false
 * if returned string is to be shown as a single line). The dataLoaded flag should be passed as
 * false if the optional data has not been loaded yet. */
char *blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, GList *columnList)
{
  char *result = NULL;
  
  if (blxSeq)
    {
      GString *resultStr = g_string_new("");
      char separator = allowNewlines ? '\n' : '\t';
      char strand = blxSeq->strand == BLXSTRAND_REVERSE ? '-' : '+';

      /* Loop through all the columns, appending them to the result string */
      GList *item = columnList;
      const int titleWidth = 17; /* if title is less than this width then it will be padded */

      for ( ; item; item = item->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
          
          const char *valueText = blxSequenceGetValueAsString(blxSeq, columnInfo->columnId);
          const char *text = valueText ? valueText : "";
              
          if (columnInfo->columnId == BLXCOL_SEQNAME)
            g_string_append_printf(resultStr, "%-*s  %s%c%c", titleWidth, columnInfo->title, text, strand, separator);
          else
            g_string_append_printf(resultStr, "%-*s  %s%c", titleWidth, columnInfo->title, text, separator);
        }
      
      /* Loop through the child features and show their coords */
      const char *title = "Coords";
      g_string_append_printf(resultStr, "%-*s  ", titleWidth, title);
      GList *mspItem = blxSeq->mspList;
      
      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);
          
          /* Don't show 'exon' msps because we show the exon's
           * individual 'cds' and 'utr' instead */
          if (msp->type != BLXMSP_EXON)
            {
              char *coordsStr = mspGetCoordsAsString(msp);
              g_string_append_printf(resultStr, "%s  ", coordsStr ? coordsStr : "");
              
              if (coordsStr)
                g_free(coordsStr);
            }
        }
      
      result = g_string_free(resultStr, FALSE);
    }
  
  return result;
}

/* Find a BlxSequence by its full name (must be an exact match but is case insensitive) */
static BlxSequence* blxSequenceFindByName(const char *name, GList *allSeqs)
{
  BlxSequence *result = NULL;
  GList *listItem = allSeqs;
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxSequence *curSeq = (BlxSequence*)(listItem->data);
      const char *curName = blxSequenceGetName(curSeq);
    
      if (name && stringsEqual(curName, name, FALSE))
        {
          result = curSeq;
          break;
        }
    }
  
  return result;
}


/* Get the "parent" sequence of the given protein variant. Assumes the variant
 * contains a dash '-' in the name followed by the variant number as a digit. 
 * e.g. SW:P51531-2.2. The function looks for a sequence with the same name but
 * with this dash and the following digit(s) (up to the end of the name or the '.'
 * if there is one) removed. Returns NULL if no parent was found. */
BlxSequence* blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs)
{
  BlxSequence *result = NULL;
  
  const char *variantName = blxSequenceGetName(variant);
  
  if (variantName)
    {
      char *parentName = g_strdup(variantName);
      char *insertPoint = strchr(parentName, '-');
      
      if (insertPoint)
        {
          /* Replace '-' by terminating char, in case there's nothing else to copy in. */
          *insertPoint = '\0';
          
          /* The insert point is where we'll copy into. Create another pointer that we'll increment
           * until we find a '.' and then we'll copy from that point. */
          char *copyPoint = insertPoint;
          ++copyPoint;
          
          gboolean foundRestartPoint = FALSE; /* set to true when we find where to start copying from again */
          
          while (copyPoint && *copyPoint != '\0')
            {
              if (foundRestartPoint)
                {
                  *insertPoint = *copyPoint;
                  ++insertPoint;
                }
              else if (*copyPoint == '.')
                {
                  foundRestartPoint = TRUE;
                  *insertPoint = *copyPoint;
                  ++insertPoint;
                }
          
              ++copyPoint;
            }
        
          *insertPoint = '\0';
        
          result = blxSequenceFindByName(parentName, allSeqs);
          g_free(parentName);
        }
    }
  
  return result;
}


/* Destroy all of the BlxSequences */
void destroyBlxSequenceList(GList **seqList)
{
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      destroyBlxSequence(blxSeq);
    }
  
  g_list_free(*seqList);
  *seqList = NULL;
}


/* Frees all memory used by a BlxSequence */
void destroyBlxSequence(BlxSequence *seq)
{
  if (seq)
    {
      if (seq->values)
        {
          /* Free all column values */
          int i = 0;      
          for ( ; i < (int)seq->values->len; ++i)
            {
              GValue *value = &g_array_index(seq->values, GValue, i);
              g_value_unset(value);
            }
        }

      if (seq->idTag)
        g_free(seq->idTag);
          
      g_free(seq);
    }
}


/* Destroy an msp and remove it from the BlxSequence and msp lists */
static void destroyMspFull(MSP *msp, BlxSequence *seq, GArray *featureLists[], MSP **lastMsp, MSP **mspList)
{
  /* Remove from the feature list */
  GArray *array = featureLists[msp->type];
  int i = 0;
  for ( ; i < array->len; ++i)
    {
      MSP *curMsp = g_array_index(array, MSP*, i);
      
      if (curMsp == msp)
        {
          array = g_array_remove_index(array, i);
          break;
        }
    }
  featureLists[msp->type] = array;
  
  /* Remove from mspList */
  MSP *curMsp = *mspList;
  MSP *prevMsp = NULL;
  for ( ; curMsp; curMsp = curMsp->next)
    {
      if (msp == curMsp)
        {
          if (!prevMsp)
            *mspList = curMsp->next; /* Remove msp from start of list */
          else 
            prevMsp->next = curMsp->next; /* Remove link to msp */
          
          if (*lastMsp == msp)
            *lastMsp = prevMsp;  /* Update pointer to last msp */

          break;
        }
      
      prevMsp = curMsp;
    }
  
  /* Remove from the BlxSequence, if given */
  if (seq && seq->mspList)
    {
      GList *item = seq->mspList;

      for ( ; item; item = item->next)
        {
          MSP *curMsp = (MSP*)(item->data);
          if (curMsp == msp)
            {
              seq->mspList = g_list_remove_link(seq->mspList, item);
              break;
            }
        }
    }

  destroyMspData(msp);
}


/* destroy a blxsequence and it's msps, and remove them from the feature lists */
static void destroyBlxSequenceFull(BlxSequence *seq, GArray *featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList)
{
  /* destroy each msp */
  GList *mspItem = seq->mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      MSP *msp = (MSP*)(mspItem->data);
      destroyMspFull(msp, NULL, featureLists, lastMsp, mspList); /* don't remove from
                                                                    seq->mspList because we
                                                                    destroy this anyway */
    }

  g_list_free(seq->mspList);

  /* Remove BlxSequence from seqList */
  *seqList = g_list_remove(*seqList, seq);

  destroyBlxSequence(seq);
}


/* Set the value for a particular column (by column name) */
void blxSequenceSetColumn(BlxSequence *seq, const char *colName, const char *value, GList *columnList)
{
  if (!colName || !value)
    return;
  
  gboolean found = FALSE;

  /* Loop through the column list and find the one with this name */
  GList *item = columnList;
  for ( ; item && !found; item = item->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
      
      if (stringsEqual(colName, columnInfo->title, FALSE))
        {
          blxSequenceSetValueFromString(seq, columnInfo->columnId, value);
          found = TRUE;
        }
    }
  
  if (!found)
    g_warning("Unrecognised column '%s'\n", colName);
}


/* Utility to create a BlxSequence with the given name. */
BlxSequence* createEmptyBlxSequence()
{
  BlxSequence *seq = (BlxSequence*)g_malloc(sizeof(BlxSequence));
  
  seq->type = BLXSEQUENCE_UNSET;
  seq->dataType = NULL;
  seq->idTag = NULL;
  seq->mspList = NULL;
  seq->sequenceReqd = FALSE;
  seq->values = NULL;
  seq->organismAbbrev = NULL;

  return seq;
}


/* Copies all fields in a sequence. Copies all MSPs apart from CDSs whose name does not match the
 * given quark */
static void copyBlxSequenceNamedCds(const BlxSequence *src, 
                                    const GQuark cdsQuark,
                                    GArray *featureLists[], 
                                    MSP **lastMsp, 
                                    MSP **mspList,
                                    GList **seqList,
                                    GList *columnList,
                                    GHashTable *lookupTable, 
                                    GError **error)
{
  GError *tmpError = NULL;
  
  /* We must give the new BlxSequence a unique id - use the cds name */
  const char *idTag = g_quark_to_string(cdsQuark);

  /* We'll copy these values directly from the source */
  const char *source = blxSequenceGetSource(src);
  const BlxStrand sStrand = src->strand;
  BlxDataType *dataType = src->dataType;
  
  /* Copy all MSPs except CDSs whose name does not match cdsQuark */
  GList *mspItem = src->mspList;
  
  for ( ; mspItem && !tmpError; mspItem = mspItem->next)
    {
      const MSP* msp = (const MSP*)(mspItem->data);
      
      if (msp->type != BLXMSP_CDS || g_quark_from_string(msp->sname) == cdsQuark)
        {
          MSP *newMsp = copyMsp(msp, featureLists, lastMsp, mspList, FALSE);

          /* Add the new msp to the new blx sequence (this creates it if it does not exist
           * i.e. the first time we get here for this idTag) */
          addBlxSequence(newMsp->sname, newMsp->sname_orig, idTag, sStrand, dataType, 
                         source, seqList, columnList, 
                         NULL, newMsp, lookupTable, &tmpError);
        }
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


BlxDataType* createBlxDataType()
{
  BlxDataType *result = (BlxDataType*)g_malloc(sizeof *result);
  
  result->name = 0;
  result->bulkFetch = NULL;
  result->userFetch = NULL;

  /* Loop through the flags, setting them all to false initially */
  int flag = MSPFLAG_MIN;
  for ( ; flag < MSPFLAG_NUM_FLAGS; ++flag)
    {
      result->flags[flag] = FALSE;
    }
  
  /* Set any specific flags that we want to be true by default */
  result->flags[MSPFLAG_SQUASH_LINKED_FEATURES] = TRUE;
  result->flags[MSPFLAG_STRAND_SPECIFIC] = TRUE;
  result->flags[MSPFLAG_SHOW_REVERSE_STRAND] = TRUE;

  return result;
}

void destroyBlxDataType(BlxDataType **blxDataType)
{
  if (!blxDataType)
    return;
  
  g_free((*blxDataType));
  *blxDataType = NULL;
}

/* The following functions get the data type values as strings */
const char* getDataTypeName(BlxDataType *blxDataType)
{
  return blxDataType ? g_quark_to_string(blxDataType->name) : NULL;
}


/* does the work for compareFuncMspPos and compareFuncMspArray */
static gint compareMsps(const MSP* const msp1, const MSP* const msp2)
{
  gint result = 0;
  
  if (result == 0)
    {
      if (msp1->qRange.min == msp2->qRange.min)
        {
          /* Sort by type. Lower type numbers should appear first. */
          result = msp2->type - msp1->type;
        }
      else 
        {
          result = msp1->qRange.min -  msp2->qRange.min;
        }
    }

  return result;
}

/* Compare the start position in the ref seq of two MSPs. Returns a negative value if a < b; zero
 * if a = b; positive value if a > b. Secondarily sorts by type in the order that types appear in 
 * the BlxMspType enum. Note that this sorts first by strand. */
gint compareFuncMspPos(gconstpointer a, gconstpointer b)
{
  gint result = 0;

  const MSP* const msp1 = (const MSP*)a;
  const MSP* const msp2 = (const MSP*)b;

  /* First, sort by strand */
  result = (int)msp1->qStrand - (int)msp2->qStrand;
  
  if (!result)
    result = compareMsps(msp1, msp2);

  return result;
}


/* Same as compareFuncMspPos but accepts pointers to MSP pointers (which is 
 * what the GArray of MSPs holds). Note that this does NOT sort first by strand,
 * unlike compareFuncMspPos. This is important for the detail-view filtering 
 * functions. */
gint compareFuncMspArray(gconstpointer a, gconstpointer b)
{
  const MSP* const msp1 = *((const MSP**)a);
  const MSP* const msp2 = *((const MSP**)b);

  return compareMsps(msp1, msp2);
}


/* returns true if the given msp should be output when piping features to dotter */
static gboolean outputMsp(const MSP* const msp, IntRange *range1, IntRange *range2)
{
  return ((msp->type == BLXMSP_FS_SEG || mspIsExon(msp) || mspIsIntron(msp) || mspIsBlastMatch(msp)) && 
          (rangesOverlap(&msp->qRange, range1) || rangesOverlap(&msp->qRange, range2))
         );
}


/* Counts the number of msps that should be output when piping to dotter */
static int countMspsToOutput(const BlxSequence* const blxSeq, IntRange *range1, IntRange *range2)
{
  int numMsps = 0;
  
  GList *mspItem = blxSeq->mspList;
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP* const msp = (const MSP*)(mspItem->data);
      
      if (outputMsp(msp, range1, range2))
        ++numMsps;
    }
  
  return numMsps;
}


/* write data from the given blxsequence to the given output pipe. if the ranges
 * are given, only outputs blxsequences that overlap either range */
void writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *range1, IntRange *range2)
{
  gboolean outputSeq = (blxSeq && (blxSeq->type == BLXSEQUENCE_TRANSCRIPT || blxSeq->type == BLXSEQUENCE_MATCH));
  
  int numMsps = countMspsToOutput(blxSeq, range1, range2);
  outputSeq &= numMsps > 0; /* only output the sequence if it has some valid msps */
  
  if (outputSeq)
    {
      fprintf(pipe, "%d %d %d",
              blxSeq->type,
              blxSeq->strand,
              numMsps); /* output number of msps so we know how many to read in */
      
      stringProtect(pipe, blxSequenceGetName(blxSeq));
      stringProtect(pipe, blxSeq->idTag);
      
      fputc('\n', pipe);
      
      /* now output the msps */
      GList *mspItem = blxSeq->mspList;
      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);
          
          if (outputMsp(msp, range1, range2))
            {
              writeMspToOutput(pipe, msp);
            }
        }
    }
}


/* increment the given char pointer until it is not whitespace */
static void nextChar(char **curChar)
{
  while (**curChar == ' ')
    ++(*curChar);
}


/* read in a blxsequence from the given text. */
BlxSequence* readBlxSequenceFromText(char *text, int *numMsps)
{
  DEBUG_ENTER("readBlxSequenceFromText(text=%s)", text);

  char *curChar = text;
  
  GError *error = NULL;
  BlxSequence *blxSeq = createEmptyBlxSequence();
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      DEBUG_EXIT("readBlxSequenceFromText returning NULL");
      return NULL;
    }
  
  nextChar(&curChar);
  
  blxSeq->type = (BlxSequenceType)strtol(curChar, &curChar, 10);
  nextChar(&curChar);

  blxSeq->strand = (BlxStrand)strtol(curChar, &curChar, 10);
  nextChar(&curChar);

  *numMsps = strtol(curChar, &curChar, 10);
  nextChar(&curChar);
  
  char *fullName = stringUnprotect(&curChar, NULL);
  blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQNAME, fullName);
  g_free(fullName);

  blxSeq->idTag = stringUnprotect(&curChar, NULL);
  
  DEBUG_EXIT("readBlxSequenceFromText returning numMsps=%d", *numMsps);
  return blxSeq;
}


/* write data from the given msp to the given output pipe. */
void writeMspToOutput(FILE *pipe, const MSP* const msp)
{
  fprintf(pipe, "%d %f %f %d %d %d %d %d %d %d", 
          msp->type,
          msp->score, 
          msp->id,
          msp->phase,
//          msp->fsColor, 
          msp->qRange.min,
          msp->qRange.max,
          msp->sRange.min,
          msp->sRange.max,
          msp->qStrand,
          msp->qFrame);
  stringProtect(pipe, msp->qname);
  stringProtect(pipe, msp->sname);
  stringProtect(pipe, msp->desc);
  fputc('\n', pipe);
}


/* Read in an msp from the given text. the text is in the format as written by writeMspToFile  */
void readMspFromText(MSP *msp, char *text)
{
  //DEBUG_ENTER("readMspFromText(text=%s)", text);

  char *curChar = text;

  nextChar(&curChar);
  msp->type = (BlxMspType)strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->score = g_strtod(curChar, &curChar);

  nextChar(&curChar);
  msp->id = g_strtod(curChar, &curChar);

  nextChar(&curChar);
  msp->phase = strtol(curChar, &curChar, 10);

//  nextChar(&curChar);
//  msp->fsColor = strtol(curChar, &curChar, 10);

  /* ref seq range */
  nextChar(&curChar);
  const int qStart = strtol(curChar, &curChar, 10); 
  nextChar(&curChar);
  const int qEnd = strtol(curChar, &curChar, 10); 
  intrangeSetValues(&msp->qRange, qStart, qEnd);

  /* match seq range */
  nextChar(&curChar);
  const int sStart = strtol(curChar, &curChar, 10); 
  nextChar(&curChar);
  const int sEnd = strtol(curChar, &curChar, 10); 
  intrangeSetValues(&msp->sRange, sStart, sEnd);
  
  nextChar(&curChar);
  msp->qStrand = (BlxStrand)strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->qFrame = strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->qname = stringUnprotect(&curChar, NULL);
  msp->sname = stringUnprotect(&curChar, NULL);
  msp->desc = stringUnprotect(&curChar, NULL);

  //DEBUG_EXIT("readMspFromText returning");
}


/* Insert the given MSP into the given list */
static void insertMsp(MSP *msp, MSP **mspList, MSP **lastMsp)
{
  /* Add it to the list */
  if (!*mspList) 
    {
      /* Nothing in list yet: make this the first entry */
      *mspList = msp;
    }
  
  if (*lastMsp)
    {
      /* Tag it on to the end of the list */
      (*lastMsp)->next = msp;
    }
  
  /* Make the 'lastMsp' pointer point to the new end of the list */
  *lastMsp = msp;
}


/* Allocate memory for an MSP and initialise all its fields to relevant 'empty' values.
 * Add it into the given MSP list and make the end pointer ('lastMsp') point to the new
 * end of the list. NB the msplist is obsolete and is being replaced by featureLists, but we don't
 * add it to a feature list yet because we don't know its type. Returns a pointer to the newly-created MSP */
MSP* createEmptyMsp(MSP **lastMsp, MSP **mspList)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  int i = 0;
  for ( ; i < BLXMODEL_NUM_MODELS; ++i)
    msp->treePaths[i] = NULL;

  msp->next = NULL;
  msp->childMsps = NULL;
  
  msp->type = BLXMSP_INVALID;
  msp->score = 0.0;
  msp->id = 0.0;
  msp->phase = 0;
  
  msp->qname = NULL;
  msp->qFrame = UNSET_INT;
  
  msp->qRange.min = 0;
  msp->qRange.max = 0;
  msp->displayRange.min = 0;
  msp->displayRange.max = 0;
  msp->fullRange.min = 0;
  msp->fullRange.max = 0;
  
  msp->sSequence = NULL;
  msp->sname = msp->sname_orig = NULL;
  
  msp->sRange.min = 0;
  msp->sRange.max = 0;
  msp->fullSRange.min = 0;
  msp->fullSRange.max = 0;
  
  msp->desc = NULL;
  
  msp->style = NULL;
  
  msp->fs = NULL;
  msp->fsColor = 0;
  msp->fsShape = BLXCURVE_BADSHAPE;
  
  msp->xy = NULL;
  msp->gaps = NULL;
  
  insertMsp(msp, mspList, lastMsp);
  
  return msp;
}


static void freeStringPointer(char **ptr)
{
  if (*ptr)
    {
      g_free(*ptr);
      *ptr = NULL;
    }
}


/* Destroy all of the MSPs in the given list */
void destroyMspList(MSP **mspList)
{
  /* Free the allocated sequences and names */
  MSP *msp = NULL;
  for (msp = *mspList; msp; msp = msp->next)
    {
      destroyMspData(msp);
    }
  
  /* Now free the MSPs themselves. */
  MSP *fmsp = NULL;
  for (msp = *mspList; msp; )
    {
      fmsp = msp;
      msp = msp->next;
      g_free(fmsp);
    }
  
  *mspList = NULL;
  
  return ;
}


/* Free all of the memory used by an MSP */
void destroyMspData(MSP *msp)
{
  freeStringPointer(&msp->qname);
  freeStringPointer(&msp->sname);
  freeStringPointer(&msp->sname_orig);
  freeStringPointer(&msp->desc);
  
  if (msp->gaps)
    {
      /* free the child msp list */
      if (msp->childMsps)
        {
          g_list_free(msp->childMsps);
          msp->childMsps = NULL;
        }

      /* free memory allocated for the gap ranges */
      GSList *item = msp->gaps;
      for ( ; item; item = item->next)
        g_free(item->data);
      
      g_slist_free(msp->gaps);
      msp->gaps = NULL;
    }
  
  if (msp->xy)
    {
      g_array_free(msp->xy, TRUE);
      msp->xy = NULL;
    }
}


/* Allocates memory for an MSP and initialise all its fields to the given values.
 * For backwards compatibility, adds it into the given MSP list and makes the end pointer
 * ('lastMsp') point to the new end of the list. We will hopefully get rid of mspList eventually
 * and replace it by featureLists. The new msp is added to the relevant list in the featureLists
 * array according to its type. Returns a pointer to the newly-created MSP. Also creates a BlxSequence
 * for this MSP's sequence name (or adds the MSP to the existing one, if it exists already), 
 * and adds that BlxSequence to the given seqList. Takes ownership of 'sequence'. */
MSP* createNewMsp(GArray* featureLists[],
                  MSP **lastMsp, 
                  MSP **mspList,
                  GList **seqList,
                  GList *columnList,
                  const BlxMspType mspType,
                  BlxDataType *dataType,
		  const char *source,
                  const gdouble score,
                  const gdouble percentId,
                  const int phase,
		  const char *idTag,
                  const char *qName,
                  const int qStart,
                  const int qEnd,
                  const BlxStrand qStrand,
                  const int qFrame,
                  const char *sName,
                  const char *sName_orig,
                  const int sStart,
                  const int sEnd,
                  BlxStrand sStrand,
                  char *sequence,
                  const GQuark filename,
                  GHashTable *lookupTable, 
                  GError **error)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = mspType;
  msp->score = score; 
  msp->id = percentId; 
  msp->phase = phase;
  msp->filename = filename;
  
  msp->qname = qName ? g_strdup(qName) : NULL;
  
  msp->qFrame = qFrame;
  msp->qStrand = qStrand;
  
  msp->sname = sName ? g_strdup(sName) : NULL;
  msp->sname_orig = sName_orig ? g_strdup(sName_orig) : NULL;

  
  intrangeSetValues(&msp->qRange, qStart, qEnd);  
  intrangeSetValues(&msp->sRange, sStart, sEnd);

  /* For exons and introns, the s strand is not applicable. We always want the exon
   * to be in the same direction as the ref sequence, so set the match seq strand to be 
   * the same as the ref seq strand */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      sStrand = qStrand;
    }
  
  /* Add it to the relevant feature list. */
  featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);

  /* For matches, exons and introns, add a new (or add to an existing) BlxSequence */
  if (typeIsExon(mspType) || typeIsIntron(mspType) || 
      typeIsMatch(mspType) || 
      typeIsVariation(mspType) || typeIsRegion(mspType))
    {
      addBlxSequences(msp->sname, msp->sname_orig, idTag, sStrand, dataType, source, 
                      featureLists, lastMsp, mspList, seqList,
                      columnList, sequence, msp, lookupTable, error);
    }

  if (error && *error)
    {
      prefixError(*error, "Error creating MSP (ref seq='%s' [%d - %d], match seq = '%s' [%d - %d]). ",
                  qName, qStart, qEnd, sName, sStart, sEnd);
    }
  
  return msp;
}


/* Make a copy of an MSP. If addToParent is true it also copies the pointer to the parent
 * BlxSequence and adds the new msp to that BlxSequence's mspList. Otherwise it's parent is null. */
MSP* copyMsp(const MSP* const src,
             GArray* featureLists[],             
             MSP **lastMsp, 
             MSP **mspList,
             const gboolean addToParent)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = src->type;
  msp->score = UNSET_INT; 
  msp->id = UNSET_INT; 
  msp->phase = src->phase;
  
  msp->qname = src->qname ? g_strdup(src->qname) : NULL;
  
  msp->qFrame = src->qFrame;
  msp->qStrand = src->qStrand;
  
  msp->sname = src->sname ? g_strdup(src->sname) : NULL;
  msp->sname_orig = src->sname_orig ? g_strdup(src->sname_orig) : NULL;
  
  intrangeSetValues(&msp->qRange, src->qRange.min, src->qRange.max);  
  intrangeSetValues(&msp->sRange, src->sRange.min, src->sRange.max);
  
  /* For matches, exons and introns, add (or add to if already exists) a BlxSequence */
  if (addToParent && src->sSequence)
    {
      src->sSequence->mspList = g_list_insert_sorted(src->sSequence->mspList, msp, compareFuncMspPos);
      msp->sSequence = src->sSequence;
    }

  /* Add it to the relevant feature list. */
  featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);
  
  return msp;
}



/* Set the given child list in the given exon. Takes ownership of the child list
 * (and frees it if exon is null).
 * Exons and UTRs don't have phase, but we want to display them in the same reading frame
 * as the CDS in the same exon, if there is one; this function copies it to its siblings. */
static void setExonChildList(MSP *exon, 
                             GList *childList,
                             GArray* featureLists[], 
                             MSP **lastMsp,
                             MSP **mspList,
                             GList **seqList,
                             GList *columnList,
                             GHashTable *lookupTable,
                             MSP **spanningCds)
{
  if (!exon)
    {
      g_list_free(childList);
      return;
    }

  /* Take ownership of the child list */
  if (exon->childMsps)
    g_list_free(exon->childMsps);
  
  exon->childMsps = childList;
  
  /* Loop through and see if there's a CDS */
  GError *tmpError = NULL;
  int frame = UNSET_INT;
  int phase = UNSET_INT;
  gboolean found = FALSE;
  
  GList *childItem = exon->childMsps;
  for ( ; childItem; childItem = childItem->next)
    {
      MSP *cds = (MSP*)(childItem->data);
      
      if (cds->type == BLXMSP_CDS)
        {
          frame = cds->qFrame;
          phase = cds->phase;
          found = TRUE;

          /* Hack to support invalid gff from zmap where a single cds object spans the entire
           * length rather than having separate cds objects for each exon. Trim the cds to the
           * exon and "save" the cds to check against later exons. We should only have one cds in
           * this case (we should verify this but don't at the moment). */
          if (rangesOverlap(&cds->qRange, &exon->qRange) &&
              (cds->qRange.min < exon->qRange.min || cds->qRange.max > exon->qRange.max))
            {
              /* Replace the original cds with a new one truncated to this exon */
              int start = max(cds->qRange.min, exon->qRange.min);
              int end = min(cds->qRange.max, exon->qRange.max);
              MSP *newCds = createMissingMsp(BLXMSP_CDS, start, end, cds->qname, cds->qFrame, cds->style, cds->sSequence, 
                                             featureLists, lastMsp, mspList, seqList, columnList, lookupTable, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);

              *spanningCds = cds;
              exon->childMsps = g_list_remove(exon->childMsps, cds);
              exon->childMsps = g_list_insert_sorted(exon->childMsps, newCds, compareFuncMspPos);
            }

          break;
        }
    }
  
  if (found)
    {
      /* Update the exon */
      exon->qFrame = frame;
      exon->phase = phase;
      
      /* Loop through and update the other msps (i.e. exon and UTR get the same frame/phase info
         as the CDS) */
      for (childItem = exon->childMsps; childItem; childItem = childItem->next)
        {
          MSP *msp = (MSP*)(childItem->data);
          msp->qFrame = frame;
          msp->phase = phase;
        }
    }
}


/* Utility to create the msp for the createMissingExonCdsUtr function. */
static MSP* createMissingMsp(const BlxMspType newType,
                             const int newStart,
                             const int newEnd,
                             const char *qname,
                             const int newFrame,
                             BlxStyle *newStyle,
                             BlxSequence *blxSeq, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             GList **seqList,
                             GList *columnList,
                             GHashTable *lookupTable,
                             GError **error)
{
  MSP *result = NULL;
  
  if (newType != BLXMSP_INVALID)
    {
      /* Create the new exon/cds/utr */
      DEBUG_OUT("Creating MSP for transcript '%s' of type %d.\n", blxSequenceGetName(blxSeq), newType);
      
      GError *tmpError = NULL;
      
      result = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, newType, NULL, blxSequenceGetSource(blxSeq),
                            UNSET_INT, UNSET_INT, UNSET_INT, blxSeq->idTag,
                            qname, newStart, newEnd, blxSeq->strand, newFrame,
                            blxSequenceGetName(blxSeq), blxSequenceGetName(blxSeq),
                            UNSET_INT, UNSET_INT, blxSeq->strand, NULL,
                            0, lookupTable, &tmpError);
      
      result->style = newStyle;
      
      if (tmpError)
        {
          prefixError(tmpError, "Error constructing missing exon/cds/utr [type='%d']", newType);
          g_propagate_error(error, tmpError);
        }
    }
  
  return result;
}


/* We have an exon and one or more child CDS/UTRs. Check whether there is a gap
 * at the startor end of the exon and, if so, construct a CDS/UTR to fill it. */
static void createMissingCdsUtr(MSP *exon,
                                GList **childList,
                                BlxSequence *blxSeq, 
                                GArray* featureLists[], 
                                MSP **lastMsp, 
                                MSP **mspList, 
                                GList **seqList, 
                                GList *columnList,
                                GHashTable *lookupTable,
                                GError **error)
{
  MSP *startMsp = (MSP*)(g_list_first(*childList)->data);
  MSP *endMsp = (MSP*)(g_list_last(*childList)->data);
  
  if (exon->qRange.min < startMsp->qRange.min)
    {
      const BlxMspType type = (startMsp->type == BLXMSP_CDS ? BLXMSP_UTR : BLXMSP_CDS);
      MSP *result = createMissingMsp(type, exon->qRange.min, startMsp->qRange.min - 1, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
      *childList = g_list_append(*childList, result);
    }

  if (exon->qRange.max > endMsp->qRange.max)
    {
      const BlxMspType type = (endMsp->type == BLXMSP_CDS ? BLXMSP_UTR : BLXMSP_CDS);
      MSP *result = createMissingMsp(type, endMsp->qRange.max + 1, exon->qRange.max, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
      *childList = g_list_append(*childList, result);
    }
}


/* Create a UTR that spans the given exon, and add it to the given childList */
static void createMissingUtr(MSP *exon,
                             GList **childList,
                             BlxSequence *blxSeq, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             GList **seqList, 
                             GList *columnList,
                             GHashTable *lookupTable,
                             GError **error)
{
  MSP *result = createMissingMsp(BLXMSP_UTR, exon->qRange.min, exon->qRange.max, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
  *childList = g_list_append(*childList, result);
}


/* Create an exon that spans the given child CDSs/UTRs. Must not be called with
 * an empty childList. */
static MSP* createMissingExon(GList *childList,
                              BlxSequence *blxSeq, 
                              GArray* featureLists[], 
                              MSP **lastMsp, 
                              MSP **mspList, 
                              GList **seqList, 
                              GList *columnList,
                              GHashTable *lookupTable,
                              GError **error)
{
  /* Get the max and min extent of the children. They should be in order
   * of increasing coords and should not overlap etc. */
  MSP *startMsp = (MSP*)(g_list_first(childList)->data);
  MSP *endMsp = (MSP*)(g_list_last(childList)->data);
  
  MSP *result = createMissingMsp(BLXMSP_EXON, startMsp->qRange.min, endMsp->qRange.max, startMsp->qname, startMsp->qFrame, startMsp->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
  return result;
}


/* Utility used by constructExonData to create a missing exon/cds/utr given 
 * two others out of the three - i.e. if we have an overlapping exon and cds we can
 * construct the corresponding utr. If created, the new msp is added to the given 
 * BlxSequence and the  MSP list. If a CDS is given and no UTR exists, assume the exon
 * spans the entire CDS (and similarly if a UTR is given but no CDS exists) */
static void createMissingExonCdsUtr(MSP **exon, 
                                    GList **childList,
                                    BlxSequence *blxSeq, 
                                    GArray* featureLists[], 
                                    MSP **lastMsp, 
                                    MSP **mspList, 
                                    GList **seqList, 
                                    GList *columnList,
                                    GHashTable *lookupTable,
                                    GError **error)
{
  if (*exon && g_list_length(*childList) > 0)
    {
      /* We have an exon and one or more CDS/UTRs. Check if there any other
       * child CDS/UTRs that we need to create. */
      createMissingCdsUtr(*exon, childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
    }
  else if (*exon)
    {
      /* We have an exon but no children. Assume the exon is a single UTR */
      createMissingUtr(*exon, childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
    }
  else if (g_list_length(*childList) > 0)
    {
      /* We have children but no exon; create the exon from the children */
      *exon = createMissingExon(*childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, error);
    }
}


/* Construct any missing exon data, i.e.
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructExonData(BlxSequence *blxSeq, 
                              GArray* featureLists[], 
                              MSP **lastMsp,
                              MSP **mspList,
                              GList **seqList,
                              GList *columnList,
                              GHashTable *lookupTable)
{
  GError *tmpError = NULL;
  
  const MSP *prevMsp = NULL;
  const MSP *prevExon = NULL;
  
  MSP *curExon = NULL;          /* the current exon we're looking at */
  GList *curChildMsps = NULL;   /* the child CDS/UTRs of the current exon */
  MSP *spanningCds = NULL;      /* hack to support invalid GFF used by zmap where a single CDS
                                 * spanning the entire range is given, rather than a separate CDS
                                 * feature for each exon */
  
  /* Loop through all MSPs on this sequence (which must be sorted by position on the
   * ref seq - createNewMsp automatically sorts them for us) and create any missing 
   * exon/cds/utr/introns. */
  GList *mspItem = blxSeq->mspList;
  gboolean finished = FALSE;
  
  while (!finished)
    {
      MSP *msp = mspItem ? (MSP*)(mspItem->data) : NULL;
      
      /* Only consider exons and introns */
      if (mspIsExon(msp) || mspIsIntron(msp) || !msp)
        {
          /* See if there was a gap between this exon and the previous one. There's a gap if
           * we have two exons with space between them, or if we're at the first or last exon 
           * and there's a gap to the end of the transcript. */
          gboolean foundGap = FALSE;

          if (msp && mspIsIntron(msp))
            {
              /* An intron IS a gap! */
              foundGap = TRUE;
            }
          else if (msp && prevMsp)
            {
              /* We have a previous msp; there's a gap if the previous msp is an exon
               * and if there is a gap between it and the current exon */
              foundGap = mspIsExon(prevMsp) && !rangesOverlap(&prevMsp->qRange, &msp->qRange) && !rangesAdjacent(&prevMsp->qRange, &msp->qRange);
            }
          else if (msp)
            {
              /* No previous msp; check if we're at the first exon and if so whether
               * there's a gap between it and the start of the transcript */
              foundGap = blxSequenceGetStart(blxSeq, blxSeq->strand) < msp->qRange.min;
            }
          
          if (foundGap || msp == NULL)
            {
              /* We've found a gap between exons, or reached the end. First, see if the current exon/cds or utr
               * is missing and construct it if possible. Also do this if we're at the last MSP. */
              createMissingExonCdsUtr(&curExon, &curChildMsps, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
              
              IntRange newRange = {UNSET_INT, UNSET_INT};
              
              if (prevExon && curExon && !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron to span the gap */
                  newRange.min = prevExon->qRange.max + 1;
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (!prevExon && curExon && blxSequenceGetStart(blxSeq, blxSeq->strand) < curExon->qRange.min && 
		       !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the start */
                  newRange.min = blxSequenceGetStart(blxSeq, blxSeq->strand);
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (msp == NULL && curExon && blxSequenceGetEnd(blxSeq, blxSeq->strand) > curExon->qRange.max &&
		       !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the end */
                  newRange.min = curExon->qRange.max + 1;
                  newRange.max = blxSequenceGetEnd(blxSeq, blxSeq->strand);
                }
              
              if (curExon && newRange.min != UNSET_INT && newRange.max != UNSET_INT)
                {
                  createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_INTRON, NULL, blxSequenceGetSource(blxSeq), 
                               curExon->score, curExon->id, 0, blxSeq->idTag, 
                               curExon->qname, newRange.min, newRange.max, blxSeq->strand, curExon->qFrame, 
                               blxSequenceGetName(blxSeq), blxSequenceGetName(blxSeq),
                               UNSET_INT, UNSET_INT, blxSeq->strand, NULL, 
                               0, lookupTable, &tmpError);
                  
                  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
                }
              
              /* We're done with this exon, so set the exon's list of child msps
               * and reset the pointers */
              setExonChildList(curExon, curChildMsps, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, &spanningCds);

              prevExon = curExon;
              curExon = NULL;
              curChildMsps = NULL;
            }
          
          if (msp && msp->type == BLXMSP_EXON)
            {
              curExon = msp;

              /* If there's a spanning cds, add it to the child list for each exon */
              if (spanningCds)
                curChildMsps = g_list_append(curChildMsps, spanningCds);
            }
          else if (msp && (msp->type == BLXMSP_CDS || msp->type == BLXMSP_UTR))
            {
              curChildMsps = g_list_append(curChildMsps, msp);
            }
          
          /* Remember the last MSP we saw */
          prevMsp = msp;
          
          /* Proceed to the next MSP. We allow an extra loop with a NULL mspItem at the end, and then finish. */
          if (mspItem)
            mspItem = mspItem->next;
          else
            finished = TRUE;
        } 
      else
        {
          /* Something that's not an exon/intron. Skip this BlxSequence. */
          finished = TRUE;
        }
    }

  if (spanningCds)
    destroyMspFull(spanningCds, blxSeq, featureLists, lastMsp, mspList);
}


/* Construct any missing transcript data, i.e.
 *   - if we have multiple CDSs, copy the transcript so we can show each variant;
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructTranscriptData(GArray* featureLists[], 
                                    MSP **lastMsp,
                                    MSP **mspList,
                                    GList **seqList,
                                    GList *columnList,
                                    GHashTable *lookupTable,
                                    GError **error)
{
  GError *tmpError = NULL;

  /* Loop through all transcripts  */
  GList *seqItem = *seqList;

  while (seqItem)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);

      if (blxSeq->type != BLXSEQUENCE_TRANSCRIPT)
        {
          seqItem = seqItem->next;
          continue;
        }

      /* Get a list of all CDS names in the transcript */
      GList *cdsList = blxSequenceConstructCdsList(blxSeq);

      const int numVariants = g_list_length(cdsList);

      if (numVariants <= 1)
        {
          /* Only one variant so process it as it is */
          findSequenceExtents(blxSeq);
          constructExonData(blxSeq, featureLists, lastMsp, mspList, seqList, columnList, lookupTable);

          seqItem = seqItem->next;
        }
      else
        {
          /* More than one variant: create copies of the transcript for each variant */
          GList *cdsItem = cdsList;
  
          for ( ; cdsItem && !tmpError; cdsItem = cdsItem->next)
            {
              GQuark cdsQuark = GPOINTER_TO_INT(cdsItem->data);

              /* The copy is added to the end of seqList so it will be processed later as we continue
               * to loop through the list */
              copyBlxSequenceNamedCds(blxSeq, cdsQuark, featureLists, lastMsp, mspList, seqList, columnList, lookupTable, &tmpError);
            }

          /* We have made copies for each CDS but we don't need the original transcript so delete
           * it. This removes the current seqItem from seqList so we have to increment the pointer
           * before we delete it. */
          seqItem = seqItem->next;
          destroyBlxSequenceFull(blxSeq, featureLists, lastMsp, mspList, seqList);
        }

      g_list_free(cdsList);
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Adjust the MSP's q coords (as parsed from the input file) by the given offset i.e. to convert
 * them to real coords */
static void adjustMspCoordsByOffset(MSP *msp, const int offset)
{
  if (offset != 0)
    {
      /* Convert the input coords (which are 1-based within the ref sequence section
       * that we're dealing with) to "real" coords (i.e. coords that the user will see). */
      msp->qRange.min += offset;
      msp->qRange.max += offset;
      
      /* Gap coords are also 1-based, so convert those too */
      GSList *rangeItem = msp->gaps;
      
      for ( ; rangeItem; rangeItem = rangeItem->next)
        {
          CoordRange *curRange = (CoordRange*)(rangeItem->data);
          curRange->qStart += offset;
          curRange->qEnd += offset;
        }
    }
}


/* Find the first/last base in an entire sequence, if not already set. For transcripts, the
 * range should already be set from the parent transcript item. For matches, get the start/end
 * of the first/last MSP in the sequence */
static void findSequenceExtents(BlxSequence *blxSeq)
{
  blxSeq->qRangeFwd.min = findMspListQExtent(blxSeq->mspList, TRUE, BLXSTRAND_FORWARD);
  blxSeq->qRangeFwd.max = findMspListQExtent(blxSeq->mspList, FALSE, BLXSTRAND_FORWARD);
  blxSeq->qRangeRev.min = findMspListQExtent(blxSeq->mspList, TRUE, BLXSTRAND_REVERSE);
  blxSeq->qRangeRev.max = findMspListQExtent(blxSeq->mspList, FALSE, BLXSTRAND_REVERSE);
}


/* Get the offset required from the given base to give the coord that is the first base in 
 * the first codon of reading frame 1 (or the last codon in reading frame 3 for the reverse strand) */
static int getOffsetToCodonStart(const int coord, const int numFrames, const BlxStrand strand)
{
  int offset = 0;
  
  if (strand == BLXSTRAND_FORWARD)
    {
      /* If the sequence is 
       *     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ...
       * then the calculated base below will be
       *     0, 1, 2, 0, 1, 2, 0, 1, 2, 0,  1, ...
       * so we have to offset the coords by the following values to get to a base1 coord
       *     0, 2, 1, 0, 2, 1, 0, 2, 1, 0,  2, ...
       */
      
      int base = ((coord - 1) % numFrames); /* 0, 1 or 2 */
      offset = numFrames - base;            /* 3, 2 or 1 */
      
      if (offset >= numFrames)
	offset -= numFrames;                /* 0, 2 or 1 */
    }
  else
    {
      /* I'm not sure if there is a convention for where the reading frame starts in the reverse
       * strand, so I've made up my own convention. It essentially means that we want the first 
       * coord in the reversed sequence to be the last base in the last reading frame, i.e. base 3 in frame 3.
       *
       * If the forward strand coords are
       *             1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 
       * then the have the following base numbers for each frame
       * frame 1:    1, 2, 3, 1, 2, 3, 1, 2, 3, 1,  2
       * frame 2:    3, 1, 2, 3, 1, 2, 3, 1, 2, 3,  1
       * frame 3:    2, 3, 1, 2, 3, 1, 2, 3, 1, 2,  3 (I've included up to coord 11 so we end at base 3 in frame 3)
       *
       * Looking at frame 1 in reverse we have
       *             11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
       *              2,  1, 3, 2, 1, 3, 2, 1, 3, 2, 1
       *
       * Base 3 in frame 3 is base 2 in frame 1, so this is where we want to start. The offsets we need
       * to get to a frame1-base2 coord are:
       *              0,  2, 1, 0, 2, 1, 0, 2, 1, 0, 2
       */
      
      int base = (coord + 1) % numFrames;   /* 0, 2 or 1 */
      offset = base;
    }
  
  return offset;
}


/* Calculate the reference sequence reading frame that the given MSP belongs to, if not
 * already set. Requires either the phase to be set, or the frame to already be set; otherwise
 * assumes a phase of 0 and gives a warning. */
static void calcReadingFrame(MSP *msp, const BlxSeqType seqType, const int numFrames, const IntRange* const refSeqRange)
{
  /* For matches and exons, calculate frame if the phase is known, because the old code that 
   * used to pass the reading frame in exblx files seemed to occasionally pass an incorrect reading frame. */
  if (!mspIsIntron(msp))
    {
      /* Get the first coord of the first complete codon. This is the start coord (5' end) of the match
       * plus (or minus) the phase (if non-zero), which initially gets stored in the qFrame field in the MSP... */
      const int direction = (msp->qStrand == BLXSTRAND_FORWARD ? 1 : -1);
      const int coord = mspGetQStart(msp) + (direction * msp->phase);
      
      /* Find the reading frame that this coord belongs in. This is the same as the base number within
       * reading frame 1. */
      int frame = UNSET_INT;
      const gboolean invertCoords = (mspGetRefStrand(msp) == BLXSTRAND_REVERSE);
      
      convertDnaIdxToDisplayIdx(coord, seqType, 1, numFrames, invertCoords, refSeqRange, &frame);
      
      if (frame > 0)
	{
	  if (msp->qFrame > 0 && mspIsExon(msp))
	    {
	      /* Frame is already set. This means we have the old exblx file format, which does not provide phase for exons
	       * so we must use the reading frame that was passed in the file. However, this used to assume that the first base
	       * in the ref seq was base 1 in reading frame 1, which is not always the case (we now calculate the correct reading frame
	       * based on mod3 of the coord). Therefore, we must offset the given frame if the first coord is not base1 in frame1. */
	      const int startCoord = (msp->qStrand == BLXSTRAND_REVERSE ? refSeqRange->max : refSeqRange->min);
	      const int offset = getOffsetToCodonStart(startCoord, numFrames, msp->qStrand);
	      msp->qFrame += offset;
	      
	      if (msp->qFrame > numFrames)
		{
		  msp->qFrame -= numFrames;
		}
	      
	      if (msp->qFrame != frame && seqType == BLXSEQ_PEPTIDE)
		{
		  g_warning("MSP '%s' (q=%d-%d; s=%d-%d) has reading frame '%d' but calculated frame was '%d'\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max, msp->qFrame, frame);
		}
	    }
	  else
	    {
	      /* We either have the new file format which provides phase or it's not an exon so phase 
	       * is not applicable, so we can trust the calculated value. */
	      msp->qFrame = frame;
	    }
	}
      
      if (msp->qFrame == UNSET_INT)
	{
	  g_warning("Reading frame could not be calculated for MSP '%s' (q=%d-%d; s=%d-%d) - setting to 1.\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
	  msp->qFrame = 1;
	}
    }
}


/* Should be called after all parsed data has been added to a BlxSequence. Calculates summary
 * data and the introns etc. */
void finaliseBlxSequences(GArray* featureLists[], 
			  MSP **mspList, 
			  GList **seqList, 
                          GList *columnList,
			  const int offset,
			  const BlxSeqType seqType, 
			  const int numFrames,
			  const IntRange* const refSeqRange,
			  const gboolean calcFrame,
                          GHashTable *lookupTable)
{
  GError *tmpError = NULL;

  /* Loop through all MSPs and adjust their coords by the offest, then calculate their reading
   * frame. Also find the last MSP in the list. */
  MSP *msp = *mspList;
  MSP *lastMsp = msp;

  while (msp)
    {
      adjustMspCoordsByOffset(msp, offset);
    
      if (calcFrame)
	calcReadingFrame(msp, seqType, numFrames, refSeqRange);
    
      msp = msp->next;
      if (msp)
        lastMsp = msp;
    }
  
  /* Loop through all BlxSequences */
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      /* So far we only have the forward strand version of each sequence. We must complement any 
       * that need the reverse strand */
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);

      if (blxSeq && 
          blxSeq->strand == BLXSTRAND_REVERSE && 
          blxSequenceGetFlag(blxSeq, MSPFLAG_STRAND_SPECIFIC) &&
          blxSequenceGetFlag(blxSeq, MSPFLAG_SHOW_REVERSE_STRAND) &&
          blxSequenceGetSequence(blxSeq))
        {
          char *sequence = g_strdup(blxSequenceGetSequence(blxSeq));
          blxComplement(sequence);
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQUENCE, sequence);
        }
    }

  /* We need to construct any missing transcript data */
  constructTranscriptData(featureLists, &lastMsp, mspList, seqList, columnList, lookupTable, &tmpError);
  prefixError(tmpError, "Error constructing transcript data: ");
  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);

  /* We need to clear GFF ID flags once we've finished importing a file so that they don't clash
   * with IDs from any other files that the user may import later. */
  for (seqItem = *seqList; seqItem; seqItem = seqItem->next)
    blxSequenceClearGFFIds((BlxSequence*)(seqItem->data));

  /* Sort msp arrays by start coord (only applicable to msp types that
   * appear in the detail-view because the order is only applicable when
   * filtering detail-view rows) */
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      if (typeShownInDetailView((BlxMspType)typeId))
        g_array_sort(featureLists[typeId], compareFuncMspArray);
    }
}


/* Given a list of msps, find the min or max q coords (according to the given flag)
 * for msps on the given q strand (or on both strands if strand is BLXSTRAND_NONE. */
int findMspListQExtent(GList *mspList, const gboolean findMin, const BlxStrand strand)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP* const msp = (const MSP*)(mspItem->data);
    
      if (msp->qStrand == strand || strand == BLXSTRAND_NONE)
	{      
	  if (first)
	    {
	      result = findMin ? msp->qRange.min : msp->qRange.max;
	      first = FALSE;
	    }
	  else if (findMin && msp->qRange.min < result)
	    {
	      result = msp->qRange.min;
	    }
	  else if (!findMin && msp->qRange.max > result)
	    {
	      result = msp->qRange.max;
	    }
	}
    }
  
  return result;
}


/* Given a list of msps, find the min or max s coords, according to the given flag. */
int findMspListSExtent(GList *mspList, const gboolean findMin)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP* const msp = (const MSP*)(mspItem->data);
      
      if (first)
	{
	  result = findMin ? msp->sRange.min : msp->sRange.max;
	  first = FALSE;
	}
      else if (findMin && msp->sRange.min < result)
	{
	  result = msp->sRange.min;
	}
      else if (!findMin && msp->sRange.max > result)
	{
	  result = msp->sRange.max;
	}
    }
  
  return result;
}


/* Get the default value for the given flag */
gboolean mspFlagGetDefault(const MspFlag flag)
{
  gboolean result = FALSE;
  
  if (flag > MSPFLAG_MIN && flag < MSPFLAG_NUM_FLAGS)
    {
      /* The defaults are populated when we create a BlxDataType,
       * so create a dummy one so that we can access them. This
       * is made a global so that we can also change the defaults. */
      if (!g_DefaultDataType)
        g_DefaultDataType = createBlxDataType();
      
      if (g_DefaultDataType)
        result = g_DefaultDataType->flags[flag];
      else
        g_critical("Program error: Failed to create default data-type struct. Default feature properties may be incorrect.\n");
    }
  else
    {
      g_critical("Program error: attempt to use unknown MSP flag '%d'\n", flag);
    }
  
  return result;
}


/* Set the default value for the given flag */
void mspFlagSetDefault(const MspFlag flag, const gboolean value)
{
  if (flag > MSPFLAG_MIN && flag < MSPFLAG_NUM_FLAGS)
    {
      /* The defaults are populated when we create a BlxDataType,
       * so create a dummy one so that we can access and set them. */
      if (!g_DefaultDataType)
        g_DefaultDataType = createBlxDataType();
      
      if (g_DefaultDataType)
        g_DefaultDataType->flags[flag] = value;
    }
  else
    {
      g_critical("Program error: attempt to use unknown MSP flag '%d'\n", flag);
    }
}


/* Returns the colinearity of the two msps */
ColinearityType mspIsColinear(const MSP* const msp1, const MSP* const msp2)
{
  ColinearityType result = COLINEAR_INVALID;

  if (msp1 && msp2 && msp1->qStrand == msp2->qStrand && 
      msp1->sSequence && msp2->sSequence && 
      msp1->sSequence->strand == msp2->sSequence->strand)
    {
      if (msp2->sRange.min < msp1->sRange.max)
        result = COLINEAR_NOT;
      else if (msp2->sRange.min == msp1->sRange.max + 1)
        result = COLINEAR_PERFECT;
      else
        result = COLINEAR_IMPERFECT;
    }

  return result;
}


/* Return the value of the given boolean flag */
gboolean dataTypeGetFlag(const BlxDataType* const dataType, const MspFlag flag)
{
  gboolean result = mspFlagGetDefault(flag);

  if (flag > MSPFLAG_MIN && flag < MSPFLAG_NUM_FLAGS)
    {  
      if (dataType)
        result = dataType->flags[flag];
    }
  else
    {
      g_critical("Program error: attempt to use unknown MSP flag '%d'\n", flag);
    }
  
  return result;
}


/* Return the value of the given boolean flag. */
gboolean blxSequenceGetFlag(const BlxSequence* const blxSeq, const MspFlag flag)
{
  gboolean result = mspFlagGetDefault(flag);
  
  if (blxSeq)
    result = dataTypeGetFlag(blxSeq->dataType, flag);
    
  return result;
}


/* Return the value of the given boolean flag */
gboolean mspGetFlag(const MSP* const msp, const MspFlag flag)
{
  gboolean result = mspFlagGetDefault(flag);
  
  if (msp)
    result = blxSequenceGetFlag(msp->sSequence, flag);
  
  return result;
}

/* Get the config key to use for the given flag */
const char* mspFlagGetConfigKey(const MspFlag flag)
{
  const char *result = NULL;
  
  if (flag > MSPFLAG_MIN && flag < MSPFLAG_NUM_FLAGS)
    {
      /* To make sure we don't access the array out of bounds, loop
       * through from the beginning checking that no values are "dummy".
       * The array is terminated with a "dummy" value, so if we find this
       * before we find our flag, there must be enum values that haven't
       * been added to the array, which is a programming error. */
      int i = MSPFLAG_MIN + 1;
      for ( ; i <= flag; ++i)
        {
          if (stringsEqual(g_MspFlagConfigKeys[i], "dummy", FALSE))
            {
              g_critical("Program error: MSP flag '%d' does not have a config key defined\n", flag);
              break;
            }
          else if (i == flag)
            {
              result = g_MspFlagConfigKeys[i];
              break;
            }          
        }
    }
  else
    {
      g_critical("Program error: Tried to use an unknown MSP flag '%d'\n", flag);
    }
  
  return result;
}


/* Determine the type of BlxSequence that a particular MSP type belongs to */
static BlxSequenceType getBlxSequenceTypeForMsp(const BlxMspType mspType)
{
  BlxSequenceType result = BLXSEQUENCE_UNSET;
  
  if (mspType == BLXMSP_MATCH)
    {
      result = BLXSEQUENCE_MATCH;
    }
  else if (typeIsExon(mspType) || typeIsIntron(mspType) || mspType == BLXMSP_TRANSCRIPT)
    {
      result = BLXSEQUENCE_TRANSCRIPT;
    }
  else if (typeIsVariation(mspType))
    {
      result = BLXSEQUENCE_VARIATION;
    }
  else if (mspType == BLXMSP_REGION)
    {
      result = BLXSEQUENCE_REGION;
    }

  return result;
}


/* Utility to get the unique key for a given text string and strand */
static GQuark getLookupKey(const char *text, const BlxStrand strand)
{
  char *keyStr = g_strdup_printf("%s%c", text, (strand == BLXSTRAND_FORWARD ? '+' : '-'));
  GQuark key = g_quark_from_string(keyStr);
  g_free(keyStr);
  
  return key;
}


/* Utility to find a blxsequence with the given name/id/strand in the 
 * given hash table. Returns null if it is not there. */
static BlxSequence* findBlxSequence(GHashTable *lookupTable,
                                    const char *name, 
                                    const char *idTag,
                                    const BlxStrand strand,
                                    const gboolean linkFeaturesByName)
{
  BlxSequence *result = NULL;

  if (idTag)
    {
      /* We compare on the id tag and strand, so combine these into a single
       * string and convert it to a quark for quicker comparisions.
       * (This is the key for the hash table.) */
      GQuark key = getLookupKey(idTag, strand);
      result = (BlxSequence*)g_hash_table_lookup(lookupTable, GINT_TO_POINTER(key));
    }

  if (!result && name && linkFeaturesByName)
    {
      /* No id tag is given but we are asked to link features with the same
       * name, so do the comparison using name and strand */
      GQuark key = getLookupKey(name, strand);
      result = (BlxSequence*)g_hash_table_lookup(lookupTable, GINT_TO_POINTER(key));
    }

  return result;
}


/* Sort two columns by index */
gint columnIdxCompareFunc(gconstpointer a, gconstpointer b)
{
  BlxColumnInfo *col1 = (BlxColumnInfo*)a;
  BlxColumnInfo *col2 = (BlxColumnInfo*)b;
  return col1->columnIdx - col2->columnIdx;
}

/* Sort two columns by ID */
static gint columnIdCompareFunc(gconstpointer a, gconstpointer b)
{
  BlxColumnInfo *col1 = (BlxColumnInfo*)a;
  BlxColumnInfo *col2 = (BlxColumnInfo*)b;
  return col1->columnId - col2->columnId;
}


/* Create a new sequence with the given name and/or ID tag */
static BlxSequence* createBlxSequence(const char *name,
                                      const char *idTag, 
                                      const BlxStrand strand,
                                      BlxDataType *dataType,
                                      const char *source,
                                      GList *columnList)
{
  BlxSequence *seq = createEmptyBlxSequence();

  /* Set the hard-coded values */
  seq->strand = strand;
  seq->dataType = dataType;

  if (idTag)
    seq->idTag = g_strdup(idTag);

  /* Allocate the list of values and initialise the value types
   * to the relevant type for that column */
  //const int numColumns = g_list_length(columnList);
  seq->values = g_array_new(FALSE, FALSE, sizeof(GValue));

  /* Sort the columns by column ID (NOT column index) so that we
   * can easily index on the ID when looking up values in the array */
  columnList = g_list_sort(columnList, columnIdCompareFunc);
  
  GList *item = columnList;
  for ( ; item; item = item->next)
    {
      /* Note that we index on columnId, not columnIdx. This is so that
       * we can quickly index on the columnId enum. The columnIdx represents
       * the display order, which is not relevant here. */
      BlxColumnInfo * columnInfo = (BlxColumnInfo*)(item->data);
      GType type = columnInfo->type;
      
      /* Bit of a hack; the sequence column type is pointer because we
       * pass a pointer to the msp to the cell renderer, but here we want to
       * store the actual sequence string */
      if (columnInfo->columnId == BLXCOL_SEQUENCE)
        type = G_TYPE_STRING;
      
      GValue value = {0};
      g_value_init(&value, type);

      g_array_append_val(seq->values, value);
    }

  /* Set the given name and source */
  blxSequenceSetValueFromString(seq, BLXCOL_SEQNAME, name ? name : idTag);
  blxSequenceSetValueFromString(seq, BLXCOL_SOURCE, source);

  /* Make sure we change back to the original sort order (sorted by index) */
  columnList = g_list_sort(columnList, columnIdxCompareFunc);

  return seq;
}


/* Wrapper for addBlxSequence to add multiple sequences. The idTag might be a comma-separated
 * list of parent IDs, in which case we need to add the msp to multiple BlxSequences (creating
 * those BlxSequences if they don't exist. */
static void addBlxSequences(const char *name, 
                            const char *name_orig, 
                            const char *idTag, 
                            BlxStrand strand,
                            BlxDataType *dataType,
                            const char *source,
                            GArray *featureLists[],
                            MSP **lastMsp,
                            MSP **mspList,
                            GList **seqList, 
                            GList *columnList,
                            char *sequence, 
                            MSP *msp_in, 
                            GHashTable *lookupTable,
                            GError **error)
{
  GError *tmpError = NULL;
  MSP *msp = msp_in;

  if (idTag)
    {
      /* For exons and introns, the ID tag we receive is the parent tag, and it may contain
       * multiple parent transcripts. In this case we need to add the msp to multiple parent
       * BlxSequences */
      char **tokens = g_strsplit_set(idTag, ",", -1);   /* -1 means do all tokens. */
      char **token = tokens;
      gboolean usedMsp = FALSE;
      
      while (token && *token && **token && !tmpError)
        {
          /* If we've already used the passed-in msp, then we need to make a copy of it to 
           * add to the next BlxSequence (because the msp points to its BlxSequence so can't
           * be added to multiple BlxSequences, at least at the moment) */
          if (usedMsp)
            msp = copyMsp(msp, featureLists, lastMsp, mspList, FALSE);

          if (!tmpError)
            addBlxSequence(msp->sname, msp->sname_orig, *token, strand,
                           dataType, source, seqList, columnList, sequence, msp, lookupTable, &tmpError);

          usedMsp = TRUE;
          ++token;
        }

      g_strfreev(tokens);
    }
  else
    {
      addBlxSequence(msp->sname, msp->sname_orig, idTag, strand, dataType, source, seqList, columnList, sequence, msp, lookupTable, &tmpError);
    }

  if (tmpError)
    g_propagate_error(error, tmpError);
}


/* Add or create a BlxSequence struct, creating the BlxSequence if one does not
 * already exist for the MSP's sequence name. Seperate BlxSequence structs are created
 * for the forward and reverse strands of the same sequence. The passed-in sequence 
 * should always be forwards, and we reverse complement it here if we need the 
 * reverse strand. Returns the new BlxSequence */
BlxSequence* addBlxSequence(const char *name, 
                            const char *name_orig, 
			    const char *idTag, 
			    BlxStrand strand,
			    BlxDataType *dataType,
                            const char *source,
			    GList **seqList, 
                            GList *columnList,
			    char *sequence, 
			    MSP *msp, 
                            GHashTable *lookupTable,
			    GError **error)
{
  BlxSequence *blxSeq = NULL;

  if (name || idTag)
    {
      /* If this is an exon or intron the match strand is not applicable. The exon should 
       * be in the same direction as the ref seq, so use the ref seq strand. */
      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
        {
          strand = msp->qStrand;
        }
    
      /* See if this sequence already exists. This matches on name (if linkFeaturesByName is
       * true) or on tag, and strand. */
      gboolean linkFeaturesByName = dataTypeGetFlag(dataType, MSPFLAG_LINK_FEATURES_BY_NAME);
      blxSeq = findBlxSequence(lookupTable, name, idTag, strand, linkFeaturesByName);
      
      if (!blxSeq)
        {
          /* Create a new BlxSequence, and take ownership of the passed in sequence (if any) */
          blxSeq = createBlxSequence(name, idTag, strand, dataType, source, columnList);
          
          /* Add it to the return sequence list (must append it because this function can be
           * called from within a loop which relies on new sequences being appended) */
          *seqList = g_list_append(*seqList, blxSeq);

          /* Add an entry to the lookup table (add an entry for both id and name, if given,
           * because the next feature may have only one or the other set.) */
          if (idTag)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(idTag, strand)), blxSeq);

          if (name && linkFeaturesByName)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(name, strand)), blxSeq);
        }
      else
        {
          if (dataType && !blxSeq->dataType)
            blxSeq->dataType = dataType;
          else if (dataType && blxSeq->dataType != dataType)
            g_warning("Duplicate sequences have different data types [name=%s, ID=%s, strand=%d, orig type=%s, new type=%s].\n", name, idTag, strand, g_quark_to_string(blxSeq->dataType->name), g_quark_to_string(dataType->name));
          
          const char *oldSource = blxSequenceGetSource(blxSeq);

          if (source && !oldSource)
            blxSequenceSetValueFromString(blxSeq, BLXCOL_SOURCE, source);
          else if (source && !stringsEqual(oldSource, source, FALSE))
            g_warning("Duplicate sequences have different sources [name=%s, ID=%s, strand=%d, orig source=%s, new source=%s].\n", name, idTag, strand, oldSource, source);
        }
      
      if (name && !blxSequenceGetName(blxSeq))
	{
	  /* It's possible that the BlxSequence was created without a name if we found an
	   * unnamed child exon before we found the parent transcript, so set the name if we have it. */
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQNAME, name);

          /* Add an entry to the lookup table keyed on name, now that we know it */
          if (linkFeaturesByName)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(name, strand)), blxSeq);
	}
      
      if (msp)
        {
          /* Add the MSP to the BlxSequence's list. Keep it sorted by position. */
          blxSeq->mspList = g_list_insert_sorted(blxSeq->mspList, msp, compareFuncMspPos);
          msp->sSequence = blxSeq;
          
          if (blxSeq->type == BLXSEQUENCE_UNSET)
            {
              blxSeq->type = getBlxSequenceTypeForMsp(msp->type);
            }
          else if (blxSeq->type != getBlxSequenceTypeForMsp(msp->type))
            {
              g_warning("Adding MSP of type %d to parent of type %d (expected parent type to be %d)\n", 
			msp->type, blxSeq->type, getBlxSequenceTypeForMsp(msp->type));
            }
        }
      
      /* Add the sequence data */
      addBlxSequenceData(blxSeq, sequence, error);
    }
  else
    {
      g_set_error(error, BLX_ERROR, 1, "Sequence name or parent ID must be set.\n");
    }
  
  return blxSeq;
}


/* Add the given sequence data to a BlxSequence. Validates that the existing sequence data is
 * either null or is the same as the new data; sets the given error if not. We claim ownership
 * of the given sequence data (either the BlxSequence owns it, or we delete it if it is not 
 * required). The given sequence should always be the forward strand; we complement it ourselves
 * here if this BlxSequence requires the reverse strand. */
void addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error)
{
  if (!sequence)
    {
      return;
    }
  
  gboolean sequenceUsed = FALSE;

  const char *oldSequence = blxSequenceGetSequence(blxSeq);
  
  if (blxSeq && blxSequenceRequiresSeqData(blxSeq))
    {
      if (!oldSequence)
        {
          /* Sequence does not yet exist, so add it */
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQUENCE, sequence);
          sequenceUsed = TRUE;
        }
      else if (error && *error)
        {
          /* Sequence already exists. Validate that it's the same as the existing one. */
          if (!stringsEqual(sequence, oldSequence, FALSE))
            {
              g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_DATA_MISMATCH, "Sequence data for '%s' does not match previously-found data.\n", blxSequenceGetName(blxSeq));
            }
        }
    }      
  
  if (!sequenceUsed)
    {
      g_free(sequence);
    }
}


/***********************************************************
 *                      Columns
 ***********************************************************/

/* Creates a data "column" from the given info and adds it to the columnList. */
void blxColumnCreate(BlxColumnId columnId, 
                     const gboolean createHeader,
                     const char *title,
                     GType type,
                     const char *propertyName,
                     const int defaultWidth,
                     const gboolean dataLoaded,
                     const gboolean showColumn,
                     const gboolean showSummary,
                     const gboolean canShowSummary,
                     const gboolean searchable,
                     const char *sortName,
                     const char *emblId,
                     const char *emblTag,
                     GList **columnList)
{
  /* Create a simple label for the header (unless told not to) */
  GtkWidget *headerWidget = NULL;
  
  if (createHeader)
    {
      headerWidget = createLabel(title, 0.0, 1.0, TRUE, TRUE, TRUE);
      gtk_widget_set_size_request(headerWidget, defaultWidth, -1);
    }
  
  /* Create the column info */
  BlxColumnInfo *columnInfo = (BlxColumnInfo*)g_malloc(sizeof(BlxColumnInfo));

  static int columnIdx = 0;  
  columnInfo->columnIdx = columnIdx;
  ++columnIdx;

  columnInfo->columnId = columnId;
  columnInfo->headerWidget = headerWidget;
  columnInfo->refreshFunc = NULL;
  columnInfo->title = title;
  columnInfo->propertyName = propertyName;
  columnInfo->width = defaultWidth;
  columnInfo->sortName = sortName;
  columnInfo->emblId = g_quark_from_string(emblId);
  columnInfo->emblTag = g_quark_from_string(emblTag);
  columnInfo->dataLoaded = TRUE;
  columnInfo->showColumn = showColumn;
  columnInfo->showSummary = showSummary;
  columnInfo->canShowSummary = canShowSummary;
  columnInfo->searchable = searchable;
  columnInfo->type = type;
  
  /* Place it in the list. List must be sorted in the same order
   * as the GtkListStore or gtk_list_store_set fails */
  *columnList = g_list_insert_sorted(*columnList, columnInfo, columnIdxCompareFunc);
}


/* Create a list of CDS names in a transcript. If there are multiple names it means we have
 * e.g. different start codons and we need to create multiple versions of the transcript */
GList* blxSequenceConstructCdsList(BlxSequence *seq)
{
  GList *result = NULL;

  /* If none of the CDSs have names, that's fine - we just assume they're all the same. It's
   * ambiguous what we should do if some have names and others don't, though. I guess we just have
   * to lump all the nameless ones together. */
  if (seq && seq->type == BLXSEQUENCE_TRANSCRIPT && seq->mspList && g_list_length(seq->mspList) > 0)
    {
      GList *mspItem = seq->mspList;

      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP *msp = (const MSP*)(mspItem->data);
          
          if (msp->type == BLXMSP_CDS)
            {
              GQuark name = 0;

              if (msp->sname)
                name = g_quark_from_string(msp->sname);

              if (!g_list_find(result, GINT_TO_POINTER(name)))
                result = g_list_append(result, GINT_TO_POINTER(name));
            }
        }
    }

  return result;
}
