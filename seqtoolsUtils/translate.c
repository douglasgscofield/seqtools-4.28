/*  File: translate.c
 *  Author: sre, 1993-01-12
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
 * Description: Functions for complementing, reversing and translating nucleic
 *              acid sequences
 *
 * Exported functions: See utilities.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/iupac.h>
#include <glib.h>
#include <ctype.h>
#include <string.h>


/* Translation errors domain */
#define SEQTOOLS_TRANSLATION_ERROR g_quark_from_string("SeqTools")

/* Error codes */
typedef enum
  {
    SEQTOOLS_ERROR_INVALID_NUCLEOTIDE      /* invalid dna/rna nucleotide */
  } SeqToolsTranslationError;




/* THIS FILE NEEDS RENAMING TO SOMETHING LIKE utils.c */

/* Function: Translate(char *seq, char **code)
 * 
 * Given a ptr to the start of a nucleic acid sequence, and a genetic code, translate the sequence
 * into amino acid sequence.
 * 
 * code is an array of 65 strings, representing the translations of the 64 codons, arranged in
 * order AAA, AAC, AAG, AAU, ..., UUA, UUC, UUG, UUU.  '*' or '***' is used to represent
 * termination codons, usually. The final string, code[64], is the code for an ambiguous amino
 * acid.
 *
 * Because of the way space is allocated for the amino acid sequence, the amino acid strings
 * cannot be longer than 3 letters each. (I don't foresee using anything but the single- and
 * triple- letter codes.)
 * 
 * Returns a ptr to the translation string on success, or NULL on failure.
 */
char *blxTranslate(const char *seq, char **code)
{
  char *aaseq = NULL ;					    /* RETURN: the translation */
  int   codon;						    /* index for codon         */
  char *aaptr;						    /* ptr into aaseq */
  int   i;
  
  if (seq && *seq)
    {
      aaseq = (char *)g_malloc(strlen(seq) + 1) ;

      for (aaptr = aaseq ; *seq != '\0' && *(seq+1) != '\0' && *(seq+2) != '\0'; seq += 3)
	{
	  /* calculate the lookup value for this codon */
	  codon = 0;
	  for (i = 0; i < 3; i++)
	    {
	      codon *= 4;

	      switch (*(seq + i))
		{
		case 'A': case 'a':             break;
		case 'C': case 'c': codon += 1; break;
		case 'G': case 'g': codon += 2; break;
		case 'T': case 't': codon += 3; break;
		case 'U': case 'u': codon += 3; break;
		default: codon = 64; break;
		}

	      if (codon == 64)
		break;
	    }

	  strcpy(aaptr, code[codon]);
	  aaptr += strlen(code[codon]);
	}
    }  

  return aaseq ;
}


/* Get the complement of the given nucleotide. Returns the original char and sets
 * the error if no valid complement exists */
char complementChar(const char inputChar, GError **error)
{
  /* Loop through each iupac code looking for this char. iupac chars are all
   * uppercase */
  char result = '\0';
  char c = toupper(inputChar);
  int idx = 0;
  
  for ( ; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);
  
  if (idx >= IUPACSYMNUM)
    {
      /* not found; return original char */
      result = inputChar;
      g_set_error(error, SEQTOOLS_TRANSLATION_ERROR, SEQTOOLS_ERROR_INVALID_NUCLEOTIDE, "Invalid nucleotide '%c'; could not find complement.\n", inputChar);
    }
  else
    {
      result = iupac[idx].symcomp;
      
      if (islower(inputChar))
        result = tolower(result);
    }
  
  return result;
}


/* All these calls need rationalising into a single function with options. */

/* revComplement.c
 * 
 * Reverse complement of a IUPAC character string
 * 
 */
/* Ratinalise this with my func. below..... */
char *revComplement(char *comp, char *seq)
{
  long  bases;
  char *bckp, *fwdp;
  long  pos;

  if (comp == NULL)
    return NULL;
  if (seq == NULL)
    return NULL;

  bases = strlen(seq);

  fwdp = comp;
  bckp = seq + bases -1;
  for (pos = 0; pos < bases; pos++)
    {
      *fwdp = complementChar(*bckp, NULL);
      fwdp++;
      bckp--;
    }
  *fwdp = '\0';

  return comp;
}
  

/* blxComplement.c
 * 
 * Just complement of a IUPAC character string
 * 
 * Note that it overwrites calling string!!!! (revcomp doesn't)
 */
void blxComplement(char *seq)
{
  char *fwdp;
  long  pos;

  if (seq == NULL)
    return ;

  fwdp = seq;
  for (pos = 0; pos < (int)strlen(seq); pos++)
    {
      *fwdp = complementChar(*fwdp, NULL);
      fwdp++;
    }

  *fwdp = '\0';

  return ;
}
  


