/*  File: dotter.h
 *  Author: Erik Sonnhammer, 1999-08-26
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
 * Description: 
 * Exported functions:
 *   Only 3 parameters are mandatory, the rest can be set to NULL.
 *   A minimal call would look like:
 *
 *   dotter(type, 0, 0, qseq, 0, 0, sseq, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
 *
 *   NOTE: qseq and sseq must be g_malloc'ed in the calling routine.  
 *   They are g_free'd by Dotter.
 *----------------------------------------------------------------------------
 */

#ifndef DEF_DOTTER_H
#define DEF_DOTTER_H

#include <seqtoolsUtils/blxmsp.h>


/* Possible graphical formats the plot can be exported as */
typedef enum _DotterExportFormat
  {
    DOTTER_EXPORT_NONE,         /* Don't export the plot to a graphical format */
    DOTTER_EXPORT_PDF,          /* PDF */
    DOTTER_EXPORT_PS,           /* Postscript */
    DOTTER_EXPORT_SVG           /* Scalable vector graphics */
  } DotterExportFormat;


/* Options specifying the initial state for dotter */
typedef struct _DotterOptions
  {
    int qoffset;              /* qoffset + 1 gives the value of the first coord on the ref seq */
    int soffset;              /* soffset + 1 gives the value of the first coord on the match seq */
    gboolean selfcall;        /* whether called internally, i.e. so that features/sequences will be piped into dotter rather than read from files */
    int qlen;                 /* length of the ref seq */
    int slen;                 /* length of the match seq */
    int dotterZoom;           /* initial zoom level */
    int install : 1;          /* whether to add -install to the dotter args (for private colormaps) */
    int pixelFacset;
    int seqInSFS;             /* whether the sequences are in the features file, i.e. there are no separate sequence files */
    
    float memoryLimit;
    
    char *savefile;           /* file to save the dot-plot to (batch mode; saves the dot-matrix so it can be loaded later and interacted with) */
    char *exportfile;         /* file to export the dot-plot to (batch mode; exports to a graphical format, e.g. pdf or ps. Default is pdf unless file extension indicates otherwise) */
    char *loadfile;           /* file to load a dot-plot from */
    char *FSfilename;         /* file containing features i.e. MSPs */
    char *mtxfile;            /* caller-supplied matrix file */
    
    char *winsize;            /* caller-supplied sliding-window size */
    
    char *qname;              /* reference (horizontal) sequence name */
    char *qseq;               /* reference (horizontal) sequence data */
    char *sname;              /* match (vertical) sequence name */
    char *sseq;               /* match (vertical) sequence data */
    
    gboolean mirrorImage;     /* display mirror image in self comparisons (i.e. so we only have to calculate half of the dot-plot) */
    gboolean watsonOnly;      /* only show the watson (forward) strand of the ref seq */
    gboolean crickOnly;       /* only show the crick (reverse) strand of the ref seq */
    gboolean hspsOnly;        /* only draw HSPs (i.e. don't calculate the dot-plot, just draw lines where we know HSPs should be) */
    int greyrampWhite;        /* greyramp white limit */
    int greyrampBlack;        /* greyramp black limit */
    gboolean swapGreyramp;    /* swap the default black/white points on the greyramp tool (inverts the colors) */
    gboolean breaklinesOn;    /* whether to enable breaklines between sequences */
    gboolean hozScaleRev;     /* revese the horizontal scale */
    gboolean vertScaleRev;    /* revese the vertical scale */
    gboolean negateCoords;    /* negate the displayed coords when the scale is reversed, i.e. so they still appear to increase from left-to-right */
    gboolean suppressScale;   /* suppress printing the coordinate scale on x- and y-axis */
    gboolean abbrevTitle;     /* abbrev window title prefix to save space */
    
    BlxMessageData msgData;   /* data to be passed to the message handlers */

    char *breaklineColor;     /* if not null, color for breaklines between multiple sequences */
    char *windowColor;        /* if not null, background color for the window */
  } DotterOptions;


void dotter(
	const BlxBlastMode blastMode, /* Mandatory, one of { BLASTP, BLASTN, BLASTX } 
			      P -> Protein-Protein
			      N -> DNA-DNA
			      X -> DNA-Protein */
	
	DotterOptions *options, /* Optional, may be NULL 
                                   Various options for display features */

	const BlxStrand refSeqStrand,   /* which strand of the reference sequence was passed */

        const BlxStrand matchSeqStrand, /* which strand of the match sequence was passed */

	int   qcenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	int   scenter,	   /* Optional, may be NULL 
			      Coordinate to centre horisontal sequence on */

	MSP *MSPs,	   /* Optional, may be NULL
			      List of MSPs containing genes and blast matches */

	GList *seqList,	   /* Optional, may be NULL
			      List of all match sequences, as BlxSequences */
	    
	int   MSPoff	   /* Optional, may be NULL
			      Coordinate offset of MSPs */
);


#endif /*  !defined DEF_DOTTER_H */
