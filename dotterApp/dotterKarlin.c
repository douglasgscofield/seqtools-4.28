/*  File: dotterKarlin.c
 *  Author: Erik Sonnhammer, 1995-08-28
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
 * Description: Karlin/Altschul statistics calculations
 *----------------------------------------------------------------------------
 */


#include <dotterApp/dotter_.h>
#include <gtk/gtk.h>
#include <math.h>
#include <stdlib.h>


#define MAXIT 20	/* Maximum number of iterations used in calculating K */
/* For greater accuracy, set SUMLIMIT to 0.00001 */
/* For higher speed, set SUMLIMIT to 0.001 */
#define SUMLIMIT 0.01


/* integer power function
 * Originally submission by John Spouge, 6/25/90
 * From blast/gish/fct/fct_powi.c
 */
double fct_powi(double x, register int n)
/* x  argument */ 
/* n  power */
{
    register int	i;
    double	y;
    
    y = 1.;
    for (i = abs(n); i > 0; i /= 2) {
	if (i&1)
	    y *= x;
	x *= x;
    }
    if (n >= 0)
	return y;
    return 1./y;
}


/* fct_gcd(a, b)
 * Return the greatest common divisor of a and b.
 * From blast/gish/fct/fct_gcd.c
 */
long fct_gcd(register long a, register long b)
{
    register long	c;
    
    b = abs(b);
    if (b > a)
	c=a, a=b, b=c;
    
    while (b != 0) {
	c = a%b;
	a = b;
	b = c;
    }
    return a;
}


/*
 * Return values accurate to approx. 16 digits for the quantity exp(x)-1
 * for all values of x, both large and small.
 * from blast/gish/fct/fct_expm.c
 */

double fct_expm1(double x)
{
    double absx = ((x < 0) ? -x : x) ;
    
    if (absx > .33)
	return exp(x) - 1.;
    
    if (absx < 1.e-16)
	return x;
    
    return x * (1. + x *
		(0.5 + x * (1./6. + x *
			    (1./24. + x * (1./120. + x *
					   (1./720. + x * (1./5040. + x *
							   (1./40320. + x * (1./362880. + x *
									     (1./3628800. + x * (1./39916800. + x *
												 (1./479001600. + x/6227020800.)
												 ))
									     ))
							   ))
					   ))
			    )));
}


/* etop() -- given an Expect value, return the associated probability 
 * From blast/blast/lib/etop.c
 */
double etop(double E)
{
        return -fct_expm1(-E);
}


/* From blast/blast/lib/karlin.c
 *
 *  long	low;	 	 Lowest score (must be negative)    
 *  long	high;		 Highest score (must be positive)   
 *  double	*pr;		 Probabilities for various scores   
 *  double	*lambda;	 Pointer to parameter lambda        
 *  double	*K;		 Pointer to parmeter K              
 *  double	*H;		 Pointer to parmeter H              
 */
double karlin(long low, long high, double *pr, double *lambda, double *K, double *H)
{

/**************** Statistical Significance Parameter Subroutine ****************

	Version 1.0	February 2, 1990

	Program by:	Stephen Altschul

	Address:	National Center for Biotechnology Information
			National Library of Medicine
			National Institutes of Health
			Bethesda, MD  20894

	Internet:	altschul@ncbi.nlm.nih.gov

See:	Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
	Significance of Molecular Sequence Features by Using General Scoring
	Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

	Computes the parameters lambda and K for use in calculating the
	statistical significance of high-scoring segments or subalignments.

	The scoring scheme must be integer valued.  A positive score must be
	possible, but the expected (mean) score must be negative.

	A program that calls this routine must provide the value of the lowest
	possible score, the value of the greatest possible score, and a pointer
	to an array of probabilities for the occurence of all scores between
	these two extreme scores.  For example, if score -2 occurs with
	probability 0.7, score 0 occurs with probability 0.1, and score 3
	occurs with probability 0.2, then the subroutine must be called with
	low = -2, high = 3, and pr pointing to the array of values
	{ 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
	pointers to lambda and K; the subroutine will then calculate the values
	of these two parameters.  In this example, lambda=0.330 and K=0.154.

	The parameters lambda and K can be used as follows.  Suppose we are
	given a length N random sequence of independent letters.  Associated
	with each letter is a score, and the probabilities of the letters
	determine the probability for each score.  Let S be the aggregate score
	of the highest scoring contiguous segment of this sequence.  Then if N
	is sufficiently large (greater than 100), the following bound on the
	probability that S is greater than or equal to x applies:
	
		P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].
	
	In other words, the p-value for this segment can be written as
	1-exp[-KN*exp(-lambda*S)].

	This formula can be applied to pairwise sequence comparison by assigning
	scores to pairs of letters (e.g. amino acids), and by replacing N in the
	formula with N*M, where N and M are the lengths of the two sequences
	being compared.

	In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
	distinct segments all with score >= S is given by:

                               2             m-1           -y
                1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

	Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
	as the previous formula.

*******************************************************************************/

    int		i, j;
    long	range, lo, hi, first, last;
    double	up, new_val, sum, Sum, av, beta, oldsum, oldsum2;
    double	*p = NULL, *P = NULL, *ptrP, *ptr1, *ptr2;
    double	ratio;

    /* Check that scores and their associated probabilities are valid     */

    if (low >= 0.) {
	g_critical("Karlin-Altschul statistics error: There must be at least one negative score in the substitution matrix.");
	return -1.0;
    }

    for (i=range=high-low; i > -low && pr[i] == 0.0; --i);
    if (i <= -low) {
	g_critical("Karlin-Altschul statistics error: A positive score is impossible in the context of the scoring scheme, the residue composition of the query sequence, and the residue composition assumed for the database.");
	return -1.0;
    }

    for (sum=i=0; i<=range; sum += pr[i++])
	if (pr[i] < 0.) {
	    g_critical("Karlin-Altschul statistics error: Negative probabilities for scores are disallowed.");
	    return -1.0;
	}

    if (sum<0.99995 || sum>1.00005)
      g_message("Score probabilities sum to %.5lf and will be normalized to 1.\n", sum);

    p = (double *)g_malloc(sizeof(*p) * (range+1));
    for (Sum=low,i=0; i<=range; ++i)
	Sum += i*(p[i]=pr[i]/sum);

    if (Sum >= 0.) {
	g_critical("Karlin/Altschul statistics failed due to non-negative expected score: %#0.3lg", Sum);
	return Sum;
    }

    /* Calculate the parameter lambda */

    up = 0.5;
    do {
	up *= 2;
	ptr1 = p;
	for (sum=0,i=low; i<=high; ++i)
	    sum += *ptr1++ * exp(up*i);
    } while (sum<1.0);
    
    /* Root solving by the bisection method */
    for (*lambda=0., j=0; j<25; ++j) {
	new_val = (*lambda+up)/2.0;
	ptr1 = p;
	for (sum=0., i=low; i <= high; ++i)
	    sum += *ptr1++ * exp(new_val*i);
	if (sum > 1.0)
	    up = new_val;
	else
	    *lambda = new_val;
    }
    beta = exp(*lambda);


    /* Calculate the relative entropy of the p's and q's, parameter H */
    ptr1 = p;
    for (av=0, i=low; i<=high; ++i)
	av += *ptr1++ *i*exp(*lambda*i);
    *H = *lambda*av;
    
    /* Calculate the parameter K */
    if (low == -1 || high == 1) {
	*K = (high == 1 ? av : Sum*Sum/av);
	*K *= 1.0 - 1./beta;
	goto OKExit;
    }
    Sum = 0.;
    lo = hi = 0;
    P = (double *)g_malloc(MAXIT* (range+1) * sizeof(*P));
    *P = sum = oldsum = oldsum2 = 1.;
    for (j=0;  j<MAXIT && sum > SUMLIMIT; oldsum = sum, Sum += sum /= ++j) {
	first = last = range;
	for (ptrP = P + (hi += high) - (lo += low); ptrP >= P; *ptrP-- =sum) {
	    ptr1 = ptrP - first;
	    ptr2 = p + first;
	    for (sum=0., i=first; i <= last; ++i)
		sum += *ptr1-- * *ptr2++;
	    if (first != 0)
		--first;
	    if (ptrP-P <= range)
		--last;
	}
	new_val = fct_powi(beta, lo-1);
	for (sum=0, i=lo; i != 0; ++i)
	    sum += *++ptrP * (new_val *= beta);
	for (; i <= hi; ++i)
	    sum += *++ptrP;
	oldsum2 = oldsum;
    }


    /* Geometric progression correction terms to accommodate fewer iterations */
    ratio = oldsum / oldsum2;
    if (ratio >= (1.0 - SUMLIMIT*0.001))
      {
	g_critical ("Value calculated for K was too high due to insufficient iterations.  "
		   "Fudging it.") ;
	*K = 0.1 ;
	goto OKExit ;
/* was:
        g_error("Value calculated for K was too high due to insufficient iterations.  "
       	      "Alternatively, the expected average score is insufficiently negative.") ;
*/
      }
    while (sum > SUMLIMIT*0.01) {
	oldsum *= ratio;
	Sum += sum = oldsum / ++j;
    }
    
    for (i=low; p[i-low] == 0.; ++i)
	;
    for (j= -i;i<high && j>1;)
	if (p[++i-low])
	    j = fct_gcd(j,i);
    
    *K = (j*exp(-2.*Sum))/(av*etop(*lambda * j));

OKExit:
	   
    g_free(p);
    g_free(P);
    return 0;		/* Parameters calculated successfully */
}


/* Adapted from blastp.c */
int winsizeFromlambdak(gint32 mtx[24][24], int *tob, int abetsize, const char *qseq, const char *sseq, 
		       double *exp_res_score, double *Lambda)
{
    gint32 
        lows=0, highs=0,
        range;
  
    int    
	i, j,
	*n1, *n2,
	qlen=0, slen=0,
	retval,
	n = 100;		/* Nominal size of dot-matrix */
    double  
	*fq1, *fq2, *prob, K, H,
	qij, exp_MSP_score, sum;
    
    
    n1 = (int *)g_malloc((abetsize+4)*sizeof(int));
    n2 = (int *)g_malloc((abetsize+4)*sizeof(int));
    fq1 = (double *)g_malloc((abetsize+4)*sizeof(double));
    fq2 = (double *)g_malloc((abetsize+4)*sizeof(double));
    

    /* Find high and lows score in score matrix */
    for (i = 0; i < abetsize; ++i)
	for (j = 0; j <  abetsize; ++j) {
	    if (mtx[i][j] < lows ) lows  = mtx[i][j];
	    if (mtx[i][j] > highs) highs = mtx[i][j];
	}


    /* Sum counts of residues */
    for (i = 0; i < abetsize; ++i)
      {
	n1[i] = 0;
      }
    for (i = 0; qseq[i]; ++i)
      {
	/* only count unambiguous letters */
	if (tob[(int)qseq[i]] < abetsize )
	  {
	    n1[tob[(int)qseq[i]]]++;
	    qlen++;
	  }
      }
    for (i = 0; i < abetsize; ++i) n2[i] = 0;
    for (i = 0; sseq[i]; ++i) {
	/* only count unambiguous letters */
	if (tob[(int)sseq[i]] != NA ) {
	    n2[tob[(int)sseq[i]]]++;
	    slen++;
	}
    }
		    

    /* Convert counts to frequencies */
    for (i = 0; i < abetsize; ++i) {
	fq1[i] = (double)n1[i] / qlen;
	fq2[i] = (double)n2[i] / slen;
    }	
    

    /* Calculate probability of each score */
    range = highs - lows;
    prob = (double *)g_malloc(sizeof(double)*(range+1));
    for (i = 0; i <= range; ++i) prob[i] = 0.0;
    
    for (i = 0; i < abetsize; ++i)
      {
	for (j = 0; j < abetsize ; ++j)
	  {
	    prob[mtx[i][j]-lows] += fq1[i] * fq2[j];
	  }
      }

    if ((*exp_res_score = karlin(lows, highs, prob, Lambda, &K, &H)))
      {
	g_critical("Setting ad hoc values to winsize=%d and expected score=%.3f", 25, *exp_res_score);
	return 25;
      }
    

    /* Calculate expected score per residue in MSP */
    *exp_res_score = sum = 0;
    for (i = 0; i < abetsize; ++i)
	for (j = 0; j < abetsize ; ++j) {
	    qij = fq1[i]*fq2[j]*exp(*Lambda*mtx[i][j]); /* Is this correct? */
	    sum += qij;
	    *exp_res_score += qij*mtx[i][j];
	}
    if (sum -1.0 > 0.0001)
	g_warning("Warning: SUM(PiPj*exp(Lambda*Sij)) = %f (Should be 1.0)\n", sum);

    exp_MSP_score = (log(n*n) + log(K)) / *Lambda;

    retval = (int) (exp_MSP_score / *exp_res_score + 0.5);

    g_message("Karlin/Altschul statistics for these sequences and score matrix:\n");
    g_message("   K      = %.3f\n", K);
    g_message("   Lambda = %.3f\n", *Lambda);
    g_message("   => Expected MSP score in a %dx%d matrix = %.3f\n", n, n, exp_MSP_score);

    g_message("   Expected residue score in MSP = %.3f\n", *exp_res_score);
    g_message("   => Expected MSP length = %d\n", retval);


    g_free(prob);
    g_free(n1);
    g_free(n2);
    g_free(fq1);
    g_free(fq2);

    return retval;
}
 
