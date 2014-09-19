/*  File: version.h
 *  Author: Ed Griffiths
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
 * Description: Macros to support version numbering of libraries and
 *              applications in SeqTools.
 *----------------------------------------------------------------------------
 */

#ifndef UT_VERSION_H
#define UT_VERSION_H


/* Tools for creating version strings in an application or library.          */
/*                                                                           */


/* This macro creates a routine that must be provided by all applications    */
/* that use the ACECB kernel code or libace. libace routines expect to be    */
/* able to query the date on which the applications main routine was         */
/* compiled so that this information can be displayed to the user.           */
/* The function must have this prototype and must return a string that gives */
/* the build date:                                                           */
/*                                                                           */
char *utAppGetCompileDate(void) ;
/* The acedb                                                                 */
/* makefile is arranged so that the main routine is recompiled every time    */
/* the application is relinked. This means that the date represents the      */
/* 'build' date of the application.                                          */
/*                                                                           */
/* Code the macro by simply putting it in the .c file that contains the      */
/* main function of the application, it's probably best to put it just       */
/* before or after the main function. Do not put a terminating ';' after     */
/* the macro, this will cause a compile error.                               */
/*                                                                           */
#define UT_COMPILE_PHRASE "compiled on:"

#define UT_MAKE_GETCOMPILEDATEROUTINE()                                      \
char *utAppGetCompileDate(void) { return UT_COMPILE_PHRASE " " __DATE__ " " __TIME__ ; }



/* These tools assume that various numbers/strings are defined, e.g.         */
/*                                                                           */
/* #define SOME_TITLE   "UT library"         (definitive name for library)   */
/* #define SOME_DESC    "brief description"  (purpose of library - one liner)*/
/* #define SOME_VERSION 1                    (major version)                 */
/* #define SOME_RELEASE 0                    (minor version)                 */
/* #define SOME_UPDATE  1                    (latest fix number)             */
/*                                                                           */


/*  Use UT_MAKESTRING to make strings out of #define'd numbers.          
 *  (required because of the way ANSI preprocessor handles strings)      
 *  e.g. UT_MAKESTRING(6)  produces "6"                                 */   
#define UT_PUTSTRING(x) #x
#define UT_MAKESTRING(x) UT_PUTSTRING(x)


/* Make a single version number out of a version, release and update number.
 * NOTE that there should be no more than 100 (i.e. 0 - 99) revisions per    
 * version, or updates per revision, otherwise version will be wrong. */
#define UT_MAKE_VERSION_NUMBER(VERSION, RELEASE, UPDATE) \
((VERSION * 10000) + (RELEASE * 100) + UPDATE)


/* Make a single version string out of the version, release and update numbers */
#define UT_MAKE_VERSION_STRING(VERSION, RELEASE, UPDATE) \
UT_MAKESTRING(VERSION) "." UT_MAKESTRING(RELEASE) "." UT_MAKESTRING(UPDATE)


/* Make a title string containing the title of the application/library and the version. */
#define UT_MAKE_TITLE_STRING(TITLE, VERSION) \
TITLE " - " VERSION


/* Make a string containing the compile time and date */
#define UT_MAKE_COMPILE_DATE() \
__TIME__ " " __DATE__


/* Make a copyright string, where the copyright for the given year(s) (passed as a string, e.g. "2009 - 2010") */
#define UT_MAKE_COPYRIGHT_STRING(YEARS_STRING) \
"Copyright (c) " YEARS_STRING ": Genome Research Ltd."


/* Make a version-info string (has the package name and version string) */
#define UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, VERSION_STRING) \
PACKAGE_NAME" - "VERSION_STRING


/* Make a licence string */
#define UT_MAKE_LICENCE_STRING(TITLE) \
TITLE" is distributed under the GNU Public License; see http://www.gnu.org/copyleft/gpl.txt"


/* Define the authors of the SeqTools package. AUTHOR_LIST is a comma-separated list of all authors
 * that should be credited. AUTHOR_TEXT is a string containing the main authors. */
#define AUTHOR_LIST	   "Gemma Barson (Sanger Institute, UK) <gb10@sanger.ac.uk>",\
                           "Erik Sonnhammer (SBC, Sweden) <Erik.Sonnhammer@sbc.su.se>",\
                           "Jean Thierry-Mieg (CRBM du CNRS, France) <mieg@kaa.crbm.cnrs-mop.fr>",\
                           "Richard Durbin (Sanger Institute, UK) <rd@sanger.ac.uk>",\
                           "Ed Griffiths (Sanger Institute, UK) <edgrif@sanger.ac.uk>",\
                           "Roy Storey (Sanger Institute, UK) <rds@sanger.ac.uk>",\
                           "Malcolm Hinsley (Sanger Institute, UK) <mh17@sanger.ac.uk>"

#define AUTHOR_TEXT        "Gemma Barson <gb10@sanger.ac.uk>\n"\
                           "Erik Sonnhammer <Erik.Sonnhammer@sbc.su.se>"
                           
#define AUTHOR_TEXT_FULL   " Written by Gemma Barson <gb10@sanger.ac.uk>\n"\
                           " Based on original code by Erik Sonnhammer <Erik.Sonnhammer@sbc.su.se>"
                           


/*    Macro for creating a standard copyright string to be inserted into     
 *    compiled applications and libraries. The macro ensures a common        
 *    format for version numbers etc.                                        
 *                                                                           
 * The macro is a statement, NOT an expression, but does NOT require a       
 * terminating semi-colon. The macro should be coded like this:              
 *                                                                           
 *    UT_COPYRIGHT_STRING(prefix, title, description, copyright_years)       
 *                                                                           
 *    where  prefix is some a string locally used to prefix variables        
 *    where  title is a string of the form   "Appname  1.0.1"                
 *    where  description is of the form  "Application to blah, blah."        
 *    and    copyright_years is of the form  "2010", "2008-2010", "2008, 2010" etc.
 */
#define UT_COPYRIGHT(YEARS_STRING)                                                   \
"@(#) "UT_MAKE_COPYRIGHT_STRING(YEARS_STRING)"\n"                                    \
"@(#) \n"                                                                            \
"@(#) This file contains the above Sanger Informatics Group library, \n"             \
"@(#) "AUTHOR_TEXT"\n\n"                                                             \
"@(#) You may redistribute this software subject to the conditions in the \n"        \
"@(#) accompanying copyright file. Anyone interested in obtaining an up to date \n"  \
"@(#) version should contact one of the authors at the above email addresses. \n"


#define UT_COPYRIGHT_STRING(TITLE, VERSION, DESCRIPTION_STRING, COPYRIGHT_YEARS)     \
"@(#) \n"                                                                            \
"@(#) --------------------------------------------------------------------------\n"  \
"@(#) Title/Version:  "UT_MAKE_TITLE_STRING(TITLE, VERSION)"\n"                      \
"@(#)      Compiled:  "UT_MAKE_COMPILE_DATE"\n"                                      \
"@(#)   Description:  "DESCRIPTION_STRING"\n\n"                                        \
UT_COPYRIGHT()                                                                       \
"@(#) --------------------------------------------------------------------------\n"  \
"@(#) \n" ;


#endif	/* UT_VERSION_H */
