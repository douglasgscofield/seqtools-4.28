/*  File: blxFetchDb.c
 *  Author: Gemma Barson, 2012-08-06
 *  Copyright (c) 2012 Genome Research Ltd
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
 * Description: Blixem functions for control of sequence fetching and
 *              display.
 *              
 *              Compiling with -DSQLITE3 includes code to issue
 *              fetch requests to an sqlite3 database.
 *              This requires sqlite3 libs to be installed.
 *----------------------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>

//#ifdef SQLITE3
#include <sqlite3.h>
//#endif

#include <seqtoolsUtils/utilities.h>
#include <blixemApp/blixem_.h>
#include <blixemApp/blxwindow.h>

/* Error codes and domain */
#define BLX_FETCH_ERROR g_quark_from_string("Blixem config")

/* Function pointer for sqlite callback functions */
typedef int (*SqliteFunc)(void*,int,char**,char**);


typedef struct _SqliteFetchData
{
  GList *seqList;
  GList *columnList;
} SqliteFetchData;



//#ifdef SQLITE3

/********************/
/* Local functions */
/********************/

/* Utility to find the given sequence in the given list */
static BlxSequence* findBlxSequence(const char *seqName, GList *seqList)
{
  BlxSequence* result = NULL;
  GList *item = seqList;

  for (; item && !result; item = item->next)
    {
      BlxSequence* blxSeq = (BlxSequence*)(item->data);
      const char *curName = blxSequenceGetName(blxSeq);
      
      if (stringsEqual(curName, seqName, FALSE))
        result = blxSeq;
    }

  return result;
}


/* Callback to print the results of an sql query */
//static int printResultsCB(void *NotUsed, int argc, char **argv, char **azColName)
//{
//  int i;
//
//  for (i = 0; i < argc; ++i)
//    {
//      g_message("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
//    }
//
//  g_message("\n");
//  return 0;
//}


/* Callback to display the results of an sql query in a pop-up dialog.
 * The parent window is passed in the user data. */
static int displayResultsCB(void *data, int argc, char **argv, char **azColName)
{
  GtkWidget *parent = NULL;

  if (data)
    parent = GTK_WIDGET(data);

  /* Compile the results into a single string */
  GString *result = g_string_new("");
  int i;

  for (i = 0; i < argc; ++i)
    {
      g_string_append_printf(result, "%s = %s\n\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }

  char *title = g_strdup_printf("%s - %s", g_get_prgname(), "SQL query");
  displayFetchResults(title, result->str, parent, NULL, NULL);
  g_free(title);

  return 0;
}


/* Callback to populate the results of an sql query into a list of sequences */
static int populateResultsListCB(void *data, int argc, char **argv, char **azColName)
{
  DEBUG_ENTER("populateResultsListCB");

  static GQuark nameCol = 0;
  if (!nameCol)
    {
      nameCol = g_quark_from_string("Name");
    }
  
  SqliteFetchData *fetchData = (SqliteFetchData*)data;

  /* Loop through the columns to find the name column,
   * and find the BlxSequence with that name */
  int i;
  BlxSequence *blxSeq = NULL;

  for (i = 0; i < argc && !blxSeq; ++i)
    {
      GQuark column = g_quark_from_string(azColName[i]);
      
      if (column == nameCol)
        {
          blxSeq = findBlxSequence(argv[i], fetchData->seqList);
        }
    }

  /* Loop again to populate the sequence data */
  if (blxSeq)
    {
      for (i = 0; i < argc; ++i)
        {
          DEBUG_OUT("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
          blxSequenceSetColumn(blxSeq, azColName[i], argv[i], fetchData->columnList);
        }
    }

  DEBUG_EXIT("populateResultsListCB");
  return 0;
}


/* Execute the given sql query and call the callback on the results. */
static void sqliteRequest(const char *database, const char *query, SqliteFunc callbackFunc, void *callbackData, GError **error)
{
  DEBUG_ENTER("sqliteRequest");

  //#ifdef SQLITE3
  sqlite3 *db;

  int rc = sqlite3_open(database, &db);
  
  if (rc)
    {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      sqlite3_close(db);
      return;
    }

  DEBUG_OUT("Database: %s\nQuery: %s\n", database, query);

  char *zErrMsg = 0;
  rc = sqlite3_exec(db, query, callbackFunc, callbackData, &zErrMsg);

  if (rc != SQLITE_OK)
    {
      fprintf(stderr, "SQL error: %s\n", zErrMsg);
      sqlite3_free(zErrMsg);
    }

  sqlite3_close(db);
//#else
//  g_set_error(error, BLX_ERROR, 1, "SQLite not available: sqlite3 may not be installed on your system, or has not been compiled into the package correctly.\n");
//#endif

  DEBUG_EXIT("sqliteRequest");
}


/* Check that the fetch method is non-null and has a db and query */
static void validateFetchMethod(const BlxFetchMethod* const fetchMethod, GError **error)
{
  if (!fetchMethod)
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_NULL_FETCH, "Program error: fetch method is null\n");
    }
  else if (!fetchMethod->location)
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_NO_EXE, "No database file specified for fetch method '%s'\n", g_quark_to_string(fetchMethod->name));
    }
  else if (!fetchMethod->args)
    {
      g_set_error(error, BLX_CONFIG_ERROR, BLX_CONFIG_ERROR_NO_ARGS, "No query specified for fetch method '%s'\n", g_quark_to_string(fetchMethod->name));
    }
}


//#endif


/********************/
/* Public functions */
/********************/


/* Fetch a single sequence using sqlite. If displayResults is true, disply the
 * results in a pop-up window. If result_out is non-null, populate it with the result. */
void sqliteFetchSequence(const BlxSequence* const blxSeq, 
                         const BlxFetchMethod* const fetchMethod,
                         const gboolean displayResults,
                         const int attempt,
                         GtkWidget *blxWindow)
{
  DEBUG_ENTER("sqliteFetchSequence");

  GError *tmpError = NULL;
  validateFetchMethod(fetchMethod, &tmpError);
    
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GString *query = NULL;
  
  if (!tmpError)
    {
      query = getFetchArgs(fetchMethod, blxSeq, NULL, 
                           bc->refSeqName, bc->refSeqOffset, &bc->refSeqRange,
                           bc->dataset, &tmpError);
    }
  
  if (query && !tmpError)
    {
      sqliteRequest(fetchMethod->location, 
                    query->str,
                    displayResultsCB,
                    blxWindow,
                    &tmpError);
    }

  g_string_free(query, TRUE);
  

  if (tmpError)
    {
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
      fetchSequence(blxSeq, displayResults, attempt + 1, blxWindow, NULL, NULL);
    }
    
  DEBUG_EXIT("sqliteFetchSequence");
}


/* Fetch multiple sequences using sqlite */
void sqliteFetchSequences(GList *seqsToFetch, 
                          const BlxFetchMethod* const fetchMethod, 
                          GList *columnList,
                          GError **error)
{
  DEBUG_ENTER("sqliteFetchSequences");

  GError *tmpError = NULL;
  validateFetchMethod(fetchMethod, &tmpError);
    
  GString *query = NULL;

  if (!tmpError)
    {
      query = getFetchArgsMultiple(fetchMethod, seqsToFetch, &tmpError);
    }

  SqliteFetchData fetchData = {seqsToFetch, columnList};
  
  if (query && !tmpError)
    {
      sqliteRequest(fetchMethod->location, 
                    query->str,
                    populateResultsListCB,
                    &fetchData,
                    &tmpError);
    }

  g_string_free(query, TRUE);
  
  if (tmpError)
    g_propagate_error(error, tmpError);

  DEBUG_EXIT("sqliteFetchSequences");
}

