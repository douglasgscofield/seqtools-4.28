/*  File: blxdotter.c
 *  Author: Gemma Barson, 2010-02-03
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
 * Description: See blxdotter.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/blxdotter.h>
#include <blixemApp/blxwindow.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/detailview.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define DEFAULT_DOTTER_RANGE_SELF	2000

typedef struct _DotterDialogData
  {
    GtkWidget *blxWindow;           /* pointer to the main blixem window */
    
    GtkWidget *autoButton;          /* the radio button on the dialog for automatic dotter parameters */
    GtkWidget *manualButton;        /* the radio button on the dialog for manual dotter parameters */

    GtkWidget *startEntry;          /* the text entry box on the dialog for the start coord */
    GtkWidget *endEntry;            /* the text entry box on the dialog for the end coord */
    GtkWidget *zoomEntry;           /* the text entry box on the dialog for the zoom value */
    
    gboolean callOnSelf;            /* whether to call dotter on the query seq versus itself */
    gboolean hspsOnly;              /* whether to call dotter on HSPs only */
    
    char *dotterSSeq;               /* the match sequence to call dotter on */
  } DotterDialogData;


/* Error codes and domain */
#define BLX_DOTTER_ERROR g_quark_from_string("Dotter")

typedef enum {
  BLX_DOTTER_ERROR_NO_SEQS,	      /* no sequences are selected */
  BLX_DOTTER_ERROR_TOO_MANY_SEQS,     /* too many sequences are selected */
  BLX_DOTTER_ERROR_INVALID_SEQ,       /* selected sequence(s) are invalid */
  BLX_DOTTER_ERROR_NOT_FOUND,	      /* failed to find the sequence */
  BLX_DOTTER_ERROR_NO_REF_SEQ,	      /* failed to find the query sequence segment */
  BLX_DOTTER_ERROR_INTERNAL_SEQ,      /* using internally-stored sequence (because fetch failed) */
  BLX_DOTTER_ERROR_NO_MATCHES,        /* there are no matches on the requested sequence */
  BLX_DOTTER_ERROR_NO_SEQ_DATA,       /* the match sequence has no sequence data (e.g. if could not pfetch it) */
  BLX_DOTTER_ERROR_INVALID_STRAND,    /* there are no match sequences on the correct strand */
  BLX_DOTTER_ERROR_NO_EXE,            /* no dotter executable found */
  BLX_DOTTER_ERROR_PIPE,	      /* error opening pipe */
  BLX_DOTTER_ERROR_FORK,	      /* error forking process */
  BLX_DOTTER_ERROR_NEGATED_COORD,     /* warning that an original coord was not in range so we negated it */
  BLX_DOTTER_ERROR_OUT_OF_RANGE       /* dotter coord was out of range */
} BlxDotterError;


/* Local function declarations */
static gboolean	      getDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, const gboolean callOnSelf, const gboolean autoDotter, int *dotterStart, int *dotterEnd, int *dotterZoom, GError **error);
static gboolean	      smartDotterRange(GtkWidget *blxWindow, const char *dotterSSeq, int *dotter_start_out, int *dotter_end_out, GError **error);
static gboolean	      smartDotterRangeSelf(GtkWidget *blxWindow, int *dotter_start_out, int *dotter_end_out, GError **error);
static gboolean	      callDotterSelf(GtkWidget *blxWindow, GError **error);


/*******************************************************************
 *                      Dotter settings dialog                     *
 *******************************************************************/

/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It sets the dotter mode according to the toggle state of the 
 * "auto" button. */
static gboolean onSaveDotterMode(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->autoDotter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the state of the 'call on self' button. */
static gboolean onSaveDotterSelf(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->dotterSelf = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the state of the 'HSPs only' button. */
static gboolean onSaveDotterHsps(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->dotterHsps = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  return TRUE;
}


/* Check whether the given dotter coord is within the blixem ref seq range.
 * The coord should be as the user sees it, i.e. negated if the 'negate coords'
 * option is enabled and the display is reversed - in this case in input arg
 * will be updated to reflect the real coord. 
 * Returns false and sets the error if not in range. Also sets the error but 
 * returns true if we succeeded but with a warning. */
static gboolean boundsCheckDotterCoord(int *coordIn, BlxViewContext *bc, GError **error)
{
  gboolean ok = TRUE;
  
  const gboolean negate = bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS];
  
  int coord = *coordIn;
  
  if (negate)
    {
      /* Display coords are negated - un-negate to get the real coord */
      coord *= -1;
    }
  
  if (!valueWithinRange(coord, &bc->refSeqRange))
    {
      /* Try negating it in case the user missed the minus sign off. */
      coord *= -1;
      
      if (valueWithinRange(coord, &bc->refSeqRange))
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NEGATED_COORD,
                      "Negated coord '%d' because original coord was not in range.\n", *coordIn);
        }
      else
        {
          ok = FALSE;
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_OUT_OF_RANGE,
                      "Coord '%d' is outside reference sequence range [%d -> %d].\n", 
                      *coordIn, 
                      (negate ? bc->refSeqRange.max * -1 : bc->refSeqRange.min),
                      (negate ? bc->refSeqRange.min * -1 : bc->refSeqRange.max));
        }
    }
  
  *coordIn = coord;
  
  return ok;
}


/* Callback to be called when the user clicks OK or Apply on the dotter
 * dialog. It saves the start parameter that was entered (if manual dotter
 * params are being used). */
static gboolean onSaveDotterStart(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Only save the parameter if we are using manual parameters */
  if (!bc->autoDotter)
    {
      int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      GError *error = NULL;
      
      result = boundsCheckDotterCoord(&newVal, bc, &error);
      
      if (result)
        {
          bc->dotterStart = newVal;
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }
      else
        {
          postfixError(error, "Value not saved.");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterEnd(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!bc->autoDotter)
    {
      int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));

      GError *error = NULL;
      
      result = boundsCheckDotterCoord(&newVal, bc, &error);
      
      if (result)
        {
          bc->dotterEnd = newVal;
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
        }
      else
        {
          postfixError(error, "Value not saved.");
          reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
        }
    }  
  
  return result;
}

static gboolean onSaveDotterZoom(GtkWidget *entry, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only save the parameter if we are using manual parameters */
  if (!blxContext->autoDotter)
    {
      const int newVal = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
      
      if (newVal < 0)
        {
          result = FALSE;
          g_critical("Zoom must be greater than zero (or equal to zero to use auto zoom)");
        }
      else
        {
          blxContext->dotterZoom = newVal;
        }
    }  
  
  return result;
}


/* Get the display coord version of the given coord (i.e. negated if the
 * 'negate coords' option is enabled and the display is reversed). Only 
 * applicable to reference sequence coords. */
static int getDisplayCoord(const int coordIn, BlxViewContext *bc)
{
  int result = coordIn;
  
  if (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS])
    {
      result *= -1;
    }
  
  return result;
}


/* Called when the 'last saved' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end coords that were last saved */
static void onLastSavedButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  char *startString = convertIntToString(getDisplayCoord(bc->dotterStart, bc));
  char *endString = convertIntToString(getDisplayCoord(bc->dotterEnd, bc));
  char *zoomString = convertIntToString(bc->dotterZoom);
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the 'full range' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end of the full ref seq range */
static void onFullRangeButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  const int startCoord = (bc->displayRev ? bc->refSeqRange.max : bc->refSeqRange.min);
  const int endCoord = (bc->displayRev ? bc->refSeqRange.min : bc->refSeqRange.max);
  
  char *startString = convertIntToString(getDisplayCoord(startCoord, bc));
  char *endString = convertIntToString(getDisplayCoord(endCoord, bc));
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the 'big picture range' button in the dotter dialog is clicked. Populates
 * the coord boxes with the start/end of the big picture's current display range */
static void onBpRangeButtonClicked(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(dialogData->blxWindow);
  const IntRange* const displayRange = bigPictureGetDisplayRange(bigPicture);

  int qStart = convertDisplayIdxToDnaIdx(displayRange->min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  int qEnd = convertDisplayIdxToDnaIdx(displayRange->max, bc->seqType, bc->numFrames, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange);

  boundsLimitValue(&qStart, &bc->refSeqRange);
  boundsLimitValue(&qEnd, &bc->refSeqRange);
  
  char *startString = convertIntToString(getDisplayCoord(qStart, bc));
  char *endString = convertIntToString(getDisplayCoord(qEnd, bc));
  
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  
  g_free(startString);
  g_free(endString);
  
  /* Change the mode to manual */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(dialogData->manualButton), TRUE);
}


/* Called when the "call on self" button is toggled on the dotter dialog */
static void onSelfButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  dialogData->callOnSelf = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  /* Recalculate auto start/end */
  int autoStart = UNSET_INT, autoEnd = UNSET_INT;
  const gboolean autoDotter = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton));

  getDotterRange(dialogData->blxWindow, dialogData->dotterSSeq, dialogData->callOnSelf, autoDotter, &autoStart, &autoEnd, NULL, NULL);

  if (autoStart == UNSET_INT)
    autoStart = bc->displayRev ? bc->refSeqRange.max : bc->refSeqRange.min;
  
  if (autoEnd == UNSET_INT)
    autoEnd = bc->displayRev ? bc->refSeqRange.min : bc->refSeqRange.max;
  
  char *startString = convertIntToString(getDisplayCoord(autoStart, bc));
  char *endString = convertIntToString(getDisplayCoord(autoEnd, bc));
  char *zoomString = convertIntToString(bc->dotterZoom);

  /* Display the new values */
  gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
  gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
  
  g_free(startString);
  g_free(endString);
  g_free(zoomString);
}


/* Called when the "call on self" button is toggled on the dotter dialog */
static void onHspsButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  dialogData->hspsOnly = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
}


/* Called when the user has hit a response button on the dotter settings dialog */
static void onResponseDotterDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  DotterDialogData *dialogData = (DotterDialogData*)(data);
  
  GError *error = NULL;
  
  switch (responseId)
    {
      case GTK_RESPONSE_ACCEPT:
        /* Only continue to call dotter if saving the values is successful */
	if (widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId)))
          {
            if (dialogData->callOnSelf)
              {
                destroy = callDotterSelf(dialogData->blxWindow, &error);
              }
            else
              {
                destroy = callDotter(dialogData->blxWindow, dialogData->hspsOnly, dialogData->dotterSSeq, &error);
              }
          }
        else
          {
            destroy = FALSE; /* there was an error, so leave the dialog open */
          }
        
	break;
	
      case GTK_RESPONSE_APPLY:
	widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
	destroy = FALSE;
	break;
	
      default:
	break;
    };

  /* If any errors were found, report them */
  if (error)
    {
      prefixError(error, "Could not start Dotter. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }

  if (destroy)
    {
      /* The dialog is persistent so hide it rather than destroying it. */
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
}


static void onDestroyDotterDialog(GtkWidget *dialog, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  
  if (dialogData)
    {
      if (dialogData->dotterSSeq)
        {
          g_free(dialogData->dotterSSeq);
          dialogData->dotterSSeq = NULL;
        }
      
      g_free(dialogData);
    }
}


/* Called when the auto/manual radio button is toggled. */
static void onRadioButtonToggled(GtkWidget *button, gpointer data)
{
  DotterDialogData *dialogData = (DotterDialogData*)data;
  BlxViewContext *bc = blxWindowGetContext(dialogData->blxWindow);
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(dialogData->autoButton)))
    {
      /* Recalculate auto start/end in case user has selected a different sequence */
      int autoStart = UNSET_INT, autoEnd = UNSET_INT;
      getDotterRange(dialogData->blxWindow, dialogData->dotterSSeq, dialogData->callOnSelf, TRUE, &autoStart, &autoEnd, NULL, NULL);

      if (autoStart == UNSET_INT && autoEnd == UNSET_INT)
        {
          autoStart = bc->displayRev ? bc->refSeqRange.max : bc->refSeqRange.min;
          autoEnd = bc->displayRev ? bc->refSeqRange.min : bc->refSeqRange.max;
        }
      
      char *startString = convertIntToString(getDisplayCoord(autoStart, bc));
      char *endString = convertIntToString(getDisplayCoord(autoEnd, bc));
      char *zoomString = convertIntToString(bc->dotterZoom);

      /* Display the new values */
      gtk_entry_set_text(GTK_ENTRY(dialogData->startEntry), startString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->endEntry), endString);
      gtk_entry_set_text(GTK_ENTRY(dialogData->zoomEntry), zoomString);
      
      g_free(startString);
      g_free(endString);
      g_free(zoomString);
      
      /* Lock out the entry boxes so they cannot be edited */
      gtk_widget_set_sensitive(dialogData->startEntry, FALSE);
      gtk_widget_set_sensitive(dialogData->endEntry, FALSE);
      gtk_widget_set_sensitive(dialogData->zoomEntry, FALSE);
    }
  else
    {
      /* Manual coords. Leave values as they are but unlock the boxes so they can be edited. */
      gtk_widget_set_sensitive(dialogData->startEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->endEntry, TRUE);
      gtk_widget_set_sensitive(dialogData->zoomEntry, TRUE);
    }
}


static GtkWidget* createTextEntry(GtkTable *table, 
				  int col, 
				  int row, 
				  const int xpad, 
				  const int ypad, 
				  const char *title,
				  BlxResponseCallback callbackFunc,
				  GtkWidget *blxWindow,
				  const int initValue)
{
  /* Create a label in the given column */
  GtkWidget *label = gtk_label_new(title);
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create the text entry widget in the next column in the same row */
  ++col;
  GtkWidget *entry = gtk_entry_new();
  gtk_table_attach(table, entry, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  /* Set the initial text */
  char *initText = convertIntToString(initValue);
  gtk_entry_set_text(GTK_ENTRY(entry), initText);
  g_free(initText);
  
  /* Add the callback data. This specifies what callback to use when the user hits OK or Apply on the dialog. */
  widgetSetCallbackData(entry, callbackFunc, blxWindow);
  
  return entry;
}


/* Utility to get the title for the dotter dialog. Uses the selected sequence name if a single
 * sequence is selected, or shows <no sequences> or <multiple sequences>. The result should be
 * free'd with g_free. */
static char* getDotterTitle(const BlxViewContext *bc)
{
  char *result = NULL;
  GString *resultStr = g_string_new(blxGetTitlePrefix(bc));
  g_string_append(resultStr, "Dotter sequence: ");

  const int numSeqs = g_list_length(bc->selectedSeqs);
  
  if (numSeqs == 1)
    {
      const BlxSequence *blxSeq = (const BlxSequence*)(bc->selectedSeqs->data);
      g_string_append(resultStr, blxSequenceGetName(blxSeq));
    }
  else if (numSeqs < 1)
    {
      g_string_append(resultStr, "<no sequences selected>");
    }
  else
    {
      g_string_append(resultStr, "<multiple sequences>");
    }
  
  result = g_string_free(resultStr, FALSE);
  
  return result;
}


/* Pop up a dialog to allow the user to edit dotter parameters and launch dotter */
void showDotterDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_DOTTER;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  static DotterDialogData *dialogData = NULL;
  
  char *title = getDotterTitle(bc);
  
  if (!dialog)
    {
      dialog = gtk_dialog_new_with_buttons(title, 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_SAVE,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_EXECUTE,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      /* Create the dialog data struct first time round, but re-populate it each time. Create 
       * a destructor function that will free the struct. */
      dialogData = (DotterDialogData*)g_malloc(sizeof(DotterDialogData));
      dialogData->blxWindow = NULL;
      dialogData->autoButton = NULL;
      dialogData->manualButton = NULL;
      dialogData->startEntry = NULL;
      dialogData->endEntry = NULL;
      dialogData->zoomEntry = NULL;
      dialogData->callOnSelf = FALSE;
      dialogData->hspsOnly = FALSE;
      dialogData->dotterSSeq = NULL;

      g_signal_connect(G_OBJECT(dialog), "destroy", G_CALLBACK(onDestroyDotterDialog), dialogData);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDotterDialog), dialogData);
    }
  else
    {
      /* Refresh the dialog by clearing its contents an re-creating it */
      dialogClearContentArea(GTK_DIALOG(dialog));
      
      gtk_window_set_title(GTK_WINDOW(dialog), title);
    }
  
  g_free(title);
  
  GtkContainer *contentArea = GTK_CONTAINER(GTK_DIALOG(dialog)->vbox);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  const int spacing = 4;

  /* Create a container for the child widgets */
  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 0));
  gtk_container_add(contentArea, GTK_WIDGET(hbox));

  /* Radio buttons for auto/manual params */
  GtkBox *vbox1 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox1), FALSE, FALSE, spacing);
  
  GtkWidget *autoButton = gtk_radio_button_new_with_mnemonic(NULL, "_Auto");
  gtk_box_pack_start(vbox1, autoButton, FALSE, FALSE, spacing);
  gtk_widget_set_tooltip_text(autoButton, "Always use the visible big picture range");
  widgetSetCallbackData(autoButton, onSaveDotterMode, blxWindow);
		   
  GtkWidget *manualButton = gtk_radio_button_new_with_mnemonic_from_widget(GTK_RADIO_BUTTON(autoButton), "_Manual");
  gtk_widget_set_tooltip_text(manualButton, "Set the coords manually. Once saved, the same coords will be used until you change them manually again.");
  gtk_box_pack_start(vbox1, manualButton, FALSE, FALSE, spacing);

  /* Buttons that the user can click to populate the parameter boxes with certain values */
  GtkWidget *lastSavedButton = gtk_button_new_with_mnemonic("_Last saved ->");
  GtkWidget *fullRangeButton = gtk_button_new_with_mnemonic("_Full range ->");
  GtkWidget *bpRangeButton = gtk_button_new_with_mnemonic("_Big picture range ->");

  gtk_widget_set_tooltip_text(lastSavedButton, "Put the last-saved coords in the manual coords boxes");
  gtk_widget_set_tooltip_text(fullRangeButton, "Put the full range in the manual coords boxes");
  gtk_widget_set_tooltip_text(bpRangeButton, "Put the big picture range in the manual coords boxes");

  gtk_box_pack_start(vbox1, lastSavedButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox1, fullRangeButton, FALSE, FALSE, spacing);
  gtk_box_pack_start(vbox1, bpRangeButton, FALSE, FALSE, spacing);

  /* Disable last-saved button if no saved values exist */
  if (bc->dotterStart == UNSET_INT)
    {
      gtk_widget_set_sensitive(lastSavedButton, FALSE);
    }

  /* Create the text entry boxes for the dotter parameters */
  GtkTable *table = GTK_TABLE(gtk_table_new(3, 2, FALSE));
  gtk_box_pack_start(hbox, GTK_WIDGET(table), FALSE, FALSE, spacing);
  int xpad = 4, ypad = 4;
  
  GtkWidget *startEntry = createTextEntry(table, 1, 1, xpad, ypad, "<i>Start:</i>", onSaveDotterStart, blxWindow, getDisplayCoord(bc->dotterStart, bc));
  GtkWidget *endEntry = createTextEntry(table, 1, 2, xpad, ypad, "<i>End:</i>", onSaveDotterEnd, blxWindow, getDisplayCoord(bc->dotterEnd, bc));
  GtkWidget *zoomEntry = createTextEntry(table, 1, 3, xpad, ypad, "<i>Zoom:</i>", onSaveDotterZoom, blxWindow, bc->dotterZoom);
  gtk_widget_set_tooltip_text(zoomEntry, "The level of zoom to open Dotter with (higher values zoom in)");

  /* Create check buttons for the various options */
  GtkBox *vbox3 = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox3), FALSE, FALSE, spacing);
  
  GtkWidget *selfButton = gtk_check_button_new_with_mnemonic("Call on _self");
  gtk_widget_set_tooltip_text(selfButton, "Call dotter on the reference sequence versus itself (useful for finding internal repeats)");
  gtk_box_pack_start(vbox3, selfButton, FALSE, FALSE, 0);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(selfButton), bc->dotterSelf);
  widgetSetCallbackData(selfButton, onSaveDotterSelf, blxWindow);
  
  GtkWidget *hspsButton = gtk_check_button_new_with_mnemonic("_HSPs only");
  gtk_widget_set_tooltip_text(hspsButton, "Start dotter showing high-scoring pairs only");
  gtk_box_pack_start(vbox3, hspsButton, FALSE, FALSE, 0);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hspsButton), bc->dotterHsps);
  widgetSetCallbackData(hspsButton, onSaveDotterHsps, blxWindow);
  
  /* Create the data struct we need to pass to the toggled callback, and connect signals */
  dialogData->blxWindow = blxWindow;
  dialogData->autoButton = autoButton;
  dialogData->manualButton = manualButton;
  dialogData->startEntry = startEntry;
  dialogData->endEntry = endEntry;
  dialogData->zoomEntry = zoomEntry;
  dialogData->callOnSelf = bc->dotterSelf;
  dialogData->hspsOnly = bc->dotterHsps;
  
  if (dialogData->dotterSSeq)
    {
      g_free(dialogData->dotterSSeq);
      dialogData->dotterSSeq = NULL;
    }

  dialogData->dotterSSeq = getDotterSSeq(blxWindow, NULL);
  
  /* There is an issue if the user selects a different sequence while the dotter dialog
   * is still open: the auto range does not update automatically for the new sequence. To 
   * mitigate this, connect the 'clicked' signal as well as the toggle signal, so that they can
   * click on the 'auto' toggle button and have it refresh, even if that button is already selected.*/
//  g_signal_connect(G_OBJECT(autoButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
//  g_signal_connect(G_OBJECT(manualButton), "toggled", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(autoButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(manualButton), "clicked", G_CALLBACK(onRadioButtonToggled), dialogData);
  
  g_signal_connect(G_OBJECT(lastSavedButton), "clicked", G_CALLBACK(onLastSavedButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(fullRangeButton), "clicked", G_CALLBACK(onFullRangeButtonClicked), dialogData);
  g_signal_connect(G_OBJECT(bpRangeButton), "clicked", G_CALLBACK(onBpRangeButtonClicked), dialogData);

  g_signal_connect(G_OBJECT(selfButton), "toggled", G_CALLBACK(onSelfButtonToggled), dialogData);
  g_signal_connect(G_OBJECT(hspsButton), "toggled", G_CALLBACK(onHspsButtonToggled), dialogData);
  
  /* Set the initial state of the toggle buttons and entry widgets */
  if (bc->autoDotter)
    {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manualButton), TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(autoButton), TRUE);
    }
  else
    {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(autoButton), TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(manualButton), TRUE);
    }
  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/*******************************************************************
 *		      Dotter utility functions                     *
 *******************************************************************/

/* Convert a blixem blast mode to a char for passing to dotter */
char getDotterMode(const BlxBlastMode blastMode)
{
  char type = ' ';
  
  if (blastMode == BLXMODE_BLASTP || blastMode == BLXMODE_TBLASTN)
    {
      type = 'P';
    }
  else if (blastMode == BLXMODE_BLASTX)
    {
      type = 'X';
    }
  else if (blastMode == BLXMODE_BLASTN || blastMode == BLXMODE_TBLASTX)
    {
      type = 'N';
    }
  
  return type;
}


/* Get the start/end coords. If the passed autoDotter flag is true, calculate coords
 * automatically - otherwise use the stored manual coords */
static gboolean getDotterRange(GtkWidget *blxWindow, 
			       const char *dotterSSeq, 
			       const gboolean callOnSelf,
                               const gboolean autoDotter,
			       int *dotterStart, 
			       int *dotterEnd, 
			       int *dotterZoom, 
			       GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */

  GError *tmpError = NULL;
  gboolean success = TRUE;
  
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (!autoDotter)
    {
      /* Use manual coords */
      if (dotterStart) *dotterStart = bc->dotterStart;
      if (dotterEnd) *dotterEnd = bc->dotterEnd;
      if (dotterZoom) *dotterZoom = bc->dotterZoom;
      
      if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
	{
	  g_debug("Manual dotter parameters were requested but one or more coord is not set (start=%d, end=%d). Calculating automatic parameters instead.\n",
		  *dotterStart, *dotterEnd);
	}
    }
  
  if ((dotterStart && *dotterStart == UNSET_INT) || (dotterEnd && *dotterEnd == UNSET_INT))
    {
      /* Calculate automatic coords */
      if (callOnSelf)
	success = smartDotterRangeSelf(blxWindow, dotterStart, dotterEnd, &tmpError);
      else
	success = smartDotterRange(blxWindow, dotterSSeq, dotterStart, dotterEnd, &tmpError);
    }

  if (success && !tmpError)
    {
      /* Check that there are valid MSPs within the dotter range. Set a warning if not. */
      gboolean found = FALSE;
      const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD);
      const int qMin = min(*dotterStart, *dotterEnd);
      const int qMax = max(*dotterStart, *dotterEnd);

      GList *seqItem = bc->selectedSeqs;
      for ( ; seqItem; seqItem = seqItem->next)
        {
          const BlxSequence *seq = (const BlxSequence*)(seqItem->data);
          GList *mspItem = seq->mspList;  
          
          for ( ; mspItem ; mspItem = mspItem->next)
            {
              const MSP *msp = (MSP*)(mspItem->data);
              
              if ((msp->qStrand == activeStrand || bc->blastMode == BLXMODE_BLASTN) &&
                  msp->qRange.min <= qMax && msp->qRange.max >= qMin)
                {
                  found = TRUE;
                  break;
                }
            }
        }
      
      if (!found)
        {
          g_set_error(&tmpError, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_MATCHES,
            "There were no matches for the selected sequence(s) within the selected dotter range.\nThis may be because the range was trimmed - try zooming in on the region of interest.\n");
        }
    }

  if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }

  return success;
}


/* Utility to fetch the selected match sequence or get it from the selected MSP.
 * This function assumes that if multiple MSPs are selected, that they are all for 
 * the same match sequence. Returns null if no MSPs are selected, with details of the error
 * in 'error'.  If the sequence was found but there were warnings, it returns non-null with
 * the warnings in 'error'. The return value should be free'd with g_free */
char* getDotterSSeq(GtkWidget *blxWindow, GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */
  
  char *dotterSSeq = NULL;
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the selected sequence name */
  if (g_list_length(bc->selectedSeqs) < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return dotterSSeq;
    }
  
  const BlxSequence *blxSeq = (const BlxSequence*)(bc->selectedSeqs->data);

  /* If we're in seqbl mode, only part of the sequence is stored
   * internally, so try to fetch the full sequence.
   * gb10: I don't think this is applicable any more (or even if it
   * will work if partial sequences are stored). If we do need to do
   * a fetch here then we will need to look for a fetch method that
   * returns the fasta sequence (rather than the embl entry). */
  const BlxBlastMode blastMode = bc->blastMode;
  /*  if (blastMode != BLXMODE_TBLASTN)
    {
      fetchSequence(blxSeq, FALSE, 0, blxWindow, &dotterSSeq);
      }
  */

  if (!dotterSSeq && blastMode != BLXMODE_TBLASTX)
    {
      /* Check if sequence is stored internally (i.e. it was passed from acedb) */
      g_debug("Looking for sequence stored internally... ");
    
      dotterSSeq = g_strdup(blxSequenceGetSequence(blxSeq));
      
      if (!dotterSSeq)
	{
	  g_debug("not found.\n");
	  g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NOT_FOUND, "Failed to find sequence for '%s'.\n", blxSequenceGetName(blxSeq));
	  return FALSE;
	}

      g_debug("found.\n");
      
      /* Dotter expects the passed sequence to be forwards and un-complemented but if this is
       * the reverse strand the sequence will be complemented, so un-complement it. */
      if (blxSeq->strand == BLXSTRAND_REVERSE)
        {
          blxComplement(dotterSSeq);
        }
    }

  if (dotterSSeq && (strchr(dotterSSeq, SEQUENCE_CHAR_PAD) || blastMode == BLXMODE_TBLASTN))
    {
      g_warning("The sequence for '%s' is incomplete.\n", blxSequenceGetName(blxSeq));
    }
  
  return dotterSSeq;
}


/* Attempts to set the range of dotter in a sensible way, when calling dotter on the reference
 * sequence versus itself. */
static gboolean smartDotterRangeSelf(GtkWidget *blxWindow,
				     int *dotter_start_out,
				     int *dotter_end_out,
				     GError **error)
{
  /* We'll just use a large-ish range centred on the current display range */
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
  
  int mid = getRangeCentre(displayRange);

  /* Convert to DNA coords */
  mid = convertDisplayIdxToDnaIdx(mid, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  
  const int offset = DEFAULT_DOTTER_RANGE_SELF / 2;

  *dotter_start_out = mid - offset;
  *dotter_end_out = *dotter_start_out + DEFAULT_DOTTER_RANGE_SELF;

  boundsLimitValue(dotter_start_out, &bc->refSeqRange);
  boundsLimitValue(dotter_end_out, &bc->refSeqRange);
  
  return TRUE;
}


/* Attempts to set the range of dotter in a sensible way, when calling dotter on the reference
 * sequence versus itself. For now it just uses the big picture range. It's different to setting
 * the big picture range using manual coords because it gets updated whenever the user scrolls,
 * whereas manual coords don't change until the user hits the button again.
 * gb10: this is a much simpler replacement for the ifdef'd out smartDotterRange. We could do
 * something much smarter but at the moment it doesn't seem to be necessary. */
static gboolean smartDotterRange(GtkWidget *blxWindow,
                                 const char *dotterSSeq, 
                                 int *dotter_start_out,
                                 int *dotter_end_out,
                                 GError **error)
{
  /* We'll just use a large-ish range centred on the current display range */
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const IntRange* const displayRange = bigPictureGetDisplayRange(bigPicture);
  
  /* Convert to DNA coords */
  int start = convertDisplayIdxToDnaIdx(displayRange->min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  int end = convertDisplayIdxToDnaIdx(displayRange->max, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
  
  if (dotter_start_out)
    *dotter_start_out = start;
  
  if (dotter_end_out)
    *dotter_end_out = end;

  boundsLimitValue(dotter_start_out, &bc->refSeqRange);
  boundsLimitValue(dotter_end_out, &bc->refSeqRange);
  
  return TRUE;
}


#ifdef USE_OLD_AUTO_DOTTER_COORDS
/* gb10: Users never liked this way of calculating auto dotter coords because it can return a range
 * that is nowhere near the visible range in the big picture. I've changed auto coords to just
 * use the visible big picture range instead for now. I've left the old code here because it's _almost_
 * good - it just needs to do something more sensible when clipping the range - so this might be
 * useful in future if we can improve it.  */

/* Attempts to set the range of dotter in some sort of sensible way. The problem is that
 * hits can occur over a much wider range than the user is looking at, so the function
 * attempts to find the range of hits that corresponds to what the user can see.
 * Returns TRUE if it managed to find sequences and set a sensible range, FALSE otherwise.
 * Returns true but sets the error if there is a warning.
 * NOTE: This function assumes that only a single sequence can be selected at any one time. */
static gboolean smartDotterRange(GtkWidget *blxWindow,
				 const char *dotterSSeq, 
				 int *dotter_start_out, 
				 int *dotter_end_out,
				 GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */

  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Check that a sequence is selected */
  GList *selectedSeqs = bc->selectedSeqs;
  if (g_list_length(selectedSeqs) < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }

  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  const IntRange* const bigPicRange = bigPictureGetDisplayRange(bigPicture);
  const BlxStrand activeStrand = (bc->displayRev ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD) ;

  /* Loop through all MSPs in the selected sequence. We'll estimate the wanted
   * query region from the extent of the HSP's that are completely within view. */
  const BlxSequence *selectedSeq = (const BlxSequence*)(selectedSeqs->data);
  int qMin = UNSET_INT, qMax = UNSET_INT;
  GList *mspListItem = selectedSeq->mspList;  
  
  for ( ; mspListItem ; mspListItem = mspListItem->next)
    {
      const MSP *msp = (MSP*)(mspListItem->data);
      
      /* Get the msp start/end in terms of display coords, and find the min/max */
      const IntRange* const mspDisplayRange = mspGetDisplayRange(msp);

      /* Check if the MSP is in a visible strand is entirely within the big picture range */
      if ((msp->qStrand == activeStrand || (bc->blastMode == BLXMODE_BLASTN)) &&
	  rangesOverlap(mspDisplayRange, bigPicRange))
	{
	  int qSeqMin = msp->qRange.min;
	  int qSeqMax = msp->qRange.max;
	  int sSeqMin = msp->sRange.min;
	  int sSeqMax = msp->sRange.max;
	  
	  /* Extrapolate qMin backwards to the start of the match sequence (i.e. where
	   * s==0) and qMax forwards to the end of the match sequence (i.e. where s==sLength). */
	  int distToSMin = sSeqMin - 1;
	  int distToSMax = 200; /* default amount if sequence not found or if mode is tblastn */
	  
	  if (bc->blastMode != BLXMODE_TBLASTN && dotterSSeq)
	    {
	      distToSMax = strlen(dotterSSeq) - sSeqMax;
	    }

	  /* If the match sequence is a peptide sequence, convert the number of peptide
	   * coords we want to traverse to the equivalent number of DNA coords */
	  if (bc->seqType == BLXSEQ_PEPTIDE)
	    {
	      distToSMin *= bc->numFrames;
	      distToSMax *= bc->numFrames;
	    }

	  /* If the strands are in opposite directions, the low end of the ref 
	   * sequence corresponds to the high of the match sequence, and vice versa. */
	  const gboolean sameDirection = (mspGetRefStrand(msp) == mspGetMatchStrand(msp));
	  qSeqMin -= sameDirection ? distToSMin : distToSMax;
	  qSeqMax += sameDirection ? distToSMax : distToSMin;

	  /* Set the results */
	  if (qMin == UNSET_INT || qSeqMin < qMin)
	    {
	      qMin = qSeqMin;
	    }
	  
	  if (qMax == UNSET_INT || qSeqMax > qMax)
	    {
	      qMax = qSeqMax;
	    }
	}
    }

  if (qMin == UNSET_INT && qMax == UNSET_INT)
    {
      /* No alignments found. Give a warning, and use the big picture range. */
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_MATCHES, 
                  "There were no matches for the selected sequence(s) within the big picture range.\nZoom out to ensure alignments lie entirely within the big picture range.");
      
      qMin = convertDisplayIdxToDnaIdx(bigPicRange->min, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
      qMax = convertDisplayIdxToDnaIdx(bigPicRange->max, bc->seqType, 1, bc->numFrames, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  /* Due to gaps, we might miss the ends - add some more */
  int extend = 0.1 * (qMax - qMin) ;
  qMin -= extend ;
  qMax += extend ;

  if (bc->blastMode == BLXMODE_BLASTX || bc->blastMode == BLXMODE_TBLASTX)
    {
      /* If sstart and send weren't in the end exons, we'll miss those - add some more */
      extend = 0.2 * (qMax - qMin) ;
      qMin -= extend ;
      qMax += extend ;
    }

  /* Keep it within bounds. */
  boundsLimitValue(&qMin, &bc->refSeqRange);
  boundsLimitValue(&qMax, &bc->refSeqRange);

  /* Apply min and max limits:  min 500 residues, max 10 Mb dots. */
  int qLen = qMax - qMin;
  int midCoord = qMin + qLen/2;

  if (qLen < 500)
    {
      qLen = 500;
    }

  const int sLen = blxSequenceGetLength(selectedSeq);
  if (qLen * sLen > 1e7)
    {
      qLen = 1e7 / sLen;
    }

  qMin = midCoord - (qLen / 2) ;
  qMax = midCoord + (qLen / 2) ;

  /* Bounds check again */
  boundsLimitValue(&qMin, &bc->refSeqRange);
  boundsLimitValue(&qMax, &bc->refSeqRange);

  /* Return the start/end. The values start low and end high in normal 
   * left-to-right display, or vice-versa if the display is reversed. */
  *dotter_start_out = bc->displayRev ? qMax : qMin;
  *dotter_end_out = bc->displayRev ? qMin : qMax;
  
  return TRUE;
}
#endif


/*******************************************************************
 *		      Functions to call dotter                     *
 *******************************************************************/


/* This actually executes the dotter child process */
static void callDotterChildProcess(GtkWidget *blxWindow,
                                   const char *dotterBinary, 
				   const int dotterZoom,
                                   const gboolean hspsOnly,
                                   const char *seq1Name,
                                   const IntRange* const seq1Range,
                                   const BlxStrand seq1Strand,
				   const gboolean seq1DisplayRev,
                                   const char *seq2Name,
                                   const IntRange* const seq2Range,
                                   const BlxStrand seq2Strand,
				   const gboolean seq2DisplayRev,
				   int *pipes, 
                                   BlxViewContext *bc)
{
  DEBUG_OUT("callDotterChildProcess\n");

  char *seq1OffsetStr = convertIntToString(seq1Range->min - 1);
  char *seq2OffsetStr = convertIntToString(seq2Range->min - 1);
  char *seq1LenStr = convertIntToString(getRangeLength(seq1Range));
  char *seq2LenStr = convertIntToString(getRangeLength(seq2Range));
  
  /* Create the argument list - start with any options we want to pass */
  GSList *argList = NULL;
  argList = g_slist_append(argList, g_strdup(dotterBinary));
  argList = g_slist_append(argList, g_strdup("-z"));
  argList = g_slist_append(argList, convertIntToString(dotterZoom));
  argList = g_slist_append(argList, g_strdup("-q"));
  argList = g_slist_append(argList, seq1OffsetStr);
  argList = g_slist_append(argList, g_strdup("-s"));
  argList = g_slist_append(argList, seq2OffsetStr);
  argList = g_slist_append(argList, g_strdup("--horizontal-type"));
  argList = g_slist_append(argList, g_strdup("d"));
  argList = g_slist_append(argList, g_strdup("--vertical-type"));

  if (bc->seqType == BLXSEQ_PEPTIDE)
    argList = g_slist_append(argList, g_strdup("p"));
  else
    argList = g_slist_append(argList, g_strdup("d"));

  if (bc->flags[BLXFLAG_ABBREV_TITLE])
    argList = g_slist_append(argList, g_strdup("--abbrev-title-on"));
  else
    argList = g_slist_append(argList, g_strdup("--abbrev-title-off"));

  if (bc->windowColor)
    argList = g_slist_append(argList, g_strdup_printf("--session_colour=%s", bc->windowColor));
  
  if (seq1Strand == BLXSTRAND_REVERSE)      argList = g_slist_append(argList, g_strdup("-r"));
  if (seq2Strand == BLXSTRAND_REVERSE)	    argList = g_slist_append(argList, g_strdup("-v"));
  if (seq1DisplayRev)			    argList = g_slist_append(argList, g_strdup("--reverse-h-display"));
  if (seq2DisplayRev)			    argList = g_slist_append(argList, g_strdup("--reverse-v-display"));
  if (hspsOnly)				    argList = g_slist_append(argList, g_strdup("-H"));
  if (bc->flags[BLXFLAG_NEGATE_COORDS])	    argList = g_slist_append(argList, g_strdup("-N"));

  /* now tell Dotter that we're calling it internally from another SeqTools
   * program, so that it knows to expect piped data */
  argList = g_slist_append(argList, g_strdup("-S"));
  
  /* Now pass the required arguments. These must be in the correct order. */
  argList = g_slist_append(argList, g_strdup(seq1Name));
  argList = g_slist_append(argList, seq1LenStr);
  argList = g_slist_append(argList, g_strdup(seq2Name));
  argList = g_slist_append(argList, seq2LenStr);
  argList = g_slist_append(argList, g_strdup(dotterBinary));

  /* Now pass the screen number - we want to start dotter on the same screen as blixem's main
   * window */
  if (blxWindow && GTK_IS_WINDOW(blxWindow))
    {
      GdkScreen *screen = gtk_window_get_screen(GTK_WINDOW(blxWindow));
      if (screen)
        {
          char *screenStr = convertIntToString(gdk_screen_get_number(screen));
          argList = g_slist_append(argList, g_strdup("--screen"));
          argList = g_slist_append(argList, screenStr);
        }
    }

  /* Terminate list with null */
  argList = g_slist_append(argList, NULL);

  /* Convert the list to an array */
  DEBUG_OUT("Args = ");
  char *args[g_slist_length(argList)];
  GSList *item = argList;
  int i = 0;
  
  for ( ; item; item = item->next)
    {
      char *arg = (char*)(item->data);
      args[i] = arg;
      ++i;
      DEBUG_OUT(", %s", arg  );
    }
  DEBUG_OUT("\n");
    
  DEBUG_OUT("Executing dotter\n");
  dup2(pipes[0], 0);
  close(pipes[0]);
  execv(args[0], args);
  
  exit(1);
}


/* Call dotter as an external process */
gboolean callDotterExternal(GtkWidget *blxWindow,
                            BlxViewContext *bc,
                            int dotterZoom, 
                            const gboolean hspsOnly,
                            const char *seq1Name,
                            IntRange *seq1Range,
                            char *seq1,
                            const BlxStrand seq1Strand,
			    const gboolean seq1DisplayRev,
                            const char *seq2Name,
                            IntRange *seq2Range,
                            char *seq2,
			    const BlxStrand seq2Strand,
			    const gboolean seq2DisplayRev,
                            GError **error)
{
#if !defined(NO_POPEN)

  boundsLimitRange(seq1Range, &bc->refSeqRange, FALSE);

  static char *dotterBinary = NULL;

  /* Open pipe to new dotterBinary */
  if (!dotterBinary) 
    { 
      g_debug("Looking for Dotter ...\n");
      dotterBinary = g_find_program_in_path("dotter");

      if (!dotterBinary)
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_EXE, "No dotter executable found in path.\n$PATH=%s\n", getenv("PATH"));
          dotterBinary = NULL;
          return FALSE;
        }
    }
  
  g_debug("Calling %s with region: %d,%d - %d,%d\n", dotterBinary, seq1Range->min, seq2Range->min, seq1Range->max, seq2Range->max);
  fflush(stdout);

  /* Pass on the features via pipes. */
  int pipes[2];
  if (pipe (pipes))
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_PIPE, "Error opening pipe to Dotter.\n");
      return FALSE;
    }

  /* Create the child process */
  pid_t pid = fork();
  bc->spawnedProcesses = g_slist_append(bc->spawnedProcesses, GINT_TO_POINTER(pid));
  
  if (pid < 0)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_FORK, "Error forking process for Dotter.\n");
      return FALSE;
    }
  else if (pid == 0)
    {
      /* Child process. Execute dotter. First, close the write end of the pipe */
      close(pipes[1]);

      DEBUG_OUT("Calling dotter child process\n");
      callDotterChildProcess(blxWindow, dotterBinary, dotterZoom, hspsOnly, 
                             seq1Name, seq1Range, seq1Strand, seq1DisplayRev,
                             seq2Name, seq2Range, seq2Strand, seq2DisplayRev,
                             pipes, bc);
    }
  else
    {
      /* Parent process. Pipe sequences and MSPs to dotter. First, close the read end of the pipe. */
      close(pipes[0]);

      g_debug("Spawned process %d\n", pid);
      DEBUG_OUT("Piping sequences to dotter...\n");
      
      fflush(stdout);
      FILE *pipe = fdopen(pipes[1], "w");

      /* Pass the sequences */
      fwrite(seq1, 1, getRangeLength(seq1Range), pipe);
      fwrite(seq2, 1, getRangeLength(seq2Range), pipe);
      
      DEBUG_OUT("...done\n");
      
      /* Pass the features */
      DEBUG_OUT("Piping features to dotter...\n");
      GList *seqItem = bc->matchSeqs;
      for ( ; seqItem; seqItem = seqItem->next) 
	{
	  BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
	  writeBlxSequenceToOutput(pipe, blxSeq, seq1Range, seq2Range);
	}
      
      fprintf(pipe, "%c\n", EOF);
      fflush(pipe);
      
      DEBUG_OUT("...done\n");
    }
#endif
  
  return TRUE;
}


/* Call dotter. Returns true if dotter was called; false if we quit trying. */
gboolean callDotter(GtkWidget *blxWindow, const gboolean hspsOnly, char *dotterSSeqIn, GError **error)
{
  g_return_val_if_fail(!error || *error == NULL, FALSE); /* if error is passed, its contents must be NULL */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const int numSeqsSelected = g_list_length(bc->selectedSeqs);
  
  if (numSeqsSelected < 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQS, "There are no sequences selected.\n");
      return FALSE;
    }
  else if (numSeqsSelected > 1)
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_TOO_MANY_SEQS, "Cannot run dotter on multiple sequences.\n");
      return FALSE;
    }
  
  /* Check this sequence is a valid blast match (just check the first MSP;
   * they must all the same type if they have the same seq name) */
  const BlxSequence *selectedSeq = (const BlxSequence*)(bc->selectedSeqs->data);
  const MSP *firstMsp = (const MSP*)(selectedSeq->mspList->data);

  if (!mspIsBlastMatch(firstMsp))
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_INVALID_SEQ, "You must select a valid match sequence first.\n");
      return FALSE;
    }
  
  /* We will display the active strand as the main strand in dotter */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);

  if (bc->seqType == BLXSEQ_PEPTIDE)
    {
      GList *mspItem = selectedSeq->mspList;
      gboolean found = FALSE;
      
      for ( ; mspItem && !found; mspItem = mspItem->next)
        {
          const MSP* const msp = (const MSP*)(mspItem->data);
          found = (msp->qStrand == qStrand);
        }
      
      if (!found)
        {
          g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_INVALID_STRAND, "You must select a match sequence on the active strand, or toggle strands first.\n");
          return FALSE;
        }
    }
  
  /* Make a copy of the match sequence, because dotter takes ownership of this. */
  char *dotterSSeq = NULL;
  if (dotterSSeqIn)
    {
      dotterSSeq = g_strdup(dotterSSeqIn);
    }
  else
    {
      g_set_error(error, BLX_DOTTER_ERROR, BLX_DOTTER_ERROR_NO_SEQ_DATA, "No sequence data for this sequence.\n");
      return FALSE;
    }
  
  /* Get the coords */
  int dotterStart = UNSET_INT, dotterEnd = UNSET_INT, dotterZoom = 0;
  GError *rangeError = NULL;
  
  gboolean ok = getDotterRange(blxWindow, dotterSSeq, FALSE, bc->autoDotter, &dotterStart, &dotterEnd, &dotterZoom, &rangeError);

  if (!ok)
    {
      prefixError(rangeError, "Error calculating dotter range. ");
      g_propagate_error(error, rangeError);
      return FALSE;
    }
  else if (ok && rangeError && error) /* if error is null, don't issue warnings */
    {
      /* There was a warning when calculating the range. Ask the user if they want to continue. */
      prefixError(rangeError, "Warning: ");
      postfixError(rangeError, "\nContinue?");

      char *title = g_strdup_printf("%sWarning", blxGetTitlePrefix(bc));
      ok = (runConfirmationBox(blxWindow, title, rangeError->message) == GTK_RESPONSE_ACCEPT);
      
      g_free(title);
      g_error_free(rangeError);
      rangeError = NULL;
      
      if (!ok)
	return FALSE;
    }
  
  /* Get the section of reference sequence that we're interested in */
  const int frame = mspGetRefFrame(firstMsp, bc->seqType);
  IntRange dotterRange;
  intrangeSetValues(&dotterRange, dotterStart, dotterEnd);
  GError *seqError = NULL;

  char *refSeqSegment = getSequenceSegment(bc->refSeq,
                                           &dotterRange,
                                           BLXSTRAND_FORWARD,   /* always pass forward strand to dotter */
                                           BLXSEQ_DNA,	  /* calculated dotter coords are always nucleotide coords */
                                           BLXSEQ_DNA,          /* required sequence is in nucleotide coords */
                                           frame,
                                           bc->numFrames,
                                           &bc->refSeqRange,
                                           bc->blastMode,
                                           bc->geneticCode,
                                           FALSE,		  /* input coords are always left-to-right, even if display reversed */
                                           FALSE,               /* always pass forward strand to dotter */
                                           FALSE,               /* always pass forward strand to dotter */
                                           &seqError);
  
  if (!refSeqSegment)
    {
      g_propagate_error(error, seqError);
      return FALSE;
    }
  else
    {
      /* If there was an error set but the sequence was still returned then it's a non-critical warning */
      reportAndClearIfError(&seqError, G_LOG_LEVEL_WARNING);
    }
  
  /* Get the match sequence name (chopping off the letters before the colon, if there is one). */
  const char *dotterSName = strchr(mspGetSName(firstMsp), ':');
  if (dotterSName)
    {
      dotterSName++;
    }
  else
    {
      dotterSName = mspGetSName(firstMsp);
    }

  IntRange sRange = {1, strlen(dotterSSeq)};
  
  const int offset = dotterRange.min - 1;
  const BlxStrand refSeqStrand = blxWindowGetActiveStrand(blxWindow);
  
  /* Get the list of all MSPs */
  g_message("Calling dotter on match '%s' with reference sequence region: %d -> %d\n", dotterSName, dotterStart, dotterEnd);
  
  g_debug("reference sequence: name =  %s, offset = %d\n"
          "    match sequence: name =  %s, offset = %d\n", 
          bc->refSeqName, offset, dotterSName, 0);
  
  const gboolean revHozScale = (refSeqStrand == BLXSTRAND_REVERSE);
  const gboolean revVertScale = FALSE; /* don't rev match seq scale, because it would show in dotter with -ve coords, but blixem always shows +ve coords */

  return callDotterExternal(blxWindow, bc, dotterZoom, hspsOnly, 
                            bc->refSeqName, &dotterRange, refSeqSegment, refSeqStrand, revHozScale,
                            dotterSName, &sRange, dotterSSeq, selectedSeq->strand, revVertScale,
                            error);
}


/* Call dotter on the reference sequence versus itself.
 *
 * The follow notes on this are from http://sonnhammer.sbc.su.se/Dotter.html:
 *
 * When looking for overlaps between many sequences, for instance when assembling contigs, it can
 * be uselful to make a dotplot of all sequences vs. each other. This way any overlap will show up 
 * as a diagonal in the corner of a subsequence dotplot. Dotter has a built-in mechanism for this. 
 * To run Dotter on many sequences at once, concatenate the sequence files (in fasta format). Then
 * run dotter on the concatenated sequence file against itself, and green partitioning lines will
 * appear between the sequences. At each partitioning line, the name of the following sequence is 
 * printed. These lines can be turned on and off with the button "Draw lines a segment ends" in
 * the "Feature series selection tool", which is launched from the main menu. 
 */
static gboolean callDotterSelf(GtkWidget *blxWindow, GError **error)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Get the auto range, if requested */
  int dotterStart = UNSET_INT;
  int dotterEnd = UNSET_INT;
  int dotterZoom = 0;
  
  GError *tmpError = NULL;
  if (!getDotterRange(blxWindow, NULL, TRUE, bc->autoDotter, &dotterStart, &dotterEnd, &dotterZoom, &tmpError))
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  
  /* Set the type */
  char type = ' ';
  
  if (bc->blastMode == BLXMODE_BLASTP || bc->blastMode == BLXMODE_TBLASTN || bc->blastMode == BLXMODE_TBLASTX)
    {
      type = 'P';
    }
  else if (bc->blastMode == BLXMODE_BLASTX || bc->blastMode == BLXMODE_BLASTN)
    {
      type = 'N';
    }

  /* Get the section of reference sequence that we're interested in */
  const BlxStrand qStrand = blxWindowGetActiveStrand(blxWindow);
  const int frame = 1;
  IntRange qRange;
  intrangeSetValues(&qRange, dotterStart, dotterEnd);
    
  char *refSeqSegment = getSequenceSegment(bc->refSeq,
                                           &qRange,
                                           BLXSTRAND_FORWARD,
                                           BLXSEQ_DNA,	  /* calculated dotter coords are always in nucleotide coords */
                                           BLXSEQ_DNA,      /* required sequence is in nucleotide coords */
                                           frame,
                                           bc->numFrames,
                                           &bc->refSeqRange,
                                           bc->blastMode,
                                           bc->geneticCode,
                                           FALSE,		  /* input coords are always left-to-right, even if display reversed */
                                           FALSE,  /* whether to reverse */
                                           FALSE,  /* whether to allow rev strands to be complemented */
                                           &tmpError);

  if (!refSeqSegment)
    {
      g_propagate_error(error, tmpError);
      return FALSE;
    }
  else
    {
      /* If there's an error but the sequence was still returned it's a non-critical warning */
      reportAndClearIfError(&tmpError, G_LOG_LEVEL_WARNING);
    }
  
  /* Make a copy of the reference sequence segment to pass as the match sequence */
  char *dotterSSeq = g_strdup(refSeqSegment);

  const gboolean revScale = (qStrand == BLXSTRAND_REVERSE);
  
  g_message("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);

  callDotterExternal(blxWindow, bc, dotterZoom, FALSE,
                     bc->refSeqName, &qRange, refSeqSegment, qStrand, revScale,
                     bc->refSeqName, &qRange, dotterSSeq, qStrand, revScale,
                     error);

  return TRUE;
}

