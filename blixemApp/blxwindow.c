/*  File: blxWindow.c
 *  Author: Gemma Barson, 2009-11-24
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
 * Description: See blxWindow.h
 *----------------------------------------------------------------------------
 */

#include <blixemApp/blxwindow.h>
#include <blixemApp/detailview.h>
#include <blixemApp/detailviewtree.h>
#include <blixemApp/bigpicture.h>
#include <blixemApp/blxdotter.h>
#include <blixemApp/exonview.h>
#include <blixemApp/coverageview.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxGff3Parser.h>
#include <seqtoolsUtils/blxmsp.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>
#include <ctype.h>

#define DEFAULT_WINDOW_BORDER_WIDTH      1    /* used to change the default border width around the blixem window */
#define DEFAULT_COVERAGE_VIEW_BORDER     12   /* size of border to allow around the coverage view */
#define DEFAULT_FONT_SIZE_ADJUSTMENT     -2   /* used to start with a smaller font than the default widget font */
#define DEFAULT_SCROLL_STEP_INCREMENT    5    /* how many bases the scrollbar scrolls by for each increment */
#define DEFAULT_WINDOW_WIDTH_FRACTION    0.9  /* what fraction of the screen size the blixem window width defaults to */
#define DEFAULT_WINDOW_HEIGHT_FRACTION   0.6  /* what fraction of the screen size the blixem window height defaults to */
#define MATCH_SET_GROUP_NAME             "Match set"
#define LOAD_DATA_TEXT                   "Load optional\ndata"
#define DEFAULT_TABLE_XPAD               2    /* default x-padding to use in tables */
#define DEFAULT_TABLE_YPAD               2    /* default y-padding to use in tables */
#define MAX_RECOMMENDED_COPY_LENGTH      100000 /* warn if about to copy text longer than this to the clipboard */


typedef enum {SORT_TYPE_COL, SORT_TEXT_COL, N_SORT_COLUMNS} SortColumns;


/* Utility struct used when comparing sequences to a search string */
typedef struct _CompareSeqData
  {
    const char *searchStr;    /* the string to search for */
    BlxColumnId searchCol;    /* the column ID, which defines what data to search e.g. Name or Tissue Type */
    BlxViewContext *bc;       /* the main context */
    GList *matchList;         /* resulting list of all BlxSequences that match */
    GError *error;
  } SeqSearchData;


/* Properties specific to the blixem window */
typedef struct _BlxWindowProperties
  {
    GtkWidget *bigPicture;          /* The top section of the view, showing a "big picture" overview of the alignments */
    GtkWidget *detailView;          /* The bottom section of the view, showing a detailed list of the alignments */
    GtkWidget *mainmenu;            /* The main menu */
    GtkActionGroup *actionGroup;    /* The action-group for the menus */

    BlxViewContext *blxContext;       /* The blixem view context */

    GtkPageSetup *pageSetup;          /* Page setup for printing */
    GtkPrintSettings *printSettings;  /* Used so that we can re-use the same print settings as a previous print */
  } BlxWindowProperties;


/* Local function declarations */
static BlxWindowProperties*       blxWindowGetProperties(GtkWidget *widget);

static void                       onHelpMenu(GtkAction *action, gpointer data);
static void                       onAboutMenu(GtkAction *action, gpointer data);
static void                       onQuit(GtkAction *action, gpointer data);
static void                       onPrintMenu(GtkAction *action, gpointer data);
static void                       onPageSetupMenu(GtkAction *action, gpointer data);
static void                       onSettingsMenu(GtkAction *action, gpointer data);
static void                       onLoadMenu(GtkAction *action, gpointer data);
static void                       onCopySeqsMenu(GtkAction *action, gpointer data);
static void                       onCopySeqDataMenu(GtkAction *action, gpointer data);
static void                       onCopyRefSeqMenu(GtkAction *action, gpointer data);
static void                       onSortMenu(GtkAction *action, gpointer data);
static void                       onZoomInMenu(GtkAction *action, gpointer data);
static void                       onZoomOutMenu(GtkAction *action, gpointer data);
static void                       onFindMenu(GtkAction *action, gpointer data);
static void                       onGoToMenu(GtkAction *action, gpointer data);
static void                       onPrevMatchMenu(GtkAction *action, gpointer data);
static void                       onNextMatchMenu(GtkAction *action, gpointer data);
static void                       onFirstMatchMenu(GtkAction *action, gpointer data);
static void                       onLastMatchMenu(GtkAction *action, gpointer data);
static void                       onPageLeftMenu(GtkAction *action, gpointer data);
static void                       onPageRightMenu(GtkAction *action, gpointer data);
static void                       onScrollLeft1Menu(GtkAction *action, gpointer data);
static void                       onScrollRight1Menu(GtkAction *action, gpointer data);
static void                       onSquashMatchesMenu(GtkAction *action, gpointer data);
static void                       onToggleStrandMenu(GtkAction *action, gpointer data);
static void                       onViewMenu(GtkAction *action, gpointer data);
static void                       onCreateGroupMenu(GtkAction *action, gpointer data);
static void                       onEditGroupsMenu(GtkAction *action, gpointer data);
static void                       onToggleMatchSet(GtkAction *action, gpointer data);
static void                       onDotterMenu(GtkAction *action, gpointer data);
static void                       onCloseAllDottersMenu(GtkAction *action, gpointer data);
static void                       onSelectFeaturesMenu(GtkAction *action, gpointer data);
static void                       onDeselectAllRows(GtkAction *action, gpointer data);
static void                       onStatisticsMenu(GtkAction *action, gpointer data);

static gboolean                   onKeyPressBlxWindow(GtkWidget *window, GdkEventKey *event, gpointer data);
static void                       onUpdateBackgroundColor(GtkWidget *blxWindow);

static void                       onDestroyBlxWindow(GtkWidget *widget);

static BlxStrand                  blxWindowGetInactiveStrand(GtkWidget *blxWindow);

static GtkComboBox*               widgetGetComboBox(GtkWidget *widget);
static BlxColumnId                getColumnFromComboBox(GtkComboBox *combo);

static void                       onButtonClickedDeleteGroup(GtkWidget *button, gpointer data);
static void                       blxWindowGroupsChanged(GtkWidget *blxWindow);
static void                       getSequencesThatMatch(gpointer listDataItem, gpointer data);
static GList*                     getSeqStructsFromText(GtkWidget *blxWindow, const char *inputText, const BlxColumnId searchCol, GError **error);

static void                       createSortBox(GtkBox *parent, GtkWidget *detailView, const BlxColumnId initSortColumn, GList *columnList, const char *labelText, const gboolean searchableOnly);
static GtkWidget*                 createCheckButton(GtkBox *box, const char *mnemonic, const gboolean isActive, GCallback callback, gpointer data);
static void                       blxWindowSetUsePrintColors(GtkWidget *blxWindow, const gboolean usePrintColors);
static gboolean                   blxWindowGetUsePrintColors(GtkWidget *blxWindow);

static void                       blxWindowFindDnaString(GtkWidget *blxWindow, const char *inputSearchStr, const int startCoord, const gboolean searchLeft, const gboolean findAgain, GError **error);
static GList*                     findSeqsFromList(GtkWidget *blxWindow, const char *inputText, const BlxColumnId inputCol, const gboolean rememberSearch, const gboolean findAgain, GError **error);
static int                        getSearchStartCoord(GtkWidget *blxWindow, const gboolean startBeginning, const gboolean searchLeft);
static GList*                     findSeqsFromColumn(GtkWidget *blxWindow, const char *inputText, const BlxColumnId searchCol, const gboolean rememberSearch, const gboolean findAgain, GError **error);
static GtkWidget*                 dialogChildGetBlxWindow(GtkWidget *child);
static void                       killAllSpawned(BlxViewContext *bc);
static void                       calculateDepth(BlxViewContext *bc);
static gdouble                    calculateMspData(MSP *mspList, BlxViewContext *bc);

static gboolean                   setFlagFromButton(GtkWidget *button, gpointer data);
static void                       copySelectedSeqDataToClipboard(GtkWidget *blxWindow);
static void                       copyRefSeqToClipboard(GtkWidget *blxWindow, const int fromIdx_in, const int toIdx_in);


/* MENU BUILDERS */

/* Standard menu entries */
static const GtkActionEntry mainMenuEntries[] = {
  { "CopyMenuAction",   NULL, "Copy"},

  { "Quit",             GTK_STOCK_QUIT,           "_Quit",                    "<control>Q",         "Quit  Ctrl+Q",                         G_CALLBACK(onQuit)},
  { "Help",             GTK_STOCK_HELP,           "_Help",                    "<control>H",         "Display help  Ctrl+H",                 G_CALLBACK(onHelpMenu)},
  { "About",            GTK_STOCK_ABOUT,          "About",                    NULL,                 "Program information",                  G_CALLBACK(onAboutMenu)},
  { "Print",            GTK_STOCK_PRINT,          "_Print...",                "<control>P",         "Print  Ctrl+P",                        G_CALLBACK(onPrintMenu)},
  { "PageSetup",        GTK_STOCK_PAGE_SETUP,     "Page set_up...",           NULL,                 "Page setup",                           G_CALLBACK(onPageSetupMenu)},
  { "Settings",         GTK_STOCK_PREFERENCES,    "_Settings...",             "<control>S",         "Settings  Ctrl+S",                     G_CALLBACK(onSettingsMenu)},
  { "Load",             GTK_STOCK_OPEN,           "_Open features file...",    NULL,                "Load additional features from file  Ctrl+O", G_CALLBACK(onLoadMenu)},

  { "CopySeqNames",     NULL,                     "Copy sequence name(s)",    "<control>C",         "Copy selected sequences name(s)  Ctrl+C", G_CALLBACK(onCopySeqsMenu)},
  { "CopySeqData",      NULL,                     "Copy sequence data",       "<shift><control>C",  "Copy selected sequences data  Shift+Ctrl+C", G_CALLBACK(onCopySeqDataMenu)},
  { "CopyRefSeq",       NULL,                     "Copy reference sequence from mark",  NULL,                 "Copy selected sequences name(s)  Alt+C", G_CALLBACK(onCopyRefSeqMenu)},

  { "Sort",             GTK_STOCK_SORT_ASCENDING, "Sort...",                  NULL,                 "Sort sequences",                       G_CALLBACK(onSortMenu)},
  { "ZoomIn",           GTK_STOCK_ZOOM_IN,        "Zoom in",                  "equal",              "Zoom in  =",                           G_CALLBACK(onZoomInMenu)},
  { "ZoomOut",          GTK_STOCK_ZOOM_OUT,       "Zoom out",                 "minus",              "Zoom out  -",                          G_CALLBACK(onZoomOutMenu)},
  { "GoTo",             GTK_STOCK_JUMP_TO,        "Go to position...",        "P",                  "Go to position  P",                    G_CALLBACK(onGoToMenu)},
  { "FirstMatch",       GTK_STOCK_GOTO_FIRST,     "First match",              "<control>Home",      "Go to first match in selection (or all, if none selected)  Ctrl+Home",    G_CALLBACK(onFirstMatchMenu)},
  { "PrevMatch",        GTK_STOCK_GO_BACK,        "Previous match",           "<control>Left",      "Go to previous match in selection (or all, if none selected)  Ctrl+Left", G_CALLBACK(onPrevMatchMenu)},
  { "NextMatch",        GTK_STOCK_GO_FORWARD,     "Next match",               "<control>Right",     "Go to next match in selection (or all, if none selected)  Ctrl+Right",    G_CALLBACK(onNextMatchMenu)},
  { "LastMatch",        GTK_STOCK_GOTO_LAST,      "Last match",               "<control>End",       "Go to last match in selection (or all, if none selected)  Ctrl+End",      G_CALLBACK(onLastMatchMenu)},
  { "BackPage",         NULL,                     "<<",                       "<control>comma",     "Scroll left one page  Ctrl+,",         G_CALLBACK(onPageLeftMenu)},
  { "BackOne",          NULL,                     "<",                        "comma",              "Scroll left one index  ,",             G_CALLBACK(onScrollLeft1Menu)},
  { "FwdOne",           NULL,                     ">",                        "period",             "Scroll right one index  .",            G_CALLBACK(onScrollRight1Menu)},
  { "FwdPage",          NULL,                     ">>",                       "<control>period",    "Scroll right one page  Ctrl+.",        G_CALLBACK(onPageRightMenu)},
  { "Find",             GTK_STOCK_FIND,           "Find...",                  "<control>F",         "Find sequences  Ctrl+F",               G_CALLBACK(onFindMenu)},
  { "ToggleStrand",     GTK_STOCK_REFRESH,        "Toggle strand",            "T",                  "Toggle the active strand  T",          G_CALLBACK(onToggleStrandMenu)},

  { "View",             GTK_STOCK_FULLSCREEN,     "_View...",                 "V",                  "Edit view settings  V",                G_CALLBACK(onViewMenu)},
  { "CreateGroup",      NULL,                     "Create Group...",          "<shift><control>G",  "Create group  Shift+Ctrl+G",           G_CALLBACK(onCreateGroupMenu)},
  { "EditGroups",       GTK_STOCK_EDIT,           "Edit _Groups...",          "<control>G",         "Edit groups  Ctrl+G",                  G_CALLBACK(onEditGroupsMenu)},
  { "ToggleMatchSet",   NULL,                     "Toggle _match set group",  "G",                  "Create/clear the match set group  G",  G_CALLBACK(onToggleMatchSet)},
  { "DeselectAllRows",  NULL,                     "Deselect _all",            "<shift><control>A",  "Deselect all  Shift+Ctrl+A",           G_CALLBACK(onDeselectAllRows)},

  { "Dotter",           NULL,                     "_Dotter...",               "<control>D",         "Start Dotter  Ctrl+D",                 G_CALLBACK(onDotterMenu)},
  { "CloseAllDotters",  GTK_STOCK_CLOSE,          "Close all Dotters",        NULL,                 "Close all Dotters",                    G_CALLBACK(onCloseAllDottersMenu)},
  { "SelectFeatures",   GTK_STOCK_SELECT_ALL,     "Feature series selection tool...",  NULL,           "Feature series selection tool",        G_CALLBACK(onSelectFeaturesMenu)},

  { "Statistics",       NULL,                     "Statistics",               NULL,                 "Show memory statistics",               G_CALLBACK(onStatisticsMenu)}
};


/* Menu entries for toggle-able actions */
static GtkToggleActionEntry toggleMenuEntries[] = {
  { "SquashMatches",    GTK_STOCK_DND_MULTIPLE,  "Squash matches",           NULL,                 "Squash matches", G_CALLBACK(onSquashMatchesMenu), FALSE} /* must be item 0 in list */
};


/* This defines the layout of the menu for a standard user */
static const char standardMenuDescription[] =
"<ui>"
"  <popup name='ContextMenu' accelerators='true'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Help'/>"
"      <menuitem action='Print'/>"
//"      <menuitem action='PageSetup'/>"
"      <menuitem action='Settings'/>"
"      <menuitem action='Load'/>"
"      <separator/>"
"      <menu action='CopyMenuAction'>"
"        <menuitem action='CopySeqNames'/>"
"        <menuitem action='CopySeqData'/>"
"        <menuitem action='CopyRefSeq'/>"
"      </menu>"
"      <menuitem action='View'/>"
"      <menuitem action='CreateGroup'/>"
"      <menuitem action='EditGroups'/>"
"      <menuitem action='ToggleMatchSet'/>"
"      <menuitem action='DeselectAllRows'/>"
"      <separator/>"
"      <menuitem action='Dotter'/>"
"      <menuitem action='CloseAllDotters'/>"
"  </popup>"
"  <toolbar name='Toolbar'>"
"    <toolitem action='Help'/>"
"    <toolitem action='About'/>"
"    <toolitem action='Settings'/>"
"    <separator/>"
"    <toolitem action='Sort'/>"
"    <toolitem action='SquashMatches'/>"
"    <toolitem action='ZoomIn'/>"
"    <toolitem action='ZoomOut'/>"
"    <separator/>"
"    <toolitem action='GoTo'/>"
"    <toolitem action='FirstMatch'/>"
"    <toolitem action='PrevMatch'/>"
"    <toolitem action='NextMatch'/>"
"    <toolitem action='LastMatch'/>"
"    <toolitem action='BackPage'/>"
"    <toolitem action='BackOne'/>"
"    <toolitem action='FwdOne'/>"
"    <toolitem action='FwdPage'/>"
"    <separator/>"
"    <toolitem action='Find'/>"
"    <toolitem action='ToggleStrand'/>"
"  </toolbar>"
"</ui>";


/* This defines the additional menu components for a developer user */
static const char developerMenuDescription[] =
"<ui>"
"  <popup name='ContextMenu' accelerators='true'>"
"      <menuitem action='SelectFeatures'/>"
"      <separator/>"
"      <menuitem action='Statistics'/>"
"  </popup>"
"</ui>";



/***********************************************************
 *                         Utilities                       *
 ***********************************************************/

/* Return true if the current user is in our list of developers. */
static gboolean userIsDeveloper()
{
  const gchar* developers[] = {"edgrif", "gb10"};

  gboolean result = FALSE;
  const gchar *user = g_get_user_name();
  int numDevelopers = sizeof(developers) / sizeof(gchar*);
  
  int i = 0;
  for (i = 0; i < numDevelopers; ++i)
    {
      if (strcmp(user, developers[i]) == 0)
        {
          result = TRUE;
          break;
        }
    }

  return result;
}


/* Returns the order number of the group that this sequence belongs to, or
 * UNSET_INT if it does not belong to a group. */
int sequenceGetGroupOrder(GtkWidget *blxWindow, const BlxSequence *seq)
{
  SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, seq);
  return group ? group->order : UNSET_INT;
}


/* Scroll the detail view left/right by 1 base (or by 1 page, if the modifier
 * is pressed) */
static void scrollDetailView(GtkWidget *window, const gboolean moveLeft, const gboolean modifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(window);
  
  if (moveLeft && modifier)
    scrollDetailViewLeftPage(detailView);
  else if (moveLeft)
    scrollDetailViewLeft1(detailView);
  else if (modifier)
    scrollDetailViewRightPage(detailView);
  else
    scrollDetailViewRight1(detailView);
}


/* Move the current row selection up/down */
static gboolean moveRowSelection(GtkWidget *blxWindow, const gboolean moveUp, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  const int activeFrame = detailViewGetActiveFrame(detailView);
  const BlxStrand activeStrand = detailViewGetSelectedStrand(detailView);
  
  GtkWidget *tree = detailViewGetTree(detailView, activeStrand, activeFrame);
  return treeMoveRowSelection(tree, moveUp, shiftModifier);
}


/* Move the selected base index 1 base to the left/right. Moves by individual
 * DNA bases (i.e. you have to move 3 bases in order to scroll a full peptide
 * if viewing protein matches). Scrolls the detail view if necessary to keep 
 * the new base in view. */
static void moveSelectedBaseIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  GtkWidget *detailView = blxWindowGetDetailView(window);
  DetailViewProperties *properties = detailViewGetProperties(detailView);

  const gboolean displayRev = detailViewGetDisplayRev(detailView);
  
  const int direction = (moveLeft == displayRev ? 1 : -1);
  const int newDnaIdx = properties->selectedDnaBaseIdx + direction;
  
  detailViewSetSelectedDnaBaseIdx(detailView, 
                                  newDnaIdx, 
                                  detailViewGetActiveFrame(detailView),
                                  TRUE, TRUE);
}


/* Called when user pressed Home/End. If the modifier is pressed, scroll to the
*  start/end of all matches in the current selection (or all matches, if no 
* selection), or to the start/end of the entire display if the modifier is not pressed. */
static void scrollToExtremity(GtkWidget *blxWindow, const gboolean moveLeft, const gboolean modifier)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  if (modifier)
    {
      GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);

      if (moveLeft)
        firstMatch(detailView, selectedSeqs);
      else
        lastMatch(detailView, selectedSeqs);
    }
  else
    {
      const BlxSeqType seqType = blxWindowGetSeqType(blxWindow);
      const IntRange* const fullRange = blxWindowGetFullRange(blxWindow);

      if (moveLeft)
        setDetailViewStartIdx(detailView, fullRange->min, seqType);
      else
        setDetailViewEndIdx(detailView, fullRange->max, seqType);
    }
}


/* Jump left or right to the next/prev nearest match. Only include matches in the
 * current selection, if any rows are selected. */
static void goToMatch(GtkWidget *blxWindow, const gboolean moveLeft)
{
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (moveLeft)
    {
      prevMatch(blxWindowGetDetailView(blxWindow), selectedSeqs);
    }
  else
    {
      nextMatch(blxWindowGetDetailView(blxWindow), selectedSeqs);      
    }
}


/* Move the selected display index 1 value to the left/right. Moves by full peptides
 * if viewing protein matches. Scrolls the detail view if necessary to keep the new 
 * index in view. */
static void moveSelectedDisplayIdxBy1(GtkWidget *window, const gboolean moveLeft)
{
  GtkWidget *detailView = blxWindowGetDetailView(window);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(detailView);
  
  if (detailViewProperties->selectedBaseIdx != UNSET_INT)
    {
      /* Decrement the index if moving left decrease and increment it if moving right, 
       * unless the display is toggled, in which case do the opposite */
      int newSelectedBaseIdx = detailViewProperties->selectedBaseIdx;
      
      if (moveLeft)
        {
          --newSelectedBaseIdx;
        }
      else
        {
          ++newSelectedBaseIdx;
        }
      
      IntRange *fullRange = blxWindowGetFullRange(window);
      boundsLimitValue(&newSelectedBaseIdx, fullRange);
      
      detailViewSetSelectedBaseIdx(detailView, newSelectedBaseIdx, detailViewProperties->selectedFrame, detailViewProperties->selectedBaseNum, TRUE, TRUE);
      detailViewRedrawAll(detailView);
    }
}


/* Zooms the display in/out. The modifiers control which section is zoomed */
static void zoomBlxWindow(GtkWidget *window, const gboolean zoomIn, const gboolean ctrl, const gboolean shift)
{
  if (ctrl)
    {
      if (shift)
        {
          zoomWholeBigPicture(blxWindowGetBigPicture(window));
        }
      else
        {
          zoomBigPicture(blxWindowGetBigPicture(window), zoomIn);
        }
    }
  else
    {
      zoomDetailView(blxWindowGetDetailView(window), zoomIn);
    }
}


/* Force a redraw of all widgets. Clears cached bitmaps etc. first */
void blxWindowRedrawAll(GtkWidget *blxWindow)
{
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  bigPictureRedrawAll(bigPicture);

  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  detailViewRefreshAllHeaders(detailView);
  
  callFuncOnAllDetailViewTrees(detailView, widgetClearCachedDrawable, NULL);
  
  gtk_widget_queue_draw(blxWindow);
}


/* Utility to create a vbox with the given border and pack it into the given box.
 * Also put a frame around it with the given label if includeFrame is true */
static GtkWidget* createVBoxWithBorder(GtkWidget *parent, 
                                       const int borderWidth,
                                       const gboolean includeFrame,
                                       const char *frameTitle)
{
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), borderWidth);

  if (includeFrame)
    {
      GtkWidget *frame = gtk_frame_new(frameTitle);
      gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);
      gtk_container_add(GTK_CONTAINER(frame), vbox);
    }
  else
    {
      gtk_box_pack_start(GTK_BOX(parent), vbox, FALSE, FALSE, 0);
    }
  
  return vbox;
}

/* Utility to create an hbox with the given border and pack it into the given container */
static GtkWidget* createHBoxWithBorder(GtkWidget *parent, const int borderWidth, const gboolean includeFrame, const char *frameTitle)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(hbox), borderWidth);
  
  if (includeFrame)
    {
      GtkWidget *frame = gtk_frame_new(frameTitle);
      gtk_container_add(GTK_CONTAINER(parent), frame);
      gtk_container_add(GTK_CONTAINER(frame), hbox);
    }
  else
    {
      gtk_container_add(GTK_CONTAINER(parent), hbox);
    }
  
  return hbox;
}


/* Utility to return true if any groups exist. Ignores the 'match set' group
 * if it doesn't have any sequences. */
static gboolean blxWindowGroupsExist(GtkWidget *blxWindow)
{
  gboolean result = FALSE;
  
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  GList *groupList = blxContext->sequenceGroups;
  
  if (g_list_length(groupList) > 1)
    {
      result = TRUE;
    }
  else if (g_list_length(groupList) == 1)
    {
      /* Only one group. If it's the match set group, check it has sequences */
      SequenceGroup *group = (SequenceGroup*)(groupList->data);
      
      if (group != blxContext->matchSetGroup || g_list_length(group->seqList) > 0)
        {
          result = TRUE;
        }
    }
  
  return result;
}


/* Utility to create a text entry widget displaying the given double value. The
 * given callback will be called when the user OK's the dialog that this widget 
 * is a child of. */
static GtkWidget* createTextEntryString(const char *value)
{
  GtkWidget *entry = gtk_entry_new();
  
  gtk_entry_set_text(GTK_ENTRY(entry), value);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(value) + 2);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  return entry;
}


/* Utility to create a text entry widget displaying the given double value. The
 * given callback will be called when the user OK's the dialog that this widget 
 * is a child of. */
static GtkWidget* createTextEntryInt(const int value)
{
  GtkWidget *entry = gtk_entry_new();
  
  char *displayText = convertIntToString(value);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  g_free(displayText);

  return entry;
}


/* This dialog is shown when the user attempts to load a file that
 * is not in a natively-supported format. It asks the user what the
 * source should be, and allows the user to edit the coordinate range
 * to fetch data for. If the user enters valid values and hits OK then
 * the return values are populated and we return TRUE; else return FALSE. */
static gboolean showNonNativeFileDialog(GtkWidget *window, 
                                        const char *filename,
                                        GString **source_out,
                                        int *start_out,
                                        int *end_out)
{
  BlxViewContext *bc = blxWindowGetContext(window);
  
  char *title = g_strdup_printf("%sLoad Non-Native File", blxGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(window),
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_CANCEL,
                                                  GTK_RESPONSE_REJECT,
                                                  GTK_STOCK_OK,
                                                  GTK_RESPONSE_ACCEPT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  GtkContainer *contentArea = GTK_CONTAINER(GTK_DIALOG(dialog)->vbox);

  char *labelStr = g_strdup_printf("\nFile '%s' is not a natively-supported file format.\n\nSpecify the Source to fetch data from this file using an external command\n(a fetch method for the Source must be specified in the config file)\n", filename);  
  GtkWidget *label = gtk_label_new(labelStr);
  g_free(labelStr);
  
  GtkWidget *sourceEntry = createTextEntryString("");
  GtkWidget *label2 = gtk_label_new("\n\nRegion to fetch data for:");
  GtkWidget *startEntry = createTextEntryInt(bc->refSeqRange.min);
  GtkWidget *endEntry = createTextEntryInt(bc->refSeqRange.max);
  
  GtkTable *table = GTK_TABLE(gtk_table_new(5, 2, FALSE));
  gtk_container_add(contentArea, GTK_WIDGET(table));

  gtk_table_attach(table, label, 0, 2, 0, 1, GTK_SHRINK, GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, gtk_label_new("Source"), 0, 1, 1, 2, GTK_SHRINK, GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, sourceEntry, 1, 2, 1, 2, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, label2, 0, 2, 2, 3, GTK_SHRINK, GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, gtk_label_new("Start"), 0, 1, 3, 4, GTK_SHRINK, GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, startEntry, 1, 2, 3, 4, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, gtk_label_new("End"), 0, 1, 4, 5, GTK_SHRINK, GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  gtk_table_attach(table, endEntry, 1, 2, 4, 5, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);

  gtk_widget_show_all(dialog);
  gint response = gtk_dialog_run(GTK_DIALOG(dialog));
  gboolean result = FALSE;
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      const gchar *source = gtk_entry_get_text(GTK_ENTRY(sourceEntry));
      
      /* source is mandatory */
      if (source && *source)
        {
          result = TRUE;
          *source_out = g_string_new(source);

          /* to do: start and end */
        }
    }

  gtk_widget_destroy(dialog);

  return result;
}


/* This function loads the contents of a non-natively supported features-
 * file into blixem, using an external script to convert the file into
 * a supported file format such as GFF. A fetch method stanza must exist in the
 * config to define the script and its parameters.
 * This function asks the user what Source the file relates to so that it can 
 * look up the fetch method that should be used. It optionally also allows the
 * user to specify a coordinate range to limit the fetch to.. */
static void loadNonNativeFile(const char *filename,
                              GtkWidget *blxWindow,
                              MSP **newMsps,
                              GList **newSeqs,
                              GHashTable *lookupTable,
                              const int refSeqOffset,
                              const IntRange* const refSeqRange,
                              GError **error)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GKeyFile *keyFile = blxGetConfig();

  GString *source = NULL;
  int start = bc->refSeqRange.min, end = bc->refSeqRange.max;
  
  if (!showNonNativeFileDialog(blxWindow, filename, &source, &start, &end))
    return;

  GError *tmp_error = NULL;
  BlxDataType *dataType = NULL;
  const BlxFetchMethod *fetchMethod = NULL;

  if (!source || !source->str)
    {
      g_set_error(&tmp_error, BLX_ERROR, 1, "No Source specified; cannot look up fetch method.\n");
    }

  if (!tmp_error)
    {
      dataType = getBlxDataType(0, source->str, keyFile, &tmp_error);

      if (!dataType && !tmp_error)
        g_set_error(&tmp_error, BLX_ERROR, 1, "No data-type found for source '%s'\n", source->str);
    }

  if (!tmp_error)
    {
      if (dataType->bulkFetch)
        {
          GQuark fetchMethodQuark = g_array_index(dataType->bulkFetch, GQuark, 0);
          fetchMethod = getFetchMethodDetails(fetchMethodQuark, bc->fetchMethods);
        }

      if (!fetchMethod)
        {
          g_set_error(&tmp_error, BLX_ERROR, 1, "No fetch method specified for data-type '%s'\n", g_quark_to_string(dataType->name));
        }

      /* The output of the fetch must be a natively supported file format (i.e. GFF) */
      if (!tmp_error && fetchMethod->outputType != BLXFETCH_OUTPUT_GFF)
        {
          g_set_error(&tmp_error, BLX_ERROR, 1, "Expected fetch method output type to be '%s' but got '%s'\n", outputTypeStr(BLXFETCH_OUTPUT_GFF), outputTypeStr(fetchMethod->outputType));
        }
    }
     
  if (!tmp_error)
    {
      MatchSequenceData match_data = {NULL, bc->refSeqName, start, end, bc->dataset, source->str, filename};
      GString *command = doGetFetchCommand(fetchMethod, &match_data, &tmp_error);

      if (!tmp_error && command && command->str)
        {
          const char *fetchName = g_quark_to_string(fetchMethod->name);
          GSList *styles = blxReadStylesFile(NULL, NULL);

          sendFetchOutputToFile(command, keyFile, &bc->blastMode, 
                                bc->featureLists, bc->supportedTypes, styles,
                                &bc->matchSeqs, &bc->mspList, 
                                fetchName, bc->saveTempFiles, newMsps, newSeqs,
                                bc->columnList, lookupTable, refSeqOffset, refSeqRange, &tmp_error);
        }
    }          

  if (tmp_error)
    g_propagate_error(error, tmp_error);
}


/* Dynamically load in additional features from a file. (should be called after
 * blixem's GUI has already started up, rather than during start-up where normal
 * feature-loading happens) */
static void dynamicLoadFeaturesFile(GtkWidget *blxWindow, const char *filename, const char *buffer, GError **error)
{
  /* Must be passed either a filename or buffer */
  if (!filename && !buffer)
    return;

  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GKeyFile *keyFile = blxGetConfig();
  
  /* We'll load the features from the file into some temporary lists */
  MSP *newMsps = NULL;
  GList *newSeqs = NULL;
  GError *tmp_error = NULL;

  /* Create a temporary lookup table for BlxSequences so we can link them on GFF ID */
  GHashTable *lookupTable = g_hash_table_new(g_direct_hash, g_direct_equal);

  /* Assume it's a natively-supported file and attempt to parse it. The first thing this
   * does is check that it's a native file and if not it sets the error */
  loadNativeFile(filename, buffer, keyFile, &bc->blastMode, bc->featureLists, bc->supportedTypes, NULL, &newMsps, &newSeqs, bc->columnList, lookupTable, bc->refSeqOffset, &bc->refSeqRange, &tmp_error);

  if (tmp_error && filename)
    {
      /* Input file is not natively supported. We can still load it if
       * there is a fetch method associated with it: ask the user what
       * the Source is so that we can find the fetch method. Probably 
       * should only get here if the input is an actual file so don't
       * support this for buffers for now. */
      g_error_free(tmp_error);
      tmp_error = NULL;
      
      loadNonNativeFile(filename, blxWindow, &newMsps, &newSeqs, lookupTable, bc->refSeqOffset, &bc->refSeqRange, &tmp_error);
    }

  if (!tmp_error)
    {
      /* Count how many features were added. (Need to do this before blxMergeFeatures because
       * once this list gets merged the count will no longer be correct.) */
      const int numAdded = g_list_length(newSeqs);

      /* Fetch any missing sequence data and finalise the new sequences */
      bulkFetchSequences(0, FALSE, bc->saveTempFiles, bc->seqType, &newSeqs, bc->columnList,
                         bc->bulkFetchDefault, bc->fetchMethods, &newMsps, &bc->blastMode,
                         bc->featureLists, bc->supportedTypes, NULL, bc->refSeqOffset,
                         &bc->refSeqRange, bc->dataset, FALSE, lookupTable);

      finaliseFetch(newSeqs, bc->columnList);

      finaliseBlxSequences(bc->featureLists, &newMsps, &newSeqs, bc->columnList, bc->refSeqOffset, bc->seqType, 
                           bc->numFrames, &bc->refSeqRange, TRUE, lookupTable);

      double lowestId = calculateMspData(newMsps, bc);
      bigPictureSetMinPercentId(blxWindowGetBigPicture(blxWindow), lowestId);

      /* Add the msps/sequences to the tree data models (must be done after finalise because
       * finalise populates the child msp lists for parent features) */
      detailViewAddMspData(blxWindowGetDetailView(blxWindow), newMsps, newSeqs);

      /* Merge the temporary lists into the main lists (takes ownership of the temp lists) */
      blxMergeFeatures(newMsps, newSeqs, &bc->mspList, &bc->matchSeqs);

      /* Cache the new msp display ranges and sort and filter the trees. */
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      cacheMspDisplayRanges(bc, detailViewGetNumUnalignedBases(detailView));
      detailViewResortTrees(detailView);
      callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);

      /* Recalculate the coverage */
      calculateDepth(bc);
      updateCoverageDepth(blxWindowGetCoverageView(blxWindow), bc);
  
      /* Re-calculate the height of the exon views */
      GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
      calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
      calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
      forceResize(bigPicture);
  
      blxWindowRedrawAll(blxWindow);
      
      if (numAdded == 0)
        g_warning("No features loaded\n");
      else if (numAdded == 1)
        g_message("Loaded %d new feature\n", numAdded);
      else
        g_message("Loaded %d new features\n", numAdded);
    }

  g_hash_table_unref(lookupTable);

  if (tmp_error)
    g_propagate_error(error, tmp_error);
}


/***********************************************************
 *                         View panes menu                 *
 ***********************************************************/

/* Toggle visibility the n'th tree. This is the active strand's frame n if displaying
 * protein matches (where we only display one strand), or the forward or reverse
 * strand tree if displaying DNA matches (where both strands are displayed). */
static void toggleTreeVisibility(GtkWidget *blxWindow, const int number)
{
  const gboolean toggled = blxWindowGetDisplayRev(blxWindow);
  const BlxStrand activeStrand = toggled ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
  
  /* For protein matches, trees are always displayed in frame order (i.e. 1, 2, 3), 
   * so just use the number pressed for the frame, and the active strand for the
   * strand. */
  int frame = number;
  BlxStrand strand = activeStrand;
  
  /* For DNA matches, the frame is always 1, but the strand depends on which number
   * was pressed: use 1 to toggle active strand, 2 for other strand */
  if (blxWindowGetSeqType(blxWindow) == BLXSEQ_DNA)
    {
      frame = 1;
      
      if (number == 1)
        {
          strand = activeStrand;
        }
      else if (number == 2)
        {
          strand = toggled ? BLXSTRAND_FORWARD : BLXSTRAND_REVERSE;
        }
    }
  
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (tree && gtk_widget_get_parent(tree))
    {
      widgetSetHidden(tree, !widgetGetHidden(tree));
    }
}


/* Toggle visibility of the active (1) or other (2) strand grid depending on the number pressed */
static void toggleGridVisibility(GtkWidget *blxWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
      const gboolean useFwdGrid = (number == 1) != blxWindowGetDisplayRev(blxWindow);
      
      GtkWidget *grid = useFwdGrid ? bigPictureGetFwdGrid(bigPicture) : bigPictureGetRevGrid(bigPicture);
      widgetSetHidden(grid, !widgetGetHidden(grid));

      /* We need to force a resize of the big picture because the size-allocate
       * signal doesn't get emitted by default when its contents shrink */
      forceResize(bigPicture);
    }
}


/* Toggle visibility of the active (1) or other (2) strand exon view depending on the number pressed */
static void toggleExonViewVisibility(GtkWidget *blxWindow, const int number)
{
  if (number == 1 || number == 2)
    {
      GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
      const gboolean useFwdExonView = (number == 1) != blxWindowGetDisplayRev(blxWindow);
      
      GtkWidget *exonView = useFwdExonView ? bigPictureGetFwdExonView(bigPicture) : bigPictureGetRevExonView(bigPicture);
      widgetSetHidden(exonView, !widgetGetHidden(exonView));

      forceResize(bigPicture);
    }
}


/* Toggle the visibility of tree/grid panes following a number key press */
static void togglePaneVisibility(GtkWidget *blxWindow, const int number, const gboolean modifier1, const gboolean modifier2)
{
  /* Affects big picture if modifier1 was pressed, the detail view otherwise */
  if (modifier1)
    {
      /* If modifier 2 was also pressed, affects the exon views; otherwise the grids */
      if (modifier2)
        {
          toggleExonViewVisibility(blxWindow, number);
        }
      else
        {
          toggleGridVisibility(blxWindow, number);
        }
    }
  else
    {
      toggleTreeVisibility(blxWindow, number);
    }
}


/* Repeat the last find operation. Searches for the next (rightwards) match unless the given
 * modifier is pressed, in which case it searches for the previous (leftwards) match */
static void findAgain(GtkWidget *blxWindow, const gboolean modifier)
{
  GError *error = NULL;

  const int startCoord = getSearchStartCoord(blxWindow, FALSE, modifier);

  /* Try the DNA search. Does nothing if last search was not a DNA search. */
  blxWindowFindDnaString(blxWindow, NULL, startCoord, modifier, TRUE, &error);

  if (error)
    {
      /* DNA search was attempted but not found. Try looping round to the beginning */
      g_error_free(error);
      error = NULL;
      const int newStart = getSearchStartCoord(blxWindow, TRUE, modifier);
      blxWindowFindDnaString(blxWindow, NULL, newStart, modifier, TRUE, &error);
    }
  
  if (!error)
    {
      /* Try the search-from-list search. Returns NULL if last search was not a list search */
      GList *seqList = findSeqsFromList(blxWindow, NULL, BLXCOL_NONE, FALSE, TRUE, &error);
      
      if (!seqList && !error)
        {
          /* Try the search-by-name search. Returns NULL if last search was not a name search. */
          seqList = findSeqsFromColumn(blxWindow, NULL, BLXCOL_NONE, FALSE, TRUE, &error);
        }
      
      /* If either the list or name search succeeded, select the prev/next MSP from the 
       * found sequence(s) depending on which direction we're searching. */
      if (seqList)
        {
          blxWindowSetSelectedSeqList(blxWindow, seqList);
          
          if (modifier)
            {
              prevMatch(blxWindowGetDetailView(blxWindow), seqList);
            }
          else
            {
              nextMatch(blxWindowGetDetailView(blxWindow), seqList);
            }
        }
    }
  
  if (error)
    {
      prefixError(error, "Find %s failed. ", (modifier ? "previous" : "next"));
      reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
    }
}


/* Called when the state of a check button is toggled */
static void onVisibilityButtonToggled(GtkWidget *button, gpointer data)
{
  GtkWidget *widgetToToggle = GTK_WIDGET(data);
  gboolean visible = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  widgetSetHidden(widgetToToggle, !visible);
}


/* Create a check button to control visibility of the given widget */
static void createVisibilityButton(GtkWidget *widgetToToggle, const char *mnemonic, GtkWidget *container)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_container_add(GTK_CONTAINER(container), button);

  /* Set the state depending on the widget's current visibility */
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), GTK_WIDGET_VISIBLE(widgetToToggle));
  
  g_signal_connect(G_OBJECT(button), "toggled", G_CALLBACK(onVisibilityButtonToggled), widgetToToggle);
}


/* Create a check button to control visibility of the given tree. */
static void createTreeVisibilityButton(GtkWidget *detailView, const BlxStrand strand, const int frame, GtkWidget *container)
{
  /* Some trees may have been removed from the blixem window if they are not on the active 
   * strand, so only show check boxes for those that are in the window (i.e. have a parent). 
   * Note that we get the tree container here, which consists of the tree itself plus any headers etc. */
  GtkWidget *tree = detailViewGetTreeContainer(detailView, strand, frame);
  
  if (gtk_widget_get_parent(tree))
    {
      const gboolean toggled = detailViewGetDisplayRev(detailView);
      gboolean isActiveStrand = ((strand == BLXSTRAND_FORWARD) != toggled);

      if (detailViewGetSeqType(detailView) == BLXSEQ_DNA)
        {
          /* We only have 1 frame, but trees are from both strands, so distinguish between strands.
           * Put each strand in its own frame. */
          char text1[] = "Show _active strand";
          char text2[] = "Show othe_r strand";

          GtkWidget *frame = gtk_frame_new(isActiveStrand ? "Active strand" : "Other strand");
          gtk_container_add(GTK_CONTAINER(container), frame);
          createVisibilityButton(tree, isActiveStrand ? text1 : text2, frame);
        }
      else
        {
          /* All the visible trees should be in the same strand, so just distinguish by frame number. */
          char formatStr[] = "Show frame _%d";
          char displayText[strlen(formatStr) + numDigitsInInt(frame) + 1];
          sprintf(displayText, formatStr, frame);

          createVisibilityButton(tree, displayText, container);
        }
    }
}


/* Callback called when the user clicks the 'bump exon view' button */
static void onBumpExonView(GtkWidget *button, gpointer data)
{
  GtkWidget *exonView = GTK_WIDGET(data);
  const gboolean expanded = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  exonViewSetExpanded(exonView, expanded);
}


/* Create the set of settings buttons to control display of an exon-view widget. */
static void createExonButtons(GtkWidget *exonView, const char *visLabel, const char *bumpLabel, GtkWidget *parent)
{
  /* Pack everything in an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(parent), hbox);
  
  /* Create a check button to control visibility of the exon view */
  createVisibilityButton(exonView, visLabel, hbox);

  /* Create a check button to control whether the exon view is expanded or compressed */
  const gboolean isBumped = exonViewGetExpanded(exonView);
  createCheckButton(GTK_BOX(hbox), bumpLabel, isBumped, G_CALLBACK(onBumpExonView), exonView);
}


/* Shows the "View panes" dialog. This dialog allows the user to show/hide certain portions of the window. */
void showViewPanesDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_VIEW;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sView panes", blxGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title, 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);

      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Clear contents and re-create */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  int borderWidth = 12;
  
  /* Big picture */
  GtkWidget *bp = blxWindowGetBigPicture(blxWindow);
  GtkWidget *bpVbox = createVBoxWithBorder(contentArea, borderWidth, TRUE, "Big picture");
  
  createVisibilityButton(bp, "Show _big picture", bpVbox);
  GtkWidget *bpSubBox = createVBoxWithBorder(bpVbox, borderWidth, FALSE, NULL);
  
  GtkWidget *bpActiveStrand = createVBoxWithBorder(bpSubBox, 0, TRUE, "Active strand");
  createVisibilityButton(bigPictureGetActiveGrid(bp), "Show _grid", bpActiveStrand);
  createExonButtons(bigPictureGetActiveExonView(bp), "Show _exons    ", "_Bump exons    ", bpActiveStrand);
  
  GtkWidget *bpOtherStrand = createVBoxWithBorder(bpSubBox, 0, TRUE, "Other strand");
  createVisibilityButton(bigPictureGetInactiveGrid(bp), "Show gr_id", bpOtherStrand);
  createExonButtons(bigPictureGetInactiveExonView(bp), "Show e_xons    ", "Bum_p exons    ", bpOtherStrand);
  
  /* Detail view */
  GtkWidget *dvVbox = createVBoxWithBorder(contentArea, borderWidth, TRUE, "Alignment lists");
  GtkWidget *dv = blxWindowGetDetailView(blxWindow);
  createVisibilityButton(dv, "Show alignment _lists", dvVbox);
  
  GtkWidget *dvSubBox = createVBoxWithBorder(dvVbox, borderWidth, FALSE, NULL);
  const int numFrames = blxWindowGetNumFrames(blxWindow);
  int frame = 1;
  for ( ; frame <= numFrames; ++frame)
    {
      createTreeVisibilityButton(dv, blxWindowGetActiveStrand(blxWindow), frame, dvSubBox);
      createTreeVisibilityButton(dv, blxWindowGetInactiveStrand(blxWindow), frame, dvSubBox);
    }
  
  /* Coverage view */
  GtkWidget *coverageView = bigPictureGetCoverageView(bp);
  GtkWidget *coverageVbox = createVBoxWithBorder(contentArea, borderWidth, TRUE, "Coverage view");
  createVisibilityButton(coverageView, "Show _coverage view", coverageVbox);

  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}



/***********************************************************
 *                          Find menu                      *
 ***********************************************************/

static GList* findSeqsFromColumn(GtkWidget *blxWindow, const char *inputText, const BlxColumnId inputCol, const gboolean rememberSearch, const gboolean findAgain, GError **error)
{
  /* Previous values (if applicable) */
  static char *prevSearchStr = NULL;
  static BlxColumnId prevSearchCol = BLXCOL_NONE;

  /* Current values */
  char *searchStr = NULL;
  BlxColumnId searchCol = BLXCOL_NONE;
  
  /* If it's a find-again, use the existing values; otherwise, use the input values */
  if (findAgain)
    {
      searchStr = prevSearchStr;
      searchCol = prevSearchCol;
    }
  else
    {
      g_free(searchStr);
      searchStr = g_strdup(inputText);
      searchCol = inputCol;
    
      if (rememberSearch)
        {
          prevSearchStr = searchStr;
          prevSearchCol = searchCol;
        }
    }
  
  if (!searchStr || searchCol == BLXCOL_NONE)
    {
      /* We will get here if we do a find-again when there wasn't a previous find */
      return NULL;
    }
  
  /* Loop through all the sequences and see if the sequence data for this column
   * matches the search string */
  GList *seqList = blxWindowGetAllMatchSeqs(blxWindow);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  SeqSearchData searchData = {searchStr, searchCol, bc, NULL, NULL};
  
  g_list_foreach(seqList, getSequencesThatMatch, &searchData);
  
  if (g_list_length(searchData.matchList) < 1)
    {
      GList *columnList = blxWindowGetColumnList(blxWindow);
      const char *columnName = getColumnTitle(columnList, searchCol);

      if (searchData.error)
        g_propagate_error(error, searchData.error);
      else
        g_set_error(error, BLX_ERROR, BLX_ERROR_STRING_NOT_FOUND, "No sequences found where column '%s' matches text '%s'.\n", columnName, searchStr);
    }
  
  return searchData.matchList;
}


/* Utility to extract the contents of a GtkTextView and return it as a string. The result is
 * owned by the GtkTextView and should not be free'd. */
static const char* getStringFromTextView(GtkTextView *textView)
{
  if (!textView || !GTK_WIDGET_SENSITIVE(GTK_WIDGET(textView)))
    {
      g_critical("Could not set search string: invalid text entry box\n");
      return NULL;
    }
  
  /* Get the input text from the text buffer and create the group */
  GtkTextBuffer *textBuffer = gtk_text_view_get_buffer(textView);
  
  GtkTextIter start, end;
  gtk_text_buffer_get_bounds(textBuffer, &start, &end);
  
  return gtk_text_buffer_get_text(textBuffer, &start, &end, TRUE);
}


/* Finds all the valid sequences blixem knows about whose column data matches
 * an item from the given input text. Returns the results as a GList of BlxSequences.
 * The input text may be a multi-line list (one search item per line). */
static GList* findSeqsFromList(GtkWidget *blxWindow,
                               const char *inputText, 
                               const BlxColumnId inputCol, 
                               const gboolean rememberSearch, 
                               const gboolean findAgain, 
                               GError **error)
{
  /* Previous values (if applicable) */
  static char *prevSearchStr = NULL;
  static BlxColumnId prevSearchCol = BLXCOL_NONE;
  
  /* Current values */
  char *searchStr = NULL;
  BlxColumnId searchCol = BLXCOL_NONE;
  
  /* If we-re doing a find-again, use the values from last time; otherwise use the input values */
  if (findAgain)
    {
      searchStr = prevSearchStr;
      searchCol = prevSearchCol;
    }
  else
    {
      g_free(searchStr);
      searchStr = g_strdup(inputText);
      searchCol = inputCol;
    
      if (rememberSearch)
        {
          prevSearchStr = searchStr;
          prevSearchCol = searchCol;
        }
    }
      
  if (!searchStr || searchCol == BLXCOL_NONE)
    {
      /* We may get here if we did a find-again when there was no previous find */
      return NULL;
    }

  
  GError *tmpError = NULL;
  GList *seqList = getSeqStructsFromText(blxWindow, searchStr, searchCol, &tmpError);

  if (g_list_length(seqList) < 1)
    {
      GList *columnList = blxWindowGetColumnList(blxWindow);
      const char *columnName = getColumnTitle(columnList, searchCol);
      
      if (tmpError)
        g_propagate_error(error, tmpError);
      else
        g_set_error(error, BLX_ERROR, BLX_ERROR_STRING_NOT_FOUND, "No sequences found where column '%s' matches text '%s'.\n", columnName, searchStr);
    }
  
  return seqList;
}


/* Callback called when requested to find sequences from a sequence name. Selects
 * the sequences and scrolls to the start of the first match in the selection */
static gboolean onFindSeqsFromName(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  const char *inputText = NULL;
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      inputText = getStringFromTextEntry(GTK_ENTRY(data));
    }
  
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  
  /* Find the combo box on the dialog (there should only be one), which tells
   * us which column to search. */
  GtkComboBox *combo = widgetGetComboBox(dialog);
  BlxColumnId searchCol = getColumnFromComboBox(combo);

  /* Find all sequences that match */
  GError *error = NULL;
  GList *seqList = findSeqsFromColumn(blxWindow, inputText, searchCol, TRUE, FALSE, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      result = FALSE;
    }
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      blxWindowSetSelectedSeqList(blxWindow, seqList);
      
      if (responseId == BLX_RESPONSE_FORWARD)
        {
          nextMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else if (responseId == BLX_RESPONSE_BACK)
        {
          prevMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else
        {
          firstMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
    }
    
  return result;
}


/* Callback called when requested to find sequences from a given list. Selects
 * the sequences ands scrolls to the start of the first match in the selection. */
static gboolean onFindSeqsFromList(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  const char *inputText = NULL;
  
  /* Nothing to do if this button is not active */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      inputText = getStringFromTextView(GTK_TEXT_VIEW(data));
    }
  
  /* Get the dialog and main window */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));

  /* Find the combo box on the dialog (there should only be one), which tells
   * us which column to search. */
  GtkComboBox *combo = widgetGetComboBox(dialog);
  BlxColumnId searchCol = getColumnFromComboBox(combo);
  
  GError *error = NULL;
  GList *seqList = findSeqsFromList(blxWindow, inputText, searchCol, TRUE, FALSE, &error);

  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      result = FALSE;
    }
  
  if (seqList)
    {
      blxWindowSetSelectedSeqList(blxWindow, seqList);
      
      if (responseId == BLX_RESPONSE_FORWARD)
        {
          nextMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else if (responseId == BLX_RESPONSE_BACK)
        {
          prevMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
      else
        {
          firstMatch(blxWindowGetDetailView(blxWindow), seqList);
        }
    }
  
  return result;
}


/* Search for the given DNA string in the reference sequence. Searches for the next (rightwards)
 * value from the given start coord, unless searchLeft is true in which case it searches leftwards.
 * If findAgain is true it repeats the last DNA search. */
static void blxWindowFindDnaString(GtkWidget *blxWindow, 
                                   const char *inputSearchStr, 
                                   const int refSeqStart,
                                   const gboolean searchLeft, 
                                   const gboolean findAgain, 
                                   GError **error)
{
  /* Remember the last input string for use with findAgain */
  static char *searchStr = NULL;
  
  if (!findAgain)
    {
      /* We must copy the input string because it may not exist if/when we come to do a 'find again' */
      g_free(searchStr);
      searchStr = g_strdup(inputSearchStr);
    }
  
  const int searchStrMax = searchStr ? strlen(searchStr) - 1 : -1;
  
  if (searchStrMax < 0)
    {
      return;
    }

  const int searchStart = searchLeft ? searchStrMax : 0;
  const int searchEnd = searchLeft ? 0 : searchStrMax;
  const int searchStrIncrement = searchLeft ? -1 : 1;
  
  /* Values increase left-to-right in normal display or right-to-left in reversed display */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const gboolean searchForward = (searchLeft == bc->displayRev);
  const int refSeqIncrement = searchForward ? 1 : -1;
  
  /* We'll need to complement ref seq bases if the active strand is the reverse strand */
  const gboolean complement = (blxWindowGetActiveStrand(blxWindow) == BLXSTRAND_REVERSE);
  
  int refSeqIdx = refSeqStart;
  int searchStrIdx = searchStart;
  int matchStart = UNSET_INT;
  
  while (refSeqIdx >= bc->refSeqRange.min && refSeqIdx <= bc->refSeqRange.max && searchStrIdx >= 0 && searchStrIdx <= searchStrMax)
    {
      const char refSeqBase = getSequenceIndex(bc->refSeq, refSeqIdx, complement, &bc->refSeqRange, BLXSEQ_DNA);
      char searchStrBase = convertBaseToCorrectCase(searchStr[searchStrIdx], BLXSEQ_DNA);      
      
      if (refSeqBase == searchStrBase)
        {
          /* The base matches. If it's the first matching base, set the match-start coord (or if we're 
           * searching leftwards, then always set the match-start coord, because the start is actually 
           * the last coord that will be found). Then proceed to the next position in the search string */
          if (matchStart == UNSET_INT)
            {
              matchStart = refSeqIdx;
            }
          
          searchStrIdx += searchStrIncrement;
          refSeqIdx += refSeqIncrement;
        }
      else if (matchStart != UNSET_INT)
        {
          /* We were in a match but this base doesn't match. Reset to the start of the 
           * search string, and start looking again from one base after the place where the last
           * match started. (We need to re-check all bases from there because we're comparing
           * against a different section of the search string now.) */
          searchStrIdx = searchStart;
          refSeqIdx = matchStart + refSeqIncrement;
          matchStart = UNSET_INT;
        }
      else
        {
          refSeqIdx += refSeqIncrement;
        }
    }
  
  /* Undo the last increment, so that we have the final coords of the matching section (if found) */
  refSeqIdx -= refSeqIncrement;
  searchStrIdx -= searchStrIncrement;
  
  /* If we reached the end of the search string, then we matched the whole lot. */
  const gboolean finished = searchStrIdx == searchEnd;
  
  if (matchStart != UNSET_INT && finished)
    {
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      const int frame = 1;
      int baseNum = UNSET_INT;
      
      int result = searchLeft ? refSeqIdx : matchStart;
      result = convertDnaIdxToDisplayIdx(result, bc->seqType, frame, bc->numFrames, bc->displayRev, &bc->refSeqRange, &baseNum);
      
      detailViewSetSelectedBaseIdx(detailView, result, frame, baseNum, TRUE, FALSE);
      detailViewRedrawAll(detailView);
    }
  else
    {
      g_set_error(error, BLX_ERROR, BLX_ERROR_STRING_NOT_FOUND, "The string '%s' was not found in the reference sequence searching to the %s from coord %d.\n", searchStr, (searchLeft ? "left" : "right"), refSeqStart);
    }
}


/* Get the start coord for a search. If startBeginning is false, this gets the currently-selected display
 * index (shifted by one base so that we don't start searching at the same position as a previous
 * find result) or, if no base index is selected, returns the start coord of the current display range. If
 * startBeginning is true, just start from the beginning of the reference sequence. The result is 
 * nucleotide coord on the ref sequence. */
static int getSearchStartCoord(GtkWidget *blxWindow, const gboolean startBeginning, const gboolean searchLeft)
{
  int result = UNSET_INT;
  
  const BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  if (startBeginning)
    {
      result = (searchLeft == bc->displayRev) ? bc->refSeqRange.min : bc->refSeqRange.max;
    }
  else  
    {
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      result = detailViewGetSelectedBaseIdx(detailView);
      
      if (result != UNSET_INT)
        {
          /* Increment by one to make sure we don't re-find a previously-found match
           * (or decrement if searching leftwards) */
           if (searchLeft)
             {
                --result;
              }
            else
              {
                ++result;
              }
        }
      else
        {
          /* The start display coord is the min coord if we're searching left and the max if searching right. */
          const IntRange* const displayRange = detailViewGetDisplayRange(detailView);
          result = searchLeft ? displayRange->max : displayRange->min;
        }

      /* Convert the display coord to a nucleotide coord */
      result = convertDisplayIdxToDnaIdx(result, bc->seqType, 1, 1, bc->numFrames, bc->displayRev, &bc->refSeqRange);
    }
  
  return result;
}


/* Callback called when requested to search for a DNA string. If found, sets the currently-
 * selected base index to the coord where the matching string starts. The text entry for the
 * search string is passed as the callback data. */
static gboolean onFindDnaString(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  /* Get the search string from the text entry. If the toggle button is not active, call
   * blxWindowFindDnaString with a NULL search string to "cancel" any previous searches
   * so that "findAgain" will not attempt to perform a DNA search). */
  const char *searchStr = NULL;

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      searchStr = getStringFromTextEntry(GTK_ENTRY(data));

      if (!searchStr || strlen(searchStr) < 1)
        {
          g_critical("DNA search failed. The search string was empty.\n");
          result = FALSE;
        }
    }

  /* Search left wrt the screen if the user hit 'back' search right if 'forward' */
  const gboolean searchLeft = (responseId == BLX_RESPONSE_BACK);

  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  const gboolean startBeginning = (responseId != BLX_RESPONSE_FORWARD && responseId != BLX_RESPONSE_BACK);
  int startCoord = getSearchStartCoord(blxWindow, startBeginning, searchLeft);
  
  GError *error = NULL;
  blxWindowFindDnaString(blxWindow, searchStr, startCoord, searchLeft, FALSE, &error);
  
  if (error)
    {
      if (!startBeginning)
        {
          /* Try looping round to the beginning */
          postfixError(error, " Trying again from the %s of the range.\n", (searchLeft ? "end" : "start"));
          reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
          
          startCoord = getSearchStartCoord(blxWindow, TRUE, searchLeft);
          blxWindowFindDnaString(blxWindow, searchStr, startCoord, searchLeft, FALSE, &error);
        }
    }
  
  if (error)
    {
      result = FALSE;
      prefixError(error, "DNA search failed. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  return result;
}


/* Utility to create a drop-down combo box for selecting the column to search by */
static void createSearchColumnCombo(GtkTable *table, const int col, const int row, GtkWidget *blxWindow)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_table_attach(table, hbox, col, col + 1, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, DEFAULT_TABLE_XPAD, DEFAULT_TABLE_YPAD);
  
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);

  createSortBox(GTK_BOX(hbox), detailView, BLXCOL_SEQNAME, columnList, "Search column: ", TRUE);  
}


static void onClearFindDialog(GtkWidget *button, gpointer data)
{
  GSList *entryList = (GSList*)data;
  
  for ( ; entryList; entryList = entryList->next)
    {
      if (GTK_IS_ENTRY(entryList->data))
        {
          GtkEntry *entry = GTK_ENTRY(entryList->data);
          gtk_entry_set_text(entry, "");
        }
      else if (GTK_IS_TEXT_VIEW(entryList->data))
        {
          GtkTextView *textView = GTK_TEXT_VIEW(entryList->data);
          gtk_text_buffer_set_text(gtk_text_view_get_buffer(textView), "", -1);
        }
      else
        {
          g_warning("onClearFindDialog: Unexpected widget type: expected GtkEntry or GtkTextView\n");
        }
    }
}


/* Clear up data created for the find dialog when it is destroyed
 * (here for completeness but not actually called because the dialog
 * is persistent). */
static void onDestroyFindDialog(GtkWidget *widget, gpointer data)
{
  GSList *entryList = (GSList*)data;

  if (entryList)
    {
      g_slist_free(entryList);
    }
}


/* Show the 'Find' dialog */
void showFindDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_FIND;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sFind sequences", blxGetTitlePrefix(bc));

      /* Note that we add some buttons here but some more at the end because
       * we want to create a custom Clear button in the middle somewhere */
      dialog = gtk_dialog_new_with_buttons(title, 
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_GO_BACK,     /* previous match */
                                           BLX_RESPONSE_BACK,
                                           GTK_STOCK_GO_FORWARD,  /* next match */
                                           BLX_RESPONSE_FORWARD,
                                           NULL);
      
      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      GtkBox *contentArea = GTK_BOX(GTK_DIALOG(dialog)->vbox);
      GtkBox *actionArea = GTK_BOX(GTK_DIALOG(dialog)->action_area);
      const int numRows = 3;
      const int numCols = 2;
      GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
      gtk_box_pack_start(contentArea, GTK_WIDGET(table), TRUE, TRUE, 0);

      /* This list will be populated with the text entry widgets. */
      GSList *entryList = NULL;

      /* Column 1: match-seq search options */
      GtkRadioButton *button1 = createRadioButton(table, 1, 1, NULL, "_Text search (wildcards * and ?)", TRUE, TRUE, FALSE, onFindSeqsFromName, blxWindow, &entryList);
      createRadioButton(table, 1, 2, button1, "_List search", FALSE, TRUE, TRUE, onFindSeqsFromList, blxWindow, &entryList);
      createSearchColumnCombo(table, 1, 3, blxWindow);
      
      /* Column 2: ref-seq search options */
      createRadioButton(table, 2, 1, button1, "_DNA search", FALSE, TRUE, FALSE, onFindDnaString, blxWindow, &entryList);
      
      /* Add a button to clear the text entry fields. It's easier to do this
       * here than in the response callback because we want to send different data */
      GtkWidget *clearButton = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
      g_signal_connect(G_OBJECT(clearButton), "clicked", G_CALLBACK(onClearFindDialog), entryList);
      gtk_box_pack_end(actionArea, clearButton, FALSE, FALSE, 0);

      /* Add remaining buttons after the Clear button */
      gtk_dialog_add_buttons(GTK_DIALOG(dialog),
                             GTK_STOCK_CLOSE,       /* close / cancel */
                             GTK_RESPONSE_REJECT,
                             GTK_STOCK_OK,          /* ok, do the search */
                             GTK_RESPONSE_ACCEPT,
                             NULL);

      gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(blxWindow));
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
      g_signal_connect(dialog, "destroy", G_CALLBACK(onDestroyFindDialog), entryList);
    }

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/* Show the 'Info' dialog, which displays info about the currently-selected sequence(s) */
void showInfoDialog(GtkWidget *blxWindow)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons("Blixem - Sequence info", 
						  NULL, 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CLOSE,
						  GTK_RESPONSE_REJECT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_REJECT);

  int width = blxWindow->allocation.width * 0.7;
  int height = blxWindow->allocation.height * 0.9;
  
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Compile the message text from the selected sequence(s) */
  GString *resultStr = g_string_new("");
  GList *seqItem = bc->selectedSeqs;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      char *seqText = blxSequenceGetInfo(blxSeq, TRUE, bc->columnList);
      g_string_append_printf(resultStr, "%s\n\n", seqText);
      g_free(seqText);
    }
  
  /* We'll use the same fixed-width font as the detail-view */
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  PangoFontDescription *fontDesc = detailViewGetFontDesc(detailView);
  
  GtkWidget *child = createScrollableTextView(resultStr->str, TRUE, fontDesc, TRUE, NULL, &height, NULL);
                             
  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  gtk_widget_show_all(dialog);
  
  g_string_free(resultStr, TRUE);
}


/* Toggle the bump state. Currently only the exon view can be bumped, so this just
 * affects that. */
static void toggleBumpState(GtkWidget *blxWindow)
{
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  exonViewToggleExpanded(bigPictureGetFwdExonView(bigPicture));
  exonViewToggleExpanded(bigPictureGetRevExonView(bigPicture));
}

/***********************************************************
 *                    Group sequences menu                 *
 ***********************************************************/

/* Utility to free the given list of  strings and (if the option is true) all
 * of its data items as well. */
static void freeStringList(GList **stringList, const gboolean freeDataItems)
{
  if (stringList && *stringList)
    {
      if (freeDataItems)
        {
          GList *item = *stringList;
          for ( ; item; item = item->next)
            {
              char *strData = (char*)(item->data);
              g_free(strData);
              item->data = NULL;
            }
        }
      
      g_list_free(*stringList);
      *stringList = NULL;
    }
}


/* Free the memory used by the given sequence group and its members. */
static void destroySequenceGroup(BlxViewContext *bc, SequenceGroup **seqGroup)
{
  if (seqGroup && *seqGroup)
    {
      /* Remove it from the list of groups */
      bc->sequenceGroups = g_list_remove(bc->sequenceGroups, *seqGroup);
      
      /* If this is pointed to by the match-set pointer, null it */
      if (*seqGroup == bc->matchSetGroup)
        {
          bc->matchSetGroup = NULL;
        }
      
      /* Free the memory used by the group name */
      if ((*seqGroup)->groupName)
        {
          g_free((*seqGroup)->groupName);
        }
      
      /* Free the list of sequences */
      if ((*seqGroup)->seqList)
        {
          freeStringList(&(*seqGroup)->seqList, (*seqGroup)->ownsSeqNames);
        }
      
      g_free(*seqGroup);
      *seqGroup = NULL;
    }
}


/* Delete a single group */
static void blxWindowDeleteSequenceGroup(GtkWidget *blxWindow, SequenceGroup *group)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext->sequenceGroups)
    {
      destroySequenceGroup(blxContext, &group);
      blxWindowGroupsChanged(blxWindow);
    }
}


/* Delete all groups */
static void blxContextDeleteAllSequenceGroups(BlxViewContext *bc)
{
  GList *groupItem = bc->sequenceGroups;
  
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      destroySequenceGroup(bc, &group);
    }
  
  g_list_free(bc->sequenceGroups);
  bc->sequenceGroups = NULL;
  
  /* Reset the hide-not-in-group flags otherwise we'll hide everything! */
  bc->flags[BLXFLAG_HIDE_UNGROUPED_SEQS] = FALSE ;
  bc->flags[BLXFLAG_HIDE_UNGROUPED_FEATURES] = FALSE ;
}


static void blxWindowDeleteAllSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  blxContextDeleteAllSequenceGroups(bc);
  blxWindowGroupsChanged(blxWindow);
}


/* Update function to be called whenever groups have been added or deleted,
 * or sequences have been added to or removed from a group */
static void blxWindowGroupsChanged(GtkWidget *blxWindow)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  /* Re-sort all trees, because grouping affects sort order */
  detailViewResortTrees(detailView);
  
  /* Refilter the trees (because groups affect whether sequences are visible) */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);

  /* Resize exon view because transcripts may have been hidden/unhidden */
  calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
  calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
  forceResize(bigPicture);

  /* Redraw all (because highlighting affects both big picture and detail view) */
  blxWindowRedrawAll(blxWindow);
}


/* Create a new sequence group from the given list of sequence names, with a
 * unique ID and name, and add it to the blxWindow's list of groups. The group 
 * should be destroyed with destroySequenceGroup. If ownSeqNames is true, the group
 * will take ownership of the sequence names and free them when it is destroyed. 
 * Caller can optionally provide the group name; if not provided, a default name
 * will be allocated. */
static SequenceGroup* createSequenceGroup(GtkWidget *blxWindow, GList *seqList, const gboolean ownSeqNames, const char *groupName)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* Create the new group */
  SequenceGroup *group = (SequenceGroup*)g_malloc(sizeof(SequenceGroup));
  
  group->seqList = seqList;
  group->ownsSeqNames = ownSeqNames;
  group->hidden = FALSE;
  
  /* Find a unique ID */
  GList *lastItem = g_list_last(bc->sequenceGroups);
  
  if (lastItem)
    {
      SequenceGroup *lastGroup = (SequenceGroup*)(lastItem->data);
      group->groupId = lastGroup->groupId + 1;
    }
  else
    {
      group->groupId = 1;
    }

  if (groupName)
    {
      group->groupName = g_strdup(groupName);
    }
  else
    {
      /* Create a default name based on the unique ID */
      char formatStr[] = "Group%d";
      const int nameLen = strlen(formatStr) + numDigitsInInt(group->groupId);
      group->groupName = (char*)g_malloc(nameLen * sizeof(*group->groupName));
      sprintf(group->groupName, formatStr, group->groupId);
    }
  
  /* Set the order number. For simplicity, set the default order to be the same
   * as the ID number, so groups are sorted in the order they were added */
  group->order = group->groupId;

  /* Set the default highlight color. */
  group->highlighted = TRUE;

  BlxColorId colorId = (groupName && !strcmp(groupName, MATCH_SET_GROUP_NAME)) ? BLXCOLOR_MATCH_SET : BLXCOLOR_GROUP;
  GdkColor *color = getGdkColor(colorId, bc->defaultColors, FALSE, bc->usePrintColors);
  group->highlightColor = *color;

  /* Add it to the list, and update */
  bc->sequenceGroups = g_list_append(bc->sequenceGroups, group);
  blxWindowGroupsChanged(blxWindow);
  
  return group;
}


/* This function sets the sequence-group-name text based on the given text entry widget */
static gboolean onGroupNameChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = TRUE;

  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newName = gtk_entry_get_text(entry);
  
  if (!newName || strlen(newName) < 1)
    {
      g_critical("Invalid group name '%s' entered; reverting to previous group name '%s'.", newName, group->groupName);
      gtk_entry_set_text(entry, group->groupName);
      result = FALSE;
    }
  else
    {
      if (group->groupName) 
        g_free(group->groupName);
      
      group->groupName = g_strdup(newName);
    }
  
  return result;
}


/* This function is called when the sequence-group-order text entry widget's
 * value has changed. It sets the new order number in the group. */
static gboolean onGroupOrderChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkEntry *entry = GTK_ENTRY(widget);
  SequenceGroup *group = (SequenceGroup*)data;
  
  const gchar *newOrder = gtk_entry_get_text(entry);
  
  if (!newOrder || strlen(newOrder) < 1)
    {
      g_critical("Invalid order number '%s' entered; reverting to previous order number '%d'.", newOrder, group->order);
      char *orderStr = convertIntToString(group->order);
      gtk_entry_set_text(entry, orderStr);
      g_free(orderStr);
      result = FALSE;
    }
  else
    {
      group->order = convertStringToInt(newOrder);
      
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(widget));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      blxWindowGroupsChanged(blxWindow);
    }
  
  return result;
}


/* This callback is called when the dialog settings are applied. It sets the hidden
 * status of the passed groupo based on the toggle button's state */
static gboolean onGroupHiddenToggled(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  SequenceGroup *group = (SequenceGroup*)data;
  group->hidden = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

  /* Refilter trees and redraw all immediately show the new status */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));

  blxWindowGroupsChanged(blxWindow);
  
  return result;
}


/* This callback is called when the toggle button for a group's "highlighted" flag is toggled.
 * It updates the group's highlighted flag according to the button's new status. */
static gboolean onGroupHighlightedToggled(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  SequenceGroup *group = (SequenceGroup*)data;
  group->highlighted = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  /* Redraw the blixem window to immediately show the new status. The toplevel
   * parent of the button is the dialog, and the blixem window is the transient
   * parent of the dialog. */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  
  blxWindowRedrawAll(blxWindow);
  
  return result;
}


/* Called when the user has changed the color of a group in the 'edit groups' dialog */
static gboolean onGroupColorChanged(GtkWidget *button, const gint responseId, gpointer data)
{  
  gboolean result = TRUE;
  
  SequenceGroup *group = (SequenceGroup*)data;
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), &group->highlightColor);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), &group->highlightColor, TRUE, TRUE);
  
  /* Redraw everything in the new colors */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(button)));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  blxWindowRedrawAll(blxWindow);
  
  return result;
}


/* This function creates a widget that allows the user to edit the group
 * pointed to by the given list item, and adds it to the table container
 * widget at the given row. */
static void createEditGroupWidget(GtkWidget *blxWindow, SequenceGroup *group, GtkTable *table, const int row, const int xpad, const int ypad)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  /* Only show the special 'match set' group if it has some sequences */
  if (group != blxContext->matchSetGroup || g_list_length(group->seqList) > 0)
    {
      /* Show the group's name in a text box that the user can edit */
      GtkWidget *nameWidget = gtk_entry_new();
      gtk_entry_set_text(GTK_ENTRY(nameWidget), group->groupName);
      gtk_entry_set_activates_default(GTK_ENTRY(nameWidget), TRUE);
      widgetSetCallbackData(nameWidget, onGroupNameChanged, group);
      
      /* Add a check box for the 'hidden' flag */
      GtkWidget *isHiddenWidget = gtk_check_button_new();
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHiddenWidget), group->hidden);
      widgetSetCallbackData(isHiddenWidget, onGroupHiddenToggled, group);

      /* Add a check box for the 'highlighted' flag */
      GtkWidget *isHighlightedWidget = gtk_check_button_new();
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(isHighlightedWidget), group->highlighted);
      widgetSetCallbackData(isHighlightedWidget, onGroupHighlightedToggled, group);
      
      /* Show the group's order number in an editable text box */
      GtkWidget *orderWidget = gtk_entry_new();

      char *orderStr = convertIntToString(group->order);
      gtk_entry_set_text(GTK_ENTRY(orderWidget), orderStr);
      g_free(orderStr);
      
      gtk_entry_set_activates_default(GTK_ENTRY(orderWidget), TRUE);
      gtk_widget_set_size_request(orderWidget, 30, -1);
      widgetSetCallbackData(orderWidget, onGroupOrderChanged, group);

      /* Show the group's highlight color in a button that will also launch a color-picker */
      GtkWidget *colorButton = gtk_color_button_new_with_color(&group->highlightColor);
      widgetSetCallbackData(colorButton, onGroupColorChanged, group);
      
      /* Create a button that will delete this group */
      GtkWidget *deleteButton = gtk_button_new_from_stock(GTK_STOCK_DELETE);
      g_signal_connect(G_OBJECT(deleteButton), "clicked", G_CALLBACK(onButtonClickedDeleteGroup), group);
      
      /* Put everything in the table */
      gtk_table_attach(table, nameWidget,               1, 2, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHiddenWidget,   2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, isHighlightedWidget,      3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, orderWidget,              4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, colorButton,              5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(table, deleteButton,             6, 7, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
    }
}


/* Like blxSequenceGetColumn but also supports the group column (which
 * needs the context for its data) */
static const char* blxSequenceGetColumnData(const BlxSequence* const blxSeq, 
                                            const BlxColumnId columnId,
                                            const BlxViewContext *bc)
{
  const char *result = NULL;

  if (columnId == BLXCOL_GROUP)
    {
      const SequenceGroup *group = blxContextGetSequenceGroup(bc, blxSeq);

      if (group)
        result = group->groupName;
    }
  else
    {
      result = blxSequenceGetColumn(blxSeq, columnId);
    }
  
  return result;
}
          

/* Called for an entry from a list of BlxSequences. Compares the relevant data
 * from the BlxSequence (as indicated by the searchCol member of the SeqSearchData
 * struct) to the search string (also specified in the SeqSearchData). If it
 * matches, it appends the BlxSequence to the result list in the SeqSearchData. */
static void getSequencesThatMatch(gpointer listDataItem, gpointer data)
{
  /* Get the BlxSequence for this list item */
  BlxSequence *seq = (BlxSequence*)listDataItem;
  
  /* Get the search data */
  SeqSearchData *searchData = (SeqSearchData*)data;
  
  if (searchData->error)
    return; /* already hit an error so don't try any more */

  /* Get the relevant data for the search column */
  const char *dataToCompare = blxSequenceGetColumnData(seq, searchData->searchCol, searchData->bc);
  
  if (!dataToCompare)
    {
      if (searchData->error)
        {
          /* Default error message is not that useful, so replace it */
          reportAndClearIfError(&searchData->error, G_LOG_LEVEL_DEBUG);
          g_set_error(&searchData->error, BLX_ERROR, BLX_ERROR_INVALID_COLUMN, "Invalid search column.\n");
        }
      
      return;
    }
  
  /* Do the search */
  gboolean found = wildcardSearch(dataToCompare, searchData->searchStr);

  if (searchData->searchCol == BLXCOL_SEQNAME)
    {
      /* Sequence names have a variant number postfix; if not found, try
       * to match the text without this postfix. */
      if (!found && dataToCompare)
        {
          char *seqName = g_strdup(dataToCompare);
          char *cutPoint = strchr(seqName, '.');
          
          if (cutPoint)
            {
              *cutPoint = '\0';
              found = wildcardSearch(seqName, searchData->searchStr);
            }
          
          g_free(seqName);
        }
    }
  
  if (found)
    {
      /* Add this BlxSequence onto the result list. */
      searchData->matchList = g_list_prepend(searchData->matchList, seq);
    }  
}


/* If the given radio button is enabled, add a group based on the curently-
 * selected sequences. */
static gboolean onAddGroupFromSelection(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {  
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      
      if (g_list_length(blxContext->selectedSeqs) > 0)
        {
          GList *list = g_list_copy(blxContext->selectedSeqs); /* group takes ownership of this */
          createSequenceGroup(blxWindow, list, FALSE, NULL);
        }
      else
        {
          result = FALSE;
          g_critical("Warning: cannot create group; no sequences are currently selected");
        }
    }
  
  return result;
}


/* If the given radio button is enabled, add a group based on the search text
 * in the given text entry. This function finds the combo box on the dialog that
 * specifies which column to search by, and searches the relevant column for
 * the given search text. */
static gboolean onAddGroupFromText(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  /* Nothing to do if this radio button is not active */
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return result;
    }
  
  /* Get the dialog and main window */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));

  /* Get the search text from the text entry widget (which is passed as user data) */
  const char *inputText = getStringFromTextEntry(GTK_ENTRY(data));

  /* Find the combo box on the dialog (there should only be one), which tells
   * us which column to search. */
  GtkComboBox *combo = widgetGetComboBox(dialog);
  BlxColumnId searchCol = getColumnFromComboBox(combo);
  
  GError *error = NULL;
  GList *seqList = findSeqsFromColumn(blxWindow, inputText, searchCol, FALSE, FALSE, &error);

  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      result = FALSE;
    }

  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      createSequenceGroup(blxWindow, seqList, FALSE, inputText);
    }
  
  return result;
}


/* Utility function to take a list of search strings (newline separated), and
 * return a list of BlxSequences whose column data matches one of those search
 * strings. */
static GList* getSeqStructsFromSearchStringList(GList *searchStringList, 
                                                GList *seqList, 
                                                BlxViewContext *bc, 
                                                const BlxColumnId searchCol,
                                                GError **error)
{
  SeqSearchData searchData = {NULL, searchCol, bc, NULL, NULL};
  
  /* Loop through all the names in the input list */
  GList *nameItem = searchStringList;
  for ( ; nameItem; nameItem = nameItem->next)
    {
      /* Compare this name to all names in the sequence list. If it matches,
       * add it to the result list. */
      searchData.searchStr = (const char*)(nameItem->data);
      g_list_foreach(seqList, getSequencesThatMatch, &searchData);
      
      if (searchData.error)
        break;
    }

  if (searchData.error)
    g_propagate_error(error, searchData.error);
  
  return searchData.matchList;
}


/* Utility to create a GList of BlxSequences from a textual list of search strings.
 * Returns only valid sequences that blixem knows about. Looks for sequences
 * whose relevent data for the given column matches an item in the input text
 * (which may be a multi-line list; one search string per line). */
static GList* getSeqStructsFromText(GtkWidget *blxWindow, const char *inputText, const BlxColumnId searchCol, GError **error)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  GError *tmpError = NULL;

  GList *searchStringList = parseMatchList(inputText);

  /* Extract the entries from the list that are sequences that blixem knows about */
  GList *matchSeqs = blxWindowGetAllMatchSeqs(blxWindow);
  GList *seqList = getSeqStructsFromSearchStringList(searchStringList, matchSeqs, bc, searchCol, &tmpError);

  if (tmpError)
    g_propagate_error(error, tmpError);
  
  /* Must free the original name list and all its data. */
  freeStringList(&searchStringList, TRUE);
  
  if (g_list_length(seqList) < 1)
    {
      g_list_free(seqList);
      seqList = NULL;
    }
  
  return seqList;
}


/* Callback function to be used when requesting text from the clipboard to be used
 * to create the 'match set' group from the paste text */
static void createMatchSetFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data)
{
  /* Get the list of sequences to include */
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = getSeqStructsFromText(blxWindow, clipboardText, BLXCOL_SEQNAME, NULL);
  
  /* If a group already exists, replace its list. Otherwise create the group. */
  if (seqList)
    {
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      
      if (!blxContext->matchSetGroup)
        {
          blxContext->matchSetGroup = createSequenceGroup(blxWindow, seqList, FALSE, MATCH_SET_GROUP_NAME);
        }
      else
        {
          if (blxContext->matchSetGroup->seqList)
            {
              g_list_free(blxContext->matchSetGroup->seqList);
            }
          
          blxContext->matchSetGroup->seqList = seqList;
        }
      
      /* Reset the highlighted/hidden properties to make sure the group is initially visible */
      blxContext->matchSetGroup->highlighted = TRUE;
      blxContext->matchSetGroup->hidden = FALSE;

      blxWindowGroupsChanged(blxWindow);
      
      /* Refresh the groups dialog, if it happens to be open */
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
}


/* Callback function to be used when requesting text from the clipboard to be used
 * to find and select sequences based on the paste text. Scrolls the first selected
 * match into view. */
void findSeqsFromClipboard(GtkClipboard *clipboard, const char *clipboardText, gpointer data)
{
  /* Get the list of sequences to include */
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = getSeqStructsFromText(blxWindow, clipboardText, BLXCOL_SEQNAME, NULL);
  
  if (seqList)
    {
      blxWindowSetSelectedSeqList(blxWindow, seqList);
      firstMatch(blxWindowGetDetailView(blxWindow), seqList);
    }
}


/* This function toggles the match set.  That is, if the match set (a special 
 * group) exists then it deletes it; if it does not exist, then it creates it
 * from the current clipboard text (which should contain valid sequence name(s)). */
static void toggleMatchSet(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  
  if (blxContext->matchSetGroup && blxContext->matchSetGroup->seqList)
    {
      /* Clear the list of names only (don't delete the group, because we want to
       * keep any changes the user made (e.g. to the group color etc.) for next time. */
      g_list_free(blxContext->matchSetGroup->seqList);
      blxContext->matchSetGroup->seqList = NULL;
      blxWindowGroupsChanged(blxWindow);

      /* Refresh the groups dialog, if it happens to be open */
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
  else
    {
      requestPrimaryClipboardText(createMatchSetFromClipboard, blxWindow);
    }
}


/* If the given radio button is enabled, add a group based on the list of sequences
 * in the given text entry. */
static gboolean onAddGroupFromList(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      return result;
    }
  
  /* Get the dialog and main window */
  GtkWidget *dialog = gtk_widget_get_toplevel(button);
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
  
  /* Find the combo box on the dialog (there should only be one), which tells
   * us which column to search. */
  GtkComboBox *combo = widgetGetComboBox(dialog);
  BlxColumnId searchCol = getColumnFromComboBox(combo);
  
  /* The text entry box was passed as the user data. We should have a (multi-line) text view */
  const char *inputText = getStringFromTextView(GTK_TEXT_VIEW(data));

  GError *error = NULL;
  GList *seqList = findSeqsFromList(blxWindow, inputText, searchCol, FALSE, FALSE, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      result = FALSE;
    }
  
  if (seqList)
    {
      GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
      GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

      createSequenceGroup(blxWindow, seqList, FALSE, NULL);
    }
  
  return result;
}


/* Called when the user has clicked the "delete all groups" button on the "group sequences" dialog. */
static void onButtonClickedDeleteAllGroups(GtkWidget *button, gpointer data)
{
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));

  /* Ask the user if they're sure */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  char *title = g_strdup_printf("%sDelete groups", blxGetTitlePrefix(bc));
  gint response = runConfirmationBox(blxWindow, title, "This will delete ALL groups. Are you sure?");
  g_free(title);

  if (response == GTK_RESPONSE_ACCEPT)
    {
      blxWindowDeleteAllSequenceGroups(blxWindow);
      
      /* Close the dialog, because there are no groups left to display. */
      gtk_widget_hide_all(GTK_WIDGET(dialogWindow));
    }
}


/* Called when the user has clicked the "delete group" button on the "group sequences" dialog. */
static void onButtonClickedDeleteGroup(GtkWidget *button, gpointer data)
{
  SequenceGroup *group = (SequenceGroup*)data;
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  
  /* Ask the user if they're sure */
  char formatStr[] = "Are you sure you wish to delete group '%s'?";
  char messageText[strlen(formatStr) + strlen(group->groupName)];
  sprintf(messageText, formatStr, group->groupName);
  
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  char *title = g_strdup_printf("%sDelete group", blxGetTitlePrefix(bc));
  gint response = runConfirmationBox(blxWindow, title, messageText);
  g_free(title);
  
  if (response == GTK_RESPONSE_ACCEPT)
    {
      blxWindowDeleteSequenceGroup(blxWindow, group);
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
}


/* Called when the user chooses a different tabe on the groups dialog */
static gboolean onSwitchPageGroupsDialog(GtkNotebook *notebook, GtkNotebookPage *page, guint pageNum, gpointer data)
{
  GtkDialog *dialog = GTK_DIALOG(data);

  if (pageNum == 0)
    {
      /* For the create-groups page, set the default response to be 'accept' */
      gtk_dialog_set_default_response(dialog, GTK_RESPONSE_ACCEPT);
    }
  else
    {
      /* For other pages, set default response to be 'apply' */
      gtk_dialog_set_default_response(dialog, GTK_RESPONSE_APPLY);
    }
    
  return FALSE;
}


/* Utility to find and return a notebook child of the given widget. Assumes there is only
 * one - if there are more it will just return the first found. Returns NULL if not found. */
static GtkNotebook* containerGetChildNotebook(GtkContainer *container)
{
  GtkNotebook *result = NULL;
  
  GList *children = gtk_container_get_children(container);
  GList *child = children;
  
  for ( ; child; child = child->next)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);
      
      if (GTK_IS_NOTEBOOK(childWidget))
        {
          result = GTK_NOTEBOOK(childWidget);
          break;
        }
      else if (GTK_IS_CONTAINER(childWidget))
       {
         /* recurse */
         containerGetChildNotebook(GTK_CONTAINER(childWidget));
       }
    }

  g_list_free(children);
 
  return result;
}


/* Callback called when user responds to groups dialog */
void onResponseGroupsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = FALSE;
  gboolean refresh = FALSE;
  
  /* If a notebook was passed, only call callbacks for widgets in the active tab */
  GtkNotebook *notebook = containerGetChildNotebook(GTK_CONTAINER(dialog->vbox));
               
  if (!notebook)
    {
      g_warning("Expected Groups dialog to contain a notebook widget. Dialog may not refresh properly.\n");
    }

  guint pageNo = notebook ? gtk_notebook_get_current_page(notebook) : 0;
  GtkWidget *page = notebook ? gtk_notebook_get_nth_page(notebook, pageNo) : NULL;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      destroy = widgetCallAllCallbacks(page, GINT_TO_POINTER(responseId));
      refresh = FALSE;
      break;

    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(page, GINT_TO_POINTER(responseId));
      destroy = FALSE;
      refresh = (pageNo == 0); /* if created a new group, Edit Groups section must be refreshed */
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      refresh = FALSE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      /* Groups dialog is persistent, so hide it rather than destroying it */
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
  else if (refresh)
    {
      GtkWidget *blxWindow = GTK_WIDGET(data);
      refreshDialog(BLXDIALOG_GROUPS, blxWindow);
    }
}


/* Callback for when the 'hide ungrouped sequences' option is changed */
static gboolean onHideUngroupedChanged(GtkWidget *button, const gint responseId, gpointer data)
{
  setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  refilterDetailView(detailView, NULL);

  calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
  calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
  forceResize(bigPicture);

  blxWindowRedrawAll(blxWindow);
  
  return TRUE;
}


/* Create the 'create group' tab of the groups dialog. Appends it to the notebook. */
static void createCreateGroupTab(GtkNotebook *notebook, BlxViewContext *bc, GtkWidget *blxWindow)
{
  const int numRows = 3;
  const int numCols = 2;
  const gboolean seqsSelected = g_list_length(bc->selectedSeqs) > 0;
  
  /* Put everything in a table */
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));

  /* Append the table as a new tab to the notebook */
  gtk_notebook_append_page(notebook, GTK_WIDGET(table), gtk_label_new("Create group"));
  
  /* Create the left-hand-side column */
  GtkRadioButton *button1 = createRadioButton(table, 1, 1, NULL, "_Text search (wildcards * and ?)", !seqsSelected, TRUE, FALSE, onAddGroupFromText, blxWindow, NULL);
  createRadioButton(table, 1, 2, button1, "_List search", FALSE, TRUE, TRUE, onAddGroupFromList, blxWindow, NULL);
  createSearchColumnCombo(table, 1, 3, blxWindow);
  
  /* Create the right-hand-side column */
  createRadioButton(table, 2, 1, button1, "Use current _selection", seqsSelected, FALSE, FALSE, onAddGroupFromSelection, blxWindow, NULL);
}


/* Create the 'edit groups' tab of the groups dialog. Appends it to the given notebook. */
static void createEditGroupsTab(GtkNotebook *notebook, BlxViewContext *bc, GtkWidget *blxWindow)
{
  const int numRows = g_list_length(bc->sequenceGroups) + 4; /* +4 for: header; delete-all button;
                                                                hide-all-seqs; hide-all-features */
  const int numCols = 6;
  const int xpad = DEFAULT_TABLE_XPAD;
  const int ypad = DEFAULT_TABLE_YPAD;
  int row = 1;
  
  /* Put everything in a table */
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  
  /* Append the table as a new tab to the notebook */
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(table), gtk_label_new("Edit groups"));
  
  /* Add check buttons to turn on the 'hide ungrouped sequences/features' options */
  GtkWidget *hideButton1 = gtk_check_button_new_with_mnemonic("_Hide all sequences not in a group");
  GtkWidget *hideButton2 = gtk_check_button_new_with_mnemonic("_Hide all features not in a group");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hideButton1), bc->flags[BLXFLAG_HIDE_UNGROUPED_SEQS]);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hideButton2), bc->flags[BLXFLAG_HIDE_UNGROUPED_FEATURES]);
  widgetSetCallbackData(hideButton1, onHideUngroupedChanged, GINT_TO_POINTER(BLXFLAG_HIDE_UNGROUPED_SEQS));
  widgetSetCallbackData(hideButton2, onHideUngroupedChanged, GINT_TO_POINTER(BLXFLAG_HIDE_UNGROUPED_FEATURES));

  gtk_table_attach(table, hideButton1, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;
  gtk_table_attach(table, hideButton2, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
  ++row;
  
  /* Add labels for each column in the table */
  gtk_table_attach(table, gtk_label_new("Group name"),    1, 2, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Hide"),          2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Highlight"),     3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Order"),         4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  ++row;
  
  /* Add a set of widgets for each group */
  GList *groupItem = blxWindowGetSequenceGroups(blxWindow);
  for ( ; groupItem; groupItem = groupItem->next)
    {
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      createEditGroupWidget(blxWindow, group, table, row, xpad, ypad);
      ++row;
    }
  
  /* Add a button to delete all groups */
  GtkWidget *deleteGroupsButton = gtk_button_new_with_label("Delete all groups");
  gtk_button_set_image(GTK_BUTTON(deleteGroupsButton), gtk_image_new_from_stock(GTK_STOCK_DELETE, GTK_ICON_SIZE_BUTTON));
  gtk_widget_set_size_request(deleteGroupsButton, -1, 30);
  g_signal_connect(G_OBJECT(deleteGroupsButton), "clicked", G_CALLBACK(onButtonClickedDeleteAllGroups), NULL);
  gtk_table_attach(table, deleteGroupsButton, numCols - 1, numCols + 1, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);
}


/* Shows the "Group sequences" dialog. This dialog allows the user to group sequences together.
 * This tabbed dialog shows both the 'create group' and 'edit groups' dialogs in one. If the
 * 'editGroups' argument is true and groups exist, the 'Edit Groups' tab is displayed by default;
 * otherwise the 'Create Groups' tab is shown. */
void showGroupsDialog(GtkWidget *blxWindow, const gboolean editGroups, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_GROUPS;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sGroups", blxGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_APPLY,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      /* Make sure we only connect the response event once */
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseGroupsDialog), blxWindow);
    }
  else
    {
      /* Refresh by deleting the dialog contents and re-creating them. */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  /* Create tabbed pages */
  GtkWidget *notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, TRUE, TRUE, 0);

  createCreateGroupTab(GTK_NOTEBOOK(notebook), bc, blxWindow);
  createEditGroupsTab(GTK_NOTEBOOK(notebook), bc, blxWindow);

  
  /* Connect signals and show */
  gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(blxWindow));
  g_signal_connect(notebook, "switch-page", G_CALLBACK(onSwitchPageGroupsDialog), dialog);
  
  gtk_widget_show_all(dialog);

  if (editGroups && notebook && blxWindowGroupsExist(blxWindow))
    {
      gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), 1); /* 'edit' page is the 2nd page */
    }
  else
    {
      gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), 0); /* 'create' page is the 1st page */
    }
  
  if (bringToFront)
    {
      /* If user has asked to edit groups (and some groups exist), make the second tab
       * the default and the 'close' button the default action. (Must do this after showing
       * the child widgets due to a GTK legacy whereby the notebook won't change tabs otherwise.) */
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/***********************************************************
 *                         Settings menu                   *
 ***********************************************************/

/* This function should be called on the child widget of a dialog box that is a transient
 * child of the main blixem window. It finds the parent dialog of the child and then finds
 * the blxWindow from the dialog. */
static GtkWidget* dialogChildGetBlxWindow(GtkWidget *child)
{
  GtkWidget *result = NULL;
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(child));
  
  if (dialogWindow)
    {
      result = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
    }
  
  return result;
}


/* Updates the given flag from the given button. The passed in widget is the toggle button and
 * the data is an enum indicating which flag was toggled. Returns the new value that was set. */
static gboolean setFlagFromButton(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  BlxFlag flag = (BlxFlag)GPOINTER_TO_INT(data);
  const gboolean newValue = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

  if (flag > BLXFLAG_MIN && flag < BLXFLAG_NUM_FLAGS)
    bc->flags[flag] = newValue;
  
  return newValue;
}


/* This callback is called when one of the boolean flags is toggled on the settings dialog.
 * This generic function sets the flag and redraws everything; if different 
 * updates are required a custom callback function can be used instead. */
static void onToggleFlag(GtkWidget *button, gpointer data)
{
  setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  blxWindowRedrawAll(blxWindow);
}


/* Callback function called when the 'squash matches' button is toggled */
static void onSquashMatches(GtkWidget *button, gpointer data)
{
  const gboolean squash = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));

  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  
  setToggleMenuStatus(properties->actionGroup, "SquashMatches", squash);
}


/* Callback function called when the 'Show Variation track' button is toggled */
static void onShowVariationTrackToggled(GtkWidget *button, gpointer data)
{
  const gboolean showSnpTrack = setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  detailViewUpdateShowSnpTrack(detailView, showSnpTrack);
}


/* Utility to create a check button with certain given properties, and to pack it into the parent.
 * Returns the new button. */
static GtkWidget* createCheckButton(GtkBox *box, 
                                    const char *mnemonic, 
                                    const gboolean isActive, 
                                    GCallback callback, 
                                    gpointer data)
{
  GtkWidget *button = gtk_check_button_new_with_mnemonic(mnemonic);
  gtk_box_pack_start(box, button, FALSE, FALSE, 0);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), isActive);
  
  g_signal_connect(G_OBJECT(button), "toggled", callback, data);
  
  return button;
}


/* Callback to be called when the user has entered a new column size */
static gboolean onColumnSizeChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkEntry *entry = GTK_ENTRY(widget);
  BlxColumnInfo *columnInfo = (BlxColumnInfo*)data;
  
  const gchar *newSizeText = gtk_entry_get_text(entry);
  const int newWidth = convertStringToInt(newSizeText);
  
  if (newWidth != columnInfo->width)
    {
      /* Check it's a sensible value. We could do with a better check really but
       * for now just check that it's less than the screen width. This at least
       * catches excessively large values, which can cause Blixem to crash. Slightly
       * too-large values may make things look odd but should be recoverable. */
      GtkWidget *blxWindow = dialogChildGetBlxWindow(widget);
      GdkScreen *screen = gtk_widget_get_screen(blxWindow);

      if (newWidth > gdk_screen_get_width(screen))
        {
          g_critical("Column width '%d' too large; not changed.\n", newWidth);
        }
      else
        {
          columnInfo->width = newWidth;

          GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
          updateDynamicColumnWidths(detailView);
        }
    }
  
  return result;
}


/* Callback to be called when the user has toggled the visibility of a column */
static gboolean onColumnVisibilityChanged(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  BlxColumnInfo *columnInfo = (BlxColumnInfo*)data;
  columnInfo->showColumn = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);

  updateDynamicColumnWidths(detailView);
  
  return result;
}


/* Callback to be called when the user has toggled which columns are included in summary info */
static gboolean onSummaryColumnsChanged(GtkWidget *button, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  BlxColumnInfo *columnInfo = (BlxColumnInfo*)data;
  columnInfo->showSummary = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);

  /* Just clear the moused-over feedback area to make sure it's not showing invalid data. The
   * user can mouse-over again to see the new data. */
  clearFeedbackArea(detailView);
  
  return result;
}


/* Just calls gtk_widget_set_sensitive, but so that we can use it as a callback */
static void widgetSetSensitive(GtkWidget *widget, gpointer data)
{
  gboolean sensitive = GPOINTER_TO_INT(data);
  gtk_widget_set_sensitive(widget, sensitive);
}


/* Callback when the user hits the 'load optional data' button on the settings dialog */
static void onButtonClickedLoadOptional(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Create a temporary lookup table for BlxSequences so we can link them on GFF ID */
  GHashTable *lookupTable = g_hash_table_new(g_direct_hash, g_direct_equal);
  
  GError *error = NULL;
  gboolean success = bulkFetchSequences(
    0, bc->external, bc->flags[BLXFLAG_SAVE_TEMP_FILES],
    bc->seqType, &bc->matchSeqs, bc->columnList, bc->optionalFetchDefault, bc->fetchMethods, &bc->mspList,
    &bc->blastMode, bc->featureLists, bc->supportedTypes, NULL, bc->refSeqOffset,
    &bc->refSeqRange, bc->dataset, TRUE, lookupTable);
  
  finaliseFetch(bc->matchSeqs, bc->columnList);

  if (error)
    {
      prefixError(error, "Error loading optional data. ");
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  
  if (success)
    {
      /* Set the flag to say that the data has now been loaded */
      bc->flags[BLXFLAG_OPTIONAL_COLUMNS] = TRUE;
      
      /* Disable the button so user can't try to load data again. */
      gtk_widget_set_sensitive(button, FALSE);
      
      /* Enable the text entry boxes for all columns. They are all in the container passed
       * as the user data. */
      GtkContainer *container = GTK_CONTAINER(data);
      gtk_container_foreach(container, widgetSetSensitive, GINT_TO_POINTER(TRUE));
      
      /* Update the flag in all columns to indicate that the data is now loaded */
      GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
      GList *listItem = detailViewGetColumnList(detailView);
      
      for ( ; listItem; listItem = listItem->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);
          columnInfo->dataLoaded = TRUE;
        }   
      
      /* Re-sort the trees, because the new data may affect the sort order. Also
       * resize them, because whether data is loaded affects whether columns are shown. */
      detailViewResortTrees(detailView);
      updateDynamicColumnWidths(detailView);

      /* Force a of resize the tree columns (updateDynamicColumnWidths won't resize them
       * because the widths and visibility-flags haven't actually changed, but visibility
       * IS affected because we've now loaded the data) */
      callFuncOnAllDetailViewTrees(detailView, resizeTreeColumns, NULL);
      resizeDetailViewHeaders(detailView);
      updateDetailViewRange(detailView);
      
      detailViewRedrawAll(detailView);
    }

  g_hash_table_unref(lookupTable);
}


/* Create a button to allow the user to load the data for optional columns, if not already loaded */
static GtkWidget* createColumnLoadDataButton(GtkTable *table, 
                                             GtkWidget *detailView,
                                             int *row,
                                             const int cols,
                                             int xpad,
                                             int ypad)
{
  BlxViewContext *bc = blxWindowGetContext(detailViewGetBlxWindow(detailView));
  const gboolean dataLoaded = bc->flags[BLXFLAG_OPTIONAL_COLUMNS];

  /* Create the button */
  GtkWidget *button = gtk_button_new_with_label(LOAD_DATA_TEXT);
  gtk_widget_set_sensitive(button, !dataLoaded); /* only enable if data not yet loaded */

  /* Add the button to the table, spanning all of the remaining columns */
  gtk_table_attach(table, button, 0, 1, *row, *row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create an explanatory label spanning the rest of the columns */
  GtkWidget *label = gtk_label_new("Fetches additional information e.g. from an\nEMBL file (as determined by the config).");
  gtk_table_attach(table, label, 1, cols, *row, *row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);

  *row += 1;

  return button;
}


/* Create the settings buttons for a single column */
static void createColumnButton(BlxColumnInfo *columnInfo, GtkTable *table, int *row)
{
  /* Create a label showing the column name */
  GtkWidget *label = gtk_label_new(columnInfo->title);
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
  
  /* Tick-box controlling whether column is displayed */
  GtkWidget *showColButton = gtk_check_button_new();
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(showColButton), columnInfo->showColumn);
  widgetSetCallbackData(showColButton, onColumnVisibilityChanged, (gpointer)columnInfo);

  /* Tick-box controlling whether column is included in summary info */
  GtkWidget *showSummaryButton = gtk_check_button_new();
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(showSummaryButton), columnInfo->showSummary);
  widgetSetCallbackData(showSummaryButton, onSummaryColumnsChanged, (gpointer)columnInfo);
  
  GtkWidget *entry = gtk_entry_new();
  
  if (columnInfo->columnId == BLXCOL_SEQUENCE)
    {
      /* The sequence column updates dynamically, so don't allow the user to edit it */
      char displayText[] = "<dynamic>";
      gtk_entry_set_text(GTK_ENTRY(entry), displayText);
      gtk_widget_set_sensitive(entry, FALSE);
      gtk_widget_set_sensitive(showSummaryButton, FALSE);
      gtk_widget_set_sensitive(showSummaryButton, FALSE);
      gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2); /* fudge up width a bit in case user enters longer text */
    }
  else
    {
      /* Grey out all the options if the column hasn't been loaded */
      if (!columnInfo->dataLoaded)
        {
          gtk_widget_set_sensitive(showColButton, FALSE);
          gtk_widget_set_sensitive(showSummaryButton, FALSE);
          gtk_widget_set_sensitive(entry, FALSE);
        }
      else if (!columnInfo->canShowSummary)
        {
          gtk_widget_set_sensitive(showSummaryButton, FALSE);          
        }

      char *displayText = convertIntToString(columnInfo->width);
      gtk_entry_set_text(GTK_ENTRY(entry), displayText);
      
      gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
      gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
      
      widgetSetCallbackData(entry, onColumnSizeChanged, (gpointer)columnInfo);
      
      g_free(displayText);
    }
  
  gtk_table_attach(table, label,             0, 1, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);
  gtk_table_attach(table, showColButton,     1, 2, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);
  gtk_table_attach(table, showSummaryButton, 2, 3, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);
  gtk_table_attach(table, entry,             3, 4, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);
  *row += 1;
}


/* Create labels for the column properties widgets created by createColumnButton */
static void createColumnButtonHeaders(GtkTable *table, int *row)
{
  GtkWidget *label = gtk_label_new("Show\ncolumn");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
  gtk_table_attach(table, label, 1, 2, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);

  label = gtk_label_new("Show mouse-\nover details");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
  gtk_table_attach(table, label, 2, 3, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);

  label = gtk_label_new("Column width");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 1.0);
  gtk_table_attach(table, label, 3, 4, *row, *row + 1, GTK_FILL, GTK_SHRINK, 4, 4);

  *row += 1;
}

/* Create a set of widgets that allow columns settings to be adjusted */
static void createColumnButtons(GtkWidget *parent, GtkWidget *detailView, const int border)
{
  /* put all the column settings in a table. Put the table in a
   * scrolled window, because there are likely to be many rows */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(parent), scrollWin);
  
  GList *columnList = detailViewGetColumnList(detailView);
  const int rows = g_list_length(columnList) + 1;
  const int cols = 4;
  int row = 1;

  GtkTable *table = GTK_TABLE(gtk_table_new(rows, cols, FALSE));  

  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Create a button to allow the user to load the full EMBL data, if not already loaded */
  GtkWidget *button = createColumnLoadDataButton(table, detailView, &row, cols, border, border);
  
  /* Loop through each column and create widgets to control the properties */
  createColumnButtonHeaders(table, &row);

  GList *listItem = columnList;
  for ( ; listItem; listItem = listItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);
      createColumnButton(columnInfo, table, &row);
    }
  
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(onButtonClickedLoadOptional), table);
}


/* Callback to be called when the user has entered a new percent-ID per cell */
static gboolean onIdPerCellChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  return bigPictureSetIdPerCell(bigPicture, newValue);
}

 
/* Callback to be called when the user has entered a new maximum percent-ID to display */
static gboolean onMaxPercentIdChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  return bigPictureSetMaxPercentId(bigPicture, newValue);
}

/* Callback to be called when the user has entered a new minimum percent-ID to display */
static gboolean onMinPercentIdChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *bigPicture = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  return bigPictureSetMinPercentId(bigPicture, newValue);
}


/* Callback to be called when the user has changed the depth-per-cell on the coverage view */
static gboolean onDepthPerCellChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *coverageView = GTK_WIDGET(data);  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  const gdouble newValue = g_strtod(text, NULL);
  return coverageViewSetDepthPerCell(coverageView, newValue);
}


///* Utility to create a text entry widget displaying the given integer value. The
// * given callback will be called when the user OK's the dialog that this widget 
// * is a child of. */
//static void createTextEntryFromInt(GtkWidget *parent, 
//                                 const char *title, 
//                                 const int value, 
//                                 BlxResponseCallback callbackFunc, 
//                                 gpointer callbackData)
//{
//  /* Pack label and text entry into a vbox */
//  GtkWidget *vbox = createVBoxWithBorder(parent, 4, FALSE, NULL);
//
//  GtkWidget *label = gtk_label_new(title);
//  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
//
//  GtkWidget *entry = gtk_entry_new();
//  gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
//
//  char *displayText = convertIntToString(value);
//  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
//
//  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
//  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
//
//  widgetSetCallbackData(entry, callbackFunc, callbackData);
//  
//  g_free(displayText);
//}


/* Utility to create a text entry widget displaying the given double value. The
 * given callback will be called when the user OK's the dialog that this widget 
 * is a child of. */
static void createTextEntry(GtkWidget *parent, 
                            const char *title, 
                            const gdouble value, 
                            BlxResponseCallback callbackFunc, 
                            gpointer callbackData)
{
  /* Pack label and text entry into a vbox */
  GtkWidget *vbox = createVBoxWithBorder(parent, 4, FALSE, NULL);
  
  GtkWidget *label = gtk_label_new(title);
  gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
  
  GtkWidget *entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);
  
  char *displayText = convertDoubleToString(value, 1);
  gtk_entry_set_text(GTK_ENTRY(entry), displayText);
  
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(displayText) + 2);
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  widgetSetCallbackData(entry, callbackFunc, callbackData);
  
  g_free(displayText);
}


/* Create a set of widgets to allow the user to edit grid properties */
static void createGridSettingsButtons(GtkWidget *parent, GtkWidget *bigPicture)
{
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Overview section");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);

  /* Arrange the widgets horizontally */
  GtkWidget *hbox = createHBoxWithBorder(frame, 12, FALSE, NULL);
  const DoubleRange* const percentIdRange = bigPictureGetPercentIdRange(bigPicture);
  
  createTextEntry(hbox, "%ID per cell", bigPictureGetIdPerCell(bigPicture), onIdPerCellChanged, bigPicture);
  createTextEntry(hbox, "Max %ID", percentIdRange->max, onMaxPercentIdChanged, bigPicture);
  createTextEntry(hbox, "Min %ID", percentIdRange->min, onMinPercentIdChanged, bigPicture);
}


/* Create a set of widgets to allow the user to edit coverage-view properties */
static void createCoverageSettingsButtons(GtkWidget *parent, GtkWidget *bigPicture)
{
  BigPictureProperties *properties = bigPictureGetProperties(bigPicture);
  
  /* Group these buttons in a frame */
  GtkWidget *frame = gtk_frame_new("Coverage section");
  gtk_box_pack_start(GTK_BOX(parent), frame, FALSE, FALSE, 0);
  
  /* Arrange the widgets horizontally */
  GtkWidget *hbox = createHBoxWithBorder(frame, 12, FALSE, NULL);
  const double rangePerCell = coverageViewGetDepthPerCell(properties->coverageView);
  
  createTextEntry(hbox, "Depth per cell", rangePerCell, onDepthPerCellChanged, properties->coverageView);
}


/* Callback called when user has changed a blixem color */
static gboolean onChangeBlxColor(GtkWidget *button, const gint responseId, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  
  /* update the color */
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), color);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), color, TRUE, TRUE);
  
  /* Redraw */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  blxWindowRedrawAll(blxWindow);
  
  return TRUE;
}


/* Callback called when user has changed the blixem background color */
static gboolean onChangeBackgroundColor(GtkWidget *button, const gint responseId, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  
  /* update the color */
  gtk_color_button_get_color(GTK_COLOR_BUTTON(button), color);
  gdk_colormap_alloc_color(gdk_colormap_get_system(), color, TRUE, TRUE);
  
  /* Update */
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(button));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  onUpdateBackgroundColor(blxWindow);
  
  return TRUE;
}


/* Create a button to allow user to change the color of the given setting */
static void createColorButton(GtkTable *table, GdkColor *color, BlxResponseCallback callbackFunc, gpointer callbackData,
                              const int row, const int column, const int xpad, const int ypad)
{
  GtkWidget *colorButton = gtk_color_button_new_with_color(color);
  widgetSetCallbackData(colorButton, callbackFunc, callbackData);

  gtk_table_attach(table, colorButton, column, column + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
}


/* Create buttons for the user to be able to change the blixem colour settings */
static void createColorButtons(GtkWidget *parent, GtkWidget *blxWindow, const int borderWidth)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  
  /* put all the colors in a table. Put the table in a scrolled window, because
   * there are likely to be many rows */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(parent), scrollWin);
  
  const int numCols = 5;
  const int numRows = BLXCOL_NUM_COLORS + 1; /* add one for header row */
  const int xpad = 2;
  const int ypad = 2;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Add a header row */
  int row = 1;
  gtk_table_attach(table, gtk_label_new("Normal   "), 2, 3, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("(selected)   "), 3, 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("Print   "), 4, 5, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_table_attach(table, gtk_label_new("(selected)   "), 5, 6, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* loop through all defined blixem colours */
  int colorId = BLXCOLOR_MIN + 1;
  
  for ( ; colorId < BLXCOL_NUM_COLORS; ++colorId)
    {
      ++row;
    
      BlxColor *blxCol = getBlxColor(bc->defaultColors, colorId);
      GtkWidget *label = gtk_label_new(blxCol->name);
      gtk_table_attach(table, label, 1, 2, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);

      /* Special callback for the background color */
      BlxResponseCallback callbackFunc = (colorId == BLXCOLOR_BACKGROUND) ? onChangeBackgroundColor : onChangeBlxColor;
      
      createColorButton(table, &blxCol->normal, callbackFunc, &blxCol->normal, row, 2, xpad, ypad);
      createColorButton(table, &blxCol->print, callbackFunc, &blxCol->print, row, 4, xpad, ypad);
      createColorButton(table, &blxCol->selected, callbackFunc, &blxCol->selected, row, 3, xpad, ypad);
      createColorButton(table, &blxCol->printSelected, callbackFunc, &blxCol->printSelected, row, 5, xpad, ypad);
    }
}


/* Called when the user toggles whether print colors should be used or not */
static void onTogglePrintColors(GtkWidget *button, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  const gboolean usePrintColors = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  blxWindowSetUsePrintColors(blxWindow, usePrintColors);
}


/* Callback function called when the parent button of a set of sub-buttons is toggled. Enables/disables
 * the child buttons depending on whether the parent is now active or not. The container widget of the
 * child buttons is passed as the user data. */
static void onParentBtnToggled(GtkWidget *button, gpointer data)
{
  const gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  GtkWidget *subComponents = GTK_WIDGET(data);

  gtk_widget_set_sensitive(subComponents, active); 
}


/* Callback function called when the 'Show unaligned bases' or 'Show polyA tails'
 * buttons are toggled. These are parent buttons so will cause child buttons 
 * to be enabled/disabled. */
static void onShowAdditionalSeqToggled(GtkWidget *button, gpointer data)
{
  onParentBtnToggled(button, data);
  
  /* Perform any required updates */
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  detailViewUpdateMspLengths(detailView, detailViewGetNumUnalignedBases(detailView));
}


/* Callback function called when the 'Limit unaligned bases' button is toggled */
static void onLimitUnalignedBasesToggled(GtkWidget *button, gpointer data)
{
  /* Get the new value */
  const gboolean limitUnalignedBases = setFlagFromButton(button, GINT_TO_POINTER(BLXFLAG_LIMIT_UNALIGNED_BASES));

  /* Enable/disable the sub-options. Their widgets are all in the container passed as the data. */
  GtkWidget *subComponents = GTK_WIDGET(data);
  gtk_widget_set_sensitive(subComponents, limitUnalignedBases); 
  
  /* Get the detail view from the main window */
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  /* Perform any required updates */
  detailViewUpdateMspLengths(detailView, detailViewGetNumUnalignedBases(detailView));
}


/* Callback called when the user has changed the number of additional bases to show when the
 * 'show unaligned bases' option is enabled. */
static gboolean onSetNumUnalignedBases(GtkWidget *entry, const gint responseId, gpointer data)
{
  const char *numStr = gtk_entry_get_text(GTK_ENTRY(entry));
  int numBases = convertStringToInt(numStr);
  
  GtkWidget *detailView = GTK_WIDGET(data);
  detailViewSetNumUnalignedBases(detailView, numBases);

  /* Perform any required updates */
  detailViewUpdateMspLengths(detailView, numBases);

  return TRUE;
}


/* Callback called when the 'selected seqs only' option of the 'show unaligned bases'
 * option is toggled. */
static void onToggleShowUnalignedSelected(GtkWidget *button, gpointer data)
{
  /* Get the new value */
  setFlagFromButton(button, GINT_TO_POINTER(BLXFLAG_SHOW_UNALIGNED_SELECTED));
 
  /* Perform any required updates */
  GtkWidget *detailView = GTK_WIDGET(data);
//  detailViewUpdateMspLengths(detailView, detailViewGetNumUnalignedBases(detailView));
  refilterDetailView(detailView, NULL);
}



/* Create the check button for the 'limit number of unaligned bases' option on the settings dialog
 * and pack it into the given container. */
static void createLimitUnalignedBasesButton(GtkContainer *parent, GtkWidget *detailView, BlxViewContext *bc)
{
  /* Create an hbox for the "limit to so-many bases" option, which has a check button, text
   * entry and some labels. Pack the hbox into the given parent. */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_add(parent, hbox);
  
  /* Create a text entry box so the user can enter the number of bases */
  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  widgetSetCallbackData(entry, onSetNumUnalignedBases, detailView);
  
  DetailViewProperties *properties = detailViewGetProperties(detailView);
  char *numStr = convertIntToString(properties->numUnalignedBases);
  
  gtk_entry_set_text(GTK_ENTRY(entry), numStr);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(numStr) + 2); /* fudge up width a bit in case user enters longer text */
  g_free(numStr);
  
  /* Check button to enable/disable setting the limit */
  GtkWidget *button = gtk_check_button_new_with_mnemonic("Li_mit to ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), bc->flags[BLXFLAG_LIMIT_UNALIGNED_BASES]);
  g_signal_connect(G_OBJECT(button), "toggled", G_CALLBACK(onLimitUnalignedBasesToggled), entry);
  
  /* Pack it all in the hbox */
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(" additional bases"), FALSE, FALSE, 0);
}


/* Create a "parent" option button that has a vbox container for "sub-components", i.e. more 
 * option buttons (or other dialog widgets) that will be enabled only when the parent option is
 * active. Returns the container for the sub-options, which should be packed with the sub-option widgets. */
static GtkContainer* createParentCheckButton(GtkWidget *parent, 
                                             GtkWidget *detailView,
                                             BlxViewContext *bc,
                                             const char *label,
                                             const BlxFlag flag,
                                             GtkWidget **buttonOut,
                                             GCallback callbackFunc)
{
  /* We'll the main button and any sub-components into a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(parent), vbox, FALSE, FALSE, 0);

  /* Create a vbox for the sub-components. Create the vbox now so we can pass it to the main toggle
   * button callback, but don't pack it in the container till we've added the main check button. The sub
   * components are only active if the main check button is active. */
  const gboolean active = bc->flags[flag];
  GtkWidget *subContainer = gtk_vbox_new(FALSE, 0);
  gtk_widget_set_sensitive(subContainer, active); 
  
  /* Main check button to enable/disable the option. This call puts it in the vbox. Set two callbacks:
   * one to update the flag, and one to enable/disable the child buttons. */
  GtkWidget *btn = gtk_check_button_new_with_mnemonic(label);
  gtk_box_pack_start(GTK_BOX(vbox), btn, FALSE, FALSE, 0);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(btn), active);

  /* Connect the toggleFlag callback first so that the flag is set correctly before the callbackFunc is called */
  g_signal_connect(G_OBJECT(btn), "toggled", G_CALLBACK(onToggleFlag), GINT_TO_POINTER(flag));
  
  if (callbackFunc)
    g_signal_connect(G_OBJECT(btn), "toggled", callbackFunc, subContainer);

  if (buttonOut)
    *buttonOut = btn;

  /* Now add the subcomponent container to the vbox. Bit of a hack - put it inside an hbox with 
   * a blank label preceeding it, so that the sub-components appear offset to the right slightly
   * from the main button. */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("   "), FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), subContainer, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  
  return GTK_CONTAINER(subContainer);
}


/* Refresh the given dialog, if it is open */
void refreshDialog(const BlxDialogId dialogId, GtkWidget *blxWindow)
{
  /* This is a bit crude but does the job: if the dialog is visible, just call its
   * 'show' function to re-create its contents. Only need to do anything for persistent
   * dialogs. Note that we don't want to bring the dialog to the front, just refresh it in
   * case the user looks at it again. */
  if (blxWindow)
    {
      BlxViewContext *bc = blxWindowGetContext(blxWindow);
      GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
      
      if (dialog && GTK_WIDGET_VISIBLE(dialog))
        {
           switch (dialogId)
             {
               case BLXDIALOG_SETTINGS:
                 showSettingsDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_SORT:
                 showSortDialog(blxWindow, FALSE);
                 break;
                 
               case BLXDIALOG_HELP:
                 showHelpDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_FIND:
                 showFindDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_VIEW:
                 showViewPanesDialog(blxWindow, FALSE);
                 break;

               case BLXDIALOG_DOTTER:
                 showDotterDialog(blxWindow, FALSE);
                 break;
                 
               case BLXDIALOG_GROUPS:
                 showGroupsDialog(blxWindow, TRUE, FALSE); /* show the 'edit' pane because we've got here by adding/deleting a group */
                 break;
                 
               default:
                 break;
             };
        }
    }
  else
    {
      g_warning("Could not refresh dialog [ID=%d]; parent window not found.\n", dialogId);
    }
}


/* Callback called when the user has responded to the font selection dialog */
static void onResponseFontSelectionDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  if (responseId == GTK_RESPONSE_ACCEPT || responseId == GTK_RESPONSE_OK || responseId == GTK_RESPONSE_APPLY)
    {
      GtkWidget *blxWindow = GTK_WIDGET(data);

      /* Check that the user selected a monospace font (unfortunately there's no easy way to get the
       * font family in older GTK versions so don't bother checking) */
      gboolean ok = TRUE;
      
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (14)
      GtkFontSelection *fontSeln = GTK_FONT_SELECTION(gtk_buildable_get_internal_child(GTK_BUILDABLE(dialog), gtk_builder_new(), "font_selection"));
      PangoFontFamily *family = gtk_font_selection_get_family(fontSeln);
      ok = pango_font_family_is_monospace(family);
#endif
      
      /* Get the selected font name */
      gchar *fontName = gtk_font_selection_dialog_get_font_name(GTK_FONT_SELECTION_DIALOG(dialog));

      if (!ok)
        {
          BlxViewContext *bc = blxWindowGetContext(blxWindow);

          char *msg = g_strdup_printf("Selected font '%s' is not a fixed-width font. Matches may not appear correctly aligned. Are you sure you want to continue?", fontName);
          char *title = g_strdup_printf("%sWarning", blxGetTitlePrefix(bc));

          gint response = runConfirmationBox(GTK_WIDGET(dialog), "Blixem - Warning", msg);

          g_free(title);
          g_free(msg);

          ok = (response == GTK_RESPONSE_ACCEPT);
        }
      
      if (ok)
        {
          GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
          DetailViewProperties *properties = detailViewGetProperties(detailView);
          
          pango_font_description_free(properties->fontDesc);
          properties->fontDesc = pango_font_description_from_string(fontName);
          
          updateDetailViewFontDesc(detailView);

          g_debug("Set font family to '%s'\n", fontName);
          blxWindowRedrawAll(blxWindow);
          
          if (responseId != GTK_RESPONSE_APPLY)
            {
              gtk_widget_destroy(GTK_WIDGET(dialog));
            }
        }
    }
  else
    {
      /* Cancelled */
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/* Callback for when the font selection button is pressed. Opens the font selection dialog */
static void onFontSelectionButtonPressed(GtkWidget *button, gpointer data)
{
  GtkWidget *dialog = gtk_font_selection_dialog_new("Select fixed-width font");
  g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(onResponseFontSelectionDialog), data);
  gtk_widget_show_all(dialog);
}


/* Create a button on the settings dialog to open a font-selection dialog */
static void createFontSelectionButton(GtkBox *parent, GtkWidget *blxWindow)
{
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(parent, hbox, FALSE, FALSE, 0);

  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("Change fixed-width font:   "), FALSE, FALSE, 0);
  
  GtkWidget *button = gtk_button_new_from_stock(GTK_STOCK_SELECT_FONT);
  g_signal_connect(G_OBJECT(button), "pressed", G_CALLBACK(onFontSelectionButtonPressed), blxWindow);
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
}


/* This function restores all settings to defaults */
static void resetSettings(GtkWidget *blxWindow)
{
  gint responseId = runConfirmationBox(blxWindow, "Reset all settings", 
    "This will reset all settings to their default values: are you sure you want to continue?");
  
  if (responseId == GTK_RESPONSE_ACCEPT)
    {
      resetColumnWidths(blxWindowGetColumnList(blxWindow));
      showSettingsDialog(blxWindow, FALSE);
    }
}


void onResponseSettingsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  if (responseId == BLX_RESPONSE_RESET)
    resetSettings(dialogChildGetBlxWindow(GTK_WIDGET(dialog)));
  else
    onResponseDialog(dialog, responseId, data); /* default handler */
}


/* Show/refresh the "Settings" dialog. */
void showSettingsDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_SETTINGS;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      /* note: reset-to-defaults option commented out because it is incomplete:
       * for now, the help page tells the user to delete the ~/.blixemrc file
       * to reset to defaults */
      char *title = g_strdup_printf("%sSettings", blxGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
//                                           "Reset to defaults",
//                                           BLX_RESPONSE_RESET,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_APPLY,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      GdkScreen *screen = gtk_widget_get_screen(dialog);
      const int width = gdk_screen_get_width(screen) * 0.33;
      const int height = gdk_screen_get_height(screen) * 0.33;
      gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

      g_signal_connect(dialog, "response", G_CALLBACK(onResponseSettingsDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Need to refresh the dialog contents, so clear and re-create content area */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);

  int borderWidth = 12;
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  GtkWidget *bigPicture = blxWindowGetBigPicture(blxWindow);
  
  /* We'll put everything into a tabbed notebook */
  GtkWidget *notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, TRUE, TRUE, 0);


  /* OPTIONS PAGE */
  GtkWidget *optionsPage = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(optionsPage), gtk_label_new_with_mnemonic("Opt_ions"));

  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(optionsPage), scrollWin);

  GtkWidget *optionsBox = gtk_vbox_new(FALSE, 0);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(optionsBox));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  GtkContainer *variationContainer = createParentCheckButton(optionsBox, detailView, bc, "Highlight _variations in reference sequence", BLXFLAG_HIGHLIGHT_VARIATIONS, NULL, G_CALLBACK(onParentBtnToggled));
  createCheckButton(GTK_BOX(variationContainer), "Show variations trac_k", bc->flags[BLXFLAG_SHOW_VARIATION_TRACK], G_CALLBACK(onShowVariationTrackToggled), GINT_TO_POINTER(BLXFLAG_SHOW_VARIATION_TRACK));

  /* show-polyA-tails option and its sub-options. Connect onToggleFlag twice to the 'when selected' button to also toggle the 'show signals when selected' button in unison. */
  GtkWidget *polyAParentBtn = NULL;
  GtkContainer *polyAContainer = createParentCheckButton(optionsBox, detailView, bc, "Show polyA _tails", BLXFLAG_SHOW_POLYA_SITE, &polyAParentBtn, G_CALLBACK(onShowAdditionalSeqToggled));
  GtkWidget *polyABtn = createCheckButton(GTK_BOX(polyAContainer), "Selected sequences only", bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED], G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_POLYA_SITE_SELECTED));
  g_signal_connect(G_OBJECT(polyAParentBtn), "toggled", G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_POLYA_SIG));
  g_signal_connect(G_OBJECT(polyABtn), "toggled", G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_POLYA_SIG_SELECTED));

  const gboolean squashMatches = (bc->modelId == BLXMODEL_SQUASHED);
  
  /* show-unaligned-sequence option and its sub-options */
  GtkContainer *unalignContainer = createParentCheckButton(optionsBox, detailView, bc, "Show _unaligned sequence", BLXFLAG_SHOW_UNALIGNED, NULL, G_CALLBACK(onShowAdditionalSeqToggled));
  createLimitUnalignedBasesButton(unalignContainer, detailView, bc);
  createCheckButton(GTK_BOX(unalignContainer), "Selected sequences only", bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED], G_CALLBACK(onToggleShowUnalignedSelected), detailView);

  /* show-colinearity-lines option and its sub-options */
  GtkContainer *colinearityContainer = createParentCheckButton(optionsBox, detailView, bc, "Show _colinearity lines", BLXFLAG_SHOW_COLINEARITY, NULL, G_CALLBACK(onParentBtnToggled));
  createCheckButton(GTK_BOX(colinearityContainer), "Selected sequences only", bc->flags[BLXFLAG_SHOW_COLINEARITY_SELECTED], G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_COLINEARITY_SELECTED));

  createCheckButton(GTK_BOX(optionsBox), "Show Sp_lice Sites for selected seqs", bc->flags[BLXFLAG_SHOW_SPLICE_SITES], G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_SHOW_SPLICE_SITES));

  createCheckButton(GTK_BOX(optionsBox), "_Highlight differences", bc->flags[BLXFLAG_HIGHLIGHT_DIFFS], G_CALLBACK(onToggleFlag), GINT_TO_POINTER(BLXFLAG_HIGHLIGHT_DIFFS));
  createCheckButton(GTK_BOX(optionsBox), "_Squash matches", squashMatches, G_CALLBACK(onSquashMatches), NULL);


  /* DISPLAY PAGE */
  GtkWidget *displayPage = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(displayPage), gtk_label_new_with_mnemonic("_Display"));

  GtkWidget *displayScrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(displayPage), displayScrollWin);

  GtkWidget *displayBox = gtk_vbox_new(FALSE, 0);
  gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(displayScrollWin), GTK_WIDGET(displayBox));
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(displayScrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  GtkWidget *settingsBox = createVBoxWithBorder(displayBox, borderWidth, TRUE, "General");
  const gboolean usePrintColours = blxWindowGetUsePrintColors(blxWindow);
  createCheckButton(GTK_BOX(settingsBox), "Use _print colours", usePrintColours, G_CALLBACK(onTogglePrintColors), blxWindow);
  createFontSelectionButton(GTK_BOX(settingsBox), blxWindow);

  createGridSettingsButtons(displayBox, bigPicture);

  createCoverageSettingsButtons(displayBox, bigPicture);


  /* COLUMNS PAGE */
  GtkWidget *columnsPage = gtk_vbox_new(FALSE, 0);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(columnsPage), gtk_label_new_with_mnemonic("Colum_ns"));

  createColumnButtons(columnsPage, detailView, borderWidth);


  /* COLOURS PAGE */
  GtkWidget *appearancePage = gtk_vbox_new(FALSE, borderWidth);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), GTK_WIDGET(appearancePage), gtk_label_new_with_mnemonic("Colou_rs"));

  createColorButtons(appearancePage, blxWindow, borderWidth);

  
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/***********************************************************
 *                       Sort menu                         *
 ***********************************************************/

/* See if this widget is a combo box, or if it has a child combo box */
static GtkComboBox* widgetGetComboBox(GtkWidget *widget)
{
  GtkComboBox *result = NULL;
  
  if (GTK_IS_COMBO_BOX(widget))
    {
      result = GTK_COMBO_BOX(widget);
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      GList *children = gtk_container_get_children(GTK_CONTAINER(widget));
      GList *childItem = children;
      
      for ( ; childItem; childItem = childItem->next)
        {
          GtkWidget *childWidget = GTK_WIDGET(childItem->data);
          result = widgetGetComboBox(childWidget);
          
          if (result)
            break;
        }

      g_list_free(children);
    }
  
  return result;
}


/* For a drop-down box that contains columns, find which column is currently
 * selected */
static BlxColumnId getColumnFromComboBox(GtkComboBox *combo)
{
  BlxColumnId result = BLXCOL_NONE;
  
  if (combo)
    {
      /* Get the combo box value */
      GtkTreeIter iter;
      
      if (gtk_combo_box_get_active_iter(combo, &iter))
        {
          GtkTreeModel *model = gtk_combo_box_get_model(combo);
          
          GValue val = {0};
          gtk_tree_model_get_value(model, &iter, SORT_TYPE_COL, &val);
          
          result = (BlxColumnId)g_value_get_int(&val);
        }
    }
  
  return result;
}


/* Callback function called when the 'invert sort order' button is toggled */
static gboolean onInvertSortChanged(GtkWidget *button, const gint responseId, gpointer data)
{
  const gboolean invert = setFlagFromButton(button, data);
  
  GtkWidget *blxWindow = dialogChildGetBlxWindow(button);
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  
  detailViewUpdateSortInverted(detailView, invert);
  
  return TRUE;
}


/* Callback called when the sort order has been changed in the drop-down box */
static gboolean onSortOrderChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *detailView = GTK_WIDGET(data);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);

  if (GTK_WIDGET_REALIZED(detailView) && GTK_IS_CONTAINER(widget))
    {
      /* Loop through each child of the given widget (assumes that each child is or
       * contains one combo box) */
      GList *children = gtk_container_get_children(GTK_CONTAINER(widget));
      GList *childItem = children;
      const int numColumns = g_list_length(columnList);
      int priority = 0;
      
      for ( ; childItem; childItem = childItem->next, ++priority)
        {
          if (priority >= numColumns)
            {
              g_critical("Exceeded max number of sort columns (%d).\n", numColumns);
              break;
            }
          
          /* See if this is a (or has a child) combo box */
          GtkWidget *childWidget = GTK_WIDGET(childItem->data);
          GtkComboBox *combo = widgetGetComboBox(childWidget);
          
          if (combo)
            {
              dvProperties->sortColumns[priority] = getColumnFromComboBox(combo);
            }
        }
      
      g_list_free(children);

      /* Re-sort trees */
      detailViewResortTrees(detailView);
    }
  
  return TRUE;
}


/* Add an option for the sorting drop-down box */
static GtkTreeIter* addSortBoxItem(GtkTreeStore *store, 
                                   GtkTreeIter *parent, 
                                   BlxColumnId sortColumn, 
                                   const char *sortName,
                                   BlxColumnId initSortColumn,
                                   GtkComboBox *combo)
{
  GtkTreeIter iter;
  gtk_tree_store_append(store, &iter, parent);
  
  gtk_tree_store_set(store, &iter, SORT_TYPE_COL, sortColumn, SORT_TEXT_COL, sortName, -1);
  
  if (sortColumn == initSortColumn)
    {
      gtk_combo_box_set_active_iter(combo, &iter);
    }
  
  return NULL;
}


/* Create the combo box used for selecting sort criteria */
static void createSortBox(GtkBox *parent, 
                          GtkWidget *detailView, 
                          const BlxColumnId initSortColumn, 
                          GList *columnList, 
                          const char *labelText,
                          const gboolean searchableOnly)
{
  /* Put the label and drop-down in a box */
  GtkWidget *box = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(parent, box, FALSE, FALSE, 0);
  
  /* Add a label, to make it obvious what the combo box is for */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
  gtk_container_add(GTK_CONTAINER(box), label);
  
  /* Create the data for the drop-down box. Use a tree so that we can sort by
   * multiple criteria. */
  GtkTreeStore *store = gtk_tree_store_new(N_SORT_COLUMNS, G_TYPE_INT, G_TYPE_STRING);
  GtkComboBox *combo = GTK_COMBO_BOX(gtk_combo_box_new_with_model(GTK_TREE_MODEL(store)));
  g_object_unref(store);
  gtk_container_add(GTK_CONTAINER(box), GTK_WIDGET(combo));

  /* Create a cell renderer to display the sort text. */
  GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
  gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, FALSE);
  gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", SORT_TEXT_COL, NULL);
  
  GtkTreeIter *iter = NULL;
  
  /* Add a blank row for the case where nothing is selected (unless we only
   * want searchable columns, because we can't search the NONE column) */
  if (!searchableOnly)
    iter = addSortBoxItem(store, iter , BLXCOL_NONE, "<select column>", initSortColumn, combo);

  /* Add a row for each column that has the 'sortName' property set. */
  GList *columnItem = columnList;

  for ( ; columnItem; columnItem = columnItem->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(columnItem->data);

      /* Only include columns that have a sort name and, if searchableOnly is
       * true, only include columns that are searchable. */
      if (columnInfo->sortName && (columnInfo->searchable || !searchableOnly))
        {
          iter = addSortBoxItem(store, iter, columnInfo->columnId, columnInfo->sortName, initSortColumn, combo);
        }
    }
}


/* Called when the user clicks the 'add new sort-by box' button on the sort dialog */
static void onAddNewSortByBox(GtkButton *button, gpointer data)
{
  /* The user-data is the container box that contains the sort-by boxes */
  GtkBox *box = GTK_BOX(data);

  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(GTK_WIDGET(button)));
  GtkWidget *blxWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);

  /* Add another sort-by box to the container */
  createSortBox(box, detailView, BLXCOL_NONE, columnList, "then by", FALSE);
  
  gtk_widget_show_all(GTK_WIDGET(box));
}


/* Show/refresh the "Sort" dialog. */
void showSortDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  const BlxDialogId dialogId = BLXDIALOG_SORT;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sSort", blxGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(blxWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_APPLY,
                                           GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Need to refresh the dialog contents, so clear and re-create content area */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);
  
  const int borderWidth = 12;
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  DetailViewProperties *dvProperties = detailViewGetProperties(detailView);
  GList *columnList = blxWindowGetColumnList(dvProperties->blxWindow);
  const int numColumns = g_list_length(columnList);

  /* Add a drop-down for each sort column that is currently specified (or just
   * the default number if none specified). */
  GtkWidget *vbox = gtk_vbox_new(FALSE, borderWidth);
  gtk_container_add(GTK_CONTAINER(contentArea), vbox);
  int sortPriority = 0;
  const int minBoxes = 3;
  
  for ( ; sortPriority < numColumns; ++sortPriority)
    {
      const BlxColumnId columnId = dvProperties->sortColumns[sortPriority];
      
      if (columnId != BLXCOL_NONE || sortPriority < minBoxes)
        {
          if (sortPriority == 0)
            createSortBox(GTK_BOX(vbox), detailView, columnId, columnList, "Sort by", FALSE);
          else
            createSortBox(GTK_BOX(vbox), detailView, columnId, columnList, "then by", FALSE);
        }
      else
        {
          break;
        }
    }
  
  widgetSetCallbackData(vbox, onSortOrderChanged, detailView);

  
  /* Add a button to add a new sort-by box */
  GtkWidget *button = gtk_button_new_from_stock(GTK_STOCK_ADD);
  gtk_box_pack_end(GTK_BOX(vbox), button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(onAddNewSortByBox), vbox);
  
  
  /* Add a toggle button for the 'invert sort order' option */
  GtkWidget *toggle = gtk_check_button_new_with_mnemonic("_Invert sort order");
  gtk_box_pack_end(GTK_BOX(vbox), toggle, FALSE, FALSE, 0);
  
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle), bc->flags[BLXFLAG_INVERT_SORT]);
  widgetSetCallbackData(toggle, onInvertSortChanged, GINT_TO_POINTER(BLXFLAG_INVERT_SORT));
  

  /* Shot the dialog */
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
}


/***********************************************************
 *                  Statistics menu                        *
 ***********************************************************/

static void getStats(GtkWidget *blxWindow, GString *result, MSP *MSPlist)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Compile data for sequences */
  int totalNumSeqs = 0;             /* total number of sequences */
  int numValidSeqs = 0;             /* how many sequences have sequence data filled in */

  gint32 seqDataSize = 0;           /* total memory used by the sequence data */
  gint32 seqStructSize = 0;         /* total memory used by the sequence structs */

  GList *seqItem = bc->matchSeqs;
  for ( ; seqItem; seqItem = seqItem->next)
    {
      ++totalNumSeqs;
      seqStructSize += sizeof(BlxSequence);
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      const char *sequence = blxSequenceGetSequence(blxSeq);

      if (sequence)
        {
          ++numValidSeqs;
          seqDataSize += strlen(sequence) * sizeof(char);
        }
    }
  

  /* Compile data for MSPs */
  int numMSPs = 0;                /* total number of MSPs */
  gint32 mspStructSize = 0;       /* total memory used by the MSP structs */

  MSP *msp = NULL;
  for (msp = MSPlist; msp ; msp = msp->next)
    {
      ++numMSPs;
      mspStructSize += sizeof(MSP);
    }

  
  /* Other data */
  int refSeqLen = strlen(blxWindowGetRefSeq(blxWindow));

  
  /* Create the text based on the results */
  g_string_printf(result, "%s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s %s%d%s",
                  "Length of reference sequence\t\t\t\t\t\t\t= ", refSeqLen, " characters\n\n",
                  "Total number of match sequences\t\t\t\t\t\t= ", totalNumSeqs, "\n",
                  "Number of match sequences containing sequence data\t= ", numValidSeqs, "\n",
                  "Total memory used by sequence data\t\t\t\t\t= ", seqDataSize, " bytes\n\n",
                  "Size of each sequence struct\t\t\t\t\t\t\t= ", (int)sizeof(BlxSequence), " bytes\n",
                  "Total memory used by sequence structs\t\t\t\t\t= ", seqStructSize, " bytes\n\n",
                  "Number of MSPs\t\t\t\t\t\t\t\t\t\t= ", numMSPs, "\n",
                  "Size of each MSP\t\t\t\t\t\t\t\t\t\t= ", (int)sizeof(MSP), " bytes\n",
                  "Total memory used by MSP structs\t\t\t\t\t\t= ", (int)sizeof(MSP) * numMSPs, " bytes");
}


static void showStatsDialog(GtkWidget *blxWindow, MSP *MSPlist)
{
  /* Create a dialog widget with an OK button */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  char *title = g_strdup_printf("%sStatistics", blxGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(blxWindow),
                                                  GTK_DIALOG_DESTROY_WITH_PARENT,
                                                  GTK_STOCK_OK,
                                                  GTK_RESPONSE_ACCEPT,
                                                  NULL);

  g_free(title);
  
  /* Ensure that the dialog box (along with any children) is destroyed when the user responds. */
  g_signal_connect (dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  /* Create a text buffer containing the required text*/
  GString *displayText = g_string_sized_new(200); /* will be extended if we need more space */
  getStats(blxWindow, displayText, MSPlist);
  GtkTextBuffer *textBuffer = gtk_text_buffer_new(gtk_text_tag_table_new());
  
  gtk_text_buffer_set_text(GTK_TEXT_BUFFER(textBuffer), displayText->str, -1);
  
  g_string_free(displayText, TRUE);
  
  /* Create a text view widget and put it in the vbox area of the dialog */
  GtkWidget *textView = gtk_text_view_new_with_buffer(GTK_TEXT_BUFFER(textBuffer));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), textView, TRUE, TRUE, 0);
  
  /* Show the dialog */
  gtk_widget_show(textView);
  gtk_widget_show(dialog);
}


/***********************************************************
 *                      About dialog                       *
 ***********************************************************/

/* A GtkAboutDialogActivateLinkFunc() called when user clicks on website link in "About" window. */
static void aboutDialogOpenLinkCB(GtkAboutDialog *about, const gchar *link, gpointer data)
{
  GError *error = NULL ;
    
  if (!seqtoolsLaunchWebBrowser(link, &error))
    g_critical("Cannot show link in web browser: \"%s\"", link) ;    
}


/* Shows the 'About' dialog */
void showAboutDialog(GtkWidget *parent)
{
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (6)
  const gchar *authors[] = {AUTHOR_LIST, NULL} ;

  gtk_about_dialog_set_url_hook(aboutDialogOpenLinkCB, NULL, NULL) ;
  
  gtk_show_about_dialog(GTK_WINDOW(parent),
                        "authors", authors,
                        "comments", blxGetCommentsString(), 
                        "copyright", blxGetCopyrightString(),
                        "license", blxGetLicenseString(),
                        "name", blxGetAppName(),
                        "version", blxGetVersionString(),
                        "website", blxGetWebSiteString(),
                        NULL) ;
#endif
}


/***********************************************************
 *                      Help menu                          *
 ***********************************************************/

void onResponseHelpDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      destroy = TRUE;
      break;
      
    case GTK_RESPONSE_HELP:
      showAboutDialog(NULL);
      destroy = FALSE;
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
      /* If it's a persistent dialog, just hide it, otherwise destroy it */
      const gboolean isPersistent = GPOINTER_TO_INT(data);
      
      if (isPersistent)
        {
          gtk_widget_hide_all(GTK_WIDGET(dialog));
        }
      else
        {
          gtk_widget_destroy(GTK_WIDGET(dialog));
        }
    }
}


void showHelpDialog(GtkWidget *blxWindow, const gboolean bringToFront)
{
  GError *error = NULL;

  /* The docs should live in /share/doc/seqtools/, in the same parent
   * directory that our executable's 'bin' directory is in. Open the 'quick
   * start' page. */
  char rel_path[100] = "../share/doc/seqtools/blixem_quick_start.html";

  /* Find the executable's path */
  char *exe = g_find_program_in_path(g_get_prgname());
  gboolean ok = (exe != NULL);

  if (ok)
    {
      /* Get the executable's directory */
      char *dir = g_path_get_dirname(exe);
      
      ok = dir != NULL;
      
      if (ok)
        {
          /* Get the path to the html page */
          char *path = g_strdup_printf("%s/%s", dir, rel_path);
          
          ok = path != NULL;
          
          if (ok)
            {
              g_message("Opening help page '%s'\n", path);
              seqtoolsLaunchWebBrowser(path, &error);
              g_free(path);
            }

          g_free(dir);
        }

      g_free(exe);
    }
  
  if (!ok)
    {
      if (error)
        reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      else
        g_critical("Could not find help documentation: %s\n", rel_path);
    }
}

/***********************************************************
 *                        Menu actions                     *
 ***********************************************************/

/* Called when the user selects the quit menu option, or hits the Quit shortcut key.
 * Pops up a dialog showing user help information */
static void onQuit(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  gtk_widget_destroy(blxWindow);
}

static void onHelpMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showHelpDialog(blxWindow, TRUE);
}

static void onAboutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showAboutDialog(blxWindow);
}

static void onFindMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showFindDialog(blxWindow, TRUE);
}

static void onGoToMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  
  /* We currently only accept input in terms of DNA coords on the ref seq */
  const BlxSeqType seqType = BLXSEQ_DNA;
  
  goToDetailViewCoord(blxWindowGetDetailView(blxWindow), seqType);
}

static void onPrevMatchMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(blxWindow);
  prevMatch(blxWindowGetDetailView(blxWindow), seqList);
}

static void onNextMatchMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(blxWindow);
  nextMatch(blxWindowGetDetailView(blxWindow), seqList);
}

static void onFirstMatchMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(blxWindow);
  firstMatch(blxWindowGetDetailView(blxWindow), seqList);
}

static void onLastMatchMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GList *seqList = blxWindowGetSelectedSeqs(blxWindow);
  lastMatch(blxWindowGetDetailView(blxWindow), seqList);
}

static void onPageLeftMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  scrollDetailViewLeftPage(blxWindowGetDetailView(blxWindow));
}

static void onPageRightMenu(GtkAction *action, gpointer data)
{  
  GtkWidget *blxWindow = GTK_WIDGET(data);
  scrollDetailViewRightPage(blxWindowGetDetailView(blxWindow));
}

static void onScrollLeft1Menu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  scrollDetailViewLeft1(blxWindowGetDetailView(blxWindow));
}

static void onScrollRight1Menu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  scrollDetailViewRight1(blxWindowGetDetailView(blxWindow));
}

static void onSquashMatchesMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  
  const gboolean squash = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
  
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  detailViewUpdateSquashMatches(detailView, squash);
  
  /* Refresh the settings dialog, if it is open */
  refreshDialog(BLXDIALOG_SETTINGS, blxWindow);
}

static void onToggleStrandMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  toggleStrand(blxWindowGetDetailView(blxWindow));
}

/* Called when the user selects the View menu option, or hits the Settings shortcut key */
static void onViewMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showViewPanesDialog(blxWindow, TRUE);
}


/* Called when the user selects the 'Group Sequences' menu option, or hits the relevant shortcut key */
static void onCreateGroupMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, FALSE, TRUE);
}


/* Called when the user selects the 'Groups' menu option, or hits the relevant shortcut key */
static void onEditGroupsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showGroupsDialog(blxWindow, TRUE, TRUE);
}

/* Called when the user selects the 'Toggle match set' option, or hits the relevant shortcut key */
static void onToggleMatchSet(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  toggleMatchSet(blxWindow);
}

/* Called when the user selects the Settings menu option, or hits the Settings shortcut key */
static void onSettingsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showSettingsDialog(blxWindow, TRUE);
}


static void onLoadMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  GError *tmp_error = NULL;

  char *filename = getLoadFileName(blxWindow, NULL, "Load file");
  dynamicLoadFeaturesFile(blxWindow, filename, NULL, &tmp_error);
  
  g_free(filename);

  reportAndClearIfError(&tmp_error, G_LOG_LEVEL_CRITICAL);
}

static void onCopySeqsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  copySelectionToClipboard(blxWindow);
}

static void onCopySeqDataMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  copySelectedSeqDataToClipboard(blxWindow);
}

static void onCopyRefSeqMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);

  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);

  /* Copy the portion of the ref seq from the selected index
   * to the clicked index */
  
  if (detailViewGetSelectedBaseSet(detailView))
    {
      const int fromIdx = detailViewGetSelectedDnaBaseIdx(detailView);
      const int toIdx = detailViewGetClickedBaseIdx(detailView);

      copyRefSeqToClipboard(blxWindow, fromIdx, toIdx);
    }
  else
    {
      g_critical("Please middle-click on a coordinate first to set the mark\n");
    }
}

static void onSortMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showSortDialog(blxWindow, TRUE);
}

static void onZoomInMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  zoomDetailView(blxWindowGetDetailView(blxWindow), TRUE);
}

static void onZoomOutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  zoomDetailView(blxWindowGetDetailView(blxWindow), FALSE);
}

/* Called when the user selects the Dotter menu option, or hits the Dotter shortcut key */
static void onDotterMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  showDotterDialog(blxWindow, TRUE);
}

/* Called when the user selects the 'Close all Dotters' menu option */
static void onCloseAllDottersMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  /* Check if there are actually any spawned processes */
  if (g_slist_length(bc->spawnedProcesses) > 0)
    {
      gint responseId = runConfirmationBox(blxWindow, "Close all Dotters", 
        "Are you sure you want to close all Dotters started from this Blixem?");

      if (responseId == GTK_RESPONSE_ACCEPT)
        {
          killAllSpawned(bc);
        }
    }
  else
    {
      g_message("No Dotters to close.\n");
    }
}

/* Called when the user selects the 'Select Features' menu option, or hits the relevant shortcut key */
static void onSelectFeaturesMenu(GtkAction *action, gpointer data)
{
//  selectFeatures();
}


/* Called when the user selects the 'Deselect all' menu option, or hits the relevant shortcut key */
static void onDeselectAllRows(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  blxWindowDeselectAllSeqs(blxWindow);
}


/* Called when the user selects the Statistics menu option, or hits the Statistics shortcut key.
 * Pops up a dialog showing memory usage statistics. */
static void onStatisticsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  MSP *mspList = blxWindowGetMspList(blxWindow);
  showStatsDialog(blxWindow, mspList);
}


/* Called when the user selects the Print menu option, or hits the Print shortcut key */
static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);

  /* We need to do some work to prepare the big picture for printing */
  bigPicturePrepareForPrinting(properties->bigPicture);
  
  blxPrintWidget(blxWindow, NULL, GTK_WINDOW(blxWindow), &properties->printSettings, &properties->pageSetup, NULL, TRUE, PRINT_FIT_BOTH);
  
  bigPictureRedrawAll(properties->bigPicture);
}


/* Called when the user selects the Page Setup menu option, or hits the Print shortcut key */
static void onPageSetupMenu(GtkAction *action, gpointer data)
{
  GtkWidget *blxWindow = GTK_WIDGET(data);
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);

  if (!properties->pageSetup)
    properties->pageSetup = gtk_page_setup_new();
  
  if (!properties->printSettings)
    properties->printSettings = gtk_print_settings_new();
  
  properties->pageSetup = gtk_print_run_page_setup_dialog(GTK_WINDOW(blxWindow), 
                                                          properties->pageSetup, 
                                                          properties->printSettings);
}  

/***********************************************************
 *                         Events                          *
 ***********************************************************/

/* Mouse button handler */
static gboolean onButtonPressBlxWindow(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  if (event->type == GDK_BUTTON_PRESS && event->button == 3)
    {
      gtk_menu_popup (GTK_MENU (data), NULL, NULL, NULL, NULL, event->button, event->time);
      return TRUE;
  }
  
  return TRUE;
}


/* Mouse button handler for the paned window containing the big picture and detail view */
static gboolean onButtonPressPanedWin(GtkWidget *panedWin, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  switch (event->button)
    {
    case 1: /* left button */
      {
        if (event->type == GDK_2BUTTON_PRESS) /* double-click */
          {
            /* When the user double-clicks the paned window separator, reset the splitter position
             * (i.e. so that gets automatically positioned based on the child widgets' size)
             * to do: This makes the splitter jump temporarily to the desired position but then it
             * immediately jumps back, so I'm leaving it out for now */
            /* gtk_paned_set_position(GTK_PANED(panedWin), 100); */
            handled = TRUE;
          }

        break;
      }
      
    default:
      break;
    };
    
  return handled;
}


/* Handlers for specific key presses */
static gboolean onKeyPressEscape(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  /* Reset the selected base index. Leave the selected frame as it is, though. */
  GtkWidget *detailView = blxWindowGetDetailView(window);
  detailViewUnsetSelectedBaseIdx(detailView);
  return TRUE;
}

static gboolean onKeyPressLeftRight(GtkWidget *window, const gboolean left, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  if (ctrlModifier)
    {
      goToMatch(window, left);
    }
  else if (shiftModifier)
    {
      moveSelectedBaseIdxBy1(window, left);
    }
  else
    {
      moveSelectedDisplayIdxBy1(window, left);
    }
  
  return TRUE;
}

static gboolean onKeyPressUpDown(GtkWidget *window, const gboolean up, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = moveRowSelection(window, up, ctrlModifier, shiftModifier);
  return result;
}

static gboolean onKeyPressHomeEnd(GtkWidget *window, const gboolean home, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollToExtremity(window, home, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressComma(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollDetailView(window, TRUE, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressPeriod(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  scrollDetailView(window, FALSE, ctrlModifier);
  return TRUE;
}

static gboolean onKeyPressPlusMinus(GtkWidget *window, const gboolean zoomIn, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  zoomBlxWindow(window, zoomIn, ctrlModifier, shiftModifier);
  return TRUE;
}

static gboolean onKeyPressV(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (ctrlModifier)
    {
      /* Paste from the default clipboard */
      requestDefaultClipboardText(findSeqsFromClipboard, window);
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressF(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  if (ctrlModifier)
    {
      showFindDialog(window, TRUE);
    }
  else
    {
      requestPrimaryClipboardText(findSeqsFromClipboard, window);
    }

  return TRUE;
}

static gboolean onKeyPressP(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (!ctrlModifier)
    {
      goToDetailViewCoord(blxWindowGetDetailView(window), BLXSEQ_DNA); /* for now, only accept input in terms of DNA seq coords */
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressG(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  gboolean result = FALSE;
  
  if (!ctrlModifier)
    {       
      toggleMatchSet(window);
      result = TRUE;
    }
  
  return result;
}

static gboolean onKeyPressB(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  toggleBumpState(window);
  
  /* Refresh the view dialog, if it is open */
  refreshDialog(BLXDIALOG_VIEW, window);
  
  return TRUE;
}

static gboolean onKeyPressT(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  toggleStrand(blxWindowGetDetailView(window));
  return TRUE;
}

static gboolean onKeyPressI(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  showInfoDialog(window);
  return TRUE;
}

static gboolean onKeyPressNumber(GtkWidget *window, const int number, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  togglePaneVisibility(window, number, ctrlModifier, shiftModifier);
  
  /* Refresh the view dialog, if it is open */
  refreshDialog(BLXDIALOG_VIEW, window);
  
  return TRUE;
}

static gboolean onKeyPressF3(GtkWidget *window, const gboolean ctrlModifier, const gboolean shiftModifier)
{
  findAgain(window, shiftModifier);
  return TRUE;
}

/* Key press handler */
static gboolean onKeyPressBlxWindow(GtkWidget *window, GdkEventKey *event, gpointer data)
{
  gboolean result = FALSE;
  
  const gboolean ctrlModifier = (event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK;
  const gboolean shiftModifier = (event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK;
  
  switch (event->keyval)
    {
      case GDK_Escape:      result = onKeyPressEscape(window, ctrlModifier, shiftModifier);           break;
      
      case GDK_Left:        result = onKeyPressLeftRight(window, TRUE, ctrlModifier, shiftModifier);  break;
      case GDK_Right:       result = onKeyPressLeftRight(window, FALSE, ctrlModifier, shiftModifier); break;
      
      case GDK_Up:          result = onKeyPressUpDown(window, TRUE, ctrlModifier, shiftModifier);     break;
      case GDK_Down:        result = onKeyPressUpDown(window, FALSE, ctrlModifier, shiftModifier);    break;
        
      case GDK_Home:        result = onKeyPressHomeEnd(window, TRUE, ctrlModifier, shiftModifier);    break;
      case GDK_End:         result = onKeyPressHomeEnd(window, FALSE, ctrlModifier, shiftModifier);   break;

      case GDK_comma:       result = onKeyPressComma(window, ctrlModifier, shiftModifier);            break;
      case GDK_period:      result = onKeyPressPeriod(window, ctrlModifier, shiftModifier);           break;

      case GDK_equal:       /* fall through */
      case GDK_plus:        result = onKeyPressPlusMinus(window, TRUE, ctrlModifier, shiftModifier);  break;
      
      case GDK_minus:       /* fall through */
      case GDK_underscore:  result = onKeyPressPlusMinus(window, FALSE, ctrlModifier, shiftModifier); break;
      
      case GDK_F3:          result = onKeyPressF3(window, ctrlModifier, shiftModifier);               break;

      case GDK_v:           /* fall through */
      case GDK_V:           result = onKeyPressV(window, ctrlModifier, shiftModifier);                break;
        
      case GDK_f:           /* fall through */
      case GDK_F:           result = onKeyPressF(window, ctrlModifier, shiftModifier);                break;

      case GDK_g:           /* fall through */
      case GDK_G:           result = onKeyPressG(window, ctrlModifier, shiftModifier);                break;

      case GDK_p:           /* fall through */
      case GDK_P:           result = onKeyPressP(window, ctrlModifier, shiftModifier);                break;
                
      case GDK_b:           /* fall through */
      case GDK_B:           result = onKeyPressB(window, ctrlModifier, shiftModifier);                break;
        
      case GDK_t:           /* fall through */
      case GDK_T:           result = onKeyPressT(window, ctrlModifier, shiftModifier);                break;

      case GDK_i:           /* fall through */
      case GDK_I:           result = onKeyPressI(window, ctrlModifier, shiftModifier);                break;

      case GDK_1:           /* fall through */
      case GDK_exclam:      result = onKeyPressNumber(window, 1, ctrlModifier, shiftModifier);        break;

      case GDK_2:           /* fall through */
      case GDK_quotedbl:    /* fall through */
      case GDK_at:          result = onKeyPressNumber(window, 2, ctrlModifier, shiftModifier);        break;

      case GDK_3:           /* fall through */
      case GDK_currency:    result = onKeyPressNumber(window, 3, ctrlModifier, shiftModifier);        break;
    };
  
  return result;
}

/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BlxWindowProperties* blxWindowGetProperties(GtkWidget *widget)
{
  /* optimisation: cache result, because we know there is only ever one main window */
  static BlxWindowProperties *properties = NULL;
  
  if (!properties && widget)
    properties = (BlxWindowProperties*)(g_object_get_data(G_OBJECT(widget), "BlxWindowProperties"));
  
  return properties;
}

BlxViewContext* blxWindowGetContext(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->blxContext : NULL;
}

GList* blxWindowGetColumnList(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  return (bc ? bc->columnList : NULL);
}

/* Kill all processes spawned from blixem */
static void killAllSpawned(BlxViewContext *bc)
{
  GSList *processes = bc->spawnedProcesses;
  
  for ( ; processes; processes = processes->next)
    {
      pid_t pid = GPOINTER_TO_INT(processes->data);
      kill(pid, 9);
    }
    
  if (bc->spawnedProcesses)
    {
      g_slist_free(bc->spawnedProcesses);
      bc->spawnedProcesses = NULL;
    }
}


/* utility to free the given pointer and set it to null */
static void freeAndNull(gpointer *ptr)
{
  if (ptr && *ptr)
    {
      g_free(*ptr);
      *ptr = NULL;
    }
}


static void destroyBlxContext(BlxViewContext **bcPtr)
{
  if (bcPtr && *bcPtr)
    {
      BlxViewContext *bc = *bcPtr;

      /* Free allocated strings */
      freeAndNull((gpointer*)(&bc->dataset));
      freeAndNull((gpointer*)(&bc->refSeqName));

      /* Free table of fetch methods and the fetch-method structs */
      /* to do */
      
      /* Free the list of selected sequence names (not the names themselves
       * because we don't own them). */
      if (bc->selectedSeqs)
        {
          g_list_free(bc->selectedSeqs);
          bc->selectedSeqs = NULL;
        }
      
      blxContextDeleteAllSequenceGroups(bc);

      /* Free the color array */
      if (bc->defaultColors)
        {
          int i = BLXCOLOR_MIN + 1;
          for (; i < BLXCOL_NUM_COLORS; ++i)
            {
              BlxColor *blxColor = &g_array_index(bc->defaultColors, BlxColor, i);
              destroyBlxColor(blxColor);
            }

          g_array_free(bc->defaultColors, TRUE);
          bc->defaultColors = NULL;
        }

      /* destroy the feature lists. note that the stored msps are owned
      * by the msplist, not by the feature lists */
      int typeId = 0;
      for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
        g_array_free(bc->featureLists[typeId], FALSE);
      
      destroyMspList(&(bc->mspList));
      destroyBlxSequenceList(&(bc->matchSeqs));
      blxDestroyGffTypeList(&(bc->supportedTypes));
      killAllSpawned(bc);
      
      /* Free the context struct itself */
      g_free(bc);
      *bcPtr = NULL;
    }
}


/* This returns the name for the given flag enum. It returns null if the name
 * has not been set. This is used to save settings in the config file; only 
 * flags whose name is set in this function will be saved. */
static const char* getFlagName(const BlxFlag flag)
{
  /* Create an array of names for all of the relevant settings */
  static const char* names[BLXFLAG_NUM_FLAGS] = {0, 0};
  
  if (!names[1]) /* only populate it once */
    {
      names[BLXFLAG_HIGHLIGHT_DIFFS] = SETTING_NAME_HIGHLIGHT_DIFFS;
      names[BLXFLAG_HIGHLIGHT_VARIATIONS] = SETTING_NAME_HIGHLIGHT_VARIATIONS;
      names[BLXFLAG_SHOW_VARIATION_TRACK] = SETTING_NAME_SHOW_VARIATION_TRACK;
      names[BLXFLAG_SHOW_UNALIGNED] = SETTING_NAME_SHOW_UNALIGNED;
      names[BLXFLAG_SHOW_UNALIGNED_SELECTED] = SETTING_NAME_SHOW_UNALIGNED_SELECTED;
      names[BLXFLAG_LIMIT_UNALIGNED_BASES] = SETTING_NAME_LIMIT_UNALIGNED_BASES;
      names[BLXFLAG_SHOW_POLYA_SITE] = SETTING_NAME_SHOW_POLYA_SITE;
      names[BLXFLAG_SHOW_POLYA_SITE_SELECTED] = SETTING_NAME_SHOW_POLYA_SITE_SELECTED;
      names[BLXFLAG_SHOW_POLYA_SIG] = SETTING_NAME_SHOW_POLYA_SIG;
      names[BLXFLAG_SHOW_POLYA_SIG_SELECTED] = SETTING_NAME_SHOW_POLYA_SIG_SELECTED;
      names[BLXFLAG_SHOW_SPLICE_SITES] = SETTING_NAME_SHOW_SPLICE_SITES;
      names[BLXFLAG_SHOW_COLINEARITY] = SETTING_NAME_SHOW_COLINEARITY;
      names[BLXFLAG_SHOW_COLINEARITY_SELECTED] = SETTING_NAME_SHOW_COLINEARITY_SELECTED;
    }

  if (flag <= BLXFLAG_MIN || flag >= BLXFLAG_NUM_FLAGS)
    return NULL;
  else
    return names[flag];
}


/* Called by saveBlixemSettings; does the work to save the boolean flags */
static void saveBlixemSettingsFlags(BlxViewContext *bc, GKeyFile *key_file)
{
  /* loop through each (save-able) setting */
  int flag = BLXFLAG_MIN + 1;
  
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      const char *flagName = getFlagName((BlxFlag)flag);

      if (flagName)
        g_key_file_set_integer(key_file, SETTINGS_GROUP, flagName, bc->flags[flag]);
    }
}


/* This function saves all of blixem's customisable settings to
 * a config file. */
static void saveBlixemSettings(GtkWidget *blxWindow)
{
  char *filename = g_strdup_printf("%s/%s", g_get_home_dir(), BLIXEM_SETTINGS_FILE);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  GKeyFile *key_file = g_key_file_new();
  GKeyFileFlags flags = G_KEY_FILE_NONE;
  GError *error = NULL;
  
  /* Load existing contents so they can be merged, if the file already exists */
  g_key_file_load_from_file(key_file, filename, flags, &error);
  g_message("Saving Blixem settings to '%s'.\n", filename);

  /* Write the settings */
  saveBlixemSettingsFlags(blxWindowGetContext(blxWindow), key_file);
  g_key_file_set_integer(key_file, SETTINGS_GROUP, SETTING_NAME_SQUASH_MATCHES, (bc->modelId == BLXMODEL_SQUASHED));
  detailViewSaveProperties(blxWindowGetDetailView(blxWindow), key_file);

  /* Output to the config file */
  gchar *file_content = g_key_file_to_data(key_file, NULL, NULL);
      
  if (!g_file_set_contents(filename, file_content, -1, NULL))
    g_warning("Error saving settings to '%s'.\n", filename);

  g_free(file_content);
  g_key_file_free(key_file);
}


static void onDestroyBlxWindow(GtkWidget *widget)
{
  BlxWindowProperties *properties = blxWindowGetProperties(widget);
  
  saveBlixemSettings(widget);

  if (properties)
    {
      destroyBlxContext(&properties->blxContext);

      if (properties->mainmenu)
        {
          gtk_widget_destroy(properties->mainmenu);
          properties->mainmenu = NULL;
        }
      
      /* Destroy the print settings */
      if (properties->printSettings)
        {
          g_object_unref(properties->printSettings);
          properties->printSettings = NULL;
        }
      
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(widget), "BlxWindowProperties", NULL);
    }

  /* Reset any globals */
  blviewResetGlobals();
  
  gtk_main_quit();
}


/* Create the colors that blixem will use for various specific purposes */
static void createBlxColors(BlxViewContext *bc, GtkWidget *widget)
{
  /* Initialise the array with empty BlxColor structs */
  bc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), BLXCOL_NUM_COLORS);
  int i = BLXCOLOR_MIN + 1;
  
  for ( ; i < BLXCOL_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = (BlxColor*)g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(bc->defaultColors, *blxColor);
    }
  
  /* Get the default background color of our widgets (i.e. that inherited from the theme).
   * Convert it to a string so we can use the same creation function as the other colors */
  char *defaultBgColorStr = convertColorToString(&widget->style->bg[GTK_STATE_NORMAL]);
  createBlxColor(bc->defaultColors, BLXCOLOR_BACKGROUND, "Background", "Background color", defaultBgColorStr, BLX_WHITE, "#bdbdbd", NULL);
  
  /* reference sequence */
  createBlxColor(bc->defaultColors, BLXCOLOR_REF_SEQ, "Reference sequence", "Default background color for the reference sequence", BLX_YELLOW, BLX_VERY_LIGHT_GREY, BLX_DARK_YELLOW, NULL);
  
  /* matches */
  createBlxColor(bc->defaultColors, BLXCOLOR_MATCH, "Exact match", "Exact match", BLX_LIGHT_CYAN, BLX_LIGHT_GREY, BLX_CYAN, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_CONS, "Conserved match", "Conserved match", BLX_VIOLET, BLX_VERY_LIGHT_GREY, BLX_DARK_VIOLET, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_MISMATCH, "Mismatch", "Mismatch", "#FFFFFF", BLX_WHITE, "#FED4EA", NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_INSERTION, "Insertion", "Insertion", "#9E00FF", BLX_VERY_DARK_GREY, NULL, NULL);
  
  /* exons */
  createBlxColor(bc->defaultColors, BLXCOLOR_EXON_START, "Exon start", "Exon start boundary", BLX_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_EXON_END, "Exon end", "Exon end boundary", BLX_DARK_BLUE, BLX_GREY, NULL, NULL);

  createBlxColor(bc->defaultColors, BLXCOLOR_EXON_FILL, "Exon fill color", "Exon fill color in big picture", BLX_YELLOW, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_EXON_LINE, "Exon line color", "Exon line color in big picture", BLX_BLUE, BLX_VERY_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_CDS_FILL, "CDS fill color", "Coding section fill color in big picture", BLX_LIGHT_GREEN, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_CDS_LINE, "CDS line color", "Coding section line color in big picture", BLX_DARK_GREEN, BLX_DARK_GREY, BLX_VERY_DARK_GREEN, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_UTR_FILL, "Exon fill color (UTR)", "Untranslated region fill color in big picture", BLX_LIGHT_RED, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_UTR_LINE, "Exon line color (UTR)", "Untranslated region line color in big picture", BLX_DARK_RED, BLX_VERY_DARK_GREY, BLX_VERY_DARK_RED, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_PARTIAL_EXON_CROSSHATCH, "Cross-hatch line color for partial exons", "Line color of cross-hatch highlighting for partial exons", BLX_GREY, BLX_GREY, NULL, NULL);
  
  /* codons */
  createBlxColor(bc->defaultColors, BLXCOLOR_CODON, "Codon nucleotides", "Codon nucleotides", BLX_SKY_BLUE, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_MET, "MET codons", "MET codons", BLX_LAWN_GREEN, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_STOP, "STOP codons", "MET codons", BLX_SALMON_RED, BLX_LIGHT_GREY, NULL, NULL);
  
  /* SNPs */
  createBlxColor(bc->defaultColors, BLXCOLOR_SNP, "SNPs", "SNPs", BLX_ORANGE, BLX_GREY, NULL, NULL);

  /* Big Picture */
  createBlxColor(bc->defaultColors, BLXCOLOR_GRID_LINE, "Grid lines", "Big Picture grid lines", BLX_YELLOW, BLX_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_GRID_TEXT, "Grid text", "Big Picture grid text", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_HIGHLIGHT_BOX, "Highlight box", "Highlight box in the big picture", BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_PREVIEW_BOX, "Preview box", "Preview box in the big picture", BLX_BLACK, BLX_BLACK, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_MSP_LINE, "Big picture match line", "Color of the lines representing matches in the Big Picture", BLX_BLACK, BLX_BLACK, BLX_CYAN, BLX_GREY);

  /* groups */
  createBlxColor(bc->defaultColors, BLXCOLOR_GROUP, "Default group color", "Default highlight color for a new group", BLX_ORANGE_RED, BLX_VERY_LIGHT_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_MATCH_SET, "Default match set color", "Default color for the match set group (applies only when it is created for the first time or after being deleted)", BLX_RED, BLX_VERY_LIGHT_GREY, NULL, NULL);

  /* colinearity */
  createBlxColor(bc->defaultColors, BLXCOLOR_COLINEAR_PERFECT, "Perfect colinearity", "Color of lines joining alignment blocks with perfect colinearity", BLX_DARK_GREEN, BLX_LIGHT_GREY, BLX_DARK_GREEN, BLX_LIGHT_GREY);
  createBlxColor(bc->defaultColors, BLXCOLOR_COLINEAR_IMPERFECT, "Imperfect colinearity", "Color of lines joining alignment blocks with imperfect colinearity", BLX_ORANGE, BLX_GREY, BLX_ORANGE, BLX_GREY);
  createBlxColor(bc->defaultColors, BLXCOLOR_COLINEAR_NOT, "Not colinear", "Color of lines joining alignment blocks that are not colinear", BLX_RED, BLX_DARK_GREY, BLX_RED, BLX_DARK_GREY);

  /* polyA features */
  createBlxColor(bc->defaultColors, BLXCOLOR_POLYA_TAIL, "polyA tail", "polyA tail", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_POLYA_SIGNAL, "polyA signal", "polyA signal", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_POLYA_SIGNAL_ANN, "Annotated polyA signal", "Annotated polyA signal", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_POLYA_SITE_ANN, "Annotated polyA site", "Annotated polyA site", BLX_RED, BLX_DARK_GREY, NULL, NULL);

  /* misc */
  createBlxColor(bc->defaultColors, BLXCOLOR_UNALIGNED_SEQ, "Unaligned sequence", "Addition sequence in the match that is not part of the alignment", "#FFC432", BLX_WHITE, "#FFE8AD", NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_CANONICAL, "Canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are canonical", BLX_GREEN, BLX_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_NON_CANONICAL, "Non-canonical intron bases", "The two bases at the start/end of the intron for the selected MSP are colored this color if they are not canonical", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_TREE_GRID_LINES, "Tree grid lines", "Tree grid lines", BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY, BLX_VERY_DARK_GREY);
  createBlxColor(bc->defaultColors, BLXCOLOR_CLIP_MARKER, "Clipped-match indicator", "Marker to indicate a match has been clipped to the display range", BLX_RED, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_COVERAGE_PLOT, "Coverage plot", "Coverage plot", BLX_ROYAL_BLUE, BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_ASSEMBLY_GAP, "Assembly gaps", "Highlight color for assembly gaps", "#D14553", BLX_DARK_GREY, NULL, NULL);
  createBlxColor(bc->defaultColors, BLXCOLOR_SELECTION, "Selection color", "Highlight color for selections", BLX_DARK_GREY, BLX_DARK_GREY, NULL, NULL);
  
  g_free(defaultBgColorStr);
}


/* Whether to include the given msp type in depth coverage calculations */
static gboolean includeTypeInCoverage(BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH);
}


/* Calculate the depth of coverage of short-reads for each reference sequence display coord.
 * depthArray must be the same length as displayRange. */
static void calculateDepth(BlxViewContext *bc)
{
  /* Allocate the depth array, if null */
  const int displayLen = getRangeLength(&bc->fullDisplayRange);
  
  if (displayLen < 1)
    return; 
  
  bc->depthArray = (int*)g_malloc(sizeof(int) * displayLen);
  
  /* Initialise each entry to zero */  
  int i = 0;
  for ( ; i < displayLen; ++i)
    {
      bc->depthArray[i] = 0;
    }
  
  /* Loop through all MSP lists */
  int mspType = 0;
  
  for ( ; mspType < BLXMSP_NUM_TYPES; ++mspType)
    {
      /* Only include MSPs of relevant types */
      if (!includeTypeInCoverage((BlxMspType)mspType))
        continue;
      
      /* Loop through all MSPs in this list */
      GArray *mspArray = bc->featureLists[mspType];
      const int fullDisplayLen = getRangeLength(&bc->fullDisplayRange);
    
      i = 0;
      const MSP *msp = mspArrayIdx(mspArray, i);
  
      for ( ; msp; msp = mspArrayIdx(mspArray, ++i))
        {
          /* For each ref-seq coord that this alignment spans, increment the depth */
          int alignIdx = msp->displayRange.min;
          for ( ; alignIdx <= msp->displayRange.max; ++alignIdx)
            {
              /* Convert the msp coord to a zero-based coord. Note that parts of the
               * msp range may be outside the ref seq range. */
              const int displayIdx = alignIdx - bc->fullDisplayRange.min;
            
              if (displayIdx >= 0 && displayIdx < fullDisplayLen)
                bc->depthArray[displayIdx] += 1;
            }
        }
    } 
  
  /* Find the max and min depth */
  bc->minDepth = bc->depthArray[0];
  bc->maxDepth = bc->depthArray[0];
  
  for (i = 1 ; i < displayLen; ++i)
    {
      if (bc->depthArray[i] < bc->minDepth)
        bc->minDepth = bc->depthArray[i];
      
      if (bc->depthArray[i] > bc->maxDepth)
        bc->maxDepth = bc->depthArray[i];
    }  
}


/* Called on startup to set the initial state of the flags. Gets the state for
 * the settings from the config file if specified, otherwises uses hard-coded
 * defaults. */
static void initialiseFlags(BlxViewContext *blxContext, CommandLineOptions *options)
{
  /* Initialise all the flags to false */
  int flag = BLXFLAG_MIN + 1;
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      blxContext->flags[flag] = FALSE;
    }
  
  /* Set any specific flags that we want initialised to TRUE */
  blxContext->flags[BLXFLAG_LIMIT_UNALIGNED_BASES] = TRUE;
  blxContext->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED] = TRUE;
  blxContext->flags[BLXFLAG_SHOW_POLYA_SIG_SELECTED] = TRUE;
  blxContext->flags[BLXFLAG_SHOW_SPLICE_SITES] = TRUE;
  blxContext->flags[BLXFLAG_NEGATE_COORDS] = options->negateCoords;
  blxContext->flags[BLXFLAG_HIGHLIGHT_DIFFS] = options->highlightDiffs;
  blxContext->flags[BLXFLAG_SAVE_TEMP_FILES] = options->saveTempFiles;
  blxContext->flags[BLXFLAG_ABBREV_TITLE] = options->abbrevTitle;
}


/* load settings from the config file */
static void loadBlixemSettings(BlxViewContext *blxContext)
{
  /* Override the defaults settings with those given in the config file, if any */
  GKeyFile *key_file = blxGetConfig();

  if (!key_file)
    return;
  
  GError *error = NULL;

  /* squash-matches */
  int squashMatches = g_key_file_get_integer(blxGetConfig(), SETTINGS_GROUP, SETTING_NAME_SQUASH_MATCHES, &error);
  
  if (error)
    {
      /* we don't care if it wasn't found; just clear the error */
      g_error_free(error);
      error = NULL;
    }
  else
    {
      blxContext->modelId = squashMatches ? BLXMODEL_SQUASHED : BLXMODEL_NORMAL;
    }

  /* loop through all the flags and see if any of them are given */
  int flag = BLXFLAG_MIN + 1;
  
  for ( ; flag < BLXFLAG_NUM_FLAGS; ++flag)
    {
      const char *flagName = getFlagName((BlxFlag)flag);
      
      if (flagName)
        {
          int result = g_key_file_get_integer(key_file, SETTINGS_GROUP, flagName, &error);
          
          if (error)
            {
              g_error_free(error);
              error = NULL;
            }
          else
            {
              blxContext->flags[flag] = result;
            }
        }
    }
}



static BlxViewContext* blxWindowCreateContext(CommandLineOptions *options,
                                              const IntRange* const refSeqRange,
                                              const IntRange* const fullDisplayRange,
                                              const char *paddingSeq,
                                              GArray* featureLists[],
                                              GList *seqList,
                                              GSList *supportedTypes,
                                              GtkWidget *widget,
                                              GtkWidget *statusBar,
                                              const gboolean External,
                                              GSList *styles)
{
  BlxViewContext *blxContext = (BlxViewContext*)g_malloc(sizeof *blxContext);
  
  blxContext->statusBar = statusBar;
  
  blxContext->refSeq = options->refSeq;
  blxContext->refSeqName = options->refSeqName ? g_strdup(options->refSeqName) : g_strdup("Blixem-seq");
  blxContext->refSeqRange.min = refSeqRange->min;
  blxContext->refSeqRange.max = refSeqRange->max;
  blxContext->fullDisplayRange.min = fullDisplayRange->min;
  blxContext->fullDisplayRange.max = fullDisplayRange->max;
  blxContext->refSeqOffset = options->refSeqOffset;
  blxContext->optionalColumns = options->optionalColumns;

  blxContext->mspList = options->mspList;
  blxContext->columnList = options->columnList;
  blxContext->styles = styles;
  
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      blxContext->featureLists[typeId] = featureLists[typeId];
    }
  
  blxContext->geneticCode = options->geneticCode;
  blxContext->blastMode = options->blastMode;
  blxContext->seqType = options->seqType;
  blxContext->numFrames = options->numFrames;
  blxContext->paddingSeq = paddingSeq;
  blxContext->bulkFetchDefault = options->bulkFetchDefault;
  blxContext->userFetchDefault = options->userFetchDefault;
  blxContext->optionalFetchDefault = options->optionalFetchDefault;
  blxContext->fetchMethods = options->fetchMethods;
  blxContext->dataset = g_strdup(options->dataset);
  blxContext->matchSeqs = seqList;
  blxContext->supportedTypes = supportedTypes;
  
  blxContext->displayRev = FALSE;
  blxContext->external = External;
  
  blxContext->selectedSeqs = NULL;
  blxContext->sequenceGroups = NULL;
  blxContext->matchSetGroup = NULL;
  
  blxContext->autoDotter = TRUE;
  blxContext->dotterSelf = FALSE;
  blxContext->dotterHsps = FALSE;
  blxContext->dotterStart = UNSET_INT;
  blxContext->dotterEnd = UNSET_INT;
  blxContext->dotterZoom = 0;
  
  blxContext->defaultColors = NULL;
  blxContext->usePrintColors = FALSE;
  blxContext->windowColor = options->windowColor;
  
  createBlxColors(blxContext, widget);
  
  initialiseFlags(blxContext, options);
    
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < BLXDIALOG_NUM_DIALOGS; ++dialogId)
    {
      blxContext->dialogList[dialogId] = NULL;
    }
    
  blxContext->spawnedProcesses = NULL;
  blxContext->depthArray = NULL;
  blxContext->minDepth = 0;
  blxContext->maxDepth = 0;
 
  loadBlixemSettings(blxContext);

  /* do this after loading settings because the passed-in squashed 
   * matches option should override the saved option in the settings */
  blxContext->modelId = options->squashMatches ? BLXMODEL_SQUASHED : BLXMODEL_NORMAL;

  return blxContext;
}


/* Create the properties struct and initialise all values. */
static void blxWindowCreateProperties(CommandLineOptions *options,
                                      BlxViewContext *blxContext,
                                      GtkWidget *widget, 
                                      GtkWidget *bigPicture, 
                                      GtkWidget *detailView,
                                      GtkWidget *mainmenu,
                                      GtkActionGroup *actionGroup,
                                      const IntRange* const refSeqRange,
                                      const IntRange* const fullDisplayRange,
                                      const char *paddingSeq)
{
  if (widget)
    {
      BlxWindowProperties *properties = (BlxWindowProperties*)g_malloc(sizeof *properties);
      
      properties->blxContext = blxContext;

      properties->bigPicture = bigPicture;
      properties->detailView = detailView;
      properties->mainmenu = mainmenu;
      properties->actionGroup = actionGroup;

      properties->pageSetup = gtk_page_setup_new();
      gtk_page_setup_set_orientation(properties->pageSetup, GTK_PAGE_ORIENTATION_LANDSCAPE);
      
      properties->printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(properties->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
      gtk_print_settings_set_quality(properties->printSettings, GTK_PRINT_QUALITY_HIGH);
      gtk_print_settings_set_resolution(properties->printSettings, DEFAULT_PRINT_RESOLUTION);
    
      g_object_set_data(G_OBJECT(widget), "BlxWindowProperties", properties);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK (onDestroyBlxWindow), NULL);
    }
}

gboolean blxWindowGetDisplayRev(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->displayRev : FALSE;
}

GtkWidget* blxWindowGetBigPicture(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->bigPicture : NULL;
}

GtkWidget* blxWindowGetDetailView(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->detailView : NULL;
}

GtkWidget* blxWindowGetCoverageView(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? bigPictureGetCoverageView(properties->bigPicture) : NULL;
}

GtkWidget* blxWindowGetMainMenu(GtkWidget *blxWindow)
{
  BlxWindowProperties *properties = blxWindowGetProperties(blxWindow);
  return properties ? properties->mainmenu : NULL;
}

BlxBlastMode blxWindowGetBlastMode(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->blastMode : (BlxBlastMode)0;
}

char * blxWindowGetRefSeq(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->refSeq : NULL;
}

const char * blxWindowGetRefSeqName(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->refSeqName : NULL;
}

char** blxWindowGetGeneticCode(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->geneticCode : NULL;
}

MSP* blxWindowGetMspList(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->mspList : NULL;
}

GList* blxWindowGetAllMatchSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->matchSeqs : NULL;
}

BlxSeqType blxWindowGetSeqType(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->seqType : BLXSEQ_NONE;
}

IntRange* blxWindowGetFullRange(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? &blxContext->fullDisplayRange : NULL;
}

IntRange* blxWindowGetRefSeqRange(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? &blxContext->refSeqRange : NULL;
}

int blxWindowGetNumFrames(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->numFrames : UNSET_INT;
}

int blxWindowGetAutoDotter(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->autoDotter : TRUE;
}

int blxWindowGetDotterStart(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterStart : UNSET_INT;
}

int blxWindowGetDotterEnd(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterEnd : UNSET_INT;
}

int blxWindowGetDotterZoom(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->dotterZoom : UNSET_INT;
}

const char* blxWindowGetPaddingSeq(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->paddingSeq : NULL;
}

/* Return the active strand - forward strand by default, reverse strand if display toggled */
BlxStrand blxWindowGetActiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? BLXSTRAND_REVERSE : BLXSTRAND_FORWARD;
}

/* Return the inactive strand - reverse strand by default, forward strand if display toggled */
static BlxStrand blxWindowGetInactiveStrand(GtkWidget *blxWindow)
{
  return blxWindowGetDisplayRev(blxWindow) ? BLXSTRAND_FORWARD : BLXSTRAND_REVERSE;
}

/* returns true if display coords should be negated */
gboolean blxWindowGetNegateCoords(GtkWidget *blxWindow)
{
  /* We negate coords (literally just stick a '-' on the front) for display purposes if the
   * display is reversed and the negate-coords option is enabled. This gives the effect that coords
   * always increase left-to-right, whereas when the display is reversed they really decrease. */
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  return (bc->displayRev && bc->flags[BLXFLAG_NEGATE_COORDS]);
}

/* Get the column info for a particular column */
BlxColumnInfo *getColumnInfo(GList *columnList, const BlxColumnId columnId)
{
  BlxColumnInfo *result = NULL;
  
  GList *listItem = columnList;
  
  for ( ; listItem; listItem = listItem->next)
  {
    BlxColumnInfo *columnInfo = (BlxColumnInfo*)(listItem->data);
    if (columnInfo && columnInfo->columnId == columnId)
      {
        result = columnInfo;
        break;
      }
  }
  
  return result;
}

/* Returns the list of all sequence groups */
GList *blxWindowGetSequenceGroups(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext->sequenceGroups;
}

/* Returns the group that the given sequence belongs to, if any (assumes the sequence
 * is only in one group; otherwise it just returns the first group it finds). */
SequenceGroup *blxContextGetSequenceGroup(const BlxViewContext *bc, const BlxSequence *seqToFind)
{
  SequenceGroup *result = NULL;
  
  if (!seqToFind)
    return result;
  
  /* Loop through all the groups until we find this sequence in one */
  GList *groupItem = bc->sequenceGroups;
  for ( ; groupItem; groupItem = groupItem->next)
    {
      /* See if our sequence struct is in this group's list */
      SequenceGroup *group = (SequenceGroup*)(groupItem->data);
      GList *foundItem = g_list_find(group->seqList, seqToFind);
      
      if (foundItem)
        {
          result = group;
          break;
        }
    }
  
  return result;
}


/* Returns the group that the given sequence belongs to, if any (assumes the sequence
 * is only in one group; otherwise it just returns the first group it finds). */
SequenceGroup *blxWindowGetSequenceGroup(GtkWidget *blxWindow, const BlxSequence *seqToFind)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  return blxContextGetSequenceGroup(bc, seqToFind);
}


static gboolean blxWindowGetUsePrintColors(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  return bc->usePrintColors;
}


/* This should be called whenever the background color has changed */
static void onUpdateBackgroundColor(GtkWidget *blxWindow)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  GdkColor *defaultBgColor = getGdkColor(BLXCOLOR_BACKGROUND, bc->defaultColors, FALSE, bc->usePrintColors);
  setWidgetBackgroundColor(blxWindow, defaultBgColor);
  
  blxWindowRedrawAll(blxWindow);
}


/* This sets the 'use print colors' flag and then updates the display */
static void blxWindowSetUsePrintColors(GtkWidget *blxWindow, const gboolean usePrintColors)
{
  BlxViewContext *bc = blxWindowGetContext(blxWindow);
  bc->usePrintColors = usePrintColors;
  onUpdateBackgroundColor(blxWindow);
}


/***********************************************************
 *                        Selections                       *
 ***********************************************************/

GList* blxWindowGetSelectedSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContext ? blxContext->selectedSeqs : NULL;
}


/* Get the selected sequences as a list of sequence names. The returned list
 * is formatted as newline-separated values. */
static GString* blxWindowGetSelectedSeqNames(GtkWidget *blxWindow)
{
  GList *listItem = blxWindowGetSelectedSeqs(blxWindow);
  GString *result = g_string_new_len(NULL, 50);
  gboolean first = TRUE;
  
  for ( ; listItem; listItem = listItem->next)
    {
      /* Add a separator before the name, unless it's the first one */
      if (!first)
        {
          g_string_append(result, "\n");
        }
      else
        {
          first = FALSE;
        }

      const BlxSequence *seq = (const BlxSequence*)(listItem->data);
      g_string_append(result, blxSequenceGetName(seq));
    }

  return result;
}


/* Get the selected sequence data. Only works for a single selection.
 * Returns null if fails. */
static const char* blxWindowGetSelectedSeqData(GtkWidget *blxWindow)
{
  GList *listItem = blxWindowGetSelectedSeqs(blxWindow);
  const char *result = NULL;

  if (g_list_length(listItem) < 1)
    g_critical("Please select a sequence.\n");
  else if (g_list_length(listItem) > 1)
    g_critical("Please select a single sequence.\n");
  else
    {
      const BlxSequence *seq = (const BlxSequence*)(listItem->data);
      result = blxSequenceGetSequence(seq);
    }

  return result;
}


/* This function copies the currently-selected sequences' names to the default
 * clipboard. */
void copySelectionToClipboard(GtkWidget *blxWindow)
{
  if (g_list_length(blxWindowGetSelectedSeqs(blxWindow)) < 1)
    {
      g_critical("Please select a sequence");
    }
  else
    {
      GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);

      if (displayText)
        {
          setDefaultClipboardText(displayText->str);
          g_string_free(displayText, TRUE);
          g_message("Copied selected sequence name(s) to clipboard\n");
        }
    }
}


/* This function copies the currently-selected sequence's data to the default
 * clipboard. */
static void copySelectedSeqDataToClipboard(GtkWidget *blxWindow)
{
  const char *displayText = blxWindowGetSelectedSeqData(blxWindow);

  if (displayText)
    {
      const int len = strlen(displayText);
      
      /* Warn user if they're about to copy a large sequence */
      if (len <= MAX_RECOMMENDED_COPY_LENGTH || 
          runConfirmationBox(blxWindow, "Copy sequence", "You are about to copy a large amount of text to the clipboard\n\nAre you sure you want to continue?") == GTK_RESPONSE_ACCEPT)
        {
          setDefaultClipboardText(displayText);
          g_message("Copied selected sequence data to clipboard\n");
        }
    }
}


/* This function copies the reference sequence, from the 
 * clicked position to the marked position, onto the clipboard. */
static void copyRefSeqToClipboard(GtkWidget *blxWindow, const int fromIdx_in, const int toIdx_in)
{
  const char *refSeq = blxWindowGetRefSeq(blxWindow);

  if (refSeq)
    {
      /* Need to get 0-based indices */
      const IntRange* const refSeqRange = blxWindowGetRefSeqRange(blxWindow);
      
      const int fromIdx = min(fromIdx_in, toIdx_in) - refSeqRange->min;
      const int toIdx = max(fromIdx_in, toIdx_in) - refSeqRange->min;
      const int len = toIdx - fromIdx + 1;

      /* Warn user if they're about to copy a large sequence */
      if (len <= MAX_RECOMMENDED_COPY_LENGTH || 
          runConfirmationBox(blxWindow, "Copy sequence", "You are about to copy a large amount of text to the clipboard\n\nAre you sure you want to continue?") == GTK_RESPONSE_ACCEPT)
        {
          const char *displayText = g_strndup(refSeq + fromIdx, toIdx - fromIdx + 1);
      
          if (displayText)
            {
              setDefaultClipboardText(displayText);
              g_message("Copied reference sequence from %d to %d\n", fromIdx_in, toIdx_in);
            }
        }
    }
}


/* Update function to be called whenever the MSP selection has changed */
void blxWindowSelectionChanged(GtkWidget *blxWindow)
{
  GtkWidget *detailView = blxWindowGetDetailView(blxWindow);
  BlxViewContext *bc = blxWindowGetContext(blxWindow);

  if ((bc->flags[BLXFLAG_SHOW_UNALIGNED] && bc->flags[BLXFLAG_SHOW_UNALIGNED_SELECTED]) ||
      (bc->flags[BLXFLAG_SHOW_POLYA_SITE] && bc->flags[BLXFLAG_SHOW_POLYA_SITE_SELECTED]))
    {
      /* When these options are enabled, the length of an MSP can depend on
       * whether it is selected or not; changing the selection can therefore
       * affect which MPSs are visible, so we need to re-filter the trees */
      refilterDetailView(detailView, NULL);
    }
  
  /* Redraw */
  updateFeedbackBox(detailView);
  blxWindowRedrawAll(blxWindow);
  
  /* Copy the selected sequence names to the PRIMARY clipboard */
  GString *displayText = blxWindowGetSelectedSeqNames(blxWindow);
  setPrimaryClipboardText(displayText->str);
  g_string_free(displayText, TRUE);
  
  /* Refresh the dotter dialog, if it happens to be open */
  refreshDialog(BLXDIALOG_DOTTER, blxWindow);
}


/* Call this function to select the given match sequence */
void blxWindowSelectSeq(GtkWidget *blxWindow, BlxSequence *seq)
{
  if (!blxWindowIsSeqSelected(blxWindow, seq))
    {
      BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
      blxContext->selectedSeqs = g_list_append(blxContext->selectedSeqs, seq);
      blxWindowSelectionChanged(blxWindow);
    }
}


/* Utility function to set the list of selected sequences */
void blxWindowSetSelectedSeqList(GtkWidget *blxWindow, GList *seqList)
{
  blxWindowDeselectAllSeqs(blxWindow);

  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  blxContext->selectedSeqs = seqList;

  blxWindowSelectionChanged(blxWindow);
}


/* Call this function to deselect the given sequence */
void blxWindowDeselectSeq(GtkWidget *blxWindow, BlxSequence *seq)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  /* See if it's in the list and, if so, get a pointer to the list element */
  GList *foundSeq = g_list_find(blxContext->selectedSeqs, seq);

  if (foundSeq)
    {
      blxContext->selectedSeqs = g_list_remove(blxContext->selectedSeqs, foundSeq->data);
      blxWindowSelectionChanged(blxWindow);
    }
}


/* Call this function to deselect all sequences */
void blxWindowDeselectAllSeqs(GtkWidget *blxWindow)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);

  if (g_list_length(blxContext->selectedSeqs) > 0)
    {
      g_list_free(blxContext->selectedSeqs);
      blxContext->selectedSeqs = NULL;
      blxWindowSelectionChanged(blxWindow);
    }
}


gboolean blxContextIsSeqSelected(const BlxViewContext* const bc, const BlxSequence *seq)
{
  GList *foundItem = NULL;
  
  if (bc && seq)
    {
      foundItem = g_list_find(bc->selectedSeqs, seq);
    }
  
  return (foundItem != NULL);
}


/* Returns true if the given sequence is selected */
gboolean blxWindowIsSeqSelected(GtkWidget *blxWindow, const BlxSequence *seq)
{
  BlxViewContext *blxContext = blxWindowGetContext(blxWindow);
  return blxContextIsSeqSelected(blxContext, seq);
}


/* Set the given sequence as selected or unselected, depending on the given argument */
void blxWindowSetSeqSelected(GtkWidget *blxWindow, BlxSequence *seq, const gboolean selected)
{
  if (selected)
    {
      blxWindowSelectSeq(blxWindow, seq);
    }
  else
    {
      blxWindowDeselectSeq(blxWindow, seq);
    }
}


BlxSequence* blxWindowGetLastSelectedSeq(GtkWidget *blxWindow)
{
  BlxSequence *result = NULL;
  
  GList *selectedSeqs = blxWindowGetSelectedSeqs(blxWindow);
  
  if (g_list_length(selectedSeqs) > 0)
    {
      /* Get the last-selected sequence */
      GList *lastItem = g_list_last(selectedSeqs);
      result = (BlxSequence*)(lastItem->data);
    }
  
  return result;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

static void onDragDataReceived(GtkWidget *widget, 
                               GdkDragContext *context, 
                               int x, 
                               int y,
                               GtkSelectionData *selectionData, 
                               guint info, 
                               guint time,
                               gpointer userdata)
{
  DEBUG_ENTER("onDragDataReceived()");

  g_return_if_fail(selectionData);

  if ((info == TARGET_STRING || info == TARGET_URL) && selectionData->data)
    {
      DEBUG_OUT("Received drag and drop text:\n%s\n", selectionData->data);
      GError *tmp_error = NULL;
      
      /* For now just assume the text contains supported file contents. The file parsing
       * will fail if it's not a supported format. */
      char *text = (char*)(gtk_selection_data_get_text(selectionData));
      dynamicLoadFeaturesFile(widget, NULL, text, &tmp_error);

      if (tmp_error)
        {
          prefixError(tmp_error, "Error processing text from drag-and-drop: ");
          reportAndClearIfError(&tmp_error, G_LOG_LEVEL_CRITICAL);
        }
    }

  DEBUG_EXIT("onDragDataReceived returning ");
}


/* Called when the user moves the cursor over the window during a drag */
static gboolean onDragMotion(GtkWidget *widget, GdkDragContext *event, gint x, gint y, guint time, gpointer data)
{
  /* Bring the blixem window to the front so the user can see where they're going to drop */
  gtk_window_present(GTK_WINDOW(widget));

  /* Return true to indicate that the whole window is a drop zone */
  return TRUE;
}


static void setDragDropProperties(GtkWidget *widget)
{
  DEBUG_ENTER("setDragDropProperties()");

  static GtkTargetEntry targetentries[] =
    {
      { "STRING",        0, TARGET_STRING },
      { "text/plain",    0, TARGET_STRING },
      { "text/uri-list", 0, TARGET_URL },
    };
  
  gtk_drag_dest_set(widget, GTK_DEST_DEFAULT_ALL, targetentries, 3,
                    GDK_ACTION_COPY|GDK_ACTION_MOVE|GDK_ACTION_LINK);
 
  g_signal_connect(widget, "drag-data-received", G_CALLBACK(onDragDataReceived), NULL);
  g_signal_connect(widget, "drag-motion", G_CALLBACK(onDragMotion),	NULL);
  
  DEBUG_EXIT("setDragDropProperties returning ")
}

/* Set various properties for the blixem window */
static void setStyleProperties(GtkWidget *widget, char *windowColor)
{
  DEBUG_ENTER("setStyleProperties()");

  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(widget);
  const int width = gdk_screen_get_width(screen) * DEFAULT_WINDOW_WIDTH_FRACTION;
  const int height = gdk_screen_get_height(screen) * DEFAULT_WINDOW_HEIGHT_FRACTION;
  
  gtk_window_set_default_size(GTK_WINDOW(widget), width, height);
  
  gtk_container_set_border_width (GTK_CONTAINER(widget), DEFAULT_WINDOW_BORDER_WIDTH); 
  gtk_window_set_mnemonic_modifier(GTK_WINDOW(widget), GDK_MOD1_MASK); /* MOD1 is ALT on most systems */
  
  /* Set the default font size to be a bit smaller than usual */
  int origSize = pango_font_description_get_size(widget->style->font_desc) / PANGO_SCALE;
  const char *origFamily = pango_font_description_get_family(widget->style->font_desc);

  char parseString[500];
  sprintf(parseString, "gtk-font-name = \"%s %d\"", origFamily, origSize + DEFAULT_FONT_SIZE_ADJUSTMENT);
  gtk_rc_parse_string(parseString);

  DEBUG_EXIT("setStyleProperties returning ");
}


/* Create the main menu */
static void createMainMenu(GtkWidget *window, BlxViewContext *bc, GtkWidget **mainmenu, GtkWidget **toolbar, GtkActionGroup **actionGroupOut)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  if (actionGroupOut)
    *actionGroupOut = action_group;

  /* Set the squash-matches toggle button depending on the initial state */
  toggleMenuEntries[0].is_active = (bc->modelId == BLXMODEL_SQUASHED);  

  gtk_action_group_add_actions (action_group, mainMenuEntries, G_N_ELEMENTS (mainMenuEntries), window);
  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS (toggleMenuEntries), window);

  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);
  
  GError *error = NULL;
  gboolean ok = gtk_ui_manager_add_ui_from_string (ui_manager, standardMenuDescription, -1, &error);
  
  if (ok && userIsDeveloper())
    ok = gtk_ui_manager_add_ui_from_string (ui_manager, developerMenuDescription, -1, &error);
    
  if (!ok)
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  *mainmenu = gtk_ui_manager_get_widget (ui_manager, "/ContextMenu");
  *toolbar = gtk_ui_manager_get_widget (ui_manager, "/Toolbar");
}


/* calcID: caculated percent identity of an MSP
 * 
 * There seems to be a general problem with this routine for protein
 * alignments, the existing code certainly does not do the right thing.
 * I have fixed this routine for gapped sequence alignments but not for
 * protein stuff at all.
 * 
 * To be honest I think this routine is a _waste_ of time, the alignment
 * programs that feed data to blixem produce an identity anyway so why
 * not use that...why reinvent the wheel......
 * 
 * */
static void calcID(MSP *msp, BlxViewContext *bc)
{
  const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
  const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
  
  if (mspIsBlastMatch(msp) && msp->id < 0) /* Only calculate if ID is not already set */
    {
      msp->id = 0.0;

      /* If there is no sequence data, leave the ID as zero */
      const char *matchSeq = mspGetMatchSeq(msp);

      if (matchSeq)
        {
          /* Note that this will reverse complement the ref seq if it is the reverse 
           * strand. This means that where there is no gaps array the comparison is trivial
           * as coordinates can be ignored and the two sequences just whipped through. */
          GError *error = NULL;
          IntRange qRange;
          intrangeSetValues(&qRange, msp->qRange.min, msp->qRange.max); /* make a copy because it will be updated */
          
          char *refSeqSegment = getSequenceSegment(bc->refSeq,
                                                   &qRange,
                                                   mspGetRefStrand(msp), 
                                                   BLXSEQ_DNA,        /* msp q coords are always nucleotide coords */
                                                   bc->seqType,       /* required seq type is the display seq type */
                                                   mspGetRefFrame(msp, bc->seqType),
                                                   bc->numFrames,
                                                   &bc->refSeqRange,
                                                   bc->blastMode,
                                                   bc->geneticCode,
                                                   bc->displayRev,
                                                   !qForward,
                                                   TRUE,
                                                   &error);
          
          if (!refSeqSegment)
            {
              prefixError(error, "Failed to calculate ID for sequence '%s' (match coords = %d - %d). ", mspGetSName(msp), msp->sRange.min, msp->sRange.max);
              reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
              return;
            }
          else
            {
              /* If there's an error but the sequence was still returned it's 
               * a non-critical warning. Only issue one warning because we can
               * get many thousands and it can fill up the terminal if we output
               * them all. */
              if (error)
                {
                  static gboolean done = FALSE;
                  
                  if (!done)
                    {
                      g_warning("There were errors calculating the percent ID for some sequences because the match extends out of the reference sequence range; some IDs may be incorrect.\n");
                      done = TRUE;
                    }
                  
                  g_error_free(error);
                  error = NULL;
                }
            }
          
          /* We need to find the number of characters that match out of the total number */
          int numMatchingChars = 0;
          int totalNumChars = 0;
          const int numGaps = msp->gaps ? g_slist_length(msp->gaps) : 0;
          
          if (numGaps == 0)
            {
              /* Ungapped alignments. */
              totalNumChars = (qRange.max - qRange.min + 1) / bc->numFrames;

              if (bc->blastMode == BLXMODE_TBLASTN || bc->blastMode == BLXMODE_TBLASTX)
                {
                  int i = 0;
                  for ( ; i < totalNumChars; i++)
                    {
                      if (toupper(matchSeq[i]) == toupper(refSeqSegment[i]))
                        {
                          numMatchingChars++;
                        }
                    }
                }
              else                                                  /* blastn, blastp & blastx */
                {
                  int i = 0;
                  for ( ; i < totalNumChars; i++)
                    {
                      int sIndex = sForward ? msp->sRange.min + i - 1 : msp->sRange.max - i - 1;
                      if (toupper(matchSeq[sIndex]) == toupper(refSeqSegment[i]))
                        {
                          numMatchingChars++;
                        }
                    }
                }
            }
          else
            {
              /* Gapped alignments. */

              /* To do tblastn and tblastx is not imposssible but would like to work from
               * examples to get it right.... */
              if (bc->blastMode == BLXMODE_TBLASTN)
                {
                  g_message("not implemented yet\n") ;
                }
              else if (bc->blastMode == BLXMODE_TBLASTX)
                {
                  g_message("not implemented yet\n") ;
                }
              else
                {
                  /* blastn and blastp remain simple but blastx is more complex since the query
                   * coords are nucleic not protein. */
                  GSList *rangeItem = msp->gaps;
                  
                  for ( ; rangeItem; rangeItem = rangeItem->next)
                    {
                      CoordRange *range = (CoordRange*)(rangeItem->data);
                      
                      int qRangeMin = 0, qRangeMax = 0, sRangeMin = 0, sRangeMax = 0;
                      getCoordRangeExtents(range, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
                      
                      totalNumChars += sRangeMax - sRangeMin + 1;
                      
                      /* Note that refSeqSegment is just the section of the ref seq relating to this msp.
                       * We need to translate the first coord in the range (which is in terms of the full
                       * reference sequence) into coords in the cut-down ref sequence. */
                      int q_start = qForward ? (qRangeMin - qRange.min) / bc->numFrames : (qRange.max - qRangeMax) / bc->numFrames;
                      const int qLen = strlen(refSeqSegment);
                      
                      /* We can index sseq directly (but we need to adjust by 1 for zero-indexing). We'll loop forwards
                       * through sseq if we have the forward strand or backwards if we have the reverse strand,
                       * so start from the lower or upper end accordingly. */
                      int s_start = sForward ? sRangeMin - 1 : sRangeMax - 1 ;

                      int sIdx = s_start, qIdx = q_start ;
                      while (((sForward && sIdx < sRangeMax) || (!sForward && sIdx >= sRangeMin - 1)) && qIdx < qLen)
                        {
                          /* Check that qIdx is not less that 0, which could happen if we have somehow got duff data. */
                          if (qIdx >= 0 && toupper(matchSeq[sIdx]) == toupper(refSeqSegment[qIdx]))
                            {
                              numMatchingChars++ ;
                            }

                          /* Move to the next base. The refSeqSegment is always forward, but we might have to
                           * traverse the s sequence in reverse. */
                          ++qIdx ;
                          if (sForward) ++sIdx ;
                          else --sIdx ;
                        }
                    }
                }
            }
          
          msp->id = (100.0 * numMatchingChars / totalNumChars);
          
          g_free(refSeqSegment);
        }
    }
  
  return ;
}


/* Calculate the ID and q frame for the given MSP and store 
 * them in the MSP struct. Returns the calculated ID (or UNSET_INT if this msp
 * is not a blast match). */
static void calcMspData(MSP *msp, BlxViewContext *bc)
{  
  /* Calculate the ID */
  if (mspIsBlastMatch(msp))
    {
      calcID(msp, bc);
    }
}


static gdouble calculateMspData(MSP *mspList, BlxViewContext *bc)
{
  MSP *msp = mspList;
  gdouble lowestId = -1.0;
  
  for ( ; msp; msp = msp->next)
    {
      calcMspData(msp, bc);
      
      if (mspIsBlastMatch(msp) && (lowestId == -1.0 || msp->id < lowestId))
        {
          lowestId = msp->id;
        }
    }
  
  return lowestId;
}


/* Calculate the reference sequence range from the range and offset given in 
 * the option. Also translate this to display coords. */
static void calculateRefSeqRange(CommandLineOptions *options,
                                 IntRange *refSeqRange,
                                 IntRange *fullDisplayRange)
{
  
  /* Offset the reference sequence range, if an offset was specified. */ 
  refSeqRange->min = options->refSeqRange.min + options->refSeqOffset;
  refSeqRange->max = options->refSeqRange.max + options->refSeqOffset;
  
  fullDisplayRange->min = refSeqRange->min;
  fullDisplayRange->max = refSeqRange->max;
  
  if (options->seqType == BLXSEQ_PEPTIDE)
    {
      /* Adjust the reference sequence reading frame so that it always starts at
        * base 1 in frame 1. This makes the display easier because we can always just
        * start drawing from the 1st base in the reference sequence. */
      int base = UNSET_INT;
      convertDnaIdxToDisplayIdx(refSeqRange->min, options->seqType, 1, options->numFrames, FALSE, refSeqRange, &base);
      
      int offset = (options->numFrames - base + 1);
      
      if (offset >= options->numFrames)
          offset -= options->numFrames;
      
      refSeqRange->min += offset;
      options->refSeq = options->refSeq + offset;
      
      /* Now do the same for when the ref seq is reversed */
      convertDnaIdxToDisplayIdx(refSeqRange->max, options->seqType, 1, options->numFrames, TRUE, refSeqRange, &base);
      offset = (options->numFrames - base + 1);
      
      if (offset >= options->numFrames) 
          offset -= options->numFrames;
      
      refSeqRange->max -= offset;
      options->refSeq[refSeqRange->max - refSeqRange->min + 1] = '\0';
      
      /* Now calculate the full display range in display coords */
      fullDisplayRange->min = convertDnaIdxToDisplayIdx(refSeqRange->min, options->seqType, 1, options->numFrames, FALSE, refSeqRange, NULL);
      fullDisplayRange->max = convertDnaIdxToDisplayIdx(refSeqRange->max, options->seqType, 1, options->numFrames, FALSE, refSeqRange, NULL);
    }  
}


/* Create the main blixem window */
GtkWidget* createBlxWindow(CommandLineOptions *options, 
                           const char *paddingSeq, 
                           GArray* featureLists[],
                           GList *seqList, 
                           GSList *supportedTypes,
                           const gboolean External,
                           GSList *styles)
{
  IntRange refSeqRange;
  IntRange fullDisplayRange;
  
  calculateRefSeqRange(options, &refSeqRange, &fullDisplayRange);
  
  g_message("Reference sequence [%d - %d], display range [%d - %d]\n", 
            refSeqRange.min, refSeqRange.max, fullDisplayRange.min, fullDisplayRange.max);
  
  /* Offset the start coords, if applicable, and convert it to display coords */
  int startCoord = options->startCoord + options->refSeqOffset;
  if (options->seqType == BLXSEQ_PEPTIDE)
    startCoord = convertDnaIdxToDisplayIdx(startCoord, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);

  if (options->bigPictRange.min != UNSET_INT && options->bigPictRange.max != UNSET_INT)
    {
      /* Apply any offset */
      options->bigPictRange.min += options->refSeqOffset;
      options->bigPictRange.max += options->refSeqOffset;

      /* Make sure the big picture range is not outside our ref seq range */
      boundsLimitRange(&options->bigPictRange, &refSeqRange, FALSE);

      if (options->seqType == BLXSEQ_PEPTIDE)
        {
          /* Convert to peptide coords */
          options->bigPictRange.min = convertDnaIdxToDisplayIdx(options->bigPictRange.min, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);
          options->bigPictRange.max = convertDnaIdxToDisplayIdx(options->bigPictRange.max, options->seqType, 1, options->numFrames, FALSE, &refSeqRange, NULL);
        }
    }
  
  
  
  /* Create the main blixem window */
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setStyleProperties(window, options->windowColor);
  setDragDropProperties(window);

  /* Create a status bar */
  GtkWidget *statusBar = gtk_statusbar_new();
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), TRUE);
  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");
  
  /* Set the window and statusbar in the message handler data, now that we know them */
  options->msgData.parent = GTK_WINDOW(window);
  options->msgData.statusBar = GTK_STATUSBAR(statusBar);
  
  /* Create a vertical box to pack everything in */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), vbox);
  
  /* Create the widgets. We need a single adjustment for the entire detail view, which will also be referenced
   * by the big picture view, so create it first. */
  GtkAdjustment *detailAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, /* initial value = 0 */
                                                                      fullDisplayRange.min, /* lower value */
                                                                      fullDisplayRange.max, /* upper value */
                                                                      DEFAULT_SCROLL_STEP_INCREMENT, /* step increment used for mouse wheel scrolling */
                                                                      0,   /* page increment dynamically set based on display range */
                                                                      0)); /* page size dunamically set based on display range */
  
  BlxViewContext *blxContext = blxWindowCreateContext(options, 
                                                      &refSeqRange, 
                                                      &fullDisplayRange, 
                                                      paddingSeq, 
                                                      featureLists,
                                                      seqList, 
                                                      supportedTypes,
                                                      window, 
                                                      statusBar,
                                                      External,
                                                      styles);

  /* Create the main menu */
  GtkWidget *mainmenu = NULL;
  GtkWidget *toolbar = NULL;
  GtkActionGroup *actionGroup = NULL;
  createMainMenu(window, blxContext, &mainmenu, &toolbar, &actionGroup);
  
  const gdouble lowestId = calculateMspData(options->mspList, blxContext);
  
  GtkWidget *fwdStrandGrid = NULL, *revStrandGrid = NULL;

  /* Create the two main sections - the big picture and detail view - in a paned window */
  GtkWidget *panedWin = gtk_vpaned_new();
  gtk_box_pack_start(GTK_BOX(vbox), panedWin, TRUE, TRUE, 0);

  GtkWidget *coverageView = createCoverageView(window, blxContext);
  
  GtkWidget *bigPicture = createBigPicture(window,
                                           GTK_CONTAINER(panedWin),
                                           coverageView,
                                           &fwdStrandGrid, 
                                           &revStrandGrid,
                                           &options->bigPictRange,
                                           &refSeqRange,
                                           options->bigPictZoom,
                                           lowestId);

  GtkWidget *detailView = createDetailView(window,
                                           GTK_CONTAINER(panedWin),
                                           toolbar,
					   detailAdjustment, 
					   fwdStrandGrid, 
					   revStrandGrid,
					   options->mspList,
                                           options->columnList,
					   options->blastMode,
					   options->seqType,
					   options->numFrames,
					   options->refSeqName,
					   startCoord,
					   options->sortInverted,
					   options->initSortColumn,
                                           options->optionalColumns,
                                           options->windowColor);

  
  /* Add the coverage view underneath the main panes */
  gtk_box_pack_start(GTK_BOX(vbox), coverageView, FALSE, FALSE, DEFAULT_COVERAGE_VIEW_BORDER);

  
  /* Create a custom scrollbar for scrolling the sequence column and put it at the bottom of the window */
  GtkWidget *scrollBar = createDetailViewScrollBar(detailAdjustment, detailView);
  gtk_box_pack_start(GTK_BOX(vbox), scrollBar, FALSE, TRUE, 0);

  
  /* Put the statusbar at the bottom */
  gtk_box_pack_start(GTK_BOX(vbox), statusBar, FALSE, TRUE, 0);

  
  /* Set required data for the blixem window */
  blxWindowCreateProperties(options,
                            blxContext,
                            window, 
                            bigPicture, 
                            detailView, 
                            mainmenu,
                            actionGroup,
                            &refSeqRange, 
                            &fullDisplayRange,
                            paddingSeq);
  
  /* Connect signals */
  g_signal_connect(G_OBJECT(panedWin), "button-press-event", G_CALLBACK(onButtonPressPanedWin), window);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressBlxWindow), mainmenu);
  g_signal_connect(G_OBJECT(window), "key-press-event", G_CALLBACK(onKeyPressBlxWindow), NULL);
  
  
  /* Add the MSP's to the trees and sort them by the initial sort mode. This must
   * be done after all widgets have been created, because it accesses their properties.*/
  detailViewAddMspData(detailView, options->mspList, seqList);
  
  /* Updated the cached display range and full extents of the MSPs */
  detailViewUpdateMspLengths(detailView, detailViewGetNumUnalignedBases(detailView));
  cacheMspDisplayRanges(blxContext, detailViewGetNumUnalignedBases(detailView));
  calculateDepth(blxContext);
  updateCoverageDepth(coverageView, blxContext);
  
  /* Set the detail view font (again, this accesses the widgets' properties). */
  updateDetailViewFontDesc(detailView);

  /* Calculate the number of vertical cells in the grids (again, requires properties) */
  calculateNumVCells(bigPicture);


  /* Realise the widgets */
  g_message_info("Starting %s\n", g_get_prgname());
  gtk_widget_show_all(window);


  /* If the options say to make the reverse strand the active strand, toggle the display now */
  if (options->activeStrand == BLXSTRAND_REVERSE)
    toggleStrand(detailView);

  /* Hide the coverage view by default (unless told to display it) */
  if (!options->coverageOn)
    widgetSetHidden(coverageView, TRUE);
  
  /* The trees use the normal model by default, so if we're starting in 
   * 'squash matches' mode we need to change the model */
  if (blxContext->modelId == BLXMODEL_SQUASHED)
    {
      callFuncOnAllDetailViewTrees(detailView, treeUpdateSquashMatches, NULL);
      detailViewRedrawAll(detailView);
    }
  
  /* If the options say to hide the inactive strand, hide it now. (This must be done
   * after showing the widgets, or it will get shown again in show_all.). To do: we just
   * hide the grid at the moment; hide the detail-view pane as well?  */
  if (options->hideInactive && options->activeStrand == BLXSTRAND_FORWARD)
    widgetSetHidden(revStrandGrid, TRUE);
  else if (options->hideInactive && options->activeStrand == BLXSTRAND_REVERSE)
    widgetSetHidden(fwdStrandGrid, TRUE);

  if (!options->bigPictON)
    widgetSetHidden(bigPicture, TRUE);

  if (options->sortInverted)
    {
      blxContext->flags[BLXFLAG_INVERT_SORT] = TRUE;
      detailViewUpdateSortInverted(detailView, options->sortInverted);
    }
  
  /* Set the initial column widths. (This must be called after the widgets are 
   * realised because it causes the scroll range to be updated, which in turn causes
   * the big picture range to be set. The widgets must be realised before this because
   * the initial big picture range depends on the detail view range, which is calculated
   * from its window's width, and this will be incorrect if it has not been realised.) */
  updateDynamicColumnWidths(detailView);
  
  /* Just once, at the start, update the visibility of all tree rows. (After this,
   * filter updates will be done on affected rows only.) */
  callFuncOnAllDetailViewTrees(detailView, refilterTree, NULL);
  detailViewResortTrees(detailView);
  
  /* Calculate initial size of the exon views (depends on big picture range) */
  calculateExonViewHeight(bigPictureGetFwdExonView(bigPicture));
  calculateExonViewHeight(bigPictureGetRevExonView(bigPicture));
  forceResize(bigPicture);
  
  return window;
}
