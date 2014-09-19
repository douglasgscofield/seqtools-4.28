/*  File: belvuWindow.c
 *  Author: Gemma Barson, 2011-04-11
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
 * Description: see belvuWindow.h
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvuWindow.h"
#include "belvuApp/belvuAlignment.h"
#include "belvuApp/belvuTree.h"
#include "belvuApp/belvuConsPlot.h"
#include "belvuApp/belvu_.h"
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <string.h>
#include <unistd.h>


#define DEFAULT_WINDOW_BORDER_WIDTH      1    /* used to change the default border width around the blixem window */
#define DEFAULT_FONT_SIZE_ADJUSTMENT     0   /* used to start with a smaller font than the default widget font */
#define MAIN_BELVU_WINDOW_NAME           "BelvuWindow"
#define WRAPPED_BELVU_WINDOW_NAME        "WrappedBelvuWindow"
#define BELVU_ORGS_WINDOW_NAME           "BelvuOrgsWindow" 
#define DEFAULT_WRAP_WINDOW_WIDTH_FRACTION      0.6    /* default width of wrap window (as fraction of screen width) */
#define DEFAULT_WRAP_WINDOW_HEIGHT_FRACTION     0.85   /* default height of wrap window (as fraction of screen height) */
#define MAX_ORGS_WINDOW_WIDTH_FRACTION          0.5    /* max width of organisms window (as fraction of screen width) */
#define MAX_ORGS_WINDOW_HEIGHT_FRACTION         0.7    /* max height of organisms window (as fraction of screen height) */
#define MAX_ANNOTATION_WINDOW_WIDTH_FRACTION    0.6    /* max width of annotation window (as fraction of screen width) */
#define MAX_ANNOTATION_WINDOW_HEIGHT_FRACTION   0.7    /* max height of annotation window (as fraction of screen height) */
#define ORGS_WINDOW_XPAD                        20     /* x padding for the organisms window */
#define ORGS_WINDOW_YPAD                        20     /* y padding for the organisms window */

/* Utility struct to pass data to the color-changed callback
 * when a color has been changed on the edit-colors dialog */
typedef struct _ColorChangeData
{
  GtkWidget *colorButton;               /* The color-picker button */
  char *residue;                        /* The residue that this color applies to */
} ColorChangeData;


/* Properties specific to the belvu window */
typedef struct _BelvuWindowProperties
  {
    BelvuContext *bc;                   /* The belvu context */
    GtkWidget *statusBar;               /* Message bar at the bottom of the main window */
    GtkWidget *feedbackBox;             /* Feedback area showing info about the current selction */
    GtkActionGroup *actionGroup;
  } BelvuWindowProperties;


/* Properties for generic windows */
typedef struct _GenericWindowProperties
  {
    BelvuContext *bc;                   /* The belvu context */
    GtkActionGroup *actionGroup;
  } GenericWindowProperties;




/* Local function declarations */
static void                      onCloseMenu(GtkAction *action, gpointer data);
static void                      onQuitMenu(GtkAction *action, gpointer data);
static void                      onHelpMenu(GtkAction *action, gpointer data);
static void                      onAboutMenu(GtkAction *action, gpointer data);
static void                      onPrintMenu(GtkAction *action, gpointer data);
static void                      onFindMenu(GtkAction *action, gpointer data);
static void                      onWrapMenu(GtkAction *action, gpointer data);
static void                      onShowTreeMenu(GtkAction *action, gpointer data);
static void                      onRecalcTreeMenu(GtkAction *action, gpointer data);
static void                      onTreeOptsMenu(GtkAction *action, gpointer data);
static void                      onConsPlotMenu(GtkAction *action, gpointer data);
static void                      onSaveMenu(GtkAction *action, gpointer data);
static void                      onSaveAsMenu(GtkAction *action, gpointer data);
static void                      onOutputMenu(GtkAction *action, gpointer data);
static void                      onFetchWWWMenu(GtkAction *action, gpointer data);
static void                      onCompareMenu(GtkAction *action, gpointer data);
static void                      onCleanUpMenu(GtkAction *action, gpointer data);

static void                      onrmPickedMenu(GtkAction *action, gpointer data);
static void                      onRemoveSeqsMenu(GtkAction *action, gpointer data);
static void                      onrmGappySeqsMenu(GtkAction *action, gpointer data);
static void                      onrmPartialSeqsMenu(GtkAction *action, gpointer data);
static void                      onrmRedundantMenu(GtkAction *action, gpointer data);
static void                      onrmOutliersMenu(GtkAction *action, gpointer data);
static void                      onrmScoreMenu(GtkAction *action, gpointer data);
static void                      onrmColumnPromptMenu(GtkAction *action, gpointer data);
static void                      onrmColumnLeftMenu(GtkAction *action, gpointer data);
static void                      onrmColumnRightMenu(GtkAction *action, gpointer data);
static void                      onrmColumnCutoffMenu(GtkAction *action, gpointer data);
static void                      onrmGappyColumnsMenu(GtkAction *action, gpointer data);
static void                      onAutoRmEmptyColumnsMenu(GtkAction *action, gpointer data);
static void                      onreadLabelsMenu(GtkAction *action, gpointer data);
static void                      onselectGapsMenu(GtkAction *action, gpointer data);
static void                      onhideMenu(GtkAction *action, gpointer data);
static void                      onunhideMenu(GtkAction *action, gpointer data);

static void                      onToggleSchemeType(GtkRadioAction *action, GtkRadioAction *current, gpointer data);
static void                      onToggleColorScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data);
static void                      onToggleSortOrder(GtkRadioAction *action, GtkRadioAction *current, gpointer data);

static void                      ontogglePaletteMenu(GtkAction *action, gpointer data);
static void                      ontoggleColorByResIdMenu(GtkAction *action, gpointer data);
static void                      oncolorByResIdMenu(GtkAction *action, gpointer data);
static void                      onsaveColorSchemeMenu(GtkAction *action, gpointer data);
static void                      onloadColorSchemeMenu(GtkAction *action, gpointer data);
static void                      onignoreGapsMenu(GtkAction *action, gpointer data);
static void                      onprintColorsMenu(GtkAction *action, gpointer data);
static void                      onexcludeHighlightedMenu(GtkAction *action, gpointer data);
static void                      ondisplayColorsMenu(GtkAction *action, gpointer data);
static void                      onlowercaseMenu(GtkAction *action, gpointer data);
static void                      oneditColorSchemeMenu(GtkAction *action, gpointer data);

static void                      onSaveTreeMenu(GtkAction *action, gpointer data);
static void                      onFindOrthogsMenu(GtkAction *action, gpointer data);
static void                      onShowOrgsMenu(GtkAction *action, gpointer data);

static void                      onZoomInMenu(GtkAction *action, gpointer data);
static void                      onZoomOutMenu(GtkAction *action, gpointer data);
static void                      onSetFontMenu(GtkAction *action, gpointer data);

static void                      showHelpDialog();
static void                      showAboutDialog(GtkWidget *parent);
static void                      showFindDialog(BelvuContext *bc, GtkWidget *window);
static void                      showWrapDialog(BelvuContext *bc, GtkWidget *belvuWindow);
static void                      createWrapWindow(GtkWidget *belvuWindow, const int linelen, const gchar *title);
static void                      widgetGetDrawing(GtkWidget *window, gpointer data);
static void	                 showMakeNonRedundantDialog(GtkWidget *belvuWindow);
static void                      showRemoveOutliersDialog(GtkWidget *belvuWindow);
static void                      showRemoveByScoreDialog(GtkWidget *belvuWindow);
static void                      startRemovingSequences(GtkWidget *belvuWindow);
static void                      endRemovingSequences(GtkWidget *belvuWindow);

static void                      showRemoveGappySeqsDialog(GtkWidget *belvuWindow);
static void                      showRemoveColumnsDialog(GtkWidget *belvuWindow);
static void                      showRemoveColumnsCutoffDialog(GtkWidget *belvuWindow);
static void                      showRemoveGappyColumnsDialog(GtkWidget *belvuWindow);

static void                      showColorByResIdDialog(GtkWidget *belvuWindow);
static void                      showEditResidueColorsDialog(GtkWidget *belvuWindow, const gboolean bringToFront);
static void                      showEditConsColorsDialog(GtkWidget *belvuWindow, const gboolean bringToFront);
static void                      saveOrResetConsColors(BelvuContext *bc, const gboolean save);
static void                      showSelectGapCharDialog(GtkWidget *belvuWindow);

static GtkWidget*                createPromptDialog(GtkWidget *window, const char *defaultResult, const char *title, const char *text1, const char *text2, GtkWidget **entry);

static gboolean                  saveAlignmentPrompt(GtkWidget *window, BelvuContext *bc);
static void                      showFontDialog(BelvuContext *bc, GtkWidget *window);

static const char*               saveFasta(BelvuContext *bc, GtkWidget *parent);
static const char*               saveMul(BelvuContext *bc, GtkWidget *parent);
static const char*               saveMsf(BelvuContext *bc, GtkWidget *parent);
static void                      showSaveAsDialog(BelvuContext *bc, GtkWidget *window);
static gboolean                  saveAlignment(BelvuContext *bc, GtkWidget *window);

static BelvuWindowProperties*    belvuWindowGetProperties(GtkWidget *widget);
static GenericWindowProperties*  windowGetProperties(GtkWidget *widget);
static BelvuContext*             windowGetContext(GtkWidget *window);
static void                      createOrganismWindow(BelvuContext *bc);
static void                      onDestroyBelvuWindow(GtkWidget *belvuWindow);


/***********************************************************
 *                      Menus and Toolbar                  *
 ***********************************************************/

#define rmPickedStr        "Remove highlighted line"
#define rmPickedDesc       "Remove highlighted line"
#define rmManyStr          "Remove many sequences..."
#define rmManyDesc         "Remove many sequences"
#define rmGappySeqsStr     "Remove gappy sequences..."
#define rmGappySeqsDesc    "Remove gappy sequences"
#define rmPartialSeqsStr   "Remove partial sequences"
#define rmPartialSeqsDesc  "Remove partial sequences"
#define rmRedundantStr     "Make non-redundant..."
#define rmRedundantDesc    "Remove sequences that are more than a given percentage identical"
#define rmOutliersStr      "Remove outliers..."
#define rmOutliersDesc     "Remove sequences that are less than a given percentage identical"
#define rmScoreStr         "Remove sequences by score..."
#define rmScoreDesc        "Remove sequences below a given score"
#define rmColumnPromptStr  "Remove columns..."
#define rmColumnPromptDesc "Remove specific columns"
#define rmColumnLeftStr    "<- Remove columns left of selection (inclusive)"
#define rmColumnLeftDesc   "Remove columns to the left of the currently-selected column (inclusive)"
#define rmColumnRightStr   "Remove columns right of selection (inclusive) ->"
#define rmColumnRightDesc  "Remove columns to the right of the currently-selected column (inclusive) -> "
#define rmColumnCutoffStr  "Remove columns by conservation..."
#define rmColumnCutoffDesc "Remove columns with conservation between specific values"
#define rmGappyColumnsStr  "Remove gappy columns..."
#define rmGappyColumnsDesc "Remove columns with more than a specified percentage of gaps"
#define readLabelsStr      "Read labels of highlighted sequence and spread them"
#define readLabelsDesc     "Read labels of highlighted sequence and spread them"
#define selectGapsStr      "Select gap character..."
#define selectGapsDesc     "Select the character to use for displaying gaps"
#define hideStr            "Hide highlighted line"
#define hideDesc           "Hide the currently-highlighted line"
#define unhideStr          "Unhide all hidden lines"
#define unhideDesc         "Unhide all hidden lines"

#define togglePaletteStr      "Toggle color schemes"
#define togglePaletteDesc     "Toggle between conservation and residue color schemes"
#define colorByResIdStr       "Set %ID threshold..."
#define colorByResIdDesc      "Set the threshold above which to color residues"
#define saveColorSchemeStr    "Save colour scheme..."
#define saveColorSchemeDesc   "Save current colour scheme"
#define loadColorSchemeStr    "Load colour scheme..."
#define loadColorSchemeDesc   "Read colour scheme from file"
#define editColorSchemeStr    "Edit current colour scheme..."
#define editColorSchemeDesc   "Edit the current colour scheme"

#define SaveTreeStr           "Save Tree"
#define SaveTreeDesc          "Save Tree in New Hampshire format"
#define FindOrthogsStr        "Find putative orthologs"
#define FindOrthogsDesc       "Find putative orthologs"
#define ShowOrgsStr           "Show organisms"
#define ShowOrgsDesc          "Show current organisms"
#define PlotOptsStr           "Plot settings..."
#define PlotOptsDesc          "Edit conservation-plot settings"

#define autoRmEmptyColumnsStr  "Automatically remove empty columns"
#define autoRmEmptyColumnsDesc "Automatically remove columns that are 100% gaps after sequence deletions"
#define excludeHighlightedStr  "Exclude highlighted from calculations"
#define excludeHighlightedDesc "Exclude highlighted from calculations"
#define lowercaseStr           "Highlight lowercase characters"
#define lowercaseDesc          "Highlight lowercase characters"

#define colorSimStr            "By average similarity by Blosum62"
#define colorIdStr             "By percent identity"
#define colorIdSimStr          "By percent identity + Blosum62"
#define displayColorsStr        "Display colors (faster without)"
#define printColorsStr         "Use gray shades (for printing)"
#define ignoreGapsStr          "Ignore gaps in conservation calculation"
#define thresholdStr           "Only colour residues above %ID threshold"

#define ConsPlotStr            "Show conservation p_lot"
#define ConsPlotDesc           "Plot conservation profile"
#define WrapStr                "_Wrap for printing..."
#define WrapDesc               "Wrap alignments for printing"
#define OutputStr              "_Output score/coords"
#define OutputDesc             "Output current alignment's score and coords"
#define CompareStr             "Compare all and output identities"
#define CompareDesc            "Compage all sequences against all others and output their identities"
#define FetchWWWStr            "Fetch sequences via WWW"
#define FetchWWWDesc           "Fetch sequences via WWW"


/* Define the menu actions for standard menu entries */
static const GtkActionEntry menuEntries[] = {
  { "FileMenuAction",  NULL, "_File"},
  { "EditMenuAction",  NULL, "_Edit"},
  { "ColorMenuAction", NULL, "_Color"},
  { "SettingsMenuAction", NULL, "_Settings"},
  { "SortMenuAction",  NULL, "S_ort"},
  { "HelpMenuAction",  NULL, "_Help"},

  { "Close",	           GTK_STOCK_CLOSE,      "_Close",             "<control>W",        "Close",                 G_CALLBACK(onCloseMenu)},
  { "Quit", 	           GTK_STOCK_QUIT,       "_Quit",              "<control>Q",        "Quit  Ctrl+Q",          G_CALLBACK(onQuitMenu)},
  { "Help",	           GTK_STOCK_HELP,       "_Help",              "<control>H",        "Display help  Ctrl+H",  G_CALLBACK(onHelpMenu)},
  { "About",	           GTK_STOCK_ABOUT,      "A_bout",             NULL,                "About",                 G_CALLBACK(onAboutMenu)},
  { "Print",	           GTK_STOCK_PRINT,      "_Print...",          "<control>P",        "Print  Ctrl+P",         G_CALLBACK(onPrintMenu)},
  { "Find",                GTK_STOCK_FIND,       "_Find...",           "<control>F",        "Find  Ctrl+F",          G_CALLBACK(onFindMenu)},
  { "Wrap", 	           NULL,                 WrapStr,              NULL,                WrapDesc,                G_CALLBACK(onWrapMenu)},
  { "ShowTree",	           NULL,                 "Show _tree",         NULL,                "Show tree",             G_CALLBACK(onShowTreeMenu)},
  { "RecalcTree",          NULL,                 "Recalculate tree",   NULL,                "Recalculate tree (e.g. after alignment has changed or to reset after swapping nodes)",      G_CALLBACK(onRecalcTreeMenu)},
  { "TreeOpts",	           GTK_STOCK_PREFERENCES,"Tree settings...",   NULL,                "Edit tree settings",    G_CALLBACK(onTreeOptsMenu)},
  { "ConsPlot",	           NULL,                 ConsPlotStr,          NULL,                ConsPlotDesc,            G_CALLBACK(onConsPlotMenu)},
  { "Save",                GTK_STOCK_SAVE,       "_Save",              "<control>S",        "Save alignment",        G_CALLBACK(onSaveMenu)},
  { "SaveAs",	           GTK_STOCK_SAVE_AS,    "Save _as...",        "<shift><control>S", "Save alignment as",     G_CALLBACK(onSaveAsMenu)},
  { "Output",	           NULL,                 OutputStr,            NULL,                OutputDesc,              G_CALLBACK(onOutputMenu)},
  { "Compare",	           NULL,                 CompareStr,           NULL,                CompareDesc,             G_CALLBACK(onCompareMenu)},
  { "CleanUp",	           GTK_STOCK_CLEAR,      "Clean _up windows",  NULL,                "Clean up windows",      G_CALLBACK(onCleanUpMenu)},

  {"rmPicked",             NULL,                 rmPickedStr,          NULL,                rmPickedDesc,            G_CALLBACK(onrmPickedMenu)},
  {"rmGappySeqs",          NULL,                 rmGappySeqsStr,       NULL,                rmGappySeqsDesc,         G_CALLBACK(onrmGappySeqsMenu)},
  {"rmPartialSeqs",        NULL,                 rmPartialSeqsStr,     "<control>T" ,       rmPartialSeqsDesc,       G_CALLBACK(onrmPartialSeqsMenu)},
  {"rmRedundant",          NULL,                 rmRedundantStr,       "<control>R",        rmRedundantDesc,         G_CALLBACK(onrmRedundantMenu)},
  {"rmOutliers",           NULL,                 rmOutliersStr,        NULL,                rmOutliersDesc,          G_CALLBACK(onrmOutliersMenu)},
  {"rmScore",              NULL,                 rmScoreStr,           NULL,                rmScoreDesc,             G_CALLBACK(onrmScoreMenu)},
  {"rmColumnPrompt",       NULL,                 rmColumnPromptStr,    NULL,                rmColumnPromptDesc,      G_CALLBACK(onrmColumnPromptMenu)},
  {"rmColumnLeft",         GTK_STOCK_GO_BACK,    rmColumnLeftStr,      NULL,                rmColumnLeftDesc,        G_CALLBACK(onrmColumnLeftMenu)},
  {"rmColumnRight",        GTK_STOCK_GO_FORWARD, rmColumnRightStr,     NULL,                rmColumnRightDesc,       G_CALLBACK(onrmColumnRightMenu)},
  {"rmColumnCutoff",       NULL,                 rmColumnCutoffStr,    NULL,                rmColumnCutoffDesc,      G_CALLBACK(onrmColumnCutoffMenu)},
  {"rmGappyColumns",       NULL,                 rmGappyColumnsStr,    NULL,                rmGappyColumnsDesc,      G_CALLBACK(onrmGappyColumnsMenu)},
  {"readLabels",           NULL,                 readLabelsStr,        NULL,                readLabelsDesc,          G_CALLBACK(onreadLabelsMenu)},
  {"selectGaps",           NULL,                 selectGapsStr,        NULL,                selectGapsDesc,          G_CALLBACK(onselectGapsMenu)},
  {"hide",                 NULL,                 hideStr,              NULL,                hideDesc,                G_CALLBACK(onhideMenu)},
  {"unhide",               NULL,                 unhideStr,            NULL,                unhideDesc,              G_CALLBACK(onunhideMenu)},

  {"togglePalette",        NULL,                 togglePaletteStr,     "T",                 togglePaletteDesc,       G_CALLBACK(ontogglePaletteMenu)},
  {"colorByResId",         NULL,                 colorByResIdStr,      NULL,                colorByResIdDesc,        G_CALLBACK(oncolorByResIdMenu)},
  {"saveColorScheme",      NULL,                 saveColorSchemeStr,   NULL,                saveColorSchemeDesc,     G_CALLBACK(onsaveColorSchemeMenu)},
  {"loadColorScheme",      NULL,                 loadColorSchemeStr,   NULL,                loadColorSchemeDesc,     G_CALLBACK(onloadColorSchemeMenu)},
  {"editColorScheme",    GTK_STOCK_SELECT_COLOR, editColorSchemeStr,   NULL,                editColorSchemeDesc,     G_CALLBACK(oneditColorSchemeMenu)},

  {"SaveTree",             GTK_STOCK_FLOPPY,     SaveTreeStr,          NULL,                SaveTreeDesc,            G_CALLBACK(onSaveTreeMenu)},
  {"ShowOrgs",             NULL,                 ShowOrgsStr,          NULL,                ShowOrgsDesc,            G_CALLBACK(onShowOrgsMenu)},
  {"PlotOpts",             GTK_STOCK_PROPERTIES, PlotOptsStr,          NULL,                PlotOptsDesc,            G_CALLBACK(onPlotOptsMenu)},

  {"ZoomIn",              GTK_STOCK_ZOOM_IN,    "Zoom in",            NULL,                 "Zoom in  =",               G_CALLBACK(onZoomInMenu)},
  {"ZoomOut",             GTK_STOCK_ZOOM_OUT,   "Zoom out",           NULL,                 "Zoom out  -",              G_CALLBACK(onZoomOutMenu)},
  {"SetFont",             NULL,                 "Set font size",      NULL,                 "Set font size",         G_CALLBACK(onSetFontMenu)}
};

/* Define the menu actions for toggle menu entries */
static const GtkToggleActionEntry toggleMenuEntries[] = {
  {"FetchWWW",               NULL, FetchWWWStr,                          NULL,                FetchWWWDesc,            G_CALLBACK(onFetchWWWMenu),           FALSE},
  {"rmMany",     GTK_STOCK_DELETE, rmManyStr,                            NULL,                rmManyDesc,              G_CALLBACK(onRemoveSeqsMenu),         FALSE},
  {"FindOrthogs",            NULL, FindOrthogsStr,                       NULL,                FindOrthogsDesc,         G_CALLBACK(onFindOrthogsMenu),        FALSE},
  {"autoRmEmptyColumns",     NULL, autoRmEmptyColumnsStr,                NULL,                autoRmEmptyColumnsDesc,  G_CALLBACK(onAutoRmEmptyColumnsMenu), TRUE}, 
  {"toggleColorByResId",     NULL, thresholdStr,                         NULL,                thresholdStr,            G_CALLBACK(ontoggleColorByResIdMenu), FALSE},
  {"ignoreGaps",             NULL, ignoreGapsStr,                        NULL,                ignoreGapsStr,           G_CALLBACK(onignoreGapsMenu),         FALSE},
  {"printColors",            NULL, printColorsStr,                       NULL,                printColorsStr,          G_CALLBACK(onprintColorsMenu),        FALSE},
  {"excludeHighlighted",     NULL, excludeHighlightedStr,                NULL,                excludeHighlightedDesc,  G_CALLBACK(onexcludeHighlightedMenu), FALSE},
  {"displayColors",          NULL, displayColorsStr,                     NULL,                displayColorsStr,        G_CALLBACK(ondisplayColorsMenu),      TRUE},
  {"lowercase",              NULL, lowercaseStr,                         NULL,                lowercaseDesc,           G_CALLBACK(onlowercaseMenu),          FALSE}
};


/* Define the menu actions for radio-button menu entries */
static const GtkRadioActionEntry schemeMenuEntries[] = {
  {"ColorByResidue",       NULL, "Color by _residue",                       NULL, "Color by residue",                  BELVU_SCHEME_TYPE_RESIDUE},
  {"ColorByCons",          NULL, "Color by _conservation",                  NULL, "Color by conservation",             BELVU_SCHEME_TYPE_CONS}
};

static const GtkRadioActionEntry colorSchemeMenuEntries[] = {
  {"colorSchemeStandard",  NULL, "By residue: Erik's",                            NULL, "Erik's",                            BELVU_SCHEME_ERIK},
  {"colorSchemeGibson",    NULL, "By residue: Toby's",                            NULL, "Toby's",                            BELVU_SCHEME_GIBSON},
  {"colorSchemeCGP",       NULL, "By residue: Cys/Gly/Pro",                       NULL, "Cys/Gly/Pro",                       BELVU_SCHEME_CGP},
  {"colorSchemeCGPH",       NULL, "By residue: Cys/Gly/Pro/His",                  NULL, "Cys/Gly/Pro/His",                   BELVU_SCHEME_CGPH},
  {"colorSchemeEmpty",     NULL, "By residue: Clean slate",                       NULL, "Clean slate",                       BELVU_SCHEME_NONE},
  {"colorSchemeCustom",    NULL, "By residue: Custom",                            NULL, "Custom",                            BELVU_SCHEME_CUSTOM},
  {"colorSim",             NULL, colorSimStr,                         NULL, colorSimStr,                         BELVU_SCHEME_BLOSUM},
  {"colorId",              NULL, colorIdStr,                          NULL, colorIdStr,                          BELVU_SCHEME_ID},
  {"colorIdSim",           NULL, colorIdSimStr,                       NULL, colorIdSimStr,                       BELVU_SCHEME_ID_BLOSUM}
};

static const GtkRadioActionEntry sortMenuEntries[] = {
  {"unsorted",             NULL,                     "unsorted",                           NULL, "Unsorted",                          BELVU_UNSORTED},
  {"consSort",             NULL,                     "by conservation",                    NULL, "Sort by conservation",              BELVU_SORT_CONS},
  {"scoreSort",            NULL,                     "by score",                           NULL, "Sort by score",                     BELVU_SORT_SCORE},
  {"alphaSort",            GTK_STOCK_SORT_ASCENDING, "alphabetically",                     NULL, "Sort alphabetically",               BELVU_SORT_ALPHA},
  {"organismSort",         NULL,                     "by swissprot organism",              NULL, "Sort by swissprot organism",        BELVU_SORT_ORGANISM},
  {"treeSort",             NULL,                     "by tree order",                      NULL, "Sort by tree order",                BELVU_SORT_TREE},
  {"simSort",              NULL,                     "by similarity to selected sequence", NULL, "Sort by similarity to highlighted sequence", BELVU_SORT_SIM},
  {"idSort",               NULL,                     "by identity to selected sequence",   NULL, "Sort by identity to highlighted sequence",   BELVU_SORT_ID}
};


/* Define the menu layout */
static const char standardMenuDescription[] =
"<ui>"
/* ACCELERATORS */
"  <accelerator action='togglePalette'/>"
"  <accelerator action='Find'/>"

/* MAIN MENU BAR */
"  <menubar name='MenuBar' accelerators='true'>"
     /* File menu */
"    <menu action='FileMenuAction'>"
"      <menuitem action='Quit'/>"
"      <menuitem action='Wrap'/>"
"      <menuitem action='Print'/>"
"      <separator/>"
"      <menuitem action='ShowTree'/>"
"      <menuitem action='TreeOpts'/>"
"      <menuitem action='RecalcTree'/>"
"      <separator/>"
"      <menuitem action='ConsPlot'/>"
"      <separator/>"
"      <menuitem action='Save'/>"
"      <menuitem action='SaveAs'/>"
"      <menuitem action='Output'/>"
"      <separator/>"
"      <menuitem action='FetchWWW'/>"
"      <menuitem action='Compare'/>"
"      <menuitem action='CleanUp'/>"
"    </menu>"
    /* Edit menu */
"    <menu action='EditMenuAction'>"
"      <menuitem action='rmPicked'/>"
"      <menuitem action='rmMany'/>"
"      <menuitem action='rmGappySeqs'/>"
"      <menuitem action='rmPartialSeqs'/>"
"      <menuitem action='rmRedundant'/>"
"      <menuitem action='rmOutliers'/>"
"      <menuitem action='rmScore'/>"
"      <separator/>"
"      <menuitem action='rmColumnPrompt'/>"
"      <menuitem action='rmColumnLeft'/>"
"      <menuitem action='rmColumnRight'/>"
"      <menuitem action='rmColumnCutoff'/>"
"      <menuitem action='rmGappyColumns'/>"
"      <menuitem action='autoRmEmptyColumns'/>"
"      <separator/>"
"      <menuitem action='readLabels'/>"
"      <menuitem action='selectGaps'/>"
"      <menuitem action='hide'/>"
"      <menuitem action='unhide'/>"
"    </menu>"
     /* Color schemes menu */
"    <menu action='ColorMenuAction'>"
"      <menuitem action='colorSchemeStandard'/>"
"      <menuitem action='colorSchemeGibson'/>"
"      <menuitem action='colorSchemeCGP'/>"
"      <menuitem action='colorSchemeCGPH'/>"
"      <menuitem action='colorSchemeEmpty'/>"
"      <menuitem action='colorSchemeCustom'/>"
"      <separator/>"
"      <menuitem action='colorSim'/>"
"      <menuitem action='colorId'/>"
"      <menuitem action='colorIdSim'/>"
"      <separator/>"
"      <menuitem action='editColorScheme'/>"
"      <menuitem action='saveColorScheme'/>"
"      <menuitem action='loadColorScheme'/>"
"    </menu>"
     /* Settings menu */
"    <menu action='SettingsMenuAction'>"
"      <menuitem action='toggleColorByResId'/>"
"      <menuitem action='colorByResId'/>"
"      <menuitem action='ignoreGaps'/>"
"      <separator/>"
"      <menuitem action='excludeHighlighted'/>"
"      <menuitem action='printColors'/>"
"      <menuitem action='displayColors'/>"
"      <menuitem action='lowercase'/>"
"    </menu>"
     /* Sort menu */
"    <menu action='SortMenuAction'>"
"      <menuitem action='scoreSort'/>"
"      <menuitem action='alphaSort'/>"
"      <menuitem action='organismSort'/>"
"      <menuitem action='treeSort'/>"
"      <menuitem action='simSort'/>"
"      <menuitem action='idSort'/>"
"    </menu>"
/* Help menu */
"    <menu action='HelpMenuAction'>"
"      <menuitem action='Help'/>"
"      <menuitem action='About'/>"
"    </menu>"
"  </menubar>"

/* CONTEXT MENUS */

/* Main context menu */
"  <popup name='ContextMenu' accelerators='true'>"
"    <menuitem action='Quit'/>"
"    <menuitem action='Wrap'/>"
"    <menuitem action='Print'/>"
"    <separator/>"
"    <menuitem action='ShowTree'/>"
"    <menuitem action='TreeOpts'/>"
"    <menuitem action='RecalcTree'/>"
"    <separator/>"
"    <menuitem action='ConsPlot'/>"
"    <separator/>"
"    <menuitem action='Save'/>"
"    <menuitem action='SaveAs'/>"
"    <menuitem action='Output'/>"
"    <separator/>"
"    <menuitem action='FetchWWW'/>"
"    <menuitem action='Compare'/>"
"    <menuitem action='CleanUp'/>"
"  </popup>"
/* Wrapped-alignments window context menu */
"  <popup name='WrapContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <menuitem action='Wrap'/>"
"  </popup>"
/* Organisms window context menu */
"  <popup name='OrgsContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"  </popup>"
/* Tree context menu */
"  <popup name='TreeContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <menuitem action='SaveTree'/>"
"    <separator/>"
"    <menuitem action='TreeOpts'/>"
"    <menuitem action='RecalcTree'/>"
"    <separator/>"
"    <menuitem action='FindOrthogs'/>"
"    <menuitem action='ShowOrgs'/>"
"  </popup>"
/* Conservation-plot context menu */
"  <popup name='PlotContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Print'/>"
"    <separator/>"
"    <menuitem action='PlotOpts'/>"
"  </popup>"

/* TOOLBAR */
"  <toolbar name='Toolbar'>"
"    <toolitem action='Help'/>"
"    <toolitem action='rmMany'/>"
"    <toolitem action='editColorScheme'/>"
"    <toolitem action='alphaSort'/>"
"    <separator/>"
"    <toolitem action='ZoomIn'/>"
"    <toolitem action='ZoomOut'/>"
"    <toolitem action='Find'/>"
"  </toolbar>"
"</ui>";



/* Utility function to create the UI manager for the menus */
GtkUIManager* createUiManager(GtkWidget *window, 
                              BelvuContext *bc, 
                              GtkActionGroup **actionGroupOut)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  gtk_action_group_add_actions(action_group, menuEntries, G_N_ELEMENTS(menuEntries), window);
  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS(toggleMenuEntries), window);

  gtk_action_group_add_radio_actions(action_group, schemeMenuEntries, G_N_ELEMENTS(schemeMenuEntries), bc->schemeType, G_CALLBACK(onToggleSchemeType), window);
  gtk_action_group_add_radio_actions(action_group, colorSchemeMenuEntries, G_N_ELEMENTS(colorSchemeMenuEntries), BELVU_SCHEME_BLOSUM, G_CALLBACK(onToggleColorScheme), window);
  gtk_action_group_add_radio_actions(action_group, sortMenuEntries, G_N_ELEMENTS(sortMenuEntries), BELVU_SORT_CONS, G_CALLBACK(onToggleSortOrder), window);

  greyOutInvalidActionsForGroup(bc, action_group);

  GtkUIManager *ui_manager = gtk_ui_manager_new();
  gtk_ui_manager_insert_action_group(ui_manager, action_group, 0);
  gtk_ui_manager_set_add_tearoffs(ui_manager, TRUE);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group(ui_manager);
  gtk_window_add_accel_group(GTK_WINDOW(window), accel_group);
  
  if (actionGroupOut)
    *actionGroupOut = action_group;
  
  return ui_manager;
}


/* Create a menu */
GtkWidget* createBelvuMenu(GtkWidget *window, 
                           const char *path, 
                           GtkUIManager *ui_manager)
{
  GError *error = NULL;
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, standardMenuDescription, -1, &error))
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  GtkWidget *menu = gtk_ui_manager_get_widget (ui_manager, path);
  
  return menu;
}


/* The following functions implement the menu actions */

/* FILE MENU ACTIONS */
static void onCloseMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  if (window == bc->consPlot) /* conservation plot is persistent, so just hide it rather than closing */
    gtk_widget_hide_all(window);
  else
    gtk_widget_destroy(window);
}

static void onQuitMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  if (stringsEqual(gtk_widget_get_name(window), MAIN_BELVU_WINDOW_NAME, TRUE))
    {
      /* If this is the main window, just call its destructor function, which
       * handles save-checking and quitting the application itself (and we don't
       * want to do that twice). */
      onDestroyBelvuWindow(window);
    }
  else
    {
      gboolean quit = TRUE;
      
      /* Check if the alignment has been save and if not give the option to cancel */
      if (!bc->saved)
        quit = saveAlignmentPrompt(window, bc);
      
      if (quit)
        gtk_main_quit();
    }
}

static void onHelpMenu(GtkAction *action, gpointer data)
{
  showHelpDialog();
}

static void onAboutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showAboutDialog(belvuWindow);
}

/* For some widget the print function prints the visible window as it is, for
 * others it prints only the cached drawables. This returns true if the latter. */
static gboolean printCachedDrawablesOnly(GtkWidget *widget)
{
  const char *name = gtk_widget_get_name(widget);
  
  return (stringsEqual(name, WRAPPED_BELVU_WINDOW_NAME, TRUE) ||
          stringsEqual(name, BELVU_TREE_WINDOW_NAME, TRUE) ||
          stringsEqual(name, BELVU_ORGS_WINDOW_NAME, TRUE));
}

static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);

  static GtkPageSetup *pageSetup = NULL;
  static GtkPrintSettings *printSettings = NULL;
  
  if (!pageSetup)
    {
      pageSetup = gtk_page_setup_new();
      gtk_page_setup_set_orientation(pageSetup, GTK_PAGE_ORIENTATION_PORTRAIT);
    }
  
  if (!printSettings)
    {
      printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(printSettings, GTK_PAGE_ORIENTATION_PORTRAIT);
      gtk_print_settings_set_quality(printSettings, GTK_PRINT_QUALITY_HIGH);
      gtk_print_settings_set_resolution(printSettings, DEFAULT_PRINT_RESOLUTION);
    }
  
  /* If we're just printing the cached drawable, get the actual drawing
   * area widget that should be drawn. (Otherwise it gets clipped to the 
   * size of the container widget) */
  const gboolean printCachedOnly = printCachedDrawablesOnly(window);
  PrintScaleType scaleType = PRINT_FIT_BOTH;

  GtkWidget *widgetToPrint = NULL;
  
  if (printCachedOnly)
    {
      widgetGetDrawing(window, &widgetToPrint);
      scaleType = PRINT_FIT_WIDTH; /* might be very tall, so allow multiple pages height-wise */
    }
  else
    {
      widgetToPrint = window;
    }
  
  blxPrintWidget(widgetToPrint, NULL, GTK_WINDOW(window), &printSettings, &pageSetup, NULL, printCachedOnly, scaleType);
}

static void onFindMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  showFindDialog(bc, window);
}

static void onWrapMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  showWrapDialog(bc, bc->belvuWindow);
}

static void onShowTreeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (!properties->bc->belvuTree)
    {
      /* If the tree exists, create a window from it. Otherwise we need to
       * create it before we show it. */
      if (properties->bc->mainTree && properties->bc->mainTree->head)
        createBelvuTreeWindow(properties->bc, properties->bc->mainTree, TRUE);
      else
        createAndShowBelvuTree(properties->bc, TRUE);
    }

  if (properties->bc->belvuTree)
    gtk_window_present(GTK_WINDOW(properties->bc->belvuTree));
}

/* Utility to extract the context from any toplevel window type (i.e. the main
 * window, the tree, or a wrapped-alignment window). */
static BelvuContext* windowGetContext(GtkWidget *window)
{
  BelvuContext *bc = NULL;
  const char *name = gtk_widget_get_name(window);
  
  if (stringsEqual(name, MAIN_BELVU_WINDOW_NAME, TRUE))
    {
      BelvuWindowProperties *properties = belvuWindowGetProperties(window);
      bc = properties->bc;
    }
  else if (stringsEqual(name, BELVU_TREE_WINDOW_NAME, TRUE))
    {
      bc = belvuTreeGetContext(window);
    }
  else if (stringsEqual(name, BELVU_CONS_PLOT_WINDOW_NAME, TRUE))
    {
      bc = consPlotGetContext(window);
    }
  else /* generic windows */
    {
      GenericWindowProperties *properties = windowGetProperties(window);
      bc = properties->bc;
    }
  
  return bc;
}

static void onRecalcTreeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  if (bc->belvuTree)
    {
      /* Update the tree window */
      belvuTreeRemakeTree(bc->belvuTree);
    }
  else
    {
      /* No tree window, but make/re-make the underlying tree structure */
      separateMarkupLines(bc);
      Tree *tree = treeMake(bc, FALSE, TRUE);
      reInsertMarkupLines(bc);
      
      belvuContextSetTree(bc, &tree);
      onTreeOrderChanged(bc);
    }
}

static void onTreeOptsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  showTreeSettingsDialog(window, bc);
}

static void onConsPlotMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  if (bc->consPlot)
    {
      gtk_widget_show_all(bc->consPlot);
      gtk_window_present(GTK_WINDOW(bc->consPlot));
    }
  else
    {
      createConsPlot(bc);
    }
}


static void onSaveMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  saveAlignment(bc, window);
}

static void onSaveAsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  showSaveAsDialog(bc, window);
}

static void onOutputMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(belvuWindow);
  
  if (!bc->selectedAln) 
    {
      g_critical("Please select a sequence first.\n");
    }
  else
    {
      g_message("%.1f %s/%d-%d\n", 
                bc->selectedAln->score,
                bc->selectedAln->name,
                bc->selectedAln->start,
                bc->selectedAln->end);
    
      fflush(stdout);
    }
}

static void onCompareMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  listIdentity(bc);
}

static void onFetchWWWMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  bc->useWWWFetch = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
}

/* This just calls gtk_widget_destroy but accepts gpointer arguments so
 * that we can call it from a foreach fucntion. */
static void destroyWidget(gpointer widget, gpointer data)
{
  gtk_widget_destroy(GTK_WIDGET(widget));
}

static void onCleanUpMenu(GtkAction *action, gpointer data)
{
  /* Close all windows that were spawned from the main window */
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  g_slist_foreach(properties->bc->spawnedWindows, destroyWidget, NULL);
}

/* EDIT MENU ACTIONS */
static void onrmPickedMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  removeSelectedSequence(properties->bc, properties->bc->belvuAlignment);
}

static void onrmGappySeqsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showRemoveGappySeqsDialog(belvuWindow);
}

static void onrmPartialSeqsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  removePartialSeqs(properties->bc, properties->bc->belvuAlignment);
}

static void onrmRedundantMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showMakeNonRedundantDialog(belvuWindow);
}

static void onrmOutliersMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showRemoveOutliersDialog(belvuWindow);
}

static void onrmScoreMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showRemoveByScoreDialog(belvuWindow);
}

static void onrmColumnPromptMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showRemoveColumnsDialog(belvuWindow);
}

static void onrmColumnLeftMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (properties->bc->selectedCol > 0)
    {
      rmColumn(properties->bc, 1, properties->bc->selectedCol);
      rmFinaliseColumnRemoval(properties->bc);

      updateOnAlignmentLenChanged(properties->bc->belvuAlignment);

      properties->bc->selectedCol = 0; /* cancel selection, because this col is deleted now */
      properties->bc->highlightedCol = 0; /* cancel selection, because this col is deleted now */
      
      onColSelectionChanged(properties->bc);
    }
  else
    {
      g_critical("Please select a column first.\n\nMiddle-click with the mouse to select a column.\n");
    }
}

static void onrmColumnRightMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (properties->bc->selectedCol > 0)
    {
      rmColumn(properties->bc, properties->bc->selectedCol, properties->bc->maxLen);
      rmFinaliseColumnRemoval(properties->bc);

      updateOnAlignmentLenChanged(properties->bc->belvuAlignment);

      properties->bc->selectedCol = 0; /* cancel selection, because this col is deleted now */
      properties->bc->highlightedCol = 0; /* cancel selection, because this col is deleted now */
      onColSelectionChanged(properties->bc);
    }
  else
    {
      g_critical("Please select a column first.\n\nMiddle-click with the mouse to select a column.\n");
    }
}

/* Remove columns based on a conservation-cutoff */
static void onrmColumnCutoffMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (!colorByConservation(properties->bc)) 
    {
      g_critical("Please select a conservation coloring scheme from the Color menu first.\n");
      return;
    }

  showRemoveColumnsCutoffDialog(belvuWindow);
}

/* Remove columns with more than a specified percentage of gaps */
static void onrmGappyColumnsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showRemoveGappyColumnsDialog(belvuWindow);
}

/* Toggle the 'remove-many-sequences' option on or off */
static void onRemoveSeqsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  
  const gboolean optionOn = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
  
  if (optionOn)
    startRemovingSequences(belvuWindow);
  else
    endRemovingSequences(belvuWindow);
}

/* Toggle the 'auto-remove empty columns' option */
static void onAutoRmEmptyColumnsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  const gboolean newVal = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
  
  if (properties->bc->rmEmptyColumnsOn != newVal)
    {
      properties->bc->rmEmptyColumnsOn = newVal;
    }
}

static void onreadLabelsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (!properties->bc->selectedAln) 
    {
      g_critical("Please select a sequence first.\n");
      return;
    }

  char *title = g_strdup_printf("Read labels of %s from file", properties->bc->selectedAln->name);

  const char *filename = getLoadFileName(belvuWindow, properties->bc->dirName, title);
  g_free(title);

  FILE *fil = fopen(filename, "r");

  if (fil)
    {
      readLabels(properties->bc, fil);
      belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
    }
}

static void onselectGapsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showSelectGapCharDialog(belvuWindow);
}

static void onhideMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (!properties->bc->selectedAln)
    {
      g_critical("Please select a sequence to hide.\n");
    }
  else
    {
      properties->bc->selectedAln->hide = TRUE;
      belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
    }
}

static void onunhideMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  /* Reset the 'hide' flag to false in all sequences */
  int i = 0;
  for (i = 0; i < (int)properties->bc->alignArr->len; ++i)
    {
      g_array_index(properties->bc->alignArr, ALN*, i)->hide = FALSE;
    }
  
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}


/***********************************************************
 *                  Save on exit dialog                    *
 ***********************************************************/

/* Pops up a dialog asking the user whether they want to save the alignment
 * or not. Deals with the save, if applicable, and returns true if the user
 * still wishes to quit (or false if the user cancelled). */
static gboolean saveAlignmentPrompt(GtkWidget *widget, BelvuContext *bc)
{
  char *title = g_strdup_printf("%sSave alignment?", belvuGetTitlePrefix(bc));
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title,
                                                  GTK_WINDOW(widget), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_YES, GTK_RESPONSE_YES,        /* yes, save the alignment and exit */
                                                  GTK_STOCK_NO, GTK_RESPONSE_NO,          /* no, don't save (but still exit) */
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,  /* don't save and don't exit */
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_YES);

  /* Put message and icon into an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 0);

  GtkWidget *image = gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start(GTK_BOX(hbox), image, TRUE, TRUE, 0);

  GtkWidget *label = gtk_label_new("Alignment was modified - save ?");
  gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
  gtk_widget_show_all(hbox);

  gint response = gtk_dialog_run(GTK_DIALOG(dialog));
  gboolean quit = FALSE;
  
  if (response == GTK_RESPONSE_YES)
    {
      quit = saveAlignment(bc, widget);
    }
  else if (response == GTK_RESPONSE_NO)
    {
      quit = TRUE;
    }

  gtk_widget_destroy(dialog);
  
  return quit;
}

/***********************************************************
 *                  File menu actions                      *
 ***********************************************************/

static const char* saveFasta(BelvuContext *bc, GtkWidget *parent)
{
  char *title = g_strdup_printf("%s", bc->saveFormat == BELVU_FILE_UNALIGNED_FASTA ? "Save as unaligned Fasta file:" : "Save as aligned Fasta file:");
  const char *filename = getSaveFileName(parent, bc->fileName, bc->dirName, NULL, title);
  g_free(title);
  
  FILE *fil = fopen(filename, "w");
  
  if (fil)
    {
      writeFasta(bc, fil);
    }
  
  return filename;
}


static const char* saveMsf(BelvuContext *bc, GtkWidget *parent)
{
  const char *filename = getSaveFileName(parent, bc->fileName, bc->dirName, NULL, "Save as MSF (/) file:");
  
  FILE *fil = fopen(filename, "w");
  
  if (fil)
    {
      writeMSF(bc, fil);
    }
  
  return filename;
}


static const char* saveMul(BelvuContext *bc, GtkWidget *parent)
{
  const char *filename = getSaveFileName(parent, bc->fileName, bc->dirName, NULL, "Save as Stockholm file:");
  FILE *fil = fopen(filename, "w");
  
  if (fil)
    {
      writeMul(bc, fil);
    }
  
  return filename;
}


/***********************************************************
 *                 Colour menu actions                     *
 ***********************************************************/

/* This function is called when the color scheme has been changed. It performs all 
 * required updates. */
static void onColorSchemeChanged(BelvuWindowProperties *properties)
{
  /* Make sure the correct scheme type is set in the menus */
  switch (properties->bc->schemeType)
    {
      case BELVU_SCHEME_TYPE_RESIDUE:
	setToggleMenuStatus(properties->actionGroup, "ColorByResidue", TRUE);
	break;

      case BELVU_SCHEME_TYPE_CONS:
	setToggleMenuStatus(properties->actionGroup, "ColorByCons", TRUE);
	break;
    
      default:
	g_warning("Program error: unrecognised color scheme type '%d'.\n", properties->bc->schemeType);
	break;
    };
  
  /* Some menu actions are enabled/disabled depending on which scheme type is selected */
  greyOutInvalidActions(properties->bc);
  
  /* Update the display */
  updateSchemeColors(properties->bc);
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}


static void onToggleSchemeType(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  BelvuSchemeType newScheme = (BelvuSchemeType)gtk_radio_action_get_current_value(current);
  
  if (newScheme != properties->bc->schemeType)
    {
      properties->bc->schemeType = newScheme;
      
      /* Toggle the actual scheme to the current default for this scheme type */
      if (newScheme == BELVU_SCHEME_TYPE_RESIDUE)
        setRadioMenuStatus(properties->actionGroup, "colorSchemeStandard", properties->bc->residueScheme);
      else
        setRadioMenuStatus(properties->actionGroup, "colorSchemeStandard", properties->bc->consScheme);
      
      onColorSchemeChanged(properties);
    }
}


static void onToggleColorScheme(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;

  /* Get the new color scheme */
  BelvuColorScheme newScheme = (BelvuColorScheme)gtk_radio_action_get_current_value(current);

  /* Determine what scheme type this puts us in */
  if (newScheme > NUM_RESIDUE_SCHEMES)
    {
      /* The new scheme is a color-by-conservation scheme type */
      if (bc->schemeType != BELVU_SCHEME_TYPE_CONS || bc->consScheme != newScheme)
        {
          bc->schemeType = BELVU_SCHEME_TYPE_CONS;
          bc->consScheme = newScheme;
          onColorSchemeChanged(properties);
        }
    }
  else
    {
      /* The new scheme is a color-by-residue scheme type */
      if (bc->schemeType != BELVU_SCHEME_TYPE_RESIDUE || bc->residueScheme != newScheme)
        {
          bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
          bc->residueScheme = newScheme;
          setResidueSchemeColors(bc);
          onColorSchemeChanged(properties);
        }
    }
}


static void onToggleSortOrder(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  properties->bc->sortType = (BelvuSortType)gtk_radio_action_get_current_value(current);
  doSort(properties->bc, properties->bc->sortType, TRUE);

  /* To do: This is a hack to overcome a bug where the sort order gets messed
   * up when switching to tree-sort from a different sort order after having
   * changed the tree order by swapping nodes. Calling doSort again seems to
   * sort it out, although obviously this is not ideal. */
  if (properties->bc->sortType == BELVU_SORT_TREE)
    doSort(properties->bc, properties->bc->sortType, TRUE);

  centerHighlighted(properties->bc, properties->bc->belvuAlignment);
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}


/* Toggle between color-by-conservation and color-by-residue modes */
static void ontogglePaletteMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  
  if (stringsEqual(gtk_widget_get_name(window), MAIN_BELVU_WINDOW_NAME, TRUE))
    {
      BelvuWindowProperties *properties = belvuWindowGetProperties(window);
      
      if (properties->bc->schemeType == BELVU_SCHEME_TYPE_CONS)
        setToggleMenuStatus(properties->actionGroup, "ColorByResidue", TRUE);
      else
        setToggleMenuStatus(properties->actionGroup, "ColorByCons", TRUE);
    }
}

/* This controls whether the color-by-residue-id option is toggle on or off */
static void ontoggleColorByResIdMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  /* Toggle the flag */
  properties->bc->colorByResIdOn = !properties->bc->colorByResIdOn;
  
  /* Update the color scheme and redraw */
  updateSchemeColors(properties->bc);
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}

/* This pops up a dialog to set the color-by-res-id threshold, and turns the option on */
static void oncolorByResIdMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  showColorByResIdDialog(belvuWindow);  
}

static void onsaveColorSchemeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  const char *filename = getSaveFileName(belvuWindow, properties->bc->fileName, properties->bc->dirName, NULL, "Save colour scheme");
  FILE *fil = fopen(filename, "w");

  if (filename)
    {
      if (properties->bc->dirName) g_free(properties->bc->dirName);
      if (properties->bc->fileName) g_free(properties->bc->fileName);
      properties->bc->dirName = g_path_get_dirname(filename);
      properties->bc->fileName = g_path_get_basename(filename);
    }
  
  saveResidueColorScheme(properties->bc, fil);
}

static void onloadColorSchemeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  const char *filename = getLoadFileName(belvuWindow, properties->bc->dirName, "Read color scheme");
  FILE *fil = fopen(filename, "r");

  if (filename)
    {
      if (properties->bc->dirName) g_free(properties->bc->dirName);
      if (properties->bc->fileName) g_free(properties->bc->fileName);
      properties->bc->dirName = g_path_get_dirname(filename);
      properties->bc->fileName = g_path_get_basename(filename);
    }
  
  readResidueColorScheme(properties->bc, fil, getColorArray(), TRUE);

  setRadioMenuStatus(properties->actionGroup, "colorSchemeStandard", BELVU_SCHEME_CUSTOM);
  setToggleMenuStatus(properties->actionGroup, "ColorByResidue", TRUE);
  onColorSchemeChanged(properties);
}

static void onignoreGapsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  /* Toggle the 'ignore gaps' option */
  properties->bc->ignoreGapsOn = !properties->bc->ignoreGapsOn;
  
  onColorSchemeChanged(properties);
}

static void onprintColorsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  /* Toggle the 'ignore gaps' option */
  properties->bc->printColorsOn = !properties->bc->printColorsOn;
  
  onColorSchemeChanged(properties);
}

static void onexcludeHighlightedMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  const gboolean exclude = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
  
  setExcludeFromConsCalc(properties->bc, exclude);
  
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}

static void ondisplayColorsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  const gboolean newVal = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
  
  if (properties->bc->displayColors != newVal)
    {
      properties->bc->displayColors = newVal;
      belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
    }
}

static void onlowercaseMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  properties->bc->lowercaseOn = !properties->bc->lowercaseOn;
  belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
}

static void oneditColorSchemeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(belvuWindow);
  
  if (colorByConservation(bc))
    showEditConsColorsDialog(belvuWindow, TRUE);
  else
    showEditResidueColorsDialog(belvuWindow, TRUE);
}

static void onSaveTreeMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  const char *filename = getSaveFileName(window, bc->fileName, bc->dirName, NULL, "Save tree in New Hampshire format");
  
  FILE *file = fopen(filename, "w");
  
  if (file)
    {
      saveTreeNH(bc->mainTree, bc->mainTree->head, file);

      /* Add a terminating line and close the file. */
      fprintf(file, ";\n");
      fclose(file);
      
      g_message_info("Tree saved to %s\n", filename);
    }
}

static void onFindOrthogsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  GtkWidget *treeWindow = NULL;
  
  /* If the current window is a tree, use that; otherwise use the main belvu tree */
  if (stringsEqual(gtk_widget_get_name(window), BELVU_TREE_WINDOW_NAME, TRUE))
    treeWindow = window;
  else
    treeWindow = bc->belvuTree;
  
  bc->highlightOrthologs = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));
      
  /* If turning the option on, print the orthologs to stdout */
  if (bc->highlightOrthologs)
    treePrintOrthologs(bc, treeWindow);

  /* Refresh the tree to show or hide the orthologs */
  belvuTreeRedrawAll(bc->belvuTree, NULL);
}

static void onShowOrgsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);
  
  if (bc->organismArr->len < 1)
    {
      g_critical("No organism details found.\n");
    }
  else if (bc->orgsWindow)
    {
      gtk_widget_show_all(bc->orgsWindow);
      gtk_window_present(GTK_WINDOW(bc->orgsWindow));
    }
  else
    {
      createOrganismWindow(bc);
    }
}

static void onZoomInMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  incrementFontSize(bc);
}

static void onZoomOutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  decrementFontSize(bc);
}

static void onSetFontMenu(GtkAction *action, gpointer data)
{
  GtkWidget *window = GTK_WIDGET(data);
  BelvuContext *bc = windowGetContext(window);

  showFontDialog(bc, window);
}


/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuWindowProperties* belvuWindowGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuWindowProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuWindowProperties")) : NULL;
}


/* Does the job of destroying the belvu window */
static void destroyBelvuWindow(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (properties)
    {
      destroyBelvuContext(&properties->bc);
  
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuWindow), "BelvuWindowProperties", NULL);
    }

  gtk_main_quit();  
}


/* Signal handler for when the main belvu window is closed (note that if the
 * alignment is not saved the user may cancel) */
static void onDestroyBelvuWindow(GtkWidget *belvuWindow)
{
  gboolean destroy = TRUE;
  
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (!properties->bc->saved)
    {
      /* The alignment has not been saved - ask the user if they want to save/cancel */
      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      destroy = saveAlignmentPrompt(belvuWindow, properties->bc);
    }
  
  if (destroy)
    {
      destroyBelvuWindow(belvuWindow);
    }
}


/* Create the properties struct and initialise all values. */
static void belvuWindowCreateProperties(GtkWidget *belvuWindow, 
					BelvuContext *bc, 
					GtkWidget *statusBar,
					GtkWidget *feedbackBox,
                                        GtkActionGroup *actionGroup)
{
  if (belvuWindow)
    {
      BelvuWindowProperties *properties = (BelvuWindowProperties*)g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->statusBar = statusBar;
      properties->feedbackBox = feedbackBox;
      properties->actionGroup = actionGroup;
      
      g_object_set_data(G_OBJECT(belvuWindow), "BelvuWindowProperties", properties);
      g_signal_connect(G_OBJECT(belvuWindow), "destroy", G_CALLBACK (onDestroyBelvuWindow), NULL);
    }
}


GtkActionGroup* belvuWindowGetActionGroup(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  return (properties ? properties->actionGroup : NULL);
}


/* Properties for generic windows */
static GenericWindowProperties* windowGetProperties(GtkWidget *widget)
{
  return widget ? (GenericWindowProperties*)(g_object_get_data(G_OBJECT(widget), "GenericWindowProperties")) : NULL;
}

static void onDestroyGenericWindow(GtkWidget *window)
{
  GenericWindowProperties *properties = windowGetProperties(window);

  /* We must remove the window from the list of spawned windows */
  properties->bc->spawnedWindows = g_slist_remove(properties->bc->spawnedWindows, window);
  
  if (properties)
    {
      /* Free the properties struct */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(window), "GenericWindowProperties", NULL);
    }
}


/* Create the properties struct and initialise all values for a generic window. */
static void genericWindowCreateProperties(GtkWidget *wrapWindow, 
                                          BelvuContext *bc,
                                          GtkActionGroup *actionGroup)
{
  if (wrapWindow)
    {
      GenericWindowProperties *properties = (GenericWindowProperties*)g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->actionGroup = actionGroup;

      g_object_set_data(G_OBJECT(wrapWindow), "GenericWindowProperties", properties);
      g_signal_connect(G_OBJECT(wrapWindow), "destroy", G_CALLBACK (onDestroyGenericWindow), NULL);
    }
}


/***********************************************************
 *                        Font size                        *
 ***********************************************************/

static void showFontDialog(BelvuContext *bc, GtkWidget *window)
{
  GtkWidget *entry = NULL;

  const int oldSize = pango_font_description_get_size(window->style->font_desc) / PANGO_SCALE;
  char *defaultText = g_strdup_printf("%d", oldSize);

  char *title = g_strdup_printf("%sFont", belvuGetTitlePrefix(bc));
  GtkWidget *dialog = createPromptDialog(window, defaultText, title, "Select font size:", "", &entry);
  g_free(title);

  gtk_widget_show_all(dialog);

  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const int newSize = convertStringToInt(gtk_entry_get_text(GTK_ENTRY(entry)));

      if (newSize >= MIN_FONT_SIZE && newSize <= MAX_FONT_SIZE)
        {
          widgetSetFontSize(bc->belvuWindow, GINT_TO_POINTER(newSize));
          widgetSetFontSize(bc->belvuTree, GINT_TO_POINTER(newSize));
          
          onBelvuAlignmentFontSizeChanged(bc->belvuAlignment);
          onBelvuTreeFontSizeChanged(bc->belvuTree);
        }
      else
        {
          g_critical("Invalid font size %d; font size must be between %d and %d.\n", newSize, MIN_FONT_SIZE, MAX_FONT_SIZE);
        }
    }

  gtk_widget_destroy(dialog);
}


void incrementFontSize(BelvuContext *bc)
{
  if (bc->belvuAlignment)
    {
      int size = pango_font_description_get_size(bc->belvuAlignment->style->font_desc) / PANGO_SCALE;
      widgetSetFontSizeAndCheck(bc->belvuAlignment, size + 1);
      onBelvuAlignmentFontSizeChanged(bc->belvuAlignment);
    }
}

void decrementFontSize(BelvuContext *bc)
{
  if (bc->belvuAlignment)
    {
      int size = pango_font_description_get_size(bc->belvuAlignment->style->font_desc) / PANGO_SCALE;
      widgetSetFontSizeAndCheck(bc->belvuAlignment, size - 1);
      onBelvuAlignmentFontSizeChanged(bc->belvuAlignment);
    }
}

/***********************************************************
 *		      Remove sequences			   *
 ***********************************************************/

/* This should be called when the 'removing sequences' option has changed */
static void updateSequenceRemovalMode(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  if (properties->bc->removingSeqs)
    {
      gdk_window_set_cursor(belvuWindow->window, properties->bc->removeSeqsCursor);
      g_message_info("Double-click on sequences to remove.  Esc or right-click to cancel.\n");
    }
  else
    {
      gdk_window_set_cursor(belvuWindow->window, properties->bc->defaultCursor);
      g_message_info("Finished removing sequences.\n");
    }

  /* Force cursor to change immediately */
    while (gtk_events_pending())
    gtk_main_iteration();
}


static void startRemovingSequences(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  if (!properties->bc->removingSeqs)
    {
      properties->bc->removingSeqs = TRUE;
      updateSequenceRemovalMode(belvuWindow);
    }
}

static void endRemovingSequences(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  properties->bc->removingSeqs = FALSE;
  updateSequenceRemovalMode(belvuWindow);
}


/***********************************************************
 *			About dialog			   *
 ***********************************************************/

/* A GtkAboutDialogActivateLinkFunc() called when user clicks on website link in "About" window. */
static void aboutDialogOpenLinkCB(GtkAboutDialog *about, const gchar *link, gpointer data)
{
  GError *error = NULL ;
    
  if (!seqtoolsLaunchWebBrowser(link, &error))
    g_critical("Cannot show link in web browser: \"%s\"", link) ;    
}


/* Shows the 'About' dialog */
static void showAboutDialog(GtkWidget *parent)
{
#if GTK_MAJOR_VERSION >= (2) && GTK_MINOR_VERSION >= (6)
  const gchar *authors[] = {AUTHOR_LIST, NULL} ;

  gtk_about_dialog_set_url_hook(aboutDialogOpenLinkCB, NULL, NULL) ;
  
  gtk_show_about_dialog(GTK_WINDOW(parent),
			"authors", authors,
			"comments", belvuGetCommentsString(), 
			"copyright", belvuGetCopyrightString(),
			"license", belvuGetLicenseString(),
			"name", belvuGetAppName(),
			"version", belvuGetVersionString(),
			"website", belvuGetWebSiteString(),
			NULL) ;
#endif
  
  return ;
}

/***********************************************************
 *                      Help dialog                        *
 ***********************************************************/

static void showHelpDialog()
{
  GError *error = NULL;

  /* The docs should live in /share/doc/seqtools/, in the same parent
   * directory that our executable's 'bin' directory is in. Open the 'quick
   * start' page. */
  char rel_path[100] = "../share/doc/seqtools/belvu_quick_start.html";

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
              g_message_info("Opening help page '%s'\n", path);
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
 *                   Save As dialog                        *
 ***********************************************************/

/* Utility to call the correct save function for the current save format.
 * Returns true if successful (false if not saved, e.g. if the user cancelled
 * the save dialog). */
static gboolean saveAlignment(BelvuContext *bc, GtkWidget *window)
{
  const char *filename = NULL;
  
  if (bc->saveFormat == BELVU_FILE_MSF)
    filename = saveMsf(bc, window);
  else if (bc->saveFormat == BELVU_FILE_ALIGNED_FASTA)
    filename = saveFasta(bc, window);
  else if (bc->saveFormat == BELVU_FILE_UNALIGNED_FASTA)
    filename = saveFasta(bc, window);
  else
    filename = saveMul(bc, window);
  
  /* Remember the last filename */
  if (filename)	
    {
      if (bc->dirName) g_free(bc->dirName);
      if (bc->fileName) g_free(bc->fileName);
      bc->dirName = g_path_get_dirname(filename);
      bc->fileName = g_path_get_basename(filename);
    }
  
  /* If the filename is null, the user must have cancelled. */
  return (filename != NULL);
}


/* This creates a drop-down box for selecting a file format */
static GtkComboBox* createFileFormatCombo(const int initFormatId)
{
  GtkComboBox *combo = createComboBox();
  GtkTreeIter *iter = NULL;
  
  int i = 0;
  for ( ; i < BELVU_NUM_FILE_FORMATS; ++i)
    addComboItem(combo, iter, i, getFileFormatString(i), initFormatId);
  
  return combo;
}


/* Callback called when the 'save coords' toggle button is toggled. It updates
 * the given widget (passed as data)  to be enabled if the toggle button is active
 * or disabled otherwise. */
static void onSaveCoordsToggled(GtkWidget *button, gpointer data)
{
  GtkWidget *otherWidget = GTK_WIDGET(data);
  
  if (otherWidget)
    {
      const gboolean isActive = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
      gtk_widget_set_sensitive(otherWidget, isActive); 
    }
}


static void showSaveAsDialog(BelvuContext *bc, GtkWidget *window)
{
  char *title = g_strdup_printf("%sSave As", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(window), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  GtkBox *contentArea = GTK_BOX(GTK_DIALOG(dialog)->vbox);
  
  /* Create a drop-down for selectin the file format */
  GtkWidget *hbox = gtk_hbox_new(FALSE, DIALOG_XPAD);
  gtk_box_pack_start(contentArea, hbox, FALSE, FALSE, DIALOG_XPAD);
  
  GtkWidget *label = gtk_label_new("Format: ");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, DIALOG_XPAD);
  
  GtkComboBox *combo = createFileFormatCombo(bc->saveFormat);
  gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo), FALSE, FALSE, DIALOG_YPAD);
  
  /* Create a tick box for enabling the 'save coords' option */
  GtkToggleButton *checkButton = GTK_TOGGLE_BUTTON(gtk_check_button_new_with_mnemonic("Save _coordinates"));
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkButton), bc->saveCoordsOn);
  gtk_box_pack_start(contentArea, GTK_WIDGET(checkButton), FALSE, FALSE, DIALOG_YPAD);
  
  /* Create selection buttons to allow user to specify separator char */
  GtkBox *hbox2 = GTK_BOX(gtk_hbox_new(FALSE, DIALOG_XPAD));
  gtk_box_pack_start(contentArea, GTK_WIDGET(hbox2), FALSE, FALSE, DIALOG_YPAD);
  
  GtkWidget *separatorLabel = gtk_label_new("Separator character between name and coords:\n(Use = for GCG)");
  gtk_misc_set_alignment(GTK_MISC(separatorLabel), 0.0, 0.0);
  gtk_box_pack_start(hbox2, separatorLabel, FALSE, FALSE, DIALOG_XPAD);
  
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, DIALOG_YPAD));
  gtk_box_pack_start(hbox2, GTK_WIDGET(vbox), FALSE, FALSE, DIALOG_XPAD);
  
  const gboolean button1Active = (bc->saveSeparator == '/');
  GtkWidget *button1 = gtk_radio_button_new_with_label_from_widget(NULL, "slash  (/) ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button1), button1Active);
  gtk_box_pack_start(vbox, button1, FALSE, FALSE, 0);
  
  GtkWidget *button2 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(button1), "equals  (=) ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button2), !button1Active);
  gtk_box_pack_start(vbox, button2, FALSE, FALSE, 0);
  
  /* The save separator is only applicable if saving coords */
  gtk_widget_set_sensitive(GTK_WIDGET(hbox2), bc->saveCoordsOn); 
  g_signal_connect(G_OBJECT(checkButton), "toggled", G_CALLBACK(onSaveCoordsToggled), hbox2);
  
  
  gtk_widget_show_all(dialog);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      bc->saveFormat = (BelvuFileFormat)gtk_combo_box_get_active(combo);
      bc->saveCoordsOn = gtk_toggle_button_get_active(checkButton);
    
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button1)))
	bc->saveSeparator = '/';
      else 
	bc->saveSeparator = '=';
    
    
      saveAlignment(bc, window);
    }
  
  gtk_widget_destroy(dialog);
}


/***********************************************************
 *                     Find dialog                         *
 ***********************************************************/

static void findSeqs(BelvuContext *bc, const char *searchStr, const gboolean findAgain, const gboolean searchBackwards)
{
  /* Loop through each alignment and see if the name matches. If so,
   * remember the index of the matching alignment. We start searching
   * from the beginning of the alignment array, unless this is a 'find again'
   * search, when we search from the next alignment after the previous result */
  static int startIdx = 0;

  if (findAgain && searchBackwards)
    --startIdx;
  else if (findAgain)
    ++startIdx;
  else
    startIdx = 0;

  const int increment = (searchBackwards ? -1 : 1);

  int i = startIdx;
  gboolean found = FALSE;

  for ( ; i >= 0 && i < (int)bc->alignArr->len; i += increment)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

      if (wildcardSearch(alnp->name, searchStr))
        {
          bc->selectedAln = alnp;
          onRowSelectionChanged(bc);

          found = TRUE;
          startIdx = i;
          break;
        }
    }

  /* If it's a find-again search and we failed to find a result, try
   * again starting from the beginning */
  if (findAgain && !found &&
      ((searchBackwards && startIdx != (int)bc->alignArr->len - 1) || (!searchBackwards && startIdx != 0)) )
    {
      startIdx = (searchBackwards ? bc->alignArr->len - 1 : 0);
      i = startIdx;

      for ( ; i >= 0 && i < (int)bc->alignArr->len; i += increment)
        {
          ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

          if (wildcardSearch(alnp->name, searchStr))
            {
              bc->selectedAln = alnp;
              onRowSelectionChanged(bc);

              found = TRUE;
              startIdx = i;
              break;
            }
        }
    }

  if (!found)
    g_critical("Alignment name '%s' not found.\n", searchStr);
}


static gboolean onFindSeqs(GtkWidget *button, const gint responseId, gpointer data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      const char *searchStr = getStringFromTextEntry(GTK_ENTRY(data));

      GtkWidget *dialog = gtk_widget_get_toplevel(button);
      GtkWidget *window = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
      BelvuContext *bc = windowGetContext(window);

      /* If the the user hit forward or back, do a find-again search in the appropriate
       * direction; otherwise do a normal search */
      if (responseId == BLX_RESPONSE_FORWARD)
        findSeqs(bc, searchStr, TRUE, FALSE);
      else if (responseId == BLX_RESPONSE_BACK)
        findSeqs(bc, searchStr, TRUE, TRUE);
      else
        findSeqs(bc, searchStr, FALSE, FALSE);
    }

  return TRUE;
}


/* To do: implement functionality to search for sequences by a pattern of residues */
/* 
static gboolean onFindResidues(GtkWidget *button, const gint responseId, gpointer data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    {
      const char *searchStr = getStringFromTextEntry(GTK_ENTRY(data));

      GtkWidget *dialog = gtk_widget_get_toplevel(button);
      GtkWidget *window = GTK_WIDGET(gtk_window_get_transient_for(GTK_WINDOW(dialog)));
      BelvuContext *bc = windowGetContext(window);

      if (responseId == BLX_RESPONSE_FORWARD)
        findResidues(bc, searchStr, TRUE, FALSE);
      else if (responseId == BLX_RESPONSE_BACK)
        findResidues(bc, searchStr, TRUE, TRUE);
      else
        findResidues(bc, searchStr, FALSE, FALSE);
    }

  return TRUE;
}
*/

static void showFindDialog(BelvuContext *bc, GtkWidget *window)
{
  const BelvuDialogId dialogId = BELDIALOG_FIND;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);

  if (!dialog)
    {
      char *title = g_strdup_printf("%sFind sequences", belvuGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(window),
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_GO_BACK,
                                           BLX_RESPONSE_BACK,
                                           GTK_STOCK_GO_FORWARD,
                                           BLX_RESPONSE_FORWARD,
                                           GTK_STOCK_CLOSE,
                                           GTK_RESPONSE_REJECT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);

      GtkBox *contentArea = GTK_BOX(GTK_DIALOG(dialog)->vbox);
      const int numRows = 2;
      const int numCols = 2;
      GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
      gtk_box_pack_start(contentArea, GTK_WIDGET(table), TRUE, TRUE, 0);


      /*GtkRadioButton *button1 =*/ createRadioButton(table, 0, 0, NULL, "_Name search (wildcards * and ?)", TRUE, TRUE, FALSE, onFindSeqs, window, NULL);
      /*createRadioButton(table, 0, 1, button1, "_Residue sequence search", FALSE, TRUE, FALSE, onFindResidues, window, NULL);*/


      gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(window));
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
    }

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
}


/***********************************************************
 *                Remove sequences dialogs                 *
 ***********************************************************/

static GtkWidget* createPromptDialog(GtkWidget *window,
                                     const char *defaultResult,
                                     const char *title,
                                     const char *text1,
                                     const char *text2,
                                     GtkWidget **entry)
{ 
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(window),
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  if (window)
    pango_font_description_set_size(dialog->style->font_desc, pango_font_description_get_size(window->style->font_desc));
  
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, FALSE, FALSE, 12);
  
  GtkWidget *label1 = gtk_label_new(text1);
  gtk_misc_set_alignment(GTK_MISC(label1), 1, 0.5);
  gtk_box_pack_start(GTK_BOX(hbox), label1, FALSE, FALSE, 0);
  
  *entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), *entry, FALSE, FALSE, 0);

  gtk_entry_set_width_chars(GTK_ENTRY(*entry), 3);
  gtk_entry_set_activates_default(GTK_ENTRY(*entry), TRUE);

  gtk_entry_set_text(GTK_ENTRY(*entry), defaultResult);
  
  GtkWidget *label2 = gtk_label_new(text2);
  gtk_misc_set_alignment(GTK_MISC(label2), 0, 0.5);
  gtk_box_pack_start(GTK_BOX(hbox), label2, FALSE, FALSE, 0);

  gtk_widget_show_all(dialog);
  
  return dialog;
}


/* Show a dialog to ask the user what threshold to use for removing
 * "gappy" sequences (i.e. sequences that are more than the given 
 * percentage of gaps). */
static void showRemoveGappySeqsDialog(GtkWidget *belvuWindow)
{
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup("50");
  
  GtkWidget *entry = NULL;
  BelvuContext *bc = windowGetContext(belvuWindow);
  
  char *title = g_strdup_printf("%sremove sequences", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Remove sequences that are ", "% or more gaps.", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      if (inputText)
	g_free(inputText);
      
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
      const gdouble cutoff = g_strtod(inputText, NULL);
      
      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      removeGappySeqs(properties->bc, properties->bc->belvuAlignment, cutoff);
    }
  
  gtk_widget_destroy(dialog);
}


static void showMakeNonRedundantDialog(GtkWidget *belvuWindow)
{
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup("80.0");
  
  GtkWidget *entry = NULL;
  BelvuContext *bc = windowGetContext(belvuWindow);
  char *title = g_strdup_printf("%sremove sequences", belvuGetTitlePrefix(bc));
  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Remove sequences that are more than ", "% identical.", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      if (inputText)
	g_free(inputText);
      
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
      const gdouble cutoff = g_strtod(inputText, NULL);
      
      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      removeRedundantSeqs(properties->bc, properties->bc->belvuAlignment, cutoff);
    }
  
  gtk_widget_destroy(dialog);
}


static void showRemoveOutliersDialog(GtkWidget *belvuWindow)
{
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup("20.0");
  
  GtkWidget *entry = NULL;
  BelvuContext *bc = windowGetContext(belvuWindow);
  char *title = g_strdup_printf("%sremove sequences", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Remove sequences that are less than ", "% identical with any other.", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      if (inputText)
	g_free(inputText);
      
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
      const gdouble cutoff = g_strtod(inputText, NULL);
      
      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      removeOutliers(properties->bc, properties->bc->belvuAlignment, cutoff);
    }
  
  gtk_widget_destroy(dialog);
}


static void showRemoveByScoreDialog(GtkWidget *belvuWindow)
{
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup("20.0");
  
  GtkWidget *entry = NULL;
  BelvuContext *bc = windowGetContext(belvuWindow);
  char *title = g_strdup_printf("%sremove sequences", belvuGetTitlePrefix(bc));
  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Remove sequences that have a score less than ", "", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      if (inputText)
	g_free(inputText);
      
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
      const gdouble cutoff = g_strtod(inputText, NULL);
      
      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      removeByScore(properties->bc, properties->bc->belvuAlignment, cutoff);
    }
  
  gtk_widget_destroy(dialog);
}


/***********************************************************
 *                     Remove columns                      *
 ***********************************************************/

static void showRemoveColumnsDialog(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  char *title = g_strdup_printf("%sRemove Columns", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(belvuWindow), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, FALSE, FALSE, 12);
  
  GtkWidget *label1 = gtk_label_new("Remove columns from");
  gtk_misc_set_alignment(GTK_MISC(label1), 1, 0.5);
  gtk_box_pack_start(GTK_BOX(hbox), label1, FALSE, FALSE, 0);
  
  char *maxLenText = g_strdup_printf("%d", bc->maxLen);
  GtkWidget *entry1 = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), entry1, FALSE, FALSE, 0);
  gtk_entry_set_width_chars(GTK_ENTRY(entry1), strlen(maxLenText) + 1);
  gtk_entry_set_activates_default(GTK_ENTRY(entry1), TRUE);
  gtk_entry_set_text(GTK_ENTRY(entry1), "1");
  
  GtkWidget *label2 = gtk_label_new("to");
  gtk_misc_set_alignment(GTK_MISC(label2), 0, 0.5);
  gtk_box_pack_start(GTK_BOX(hbox), label2, FALSE, FALSE, 0);

  GtkWidget *entry2 = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), entry2, FALSE, FALSE, 0);
  gtk_entry_set_width_chars(GTK_ENTRY(entry2), strlen(maxLenText) + 3);
  gtk_entry_set_activates_default(GTK_ENTRY(entry2), TRUE);
  gtk_entry_set_text(GTK_ENTRY(entry2), maxLenText);
  g_free(maxLenText);
  maxLenText = NULL;
  
  gtk_widget_show_all(dialog);

  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const char *inputText1 = gtk_entry_get_text(GTK_ENTRY(entry1));
      const int fromVal = convertStringToInt(inputText1);
      
      const char *inputText2 = gtk_entry_get_text(GTK_ENTRY(entry2));
      const int toVal = convertStringToInt(inputText2);

      rmColumn(bc, fromVal, toVal);
      
      rmFinaliseColumnRemoval(bc);
      updateOnAlignmentLenChanged(bc->belvuAlignment);
    }
  
  gtk_widget_destroy(dialog);
}


/* Remove columns with conservation below a given cutoff */
static void showRemoveColumnsCutoffDialog(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  static char *fromText = NULL;
  static char *toText = NULL;

  if (!fromText)
    fromText = g_strdup_printf("%.2f", -1.0);
  
  if (!toText)
    toText = g_strdup_printf("%.2f", 0.9);
  
  char *title = g_strdup_printf("%sRemove Columns", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(belvuWindow), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, DIALOG_YPAD);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), vbox, FALSE, FALSE, DIALOG_YPAD);
  
  /* Place a label at the top */
  GtkWidget *label1 = gtk_label_new("Remove columns with a (maximum) conservation between: ");
  gtk_misc_set_alignment(GTK_MISC(label1), 0.0, 0.5);
  gtk_box_pack_start(GTK_BOX(vbox), label1, FALSE, FALSE, 0);
  
  /* Place the text entry boxes in an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, DIALOG_YPAD);

  gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(""), FALSE, FALSE, DIALOG_XPAD);
  
  GtkWidget *entry1 = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), entry1, FALSE, FALSE, DIALOG_XPAD);
  gtk_entry_set_width_chars(GTK_ENTRY(entry1), strlen(fromText) + 1);
  gtk_entry_set_activates_default(GTK_ENTRY(entry1), TRUE);
  gtk_entry_set_text(GTK_ENTRY(entry1), fromText);
  
  GtkWidget *label2 = gtk_label_new("<   conservation   <=");
  gtk_misc_set_alignment(GTK_MISC(label2), 0, 0.5);
  gtk_box_pack_start(GTK_BOX(hbox), label2, FALSE, FALSE, 0);
  
  GtkWidget *entry2 = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), entry2, FALSE, FALSE, DIALOG_XPAD);
  gtk_entry_set_width_chars(GTK_ENTRY(entry2), strlen(toText) + 3);
  gtk_entry_set_activates_default(GTK_ENTRY(entry2), TRUE);
  gtk_entry_set_text(GTK_ENTRY(entry2), toText);
  
  gtk_widget_show_all(dialog);
  
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      g_free(fromText);
      fromText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry1)));
      const double fromVal = g_strtod(fromText, NULL);

      g_free(toText);
      toText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry2)));
      const double toVal = g_strtod(toText, NULL);
    
      rmColumnCutoff(bc, fromVal, toVal);
      updateOnAlignmentLenChanged(bc->belvuAlignment);
    }
  
  gtk_widget_destroy(dialog);
}


/* Remove columns with a higher fraction of gaps than a specified cutoff */
static void showRemoveGappyColumnsDialog(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup_printf("%.0f", 50.0);
  
  GtkWidget *entry = NULL;
  char *title = g_strdup_printf("%sRemove Columns", belvuGetTitlePrefix(bc));
  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Remove columns with more than ", " % gaps", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      g_free(inputText);
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));
      const double cutoff = g_strtod(inputText, NULL);
      
      rmEmptyColumns(bc, cutoff/100.0);
      
      rmFinaliseColumnRemoval(bc);
      updateOnAlignmentLenChanged(bc->belvuAlignment);
    }
  
  gtk_widget_destroy(dialog);
}

/***********************************************************
 *               Select gap character dialog                *
 ***********************************************************/

static void showSelectGapCharDialog(GtkWidget *belvuWindow)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;

  char *title = g_strdup_printf("%sGap Character", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title,
                                                  GTK_WINDOW(belvuWindow), 
                                                  (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_MODAL),
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  NULL);

  g_free(title);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), belvuWindow);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  GtkBox *hbox = GTK_BOX(gtk_hbox_new(FALSE, 12));
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), GTK_WIDGET(hbox), FALSE, FALSE, 12);

  GtkWidget *label = gtk_label_new("Select gap character:");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_box_pack_start(hbox, label, FALSE, FALSE, 12);
  
  /* Create radio buttons for each gap character (only dot or dash at the moment) */
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 12));
  gtk_box_pack_start(hbox, GTK_WIDGET(vbox), TRUE, TRUE, 0);
  
  const gboolean button1Active = (bc->gapChar == '.');
  
  GtkWidget *button1 = gtk_radio_button_new_with_label_from_widget(NULL, "dot  (.) ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button1), button1Active);
  gtk_box_pack_start(vbox, button1, FALSE, FALSE, 0);

  GtkWidget *button2 = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(button1), "dash  (-) ");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button2), !button1Active);
  gtk_box_pack_start(vbox, button2, FALSE, FALSE, 0);
  
  gtk_window_set_default_size(GTK_WINDOW(dialog), 300, -1);
  gtk_widget_show_all(dialog);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      char newChar;
      
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button1)))
        newChar = '.';
      else 
        newChar = '-';
      
      if (newChar != bc->gapChar)
        {
          bc->gapChar = newChar;
          
          /* Change all the gaps in the sequences to use the new gap char */
          int i,j;
          for (i = 0; i < (int)bc->alignArr->len; ++i) 
            {
              ALN *alnp = g_array_index(bc->alignArr, ALN*, i);
              char *alnpSeq = alnGetSeq(alnp);
              
              for (j = 0; j < bc->maxLen; ++j) 
                {
                  if (isGap(alnpSeq[j])) 
                    alnpSeq[j] = bc->gapChar;
                }
            }
          
          belvuAlignmentRedrawAll(bc->belvuAlignment);
        }
    }
  
  gtk_widget_destroy(dialog);
}


/***********************************************************
 *                Color by residue ID dialog               *
 ***********************************************************/

/* Dialog to prompt the user to enter a threshold for coloring residues by ID */
static void showColorByResIdDialog(GtkWidget *belvuWindow)
{
  static char *inputText = NULL;
  
  if (!inputText)
    inputText = g_strdup("20.0");
  
  GtkWidget *entry = NULL;
  BelvuContext *bc = windowGetContext(belvuWindow);
  char *title = g_strdup_printf("%sColor by Residue ID", belvuGetTitlePrefix(bc));
  GtkWidget *dialog = createPromptDialog(belvuWindow, inputText, title, "Only colour residues above ", "% identity", &entry);
  g_free(title);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      if (inputText)
	g_free(inputText);
      
      inputText = g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));

      BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
      properties->bc->colorByResIdCutoff = g_strtod(inputText, NULL);
    
      /* This sets the flag and also updates the associated 'toggle' menu item */
      setToggleMenuStatus(properties->actionGroup, "toggleColorByResId", TRUE);
    
      /* Update the color scheme and redraw */
      updateSchemeColors(properties->bc);
      belvuAlignmentRedrawAll(properties->bc->belvuAlignment);
    }
  
  gtk_widget_destroy(dialog);  
}


/***********************************************************
 *                    Edit Colors dialog                   *
 ***********************************************************/

/* This creates a drop-down box for selecting a color */
static GtkComboBox* createColorCombo(const int colorNum)
{
  GtkComboBox *combo = createComboBox();
  GtkTreeIter *iter = NULL;
  
  int i = 0;
  for (i = 0; i < NUM_TRUECOLORS; ++i)
    addComboItem(combo, iter, i, getColorNumName(i), colorNum);
  
  return combo;
}


/* This is called when a color combo box has been changed and it updates
 * the given color button (passed as the user data) to be the same color */
static void updateColorButton(GtkWidget *combo, gpointer data)
{
  GtkWidget *colorButton = GTK_WIDGET(data);
  
  const int colorNum = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
  
  GdkColor color;
  convertColorNumToGdkColor(colorNum, FALSE, &color);
  
  gtk_widget_modify_bg(colorButton, GTK_STATE_NORMAL, &color);
}


/* This is called when a color combo box has been changed and it updates the
 * given residue (passed as the user data) */
static void updateColorResidue(GtkWidget *combo, gpointer data)
{
  const int colorNum = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
  const char *residue = (const char*)data;
  
  setColor(*residue, colorNum);
  
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(combo));
  GtkWidget *belvuWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  onColorSchemeChanged(properties);
}


/* Create a color chooser button. This is a bit of a hack because it has to 
 * deal with the old-style acedb colors, whereas ideally we would get rid of
 * those and just use a GTK color-button widget. */
static void createColorButton(const int colorNum, 
                              GtkTable *table, 
                              const int row, 
                              const int col,
                              const int xpad,
                              const int ypad,
                              GCallback callbackFunc,
                              gpointer data)
{
  /* Create the color chooser. First get the current color for this residue */
  GdkColor color;
  convertColorNumToGdkColor(colorNum, FALSE, &color);
  
  /* Create a color box to display the currently-selected color. NB ideally we would get
   * rid of the drop-down box and have a color-picker button instead, but this requires
   * changes to the way colors are stored and saved in order to support more colors than
   * just the old-style acedb list. */
  GtkWidget *colorButton = gtk_event_box_new();
  gtk_table_attach(table, colorButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  gtk_widget_set_size_request(colorButton, 40, 20);
  gtk_widget_modify_bg(colorButton, GTK_STATE_NORMAL, &color);
  
  /* Create a drop-down box to allow the user to select one of the old-style acedb colors. */
  GtkWidget *combo = GTK_WIDGET(createColorCombo(colorNum));
  gtk_table_attach(table, combo, col + 1, col + 2, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);

  /* This just updates the gtk color button with the new color when the combo box value has changed. */
  g_signal_connect(G_OBJECT(combo), "changed", G_CALLBACK(updateColorButton), colorButton);

  /* This attaches the 'proper' user callback function */
  g_signal_connect(G_OBJECT(combo), "changed", callbackFunc, data);
}


/* This creates a single color item on the edit-residue-colors dialog */
static void createResidueColorBox(GtkTable *table, 
                                  char *residue,
                                  GString *groups[],
                                  int *row, 
                                  int *col,
                                  const int xpad,
                                  const int ypad)
{
  /* Create the label */
  GtkWidget *label = gtk_label_new(residue);
  gtk_table_attach(table, label, *col, *col + 1, *row, *row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create the color chooser. First get the current color for this residue */
  const int colorNum = getColor(*residue);
  createColorButton(colorNum, table, *row, *col + 1, xpad, ypad, G_CALLBACK(updateColorResidue), residue);
  
  /* Append this residue to the group for this color */
//to do:  g_string_append(groups[colorNum], residue);
  
  *row += 1;
}


/* Utility function to create the list of residues that belvu knows about */
static GSList* createResidueList()
{
  GSList *list = NULL;
  
  list = g_slist_prepend(list, g_strdup("V"));
  list = g_slist_prepend(list, g_strdup("Y"));
  list = g_slist_prepend(list, g_strdup("W"));
  list = g_slist_prepend(list, g_strdup("T"));
  list = g_slist_prepend(list, g_strdup("S"));
  list = g_slist_prepend(list, g_strdup("P"));
  list = g_slist_prepend(list, g_strdup("F"));
  list = g_slist_prepend(list, g_strdup("M"));
  list = g_slist_prepend(list, g_strdup("K"));
  list = g_slist_prepend(list, g_strdup("L"));
  list = g_slist_prepend(list, g_strdup("I"));
  list = g_slist_prepend(list, g_strdup("H"));
  list = g_slist_prepend(list, g_strdup("G"));
  list = g_slist_prepend(list, g_strdup("E"));
  list = g_slist_prepend(list, g_strdup("Q"));
  list = g_slist_prepend(list, g_strdup("C"));
  list = g_slist_prepend(list, g_strdup("D"));
  list = g_slist_prepend(list, g_strdup("N"));
  list = g_slist_prepend(list, g_strdup("R"));
  list = g_slist_prepend(list, g_strdup("A"));
  
  return list;
}


/* Create the colour boxes for the edit-residue-colors dialog */
static void createResidueColorBoxes(GtkBox *box, GSList *residueList, GString *groups[])
{
  const int numRows = 10;
  const int numCols = 7;
  int xpad = 4;
  int ypad = TABLE_YPAD;
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_box_pack_start(box, GTK_WIDGET(table), FALSE, FALSE, 0);
  
  int row = 0;
  int col = 0;
  
  int i = 0;
  GSList *residueItem = residueList;
  const int numResidues = g_slist_length(residueList);
  
  for ( ; residueItem && i < numResidues / 2; ++i, residueItem = residueItem->next)
    {
      createResidueColorBox(table, (char*)(residueItem->data), groups, &row, &col, xpad, ypad);
    }

  GtkWidget *spacer = gtk_label_new(" ");
  gtk_table_attach(table, spacer, col + 3, col + 4, row, row + 1, GTK_SHRINK, GTK_SHRINK, 30, ypad);

  col += 4;
  row = 0;

  for ( ; residueItem; residueItem = residueItem->next)
    {
      createResidueColorBox(table, (char*)(residueItem->data), groups, &row, &col, xpad, ypad);
    }
}


/* Create the content for the edit-residue-colors dialog */
static void createEditResidueContent(GtkBox *box)
{
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(box, GTK_WIDGET(vbox), TRUE, TRUE, 0);
  
  static GSList *residueList = NULL;
  
  if (!residueList)
    residueList = createResidueList();
  
//  GString* groups[NUM_TRUECOLORS];
//  int i = 0;
//  for (i = 0; i < NUM_TRUECOLORS; ++i)
//    groups[i] = g_string_new("");
  
  createResidueColorBoxes(vbox, residueList, NULL);
  
  
//  GtkTable *table = GTK_TABLE(gtk_table_new(g_slist_length(residueList) + 1, 3, FALSE));
//  gtk_box_pack_start(vbox, GTK_WIDGET(table), TRUE, TRUE, 0);
//  int row = 0;
//  const int xpad = TABLE_XPAD;
//  const int ypad = TABLE_YPAD;
//  
//  GtkWidget *header = gtk_label_new("Groups:");
//  gtk_table_attach(table, header, 0, 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
//  ++row;
//  
//  for (i = 0; i < NUM_TRUECOLORS; ++i)
//    {
//      if (groups[i]->len > 0 && groups[i]->str)
//        {
//          GdkColor color;
//          convertColorNumToGdkColor(i, FALSE, &color);
//          
//          GtkWidget *eventBox = gtk_event_box_new();
//          gtk_widget_modify_bg(eventBox, GTK_STATE_NORMAL, &color);
//          gtk_table_attach(table, eventBox, 0, 1, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
//          
//          GtkWidget *label = gtk_label_new(getColorNumName(i));
//          gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
//          gtk_container_add(GTK_CONTAINER(eventBox), label);
//          
//          gtk_table_attach(table, gtk_label_new(":"), 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
//
//          label = gtk_label_new(groups[i]->str);
//          gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
//          gtk_table_attach(table, label, 2, 3, row, row + 1, GTK_FILL, GTK_SHRINK, xpad, ypad);
//          
//          ++row;
//        }
//    }
// 
}


/* Called when the user responds to the edit-residue-colors dialog */
void onResponseEditResidueColorsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      if (widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId)))
        {
          /* Save the current color set and set the color scheme to 'custom' */
          destroy = TRUE;
          saveCustomColors(bc);
          setToggleMenuStatus(properties->actionGroup, "colorSchemeCustom", TRUE);
        }
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      /* Reset the color scheme, refresh, and close the dialog. */
      destroy = TRUE;
      setResidueSchemeColors(bc);
      updateSchemeColors(bc);
      belvuAlignmentRedrawAll(bc->belvuAlignment);
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
}


/* Show a dialog to allow the user to edit the residue colors */
static void showEditResidueColorsDialog(GtkWidget *belvuWindow, const gboolean bringToFront)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  const BelvuDialogId dialogId = BELDIALOG_EDIT_RESIDUE_COLORS;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sEdit Residue Colors", belvuGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(belvuWindow), 
                                           (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_MODAL),
                                           GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                           GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseEditResidueColorsDialog), belvuWindow);
    }
  else
    {
      /* Need to refresh the dialog contents, so clear and re-create content area */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkBox *vbox = GTK_BOX(GTK_DIALOG(dialog)->vbox);
  createEditResidueContent(vbox);

  /* Show / bring to front */
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
  
}



/* Called when the user responds to the edit-conservation-colors dialog */
void onResponseConsColorsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  GtkWidget *belvuWindow = GTK_WIDGET(data);
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Close dialog if successful */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      setToggleMenuStatus(properties->actionGroup, "ColorByCons", TRUE);
      break;
      
    case GTK_RESPONSE_APPLY:
      /* Never close */
      destroy = FALSE;
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      setToggleMenuStatus(properties->actionGroup, "ColorByCons", TRUE);
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      /* Reset the color scheme. */
      destroy = TRUE;
      saveOrResetConsColors(bc, FALSE); /* restores old values */
      break;
      
    default:
      break;
  };
  
  onColorSchemeChanged(properties);

  if (destroy)
    {
      gtk_widget_hide_all(GTK_WIDGET(dialog));
    }
}


/* Callback called when a foreground conservation color has been changed */
static void updateConsFgColor(GtkWidget *combo, gpointer data)
{
  int *colorNum = (int*)data;
  *colorNum = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(combo));
  GtkWidget *belvuWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  onColorSchemeChanged(properties);
}


/* Callback called when a background conservation color has been changed */
static void updateConsBgColor(GtkWidget *combo, gpointer data)
{
  int *colorNum = (int*)data;
  *colorNum = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
  
  GtkWindow *dialogWindow = GTK_WINDOW(gtk_widget_get_toplevel(combo));
  GtkWidget *belvuWindow = GTK_WIDGET(gtk_window_get_transient_for(dialogWindow));
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  
  onColorSchemeChanged(properties);
}

/* Callback called when the user 'ok's the edit-conservation-colors dialog
 * to set the new threshold values. The given widget must be a GtkEntry and the
 * user data a pointer to a double. */
static gboolean onConsThresholdChanged(GtkWidget *widget, gint responseId, gpointer data)
{
  GtkEntry *entry = GTK_ENTRY(widget);
  double *val = (double*)data;
  
  const gchar *inputText = gtk_entry_get_text(entry);
  *val = g_strtod(inputText, NULL);
  
  return TRUE;
}


/* Add a single line in the edit-cons-colors dialog */
static void addConsColorLine(BelvuContext *bc, const char *labelText, const BelvuConsLevel consLevel, double *cutoff, GtkTable *table, int *row)
{
  /* Label */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.0);
  gtk_table_attach(table, label, 0, 1, *row, *row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  /* Threshold entry box */
  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
  
  char *defaultInput = g_strdup_printf("%.1f", *cutoff);
  gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);
  gtk_entry_set_width_chars(GTK_ENTRY(entry), strlen(defaultInput) + 2);
  g_free(defaultInput);
  defaultInput = NULL;
  
  gtk_table_attach(table, entry, 1, 2, *row, *row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  widgetSetCallbackData(entry, onConsThresholdChanged, cutoff);
  
  /* Text color chooser */
  int *fgColorNum = getConsColor(bc, consLevel, TRUE);
  createColorButton(*fgColorNum, table, *row, 2, 2, TABLE_YPAD, G_CALLBACK(updateConsFgColor), fgColorNum);

  /* Background color chooser */
  int *bgColorNum = getConsColor(bc, consLevel, FALSE);
  createColorButton(*bgColorNum, table, *row, 4, 2, TABLE_YPAD, G_CALLBACK(updateConsBgColor), bgColorNum);

  *row += 1;
}


/* Hacky function to save the current conservation colors or reset them to
 * previously saved values (if 'save' is false). Ideally, the drawing 
 * functions would not use the colors in the BelvuContext directly and we'd
 * therefore be able to show the user different colors without having to change
 * the BelvuContext, and hence without the need for this hacky undo function. */
static void saveOrResetConsColors(BelvuContext *bc, const gboolean save)
{
  /* These are all the values we need to be able to restore if the dialog
   * is cancelled... */
  static double lowIdCutoff = 0;
  static double midIdCutoff = 0;
  static double maxIdCutoff = 0;
  static double lowSimCutoff = 0;
  static double midSimCutoff = 0;
  static double maxSimCutoff = 0;
  
  static int maxfgColor = 0;
  static int midfgColor = 0;
  static int lowfgColor = 0;
  static int maxbgColor = 0;
  static int midbgColor = 0;
  static int lowbgColor = 0;
  static int maxfgPrintColor = 0;
  static int midfgPrintColor = 0;
  static int lowfgPrintColor = 0;
  static int maxbgPrintColor = 0;
  static int midbgPrintColor = 0;
  static int lowbgPrintColor = 0;  
  if (save)
    {
      /* Remember the current conservation colors */
      lowIdCutoff = bc->lowIdCutoff;
      midIdCutoff = bc->midIdCutoff;
      maxIdCutoff = bc->maxIdCutoff;
      lowSimCutoff = bc->lowSimCutoff;
      midSimCutoff = bc->midSimCutoff;
      maxSimCutoff = bc->maxSimCutoff;

      maxfgColor = bc->maxfgColor;
      midfgColor = bc->midfgColor;
      lowfgColor = bc->lowfgColor;
      maxbgColor = bc->maxbgColor;
      midbgColor = bc->midbgColor;
      lowbgColor = bc->lowbgColor;
      maxfgPrintColor = bc->maxfgPrintColor;
      midfgPrintColor = bc->midfgPrintColor;
      lowfgPrintColor = bc->lowfgPrintColor;
      maxbgPrintColor = bc->maxbgPrintColor;
      midbgPrintColor = bc->midbgPrintColor;
      lowbgPrintColor = bc->lowbgPrintColor;
    }
  else
    {
      /* Reset to the previously-saved values */
      bc->lowIdCutoff = lowIdCutoff;
      bc->midIdCutoff = midIdCutoff;
      bc->maxIdCutoff = maxIdCutoff;
      bc->lowSimCutoff = lowSimCutoff;
      bc->midSimCutoff = midSimCutoff;
      bc->maxSimCutoff = maxSimCutoff;

      bc->maxfgColor = maxfgColor;
      bc->midfgColor = midfgColor;
      bc->lowfgColor = lowfgColor;
      bc->maxbgColor = maxbgColor;
      bc->midbgColor = midbgColor;
      bc->lowbgColor = lowbgColor;
      bc->maxfgPrintColor = maxfgPrintColor;
      bc->midfgPrintColor = midfgPrintColor;
      bc->lowfgPrintColor = lowfgPrintColor;
      bc->maxbgPrintColor = maxbgPrintColor;
      bc->midbgPrintColor = midbgPrintColor;
      bc->lowbgPrintColor = lowbgPrintColor;
    }
}


/* Create the content for the edit-conservation-color-scheme dialog */
static void createEditConsColorsContent(GtkBox *box, BelvuContext *bc)
{
  /* Save current values so that we can restore them later if we cancel */
  saveOrResetConsColors(bc, TRUE);
  
  /* We'll put everything in a vbox */
  GtkBox *vbox = GTK_BOX(gtk_vbox_new(FALSE, 0));
  gtk_box_pack_start(box, GTK_WIDGET(vbox), FALSE, FALSE, 0);
  
  /* Create the labels */
  char *tmpStr = g_strdup_printf("Coloring by %s", bc->consScheme == BELVU_SCHEME_BLOSUM ? "average BLOSUM62 score" : "% identity");
  
  GtkWidget *label = gtk_label_new(tmpStr);
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_box_pack_start(vbox, label, FALSE, FALSE, DIALOG_YPAD);
  
  g_free(tmpStr);
  tmpStr = NULL;

  const gboolean colorById = (bc->consScheme == BELVU_SCHEME_ID || bc->consScheme == BELVU_SCHEME_ID_BLOSUM);
  const gboolean colorBySim = (bc->consScheme == BELVU_SCHEME_BLOSUM || bc->consScheme == BELVU_SCHEME_ID_BLOSUM);
  
  if (colorBySim)
    {
      label = gtk_label_new("Similar residues according to BLOSUM62 are coloured as the most conserved one.");
      gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
      gtk_box_pack_start(vbox, label, FALSE, FALSE, DIALOG_YPAD);
    }
  
  /* Place the main widgets in a table */
  GtkTable *table = GTK_TABLE(gtk_table_new(4, 6, FALSE));
  gtk_box_pack_start(vbox, GTK_WIDGET(table), TRUE, TRUE, DIALOG_YPAD);
  
  /* 1st row contains labels */
  int row = 0;
  label = gtk_label_new("Threshold");\
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_table_attach(table, label, 1, 2, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);

  label = gtk_label_new("Text colour");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_table_attach(table, label, 2, 4, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);

  label = gtk_label_new("Background colour");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_table_attach(table, label, 4, 6, row, row + 1, GTK_FILL, GTK_SHRINK, TABLE_XPAD, TABLE_YPAD);
  
  ++row;
  addConsColorLine(bc, "Max:", CONS_LEVEL_MAX, colorById ? &bc->maxIdCutoff : &bc->maxSimCutoff, table, &row);
  addConsColorLine(bc, "Mid:", CONS_LEVEL_MID, colorById ? &bc->midIdCutoff : &bc->midSimCutoff, table, &row);
  addConsColorLine(bc, "Low:", CONS_LEVEL_LOW, colorById ? &bc->lowIdCutoff : &bc->lowSimCutoff, table, &row);
  
  label = gtk_label_new("Press Enter or click Add to update the display after changing threshold values.");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_box_pack_start(vbox, label, FALSE, FALSE, DIALOG_YPAD);

  label = gtk_label_new("Click OK to save changes.");
  gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
  gtk_box_pack_start(vbox, label, FALSE, FALSE, DIALOG_YPAD);
}


/* Show a dialog to allow the user to edit the conservation color scheme */
static void showEditConsColorsDialog(GtkWidget *belvuWindow, const gboolean bringToFront)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);
  BelvuContext *bc = properties->bc;
  
  const BelvuDialogId dialogId = BELDIALOG_EDIT_CONS_COLORS;
  GtkWidget *dialog = getPersistentDialog(bc->dialogList, dialogId);
  
  if (!dialog)
    {
      char *title = g_strdup_printf("%sEdit Conservation Colors", belvuGetTitlePrefix(bc));

      dialog = gtk_dialog_new_with_buttons(title, 
                                           GTK_WINDOW(belvuWindow), 
                                           (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT | GTK_DIALOG_MODAL),
                                           GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                           GTK_STOCK_ADD, GTK_RESPONSE_APPLY,
                                           GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                           NULL);

      g_free(title);
      
      /* These calls are required to make the dialog persistent... */
      addPersistentDialog(bc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseConsColorsDialog), belvuWindow);
    }
  else
    {
      /* Need to refresh the dialog contents, so clear and re-create content area */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_APPLY);
  
  GtkBox *vbox = GTK_BOX(GTK_DIALOG(dialog)->vbox);
  createEditConsColorsContent(vbox, bc);
  
  /* Show / bring to front */
  gtk_widget_show_all(dialog);
  
  if (bringToFront)
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
  
}


/***********************************************************
 *                         Wrap window                     *
 ***********************************************************/

/* Create a text entry with a label.  'labelText' gives the label text and 
 * defaultInput gives the default input to show in the text entry (may be null).
 * Adds the result to table, if given, and returns the text entry widget */
static GtkWidget* createTextEntryWithLabel(const char *labelText, 
                                           const char *defaultInput, 
                                           GtkTable *table,
                                           const int col,
                                           const int row)
{
  const int xpad = 2;
  const int ypad = 2;
  
  /* Create the label in the given column */
  GtkWidget *label = gtk_label_new(labelText);
  gtk_misc_set_alignment(GTK_MISC(label), 1, 0);
  gtk_table_attach(table, label, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, xpad, ypad);
  
  /* Create the entry in the next column */
  GtkWidget *entry = gtk_entry_new();
  gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

  gtk_table_attach(table, entry, col + 1, col + 2, row, row + 1, (GtkAttachOptions)(GTK_EXPAND | GTK_FILL), GTK_SHRINK, xpad, ypad);
  
  if (defaultInput)
    {
      gtk_entry_set_text(GTK_ENTRY(entry), defaultInput);
      const int defaultLen = min(strlen(defaultInput) * 8, 500);
      gtk_widget_set_size_request(entry, defaultLen, -1);
    }
  
  return entry;
}


/* This shows a dialog that asks the user for settings for the wrap-alignment
 * view and then opens the wrap-alignment window on ok. */
static void showWrapDialog(BelvuContext *bc, GtkWidget *belvuWindow)
{
  char *title = g_strdup_printf("%swrap alignment", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  GTK_WINDOW(belvuWindow), 
                                                  (GtkDialogFlags)(GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT),
                                                  GTK_STOCK_OK, GTK_RESPONSE_ACCEPT,
                                                  GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;

  /* Create a text entry for the line width and title */
  GtkWidget *table = gtk_table_new(2, 2, FALSE);
  gtk_box_pack_start(GTK_BOX(contentArea), table, TRUE, TRUE, 0);
  
  GtkWidget *widthEntry = createTextEntryWithLabel("Line width", "80", GTK_TABLE(table), 0, 0);
  GtkWidget *titleEntry = createTextEntryWithLabel("Title", bc->Title, GTK_TABLE(table), 0, 1);
  
  gtk_widget_show_all(dialog);
  
  if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
      const gchar *inputText = gtk_entry_get_text(GTK_ENTRY(widthEntry));
      const int linelen = convertStringToInt(inputText);
      
      const gchar *title = gtk_entry_get_text(GTK_ENTRY(titleEntry));
      
      createWrapWindow(belvuWindow, linelen, title);
    }
  
  gtk_widget_destroy(dialog);
}


/* Utility function to return the child widget of the given container that
 * has a cached drawable. Returns the first found (i.e. expects there to be
 * only one).
 * Sets the result in the user data, if found. The value of *result must be
 * null on the first call. */
static void widgetGetDrawing(GtkWidget *widget, gpointer data)
{
  GtkWidget **result = (GtkWidget**)data;

  if (*result != NULL)
    {
      return;
    }
  else if (widgetGetDrawable(widget))
    {
      *result = widget;
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), widgetGetDrawing, result);
    }
}


static void setWrapWindowStyleProperties(GtkWidget *window)
{
  gtk_widget_set_name(window, WRAPPED_BELVU_WINDOW_NAME);
  
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);

  const int width = screenWidth * DEFAULT_WRAP_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_WRAP_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


static void createWrapWindow(GtkWidget *belvuWindow, const int linelen, const gchar *title)
{
  BelvuWindowProperties *properties = belvuWindowGetProperties(belvuWindow);

  /* Create the window */
  GtkWidget *wrapWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setWrapWindowStyleProperties(wrapWindow);
  
  BelvuContext *bc = windowGetContext(belvuWindow);
  char *windowTitle = g_strdup_printf("%s%s", belvuGetTitlePrefix(bc), properties->bc->Title);
  gtk_window_set_title(GTK_WINDOW(wrapWindow), windowTitle);
  g_free(windowTitle);
  
  /* We must add all toplevel windows to the list of spawned windows */
  properties->bc->spawnedWindows = g_slist_prepend(properties->bc->spawnedWindows, wrapWindow);
  
  /* Create the context menu and set a callback to show it */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(wrapWindow, properties->bc, &actionGroup);
  GtkWidget *contextmenu = createBelvuMenu(wrapWindow, "/WrapContextMenu", uiManager);
  
  gtk_widget_add_events(wrapWindow, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(wrapWindow), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(wrapWindow), vbox);
  
  /* Add the alignment section */
  GtkWidget *wrappedAlignment = createBelvuAlignment(properties->bc, title, linelen);
  gtk_box_pack_start(GTK_BOX(vbox), wrappedAlignment, TRUE, TRUE, 0);
  
  /* Set properties */
  genericWindowCreateProperties(wrapWindow, properties->bc, actionGroup);

  /* Show window */
  gtk_widget_show_all(wrapWindow);
  gtk_window_present(GTK_WINDOW(wrapWindow));
  
  /* Make sure the font size is up to date */
  onBelvuAlignmentFontSizeChanged(wrappedAlignment);
}


/***********************************************************
 *                      Organisms window                   *
 ***********************************************************/

static void setOrgsWindowStyleProperties(GtkWidget *window, BelvuContext *bc)
{
  gtk_widget_set_name(window, BELVU_ORGS_WINDOW_NAME);

  /* Set default size based on number of alignments and max name width */
  gdouble charWidth = 0, charHeight = 0;
  getFontCharSize(window, window->style->font_desc, &charWidth, &charHeight);

  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int maxWidth = screenWidth * MAX_ORGS_WINDOW_WIDTH_FRACTION;
  const int maxHeight = screenHeight * MAX_ORGS_WINDOW_HEIGHT_FRACTION;
  
  int width = min(maxWidth, charWidth * bc->maxNameLen + ORGS_WINDOW_XPAD * 2);
  int height = min(maxHeight, charHeight * bc->organismArr->len + ORGS_WINDOW_YPAD * 2);
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


/* This does the work to draw the organisms. */
static void drawOrganisms(GtkWidget *widget, GdkDrawable *drawable, BelvuContext *bc)
{
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor color;

  gdouble charHeight = 0;
  getFontCharSize(widget, widget->style->font_desc, NULL, &charHeight);

  int y = ORGS_WINDOW_YPAD;
  const int x = ORGS_WINDOW_XPAD;
  int i = 0;

  for ( ; i < (int)bc->organismArr->len; ++i) 
    {
      ALN *alnp = g_array_index(bc->organismArr, ALN*, i);
    
      convertColorNumToGdkColor(alnp->color, FALSE, &color);
      gdk_gc_set_foreground(gc, &color);
    
      drawText(widget, drawable, gc, x, y, alnp->organism, NULL, NULL);
    
      y += charHeight;
    }
 
  g_message_info("%d organisms found\n", bc->organismArr->len);
  g_object_unref(gc);
}


/* Expose handler for the organisms view */
static gboolean onExposeOrganismsView(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  BelvuContext *bc = (BelvuContext*)data;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
	{
	  /* There isn't a bitmap yet. Create it now. */
          guint width = 0, height = 0;
          gtk_layout_get_size(GTK_LAYOUT(widget), &width, &height);
	  bitmap = createBlankSizedPixmap(widget, window, width, height);
	  drawOrganisms(widget, bitmap, bc);
	}
      
      if (bitmap)
	{
	  /* Push the bitmap onto the window */
	  GdkGC *gc = gdk_gc_new(window);
	  gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
	  g_object_unref(gc);
	}
      else
	{
	  g_warning("Failed to draw Organisms view [%p] - could not create bitmap.\n", widget);
	}
    }
  
  return TRUE;
}


static void createOrganismWindow(BelvuContext *bc)
{
  /* Create the window */
  bc->orgsWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setOrgsWindowStyleProperties(bc->orgsWindow, bc);
  g_signal_connect(bc->orgsWindow, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
  
  char *title = g_strdup_printf("%sOrganisms", belvuGetTitlePrefix(bc));
  gtk_window_set_title(GTK_WINDOW(bc->orgsWindow), title);
  g_free(title);
  
  /* We must add all toplevel windows to the list of spawned windows */
  bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, bc->orgsWindow);
  
  /* Create the context menu and set a callback to show it */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(bc->orgsWindow, bc, &actionGroup);
  GtkWidget *contextmenu = createBelvuMenu(bc->orgsWindow, "/OrgsContextMenu", uiManager);
  
  gtk_widget_add_events(bc->orgsWindow, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(bc->orgsWindow), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* Create a scrollable drawing area */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(bc->orgsWindow), scrollWin);
  
  GtkWidget *drawing = gtk_layout_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrollWin), drawing);
  g_signal_connect(G_OBJECT(drawing), "expose-event", G_CALLBACK(onExposeOrganismsView), bc);

  /* Set layout size big enough to fit all organisms */
  gdouble width = 0, height = 0;
  getFontCharSize(drawing, drawing->style->font_desc, &width, &height);
  width *= (bc->maxNameLen + 1);
  height *= (bc->organismArr->len + 1);
  gtk_layout_set_size(GTK_LAYOUT(drawing), width, height);
  
  /* Set default background color */
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);
  gtk_widget_modify_bg(drawing, GTK_STATE_NORMAL, bgColor);
  
  /* Set properties */
  genericWindowCreateProperties(bc->orgsWindow, bc, actionGroup);
  
  gtk_widget_show_all(bc->orgsWindow);
  gtk_window_present(GTK_WINDOW(bc->orgsWindow));
}


/***********************************************************
 *                           Updates                       *
 ***********************************************************/

void showAnnotationWindow(BelvuContext *bc)
{
  /* If there are no annotations, there's nothing to do */
  if (g_slist_length(bc->annotationList) < 1) 
    return;
  
  /* Loop through each annotation line */
  GString *resultStr = g_string_new("");
  int maxLen = 0;
  GSList *annItem = bc->annotationList;
  
  for ( ; annItem; annItem = annItem->next)
    {
      /* Calculate the max line length */
      const char *cp = (const char*)(annItem->data);
      
      if ((int)strlen(cp) > maxLen) 
        maxLen = strlen(cp);
      
      /* Append this text to the string */
      g_string_append_printf(resultStr, "%s\n", cp);
    }
  
  /* Create the dialog */
  char *title = g_strdup_printf("%sAnnotations", belvuGetTitlePrefix(bc));

  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
                                                  NULL, 
                                                  GTK_DIALOG_DESTROY_WITH_PARENT,
                                                  GTK_STOCK_CLOSE, GTK_RESPONSE_ACCEPT,
                                                  NULL);

  g_free(title);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_CLOSE);
  g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(gtk_widget_destroy), NULL);
  
  /* Use a fixed-width font */
  const char *fontFamily = findFixedWidthFont(dialog);
  PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
  pango_font_description_set_size(fontDesc, pango_font_description_get_size(dialog->style->font_desc));
  
  GtkWidget *textView = createScrollableTextView(resultStr->str, FALSE, fontDesc, FALSE, NULL, NULL, NULL);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), textView, TRUE, TRUE, 0);

  const gchar *env = g_getenv(FONT_SIZE_ENV_VAR);
  if (env)
    widgetSetFontSizeAndCheck(dialog, convertStringToInt(env));


  /* Set the initial size */
  double charWidth, charHeight;
  getFontCharSize(dialog, fontDesc, &charWidth, &charHeight);
  int width = ((maxLen + 1) * charWidth) + scrollBarWidth();
  int height = ((g_slist_length(bc->annotationList) + 1) * charHeight) + scrollBarWidth() + 40; /* extra fudge to allow space for buttons */

  GdkScreen *screen = gtk_widget_get_screen(dialog);
  const int maxWidth = gdk_screen_get_width(screen) * MAX_ANNOTATION_WINDOW_WIDTH_FRACTION;
  const int maxHeight = gdk_screen_get_height(screen) * MAX_ANNOTATION_WINDOW_HEIGHT_FRACTION;
  width = min(width, maxWidth);
  height = min(height, maxHeight);
  
  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);

  gtk_widget_show_all(dialog);
}


/***********************************************************
 *                           Updates                       *
 ***********************************************************/


/* Update the feedback box to show info about the currently-selected items */
static void updateFeedbackBox(BelvuContext *bc, GtkWidget *feedbackBox)
{
  GString *resultStr = g_string_new("");
  
  char *tmpStr = NULL;
  
  /* If a column is selected, display the column number */
  if (bc->selectedCol > 0)
    {
      tmpStr = g_strdup_printf("Column %d: ", bc->selectedCol);
      g_string_append(resultStr, tmpStr);
      g_free(tmpStr);
    }
  
  /* If an alignment is selected, display info about it */
  if (bc->selectedAln && alnGetSeq(bc->selectedAln))
    {
      char *selectedSeq = alnGetSeq(bc->selectedAln);
      
      tmpStr = g_strdup_printf("%s/%d-%d", bc->selectedAln->name, bc->selectedAln->start, bc->selectedAln->end);
      g_string_append(resultStr, tmpStr);
      g_free(tmpStr);
      
      /* If a column is selected, display info about the selected alignment's
       * coord at that column position. */
      if (bc->selectedCol > 0)
        {
          /* Print the char of the current sequence at the selected column */
          tmpStr = g_strdup_printf("  %c = ", selectedSeq[bc->selectedCol - 1]);
          g_string_append(resultStr, tmpStr);
          g_free(tmpStr);
          
          /* Loop through each column before the selected column and calculate the
           * number of gaps. Also note whether we see an asterisk in the sequence */
          gboolean hasAsterisk = FALSE;
          int numGaps = 0;
          int colIdx = 0;
          
          for ( ; colIdx < bc->selectedCol; colIdx++)
            {
              if (isGap(selectedSeq[colIdx])) 
                numGaps++;
              else if (selectedSeq[colIdx] == '*') 
                hasAsterisk = TRUE;
            }
          
          if (hasAsterisk)
            {
              g_string_append(resultStr, "(unknown position due to insertion)");
            }
          else
            {
              tmpStr = g_strdup_printf("%d", bc->selectedCol - 1 + bc->selectedAln->start - numGaps);
              g_string_append(resultStr, tmpStr);
              g_free(tmpStr);
            }
        }
      
      /* Display the total number of highlighted alignments */
      const int numHighlighted = g_slist_length(bc->highlightedAlns);
      
      tmpStr = g_strdup_printf(" (%d match", numHighlighted);
      g_string_append(resultStr, tmpStr);
      g_free(tmpStr);
      
      if (numHighlighted != 1)
        g_string_append(resultStr, "es");
      
      g_string_append(resultStr, ")");
    }
  
  gtk_entry_set_text(GTK_ENTRY(feedbackBox), resultStr->str);
  
  g_string_free(resultStr, TRUE);
}


/* This should be called whenever the selected sequence has changed */
void onRowSelectionChanged(BelvuContext *bc)
{
  /* Redraw the alignment widget */
  belvuAlignmentRedrawAll(bc->belvuAlignment);
  centerHighlighted(bc, bc->belvuAlignment);
  
  /* Redraw all of the trees */
  g_slist_foreach(bc->spawnedWindows, belvuTreeRedrawAll, NULL);
  
  /* Set the status of the 'exclude highlighted' toggle menu option
   * depending on whether the newly-selected sequence is selected or not
   * (or grey it out if nothing is selected) */
  BelvuWindowProperties *properties = belvuWindowGetProperties(bc->belvuWindow);
  enableMenuAction(properties->actionGroup, "excludeHighlighted", bc->selectedAln != NULL);
  
  if (bc->selectedAln)
    {
      setToggleMenuStatus(properties->actionGroup, "excludeHighlighted", bc->selectedAln->nocolor);

      /* Copy the selected sequence name to the PRIMARY clipboard */
      setPrimaryClipboardText(bc->selectedAln->name);
    }

  /* Highlight any alignments that have the same name as the selected alignment */
  g_slist_free(bc->highlightedAlns); /* clear current list */
  bc->highlightedAlns = NULL;

  int i = 0;
  for (i = 0; i < (int)bc->alignArr->len; ++i)
    {
      ALN *alnp = g_array_index(bc->alignArr, ALN*, i);

      if (alignmentHighlighted(bc, alnp))
        bc->highlightedAlns = g_slist_prepend(bc->highlightedAlns, alnp);
    }

  /* Update the feedback box */
  updateFeedbackBox(properties->bc, properties->feedbackBox);
  
  /* If the current sort method is by similarity/id to the selected sequence,
   * this has effectively been 'invalidated' by the fact that the selected
   * sequence has changed. We don't undo the sort, but we must set the sort
   * method to 'unsorted' so that the user can re-select the sort-by-sim/sort-by-id
   * radio button if they want to sort by similarity/id to the newly-selected sequence */
  setRadioMenuStatus(properties->actionGroup, "unsorted", BELVU_UNSORTED);
}


/* This should be called whenever the selected column has changed */
void onColSelectionChanged(BelvuContext *bc)
{
  /* Update the feedback box */
  BelvuWindowProperties *properties = belvuWindowGetProperties(bc->belvuWindow);
  updateFeedbackBox(properties->bc, properties->feedbackBox);
  
  /* Refresh the alignment widget */
  belvuAlignmentRefreshAll(bc->belvuAlignment);
}


/* This is called after a tree has changed */
void onTreeOrderChanged(BelvuContext *bc)
{
  /* If sorting by tree order, we need to refresh the sort order */
  if (bc->sortType == BELVU_SORT_TREE)
    doSort(bc, bc->sortType, FALSE);

  /* Recenter on the highlighted alignment */
  centerHighlighted(bc, bc->belvuAlignment);

  /* Redraw the tree and the alignment list */
  belvuTreeRedrawAll(bc->belvuTree, NULL);
  belvuAlignmentRedrawAll(bc->belvuAlignment);
}


/***********************************************************
 *                         Key handlers                    *
 ***********************************************************/

static gboolean onKeyPressEscape(BelvuContext *bc)
{
  if (bc->removingSeqs)
    {
      /* Cancel 'removing sequences' mode */
      BelvuWindowProperties *properties = belvuWindowGetProperties(bc->belvuWindow);
      setToggleMenuStatus(properties->actionGroup, "rmMany", !bc->removingSeqs);
    }
  
  return TRUE;
}

static gboolean onKeyPressHomeEnd(BelvuContext *bc, const gboolean home, const gboolean ctrl, const gboolean shift)
{
  /* Scroll to the top/bottom of the alignment list */
  vScrollStartEnd(bc->belvuAlignment, home);
  return TRUE;
}

static gboolean onKeyPressPageUpDown(BelvuContext *bc, const gboolean up, const gboolean ctrl, const gboolean shift)
{
  /* Scroll to the top/bottom of the alignment list */
  vScrollPageUpDown(bc->belvuAlignment, up);
  return TRUE;
}

static gboolean onKeyPressLeftRight(BelvuContext *bc, const gboolean left, const gboolean ctrl, const gboolean shift)
{
  if (ctrl)
    {
      /* Scroll one character left/right */
      hScrollLeftRight(bc->belvuAlignment, left, 1);
    }
  else
    {
      /* Scroll to the top/bottom of the alignment list */
      hScrollPageLeftRight(bc->belvuAlignment, left);
    }
  
  return TRUE;
}

static gboolean onKeyPressUpDown(BelvuContext *bc, const gboolean up, const gboolean ctrl, const gboolean shift)
{
  if (ctrl)
    {
      /* Scroll one row up/down */
      vScrollUpDown(bc->belvuAlignment, up, 1);
    }
  else
    {
      /* Scroll one page up/down */
      vScrollPageUpDown(bc->belvuAlignment, up);
    }
  
  return TRUE;
}

static gboolean onKeyPressCommaPeriod(BelvuContext *bc, const gboolean comma, const gboolean ctrl, const gboolean shift)
{
  if (ctrl)
    {
      if (shift) /* scroll to very start/end */
        hScrollStartEnd(bc->belvuAlignment, comma);
      else
        hScrollPageLeftRight(bc->belvuAlignment, comma); /* Scroll left/right by one page */
    }
  else
    {
      /* Scroll left/right by one character */
      hScrollLeftRight(bc->belvuAlignment, comma, 1);
    }
  
  return TRUE;
}

static gboolean onKeyPressInsDel(BelvuContext *bc, const gboolean insert, const gboolean ctrl, const gboolean shift)
{
  /* Scroll to the leftmost/rightmost extent of the display */
  hScrollStartEnd(bc->belvuAlignment, insert);
  return TRUE;
}

static gboolean onKeyPressPlusMinus(BelvuContext *bc, const gboolean plus, const gboolean ctrl, const gboolean shift)
{
  /* Zoom in/out */
  if (plus)
    incrementFontSize(bc);
  else
    decrementFontSize(bc);
    
  return TRUE;
}


/***********************************************************
 *                           Events                        *
 ***********************************************************/

/* Key press handler */
gboolean onKeyPressBelvu(GtkWidget *window, GdkEventKey *event, gpointer data)
{
  gboolean handled = FALSE;
  
  BelvuContext *bc = (BelvuContext*)data;

  const gboolean ctrl = (event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK;	
  const gboolean shift = (event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK;	
  
  switch (event->keyval)
    {
      case GDK_Escape:    handled = onKeyPressEscape(bc);                          break;
        
      case GDK_Home:      handled = onKeyPressHomeEnd(bc, TRUE, ctrl, shift);      break;
      case GDK_End:       handled = onKeyPressHomeEnd(bc, FALSE, ctrl, shift);     break;
      case GDK_Page_Up:   handled = onKeyPressPageUpDown(bc, TRUE, ctrl, shift);   break;
      case GDK_Page_Down: handled = onKeyPressPageUpDown(bc, FALSE, ctrl, shift);  break;
        
      case GDK_Left:      handled = onKeyPressLeftRight(bc, TRUE, ctrl, shift);    break;
      case GDK_Right:     handled = onKeyPressLeftRight(bc, FALSE, ctrl, shift);   break;
      case GDK_Up:        handled = onKeyPressUpDown(bc, TRUE, ctrl, shift);       break;
      case GDK_Down:      handled = onKeyPressUpDown(bc, FALSE, ctrl, shift);      break;
        
      case GDK_less:      /* fall through */
      case GDK_comma:     handled = onKeyPressCommaPeriod(bc, TRUE, ctrl, shift);  break;
      case GDK_greater:   /* fall through */
      case GDK_period:    handled = onKeyPressCommaPeriod(bc, FALSE, ctrl, shift); break;
        
      case GDK_Insert:    handled = onKeyPressInsDel(bc, TRUE, ctrl, shift);       break;
      case GDK_Delete:    handled = onKeyPressInsDel(bc, FALSE, ctrl, shift);      break;

      case GDK_plus:      /* fall through */
      case GDK_equal:     handled = onKeyPressPlusMinus(bc, TRUE, ctrl, shift);    break;
      case GDK_minus:     /* fall through */
      case GDK_underscore: handled = onKeyPressPlusMinus(bc, FALSE, ctrl, shift);  break;

      default: break;
    };
  
  return handled;
}


/* Mouse button handler */
gboolean onButtonPressBelvu(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      /* For the main window, if we're removing sequences, then just cancel that mode. 
       * Otherwise (and for any other window type) show the context menu. */
      if (stringsEqual(gtk_widget_get_name(window), MAIN_BELVU_WINDOW_NAME, TRUE))
	{
	  BelvuWindowProperties *properties = belvuWindowGetProperties(window);
	  if (properties->bc->removingSeqs)
	    {
              setToggleMenuStatus(properties->actionGroup, "rmMany", !properties->bc->removingSeqs);
	      handled = TRUE;
	    }
	}

      if (!handled)
	{
	  GtkMenu *contextMenu = GTK_MENU(data);
	  gtk_menu_popup (contextMenu, NULL, NULL, NULL, NULL, event->button, event->time);
	  handled = TRUE;
	}
    }
  
  return handled;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

static GtkWidget* createFeedbackBox(GtkToolbar *toolbar)
{
  /* Bit of a hack, but add some space before the feedback box, just
   * to avoid the toolbar being too cluttered. */
  addToolbarWidget(toolbar, gtk_label_new("    "), -1);
  
  GtkWidget *feedbackBox = gtk_entry_new();
  
  /* User can copy text out but not edit contents */
  gtk_editable_set_editable(GTK_EDITABLE(feedbackBox), FALSE);
  
  GtkToolItem *item = addToolbarWidget(toolbar, feedbackBox, -1);
  gtk_tool_item_set_expand(item, TRUE); 
  
  /* We want the box to be printed, so connect the expose function that will 
   * draw to a pixmap for printing */
  g_signal_connect(G_OBJECT(feedbackBox), "expose-event", G_CALLBACK(onExposePrintable), NULL);
  
  return feedbackBox;
}


/* Set various properties for the main belvu window components */
static void setStyleProperties(GtkWidget *window, GtkToolbar *toolbar)
{
  /* Set the initial window size based on some fraction of the screen size */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_BELVU_WINDOW_WIDTH_FRACTION;
  const int height = screenHeight * DEFAULT_BELVU_WINDOW_HEIGHT_FRACTION;
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  gtk_container_set_border_width (GTK_CONTAINER(window), DEFAULT_WINDOW_BORDER_WIDTH); 
  gtk_window_set_mnemonic_modifier(GTK_WINDOW(window), GDK_MOD1_MASK); /* MOD1 is ALT on most systems */
  
    /* Set toolbar style properties */
  gtk_toolbar_set_style(toolbar, GTK_TOOLBAR_ICONS);
  gtk_toolbar_set_icon_size(toolbar, GTK_ICON_SIZE_SMALL_TOOLBAR);
}


gboolean createBelvuWindow(BelvuContext *bc, BlxMessageData *msgData)
{
  gboolean ok = TRUE;
  
  GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_widget_set_name(window, MAIN_BELVU_WINDOW_NAME);

  /* Set a pointer to the main window in the context */
  bc->belvuWindow = window;

  /* Set the title */
  char *title = g_strdup_printf("%s%s", belvuGetTitlePrefix(bc), bc->Title);
  gtk_window_set_title(GTK_WINDOW(window), title);
  g_free(title);
  
  /* Create the status bar */
  GtkWidget *statusBar = gtk_statusbar_new();
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), TRUE);
  setStatusBarShadowStyle(statusBar, "GTK_SHADOW_NONE");
  gtk_statusbar_set_has_resize_grip(GTK_STATUSBAR(statusBar), FALSE);

  /* Set the window and statusbar in the message handler data, now that we know them */
  msgData->parent = GTK_WINDOW(window);
  msgData->statusBar = GTK_STATUSBAR(statusBar);

  /* Create the menu and toolbar */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(window, bc, &actionGroup);
  
  GtkWidget *menubar = createBelvuMenu(window, "/MenuBar", uiManager);
  GtkWidget *contextmenu = createBelvuMenu(window, "/ContextMenu", uiManager);
  GtkWidget *toolbar = createBelvuMenu(window, "/Toolbar", uiManager);
  
  addToolbarWidget(GTK_TOOLBAR(toolbar), gtk_label_new("    "), 0); /* hacky way to add some space at start of toolbar */

  /* Create the feedback box on the toolbar */
  GtkWidget *feedbackBox = createFeedbackBox(GTK_TOOLBAR(toolbar));
  
  /* Set the style properties */
  setStyleProperties(window, GTK_TOOLBAR(toolbar));

  /* Create the alignment section. Store it in the context so that we can update it. */
  bc->belvuAlignment = createBelvuAlignment(bc, NULL, UNSET_INT);
  
  /* We'll put everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(window), GTK_WIDGET(vbox));

  /* Put the menu and toolbar on the same row using an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), menubar, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(hbox), toolbar, TRUE, TRUE, 0);
  
  gtk_box_pack_start(GTK_BOX(vbox), bc->belvuAlignment, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), statusBar, FALSE, FALSE, 0);

  /* Connect signals */
  gtk_widget_add_events(window, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(window, GDK_KEY_PRESS_MASK);
  g_signal_connect(G_OBJECT(window), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  g_signal_connect(G_OBJECT(window), "key-press-event", G_CALLBACK(onKeyPressBelvu), bc);

  belvuWindowCreateProperties(window, bc, statusBar, feedbackBox, actionGroup);

  /* Set the default cursor (can only get the window's cursor after window is shown) */
  bc->defaultCursor = NULL;

  /* Show the main window (unless we only want the tree) */
  if (!bc->onlyTree)
    gtk_widget_show_all(window);
  
  /* If the BELVU_FONT_SIZE environment variable is set, we'll use it to set the
   * default font size for all the widgets. */
  const gchar *env = g_getenv(FONT_SIZE_ENV_VAR);
  
  if (env)
    {
      if (bc->belvuWindow)
        widgetSetFontSizeAndCheck(bc->belvuWindow, convertStringToInt(env));

      if (bc->belvuTree)
        widgetSetFontSizeAndCheck(bc->belvuTree, convertStringToInt(env));
    }

  env = g_getenv(STATUSBAR_SIZE_ENV_VAR);
  if (env)
    {
      const int height = convertStringToInt(env);
      
      /* If too small, hide the statusbar */
      if (height < MIN_FONT_SIZE)
        gtk_widget_hide_all(statusBar);
      else
        widgetSetFontSizeAndCheck(statusBar, height);
    }

  /* Make sure the alignment font size isup to date. Note: do this before
   * creating the tree, because the tree flushes all pending gtk calls and 
   * we must update the font size before the alignment is realised. */
  onBelvuAlignmentFontSizeChanged(bc->belvuAlignment);

  /* Show the tree on startup, if requested */
  if (bc->initTree)
    {
      createAndShowBelvuTree(bc, TRUE);
      onBelvuTreeFontSizeChanged(bc->belvuTree);
    }
  
  if (!bc->onlyTree)
    {
      gtk_window_present(GTK_WINDOW(window));
      
      if (bc->sortType)
        setRadioMenuStatus(actionGroup, "unsorted", bc->sortType);
      
      setToggleMenuStatus(actionGroup, "displayColors", bc->displayColors);

      if (bc->schemeType == BELVU_SCHEME_TYPE_RESIDUE)
        setRadioMenuStatus(actionGroup, "colorSchemeStandard", bc->residueScheme);
      else
        setRadioMenuStatus(actionGroup, "colorSchemeStandard", bc->consScheme);
      
      if (bc->initTree)
        belvuAlignmentRedrawAll(bc->belvuAlignment); /* redraw, because tree creation removes markup which can mess this up */
    }
  
  return ok;
}
