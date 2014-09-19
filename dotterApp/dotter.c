/*  File: dotplot.c
 *  Author: Erik Sonnhammer, 1993-09-04
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
 * Description: This file contains functions that initialise Dotter and create
 *              the main Dotter window. 
 *----------------------------------------------------------------------------
 */


/*
   DOTTER - Sequence-sequence dotplot in a pixel background image

   Memory requirements:
   The pixelmap + (alphabetsize+4) x qlen + 2 x slen

          Pending (this is an old list - may be out of date):

          change filqueryopen to filqueryeditopen (.dotter) when getting new code.

          Fix revcomp when zooming in.

          Add to alignment tool the extent and score of the current window.

          (Raise initial cutoffs when compressed. ?)

          HSP drawing bugged for reverse matches and for reversed scale
          (Huge job to fix ... do only if really necessary)

          Score matrix for DNA w/ transversions/transitions...? literature?

-------------------------------------------------------------------------------- */

/* CPU usage profiling on SGI:

   cc dotter.c, no optimization
   pixie -o dotter.pixie dotter
   dotter.pixie -wD -b t seq seq
   prof -pixie -h -only calcWindow dotter dotter.Addrs dotter.Counts
*/

#include <seqtoolsUtils/version.h>
#include <seqtoolsUtils/utilities.h>
#include <dotterApp/dotter_.h>
#include <dotterApp/seqtoolsExonView.h>
#include <gdk/gdkkeysyms.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* tint stuff used to be in graph.h, now local - rd 960524
   NB colours as Erik liked them before Jean tinkered!
   could rename #define's more sensibly now
*/

//#define TINT_WHITE      0x00
//#define TINT_HIGHLIGHT1 0x01  /* highest priority, dna highlighting */
//#define TINT_HIGHLIGHT2 0x02  /* highlight friends */
//#define TINT_RED        0x04 
//#define TINT_LIGHTGRAY  0x08 
//#define TINT_MAGENTA    0x10 
//#define TINT_CYAN       0x20 
//#define TINT_LIGHTGREEN 0x40 
//#define TINT_YELLOW     0x80 

//static int tints[8] = { LIGHTRED, MIDBLUE, RED, LIGHTGRAY, 
//                      MAGENTA, CYAN, LIGHTGREEN, YELLOW } ;

#define MAX_WINDOW_WIDTH_FRACTION             0.8 /* max init width of dotter window as fraction of screen size */
#define MAX_WINDOW_HEIGHT_FRACTION            0.8 /* max init height of dotter window as fraction of screen size */
#define MAIN_WINDOW_NAME                      "DotterMainWindow"

#define DOCK_WINDOWS_DEFAULT TRUE

//typedef struct
//{
//  const char *name ;
//  int start, end ;
//  char strand ;
//  MSP *msp_start, *msp_end ;
//  float y_pos ;
//} GeneDataStruct, *GeneData ;
//
//typedef struct
//{
//  GList *forward_genes ;
//  GList *reverse_genes ;
//  char strand  ;
//} GeneStrandStruct, *GeneStrand ;
//


/* Struct containing properties for a dotter window */
typedef struct _DotterProperties
{
  GtkWidget *greyrampTool;                  /* the greyramp tool */
  GtkWidget *greyrampWindow;                /* the window containing the greyramp too when undocked */
  GtkWidget *greyrampContainer;             /* the container containing the greyramp tool when docked */
  GtkWidget *alignmentTool;                 /* the alignment tool */
  GtkWidget *alignmentWindow;               /* the window containing the alignment tool when undocked */
  GtkWidget *alignmentContainer;            /* the container containing the alignment tool when docked */
  GtkWidget *dotplot;                       /* the dotplot drawing area */
  
  gboolean windowsDocked;                   /* if true, all tools are docked into a single window */
  DotterWindowContext *dotterWinCtx;
  const char *exportFileName;
} DotterProperties;


/* Local function declarations */

static void showSettingsDialog(GtkWidget *dotterWindow);
static void readmtx(int mtx[24][24], char *mtxfile);
static void mtxcpy(int mtx[24][24], int BLOSUM62[24][24]);
static void DNAmatrix(int mtx[24][24]);
//static void drawAllFeatures(MSP *msp) ;
//static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth) ;
//gint compareMSPs(gconstpointer a, gconstpointer b) ;
//static void getGenePositionsCB(gpointer data, gpointer user_data) ;
//gint compareGenes(gconstpointer a, gconstpointer b) ;
//static void setYoffsets(GList *first, float min_offset, float max_offset) ;
//static void drawGenesCB(gpointer data, gpointer user_data) ;
//static void drawMSPGene(MSP *msp, float y_offset) ;
//static int gArrayGetLen(GArray *array);

static void                   showHideGreyrampTool(GtkWidget *dotterWindow, const gboolean show);
static void                   showHideAlignmentTool(GtkWidget *dotterWindow, const gboolean show);
static GtkWidget*             createDotterWindow(DotterContext *dc, DotterWindowContext *dwc, const DotterHspMode hspMode, GtkWidget *dotplot, GtkWidget *dotplotContainer, GtkWidget *greyrampContainer, GtkWidget *alignmentContainer, const char *exportFileName, char *windowColor);
static DotterContext*         dotterGetContext(GtkWidget *dotterWindow);
static void                   redrawAll(GtkWidget *dotterWindow, gpointer data);
static void                   refreshAll(GtkWidget *dotterWindow, gpointer data);
static gboolean               onKeyPressDotter(GtkWidget *widget, GdkEventKey *event, gpointer data);
static gboolean               onKeyPressDotterCoords(GtkWidget *widget, GdkEventKey *event, gpointer data);
static gboolean               negateDisplayCoord(DotterContext *dc, const gboolean horizontal);
static gboolean               setStartCoord(GtkWidget *dotterWindow, DotterWindowContext *dwc, const gboolean horizontal, const int newValue);
static gboolean               setEndCoord(GtkWidget *dotterWindow, DotterWindowContext *dwc, const gboolean horizontal, const int newValue);
static void                   printDotterWindow(GtkWidget *dotterWindow);
static DotterProperties*      dotterGetProperties(GtkWidget *dotterWindow);

static GtkWidget* createDotterInstance(DotterContext *dotterCtx,
                                       DotterWindowContext *dotterWinCtx,
                                       const char *loadFileName,
                                       const char *saveFileName,
                                       const char *exportFileName,
                                       const gboolean hspsOn,
                                       const gboolean breaklinesOn,
                                       const char* winsizeIn,
                                       const int pixelFacIn,
                                       const int zoomFacIn,
                                       const int qcenter,
                                       const int scenter,
                                       const gboolean greyramSwap,
                                       char *windowColor);

static void                       onQuitMenu(GtkAction *action, gpointer data);
static void                       onCloseMenu(GtkAction *action, gpointer data);
static void                       onSavePlotMenu(GtkAction *action, gpointer data);
static void                       onExportPlotMenu(GtkAction *action, gpointer data);
static void                       onPrintMenu(GtkAction *action, gpointer data);
static void                       onSettingsMenu(GtkAction *action, gpointer data);
static void                       onShowHideGreyrampMenu(GtkAction *action, gpointer data);
static void                       onShowHideAlignmentMenu(GtkAction *action, gpointer data);
static void                       onToggleCrosshairMenu(GtkAction *action, gpointer data);
static void                       onToggleCoordsMenu(GtkAction *action, gpointer data);
static void                       onToggleFullscreenMenu(GtkAction *action, gpointer data);
static void                       onTogglePixelmapMenu(GtkAction *action, gpointer data);
static void                       onToggleGridMenu(GtkAction *action, gpointer data);
static void                       onHelpMenu(GtkAction *action, gpointer data);
static void                       onAboutMenu(GtkAction *action, gpointer data);
static void                       onCopyHCoordMenu(GtkAction *action, gpointer data);
static void                       onCopyVCoordMenu(GtkAction *action, gpointer data);
static void                       onToggleUsePrintColorsMenu(GtkAction *action, gpointer data);
static void                       onToggleBumpExonsMenu(GtkAction *action, gpointer data);
static void                       onToggleDockWindowsMenu(GtkAction *action, gpointer data);
static void                       onPrintColorsChanged(GtkWidget *dotterWindow);


/* Menu builders: the action entry list lists menu actions for all menus */
static const GtkActionEntry menuEntries[] = {
{ "FileMenuAction",   NULL, "_File"},
{ "EditMenuAction",   NULL, "_Edit"},
{ "ViewMenuAction",   NULL, "_View"},
{ "HelpMenuAction",   NULL, "_Help"},

{ "Quit",           GTK_STOCK_QUIT,         "_Quit all dotters",      "<control>Q", "Quit all dotters",           G_CALLBACK(onQuitMenu)},
{ "Close",          GTK_STOCK_CLOSE,        "_Close this dotter",     "<control>W", "Close this dotter",          G_CALLBACK(onCloseMenu)},
{ "SavePlot",       GTK_STOCK_SAVE,         "_Save plot",             NULL,         "Save plot",                  G_CALLBACK(onSavePlotMenu)},
{ "ExportPlot",     NULL,                   "_Export plot",           NULL,         "Export plot",                G_CALLBACK(onExportPlotMenu)},
{ "Print",          GTK_STOCK_PRINT,        "_Print...",              "<control>P", "Print",                      G_CALLBACK(onPrintMenu)},
{ "Settings",       GTK_STOCK_PREFERENCES,  "Settings",               "<control>S", "Set dotter parameters",      G_CALLBACK(onSettingsMenu)},
{ "Help",           GTK_STOCK_HELP,         "_Help",                  "<control>H", "Dotter Help",                G_CALLBACK(onHelpMenu)},
{ "About",          GTK_STOCK_ABOUT,        "_About",                 NULL,         "About Dotter",               G_CALLBACK(onAboutMenu)},
{ "CopyHCoord",     NULL,                   "Copy _horizontal coord", NULL,         "Copy the current horizontal sequence coord to the clipboard", G_CALLBACK(onCopyHCoordMenu)},
{ "CopyVCoord",     NULL,                   "Copy _vertical coord",   NULL,         "Copy the current vertical sequence coord to the clipboard", G_CALLBACK(onCopyVCoordMenu)}
};

/* Toggle-able menu entries are listed here: */
static GtkToggleActionEntry toggleMenuEntries[] = {
{ "TogglePixmap",     NULL, "Pixelmap",              NULL,  "Show the pixelmap",              G_CALLBACK(onTogglePixelmapMenu),        TRUE}, /* must be item 0 in list */
{ "ToggleGrid",       NULL, "Gridlines",             NULL,  "Show grid lines",                G_CALLBACK(onToggleGridMenu),            FALSE},
{ "ToggleCrosshair",  NULL, "Crosshair",             NULL,  "Show the crosshair",             G_CALLBACK(onToggleCrosshairMenu),       TRUE},
{ "ToggleCoords",     NULL, "Crosshair label",       NULL,  "Show the crosshair label",       G_CALLBACK(onToggleCoordsMenu),          TRUE},
{ "ToggleFullscreen", NULL, "Crosshair fullscreen",  NULL,  "Show the crosshair full screen", G_CALLBACK(onToggleFullscreenMenu),      TRUE},
{ "TogglePrintColors",NULL, "Use print colors",      NULL,  "Use print _colors",              G_CALLBACK(onToggleUsePrintColorsMenu),  FALSE},
{ "ToggleBumpExons",  NULL, "Bump exons",            "B",   "_Bump exons",                    G_CALLBACK(onToggleBumpExonsMenu),       FALSE},
{ "ToggleGreyramp",   NULL, "_Greyramp tool",        "<control>G", "Show/hide the greyramp tool",G_CALLBACK(onShowHideGreyrampMenu),   TRUE},
{ "ToggleAlignment",  NULL, "_Alignment tool",       "<control>A", "Show/hide the alignment tool",G_CALLBACK(onShowHideAlignmentMenu), TRUE},
{ "DockWindows",      NULL, "Dock windows",          "<control>K", "_Dock windows",           G_CALLBACK(onToggleDockWindowsMenu),     DOCK_WINDOWS_DEFAULT}
};

/* Radio-button menu entries are listed here: */
static const GtkRadioActionEntry radioMenuEntries[] = {
{ "HspsOff",    NULL, "HSPs off",                     NULL,  "Hide Blast HSPs",                                    DOTTER_HSPS_OFF},
{ "HspsGrey",   NULL, "Draw HSPs (greyramp)",         NULL,  "Draw Blast HSPs as greyramp",                        DOTTER_HSPS_GREYSCALE},
{ "HspsLine",   NULL, "Draw HSPs (red lines)",        NULL,  "Draw Blast HSPs as solid red lines",                 DOTTER_HSPS_LINE},
{ "HspsFunc",   NULL, "Draw HSPs (color = f(score))", NULL,  "Draw Blast HSPs with color as a function of score",  DOTTER_HSPS_FUNC}
};



/* Menu UI descriptions */
static const char mainMenuDescription[] =
"<ui>"
"  <menubar name='MenuBar'>"
"    <menu action='FileMenuAction'>"
"      <menuitem action='SavePlot'/>"
"      <menuitem action='ExportPlot'/>"
"      <separator/>"
"      <menuitem action='Print'/>"
"      <menuitem action='TogglePrintColors'/>"
"      <separator/>"
"      <menuitem action='Close'/>"
"      <menuitem action='Quit'/>"
"    </menu>"
"    <menu action='EditMenuAction'>"
"      <menuitem action='CopyHCoord'/>"
"      <menuitem action='CopyVCoord'/>"
"      <separator/>"
"      <menuitem action='Settings'/>"
"    </menu>"
"    <menu action='ViewMenuAction'>"
"      <menuitem action='ToggleGreyramp'/>"
"      <menuitem action='ToggleAlignment'/>"
"      <menuitem action='DockWindows'/>"
"      <separator/>"
"      <menuitem action='ToggleCrosshair'/>"
"      <menuitem action='ToggleCoords'/>"
"      <menuitem action='ToggleFullscreen'/>"
"      <separator/>"
"      <menuitem action='TogglePixmap'/>"
"      <menuitem action='ToggleGrid'/>"
"      <separator/>"
"      <menuitem action='HspsOff'/>"
"      <menuitem action='HspsGrey'/>"
"      <menuitem action='HspsLine'/>"
"      <menuitem action='HspsFunc'/>"
"     <separator/>"
"      <menuitem action='ToggleBumpExons'/>"
"    </menu>"
"    <menu action='HelpMenuAction'>"
"      <menuitem action='Help'/>"
"      <menuitem action='About'/>"
"    </menu>"
"  </menubar>"
"  <popup name='ContextMenu' accelerators='true'>"
"    <menuitem action='Close'/>"
"    <menuitem action='Quit'/>"
"    <menuitem action='Help'/>"
"    <menuitem action='SavePlot'/>"
"    <menuitem action='Print'/>"
"    <separator/>"
"    <menuitem action='Settings'/>"
"    <separator/>"
"    <menuitem action='ToggleGreyramp'/>"
"    <menuitem action='ToggleAlignment'/>"
"    <menuitem action='DockWindows'/>"
"    <separator/>"
"    <menu action='ViewMenuAction'>"
"      <separator/>"
"      <menuitem action='ToggleCrosshair'/>"
"      <menuitem action='ToggleCoords'/>"
"      <menuitem action='ToggleFullscreen'/>"
"      <separator/>"
"      <menuitem action='TogglePixmap'/>"
"      <menuitem action='ToggleGrid'/>"
"      <separator/>"
"      <menuitem action='HspsOff'/>"
"      <menuitem action='HspsGrey'/>"
"      <menuitem action='HspsLine'/>"
"      <menuitem action='HspsFunc'/>"
"    </menu>"
"  </popup>"
"</ui>";



/* Global variables */

static int    MATRIX[24][24];
//              MSPoffset,      /* Difference between real MSP coord and coord stored in MSP */
//              HSPgaps = 0.
//              fsBoxStart,
//              fsRightOn = 1,
//              fsBottomOn = 1,
//              fsAnnRightOn = 1,
//              fsAnnBottomOn = 1,
//              fsEndLinesOn = 0;

//static Graph  fsGraph=0;

#define MAXALIGNLEN 501

//static char fsPlotHeighttx[10];

       float fsPlotHeight=2.0;  /* The height of feature series XY plots */
//static MSP   *msp;

/*  BLOSUM62 930809

     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
int BLOSUM62[24][24] = {
  {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4 },
  { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4 },
  { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4 },
  { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
  {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 },
  { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4 },
  { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4 },
  { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4 },
  { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4 },
  { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4 },
  { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4 },
  { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4 },
  { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4 },
  { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4 },
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4 },
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4 },
  { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4 },
  { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4 },
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4 },
  { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
  { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
  { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4 },
  { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1 }
 };


//GArray *fsArr = NULL;  /* Stores Feature Series - t he actual segments are stored
//                         as MSPs, using these fields:
//                         msp->sframe  = [1..2] sequence
//                         msp->qstart = segment start
//                         msp->qend   = segment end
//                         msp->fs     = Ordinal number of series that this MSP belongs to.
//                         msp->fsColor  = color
//                         msp->desc   = annotation
//                         */


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

/* free memory used by color structs */
static void destroyDotterColors(DotterContext *dc)
{
  if (!dc || !dc->defaultColors)
    return;

  int i = DOTCOLOR_MIN + 1;
  
  for ( ; i < DOTCOLOR_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = &g_array_index(dc->defaultColors, BlxColor, i);
      destroyBlxColor(blxColor);
    }

  g_array_free(dc->defaultColors, FALSE);
  dc->defaultColors = NULL;
}


/* Create the colors that Dotter will use for various specific purposes */
static void createDotterColors(DotterContext *dc,
                               DotterOptions *options)
{
  /* Initialise the array with empty BlxColor structs */
  dc->defaultColors = g_array_sized_new(FALSE, FALSE, sizeof(BlxColor), DOTCOLOR_NUM_COLORS);
  int i = DOTCOLOR_MIN + 1;

  for ( ; i < DOTCOLOR_NUM_COLORS; ++i)
    {
      BlxColor *blxColor = (BlxColor*)g_malloc(sizeof(BlxColor));
      blxColor->name = NULL;
      blxColor->desc = NULL;
      g_array_append_val(dc->defaultColors, *blxColor);
    }

  /* Set default colors that can be set by user options */
  const char* breaklineColor = options->breaklineColor ? options->breaklineColor : BLX_GREEN;
  
  /* Get the default background color of our widgets (i.e. that inherited from the theme).
   * Convert it to a string so we can use the same creation function as the other colors */
  GtkWidget *tmp = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  char *defaultBgColorStr = convertColorToString(&tmp->style->bg[GTK_STATE_NORMAL]);
  createBlxColor(dc->defaultColors, DOTCOLOR_BACKGROUND, "Background", "Background color", defaultBgColorStr, BLX_WHITE, "#bdbdbd", NULL);
  g_free(defaultBgColorStr);
  gtk_widget_destroy(tmp);
  
  /* matches */
  createBlxColor(dc->defaultColors, DOTCOLOR_MATCH, "Exact match", "Exact match", BLX_LIGHT_CYAN, BLX_LIGHT_CYAN, BLX_CYAN, BLX_CYAN);
  createBlxColor(dc->defaultColors, DOTCOLOR_CONS, "Conserved match", "Conserved match", BLX_VIOLET, BLX_VIOLET, BLX_DARK_VIOLET, BLX_DARK_VIOLET);
  createBlxColor(dc->defaultColors, DOTCOLOR_MISMATCH, "Mismatch", "Mismatch", "#cacaca", BLX_WHITE, "#cacaca", BLX_WHITE);
  
  /* exons */
  createBlxColor(dc->defaultColors, DOTCOLOR_EXON_FILL, "Exon fill color", "Exon fill color", BLX_YELLOW, BLX_YELLOW, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_EXON_LINE, "Exon line color", "Exon outline color", BLX_BLUE, BLX_BLUE, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CDS_FILL, "CDS fill color", "Coding section fill color", BLX_LIGHT_GREEN, BLX_LIGHT_GREEN, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CDS_LINE, "CDS line color", "Coding section outline color", BLX_DARK_GREEN, BLX_DARK_GREEN, BLX_VERY_DARK_GREEN, BLX_VERY_DARK_GREEN);
  createBlxColor(dc->defaultColors, DOTCOLOR_UTR_FILL, "UTR fill color", "Untranslated region fill color", BLX_LIGHT_RED, BLX_LIGHT_RED, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_UTR_LINE, "UTR line color", "Untranslated region outline color", BLX_DARK_RED, BLX_DARK_RED, BLX_VERY_DARK_RED, BLX_VERY_DARK_RED);

  /* dot plot */
  createBlxColor(dc->defaultColors, DOTCOLOR_CROSSHAIR, "Crosshair", "Color of the crosshair on the dot plot", BLX_BLUE, BLX_BLUE, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_GRID, "Grid", "Line color of the grid on the dot plot", BLX_LIGHT_RED, BLX_LIGHT_RED, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_BORDER, "Grid", "Highlight color for the border where the alignment cannot be calculated", "#ffeeee", "#bbbbbb", NULL, NULL);

  /* greyramp */
  createBlxColor(dc->defaultColors, DOTCOLOR_THRESHOLD_MARKER, "Greyramp threshold marker color", "Outline color of the threshold marker on the greyramp tool", BLX_RED, BLX_RED, BLX_GREEN, BLX_GREEN);
  createBlxColor(dc->defaultColors, DOTCOLOR_MARKER_LINE, "Greyramp marker outline color", "Outline color of the triangle markers on the greyramp tool", BLX_BLACK, BLX_BLACK, BLX_GREEN, BLX_GREEN);
  createBlxColor(dc->defaultColors, DOTCOLOR_MARKER_FILL, "Greyramp marker fill color", "Fill color of the triangle markers on the greyramp tool", BLX_WHITE, BLX_WHITE, NULL, NULL);

  /* misc */
  createBlxColor(dc->defaultColors, DOTCOLOR_BREAKLINE, "Breakline color", "Color of the separator lines between sequences, if there were multiple sequences in the input file", breaklineColor, breaklineColor, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_CANONICAL, "Canonical", "Canonical splice sites", BLX_GREEN, BLX_GREEN, NULL, NULL);
  createBlxColor(dc->defaultColors, DOTCOLOR_NON_CANONICAL, "Non-canonical", "Non-canonical splice sites", BLX_RED, BLX_RED, NULL, NULL);
}


/* Get the initial zoom factor, calculating it if the passed-in zoom is 0 */
static gdouble getInitZoomFactor(DotterContext *dc, const gdouble zoomFacIn, const int qLen, const int sLen)
{
  DEBUG_ENTER("getInitZoomFactor");

  gdouble result = zoomFacIn;
  
  if (result <= 0)
    {
      if (dc->memoryLimit)
        {
          result = (int)sqrt((qLen / dc->numFrames / 1e6 * sLen - 1e-6) / dc->memoryLimit) + 1;
        }
      else
        {
          g_error("Cannot calculate zoom; division by 0 (memory limit is 0).\n");
        }
    }
  
  DEBUG_EXIT("getInitZoomFactor returning %f", result);
  return result;
}


static DotterContext* createDotterContext(DotterOptions *options,
                                          BlxBlastMode blastMode, 
                                          const gboolean showWindow,
                                          const BlxStrand refSeqStrand,
                                          const BlxStrand matchSeqStrand,
                                          MSP *mspList,
                                          GList *seqList,
                                          int matrix[24][24],
                                          char *matrixName)
{
  DEBUG_ENTER("createDotterContext");

  DotterContext *result = (DotterContext*)g_malloc(sizeof *result);
  
  result->blastMode = blastMode;
  result->displaySeqType = (blastMode == BLXMODE_BLASTN) ? BLXSEQ_DNA : BLXSEQ_PEPTIDE;
  result->numFrames = (blastMode == BLXMODE_BLASTX) ? 3 : 1;
  result->geneticCode = stdcode1;
  mtxcpy(result->matrix, matrix);
  result->matrixName = matrixName;
  result->mspList = mspList;
  result->seqList = seqList;
  result->windowList = NULL;
  result->abbrevTitle = options->abbrevTitle;
  
  result->watsonOnly = options->watsonOnly;
  result->crickOnly = options->crickOnly;
  
  /* Set the fixed-width font (not applicable in batch mode) */
  if (showWindow)
    {
      GtkWidget *tmp = gtk_window_new(GTK_WINDOW_TOPLEVEL);
      const char *fontFamily = findFixedWidthFont(tmp);
      result->fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(result->fontDesc, pango_font_description_get_size(tmp->style->font_desc));
      getFontCharSize(tmp, result->fontDesc, &result->charWidth, &result->charHeight);
      gtk_widget_destroy(tmp);
    }
  else
    {
      result->fontDesc = NULL;
    }
  
  result->refSeqName = g_strdup(options->qname);
  result->refSeq = options->qseq; /* take ownership of passed-in seq */
  result->refSeqRev = NULL;
  result->refSeqType = (blastMode == BLXMODE_BLASTP ? BLXSEQ_PEPTIDE : BLXSEQ_DNA);
  
  /* for dna ref sequences, reverse-complement the ref seq */
  if (result->refSeqType == BLXSEQ_DNA && result->refSeq)
    {
      result->refSeqRev = (char*)g_malloc(strlen(result->refSeq) + 1);
      revComplement(result->refSeqRev, result->refSeq);
    }
  else if (result->refSeq)
    {
      /* Just reverse it */
      result->refSeqRev = g_strdup(result->refSeq);
      g_strreverse(result->refSeqRev);
    }
  
  result->refSeqStrand = refSeqStrand;

  result->matchSeqName = g_strdup(options->sname);
  result->matchSeq = options->sseq; /* take ownership of passed-in seq */
  result->matchSeqRev = NULL;
  result->matchSeqType = (blastMode == BLXMODE_BLASTN ? BLXSEQ_DNA : BLXSEQ_PEPTIDE);
  result->matchSeqStrand = matchSeqStrand;
  
  result->refSeqFullRange.min = options->qoffset + 1;
  result->refSeqFullRange.max = options->qoffset + strlen(options->qseq);
  result->matchSeqFullRange.min = options->soffset + 1;
  result->matchSeqFullRange.max = options->soffset + strlen(options->sseq);
  
  result->hozScaleRev = options->hozScaleRev;
  result->vertScaleRev = options->vertScaleRev;
  result->negateCoords = options->negateCoords;
  
  result->displayMirror = options->mirrorImage;
  
  result->memoryLimit = options->memoryLimit;
  
  result->defaultColors = NULL;
  
  
  /* Reverse/comp match seq if applicable */
  if (result->matchSeqType == BLXSEQ_DNA && result->matchSeqStrand == BLXSTRAND_REVERSE && result->matchSeq)
    {
      result->matchSeqRev = (char*)g_malloc(strlen(result->matchSeq) + 1);
      revComplement(result->matchSeqRev, result->matchSeq);
    }
  else if (result->matchSeqStrand == BLXSTRAND_REVERSE && result->matchSeq)
    {
      /* Peptide sequence. Just reverse */
      result->matchSeqRev = g_strdup(result->matchSeq);
      g_strreverse(result->matchSeqRev);
    }
  
  int i = 0;
  for ( ; i < result->numFrames; ++i)
    {
      result->peptideSeqs[i] = NULL;
    }
  
  if (result->blastMode == BLXMODE_BLASTX) 
    {
      /* Create the 3 frame translations (for the strand we're interested in only). */
      const gboolean rev = (result->hozScaleRev);
      char *refSeqToUse = (rev ? result->refSeqRev : result->refSeq);
      
      int i = 0;
      for (i = 0; i < result->numFrames; i++)
        {
          /* Get the start coord at this index and calculate which reading frame it really is
           * (because the first coord in the sequence might not be base 1 in frame 1). */
          const int startCoord = rev ? result->refSeqFullRange.max - i : result->refSeqFullRange.min + i;
        
          int frame = UNSET_INT;
          convertToDisplayIdx(startCoord, TRUE, result, 1, &frame);

          result->peptideSeqs[frame - 1] = blxTranslate(refSeqToUse + i, result->geneticCode);

          DEBUG_OUT("Frame %d starts at coord %d for hoz seq strand = %d.\n", frame, startCoord, result->refSeqStrand);
        }
      
      /* Check all of the frames got set */
      for (i = 0; i < result->numFrames; ++i)
        {
          if (result->peptideSeqs[i] == NULL)
            {
              g_error("Error calculating translated sequence for reading frame %d.\n", i + 1);
            }
        }
    }
  
  if (showWindow)
    { 
      createDotterColors(result, options);
    }
  
  /* Calculate the height and width of the horizontal and vertical scales */
  const int leftBorderChars = max(numDigitsInInt(result->matchSeqFullRange.min), numDigitsInInt(result->matchSeqFullRange.max)) + 1;
  result->scaleWidth = DEFAULT_MAJOR_TICK_HEIGHT * 2 + (roundNearest)((gdouble)leftBorderChars * result->charWidth) + SCALE_LINE_WIDTH;
  result->scaleHeight = DEFAULT_MAJOR_TICK_HEIGHT + roundNearest(result->charHeight) + SCALE_LINE_WIDTH;
  
  result->msgData = &options->msgData;
  
  DEBUG_EXIT("createDotterContext returning");
  return result;
}


static void destroyDotterContext(DotterContext **dc)
{
  DEBUG_ENTER("destroyDotterContext");

  if ((*dc))
    {
    if ((*dc)->blastMode == BLXMODE_BLASTX)
      {
        int i = 0;
        for ( ; i < (*dc)->numFrames; i++)
          {
            if ((*dc)->peptideSeqs && (*dc)->peptideSeqs[i])
              {
                g_free((*dc)->peptideSeqs[i]);
                (*dc)->peptideSeqs[i] = NULL;
              }
          }
      }
    
    /* destroy the msps, sequence structs and colors */
    destroyMspList(&(*dc)->mspList);
    destroyBlxSequenceList(&(*dc)->seqList);
    destroyDotterColors((*dc));

    /* Free stuff g_malloc'ed in calling routine (usually dotterMain) */
    g_free((*dc)->refSeq);
    (*dc)->refSeq = NULL;
    
    g_free((*dc)->refSeqName);
    (*dc)->refSeqName = NULL;
    
    g_free((*dc)->matchSeq);
    (*dc)->matchSeq = NULL;
    
    g_free((*dc)->matchSeqName);
    (*dc)->matchSeqName = NULL;

      if ((*dc)->matrixName)
        {
          g_free((*dc)->matrixName);
          (*dc)->matrixName = NULL;
        }
    }  

  /* free the context struct itself */
  g_free(*dc);
  *dc = NULL;
  
  DEBUG_EXIT("destroyDotterContext returning ");
}

static void destroyDotterWindowContext(DotterWindowContext **dwc)
{
  DEBUG_ENTER("destroyDotterWindowContext");

  if ((*dwc)->printSettings)
    g_object_unref((*dwc)->printSettings);
  
  g_free(*dwc);
  *dwc = NULL;
  
  DEBUG_EXIT("destroyDotterWindowContext returning ");
}

static DotterProperties* dotterGetProperties(GtkWidget *dotterWindow)
{
  return dotterWindow ? (DotterProperties*)(g_object_get_data(G_OBJECT(dotterWindow), "DotterProperties")) : NULL;
}

static void onDestroyDotterWindow(GtkWidget *dotterWindow)
{
  DEBUG_ENTER("onDestroyDotterWindow");

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties)
    {
      if (properties->dotterWinCtx)
        {
          DotterContext *dc = properties->dotterWinCtx->dotterCtx;

          g_object_unref(properties->dotterWinCtx->uiManager);
          properties->dotterWinCtx->uiManager = NULL;

          /* free the context for this window */
          destroyDotterWindowContext(&properties->dotterWinCtx);

          /* if it's the last window, then we also need to destroy the main context
           * and quit the program */
          if (dc && dc->windowList)
            {
              dc->windowList = g_slist_remove(dc->windowList, dotterWindow);

              if (g_slist_length(dc->windowList) < 1)
                {
                  g_slist_free(dc->windowList);
                  dc->windowList = NULL;
                  destroyDotterContext(&dc);
                  gtk_main_quit();
                }
            }
        }

      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(dotterWindow), "DotterProperties", NULL);
    }
  
  DEBUG_EXIT("onDestroyDotterWindow returning ");
}


/* Close all windows associated with the given dotter context */
static void dotterContextCloseAllWindows(DotterContext *dc)
{
  if (!dc || !dc->windowList)
    return;

  /* Can't loop through the list because pointers are removed from it when we destroy the
   * windows. Also, when the last window is destroyed, dc will be destroyed. So: loop while we
   * still have more than one entry in the list. We know we must quit after the last entry. */
  int len = g_slist_length(dc->windowList);
  
  while (len > 0)
    {
      GtkWidget *window = GTK_WIDGET(dc->windowList->data);
      gtk_widget_destroy(window);
      
      if (len == 1) /* just processed the last one so break */
        break;
      else
        len = g_slist_length(dc->windowList);
    }
}

/* "Close" the given window. Only really closes it if it's the main window;
 * for other windows, return false to indicate we haven't handled it. */
static gboolean closeWindow(GtkWidget *widget)
{
  gboolean handled = FALSE;
  const char *name = gtk_widget_get_name(widget);
  
  if (name && strcmp(name, MAIN_WINDOW_NAME) == 0)
    {
      handled = TRUE;
      gtk_widget_destroy(widget);
    }

  return handled;
}


/* Convert a dna index to display (dna or peptide) index, if applicable. If horizontal is true
 * we have the horizontal (reference) sequence, otherwise the vertical (match) sequence. */
int convertToDisplayIdx(const int dnaIdx, const gboolean horizontal, DotterContext *dc, const int frameIn, int *baseNum)
{
  //DEBUG_ENTER("convertToDisplayIdx(dnaIdx=%d, hoz=%d, frameIn=%d)", dnaIdx, horizontal, frameIn);

  int result = dnaIdx;
  
  if (baseNum)
    *baseNum = 1;
  
  /* Match seq coords are always in display coords already, so only do anything if this is the
   * ref seq. Also, we only need to convert if it's peptide-nucleotide match. */
  if (horizontal && dc->blastMode == BLXMODE_BLASTX)
    {
      /* If the strand is reversed, the frame will be inverted; un-invert it first */
      const gboolean rev = (horizontal && dc->hozScaleRev) || (!horizontal && dc->vertScaleRev);
      int frame = frameIn;
      
      double fraction = 0.0;
      
      if (rev)
        fraction = (double)(dnaIdx + frame - 1) / (double)dc->numFrames;
      else
        fraction = (double)(dnaIdx - frame + 1) / (double)dc->numFrames;
      
      DEBUG_OUT("frame=%d, fraction=%f\n", frame, fraction);

      /* Round to the higher value so that 0.3, 0.6 and 1.0 all round
       * to the same value. If values are negative this still works; 0, -0.3
       * and -0.6 will all round to 0.
       *
       * To illustrate (the top section shows nucleotide coords, reading 
       * from top-to-bottom then left-to-right): 
       * 
       *              Display forwards        Display reversed
       * Frame1:         -5 -2  1  4             6  3  0 -3
       * Frame2:         -4 -1  2  5             5  2 -1 -4
       * Frame3:         -3  0  3  6             4  1 -2 -5
       * 
       * Peptide coord:  -1  0  1  2             2  1  0 -1
       *
       * */
      result = ceil(fraction);
    
      /* We want base 1 in the requested reading frame. */
      if (baseNum)
        {
          *baseNum = (dnaIdx - frame + 1) % dc->numFrames;
        
          /* If we have a negative base number, adding 3 shifts it to the
           * correct base number. (This also fixes the fact that mod3 gives 0
           * for base 3.) */
          if (*baseNum < 1)
            *baseNum += dc->numFrames;

          /* invert base number order when the sequence is reversed, i.e. if the mod3 of
           * the coords gives base numbers 123123123 etc then change these to 321321321 etc */
          if (rev)
            {
              *baseNum = dc->numFrames - *baseNum + 1;
            }
        }
    }
  
  //DEBUG_EXIT("convertToDisplayIdx returning %d", result);
  return result;
}


/* Convert a display index (dna or peptide) to a dna index (original sequence coord),
 * if applicable. If horizontal is true we have the horizontal (reference) sequence, 
 * otherwise the vertical (match) sequence. */
int convertToDnaIdx(const int displayIdx, 
                    const gboolean horizontal, 
                    DotterContext *dc, 
                    const int frame, 
                    int baseNum)
{
  DEBUG_ENTER("convertToDnaIdx(displayIdx=%d, hoz=%d, frame=%d, base=%d)", displayIdx, horizontal, frame, baseNum);

  int result = displayIdx;
  
  /* Match seq coords are always in display coords already, so only do anything if this is the
   * ref seq. Also, we only need to convert if it's peptide-nucleotide match. */
  if (horizontal && dc->blastMode == BLXMODE_BLASTX)
    {
      if (dc->hozScaleRev)
        result = (displayIdx * dc->numFrames) - dc->numFrames - frame + baseNum + 1;
      else
        result = (displayIdx * dc->numFrames) - dc->numFrames + frame + baseNum - 1;
      
      DEBUG_OUT("dnaIdx=%d\n", result);
    }
  
  DEBUG_EXIT("convertToDnaIdx returning %d", result);
  return result;
}


static DotterWindowContext* createDotterWindowContext(DotterContext *dotterCtx,
                                                      const IntRange* const refSeqRange,
                                                      const IntRange* const matchSeqRange,
                                                      const gdouble zoomFacIn,
						      const gboolean showWindow)
{
  DEBUG_ENTER("createDotterWindowContext");

  DotterWindowContext *result = (DotterWindowContext*)g_malloc(sizeof *result);
  
  result->dotterCtx = dotterCtx;
  
  result->refSeqRange.min = refSeqRange->min;
  result->refSeqRange.max = refSeqRange->max;
  result->matchSeqRange.min = matchSeqRange->min;
  result->matchSeqRange.max = matchSeqRange->max;

  result->refCoord = UNSET_INT;
  result->matchCoord = UNSET_INT;
  
  result->zoomFactor = getInitZoomFactor(dotterCtx, zoomFacIn, getRangeLength(refSeqRange), getRangeLength(matchSeqRange));

  /* See if we're comparing the same portion of sequence against itself */
  result->selfComp = (refSeqRange->min == matchSeqRange->min && 
                      refSeqRange->max == matchSeqRange->max &&
                      stringsEqual(dotterCtx->refSeq, dotterCtx->matchSeq, FALSE));

  result->usePrintColors = FALSE;

  if (showWindow) 
    {
      result->pageSetup = gtk_page_setup_new();
      gtk_page_setup_set_orientation(result->pageSetup, GTK_PAGE_ORIENTATION_LANDSCAPE);
      
      result->printSettings = gtk_print_settings_new();
      gtk_print_settings_set_orientation(result->printSettings, GTK_PAGE_ORIENTATION_LANDSCAPE);
      gtk_print_settings_set_quality(result->printSettings, GTK_PRINT_QUALITY_HIGH);
      gtk_print_settings_set_resolution(result->printSettings, DEFAULT_PRINT_RESOLUTION);
    }
  else 
    {
      result->pageSetup = NULL;
      result->printSettings = NULL;
    }

  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < DOTDIALOG_NUM_DIALOGS; ++dialogId)
    {
      result->dialogList[dialogId] = NULL;
    }
  
  DEBUG_EXIT("createDotterWindowContext returning");
  return result;
}

/* properties specific to a particular dotter window */
static void dotterCreateProperties(GtkWidget *dotterWindow, 
                                   GtkWidget *greyrampTool, 
                                   GtkWidget *greyrampWindow,
                                   GtkWidget *greyrampContainer,
                                   GtkWidget *alignmentTool,
                                   GtkWidget *alignmentWindow,
                                   GtkWidget *alignmentContainer,
                                   GtkWidget *dotplot,
                                   DotterWindowContext *dotterWinCtx,
                                   const char *exportFileName)
{
  DEBUG_ENTER("dotterCreateProperties");

  if (dotterWindow)
    {
      DotterProperties *properties = (DotterProperties*)g_malloc(sizeof *properties);

      properties->greyrampTool = greyrampTool;
      properties->greyrampWindow = greyrampWindow;
      properties->greyrampContainer = greyrampContainer;
      properties->alignmentTool = alignmentTool;
      properties->alignmentWindow = alignmentWindow;
      properties->alignmentContainer = alignmentContainer;
      properties->dotplot = dotplot;
      properties->windowsDocked = DOCK_WINDOWS_DEFAULT;
      properties->dotterWinCtx = dotterWinCtx;
      properties->exportFileName = exportFileName;
      
      g_object_set_data(G_OBJECT(dotterWindow), "DotterProperties", properties);
      g_signal_connect(G_OBJECT(dotterWindow), "destroy", G_CALLBACK(onDestroyDotterWindow), NULL); 
    }
  
  DEBUG_EXIT("dotterCreateProperties returning ");
}


/* Utility to extract the dotter context from the dotter window properties */
static DotterContext* dotterGetContext(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  return properties ? properties->dotterWinCtx->dotterCtx : NULL;
}


/* Perform required updates following a change to the currently-selected coords */
static void updateOnSelectedCoordsChanged(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  /* Make sure the coords are in range */
  boundsLimitValue(&properties->dotterWinCtx->refCoord, &properties->dotterWinCtx->refSeqRange);
  boundsLimitValue(&properties->dotterWinCtx->matchCoord, &properties->dotterWinCtx->matchSeqRange);
  
  /* Update the alignment view and dotplot */
  updateAlignmentRange(properties->alignmentTool, properties->dotterWinCtx);
  alignmentToolRedrawAll(properties->alignmentTool);
  
  /* Need to clear cached drawables for the alignment tool but can just refresh the dotplot */
  widgetClearCachedDrawable(properties->alignmentTool, NULL);
  refreshDotplot(properties->dotplot);
}


/* Set the initial currently-selected ref seq and match seq coords (i.e. the coords 
 * where the crosshair is centred) */
static void setInitSelectedCoords(GtkWidget *dotterWindow, const int refCoord, const int matchCoord)
{
  DEBUG_ENTER("setInitSelectedCoords");

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  if (valueWithinRange(refCoord, &dwc->refSeqRange))
    dwc->refCoord = refCoord;
  else
    dwc->refCoord = getRangeCentre(&dwc->refSeqRange);

  if (valueWithinRange(matchCoord, &dwc->matchSeqRange))
    dwc->matchCoord = matchCoord;
  else
    dwc->matchCoord = getRangeCentre(&dwc->matchSeqRange);
  
  updateOnSelectedCoordsChanged(dotterWindow);
  
  DEBUG_EXIT("setInitSelectedCoords returning ");
}


/***********************************************************
 *                    External routines                    *
 ***********************************************************/



void dotter (const BlxBlastMode blastMode,
             DotterOptions *options,
             const BlxStrand refSeqStrand,
             const BlxStrand matchSeqStrand,
              int   qcenter,
              int   scenter,
              MSP  *mspList,
              GList *seqList,
             int   MSPoff)
{
  DEBUG_ENTER("dotter(mode=%d, qname=%s, qoff=%d, qstrand=%d, sname=%s, soff=%d, sstrand=%d)",
              blastMode, options->qname, options->qoffset, refSeqStrand, options->sname, options->soffset, matchSeqStrand);
  
  const int qlen = strlen(options->qseq);
  const int slen = strlen(options->sseq);
  
  if (qlen < 1) g_error("queryseq is empty\n");
  if (slen < 1) g_error("subjectseq is empty\n");

  int i = 0;
  for (i = 0; i < qlen; i++) options->qseq[i] = toupper(options->qseq[i]);
  for (i = 0; i < slen; i++) options->sseq[i] = toupper(options->sseq[i]);

  if (!options->memoryLimit) 
    {
      options->memoryLimit = 0.5; /* Mb */
    }

  /* Get score matrix */
  char *matrixName = (char*)g_malloc((MAX_MATRIX_NAME_LENGTH + 1) * sizeof(char));
  
  if (options->mtxfile) 
    {
      readmtx(MATRIX, options->mtxfile);
      strncpy(matrixName, options->mtxfile, MAX_MATRIX_NAME_LENGTH);
    }
  else if (blastMode == BLXMODE_BLASTN) 
    {
      DNAmatrix(MATRIX);
      strcpy(matrixName, "DNA+5/-4");
    }
  else
    {
      mtxcpy(MATRIX, BLOSUM62);
      strcpy(matrixName, "BLOSUM62");
    }

  /* If a save/export file was given, that implies we're in batch mode.
   * We don't display the window in batch mode, unless we're exporting to PDF, in
   * which case we need to create the window so that we have something to print  */
  const gboolean batchMode = (options->savefile || options->exportfile);
  const gboolean createWindow = options->exportfile || !batchMode;
  
  if (batchMode) 
    {
      /* Don't do batch processing if output file can't be opened */
      if (options->savefile && !fopen (options->savefile, "wb"))
        g_error("Failed to open %s\n", options->savefile);

      if (options->exportfile && !fopen (options->exportfile, "wb"))
        g_error("Failed to open %s\n", options->exportfile);
    }
  
  /* Create the main dotter context (shared properties for all dotter windows in this process) */
  DotterContext *dotterCtx = createDotterContext(options, blastMode, createWindow, refSeqStrand, matchSeqStrand, mspList, seqList, MATRIX, matrixName);

  /* Create a context specific to the initial dotter window */
  DotterWindowContext *dotterWinCtx = createDotterWindowContext(dotterCtx, &dotterCtx->refSeqFullRange, &dotterCtx->matchSeqFullRange, options->dotterZoom, createWindow);

  /* Create the widgets */
  createDotterInstance(dotterCtx, 
                       dotterWinCtx,
                       options->loadfile,
                       options->savefile,
                       options->exportfile,
                       options->hspsOnly,
                       options->breaklinesOn,
                       options->winsize,
                       options->pixelFacset,
                       options->dotterZoom,
                       qcenter,
                       scenter,
                       options->swapGreyramp,
                       options->windowColor);

  DEBUG_EXIT("dotter returning");
  return ;
}


/* Create all the widgets for a dotter instance. Uses the existing dotter context. Multiple
 * instances (i.e. multiple dotter windows) can exist that share the same main context but display
 * a different range of coords etc,. This creates the widgets and shows them. */
static GtkWidget* createDotterInstance(DotterContext *dotterCtx,
                                       DotterWindowContext *dotterWinCtx,
                                       const char *loadFileName,
                                       const char *saveFileName,
                                       const char *exportFileName,
                                       const gboolean hspsOn,
                                       const gboolean breaklinesOn,
                                       const char* winsizeIn,
                                       const int pixelFacIn,
                                       const int zoomFacIn,
                                       const int qcenter,
                                       const int scenter,
                                       const gboolean greyrampSwap,
                                       char *windowColor)
{
  DEBUG_ENTER("createDotterInstance");

  GtkWidget *dotterWindow = NULL;

  GtkWidget *dotplot = NULL;
  GtkWidget *dotplotWidget = createDotplot(dotterWinCtx, 
                                           loadFileName,
                                           saveFileName,
                                           exportFileName,
                                           hspsOn,
                                           breaklinesOn,
                                           winsizeIn,
                                           pixelFacIn,
                                           zoomFacIn,
                                           qcenter,
                                           scenter,
                                           &dotplot);
  
  /* Only create the graphical elements if there is a graphical dotplot widget */
  if (dotplotWidget)
    {
      GtkWidget *greyrampContainer = gtk_frame_new(NULL); /* container for when docked */
      gtk_frame_set_shadow_type(GTK_FRAME(greyrampContainer), GTK_SHADOW_NONE);
      GtkWidget *greyrampWindow = NULL; /* container for when undocked */
      GtkWidget *greyrampTool = createGreyrampTool(dotterWinCtx, 40, 100, greyrampSwap, &greyrampWindow);
      registerGreyrampCallback(greyrampTool, dotplot, dotplotUpdateGreymap);
      blxSetWidgetColor(greyrampWindow, windowColor);

      GtkWidget *alignmentContainer = gtk_frame_new(NULL); /* container for when docked */
      gtk_frame_set_shadow_type(GTK_FRAME(alignmentContainer), GTK_SHADOW_NONE);
      GtkWidget *alignmentWindow = NULL; /* container for when undocked */
      GtkWidget *alignmentTool = createAlignmentTool(dotterWinCtx, &alignmentWindow);
      blxSetWidgetColor(alignmentWindow, windowColor);

      if (DOCK_WINDOWS_DEFAULT)
        {
          /* Dock into containers that will go into the main window */
          gtk_container_add(GTK_CONTAINER(alignmentContainer), alignmentTool);
          gtk_container_add(GTK_CONTAINER(greyrampContainer), greyrampTool);
        }
      else
        {
          /* Add to the separate toplevel windows */
          gtk_container_add(GTK_CONTAINER(alignmentWindow), alignmentTool);
          gtk_container_add(GTK_CONTAINER(greyrampWindow), greyrampTool);
        }
  
      const DotterHspMode hspMode = dotplotGetHspMode(dotplot);
      dotterWindow = createDotterWindow(dotterCtx, dotterWinCtx, hspMode, 
                                        dotplot, dotplotWidget, greyrampContainer, alignmentContainer, 
                                        exportFileName, windowColor);
      
      /* Set the tool windows as transient for the main window and clear them up when the
       * main window is destroyed */
      gtk_window_set_transient_for(GTK_WINDOW(greyrampWindow), GTK_WINDOW(dotterWindow));
      gtk_window_set_transient_for(GTK_WINDOW(alignmentWindow), GTK_WINDOW(dotterWindow));
      gtk_window_set_destroy_with_parent(GTK_WINDOW(greyrampWindow), TRUE);
      gtk_window_set_destroy_with_parent(GTK_WINDOW(alignmentWindow), TRUE);

      /* Set the handlers for the alignment and greyramp tools. Connect them here so we can pass
       * the main window as data. */
      gtk_widget_add_events(alignmentWindow, GDK_KEY_PRESS_MASK);
      gtk_widget_add_events(greyrampWindow, GDK_KEY_PRESS_MASK);
      g_signal_connect(G_OBJECT(alignmentWindow), "key-press-event", G_CALLBACK(onKeyPressDotterCoords), dotterWindow);
      g_signal_connect(G_OBJECT(greyrampWindow), "key-press-event", G_CALLBACK(onKeyPressDotter), dotterWindow);

      /* Keep track of all the windows we create, so that we can destroy the DotterContext when
       * the last one is closed */
      dotterCtx->windowList = g_slist_append(dotterCtx->windowList, dotterWindow);
      
      dotterCreateProperties(dotterWindow, 
                             greyrampTool, greyrampWindow, greyrampContainer,
                             alignmentTool, alignmentWindow, alignmentContainer,
                             dotplot, dotterWinCtx, exportFileName);
      DotterProperties *properties = dotterGetProperties(dotterWindow);
      
      setInitSelectedCoords(dotterWindow, qcenter, scenter);
      
      updateGreyMap(properties->greyrampTool);
    }
  
  DEBUG_EXIT("createDotterInstance returning ");
  return dotterWindow;
}


/* Open another dotter window, internal to the existing process (i.e. using the same sequences
 * etc. but just displaying a different range). */
void callDotterInternal(DotterContext *dc, 
                        const IntRange* const refSeqRange,
                        const IntRange* const matchSeqRange,
                        const gdouble zoomFactor,
                        const gboolean breaklinesOn)
{
  DotterWindowContext *dwc = createDotterWindowContext(dc, refSeqRange, matchSeqRange, zoomFactor, TRUE);
  createDotterInstance(dc, dwc, NULL, NULL, NULL, FALSE, breaklinesOn, NULL, 0, 0, 0, 0, FALSE, NULL);
}



/***********************************************************
 *                    External routines                    *
 ***********************************************************/

/* RMEXP will subtract the expected score to remove background noise
   due to biased composition. Gos's idea - not explored yet.

static void rmExp(void){}
*/


///* Utility to return the Feature Series that the given MSP belongs to */
//static FeatureSeries* mspGetFeatureSeries(const MSP *msp)
//{
//  return msp->fs;
//}


/* Returns true if this MSP is part of a Feature Series and that feature series is currently displayed */
//static gboolean mspShowFs(const MSP *msp)
//{
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//  return (fs && fs->on);
//}


/* Return the y position of the lower edge of a feature-series MSP, given its height.
   Calculate it if it is not set yet. Updates the max y coord to be the lower edge of
   this MSP. Returns 0 if the MSP is not in a Feature Series */
//float mspGetFsBottomEdge(MSP *msp, float *maxy, float height)
//{
//  float result = 0;
//  
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//  
//  if (fs)
//    {
//      result = fs->y;
//
//      if (!result)
//        {
//          *maxy += height;
//          fs->y = *maxy;
//          result = *maxy;
//        }
//    }
//  
//  return result;
//}


/* Return the x position of the rightmost edge of a feature-series MSP, given its height.
   Calculae it if it is not yet set. Updates the max x coord to be the rightmost edge of
   this MSP. Returns 0 if the MSP is not in a FeatureSeries */
//float mspGetFsRightEdge(MSP *msp, float *maxx, float height)
//{
//  float result = 0;
//
//  FeatureSeries *fs = mspGetFeatureSeries(msp);
//
//  if (fs)
//    {
//      result = fs->x;
//      
//      if (!result)
//        {
//          *maxx += height;
//          fs->x = *maxx;
//          result = *maxx;
//        }
//    }
//  
//  return result;
//}


//static float seq2graphX(int pos)
//{
//    float x;
//
//    x = ceil((float)(pos + MSPoffset - qoffset)/zoom/resfac);
//
//    if (reversedScale)
//      x = RightBorder - x;
//    else
//      x += LeftBorder-1;
//
//    return x;
//}
//
//static float seq2graphY(int pos)
//{
//    float y;
//
//    y = ceil((float)(pos - soffset)/zoom);
//    
//    y += TopBorder-1;
//
//    return y;
//}


//static void XdrawSEG(MSP *msp, float offset)
//{
//    /* Horizontal axis */
//
//    float  
//      sx = seq2graphX(mspGetQStart(msp)), 
//      ex = seq2graphX(mspGetQEnd(msp));
//    
//    offset += TopBorder + slen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    if (fsEndLinesOn) {
//      graphLine(sx, TopBorder, sx, TopBorder+slen4-2);
//      graphLine(ex, TopBorder, ex, TopBorder+slen4-2);
//    }
//
//    graphFillRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);
//    graphColor(BLACK);
//    graphRectangle(sx, offset, ex, offset - msp->score/100.0 * fonth);
//
//    graphColor(BLACK);
//    if (fsAnnBottomOn && msp->desc) {
//      graphText(msp->desc, sx, offset);
//    }    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


//static void XdrawSEGxy(MSP *msp, float offset)
//{
//    int i, inNotFilled=0, descShown=0;
//    float  
//      x, y, 
//      xold=0, yold=0;
//
//    offset += TopBorder + slen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    for (i = 0; i < qlen; i++)
//      {
//      const int xyVal = g_array_index(msp->xy, int, i);
//
//      if (xyVal == XY_NOT_FILLED)
//        {
//          inNotFilled = 1;
//        }
//      else
//        {
//          x = seq2graphX(i);
//          y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
//          
//          if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE))
//            {
//              graphLine(xold, yold, x, y);
//              
//              if (fsAnnBottomOn && msp->desc && !descShown)
//                {
//                  int linecolor = graphColor(BLACK);
//                  graphText(msp->desc, xold, offset);
//                  graphColor(linecolor);
//                  descShown = 1;
//                }
//            }
//              
//          xold = x;
//          yold = y;
//          inNotFilled = 0;
//        }
//      }
//    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}
//
//
//static void YdrawSEG(MSP *msp, float offset)
//{
//    /* Vertical axis */
//
//    float  
//      sx = seq2graphY(mspGetQStart(msp)), 
//      ex = seq2graphY(mspGetQEnd(msp));
//    
//    offset += LeftBorder + qlen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    if (fsEndLinesOn) {
//      graphLine(LeftBorder, sx, LeftBorder+qlen4-2, sx);
//      graphLine(LeftBorder, ex, LeftBorder+qlen4-2, ex);
//    }
//
//    graphFillRectangle(offset, sx, offset - msp->score/100.0 * fonth, ex);
//    graphColor(BLACK);
//    graphRectangle    (offset, sx, offset - msp->score/100.0 * fonth, ex);
//
//    graphColor(BLACK);
//    if (fsAnnRightOn && msp->desc) {
//      graphText(msp->desc, offset, sx);
//    }    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}
//
//
//static void YdrawSEGxy(MSP *msp, float offset)
//{
//    int i, inNotFilled=0, descShown=0;
//    float  
//      x, y, 
//      xold=0, yold=0;
//
//    offset += LeftBorder + qlen4 -1;
//    oldcolor = graphColor(msp->fsColor); oldLinew = graphLinewidth(.25);
//
//    for (i = 0; i < qlen; i++)
//      {
//      const int xyVal = g_array_index(msp->xy, int, i);
//      
//      if (xyVal == XY_NOT_FILLED)
//        {
//          inNotFilled = 1;
//        }
//      else
//        {
//          x = seq2graphY(i);
//          y = offset-1 - (float)xyVal / 100 * fsPlotHeight * fonth;
//          
//          if (xold && (x != xold || y != yold) && (!inNotFilled || msp->fsShape == BLXCURVE_INTERPOLATE)) 
//            {
//              graphLine(yold, xold, y, x);
//              
//              if (fsAnnRightOn && msp->desc && !descShown) 
//                {
//                  int linecolor = graphColor(BLACK);
//                  graphText(msp->desc, offset, xold);
//                  graphColor(linecolor);
//                  descShown = 1;
//                }
//            }
//
//          xold = x;
//          yold = y;
//          inNotFilled = 0;
//        }
//      }
//    
//    graphColor(oldcolor); graphLinewidth(oldLinew);
//}


//static int isHorizontalMSP(MSP *msp) {
//
//    if (!msp->qname)
//      g_error("No qname set in MSP - I need this to associate it with one of the sequences");
//
//    if (!strcasecmp(msp->qname, qname) || !strcmp(msp->qname, "@1"))
//        return 1;
//    else
//        return 0;
//}
//
//
//static int isVerticalMSP(MSP *msp) {
//
//    if (!msp->qname) 
//      g_error("No qname set in MSP - I need this to associate it with one of the sequences");
//
//    if (!strcasecmp(msp->qname, sname) || !strcmp(msp->qname, "@2"))
//        return 1;
//    else
//        return 0;
//}


//static void drawAllFeatures(MSP *msp)
//{
//  float sx, ex, 
//    y, 
//    height,           /* box height/width */
//    boxHeight=10,
//    textHeight,
//    oldLinew;
//  int i;
//  float posx=0, posy=0;               /* Next series to be drawn */
//
//  int top, bottom ;
//  float feature_boundary, feature_top, feature_bottom, feature_strand ;
//  float forward_y, reverse_y, depth ;
//  float old_line_width ;
//
//
//  /* Set forward/reverse strand gene drawing positions. */
//  graphGetBounds(&top, &bottom) ;
//  feature_boundary = 3.0 ;                                /* Allows for line thickness etc... */
//  feature_top = TopBorder + slen4 ;
//  feature_bottom = bottom ;
//  feature_strand = feature_top + ((feature_bottom - feature_top + 1) / 2) ;
//
//  forward_y = feature_top + feature_boundary ;
//  reverse_y = feature_strand + feature_boundary ;
//  depth = (feature_strand - feature_boundary) - (feature_top + feature_boundary) - fonth ;
//
//
//  /* Draw a strand separator. */
//  old_line_width = graphLinewidth(2) ;
//  graphLine(LeftBorder - 1, feature_strand, LeftBorder - 1 + qlen4, feature_strand) ;
//  graphLinewidth(old_line_width) ;
//
//
//  /* Now draw the genes. */
//  drawGenes(msp, forward_y, reverse_y, depth) ;
//
//
//  if (fsArr)
//    {
//      /* Loop through each feature-series in the feature-series array and set coords to 0 */
//      for (i = 0; i < gArrayGetLen(fsArr); i++)
//      {
//        FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//        fs->y = fs->x = 0;
//      }
//    }
//
//  textHeight = fonth;
//  boxHeight = fonth;
//
//  oldLinew = graphLinewidth(1);
//
//  BlxStrand strand = BLXSTRAND_NONE;
//  
//  if (selfcomp || !reversedScale) 
//    strand = BLXSTRAND_FORWARD;
//  else 
//    strand = BLXSTRAND_REVERSE;
//
//
//  for (; msp; msp = msp->next)
//    {    
//      height = boxHeight;
//
//      if (mspHasFs(msp))
//      {
//        sx = seq2graphX(mspGetQStart(msp));
//        ex = seq2graphX(mspGetQEnd(msp));
//
//        if (msp->qStrand != strand)
//          y = reverse_y ;
//        else
//          y = forward_y ;
//
//          if (!mspShowFs(msp))
//          continue;
//
//        /* Adjust height to score */
//        if (msp->score > 0)
//          {
//            height = boxHeight * msp->score / 100.0;
//          }
//
//        if (fsBottomOn && (isHorizontalMSP(msp) || (selfcomp && isVerticalMSP(msp))))
//          {
//            /* HORIZONTAL AXIS (X) */
//                  
//            if (msp->type == BLXMSP_XY_PLOT)
//              {
//                XdrawSEGxy(msp, mspGetFsRightEdge(msp, &posx, fonth*(fsPlotHeight+1)));
//              }
//            else if (msp->type == BLXMSP_FS_SEG)
//              {
//                XdrawSEG(msp, mspGetFsRightEdge(msp, &posx, boxHeight+textHeight));
//              }
//          }
//
//        if (fsRightOn &&  (isVerticalMSP(msp) || (selfcomp && isHorizontalMSP(msp))))
//          {
//            /* VERTICAL AXIS (Y) */
//
//            if (msp->type == BLXMSP_XY_PLOT)
//              {
//                YdrawSEGxy(msp, mspGetFsBottomEdge(msp, &posy, fonth*(fsPlotHeight+1)));
//              }
//            else if (msp->type == BLXMSP_FS_SEG)
//              {
//                YdrawSEG(msp, mspGetFsBottomEdge(msp, &posy, boxHeight+textHeight));
//              }
//          }
//
//        graphColor(oldcolor); 
//        graphLinewidth(oldLinew);
//      }
//    }
//
//  return ;
//}


//static void drawGenes(MSP *msp, float forward_y, float reverse_y, float depth)
//{
//  gboolean bump_genes = FALSE ;                                   /* Make this user settable... */
//
//  float height,               /* box height/width */
//    sy, ey, midy, x, y, oldLinew ;
//
//  height = fonth ;
//
//  oldLinew = graphLinewidth(1);
//
//  BlxStrand strand = BLXSTRAND_NONE;
//  if (selfcomp || !reversedScale) 
//    strand = BLXSTRAND_FORWARD;
//  else 
//    strand = BLXSTRAND_REVERSE;
//
//
//  if (bump_genes)
//    {
//      MSP *gene_msp, *tmp ;
//      GList *exon_intron_list = NULL ;
//      GeneStrandStruct strand_genes = {NULL} ;
//
//      /* MSP's are in any old order and only some are exons/introns, here we make a list of
//       * intron/exon MSP's and then sort that list into transcripts for drawing... */
//
//      /* Sort all introns/exons by gene and position within gene. */
//      gene_msp = tmp = msp ;
//      for ( ; tmp ; tmp = tmp->next)
//      {
//        if (tmp->score < 0)
//          exon_intron_list = g_list_append(exon_intron_list, tmp) ;
//      }
//      exon_intron_list = g_list_sort(exon_intron_list, compareMSPs) ;
//
//
//#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
//      g_list_foreach(exon_intron_list, printMSP, NULL) ;
//#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
//
//
//      /* Produce two lists, one for each strand, of genes sorted by position from the sorted intron/exon list. */
//      strand_genes.strand = '+' ;
//      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
//      strand_genes.forward_genes = g_list_sort(strand_genes.forward_genes, compareGenes) ;
//      strand_genes.forward_genes = g_list_first(strand_genes.forward_genes) ;
//
//      strand_genes.strand = '-' ;
//      g_list_foreach(exon_intron_list, getGenePositionsCB, &strand_genes) ;
//      strand_genes.reverse_genes = g_list_sort(strand_genes.reverse_genes, compareGenes) ;
//      strand_genes.reverse_genes = g_list_first(strand_genes.reverse_genes) ;
//
//      /* Set y offsets for forward and reverse stranded genes. */
//      setYoffsets(strand_genes.forward_genes, forward_y, forward_y + depth) ;
//      setYoffsets(strand_genes.reverse_genes, reverse_y, reverse_y + depth) ;
//
//
//#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
//      printf("Forward:\n") ;
//      g_list_foreach(strand_genes.forward_genes, printGene, NULL) ;
//      printf("Reverse:\n") ;
//      g_list_foreach(strand_genes.reverse_genes, printGene, NULL) ;
//      printf("\n") ;
//#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */
//
//      /* Draw the strands. */
//      oldcolor = graphColor(BLUE); 
//      g_list_foreach(strand_genes.forward_genes, drawGenesCB, &strand_genes) ;
//
//      oldcolor = graphColor(MAGENTA); 
//      g_list_foreach(strand_genes.reverse_genes, drawGenesCB, &strand_genes) ;
//
//      graphColor(oldcolor); 
//
//    }
//  else
//    {
//      for (; msp; msp = msp->next)
//      {    
//        if (msp->score < 0)
//          {
//            if (msp->qStrand != strand)
//              y = reverse_y ;
//            else
//              y = forward_y ;
//
//            oldcolor = graphColor(BLUE); 
//            drawMSPGene(msp, y) ;
//
//            if (selfcomp) /* Draw the boxes on vertical axes too */
//              {
//                sy = ceil((float)(mspGetQStart(msp)+MSPoffset - qoffset)/zoom);
//                ey = ceil((float)(mspGetQEnd(msp)+MSPoffset - qoffset)/zoom);
//              
//                sy += TopBorder-1;
//                ey += TopBorder-1;
//              
//                x = LeftBorder + qlen4 + 10;
//                if (msp->qStrand != strand) x += 20;
//              
//                if (msp->score == -1.0) /* EXON */
//                  {
//                    oldcolor = graphColor(BLUE); 
//                    graphRectangle(x, sy, x + height, ey);
//                  }
//                else if (msp->score == -2.0) /* INTRON */
//                  {
//                    oldcolor = graphColor(BLUE); 
//                    midy = 0.5 * (sy + ey) ;
//                    graphLine (x + height/2, sy, x, midy) ;
//                    graphLine (x + height/2, ey, x, midy) ;
//                  }
//              }
//
//            graphColor(oldcolor) ;
//          }
//      }
//    }
//
//  graphLinewidth(oldLinew) ;
//
//  return ;
//}


/* A GCompareFunc() to compare gene MSPs.... */
//gint compareMSPs(gconstpointer a, gconstpointer b)
//{
//  int result = 0 ;
//  MSP *msp_a = (MSP *)a, *msp_b = (MSP *)b ;
//
//  if (!(result = strcmp(msp_a->sname, msp_b->sname)))
//    {
//      if (mspGetSStart(msp_a) < mspGetSStart(msp_b))
//      result = -1 ;
//      else if (mspGetSStart(msp_a) > mspGetSStart(msp_b))
//      result = 1 ;
//      else
//      {
//        /* I actually don't think this should ever happen as it means either there are
//         * duplicates or worse there are overlapping introns/exons within a gene.... */
//        if (mspGetSEnd(msp_a) < mspGetSEnd(msp_b))
//          result = -1 ;
//        else if (mspGetSEnd(msp_a) > mspGetSEnd(msp_b))
//          result = 1 ;
//        else
//          result = 0 ;
//      }
//    }
//
//  return result ;
//}


///* A GFunc() to record start/end of genes. */
//static void getGenePositionsCB(gpointer data, gpointer user_data)
//{
//  MSP *msp = (MSP *)data ;
//  GeneStrand strand_genes = (GeneStrand)user_data ;
//  GList *gene_list ;
//  GeneData gene ;
//
//  /* We need genes added according to strand... */
//  if (msp->qframe[1] == strand_genes->strand)
//    {
//      const char *curr_name = NULL ;
//      int curr_length ;
//
//      if (strand_genes->strand == '+')
//      gene_list = strand_genes->forward_genes ;
//      else
//      gene_list = strand_genes->reverse_genes ;
//
//      if (gene_list)
//      {
//        /* Look at last element, this is the last gene we added. */
//        gene_list = g_list_last(gene_list) ;
//
//        curr_name = (((GeneData)(gene_list->data))->name) ;
//        curr_length = strlen(curr_name) - 1 ;
//      }
//
//      /* If there's no gene or we are starting a new gene then just add to the list,
//       * otherwise record latest msp end position. */
//      if (!gene_list || strncmp(mspGetSName(msp), curr_name, curr_length) != 0)
//      {
//        gene = g_new0(GeneDataStruct, 1) ;
//        gene->name = mspGetSName(msp) ;
//        gene->start = mspGetSStart(msp) ;
//        gene->end = mspGetSEnd(msp) ;
//        gene->strand = msp->qframe[1] ;
//        gene->msp_start = msp ;
//
//        gene_list = g_list_append(gene_list, gene) ;
//      }
//      else
//      {
//        gene = (GeneData)(gene_list->data) ;
//
//        gene->end = mspGetSEnd(msp) ;
//        gene->msp_end = msp ;
//      }
//
//      if (strand_genes->strand == '+')
//      strand_genes->forward_genes = gene_list ;
//      else
//      strand_genes->reverse_genes = gene_list ;
//    }
//
//  return ;
//}
//
//

/* A GCompareFunc() to compare genes for position.... */
//gint compareGenes(gconstpointer a, gconstpointer b)
//{
//  int result = 0 ;
//  GeneData gene_a = (GeneData)a, gene_b = (GeneData)b ;
//
//  if (gene_a->strand == '+' && gene_b->strand == '-')
//    result = -1 ;
//  else if (gene_a->strand == '-' && gene_b->strand == '+')
//    result = 1 ;
//  else
//    {
//      if (gene_a->start < gene_b->start)
//      result = -1 ;
//      else if (gene_a->start > gene_b->start)
//      result = 1 ;
//      else
//      result = 0 ;
//    }
//
//  return result ;
//}



//static void setYoffsets(GList *gene_list, float min_offset, float max_offset)
//{
//  GList *curr_ptr, *next_ptr ;
//  float curr_y = 0.0, bump_incr = 0.5 ;
//
//
//  /* Go through the gene list comparing, reordering and assigning y offsets to the genes
//   * so they can be drawn without overlapping. */
//  curr_ptr = gene_list ;
//  next_ptr = curr_ptr->next ;
//  curr_y = min_offset - bump_incr ;
//
//  /* Go through all genes. */
//  while (curr_ptr)
//    {
//      GList *list_ptr = g_list_first(gene_list) ;
//
//      /* For each gene look at all the preceding genes and to decide its y offset so it is
//       * bumped correctly. */
//      while (list_ptr)
//      {
//        GeneData list_gene = (GeneData)list_ptr->data ;
//        GeneData curr_gene = (GeneData)curr_ptr->data ;
//
//        if (list_ptr == curr_ptr)
//          {
//            /* This gene overlaps all previous ones and needs a new offset. */
//            curr_y += bump_incr ;
//            if (curr_y > max_offset)
//              curr_y = max_offset ;
//
//            curr_gene->y_pos = curr_y ;
//            break ;
//          }
//        else if (curr_gene->start > list_gene->end)
//          {
//            /* This gene coes not overlap the list gene so give is the same offset. */
//            curr_gene->y_pos = list_gene->y_pos ;
//
//            if (list_ptr->next != curr_ptr)
//              {
//                list_ptr = g_list_remove(list_ptr, curr_gene) ;
//                list_ptr = g_list_insert_before(gene_list, list_ptr->next, curr_gene) ;
//              }
//
//            break ;
//          }
//        else
//          {
//            /* This gene overlaps the list gene so move on to the next one. */
//            list_ptr = g_list_next(list_ptr) ;
//          }
//      }
//
//      /* update curr/next until we get to the end of the list... */
//      if ((curr_ptr = next_ptr))
//      next_ptr = curr_ptr->next ;
//    }
//
//  return ;
//}
//

///* A GFunc() to record start/end of genes. */
//static void drawGenesCB(gpointer data, gpointer user_data_unused)
//{
//  GeneData gene = (GeneData)data ;
//  MSP *msp ;
//
//  msp = gene->msp_start ;
//  do
//    {
//      drawMSPGene(msp, gene->y_pos) ;
//
//      msp = msp->next ;
//    } while (msp && msp != gene->msp_end) ;
//
//  return ;
//}
//
//
//static void drawMSPGene(MSP *msp, float y_offset)
//{
//  float height,               /* box height/width */
//    sx, ex, midx ;
//
//  height = fonth ;
//
//  sx = seq2graphX(mspGetQStart(msp));
//  ex = seq2graphX(mspGetQEnd(msp));
//
//  if (msp->score == -1.0) /* EXON */
//    {
//      graphRectangle(sx, y_offset, ex, y_offset + height);
//    }
//  else if (msp->score == -2.0) /* INTRON */
//    {
//
//      midx = 0.5 * (sx + ex) ;
//      graphLine (sx, y_offset + height/2, midx, y_offset) ;
//      graphLine (ex, y_offset + height/2, midx, y_offset) ;
//    }
//
//  return ;
//}
//

//static void printMSP(gpointer data, gpointer user_data)
//{
//  MSP *msp = (MSP *)data ;
//  
//  printf("%s:\t%d,%d\t->\t%d,%d\n",
//       msp->sname, msp->sstart, msp->send, msp->qstart, msp->qend) ;
//  
//
//  return ;
//}
//

//static void printGene(gpointer data, gpointer user_data)
//{
//  GeneData gene = (GeneData)data ;
//  
//  printf("%s: '%c' %d -> %d   is at position: %f\n",
//       gene->name, gene->strand, gene->start, gene->end, gene->y_pos) ;
//  
//
//  return ;
//}


//static void fsSelAll(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      fs->on = 1;
//      graphBoxDraw(fsBoxStart+i, WHITE, BLACK);
//    }
//}
//
//static void fsSelNone(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      fs->on = 0;
//      graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//    }
//}
//
//static void fsSelNoCurves(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      
//      if (fs->xy) 
//      {
//        fs->on = 0;
//        graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//      }
//    }
//}
//
//static void fsSelNoSegments(void)
//{
//  int i;
//  graphActivate(fsGraph);
//  
//  for (i = 0; i < gArrayGetLen(fsArr); i++) 
//    {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      
//      if (!fs->xy) 
//      {
//        fs->on = 0;
//        graphBoxDraw(fsBoxStart+i, BLACK, WHITE);
//      }
//    }
//}
//
//static void fsToggleBottom(void)
//{
//    fsBottomOn = !fsBottomOn;
//    selectFeatures();
//}
//static void fsToggleRight(void)
//{
//    fsRightOn = !fsRightOn;
//    selectFeatures();
//}
//static void fsToggleAnnBottom(void)
//{
//    fsAnnBottomOn = !fsAnnBottomOn;
//    selectFeatures();
//}
//static void fsToggleAnnRight(void)
//{
//    fsAnnRightOn = !fsAnnRightOn;
//    selectFeatures();
//}
//static void fsToggleEndLines(void)
//{
//    fsEndLinesOn = !fsEndLinesOn;
//    selectFeatures();
//}
//static void fsSel(int box, double x_unused, double y_unused, int modifier_unused)
//{
//    if (box-fsBoxStart < 0 || box-fsBoxStart > gArrayGetLen(fsArr))
//      return;
//    
//    FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, box - fsBoxStart);
//    int *on = &fs->on;
//
//    graphActivate(fsGraph);
//
//
//    if (*on) {
//      *on = 0;
//      graphBoxDraw(box, BLACK, WHITE);
//    }
//    else {
//      *on = 1;
//      graphBoxDraw(box, WHITE, BLACK);
//    }
//
//}
//static void setfsPlotHeight(char *cp)
//{
//    fsPlotHeight = atof(cp);
//}
//    
//
//void selectFeatures(void)
//{
//  int i, box;
//  float y=1.0;
//
//  if (!graphActivate(fsGraph))
//    {
//      fsGraph = graphCreate (TEXT_SCROLL, "Feature Series Selection Tool", 0, 0, 0.4, 0.6);
//    }
//  graphPop();
//  graphRegister(PICK, fsSel);
//
//  if (1 /* dotterWindow */) {
//
//    box = graphButton("Show bottom series", fsToggleBottom, 1, y);
//    if (!fsBottomOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show right series", fsToggleRight, 1, y);
//    if (!fsRightOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show bottom annotations", fsToggleAnnBottom, 1, y);
//    if (!fsAnnBottomOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Show right annotations", fsToggleAnnRight, 1, y);
//    if (!fsAnnRightOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 1.5;
//
//    box = graphButton("Draw lines at segment ends", fsToggleEndLines, 1, y);
//    if (!fsEndLinesOn) graphBoxDraw(box, BLACK, WHITE);
//    else graphBoxDraw(box, WHITE, BLACK);
//    y += 2;
//  }
//
//  sprintf(fsPlotHeighttx, "%.1f", fsPlotHeight);
//  graphText("Height of curves (XY series):", 1, y);
//  graphTextScrollEntry (fsPlotHeighttx, 6, 4, 31, y, setfsPlotHeight);
//  y += 2;
//
//  graphLinewidth(0.1);
//  graphLine(0, y, 2000, y);
//  y += 0.5;
//
//  graphText("Pick to select/unselect series", 1, y);
//  y += 1.5;
//  graphButton("All", fsSelAll, 1, y);
//  graphButton("No curves", fsSelNoCurves, 13, y);
//  graphButton("No segments", fsSelNoSegments, 24, y);
//  fsBoxStart = 1+graphButton("None", fsSelNone, 6, y);
//  y += 2;
//
//  graphTextBounds(50, gArrayGetLen(fsArr) * 1.5 + y + 5);
//
//  for (i = 0; i < gArrayGetLen(fsArr); i++)
//    {
//      float  margin = 0.1;
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//
//      box = graphBoxStart();
//      graphText(fs->name, 1, y);      
//      graphRectangle(1 - margin, y - margin, 1 + margin + strlen(fs->name), y + 1 + margin);
//      graphBoxEnd();
//
//      if (!fs->on) 
//      {
//        graphBoxDraw(box, BLACK, WHITE);
//      }
//      else 
//      {
//        graphBoxDraw(box, WHITE, BLACK);
//      }
//      
//      y += 1.5;
//    }
//
//  graphRedraw();
//
//  return ;
//}
//
//
//float fsTotalHeight(MSP *msplist)
//{
//    int i;
//    float maxy = 0;
//      
//    if (!fsArr || !gArrayGetLen(fsArr))
//      {
//      return 0.0;
//      }
//
//    for (i = 0; i < gArrayGetLen(fsArr); i++) 
//      {
//      FeatureSeries *fs = &g_array_index(fsArr, FeatureSeries, i);
//      fs->y = fs->x = 0;
//      }
//      
//    for (msp = msplist; msp; msp = msp->next) 
//      {
//        if (mspShowFs(msp))
//        {
//          if (msp->type == BLXMSP_XY_PLOT) 
//            {
//              mspGetFsBottomEdge(msp, &maxy, fsPlotHeight+1);
//            }
//          else if (msp->type == BLXMSP_FS_SEG) 
//            {
//              mspGetFsBottomEdge(msp, &maxy, 1+1);
//            }
//        }
//      }
//    
//    return maxy + 2;
//}


/* Callbacks to be called when the dotter parameters have changed. */
static gboolean onZoomFactorChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = TRUE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  gdouble newValue = g_strtod(text, NULL);
  
  if (newValue <= 0)
    {
      g_critical("Zoom factor must be greater than zero.\n");
      result = FALSE;
    }
  else if (newValue != properties->dotterWinCtx->zoomFactor)
    {
      properties->dotterWinCtx->zoomFactor = newValue;
      redrawAll(dotterWindow, NULL);
    }
  
  return result;
}

static gboolean onQStartChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  /* If display coords are negated, we must un-negate it before we use it */
  if (negateDisplayCoord(dwc->dotterCtx, TRUE))
    newValue *= -1;
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->refSeqFullRange))
    g_warning("Limiting reference sequence start to range %d -> %d.\n", dwc->dotterCtx->refSeqFullRange.min, dwc->dotterCtx->refSeqFullRange.max);
  
  boundsLimitValue(&newValue, &dwc->dotterCtx->refSeqFullRange);
  
  gboolean changed = setStartCoord(dotterWindow, dwc, TRUE, newValue);
  
  /* If it's a self comparison, also update the vertical range. */
  if (dwc->selfComp)
    changed = setStartCoord(dotterWindow, dwc, FALSE, newValue) || changed;

  /* Check the crosshair is still in range and if not clip it */
  updateOnSelectedCoordsChanged(dotterWindow);

  if (changed)
    redrawAll(dotterWindow, NULL);
  
  return TRUE;
}

static gboolean onQEndChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);

  /* If display coords are negated, we must un-negate it before we use it */
  if (negateDisplayCoord(dwc->dotterCtx, TRUE))
    newValue *= -1;
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->refSeqFullRange))
    g_warning("Limiting reference sequence end to range %d -> %d.\n", dwc->dotterCtx->refSeqFullRange.min, dwc->dotterCtx->refSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->refSeqFullRange);

  gboolean changed = setEndCoord(dotterWindow, dwc, TRUE, newValue);
  
  /* If it's a self comparison, also update the vertical range. */
  if (dwc->selfComp)
    changed = setEndCoord(dotterWindow, dwc, FALSE, newValue) || changed;
  
  /* Check the crosshair is still in range and if not clip it */
  updateOnSelectedCoordsChanged(dotterWindow);

  if (changed)
    redrawAll(dotterWindow, NULL);

  return TRUE;
}

static gboolean onSStartChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  /* If display coords are negated, we must un-negate it before we use it */
  if (negateDisplayCoord(dwc->dotterCtx, FALSE))
    newValue *= -1;
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->matchSeqFullRange))
    g_warning("Limiting vertical sequence start to range %d -> %d.\n", dwc->dotterCtx->matchSeqFullRange.min, dwc->dotterCtx->matchSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->matchSeqFullRange);
  
  gboolean changed = setStartCoord(dotterWindow, dwc, FALSE, newValue);

  /* Check the crosshair is still in range and if not clip it */
  updateOnSelectedCoordsChanged(dotterWindow);

  if (changed)
    redrawAll(dotterWindow, NULL);

  return TRUE;
}

static gboolean onSEndChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);
  
  /* If display coords are negated, we must un-negate it before we use it */
  if (negateDisplayCoord(dwc->dotterCtx, FALSE))
    newValue *= -1;
  
  if (!valueWithinRange(newValue, &dwc->dotterCtx->matchSeqFullRange))
    g_warning("Limiting vertical sequence end to range %d -> %d.\n", dwc->dotterCtx->matchSeqFullRange.min, dwc->dotterCtx->matchSeqFullRange.max);

  boundsLimitValue(&newValue, &dwc->dotterCtx->matchSeqFullRange);
  
  gboolean changed = setEndCoord(dotterWindow, dwc, FALSE, newValue);
  
  /* Check the crosshair is still in range and if not clip it */
  updateOnSelectedCoordsChanged(dotterWindow);

  if (changed)
    redrawAll(dotterWindow, NULL);
  
  return TRUE;
}

static gboolean onSlidingWinSizeChanged(GtkWidget *widget, const gint responseId, gpointer data)
{
  gboolean result = FALSE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const char *text = gtk_entry_get_text(GTK_ENTRY(widget));
  int newValue = convertStringToInt(text);

  GError *error = NULL;
  gboolean changed = dotplotSetSlidingWinSize(properties->dotplot, newValue, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
    }
  else
    {
      result = TRUE;
    }
  
  if (changed)
    redrawAll(dotterWindow, NULL);
  
  return result;
}


/* Callback called when the user has changed the 'splice sites on' option */
static gboolean onSetSpliceSitesOn(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const gboolean spliceSitesOn = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  alignmentToolSetSpliceSitesOn(properties->alignmentTool, spliceSitesOn);
  
  return TRUE;
}

/* Callback called when the user has changed the 'breaklines on' option */
static gboolean onSetBreaklinesOn(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const gboolean breaklinesOn = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  dotplotSetBreaklinesOn(properties->dotplot, breaklinesOn);
  
  return TRUE;
}

/* Callback called when the user has changed the 'horizontal labels on' option */
static gboolean onSetHozLabelsOn(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const gboolean labelsOn = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  dotplotSetHozLabelsOn(properties->dotplot, labelsOn);
  
  return TRUE;
}

/* Callback called when the user has changed the 'vertical labels on' option */
static gboolean onSetVertLabelsOn(GtkWidget *button, const gint responseId, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const gboolean labelsOn = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  dotplotSetVertLabelsOn(properties->dotplot, labelsOn);
  
  return TRUE;
}


/* Callback when we receive a response for the settings dialog */
static void onResponseSettingsDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Destroy if successful */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      break;
      
    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
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


/* Create the paraemter control widgets for the settings dialog */
static void settingsDialogParamControls(GtkWidget *dialog, GtkWidget *dotterWindow, const int border)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  DotterContext *dc = dwc->dotterCtx;
  
  /* Put everything in a frame */
  GtkWidget *frame = gtk_frame_new("Parameters");
  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), frame);
  gtk_container_set_border_width(GTK_CONTAINER(frame), border); 
  
  /* Create a table to lay out the widgets */
  const int numRows = 4;
  const int numCols = 3;
  const int xpad = 2;
  const int ypad = 2;
  
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(table));
  
  /* Get the start and end values of each range, and negate them for display if necessary */
  const int qStart = getDisplayCoord(getStartCoord(dwc, TRUE), dc, TRUE);
  const int qEnd = getDisplayCoord(getEndCoord(dwc, TRUE), dc, TRUE);
  const int sStart = getDisplayCoord(getStartCoord(dwc, FALSE), dc, FALSE);
  const int sEnd = getDisplayCoord(getEndCoord(dwc, FALSE), dc, FALSE);
  
  GtkWidget *zoomEntry = createTextEntryFromDouble(dotterWindow, table, 1, 2, xpad, ypad, "_Zoom: ", dwc->zoomFactor, onZoomFactorChanged);
  gtk_widget_set_tooltip_text(zoomEntry, "Zoom out by this factor, e.g. a zoom factor of 3 will shrink the window to 1/3 of its full size");

  GtkWidget *windowEntry = NULL;
  
  /* Create the boxes for the sequence ranges. If it's a self comparison, we only really have one range. */
  if (dwc->selfComp)
    {
      createTextEntryFromInt(dotterWindow, table, 2, 2, xpad, ypad, "Range: ", qStart, onQStartChanged);
      createTextEntryFromInt(dotterWindow, table, 2, 3, xpad, ypad, NULL, qEnd, onQEndChanged);
      windowEntry = createTextEntryFromInt(dotterWindow, table, 3, 2, xpad, ypad, "Sliding _window size: ", dotplotGetSlidingWinSize(properties->dotplot), onSlidingWinSizeChanged);
    }
  else
    {
      createTextEntryFromInt(dotterWindow, table, 2, 2, xpad, ypad, "_Horizontal range: ", qStart, onQStartChanged);
      createTextEntryFromInt(dotterWindow, table, 2, 3, xpad, ypad, NULL, qEnd, onQEndChanged);
      createTextEntryFromInt(dotterWindow, table, 3, 2, xpad, ypad, "_Vertical range: ", sStart, onSStartChanged);
      createTextEntryFromInt(dotterWindow, table, 3, 3, xpad, ypad, NULL, sEnd, onSEndChanged);
      windowEntry = createTextEntryFromInt(dotterWindow, table, 4, 2, xpad, ypad, "Sliding _window size: ", dotplotGetSlidingWinSize(properties->dotplot), onSlidingWinSizeChanged);
    }

  if (windowEntry)
    gtk_widget_set_tooltip_text(windowEntry, "The size of the sliding window used to average pairwise scores. Note that this causes the matrix to be recalculated, which may time a long time for a large plot.");
}


/* Create the display control widgets for the settings dialog */
static void settingsDialogDisplayControls(GtkWidget *dialog, GtkWidget *dotterWindow, const int border)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotplotProperties *dotplotProperties = dotplotGetProperties(properties->dotplot);
  
  /* Put everything in a vbox inside a frame */
  GtkWidget *frame = gtk_frame_new("Display");
  gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), frame);
  gtk_container_set_border_width(GTK_CONTAINER(frame), border); 
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);

  /* Create a check box for toggling splice sites on and off */
  GtkWidget *splicesBtn = gtk_check_button_new_with_mnemonic("Highlight _splice sites");
  gtk_widget_set_tooltip_text(splicesBtn, "For known high-scoring pairs, highlight splice-sites in the alignment tool");
  gtk_container_add(GTK_CONTAINER(vbox), splicesBtn);
  gboolean spliceSitesOn = alignmentToolGetSpliceSitesOn(properties->alignmentTool);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(splicesBtn), spliceSitesOn);
  widgetSetCallbackData(splicesBtn, onSetSpliceSitesOn, dotterWindow);
  
  /* Create a check box for toggling breaklines on and off. If breaklines are
   * off at startup then it means that there are not multiple sequences, so
   * the option is not applicable. */
  static int disableBreaklines = -1; /* -1 for unset; 0 for false; 1 for true */

  if (disableBreaklines == -1)
    disableBreaklines = !dotplotProperties->breaklinesOn;
  
  GtkWidget *breaklinesBtn = gtk_check_button_new_with_mnemonic("Show _breaklines");
  gtk_container_add(GTK_CONTAINER(vbox), breaklinesBtn);
  gtk_widget_set_tooltip_text(breaklinesBtn, "Show breaklines between sequences when dottering multiple sequences that have been concatenated");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(breaklinesBtn), dotplotProperties->breaklinesOn);
  
  if (disableBreaklines)
    gtk_widget_set_sensitive(breaklinesBtn, FALSE);
  
  widgetSetCallbackData(breaklinesBtn, onSetBreaklinesOn, dotterWindow);
  
  /* Add buttons to allow the user to turn off hoz/vert annotation labels */
  GtkWidget *hozBtn = gtk_check_button_new_with_mnemonic("Show _horizontal sequence labels");
  gtk_widget_set_tooltip_text(hozBtn, "Show labels for each breakline between multiple sequences on the horizontal axis");
  gtk_container_add(GTK_CONTAINER(vbox), hozBtn);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hozBtn), dotplotProperties->hozLabelsOn);
  widgetSetCallbackData(hozBtn, onSetHozLabelsOn, dotterWindow);

  GtkWidget *vertBtn = gtk_check_button_new_with_mnemonic("Show _vertical sequence labels");
  gtk_widget_set_tooltip_text(vertBtn, "Show labels for each breakline between multiple sequences on the vertical axis");
  gtk_container_add(GTK_CONTAINER(vbox), vertBtn);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vertBtn), dotplotProperties->vertLabelsOn);
  widgetSetCallbackData(vertBtn, onSetVertLabelsOn, dotterWindow);
  
}


/* Pop up a dialog to allow the user to set the dotter parameters */
static void showSettingsDialog(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  const DotterDialogId dialogId = DOTDIALOG_SETTINGS;
  GtkWidget *dialog = getPersistentDialog(dwc->dialogList, dialogId);
  
  if (!dialog)
    {
      /* Create the dialog */
      char *title = g_strdup_printf("%sSettings", dotterGetTitlePrefix(dwc->dotterCtx));

      dialog = gtk_dialog_new_with_buttons(title,
                                           GTK_WINDOW(dotterWindow), 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           GTK_STOCK_CANCEL,
                                           GTK_RESPONSE_REJECT,
                                           NULL);

      g_free(title);
      
      /* These 2 calls are required to make the dialog persistent... */
      addPersistentDialog(dwc->dialogList, dialogId, dialog);
      g_signal_connect(dialog, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseSettingsDialog), GINT_TO_POINTER(TRUE));
    }
  else
    {
      /* Clear contents and re-create */
      dialogClearContentArea(GTK_DIALOG(dialog));
    }

  /* Create the contents */
  const int border = 12; /* border around individual sections */
  settingsDialogParamControls(dialog, dotterWindow, border);
  settingsDialogDisplayControls(dialog, dotterWindow, border);

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  gtk_widget_show_all(dialog);
  gtk_window_present(GTK_WINDOW(dialog));
  
  return;
}


/* Redraw the main dotter window, the alignment tool and the greyramp tool. Re-calculates borders etc. */
static void redrawAll(GtkWidget *dotterWindow, gpointer data)
{
  gtk_widget_queue_draw(dotterWindow);
  
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  /* Check the range values are the correct way round. */
  sortValues(&dwc->refSeqRange.min, &dwc->refSeqRange.max, TRUE);
  sortValues(&dwc->matchSeqRange.min, &dwc->matchSeqRange.max, TRUE);
  
  if (properties)
    {
      gtk_widget_queue_draw(properties->greyrampTool);
      gtk_widget_queue_draw(properties->alignmentTool);
      recalcDotplot(properties->dotplot);
    }
}


/* Refresh the main dotter window, the alignment tool and the greyramp tool. Clears any cached
 * pixmaps but does not recalculate borders etc. */
static void refreshAll(GtkWidget *dotterWindow, gpointer data)
{
  callFuncOnAllChildWidgets(dotterWindow, (gpointer)widgetClearCachedDrawable);
  gtk_widget_queue_draw(dotterWindow);

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  if (properties)
    {
      gtk_widget_queue_draw(properties->greyrampTool);
      callFuncOnAllChildWidgets(properties->alignmentTool, (gpointer)widgetClearCachedDrawable);
      callFuncOnAllChildWidgets(properties->dotplot, (gpointer)widgetClearCachedDrawable);
    }
}


static void readmtx(int MATRIX[24][24], char *mtxfile)
{
    FILE *fil;
    int row, col;
    char line[1025] = "#", *p;
    
    char *mtxfileText = g_strdup_printf("%s/%s", getenv("BLASTMAT"), mtxfile);
  
    if (!(fil = fopen(mtxfile, "r")) &&
        !(fil = fopen(mtxfileText, "r")))
      {
        char *msg = g_strdup_printf("Failed to open score matrix file %s - not found in ./ or $BLASTMAT/.\n", mtxfile);
        g_error(msg, mtxfile);
        g_free(msg);
      }
  
    g_free(mtxfileText);
    mtxfileText = NULL;
    
    /* Ignore header ... */
    while (!feof(fil) && *line == '#') 
        fgets(line, 1024, fil);

    /* Read in the pairwise scores */
    for (row = 0; row < 24; row++)
    {
        if (!(fgets(line, 1024, fil)))
            g_error("Wrong number of rows in matrix file: %d (should be 24).\n", row);

        p = strtok(line, " \t\n");
        for (col = 0; col < 24; col++) 
        {
            while (*p == '*' || isalpha((int) *p))
                p = strtok(NULL, " \t\n");
          
            if (!p) 
              g_error("Error on row %d in matrix file.\n", row);

            MATRIX[row][col] = atoi(p);

            p = strtok(NULL, " \t\n");
        }
    }

    g_message("I read your score matrix %s.\n", mtxfile);
    fclose(fil);
}

static void mtxcpy(int dest[24][24], int src[24][24])
{
    int i, j;

    for (i = 0 ; i < 24 ; i++)
        for (j = 0 ; j < 24 ; j++)
            dest[i][j] = src[i][j];
}


static void DNAmatrix(int mtx[24][24])
{
    int i, j;

    for (i = 0 ; i < 6 ; i++)
        for (j = 0 ; j < 6 ; j++) {
            if ( i < 4 && j < 4) 
                mtx[i][j] = (i == j ? 5 : -4);
            else 
                mtx[i][j] = -4;
        }
}


/* Utility to get the length of the given GArray. Returns 0 if array is null. */
//static int gArrayGetLen(GArray *array)
//{
//  return (array ? array->len : 0);
//}


/***********************************************************
 *               Show/hide parts of the view               *
 ***********************************************************/

/* Show/hide the greyramp tool, bringing it to the front */
static void showHideGreyrampTool(GtkWidget *dotterWindow, const gboolean show)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  if (properties && properties->greyrampWindow && GTK_IS_WIDGET(properties->greyrampWindow))
    {
      /* Get the parent widget for the greyramp too. This is the container if docked or the
       * window otherwise */
      GtkWidget *parent = properties->windowsDocked ? properties->greyrampContainer : properties->greyrampWindow;

      if (parent && show)
        {
          /* Show it, and bring it to the front if it's a toplevel window */
          gtk_widget_show_all(parent);
      
          if (GTK_IS_WINDOW(parent))
            gtk_window_present(GTK_WINDOW(parent));
        }
      else if (parent)
        {
          /* Hide it */
          gtk_widget_hide(parent);
        }
    }
}

/* Show/hide the alignment tool, bringing it to the front */
static void showHideAlignmentTool(GtkWidget *dotterWindow, const gboolean show)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  if (properties && properties->alignmentWindow && GTK_IS_WIDGET(properties->alignmentWindow))
    {
      /* Get the parent widget for the alignment too. This is the container if docked or the
       * window otherwise */
      GtkWidget *parent = properties->windowsDocked ? properties->alignmentContainer : properties->alignmentWindow;

      if (parent && show)
        {
          /* Show it, and bring it to the front if it's a toplevel window */
          gtk_widget_show_all(parent);
      
          if (GTK_IS_WINDOW(parent))
            gtk_window_present(GTK_WINDOW(parent));
        }
      else if (parent)
        {
          /* Hide it */
          gtk_widget_hide(parent);
        }

    }
}

/* Bring the main dotter window to the front */
static void showDotterWindow(GtkWidget *dotterWindow)
{
  gtk_widget_show_all(dotterWindow);
  gtk_window_present(GTK_WINDOW(dotterWindow));
}


/* Move the given widget from source to dest */
static void reparentWidget(GtkWidget *widget, GtkContainer *source, GtkContainer *dest, const gboolean show)
{
  g_return_if_fail(widget && source && dest);

  /* If the container ref to the widget is removed it might be destroyed, so make sure there is a
   * ref to it. */
  g_object_ref(widget);
  gtk_container_remove(source, widget);
  gtk_container_add(dest, widget);
  g_object_unref(widget);

  /* Hide the old widget */
  gtk_widget_hide(GTK_WIDGET(source));

  /* Show the new widget (only if it's toggled on) */
  if (show)
    gtk_widget_show_all(GTK_WIDGET(dest));
}


static void dotterToggleDockWindows(GtkWidget *dotterWindow)
{
  g_return_if_fail(dotterWindow);

  DotterProperties *properties = dotterGetProperties(dotterWindow);
  gboolean greyrampVisible = getToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleGreyramp");
  gboolean alignmentVisible = getToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleAlignment");

  if (properties->windowsDocked)
    {
      reparentWidget(properties->alignmentTool, GTK_CONTAINER(properties->alignmentContainer), GTK_CONTAINER(properties->alignmentWindow), alignmentVisible);
      reparentWidget(properties->greyrampTool, GTK_CONTAINER(properties->greyrampContainer), GTK_CONTAINER(properties->greyrampWindow), greyrampVisible);
    }
  else
    {
      reparentWidget(properties->alignmentTool, GTK_CONTAINER(properties->alignmentWindow), GTK_CONTAINER(properties->alignmentContainer), alignmentVisible);
      reparentWidget(properties->greyrampTool, GTK_CONTAINER(properties->greyrampWindow), GTK_CONTAINER(properties->greyrampContainer), greyrampVisible);
    }

  properties->windowsDocked = !properties->windowsDocked;
}

/***********************************************************
 *                       Help Dialog                       *
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
                        "comments", dotterGetCommentsString(), 
                        "copyright", dotterGetCopyrightString(),
                        "license", dotterGetLicenseString(),
                        "name", dotterGetAppName(),
                        "version", dotterGetVersionString(),
                        "website", dotterGetWebSiteString(),
                        NULL) ;
#endif
  
  return ;
}


static void showHelpDialog(GtkWidget *dotterWindow)
{
  GError *error = NULL;
  
  /* The docs should live in /share/doc/seqtools/, in the same parent
   * directory that our executable's 'bin' directory is in. Open the 'quick
   * start' page. */
  char rel_path[100] = "../share/doc/seqtools/dotter_quick_start.html";
  
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
 *                          Events                         *
 ***********************************************************/

static void onQuitMenu(GtkAction *action, gpointer data)
{
  /* Close all associated windows */
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterContext *dc = dotterGetContext(dotterWindow);
  dotterContextCloseAllWindows(dc);
}

static void onCloseMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  closeWindow(dotterWindow);
}

static void onSavePlotMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  GError *error = NULL;
  savePlot(properties->dotplot, NULL, NULL, &error);
  
  prefixError(error, "Error saving plot. ");
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
}

static void onExportPlotMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  /* Set the background colour to something sensible for printing */
  GdkColor *defaultBgColor = getGdkColor(DOTCOLOR_BACKGROUND, properties->dotterWinCtx->dotterCtx->defaultColors, FALSE, TRUE);
  setWidgetBackgroundColor(dotterWindow, defaultBgColor);
  redrawAll(dotterWindow, NULL);
  
  GError *error = NULL;
  exportPlot(properties->dotplot, GTK_WINDOW(dotterWindow), NULL, &error);
  
  prefixError(error, "Error exporting plot. ");
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

  /* Revert the background colour */
  onPrintColorsChanged(dotterWindow);

  redrawAll(dotterWindow, NULL);
}

static void onPrintMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  printDotterWindow(dotterWindow);
}

static void onSettingsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showSettingsDialog(dotterWindow);
}

static void onHelpMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showHelpDialog(dotterWindow);
}

static void onAboutMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  showAboutDialog(dotterWindow);
}

/* Callback called when the user selects the 'copy horizontal coord' menu option */
static void onCopyHCoordMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  copyIntToDefaultClipboard(properties->dotterWinCtx->refCoord);
}

/* Callback called when the user selects the 'copy vertical coord' menu option */
static void onCopyVCoordMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  copyIntToDefaultClipboard(properties->dotterWinCtx->matchCoord);
}

static void onShowHideGreyrampMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);

  gboolean show = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));

  showHideGreyrampTool(dotterWindow, show);
}

static void onShowHideAlignmentMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);

  gboolean show = gtk_toggle_action_get_active(GTK_TOGGLE_ACTION(action));

  showHideAlignmentTool(dotterWindow, show);
}

static void onToggleCrosshairMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairOn(properties->dotplot);
}

static void onToggleCoordsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairCoordsOn(properties->dotplot);
}

static void onToggleFullscreenMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  toggleCrosshairFullscreen(properties->dotplot);
}

static void onTogglePixelmapMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  dotplotTogglePixelmap(properties->dotplot);
}

static void onToggleGridMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  dotplotToggleGrid(properties->dotplot);
}

static void onToggleHspMode(GtkRadioAction *action, GtkRadioAction *current, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  
  const DotterHspMode hspMode = (const DotterHspMode)gtk_radio_action_get_current_value(current);
  setHspMode(properties->dotplot, hspMode);
}

static void onPrintColorsChanged(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  /* Refresh the background colors for both the main window and the alignment tool */
  GdkColor *defaultBgColor = getGdkColor(DOTCOLOR_BACKGROUND, dwc->dotterCtx->defaultColors, FALSE, dwc->usePrintColors);
  
  setWidgetBackgroundColor(dotterWindow, defaultBgColor);
  setWidgetBackgroundColor(properties->alignmentTool, defaultBgColor);
  
  /* Redraw everything */
  refreshAll(dotterWindow, NULL);
}

static void onToggleUsePrintColorsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  /* Toggle the flag*/
  dwc->usePrintColors = !dwc->usePrintColors;
  onPrintColorsChanged(dotterWindow);
}

static void onToggleBumpExonsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  dotplotToggleBumpExons(properties->dotplot);
}

static void onToggleDockWindowsMenu(GtkAction *action, gpointer data)
{
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  dotterToggleDockWindows(dotterWindow);
}


/* Mouse button handler */
static gboolean onButtonPressDotter(GtkWidget *window, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 3) /* right click */
    {
      GtkMenu *contextMenu = GTK_MENU(data);
      gtk_menu_popup (contextMenu, NULL, NULL, NULL, NULL, event->button, event->time);
    }
  else if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      /* If the dot-plot was clicked the selected coords will have changed. Perform required updates. */
      updateOnSelectedCoordsChanged(window);
      handled = TRUE;
    }
  
  return handled;
}


/* Mouse move handler */
static gboolean onMouseMoveDotter(GtkWidget *window, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  if (event->state & GDK_BUTTON1_MASK)  /* left-drag */
    {
      /* If the dot-plot was clicked the selected coords will have changed. Perform required updates. */
      updateOnSelectedCoordsChanged(window);
      handled = TRUE;
    }

  return handled;
}


/* Get the factor to convert between nucleotide and display (e.g. peptide) coords, or 1 if no
 * conversion is necessary */
int getResFactor(DotterContext *dc, const gboolean horizontal)
{
  return horizontal && dc->blastMode == BLXMODE_BLASTX ? dc->numFrames : 1;
}


/* Returns true if display coordinates should be shown negated for the given scale. */
static gboolean negateDisplayCoord(DotterContext *dc, const gboolean horizontal)
{
  gboolean result = FALSE;
  
  if (dc->negateCoords)
    {
    if (horizontal)
      result = dc->hozScaleRev;
    else 
      result = dc->vertScaleRev;
    }
  
  return result;
}


/* Convert the given coord to a display coord (just negates it for the display if necessary) */
int getDisplayCoord(const int coordIn, DotterContext *dc, const gboolean horizontal)
{
  int result = coordIn;
  
  if (negateDisplayCoord(dc, horizontal))
    result *= -1;
  
  return result;
}


/* Get the currently-selected (i.e. crosshair) coord for the horizontal or 
 * vertical sequence, as indicated by the bool */
int getSelectedCoord(DotterWindowContext *dwc, const gboolean horizontal)
{
  return (horizontal ? dwc->refCoord : dwc->matchCoord);
}


/* Get the start coord of the display range for the given sequence */
int getStartCoord(DotterWindowContext *dwc, const gboolean horizontal)
{
  int result = UNSET_INT;
  
  if (horizontal)
    result = dwc->dotterCtx->hozScaleRev ? dwc->refSeqRange.max : dwc->refSeqRange.min;
  else
    result = dwc->dotterCtx->vertScaleRev ? dwc->matchSeqRange.max : dwc->matchSeqRange.min;
  
  return result;
}

/* Get the end coord of the display range for the given sequence */
int getEndCoord(DotterWindowContext *dwc, const gboolean horizontal)
{
  int result = UNSET_INT;
  
  if (horizontal)
    result = dwc->dotterCtx->hozScaleRev ? dwc->refSeqRange.min : dwc->refSeqRange.max;
  else
    result = dwc->dotterCtx->vertScaleRev ? dwc->matchSeqRange.min : dwc->matchSeqRange.max;
  
  return result;
}



/* Set the start coord of the display range for the given sequence */
static gboolean setStartCoord(GtkWidget *dotterWindow, DotterWindowContext *dwc, const gboolean horizontal, const int newValue)
{
  gboolean changed = FALSE; 
  int *valueToUpdate = NULL;
  
  if (horizontal)
    valueToUpdate = (dwc->dotterCtx->hozScaleRev ? &dwc->refSeqRange.max : &dwc->refSeqRange.min);
  else
    valueToUpdate = (dwc->dotterCtx->vertScaleRev ? &dwc->matchSeqRange.max : &dwc->matchSeqRange.min);

  if (valueToUpdate && *valueToUpdate != newValue)
    {
      *valueToUpdate = newValue;
      changed = TRUE;
    }
  
  return changed;
}

/* Set the end coord of the display range for the given sequence */
static gboolean setEndCoord(GtkWidget *dotterWindow, DotterWindowContext *dwc, const gboolean horizontal, const int newValue)
{
  gboolean changed = FALSE;
  int *valueToUpdate = NULL;
  
  if (horizontal)
    valueToUpdate = (dwc->dotterCtx->hozScaleRev ? &dwc->refSeqRange.min : &dwc->refSeqRange.max);
  else
    valueToUpdate = (dwc->dotterCtx->vertScaleRev ? &dwc->matchSeqRange.min : &dwc->matchSeqRange.max);
  
  if (valueToUpdate && *valueToUpdate != newValue)
    {
      *valueToUpdate = newValue;
      changed = TRUE;
    }
  
  return changed;
}

/* Move the given sequence coord by the given number of coords (which can be negative to move 
 * in the decreasing direction. 'horizontal' indicates whether it's the horizontal or vertical
 * sequence that we're modifying and 'reverse' indicates whether that sequence's scale is shown reversed. */
static void incrementCoord(GtkWidget *dotterWindow, 
                           DotterContext *dc,
                           int *coord, 
                           const gboolean reverse, 
                           const gboolean horizontal, 
                           const gboolean convertCoords,
                           const int numCoords)
{
  const int incValue = convertCoords ? numCoords * getResFactor(dc, horizontal) : numCoords;

  if (reverse)
    {
      *coord -= incValue;
    }
  else
    {
      *coord += incValue;
    }
  
  updateOnSelectedCoordsChanged(dotterWindow);
}


/* Handle Q key press (Ctrl-Q => close all windows) */
static gboolean onKeyPressQ(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    dotterContextCloseAllWindows(dotterGetContext(dotterWindow));
  
  return ctrlModifier;
}

/* Handle W key press (Ctrl-W => close window) */
static gboolean onKeyPressW(GtkWidget *widget, const gboolean ctrlModifier)
{
  gboolean handled = FALSE;

  if (ctrlModifier)
    {
      handled = closeWindow(widget);
    }
  
  return handled;
}

/* Handle H key press (Ctrl-H => show help dialog) */
static gboolean onKeyPressH(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showHelpDialog(dotterWindow);
  
  return ctrlModifier;
}


/* Handle S key press (Ctrl-S => show settings dialog) */
static gboolean onKeyPressS(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showSettingsDialog(dotterWindow);
    
  return ctrlModifier;
}

/* Handle G key press (Ctrl-G => show/hide greyramp tool) */
static gboolean onKeyPressG(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    {
      DotterProperties *properties = dotterGetProperties(dotterWindow);
      gboolean active = getToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleGreyramp");

      /* Toggle the visiblity */
      setToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleGreyramp", !active);
    }
  
  return ctrlModifier;
}

/* Handle A key press (Ctrl-A => show/hide alignment tool) */
static gboolean onKeyPressA(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    {
      DotterProperties *properties = dotterGetProperties(dotterWindow);
      gboolean active = getToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleAlignment");

      /* Toggle the visiblity */
      setToggleMenuStatus(properties->dotterWinCtx->actionGroup, "ToggleAlignment", !active);
    }
  
  return ctrlModifier;
}

/* Handle D key press (Ctrl-D => show dotplot) */
static gboolean onKeyPressD(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    showDotterWindow(dotterWindow);
  
  return ctrlModifier;
}

/* Handle K key press (Ctrl-K => dock/undock windows) */
static gboolean onKeyPressK(GtkWidget *dotterWindow, const gboolean ctrlModifier)
{
  if (ctrlModifier)
    {
      DotterProperties *properties = dotterGetProperties(dotterWindow);
      setToggleMenuStatus(properties->dotterWinCtx->actionGroup, "DockWindows", !properties->windowsDocked);
    }
  
  return ctrlModifier;
}

/* Handle up/down key presses */
static gboolean onKeyPressUpDown(GtkWidget *dotterWindow, const gboolean isUp, const gboolean modifier)
{
  /* Increment/decrement the vertical (i.e. match) sequence coord */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;
  
  const int numCoords = isUp ? -1 : 1;
  
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, numCoords);  
  return TRUE;
}

/* Handle left/right key presses */
static gboolean onKeyPressLeftRight(GtkWidget *dotterWindow, const gboolean isLeft, const gboolean modifier)
{
  /* Increment/decrement the horizontal (i.e. reference) sequence coord */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  const int numCoords = isLeft ? -1 : 1;
  
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, numCoords);
  return TRUE;
}

/* Handle comma/period key presses */
static gboolean onKeyPressCommaPeriod(GtkWidget *dotterWindow, const gboolean isComma, const gboolean modifier)
{
  /* Increment/decrement both sequence coords */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, isComma ? -1 : 1);
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, isComma ? -1 : 1);
  return TRUE;
}

/* Handle left/right square bracket presses */
static gboolean onKeyPressLeftRightBracket(GtkWidget *dotterWindow, const gboolean isLeft, const gboolean modifier)
{
  /* Increment the horizontal and decrement the vertical, or vice versa */
  DotterProperties *properties = dotterGetProperties(dotterWindow);
  DotterWindowContext *dwc = properties->dotterWinCtx;

  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->refCoord, dwc->dotterCtx->hozScaleRev, TRUE, !modifier, isLeft ? -1 : 1);
  incrementCoord(dotterWindow, dwc->dotterCtx, &dwc->matchCoord, dwc->dotterCtx->vertScaleRev, FALSE, !modifier, isLeft ? 1 : -1);
  return TRUE;
}


/* Key presses applicable to all windows */
gboolean onKeyPressDotter(GtkWidget *widget, GdkEventKey *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *dotterWindow = GTK_WIDGET(data);
  
  const gboolean ctrlModifier = (event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK;  
  
  switch (event->keyval)
    {
      case GDK_Q:   /* fall through */
      case GDK_q:   handled = onKeyPressQ(dotterWindow, ctrlModifier);              break;

      case GDK_W:   /* fall through */
      case GDK_w:   handled = onKeyPressW(widget, ctrlModifier);                    break;
        
      case GDK_H:   /* fall through */
      case GDK_h:   handled = onKeyPressH(dotterWindow, ctrlModifier);              break;

      case GDK_S:   /* fall through */
      case GDK_s:   handled = onKeyPressS(dotterWindow, ctrlModifier);              break;

      case GDK_G:   /* fall through */
      case GDK_g:   handled = onKeyPressG(dotterWindow, ctrlModifier);              break;

      case GDK_A:   /* fall through */
      case GDK_a:   handled = onKeyPressA(dotterWindow, ctrlModifier);              break;

      case GDK_D:   /* fall through */
      case GDK_d:   handled = onKeyPressD(dotterWindow, ctrlModifier);              break;

      case GDK_K:   /* fall through */
      case GDK_k:   handled = onKeyPressK(dotterWindow, ctrlModifier);              break;

      default: break;
  }
  
  return handled;
}


/* Key press handler for dotter windows that move the selected coords (including main window and
 *  alignment tool.). */
gboolean onKeyPressDotterCoords(GtkWidget *widget, GdkEventKey *event, gpointer data)
{
  gboolean handled = onKeyPressDotter(widget, event, data);
  
  if (!handled)
    {
      GtkWidget *dotterWindow = GTK_WIDGET(data);

      const gboolean shiftModifier = (event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK;
      //  const gboolean altModifier = (event->state & GDK_MOD1_MASK) == GDK_MOD1_MASK;
      
      switch (event->keyval)
        {
          case GDK_Up:            handled = onKeyPressUpDown(dotterWindow, TRUE, shiftModifier);        break;
          case GDK_Down:          handled = onKeyPressUpDown(dotterWindow, FALSE, shiftModifier);       break;
            
          case GDK_Left:          handled = onKeyPressLeftRight(dotterWindow, TRUE, shiftModifier);     break;
          case GDK_Right:         handled = onKeyPressLeftRight(dotterWindow, FALSE, shiftModifier);    break;
            
          case GDK_comma:         handled = onKeyPressCommaPeriod(dotterWindow, TRUE, shiftModifier);   break;
          case GDK_less:          handled = onKeyPressCommaPeriod(dotterWindow, TRUE, shiftModifier);   break;
          case GDK_period:        handled = onKeyPressCommaPeriod(dotterWindow, FALSE, shiftModifier);  break;
          case GDK_greater:       handled = onKeyPressCommaPeriod(dotterWindow, FALSE, shiftModifier);  break;
            
          case GDK_bracketleft:   handled = onKeyPressLeftRightBracket(dotterWindow, TRUE, shiftModifier);        break;
          case GDK_braceleft:     handled = onKeyPressLeftRightBracket(dotterWindow, TRUE, shiftModifier);        break;
          case GDK_bracketright:  handled = onKeyPressLeftRightBracket(dotterWindow, FALSE, shiftModifier);       break;
          case GDK_braceright:    handled = onKeyPressLeftRightBracket(dotterWindow, FALSE, shiftModifier);       break;
            
          default: break;
        }
    }
  
  return handled;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* This creates BlxColumnInfo entries for each "column" of data (a misnomer in dotter
 * because we don't currently display this data in columns! The format comes from blixem which does
 * display the data in columns, and we use the same parser, and the same data structures, for dotter). */
GList* dotterCreateColumns()
{
  GList *columnList = NULL;
  
  /* Create the columns' data structs. The columns appear in the order
   * that they are added here. */
  blxColumnCreate(BLXCOL_SEQNAME, FALSE, "Name", G_TYPE_STRING, NULL, 0, TRUE, TRUE, FALSE, FALSE, FALSE, "Name", NULL, NULL, &columnList);

  return columnList;
}


/* Create the UI manager for the menus */
static GtkUIManager* createUiManager(GtkWidget *window, const DotterHspMode hspMode, GtkActionGroup **actionGroup_out)
{
  GtkActionGroup *action_group = gtk_action_group_new ("MenuActions");
  
  /* If the HSPs are on initially then we're in hspOnly mode, so we don't show the pixmap; therefore,
   * set the toggled state of the pixelmap menu option (first in the toggleMenuEntries list) to be false
   * if HSPs are on. */
  toggleMenuEntries[0].is_active = (hspMode == 0);  
  
  gtk_action_group_add_actions(action_group, menuEntries, G_N_ELEMENTS (menuEntries), window);
  gtk_action_group_add_toggle_actions(action_group, toggleMenuEntries, G_N_ELEMENTS (toggleMenuEntries), window);
  gtk_action_group_add_radio_actions(action_group, radioMenuEntries, G_N_ELEMENTS (radioMenuEntries), hspMode, G_CALLBACK(onToggleHspMode), window);
  
  GtkUIManager *ui_manager = gtk_ui_manager_new ();
  gtk_ui_manager_insert_action_group (ui_manager, action_group, 0);
  gtk_ui_manager_set_add_tearoffs(ui_manager, TRUE);
  
  GtkAccelGroup *accel_group = gtk_ui_manager_get_accel_group (ui_manager);
  gtk_window_add_accel_group (GTK_WINDOW (window), accel_group);

  if (actionGroup_out)
    *actionGroup_out = action_group;

  return ui_manager;
}


/* Create a menu. Optionally add it to the given menu bar, if menuBar is not null, with the given label */
static GtkWidget* createDotterMenu(GtkWidget *window, 
                                   const char *menuDescription, 
                                   const char *path, 
                                   GtkUIManager *ui_manager)
{
  GError *error = NULL;
  if (!gtk_ui_manager_add_ui_from_string (ui_manager, menuDescription, -1, &error))
    {
      prefixError(error, "Building menus failed: ");
      reportAndClearIfError(&error, G_LOG_LEVEL_ERROR);
    }
  
  GtkWidget *menu = gtk_ui_manager_get_widget (ui_manager, path);
  
  return menu;
}


static GtkWidget* createDotterWindow(DotterContext *dc, 
                                     DotterWindowContext *dwc,
                                     const DotterHspMode hspMode, 
                                     GtkWidget *dotplot,
                                     GtkWidget *dotplotContainer, 
                                     GtkWidget *greyrampContainer,
                                     GtkWidget *alignmentContainer,
                                     const char *exportFileName,
                                     char *windowColor)
{ 
  DEBUG_ENTER("createDotterWindow");

  GtkWidget *dotterWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_widget_set_name(dotterWindow, MAIN_WINDOW_NAME);
  
  char *title = g_strdup_printf("%s%s vs. %s", dotterGetTitlePrefix(dc), dc->refSeqName, dc->matchSeqName);
  gtk_window_set_title(GTK_WINDOW(dotterWindow), title);
  g_free(title);
  
  /* Set the parent window in the message handlers' data, now we know it */
  dc->msgData->parent = GTK_WINDOW(dotterWindow);
  
  /* Create the menu bar, and a right-click context menu */
  dwc->uiManager = createUiManager(dotterWindow, hspMode, &dwc->actionGroup);
  GtkWidget *menuBar = createDotterMenu(dotterWindow, mainMenuDescription, "/MenuBar", dwc->uiManager);
  GtkWidget *contextMenu = createDotterMenu(dotterWindow, mainMenuDescription, "/ContextMenu", dwc->uiManager);
  
  blxSetWidgetColor(menuBar, windowColor);

  /* We'll set the default window size based on the dotplot/exon widget size, up to a 
   * max based on screen size. */
  GdkScreen *screen = gtk_widget_get_screen(dotterWindow);
  const int maxWidth = gdk_screen_get_width(screen) * MAX_WINDOW_WIDTH_FRACTION;
  const int maxHeight = gdk_screen_get_height(screen) * MAX_WINDOW_HEIGHT_FRACTION;

  const int exonViewHeight = 2 * (DEFAULT_EXON_HEIGHT + (2 * DEFAULT_EXON_YPAD));
  DotplotProperties *dotplotProperties = dotplotGetProperties(dotplot);
  const int dotplotWidth = getDotplotWidth(dotplot, dotplotProperties);
  int greyrampWidth = 400; /* roughly */
  const int alignmentToolHeight = 300; /* roughly */

  /* We'll base the layout on the relative size of the dotplot to the window size - we'll place
   * the greyramp tool on the same row as the dotplot if it'll fit and we can still display the
   * entire dotplot (not worrying about exons for now). If it won't fit, we'll place the greyramp
   * tool on the row below, adjacent to the alignment tool. */
  /*! \todo Ideally we'd adjust the layout after the user changes the settings, i.e. the zoom or
   * the range of sequence displayed */
  gboolean maximise_dotplot = FALSE;

  if (dotplotWidth > maxWidth - greyrampWidth)
    {
      maximise_dotplot = TRUE;
      greyrampWidth = 0; /* on a different row so don't include it in the width calculation */
    }

  int width = dotplotWidth + exonViewHeight + greyrampWidth;
  int height = getDotplotHeight(dotplot, dotplotProperties) + exonViewHeight + alignmentToolHeight;
  width = min(width, maxWidth);
  height = min(height, maxHeight);
  
  gtk_window_set_default_size(GTK_WINDOW(dotterWindow), width, height);

  /* Put the widgets in a table */
  const int numRows = 3;
  const int numCols = 2;
  int padding = 0;
  int row = 0;

  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_container_add(GTK_CONTAINER(dotterWindow), GTK_WIDGET(table));

  gtk_table_attach(table, menuBar, 0, numCols, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_SHRINK, padding, padding);
  ++row;

  if (maximise_dotplot)
    {
      /* dotplot spans all columns; alignment tool + greyramp on same row. */
      gtk_table_attach(table, dotplotContainer, 0, numCols, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, padding, padding);
      ++row;
      gtk_table_attach(table, alignmentContainer, 0, 1, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_SHRINK, padding, padding);
      gtk_table_attach(table, greyrampContainer, 1, 2, row, row + 1, GTK_SHRINK, GTK_FILL, padding, padding);
    }
  else
    {
      /* dotplot and greyramp on same row; alignment tool spans all columns */
      gtk_table_attach(table, dotplotContainer, 0, 1, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND, padding, padding);
      gtk_table_attach(table, greyrampContainer, 1, 2, row, row + 1, GTK_SHRINK, GTK_FILL, padding, padding);
      ++row;
      gtk_table_attach(table, alignmentContainer, 0, numCols, row, row + 1, GTK_FILL | GTK_EXPAND, GTK_SHRINK, padding, padding);
    }

  gtk_widget_add_events(dotterWindow, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(dotterWindow, GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(dotterWindow), "key-press-event", G_CALLBACK(onKeyPressDotterCoords), dotterWindow);
  g_signal_connect(G_OBJECT(dotterWindow), "button-press-event", G_CALLBACK(onButtonPressDotter), contextMenu);
  g_signal_connect(G_OBJECT(dotterWindow), "motion-notify-event", G_CALLBACK(onMouseMoveDotter), NULL);
  
  gtk_widget_show_all(dotterWindow);
  
  DEBUG_EXIT("createDotterWindow returning ");
  return dotterWindow;
}


/***********************************************************
 *                       Utilities                         *
 ***********************************************************/


/* Returns a string which is the name of the Dotter application. */
const char *dotterGetAppName(void)
{
  return DOTTER_TITLE ;
}

/* Returns a string which is the prefix to window titles. */
const char *dotterGetTitlePrefix(DotterContext *dc)
{
  return dc->abbrevTitle ? DOTTER_PREFIX_ABBREV : DOTTER_PREFIX ;
}

/* Returns a copyright string for the Dotter application. */
const char *dotterGetCopyrightString(void)
{
  return DOTTER_COPYRIGHT_STRING ;
}

/* Returns the Dotter website URL. */
const char *dotterGetWebSiteString(void)
{
  return DOTTER_WEBSITE_STRING ;
}

/* Returns a comments string for the Dotter application. */
const char *dotterGetCommentsString(void)
{
  return DOTTER_COMMENTS_STRING() ;
}

/* Returns a license string for the dotter application. */
const char *dotterGetLicenseString(void)
{
  return DOTTER_LICENSE_STRING ;
}

/* Returns a string representing the Version/Release/Update of the Dotter code. */
const char *dotterGetVersionString(void)
{
  return DOTTER_VERSION_STRING ;
}

/* Utility to copy an integer value as a string to the default clipboard */
void copyIntToDefaultClipboard(const int val)
{
  char *displayText = convertIntToString(val);
  setDefaultClipboardText(displayText);
  g_free(displayText); 
}


/* Print the main dotter window */
static void printDotterWindow(GtkWidget *dotterWindow)
{
  DotterProperties *properties = dotterGetProperties(dotterWindow);

  /* Set the background colour to something sensible for printing */
  GdkColor *defaultBgColor = getGdkColor(DOTCOLOR_BACKGROUND, properties->dotterWinCtx->dotterCtx->defaultColors, FALSE, TRUE);
  setWidgetBackgroundColor(dotterWindow, defaultBgColor);
  redrawAll(dotterWindow, NULL);

  /* The crosshair on the dotplot does not get cached in the dotplot's drawable,
   * but we want it to show in the print, so draw it on now. */
  GtkWidget *dotplot = properties->dotplot;
  dotplotPrepareForPrinting(dotplot);
  
  /* Print the parent of the dotplot, because this contains the exon views as well.
   * Note that we don't want to print the scrolled window, because that will chop off
   * parts of the plot that are not currently visible; we want to print the whole plot) */
  GtkWidget *parent = gtk_widget_get_parent(dotplot);

  /* Do the print */
  DotterWindowContext *dwc = properties->dotterWinCtx;
  blxPrintWidget(parent, NULL, GTK_WINDOW(dotterWindow), &dwc->printSettings, &dwc->pageSetup, NULL, TRUE, PRINT_FIT_BOTH);
  
  /* Revert the background colour */
  onPrintColorsChanged(dotterWindow);

  /* Redraw the entire dotplot to make sure the crosshair we added gets cleared */
  redrawAll(dotterWindow, NULL);
}


/**************************** eof ***************************/
