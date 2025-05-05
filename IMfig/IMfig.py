"""
    5/5/2025  replaced global gv dictionary with args passed to functions 
    
    Copyright 2009-2025  Jody Hey
    IMfig makes an image file containing a figure of the population phylogeny
        in an Isolation-with-Migration framework
    To see the helpscreen run the following at a command prompt:
        python IMfig.py
    If desired IMfig can be run from within an editing environment
        by setting 'cmdstr' (see bottom of this file)
    Check releasedate.
    Last tested on python 3.10  but used to run fine on 3.6
    Read the documentation for details on running the program

usage: temp.py [-h] -i IMFILENAME [-a] [-b POPBOXSPACEADJ] [-c {j,p,n}] [-d] [-e] [-f FONT] [-g GLOBALSCALE]
               [-j ARROWHEIGHTADJ] [-k] [-l HEIGHTSCALE] [-m MOPTION] [-n ALTNAMEFILENAME]
               [-o OUTPUTFILENAME] [-p FONTSIZE] [-q] [-r] [-s] [-t LASTT_LOWER_Y] [-u] [-v]
               [-w WIDTHSCALAR] [-x XADJUST] [-y LOCALYSCALE] [-z]

IMfig program. Copyright 2009-2025 Jody Hey. Release Date Apr 18, 2025

options:
  -h, --help            show this help message and exit
  -i IMFILENAME, --input IMFILENAME
                        input file name
  -a, --label-ancestor-pops
                        include ancestral population #'s in plot
  -b POPBOXSPACEADJ, --box-spacing POPBOXSPACEADJ
                        adjust width spacing of population boxes, values > 0, default = 1
  -c {j,p,n}, --convert {j,p,n}
                        output format, default is eps, see also -w
                        -c j : make a jpeg file
                        -c p : make a pdf file
                        -c n : make a png file
  -d, --no-demographic-scale
                        do not use demographic scale information even if in input file
  -e, --even-split-times
                        space split times evenly (not proportional to time, no confidence intervals shown)
  -f FONT, --font FONT  font. Default=Arial. Use postscript fonts available on the computer
                        e.g. Arial, Helvetica, Times-roman, Courier
  -g GLOBALSCALE, --global-scale GLOBALSCALE
                        global plot scale sets the size of the plot, max = 1, default = 1
  -j ARROWHEIGHTADJ, --arrow-width ARROWHEIGHTADJ
                        arrow width, default = 1
  -k, --angled-names    print population names on an angle
  -l HEIGHTSCALE, --height-scale HEIGHTSCALE
                        expand/shrink height by a positive scalar, >1 means taller, <1 means shorter
  -m MOPTION, --migration MOPTION
                        options for printing of arrows and 2Nm values for migration:
                        -m x : do not print migration arrows
                        -m a : 2Nm migration arrows for all cases when both m > 0 and 2Nm > 0
                        -m s : 2Nm migration arrows only if m is statistically significant p <= 0.05 (default)
                        -m S : 2Nm migration arrows only if m is statistically significant p <= 0.01
                        -m # : "#" is a number, migration arrows appear when 2NM >= # (e.g. -m0.1)
  -n ALTNAMEFILENAME, --alt-names ALTNAMEFILENAME
                        file with alternative species names
  -o OUTPUTFILENAME, --output OUTPUTFILENAME
                        output file name, default is imfig_output.eps
  -p FONTSIZE, --font-size FONTSIZE
                        fontsize (default is 14 for full scale, default follows global scale)
  -q, --no-confidence-interval-boxes
                        no confidence interval boxes for population boxes printed
  -r, --no-confidence-interval-arrows
                        no confidence interval arrows for population boxes printed
  -s, --square          print square, rather than landscape
  -t LASTT_LOWER_Y, --time-height LASTT_LOWER_Y
                        relative height of oldest time point, values between 0 and 1
                        default value = 1/(# sampled populations+1)
  -u, --simple-colors   simple colors, blue for population boxes, red arrows for migration (default grayscale)
  -v, --color           multiple colors for population boxes, red arrows for migration (default grayscale)
  -w WIDTHSCALAR, --image-width WIDTHSCALAR
                        file image width, integer multiple of 720 pixels (only if using -c)
  -x XADJUST, --width-adjust XADJUST
                        expand/shrink width of plot by a positive scalar, >1 means wider, <1 means narrower
  -y LOCALYSCALE, --height-adjust LOCALYSCALE
                        adjust height of splittimes, relative to bottom of figure, max = 1.
  -z, --exclude-ghost   exclude the ghost population from the figure
"""

import math
import sys
import os
import string
import traceback 
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Union, TextIO
# some users won't have colormath
try:
    from colormath.color_objects import LabColor, sRGBColor
    from colormath.color_conversions import convert_color
    from colormath.color_objects import SpectralColor
    check_colormath = True
except ImportError:
    check_colormath = False

## some users won't have PIL
check_PIL = False
try:
    from PIL import Image
    check_PIL = True
except ImportError:
    check_PIL = False

# These constants are still at the top level of the module since they aren't 
# configuration parameters that change during execution
RELEASE_YEAR = "2025"
RELEASE_DATE = "May 5, 2025"

# Constants related to chart generation - these could potentially move to the Config class
# if they might need to be adjusted during execution
ARROWHEAD_WIDTH_DEFAULT = 0.01  # arrow size
POPBOX_SPACE_DEFAULT = 0.1  # spacing between population boxes
CURVEHEIGHT_DEFAULT = 0.03  # curvature of migration arrows
TFACTOR = 1.0  # a fudge factor for moving things to the right of splittime arrows

# The number of populations will still be a global variable as it's set during 
# file parsing and is widely used
numpops = 0

@dataclass
class Config:
    """
    Class to hold all configuration parameters for IMfig.
    Replaces the global variables dictionary approach.
    """
    # Basic configuration
    label_a_pops: bool = False
    simplecolor: bool = False
    rgbcolor: bool = False
    dosquare: bool = False
    eventimes: bool = False
    popboxcintervalboxes: bool = True
    popboxcintervalarrows: bool = True
    imfilename: str = "im_eps.txt"
    outputfilename: str = "imfig_output.eps"
    globalscale: float = 1.0
    font: str = "Arial"
    bifont: str = "Arial-BoldItalic"
    fontsize: int = 14
    fontfixed: bool = False
    line0y: float = 0.95
    lineINFy: float = 0.1
    localxscale: float = -1
    localyscale: float = -1
    arrowheightadj: float = 1.0
    maximumxpoint: float = 756.1
    maximumypoint: float = 576.1
    lastt_lower_y: float = -1
    set_lastt_lower_y: bool = True
    xadjust: Optional[float] = None
    heightscale: Optional[float] = None
    
    # Color settings
    blue: List[float] = field(default_factory=lambda: [0, 0, 1])
    red: List[float] = field(default_factory=lambda: [1, 0, 0])
    black: List[float] = field(default_factory=lambda: [0, 0, 0])
    darkgreen: List[float] = field(default_factory=lambda: [0, 0.58823, 0.19607])
    graylevel: float = 0.6
    
    # Layout settings
    popboxspaceadj: float = 1.0
    moption: Union[str, float] = 's'
    skipdemographicscaling: bool = False
    anglenames: bool = False
    dashinterval: int = 3
    
    # Names and file handling
    usealtnames: bool = False
    imagefileextension: str = ""
    altnamefilename: str = ""
    altpopnames: List[str] = field(default_factory=list)
    widthscalar: int = -1
    
    # Ghost population
    excludeghost: bool = False
    useghost: bool = False
    
    # IMa version
    imaversion: int = 3
    newercode: bool = False
    
    # Fixed dimensions
    fixedLL: List[float] = field(default_factory=lambda: [36.1, 36.1])
    fixedUR: List[float] = field(default_factory=lambda: [756.1, 576.1])
    
    # Output file handle - not initialized in the constructor
    epsf: Optional[TextIO] = None
    
    def adjust_scales(self):
        """
        Set the scale values based on the current configuration.
        This replaces the setbasexyscale function.
        """
        # Set fixed dimensions
        minimumpoint = 36.1  # Minimum x and y points
        self.fixedLL = [minimumpoint, minimumpoint]
        
        # Adjust x scale if needed
        if self.xadjust is not None:
            if self.xadjust > 1.0:
                self.maximumxpoint = 756.1 * self.xadjust
                self.localxscale = -1
            else:
                self.localxscale = self.xadjust
                self.maximumxpoint = 756.1
        
        # Adjust for square format
        if self.dosquare:
            self.maximumxpoint = 576.1
        
        # Adjust for height scale
        if self.heightscale is not None:
            self.maximumypoint = 576.1 * self.heightscale
        
        # Set the upper right point
        self.fixedUR = [self.maximumxpoint, self.maximumypoint]
    
    def update_from_args(self, args):
        """
        Update configuration from parsed command line arguments.
        This replaces the global update_globals_from_args function.
        """
        # Handle special cases before general attribute transfer
        self.fontfixed = False
        if hasattr(args, 'fontsize') and args.fontsize is not None:
            self.fontfixed = True
        
        self.usealtnames = False
        if hasattr(args, 'altnamefilename') and args.altnamefilename is not None:
            self.usealtnames = True
            self.altnamefilename = args.altnamefilename
        
        self.set_lastt_lower_y = True
        if hasattr(args, 'lastt_lower_y') and args.lastt_lower_y is not None:
            self.set_lastt_lower_y = False
            self.lastt_lower_y = args.lastt_lower_y
        
        self.popboxcintervalboxes = True
        if hasattr(args, 'no_popboxcintervalboxes') and args.no_popboxcintervalboxes:
            self.popboxcintervalboxes = False
        
        self.popboxcintervalarrows = True
        if hasattr(args, 'no_popboxcintervalarrows') and args.no_popboxcintervalarrows:
            self.popboxcintervalarrows = False
        
        # Process image format
        self.imagefileextension = ""
        if hasattr(args, 'imageformat') and args.imageformat is not None:
            if args.imageformat == 'j':
                self.imagefileextension = ".jpg"
            elif args.imageformat == 'p':
                self.imagefileextension = ".pdf"
            elif args.imageformat == 'n':
                self.imagefileextension = ".png"
        
        # Transfer general attributes 
        for key, value in vars(args).items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        
        # Post-processing
        if self.anglenames:
            self.line0y = 0.88
        
        # Set bifont based on font
        self.bifont = self.font + "-BoldItalic"
        
        # Convert moption to float if it's a number
        if self.moption not in ['x', 'a', 's', 'S']:
            try:
                self.moption = float(self.moption)
            except ValueError:
                self.moption = 's'  # Default if conversion fails
        
        # Ensure output has .eps extension
        if not self.outputfilename.lower().endswith('.eps'):
            self.outputfilename += '.eps'

##***********************************************************************************
##////////////// FUNCTIONS FOR GENERATING EPS FILE   ////////////////////////////////
##***********************************************************************************

## the coordinates of the figure assume that the origin is at the lower left of the figure
## increasing x values move to the right of the origin
## increasing y values more up from the origin
## for a point, given as a list of two values, the x value is in position 0, the y value in position 1
## a box is defined as two points, the lower left and the upper right
## almost all of the code assumes a coordinate system of 0 to 1 on both the x and y axes
## all conversions to the absolute scale are handled by going thru apoint()


def w(args, s):
    """
    Simple function to make it easier to write to the eps file without having to repeat code text
    """
    args.epsf.write(s + "\n")


def apoint(args, rpoint):
    """
    rpoint is a list of length 2, convert a relative point to an absolute point

    x value is in position 0, y value in position 1

    this function handles all conversions from relative to absolute scales
    all the rest of the code assumes coordinates from 0 to 1, on both x and y axes

    some other configuration parameters that are used here:
        fixedLL is the lower left point of the plot - no values to left or below this.
        fixedUR is the upper right point - not values above or to the right of this.
        globalscale is an overall scalar of plot size
        localxscale is an x dimensional scalar of plot size
            so the x dimension can be changed without affecting the y dimension
    """
    tempy = args.fixedLL[1] + args.globalscale * rpoint[1] * (args.fixedUR[1] - args.fixedLL[1])
    if args.localxscale != -1:
        tempx = args.fixedLL[0] + args.localxscale * args.globalscale * rpoint[0] * (args.fixedUR[0] - args.fixedLL[0])
    else:
        tempx = args.fixedLL[0] + args.globalscale * rpoint[0] * (args.fixedUR[0] - args.fixedLL[0])
    if tempx - args.fixedUR[0] > 0 and tempx - args.fixedUR[0] < 1e-7:
        tempx = args.fixedUR[0]
    if tempx > args.fixedUR[0]:
        print("problem x value : ", tempx, " max x allowed : ", args.fixedUR[0])
    return [tempx, tempy]


def rapoint(args, rpoint):
    """ relative point
        this is called from a function where the scale has been reset
    """
    return [rpoint[0] * args.globalscale * (args.fixedUR[0] - args.fixedLL[0]),
            rpoint[1] * args.globalscale * (args.fixedUR[1] - args.fixedLL[1])]


def textwide(args, s, tf):
    """ approx width of text """
    width = 350  # default ok for Arial or Helvetica
    if args.font == "Times-roman":
        width = 330
    if args.font == "Courier":
        width = 390
    if args.fontfixed is False:
        localfontsize = int(args.fontsize * args.globalscale)
    else:
        localfontsize = int(args.fontsize)
    return tf * localfontsize * len(s) * width / (1000 * (args.fixedUR[0] - args.fixedLL[0]))


def dotext(args, rpoint, text, angle, bi):
    """
    Print text beginning at rpoint at angle
    bi is a boolean, True indicates bold italic font
    font and bifont are global
    """
    if bi:
        w(args, "/%s findfont" % args.bifont)
    else:
        w(args, "/%s findfont" % args.font)
    if args.fontfixed is False:
        localfontsize = int(args.fontsize * args.globalscale)
    else:
        localfontsize = int(args.fontsize)
    w(args, "%d scalefont" % localfontsize)
    w(args, "setfont")
    w(args, "newpath")
    p = apoint(args, rpoint)
    if angle != 0:
        w(args, "gsave")
        w(args, "%d %d translate" % (p[0], p[1]))
        w(args, "%d  rotate" % angle)
        w(args, "0  0 moveto")
        w(args, "(" + text + ") show")
        w(args, "grestore")
    else:
        w(args, "%d %d moveto" % (p[0], p[1]))
        w(args, "(" + text + ") show")


def curvecontrol(p1, p2, u_or_d):
    """ returns two control points to draw a curve between two points
        that are the corners of a box
        u_or_d is 1 to draw the curve above the line between the two point
        or 0 to draw the curve below the line between the two points
    """
##    four possibile orders:
##      A  p1 lower and to left of p2
##      B  p1 lower and to right of p2
##      C  p1 higher and to left of p2
##      D  p1 higher and to right of p2
##    B and C are reverse of each other
##    A and D are reverse of each other
##    so only 2 types of pairs really
##    each has a curve up or curve down possibility
##    start by converting D to A, and C to B
    e1 = 0.0001
    e2 = 0.9
    e1c = 1 - e1
    e2c = 0.5
    cp1 = []
    cp2 = []
    if p2[1] < p1[1]:
        resort = True
        ptemp = p2
        p2 = p1
        p1 = ptemp
    else:
        resort = False
    if p1[0] < p2[0]:   # type A
        if u_or_d:   # curve up
            cp1.append(((p2[0] - p1[0]) * e1) + p1[0])
            cp1.append(((p2[1] - p1[1]) * e2) + p1[1])
            cp2.append(((p2[0] - p1[0]) * e2c) + p1[0])
            cp2.append(((p2[1] - p1[1]) * e1c) + p1[1])
        else:
            cp1.append(((p2[0] - p1[0]) * e2) + p1[0])
            cp1.append(((p2[1] - p1[1]) * e1) + p1[1])
            cp2.append(((p2[0] - p1[0]) * e1c) + p1[0])
            cp2.append(((p2[1] - p1[1]) * e2c) + p1[1])
    else:  # type B
        if u_or_d:   # curve up
            cp1.append(p1[0] - ((p1[0] - p2[0]) * e1))
            cp1.append(((p2[1] - p1[1]) * e2) + p1[1])
            cp2.append(p1[0] - ((p1[0] - p2[0]) * e2c))
            cp2.append(((p2[1] - p1[1]) * e1c) + p1[1])
        else:
            cp1.append(p1[0] - ((p1[0] - p2[0]) * e2))
            cp1.append(((p2[1] - p1[1]) * e1) + p1[1])
            cp2.append(p1[0] - ((p1[0] - p2[0]) * e1c))
            cp2.append(((p2[1] - p1[1]) * e2c) + p1[1])
    if resort:
        ptemp = cp2
        cp2 = cp1
        cp1 = ptemp
    return cp1, cp2


def calccdim(args, cdimval, cbox):
    ll = cbox[0]
    ur = cbox[1]
    curvesizedefine = 0.02
    if cdimval == -1:
        cdimval = (args.fixedUR[0] - args.fixedLL[0]) * curvesizedefine
    lla = apoint(args, ll)
    ura = apoint(args, ur)
    if ura[0] - lla[0] < 2 * cdimval or ura[1] - lla[1] < 2 * cdimval:
        cdimval = min(ura[0] - lla[0], ura[1] - lla[1]) / 2.0
    return cdimval


def curvebox(args, cdim, cbox, width, color, grayamount, popnum, dash, poptree):
    """
        creates a box with curved corners, size of the curve set by curvesize
        if dash==0 and rgbcolor is True, fills the box with a lighter version of
        the rgbcolor for the population the box is for
        returns cdim which has something to do with the size of the box
    """
    if dash > 0:
        w(args, "[%d %d] 0 setdash" % (dash, dash))
    ll = cbox[0]
    ur = cbox[1]
    curvesizedefine = 0.02
    if cdim == -1:
        cdim = (args.fixedUR[0] - args.fixedLL[0]) * curvesizedefine
    lla = apoint(args, ll)
    ura = apoint(args, ur)
    if ura[0] - lla[0] < 2 * cdim or ura[1] - lla[1] < 2 * cdim:
        cdim = min(ura[0] - lla[0], ura[1] - lla[1]) / 2.0
    ula = [lla[0], ura[1]]
    lra = [ura[0], lla[1]]
    if args.rgbcolor:
        boxcolorstring = ("%f %f %f setrgbcolor" %
                          (poptree[popnum][5][0], poptree[popnum][5][1], poptree[popnum][5][2]))
        lightcolor = []
        for ii in range(3):
            # lighten to 10% of color
            lightcolor.append(1.0 - (0.1 * (1.0 - poptree[popnum][5][ii])))
        lightboxfillcolorstring = "%f %f %f setrgbcolor" % (lightcolor[0], lightcolor[1], lightcolor[2])
    else:
        if color != args.black:
            color = args.blue
            gcolor = []
            for i in range(3):
                if color[i] == 0:
                    gcolor.append(grayamount)
                else:
                    gcolor.append(color[i])
            boxcolorstring = "%f %f %f setrgbcolor" % (gcolor[0], gcolor[1], gcolor[2])
        else:
            boxcolorstring = "%f setgray" % grayamount

    w(args, "newpath")
    w(args, "%d  %d  moveto" % (lla[0] + cdim, lla[1]))
    cp1 = [lra[0] - cdim, lra[1]]
    cp2 = [lra[0], lra[1] + cdim]
    w(args, "%d  %d  lineto" % (cp1[0], cp1[1]))
    ccpoints = curvecontrol(cp1, cp2, 0)
    w(args, "%d %d %d %d %d  %d  curveto" % (ccpoints[0][0], ccpoints[0][1], ccpoints[1][0], ccpoints[1][1], cp2[0], cp2[1]))
    cp1 = [ura[0], ura[1] - cdim]
    cp2 = [ura[0] - cdim, ura[1]]
    w(args, "%d  %d  lineto" % (cp1[0], cp1[1]))
    ccpoints = curvecontrol(cp1, cp2, 1)
    w(args, "%d %d %d %d %d  %d  curveto" % (ccpoints[0][0], ccpoints[0][1], ccpoints[1][0], ccpoints[1][1], cp2[0], cp2[1]))
    cp1 = [ula[0] + cdim, ula[1]]
    cp2 = [ula[0], ula[1] - cdim]
    w(args, "%d  %d  lineto" % (cp1[0], cp1[1]))
    ccpoints = curvecontrol(cp1, cp2, 1)
    w(args, "%d %d %d %d %d  %d  curveto" % (ccpoints[0][0], ccpoints[0][1], ccpoints[1][0], ccpoints[1][1], cp2[0], cp2[1]))
    cp1 = [lla[0], lla[1] + cdim]
    cp2 = [lla[0] + cdim, lla[1]]
    w(args, "%d  %d  lineto" % (cp1[0], cp1[1]))
    ccpoints = curvecontrol(cp1, cp2, 0)
    w(args, "%d %d %d %d %d  %d  curveto" % (ccpoints[0][0], ccpoints[0][1], ccpoints[1][0], ccpoints[1][1], cp2[0], cp2[1]))
    w(args, "closepath")
    w(args, "gsave")
    if args.rgbcolor and dash == 0:  # fill the box with a lighter version of the color used for the lines of the box
        w(args, lightboxfillcolorstring)
        w(args, "fill")
        w(args, "grestore")
    width = float(width)
    w(args, "%f setlinewidth" % (width * args.globalscale))
    w(args, boxcolorstring)
    w(args, "stroke")
    if args.simplecolor or args.rgbcolor:
        w(args, "0 0 0  setrgbcolor")
    else:
        w(args, "0 setgray")
    if dash > 0:
        w(args, "[] 0 setdash")
    return cdim


def aline(args, p, width, dash, grayamount):
    """ p is a list of points in relative space (0-1)
        dash is the spacing (in point scale) of dashes in the line
        if dash is zero there is no dashing """
    if grayamount > 0:
        w(args, "%f setgray" % grayamount)
    ap = []
    for i in range(len(p)):
        ap.append(apoint(args, p[i]))
    if dash > 0:
        w(args, "[%d %d] 0 setdash" % (dash, dash))

    w(args, "%d %d moveto" % (ap[0][0], ap[0][1]))
    for j in range(1, len(p)):
        w(args, "%d %d lineto" % (ap[j][0], ap[j][1]))
    width *= args.globalscale
    w(args, "%f setlinewidth" % width)
    w(args, "stroke")
    w(args, "[ ] 0 setdash")
    if grayamount > 0:
        w(args, "0 setgray")


def arrowhead(args, head, headwidth, angle):
    """ draw arrowhead width on the same scale as points in head
        head is the center of the arrowhead
        angle = 0 has the arrow pointing to the right
    """
    w(args, "%% begin arrowhead")
    holdhead = apoint(args, head)
    head = [0, 0]
    tip = rapoint(args, [head[0] + headwidth, head[1]])
    p1 = rapoint(args, [head[0] - headwidth, head[1] + headwidth])
    p2 = rapoint(args, [head[0] - headwidth, head[1] - headwidth])
    c1 = rapoint(args, [head[0], head[1] - headwidth / 2])
    c2 = rapoint(args, [head[0], head[1] + headwidth / 2])
    w(args, "gsave")
    w(args, "%d %d translate" % (holdhead[0], holdhead[1]))
    w(args, "%d  rotate" % angle)
    w(args, "%d %d moveto" % (p1[0], p1[1]))
    w(args, "%d %d lineto" % (tip[0], tip[1]))
    w(args, "%d %d lineto" % (p2[0], p2[1]))
    w(args, "%d %d %d %d %d %d curveto" % (c1[0], c1[1], c2[0], c2[1], p1[0], p1[1]))
    w(args, "closepath")
    w(args, "fill")
    w(args, "grestore")
    w(args, "%% end arrowhead")


def arrow(args, head, tail, direc, color):
    """
        draw an arrow. head and tail are points, width is on the same scale
        direc = 0 right, 1 up, 2 left, 3 down
        arrow is gray
    """
    headwidth = ARROWHEAD_WIDTH_DEFAULT * args.arrowheightadj
    if (direc == 0):
        headadj = [head[0] - headwidth, head[1]]
    if (direc == 1):
        headadj = [head[0], head[1] - headwidth]
    if (direc == 2):
        headadj = [head[0] + headwidth, head[1]]
    if (direc == 3):
        headadj = [head[0], head[1] + headwidth]
    if color != args.black:
        color = args.blue
        gcolor = []
        for i in range(3):
            if color[i] == 0:
                gcolor.append(args.graylevel)
            else:
                gcolor.append(color[i])
        w(args, "%f %f %f setrgbcolor" % (gcolor[0], gcolor[1], gcolor[2]))
    else:
        w(args, "%f setgray" % args.graylevel)
    arrowhead(args, headadj, headwidth, direc * 90)
    ahead = apoint(args, headadj)
    atail = apoint(args, tail)
    w(args, "%d %d moveto" % (ahead[0], ahead[1]))
    w(args, "%d %d lineto" % (atail[0], atail[1]))
    w(args, "%f setlinewidth" % (2 * args.globalscale))
    w(args, "stroke")
    if args.simplecolor or args.rgbcolor:
        w(args, "0 0 0  setrgbcolor")
    else:
        w(args, "0 setgray")


def migrationstraightarrow(args, val2NM, head, tail, direc, color):
    """
        draw an arrow. head and tail are points, width is on the same scale
        direc = 0 right, 1 up, 2 left, 3 down
        arrow is gray
    """
    headwidth = ARROWHEAD_WIDTH_DEFAULT * 1.5 * args.arrowheightadj
    headwidth = ARROWHEAD_WIDTH_DEFAULT * args.arrowheightadj
    if (direc == 0):
        headadj = [head[0] - headwidth, head[1]]
    if (direc == 1):
        headadj = [head[0], head[1] - headwidth]
    if (direc == 2):
        headadj = [head[0] + headwidth, head[1]]
    if (direc == 3):
        headadj = [head[0], head[1] + headwidth]
    if args.simplecolor or args.rgbcolor:
        w(args, "%f %f %f setrgbcolor" % (color[0], color[1], color[2]))
    else:
        w(args, "0 0 0 setrgbcolor")
    arrowhead(args, headadj, headwidth, direc * 90)
    ahead = apoint(args, headadj)
    atail = apoint(args, tail)
    w(args, "%d %d moveto" % (ahead[0], ahead[1]))
    w(args, "%d %d lineto" % (atail[0], atail[1]))
    w(args, "%f setlinewidth" % (2 * args.globalscale))
    w(args, "stroke")
    if args.simplecolor or args.rgbcolor:
        w(args, "0 0 0  setrgbcolor")
    text2NMwidth = textwide(args, val2NM, 1.5)
    if (direc == 0):  # arrow to the right, line is shifted up, text is below line
        if text2NMwidth > abs(tail[0] - headadj[0]):
            textpoint = tail
        else:
            textpoint = [(headadj[0] + tail[0]) / 2, tail[1]]
    if (direc == 2):
        if text2NMwidth > abs(tail[0] - headadj[0]):
            textpoint = headadj
        else:
            textpoint = [(headadj[0] + tail[0]) / 2, headadj[1]]
    dotext(args, textpoint, val2NM, 0, True)


## stopped using this 7/9/2018  arrows were taking up too much space and causing headaches over location
## switched to migrationstraightarrow()
# def migrationcurvearrow(args, val2NM, head, tail, direc, color):
#     """ direct can be 0 or 2 (right or left)  if 0 curveheight is positive and curve goes up from
#         the tail and then down to the head
#         if direc is 2  then curve is interpreted to be negative and curve goes down from the tail
#         and then up to the head """
#     w(args, "%% BEGIN MIGRATION ARROW: %s" % val2NM)
#     curveheight = CURVEHEIGHT_DEFAULT
#     c2height = ARROWHEAD_WIDTH_DEFAULT
#     headwidth = c2height * 1.5 * args.arrowheightadj
#     width = 3.5
#     if (direc == 0):  # arrow to the right, line is shifted up, text is below line
#         textpoint = [tail[0], tail[1] - curveheight]
#         cheadadj = [head[0] - headwidth, head[1] + c2height]
#         ctail = [tail[0], tail[1] + c2height]
#         arrowheadpoint = [cheadadj[0], head[1] + c2height / 1.2]
#         if args.simplecolor or args.rgbcolor:
#             w(args, "%f %f %f setrgbcolor" % (color[0], color[1], color[2]))
#         arrowhead(args, arrowheadpoint, headwidth, 330)        # head tilted down to the right
#         if args.simplecolor or args.rgbcolor:
#             w(args, "0 0 0 setrgbcolor")
#         if abs(cheadadj[0] - ctail[0]) > 0:
#             curveheightmultiplier = math.pow(abs(cheadadj[0] - ctail[0]) / 0.15, 0.1)
#         else:
#             curveheightmultiplier = 1
#         cp1 = [ctail[0] + (cheadadj[0] - ctail[0]) * 0.8, cheadadj[1] + curveheight * curveheightmultiplier]
#         cp2 = [ctail[0] + (cheadadj[0] - ctail[0]) * 0.2, cheadadj[1] + curveheight * curveheightmultiplier]
#         textpoint = [cp2[0], cheadadj[1] - curveheight / 3]
#     if (direc == 2):   # arrow to the left, line is shifted down, text is above line
#         cheadadj = [head[0] + headwidth, head[1]]
#         textpoint = [cheadadj[0] + c2height, cheadadj[1]]
#         ctail = tail
#         arrowheadpoint = [cheadadj[0], cheadadj[1] + c2height / 3.5]
#         if args.simplecolor or args.rgbcolor:
#             w(args, "%f %f %f setrgbcolor" % (color[0], color[1], color[2]))
#         arrowhead(args, arrowheadpoint, headwidth, 150)       # head tilted up to the left
#         if args.simplecolor or args.rgbcolor:
#             w(args, "0 0 0 setrgbcolor")
#         if abs(cheadadj[0] - ctail[0]) > 0:
#             curveheightmultiplier = math.pow(abs(cheadadj[0] - ctail[0]) / 0.15, 0.1)
#         else:
#             curveheightmultiplier = 1

#         cp1 = [cheadadj[0] + (ctail[0] - cheadadj[0]) * 0.2, cheadadj[1] - curveheight * curveheightmultiplier]
#         cp2 = [cheadadj[0] + (ctail[0] - cheadadj[0]) * 0.8, cheadadj[1] - curveheight * curveheightmultiplier]
#         textpoint = [cp1[0], cheadadj[1] - curveheight / 3]

#     ahead = apoint(args, cheadadj)
#     atail = apoint(args, ctail)
#     acp1 = apoint(args, cp1)
#     acp2 = apoint(args, cp2)
#     if width > 0:
#         if args.simplecolor or args.rgbcolor:
#             w(args, "%f %f %f setrgbcolor" % (color[0], color[1], color[2]))
#         w(args, "%f setlinewidth" % (width * args.globalscale))
#         w(args, "%d %d moveto" % (ahead[0], ahead[1]))
#         w(args, "%d %d  %d  %d  %d  %d curveto" % (acp1[0], acp1[1], acp2[0], acp2[1], atail[0], atail[1]))
#         w(args, "stroke")
#         if args.simplecolor or args.rgbcolor:
#             w(args, "0 0 0 setrgbcolor")
#         dotext(args, textpoint, val2NM, 0, True)
#         if args.simplecolor or args.rgbcolor:
#             w(args, "0 0 0 setrgbcolor")
#     w(args, "%% END MIGRATION ARROW")



##***********************************************************************************
##////////////// FUNCTIONS FOR GETTING VALUES OUT OF INPUT FILE        //////////////
##***********************************************************************************
##    These functions put information in slist, a 2D global list of lists array
##    The main function here readimfile() which builds slist[][]
##    Each list in slist contains the details regarding a particular type of info to be
##    obtained from the input file
##    slist[i][0] - a brief text description about the category of information
##    slist[i][1] - a boolean value that is initialized as False, but becomes True after the info is obtained
##    slist[i][2] - the name of the function that reads the information of that type
##    slist[i][3] - the string used to search the input file, when it is found the function is called
##    slist[i][ > 3 ] - the actual information, the types and number of values vary depending on the category of information
##      all of the functions (names in slist[i][2] are called with
##            slist[i][2](args,imfile, imfileline, slist[i][3], numpops)
##        what the function returns is appended to slist[i]


def get_input_file_name(args,f, a, s):
    # f not used but needed for function to match general function format
    return a[len(s):len(a)].strip()


def check_ghost_status(args,f, a, s):
    # f, s not used but needed for function to match general function format
    if args.newercode:
        args.useghost = (a.find("-j") >= 0) # don't think this is needed, 1 is the default ' and ("1" in a[a.find("-j") + 1:])   # should only be true if -j is there with a 1'
    else:
        args.useghost = (a.find("-j") >= 0) and ("4" in a[a.find("-j") + 1:])  # should only be true if -j is there with a 4


def get_population_names(args,f, a, s):
    # a, s not used but needed for function to match general function format
    """
        usealtnames and altnamefilename defined previously
    """

    aa = f.readline().strip()
    popnamelist = []
    i = 0
    foundghost = False
    while aa.find("Population") >= 0:
        popname = aa.split()[3]
        if popname.upper() == "GHOST":
            foundghost = True
        popnamelist.append(popname)
        i += 1
        aa = f.readline().strip()
    if args.useghost is True and foundghost is False:   # for compatibility with older output files
        popnamelist.append('ghost')
    anames = []
    if args.usealtnames:
        for line in open(args.altnamefilename, "r"):
            temp = line.strip()
            if len(temp) > 0:
                anames.append(temp)
        anames = anames[0:len(popnamelist)]
    args.altpopnames = list(anames)
    return popnamelist


def get_population_tree(args,f, a, s):
    # s variables is not used but needed for function to match general function format
    """
         a couple possible things to read here
    """
    tempstring = a.split()[-1]
    while a.find("Population Tree") >= 0:
        if a.find("standard ordering") >= 0:
            tempstring = a.split()[-1]
        if a.find("Ghost Population") >= 0 and args.excludeghost is False:
            assert args.useghost
            tempstring = a.split()[-1]
        a = f.readline().strip()
    return tempstring


def get_popsize_param(args,f, a, s):
    # a, s not used but needed for function to match general function format
    """ read the histogram table of marginal distributions for population sizes:
    For each population it reads:
    the label of the parameter
    the HiSmth value
    the HPD95Lo value
    the HPD95Hi value """
    psp = []
    for i in range(4):
        aa = f.readline().split()
    for i in range(2 * numpops - 1):
        psp.append([])
        psp[i].append(aa[i + 1])
    aa = f.readline().split()
    while len(aa) > 0:
        while aa.count("#"):
            aa.remove("#")
        while aa.count("#?"):
            aa.remove("#?")
        while aa.count("?"):
            aa.remove("?")
        found = False
        if aa[0] == "HiPt" or aa[0] == "HPD95Lo" or aa[0] == "HPD95Hi":
            found = True
        if found:
            for i in range(2 * numpops - 1):
                psp[i].append(float(aa[i + 1].strip('?#')))
        aa = f.readline().split()
    return psp


def get_t_param(args,f, a, s):
    # a, s not used but needed for function to match general function format
    """ read the table of marginal distributions for splitting times:
        For each splittingtime it reads:
            the label of the parameter
            the HiSmth value
            the HPD95Lo value
            the HPD95Hi value """
    psp = []
    aa = f.readline().split()
    while aa[0] != "Value":
        aa = f.readline().split()
    for i in range(numpops - 1):
        psp.append([])
        psp[i].append(aa[i + 1])
    aa = f.readline().split()
    while len(aa) > 0:
        while aa.count("?"):
            aa.remove("?")
        while aa.count("#"):
            aa.remove("#")
        found = False
        if aa[0] == "HiSmth" or aa[0] == "HPD95Lo" or aa[0] == "HPD95Hi":
            found = True
        if found:
            for i in range(numpops - 1):
                psp[i].append(float(aa[i + 1].strip('?#')))
        aa = f.readline().split()
    return psp


def msigvals(args,ss):
    """
        get the significance levels for the migration(m) rates
        make a list, each element is a list
            the migration parameter
            the significance level
    """
    si = 0
    msiglist = []
    nummp = 0
    while si < len(ss):
        if ss[si].find("Migration Rate Parameters") == 0:
            si += 1
            aa = ss[si].split()
            ainc = 1 if args.newercode else 2
            for i in range(1, len(aa), ainc):
                msiglist.append([aa[i]])
            si += 3
            aa = ss[si].split()[1:]
            ii = 0
            for i in range(0, len(aa), ainc):
                temp = aa[i]
                try:
                    mtemp = float(temp)
                except ValueError:
                    mtemp = 0.0
                msiglist[nummp + ii].append(mtemp)
                ii += 1
            si += 4 if args.newercode else 3
            aa = ss[si].split()[1:]
            for i, temp in enumerate(aa):
                if (args.moption == 's' and temp.count('*') == 0) or (args.moption == 'S' and temp.count('*') <= 1):
                    msiglist[nummp + i].append('ns')
                else:
                    msiglist[nummp + i].append('*' * temp.count('*'))
            nummp = len(msiglist)
        si += 1
    return msiglist


def get_2NM(args,f, a, s):
    """
        replaces older code as of 9/18/2017
        reads a chunk of the file from "Marginal Peak Locations and Probabilities"  to "HISTOGRAMS"
        then scans this for the significance values for the migration (m) rates
        then gets the 2Nm numbers
        and then matches the m values with the 2Nm numbers
        returns a list, each element has 3 items:
            the name, e.g. 2N0m0>1
            the estimated 2NM value
            a string with the significance level
    """
    ss = []
    while True:
        ss.append(f.readline().strip())
        if ss[-1].upper().find("HISTOGRAMS") == 0:
            break
    msiglist = msigvals(args,ss)
    si = 0
    psp = []
    nummp = 0
    ainc = 1 if args.newercode else 2
    while si < len(ss):
        if ss[si].upper().find("POPULATION MIGRATION (2NM) TERMS") == 0:
            si += 1
            aa = ss[si].split()
            for i in range(1, len(aa), ainc):
                psp.append([aa[i]])
            si += 3
            aa = ss[si].split()
            ii = 0
            for i in range(1, len(aa), ainc):
                psp[nummp + ii].append(float(aa[i]))
                ii += 1
            nummp = len(psp)
        si += 1

    for pi, p in enumerate(psp):
        mn = p[0][p[0].upper().find('M') + 1:]
        i = 0
        while True:
            assert i < len(msiglist)
            if msiglist[i][0][1:] == mn:
                psp[pi].append(msiglist[i][2])
                psp[pi].append(msiglist[i][1])
                break
            i += 1
    return psp


def get_demog_scales(args,f, a, s):
    psp = [0, 0, 0]
    for i in range(10):  # go down several lines and look for the necessary information, very crude and
        aa = f.readline().split()
        if aa[0] == "Generation" and aa[1] == "time":
            psp[0] = float(aa[len(aa) - 1])
        if aa[0] == "Geometric" and aa[3] == "mutation":
            psp[1] = float(aa[len(aa) - 1])
        if aa[0] == "Geometric" and aa[3] == "ML":
            psp[2] = float(aa[len(aa) - 1])
    return psp


def get_parameter_priors(f, a, s):
    psp = [["population size", "uniform"], ["migration", "uniform"], ["splittime", "uniform"]]
    aa = f.readline()
    for i in range(3):
        aa = f.readline().split()
        psp[i].append(float(aa[len(aa) - 1]))
        if aa.count("exponential") > 0:
            psp[i][1] = "exponential"
    return psp

def calc_scaledvals(slist):
    gentime = slist[7][4][0]
    timeumean = slist[7][4][1]
    scaleumean = slist[7][4][2]
    scaledpop = []
    for i in range(2 * numpops - 1):
        scaledpop.append(slist[4][4][i][1] / (4.0 * timeumean * gentime))
    scaledtime = []
    for i in range(numpops - 1):
        scaledtime.append(slist[5][4][i][1] * (scaleumean / timeumean))
    return scaledpop, scaledtime


def mysplit(s, delims):
    """
        a simple function to split a string s at the occurence of all instances of each symbol in delims
        works simply by replacing every instance of each of those symbols  with a space
        and then calling regular split
        does not work with more than one symbol at a time
    """
    for c in delims:
        s = s.replace(c, ' ')
    return s.split()


def removeghost(slist, scaledpop, scaledtime):
    # remove 'ghost' from names
    global numpops
    slist[2][4] = slist[2][4][:-1]
    npops = len(slist[2][4])
    slist[4][4].pop(npops)  # remove ghost
    slist[4][4].pop(-1)  # remove last ancestor
    if scaledpop != []:   # do same for scaled pops
        scaledpop.pop(npops)
        scaledpop.pop(-1)

    # renumber ancestors
    for i in range(npops, 2 * npops - 1):
        assert slist[4][4][i][0] == 'q' + str(i + 1)
        slist[4][4][i][0] = 'q' + str(i)
    slist[5][4].pop(-1)  # remove last splitting time parameter
    if scaledtime != []:
        scaledtime.pop(-1)  # do the same for scaledtime
    #remove migration parameters associated with ghost
    ghoststr = str(npops)
    numpops = npops
    assert len(ghoststr) == 1
    rlist = []
    for i, m in enumerate(slist[6][4]):
        if m[0][2] == ghoststr or m[0][-1] == ghoststr:
            rlist.append(i)
    rlist.reverse()
    for r in rlist:
        slist[6][4].pop(r)
    #renumber ancestors in migration parameters
    for i, minfo in enumerate(slist[6][4]):
        j = 1
        m = minfo[0]
        newm = m[0]  # get the 2 in 2Nm
        while True:
            if m[j] not in string.digits:
                newm += m[j]
                j += 1
            else:
                d = m[j]
                j += 1
                if j < len(m) and m[j] in string.digits:
                    d += m[j]
                    j += 1
                if int(d) >= int(ghoststr):
                    newm += str(int(d) - 1)
                else:
                    newm += str(int(d))
            if j >= len(m):
                break
        slist[6][4][i][0] = newm
    return slist, scaledpop, scaledtime


def readimfile(args):
    """
        gets info from the input file
        returns all information in slist
        and scale information in scaledpop and scaledtime
        very kludgy, all sorts of awkward things to deal with IMa2 vs IMa3 differences and to catch wrong kinds of input files
    """
    global numpops
    imfile = open(Path(args.imfilename), "r")
    args.useghost = False
    imfileline = imfile.readline()
    if "IMa3" == imfileline[0:4]:
        args.imaversion = 3
    else:
        args.imaversion = 2
    args.newercode = False
    while imfileline != '':
        if imfileline.upper().find("IMa3 program compiled on".upper()) >= 0:
            from datetime import datetime
            linesplit = mysplit(imfileline.strip(), ", ")
            date = " ".join(linesplit[4:7])
            # changed ima3 format afer sep 12 2017
            newcodedatetime = datetime.strptime("sep 12 2017", '%b %d %Y')
            filedatetime = datetime.strptime(date, '%b %d %Y')
            args.newercode = filedatetime >= newcodedatetime
        if imfileline.upper().find("Command line string :".upper()) >= 0:
            break
        imfileline = imfile.readline()
    slist = [["ghost status", False, check_ghost_status, "Model options on command line"],
             ["inputfile", False, get_input_file_name, "Text from input file:"],
             ["pop names", False, get_population_names, "Population Names"],
             ["pop tree", False, get_population_tree, "Population Tree :"],
             ["population size parameter info", False, get_popsize_param, "MARGINAL DISTRIBUTION VALUES AND HISTOGRAMS OF POPULATION SIZE AND MIGRATION PARAMETERS"],
             ["splitting time parameter info", False, get_t_param, "MARGINAL DISTRIBUTION VALUES AND HISTOGRAMS OF PARAMETERS IN MCMC"],
             ["migration parameter info", False, get_2NM, "Marginal Peak Locations and Probabilities"],
             ["demographic scales", args.skipdemographicscaling, get_demog_scales, "MARGINAL DISTRIBUTION VALUES IN DEMOGRAPHIC UNITS"]  #,
##             ["parameter priors", False, get_parameter_priors, "Parameter Priors"] \  ignore this I think
             ]
    while imfileline != '':
        if imfileline.upper().find("LOCATIONS OF PARAMETER ESTIMATES IN THIS FILE".upper()) >= 0:
            while True:
                imfileline = imfile.readline()
                if imfileline.upper().find("Hyperparameter".upper()) >= 0:
                    print("**IMfig error - input while was run using hyperparameters")
                    quit()
                if imfileline.upper().find("ESTIMATED POSTERIOR PROBABILITIES OF POPULATION TREE TOPOLOGIES".upper()) >= 0:
                    print("**IMfig error - input while was generated to estimate phyhlogeny")
                    quit()
                if imfileline.upper().find("INPUT AND STARTING INFORMATION".upper()) >= 0:
                    break
        checkdone = True
        for i in range(len(slist)):
            checkdone = checkdone and slist[i][1]
            if slist[i][1] is False and imfileline.upper().find(slist[i][3].upper()) >= 0:
                if slist[i][0] == "ghost status":
                    slist[i][2](args,imfile, imfileline, slist[i][3])
                else:
                    slist[i].append(slist[i][2](args,imfile, imfileline, slist[i][3]))
                slist[i][1] = True
                if slist[i][0] == "pop names":
                    numpops = len(slist[i][4])
        if checkdone:
            break
        imfileline = imfile.readline()
        if "**NO DATA **" in imfileline:
            print("**IMfig error - input while was run without data")
            quit()
    imfile.close()
    (scaledpop, scaledtime) = ([], [])
    if args.skipdemographicscaling:
        slist[7][1] = False
    else:
        if len(slist[7]) == 4:
            print("**IMfig error - Information in demographic units not found, use -d option")
##            printcommandset()
            quit()
        if len(slist[7][4]) == 3:
            (scaledpop, scaledtime) = calc_scaledvals(slist)
    if args.excludeghost and args.useghost:
        slist, scaledpop, scaledtime = removeghost(slist, scaledpop, scaledtime)
    return args,slist, scaledpop, scaledtime



##***********************************************************************************
##////////////// FUNCTIONS FOR READING THE POPULATION TREE STRING  //////////////////
##***********************************************************************************


def parenth(tempcurrent, poptree, poptreestring, stringspot, ancestralpopnums, rootpop, nextnode, periodi):
    current = ancestralpopnums[tempcurrent]
    stringspot += 1
    while poptreestring[stringspot].isspace():
        stringspot += 1
    while True:
        if poptreestring[stringspot].isdigit():
            if stringspot <= len(poptreestring) - 2 and poptreestring[stringspot + 1].isdigit():
                ts = poptreestring[stringspot] + poptreestring[stringspot + 1]
                itemp = int(ts)
            else:
                itemp = int(poptreestring[stringspot])
            stringspot += 1
            if poptree[current][2] == -1:
                poptree[current][2] = itemp
            else:
                poptree[current][3] = itemp
            poptree[itemp][4] = current
        if poptreestring[stringspot] == ',':
            stringspot += 1
        if poptreestring[stringspot] == '(':
            if nextnode == -1:
                nextnode = numpops + 1
            else:
                nextnode += 1
            poptree[ancestralpopnums[nextnode]][4] = current
            if poptree[current][2] == -1:
                poptree[current][2] = ancestralpopnums[nextnode]
            else:
                poptree[current][3] = ancestralpopnums[nextnode]
            (poptree, rootpop, stringspot, periodi, nextnode) = parenth(nextnode, poptree, poptreestring, stringspot, ancestralpopnums, rootpop, nextnode, periodi)
        if poptreestring[stringspot] == ')':
            break
    stringspot += 1
    if poptreestring[stringspot] == ':':
        stringspot += 1
        if stringspot <= len(poptreestring) - 2 and poptreestring[stringspot + 1].isdigit():
            ts = poptreestring[stringspot] + poptreestring[stringspot + 1]
            i = int(ts)
        else:
            i = int(poptreestring[stringspot])
        if i < numpops:
            print(" wrong number of ancestral populations indicated. string %c " % poptreestring[stringspot])
        periodi = i - numpops
        poptree[current][0] = periodi + 1
        poptree[poptree[current][2]][1] = periodi + 1
        poptree[poptree[current][3]][1] = periodi + 1
        if i >= 10:
            stringspot += 2
        else:
            stringspot += 1
    else:
        poptree[current][0] = periodi + 1
        poptree[poptree[current][2]][1] = periodi + 1
        poptree[poptree[current][3]][1] = periodi + 1
        periodi += 1
    if poptree[current][4] != -1:
        current = poptree[current][4]
    else:
        periodi += 1
        poptree[current][1] = - 1
        rootpop = current

    return poptree, rootpop, stringspot, periodi, nextnode


def parenth0(current, poptree, poptreestring, stringspot, ancestralpopnums):
    nextlistspot = 0
    ne = stringspot
    popennum = 0
    psetlist = []
    for i in range(current):
        psetlist.append(-1)
    while ne < len(poptreestring):
        if poptreestring[ne] == '(':
            psetlist[nextlistspot] = popennum
            nextlistspot += 1
            popennum += 1
            ne += 1
        else:
            if poptreestring[ne] == ')':
                ne += 2
                if ne <= len(poptreestring) - 2 and poptreestring[ne + 1].isdigit():
                    ts = poptreestring[ne] + poptreestring[ne + 1]
                    itemp = int(ts)
                else:
                    itemp = int(poptreestring[ne])
                ancestralpopnums[current + psetlist[nextlistspot - 1]] = itemp
                nextlistspot -= 1
            else:
                ne += 1
    return poptree, ancestralpopnums


def set0(strlist, pos):
    """ removes elements of a list from pos to the end, save these as a separate list """
    hold = []
    while len(strlist) > pos:
        hold.append(strlist.pop(len(strlist) - 1))
    hold.reverse()
    return strlist, hold


def strlistadd(strlist, pos, c):
    if pos > (len(strlist) - 1):
        strlist.append(c)
    else:
        strlist[pos] = c
    return strlist


def joinlist(list1, list2):
    for i in range(len(list2)):
        list1.append(list2[i])
    return list1


def rewrite(substr):
    """    rewrite() rewrites the treestring in a standard order
        swivels nodes, if both have node sequence values, the one with the lower node sequence value (periodi[]) goes on the left
        if only one has a node sequence value, it goes on the right
            when neither has a node sequence value, the one with the lowest node number go on the left
       uses simple sorting for a pair.  To handle multifurcations, must put in proper sorting
       based on code in imamp  3_9_09
       works recursively  """

    slengths = [0] * (2 * numpops - 1)
    firstint = [0] * (2 * numpops - 1)
    holdsubs = [[]] * (2 * numpops - 1)
    periodi = [0] * (2 * numpops - 1)
    pos = 1
    subpos = pos
    subcount = 0
    pcount = 0
    slengths[subcount] = 0
    while 1:
        if substr[pos] == '(':
            pcount += 1
        if substr[pos] == ')':
            pcount -= 1
        pos += 1
        slengths[subcount] += 1
        if (pcount == 0):
            if (slengths[subcount] > 1):
                pos += 1
                i = int(substr[pos])
                if pos <= len(substr) - 2 and substr[pos + 1].isdigit():
                    ts = substr[pos] + substr[pos + 1]
                    i = int(ts)
                else:
                    i = int(substr[pos])
                periodi[subcount] = i
                if (i >= 10):
                    pos += 2
                    slengths[subcount] += 3
                else:
                    pos += 1
                    slengths[subcount] += 2
            else:
                periodi[subcount] = -1
            holdsubs[subcount] = substr[subpos:pos]
            (holdsubs[subcount], hold) = set0(holdsubs[subcount], slengths[subcount])
            i = 0
            while (holdsubs[subcount][i].isdigit() is False):
                i += 1
            firstint[subcount] = int(holdsubs[subcount][i])
            subcount += 1
            slengths[subcount] = 0
            if (substr[pos] == ','):
                pos += 1
            subpos = pos
        if pos >= len(substr):
            break
    if ((periodi[0] > periodi[1] and periodi[0] >= 0 and periodi[1] >= 0) or (periodi[0] >= 0 and periodi[1] < 0)):
        substr = strlistadd(substr, 0, '(')
        j = slengths[1]
        k = 0
        i = 1
        while i <= j:
            substr = strlistadd(substr, i, holdsubs[1][k])
            k += 1
            i += 1
        subpos = 1
        if (slengths[1] > 2):
            (substr, hold) = set0(substr, i)
            substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
            substr = joinlist(substr, hold)
        substr = strlistadd(substr, i, ',')
        i += 1
        subpos = i
        j += 1 + slengths[0]
        k = 0
        while i <= j:
            substr = strlistadd(substr, i, holdsubs[0][k])
            i += 1
            k += 1
        if (slengths[0] > 2):
            (substr, hold) = set0(substr, i)
            substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
            substr = joinlist(substr, hold)
        substr = strlistadd(substr, i, ')')
    else:
        if (firstint[0] > firstint[1] and periodi[0] < 0 and periodi[1] < 0):
            substr = strlistadd(substr, 0, '(')
            j = slengths[1]
            k = 0
            i = 1
            while  i <= j:
                substr = strlistadd(substr, i, holdsubs[1][k])
                k += 1
                i += 1
            subpos = 1
            if (slengths[1] > 2):
                substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
            substr = strlistadd(substr, i, ',')
            i += 1
            subpos = i
            j += 1 + slengths[0]
            k = 0
            while i <= j:
                substr = strlistadd(substr, i, holdsubs[0][k])
                i += 1
                k += 1
            if (slengths[0] > 2):
                substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
            substr = strlistadd(substr, i, ')')
        else:
            substr = strlistadd(substr, 0, '(')
            subpos = 1
            if (slengths[0] > 2):
                (substr, hold) = set0(substr, slengths[0] + 1)
                substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
                substr = joinlist(substr, hold)
            substr = strlistadd(substr, slengths[0] + 1, ',')
            subpos = slengths[0] + 2
            if (slengths[1] > 2):
                (substr, hold) = set0(substr, slengths[0] + slengths[1] + 2)
                substr[subpos:len(substr)] = rewrite(substr[subpos:len(substr)])
                substr = joinlist(substr, hold)
            substr = strlistadd(substr, slengths[0] + slengths[1] + 2, ')')
    return substr


def plistbyperiod(poptreestring, poptree):
    """ generate a list, for each period this a list of the populations in that period,
         by their number in order from from left to right as they appear in the plot"""
    plist = [[]]
    for i in range(1, len(poptreestring)):
        if poptreestring[i - 1] != ":" and (not poptreestring[i - 1].isdigit()) and poptreestring[i].isdigit():
            plist[0].append(int(poptreestring[i]))
    droppops = [[-1, -1]]
    addpop = [-1]
    numtreepops = 2 * numpops - 1
    for pi in range(1, numpops):
        droppops.append([])
        k = 0
        for j in range(numtreepops):
            if poptree[j][1] == pi:
                droppops[pi].append(j)
                k += 1
                if k > 2:
                    print("droppop problem ")
                    break
            if poptree[j][0] == pi:
                addpop.append(j)
        tplist1 = plist[pi - 1]
        tplist2 = []
        added = False
        j = 0
        while j < len(tplist1):
            if tplist1[j] == droppops[pi][0] or tplist1[j] == droppops[pi][1]:
                if added is False:
                    added = True
                    tplist2.append(addpop[pi])
            else:
                tplist2.append(tplist1[j])
            j += 1
        plist.append(tplist2)
    return plist, droppops, addpop


def poptreeread(poptreestring):
    """ copy of the function in imamp
         use a list of lists to hld poptree
         poptree[i] is the info for population [i]
         poptree[i][0] is the period number population [i] starts in
         poptree[i][1] is the period it ends in
         poptree[i][2] is the left up pop
         poptree[i][3] is the right up pop
         poptree[i][4] is the downpop
         examples
        (poptree, rootpop, poptreestring, plist, droppops, addpop) = poptreeread("(((5, 6):12, 7):13, (4, ((3, 1):9, (2, 0):8):10):11):14", 8)
        (poptree, rootpop, poptreestring, plist) = poptreeread("(4, ((3, 1):6, (2, 0):5):7):8", 5)
         """

    if ':' not in poptreestring:    # deal with change in treestring format
        newpst = ''
        for c in poptreestring:
            if c == ')':
                newpst += '):'
            else:
                newpst += c
        poptreestring = newpst
    poptree = []
    for i in range(numpops):
        poptree.append([-1, -1, -1, -1, -1])
        poptree[i][0] = 0
    numtreepops = 2 * numpops - 1
    for i in range(numpops, numtreepops):
        poptree.append([-1, -1, -1, -1, -1])
    poptreelist = []
    for i in range(len(poptreestring)):
        poptreelist.append(poptreestring[i])
    poptreelist = rewrite(poptreelist)
    newpoptreestring = ''
    for i in range(len(poptreelist)):
        newpoptreestring += poptreelist[i]
    stringspot = 0
    ancestralpopnums = []
    for i in range(2 * numpops - 1):
        ancestralpopnums.append(0)
    (poptree, ancestralpopnums) = parenth0(numpops, poptree, newpoptreestring, stringspot, ancestralpopnums)
    (poptree, rootpop, stringspot, periodi, nextnode) = parenth(numpops, poptree, newpoptreestring, stringspot, ancestralpopnums, -1, -1, 0)
    (plist, droppops, addpop) = plistbyperiod(newpoptreestring, poptree)
    return poptree, rootpop, newpoptreestring, plist, droppops, addpop


def meanrgb(color1, color2):
    """
        generates a 'mean' color
        check_colormath is global
        if colormath is available  this will do take the average in LAB space
            see colormath documentation
        otherwise there is a crude function for averaging rgb values
    """
    if check_colormath:
        srgb1 = sRGBColor(color1[0], color1[1], color1[2])
        srgb2 = sRGBColor(color2[0], color2[1], color2[2])

        lab1 = convert_color(srgb1, LabColor)
        lab2 = convert_color(srgb2, LabColor)
        lab1tuple = SpectralColor.get_value_tuple(lab1)
        lab2tuple = SpectralColor.get_value_tuple(lab2)
        labAtuple = ((lab1tuple[0] + lab2tuple[0]) / 2.0, (lab1tuple[1] + lab2tuple[1]) / 2.0,
                (lab1tuple[2] + lab2tuple[2]) / 2.0)
        labA = LabColor(labAtuple[0], labAtuple[1], labAtuple[2])
        rgbA = convert_color(labA, sRGBColor)
        rgbAtuple = SpectralColor.get_value_tuple(rgbA)
        return list(rgbAtuple)
    else:
        acolor = [0, 0, 0]
        for j in range(3):
            # this seems to give a useful average color
            meancolor = (color1[j] + color2[j]) / 2.0
            # now lighten it a bit
            acolor[j] = (1.0 - (0.8 * (1.0 - meancolor)))
        return acolor


def addcolors(poptree):
    """
        add colors to poptree.  set ancestors to average of descendant populations
    """
    rgbset = ([[0.4, 0.4, 0.4],
                [0.650980392, 0.462745098, 0.11372549],
                [0.4, 0.650980392, 0.117647059],
                [0.121568627, 0.470588235, 0.705882353],
                [0.905882353, 0.160784314, 0.541176471],
                [0.458823529, 0.439215686, 0.701960784],
                [0.850980392, 0.37254902, 0.007843137],
                [0.105882353, 0.619607843, 0.466666667],
                [0.901960784, 0.670588235, 0.007843137],
                [0.984313725, 0.603921569, 0.6]])
    for i in range(numpops):
        poptree[i].append(rgbset[i])
    for i in range(numpops, 2 * numpops - 1):
        poptree[i].append([])
    while True:
        notdone = False
        for i in range(numpops, 2 * numpops - 1):
            if poptree[i][5] == []:
                ld = poptree[i][2]
                rd = poptree[i][3]
                if poptree[ld][5] != [] and poptree[rd][5] != []:
                    acolor = [0, 0, 0]
                    acolor = meanrgb(poptree[ld][5], poptree[rd][5])
                    poptree[i][5] = acolor
                else:
                    notdone = True
        if notdone is False:
            break
    return poptree


def yline(args,y, farright, width, dash, grayamount):
    """ draw a line at a specific height in relative terms """
    aline(args,[[0, y], [farright, y]], width, dash, grayamount)


def centerbox(args,pop, leftpoint, rightpoint, poptree, popxvals):
    """
        centerbox is a recursive function to find locations of population boxes on the x axis
        pop is the population for which we are finding the left and right sides of the box

        centerbox returns:
            width - the difference between right and left side of population it was called with
            center - the center of the population it was called with
            new value of popxvals
                popxvals[pop] is set, they start out with values of 0 and box width, but get new values
            leftpoint
                leftpoint is the point less than which we cannot find values because a box can't be drawn there
                leftpoint is partly determined by the descendant population and partly by the population
            rightpoint
                the rightmost point of any descendant population

        Start at the bottom, go up the left side then the right, recursively
        take the width of the left side and the width of the right side,
        put a spacer between, add them, and find the center
        the width and center gets returned to the basal population

        receives the population #, leftmost side of the ancestor population, rightside, the tree and popxvals (which this modifies)

        leftpoint starts at 0 for left branch up from root
            starts at right side of what is returned from the left branch up from the root
          if terminal pop
            returns the width and center location for that terminal pop
        else
            makes recursive calls for both left and right descendants
                each returns a width and center
            overall width left width + spacer + right width
            overall center   left center + (right center - left center)/2
            returns the width and center location based on both descendants of that pop
          center goes in the middle
    """

    if poptree[pop][2] == -1:  # pop is a terminal population
        # at this point popxvals[pop] holds just the width of the box (i.e. popxvals[pop][0] is 0)
        popxvals[pop][1] = popxvals[pop][1] - popxvals[pop][0] + leftpoint
        popxvals[pop][0] = leftpoint
        return popxvals[pop][1] - popxvals[pop][0], leftpoint + (popxvals[pop][1] - popxvals[pop][0]) / 2.0, popxvals, leftpoint, popxvals[pop][1]
    else:
        popspacer = args.popboxspaceadj * POPBOX_SPACE_DEFAULT
        (lw, lc, popxvals, leftpoint, rightpoint) = centerbox(args,poptree[pop][2], leftpoint, rightpoint, poptree, popxvals)
        rleftpoint = rightpoint + popspacer
        (rwidth, rcenter, popxvals, rleftpoint, rightpoint) = centerbox(args,poptree[pop][3], rleftpoint, rightpoint, poptree, popxvals)
        newwidth = lw + popspacer + rwidth

        newwidth = popxvals[pop][1] - popxvals[pop][0]
        newcenter = lc + (rcenter - lc) / 2.0
        if newcenter - (newwidth / 2.0) < leftpoint:
            newcenter += leftpoint - (newcenter - (newwidth / 2.0))
        templeft = newcenter - newwidth / 2.0
        popxvals[pop][0] = templeft
        popxvals[pop][1] = templeft + newwidth
        return newwidth, newcenter, popxvals, leftpoint, rightpoint


def fround(val):
    """ crude rounding to a couple decimal points for a positive val"""
    if val == 0:
        return "0.0"
    lval = math.log10(val)
    if lval < 0:
        lval -= 1
    rval = -int(lval) + 2
    if lval > 3:
        return str(int(round(val, rval)))
    return str(round(val, rval))


def popadjustx(popxvals, minx_popbox):
    " shift box locations to the left, as needed to fit with minx_popbox"
    minx = popxvals[0][0]
    for i in range(1, len(popxvals)):
        if minx > popxvals[i][0]:
            minx = popxvals[i][0]
    for i in range(len(popxvals)):
        popxvals[i][0] -= (minx - minx_popbox)
        popxvals[i][1] -= (minx - minx_popbox)
    return popxvals


def setpopbox(args,ty, slist, scaledtime, rootpop, poptree):
    """popbox[i][0] is the lowerleft point of the box
            popbox[i][0][0] contains the xdimension for the left side of the box
            popbox[i][0][1] contains the y dimension for the bottom of the box
        popbox[i][1] is the upper right
            popbox[i][1][0] contains the xdimension for the right side of the box
            popbox[i][1][1] contains the y dimension for the top of the box
        slist[4] holds population size info
            slist[4][4] holds actual parameter names and values
                slist[4][4][i] holds actual parameter names and values for popsize i
                    slist[4][4][i][0] is the name
                    slist[4][4][i][1] is the estimate
                    slist[4][4][i][2] is lower 95%
                    slist[4][4][i][3] is upper  95%

    """
    wadjust = ""
    for i in range(numpops - 1):
        wadjust += "00"
    if(scaledtime != []):
        minx_popbox = textwide(args,wadjust + "0.00 MYR", TFACTOR)
    else:
        minx_popbox = textwide(args,wadjust + "0.00 tu", TFACTOR)
    minx_popbox /= args.globalscale
    if args.localxscale > 0:
        minx_popbox /= args.localxscale

    popxvals = []
    for i in range(2 * numpops - 1):   ## left side temporarily at zero, right side temporarily at upper confidence interval
        popxvals.append([0, slist[4][4][i][1]])
    (width, c, popxvals, leftpoint, rightpoint) = centerbox(args,rootpop, 0, popxvals[rootpop][1], poptree, popxvals)
    popxvals = popadjustx(popxvals, minx_popbox)
    popbox = []

    # maxwide will be used to adjust the width as a scaler  so the part furthest to the right is not too far out
    maxwide = 0
    for i in range(2 * numpops - 1):
        if maxwide < (popxvals[i][1] + (slist[4][4][i][3] - slist[4][4][i][1])):
            maxwide = (popxvals[i][1] + (slist[4][4][i][3] - slist[4][4][i][1]))
    maxwide = maxwide / (1.0 - minx_popbox)

    if args.localxscale > 0:
        maxwide *= args.localxscale

    farright = 0
    confint = []
    for i in range(2 * numpops - 1):
        confint.append([])
        confint[i].append(minx_popbox + ((popxvals[i][1] - (slist[4][4][i][1] - slist[4][4][i][2])) / maxwide))
        confint[i].append(minx_popbox + ((popxvals[i][1] + (slist[4][4][i][3] - slist[4][4][i][1])) / maxwide))
        if confint[i][1] > farright:
            farright = confint[i][1]
        popbox.append([[], []])
        popbox[i][0].append(minx_popbox + popxvals[i][0] / maxwide)
        popbox[i][1].append(minx_popbox + popxvals[i][1] / maxwide)
        if poptree[i][1] == -1:
            popbox[i][0].append(args.lineINFy)
        else:
            popbox[i][0].append(ty[poptree[i][1] - 1][0])
        if poptree[i][0] == 0:
            popbox[i][1].append(args.line0y)
        else:
            popbox[i][1].append(ty[poptree[i][0] - 1][0])
    return popbox, maxwide, confint, farright


def printpopbox(args,popbox, maxwide, confint, slist, plist, rootpop, poptree, ty, scaledpop, droppops):
    """
        print popbox representing populations in different time periods

        popbox contains the corner locations of all the population boxes (see setpopbox())
        maxwide ?  maybe the largest x value for a populations's  upper confidence interval box
        confint contains the x value locations of the two confidence interval boxes for a population
        slist is the large information array built by reading the IMa3 output file
        plist is 2D list of population numbers plist[i][j] is the number of the jth population in interval i
        rootpop is the number of the pop that is ancestral to all
        poptree is the array that holds the tree topology (and population color)  info
        ty is the array of y axis values associated with the splitting times
        scaledpop contains the Ne values if they were available
        droppops[i] contains the numbers of the two populations that join into an ancestor after interval i

        graylevel is global
        dashinterval is global
    """
    if args.simplecolor:
        color = args.blue
    else:
        color = args.black
    cdim = []
    for i in range(2 * numpops - 1):
        tempbox = [row[:] for row in popbox[i]] # copy 2d list
        tempbox[1][0] = popbox[i][1][0] - (slist[4][4][i][1] - slist[4][4][i][2]) / maxwide
        cdim.append(calccdim(args,-1, tempbox))
        w(args,"%%begin box %d" % i)
        cdimtemp = curvebox(args,cdim[i], popbox[i], 2.5, color, 0, i, 0, poptree)
        w(args,"%%done box %d" % i)
    cdim = []
    if args.popboxcintervalboxes: # print confidence interval boxes
        for i in range(2 * numpops - 1):
            tempbox = [row[:] for row in popbox[i]] # copy 2d list
            tempbox[1][0] = popbox[i][1][0] - (slist[4][4][i][1] - slist[4][4][i][2]) / maxwide
            w(args,"%%begin left confidence for box %d" % i)
            cdim.append(calccdim(args,-1, tempbox))
            cdimtemp = curvebox(args,cdim[i], tempbox, 1.5, color, args.graylevel, i, args.dashinterval, poptree)
            w(args,"%%done left confidence for box %d" % i)
            tempbox = [row[:] for row in popbox[i]] # copy 2d list
            tempbox[1][0] = popbox[i][1][0] + (slist[4][4][i][3] - slist[4][4][i][1]) / maxwide
            w(args,"%%begin right confidence for box %d" % i)
            cdimtemp =curvebox(args,cdim[i], tempbox, 1.5, color, args.graylevel, i, args.dashinterval, poptree)
            w(args,"%%done right confidence for box %d" % i)
    popprintinc = 0.01
    if args.usealtnames:
        namelist = args.altpopnames
##        print(namelist)
        if args.useghost and args.excludeghost is False and namelist[-1].upper() != "GHOST":
            assert len(namelist) < numpops
            namelist.append("Ghost")
    else:
        namelist = slist[2][4]
    if args.anglenames:
        angle = 30
    else:
        angle = 0
    for i in range(numpops): # population names are in slist[2][4]
        if poptree[i][1] == 1 and i== droppops[1][1]: # right side of most recent split
            dotext(args,[popbox[i][0][0] + (popbox[i][1][0] - popbox[i][0][0]) / 2, popbox[i][1][1] + popprintinc], namelist[i], angle, False)
        else:
            if (popbox[i][1][0] - popbox[i][0][0] > 0.15):  # if a wide box move the text in a bit
                dotext(args,[popbox[i][0][0] + (popbox[i][1][0] - popbox[i][0][0]) / 4, popbox[i][1][1] + popprintinc], namelist[i], angle, False)
            else:
                dotext(args,[popbox[i][0][0], popbox[i][1][1] + popprintinc], namelist[i], angle, False)
    popprintinc = 0.025
    if args.label_a_pops:
        for i in range(numpops, 2 * numpops - 1):
            dotext(args,[max(popbox[i][0][0], popbox[i][0][0] + (popbox[i][1][0] - popbox[i][0][0]) / 2.0 - popprintinc),
                popbox[i][0][1] + (popbox[i][1][1] - popbox[i][0][1]) / 2.0], "pop #" + str(i), 0, False)
## plot the confidence arrows for population boxes
    lastperiod = [0] * (2 * numpops - 1)
    for i in range(2 * numpops - 1):
        for j in range(len(plist)):
            for k in range(len(plist[j])):
                if plist[j][k] == i and j > lastperiod[i]:
                    lastperiod[i] = j
    periodposcount = [0] * numpops
    arrowheightinc = 0.006
    arrowheights = []
    for i in range(numpops):
        if i == 0:
            top = args.line0y
            bot = ty[i][0]
        else:
            top = ty[i - 1][0]
            if i == numpops - 1:
                bot = args.lineINFy
            else:
                bot = ty[i][0]

        if top - bot < 0.1:
            frac = 0.5
        else:
            frac = 0.8
        arrowheights.append(top - (top - bot) * frac)
    if args.popboxcintervalarrows:  # print confidence interval arrows
        for i in range(2 * numpops - 1):
            period = lastperiod[i]
            arrowheight = max(popbox[i][0][1], arrowheights[period] - periodposcount[period] * 2 * arrowheightinc)
            head = [confint[i][0], arrowheight]
            tail = [popbox[i][1][0], arrowheight]
            # head is tip of arrow to lower bound of confidence interval
            # if there is not room for the arrowhead, don't print it
            if tail[0] - head[0] > ARROWHEAD_WIDTH_DEFAULT:
                arrow(args,head, tail, 2, color)
            head = [confint[i][1], arrowheight]
            tail = [popbox[i][1][0], arrowheight]
            arrow(args,head, tail, 0, color)
            periodposcount[period] += 1
    if scaledpop != []:
        ane = scaledpop[rootpop] / 1000
        anes = fround(ane)
        dotext(args,[0.15, 0.05], " Ancestral Ne (thousands): " + anes, 0, False)
    else:
        dotext(args,[0.15, 0.05], " Ancestral 4Nu: " + str(slist[4][4][rootpop][1]), 0, False)

    if args.simplecolor:
        w(args,"0 0 0  setrgbcolor")
    return popbox


def set_tlines(args,slist):
    """
        line0y - default relative height of time 0
        eventimes - if True, space split times evenly
        lastt_lower_y - height of oldest split time, by default is 1 / (numpops + 1), else can be set by user
    """
    t = []
    for i in range(numpops - 1):
        t.append([slist[5][4][i][1], slist[5][4][i][2], slist[5][4][i][3]])  # [time, upper ci, lower ci]
    ty = []
    if args.localyscale == -1:
        yint = args.line0y - args.lastt_lower_y
        for i in range(numpops - 1):
            ty.append([])
            if args.eventimes == False:
                tmax = slist[5][4][numpops - 2][3] # bottom of confidence interval of largest(oldest) t
                for j in range(3):
                    ty[i].append(args.line0y - (t[i][j] * yint) / tmax)
            else:
##                ty[i].append(args.line0y - ((i + 1) / float(numpops + 1) * yint) / tmax)
                ty[i].append(args.line0y - yint * (i + 1) / float(numpops))
    else:
        timeumean = slist[7][4][1]
        scaleumean = slist[7][4][2]
        for i in range(numpops - 1):
            ty.append([])
            for j in range(3):
                ty[i].append(args.line0y - (t[i][j] * (scaleumean / timeumean / 1e6) * args.localyscale))
                if ty[i][j] < args.lineINFy:
                    print(" time line too low in graph, reduce local y scale (-y value) ")
        args.lastt_lower_y = ty[numpops - 2][2]
    return ty


def print_tlines(args,ty, slist, scaledtime, farright):
    """
        print the split time lines and confidence interval lines
        graylevel is global
    """
    xinc = 0.005
    yinc = 0.002
    if(scaledtime != []):
        if max(scaledtime) / 1e6 < 1.0:
            yearscaler = 1e3
            yearscalestring = " KYR"
        else:
            yearscaler = 1e6
            yearscalestring = " MYR"
    if args.eventimes is False:
        for i in range(numpops - 1):
            if (ty[i][1] > ty[i][0]):
                yline(args,ty[i][1], farright, 1, 2, args.graylevel)
            yline(args,ty[i][0], farright, 0.5, 0, 0)
            if (ty[i][2] < ty[i][0]):
                yline(args,ty[i][2], farright, 1, 2, args.graylevel)
            if(scaledtime != []):
                scaledtime[i] /= yearscaler
                mtime = round(scaledtime[i], -int(math.log10(scaledtime[i]) - 2))
                nstr = str(mtime) + yearscalestring
                dotext(args,[xinc * (i + 2), ty[i][0] + yinc], nstr, 0, False)
            else:
                nstr = fround(slist[5][4][i][1]) + "tu"
                dotext(args,[xinc * (i + 2), ty[i][0] + yinc], nstr, 0, False)
            if (ty[i][1] > ty[i][0]):
                arrow(args,[xinc * (i + 1), ty[i][1]], [xinc * (i + 1), ty[i][0]], 1, args.black)
            if (ty[i][2] < ty[i][0]):
                arrow(args,[xinc * (i + 1), ty[i][2]], [xinc * (i + 1), ty[i][0]], 3, args.black)
    else:
        for i in range(numpops - 1):
            yline(args,ty[i][0], farright, 0.5, 0, 0)
            if(scaledtime != []):
                scaledtime[i] /= yearscaler
                mtime = round(scaledtime[i], -int(math.log10(scaledtime[i]) - 2))
                nstr = str(mtime) + yearscalestring
                dotext(args,[xinc * (i + 2), ty[i][0] + yinc], nstr, 0, False)
            else:
                nstr = fround(slist[5][4][i][1]) + "tu"
                dotext(args,[xinc * (i + 2), ty[i][0] + yinc], nstr, 0, False)
    return ty


def print_mcurves(args,slist, popbox, plist):
    """migration arrows:
    note - migration arrows are drawn in the forward direction!!
    likelihood ratio=ratio of the highest probability to the probability at 2NM = 0
    Sinficant likelihood ratios:
    2.70554  at p=0.05   The ratio of probabilities (as opposed to twice the log ratio) is 3.86813
    5.41189	  at p = 0.01  the ratio of prbabilities is 14.9685
    9.54954	 at p = 0.001  the ration of probabilities is 118.483
    3.86813 <= ratio <= 14.9685 upper arrow is a dash  (0.95 on chi square 50% 0.0 and 50% 1df)
    14.9685 <= ratio <= 118.483  upper arrow is a dotted (0.99 on chi square 50% 0.0 and 50% 1df)
    118.483 <= ratio upper arrow is a solid line       (0.999 on chi square 50% 0.0 and 50% 1df)

    list of things in miginfo[i]
    0 topop
    1 frompop
    2 direction
    3 period
    4 the number in this period
    5 2NM est
    6 log likelihood ratio stat
    also save # events to print in the period"""
    def checkm(val2NM, valm, llr):
##        return  (args.moption == 'a' and val2NM > min2NM) or (args.moption == 's' and llr >= 2.74)  or val2NM > args.moption
##  can happen that singificant m  is associated with 2NM = 0 (or near).  Catch these cases and do not print
        minsig2NM = 0.001
        if type(args.moption) is float:
            return val2NM >= args.moption
        else:
            returnval = ((args.moption == 'a' and valm > 0 and val2NM > 0)
                    or (args.moption == 's' and llr.count('*') > 0 and val2NM >= minsig2NM)
                    or (args.moption == 'S' and llr.count('*') > 1) and val2NM >= minsig2NM)
            return returnval
    if args.moption == 'x':
        return
    mperiodnum = [0] * (numpops - 1)
    if len(slist[6]) > 4:
        sml = slist[6][4]
        miginfo = []
        mi = 0
        for i in range(len(sml)):
            llr = sml[i][2]

            usem = checkm(sml[i][1], sml[i][3], llr)
            if usem:
                miginfo.append([])
                c1 = max(sml[i][0].find("M"), sml[i][0].find("m"))  # either upper of lower case
                c2 = sml[i][0].find(">")
                miginfo[mi].append(int(sml[i][0][c2 + 1:len(sml[i][0])]))
                miginfo[mi].append(int(sml[i][0][c1 + 1:c2]))
                pos1 = -1
                pos2 = -1
                period = 0
                while 1:
                    for j in range(len(plist[period])):
                        if plist[period][j] == miginfo[mi][0]:
                            pos1 = j
                        if plist[period][j] == miginfo[mi][1]:
                            pos2 = j
                    if pos1 >= 0 and pos2 >= 0:
                        if pos1 < pos2:
                            direction = 0    # arrow points to right
                        else:
                            direction = 2     # arrow points to left
                        break
                    else:
                        period += 1
                        pos1 = -1
                        pos2 = -1
                miginfo[mi].append(direction)
                miginfo[mi].append(period)
                miginfo[mi].append(mperiodnum[period])
                mperiodnum[period] += 1
                miginfo[mi].append(sml[i][1])
                miginfo[mi].append(llr)
                mi += 1
            else:
                continue

        wideboxfrac = 0.4
        narrowboxfrac = 0.8
        # set height of curves
        y = []
        for i in range(len(miginfo)):
            frompop = miginfo[i][0]
            period = miginfo[i][3]
            hi = popbox[frompop][1][1]
            for j in range(len(plist[period])):
                if hi > popbox[plist[period][j]][1][1]:
                    hi = popbox[plist[period][j]][1][1]
            lo = 0
            for j in range(len(plist[period])):
                if lo < popbox[plist[period][j]][0][1]:
                    lo = popbox[plist[period][j]][0][1]
            y.append(hi - (hi - lo) * (miginfo[i][4] + 1) / (mperiodnum[miginfo[i][3]] + 1))
        for i in range(len(miginfo)):
            frompop = miginfo[i][0]
            topop = miginfo[i][1]
            period = miginfo[i][3]
            direc = miginfo[i][2]
            val2NM = fround(miginfo[i][5])
            if (miginfo[i][6] != 'ns'):
                val2NM += miginfo[i][6]
            text2NMwidth = textwide(args,val2NM, 1.5)
            if direc == 0:
                tailx = popbox[frompop][1][0] - (popbox[frompop][1][0] - popbox[frompop][0][0]) * wideboxfrac
                headx = popbox[topop][0][0] + (popbox[topop][1][0] - popbox[topop][0][0]) * wideboxfrac
                if (text2NMwidth > abs(tailx - headx)):
                    tailx = popbox[frompop][1][0] - (popbox[frompop][1][0] - popbox[frompop][0][0]) * narrowboxfrac
                    headx = popbox[topop][0][0] + (popbox[topop][1][0] - popbox[topop][0][0]) * narrowboxfrac
            if direc == 2:
                tailx = popbox[frompop][0][0] + (popbox[frompop][1][0] - popbox[frompop][0][0]) * wideboxfrac
                headx = popbox[topop][1][0] - (popbox[topop][1][0] - popbox[topop][0][0]) * wideboxfrac
                if (text2NMwidth > abs(tailx - headx)):
                    tailx = popbox[frompop][0][0] + (popbox[frompop][1][0] - popbox[frompop][0][0]) * narrowboxfrac
                    headx = popbox[topop][1][0] - (popbox[topop][1][0] - popbox[topop][0][0]) * narrowboxfrac
            if args.rgbcolor:
                migrationstraightarrow(args,val2NM, [headx, y[i]], [tailx, y[i]], direc, args.darkgreen)
            else:
                migrationstraightarrow(args,val2NM, [headx, y[i]], [tailx, y[i]], direc, args.red)


##***********************************************************************************
##////////////// PARSING and MAIN PROGRAM ///////////////////////////////////////////
##***********************************************************************************

def parse_arguments(args=None):
    """
    Parse command line arguments using argparse
    Returns parsed arguments and command string
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description=f"IMfig program. Copyright 2009-{RELEASE_YEAR} Jody Hey. Release Date {RELEASE_DATE}",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Required argument
    parser.add_argument('-i', '--input', required=True, dest='imfilename', help='input file name')
    # Optional arguments
    parser.add_argument('-a', '--label-ancestor-pops', action='store_true', dest='label_a_pops',
                        help='include ancestral population #\'s in plot')
    parser.add_argument('-b', '--box-spacing', type=float, default=1.0, dest='popboxspaceadj',
                        help='adjust width spacing of population boxes, values > 0, default = 1')
    # Add check for PIL before adding image format options
    if check_PIL:
        parser.add_argument('-c', '--convert', choices=['j', 'p', 'n'], dest='imageformat',
                           help='output format, default is eps, see also -w\n'
                                '-c j : make a jpeg file\n'
                                '-c p : make a pdf file\n'
                                '-c n : make a png file')
    
    parser.add_argument('-d', '--no-demographic-scale', action='store_true', dest='skipdemographicscaling',
                        help='do not use demographic scale information even if in input file')
    parser.add_argument('-e', '--even-split-times', action='store_true', dest='eventimes',
                        help='space split times evenly (not proportional to time, no confidence intervals shown)')
    parser.add_argument('-f', '--font', default='Arial', dest='font',
                        help='font. Default=Arial. Use postscript fonts available on the computer\n'
                             'e.g. Arial, Helvetica, Times-roman, Courier')
    parser.add_argument('-g', '--global-scale', type=float, default=1.0, dest='globalscale',
                        help='global plot scale sets the size of the plot, max = 1, default = 1')
    parser.add_argument('-j', '--arrow-width', type=float, default=1.0, dest='arrowheightadj',
                        help='arrow width, default = 1')
    parser.add_argument('-k', '--angled-names', action='store_true', dest='anglenames',
                        help='print population names on an angle')
    parser.add_argument('-l', '--height-scale', type=float, dest='heightscale',
                        help='expand/shrink height by a positive scalar, >1 means taller, <1 means shorter')
    parser.add_argument('-m', '--migration', default='s', dest='moption',
                        help='options for printing of arrows and 2Nm values for migration:\n'
                             '-m x : do not print migration arrows\n'
                             '-m a : 2Nm migration arrows for all cases when both m > 0 and 2Nm > 0\n'
                             '-m s : 2Nm migration arrows only if m is statistically significant p <= 0.05 (default)\n'
                             '-m S : 2Nm migration arrows only if m is statistically significant p <= 0.01\n'
                             '-m # : "#" is a number, migration arrows appear when 2NM >= # (e.g. -m0.1)')
    parser.add_argument('-n', '--alt-names', dest='altnamefilename',
                        help='file with alternative species names')
    parser.add_argument('-o', '--output', default='imfig_output.eps', dest='outputfilename',
                        help='output file name, default is imfig_output.eps')
    parser.add_argument('-p', '--font-size', type=float, dest='fontsize',
                        help='fontsize (default is 14 for full scale, default follows global scale)')
    parser.add_argument('-q', '--no-confidence-interval-boxes', action='store_true', dest='no_popboxcintervalboxes',
                        help='no confidence interval boxes for population boxes printed')
    parser.add_argument('-r', '--no-confidence-interval-arrows', action='store_true', dest='no_popboxcintervalarrows',
                        help='no confidence interval arrows for population boxes printed')
    parser.add_argument('-s', '--square', action='store_true', dest='dosquare',
                        help='print square, rather than landscape')
    parser.add_argument('-t', '--time-height', type=float, dest='lastt_lower_y',
                        help='relative height of oldest time point, values between 0 and 1\n'
                             'default value = 1/(# sampled populations+1)')
    parser.add_argument('-u', '--simple-colors', action='store_true', dest='simplecolor',
                        help='simple colors, blue for population boxes, red arrows for migration (default grayscale)')
    parser.add_argument('-v', '--color', action='store_true', dest='rgbcolor',
                        help='multiple colors for population boxes, red arrows for migration (default grayscale)')
    if check_PIL:
        parser.add_argument('-w', '--image-width', type=int, dest='widthscalar',
                           help='file image width, integer multiple of 720 pixels (only if using -c)')
    parser.add_argument('-x', '--width-adjust', type=float, dest='xadjust',
                        help='expand/shrink width of plot by a positive scalar, >1 means wider, <1 means narrower')
    parser.add_argument('-y', '--height-adjust', type=float, default=-1, dest='localyscale',
                        help='adjust height of splittimes, relative to bottom of figure, max = 1.')
    parser.add_argument('-z', '--exclude-ghost', action='store_true', dest='excludeghost',
                        help='exclude the ghost population from the figure')
    
    # Parse arguments
    if args is None:
        args = sys.argv[1:]
    
    parsed_args = parser.parse_args(args)
    cmdstr = " ".join(["IMfig.py"] + args)
    
    return parsed_args, cmdstr


def run_imfig(args, cmdstr):
    """
    Main processing function using the configuration object.
    Replaces the dostuff() function.
    """
    global numpops
    
    # Get info from the input file
    print(f"input file: {args.imfilename}")
    if args.imagefileextension != "":
        tempname = args.outputfilename[0:-4] + args.imagefileextension
    else:
        tempname = args.outputfilename
    print(f"output file: {tempname}")
    print("Reading input file...")
    
    # Read the input file
    args,slist, scaledpop, scaledtime = readimfile(args)
    
    # Read the tree, set up plist
    poptree, rootpop, poptreestring, plist, droppops, addpop = poptreeread(slist[3][4])
    if args.rgbcolor:
        poptree = addcolors(poptree)
    
    # Set scales
    print("Setting scales...")
    args.adjust_scales()
    if args.set_lastt_lower_y:
        args.lastt_lower_y = 1.0 / (numpops + 1)
    
    # Set up time lines
    ty = set_tlines(args, slist)
    
    # Set up population boxes
    popbox, maxwide, confint, farright = setpopbox(args, ty, slist, scaledtime, rootpop, poptree)
    
    # Write the EPS output file
    args.epsf = open(args.outputfilename, "w")
    w(args, "%!PS-Adobe-3.0 EPSF-3.0")
    w(args, "%%legal size in landscape is 792x612 set bounding box with 0.5inch margins")
    w(args, "%%the lower corner is at 36 36, x dim is 720 wide, y dim is 540 hi")
    w(args, f"%%BoundingBox: {int(args.fixedLL[0])} {int(args.fixedLL[1])} {int(args.fixedUR[0])} {int(args.fixedUR[1])}")
    w(args, f"%%IMfig program author: Jody Hey   Copyright 2009-{RELEASE_YEAR}")    
    w(args, f"%%Command line for IMfig program that generated this file: {cmdstr}")
    
    # Generate the figure
    print("Creating figure...")
    print("Plotting splitting times...")
    ty = print_tlines(args, ty, slist, scaledtime, farright)
    
    print("Plotting population boxes...")
    popbox = printpopbox(args, popbox, maxwide, confint, slist, plist, rootpop, poptree, ty, scaledpop, droppops)
    
    print("Plotting migration arrows...")
    print_mcurves(args, slist, popbox, plist)
    
    # Close output file
    args.epsf.close()
    print("Plot completed")
    
    # Convert to image format if requested
    if args.imagefileextension != "":
        success = write_image_file(args)
        if success:
            print("Image file created")


def write_image_file(args):
    """Convert EPS file to the specified image format using PIL."""
    fn, tempext = os.path.splitext(args.outputfilename)
    if args.widthscalar < 1:
        args.widthscalar = 1
    try:
        im = Image.open(args.outputfilename)
        im.load(scale=args.widthscalar)
    except Exception as e:
        print("--- An exception occurred ---")
        print(f"Error reading {args.outputfilename}: {e}")
        print(f"{args.imagefileextension} file not written")
        print(f"Exception Type: {type(e)}")
        print(f"Exception Message: {e}")        
        # This prints the exception type, message, and the call stack
        traceback.print_exc(file=sys.stdout) # Print to standard output
        print("---------------------------")        
        return False
    
    outfn = fn + args.imagefileextension
    try:
        im.save(outfn)
    except Exception as e:
        print(f"Error converting {args.outputfilename} to {args.imagefileextension}: {e}")
        return False
    
    try:
        im.close()
        os.remove(args.outputfilename)
    except Exception as e:
        print(f"Warning: Could not clean up temporary file: {e}")
    
    return True


def write_caption(args, cmdstr):
    """Write a caption file for the generated figure."""
    fn = args.outputfilename[0:-4] + "_caption.txt"
    with open(fn, 'w') as f:
        f.write(f"IMfig program Copyright 2009-{RELEASE_YEAR} Jody Hey\n")
        f.write(f"command line string: {cmdstr}\n\n")
        s = ""
        if args.imaversion == 3:
            s += "Figure ?. A representation of an estimated Isolation with Migration model generated by IMa3 and the IMfig program (Hey et al., 2018). "
        else:
            s += "Figure ?. A representation of an estimated Isolation with Migration model generated by IMa2 and the IMfig program (Hey 2010). "
        s += "The phylogeny is depicted as a series of boxes organized hierarchically, with ancestor boxes positioned in between the corresponding descendants, and the width of boxes proportional to estimated Ne."
        if args.popboxcintervalboxes:
            s += " 95% confidence intervals for each Ne value are shown as dashed lines to the right of the left side of the corresponding population box. "
            if args.popboxcintervalarrows:
                s += " Gray arrows to the 95% Ne intervals are also shown extending to the left and right of the right boundary of each population box.  "
        elif args.popboxcintervalarrows:
            s += " 95% confidence intervals for Ne values are shown as gray arrows extending to the left and right of the right boundary of each population box.  "
        s += "Splitting times"
        if args.eventimes:
            s += ", positioned at even intervals, "
        s += " are depicted as solid horizontal lines, with text values on the left. "
        if not args.eventimes:
            s += "Confidence intervals for splitting times are shown as vertical gray arrows on the left, and parallel dashed lines. "
        if args.moption != 'x':
            s += "Migration arrows (if shown) indicate estimated 2Nm values from one population to another over the time interval when both populations exist. "
            s += "Arrows are shown only for estimated migration rates  "
            if args.moption == 's':
                s += "that are statistically significant (Nielsen and Wakeley, 2001) at or above the 0.05 level (* p < 0.05, ** p< 0.01, *** p < 0.001). "
            else:
                if args.moption == 'S':
                    s += " that are statistically significant (Nielsen and Wakeley, 2001) at or above the 0.01 level (** p< 0.01, *** p < 0.001). "
                else:
                    if args.moption != 'a':
                        s += f" that are above {args.moption:.3f}. "
        if args.useghost and args.excludeghost:
            s += " The ghost population is not shown in this figure. "
        if args.label_a_pops:
            s += "Ancestral population numbers are shown in ancestral boxes. "
        if args.skipdemographicscaling:
            s += " Population size (Ne) and splitting times values are scaled by the geometric mean of the mutation rates of the loci used for the analysis. "
        if args.imaversion == 3:
            s += "\n\nHey J, Chung Y, Sethuraman A, Lachance J, Tishkoff SA, Soudsa VC, Wang Y. 2018. Phylogeny Estimation by Integration over Isolation with Migration Models. Mol Biol Evol in press.\n"
        else:
            s += "\n\nHey J. 2010. The Divergence of Chimpanzee Species and Subspecies as Revealed in Multipopulation Isolation-with-Migration Analyses. Mol Biol Evol 27:921-933.\n"
        if args.moption.upper() == 'S':
            s += "Nielsen R, Wakeley J. 2001. Distinguishing migration from isolation. A Markov chain Monte Carlo approach. Genetics 158:885-896.\n"
        f.write(s)
    print("Caption file written")



if __name__ == "__main__":    
    """Main entry point for the IMfig program."""
    print(f"IMfig program. Copyright 2009-{RELEASE_YEAR} Jody Hey Release Date {RELEASE_DATE}")
    
    parsed_args, cmdstr = parse_arguments(sys.argv[1:])
    
    # Create configuration object
    args = Config()
    args.update_from_args(parsed_args)
    
    # Process the input file and generate the figure
    run_imfig(args, cmdstr)
    
    # Write caption
    write_caption(args, cmdstr)


