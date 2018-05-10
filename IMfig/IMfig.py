"""
    Copyright 2009-2018  Jody Hey
    IMfig makes an eps file containing a figure of the population phylogeny
        in an Isolation-with-Migration framework
    To see the helpscreen run the following at a command prompt:
        python IMfig.py
    If desired IMfig can be run from within an editing environment
        by setting 'cmdstr' (see bottom of this file)
    Check releasedate.
    Tested using python 3.6  but probably runs ok in python 2.7
    Read the documentation for details on running the program
"""

import math
import sys
import os
## some users won't have colormath
try:
    from colormath.color_objects import LabColor, sRGBColor
    from colormath.color_conversions import convert_color
    from colormath.color_objects import SpectralColor
    check_colormath = True
except ImportError:
    check_colormath = False

releasedate = "Jan 2, 2018"

# global variables.  gv is a dictionary that holds nearly all of them
gv = {}  ## dictionary to hold the many global constants, 'gv' for global variables
numpops = 0  ## widely used global, gets reset when reading file

# constants that don't change (some can be adjusted in effect by multiplying by scalars given by user)
arrowheadwidthdefault = 0.01  ## arrow size
popboxspacedefault = 0.1  ## spacing between population boxes
curveheightdefault = 0.03  ## curvature of migration arrows
tfactor = 1.0  ## a fudge factor for moving things to the right of splittime arrows
min2NM = 0.0  ## not sure why this was 0.0005 ## the smallest value plotted by IM programs for 2NM


##***********************************************************************************
##////////////// FUNCTIONS FOR GENERATING EPS FILE   ////////////////////////////////
##***********************************************************************************

def w(s):
    """
        simple function to make it easier to print to the eps file without having to repeat code text
    """
    gv["epsf"].write(s + "\n")

def apoint(rpoint):
    """ rpoint is a list of length 2, convert a relative point to an absolute point
        x value is in position 0,  y value in position 1
        some other global variables that are used here:
            fixedLL is the lower left point of the plot - no values to left or below this.
            fixedUR is the upper right point - not values above or to the right of this.
            globalscale is an overall scalar of plot size
            localxscale is an x dimensional scalar of plot size
                so the x dimension can be changed without affecting the y dimension
    """
    tempy = gv["fixedLL"][1] + gv["globalscale"]*rpoint[1]*(gv["fixedUR"][1]-gv["fixedLL"][1])
    if gv["localxscale"] != -1:
        tempx = gv["fixedLL"][0] + gv["localxscale"]*gv["globalscale"]*rpoint[0]*(gv["fixedUR"][0]-gv["fixedLL"][0])
    else:
        tempx = gv["fixedLL"][0] + gv["globalscale"]*rpoint[0]*(gv["fixedUR"][0]-gv["fixedLL"][0])
    if tempx - gv["fixedUR"][0] > 0 and tempx - gv["fixedUR"][0] < 1e-7:
        tempx = gv["fixedUR"][0]
    if tempx > gv["fixedUR"][0]:
        print ( "problem x value : ",tempx,   " max x allowed : ",gv["fixedUR"][0])
    return [tempx,tempy]

def rapoint(rpoint):
    """ relative point
        this is called from a function where the scale has been reset
    """
    return [rpoint[0]*gv["globalscale"]*(gv["fixedUR"][0]-gv["fixedLL"][0]),
            rpoint[1]*gv["globalscale"]*(gv["fixedUR"][1]-gv["fixedLL"][1])]

def textwide(s, tf):
    """ approx width of text """
    width = 350  ## default ok for Arial or Helvetica
    if gv["font"] == "Times-roman":
        width = 330
    if gv["font"] == "Courier":
        width = 390
    if gv["fontfixed"] is False:
        localfontsize = int(gv["fontsize"]*gv["globalscale"])
    else:
        localfontsize = int(gv["fontsize"])
    return tf*localfontsize * len(s)*width/(1000*(gv["fixedUR"][0] - gv["fixedLL"][0]))

def dotext(rpoint, text, angle, bi):
    """
        print text beginning at rpoint at angle
        font and bifont are global
    """
##    w("/Arial findfont")
    if bi:
        w("/%s findfont" % gv["bifont"])
    else:
        w("/%s findfont" % gv["font"])
    if gv["fontfixed"] is False:
        localfontsize = int(gv["fontsize"]*gv["globalscale"])
    else:
        localfontsize = int(gv["fontsize"])
    w("%d scalefont" % localfontsize)
    w("setfont")
    w("newpath")
    p = apoint(rpoint)
    if angle != 0:
        w("gsave")
        w("%d %d translate" % (p[0], p[1]))
        w("%d  rotate" % angle)
        w("0  0 moveto")
        w("(" + text + ") show")
        w("grestore")
    else:
        w("%d %d moveto" % (p[0],p[1]))
        w("(" + text + ") show")

def curvecontrol(p1,p2, u_or_d):
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
##    start by converting D to A,  and C to B
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
    if p1[0] < p2[0]:   ## type A
        if u_or_d:   ## curve up
            cp1.append( ((p2[0]-p1[0]) * e1) + p1[0])
            cp1.append( ((p2[1]-p1[1]) * e2) + p1[1])
            cp2.append( ((p2[0]-p1[0]) * e2c) + p1[0])
            cp2.append( ((p2[1]-p1[1]) * e1c) + p1[1])
        else:
            cp1.append( ((p2[0]-p1[0]) * e2) + p1[0])
            cp1.append( ((p2[1]-p1[1]) * e1) + p1[1])
            cp2.append( ((p2[0]-p1[0]) * e1c) + p1[0])
            cp2.append( ((p2[1]-p1[1]) * e2c) + p1[1])
    else:  ## type B
        if u_or_d:   ## curve up
            cp1.append( p1[0]-((p1[0]-p2[0]) * e1))
            cp1.append( ((p2[1]-p1[1]) * e2) + p1[1])
            cp2.append( p1[0] - ((p1[0]-p2[0]) * e2c))
            cp2.append( ((p2[1]-p1[1]) * e1c) + p1[1])
        else:
            cp1.append( p1[0]-((p1[0]-p2[0]) * e2))
            cp1.append( ((p2[1]-p1[1]) * e1) + p1[1])
            cp2.append( p1[0]-((p1[0]-p2[0]) * e1c))
            cp2.append( ((p2[1]-p1[1]) * e2c) + p1[1])
    if resort:
        ptemp = cp2
        cp2 = cp1
        cp1 = ptemp
    return cp1,cp2

def calccdim(cdimval,cbox):
    ll = cbox[0]
    ur = cbox[1]
    curvesizedefine = 0.02
    if cdimval == -1:
        cdimval = (gv["fixedUR"][0]-gv["fixedLL"][0]) *curvesizedefine
    lla = apoint(ll)
    ura = apoint(ur)
    if ura[0]-lla[0] < 2*cdimval or  ura[1]-lla[1] < 2* cdimval:
        cdimval = min(ura[0]-lla[0],ura[1]-lla[1])/2.0
    return cdimval

def curvebox(cdim, cbox, width,color,grayamount, popnum,dash,poptree):
    """
        creates a box with curved corners, size of the curve set by curvesize
        if dash==0 and rgbcolor == True, fills the box with a lighter version of
        the rgbcolor for the population the box is for
        returns cdim which has something to do with the size of the box
    """
    if dash > 0:
        w("[%d %d] 0 setdash" % (dash,dash))
    ll = cbox[0]
    ur = cbox[1]
    curvesizedefine = 0.02
    if cdim == -1:
        cdim = (gv["fixedUR"][0]-gv["fixedLL"][0]) *curvesizedefine
    lla = apoint(ll)
    ura = apoint(ur)
    if ura[0]-lla[0] < 2*cdim or  ura[1]-lla[1] < 2* cdim:
        cdim = min(ura[0]-lla[0],ura[1]-lla[1])/2.0
    ula = [lla[0],ura[1]]
    lra = [ura[0],lla[1]]
    if gv["rgbcolor"]:
        boxcolorstring = ( "%f %f %f setrgbcolor" %
                (poptree[popnum][5][0],poptree[popnum][5][1],poptree[popnum][5][2]))
        lightcolor = []
        for ii in range(3):
            ## lighten to 10% of color
            lightcolor.append(1.0 - (0.1* (1.0 -poptree[popnum][5][ii])))
        lightboxfillcolorstring = "%f %f %f setrgbcolor" % (lightcolor[0],lightcolor[1],lightcolor[2])
    else:
        if color != gv["black"]:
            color = gv["blue"]
            gcolor = []
            for i in range(3):
                if color[i] == 0:
                    gcolor.append(grayamount)
                else:
                    gcolor.append(color[i])
            boxcolorstring = "%f %f %f setrgbcolor" % (gcolor[0],gcolor[1],gcolor[2])
        else:
            boxcolorstring = "%f setgray" % grayamount

    w("newpath")
    w("%d  %d  moveto" %(lla[0]+cdim,lla[1]))
    cp1 = [lra[0]-cdim,lra[1]]
    cp2 = [lra[0],lra[1]+cdim]
    w("%d  %d  lineto" %(cp1[0],cp1[1]))
    ccpoints = curvecontrol(cp1,cp2,0)
    w("%d %d %d %d %d  %d  curveto" %(ccpoints[0][0],ccpoints[0][1],ccpoints[1][0],ccpoints[1][1],cp2[0],cp2[1]))
    cp1=[ura[0],ura[1]-cdim]
    cp2 = [ura[0]-cdim,ura[1]]
    w("%d  %d  lineto" %(cp1[0],cp1[1]))
    ccpoints = curvecontrol(cp1,cp2,1)
    w("%d %d %d %d %d  %d  curveto" %(ccpoints[0][0],ccpoints[0][1],ccpoints[1][0],ccpoints[1][1],cp2[0],cp2[1]))
    cp1 = [ula[0]+cdim,ula[1]]
    cp2 = [ula[0],ula[1]-cdim]
    w("%d  %d  lineto" %(cp1[0],cp1[1]))
    ccpoints = curvecontrol(cp1,cp2,1)
    w("%d %d %d %d %d  %d  curveto" %(ccpoints[0][0],ccpoints[0][1],ccpoints[1][0],ccpoints[1][1],cp2[0],cp2[1]))
    cp1 = [lla[0],lla[1]+cdim]
    cp2 = [lla[0]+cdim,lla[1]]
    w("%d  %d  lineto" %(cp1[0],cp1[1]))
    ccpoints = curvecontrol(cp1,cp2,0)
    w("%d %d %d %d %d  %d  curveto" %(ccpoints[0][0],ccpoints[0][1],ccpoints[1][0],ccpoints[1][1],cp2[0],cp2[1]))
    w("closepath")
    w("gsave")
    if gv["rgbcolor"] and dash == 0: ## fill the box with a lighter version of the color used for the lines of the box
        w(lightboxfillcolorstring)
        w("fill")
        w("grestore")
    width = float(width)
    w("%f setlinewidth" % (width*gv["globalscale"]))
    w(boxcolorstring)
    w("stroke")
    if gv["simplecolor"] or gv["rgbcolor"]:
        w("0 0 0  setrgbcolor")
    else:
        w("0 setgray")
    if dash > 0:
        w("[] 0 setdash")
    return cdim

def aline(p, width, dash, grayamount):
    """ p is a list of points in relative space (0-1)
        dash is the spacing (in point scale) of dashes in the line
        if dash is zero there is no dashing """
    if grayamount > 0:
        w("%f setgray" %grayamount)
    ap = []
    for i in range(len(p)):
        ap.append(apoint(p[i]))
    if dash > 0:
        w("[%d %d] 0 setdash" % (dash,dash))

    w("%d %d moveto" % (ap[0][0],ap[0][1]))
    for j in range(1,len(p)):
        w("%d %d lineto" % (ap[j][0],ap[j][1]))
    width*= gv["globalscale"]
    w("%f setlinewidth" % width)
    w("stroke")
    w("[ ] 0 setdash")
    if grayamount > 0:
        w("0 setgray")

def arrowhead(head,headwidth,angle):
    """ draw arrowhead width on the same scale as points in head
        head is the center of the arrowhead
        angle = 0 has the arrow pointing to the right
    """
    w("%% begin arrowhead")
    holdhead = apoint(head)
    head = [0,0]
    tip = rapoint([head[0] + headwidth,head[1]])
    p1 = rapoint([head[0] - headwidth,head[1] + headwidth])
    p2 = rapoint([head[0] - headwidth,head[1] - headwidth])
    c1 = rapoint([head[0],head[1]-headwidth/2])
    c2 = rapoint([head[0],head[1]+headwidth/2])
    w("gsave")
    w("%d %d translate" % (holdhead[0],holdhead[1]))
    w("%d  rotate" % angle)
    w("%d %d moveto" % (p1[0],p1[1]))
    w("%d %d lineto" % (tip[0],tip[1]))
    w("%d %d lineto" % (p2[0],p2[1]))
    w("%d %d %d %d %d %d curveto"% (c1[0],c1[1],c2[0],c2[1],p1[0],p1[1]))
    w("closepath")
    w("fill")
    w("grestore")
    w("%% end arrowhead")

def arrow(head,tail,direc, color):
    """
        draw an arrow. head and tail are points, width is on the same scale
        direc = 0 right, 1 up, 2 left, 3 down
        arrow is gray
    """
    headwidth = arrowheadwidthdefault*gv["arrowheightadj"]
    if (direc == 0):
        headadj = [head[0]-headwidth,head[1]]
    if (direc == 1):
        headadj = [head[0],head[1]-headwidth]
    if (direc == 2):
        headadj = [head[0]+headwidth,head[1]]
    if (direc == 3):
        headadj = [head[0],head[1]+headwidth]
    if color != gv["black"]:
        color = gv["blue"]
        gcolor = []
        for i in range(3):
            if color[i] == 0:
                gcolor.append(gv["graylevel"])
            else:
                gcolor.append(color[i])
        w("%f %f %f setrgbcolor" % (gcolor[0],gcolor[1],gcolor[2]))
    else:
        w("%f setgray" % gv["graylevel"])
    arrowhead(headadj,headwidth,direc*90)
    ahead = apoint(headadj)
    atail = apoint(tail)
    w("%d %d moveto" % (ahead[0],ahead[1]))
    w("%d %d lineto" % (atail[0],atail[1]))
    w("%f setlinewidth" % (2*gv["globalscale"]))
    w("stroke")
    if gv["simplecolor"] or gv["rgbcolor"]:
        w("0 0 0  setrgbcolor")
    else:
        w("0 setgray")

def migrationcurvearrow(val2NM,head,tail,direc, color):
    """ direct can be 0 or 2 (right or left)  if 0 curveheight is positive and curve goes up from
        the tail and then down to the head
        if direc is 2  then curve is interpreted to be negative and curve goes down from the tail
        and then up to the head """
    w("%% BEGIN MIGRATION ARROW: %s"%val2NM)
    curveheight = curveheightdefault
    c2height = arrowheadwidthdefault
    headwidth = c2height*1.5*gv["arrowheightadj"]
    width = 3.5
    if (direc == 0): ## arrow to the right,  line is shifted up, text is below line
        textpoint=[tail[0],tail[1]-curveheight]
        cheadadj = [head[0]-headwidth,head[1]+c2height]
        ctail =  [tail[0],tail[1]+c2height]
        arrowheadpoint = [cheadadj[0], head[1] + c2height/1.2]
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("%f %f %f setrgbcolor" % (color[0],color[1],color[2]))
        arrowhead(arrowheadpoint,headwidth,330)        ## head tilted down to the right
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("0 0 0 setrgbcolor")
        if abs(cheadadj[0] - ctail[0]) > 0:
            curveheightmultiplier =math.pow(abs(cheadadj[0] - ctail[0])/0.15,0.1)
        else:
            curveheightmultiplier = 1
        cp1 = [ctail[0] + (cheadadj[0] - ctail[0])*0.8,cheadadj[1] + curveheight*curveheightmultiplier]
        cp2 = [ctail[0] + (cheadadj[0] - ctail[0])*0.2,cheadadj[1] + curveheight*curveheightmultiplier]
        textpoint = [cp2[0],cheadadj[1]-curveheight/3]
    if (direc == 2): ## arrow to the left, line is shifted down, text is above line
        cheadadj = [head[0]+headwidth,head[1]]
        textpoint = [cheadadj[0]+c2height,cheadadj[1]]
        ctail = tail
        arrowheadpoint = [cheadadj[0], cheadadj[1] + c2height/3.5]
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("%f %f %f setrgbcolor" % (color[0],color[1],color[2]))
        arrowhead(arrowheadpoint,headwidth,150)       ## head tilted up to the left
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("0 0 0 setrgbcolor")
        if abs(cheadadj[0] - ctail[0]) > 0:
            curveheightmultiplier = math.pow(abs(cheadadj[0] - ctail[0])/0.15,0.1)
        else:
            curveheightmultiplier = 1

        cp1 = [cheadadj[0] + (ctail[0] - cheadadj[0])*0.2,cheadadj[1] - curveheight*curveheightmultiplier]
        cp2 = [cheadadj[0] + (ctail[0] - cheadadj[0])*0.8,cheadadj[1] - curveheight*curveheightmultiplier]
        textpoint = [cp1[0],cheadadj[1]-curveheight/3]

    ahead = apoint(cheadadj)
    atail = apoint(ctail)
    acp1 = apoint(cp1)
    acp2 = apoint(cp2)
    if width > 0:
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("%f %f %f setrgbcolor" % (color[0],color[1],color[2]))
        w("%f setlinewidth" % (width*gv["globalscale"]))
        w("%d %d moveto" % (ahead[0],ahead[1]))
        w("%d %d  %d  %d  %d  %d curveto" % (acp1[0],acp1[1],acp2[0],acp2[1],atail[0],atail[1]))
        w("stroke")
        ## stopped using the white line
        ## put a white line down middle of the migration arrow
##        if gv["simplecolor"] or gv["rgbcolor"]:
##            w("%f %f %f setrgbcolor" % (255,255,255))#0,0,0))
##        w("%f setlinewidth" % 0.5)
##        w("%d %d moveto" % (ahead[0],ahead[1]))
##        w("%d %d  %d  %d  %d  %d curveto" % (acp1[0],acp1[1],acp2[0],acp2[1],atail[0],atail[1]))
##        w("stroke")
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("0 0 0 setrgbcolor")
        dotext(textpoint,val2NM,0, True)
        if gv["simplecolor"] or gv["rgbcolor"]:
            w("0 0 0 setrgbcolor")
    w("%% END MIGRATION ARROW")


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
##    slist[i][3] - the string used to search the input file,  when it is found the function is called
##    slist[i][ > 3 ] - the actual information, the types and number of values vary depending on the category of information
##      all of the functions (names in slist[i][2] are called with
##            slist[i][2](imfile,imfileline,slist[i][3],numpops)
##        what the function returns is appended to slist[i]

def get_input_file_name (f, a,s):
    ## f not used but needed for function to match general function format
    return a[len(s):len(a)].strip()

def check_ghost_status(f,a,s):
    ## f, s not used but needed for function to match general function format
    global gv
    if gv["newercode"]:
        gv["useghost"] =  (a.find("-j") >= 0) and ( "1" in  a[a.find("-j")+1:])  ## should only be true if -j is there with a 1
    else:
        gv["useghost"] =  (a.find("-j") >= 0) and ( "4" in  a[a.find("-j")+1:])  ## should only be true if -j is there with a 4

def get_population_names (f,a,s):
    ## a, s not used but needed for function to match general function format
    """
        usealtnames and altnamefilename defined previously
    """
    global gv
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
    if gv["useghost"] == True and foundghost == False:  # for compatibility with older output files
        popnamelist.append('ghost')
    anames = []
    if gv["usealtnames"]:
        for line in open(gv["altnamefilename"],"r"):
            temp = line.strip()
            anames.append(temp)
        anames = anames[0:len(popnamelist)]
    gv["altpopnames"] = list(anames)
    return popnamelist

def get_population_tree (f,a,s):
    ## s not used but needed for function to match general function format
    """
         a couple possible things to read here
    """
    tempstring = a.split()[-1]
    while a.find("Population Tree") >= 0:
        if a.find("standard ordering") >= 0:
            tempstring = a.split()[-1]
        if a.find("Ghost Population") >= 0:
            assert gv["useghost"]
            tempstring = a.split()[-1]
        a = f.readline().strip()
    return tempstring

def get_popsize_param (f,a,s):
    ## a, s not used but needed for function to match general function format
    """ read the histogram table of marginal distributions for population sizes:
    For each population it reads:
    the label of the parameter
    the HiSmth value
    the HPD95Lo value
    the HPD95Hi value """
    psp = []
    for i in range(4):
        aa = f.readline().split()
    for i in range(2*numpops - 1):
        psp.append([])
        psp[i].append(aa[i+1])
    aa = f.readline().split()
    while len(aa) > 0:
        while aa.count("#"):
            aa.remove("#")
        while aa.count("#?"):
            aa.remove("#?")
        while aa.count("?"):
            aa.remove("?")
        found = False
        if aa[0]=="HiPt" or aa[0]=="HPD95Lo"  or aa[0]=="HPD95Hi":
            found = True
        if found:
            for i in range(2*numpops - 1):
                psp[i].append(float(aa[i+1].strip('?#')))
        aa = f.readline().split()
    return psp

def get_t_param (f,a,s):
    ## a, s not used but needed for function to match general function format
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
    for i in range(numpops-1):
        psp.append([])
        psp[i].append(aa[i+1])
    aa = f.readline().split()
    while len(aa)> 0:
        while aa.count("?"):
            aa.remove("?")
        while aa.count("#"):
            aa.remove("#")
        found = False
        if aa[0]=="HiSmth" or aa[0]=="HPD95Lo"  or aa[0]=="HPD95Hi":
            found = True
        if found:
            for i in range(numpops-1):
                psp[i].append(float(aa[i+1].strip('?#')))
        aa = f.readline().split()
    return psp

def msigvals(ss):
    """
        get the significance levels for the migration(m) rates
        make a list,  each element is a list
            the migration parameter
            the significance level
    """
    si = 0
    mlist = []
    nummp = 0
    while si < len(ss):
        if ss[si].find("Migration Rate Parameters") == 0:
            si += 1
            aa = ss[si].split()
            ainc = 1 if gv["newercode"] else 2
            for i in range(1,len(aa),ainc):
                mlist.append([aa[i]])
            si += 7 if gv["newercode"] else 6
            aa = ss[si].split()[1:]
            for i,temp in enumerate(aa):
                if (gv["moption"] == 's' and temp.count('*') == 0 ) or (gv["moption"] == 'S' and  temp.count('*') <= 1 ):
                    mlist[nummp+i].append('ns')
                else:
                    mlist[nummp+i].append( '*' * temp.count('*'))
            nummp = len(mlist)
        si += 1
    return mlist

def get_2NM (f,a,s):
    """
        replaces older code as of 9/18/2017
        reads a chunk of the file from "Marginal Peak Locations and Probabilities"  to "HISTOGRAMS"
        then scans this for the significance values for the migration (m) rates
        then gets the 2Nm numbers
        and then matches the m values with the 2Nm numbers
        returns a list, each element has 3 items:
            the name,  e.g. 2N0m0>1
            the estimated 2NM value
            a string with the significance level
    """
    ss = []
    while True:
        ss.append(f.readline().strip())
        if ss[-1].upper().find("HISTOGRAMS")== 0:
            break
    mlist = msigvals(ss)
    si = 0
    psp = []
    nummp = 0
    ainc = 1 if gv["newercode"] else 2
    while si < len(ss):
        if ss[si].upper().find("POPULATION MIGRATION (2NM) TERMS") == 0:
            si += 1
            aa = ss[si].split()
            for i in range(1,len(aa),ainc):
                psp.append([aa[i]])
            si += 3
            aa = ss[si].split()
            ii = 0
            for i in range(1,len(aa),ainc):
                psp[nummp + ii].append(float(aa[i]))
                ii += 1
            nummp = len(psp)
        si += 1

    for pi,p in enumerate(psp):
        mn = p[0][p[0].upper().find('M')+1:]
        i = 0
        while True:
            assert i < len(mlist)
            if mlist[i][0][1:] == mn:
                psp[pi].append(mlist[i][1])
                break
            i += 1
    return psp

def get_demog_scales (f,a,s):
    psp = [0,0,0]
    for i in range(10): ## go down several lines and look for the necessary information,  very crude and
        aa = f.readline().split()
        if aa[0]=="Generation" and aa[1]=="time":
            psp[0] = float(aa[len(aa)-1])
        if aa[0]=="Geometric" and aa[3]=="mutation":
            psp[1] = float(aa[len(aa)-1])
        if aa[0]=="Geometric" and aa[3]=="ML":
            psp[2] = float(aa[len(aa)-1])
    return psp

def get_parameter_priors (f,a,s):
    psp = [["population size","uniform"],["migration","uniform"],["splittime","uniform"]]
    aa = f.readline()
    for i in range(3):
        aa = f.readline().split()
        psp[i].append(float(aa[len(aa)-1]))
        if aa.count("exponential") > 0:
            psp[i][1] = "exponential"
    return psp


## localscale is 1 if (gv["fixedUR"][0]-gv["fixedLL"][0]) corresponds to Ne of 1e6
def calc_scaledvals(slist):
    gentime= slist[7][4][0]
    timeumean = slist[7][4][1]
    scaleumean = slist[7][4][2]
    scaledpop = []
    for i in range(2*numpops-1):
        scaledpop.append(slist[4][4][i][1]/(4.0 * timeumean*gentime))
    scaledtime = []
    for i in range(numpops-1):
        scaledtime.append(slist[5][4][i][1] * (scaleumean/timeumean))
    return scaledpop, scaledtime


def checkimcommandline(line): ## not sure what this is
    s = line.split('-')
    return

def mysplit(s,delims):
    """
        a simple function to split a string s at the occurence of all instances of each symbol in delims
        works simply by replacing every instance of each of those symbols  with a space
        and then calling regular split
        does not work with more than one symbol at a time
    """
    for c in delims:
        s = s.replace(c,' ')
    return s.split()

def readimfile():
    """
        gets info from the input file
        returns all information in slist
        and scale information in scaledpop and scaledtime
    """
    global numpops
    global gv
    imfile = open(gv["imfilename"],"r")
    gv["useghost"] = False
    imfileline  = imfile.readline()
    gv["newercode"] = False
    while imfileline != '':
        if imfileline.upper().find("IMa3 program compiled on".upper()) >= 0:
            import re
            from datetime import datetime
            linesplit = mysplit(imfileline.strip(),",")
            date = " ".join(linesplit[4:7])
            ## changed ima3 format afer sep 12 2017
            newcodedatetime = datetime.strptime("sep 12 2017", '%b %d %Y')
            filedatetime = datetime.strptime(date, '%b %d %Y')
            gv["newercode"] = filedatetime >= newcodedatetime
        if imfileline.upper().find("Command line string :".upper()) >= 0:
            checkimcommandline(imfileline[22:])
            break
        imfileline  = imfile.readline()
    slist = [["ghost status",False,check_ghost_status,"Model options on command line"],\
             ["inputfile",False,get_input_file_name,"Text from input file:"],\
             ["pop names",False,get_population_names,"Population Names"],\
             ["pop tree",False,get_population_tree,"Population Tree :"],\
             ["population size parameter info",False,get_popsize_param,"MARGINAL DISTRIBUTION VALUES AND HISTOGRAMS OF POPULATION SIZE AND MIGRATION PARAMETERS"],\
             ["splitting time parameter info",False,get_t_param,"MARGINAL DISTRIBUTION VALUES AND HISTOGRAMS OF PARAMETERS IN MCMC"],\
             ["migration parameter info",False,get_2NM,"Marginal Peak Locations and Probabilities"],\
             ["demographic scales",gv["skipdemographicscaling"],get_demog_scales,"MARGINAL DISTRIBUTION VALUES IN DEMOGRAPHIC UNITS"] #,\
##             ["parameter priors",False,get_parameter_priors,"Parameter Priors"] \  ignore this I think
             ]
    while imfileline != '':
        checkdone = True
        for i in range(len(slist)):
            checkdone =  checkdone and slist[i][1]
            if slist[i][1] == False and imfileline.upper().find(slist[i][3].upper()) >= 0:
                if slist[i][0] == "ghost status":
                    slist[i][2](imfile,imfileline,slist[i][3])
                else:
                    slist[i].append(slist[i][2](imfile,imfileline,slist[i][3]))
                slist[i][1] = True
                if slist[i][0] == "pop names":
                    numpops = len(slist[i][4])
        if checkdone:
            break
        imfileline  = imfile.readline()
    imfile.close()
    (scaledpop,scaledtime) = ([],[])
    if gv["skipdemographicscaling"]:
        slist[7][1] = False
    else:
        if len(slist[7]) == 4:
            print (  "**IMfig error - Information in demographic units not found, use -d option")
            printcommandset()
            quit()
        if len(slist[7][4])==3:
            (scaledpop, scaledtime) = calc_scaledvals(slist)
    return slist, scaledpop, scaledtime

##***********************************************************************************
##////////////// FUNCTIONS FOR READING THE POPULATION TREE STRING  //////////////////
##***********************************************************************************
def parenth(tempcurrent,poptree,poptreestring,stringspot,ancestralpopnums,rootpop,nextnode,periodi):
    current = ancestralpopnums[tempcurrent]
    stringspot += 1
    while poptreestring[stringspot].isspace():
        stringspot += 1
    while True:
        if poptreestring[stringspot].isdigit():
            if stringspot <= len(poptreestring)-2 and poptreestring[stringspot+1].isdigit():
                ts = poptreestring[stringspot] + poptreestring[stringspot+1]
                itemp = int(ts)
            else:
                itemp = int(poptreestring[stringspot])
            stringspot += 1
            if  poptree[current][2] == -1:
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
            if  poptree[current][2] == -1:
                poptree[current][2] =   ancestralpopnums[nextnode]
            else:
                poptree[current][3] =   ancestralpopnums[nextnode]
            (poptree,rootpop,stringspot,periodi,nextnode) = parenth(nextnode,poptree,poptreestring,stringspot,ancestralpopnums,rootpop,nextnode,periodi)
        if poptreestring[stringspot] == ')':
            break
    stringspot += 1
    if poptreestring[stringspot] == ':':
        stringspot += 1
        if stringspot <= len(poptreestring)-2 and poptreestring[stringspot+1].isdigit():
            ts = poptreestring[stringspot] + poptreestring[stringspot+1]
            i = int(ts)
        else:
            i = int(poptreestring[stringspot])
        if i < numpops:
            print ( " wrong number of ancestral populations indicated. string %c " % poptreestring[stringspot])
        periodi = i- numpops
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
        poptree[current][1] = -1
        rootpop = current

    return poptree,rootpop,stringspot,periodi,nextnode

def parenth0(current,poptree,poptreestring,stringspot,ancestralpopnums):
    nextlistspot = 0
    ne = stringspot
    popennum = 0
    psetlist = []
    for i in range(current):
        psetlist.append(-1)
    while ne < len(poptreestring):
        if poptreestring[ne]=='(':
            psetlist[nextlistspot] = popennum
            nextlistspot += 1
            popennum += 1
            ne += 1
        else:
            if poptreestring[ne]==')':
                ne += 2
                if ne <= len(poptreestring)-2 and poptreestring[ne+1].isdigit():
                    ts = poptreestring[ne] + poptreestring[ne+1]
                    itemp = int(ts)
                else:
                    itemp = int(poptreestring[ne])
                ancestralpopnums[current + psetlist[nextlistspot - 1]] = itemp
                nextlistspot -= 1
            else:
                ne += 1
    return poptree, ancestralpopnums

def set0 (strlist,pos):
    """ removes elements of a list from pos to the end, save these as a separate list """
    hold = []
    while len(strlist) > pos:
        hold.append(strlist.pop(len(strlist)-1))
    hold.reverse()
    return strlist,hold

def strlistadd(strlist,pos,c):
    if pos > (len(strlist)-1):
        strlist.append(c)
    else:
        strlist[pos] = c
    return strlist

def joinlist(list1, list2):
    for i in range(len(list2)):
        list1.append(list2[i])
    return list1


def rewrite (substr):
    """    rewrite() rewrites the treestring in a standard order
        swivels nodes,  if both have node sequence values, the one with the lower node sequence value (periodi[]) goes on the left
        if only one has a node sequence value,  it goes on the right
            when neither has a node sequence value, the one with the lowest node number go on the left
       uses simple sorting for a pair.  To handle multifurcations, must put in proper sorting
       based on code in imamp  3_9_09
       works recursively  """

    slengths = [0]* (2*numpops-1)
    firstint = [0] * (2*numpops - 1)
    holdsubs = [[]] * (2*numpops - 1)
    periodi = [0] * (2*numpops - 1)
    pos = 1
    subpos = pos
    subcount = 0
    pcount = 0
    slengths[subcount] = 0
    while 1:
        if substr[pos] == '(':
            pcount += 1
        if substr[pos] == ')':
            pcount-= 1
        pos+= 1
        slengths[subcount]+= 1
        if (pcount == 0):
            if (slengths[subcount] > 1):
                pos+= 1
                i = int(substr[pos])
                if pos <= len(substr)-2 and substr[pos+1].isdigit():
                    ts = substr[pos] + substr[pos+1]
                    i = int(ts)
                else:
                    i = int(substr[pos])
                periodi[subcount] = i
                if (i >= 10):
                    pos += 2
                    slengths[subcount] += 3
                else:
                    pos+= 1
                    slengths[subcount] += 2
            else:
                periodi[subcount] = -1
            holdsubs[subcount] = substr[subpos:pos]
            (holdsubs[subcount],hold) = set0(holdsubs[subcount],slengths[subcount])
            i = 0
            while (holdsubs[subcount][i].isdigit() == False):
                i+= 1
            firstint[subcount] = int(holdsubs[subcount][i])
            subcount+= 1
            slengths[subcount] = 0
            if (substr[pos] == ','):
                pos+= 1
            subpos = pos
        if pos >= len(substr):
            break
    if ((periodi[0] > periodi[1] and periodi[0] >= 0 and periodi[1] >= 0) or (periodi[0] >= 0 and periodi[1] < 0)):
        substr = strlistadd(substr,0,'(')
        j = slengths[1]
        k = 0
        i = 1
        while i<= j:
            substr = strlistadd(substr,i,holdsubs[1][k])
            k += 1
            i += 1
        subpos = 1
        if (slengths[1] > 2):
            (substr,hold) = set0(substr,i)
            substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
            substr = joinlist(substr,hold)
        substr = strlistadd(substr,i,',')
        i+= 1
        subpos = i
        j += 1 + slengths[0]
        k = 0
        while i <= j:
            substr = strlistadd(substr,i,holdsubs[0][k])
            i += 1
            k += 1
        if (slengths[0] > 2):
            (substr,hold) = set0(substr,i)
            substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
            substr = joinlist(substr,hold)
        substr = strlistadd(substr,i,')')
    else:
        if (firstint[0] > firstint[1] and periodi[0] < 0 and periodi[1] < 0):
            substr = strlistadd(substr,0,'(')
            j = slengths[1]
            k = 0
            i = 1
            while  i<= j:
                substr = strlistadd(substr,i,holdsubs[1][k])
                k += 1
                i += 1
            subpos = 1
            if (slengths[1] > 2):
                substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
            substr = strlistadd(substr,i,',')
            i+= 1
            subpos = i
            j += 1 + slengths[0]
            k = 0
            while i <= j:
                substr = strlistadd(substr,i,holdsubs[0][k])
                i+= 1
                k+= 1
            if (slengths[0] > 2):
                substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
            substr = strlistadd(substr,i,')')
        else:
            substr = strlistadd(substr,0,'(')
            subpos = 1
            if (slengths[0] > 2):
                (substr,hold) = set0(substr,slengths[0] + 1)
                substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
                substr = joinlist(substr,hold)
            substr = strlistadd(substr,slengths[0] + 1,',')
            subpos = slengths[0] + 2
            if (slengths[1] > 2):
                (substr,hold) = set0(substr,slengths[0] + slengths[1] + 2)
                substr[subpos:len(substr)] = rewrite (substr[subpos:len(substr)])
                substr = joinlist(substr,hold)
            substr = strlistadd(substr,slengths[0] + slengths[1] + 2,')')
    return substr


def plistbyperiod(poptreestring,poptree):
    """ generate a list, for each period this a list of the populations in that period,
         by their number in order from from left to right as they appear in the plot"""
    plist = [[]]
    for i in range(1,len(poptreestring)):
        if poptreestring[i-1] != ":" and (not poptreestring[i-1].isdigit()) and poptreestring[i].isdigit():
            plist[0].append(int(poptreestring[i]))
    droppops = [[-1,-1]]
    addpop = [-1]
    numtreepops = 2*numpops - 1
    for pi in range(1,numpops):
        droppops.append([])
        k=0
        for j in range(numtreepops):
            if poptree[j][1] == pi:
                droppops[pi].append(j)
                k += 1
                if k > 2:
                    print ( "droppop problem ")
                    break
            if poptree[j][0] == pi:
                addpop.append(j)
        tplist1 = plist[pi-1]
        tplist2 = []
        added = False
        j = 0
        while j < len(tplist1):
            if tplist1[j] == droppops[pi][0] or tplist1[j] == droppops[pi][1]:
                if added == False:
                    added = True
                    tplist2.append(addpop[pi])
            else:
                tplist2.append(tplist1[j])
            j += 1
        plist.append(tplist2)
    return plist,droppops,addpop

def poptreeread (poptreestring):
    """ copy of the function in imamp
         use a list of lists to hld poptree
         poptree[i] is the info for population [i]
         poptree[i][0] is the period number population [i] starts in
         poptree[i][1] is the period it ends in
         poptree[i][2] is the left up pop
         poptree[i][3] is the right up pop
         poptree[i][4] is the downpop
         examples
        (poptree,rootpop,poptreestring,plist, droppops,addpop) = poptreeread("(((5,6):12,7):13,(4,((3,1):9,(2,0):8):10):11):14",8)
        (poptree,rootpop,poptreestring,plist) = poptreeread("(4,((3,1):6,(2,0):5):7):8",5)
         """

    poptree = []
    for i in range(numpops):
        poptree.append([-1,-1,-1,-1,-1])
        poptree[i][0] = 0
    numtreepops = 2*numpops - 1
    for i in range(numpops, numtreepops):
        poptree.append([-1,-1,-1,-1,-1])
    poptreelist = []
    for i in range(len(poptreestring)):
        poptreelist.append(poptreestring[i])
    poptreelist = rewrite(poptreelist)
    newpoptreestring = ''
    for i in range(len(poptreelist)):
        newpoptreestring += poptreelist[i]
    stringspot = 0
    ancestralpopnums = []
    for i in range(2*numpops - 1):
        ancestralpopnums.append(0)
    (poptree, ancestralpopnums) = parenth0(numpops,poptree,newpoptreestring,stringspot,ancestralpopnums)
    (poptree,rootpop,stringspot,periodi,nextnode) = parenth(numpops,poptree,newpoptreestring,stringspot,ancestralpopnums,-1,-1,0)
    (plist, droppops,addpop) = plistbyperiod(newpoptreestring,poptree)
    return poptree, rootpop, newpoptreestring, plist, droppops,addpop

def meanrgb(color1,color2):
    """
        generates a 'mean' color
        check_colormath is global
        if colormath is available  this will do take the average in LAB space
            see colormath documentation
        otherwise there is a crude function for averaging rgb values
    """
    if check_colormath:
        srgb1 = sRGBColor(color1[0],color1[1],color1[2])
        srgb2 = sRGBColor(color2[0],color2[1],color2[2])

        lab1 = convert_color (srgb1,LabColor)
        lab2 = convert_color (srgb2,LabColor)
        lab1tuple = SpectralColor.get_value_tuple(lab1)
        lab2tuple = SpectralColor.get_value_tuple(lab2)
        labAtuple = ( (lab1tuple[0] + lab2tuple[0])/2.0 , (lab1tuple[1] + lab2tuple[1])/2.0,
                (lab1tuple[2] + lab2tuple[2])/2.0 )
        labA = LabColor(labAtuple[0],labAtuple[1],labAtuple[2])
        rgbA = convert_color(labA,sRGBColor)
        rgbAtuple = SpectralColor.get_value_tuple(rgbA)
        return list(rgbAtuple)
    else:
        acolor = [0,0,0]
        for j in range(3):
            ## this seems to give a useful average color
            meancolor = (color1[j] + color2[j])/2.0
            # now lighten it a bit
            acolor[j] = (1.0 - (0.8 * (1.0 -meancolor )))
        return acolor


def addcolors(poptree):
    """
        add colors to poptree.  set ancestors to average of descendant populations
    """
    rgbset = ([[0.4,0.4,0.4],
                [0.650980392,0.462745098,0.11372549],
                [0.4,0.650980392,0.117647059],
                [0.121568627,0.470588235,0.705882353],
                [0.905882353,0.160784314,0.541176471],
                [0.458823529,0.439215686,0.701960784],
                [0.850980392,0.37254902,0.007843137],
                [0.105882353,0.619607843,0.466666667],
                [0.901960784,0.670588235,0.007843137],
                [0.984313725,0.603921569,0.6]])
    for i in range(numpops):
        poptree[i].append(rgbset[i])
    for i in range(numpops,2*numpops-1):
        poptree[i].append([])
    while True:
        notdone = False
        for i in range(numpops,2*numpops-1):
            if poptree[i][5] == []:
                ld = poptree[i][2]
                rd = poptree[i][3]
                if poptree[ld][5] != [] and poptree[rd][5] != []:
                    acolor = [0,0,0]
                    acolor = meanrgb(poptree[ld][5],poptree[rd][5])
                    poptree[i][5] = acolor
                else:
                    notdone = True
        if notdone == False:
            break
    return poptree


def yline(y,farright, width, dash, grayamount):
    """ draw a line at a specific height in relative terms """
    aline([[0,y],[farright,y]],width, dash, grayamount)


def centerbox(pop,leftpoint,rightpoint,poptree,popxvals):
    """
        centerbox is a recursive function to find locations of population boxes on the x axis
        pop is the population for which we are finding the left and right sides of the box

        centerbox returns:
            width - the difference between right and left side of population it was called with
            center - the center of the population it was called with
            new value of popxvals
                popxvals[pop] is set,  they start out with values of 0 and box width, but get new values
            leftpoint
                leftpoint is the point less than which we cannot find values because a box can't be drawn there
                leftpoint is partly determined by the descendant population and partly by the population
            rightpoint
                the rightmost point of any descendant population

        Start at the bottom, go up the left side then the right, recursively
        take the width of the left side and the width of the right side,
        put a spacer between, add them, and find the center
        the width and center gets returned to the basal population

        receives the population #,  leftmost side of the ancestor population, rightside,   the tree and popxvals (which this modifies)

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

    if poptree[pop][2] == -1:  ## pop is a terminal population
        ## at this point popxvals[pop] holds just the width of the box (i.e. popxvals[pop][0] is 0)
        popxvals[pop][1] = popxvals[pop][1] - popxvals[pop][0] + leftpoint
        popxvals[pop][0] = leftpoint
        return popxvals[pop][1]-popxvals[pop][0],leftpoint + (popxvals[pop][1]-popxvals[pop][0])/ 2.0,popxvals, leftpoint, popxvals[pop][1]
    else:
        popspacer = gv["popboxspaceadj"] * popboxspacedefault
        (lw,lc, popxvals, leftpoint,rightpoint) = centerbox(poptree[pop][2],leftpoint,rightpoint, poptree,popxvals)
        rleftpoint = rightpoint + popspacer
        (rwidth,rcenter, popxvals, rleftpoint,rightpoint) = centerbox(poptree[pop][3],rleftpoint,rightpoint, poptree,popxvals)
        newwidth = lw + popspacer + rwidth

        newwidth = popxvals[pop][1] - popxvals[pop][0]
        newcenter = lc + (rcenter - lc)/2.0
        if newcenter - (newwidth/2.0) < leftpoint :
            newcenter += leftpoint - (newcenter - (newwidth/2.0))
        templeft = newcenter - newwidth/2.0
        popxvals[pop][0] = templeft
        popxvals[pop][1]  = templeft + newwidth
        return newwidth, newcenter,popxvals,leftpoint,rightpoint

def fround(val):
    """ crude rounding to a couple decimal points for a positive val"""
    if val==0:
        return "0.0"
    lval = math.log10(val)
    if lval < 0:
        lval -= 1
    rval = -int(lval) + 2
    if lval > 3:
        return str(int(round(val, rval)))
    return str(round(val, rval))

def popadjustx(popxvals,minx_popbox):
    " shift box locations to the left, as needed to fit with minx_popbox"
    minx = popxvals[0][0]
    for i in range(1,len(popxvals)):
        if minx > popxvals[i][0]:
            minx = popxvals[i][0]
    for i in range(len(popxvals)):
        popxvals[i][0] -= (minx - minx_popbox)
        popxvals[i][1] -= (minx - minx_popbox)
    return popxvals

def setpopbox(ty,slist,scaledtime,rootpop,poptree):
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
    for i in range(numpops-1):
        wadjust += "00"
    if(scaledtime != []):
        minx_popbox = textwide(wadjust+"0.00 MYR", tfactor)
    else:
        minx_popbox = textwide(wadjust+"0.00 tu", tfactor)
    minx_popbox /= gv["globalscale"]
    if gv["localxscale"] > 0:
        minx_popbox /= gv["localxscale"]

    popxvals = []
## if scaledpop == [] then no text is written on time split line and there is more width to work with
    for i in range(2*numpops - 1):
## left side temporarily at zero, right side temporarily at upper confidence interval
        popxvals.append( [0,slist[4][4][i][1]])
    (width,c,popxvals, leftpoint,rightpoint) = centerbox(rootpop,0,popxvals[rootpop][1],poptree,popxvals)
    popxvals = popadjustx(popxvals,minx_popbox)
    popbox = []

    ## maxwide will be used to adjust the width as a scaler  so the part furthest to the right is not too far out
    maxwide = 0
    for i in range(2*numpops-1):
        if maxwide < (popxvals[i][1] + (slist[4][4][i][3]-slist[4][4][i][1])):
            maxwide = (popxvals[i][1] + (slist[4][4][i][3]-slist[4][4][i][1]))
    maxwide = maxwide/(1.0-minx_popbox)

    if gv["localxscale"] > 0:
        maxwide *= gv["localxscale"]

    farright = 0
    confint = []
    for i in range(2*numpops-1):
        confint.append([])
        confint[i].append(minx_popbox + ((popxvals[i][1] - (slist[4][4][i][1]-slist[4][4][i][2]))/maxwide))
        confint[i].append(minx_popbox + ((popxvals[i][1] + (slist[4][4][i][3]-slist[4][4][i][1]))/maxwide))
        if confint[i][1] > farright:
            farright = confint[i][1]
        popbox.append([[],[]])
        popbox[i][0].append(minx_popbox + popxvals[i][0]/maxwide)
        popbox[i][1].append(minx_popbox + popxvals[i][1]/maxwide)
        if poptree[i][1] == -1:
            popbox[i][0].append(gv["lineINFy"])
        else:
            popbox[i][0].append(ty[poptree[i][1]-1][0])
        if poptree[i][0] == 0:
            popbox[i][1].append(gv["line0y"])
        else:
            popbox[i][1].append(ty[poptree[i][0]-1][0])
    return popbox,maxwide, confint, farright

def printpopbox(popbox,maxwide,confint,slist,plist,rootpop, poptree, ty,scaledpop,droppops):
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
    if gv["simplecolor"]:
        color = gv["blue"]
    else:
        color = gv["black"]
    cdim = []
    for i in range(2*numpops-1):
        tempbox = [row[:] for row in popbox[i]] # copy 2d list
        tempbox[1][0] = popbox[i][1][0] - (slist[4][4][i][1]-slist[4][4][i][2])/maxwide
        # cdim.append(curvebox(-1,tempbox,1.5,color,gv["graylevel"],i,gv["dashinterval"],poptree))
        cdim.append(calccdim(-1,tempbox))
        w("%%begin box %d" % i)
        cdimtemp = curvebox(cdim[i],popbox[i],2.5,color,0,i,0,poptree)
        w("%%done box %d" % i)
    cdim = []
    if gv["popboxcintervalboxes"]: ## print confidence interval boxes
        for i in range(2*numpops-1):
            tempbox = [row[:] for row in popbox[i]] # copy 2d list
            tempbox[1][0] = popbox[i][1][0] - (slist[4][4][i][1]-slist[4][4][i][2])/maxwide
            w("%%begin left confidence for box %d" % i)
            cdim.append(calccdim(-1,tempbox))
            cdimtemp = curvebox(cdim[i],tempbox,1.5,color,gv["graylevel"],i,gv["dashinterval"],poptree)
            w("%%done left confidence for box %d" % i)
            tempbox = [row[:] for row in popbox[i]] # copy 2d list
            tempbox[1][0] = popbox[i][1][0] + (slist[4][4][i][3]-slist[4][4][i][1])/maxwide
            w("%%begin right confidence for box %d" % i)
            cdimtemp =curvebox(cdim[i],tempbox,1.5,color, gv["graylevel"],i,gv["dashinterval"],poptree)
            w("%%done right confidence for box %d" % i)
    popprintinc = 0.01
    if gv["usealtnames"]:
        namelist = gv["altpopnames"]
    else:
        namelist = slist[2][4]
    if gv["anglenames"]:
        angle = 30
    else:
        angle = 0
    for i in range(numpops): ## population names are in slist[2][4]
        if poptree[i][1] == 1 and i== droppops[1][1]: ## right side of most recent split
            dotext([popbox[i][0][0] + (popbox[i][1][0]-popbox[i][0][0])/2,popbox[i][1][1]+popprintinc],namelist[i],angle, False)
        else:
            dotext([popbox[i][0][0],popbox[i][1][1]+popprintinc],namelist[i],angle, False)
    popprintinc = 0.025
    if gv["label_a_pops"]:
        for i in range(numpops,2*numpops-1):
            dotext([max(popbox[i][0][0],popbox[i][0][0] + (popbox[i][1][0] - popbox[i][0][0])/2.0 - popprintinc),\
                popbox[i][0][1] + (popbox[i][1][1] - popbox[i][0][1])/2.0],"pop #"+str(i),0, False)
## plot the confidence arrows for population boxes
    lastperiod = [0]*(2*numpops-1)
    for i in range(2*numpops-1):
        for j in range(len(plist)):
            for k in range(len(plist[j])):
                if plist[j][k] == i and j > lastperiod[i]:
                    lastperiod[i] = j
    periodposcount = [0]*numpops
    arrowheightinc = 0.006
    arrowheights = []
    for i in range(numpops):
        if i==0:
            top = gv["line0y"]
            bot = ty[i][0]
        else:
            top = ty[i-1][0]
            if i== numpops - 1:
                bot = gv["lineINFy"]
            else:
                bot = ty[i][0]

        if top-bot < 0.1:
            frac = 0.5
        else:
            frac = 0.8
        arrowheights.append(top - (top-bot)*frac)
    if gv["popboxcintervalarrows"]: ## print confidence interval arrows
        for i in range(2*numpops-1):
            period = lastperiod[i]
            arrowheight = max(popbox[i][0][1],arrowheights[period] -periodposcount[period]*2*arrowheightinc)
            head = [confint[i][0],arrowheight]
            tail = [popbox[i][1][0],arrowheight]
            # head is tip of arrow to lower bound of confidence interval
            # if there is not room for the arrowhead, don't print it
            if tail[0] - head[0] > arrowheadwidthdefault:
                arrow(head,tail,2,color)
            head = [confint[i][1],arrowheight]
            tail = [popbox[i][1][0],arrowheight]
            arrow(head,tail,0, color)
            periodposcount[period] += 1
    if scaledpop != []:
        ane = scaledpop[rootpop]/1000
        anes = fround(ane)
        dotext([0.15,0.05]," Ancestral Ne (thousands): " + anes,0, False)
    else:
        dotext([0.15,0.05]," Ancestral 4Nu: " + str(slist[4][4][rootpop][1]),0, False)

    if gv["simplecolor"]:
        w("0 0 0  setrgbcolor")
    return popbox

def set_tlines(ty,slist):
    """
        line0y - default relative height of time 0
        eventimes - if True,   space split times evenly
        lastt_lower_y - height of oldest split time,  by default is 1/(numpops+1),   else can be set by user
    """
    t = []
    for i in range(numpops-1):
        t.append([slist[5][4][i][1],slist[5][4][i][2],slist[5][4][i][3]])  ## [time,  upper ci,  lower ci]
    ty = []
    if gv["localyscale"] == -1:
        yint = gv["line0y"] - gv["lastt_lower_y"]
        for i in range(numpops-1):
            ty.append([])
            if gv["eventimes"] ==  False:
                tmax = slist[5][4][numpops-2][3] ## bottom of confidence interval of largest(oldest) t
                for j in range(3):
                    ty[i].append(gv["line0y"] - (t[i][j]*yint)/tmax)
            else:
##                ty[i].append(gv["line0y"] - ((i+1)/float(numpops+1)*yint)/tmax)
                ty[i].append(gv["line0y"] - yint * (i+1)/float(numpops) )
    else:
        timeumean = slist[7][4][1]
        scaleumean = slist[7][4][2]
        for i in range(numpops-1):
            ty.append([])
            for j in range(3):
                ty[i].append(gv["line0y"] - (t[i][j] * (scaleumean/timeumean/1e6)* gv["localyscale"]))
                if ty[i][j] < gv["lineINFy"]:
                    print ( " time line too low in graph,  reduce local y scale (-y value) ")
        gv["lastt_lower_y"] = ty[numpops-2][2]
##    print "ty : ",ty
    return ty


def print_tlines(ty,slist,scaledtime, farright):
    """
        print the split time lines and confidence interval lines
        graylevel is global
    """
    xinc = 0.005
    yinc = 0.002
    if(scaledtime != []):
        if max(scaledtime)/1e6  < 1.0:
            yearscaler = 1e3
            yearscalestring = " KYR"
        else:
            yearscaler = 1e6
            yearscalestring = " MYR"
    if gv["eventimes"] == False:
        for i in range(numpops-1):
            if (ty[i][1] > ty[i][0]):
                yline(ty[i][1],farright,1,2,gv["graylevel"])
            yline(ty[i][0],farright,0.5,0,0)
            if (ty[i][2] < ty[i][0]):
                yline(ty[i][2],farright,1,2,gv["graylevel"])
            if(scaledtime != []):
                scaledtime[i] /= yearscaler
                mtime = round(scaledtime[i],-int(math.log10(scaledtime[i])-2))
                nstr = str(mtime) + yearscalestring
    ##            str(int(round(scaledtime[i],-int(math.log10(scaledtime[i])-2)))) + " yrs"
                dotext([xinc*(i+2),ty[i][0]+yinc],nstr,0, False)
            else:
                nstr = fround(slist[5][4][i][1]) + "tu"
                dotext([xinc*(i+2),ty[i][0]+yinc],nstr,0, False)
            if (ty[i][1] > ty[i][0]):
                arrow([xinc*(i+1),ty[i][1]],[xinc*(i+1),ty[i][0]],1, gv["black"])
            if (ty[i][2] < ty[i][0]):
                arrow([xinc*(i+1),ty[i][2]],[xinc*(i+1),ty[i][0]],3, gv["black"])
    else:
        for i in range(numpops-1):
            yline(ty[i][0],farright,0.5,0,0)
            if(scaledtime != []):
                scaledtime[i] /= yearscaler
                mtime = round(scaledtime[i],-int(math.log10(scaledtime[i])-2))
                nstr = str(mtime) + yearscalestring
    ##            str(int(round(scaledtime[i],-int(math.log10(scaledtime[i])-2)))) + " yrs"
                dotext([xinc*(i+2),ty[i][0]+yinc],nstr,0, False)
            else:
                nstr = fround(slist[5][4][i][1]) + "tu"
                dotext([xinc*(i+2),ty[i][0]+yinc],nstr,0, False)
    return ty

def print_mcurves(slist, popbox, plist):
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
    def checkm(val2NM, llr):
##        return  (gv["moption"] == 'a' and val2NM > min2NM) or  (gv["moption"] == 's' and llr >= 2.74)  or  val2NM > gv["moption"]
##  can happen that singificant m  is associated with 2NM = 0 (or near).  Catch these cases and do not print
        minsig2NM = 0.001
        if type(gv["moption"]) is float:
            return val2NM >= gv["moption"]
        else:
            returnval =  ((gv["moption"] == 'a' and val2NM >= min2NM)
                    or  (gv["moption"] == 's' and llr.count('*') > 0 and val2NM >=minsig2NM )
                    or  (gv["moption"] == 'S' and llr.count('*') > 1) and val2NM >=minsig2NM )
            return  returnval
    if gv["moption"] == 'x':
        return
    mperiodnum = [0]*(numpops-1)
    if len(slist[6]) > 4:
        sml = slist[6][4]
        miginfo = []
        mi = 0
        for i in range(len(sml)):
##            pratio =  sml[i][3]/sml[i][2]
##            llr = 2*math.log(pratio)
## alternate code to get values from Marginal peak location tables
            llr = sml[i][2]

            usem = checkm(sml[i][1],llr)
            if usem:
                miginfo.append([])
                c1 = max(sml[i][0].find("M"),sml[i][0].find("m")) ## either upper of lower case
                c2 = sml[i][0].find(">")
                miginfo[mi].append(int(sml[i][0][c2+1:len(sml[i][0])]))
                miginfo[mi].append(int(sml[i][0][c1+1:c2]))
                pos1 = -1
                pos2 = -1
                p  = 0
                while 1:
                    for j in range(len(plist[p])):
                        if plist[p][j] == miginfo[mi][0]:
                            pos1 = j
                        if  plist[p][j] == miginfo[mi][1]:
                            pos2 = j
                    if pos1 >= 0 and pos2 >=0:
                        if pos1 < pos2:
                            direction = 0
                        else:
                            direction = 2
                        break
                    else:
                        p += 1
                        pos1 = -1
                        pos2 = -1
                miginfo[mi].append(direction)
                miginfo[mi].append(p)
                miginfo[mi].append(mperiodnum[p])
                mperiodnum[p] += 1
                miginfo[mi].append(sml[i][1])
                miginfo[mi].append(llr)
                mi += 1
            else:
                continue

        mboxfrac = 0.3
        ## set height of curves
        y = []
        for i in range(len(miginfo)):
            frompop = miginfo[i][0]
            period = miginfo[i][3]
            hi = popbox[frompop][1][1]
            for j in range (len(plist[period])):
                if hi > popbox[plist[period][j]][1][1]:
                    hi = popbox[plist[period][j]][1][1]
            lo = 0
            for j in range (len(plist[period])):
                if lo < popbox[plist[period][j]][0][1]:
                    lo = popbox[plist[period][j]][0][1]
            y.append(hi - (hi - lo)*(miginfo[i][4]+1)/(mperiodnum[miginfo[i][3]]+1))
        for i in range(len(miginfo)):
            frompop = miginfo[i][0]
            topop = miginfo[i][1]
            period = miginfo[i][3]
            direc = miginfo[i][2]
            val2NM = fround(miginfo[i][5])
            if (miginfo[i][6] != 'ns'):
                val2NM += miginfo[i][6]
            text2NMwidth = textwide(val2NM,2.5)
            if direc == 0:
                tailx =  popbox[frompop][1][0] - (popbox[frompop][1][0]-popbox[frompop][0][0])*mboxfrac
                headx =  popbox[topop][0][0] + (popbox[topop][1][0] - popbox[topop][0][0]) * mboxfrac
                if (text2NMwidth > abs(tailx-headx)):
                    tailx -= (text2NMwidth - abs(tailx-headx))/2
                    headx += (text2NMwidth - abs(tailx-headx))/2
            if direc == 2:
                tailx =  popbox[frompop][0][0] + (popbox[frompop][1][0] - popbox[frompop][0][0]) * mboxfrac
                headx =  popbox[topop][1][0] - (popbox[topop][1][0]-popbox[topop][0][0])* mboxfrac
                if (text2NMwidth > abs(tailx-headx)):
                    tailx += (text2NMwidth - abs(tailx-headx))/2
                    headx -= (text2NMwidth - abs(tailx-headx))/2
            migrationcurvearrow(val2NM,[headx,y[i]],[tailx,y[i]],direc,gv["red"])
            if gv["rgbcolor"]:
                migrationcurvearrow(val2NM,[headx,y[i]],[tailx,y[i]],direc,gv["darkgreen"])
            else:
                migrationcurvearrow(val2NM,[headx,y[i]],[tailx,y[i]],direc,gv["red"])


##***********************************************************************************
##////////////// Command line use ///////////////////////////////////////////////////
##***********************************************************************************

def scancommandline(args):
    """ command line consists of flags, each with a dash, '-', followed immediately by a letter
        some flags should be followed by a value, depending on the flag.  The value can be placed
        immediately after the flag or spaces can be inserted """
    global gv
    def aflag ():
        gv["label_a_pops"] = True
    def bflag (tempval):
        gv["popboxspaceadj"] = float(tempval)

    def dflag ():
        gv["skipdemographicscaling"]  =  True
    def eflag():
        gv["eventimes"] = True
    def iflag (tempname):
        gv["imfilename"] = tempname.strip()
    def oflag (tempname):
        gv["outputfilename"]= tempname.strip()
    def gflag (tempval):
        gv["globalscale"] = float(tempval)
    def xflag (tempval):
        ##  edited 9/1/2017,  this seemed to work better.  use maximumxpoint for making plot wider,  and use localxscale for makeing it narrower
        f = float(tempval)
        if f > 1.0:
            gv["maximumxpoint"]  = gv["maximumxpoint"] * f
        else:
            gv["localxscale"] = f
    def yflag (tempval):
        gv["localyscale"] = float(tempval)
    def hflag (tempval):
        gv["arrowheightadj"] = float(tempval)
    def fflag(tempval):
        gv["font"] = tempval
        gv["bifont"] = gv["font"] + "-BoldItalic"
    def kflag ():
        gv["line0y"] = 0.88 ## a tradeof, between need to make room and not wanting to squash figure
        gv["anglenames"] = True
    def mflag(tempval):
        if tempval[0].isdigit():
            gv["moption"] = float(tempval)
        else:
            gv["moption"] = tempval
    def nflag (tempname):
        gv["usealtnames"] = True
        gv["altnamefilename"] = tempname.strip()
    def qflag ():
        gv["popboxcintervalboxes"]  =  False
    def rflag ():
        gv["popboxcintervalarrows"]  =  False
    def pflag(tempval):
        gv["fontsize"] = float(tempval)
        gv["fontfixed"] = True
    def tflag(tempval):
        gv["lastt_lower_y"] = float(tempval)
        gv["set_lastt_lower_y"] = False
    def sflag ():
        gv["dosquare"] = True
        gv["maximumxpoint"] = 576.1
    def uflag ():
        gv["simplecolor"]  = True
    def vflag ():
        gv["rgbcolor"] = True

    def removewhitespace(temps):
        return "".join(temps.split())

    def cleanarglist(arglist,flags_with_values):
        """
        """
        newarg = []
        if arglist[0][0] != "-":  # skip program name at beginning of list
            arglist = arglist[1:]
        ai  = 0
        while ai < len(arglist):
            if removewhitespace(arglist[ai]) != "":
                arglist[ai] = removewhitespace(arglist[ai])
            else:
                print ( "bad whitespace in command line: ",repr(" ".join(arglist)))
                sys.exit(1)
            if arglist[ai][0] == '-':
                if arglist[ai][1] in flags_with_values  and len(arglist[ai])==2:  ## found a space in the command line
                    arglist[ai] = arglist[ai] + arglist[ai+1]
                    newarg.append(arglist[ai])
                    ai += 1
                else:
                    newarg.append(arglist[ai])
            else:
                print ( "error on command line,  \"-\" not found:",arglist[ai])
                printcommandset()
                sys.exit(1)
            ai += 1

        return newarg

    def checkallflags(flags_with_values,flags_withoutvalues,cldic):
        """
            checks that flags that must be used are used
            checks that flags_with_values,flags_withoutvalues and cldic all make use of the appropriate flags
        """
        if len(set(flags_with_values).intersection(set(flags_without_values))) > 0:
            print ( "error some flags appear in two lists of flags,  with and without required values:",set(flags_with_values).intersection(set(flags_without_values)))
            printcommandset()
            sys.exit(1)
        for flag in set(flags_with_values).union(set(flags_withoutvalues)):
            if flag not in cldic:
                print ( "error some flag mismatch between strings of flags and dictionary of flags:",flag)
                printcommandset()
                sys.exit(1)
        return
    def check_flags_used(flagsused, flags_must_use):
        for f in flags_must_use:
            if f not in flagsused:
                print("-%c missing from command line. Run without any commands to get the help screen."%f)
                sys.exit(1)
        return

    cldic = {'a':aflag,'b':bflag,'d':dflag,'e':eflag,'f':fflag,\
             'g':gflag,'h':hflag,'i':iflag,'k':kflag,'m':mflag,'n':nflag,'o':oflag,\
             'p':pflag, 'q':qflag,'r':rflag,'s':sflag, 't':tflag,'u':uflag,'v':vflag,\
             'x':xflag,'y':yflag}
    flags_must_use = 'i'
    flags_with_values =  "bfghimoptxyn"
    flags_without_values = "adesuvkqr"
    cmdstr = " ".join(args)
    checkallflags(flags_with_values,flags_without_values,cldic)
    argv = cleanarglist(args,flags_with_values)
    flagsused = ''
    for i in range(0,len(argv)):
        if argv[i][0] == '-':
            flaglet = argv[i][1].lower()
            flagsused += flaglet
##            print ( i, flaglet)
            if len(argv[i]) == 2:
                if i == (len(argv)-1):
                    cldic[flaglet]()
                else:
                    if  argv[i+1][0] == '-':
                        cldic[flaglet]()
                    else:
                        cldic[flaglet](argv[i+1])
                        i += 1
            else:
                if (len(argv[i]) < 2):
                    print ( "problem on command line ")
                    exit()
                cldic[flaglet](argv[i][2:len(argv[i])])
        else:
            print ( "error on command line,  \"-\" not found:",argv[i])
            printcommandset()
            sys.exit(1)
    check_flags_used(flagsused, flags_must_use)
    return cmdstr

def printcommandset():
    print ( "IMfig command line terms  (-i is required):")
    print ( "help or h causes this screen to be written")
    print ( "-a : include ancestral population #'s in plot")
    print ( "-b : adjust width spacing of population boxes, values > 0, default = 1")
    print ( "-d : do not use demographic scale information even if in input file")
    print ( "-e : space split times evenly  (not proportional to time,  no confidence intervals shown)")
    print ( "-f : font.  Default=Arial. Use postscript fonts available on the computer")
    print ( "     e.g. Arial, Helvetica, Times-roman, Courier")
    print ( "-g : global plot scale sets the size of the plot, max = 1, default = 1")
    print ( "-h : arrow width, default = 1")
    print ( "-i : input file name")
    print ( "-k : print population names on an angle")
    print ( "-m : options for printing of arrows and 2Nm values for migration :")
    print ( "      -m x :  do not print migration arrows")
    print ( "      -m a : 2Nm migration arrows for all 2NM > 0")
    print ( "      -m s : 2Nm migration arrows only if m is statistically significant p <= 0.05 (default)")
    print ( "      -m S : 2Nm migration arrows only if m is statistically significant p <= 0.01")
    print ( "      -m # : '#' is a number, migration arrows appear when 2NM >= # (e.g. -m0.1)")
    print ( "-n : file with alternative species names")
    print ( "-o : output file name (with 'eps' extension) e.g. -o myoutputfile.eps, default=im_eps_file.eps")
    print ( "-p : fontsize (default is 14 for full scale, default follows global scale)")
    print ( "-q : no confidence interval boxes for population boxes printed")
    print ( "-r : no confidence interval arrows for population boxes printed")
    print ( "-s : print square, rather than landscape")
    print ( "-t : relative height of oldest time point, values between 0 and 1")
    print ( "     default value = 1/(# sampled populations+1)")
    print ( "-u : simple colors, blue for population boxes, red arrows for migration (default grayscale)")
    print ( "-v : multiple colors for population boxes, red arrows for migration (default grayscale)")
    print ( "-x : adjust width of plot,  >1 means wider, <1 means narrower")
    print ( "-y : adjust height of splittimes, relative to bottom of figure, max = 1.")  ## not clear what this does  5/12/2016



##*************************************************************
##///////////// default values, basic scale////////////////////
##*************************************************************


def setbasexyscale():
    global gv
    minimumxpoint = minimumypoint = 36.1   #forgot where this came from
    gv["fixedLL"] = [minimumxpoint,minimumypoint]
    gv["fixedUR"] = [gv["maximumxpoint"],gv["maximumypoint"]]

def setdefaults():
    global gv
    gv = {"label_a_pops":False,"simplecolor":False,"dosquare":False,
        "eventimes":False,"popboxcintervalboxes":True,"popboxcintervalarrows":True,"imfilename":"im_eps.txt","outputfilename":"im_eps_file.eps","globalscale":1,
         "font":"Arial","bifont":"Arial-BoldItalic","fontsize":14,"fontfixed":False,
         "line0y":0.95,"lineINFy":0.1,"localxscale":-1,"localyscale":-1,"arrowheightadj":1,
        "maximumxpoint":756.1,"maximumypoint":576.1,"lastt_lower_y":-1,"set_lastt_lower_y":True,
        "blue":[0,0,1],"red":[1,0,0],"black":[0,0,0],"darkgreen":[0,0.58823,0.19607],"graylevel":0.3,
        "popboxspaceadj":1.0,"moption":'s',"skipdemographicscaling":False,"rgbcolor":False,
        "anglenames":False,"dashinterval":3,"usealtnames":False,"altnamefilename":"","altpopnames":[]
         }


def dostuff(args):
    global gv
    setdefaults()

    ##////////////// get info from the command line ///////////////////
    cmdstr = scancommandline(args)

    ##////////////// get info from the input file (i.e. the IM results files) ///////////////////
    print ( "input file: %s\noutput file %s" % (gv["imfilename"], gv["outputfilename"]))
    print ( "read inputfile")
    (slist,scaledpop,scaledtime) = readimfile()

    ##////////////// read the tree, set up plist ///////////////////
    ## plist has the population numbers in each period, in order from left to right as they appear in the plot
    (poptree,rootpop,poptreestring,plist, droppops,addpop) = poptreeread(slist[3][4])
    if gv["rgbcolor"]:
        poptree = addcolors(poptree)

    ##//////////////// set scales ///////////////////////////
    print ( "set scales")
    setbasexyscale()
    if gv["set_lastt_lower_y"]:
        gv["lastt_lower_y"] = 1.0/(numpops + 1)
    ty = set_tlines(slist, slist)
    (popbox,maxwide,confint, farright) = setpopbox(ty,slist,scaledtime,rootpop,poptree)

    ##//////////////// write the output file ///////////////////////////
    gv["epsf"] = open(gv["outputfilename"],"w")
    w("%!PS-Adobe-3.0 EPSF-3.0")
    w("%%legal size in landscape is 792x612 set bounding box with 0.5inch margins")
    w("%%the lower corner is at 36 36, x dim is 720 wide,  y dim is 540 hi")
    w("%%%%BoundingBox: %d %d  %d  %d" % (int(gv["fixedLL"][0]),int(gv["fixedLL"][1]),int(gv["fixedUR"][0]),int(gv["fixedUR"][1])))
    w("%%%%IMfig program author: Jody Hey   Copyright 2009-2016")
    w("%%%%Command line for IMfig program that generated this file: %s"%cmdstr)

    #### useful for debugging, include this DrawAnX function in the code
    ##w("/DrawAnX")
    ##w("{ 3 3 rmoveto -6 -6 rlineto")
    ##w("0 6 rmoveto 6 -6 rlineto")
    ##w("0.01 setlinewidth")
    ##w("stroke } def")
    #### use by calling and passing x and y values
    ####e.g. w("%f %f moveto DrawAnX" %(point[0],point[1]))

    print ( "make figure")
    print ( "splitting times")
    ty = print_tlines(ty,slist, scaledtime, farright)
    print ( "population boxes")
    popbox = printpopbox(popbox,maxwide,confint,slist,plist,rootpop, poptree, ty, scaledpop,droppops)
    print ( "migration arrows")
    print_mcurves(slist,popbox, plist)
    gv["epsf"].close()
    print ( "plot completed")
    return


##***********************************************************************************
##////////////// MAIN PROGRAM ///////////////////////////////////////////////////////
##***********************************************************************************

##This program can be run from the command line (see user manual)
##Alternatively, if desired the user can uncomment this block and define the working directory and 'cmdstr':

## -- comment this block out when running from the command line --
##import os
##os.chdir(r"E:\genemod\ML-MCMC\SEAI\IMa3_work\IMfig\testbedwork")  ## directory where the input file is located
##cmdstr = r"IMfig.py -iBaHzSwYr_200loci_mode3ghost_j1_mp5.out -odebugBaHzSwYr_200loci_mode3ghost_j1_mp5.eps -ms -v -e"
##sys.argv = cmdstr.split()
## -----------

print ( "IMfig program. Copyright 2009-2018  Jody Hey  Release Date %s"%releasedate)
if len(sys.argv) <= 1 or sys.argv[1].upper() == "HELP" or sys.argv[1].upper() == 'H' or sys.argv[1].upper() == "-HELP" or sys.argv[1].upper() == '-H':  ## no arguments or only help or h
    printcommandset()
    sys.exit()
else:
    dostuff(sys.argv[1:])
    sys.exit()

# when debugging can define cmdstr (or edit args in launch.json)
# dostuff(cmdstr.split()[1:])  ## run using cmdstr
# sys.exit()
