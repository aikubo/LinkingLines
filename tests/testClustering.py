# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Tue Feb  1 16:27:58 2022

# @author: akh
# """

import numpy as np
import pandas as pd
import linkinglines as ll
import pytest

def testValidity(lines, dtheta, drho):
    
    if any(lines['ThetaRange'] > dtheta):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")
        
    if any(lines['RhoRange'] > drho):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")

def Run3Times(dikeset, dtheta, drho, **kwargs):
    xc1,yc1=ll.HT_center(dikeset)
    # run the first time 
    dikeset, Z1=ll.AggCluster(dikeset, dtheta, drho, **kwargs)
    lines1, IC1=ll.examineClusters(dikeset)
    
    # rotate 45 deg 
    dikeset2=ll.rotateData(dikeset, 45)
    dikeset2, xc2, yc2=ll.HT(dikeset2)
    dikeset2, Z2=ll.AggCluster(dikeset2, dtheta, drho, **kwargs)
    dikeset2=ll.rotateData(dikeset2,-45)
    lines2, __=ll.examineClusters(dikeset2)
    
    #move to lower left
    dx=np.min([dikeset['Xstart'].min(), dikeset['Xend'].min()])
    dy=np.min([dikeset['Ystart'].min(), dikeset['Yend'].min()])
    dikeset3, xc3, yc3=ll.HT(dikeset, xc=dx,yc=dy)

    dikeset3, Z3=ll.AggCluster(dikeset3, dtheta, drho, **kwargs)
    
    lines3, IC3=ll.examineClusters(dikeset3)
    
    #check the changes 
    eq12, diff12=ll.checkClusterChange(lines1, lines2)
    eq13, diff13=ll.checkClusterChange(lines1, lines3)
    eq32, diff32=ll.checkClusterChange(lines3, lines2)
    
    Flines1=ll.FilterLines(lines1)
    Flines2=ll.FilterLines(lines2)
    Flines3=ll.FilterLines(lines3)
    
    print("Comparing the filtered lines")
    eq12, diff12=ll.checkClusterChange(Flines1, Flines2)
    eq13, diff13=ll.checkClusterChange(Flines1, Flines3)
    eq32, diff32=ll.checkClusterChange(Flines3, Flines2)
    
        
    return lines1, lines2, lines3

def RotateAndCluster(dikeset,dtheta, drho,**kwargs):
    dikesetOrig, ZOrig=ll.AggCluster(dikeset, dtheta, drho)
    
    dikeset45=ll.rotateData(dikeset, 45)
    dikeset45=ll.DikesetReProcess(dikeset45, HTredo=True)
    dikeset45, Z45=ll.AggCluster(dikeset45, dtheta, drho)
    linesOrig=ll.examineClusterShort(dikesetOrig)
    lines45=ll.examineClusterShort(dikeset45)
    
    changedect=ll.checkAllClusterChange(linesOrig, lines45)
    if changedect:
        return dikesetOrig
    if ~changedect:
        eqLabels, diffLabels=ll.checkIndividualClusterChange(linesOrig, lines45)
        dikesetFinal=dikesetOrig.iloc[eqLabels[:len(dikesetOrig)]]

        return dikesetFinal
    
