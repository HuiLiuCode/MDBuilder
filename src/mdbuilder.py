# ----------------------------------------------------------------------
# Copyright (c) 2016-2017
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------

from __future__ import print_function, division

import os
import sys
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from Tkinter import *

import tkMessageBox
import tkFileDialog
try:
    from math import fsum as sum
except ImportError:
    pass

import Pmw
from pymol import cmd, util
from pymol.cgo import *

try:
    import numpy
except ImportError:
    _HAS_LIB = 0
else:
    _HAS_LIB = 1

#import traceback
#import getopt
import time
import math
import operator
from copy import deepcopy
from struct import calcsize, unpack
from collections import OrderedDict, defaultdict
from math import sqrt, cos, sin, radians
from itertools import combinations, product, repeat, chain
from timeit import default_timer as clock
from getpass import getuser

try:
    from __builtin__ import xrange as range
except ImportError:
    pass
try:
    from functools import reduce
except ImportError:
    pass

__home_dir__ = os.path.split(__file__)[0]
__program__ = 'MDBuilder'
__version__ = '1.0'
__author__ = 'Hui Liu, Ye Jin'
__email__ = 'teffliu@hotmail.com, jinyecpu@163.com'
__url__ = 'https://github.com/huiliucode/mdbuilder'
__desc__ = 'a program for the preparation of biomolecular simulations'
__doc__ = """
The PyMOL plugin helps to setup MDBuilder runs and provides visual support.

MDBuilder is a program for the preparation of biomolecular simulations.

Do your job as follows,

1) Make sure that the MDBuilder package is installed correctly.
2) Launch the plugin.
3) Set up the configuration in the "I/O" page.
4) Set up the configuration in the "Preparation" page.
    Click the "Execute" button and wait.
5) Set up the configuration in the "Solvation" page.
    Click the "Execute" button and wait.
6) Set up the configuration in the "Ionization" page.
    Click the "Execute" button and wait.
7) Click the "Output" button to generate the files.

URL: %s
Author: %s <%s>
""" % (__url__, __author__, __email__)

_debug = False

IONS = {
    'Na+': 'SOD',
    'K+': 'POT',
    'Mg2+': 'MG',
    'Ca2+': 'CAL',
    'Zn2+': 'ZN2',
    'Cl-': 'CLA'
    }

def __init__(self):
    """Register function for the plugin."""
    self.menuBar.addmenuitem('Plugin', 'command',
                             '%s %s'%(__program__, __version__),
                             label='%s %s'%(__program__, __version__),
                             command=lambda x=self: MDBuilderGui(x))


def perr(*errMsg):
    print("Error:", *errMsg, file=sys.stderr)
    #sys.exit(1)


def wopen(filename, mode='wb'):
    assert mode[0] == 'w'
    if os.path.isfile(filename):
        copy = filename + '.copy'
        if os.path.isfile(copy):
            try:
                os.remove(copy)
                print("Romove", copy)
            except OSError as err:
                print(err)
                perr("Cannot remove", copy)
        try:
            os.rename(filename, copy)
            print("Rename", filename, "to", copy)
        except OSError as err:
            print(err)
            perr("Cannot rename", filename)

    return open(filename, mode)


def read_pdb(filename, epreInpDict=None):
    aliasAtomList = epreInpDict['ALIASATOM'] if epreInpDict is not None else []
    aliasResiList = epreInpDict['ALIASRES'] if epreInpDict is not None else []
    pdbDataList = []
    newSegList = []
    pdbResCnt = 0
    lastResNum = 0
    lastResName = None

    try:
        pdbFile = open(filename)
    except IOError:
        perr("Cannot open", filename)

    print("Read PDB file", filename)

    fmt = '6x 5s x 4s x 4s s 4s 4x 8s 8s 8s 13x 8s'
    size = calcsize(fmt)

    for line in pdbFile:
        pdbHead = line[:6].rstrip()

        if pdbHead in ['ATOM', 'HETATM']:
            (pdbAtomNum, pdbAtomName, pdbResName,
             pdbChainID, pdbResNum, pdbPosX, pdbPosY,
             pdbPosZ, pdbSegName) = unpack(fmt, line[:size])

            pdbAtomNum = int(pdbAtomNum)
            pdbAtomName = pdbAtomName.strip()
            pdbResName = pdbResName.strip()
            pdbResNum = int(pdbResNum)
            pdbPosX = float(pdbPosX)
            pdbPosY = float(pdbPosY)
            pdbPosZ = float(pdbPosZ)
            pdbSegName = pdbSegName.strip()

            if aliasAtomList:
                for aliasAtom in aliasAtomList:
                    if (aliasAtom[0][1] == pdbAtomName and
                                aliasAtom[0][0] in [pdbResName, '']):
                        pdbAtomName = aliasAtom[1]

            newAtom = [pdbAtomNum, pdbAtomName, pdbPosX, pdbPosY, pdbPosZ]

            if aliasResiList:
                for aliasResi in aliasResiList:
                    if pdbResName == aliasResi[0]:
                        pdbResName = aliasResi[1]

            if pdbResNum == lastResNum and pdbResName == lastResName:
                newSegList[pdbResCnt - 1][2].append(newAtom)
            else:
                newRes = [pdbResNum, pdbResName, [newAtom]]
                newSegList.append(newRes)
                pdbResCnt += 1

            lastResNum = pdbResNum
            lastResName = pdbResName

        elif not pdbHead or pdbHead in ['TER', 'END']:
            if newSegList:
                pdbDataList.append([None, newSegList])
            newSegList = []
            pdbResCnt = 0
            lastResNum = 0
            lastResName = None

    if newSegList:
        pdbDataList.append([None, newSegList])

    inpSegList = epreInpDict['SEGMENT'] if epreInpDict is not None else []
    inpSegNum = len(inpSegList)
    pdbSegNum = len(pdbDataList)

    if inpSegList and inpSegNum == pdbSegNum:
        for i, seg in enumerate(pdbDataList):
            seg[0] = inpSegList[i]
    else:
        print("Use default segment names: 'S1' for segment 1")
        for i, seg in enumerate(pdbDataList):
            seg[0] = 'S%d' % (i + 1)

    pdbFile.close()

    global _debug
    if _debug:
        print(pdbDataList)

    return pdbDataList


def write_pdb(filename, pdbDataList, boxInfo=None):
    pdbFile = wopen(filename)

    print("Write PDB file", filename)

    date = time.strftime("%d-%b-%y", time.localtime()).upper()
    pdbFile.write('HEADER    %-40s%-30s\n' % ('GENERATED BY MDBuilder 1.0', date))
    pdbFile.write('AUTHOR    %-70s\n' % getuser().upper())

    if boxInfo is not None:
        a, b, c = boxInfo[2]
        pdbFile.write('%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d%10s\n' %
                      ('CRYST1', a, b, c, 90.0, 90.0, 90.0, 'P 1', 1, ' '))

    fmt = '%-6s%5d %-4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s  \n'

    atomCnt = 0
    for segName, segDataList in pdbDataList:
        chainID = segName[0]

        for resNum, resName, resAtom in segDataList:

            for atomNum, atomName, atomPosX, atomPosY, atomPosZ in resAtom:
                atomCnt += 1
                atomNum = atomCnt
                pdbFile.write(fmt % ('ATOM', atomNum, atomName.ljust(3).rjust(4),
                                     resName, chainID, resNum, atomPosX, atomPosY,
                                     atomPosZ, 1.0, 0.0, segName, atomName[0]))

    nxform = 0 if boxInfo is None else 1
    natom = atomCnt
    pdbFile.write('%6s%4s%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%10s\n' %
                  ('MASTER', ' ', 0, 0, 0, 0, 0, 0, 0, nxform, natom, 0, 0, 0, ' '))
    pdbFile.write('%-80s\n' % 'END')

    pdbFile.close()


def write_tmppdb(filename, pdbDataList, boxInfo=None):
    if hasattr(filename, 'seek'):      # file object
        if hasattr(filename, 'name'):
            fname = os.path.basename(filename.name)
        else:
            fname = ''
    if fname:
        print("ERROR")
    else:
        print("Write PDB file from string")

    date = time.strftime("%d-%b-%y", time.localtime()).upper()
    filename.write('HEADER    %-40s%-30s\n' % ('GENERATED BY MDBuilder 1.0', date))
    filename.write('AUTHOR    %-70s\n' % getuser().upper())

    if boxInfo is not None:
        a, b, c = boxInfo[2]
        filename.write('%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d%10s\n' %
                      ('CRYST1', a, b, c, 90.0, 90.0, 90.0, 'P 1', 1, ' '))

    fmt = '%-6s%5d %-4s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s  \n'

    atomCnt = 0
    for segName, segDataList in pdbDataList:
        chainID = segName[0]

        for resNum, resName, resAtom in segDataList:

            for atomNum, atomName, atomPosX, atomPosY, atomPosZ in resAtom:
                atomCnt += 1
                atomNum = atomCnt
                filename.write(fmt % ('ATOM', atomNum, atomName.ljust(3).rjust(4),
                                     resName, chainID, resNum, atomPosX, atomPosY,
                                     atomPosZ, 1.0, 0.0, segName, atomName[0]))

    nxform = 0 if boxInfo is None else 1
    natom = atomCnt
    filename.write('%6s%4s%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%10s\n' %
                  ('MASTER', ' ', 0, 0, 0, 0, 0, 0, 0, nxform, natom, 0, 0, 0, ' '))
    filename.write('%-80s\n' % 'END')


def download_pdb(filename):
    """Download a structure file from RCSB Protein Date Bank."""
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen
        from urllib.error import URLError

    filename = filename.upper()
    if filename.endswith('.PDB'):
        pdbCode = filename[:-4]
    else:
        pdbCode = filename
    if len(pdbCode) != 4:
        raise ValueError("PDB ID should be a 4-character string")
    pdbFile = pdbCode + '.pdb'

    print("Download", pdbFile, "from Protein Data Bank")
    sys.stdout.flush()

    try:
        pdb = urlopen('http://www.rcsb.org/pdb/download/downloadFile.do?' +
                      'fileFormat=pdb&compression=NO&structureId=' + pdbCode)
    except URLError:
        #perr("Cannot download the file")
        print("Cannot download the file")

    wopen(pdbFile, 'wb').write(pdb.read())

    if os.path.getsize(pdbFile) < 1024:
        #perr("Invalid pdb code", pdbCode)
        raise PDBDownloadError('Download is failed')
    print('Download is successfully completed')


def read_charmm_top_title(ctopLines):
    ctopTitleList = []

    for i, line in enumerate(ctopLines):
        ctopVer = '0 0'

        if line[0] == '*':
            ctopTitle = line
            ctopTitleList.append(ctopTitle)
        elif not line.strip():
            continue
        elif line[:4].strip().upper() in ['MASS', 'END']:
            break
        elif ctopLines[i - 1][0] == '*':
            ctopVer = ' '.join(line.strip().split())
            break

    ctopTitleList.append(ctopVer)

    return ctopTitleList


def read_charmm_top_mass(ctopLines):
    ctopAtomMassDict = {}

    for line in ctopLines:

        if line.startswith('MASS'):
            splitLine = line.split()
            ctopAtomType = splitLine[2]
            ctopAtomMass = float(splitLine[3])
            ctopAtomMassDict[ctopAtomType] = ctopAtomMass
        elif line.startswith(('DECL', 'RESI', 'PRES', 'END')):
            break

    return ctopAtomMassDict


def read_charmm_top_defa(ctopLines):
    ctopDefaDict = {}

    for line in ctopLines:

        if line.startswith('DEFA'):
            splitLine = line.split()
            for i in range(1, len(splitLine), 2):
                ctopPatch = splitLine[i][:4]
                ctopPatchName = splitLine[i + 1]
                if ctopPatch in ['FIRS', 'LAST']:
                    ctopDefaDict[ctopPatch] = ctopPatchName
                else:
                    perr("Unknown argument", ctopPatch)
        elif line.startswith(('RESI', 'PRES', 'END')):
            break

    return ctopDefaDict


def read_charmm_top_resi(ctopLines, ctopDefaDict, noCmap):
    ctopResiDataDict = {}
    lastResiName = None

    for line in ctopLines:
        ctopHead = line[:4]

        if ctopHead == 'RESI':
            splitLine = line.split()
            ctopResiName = splitLine[1]
            lastResiName = ctopResiName
            ctopResiDataDict[ctopResiName] = [OrderedDict(), [], [], [], [], {}]
            ctopResiPatchDict = ctopResiDataDict[ctopResiName][5]
            ctopResiPatchDict['FIRS'] = ctopDefaDict['FIRS']
            ctopResiPatchDict['LAST'] = ctopDefaDict['LAST']
        elif ctopHead == 'ATOM':
            if lastResiName is None: continue
            (ctopResiAtomName, ctopResiAtomType,
             ctopResiAtomChrg) = line.split()[1:]
            ctopResiAtomChrg = float(ctopResiAtomChrg)
            ctopResiAtomDict = ctopResiDataDict[lastResiName][0]
            ctopResiAtomDict[ctopResiAtomName] = [ctopResiAtomType,
                                                  ctopResiAtomChrg]
        elif ctopHead in ['BOND', 'DOUB', 'TRIP', 'AROM']:
            if lastResiName is None: continue
            splitLine = line.split()
            ctopResiBondList = ctopResiDataDict[lastResiName][1]
            for atomA in range(1, len(splitLine) - 1, 2):
                atomB = atomA + 1
                newBond = [splitLine[atomA], splitLine[atomB]]
                ctopResiBondList.append(newBond)
        elif ctopHead == 'IMPR':
            if lastResiName is None: continue
            splitLine = line.split()
            for atomA in range(1, len(splitLine) - 3, 4):
                atomB = atomA + 1
                atomC = atomA + 2
                atomD = atomA + 3
                newImpr = [splitLine[atomA], splitLine[atomB],
                           splitLine[atomC], splitLine[atomD]]
                ctopResiImprList = ctopResiDataDict[lastResiName][2]
                ctopResiImprList.append(newImpr)
        elif ctopHead == 'CMAP':
            if noCmap or not lastResiName: continue
            splitLine = line.split()
            ctopResiCmapList = ctopResiDataDict[lastResiName][3]
            ctopResiCmapList.append(splitLine[1:9])
        elif line[:2] == 'IC':
            if lastResiName is None: continue
            splitLine = line.split()
            ic = splitLine[1:]
            ic[4] = float(ic[4])
            ic[5] = radians(float(ic[5]))
            ic[6] = radians(float(ic[6]))
            ic[7] = radians(float(ic[7]))
            ic[8] = float(ic[8])
            ctopResiICList = ctopResiDataDict[lastResiName][4]
            ctopResiICList.append(ic)
        elif ctopHead == 'PATC':
            if lastResiName is None: continue
            splitLine = line.split()
            splitLineLen = len(splitLine)
            ctopResiPatchDict = ctopResiDataDict[lastResiName][5]
            if splitLineLen == 3:
                firstPatch = splitLine[1][:4]
                firstPres = splitLine[2]
                if firstPatch in ctopResiPatchDict:
                    ctopResiPatchDict[firstPatch] = firstPres
                else:
                    perr("Unknown argument", firstPatch)
            elif splitLineLen == 5:
                firstPatch = splitLine[1][:4]
                firstPres = splitLine[2]
                secondPatch = splitLine[3][:4]
                secondPres = splitLine[4]
                if firstPatch in ctopResiPatchDict:
                    ctopResiPatchDict[firstPatch] = firstPres
                else:
                    perr("Unknown argument", firstPatch)
                if secondPatch in ctopResiPatchDict:
                    ctopResiPatchDict[secondPatch] = secondPres
                else:
                    perr("Unknown argument", secondPatch)
            else:
                perr("Unknown patching format")
        elif ctopHead == 'PRES' or ctopHead.strip() == 'END':
            break

    return ctopResiDataDict


def read_charmm_top_pres(ctopLines, ctopDefaDict, noCmap):
    ctopPresDataDict = {}
    lastPresName = None

    for line in ctopLines:
        ctopHead = line[:4]

        if ctopHead == 'PRES':
            splitLine = line.split()
            ctopPresName = splitLine[1]
            lastPresName = ctopPresName
            ctopPresDataDict[ctopPresName] = [[], OrderedDict(), [], [], [], []]
        elif ctopHead == 'DELE':
            if lastPresName is None: continue
            splitLine = line.split()
            if splitLine[1] == 'ATOM':
                ctopPresDelAtomName = splitLine[2]
                ctopPresDelAtomList = ctopPresDataDict[lastPresName][0]
                ctopPresDelAtomList.append(ctopPresDelAtomName)
        elif ctopHead == 'ATOM':
            if lastPresName is None: continue
            (ctopPresAtomName, ctopPresAtomType,
             ctopPresAtomChrg) = line.split()[1:]
            ctopPresAtomChrg = float(ctopPresAtomChrg)
            ctopPresAtomDict = ctopPresDataDict[lastPresName][1]
            ctopPresAtomDict[ctopPresAtomName] = [ctopPresAtomType,
                                                  ctopPresAtomChrg]
        elif ctopHead in ['BOND', 'DOUB', 'TRIP', 'AROM']:
            if lastPresName is None: continue
            splitLine = line.split()
            ctopPresBondList = ctopPresDataDict[lastPresName][2]
            for atomA in range(1, len(splitLine) - 1, 2):
                atomB = atomA + 1
                newBond = [splitLine[atomA], splitLine[atomB]]
                ctopPresBondList.append(newBond)
        elif ctopHead == 'IMPR':
            if lastPresName is None: continue
            splitLine = line.split()
            for atomA in range(1, len(splitLine) - 3, 4):
                atomB = atomA + 1
                atomC = atomA + 2
                atomD = atomA + 3
                newImpr = [splitLine[atomA], splitLine[atomB],
                           splitLine[atomC], splitLine[atomD]]
                ctopPresImprList = ctopPresDataDict[lastPresName][3]
                ctopPresImprList.append(newImpr)
        elif ctopHead == 'CMAP':
            if noCmap or not lastPresName: continue
            splitLine = line.split()
            ctopPresCmapList = ctopPresDataDict[lastPresName][4]
            ctopPresCmapList.append(splitLine[1:9])
        elif line[:2] == 'IC':
            if lastPresName is None: continue
            splitLine = line.split()
            ic = splitLine[1:]
            ic[4] = float(ic[4])
            ic[5] = radians(float(ic[5]))
            ic[6] = radians(float(ic[6]))
            ic[7] = radians(float(ic[7]))
            ic[8] = float(ic[8])
            ctopPresICList = ctopPresDataDict[lastPresName][5]
            ctopPresICList.append(ic)
        elif ctopHead.strip() == 'END':
            break

    return ctopPresDataDict


def read_charmm_top(filename, epreInpDict):
    noCmap = epreInpDict['NOCMAP']

    try:
        ctopFile = open(filename)
    except IOError:
        perr("Cannot open", filename)

    ctopLines = ctopFile.read()

    print("Read CHARMM topology file", filename, end="\n\n")

    ctopTitleLines = ctopLines.splitlines()
    ctopTitleList = read_charmm_top_title(ctopTitleLines)
    for ctopTitle in ctopTitleList[:-1]:
        print(ctopTitle)
    ctopVer = ctopTitleList[-1]
    print('Developed under CHARMM version ' + ctopVer)

    ctopLines = ctopLines.upper()

    massStart = ctopLines.find('\nMASS')
    declStart = ctopLines.find('\nDECL')
    resiStart = ctopLines.find('\nRESI')
    presStart = ctopLines.find('\nPRES')

    if massStart != -1:
        hasMass = True
        if declStart != -1:
            ctopMassLines = ctopLines[massStart:declStart]
        elif resiStart != -1:
            ctopMassLines = ctopLines[massStart:resiStart]
        elif presStart != -1:
            ctopMassLines = ctopLines[massStart:presStart]
        else:
            ctopMassLines = ctopLines[massStart:]
        ctopMassLines = [line.split('!')[0] for line in ctopMassLines.
            splitlines() if line]
    else:
        hasMass = False

    if declStart != -1:
        hasDecl = True
        if resiStart != -1:
            ctopDeclLines = ctopLines[declStart:resiStart]
        elif presStart != -1:
            ctopDeclLines = ctopLines[declStart:presStart]
        else:
            ctopDeclLines = ctopLines[declStart:]
        ctopDeclLines = [line.split('!')[0] for line in ctopDeclLines.
            splitlines() if line]
    else:
        hasDecl = False

    if resiStart != -1:
        hasResi = True
        if presStart != -1:
            ctopResiLines = ctopLines[resiStart:presStart]
        else:
            ctopResiLines = ctopLines[resiStart:]
        ctopResiLines = [line.split('!')[0] for line in ctopResiLines.
            splitlines() if line]
    else:
        hasResi = False

    if presStart != -1:
        hasPres = True
        ctopPresLines = [line.split('!')[0] for line in ctopLines[resiStart:].
            splitlines() if line]
    else:
        hasPres = False

    if hasMass:
        print("\nLoad mass data")
        ctopAtomMassDict = read_charmm_top_mass(ctopMassLines)
    else:
        ctopAtomMassDict = {}

    if hasDecl:
        print("Load decl data")
        ctopDefaDict = read_charmm_top_defa(ctopDeclLines)
    else:
        ctopDefaDict = {}

    if hasResi:
        print("Load resi data")
        ctopResiDataDict = read_charmm_top_resi(ctopResiLines,
                                                ctopDefaDict, noCmap)
    else:
        ctopResiDataDict = {}

    if hasPres:
        print("Load pres data")
        ctopPresDataDict = read_charmm_top_pres(ctopPresLines,
                                                ctopDefaDict, noCmap)
    else:
        ctopPresDataDict = {}

    global _debug
    if _debug:
        print(ctopAtomMassDict)
        print(ctopResiDataDict)
        print(ctopPresDataDict)

    return ctopTitleList, ctopAtomMassDict, ctopResiDataDict, ctopPresDataDict


def read_charmm_prm(filename):
    try:
        cprmFile = open(filename)
    except IOError:
        perr("Cannot open", filename)

    print("Read CHARMM parameter file", filename)

    cprmBondDict = {}
    cprmAnglDict = {}
    cprmUBDict = {}
    cprmDiheDict = {}
    cprmImprDict = {}
    cprmCmapDict = {}
    cprmNbndDict = {}
    cprmNb14Dict = {}

    while True:
        line = cprmFile.next()
        if line.startswith('BOND'): break

    while True:
        line = cprmFile.next()
        if line.startswith('ANGLES'): break
        splitLine = line.split('!')[0].split()
        if splitLine:
            cprmBondDict[(splitLine[0], splitLine[1])] = (float(splitLine[2]), float(splitLine[3]))
            cprmBondDict[(splitLine[1], splitLine[0])] = (float(splitLine[2]), float(splitLine[3]))

    while True:
        line = cprmFile.next()
        if line.startswith('DIHEDRALS'): break
        splitLine = line.split('!')[0].split()
        if splitLine:
            cprmAnglDict[(splitLine[0], splitLine[1], splitLine[2])] = (
            float(splitLine[3]), radians(float(splitLine[4])))
            cprmAnglDict[(splitLine[2], splitLine[1], splitLine[0])] = (
            float(splitLine[3]), radians(float(splitLine[4])))
        if len(splitLine) == 7:
            cprmUBDict[(splitLine[0], splitLine[1], splitLine[2])] = (float(splitLine[5]), float(splitLine[6]))
            cprmUBDict[(splitLine[2], splitLine[1], splitLine[0])] = (float(splitLine[5]), float(splitLine[6]))

    while True:
        line = cprmFile.next()
        if line.startswith('IMPROPER'): break
        splitLine = line.split('!')[0].split()
        if splitLine:
            key = (splitLine[0], splitLine[1], splitLine[2], splitLine[3])
            if key in cprmDiheDict:
                cprmDiheDict[key].append((float(splitLine[4]), int(splitLine[5]), radians(float(splitLine[6]))))
            else:
                cprmDiheDict[key] = [(float(splitLine[4]), int(splitLine[5]), radians(float(splitLine[6])))]

    while True:
        line = cprmFile.next()
        if line.startswith('CMAP') or line.startswith('NONBONDED'): break
        splitLine = line.split('!')[0].split()
        if splitLine:
            cprmImprDict[(splitLine[0], splitLine[1], splitLine[2], splitLine[3])] = float(splitLine[4]), float(
                splitLine[6])

    if line.startswith('CMAP'):
        while True:
            line = cprmFile.next()
            if line.startswith('NONBONDED'): break
            splitLine = line.split('!')[0].split()
            if splitLine and len(splitLine) == 9:
                key = tuple(splitLine[:8])
                resolution = int(splitLine[8])
                cprmFile.next()
                val = []
                for i in range(resolution):
                    aval = []
                    i = 0
                    while i < int(math.ceil(resolution / 5)):
                        line = cprmFile.next()
                        splitLine = line.split('!')[0].split()
                        if splitLine:
                            aval.extend(map(float, splitLine))
                            i += 1
                    cprmFile.next()
                    assert len(aval) == 24
                    val.append(aval)
                assert len(val) == 24
                cprmCmapDict[key] = val

    line = cprmFile.next()
    while True:
        line = cprmFile.next()
        if line.startswith('HBOND') or line.startswith('END'): break
        splitLine = line.split('!')[0].split()
        if splitLine:
            cprmNbndDict[splitLine[0]] = (float(splitLine[2]), float(splitLine[3]))
            if len(splitLine) == 7:
                cprmNb14Dict[splitLine[0]] = (float(splitLine[5]), float(splitLine[6]))
            else:
                cprmNb14Dict[splitLine[0]] = (float(splitLine[2]), float(splitLine[3]))

    cprmFile.close()

    global _debug
    if _debug:
        print(cprmBondDict)
        print(cprmAnglDict)

    return cprmBondDict, \
           cprmAnglDict, \
           cprmUBDict, \
           cprmDiheDict, \
           cprmImprDict, \
           cprmCmapDict, \
           cprmNbndDict, \
           cprmNb14Dict


def read_rename_rule(filename, epreInpDict=None):
    try:
        ruleFile = open(filename).read().splitlines()
    except IOError:
        perr("Cannot open", filename)

    print("Read rename rule file", filename)

    i = 0
    while i < len(ruleFile):
        splitLine = ruleFile[i].split('#')[0].split()
        if not splitLine:
            i += 1
            continue
        splitLineLen = len(splitLine)
        arg = splitLine[0].upper()

        if arg == 'ALIASRES':
            if splitLineLen == 3:
                epreInpDict['ALIASRES'].append([splitLine[1], splitLine[2]])
            else:
                perr("Missing or too many arguments", splitLine[0])
        elif arg == 'ALIASATOM':
            if splitLineLen == 3:
                epreInpDict['ALIASATOM'].append([splitLine[1].split(':'),
                                                 splitLine[2]])
            else:
                perr("Missing or too many arguments", splitLine[0])
        else:
            perr("Unknown argument", splitLine[0])
        i += 1

    return epreInpDict


def read_bond_file(filename, epreInpDict=None):
    try:
        bondFile = open(filename).read().splitlines()
    except IOError:
        perr("Cannot open", filename)

    print("Read bond file", filename)

    i = 0
    while i < len(bondFile):
        splitLine = bondFile[i].split('#')[0].split()
        if not splitLine:
            i += 1
            continue
        splitLineLen = len(splitLine)
        arg = splitLine[0].upper()

        if arg == 'DISUBOND':
            epreInpDict['DISUBOND']['DODISU'] = True
            if splitLineLen == 3:
                epreInpDict['DISUBOND']['DISULIST'].append([
                    [int(tmp) for tmp in splitLine[1].split(':')],
                    [int(tmp) for tmp in splitLine[2].split(':')]])
            else:
                perr("Missing or too many arguments", splitLine[0])
        else:
            perr("Unknown argument", splitLine[0])
        i += 1

    return epreInpDict


def do_disu(epreInpDict, pdbDataList, ctopResiDataDict, ctopPresDataDict):
    doAuto = epreInpDict['DISUBOND']['AUTO']
    cut = epreInpDict['DISUBOND']['CUT']
    if doAuto:
        print("Auto-detect disulfide bonds")
        cysList = []
        for i, (_, segDataList) in enumerate(pdbDataList):
            for resNum, resName, resAtom in segDataList:
                if resName == 'CYS':
                    for atom in resAtom:
                        if atom[1] == 'SG':
                            cysList.append([i + 1, resNum, atom[2:]])

        numCys = len(cysList)
        print("Found %d cys residues" % numCys)
        if numCys < 2:
            print("No cys pairs found.")
            return
        disuList = do_disu_detect(cysList, cut ** 2)
        print("Found %d possible disufide bonds with a %f A cutoff" %
              (len(disuList), cut))
        epreInpDict['DISUBOND']["DISULIST"] = disuList
        if not disuList: return
    else:
        disuList = epreInpDict['DISUBOND']["DISULIST"]

    for si, sj in disuList:

        for i, res in enumerate(pdbDataList[si[0] - 1][1]):
            if res[0] == si[1]:
                if res[1] == 'CYS':
                    pdbDataList[si[0] - 1][1][i][1] = '_EPRE_CYM1'
                    break
                else:
                    perr("Segment %d residue %d is not cys" % tuple(si))
        else:
            perr("Cannot find segment %d residue %d" % tuple(si))

        for i, res in enumerate(pdbDataList[sj[0] - 1][1]):
            if res[0] == sj[1]:
                if res[1] == 'CYS':
                    pdbDataList[sj[0] - 1][1][i][1] = '_EPRE_CYM2'
                    break
                else:
                    perr("Segment %d residue %d is not cys" % tuple(sj))
        else:
            perr("Cannot find segment %d residue %d" % tuple(sj))

    cysResi = deepcopy(ctopResiDataDict['CYS'])
    patchResi = 'DISU'
    patchList = ctopPresDataDict[patchResi]
    deleList = [i[1:] for i in patchList[0] if i[0] == '1']
    patomDict = OrderedDict()
    for i, j in patchList[1].items():
        if i[0] == '1': patomDict[i[1:]] = j
    atomDict = cysResi[0]
    bondList = cysResi[1]
    imprList = cysResi[2]
    cmapList = cysResi[3]
    icList = cysResi[4]
    atomDict.update(patomDict)
    for atom in deleList: del atomDict[atom]
    deleList = set(deleList)
    cysResi[0] = atomDict
    cysResi[1] = [bond for bond in bondList if not set(bond) & deleList]
    cysResi[2] = [impr for impr in imprList if not set(impr) & deleList]
    cysResi[3] = [cmap for cmap in cmapList if not set(cmap) & deleList]
    cysResi[4] = [ic for ic in icList if not set(ic[:4]) & deleList]
    ctopResiDataDict['_EPRE_CYM1'] = cysResi

    cysResi = deepcopy(ctopResiDataDict['CYS'])
    patchResi = 'DISU'
    patchList = ctopPresDataDict[patchResi]
    deleList = [i[1:] for i in patchList[0] if i[0] == '2']
    patomDict = OrderedDict()
    for i, j in patchList[1].items():
        if i[0] == '2': patomDict[i[1:]] = j
    atomDict = cysResi[0]
    bondList = cysResi[1]
    imprList = cysResi[2]
    cmapList = cysResi[3]
    icList = cysResi[4]
    atomDict.update(patomDict)
    for atom in deleList: del atomDict[atom]
    deleList = set(deleList)
    cysResi[0] = atomDict
    cysResi[1] = [bond for bond in bondList if not set(bond) & deleList]
    cysResi[2] = [impr for impr in imprList if not set(impr) & deleList]
    cysResi[3] = [cmap for cmap in cmapList if not set(cmap) & deleList]
    cysResi[4] = [ic for ic in icList if not set(ic[:4]) & deleList]
    ctopResiDataDict['_EPRE_CYM2'] = cysResi


def do_disu_detect(cysList, cut2):
    disuDict = defaultdict()
    for i, j in combinations(cysList, 2):
        pi = i[2]
        pj = j[2]
        bondLen2 = (pi[0] - pj[0]) ** 2 + (pi[1] - pj[1]) ** 2 + (pi[2] - pj[2]) ** 2
        if bondLen2 <= cut2:
            disuDict[(i[0], i[1])].append(((j[0], j[1]), sqrt(bondLen2)))

    disuList = []
    for key, val in disuDict.items():
        if len(val) > 1:
            perr("Found more than 1 disufide bonds both share a atom")
        disuList.append([[key[0], key[1]], [val[0][0], val[0][1]]])

    return disuList


def do_disu_connect(disuList, pdbDataList, bondList):
    for si, sj in disuList:

        for i, res in enumerate(pdbDataList[si[0] - 1][1]):
            if res[0] == si[1]:
                if res[1] == '_EPRE_CYM1':
                    for atom in res[2]:
                        if atom[1] == 'SG':
                            b1 = atom[0]
                            break
                    else:
                        perr("Cannot find SG atom")
                    break
                else:
                    perr("Segment %d residue %d is not cys" % tuple(si))
        else:
            perr("Cannot find segment %d residue %d" % tuple(si))

        for i, res in enumerate(pdbDataList[sj[0] - 1][1]):
            if res[0] == sj[1]:
                if res[1] == '_EPRE_CYM2':
                    for atom in res[2]:
                        if atom[1] == 'SG':
                            b2 = atom[0]
                            break
                    else:
                        perr("Cannot find SG atom")
                    break
                else:
                    perr("Segment %d residue %d is not cys" % tuple(sj))
        else:
            perr("Cannot find segment %d residue %d" % tuple(sj))

        bondList.append([b1, b2])

        print("Connect atoms %d and %d" % (b1, b2))


def do_disu_rename(pdbDataList, atomList):
    for _, segDataList in pdbDataList:
        for i, (_, resName, _) in enumerate(segDataList):
            if resName in ('_EPRE_CYM1', '_EPRE_CYM2'):
                segDataList[i][1] = 'CYS'

    for atom in atomList:
        if atom[3] in ('_EPRE_CYM1', '_EPRE_CYM2'):
            atom[3] = 'CYS'


def build_struct(ctopDataList, pdbDataList, epreInpDict, cprmDataList):
    ctopAtomMassDict, ctopResiDataDict, ctopPresDataDict = ctopDataList
    cprmBondDict, cprmAnglDict = cprmDataList[:2]

    doDisu = epreInpDict['DISUBOND']['DODISU']
    if doDisu:
        do_disu(epreInpDict, pdbDataList, ctopResiDataDict, ctopPresDataDict)

    newDataList = []
    atomCnt = 0
    for segCnt, [segName, segDataList] in enumerate(pdbDataList):
        newDataList.append([segName, []])

        firstResi = deepcopy(ctopResiDataDict[segDataList[0][1]])
        patchResi = firstResi[5]['FIRS']
        if patchResi != 'NONE':
            patchList = ctopPresDataDict[patchResi]
            deleList = patchList[0]
            patomDict = deepcopy(patchList[1])
            atomDict = firstResi[0]
            bondList = firstResi[1]
            imprList = firstResi[2]
            cmapList = firstResi[3]
            icList = firstResi[4]
            atomDict.update(patomDict)
            patomDict.update(atomDict)
            atomDict = patomDict
            for atom in deleList:
                del atomDict[atom]
            deleList = set(deleList)
            firstResi[0] = atomDict
            firstResi[1] = [bond for bond in bondList
                            if not set(bond).intersection(deleList)]
            firstResi[2] = [impr for impr in imprList
                            if not set(impr).intersection(deleList)]
            firstResi[3] = [cmap for cmap in cmapList
                            if not set(cmap).intersection(deleList)]
            firstResi[4] = [ic for ic in icList
                            if not set(ic[:4]).intersection(deleList)]
            for bond in reversed(patchList[2]):
                firstResi[1].insert(0, bond)
            for impr in reversed(patchList[3]):
                firstResi[2].insert(0, impr)
            for cmap in reversed(patchList[4]):
                firstResi[3].insert(0, cmap)
            for ic in reversed(patchList[5]):
                firstResi[4].insert(0, ic)
        ctopResiDataDict[segName + 'FIRST'] = firstResi

        lastResi = deepcopy(ctopResiDataDict[segDataList[-1][1]])
        patchResi = lastResi[5]['LAST']
        if patchResi != 'NONE':
            patchList = ctopPresDataDict[patchResi]
            deleList = patchList[0]
            patomDict = deepcopy(patchList[1])
            atomDict = lastResi[0]
            bondList = lastResi[1]
            imprList = lastResi[2]
            cmapList = lastResi[3]
            icList = lastResi[4]
            for atom in deleList:
                del atomDict[atom]
            atomDict.update(patomDict)
            deleList = set(deleList)
            lastResi[0] = atomDict
            lastResi[1] = [bond for bond in bondList
                           if not set(bond).intersection(deleList)]
            lastResi[2] = [impr for impr in imprList
                           if not set(impr).intersection(deleList)]
            lastResi[3] = [cmap for cmap in cmapList
                           if not set(cmap).intersection(deleList)]
            lastResi[4] = [ic for ic in icList
                           if not set(ic[:4]).intersection(deleList)]
            for bond in reversed(patchList[2]):
                lastResi[1].insert(0, bond)
            for impr in reversed(patchList[3]):
                lastResi[2].insert(0, impr)
            for cmap in reversed(patchList[4]):
                lastResi[3].insert(0, cmap)
            for ic in reversed(patchList[5]):
                lastResi[4].insert(0, ic)
        ctopResiDataDict[segName + 'LAST'] = lastResi

        singleResi = deepcopy(ctopResiDataDict[segName + 'FIRST'])
        patchResi = singleResi[5]['LAST']
        if patchResi != 'NONE':
            patchList = ctopPresDataDict[patchResi]
            deleList = patchList[0]
            patomDict = deepcopy(patchList[1])
            atomDict = singleResi[0]
            bondList = singleResi[1]
            imprList = singleResi[2]
            cmapList = singleResi[3]
            icList = singleResi[4]
            for atom in deleList:
                del atomDict[atom]
            atomDict.update(patomDict)
            deleList = set(deleList)
            singleResi[0] = atomDict
            singleResi[1] = [bond for bond in bondList
                             if not set(bond).intersection(deleList)]
            singleResi[2] = [impr for impr in imprList
                             if not set(impr).intersection(deleList)]
            singleResi[3] = [cmap for cmap in cmapList
                             if not set(cmap).intersection(deleList)]
            singleResi[4] = [ic for ic in icList
                             if not set(ic[:4]).intersection(deleList)]
            for bond in reversed(patchList[2]):
                singleResi[1].insert(0, bond)
            for impr in reversed(patchList[3]):
                singleResi[2].insert(0, impr)
            for cmap in reversed(patchList[4]):
                singleResi[3].insert(0, cmap)
            for ic in reversed(patchList[5]):
                singleResi[4].insert(0, ic)
        ctopResiDataDict[segName + 'SINGLE'] = singleResi

        nres = len(segDataList)
        last = nres - 1
        for resCnt, [resiNum, resiName, resiAtom] in enumerate(segDataList):
            newDataList[segCnt][1].append([resiNum, resiName, []])

            if nres == 1:
                atomDict = ctopResiDataDict[segName + 'SINGLE'][0]
            elif resCnt == 0:
                atomDict = ctopResiDataDict[segName + 'FIRST'][0]
            elif resCnt == last:
                atomDict = ctopResiDataDict[segName + 'LAST'][0]
            else:
                atomDict = ctopResiDataDict[resiName][0]

            for atomName in atomDict:
                atomCnt += 1
                for pdbAtomList in resiAtom:
                    if atomName == pdbAtomList[1]:
                        newAtom = pdbAtomList
                        newAtom[0] = atomCnt
                        break
                else:
                    newAtom = [atomCnt, atomName, None, None, None]
                newDataList[segCnt][1][resCnt][2].append(newAtom)

    pdbDataList = newDataList

    atomList = []
    segNameNumDict = {}
    for segCnt, [segName, segDataList] in enumerate(pdbDataList):
        segNameNumDict[segName] = {}
        resiNameNumDict = segNameNumDict[segName]

        nres = len(segDataList)
        last = nres - 1
        for resCnt, [resiNum, resiName, resiAtom] in enumerate(segDataList):
            resiNameNumDict[resiNum] = {}

            if nres == 1:
                resList = ctopResiDataDict[segName + 'SINGLE']
            elif resCnt == 0:
                resList = ctopResiDataDict[segName + 'FIRST']
            elif resCnt == last:
                resList = ctopResiDataDict[segName + 'LAST']
            else:
                resList = ctopResiDataDict[resiName]

            atomDict = resList[0]
            icList = resList[4]

            if len(resiAtom) == 3 and resiName in ('WAT', 'TIP3', 'HOH'):
                pos = []
                h1Cnt = h2Cnt = None
                for atomCnt, atom in enumerate(resiAtom):
                    if atom[1] == 'OH2':
                        pos = atom[2:]
                    elif atom[1] == 'H1':
                        h1Cnt = atomCnt
                        if atom[2] is not None: break
                    elif atom[1] == 'H2':
                        h2Cnt = atomCnt
                        if atom[2] is not None: break
                    else:
                        print("Unknown atom %s in water" % atom[1])
                else:
                    if pos and h1Cnt is not None and h2Cnt is not None:
                        predH1Pos, predH2Pos = fix_cryst_wat(pos)
                        (resiAtom[h1Cnt][2], resiAtom[h1Cnt][3],
                         resiAtom[h1Cnt][4]) = predH1Pos
                        (resiAtom[h2Cnt][2], resiAtom[h2Cnt][3],
                         resiAtom[h2Cnt][4]) = predH2Pos

            numMiss = lastNumMiss = 0
            missAtoms = []

            while True:
                for atomCnt, atom in enumerate(resiAtom):
                    atomNum, atomName, atomPosX, atomPosY, atomPosZ = atom

                    if lastNumMiss == 0:
                        atomType, atomChrg = atomDict[atomName]
                        atomMass = ctopAtomMassDict[atomType]
                        resiNameNumDict[resiNum][atomName] = atomNum
                        atomList.append([atomNum, segName, resiNum, resiName,
                                         atomName, atomType, atomChrg, atomMass])

                    if atomPosX is not None: continue

                    for ic in icList:
                        pos = []

                        if atomName == ic[3]:
                            icNames = ic[:4]
                            if icNames[2][0] == '*':
                                icNames[2] = icNames[2][1:]

                            for index, icName in enumerate(icNames[:3]):
                                if icName[0] == '-':
                                    if resCnt == 0: break
                                    icName = icName[1:]
                                    icNames[index] = icName
                                    xres = pdbDataList[segCnt][1][resCnt - 1][2]
                                elif icName[0] == '+':
                                    if resCnt == last: break
                                    icName = icName[1:]
                                    icNames[index] = icName
                                    xres = pdbDataList[segCnt][1][resCnt + 1][2]
                                else:
                                    xres = pdbDataList[segCnt][1][resCnt][2]

                                for [num, name, posX, posY, posZ] in xres:
                                    if icName == name and posX is not None:
                                        pos.append([posX, posY, posZ])
                                        break
                                else:
                                    break
                            else:
                                if ic[7] == 0.0:
                                    anglKey = tuple(map(lambda x: atomDict[x][0],
                                                        icNames[1:]))
                                    if anglKey in cprmAnglDict:
                                        ic[7] = cprmAnglDict[anglKey][1]
                                    else:
                                        print("Unknown %s-%s-%s, use 109 deg" %
                                              anglKey, file=sys.stderr)
                                        ic[7] = radians(109.0)
                                if ic[8] == 0.0:
                                    bondKey = tuple(map(lambda x: atomDict[x][0],
                                                        icNames[2:]))
                                    if bondKey in cprmBondDict:
                                        ic[8] = cprmBondDict[bondKey][1]
                                    else:
                                        print("Unknown %s-%s, use 1 A" %
                                              bondKey, file=sys.stderr)
                                        ic[8] = 1.0
                                val = (ic[6], ic[7], ic[8])
                                break

                        elif atomName == ic[0]:
                            isImpr = False
                            if ic[2][0] == '*':
                                isImpr = True
                                icNames = [ic[3], ic[1], ic[2][1:], ic[0]]
                            else:
                                icNames = ic[3::-1]

                            for index, icName in enumerate(icNames[:3]):
                                if icName[0] == '-':
                                    if resCnt == 0: break
                                    icName = icName[1:]
                                    icNames[index] = icName
                                    xres = pdbDataList[segCnt][1][resCnt - 1][2]
                                elif icName[0] == '+':
                                    if resCnt == last: break
                                    icName = icName[1:]
                                    icNames[index] = icName
                                    xres = pdbDataList[segCnt][1][resCnt + 1][2]
                                else:
                                    xres = pdbDataList[segCnt][1][resCnt][2]

                                for [num, name, posX, posY, posZ] in xres:
                                    if icName == name and posX is not None:
                                        pos.append([posX, posY, posZ])
                                        break
                                else:
                                    break
                            else:
                                if ic[5] == 0.0:
                                    anglKey = tuple(map(lambda x: atomDict[x][0],
                                                        icNames[1:]))
                                    if anglKey in cprmAnglDict:
                                        ic[5] = cprmAnglDict[anglKey][1]
                                    else:
                                        print("Unknown %s-%s-%s, use 109 deg" %
                                              anglKey, file=sys.stderr)
                                        ic[5] = radians(109.0)
                                if ic[4] == 0.0:
                                    bondKey = tuple(map(lambda x: atomDict[x][0],
                                                        icNames[2:]))
                                    if bondKey in cprmBondDict:
                                        ic[4] = cprmBondDict[bondKey][1]
                                    else:
                                        print("Unknown %s-%s, use 1 A" %
                                              bondKey, file=sys.stderr)
                                        ic[4] = 1.0
                                if isImpr:
                                    val = (-ic[6], ic[5], ic[4])
                                else:
                                    val = (ic[6], ic[5], ic[4])
                                break
                    else:
                        numMiss += 1
                        missAtoms.append((resiNum, resiName, atomNum, atomName))
                        continue

                    predPos = fix_miss_atom(pos, val)
                    (pdbDataList[segCnt][1][resCnt][2][atomCnt][2],
                     pdbDataList[segCnt][1][resCnt][2][atomCnt][3],
                     pdbDataList[segCnt][1][resCnt][2][atomCnt][4]) = predPos

                if numMiss == 0 or numMiss == lastNumMiss:
                    break
                else:
                    lastNumMiss = numMiss
                    numMiss = 0
                    missAtoms = []

            if missAtoms:
                for atom in missAtoms:
                    print("Cannot predict coordinates for %d %s %d %s" % atom,
                          file=sys.stderr)
                sys.exit(1)

    bondList = []
    if doDisu:
        do_disu_connect(epreInpDict['DISUBOND']["DISULIST"],
                        pdbDataList, bondList)

    for segName, segDataList in pdbDataList:
        resiNameNumDict = segNameNumDict[segName]
        nres = len(segDataList)
        last = nres - 1
        for i, resi in enumerate(segDataList):
            resiNum = resi[0]
            resiName = resi[1]

            if nres == 1:
                ctopBondList = ctopResiDataDict[segName + 'SINGLE'][1]
            elif i == 0:
                ctopBondList = ctopResiDataDict[segName + 'FIRST'][1]
            elif i == last:
                ctopBondList = ctopResiDataDict[segName + 'LAST'][1]
            else:
                ctopBondList = ctopResiDataDict[resiName][1]

            for bond in ctopBondList:
                if bond in [['H1', 'H2'], ['H2', 'H1']]:
                    print("omit the 'HT-HT' bond")
                    continue

                if i != 0:
                    prevNum = segDataList[i - 1][0]
                if i != last:
                    nextNum = segDataList[i + 1][0]

                realBond = []

                for name in bond:
                    if name[0] == '-':
                        if i == 0: break
                        realBond.append(resiNameNumDict[prevNum][name[1:]])
                    elif name[0] == '+':
                        if i == last: break
                        realBond.append(resiNameNumDict[nextNum][name[1:]])
                    else:
                        realBond.append(resiNameNumDict[resiNum][name])
                else:
                    bondList.append(realBond)

    bondList.sort(key=lambda x: x[0])

    atomAllBondList = [[] for _ in repeat(None, len(atomList) + 1)]
    for bond in bondList:
        atomAllBondList[bond[0]].append(bond[1])
        atomAllBondList[bond[1]].append(bond[0])

    anglList = []
    for atomB, atomBList in enumerate(atomAllBondList):
        if atomBList:
            for atomA, atomC in combinations(atomBList, 2):
                if atomA < atomC:
                    anglList.append([atomA, atomB, atomC])
                else:
                    anglList.append([atomC, atomB, atomA])

    anglList.sort(key=lambda x: x[0])

    diheList = []
    for atomB, atomC in bondList:
        atomBList = atomAllBondList[atomB]
        atomCList = atomAllBondList[atomC]
        if atomBList and atomCList:
            for atomA, atomD in product(atomBList, atomCList):
                if atomA != atomC and atomB != atomD and atomA != atomD:
                    if atomA < atomD:
                        diheList.append([atomA, atomB, atomC, atomD])
                    else:
                        diheList.append([atomD, atomC, atomB, atomA])

    diheList.sort(key=lambda x: x[0])

    imprList = []
    for segName, segDataList in pdbDataList:
        resiNameNumDict = segNameNumDict[segName]
        nres = len(segDataList)
        last = nres - 1
        for i, resi in enumerate(segDataList):
            resiNum = resi[0]
            resiName = resi[1]

            if nres == 1:
                ctopImprList = ctopResiDataDict[segName + 'SINGLE'][2]
            elif i == 0:
                ctopImprList = ctopResiDataDict[segName + 'FIRST'][2]
            elif i == last:
                ctopImprList = ctopResiDataDict[segName + 'LAST'][2]
            else:
                ctopImprList = ctopResiDataDict[resiName][2]

            for impr in ctopImprList:

                if i != 0:
                    prevNum = segDataList[i - 1][0]
                if i != last:
                    nextNum = segDataList[i + 1][0]

                realImpr = []

                for name in impr:
                    if name[0] == '-':
                        if i == 0: break
                        realImpr.append(resiNameNumDict[prevNum][name[1:]])
                    elif name[0] == '+':
                        if i == last: break
                        realImpr.append(resiNameNumDict[nextNum][name[1:]])
                    else:
                        realImpr.append(resiNameNumDict[resiNum][name])
                else:
                    imprList.append(realImpr)

    imprList.sort(key=lambda x: x[0])

    cmapList = []
    noCmap = epreInpDict['NOCMAP']
    for segName, segDataList in pdbDataList:
        if noCmap: break
        resiNameNumDict = segNameNumDict[segName]
        nres = len(segDataList)
        last = nres - 1
        for i, resi in enumerate(segDataList):
            resiNum = resi[0]
            resiName = resi[1]

            if nres == 1:
                ctopCmapList = ctopResiDataDict[segName + 'SINGLE'][3]
            elif i == 0:
                ctopCmapList = ctopResiDataDict[segName + 'FIRST'][3]
            elif i == last:
                ctopCmapList = ctopResiDataDict[segName + 'LAST'][3]
            else:
                ctopCmapList = ctopResiDataDict[resiName][3]

            for cmap in ctopCmapList:

                if i != 0:
                    prevNum = segDataList[i - 1][0]
                if i != last:
                    nextNum = segDataList[i + 1][0]

                realCmap = []

                for name in cmap:
                    if name[0] == '-':
                        if i == 0: break
                        realCmap.append(resiNameNumDict[prevNum][name[1:]])
                    elif name[0] == '+':
                        if i == last: break
                        realCmap.append(resiNameNumDict[nextNum][name[1:]])
                    else:
                        realCmap.append(resiNameNumDict[resiNum][name])
                else:
                    cmapList.append(realCmap)

    cmapList.sort(key=lambda x: x[0])

    if doDisu:
        do_disu_rename(pdbDataList, atomList)

    topList = [atomList, bondList, anglList, diheList, imprList, cmapList]

    return pdbDataList, topList


def fix_cryst_wat(pos, val=(radians(104.52), 0.9572)):
    ox, oy, oz = pos
    angl_hoh, bond_oh = val
    h1x = bond_oh + ox
    h1y = oy
    h1z = oz
    h2x = bond_oh * cos(angl_hoh) + ox
    h2y = bond_oh * sin(angl_hoh) + oy
    h2z = oz

    return (h1x, h1y, h1z), (h2x, h2y, h2z)


def fix_miss_atom(pos, val):
    ri, rj, rk = pos
    tors_ijkl, angl_jkl, bond_kl = val
    rjk = numpy.subtract(rk, rj)
    lenXpi = 1.0 / numpy.linalg.norm(rjk)
    xp = rjk * lenXpi
    rji = numpy.subtract(ri, rj)
    djii = 1.0 / numpy.linalg.norm(rji)
    rji *= djii
    zp = numpy.cross(xp, rji)
    lenZpi = 1.0 / numpy.linalg.norm(zp)
    zp *= lenZpi
    yp = numpy.cross(zp, xp)
    rlpx = -bond_kl * cos(angl_jkl)
    tmp = bond_kl * sin(angl_jkl)
    rlpy = tmp * cos(tors_ijkl)
    rlpz = tmp * sin(tors_ijkl)
    rl = xp * rlpx
    rl += yp * rlpy
    rl += zp * rlpz
    rl += rk

    return rl.tolist()


def add_wat(epreInpDict, sluPos):
    try:
        import cPickle as p
    except ImportError:
        import pickle as p

    watInfoDict = {
        'TIP3P': {
            'RESNAME': 'WAT',
            'ATOM': [
                ['OH2', 'OT', -0.834, 15.99940],
                ['H1', 'HT', 0.417, 1.00800],
                ['H2', 'HT', 0.417, 1.00800]
            ],
            'BOND': [[1, 2], [1, 3]],
            'ANGLE': [[2, 1, 3]],
            'COORDINATE': 'tip3p.crd',
            'SIZE': 65.4195
        }
    }

    inpInfo = epreInpDict['ADDWAT']
    watInfo = watInfoDict[inpInfo['MODEL']]
    segName = inpInfo['SEGNAME']
    cut = inpInfo['CUT']
    pad = inpInfo['PAD']
    resName = watInfo['RESNAME']
    atomParmList = watInfo['ATOM']
    bondInfo = watInfo['BOND']
    anglInfo = watInfo['ANGLE']
    #watPosFile = watInfo['COORDINATE']
    watPosFile = inpInfo['COORDINATE']
    watBoxLen = watInfo['SIZE']

    if cut <= 0.0:
        perr("Minimum distance cutoff should be > 0.0 A")
    if pad <= 0.0:
        perr("Minimum buffering distance should be > 0.0 A")

    # Add by YeJin
    #fn = os.path.join(__home_dir__, watPosFile)

    # with open(fn, 'rb') as inf:
    #     watPos = p.load(inf)

    with open(watPosFile, 'rb') as inf:
        watPos = p.load(inf)
    print("Read solvent coordinate file", watPosFile)


    print("Add explicit solvents")

    sluMax = map(max, zip(*sluPos))
    sluMin = map(min, zip(*sluPos))

    boxMax = map(lambda x: x + pad, sluMax)
    boxMin = map(lambda x: x - pad, sluMin)

    sluMax = map(lambda x: x + cut, sluMax)
    sluMin = map(lambda x: x - cut, sluMin)
    sluMaxX, sluMaxY, sluMaxZ = sluMax
    sluMinX, sluMinY, sluMinZ = sluMin
    watBoxLeni = 1.0 / watBoxLen
    boxLen = map(operator.sub, boxMax, boxMin)
    nbox = map(lambda x: int(math.ceil(x * watBoxLeni)), boxLen)

    watBox = []
    watPosArray = numpy.array(watPos)
    checkOutsideList = []
    checkOverlapList = []
    nboxX, nboxY, nboxZ = nbox
    maxX, maxY, maxZ = boxMax
    minX, minY, minZ = boxMin
    cnt = 0
    for i in range(nboxX):
        movX = minX + i * watBoxLen
        endX = movX + watBoxLen
        isInsideX = movX > minX and endX < maxX
        isNotOverlapX = movX > sluMaxX or endX < sluMinX
        for j in range(nboxY):
            movY = minY + j * watBoxLen
            endY = movY + watBoxLen
            isInsideY = movY > minY and endY < maxY
            isNotOverlapY = movY > sluMaxY or endY < sluMinY
            for k in range(nboxZ):
                movZ = minZ + k * watBoxLen
                endZ = movZ + watBoxLen
                isInsideZ = movZ > minZ and endZ < maxZ
                isNotOverlapZ = movZ > sluMaxZ or endZ < sluMinZ
                newBox = (watPosArray + [movX, movY, movZ]).tolist()
                watBox.append(newBox)
                isInside = isInsideX and isInsideY and isInsideZ
                if not isInside:
                    checkOutsideList.append(newBox)
                isNotOverlap = isNotOverlapX or isNotOverlapY or isNotOverlapZ
                if not isNotOverlap:
                    checkOverlapList.append(newBox)

    for box in checkOutsideList:
        i = 0
        while i < len(box):
            for x, y, z in box[i]:
                if x > maxX or x < minX or y > maxY or y < minY or z > maxZ or z < minZ:
                    del box[i]
                    break
            else:
                i += 1

    cuti = 1.0 / cut
    ncellX, ncellY, ncellZ = map(lambda x: int(x * cuti) + 1, boxLen)
    ncellXY = ncellX * ncellY
    ncellXYZ = ncellXY * ncellZ

    watCellList = [[] for _ in repeat(None, ncellXYZ)]
    for boxNum, box in enumerate(checkOverlapList):
        for molNum, molPos in enumerate(box):
            for atomPos in molPos:
                ix, iy, iz = map(lambda x, y: int((x - y) * cuti), atomPos, boxMin)
                cid = ix + ncellX * iy + ncellXY * iz
                watCellList[cid].append([boxNum, molNum, atomPos])

    sluCellList = [[] for _ in repeat(None, ncellXYZ)]
    for atomPos in sluPos:
        ix, iy, iz = map(lambda x, y: int((x - y) * cuti), atomPos, boxMin)
        cid = ix + ncellX * iy + ncellXY * iz
        sluCellList[cid].append(atomPos)

    cid = 0
    sluCellNb = [[] for _ in repeat(None, ncellXYZ)]
    for iz, iy, ix in product(range(ncellZ), range(ncellY), range(ncellX)):

        if not sluCellList[cid]:
            cid += 1
            continue

        newNb = watCellList[cid]
        if newNb:
            sluCellNb[cid].append(newNb)

        if ix < ncellX - 1:
            newNb = watCellList[cid + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy < ncellY - 1:
            newNb = watCellList[cid + ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy < ncellY - 1:
            newNb = watCellList[cid + ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy < ncellY - 1 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY + ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy > 0:
            newNb = watCellList[cid - ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy > 0 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY - ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy < ncellY - 1 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY + ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy < ncellY - 1 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY + ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy > 0 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY - ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy > 0 and iz < ncellZ - 1:
            newNb = watCellList[cid + ncellXY - ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0:
            newNb = watCellList[cid - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy > 0:
            newNb = watCellList[cid - ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iz > 0:
            newNb = watCellList[cid - ncellXY]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy > 0:
            newNb = watCellList[cid - ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iz > 0:
            newNb = watCellList[cid - ncellXY - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy > 0 and iz > 0:
            newNb = watCellList[cid - ncellXY - ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy < ncellY - 1:
            newNb = watCellList[cid + ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iz > 0:
            newNb = watCellList[cid - ncellXY + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if iy < ncellY - 1 and iz > 0:
            newNb = watCellList[cid - ncellXY + ncellX]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy > 0 and iz > 0:
            newNb = watCellList[cid - ncellXY - ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy > 0 and iz > 0:
            newNb = watCellList[cid - ncellXY - ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix > 0 and iy < ncellY - 1 and iz > 0:
            newNb = watCellList[cid - ncellXY + ncellX - 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        if ix < ncellX - 1 and iy < ncellY - 1 and iz > 0:
            newNb = watCellList[cid - ncellXY + ncellX + 1]
            if newNb:
                sluCellNb[cid].append(newNb)

        cid += 1

    cut2 = cut ** 2
    for cid, watMols in enumerate(sluCellNb):
        if not watMols or not sluCellList[cid]: continue
        for slu, wat in product(sluCellList[cid],
                                [atom for mol in watMols for atom in mol]):
            if (checkOverlapList[wat[0]][wat[1]] and
                                    (slu[0] - wat[2][0]) ** 2 + (slu[1] - wat[2][1]) ** 2 +
                                (slu[2] - wat[2][2]) ** 2 <= cut2):
                checkOverlapList[wat[0]][wat[1]] = []

    watBox = [filter(None, box) for box in watBox]
    watBox = filter(None, watBox)
    boxInfo = boxMax, boxMin, boxLen

    return watBox, sluPos, boxInfo


def add_ion(epreInpDict, atomList, sluPos, watPos, boxInfo):
    ionInfoDict = {
        'SOD': ['SOD', 'SOD', 1, 22.989770],
        'MG': ['MG', 'MG', 2, 24.305000],
        'POT': ['POT', 'POT', 1, 39.102000],
        'CAL': ['CAL', 'CAL', 2, 40.080000],
        'ZN2': ['ZN', 'ZN', 2, 65.370000],
        'CLA': ['CLA', 'CLA', -1, 35.450000]
    }

    atomChrgList = [atom[6] for atom in atomList]
    preTotChrg = math.fsum(atomChrgList)
    print("Total charge %f e" % preTotChrg)
    intTotChrg = int(round(preTotChrg))
    deltaChrg = preTotChrg - intTotChrg
    if math.fabs(deltaChrg) > 1.0e-3:
        perr("Non-integer charge %f e" % deltaChrg)

    if math.fabs(intTotChrg) <= 1.0e-3:
        print("The system is neutral. No extra counter-ions is needed.")

    print("Add explicit ions")

    cationList = epreInpDict['ADDION']['CATION']
    anionList = epreInpDict['ADDION']['ANION']

    nameCation = cationList[0]
    nameAnion = anionList[0]
    numCation = cationList[1]
    numAnion = anionList[1]
    chrgCation = ionInfoDict[nameCation][2]
    chrgAnion = ionInfoDict[nameAnion][2]

    if numCation == 0 and numAnion == 0:
        if intTotChrg < 0:
            if chrgCation == 1:
                numCation = -intTotChrg
            elif chrgCation == 2:
                numCation = -(intTotChrg >> 1)
                numAnion = intTotChrg & 1
        else:
            numAnion = intTotChrg

        saltCon = epreInpDict['ADDION']['SALTCON']
        if saltCon > 0.0:
            boxLen = boxInfo[2]
            # boxlen -> boxLen
            vol = reduce(operator.mul, boxLen)
            numSalt = int(round(saltCon * vol * 6.022e-4))
            numCation += numSalt
            if chrgCation == 1:
                numAnion += numSalt
            elif chrgCation == 2:
                numAnion += (numSalt << 1)

    if numCation == 0 and numAnion == 0:
        print("No ions to be added")
        return watPos, []

    ionToSolute = epreInpDict['ADDION']['ION_SOLUTE']
    ionToIon = epreInpDict['ADDION']['ION_ION']
    ionToIon2 = ionToIon ** 2
    boxMax, boxMin, boxLen = boxInfo
    ionToSolutei = 1.0 / ionToSolute
    ncellX, ncellY, ncellZ = map(lambda x: int(x * ionToSolutei) + 1, boxLen)
    ncellXY = ncellX * ncellY
    ncellXYZ = ncellXY * ncellZ

    watCellList = [[] for _ in repeat(None, ncellXYZ)]
    for boxNum, box in enumerate(watPos):
        for molNum, molPos in enumerate(box):
            for atomPos in molPos:
                ix, iy, iz = map(lambda x, y: int((x - y) * ionToSolutei),
                                 atomPos, boxMin)
                cid = ix + ncellX * iy + ncellXY * iz
                watCellList[cid].append([boxNum, molNum, atomPos])

    ionToSolute2 = ionToSolute ** 2
    sluCellNb = [None] * ncellXYZ
    delWatDict = {}
    for atomPos in sluPos:
        ix, iy, iz = map(lambda x, y: int((x - y) * ionToSolutei), atomPos, boxMin)
        cid = ix + ncellX * iy + ncellXY * iz

        if sluCellNb[cid] is None:
            sluCellNb[cid] = []

            newNb = watCellList[cid]
            if newNb:
                sluCellNb[cid].append(newNb)

            if ix < ncellX - 1:
                newNb = watCellList[cid + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy < ncellY - 1:
                newNb = watCellList[cid + ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1:
                newNb = watCellList[cid + ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY + ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0:
                newNb = watCellList[cid - ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy > 0 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY - ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY + ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY + ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY - ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz < ncellZ - 1:
                newNb = watCellList[cid + ncellXY - ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0:
                newNb = watCellList[cid - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy > 0:
                newNb = watCellList[cid - ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iz > 0:
                newNb = watCellList[cid - ncellXY]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy > 0:
                newNb = watCellList[cid - ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iz > 0:
                newNb = watCellList[cid - ncellXY - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy > 0 and iz > 0:
                newNb = watCellList[cid - ncellXY - ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1:
                newNb = watCellList[cid + ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz > 0:
                newNb = watCellList[cid - ncellXY + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz > 0:
                newNb = watCellList[cid - ncellXY + ncellX]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz > 0:
                newNb = watCellList[cid - ncellXY - ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz > 0:
                newNb = watCellList[cid - ncellXY - ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz > 0:
                newNb = watCellList[cid - ncellXY + ncellX - 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz > 0:
                newNb = watCellList[cid - ncellXY + ncellX + 1]
                if newNb:
                    sluCellNb[cid].append(newNb)

        for watMols in sluCellNb[cid]:
            if not watMols: continue
            for wat in watMols:
                boxNum = wat[0]
                molNum = wat[1]
                if ((boxNum, molNum) not in delWatDict and
                                        (atomPos[0] - wat[2][0]) ** 2 + (atomPos[1] - wat[2][1]) ** 2 +
                                    (atomPos[2] - wat[2][2]) ** 2 <= ionToSolute2):
                    delWatDict[(boxNum, molNum)] = True

    oxyList = [
        (boxNum, molNum, mol[0])
        for boxNum, box in enumerate(watPos)
        for molNum, mol in enumerate(box)
        if (boxNum, molNum) not in delWatDict
        ]

    print("Compute electrostatic potential ...")

    sluPos = numpy.array(sluPos)
    chrg = numpy.array(atomChrgList)
    noxy = len(oxyList)
    elecPotList = [None] * noxy
    for i, oxy in enumerate(oxyList):
        tmp = sluPos - oxy[2]
        tmp *= tmp
        dist2 = numpy.add.reduce(tmp, axis=1)
        tmp = chrg / dist2
        elecPotList[i] = numpy.add.reduce(tmp, axis=0).item()

        sys.stdout.flush()

    print( "Done. ")


    ionToIoni = 1.0 / ionToIon
    ncellX, ncellY, ncellZ = map(lambda x: int(x * ionToIoni) + 1, boxLen)
    ncellXY = ncellX * ncellY
    ncellXYZ = ncellXY * ncellZ

    oxyCellList = [[] for _ in repeat(None, ncellXYZ)]
    for oxy in oxyList:
        ix, iy, iz = map(lambda x, y: int((x - y) * ionToIoni), oxy[2], boxMin)
        cid = ix + ncellX * iy + ncellXY * iz
        oxyCellList[cid].append(oxy)

    oxyCellNb = [None] * ncellXYZ

    addCationList = []
    for i in range(numCation):
        minElecID = elecPotList.index(min(elecPotList))
        newIon = deepcopy(oxyList[minElecID])
        addCationList.append(newIon)

        print("Place cation %d: %s at (%.3f, %.3f, %.3f)" %
              (i + 1, nameCation, newIon[2][0], newIon[2][1], newIon[2][2]))
        sys.stdout.flush()

        boxNum = newIon[0]
        molNum = newIon[1]
        delWatDict[(boxNum, molNum)] = True
        elecPotList[minElecID] = None
        oxyList[minElecID] = None

        ionPos = newIon[2]
        ix, iy, iz = map(lambda x, y: int((x - y) * ionToIoni), ionPos, boxMin)
        cid = ix + ncellX * iy + ncellXY * iz

        if oxyCellNb[cid] is None:
            oxyCellNb[cid] = []

            newNb = oxyCellList[cid]
            if newNb:
                oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1:
                newNb = oxyCellList[cid + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0:
                newNb = oxyCellList[cid - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0:
                newNb = oxyCellList[cid - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0:
                newNb = oxyCellList[cid - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iz > 0:
                newNb = oxyCellList[cid - ncellXY]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0:
                newNb = oxyCellList[cid - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

        for oxys in oxyCellNb[cid]:
            if not oxys: continue
            for oxy in oxys:
                boxNum = oxy[0]
                molNum = oxy[1]
                oxyPos = oxy[2]
                if ((boxNum, molNum) not in delWatDict and
                                        (ionPos[0] - oxyPos[0]) ** 2 + (ionPos[1] - oxyPos[1]) ** 2 +
                                    (ionPos[2] - oxyPos[2]) ** 2 <= ionToIon2):
                    id = oxyList.index(oxy)
                    oxyList[id] = None
                    elecPotList[id] = None
                    delWatDict[(boxNum, molNum)] = True

        oxyList = [i for i in oxyList if i is not None]
        elecPotList = [i for i in elecPotList if i is not None]

        for i, oxyInfo in enumerate(oxyList):
            oxyPos = oxyInfo[2]
            r2 = ((oxyPos[0] - ionPos[0]) ** 2 + (oxyPos[1] - ionPos[1]) ** 2 +
                  (oxyPos[2] - ionPos[2]) ** 2)
            elecPotList[i] += chrgCation / r2

    addAnionList = []
    for i in range(numAnion):
        maxElecID = elecPotList.index(max(elecPotList))
        newIon = deepcopy(oxyList[maxElecID])
        addAnionList.append(newIon)

        print("Place anion %d: %s at (%.3f, %.3f, %.3f)" %
              (i + 1, nameAnion, newIon[2][0], newIon[2][1], newIon[2][2]))
        sys.stdout.flush()

        boxNum = newIon[0]
        molNum = newIon[1]
        delWatDict[(boxNum, molNum)] = True
        elecPotList[maxElecID] = None
        oxyList[maxElecID] = None

        ionPos = newIon[2]
        ix, iy, iz = map(lambda x, y: int((x - y) * ionToIoni), ionPos, boxMin)
        cid = ix + ncellX * iy + ncellXY * iz

        if oxyCellNb[cid] is None:
            oxyCellNb[cid] = []

            newNb = oxyCellList[cid]
            if newNb:
                oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1:
                newNb = oxyCellList[cid + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0:
                newNb = oxyCellList[cid - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz < ncellZ - 1:
                newNb = oxyCellList[cid + ncellXY - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0:
                newNb = oxyCellList[cid - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0:
                newNb = oxyCellList[cid - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iz > 0:
                newNb = oxyCellList[cid - ncellXY]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0:
                newNb = oxyCellList[cid - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1:
                newNb = oxyCellList[cid + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy > 0 and iz > 0:
                newNb = oxyCellList[cid - ncellXY - ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix > 0 and iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX - 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

            if ix < ncellX - 1 and iy < ncellY - 1 and iz > 0:
                newNb = oxyCellList[cid - ncellXY + ncellX + 1]
                if newNb:
                    oxyCellNb[cid].append(newNb)

        for oxys in oxyCellNb[cid]:
            if not oxys: continue
            for oxy in oxys:
                boxNum = oxy[0]
                molNum = oxy[1]
                oxyPos = oxy[2]
                if ((boxNum, molNum) not in delWatDict and
                                        (ionPos[0] - oxyPos[0]) ** 2 + (ionPos[1] - oxyPos[1]) ** 2 +
                                    (ionPos[2] - oxyPos[2]) ** 2 <= ionToIon2):
                    id = oxyList.index(oxy)
                    oxyList[id] = None
                    elecPotList[id] = None
                    delWatDict[(boxNum, molNum)] = True

        oxyList = [i for i in oxyList if i is not None]
        elecPotList = [i for i in elecPotList if i is not None]

        for i, oxyInfo in enumerate(oxyList):
            oxyPos = oxyInfo[2]
            r2 = ((oxyPos[0] - ionPos[0]) ** 2 + (oxyPos[1] - ionPos[1]) ** 2 +
                  (oxyPos[2] - ionPos[2]) ** 2)
            elecPotList[i] += chrgAnion / r2

    ionPosList = [[nameCation, ion[2]] for ion in addCationList]
    ionPosList += [[nameAnion, ion[2]] for ion in addAnionList]

    for ion in chain(addCationList, addAnionList):
        watPos[ion[0]][ion[1]] = []
    watPos = [filter(None, box) for box in watPos]
    watPos = filter(None, watPos)

    aftTotChrg = preTotChrg + numCation * chrgCation + numAnion * chrgAnion
    if numCation > 0:
        print("Added %d %s cations" % (numCation, nameCation))
    if numAnion > 0:
        print("Added %d %s anions" % (numAnion, nameAnion))
    print("Total charge %f e after ions added" % aftTotChrg)

    return watPos, ionPosList


def build_solv_top(epreInpDict, watPos, ionPosList, atomNumAdd, resNumAdd):
    watInfoDict = {
        'TIP3P': {
            'RESNAME': 'WAT',
            'ATOM': [
                ['OH2', 'OT', -0.834, 15.99940],
                ['H1', 'HT', 0.417, 1.00800],
                ['H2', 'HT', 0.417, 1.00800]
            ],
            'BOND': [[1, 2], [1, 3]],
            'ANGLE': [[2, 1, 3]],
            'COORDINATE': 'tip3p.crd',
            'SIZE': 65.4195
        }
    }

    ionInfoDict = {
        'SOD': ['SOD', 'SOD', 1, 22.989770],
        'MG': ['MG', 'MG', 2, 24.305000],
        'POT': ['POT', 'POT', 1, 39.102000],
        'CAL': ['CAL', 'CAL', 2, 40.080000],
        'ZN2': ['ZN', 'ZN', 2, 65.370000],
        'CLA': ['CLA', 'CLA', -1, 35.450000]
    }

    inpInfo = epreInpDict['ADDWAT']
    watInfo = watInfoDict[inpInfo['MODEL']]
    segNameWat = inpInfo['SEGNAME']
    resName = watInfo['RESNAME']
    atomParmList = watInfo['ATOM']
    bondInfo = watInfo['BOND']
    anglInfo = watInfo['ANGLE']

    atomNum = atomNumAdd
    pdbDataList = []
    atomList = []
    nmol = 0

    if watPos is None or watPos == []:
        bondList = []
        anglList = []
    else:
        for i, box in enumerate(watPos):
            if len(watPos) == 1:
                segName = segNameWat
            else:
                segName = segNameWat + str(i + 1)
            newSeg = [segName, []]
            res = newSeg[1]
            resNum = 0
            for mol in box:
                nmol += 1
                resNum += 1
                newRes = [resNum, resName, []]
                res.append(newRes)
                atom = newRes[2]
                for atomParm, atomPos in zip(atomParmList, mol):
                    atomNum += 1
                    atomName, atomType, atomChrg, atomMass = atomParm
                    newAtom = [atomNum, atomName, atomPos[0],
                               atomPos[1], atomPos[2]]
                    atom.append(newAtom)
                    atomList.append([atomNum, segName, resNum, resName,
                                     atomName, atomType, atomChrg, atomMass])
            pdbDataList.append(newSeg)

        molLen = len(atomParmList)
        bondList = [map(lambda x: x + atomNumAdd + i * molLen, bond)
                    for i in range(nmol) for bond in bondInfo]
        anglList = [map(lambda x: x + atomNumAdd + i * molLen, angl)
                    for i in range(nmol) for angl in anglInfo]


    segNameIon = epreInpDict['ADDION']['SEGNAME']
    resNum = 0
    segName = segNameIon
    newSeg = [segName, []]

    if ionPosList is None or ionPosList == []:
        pass
    else:
        for ion in ionPosList:
            res = newSeg[1]
            resNum += 1
            resName = ion[0]
            newRes = [resNum, resName, []]
            res.append(newRes)
            atom = newRes[2]
            atomNum += 1
            atomName, atomType, atomChrg, atomMass = ionInfoDict[ion[0]]
            atomPos = ion[1]
            newAtom = [atomNum, atomName, atomPos[0], atomPos[1], atomPos[2]]
            atom.append(newAtom)
            atomList.append([atomNum, segName, resNum, resName,
                             atomName, atomType, atomChrg, atomMass])
        pdbDataList.append(newSeg)


    return pdbDataList, atomList, bondList, anglList


def write_psf(filename, topList, epreInpDict):
    atomList, bondList, anglList, diheList, imprList, cmapList = topList
    psfFile = wopen(filename)

    print("Write X-PLOR/NAMD formatted PSF file", filename)

    noCmap = epreInpDict['NOCMAP']
    hasCmap = True if cmapList else False
    psfFile.write('PSF\n\n' if noCmap or not hasCmap else 'PSF CMAP\n\n')
    psfFile.write('%8d !NTITLE\n' % 1)
    psfFile.write(' REMARKS PSF file generated by MDBuilder v1.0 on %s\n\n' %
                  time.ctime())

    psfFile.write('%8d !NATOM\n' % len(atomList))
    fmt = '%8d %-4s %-4s %-4s %-4s %-4s %10.6f     %9.4f  %10d\n'
    for atom in atomList:
        (atomNum, segName, resNum, resName,
         atomName, atomType, atomChrg, atomMass) = atom
        psfFile.write(fmt % (atomNum, segName, resNum, resName,
                             atomName, atomType, atomChrg, atomMass, 0))
    psfFile.write('\n')

    nbond = len(bondList)
    psfFile.write('%8d !NBOND: bonds\n' % nbond)
    for i, bond in enumerate(bondList):
        bond = ' %7d %7d' % tuple(bond)
        if (i + 1) % 4 == 0 and i != nbond - 1:
            psfFile.write(bond + '\n')
        else:
            psfFile.write(bond)
    psfFile.write('\n\n')

    nangl = len(anglList)
    psfFile.write('%8d !NTHETA: angles\n' % nangl)
    for i, angl in enumerate(anglList):
        angl = ' %7d %7d %7d' % tuple(angl)
        if (i + 1) % 3 == 0 and i != nangl - 1:
            psfFile.write(angl + '\n')
        else:
            psfFile.write(angl)
    psfFile.write('\n\n')

    ndihe = len(diheList)
    psfFile.write('%8d !NPHI: dihedrals\n' % ndihe)
    for i, dihe in enumerate(diheList):
        dihe = ' %7d %7d %7d %7d' % tuple(dihe)
        if (i + 1) % 2 == 0 and i != ndihe - 1:
            psfFile.write(dihe + '\n')
        else:
            psfFile.write(dihe)
    psfFile.write('\n\n')

    nimpr = len(imprList)
    psfFile.write('%8d !NIMPHI: impropers\n' % nimpr)
    for i, impr in enumerate(imprList):
        impr = ' %7d %7d %7d %7d' % tuple(impr)
        if (i + 1) % 2 == 0 and i != nimpr - 1:
            psfFile.write(impr + '\n')
        else:
            psfFile.write(impr)
    psfFile.write('\n\n')

    psfFile.write('%8d !NDON: donors\n' % 0)
    psfFile.write('\n\n')

    psfFile.write('%8d !NACC: acceptors\n' % 0)
    psfFile.write('\n\n')

    natom = len(atomList)
    psfFile.write('%8d !NNB\n\n' % 0)
    for i in range(natom):
        nb = '%8d' % 0
        if (i + 1) % 8 == 0 and i != natom - 1:
            psfFile.write(nb + '\n')
        else:
            psfFile.write(nb)
    psfFile.write('\n\n')

    psfFile.write('%8d %7d !NGRP\n' % (1, 0))
    psfFile.write('%8d%8d%8d\n' % (0, 0, 0))
    psfFile.write('\n')

    if not noCmap and hasCmap:
        ncmap = len(cmapList)
        psfFile.write('%8d !NCRTERM: cross-terms\n' % ncmap)
        for i, cmap in enumerate(cmapList):
            psfFile.write(' %7d %7d %7d %7d %7d %7d %7d %7d\n' % tuple(cmap))
        psfFile.write('\n')

    psfFile.close()


def write_prmtop(filename, topTitles, topList, prmList, title=None,
                 isNPT=False, version='12', format='unix'):
    atomList, bondList, anglList, diheList, imprList, cmapList = topList
    bondPrm, anglPrm, ubPrm, dihePrm, imprPrm, cmapPrm, nbndPrm, nb14Prm = prmList
    prmtopFile = wopen(filename, 'wb')

    print("Write AMBER parameter/topology file", filename)

    try:
        eol = {'unix': '\n', 'dos': '\r\n', 'sys': os.linesep}[format]
    except KeyError:
        perr("Unknown file format", format)

    (natom, ntype, nbonh, mbona, nangh, manga, ndihh, mdiha, nhparm, nparm,
     nnb, nres, nbona, nanga, ndiha, nubnd, nuang, nudih, natyp, nphb,
     ifper, nbper, ngper, ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap,
     nextr) = [0] * 31

    wflg = lambda s: prmtopFile.write("%%FLAG %s%s" % (s, eol))
    wfmt = lambda s: prmtopFile.write("%%FORMAT(%s)%s" % (s, eol))
    wcmt = lambda s: prmtopFile.write("%%COMMENT  %s%s" % (s, eol))
    wdat = prmtopFile.write
    (atomNums, segNames, resNums, resNames, atomNames,
     atomTypes, atomChrgs, atomMasses) = zip(*atomList)

    isH = [0.1 < i < 2.0 for i in atomMasses]

    natom = len(atomList)
    ifbox = 1 if isNPT else 0

    resList = OrderedDict()
    lastResNum, lastResName, lastSegName = None, None, None
    for i in range(natom):
        resNum, resName, atomNum, segName = (resNums[i],
                                             resNames[i], atomNums[i], segNames[i])
        if (lastResNum is None and lastResName is None and
                    lastSegName is None) or (resNum != lastResNum
                                             or resName != lastResName or segName != lastSegName):
            key = (resNum, resName, segName)
            resList[key] = atomNum
            lastResNum, lastResName, lastSegName = resNum, resName, segName
    nres = len(resList)

    exclList = OrderedDict()
    for i in range(1, natom + 1):
        exclList[i] = [0]
    for i in bondList:
        exclList[min(i)].append(max(i))
    for i in anglList:
        exclList[min(i[::2])].append(max(i[::2]))
    for i in diheList:
        exclList[min(i[::3])].append(max(i[::3]))
    for i, j in exclList.items():
        exclList[i] = sorted(set(j))

    nnb = sum(len(i) - 1 if len(i) > 1 else 1 for _, i in exclList.items())

    startAtomNum = [i for _, i in resList.items()]
    stopAtomNum = startAtomNum[1:] + [natom]
    nmxrs = max(map(operator.sub, stopAtomNum, startAtomNum))

    utypes = OrderedDict()
    for i, j in enumerate(atomTypes):
        if j in utypes:
            utypes[j].append(i)
            continue
        utypes[j] = [i]
    ntype = len(utypes)
    natyp = ntype

    ubonds = OrderedDict()
    for i, j in enumerate(bondList):
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in ubonds:
            ubonds[key].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in ubonds:
            ubonds[key].append(i)
            continue
        ubonds[key] = [i]
    nubnd = len(ubonds)

    nbond = len(bondList)
    ubondList = [None] * nbond
    for i, j in enumerate(ubonds.items()):
        for k in j[1]:
            ubondList[k] = i + 1
    bonhaList = [list(map(lambda x: (x - 1) * 3, i)) for i in bondList]
    bonhaList = [j + [ubondList[i]] for i, j in enumerate(bonhaList)]
    bonhList = [i for i in bonhaList if any(isH[j // 3] for j in i[:2])]
    nbonh = len(bonhList)
    bonaList = [i for i in bonhaList if not any(isH[j // 3] for j in i[:2])]
    mbona = len(bonaList)
    nbona = mbona

    uangls = OrderedDict()
    for i, j in enumerate(anglList):
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in uangls:
            uangls[key].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in uangls:
            uangls[key].append(i)
            continue
        uangls[key] = [i]
    nuang = len(uangls)

    nangl = len(anglList)
    uanglList = [None] * nangl
    for i, j in enumerate(uangls.items()):
        for k in j[1]:
            uanglList[k] = i + 1
    anghaList = [list(map(lambda x: (x - 1) * 3, i)) for i in anglList]
    anghaList = [j + [uanglList[i]] for i, j in enumerate(anghaList)]
    anghList = [i for i in anghaList if any(isH[j // 3] for j in i[:3])]
    nangh = len(anghList)
    angaList = [i for i in anghaList if not any(isH[j // 3] for j in i[:3])]
    manga = len(angaList)
    nanga = manga

    ubList = []
    for i in anglList:
        atoms = list(map(lambda x: x - 1, i))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in ubPrm:
            ubList.append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in ubPrm:
            ubList.append(i)
            continue

    uubs = OrderedDict()
    for i, j in enumerate(ubList):
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in uubs:
            uubs[key].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in uubs:
            uubs[key].append(i)
            continue
        uubs[key] = [i]

    ndihe = len(diheList)
    i = 0
    while i < ndihe:
        j = diheList[i]
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in dihePrm:
            key2 = key
        elif key[::-1] in dihePrm:
            key2 = key[::-1]
        elif ('X', key[1], key[2], 'X') in dihePrm:
            key2 = ('X', key[1], key[2], 'X')
        elif ('X', key[2], key[1], 'X') in dihePrm:
            key2 = ('X', key[2], key[1], 'X')
        else:
            perr("Unknown dihedral %s" % key)
        inc = len(dihePrm[key2]) - 1
        if inc > 0:
            ndihe += inc
            for k in range(inc):
                i += 1
                if j[2] != 0:
                    diheList.insert(i, [j[0], j[1], -j[2], j[3]])
                else:
                    diheList.insert(i, [j[0], -j[1], j[2], j[3]])
        i += 1

    udihes = OrderedDict()
    for i, j in enumerate(diheList):
        atoms = list(map(lambda x: int(math.copysign(abs(x) - 1, x)), j))
        key = tuple(map(lambda x: atomTypes[abs(x)], atoms))
        if key in udihes:
            udihes[key][0].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[abs(x)], atoms[::-1]))
        if key in udihes:
            udihes[key][0].append(i)
            continue
        udihes[key] = [[i], None]

    for key, val in udihes.items():
        if key in dihePrm:
            key2 = key
        elif key[::-1] in dihePrm:
            key2 = key[::-1]
        elif ('X', key[1], key[2], 'X') in dihePrm:
            key2 = ('X', key[1], key[2], 'X')
        elif ('X', key[2], key[1], 'X') in dihePrm:
            key2 = ('X', key[2], key[1], 'X')
        else:
            perr("Unknown dihedral %s" % key)
        val[1] = len(dihePrm[key2])

    nudih = sum(val[1] for key, val in udihes.items())

    udiheList = [None] * ndihe
    cnt = 1
    for i, j in udihes.items():
        if j[1] == 1:
            for k in j[0]:
                udiheList[k] = cnt
        else:
            for k in j[0][::j[1]]:
                for l in range(j[1]):
                    udiheList[k + l] = cnt + l
        cnt += j[1]

    dihhaList = [list(map(lambda x: int(math.copysign((abs(x) - 1) * 3, x)), i)) for i in diheList]
    dihhaList = [j + [udiheList[i]] for i, j in enumerate(dihhaList)]
    dihhList = [i for i in dihhaList if any(isH[abs(j) // 3] for j in i[:4])]
    ndihh = len(dihhList)
    dihaList = [i for i in dihhaList if not any(isH[abs(j) // 3] for j in i[:4])]
    mdiha = len(dihaList)
    ndiha = mdiha

    uimprs = OrderedDict()
    for i, j in enumerate(imprList):
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in uimprs:
            uimprs[key].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in uimprs:
            uimprs[key].append(i)
            continue
        uimprs[key] = [i]

    ucmaps = OrderedDict()
    for i, j in enumerate(cmapList):
        atoms = list(map(lambda x: x - 1, j))
        key = tuple(map(lambda x: atomTypes[x], atoms))
        if key in ucmaps:
            ucmaps[key].append(i)
            continue
        key = tuple(map(lambda x: atomTypes[x], atoms[::-1]))
        if key in ucmaps:
            ucmaps[key].append(i)
            continue
        ucmaps[key] = [i]

    date = time.strftime("%m/%d/%y  %H:%M:%S", time.localtime())
    wdat("%%VERSION  VERSION_STAMP = V0001.000  DATE = %s%s" % (date, eol))

    wflg("CTITLE")
    wfmt("A80")
    if title is None:
        title = "GENERATED BY MDBuilder 1.0"
    wdat("%-80s%s" % (title, eol))

    wflg("POINTERS")
    wfmt("10I8")

    fmt = "%8d" * 10 + eol
    wdat(fmt % (natom, ntype, nbonh, mbona, nangh, manga, ndihh, mdiha, nhparm, nparm))
    wdat(fmt % (nnb, nres, nbona, nanga, ndiha, nubnd, nuang, nudih, natyp, nphb))
    wdat(fmt % (ifper, nbper, ngper, ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap))
    fmt = "%8d" + eol
    wdat(fmt % nextr)

    wflg("FORCE_FIELD_TYPE")
    wfmt("I2,A78")
    topVer = int(topTitles[-1].split()[0])
    topTitle = topTitles[0]
    wdat(" 1 CHARMM%4d  %65s%s" % (topVer, topTitle, eol))

    wflg("ATOM_NAME")
    wfmt("20A4")
    for i, atom in enumerate(atomNames):
        wdat("%-4s" % atom)
        if not (i + 1) % 20 or i == natom - 1: wdat(eol)

    wflg("CHARGE")
    wcmt("Atomic charge multiplied by sqrt(332.0716D0) (CCELEC) ")
    wfmt("3E24.16")
    ccelec = sqrt(332.0716)
    for i, atom in enumerate(atomChrgs):
        wdat("%24.16E" % (atom * ccelec))
        if not (i + 1) % 3 or i == natom - 1: wdat(eol)

    wflg("MASS")
    wfmt("5E16.8")
    for i, atom in enumerate(atomMasses):
        wdat("%16.8E" % atom)
        if not (i + 1) % 5 or i == natom - 1: wdat(eol)

    wflg("ATOM_TYPE_INDEX")
    wfmt("10I8")
    utypeList = [None] * natom
    for i, j in enumerate(utypes.items()):
        for k in j[1]:
            utypeList[k] = i + 1
    for i, atom in enumerate(utypeList):
        wdat("%8d" % atom)
        if not (i + 1) % 10 or i == natom - 1: wdat(eol)

    wflg("NUMBER_EXCLUDED_ATOMS")
    wfmt("10I8")
    for i, j in enumerate(exclList):
        num = len(exclList[j])
        num = num - 1 if num > 1 else 1
        wdat("%8d" % num)
        if not (i + 1) % 10 or i == natom - 1: wdat(eol)

    wflg("EXCLUDED_ATOMS_LIST")
    wfmt("10I8")
    cnt = 0
    for _, i in exclList.items():
        for j in (i[1:] if len(i) > 1 else i):
            wdat("%8d" % j)
            cnt += 1
            if not cnt % 10 or cnt == nnb: wdat(eol)

    wflg("NONBONDED_PARM_INDEX")
    wfmt("10I8")
    ijid = [[None] * ntype for _ in range(ntype)]
    cnt = 0
    for i in range(ntype):
        for j in range(i + 1):
            cnt += 1
            ijid[i][j] = cnt
            ijid[j][i] = cnt
    ijid = [j for i in ijid for j in i]
    for id, i in enumerate(ijid):
        wdat("%8d" % i)
        if not (id + 1) % 10 or id == len(ijid) - 1: wdat(eol)

    wflg("RESIDUE_LABEL")
    wfmt("20A4")
    for i, j in enumerate(resList):
        wdat("%-4s" % j[1])
        if not (i + 1) % 20 or i == nres - 1: wdat(eol)

    wflg("RESIDUE_POINTER")
    wfmt("10I8")
    for i, j in enumerate(resList.items()):
        wdat("%8d" % j[1])
        if not (i + 1) % 10 or i == nres - 1: wdat(eol)

    wflg("BOND_FORCE_CONSTANT")
    wfmt("5E16.8")
    for i, key in enumerate(ubonds):
        if key in bondPrm:
            wdat("%16.8E" % bondPrm[key][0])
        elif key[::-1] in bondPrm:
            wdat("%16.8E" % bondPrm[key[::-1]][0])
        else:
            perr("Unknown bond", key)
        if not (i + 1) % 5 or i == nubnd - 1: wdat(eol)

    wflg("BOND_EQUIL_VALUE")
    wfmt("5E16.8")
    for i, key in enumerate(ubonds):
        if key in bondPrm:
            wdat("%16.8E" % bondPrm[key][1])
        elif key[::-1] in bondPrm:
            wdat("%16.8E" % bondPrm[key[::-1]][1])
        else:
            perr("Unknown bond", key)
        if not (i + 1) % 5 or i == nubnd - 1: wdat(eol)

    wflg("ANGLE_FORCE_CONSTANT")
    wfmt("5E16.8")
    for i, key in enumerate(uangls):
        if key in anglPrm:
            wdat("%16.8E" % anglPrm[key][0])
        elif key[::-1] in anglPrm:
            wdat("%16.8E" % anglPrm[key[::-1]][0])
        else:
            perr("Unknown angle", key)
        if not (i + 1) % 5 or i == nuang - 1: wdat(eol)

    wflg("ANGLE_EQUIL_VALUE")
    wfmt("3E25.17")
    for i, key in enumerate(uangls):
        if key in anglPrm:
            wdat("%25.17E" % anglPrm[key][1])
        elif key[::-1] in anglPrm:
            wdat("%25.17E" % anglPrm[key[::-1]][1])
        else:
            perr("Unknown angle", key)
        if not (i + 1) % 3 or i == nuang - 1: wdat(eol)

    wflg("CHARMM_UREY_BRADLEY_COUNT")
    wcmt("V(ub) = K_ub(r_ik - R_ub)**2")
    wcmt("Number of Urey Bradley terms and types")
    wfmt("2I8")
    nub = len(ubList)
    nuub = len(uubs)
    wdat("%8d%8d%s" % (nub, nuub, eol))

    wflg("CHARMM_UREY_BRADLEY")
    wcmt("List of the two atoms and its parameter index")
    wcmt("in each UB term: i,k,index")
    wfmt("10I8")

    uubList = [None] * nub
    for i, j in enumerate(uubs.items()):
        for k in j[1]:
            uubList[k] = i + 1
    ubhaList = [j + [uubList[i]] for i, j in enumerate(ubList)]

    cnt = 0
    for atom in ubhaList:
        for i in (atom[0], atom[2], atom[3]):
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 3 * nub: wdat(eol)

    wflg("CHARMM_UREY_BRADLEY_FORCE_CONSTANT")
    wcmt("K_ub: kcal/mole/A**2")
    wfmt("5E16.8")
    for i, key in enumerate(uubs):
        if key in ubPrm:
            wdat("%16.8E" % ubPrm[key][0])
        elif key[::-1] in ubPrm:
            wdat("%16.8E" % ubPrm[key[::-1]][0])
        else:
            perr("Unknown UB angle", key)
        if not (i + 1) % 5 or i == nuub - 1: wdat(eol)

    wflg("CHARMM_UREY_BRADLEY_EQUIL_VALUE")
    wcmt("r_ub: A")
    wfmt("5E16.8")
    for i, key in enumerate(uubs):
        if key in ubPrm:
            wdat("%16.8E" % ubPrm[key][1])
        elif key[::-1] in ubPrm:
            wdat("%16.8E" % ubPrm[key[::-1]][1])
        else:
            perr("Unknown UB angle", key)
        if not (i + 1) % 5 or i == nuub - 1: wdat(eol)

    wflg("DIHEDRAL_FORCE_CONSTANT")
    wfmt("5E16.8")
    cnt = 0
    for key in udihes:
        if key in dihePrm:
            key2 = key
        elif key[::-1] in dihePrm:
            key2 = key[::-1]
        elif ('X', key[1], key[2], 'X') in dihePrm:
            key2 = ('X', key[1], key[2], 'X')
        elif ('X', key[2], key[1], 'X') in dihePrm:
            key2 = ('X', key[2], key[1], 'X')
        else:
            perr("Unknown dihedral %s" % key)
        for i in dihePrm[key2]:
            cnt += 1
            wdat("%16.8E" % i[0])
        if not cnt % 5 or cnt == nudih: wdat(eol)

    wflg("DIHEDRAL_PERIODICITY")
    wfmt("5E16.8")
    cnt = 0
    for key in udihes:
        if key in dihePrm:
            key2 = key
        elif key[::-1] in dihePrm:
            key2 = key[::-1]
        elif ('X', key[1], key[2], 'X') in dihePrm:
            key2 = ('X', key[1], key[2], 'X')
        elif ('X', key[2], key[1], 'X') in dihePrm:
            key2 = ('X', key[2], key[1], 'X')
        else:
            perr("Unknown dihedral %s" % key)
        for i in dihePrm[key2]:
            cnt += 1
            wdat("%16.8E" % i[1])
        if not cnt % 5 or cnt == nudih: wdat(eol)

    wflg("DIHEDRAL_PHASE")
    wfmt("5E16.8")
    cnt = 0
    for key in udihes:
        if key in dihePrm:
            key2 = key
        elif key[::-1] in dihePrm:
            key2 = key[::-1]
        elif ('X', key[1], key[2], 'X') in dihePrm:
            key2 = ('X', key[1], key[2], 'X')
        elif ('X', key[2], key[1], 'X') in dihePrm:
            key2 = ('X', key[2], key[1], 'X')
        else:
            perr("Unknown dihedral %s" % key)
        for i in dihePrm[key2]:
            cnt += 1
            wdat("%16.8E" % i[2])
        if not cnt % 5 or cnt == nudih: wdat(eol)

    wflg("SCEE_SCALE_FACTOR")
    wfmt("5E16.8")
    for i in range(natom):
        wdat("%16.8E" % 1.0)
        if not (i + 1) % 5 or i == natom - 1: wdat(eol)

    wflg("SCNB_SCALE_FACTOR")
    wfmt("5E16.8")
    for i in range(natom):
        wdat("%16.8E" % 1.0)
        if not (i + 1) % 5 or i == natom - 1: wdat(eol)

    wflg("CHARMM_NUM_IMPROPERS")
    wcmt("Number of terms contributing to the")
    wcmt("quadratic four atom improper energy term:")
    wcmt("V(improper) = K_psi(psi - psi_0)**2")
    wfmt("10I8")
    nimpr = len(imprList)
    wdat("%8d%s" % (nimpr, eol))

    wflg("CHARMM_IMPROPERS")
    wcmt("List of the four atoms in each improper term")
    wcmt("i,j,k,l,index  i,j,k,l,index")
    wcmt("where index is into the following two lists:")
    wcmt("CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}")
    wfmt("10I8")

    uimprList = [None] * nimpr
    for i, j in enumerate(uimprs.items()):
        for k in j[1]:
            uimprList[k] = i + 1
    imphaList = [j + [uimprList[i]] for i, j in enumerate(imprList)]

    cnt = 0
    for atom in imphaList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 5 * nimpr: wdat(eol)

    wflg("CHARMM_NUM_IMPR_TYPES")
    wcmt("Number of unique parameters contributing to the")
    wcmt("quadratic four atom improper energy term")
    wfmt("I8")
    nuimpr = len(uimprs)
    wdat("%8d%s" % (nuimpr, eol))

    wflg("CHARMM_IMPROPER_FORCE_CONSTANT")
    wcmt("K_psi: kcal/mole/rad**2")
    wfmt("5E16.8")
    for i, key in enumerate(uimprs):
        if key in imprPrm:
            wdat("%16.8E" % imprPrm[key][0])
        elif key[::-1] in imprPrm:
            wdat("%16.8E" % imprPrm[key[::-1]][0])
        elif (key[0], 'X', 'X', key[3]) in imprPrm:
            wdat("%16.8E" % imprPrm[(key[0], 'X', 'X', key[3])][0])
        elif (key[3], 'X', 'X', key[0]) in imprPrm:
            wdat("%16.8E" % imprPrm[(key[3], 'X', 'X', key[0])][0])
        else:
            perr("Unknown improper dihedral", key)
        if not (i + 1) % 5 or i == nuimpr - 1: wdat(eol)

    wflg("CHARMM_IMPROPER_PHASE")
    wcmt("psi: degrees")
    wfmt("5E16.8")
    for i, key in enumerate(uimprs):
        if key in imprPrm:
            wdat("%16.8E" % imprPrm[key][1])
        elif key[::-1] in imprPrm:
            wdat("%16.8E" % imprPrm[key[::-1]][1])
        elif (key[0], 'X', 'X', key[3]) in imprPrm:
            wdat("%16.8E" % imprPrm[(key[0], 'X', 'X', key[3])][1])
        elif (key[3], 'X', 'X', key[0]) in imprPrm:
            wdat("%16.8E" % imprPrm[(key[3], 'X', 'X', key[0])][1])
        else:
            perr("Unknown improper dihedral", key)
        if not (i + 1) % 5 or i == nuimpr - 1: wdat(eol)

    wflg("SOLTY")
    wfmt("5E16.8")
    for i in range(natyp):
        wdat("%16.8E" % 0.0)
        if not (i + 1) % 5 or i == natyp - 1: wdat(eol)

    wflg("LENNARD_JONES_ACOEF")
    wfmt("3E24.16")
    lja = []
    ljb = []
    lj14a = []
    lj14b = []
    for id, i in enumerate(utypes):
        for j in utypes.keys()[:id + 1]:
            r = nbndPrm[i][1] + nbndPrm[j][1]
            r6 = r ** 6
            e = sqrt(nbndPrm[i][0] * nbndPrm[j][0])
            b = e * r6
            a = b * r6
            b += b
            lja.append(a)
            ljb.append(b)
            r = nb14Prm[i][1] + nb14Prm[j][1]
            r6 = r ** 6
            e = sqrt(nb14Prm[i][0] * nb14Prm[j][0])
            b = e * r6
            a = b * r6
            b += b
            lj14a.append(a)
            lj14b.append(b)
            cnt += 1
    for id, i in enumerate(lja):
        wdat("%24.16E" % i)
        if not (id + 1) % 3 or id == len(lja) - 1: wdat(eol)

    wflg("LENNARD_JONES_BCOEF")
    wfmt("3E24.16")
    for id, i in enumerate(ljb):
        wdat("%24.16E" % i)
        if not (id + 1) % 3 or id == len(ljb) - 1: wdat(eol)

    wflg("LENNARD_JONES_14_ACOEF")
    wfmt("3E24.16")
    for id, i in enumerate(lj14a):
        wdat("%24.16E" % i)
        if not (id + 1) % 3 or id == len(lj14a) - 1: wdat(eol)

    wflg("LENNARD_JONES_14_BCOEF")
    wfmt("3E24.16")
    for id, i in enumerate(lj14b):
        wdat("%24.16E" % i)
        if not (id + 1) % 3 or id == len(lj14b) - 1: wdat(eol)

    wflg("BONDS_INC_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in bonhList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 3 * nbonh: wdat(eol)

    wflg("BONDS_WITHOUT_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in bonaList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 3 * mbona: wdat(eol)

    wflg("ANGLES_INC_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in anghList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 4 * nangh: wdat(eol)

    wflg("ANGLES_WITHOUT_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in angaList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 4 * manga: wdat(eol)

    wflg("DIHEDRALS_INC_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in dihhList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 5 * ndihh: wdat(eol)

    wflg("DIHEDRALS_WITHOUT_HYDROGEN")
    wfmt("10I8")
    cnt = 0
    for atom in dihaList:
        for i in atom:
            wdat("%8d" % i)
            cnt += 1
            if not cnt % 10 or cnt == 5 * mdiha: wdat(eol)

    wflg("HBOND_ACOEF")
    wfmt("5E16.8)")
    wdat(eol)

    wflg("HBOND_BCOEF")
    wfmt("5E16.8)")
    wdat(eol)

    wflg("HBCUT")
    wfmt("5E16.8)")
    wdat(eol)

    wflg("AMBER_ATOM_TYPE")
    wfmt("20A4")
    for i, atom in enumerate(atomTypes):
        wdat("%-4s" % atom)
        if not (i + 1) % 20 or i == natom - 1: wdat(eol)

    wflg("TREE_CHAIN_CLASSIFICATION")
    wfmt("20A4")
    for i in range(natom):
        wdat("%-4s" % "BLA")
        if not (i + 1) % 20 or i == natom - 1: wdat(eol)

    wflg("JOIN_ARRAY")
    wfmt("10I8")
    for i in range(natom):
        wdat("%8d" % 0)
        if not (i + 1) % 10 or i == natom - 1: wdat(eol)

    wflg("IROTAT")
    wfmt("10I8")
    for i in range(natom):
        wdat("%8d" % 0)
        if not (i + 1) % 10 or i == natom - 1: wdat(eol)

    wflg("RADIUS_SET")
    wfmt("1A80")
    wdat("%-80s%s" % ("modified Bondi radii (mbondi)", eol))

    wflg("RADII")
    wfmt("5E16.8")
    bonhDict = {}
    for i in bonhList:
        if isH[i[0] // 3]:
            bonhDict[i[0] // 3] = i[1] // 3
        elif isH[i[1] // 3]:
            bonhDict[i[1] // 3] = i[0] // 3
        else:
            perr()
    radiiList = []
    screenList = []
    for i in range(natom):
        m = int(round(atomMasses[i]))
        if m == 1:
            a = bonhDict[i]
            ma = int(round(atomMasses[a]))
            if ma == 12:
                radiiList.append(1.3)
            elif ma == 14:
                radiiList.append(1.3)
            elif ma == 16:
                radiiList.append(0.8)
            elif ma == 32:
                radiiList.append(0.8)
            else:
                radiiList.append(1.2)
            screenList.append(0.85)
        elif m == 12:
            radiiList.append(1.7)
            screenList.append(0.72)
        elif m == 14:
            radiiList.append(1.55)
            screenList.append(0.79)
        elif m == 16:
            radiiList.append(1.5)
            screenList.append(0.85)
        elif m == 19:
            radiiList.append(1.5)
            screenList.append(0.88)
        elif m == 31:
            radiiList.append(1.85)
            screenList.append(0.86)
        elif m == 32:
            radiiList.append(1.8)
            screenList.append(0.96)
        elif m == 35:
            radiiList.append(1.7)
            screenList.append(0.8)
        else:
            radiiList.append(1.5)
            screenList.append(0.8)

    for i, j in enumerate(radiiList):
        wdat("%16.8E" % j)
        if not (i + 1) % 5 or i == natom - 1: wdat(eol)

    wflg("SCREEN")
    wfmt("5E16.8")
    for i, j in enumerate(screenList):
        wdat("%16.8E" % j)
        if not (i + 1) % 5 or i == natom - 1: wdat(eol)

    wflg("CHARMM_CMAP_COUNT")
    wcmt("Number of CMAP terms, number of unique CMAP parameters")
    wfmt("2I8")

    ncmap = len(cmapList)
    nucmap = len(ucmaps)
    ucmapList = [None] * ncmap
    for i, j in enumerate(ucmaps.items()):
        for k in j[1]:
            ucmapList[k] = i + 1
    cmahaList = [j[:4] + [j[-1], ucmapList[i]] for i, j in enumerate(cmapList)]

    wdat("%8d%8d%s" % (ncmap, nucmap, eol))

    wflg("CHARMM_CMAP_RESOLUTION")
    wcmt("Number of steps along each phi/psi CMAP axis")
    wcmt("for each CMAP_PARAMETER grid")
    wfmt("20I4")
    for i, key in enumerate(ucmaps):
        if key in cmapPrm:
            wdat("%4d" % len(cmapPrm[key]))
        elif key[::-1] in cmapPrm:
            wdat("%4d" % len(cmapPrm[key[::-1]]))
        else:
            perr("Uknown CMAP", key)
        if not (i + 1) % 20 or i == nucmap - 1: wdat(eol)

    for i, key in enumerate(ucmaps):
        wflg("CHARMM_CMAP_PARAMETER_%02d" % (i + 1))
        if key in cmapPrm:
            key2 = key
        else:
            key2 = key[::-1]
        title = " ".join(map(lambda x: x.ljust(4), key2))
        wcmt("     " + title)
        wfmt("8(F9.5)")
        nval = len(cmapPrm[key2]) ** 2
        cnt = 0
        for j in cmapPrm[key2]:
            for k in j:
                wdat("%9.5f" % k)
                cnt += 1
                if not cnt % 8 or cnt == nval: wdat(eol)

    wflg("CHARMM_CMAP_INDEX")
    wcmt("Atom index i,j,k,l,m of the cross term")
    wcmt("and then pointer to CHARMM_CMAP_PARAMETER_n")
    wfmt("6I8")
    fmt = "%8d" * 6 + eol
    for atom in cmahaList:
        wdat(fmt % tuple(atom))

    if ifbox > 0:
        bondingList = [[] for _ in range(natom + 1)]
        for i in bondList:
            bondingList[min(i)].append(max(i))
        molList = []
        leftAtoms = atomNums
        while leftAtoms:
            rootAtom = min(leftAtoms)
            newMol = [rootAtom]
            leafAtom = bondingList[rootAtom]
            while leafAtom:
                newMol.extend(leafAtom)
                leafAtom = set(nextAtom
                               for atom in leafAtom
                               if bondingList[atom]
                               for nextAtom in bondingList[atom])
            newMol = set(newMol)
            leftAtoms = set(leftAtoms) - newMol
            molList.append(sorted(newMol))
        molList.sort(key=lambda x: x[0])

        wflg("SOLVENT_POINTERS")
        wcmt("TODO")
        wfmt("3I8")
        nmol = len(molList)
        wdat("%8d%8d%8d%s" % (0, nmol, 0, eol))

        wflg("ATOMS_PER_MOLECULE")
        wfmt("12I6")
        for i, mol in enumerate(molList):
            wdat("%6d" % len(mol))
            if not (i + 1) % 6 or i == nmol - 1: wdat(eol)

    prmtopFile.close()


def write_inpcrd(filename, crd, boxInfo=None, title=None):
    inpcrdFile = wopen(filename, 'wb')

    print("Write AMBER input coordinate file", filename)

    natom = len(crd)

    if title is None:
        title = "GENERATED BY MDBuilder 1.0 ON %s" % time.ctime().upper()
    inpcrdFile.write("%-80s\n%5d\n" % (title, natom))

    for i, xyz in enumerate(crd):
        inpcrdFile.write("%12.7f%12.7f%12.7f" % xyz)
        if i & 1:
            inpcrdFile.write("\n")

    if natom & 1:
        inpcrdFile.write("\n")

    if boxInfo is not None:
        fmt = "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n"
        inpcrdFile.write(fmt % tuple(boxInfo))

    inpcrdFile.close()


def print_logo():
    print("\n----------------------------------------------------------------------------\n"
          "----------------------------- %s %s -----------------------------\n"
          "----------------------------------------------------------------------------\n"
          % (__program__, __version__)
          )
    sys.stdout.flush()


class PDBDownloadError(Exception):
    pass


class CleanableEntryField(Pmw.EntryField):
    def __init__(self, *args, **kwargs):
        Pmw.EntryField.__init__(self, *args, **kwargs)
        self.component('entry').bind('<Escape>',
                                     func=lambda x: self.setvalue(''))


class StdoutRedirector:
    def __init__(self, widget):
        self.widget = widget

    def write(self, content):
        self.widget.configure(state='normal')
        self.widget.insert('end', content)
        self.widget.configure(state='disabled')
        self.widget.see('end')

    def flush(self):
        pass


class MDBuilderGui:
    """MDBuilder GUI plugin."""

    def __init__(self, app):
        self.parent = app.root
        self.mod = None
        self.top = None
        self.topList = None
        self.prm = None
        self.sluPos = None
        self.watPos = None
        self.sluPos = None
        self.boxInfo = None

        self.epreInpDict = {
            'TOPOLOGY': None,
            'FORCEFIELD': None,
            'NOCMAP': False,
            'ALIASRES': [],
            'ALIASATOM': [],
            'PATCH': [],
            'SEGMENT': [],
            'COORDPDB': None,
            'DOWNLOAD': False,
            'NOGUESSCOORD': False,
            'DISUBOND': {
                'DODISU': False,
                'AUTO': False,
                'CUT': 2.1,
                'DISULIST': []
            },
            'ADDWAT': {
                'DOADDWAT': False,
                'MODEL': 'TIP3P',
                'SEGNAME': 'WT',
                'CUT': 2.4,
                'PAD': 9.0,
                'COORDINATE': 'tip3p.crd'
            },
            'ADDION': {
                'DOADDION': False,
                'METHOD': 'RANDOM',
                'SEGNAME': 'ION',
                'CATION': ['SOD', 0],
                'ANION': ['CLA', 0],
                'ION_SOLUTE': 5.0,
                'ION_ION': 5.0,
                'SALTCON': 0.0
            },
            'WRITEPSF': None,
            'WRITEPDB': None,
            'WRITEPRMTOP': None,
            'WRITEINPCRD': None
        }

        self.pmobj = []
        self.original_stdout = sys.stdout
        self.create_widgets()
        cmd.hide('everything', 'all')

    def __del__(self):
        sys.stdout = self.original_stdout

    def create_widgets(self):
        # dialog window
        # ----------
        self.create_dialog()

        # title label
        # ----------
        self.create_title()

        # paned window
        # ----------
        self.create_panedwindow()

        # multiple pages
        # ----------
        self.create_notebook()

        # console
        # ----------
        self.create_console()

        self.show_dialog()

    def create_dialog(self):
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Execute', 'Output', 'Console',
                                          'Quit', 'About'),
                                 title='%s %s'%(__program__, __version__),
                                 command=self.on_dialog_button_clicked)
        w = self.dialog.component('buttonbox')
        for i in range(w.numbuttons()):
            w.button(i).configure(width=10, foreground='#808080', activeforeground='#336699',
                                  relief=FLAT, font = ('Helvetica', 10, 'bold'))
            w.button(i).bind("<Enter>", self.rolloverEnter)
            w.button(i).bind("<Leave>", self.rolloverLeave)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

    def rolloverEnter(self, event):
        event.widget.config(foreground='#336699',relief = GROOVE)

    def rolloverLeave(self, event):
        event.widget.config(foreground='#808080',relief = FLAT)

    def create_title(self):
        w = Frame(self.dialog.interior(), relief='raised', borderwidth=0)
        w.pack(fill='both', expand=0, padx=10, pady=4)
        Label(w,
             text='%s %s\n%s'%(__program__, __version__, __desc__),
             font=Pmw.logicalfont('Helvetica', 6, weight='bold'),
             bg='#6699CC',
             fg='white'
             ).pack(fill='both', expand=0, ipadx=10, ipady=10)

    def create_panedwindow(self):
        self.panedwin = PanedWindow(self.dialog.interior(),
                                    orient='vertical', sashrelief='sunken')
        self.panedwin.pack(fill='both', expand=1)

    def create_console(self):
        self.console_frame = Frame(self.panedwin)
        self.console = Text(self.console_frame,
              height=10,
              width=80,
              bg='#DCDCDC',
              fg='#696969',
              relief='sunken')
        sbar = Scrollbar(self.console_frame, orient='vertical',
                         command=self.console.yview)
        self.console.config(yscrollcommand=sbar.set)
        self.console.pack(side='left', fill='both', expand=1)
        sbar.pack(side='right', fill='y')
        sys.stdout = StdoutRedirector(self.console)
        self.console_shown = 0

    def create_notebook(self):
        self.notebook = Pmw.NoteBook(self.panedwin)
        self.panedwin.add(self.notebook)
        self.panedwin.paneconfigure(self.notebook, padx=10, pady=10)

        # "I/O" page
        # ==========
        self.create_io_page()

        # "Preparation" page
        # ==========
        self.create_prep_page()

        # "Solvation" page
        # ==========
        self.create_sol_page()

        # "Ionization" page
        # ==========
        self.create_ion_page()

        # make tab sizes the same
        for w in self.notebook._pageAttrs.values():
            # modify by Yejin
            #w['tabreqwidth'] = 80
            w['tabreqwidth'] = 90

        self.notebook.setnaturalsize()

    def create_io_page(self):
        page = self.notebook.add('I/O')
        self.notebook.tab(0).focus_set()
        self.notebook.component('I/O-tab').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        self.nocmap = IntVar()
        self.nocmap.set(0)

        grp_opt = {'fill': 'both', 'expand': 1, 'padx': 10, 'pady': 5}
        frm_opt = {'fill': 'both', 'expand': 1}
        ent_opt = {'side': 'left', 'fill': 'both',
                   'expand': 1, 'padx': 10, 'pady': 5}
        btn_opt = {'side': 'right', 'fill': 'x',
                   'expand': 0, 'padx': 10, 'pady': 5}

        # "Input" group
        # **********
        group = Pmw.Group(page, tag_text='Input Files')
        group.pack(**grp_opt)
        group.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.pdbloc = CleanableEntryField(
                frame,
                labelpos='w',
                label_text='PDB File:',
                command=self.on_pdbentry_pressed)

        self.openpdbbtn = Button(
                frame,
                command=self.on_openpdb_clicked,
                text='Browse',
                width=10)

        self.downloadbtn = Button(
                frame,
                command=self.on_download_clicked,
                text='Download',
                width=10)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.ffloc = CleanableEntryField(
                frame,
                labelpos='w',
                label_text='Topology File:')
                #label_text='Forcefield File:'

        self.cmapckb = Checkbutton(frame, text='No CMAP',
                    variable=self.nocmap)

        self.openffbtn = Button(
                frame,
                command=self.on_openff_clicked,
                text='Browse',
                width=10)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.parloc = CleanableEntryField(
                frame,
                labelpos='w',
                label_text='Parameter File:')

        self.openparbtn = Button(
                frame,
                command=self.on_openpar_clicked,
                text='Browse',
                width=10)

        # "Output" group
        # **********
        group = Pmw.Group(page, tag_text='Output Files')
        group.pack(**grp_opt)
        group.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.toploc = CleanableEntryField(
                frame,
                labelpos='w',
                label_text='Topology File:')

        self.savetopbtn = Button(
                frame,
                command=self.on_savetop_clicked,
                text='SaveAs',
                width=10)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.crdloc = CleanableEntryField(
                frame,
                labelpos='w',
                label_text='Coordinate File:')

        self.savecrdbtn = Button(
                frame,
                command=self.on_savecrd_clicked,
                text='SaveAs',
                width=10)

        entries = [self.pdbloc, self.ffloc, self.parloc,
                   self.toploc, self.crdloc]
        buttons = [self.openpdbbtn, self.openffbtn, self.openparbtn,
                   self.savetopbtn, self.savecrdbtn]
        for e, b in zip(entries, buttons):
            e.pack(**ent_opt)
            b.pack(**btn_opt)
        self.cmapckb.pack(**btn_opt)
        self.downloadbtn.pack(**btn_opt)
        Pmw.alignlabels(entries)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.topfmt = Pmw.OptionMenu(
                frame,
                labelpos='w',
                label_text='Topology Format:',
                items=('NAMD psf','AMBER prmtop'),
                menubutton_width=14)
        self.topfmt.pack(side='left', anchor='w', padx=10, pady=5)

        self.crdfmt = Pmw.OptionMenu(
                frame,
                labelpos='w',
                label_text='Coordinate Format:',
                items=('pdb', 'AMBER inpcrd'),
                menubutton_width=14)
        self.crdfmt.pack(side='right', anchor='w', padx=10, pady=5)

    def create_prep_page(self):
        page = self.notebook.add('Preparation')
        self.notebook.component('Preparation-tab').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        self.use_defrule = IntVar()
        self.use_userrule = IntVar()
        self.autodisu = IntVar()
        self.userdisu = IntVar()

        self.use_defrule.set(1)
        self.use_userrule.set(0)
        self.autodisu.set(0)
        self.userdisu.set(0)

        grp_opt = {'fill': 'both', 'expand': 1, 'padx': 10, 'pady': 5}
        #frm_opt = {'fill': 'both', 'expand': 1}
        frm_opt = {'fill': 'both', 'expand': 1, 'padx': 5}
        chk_opt = {'side': 'left', 'fill': 'x', 'expand': 0}
        btn_opt = {'side': 'left', 'fill': 'x',
                   'expand': 0, 'padx': 10, 'pady': 5}

        # "Rename"
        # **********
        group = Pmw.Group(page, tag_text='Rename')
        group.pack(**grp_opt)
        group.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        self.segname = CleanableEntryField(
            frame,
            labelpos='w',
            #validate={'validator': 'alphanumeric'},
            value='S1',
            label_text='Segment names:')
        self.segname.pack( anchor='w', fill='both', expand=1, padx=5, pady=5)

        Checkbutton(frame, text='Use default rules',
                    variable=self.use_defrule,
                    command=self.on_defrule_clicked).pack(**chk_opt)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        Checkbutton(frame, text='Specify a rule file:',
                    variable=self.use_userrule,
                    command=self.toggle_renloc_entry).pack(**chk_opt)
        self.renloc = CleanableEntryField(frame, entry_state='disabled')
        self.renloc.pack(side='left', fill='x', expand=1, padx=0, pady=5)
        self.openrenbtn = Button(
                frame,
                command=self.on_openren_clicked,
                text='Browse',
                width=10,
                state='disabled')
        self.openrenbtn.pack(**btn_opt)

        # "Disulfide Bond"
        # **********
        group = Pmw.Group(page, tag_text='Disulfide Bond')
        group.pack(**grp_opt)
        group.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        chk_opt['expand'] = 0
        Checkbutton(frame, text='Automatically detect with a cutoff of',
                    variable=self.autodisu,
                    command=self.toggle_disucut_entry).pack(**chk_opt)

        self.disucut = CleanableEntryField(
                frame,
                labelpos='e',
                validate={'validator': 'real', 'min': 0.1},
                value=2.1,
                label_text=u'\xc5',
                entry_state='disabled')
        self.disucut.pack(side='left', anchor='w', padx=0, pady=5)

        frame = Frame(group.interior())
        frame.pack(**frm_opt)

        Checkbutton(frame, text='Specify a bond file:',
                    variable=self.userdisu,
                    command=self.toggle_disuloc_entry).pack(**chk_opt)
        self.disuloc = CleanableEntryField(frame, entry_state='disabled')
        self.disuloc.pack(side='left', fill='x', expand=1, padx=0, pady=5)
        self.opendisubtn = Button(
                frame,
                command=self.on_opendisu_clicked,
                text='Browse',
                width=10,
                state='disabled')
        self.opendisubtn.pack(**btn_opt)

    def toggle_state(self, w):
        if w['state'] == 'normal':
            w.configure(state='disabled')
        else:
            w.configure(state='normal')
        w.update()

    def toggle_renloc_entry(self):
        self.renloc.setvalue('')
        self.toggle_state(self.renloc.component('entry'))
        self.toggle_state(self.openrenbtn)
        self.use_defrule.set(0)

    def on_defrule_clicked(self):
        self.use_userrule.set(0)
        self.renloc['entry_state'] = 'disabled'
        self.openrenbtn['state'] = 'disabled'

    def toggle_disucut_entry(self):
        self.disucut.setvalue(2.1)
        self.toggle_state(self.disucut.component('entry'))
        self.userdisu.set(0)
        self.disuloc['entry_state'] = 'disabled'
        self.opendisubtn['state'] = 'disabled'

    def toggle_disuloc_entry(self):
        self.disuloc.setvalue('')
        self.toggle_state(self.disuloc.component('entry'))
        self.toggle_state(self.opendisubtn)
        self.autodisu.set(0)
        self.disucut['entry_state'] = 'disabled'

    def create_sol_page(self):
        page = self.notebook.add('Solvation')
        self.notebook.component('Solvation-tab').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        grp_opt = {'fill': 'both', 'expand': 1, 'padx': 10, 'pady': 5}
        frm_opt = {'fill': 'both', 'expand': 1}
        btn_opt = {'side': 'left', 'fill': 'x',
                   'expand': 0, 'padx': 10, 'pady': 5}

        frame = Frame(page)
        frame.pack(side='left', **frm_opt)

        # "Solvents"
        # **********
        igroup = Pmw.Group(frame, tag_text='Solvents')
        igroup.pack(**grp_opt)
        igroup.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        w = Frame(igroup.interior())
        w.pack(fill='both', expand=1)

        self.watmod = Pmw.OptionMenu(
                w,
                labelpos='w',
                label_text='Solvent Model:',
                items=('TIP3P',),
                #items=('TIP3P', 'TIPS3P'),
                menubutton_width=12)
        self.watmod.pack(anchor='w', expand=1, padx=10, pady=5)

        w = Frame(igroup.interior())
        w.pack(fill='both', expand=1)

        self.watseg = CleanableEntryField(
                w,
                labelpos='w',
                validate={'validator': 'alphanumeric'},
                value='WAT',
                label_text='Segment Name:')
        self.watseg.pack(anchor='w', expand=1, padx=10, pady=5)


        w = Frame(igroup.interior())
        w.pack(fill='both', expand=1)

        self.cfloc = CleanableEntryField(
            w,
            labelpos='w',
            label_text='Coordinate File:')
        self.cfloc.pack(side='left', fill='x', expand=1, padx=10, pady=5)

        self.opencfbtn = Button(
            w,
            command=self.on_opencf_clicked,
            text='Browse',
            width=10).pack(**btn_opt)

        # " Box Parameters"
        # **********
        igroup = Pmw.Group(frame, tag_text='Box Parameters')
        igroup.pack(**grp_opt)
        igroup.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        w = Frame(igroup.interior())
        w.pack(fill='both', expand=1)

        self.pad = CleanableEntryField(
            w,
            labelpos='w',
            validate={'validator': 'real', 'min': 0.0},
            value=9.0,
            label_text=u'Buffer Distance (\xc5):')
        self.pad.pack(anchor='w', expand=1, padx=10, pady=5)

        w = Frame(igroup.interior())
        w.pack(fill='both', expand=1)

        self.cut = CleanableEntryField(
            w,
            labelpos='w',
            validate={'validator': 'real', 'min': 0.1},
            value=2.4,
            label_text=u'Overlap Cutoff (\xc5):')
        self.cut.pack(anchor='w', expand=1, padx=10, pady=5)

    def on_opencf_clicked(self):
        self.cfloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.crd ',
                    filetypes=[('Solvent Coordinate File', '.crd'),
                               ('All Files', '.*')]))

    def create_ion_page(self):
        page = self.notebook.add('Ionization')
        self.notebook.component('Ionization-tab').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        frm_opt = {'fill': 'both', 'expand': 1}
        ent_opt = {'anchor': 'w', 'padx': 10, 'pady': 5,
                   'fill': 'both', 'expand': 1}

        # "Ions"
        # **********
        igroup = Pmw.Group(page, tag_text='Ions')
        igroup.pack(side='left', fill='both', expand=1, padx=10, pady=5)
        igroup.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        self.catmod = Pmw.OptionMenu(
                igroup.interior(),
                labelpos='w',
                label_text='Cation Model:',
                menubutton_width=4,
                items=('Na+', 'K+', 'Mg2+', 'Ca2+', 'Zn2+'))
        self.catmod.pack(**ent_opt)

        self.catnum = Pmw.Counter(
                igroup.interior(),
                labelpos='w',
                label_text='Cation Number:',
                entry_width=4,
                entryfield_value=0,
                entry_state='disabled',
                datatype = {'counter': 'integer'},
                entryfield_validate={'validator': 'integer', 'min': '0'})
        self.catnum.pack(**ent_opt)

        self.animod = Pmw.OptionMenu(
                igroup.interior(),
                labelpos='w',
                label_text='Anion Model:',
                menubutton_width=4,
                items=('Cl-', ))
        self.animod.pack(**ent_opt)

        self.aninum = Pmw.Counter(
                igroup.interior(),
                labelpos='w',
                label_text='Anion Number:',
                entry_width=4,
                entryfield_value=0,
                entry_state='disabled',
                datatype = {'counter': 'integer'},
                entryfield_validate={'validator': 'integer', 'min': '0'})
        self.aninum.pack(**ent_opt)

        self.ionseg = CleanableEntryField(
                igroup.interior(),
                labelpos='w',
                label_text='Segment Name:',
                entry_width=10,
                value='ION',
                validate={'validator': 'alphanumeric'})
        self.ionseg.pack(**ent_opt)

        self.do_neutral = IntVar()
        self.do_neutral.set(1)

        self.calcqbtn = Button(
                igroup.interior(),
                command=self.on_calcq_clicked,
                text='Calculate the total charge')
        self.calcqbtn.pack(fill='both', expand=0, padx=10, pady=5)

        Checkbutton(igroup.interior(),
                    text='Automatically neutralize',
                    variable=self.do_neutral,
                    command=self.toggle_nions_salcon
                    ).pack(fill='both', expand=1, padx=10, pady=5)

        Pmw.alignlabels([self.catmod, self.catnum, self.animod, self.aninum,
                         self.ionseg])

        frame = Frame(page)
        frame.pack(side='right', **frm_opt)

        # "Methods"
        # **********
        igroup = Pmw.Group(frame, tag_text='Choose a method to place the ions')
        igroup.pack(fill='both', expand=1, padx=10, pady=5)
        igroup.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        self.ionmeth = IntVar()
        self.ionmeth.set(1)
        meths = [('by electrostatic potential', 1),
                 ('randomly', 2)]
        for txt, val in meths:
            Radiobutton(igroup.interior(),
                        text=txt,
                        padx=20,
                        variable=self.ionmeth,
                        value=val).pack(anchor='w', fill='y', expand=1)

        # "Parameters"
        # **********
        igroup = Pmw.Group(frame, tag_text='Parameters')
        igroup.pack(fill='both', expand=1, padx=10, pady=5)
        igroup.component('tag').configure(font=('Helvetica', 10, 'bold'), foreground='#336699')

        self.ionion = CleanableEntryField(
                igroup.interior(),
                labelpos='w',
                entry_width=5,
                validate={'validator': 'real', 'min': 0.1},
                value=5.0,
                label_text=u'Ion-to-Ion Cutoff (\xc5):')
        self.ionion.pack(**ent_opt)

        self.ionsol = CleanableEntryField(
                igroup.interior(),
                labelpos='w',
                entry_width=5,
                validate={'validator': 'real', 'min': 0.1},
                value=5.0,
                label_text=u'Ion-to-Solvent Cutoff (\xc5):')
        self.ionsol.pack(**ent_opt)

        self.salcon = CleanableEntryField(
                igroup.interior(),
                labelpos='w',
                entry_width=5,
                validate={'validator': 'real', 'min': 0.0},
                value=0.0,
                label_text='Salt Concentration (mol/L):')
        self.salcon.pack(**ent_opt)

        Pmw.alignlabels([self.ionion, self.ionsol, self.salcon])

    def on_calcq_clicked(self):
        if not _HAS_LIB or self.mod is None:
            return
        #q = fsum([atom.charge for atom in self.mod.atoms])
        #q = sum([atom.charge for atom in self.mod.atoms])
        atomList = self.topList[0]
        atomChrgList = [atom[6] for atom in atomList]
        preTotChrg = math.fsum(atomChrgList)
        tkMessageBox.showinfo('INFO', 'Total charge is %+f e'%preTotChrg,
                              parent=self.parent)

    def toggle_nions_salcon(self):
        for w in self.catnum, self.aninum, self.salcon:
            self.toggle_state(w.component('entry'))

    def on_dialog_button_clicked(self, result):
        if result == 'Execute':
            self.on_execute_button_clicked()
        elif result == 'Output':
            self.on_output_button_clicked()
        elif result == 'Console':
            self.on_console_button_clicked()
        elif result == 'About':
            self.on_about_button_clicked()
        else:
            self.on_quit_button_clicked()

    def load_input(self):
        # check
        if not _HAS_LIB:
            tkMessageBox.showerror(
                'ERROR',
                'Please install NumPy correcttly before launch the plugin',
                parent=self.parent)
            return 1

        if not self.pdbloc.getvalue():
            tkMessageBox.showerror('ERROR', 'Please specify a pdb file',
                                   parent=self.parent)
            return 1

        if not self.ffloc.getvalue():
            tkMessageBox.showerror('ERROR', 'Please specify a forcefield topology file',
                                   parent=self.parent)
            return 1

        if not self.parloc.getvalue():
            tkMessageBox.showerror('ERROR', 'Please specify a forcefield parameter file',
                                   parent=self.parent)
            return 1

        # check nocmap
        if self.nocmap.get() == 0:
            self.epreInpDict['NOCMAP'] = False
        elif self.nocmap.get() == 1:
            self.epreInpDict['NOCMAP'] = True

        print_logo()
        print("Job started on", time.ctime())
        self.startTime = clock()

        self.epreInpDict['TOPOLOGY'] = self.ffloc.getvalue()
        self.epreInpDict['FORCEFIELD'] = self.parloc.getvalue()
        self.epreInpDict['COORDPDB'] = self.pdbloc.getvalue()

        self.mod = read_pdb(self.pdbloc.getvalue(), self.epreInpDict)
        self.top = read_charmm_top(self.ffloc.getvalue(), self.epreInpDict)
        self.prm = read_charmm_prm(self.parloc.getvalue())
        return 0

    def on_execute_button_clicked(self):
        if not _HAS_LIB:
            tkMessageBox.showerror(
                'ERROR',
                'Please install NumPy correcttly before launch the plugin',
                parent=self.parent)
            return

        if self.notebook.getcurselection() == 'Preparation':
            if self.mod is None or self.top is None or self.prm is None:
                failed = self.load_input()
                if failed:
                    return
            try:
                # TODO
                splitLine = self.segname.getvalue().split()
                splitLineLen = len(splitLine)
                if splitLineLen != 0:
                    self.epreInpDict['SEGMENT'].extend(splitLine)

                if self.renloc['entry_state'] == 'normal':
                    read_rename_rule(self.renloc.getvalue(), self.epreInpDict)

                if self.disucut['entry_state'] == 'normal':
                    self.epreInpDict['DISUBOND']['DODISU'] = True
                    self.epreInpDict['DISUBOND']['AUTO'] = True
                    self.epreInpDict['DISUBOND']['CUT'] = float(self.disucut.getvalue())

                if  self.disuloc['entry_state'] == 'normal':
                    #self.epreInpDict['DISUBOND']['DODISU'] = True
                    read_bond_file(self.disuloc.getvalue(), self.epreInpDict)

                topDataList = self.top[1:]
                self.mod, self.topList = build_struct(topDataList, self.mod,
                                                      self.epreInpDict, self.prm)

            except Exception:
                tkMessageBox.showerror(
                    'ERROR',
                    'Failed',
                    parent=self.parent)
            else:
                tkMessageBox.showinfo(
                    'INFO',
                    'Successfully completed',
                    parent=self.parent)

        elif self.notebook.getcurselection() == 'Solvation':
            if self.mod is None or self.top is None or self.prm is None:
                failed = self.load_input()
                if failed:
                    return

            if not self.cfloc.getvalue():
                tkMessageBox.showerror('ERROR', 'Please specify a solvent coordinate file',
                                       parent=self.parent)
                return

            try:
                self.epreInpDict['ADDWAT']['DOADDWAT'] = True
                self.epreInpDict['ADDWAT']['MODEL'] = self.watmod.getvalue().upper()
                self.epreInpDict['ADDWAT']['SEGNAME'] = self.watseg.getvalue().upper()
                self.epreInpDict['ADDWAT']['COORDINATE'] = self.cfloc.getvalue()
                self.epreInpDict['ADDWAT']['CUT'] = float(self.cut.getvalue())
                self.epreInpDict['ADDWAT']['PAD'] = float(self.pad.getvalue())

                # if doAddWat
                self.sluPos = [
                    (atomPosX, atomPosY, atomPosZ)
                    for _, segDataList in self.mod
                    for _, _, resAtom in segDataList
                    for _, _, atomPosX, atomPosY, atomPosZ in resAtom
                    ]

                self.watPos, self.sluPos, self.boxInfo = add_wat(self.epreInpDict, self.sluPos)

                if sum(len(j[2]) for i in self.mod for j in i[1]) != len(self.sluPos):
                    perr("Numbers of atoms do not match")

                itPos = iter(self.sluPos)
                for _, segDataList in self.mod:
                    for _, _, resAtom in segDataList:
                        for atom in resAtom:
                            atom[2:] = itPos.next()

                # if doAddWat or doAddIon
                atomNumAdd = self.mod[-1][-1][-1][-1][-1][0]
                resNumAdd = self.mod[-1][-1][-1][0]
                ionPosList = []
                watDataList = build_solv_top(self.epreInpDict, self.watPos, ionPosList,
                                             atomNumAdd, resNumAdd)

                self.tmpmod = deepcopy(self.mod)
                self.tmptopList = deepcopy(self.topList)

                #self.mod += watDataList[0]
                #self.topList[0] += watDataList[1]
                #self.topList[1] += watDataList[2]
                #self.topList[2] += watDataList[3]
                self.mod = self.mod if watDataList[0] == [] else self.mod + watDataList[0]
                self.topList[0] = self.topList[0] if watDataList[1] == [] else self.topList[0] + watDataList[1]
                self.topList[1] = self.topList[1] if watDataList[2] == [] else self.topList[1] + watDataList[2]
                self.topList[2] = self.topList[2] if watDataList[3] == [] else self.topList[2] + watDataList[3]

            except Exception:
                tkMessageBox.showerror(
                    'ERROR',
                    'Failed',
                    parent=self.parent)
            else:
                tkMessageBox.showinfo(
                    'INFO',
                    'Successfully completed',
                    parent=self.parent)

        elif self.notebook.getcurselection() == 'Ionization':
            if self.mod is None or self.top is None or self.prm is None:
                failed = self.load_input()
                if failed:
                    return

            if self.do_neutral.get():
                catnum = aninum = 0
            else:
                catnum = int(self.catnum.getvalue())
                aninum = int(self.aninum.getvalue())
                if catnum == aninum == 0:
                    return

            self.epreInpDict['ADDION']['DOADDION'] = True
            catmod = IONS[self.catmod.getvalue()]
            animod = IONS[self.animod.getvalue()]
            self.epreInpDict['ADDION']['CATION'] = [catmod,catnum]
            self.epreInpDict['ADDION']['ANION'] = [animod,aninum]
            self.epreInpDict['ADDION']['SEGNAME'] = (self.ionseg.getvalue()).upper()
            self.epreInpDict['ADDION']['ION_SOLUTE'] = float(self.ionsol.getvalue())
            self.epreInpDict['ADDION']['ION_ION'] = float(self.ionion.getvalue())
            self.epreInpDict['ADDION']['SALTCON'] = float(self.salcon.getvalue())
            if self.ionmeth.get() == 1:
                self.epreInpDict['ADDION']['METHOD'] ='NORANDOM'
            elif self.ionmeth.get() == 2:
                self.epreInpDict['ADDION']['METHOD'] ='RANDOM'
                #raise NotImplementedError("This method hasn't been implemented yet.")
                tkMessageBox.showerror('ERROR', "This 'randomly' method hasn't been implemented yet.",
                                       parent=self.parent)
                return


            try:
                self.topList = self.tmptopList
                self.mod = self.tmpmod

                self.watPos, ionPosList = add_ion(self.epreInpDict, self.topList[0],
                                                  self.sluPos, self.watPos, self.boxInfo)

                # if doAddWat or doAddIon:
                atomNumAdd = self.mod[-1][-1][-1][-1][-1][0]
                resNumAdd = self.mod[-1][-1][-1][0]
                watDataList = build_solv_top(self.epreInpDict, self.watPos, ionPosList,
                                              atomNumAdd, resNumAdd)

                # self.mod += watDataList[0]
                # self.topList[0] += watDataList[1]
                # self.topList[1] += watDataList[2]
                # self.topList[2] += watDataList[3]
                self.mod = self.mod if watDataList[0] == [] else self.mod + watDataList[0]
                self.topList[0] = self.topList[0] if watDataList[1] == [] else self.topList[0] + watDataList[1]
                self.topList[1] = self.topList[1] if watDataList[2] == [] else self.topList[1] + watDataList[2]
                self.topList[2] = self.topList[2] if watDataList[3] == [] else self.topList[2] + watDataList[3]

            except Exception:
                tkMessageBox.showerror(
                    'ERROR',
                    'Failed',
                    parent=self.parent)
            else:
                tkMessageBox.showinfo(
                    'INFO',
                    'Successfully completed',
                    parent=self.parent)

        else:
            if self.mod is None or self.top is None or self.prm is None:
                self.load_input()
            return

        # generate a tmp PDB file for view
        objname = {'Preparation': 'modified', 'Solvation': 'solvated',
                   'Ionization': 'ionized'}[self.notebook.getcurselection()]
        tmpfp = StringIO()

        #write_tmppdb(tmpfp, self.mod, self.boxInfo)
        if self.boxInfo is not None:
            write_tmppdb(tmpfp, self.mod, self.boxInfo)
        else:
            write_tmppdb(tmpfp, self.mod)

        cmd.bg_color('white')
        for obj in self.pmobj:
            cmd.delete(obj)
        cmd.read_pdbstr(tmpfp.getvalue(), objname)
        util.cbag()
        self.pmobj = [objname]
        if self.notebook.getcurselection() == 'Ionization':
            cmd.show('spheres', 'segi %s'%self.epreInpDict['ADDION']['SEGNAME'])
            cmd.color('skyblue','segi %s'%self.epreInpDict['ADDION']['SEGNAME'])
        tmpfp.close()

    def on_output_button_clicked(self):
        # check
        if not _HAS_LIB:
            tkMessageBox.showerror(
                'ERROR',
                'Please install NumPy correcttly before launch the plugin',
                parent=self.parent)
            return

        if self.mod is None or self.top is None or self.prm is None:
            failed = self.load_input()
            if failed:
                return

        if not self.toploc.getvalue():
            tkMessageBox.showerror('ERROR', 'Please specify a topology file',
                                   parent=self.parent)
            return

        if not self.toploc.getvalue():
            tkMessageBox.showerror('ERROR', 'Please specify a coordinate file',
                                   parent=self.parent)
            return

        topTitles = self.top[0]
        save_files(self.epreInpDict, self.mod, self.topList, self.prm, self.boxInfo,
                   self.topfmt.getvalue(),self.toploc.getvalue(),
                   self.crdfmt.getvalue(),self.crdloc.getvalue(), topTitles)

        tkMessageBox.showinfo('INFO', '2 files were generated',
                              parent=self.parent)

        # if doAddWat or doAddIon
        mass = sum(atom[7] for atom in self.topList[0])
        boxMax, boxMin, boxLen = self.boxInfo
        volume = reduce(operator.mul, boxLen)
        center = map(lambda x, y: 0.5 * (x + y), boxMax, boxMin)
        conv = 10 / 6.0220449
        density = mass / volume * conv
        print("Total mass: %.3f amu" % mass)
        print("Unit cell volume: %.3f A^3" % volume)
        print("Unit cell center: (%.3f, %.3f, %.3f)" % tuple(center))
        print("Density: %.3f g/cm^3" % density)

        print("Time elapsed: %.3f s" % (clock() - self.startTime))
        print("Job finished on", time.ctime())

    def on_console_button_clicked(self):
        if self.console_shown:
            self.panedwin.forget(self.console_frame)
            self.console_shown = 0
        else:
            self.panedwin.add(self.console_frame)
            self.panedwin.paneconfigure(self.console_frame, padx=10, pady=10)
            self.console_shown = 1

    def on_about_button_clicked(self):
        about = Pmw.MessageDialog(self.parent, title='About the plugin',
                                  buttons=('Close',), defaultbutton=0,
                                  message_text=__doc__, message_justify='left')
        about.component('buttonbox').button(0).configure(width=10)
        about.activate(geometry='centerscreenfirst')

    def on_quit_button_clicked(self):
        for obj in self.pmobj:
            cmd.delete(obj)
        self.dialog.withdraw()

    def on_download_clicked(self):
        pdb = self.pdbloc.getvalue()
        try:
            download_pdb(pdb)
        except Exception:
            tkMessageBox.showerror('ERROR', 'Failed to download "%s"'%pdb,
                                   parent=self.parent)
        else:
            self.pdbloc.setvalue(os.path.join(os.getcwd(), pdb.upper()+'.pdb'))
            self.on_pdbentry_pressed()

    def on_pdbentry_pressed(self):
        pdb = self.pdbloc.getvalue()
        if self.check_exist(pdb) == Pmw.OK:
            cmd.load(pdb, 'original', format='pdb', quiet=0)
            util.cbag()
            self.pmobj.append('original')

    def on_openpdb_clicked(self, event=None):
        self.pdbloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.pdb .ent',
                    filetypes=[('PDB File', '.pdb .ent'),
                               ('All Files', '.*')]))
        self.on_pdbentry_pressed()

    def on_openff_clicked(self, event=None):
        self.ffloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.inp .top',
                    filetypes=[('Topology File', '.inp .top'),
                               ('All Files', '.*')]))
                    #filetypes=[('Forcefield File', '.inp .top'),
                    #          ('All Files', '.*')]))

    def on_openpar_clicked(self, event=None):
        self.parloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.prm',
                    filetypes=[('Parameter File', '.prm'),
                               ('All Files', '.*')]))

    def on_savetop_clicked(self, event=None):
        tops = {'AMBER prmtop': '.prmtop',
                'CHAMBER prmtop': '.prmtop',
                'GROMACS top': '.top',
                'NAMD psf': '.psf'}
        top = tops[self.topfmt.getvalue()]
        self.toploc.setvalue(
                tkFileDialog.asksaveasfilename(
                    defaultextension=top,
                    filetypes=[('Topology File', top),
                               ('All Files', '.*')]))

    def on_savecrd_clicked(self, event=None):
        crds = {'AMBER inpcrd': '.inpcrd',
                'GROMACS g96': '.g96',
                'GROMACS gro': '.gro',
                'NAMD bin': '.bin',
                'pdb': '.pdb'}
        crd = crds[self.crdfmt.getvalue()]
        self.crdloc.setvalue(
                tkFileDialog.asksaveasfilename(
                    defaultextension=crd,
                    filetypes=[('Coordinate File', crd),
                               ('All Files', '.*')]))

    def on_openren_clicked(self, event=None):
        self.renloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.in ',
                    filetypes=[('Rename Rule File', '.in'),
                               ('All Files', '.*')]))

    def on_opendisu_clicked(self, event=None):
        self.disuloc.setvalue(
                tkFileDialog.askopenfilename(
                    defaultextension='.in',
                    filetypes=[('Bond File', '.in'),
                               ('All Files', '.*')]))

    def check_exist(self, s):
        if not s:
            return Pmw.PARTIAL
        elif os.path.isfile(s):
            return Pmw.OK
        elif os.path.exists(s):
            return Pmw.PARTIAL
        else:
            return Pmw.PARTIAL

    def show_dialog(self):
        self.dialog.show()


def save_files(epreInpDict, crdDataList, topList, prmDataList, boxInfo, topfmt, topfile, crdfmt, crdfile, topTitles):
    if topfmt == 'NAMD psf':
        write_psf(topfile, topList, epreInpDict)
    elif topfmt == 'AMBER prmtop':
        crd = [
            (posX, posY, posZ)
            for seg in crdDataList
            for res in seg[1]
            for atomNum, _, posX, posY, posZ in res[2]
            ]
        box = boxInfo[2] + [90.0, 90.0, 90.0]
        write_inpcrd(topfile, crd, box)
    else:
        raise ValueError('Unsupported topology format' % topfmt)

    if crdfmt == 'pdb':
        write_pdb(crdfile, crdDataList, boxInfo)
    elif crdfmt == 'AMBER inpcrd':
        doAddWat = epreInpDict['ADDWAT']['DOADDWAT']
        write_prmtop(crdfile, topTitles, topList, prmDataList, isNPT=doAddWat)
    else:
        raise ValueError('Unsupported coordinate format' % crdfmt)

