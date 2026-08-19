# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, re, subprocess
import pwem.objects.data as data
from pwchem.constants import FUNCTION_BOUNDING_BOX, PML_BBOX_STR_EACH, PML_BBOX_STR
from pyworkflow.object import (Float, Integer, List, String)
from . import Plugin as schrodingerPlugin

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')

class SchrodingerAtomStruct(data.AtomStruct):
    """An AtomStruct in the file format of Maestro"""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def getExtension(self):
        return os.path.splitext(self.getFileName())[1]

    def convert2(self, oFile=None, cwd=None, ext='.cif'):
        if not oFile:
            oFile = os.path.abspath(self.getFileName().replace(self.getExtension(), ext))
        command = '{} {} {}'.format(structConvertProg, os.path.abspath(self.getFileName()), oFile)
        if oFile.endswith('.cif'):
            command += ' -PDBx'

        subprocess.check_call(command, shell=True, cwd=cwd)
        return oFile

    def convert2PDB(self, outPDB=None, cwd=None):
        outPDB = self.convert2(outPDB, cwd, ext='.pdb')
        return outPDB


class SchrodingerSystem(data.EMFile):
    """An system atom structure (prepared for MD) in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

    def getExtension(self):
        return os.path.splitext(self.getFileName())[1]

    def getDirName(self):
        return os.path.dirname(self.getFileName())

    def getBaseName(self):
        return os.path.basename(os.path.splitext(self.getFileName())[0])

    def getTrajectoryDirName(self):
        with open(self.getFileName()) as fCMS:
            cmsSTR = fCMS.read()
            trDirs = re.findall(r'[\w-]*_trj', cmsSTR)
        return trDirs[0]

    def getCMSFileName(self):
        with open(self.getFileName()) as fCMS:
            cmsSTR = fCMS.read()
            fName = re.findall(r'[\w-]*\.cms', cmsSTR)
        return fName[0]

    def changeTrajectoryDirName(self, newDirPath, trjPath=None):
        '''Change the name of the trajectory directory specified in the CMS file'''
        dirName = self.getTrajectoryDirName()
        if trjPath==None:
            trjPath = self.getDirName()
        dirPath = os.path.join(trjPath, dirName)

        newDirName = newDirPath.split('/')[-1]
        with open(self.getFileName()) as fCMS:
            cmsSTR = fCMS.read()
        cmsSTR = cmsSTR.replace(dirName, newDirName)

        with open(self.getFileName(), 'w') as fCMS:
            fCMS.write(cmsSTR)

        os.rename(dirPath, newDirPath)

    def changeCMSFileName(self, newCMSPath):
        '''The name of the CMS is found inside itself. Change it in the file along with the file itself'''
        newCMSName = newCMSPath.split('/')[-1]
        fName = self.getCMSFileName()
        with open(self.getFileName()) as fCMS:
            cmsSTR = fCMS.read()
        cmsSTR = cmsSTR.replace(fName, newCMSName)

        with open(self.getFileName(), 'w') as fCMS:
            fCMS.write(cmsSTR)

        os.rename(self.getFileName(), newCMSPath)
        self.setFileName(newCMSPath)



class SchrodingerGrid(data.EMFile):
    """A search grid in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)
        self._centerX = Float(kwargs.get('centerX', None))
        self._centerY = Float(kwargs.get('centerY', None))
        self._centerZ = Float(kwargs.get('centerZ', None))
        self._innerX = Integer(kwargs.get('innerX', None))
        self._innerY = Integer(kwargs.get('innerY', None))
        self._innerZ = Integer(kwargs.get('innerZ', None))
        self._outerX = Integer(kwargs.get('outerX', None))
        self._outerY = Integer(kwargs.get('outerY', None))
        self._outerZ = Integer(kwargs.get('outerZ', None))

        self._forceField = String(kwargs.get('forceField', None))
        self._proteinFile = String(kwargs.get('proteinFile', None))

    def __str__(self):
      s = '{} (Center: {})'.format(self.getClassName(), self.getCenter())
      return s

    def getCenter(self):
        return self._centerX.get(), self._centerY.get(), self._centerZ.get()

    def getInnerBox(self):
        return self._innerX.get(), self._innerY.get(), self._innerZ.get()

    def getOuterBox(self):
        return self._outerX.get(), self._outerY.get(), self._outerZ.get()

    def getProteinFile(self):
        return self._proteinFile.get()


class SetOfSchrodingerGrids(data.EMSet):
    ITEM_TYPE = SchrodingerGrid

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)

    def __str__(self):
      s = '{} ({} items)'.format(self.getClassName(), self.getSize())
      return s

    def getSetPath(self):
        return os.path.abspath(self._mapperPath[0])

    def getSetDir(self):
        return '/'.join(self.getSetPath().split('/')[:-1])

    def getProteinFile(self):
        return self.getFirstItem().getProteinFile()

    def getBBoxPml(self):
        return self.getSetDir() + '/BBoxes.pml'

    def buildBBoxesPML(self):
        pmlFile = self.getBBoxPml()
        toWrite = FUNCTION_BOUNDING_BOX
        for grid in self:
            toWrite += PML_BBOX_STR_EACH.format([0, 1, 0], grid.getCenter(), grid.getInnerBox(),
                                                'InnerBox_'+str(grid.getObjId()))
            toWrite += PML_BBOX_STR_EACH.format([1, 0, 1], grid.getCenter(), grid.getOuterBox(),
                                                'OuterBox_' + str(grid.getObjId()))
        with open(pmlFile, 'w') as f:
            f.write(PML_BBOX_STR.format(os.path.abspath(self.getProteinFile()), toWrite))



class SchrodingerBindingSites(data.EMFile):
    """A set of binding sites in the file format of Maestro"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)


class SchrodingerFEPResult(data.EMFile):
    """Wraps the raw output of an alchemical free-energy (FEP+) job - the ".fmp"/pose-
    viewer/log file(s) fep_plus (RBFE) or the ABFE driver produce for one ligand pair
    (RBFE edge) or one ligand (ABFE) - plus, if parsing succeeded, the resulting free
    energy. Modelled on gromacs.objects.GromacsSystem's own
    setFreeEnergy/setFreeEnergyFile pair, so downstream code that already knows that
    convention (e.g. a viewer, or a follow-up protocol) can read either object the same
    way. See claude/decisions/schrodinger/fep_plus.md for why the free energy may be
    None here far more often than in the GROMACS case: this plugin's FEP+ log parser is
    unverified (never checked against a real FEP+ log), so the raw file(s) - not the
    parsed number - are the trustworthy part of this object."""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)
        self._freeEnergy = Float(kwargs.get('freeEnergy', None))
        self._freeEnergyFile = String(kwargs.get('freeEnergyFile', None))
        self._ligandA = String(kwargs.get('ligandA', None))
        self._ligandB = String(kwargs.get('ligandB', None))  # unset for ABFE (one ligand only)

    def __str__(self):
        label = f'{self._ligandA.get()}'
        if self._ligandB.get():
            label += f' -> {self._ligandB.get()}'
        dG = self.getFreeEnergy()
        return f'{self.getClassName()} ({label}{f", dG={dG:.2f} kcal/mol" if dG is not None else ""})'

    def setFreeEnergy(self, value):
        self._freeEnergy.set(value)

    def getFreeEnergy(self):
        return self._freeEnergy.get()

    def setFreeEnergyFile(self, fn):
        self._freeEnergyFile.set(fn)

    def getFreeEnergyFile(self):
        return self._freeEnergyFile.get()

    def setLigands(self, ligandA, ligandB=None):
        self._ligandA.set(ligandA)
        if ligandB is not None:
            self._ligandB.set(ligandB)

    def getLigands(self):
        return self._ligandA.get(), self._ligandB.get()


class SchrodingerQSARModel(data.EMObject):
    """Object to store a Phase QSAR model and its metadata."""

    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self.qsarModel = String()
        self.projectPath = String() #used only with pharmacophore models
        self.hypoFile = String() #used only with pharmacophore models
        self.modelFile = String()
        self.summaryFile = String()
        self.predictionsFile = String()
        self.molFile = String()
        self.fieldFile = String()

        self.style = String()
        self.forceField = String()
        self.trainFraction = Float()
        self.lno = Integer()

    def setModelFile(self, fn):
        self.modelFile.set(fn)

    def getModelFile(self):
        return self.modelFile.get()



