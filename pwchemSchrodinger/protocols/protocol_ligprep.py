# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# General imports
import os, re, time, glob

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, BooleanParam, IntParam, PathParam,
                                        STEPS_PARALLEL)
from pyworkflow.utils.path import copyFile, moveFile, cleanPath
from pwem.protocols import EMProtocol

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseName, getBaseFileName

# Plugin imports
from .. import Plugin
from ..utils import saveMolecule

progLigPrep=Plugin.getHome('ligprep')
progStructConvert=Plugin.getHome('utilities/structconvert')

# Output attribute names
OUTPUTATTRIBUTE = "outputSmallMolecules"
OUTPUTATTRIBUTEDROPPED = "outputSmallMoleculesDropped"

# Input modes
INPUTSET, INPUTSMI = 0, 1
# Base name of the single output LigPrep writes when the whole SMILES file is prepared at once
SMIOUTPUTNAME = "ligPrepSmiles"

class ProtSchrodingerLigPrep(EMProtocol):
    """Schrodinger's LigPrep is a program to prepare ligand libraries
    """

    _label = 'ligand preparation (ligprep)'
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfSmallMolecules, OUTPUTATTRIBUTEDROPPED: SetOfSmallMolecules}
    stepsExecutionMode = STEPS_PARALLEL
    saving = False

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputFrom', EnumParam, default=INPUTSET,
                      choices=['Set of small molecules', 'SMILES file'],
                      label='Input ligands from:',
                      help='Prepare the molecules of a set coming from another protocol, or read them '
                           'directly from a SMILES file, which is handed to LigPrep with its -ismi option.')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set of small molecules:', allowsNull=True,
                       condition=f'inputFrom=={INPUTSET}')
        form.addParam('smiFile', PathParam, label='SMILES file:', condition=f'inputFrom=={INPUTSMI}',
                      help='File (.smi) with a SMILES string per line, optionally followed by the identifier '
                           'of the molecule, which LigPrep keeps as its name. The whole file is prepared by a '
                           'single LigPrep job, so no previous conversion of the molecules is needed.')
        group = form.addGroup('Ionization')
        group.addParam('ionization', EnumParam, default=0,
                       choices=["None",'Epik (recommended)','Do not neutralize or ionize',
                                'Neutralize only', 'Neutralize and ionize'],
                       label='Ionization')
        group.addParam('pH', FloatParam, default=7.0, condition='ionization!=0',
                       label='pH')
        group.addParam('pHrange', FloatParam, default=2.0, condition='ionization!=0',
                       label='pH range')
        group.addParam('emb', BooleanParam, default=True, condition='ionization==1',
                       label='Epik metal binding',
                       help='Run Epik with the metal_binding option so that states'
                            'appropriate for interactions with metal ions in'
                            'protein binding pockets are also generated.')

        group = form.addGroup('Stereoisomers')
        group.addParam('stereoisomers', BooleanParam, default=False,
                       label='Respect chirality in the input',
                       help='Do not respect existing chirality properties and do'
                            'not respect chiralities from the input geometry.'
                            'Generate stereoisomers for all chiral centers up to'
                            'the number permitted (specified using the -s option).'
                            'This is equivalent to "Generate all combinations" in'
                            'the Ligand Preparation user interface. Default'
                            'behavior is to respect only explicitly indicated'
                            'chiralities.')
        group.addParam('Niso', IntParam, default=32, condition='not stereoisomers',
                       label='Generate up to # isomers per input structure:')

        group = form.addGroup('Optimization')
        group.addParam('optimization', EnumParam, default=2,
                       choices=["None",'OPLS 2005','OPLS3e (recommended)'],
                       label='Force-field for final optimization')

        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        prepSteps = []
        if self.inputFrom.get() == INPUTSMI:
            prepSteps.append(self._insertFunctionStep(self.ligPrepSmilesStep, prerequisites=[]))
        else:
            for mol in self.inputSmallMolecules.get():
                pStep = self._insertFunctionStep(self.ligPrepStep, mol.clone(), prerequisites=[])
                prepSteps.append(pStep)
        self._insertFunctionStep(self.createOutputStep, prerequisites=prepSteps)


    def ligPrepStep(self, mol):
        args = self.getLigPrepArgs(mol)
        self.runJob(progLigPrep, args, cwd=self._getExtraPath())

    def ligPrepSmilesStep(self):
        args = self.getLigPrepSmilesArgs()
        self.runJob(progLigPrep, args, cwd=self._getExtraPath())

        molFile = self._getExtraPath(f"{getBaseName(mol.getFileName())}.sdf")
        if os.path.exists(molFile):
            for oFile in self.splitIsomers(molFile, self.getOutputMolsDir()):
                self.renameSDFTitle(oFile)

    def getSplitMolFiles(self, oDir):
        """Files split by the ligPrepSteps, indexed by molecule name. structconvert names them
        "<molName>.sdf" when a single structure was prepared and "<molName>-<confId>.sdf" if more.
        The whole directory is read once, it holds a file per structure of the input set."""
        molFiles = {}
        for oFile in sorted(os.listdir(oDir)):
            if oFile.endswith('.sdf'):
                molName = re.sub(r"-\d+$", '', oFile[:-len('.sdf')])
                molFiles.setdefault(molName, []).append(os.path.join(oDir, oFile))
        return molFiles

    def getOutputMolsDir(self):
        oDir = self._getPath('outputMols')
        # exist_ok instead of makePath, several ligPrepSteps may create it at the same time
        os.makedirs(oDir, exist_ok=True)
        return oDir

    def createOutputStep(self):
        splitMolFiles = self.getSplitMolFiles(self.getOutputMolsDir())

        self.outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
        if self.inputFrom.get() == INPUTSMI:
            self.createOutputFromSmiles(oDir)
            return

        outputSmallMoleculesDropped = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='Dropped')
        for mol in self.inputSmallMolecules.get():
            molName = getBaseName(mol.getFileName())
            molFile = self._getExtraPath(f"{molName}.sdf")

            if not os.path.exists(molFile):
                droppedFile = self._getExtraPath(f"{molName}-dropped.sdf")
                if os.path.exists(droppedFile):
                    # LigPrep ran, but discarded the molecule at one of its processing stages.
                    print(f"LigPrep dropped molecule '{molName}', "
                          f"see '{molName}.log' to see why.")
                    # saveMolecule is not used here: it reads a conformer id from the file name,
                    # and the "-dropped" suffix is not one.
                    droppedMol = SmallMolecule()
                    droppedMol.copy(mol, copyId=False)
                    droppedMol.setFileName(droppedFile)
                    outputSmallMoleculesDropped.append(droppedMol.clone())
                else:
                    # No output at all: the LigPrep job failed
                    print(f"LigPrep failed for molecule '{molName}'")
                continue

            # Already split and renamed by its ligPrepStep
            for oFile in splitMolFiles.get(molName, []):
                saveMolecule(self, oFile, self.outputSmallMolecules, mol)

        outputsDict = {OUTPUTATTRIBUTE: self.outputSmallMolecules}
        if len(outputSmallMoleculesDropped) > 0:
            outputsDict[OUTPUTATTRIBUTEDROPPED] = outputSmallMoleculesDropped

        self._defineOutputs(**outputsDict)
        self._defineSourceRelation(self.inputSmallMolecules, self.outputSmallMolecules)
        if len(outputSmallMoleculesDropped) > 0:
            self._defineSourceRelation(self.inputSmallMolecules, outputSmallMoleculesDropped)


    def getLigPrepArgs(self, mol):
        fnSmall = os.path.abspath(mol.getFileName())
        fnBase = getBaseFileName(fnSmall)
        fnRoot = getBaseName(fnSmall)

        args = '-WAIT -LOCAL'
        if self.ionization.get() != 0:
            if self.ionization.get() == 1:
                args += " -epik"
                if self.emb.get():
                    args += " -epik_metal_binding"
            else:
                args += " -i %d" % self.ionization.get() - 2
            args += " -ph %f -pht %f" % (self.pH.get(), self.pHrange.get())

        if self.stereoisomers.get():
            args += " -g"
        else:
            args += " -ac -s %d" % self.Niso.get()

        if self.optimization.get() == 1:
            args += " -bff 14"
        elif self.optimization.get() == 2:
            args += " -bff 16"

        if fnBase.endswith('.smi'):
            args += f" -ismi {fnSmall}"
        elif fnBase.endswith('.mae') or fnBase.endswith('.maegz'):
            args += f" -imae {fnSmall}"
        elif fnBase.endswith('.sdf'):
            args += f" -isd {fnSmall}"
        else:
            fnSDF = self._getTmpPath(fnRoot + '.sdf')
            self.runJob(progStructConvert, '{} {}'.format(fnSmall, fnSDF))
            mol.setFileName(fnSDF)

            args += " -isd ../tmp/%s" % (fnRoot + '.sdf')

        fnSDF = "%s.sdf" % fnRoot
        args += " -osd %s" % fnSDF
        return args

    def renameSDFTitle(self, sdfFile):
        tmpFile = self._getTmpPath(os.path.basename(sdfFile))
        with open(sdfFile) as fIn:
            fIn.readline()
            with open(tmpFile, 'w') as fOut:
                fOut.write(f'{getBaseName(sdfFile)}\n')
                for line in fIn:
                    fOut.write(line)
        os.rename(tmpFile, sdfFile)

    def splitIsomers(self, molFile, oDir):
        prefix = 'ligPrepConfs_'

        fnRoot = getBaseName(molFile)
        fnOsdf = f"{prefix}{fnRoot}.sdf"
        args = "%s %s -split-nstructures 1" % (getBaseFileName(molFile), fnOsdf)
        self.runJob(progStructConvert, args, cwd=self._getExtraPath())

        outMols = []
        for fn in glob.glob(self._getExtraPath(f"{prefix}{fnRoot}*.sdf")):
            fnOut = os.path.join(oDir, getBaseFileName(fn).replace(prefix, ''))
            moveFile(fn, fnOut)
            outMols.append(fnOut)

        return outMols

