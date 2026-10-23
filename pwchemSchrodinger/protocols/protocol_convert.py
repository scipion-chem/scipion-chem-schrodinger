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
import os

# Scipion em imports
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, STEPS_PARALLEL

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import getBaseName

# Plugin imports
from .. import Plugin
from ..objects import SchrodingerAtomStruct
from ..utils.utils import putMolFileTitle, saveMolecule

SMALLMOL, TARGET = 0, 1
molChoices = {"Maestro": 'maegz', 'PDB': 'pdb', "Sybyl Mol2": 'mol2', "Smiles": 'smi', 'V2000 SD': 'sdf'}
targetChoices = {"Maestro": 'maegz', 'PDB': 'pdb'}

class ProtSchrodingerConvert(EMProtocol):
    """Convert a set of input ligands or a receptor structure to a specific file format"""
    _label = 'convert'

    saving = False

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        # Input comparison variable
        inputTypeComparison = 'inputType=='

        form.addSection(label='Input')
        form.addParam('inputType', EnumParam, default=0, choices=["Small molecules", 'Target structure'],
                      label='Input type', help='Type of input you want to convert')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      condition='{}{}'.format(inputTypeComparison, SMALLMOL), label='Input small molecules:',
                      help='Input small molecules to convert')
        form.addParam('outputFormatSmall', EnumParam, default=0, condition='{}{}'.format(inputTypeComparison, SMALLMOL),
                      choices=list(molChoices.keys()), label='Output format',
                      help='Output format for the small molecules')
        form.addParam('usePose', BooleanParam, default=False,
                      condition='{}{}'.format(inputTypeComparison, SMALLMOL), label='Use the docked ligands: ',
                      help='Convert the docked poses instead of the original ligand files. Their receptors are '
                           'also converted to PDB (one per pose if each pose has its own, as in IFD)')
        form.addParam('inputStructure', PointerParam, pointerClass="SchrodingerAtomStruct, AtomStruct",
                      condition='{}{}'.format(inputTypeComparison, TARGET), label='Input structure:',
                      help='Input atomic structure to convert')
        form.addParam('outputFormatTarget', EnumParam, default=0, condition='{}{}'.format(inputTypeComparison, TARGET),
                      choices=list(targetChoices.keys()), label='Output format',
                      help='Output format for the atomic structure')

        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        convSteps = []
        if self.inputType.get() == SMALLMOL:
            for mol in self.inputSmallMolecules.get():
                cStep = self._insertFunctionStep('convertMolStep', mol.clone(), prerequisites=[])
                convSteps.append(cStep)
            # Docked sets: receptors go to PDB, one step per distinct file (IFD has a receptor per pose)
            for recFile in self.getReceptorFiles():
                convSteps.append(self._insertFunctionStep('convertReceptorStep', recFile, prerequisites=[]))

        elif self.inputType.get() == TARGET:
            cStep = self._insertFunctionStep('convertTargetStep', prerequisites=[])
            convSteps.append(cStep)

        self._insertFunctionStep('createOutputStep', prerequisites=convSteps)

    def usesPoses(self):
        return self.inputType.get() == SMALLMOL and self.usePose.get()

    def getMolInFile(self, mol):
        '''With usePose the docked pose is converted, not the original input file'''
        return mol.getPoseFile() if self.usesPoses() else mol.smallMoleculeFile.get()

    def getOutFile(self, inFile, outFormat):
        '''Relative to the project, as Scipion stores file paths'''
        return self._getExtraPath(f'{getBaseName(inFile)}.{outFormat}')

    def getReceptorFiles(self):
        '''Distinct receptor files of a docked set that are not PDB yet: the set one and, if any, the per pose ones'''
        if not self.usesPoses():
            return []
        inSet = self.inputSmallMolecules.get()
        recFiles = {os.path.abspath(inSet.getOriginalReceptorFile())}
        recFiles |= {os.path.abspath(mol.getProteinFile()) for mol in inSet if mol.getProteinFile()}
        return sorted(f for f in recFiles if not f.endswith('.pdb'))

    def getPDBReceptor(self, recFile):
        return os.path.relpath(recFile) if recFile.endswith('.pdb') else self.getOutFile(recFile, 'pdb')

    def convertMolStep(self, mol):
        progStructConvert=Plugin.getHome('utilities/structconvert')

        fnSmall = self.getMolInFile(mol)
        outFormat = molChoices[self.getEnumText('outputFormatSmall')]
        fnOut = self.getOutFile(fnSmall, outFormat)

        self.runJob(progStructConvert, '{} {}'.format(os.path.abspath(fnSmall), os.path.abspath(fnOut)))
        if outFormat == 'mol2':
            putMolFileTitle(fnOut, ext=outFormat)

    def convertReceptorStep(self, recFile):
        progStructConvert = Plugin.getHome('utilities/structconvert')
        self.runJob(progStructConvert, '{} {}'.format(recFile, os.path.abspath(self.getPDBReceptor(recFile))))

    def convertTargetStep(self):
        progStructConvert = Plugin.getHome('utilities/structconvert')
        fnStructure = self.inputStructure.get().getFileName()
        fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]

        outFormat = targetChoices[self.getEnumText('outputFormatTarget')]
        fnOut = self._getExtraPath("{}.{}".format(fnRoot, outFormat))

        args = '{} {}'.format(fnStructure, fnOut)

        self.runJob(progStructConvert, args)

        if fnOut.endswith('.maegz'):
            self.target = SchrodingerAtomStruct(filename=fnOut)
        else:
            self.target = AtomStruct(filename=fnOut)

    def createOutputStep(self):
        if self.inputType.get() == SMALLMOL:
            inSet = self.inputSmallMolecules.get()
            outFormat = molChoices[self.getEnumText('outputFormatSmall')]
            outSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')

            if self.usesPoses():
                outSet.setDocked(True)
                outSet.setProteinFile(self.getPDBReceptor(inSet.getOriginalReceptorFile()))
                for mol in inSet:
                    newMol = mol.clone()
                    newMol.setPoseFile(self.getOutFile(mol.getPoseFile(), outFormat))
                    if mol.getProteinFile():
                        newMol.setProteinFile(self.getPDBReceptor(mol.getProteinFile()))
                    if outFormat != 'maegz':
                        # Not Maestro anymore, so viewers other than ChimeraX (e.g. PyMOL) can open it
                        newMol.setMolClass('Standard')
                    outSet.append(newMol)
                outSet.updateMolClass()
            else:
                for mol in inSet:
                    saveMolecule(self, self.getOutFile(mol.smallMoleculeFile.get(), outFormat), outSet, mol)

            self._defineOutputs(outputSmallMolecules=outSet)
            self._defineSourceRelation(self.inputSmallMolecules, outSet)
        else:
            self._defineOutputs(outputStructure=self.target)

    def _validate(self):
        errors = []
        if self.usesPoses() and not self.inputSmallMolecules.get().isDocked():
            errors.append('"Use the docked ligands" needs a docked set of small molecules as input')
        return errors

    def _summary(self):
        summary=[]
        if self.inputType.get()==0:
            if self.outputFormatSmall.get() == 0:
                summary.append('Converted to Maestro')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to PDB')
            elif self.outputFormatSmall.get() == 2:
                summary.append('Converted to Sybyl Mol2')
            elif self.outputFormatSmall.get() == 3:
                summary.append('Converted to Smiles')
            elif self.outputFormatSmall.get() == 4:
                summary.append('Converted to V2000 SD')
        else:
            if self.outputFormatTarget.get() == 0:
                summary.append('Converted to Maestro')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to PDB')
        return summary