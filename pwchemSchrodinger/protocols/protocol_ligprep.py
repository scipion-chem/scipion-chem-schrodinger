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
import os, time, glob

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, FloatParam, BooleanParam, IntParam, STEPS_PARALLEL
from pyworkflow.utils.path import copyFile, moveFile, cleanPath
from pwem.protocols import EMProtocol

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseName

# Plugin imports
from .. import Plugin
from ..utils import saveMolecule

progLigPrep=Plugin.getHome('ligprep')
progStructConvert=Plugin.getHome('utilities/structconvert')

# Output attribute names
OUTPUTATTRIBUTE = "outputSmallMolecules"
OUTPUTATTRIBUTEDROPPED = "outputSmallMoleculesDropped"

class ProtSchrodingerLigPrep(EMProtocol):
    """Schrodinger's LigPrep is a program to prepare ligand libraries"""
    _label = 'ligand preparation (ligprep)'
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfSmallMolecules, OUTPUTATTRIBUTEDROPPED: SetOfSmallMolecules}
    saving = False

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set of small molecules:', allowsNull=False)
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
        iniStep = self._insertFunctionStep('initializeStep')
        prepSteps=[]
        for mol in self.inputSmallMolecules.get():
            pStep = self._insertFunctionStep('ligPrepStep', mol.clone(), prerequisites=[iniStep])
            prepSteps.append(pStep)
        self._insertFunctionStep('createOutputStep', prerequisites=prepSteps)

    def initializeStep(self):
        self.outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
        self.outputSmallMoleculesDropped = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMolsDropped')

    def ligPrepStep(self, mol):
        fnSmall = mol.smallMoleculeFile.get()
        fnMol = os.path.split(fnSmall)[1]
        fnRoot = os.path.splitext(fnMol)[0]

        existingFiles = glob.glob(self._getExtraPath(fnRoot+"*"))
        if len(existingFiles) == 0:
            fnSmallExtra = self._getTmpPath(fnMol)
            copyFile(fnSmall, fnSmallExtra)

            args='-WAIT -LOCAL'
            if self.ionization.get() != 0:
                if self.ionization.get() == 1:
                    args+=" -epik"
                    if self.emb.get():
                        args+=" -epik_metal_binding"
                else:
                    args+=" -i %d"%self.ionization.get()-2
                args+=" -ph %f -pht %f"%(self.pH.get(),self.pHrange.get())

            if self.stereoisomers.get():
                args+=" -g"
            else:
                args+=" -ac -s %d"%self.Niso.get()

            if self.optimization.get() == 1:
                args+=" -bff 14"
            elif self.optimization.get() == 2:
                args+=" -bff 16"

            if fnMol.endswith('.smi'):
                args+=" -ismi tmp/%s" % (fnMol)
            elif fnMol.endswith('.mae') or fnMol.endswith('.maegz'):
                args += " -imae tmp/%s" % (fnMol)
            elif fnMol.endswith('.sdf'):
                args += " -isd tmp/%s" % (fnMol)
            else:
                fnSDF = self._getTmpPath(fnRoot + '.sdf')
                self.runJob(progStructConvert, '{} {}'.format(fnSmall, fnSDF))
                mol.smallMoleculeFile.set(fnSDF)

                args += " -isd tmp/%s" % (fnRoot + '.sdf')

            fnSDF = "extra/%s.sdf" % fnRoot
            if not os.path.exists(fnSDF):
                args+=" -osd %s"%fnSDF
                self.runJob(progLigPrep, args, cwd=self._getPath())

            if os.path.exists(self._getPath(fnSDF)):
                fnOsdf="extra/o%s.sdf"%fnRoot
                args = "%s %s -split-nstructures 1" % (fnSDF, fnOsdf)
                self.runJob(progStructConvert, args, cwd=self._getPath())
                for fn in glob.glob(self._getExtraPath("o%s*.sdf" % fnRoot)):
                    fnOut = os.path.split(fn)[1]
                    fnOut = self._getExtraPath(fnOut[1:])
                    moveFile(fn, fnOut)
                    saveMolecule(self, fnOut, self.outputSmallMolecules, mol)
                if len(glob.glob(self._getExtraPath("%s-*.sdf" % fnRoot))) > 0:
                    cleanPath(self._getPath(fnSDF))
            else:
                saveMolecule(self, fnSmall, self.outputSmallMoleculesDropped, mol)
        else:
            for fn in glob.glob(self._getExtraPath("%s*.sdf" % fnRoot)):
                print("Reading %s" % fn)
                saveMolecule(self, fn, self.outputSmallMolecules, mol)

    def createOutputStep(self):
        if len(self.outputSmallMolecules)>0:
            self._defineOutputs(**{OUTPUTATTRIBUTE: self.outputSmallMolecules})
            self._defineSourceRelation(self.inputSmallMolecules, self.outputSmallMolecules)
        if len(self.outputSmallMoleculesDropped)>0:
            self._defineOutputs(**{OUTPUTATTRIBUTEDROPPED: self.outputSmallMoleculesDropped})
            self._defineSourceRelation(self.inputSmallMolecules, self.outputSmallMoleculesDropped)

    def renameSDFTitle(self, sdfFile):
        tmpFile = self._getTmpPath(os.path.basename(sdfFile))
        with open(sdfFile) as fIn:
            fIn.readline()
            with open(tmpFile, 'w') as fOut:
                fOut.write(f'{getBaseName(sdfFile)}\n')
                for line in fIn:
                    fOut.write(line)
        os.rename(tmpFile, sdfFile)
