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
from pwchem.utils import getBaseName, getBaseFileName

# Plugin imports
from .. import Plugin
from ..utils import saveMolecule

progLigPrep=Plugin.getHome('ligprep')
progStructConvert=Plugin.getHome('utilities/structconvert')

# Output attribute names
OUTPUTATTRIBUTE = "outputSmallMolecules"
OUTPUTATTRIBUTEDROPPED = "outputSmallMoleculesDropped"

class ProtSchrodingerLigPrep(EMProtocol):
    """Schrodinger's LigPrep is a program to prepare ligand libraries

    User documentation(AI GENERATED):
    The Ligand Preparation (LigPrep) protocol in Schrodinger is designed to prepare small molecules for docking simulations.
    This process involves several key steps, including ionization, stereoisomer generation, and force field optimization.
    The protocol takes as input a set of small molecules and prepares them by performing necessary transformations and optimizations
    to ensure that they are suitable for use in molecular docking or other simulations.
    The first step involves defining the parameters for ionization. The user can choose from several ionization methods, such as neutralization,
    Epik ionization (which is recommended), or no ionization. The pH value and pH range can also be specified, with default values typically used
    for standard conditions. If Epik ionization is selected, users have the option to enable metal binding, which is important for simulating interactions
    with metal ions in protein binding pockets.

    Another key feature of this protocol is the stereoisomer generation. The user can choose whether to respect the chirality of the input molecule or to generate all
    possible stereoisomers for the chiral centers in the molecule. If stereoisomers are not generated, the protocol respects the existing chirality in the input molecule.
    The user can also specify the number of isomers to generate if chirality is not respected.

    For optimization, the protocol offers different force fields to choose from. The user can select between no optimization, the OPLS 2005 force field, or the recommended OPLS3e
    force field, which is used for final optimization to ensure that the ligand is in its lowest energy conformation.
    Once the parameters are defined, the LigPrep protocol prepares the input molecules by applying these settings. For each molecule, the preparation process involves checking
    if the molecule has already been processed. If not, the ligand is ionized (if applicable), stereoisomers are generated, and the molecule is optimized using the chosen force field.
    The molecules are then converted into a compatible structure format, such as .sdf, and saved for further analysis.

    The prepared molecules are saved into two categories: the molecules that were successfully prepared and those that were dropped during the process.
    If molecules are successfully prepared, they are saved in the output set, ready for further analysis or docking. The output also includes a summary of the ligand preparation steps,
    allowing users to review the parameters used during the process. If any molecules fail during preparation, they are saved separately, and the user can review why they were dropped.
    In summary, the LigPrep protocol is a comprehensive tool for preparing small molecules for docking simulations, ensuring that they are in an appropriate form for accurate and reliable
    docking results. The ability to customize ionization, stereoisomer generation, and optimization parameters makes it a powerful tool for a variety of molecular modeling tasks in drug
    discovery.
    """

    _label = 'ligand preparation (ligprep)'
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfSmallMolecules, OUTPUTATTRIBUTEDROPPED: SetOfSmallMolecules}
    stepsExecutionMode = STEPS_PARALLEL
    saving = False

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
        prepSteps = []
        for mol in self.inputSmallMolecules.get():
            pStep = self._insertFunctionStep(self.ligPrepStep, mol.clone(), prerequisites=[])
            prepSteps.append(pStep)
        self._insertFunctionStep(self.createOutputStep, prerequisites=prepSteps)


    def ligPrepStep(self, mol):
        args = self.getLigPrepArgs(mol)
        self.runJob(progLigPrep, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        oDir = self._getPath('outputMols')
        if not os.path.exists(oDir):
            os.mkdir(oDir)

        self.outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
        for mol in self.inputSmallMolecules.get():
            molFile = mol.getFileName()
            molFile = self._getExtraPath(f"{getBaseName(molFile)}.sdf")
            outMolFiles = self.splitIsomers(molFile, oDir)

            for oFile in outMolFiles:
                self.renameSDFTitle(oFile)
                saveMolecule(self, oFile, self.outputSmallMolecules, mol)

        self._defineOutputs(**{OUTPUTATTRIBUTE: self.outputSmallMolecules})
        self._defineSourceRelation(self.inputSmallMolecules, self.outputSmallMolecules)


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

