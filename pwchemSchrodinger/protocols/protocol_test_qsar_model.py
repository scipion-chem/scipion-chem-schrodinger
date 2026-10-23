# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import csv
import requests, pandas as pd
import math
import time
import os, subprocess, shutil

# Scipion em imports
from pwem.objects.data import AtomStruct
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL, FloatParam, FileParam, StringParam, IntParam, BooleanParam

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pyworkflow.object import Float

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin
from ..utils import createSdf


class ProtSchrodingerQSARTest(EMProtocol):
    """Test QSAR models"""
    _label = 'test QSAR model'
    filters = ['stereo', 'connect', 'distinct']
    scriptsDir = os.path.join(os.path.dirname(__file__), "../scripts")
    pharmCsv = 'qsar_pred.csv'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('model', PointerParam, pointerClass="SchrodingerQSARModel",
                      label='Model: ',
                      help='Model to use for testing.')
        form.addParam('pharmModel', BooleanParam, label='Pharmacophore QSAR model: ',
                       default=False,
                       help='Whether the input model is pharmacophore-based.')

        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Input small molecules:',
                      help='Input small molecules to convert.')

        group = form.addGroup('Filtering params', condition='pharmModel')
        group.addParam('filter', EnumParam, label='Molecule grouping: ', choices=self.filters,
                       default=0,
                       help='How to group molecules: \n'
                            '-distinct             Treat each structure as a distinct molecule (i.e., one conformer only). By default, consecutive structures with identical titles and '
                        'connectivity are treated as conformers of a single molecule. '
                        '-connect              Consider connectivities only (not titles) when perceiving conformers.'
                        '-stereo               Consider stereochemistry when perceiving conformers. Consecutive structures with the same connectivity will be treated as conformers of a'
                        'single molecule if and only if they have the same stereochemistry. Titles are ignored.')

        group.addParam('match', StringParam, label='Minimum matching features: ', default='all',
                       help='Minimum number of hypothesis sites to match. The default is all sites.')

        group = form.addGroup('Scoring params', condition='pharmModel', expertLevel=LEVEL_ADVANCED)
        group.addParam('aw', FloatParam, label='Alignment weight: ', default=1.0,
                       help='Alignment score weight. Must be >= 0.')
        group.addParam('vw', FloatParam, label='Vector weight: ', default=1.0,
                       help='Vector score weight. Must be >= 0.')
        group.addParam('volw', FloatParam, label='Volume weight: ', default=1.0,
                       help='Volume score weight. Must be >= 0.')

        group.addParam('ac', FloatParam, label='Alignment cutoff: ', default=1.2,
                       help='Alignment score cutoff. Must be >= 0.')
        group.addParam('vc', FloatParam, label='Vector cutoff: ', default=-1.0,
                       help='Eliminate hits with vector scores below this value. Must lie on [-1, 1].')
        group.addParam('volc', FloatParam, label='Volume cutoff: ', default=0.0,
                       help='Eliminate hits with volume scores below this value. Must lie on [0, 1].')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.getSmilesStep)
        self._insertFunctionStep(self.createSdfStep)
        if getattr(self.model.get(), 'qsarModel') == 'Field':
            self._insertFunctionStep(self.runPhaseQSARStep)
        else:
            self._insertFunctionStep(self.runPhaseQSARPharmStep)
            self._insertFunctionStep(self.convertOutputFilesStep)
        self._insertFunctionStep(self.createOutputStep)

    def getSmilesStep(self):
        inputSet = self.inputSmallMolecules.get()
        txtFile = self._getExtraPath("sdf_list.txt")

        with open(txtFile, "w") as f:
            for mol in inputSet:
                f.write(f"{os.path.abspath(mol.getFileName())}\n")

        smilesCsv = self._getExtraPath("qsar_dataset.csv")

        script = "rdkit_IO.py"
        args = [
            "-i", os.path.abspath(txtFile),
            "-of", "smiles_csv",
            "-o", os.path.splitext(os.path.basename(smilesCsv))[0],
            "-od", os.path.abspath(self._getExtraPath())
        ]
        pwchemPlugin.runScript(
            self,
            script,
            args,
            env=RDKIT_DIC,
            cwd=self._getPath()
        )

    def createSdfStep(self):
        csvFile = self._getExtraPath("qsar_dataset.csv")
        sdfFile = self.getSDFFile()
        createSdf(self, csvFile, sdfFile)

    def runPhaseQSARStep(self):
        """Run Schrödinger Phase field-based QSAR directly from SDF input."""

        inputSdf = self.getSDFFile()

        outDir = self.getOutDir()
        os.makedirs(outDir, exist_ok=True)

        outputSdf = os.path.abspath(os.path.join(outDir, "qsar_results.sdf"))
        predictionsCsv = outputSdf.replace(".sdf", "_pred.csv")
        outputField = os.path.abspath(os.path.join(outDir, "qsar_field.csv"))
        sumFile = os.path.abspath(os.path.join(outDir, "qsar_summary.txt"))
        modelFile = os.path.abspath(self.model.get().getModelFile())
        args = [
            inputSdf,
            outputSdf,
            "pIC50",
            "-test",
            "-imod", modelFile,
            "-osum", sumFile,
            "-opred", predictionsCsv,
            "-ofield", outputField
        ]

        prog = schrodingerPlugin.getHome("phase_fqsar")
        self.runJob(prog, args, cwd=self._getExtraPath())

        if os.path.exists(sumFile):
            print("\n===== QSAR SUMMARY =====\n")
            with open(sumFile, 'r') as f:
                print(f.read())

    def runPhaseQSARPharmStep(self):
        inputModel = self.model.get()
        hypothesis = os.path.abspath(inputModel.hypoFile.get())
        sdfFileMols = self.getSDFFile()
        jobName = 'pred'

        args = [
            sdfFileMols,
            hypothesis,
            jobName,
            f"-{self.filters[self.filter.get()]}",
            "-aw", self.aw.get(),
            "-vw", self.vw.get(),
            "-volw", self.volw.get(),
            "-ac", self.ac.get(),
            "-vc", self.vc.get(),
            "-volc", self.volc.get(),
        ]

        if self.match.get() != 'all':
            args.append("-match")
            args.append(self.match.get())

        prog = schrodingerPlugin.getHome("phase_screen")
        self.runJob(prog, args, cwd=self._getExtraPath())

        timeout = 3600
        waited = 0
        interval = 2

        while not any(f.endswith(".maegz") for f in os.listdir(self._getExtraPath())):
            if waited >= timeout:
                raise RuntimeError("No .maegz file was created in output directory")

            time.sleep(interval)
            waited += interval

    def convertOutputFilesStep(self):
        extraPath = self._getExtraPath()

        maegzFiles = [f for f in os.listdir(extraPath) if f.endswith(".maegz")]
        maegzFile = os.path.join(extraPath, maegzFiles[0])
        sdfFile = maegzFile.replace(".maegz", ".sdf")
        prog = schrodingerPlugin.getHome("utilities/structcat")
        args = [
            "-imae", os.path.abspath(maegzFile),
            "-osd", os.path.abspath(sdfFile)
        ]
        self.runJob(prog, args, cwd=extraPath)

        script = "sdfToCsv.py"
        csvFile = self._getExtraPath(self.pharmCsv)
        args = [os.path.abspath(sdfFile), os.path.abspath(csvFile)]

        pwchemPlugin.runScript(
            self,
            script,
            args,
            env=RDKIT_DIC,
            cwd=self._getPath(),
            scriptDir=self.scriptsDir
        )

    def createOutputStep(self):
        if not self.pharmModel.get():
            outDir = self.getOutDir()
            predFile = os.path.join(outDir, "qsar_results_pred.csv")
            df = pd.read_csv(predFile)
            df["Title"] = df["Title"].astype(str).str.strip()

            predCols = [c for c in df.columns if c.startswith("Pred(")]
            df["Pred_mean"] = df[predCols].mean(axis=1)
            predMap = dict(zip(df["Title"], df["Pred_mean"]))
        else:
            predFile = self._getExtraPath(self.pharmCsv)
            df = pd.read_csv(predFile)
            df["name"] = df["name"].astype(str).str.strip()

            df["predicted_activity"] = df["predicted_activity"].astype(str).str.strip()
            predMap = dict(zip(df["name"], df["predicted_activity"]))

        outMols = SetOfSmallMolecules().create(outputPath=self._getPath())

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            newMol = SmallMolecule()
            newMol.copy(mol)
            molName = str(mol.getMolName()).strip()
            newMol.predictedActivity = Float()

            if molName in predMap:
                newMol.setAttributeValue('predictedActivity', predMap[molName])
            outMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=outMols)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary=[]
        if self.pharmModel.get():
            outputFile = self._getExtraPath(self.pharmCsv)
        else:
            outputFile = self._getPath("qsar_output/qsar_results_pred.csv")

        if os.path.exists(outputFile):
            with open(outputFile, 'r') as f:
                nLines = sum(1 for _ in f)
            nData = max(0, nLines - 1)

            summary.append(f"Predicted activity for {nData} molecules")
        else:
            summary.append("No prediction file found.")

        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        if (self.aw.get() or self.vw.get() or self.volw.get() or self.ac.get()) < 0:
            validations.append("Alignment, vector and volume weights, as well as alignment cutoff must be bigger than 0.")
        if self.vc.get() < -1 or self.vc.get() > 1:
            validations.append('Vector cutoff must be [-1,1].')
        if self.volc.get() < -1 or self.volc.get() > 1:
            validations.append('Volume cutoff must be [0,1].')
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    def getSDFFile(self):
        return os.path.abspath(self._getExtraPath("qsar_dataset.sdf"))

    def getOutDir(self):
        return os.path.abspath(self._getPath("qsar_output"))