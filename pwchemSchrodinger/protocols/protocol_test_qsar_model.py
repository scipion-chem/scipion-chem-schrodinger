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
from pwchemSchrodinger.objects import SchrodingerQSARModel

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin

class ProtSchrodingerQSARTest(EMProtocol):
    """Test QSAR models"""
    _label = 'test QSAR model'


    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('model', PointerParam, pointerClass="SchrodingerQSARModel",
                      label='Model: ',
                      help='Model to use for testing.')

        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Input small molecules:',
                      help='Input small molecules to convert.')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractpIC50Step')
        self._insertFunctionStep('createSdfStep')
        self._insertFunctionStep('runPhaseQSARStep')
        self._insertFunctionStep('createOutputStep')

    def extractpIC50Step(self):
        inputSet = self.inputSmallMolecules.get()
        txtFile = self._getExtraPath("sdf_list.txt")

        with open(txtFile, "w") as f:
            for mol in inputSet:
                f.write(f"{os.path.abspath(mol.getFileName())}\n")

        smilesCsv = self._getExtraPath("qsar_dataset.csv")

        script = "extractSmiles.py"
        args = [
            os.path.abspath(txtFile),
            os.path.abspath(smilesCsv)
        ]
        pwchemPlugin.runScript(
            self,
            script,
            args,
            env=RDKIT_DIC,
            cwd=self._getPath(),
            scriptDir=os.path.join(os.path.dirname(__file__), "../scripts")
        )

        updatedRows = []
        with open(smilesCsv, "r") as f:
            reader = csv.DictReader(f)

            fieldnames = reader.fieldnames

            for row in reader:
                name = row["name"]
                smiles = row["smiles"]
                url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

                params = {
                    "canonical_smiles": smiles,
                    "standard_type": "IC50",
                    "standard_relation": "=",
                    "limit": 100
                }

                pIC50_values = []

                try:
                    for attempt in range(5):
                        try:
                            r = requests.get(url, params=params, timeout=30)
                            if r.status_code == 200:
                                break
                        except:
                            time.sleep(2 ** attempt)
                    else:
                        continue

                    data = r.json()

                    for act in data.get("activities", []):

                        value = act.get("standard_value")
                        units = act.get("standard_units")

                        if not value or not units:
                            continue

                        try:
                            value = float(value)
                        except:
                            continue

                        if units == "nM":
                            ic50_m = value * 1e-9
                        elif units == "uM":
                            ic50_m = value * 1e-6
                        elif units == "mM":
                            ic50_m = value * 1e-3
                        elif units == "pM":
                            ic50_m = value * 1e-12
                        elif units == "fM":
                            ic50_m = value * 1e-15
                        else:
                            continue

                        pIC50 = -math.log10(ic50_m)
                        pIC50_values.append(pIC50)

                except Exception as e:
                    print(f"ChEMBL error for {name}: {e}")
                    continue

                if len(pIC50_values) == 0:
                    continue

                row["pIC50"] = sum(pIC50_values) / len(pIC50_values)
                updatedRows.append(row)

        if "pIC50" not in fieldnames:
            fieldnames = list(fieldnames) + ["pIC50"]
        with open(smilesCsv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(updatedRows)

    def createSdfStep(self):
        script = "csvToSDF.py"
        csvFile = self._getExtraPath("qsar_dataset.csv")
        sdfFile = self._getExtraPath("qsar_dataset.sdf")
        args = [os.path.abspath(csvFile), os.path.abspath(sdfFile)]

        pwchemPlugin.runScript(
            self,
            script,
            args,
            env=RDKIT_DIC,
            cwd=self._getPath(),
            scriptDir=os.path.join(os.path.dirname(__file__), "../scripts"))

    def runPhaseQSARStep(self):
        """Run Schrödinger Phase field-based QSAR directly from SDF input."""

        inputSdf = os.path.abspath(self._getExtraPath("qsar_dataset.sdf"))

        outDir = self._getPath("qsar_output")
        os.makedirs(outDir, exist_ok=True)

        outputSdf = os.path.abspath(os.path.join(outDir, "qsar_results.sdf"))
        predictionsCsv = outputSdf.replace(".sdf", "_pred.csv")
        outputField = os.path.abspath(os.path.join(outDir, "qsar_field.csv"))
        sumFile = os.path.abspath(os.path.join(outDir, "qsar_summary.txt"))
        print(self.model.get())
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

        print(f"Training/testing SDF saved to: {outputSdf}")
        print(f"Predictions CSV saved to: {predictionsCsv}")


    def createOutputStep(self):
        outDir = self._getPath("qsar_output")

        predFile = os.path.join(outDir,"qsar_results_pred.csv")
        df = pd.read_csv(predFile)
        df["Title"] = df["Title"].astype(str).str.strip()

        activityMap = dict(zip(df["Title"], df["Activity"]))
        print(activityMap)

        outMols = SetOfSmallMolecules().create(outputPath=self._getPath())

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            newMol = SmallMolecule()
            newMol.copy(mol)
            molName = str(mol.getMolName()).strip()
            newMol.predictedActivity = Float()

            print(f'---molName: {molName}')

            if molName in activityMap:
                print('TRUE')
                newMol.setAttributeValue('predictedActivity', activityMap[molName])
            outMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=outMols)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary=[]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings