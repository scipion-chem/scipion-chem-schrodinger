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
import requests
import math
import time
import os, subprocess

# Scipion em imports
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin
from pwchemSchrodinger import SCHRODINGER_DIC

class ProtSchrodingerQSAR(EMProtocol):
    """Create field-based QSAR model"""
    _label = 'field-based QSAR model'

    kinases = ['CHEMBL203', 'CHEMBL1862', 'CHEMBL2971', 'CHEMBL279', 'CHEMBL240']
    GPCRs = ['CHEMBL251', 'CHEMBL210', 'CHEMBL228']
    enzymes = ['CHEMBL204', 'CHEMBL325', 'CHEMBL3927']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        #form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
        #              label='Input small molecules:',
        #              help='Input small molecules to convert')

        form.addParam('type', EnumParam, label='Target: ', default=0,
                       choices=['Kinases', 'GPCRs', 'Enzymes'],
                       help='Target type to build QSAR model.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('getpIC50Step')
        self._insertFunctionStep('runPhaseQSARStep')

    def getpIC50Step(self):
        if self.type.get() == 0:
            targets = self.kinases
        elif self.type.get() == 1:
            targets = self.GPCRs
        else:
            targets = self.enzymes

        allResults = []
        for targetId in targets:
            print(f"Fetching data for {targetId}")
            url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
            params = {
                "target_chemblId": targetId,
                "standard_type": "IC50",
                "assay_type": "B",
                "limit": 1000
            }
            while url:
                try:
                    r = requests.get(url, params=params, timeout=20)
                    if r.status_code != 200:
                        print(f"Failed to fetch {targetId}")
                        break
                    data = r.json()
                    for act in data.get("activities", []):
                        smiles = act.get("canonical_smiles")
                        value = act.get("standard_value")
                        units = act.get("standard_units")
                        chemblId = act.get("molecule_chembl_id")

                        if not smiles or not value or not chemblId:
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
                        else:
                            continue

                        pIC50 = -math.log10(ic50_m)
                        allResults.append({
                            "name": chemblId,
                            "smiles": smiles,
                            "pIC50": pIC50
                        })
                    url = data.get("next")
                    params = None
                except Exception as e:
                    print(f"Error fetching target {targetId}: {e}")
                    break
        grouped = {}
        for row in allResults:
            key = row["name"]
            if key not in grouped:
                grouped[key] = {"smiles": row["smiles"], "values": []}
            grouped[key]["values"].append(row["pIC50"])

        finalData = []
        for name, data in grouped.items():
            avg_pIC50 = sum(data["values"]) / len(data["values"])
            finalData.append({
                "name": name,
                "smiles": data["smiles"],
                "pIC50": avg_pIC50
            })
        csvFile = self._getPath("qsar_dataset.csv")
        with open(csvFile, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["name", "smiles", "pIC50"])
            writer.writeheader()
            writer.writerows(finalData)

        script = "csvToSDF.py"
        sdfFile = self._getPath("qsar_dataset.sdf")
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

        inputSdf = os.path.abspath(self._getPath("qsar_dataset.sdf"))

        outDir = self._getPath("qsar_output")
        os.makedirs(outDir, exist_ok=True)

        outputSdf = os.path.abspath(os.path.join(outDir, "qsar_results.sdf"))
        predictionsCsv = outputSdf.replace(".sdf", "_pred.csv")

        args = [
            inputSdf,  # input SDF
            outputSdf,  # output SDF
            "pIC50",  # activity property (positional)
            "-build",
            "-style", "ff",
            "-opred", predictionsCsv
        ]

        prog = schrodingerPlugin.getHome("phase_fqsar")
        self.runJob(prog, args, cwd=self._getExtraPath())

        print(f"QSAR SDF saved to: {outputSdf}")
        print(f"Predictions CSV saved to: {predictionsCsv}")


    def createOutputStep(self):
        if self.inputType.get() == 0:
            self._defineOutputs(outputSmallMolecules=self.outputSmallMolecules)
        else:
            self._defineOutputs(outputStructure=self.target)

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