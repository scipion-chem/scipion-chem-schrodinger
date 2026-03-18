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
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL, FloatParam, StringParam, IntParam

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin

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

        form.addParam('type', EnumParam, label='Target type: ', default=0,
                       choices=['Kinases', 'GPCRs', 'Enzymes'],
                       help='Target type to build QSAR model.')
        form.addParam('style', StringParam, label='Fields: ', default='ff',
                      help='Fields to include (can be more than one): \n'
                           '- ff: all force fields\n'
                           '- ff_s: force field steric (Lennar-Jones)\n'
                           '- ff_e: force field electrostatic (q/r)\n'
                           '- qm_e: electrostatic field precomputed by Jaguar\n'
                           '- gauss_s: gaussian steric\n'
                           '- gauss_e: gaussian electrostatic\n'
                           '- gauss_h: gaussian hydrophobic\n'
                           '- gauss_a: gaussian H bond acceptor\n'
                           '- gauss_d: gaussian H bond donor\n'
                           '- gauss: all above gaussian fields\n'
                           '- gauss_r: gaussian aromatic ring\n'
                           '- gauss_ext: gauss + gauss_r')
        form.addParam('forceField', EnumParam, label='Force field: ', default=1,
                      choices=['OPLS_2005', 'OPLS4'],
                      help='Force field from which to draw atom based parameters.')
        form.addParam('train', FloatParam, label='Training partition: ', default=0.8,
                      help='Partition of train set.')
        form.addParam('lno', IntParam, label='Leave-n-out cross-validation: ', default=10,
                      help='Number of training set observations to exclude for cross-validation.\n'
                           'Guidelines:\n'
                           '- small datasets (<20 mols): 1\n'
                           '- medium datasets (20-100 mols): 5-10\n'
                           '- large datasets (>100 mols): 10')

        group = form.addGroup('Grid and FF params', expertLevel=LEVEL_ADVANCED)
        group.addParam('grid', FloatParam, label='Grid spacing (Å):', default=1.0,
                       help='Spacing of field points in angstroms (0.5–4.0)')
        group.addParam('extend', FloatParam, label='Grid extension (Å):', default=3.0,
                       help='Distance to extend grid beyond training set limits')
        group.addParam('buff', FloatParam, label='Field exclusion radius (Å):', default=2.0,
                       help='Ignore force field at grid points within this distance from any atom')
        group.addParam('scut', FloatParam, label='Steric cutoff (kcal/mol):', default=30.0,
                       help='Truncate steric force fields at this value')
        group.addParam('ecut', FloatParam, label='Electrostatic cutoff (kcal/mol):', default=30.0,
                       help='Truncate electrostatic fields at this value')
        group.addParam('sd', FloatParam, label='Min field stddev:', default=0.01,
                       help='Ignore fields if standard deviation over training set is less than this')

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
        csvFile = self._getExtraPath("qsar_dataset.csv")
        with open(csvFile, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["name", "smiles", "pIC50"])
            writer.writeheader()
            writer.writerows(finalData)

        script = "csvToSDF.py"
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

        inputSdf = os.path.abspath(self._getPath("qsar_dataset.sdf"))

        outDir = self._getPath("qsar_output")
        os.makedirs(outDir, exist_ok=True)

        outputSdf = os.path.abspath(os.path.join(outDir, "qsar_results.sdf"))
        predictionsCsv = outputSdf.replace(".sdf", "_pred.csv")
        sumFile = os.path.abspath(os.path.join(outDir, "summary.txt"))
        outputField = os.path.abspath(os.path.join(outDir, "qsar_field.csv"))
        modelFile = os.path.abspath(os.path.join(outDir,"qsar_model.pharm"))
        style = self.style.get()
        ff = self.forceField.choices[self.forceField.get()]

        args = [
            inputSdf,
            outputSdf,
            "pIC50",
            "-build",
            "-style", style,
            "-force_field", ff,
            "-pt", self.train.get(),
            "-LNO", self.lno.get(),
            "-grid", self.grid.get(),
            "-extend", self.extend.get(),
            "-buff", self.buff.get(),
            "-scut", self.scut.get(),
            "-ecut", self.ecut.get(),
            "-sd", self.sd.get(),
            "-omod", modelFile,
            "-opred", predictionsCsv,
            "-osum", sumFile,
            "-ofield", outputField
        ]

        prog = schrodingerPlugin.getHome("phase_fqsar")
        self.runJob(prog, args, cwd=self._getExtraPath())

        print(f"QSAR model saved to: {modelFile}")
        print(f"Training/testing SDF saved to: {outputSdf}")
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
        if self.grid.get() < 0.5 or self.grid.get() > 4.0:
            validations.append('Grid parameter for field must be between 0.5 and 4.0')
        return validations

    def _warnings(self):
        warnings = []
        return warnings