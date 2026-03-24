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
import os, subprocess, shutil

# Scipion em imports
from pwem.objects.data import AtomStruct
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL, FloatParam, StringParam, IntParam, FileParam

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules
from pwchemSchrodinger.objects import SchrodingerQSARModel

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin

class ProtSchrodingerQSAR(EMProtocol):
    """Create field-based QSAR model"""
    _label = 'field-based QSAR model'

    kinases = ['CHEMBL203', 'CHEMBL1862', 'CHEMBL2971', 'CHEMBL279', 'CHEMBL240']
    GPCRs = ['CHEMBL251', 'CHEMBL210', 'CHEMBL228']
    enzymes = ['CHEMBL204', 'CHEMBL325', 'CHEMBL3927']


    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('input', EnumParam, label='Input source: ', default=0,
                      choices=['ChEMBL', 'SetOfSmallMolecules', 'CSV activity file'],
                      help='Choose whether to obtain molecules directly from ChEMBL or from a set.')
        form.addParam('chemblInput', BooleanParam, label='Input IDs: ',
                      default=True,
                      help='INput specific CHEMBL IDs or select target type.')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Input small molecules:', condition='input==1',
                      help='Input small molecules to convert.')
        form.addParam('type', EnumParam, label='Target type: ', default=0, condition='input==0 and chemblInput',
                      choices=['Kinases', 'GPCRs', 'Enzymes'],
                      help='Target type to build QSAR model.')
        form.addParam('ids', StringParam, label='CHEMBL IDs: ', default='',
                      help='Input full CHEMBL IDs separated by commas (eg. CHEMBL3927, CHEMBL325.')
        form.addParam('inputFile', FileParam, label="Activity file: ", condition='input==2',
                      help='CSV file with activity info. Each row should be a molecule with a column containing IC50 or pIC50 activity values in nM.')

        form.addParam('actFilter', FloatParam, label='Activity threshold: ', default=0.0,
                      help='Keep molecules with activity over threshold value.')
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

        group = form.addGroup('Grid and FF params')
        group.addParam('grid', FloatParam, label='Grid spacing (Å):', default=1.0,
                       help='Spacing of field points in angstroms (0.5–4.0)')
        group.addParam('extend', FloatParam, label='Grid extension (Å):', default=3.0,
                       help='Distance to extend grid beyond training set limits')
        group.addParam('buff', FloatParam, label='Field exclusion radius (Å):', default=2.0,
                       help='Ignore force field at grid points within this distance from any atom')
        group.addParam('scut', FloatParam, label='Steric cutoff (kcal/mol):', default=30.0,
                       expertLevel=LEVEL_ADVANCED,
                       help='Truncate steric force fields at this value')
        group.addParam('ecut', FloatParam, label='Electrostatic cutoff (kcal/mol):', default=30.0,
                       expertLevel=LEVEL_ADVANCED,
                       help='Truncate electrostatic fields at this value')
        group.addParam('sd', FloatParam, label='Min field stddev:', default=0.01,
                       expertLevel=LEVEL_ADVANCED,
                       help='Ignore fields if standard deviation over training set is less than this')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.input.get() == 0:
            self._insertFunctionStep('getpIC50Step')
        elif self.input.get() == 1:
            self._insertFunctionStep('extractpIC50Step')
        else:
            self._insertFunctionStep('checkActivityStep')
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

    def checkActivityStep(self):
        csvFile = self.inputFile.get()
        extraPath = self._getExtraPath()
        destFile = os.path.join(extraPath, os.path.basename(csvFile))

        shutil.copy(csvFile, destFile)

        with open(csvFile, "r", newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
        if not fieldnames:
            return
        if "pIC50" in fieldnames:
            print("pIC50 already present. No changes made.")
            return
        if "IC50" not in fieldnames:
            print("No IC50 or pIC50 column found.")
            return
        rows = []
        for row in reader:
            value = row.get("IC50")

            if not value:
                continue
            try:
                value = float(value)
            except:
                continue

            if value <= 0:
                continue
            pIC50 = 9 - math.log10(value)

            row["pIC50"] = pIC50
            del row["IC50"]

            rows.append(row)

        newFieldnames = [f for f in fieldnames if f != "IC50"] + ["pIC50"]

        with open(csvFile, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=newFieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def getpIC50Step(self):
        if self.chemblInput.get():
            targets = [t.strip() for t in self.ids.get().split(",")]
        else:
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
                    for attempt in range(5):
                        try:
                            r = requests.get(url, params=params, timeout=60)
                            if r.status_code == 200:
                                break
                            else:
                                print(f"Failed to fetch {targetId} (status {r.status_code})")
                        except requests.exceptions.Timeout:
                            print(f"Timeout fetching {targetId} (attempt {attempt + 1}/5)")
                            time.sleep(2 ** attempt)
                        except Exception as e:
                            print(f"Error fetching {targetId}: {e}")
                            time.sleep(2 ** attempt)
                    else:
                        print(f"Skipping {targetId} after multiple failed attempts")
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
                        elif units == "pM":
                            ic50_m = value * 1e-12
                        elif units == "fM":
                            ic50_m = value * 1e-15
                        else:
                            continue

                        pIC50 = -math.log10(ic50_m)

                        if pIC50 < self.actFilter.get():
                            continue

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
        modelFile = os.path.abspath(os.path.join(outDir,"qsar_model.pharm"))
        style = self.style.get()
        ffNum = self.forceField.get()
        if ffNum == 0: ff = 'OPLS_2005'
        else: ff = 'OPLS4'

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

        print(f"QSAR model saved to: {modelFile}")
        print(f"Training/testing SDF saved to: {outputSdf}")
        print(f"Predictions CSV saved to: {predictionsCsv}")


    def createOutputStep(self):
        outDir = self._getPath("qsar_output")

        model = SchrodingerQSARModel(modelFile = os.path.join(outDir,"qsar_model.pharm"))
        model.summaryFile.set(os.path.join(outDir, "qsar_summary.txt"))
        model.predictionsFile.set(os.path.join(outDir, "qsar_results_pred.csv"))
        model.sdfFile.set(os.path.join(outDir, "qsar_results.sdf"))
        model.fieldFile.set(os.path.join(outDir, "qsar_field.csv"))

        style = self.style.get()
        ffNum = self.forceField.get()
        if ffNum == 0:
            ff = 'OPLS_2005'
        else:
            ff = 'OPLS4'
        model.style.set(style)
        model.forceField.set(ff)
        model.trainFraction.set(self.train.get())
        model.lno.set(self.lno.get())

        self._defineOutputs(SchrodingerQSARModel=model)


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