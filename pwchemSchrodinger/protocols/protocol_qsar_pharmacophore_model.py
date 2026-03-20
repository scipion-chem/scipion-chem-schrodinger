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
from pyworkflow.protocol.params import PointerParam, EnumParam, STEPS_PARALLEL, FloatParam, FileParam, StringParam, IntParam, BooleanParam

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules
from pwchemSchrodinger.objects import SchrodingerQSARModel

from pwchem.constants import RDKIT_DIC
from pwchem import Plugin as pwchemPlugin

from .. import Plugin as schrodingerPlugin

class ProtSchrodingerQSARPharmacophore(EMProtocol):
    """Create pharmacophore-based QSAR model"""
    _label = 'pharmacophore-based QSAR model'

    kinases = ['CHEMBL203', 'CHEMBL1862', 'CHEMBL2971', 'CHEMBL279', 'CHEMBL240']
    GPCRs = ['CHEMBL251', 'CHEMBL210', 'CHEMBL228']
    enzymes = ['CHEMBL204', 'CHEMBL325', 'CHEMBL3927']

    stepsExecutionMode = STEPS_PARALLEL


    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('input', EnumParam, label='Input source: ', default=0,
                      choices=['ChEMBL', 'SetOfSmallMolecules', 'CSV activity file'],
                      help='Choose whether to obtain molecules directly from ChEMBL or from a set.')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Input small molecules:', condition='input==1',
                      help='Input small molecules to convert.')
        form.addParam('type', EnumParam, label='Target type: ', default=0, condition='input==0',
                      choices=['Kinases', 'GPCRs', 'Enzymes'],
                      help='Target type to build QSAR model.')
        form.addParam('inputFile', FileParam, label="Activity file: ", condition='input==2',
                      help='CSV file with activity info. Each row should be a molecule with a column containing IC50 or pIC50 activity values in nM.')

        group = form.addGroup('Molecule preparation params')
        group.addParam('epik', EnumParam, label='Epik version: ', choices=['Classic', 'Modern'],
                       default=1,
                       help='Epik version to use for ionization.')
        group.addParam('ph', FloatParam, label='Target pH: ', default=7.4,
                      help='Effective/target pH.')
        group.addParam('phTolerance', FloatParam, label='pH tolerance: ', default=1.0,
                      help='pH tolerance for generated structures.')
        group.addParam('emb', BooleanParam, label='Epik metal binding: ', expertLevel=LEVEL_ADVANCED,
                       default=False,
                       help='Run Epik with the metal_binding option so that states appropriate for interactions with metal ions in protein binding pockets are also generated.')
        group.addParam('stereoisomers', IntParam, label='Number of stereoisomers: ', default=32,
                      help='Generate up to this many stereoisomers per input structure.')
        group.addParam('chirality', BooleanParam, label='Respect chirality: ',
                       default=True,
                       help=' Respect chiralities from input geometry when generating stereoisomers.')
        group.addParam('forceFieldLig', EnumParam, label='Force field: ', default=0,
                      choices=['OPLS_2005', 'S-OPLS'], #0=14, 1=16
                      help=' Force-field to be used for the final geometry optimization.')

        group = form.addGroup('Phase project params')
        group.addParam('active', FloatParam, label='Active threshold: ', default=7.0,
                       help='Good binders with pIC50>=*threshold*.')
        group.addParam('inactive', FloatParam, label='Inactive threshold: ', default=5.5,
                       help='Bad binders with pIC50<=*threshold*.')
        group.addParam('repr', IntParam, label='Number of representatives: ', default=30,
                       help='How many actives to keep as representatives.')

        group = form.addGroup('Pharmacophore discovery params')
        group.addParam('sites', StringParam, label='Pharmacophore size: ', default='4:6',
                       help='Search each reference ligand for common pharmacophores containing between <min> and <max> sites. The legal range is 3:7. The actual searchproceeds from <max> down to <min>, and halts before reaching <min> if common pharmacophores containing more than <min> sites are found. Use -ex to force the full range to be considered. This procedure is followed independently for each reference ligand conformer, so it is still possible to obtain common pharmacophores that contain different numbers of sites even if -ex is not used.')
        group.addParam('miss', IntParam, label='Maximum misses (flexibility): ', default=1,
                        help=('Maximum number of actives that can be missed during pharmacophore generation.\n'
                            'The algorithm starts with strict matching (no misses) and gradually relaxes the constraint\n'
                            'until this value is reached.\n\n'
                            'Higher values = more flexible models (less strict).'))
        group.addParam('keep', IntParam, label='Hypotheses per site: ', default=10,
                       help='Maximum number of hypotheses to retain for each number of sites.')

        group.addParam('redun', FloatParam, label='Redundancy: ', default=0.25, expertLevel=LEVEL_ADVANCED,
                       help='Site point positional difference for elimination of redundant pharmacophores.')
        group.addParam('site', FloatParam, label='Site score weight: ', default=1, expertLevel=LEVEL_ADVANCED,
                       help='Site score weight to use when computing Survival score.')
        group.addParam('vect', FloatParam, label='Vector score weight: ', default=1, expertLevel=LEVEL_ADVANCED,
                       help='Vector score weight to use when computing Survival score.')
        group.addParam('vol', FloatParam, label='Volume score weight: ', default=1, expertLevel=LEVEL_ADVANCED,
                       help='Volume score weight to use when computing Survival score.')
        group.addParam('select', FloatParam, label='Selectivity score weight: ', default=1, expertLevel=LEVEL_ADVANCED,
                       help='Selectivity score weight to use when computing Survival score.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.input.get() == 0:
            self._insertFunctionStep('getpIC50Step')
        elif self.input.get() == 1:
            self._insertFunctionStep('extractpIC50Step')
        else:
            self._insertFunctionStep('checkActivityStep')
        self._insertFunctionStep('runLigPrepStep')
        self._insertFunctionStep('createPhaseProjectStep')
        self._insertFunctionStep('createPharmacophoreDiscStep')
        #self._insertFunctionStep('runPhaseQSARStep')
        #self._insertFunctionStep('createOutputStep')

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

    def runLigPrepStep(self):
        inputSdf = self._getExtraPath("qsar_dataset.sdf")
        outputFile = ("ligprep.maegz")

        prog = schrodingerPlugin.getHome("ligprep")

        args = [
            "-isd", os.path.abspath(inputSdf),
            "-omae", (outputFile)
        ]

        if self.epik.get() == 0: epik = '-epik'
        else: epik = '-epikx'
        args.append(epik)
        args.extend(['-ph', self.ph.get(), '-pht', self.phTolerance.get()])
        if self.emb.get(): args.append('-emb')
        if self.chirality.get(): args.append('-g')
        args.extend(['-s', self.stereoisomers.get()])
        if self.forceFieldLig.get() == 0: ff = 14
        else: ff = 16
        args.extend(['-bff', ff])

        args.extend(['-NJOBS', self.numberOfThreads.get()])
        args.extend(['-HOST', 'localhost'])

        self.runJob(prog, args, cwd=self._getExtraPath())

        timeout = 800
        waited = 0
        interval = 2

        while not os.path.exists(os.path.join(self._getExtraPath(), outputFile)):
            if waited >= timeout:
                raise RuntimeError(f"LigPrep did not create expected file: {outputFile}")
            time.sleep(interval)
            waited += interval

    def createPhaseProjectStep(self):
        ligFile = self._getExtraPath("ligprep.maegz")
        projectFile = self._getExtraPath("phaseProject.phprj")

        prog = schrodingerPlugin.getHome("utilities/phase_project")
        # import ligands
        args = [
            os.path.abspath(projectFile),
            'import',
            '-i', os.path.abspath(ligFile),
            '-new',
            '-act', 'pIC50'
        ]
        self.runJob(prog, args, cwd=self._getExtraPath())

        # define active/inactive
        args = [
            os.path.abspath(projectFile),
            'revise',
            '-active', self.active.get(),
            '-inactive', self.inactive.get(),
            '-commit'
        ]
        self.runJob(prog, args, cwd=self._getExtraPath())

        # reduce actives
        args = [
            os.path.abspath(projectFile),
            'revise',
            '-repr', self.repr.get(),
            '-commit'
        ]
        self.runJob(prog, args, cwd=self._getExtraPath())

        # create pharmacophore sites
        args = [
            os.path.abspath(projectFile),
            'revise',
            '-sites'
        ]
        self.runJob(prog, args, cwd=self._getExtraPath())

        # export to phzip
        args = [os.path.abspath(projectFile), 'archive']
        self.runJob(prog, args, cwd=self._getExtraPath())

    def createPharmacophoreDiscStep(self):
        projectFile = self._getExtraPath("phaseProject.phzip")
        baseName = os.path.splitext(os.path.basename(projectFile))[0]
        outputFile = self._getExtraPath(f"{baseName}_find_common.zip")

        prog = schrodingerPlugin.getHome("phase_find_common")
        args = [
            os.path.abspath(projectFile),
            '-LOCAL',
            '-HOST', 'localhost',
            '-sites', self.sites.get(),
            '-miss', self.miss.get(),
            '-keep', self.keep.get(),
            '-redun', self.redun.get(),
            '-site', self.site.get(),
            '-vect', self.vect.get(),
            '-vol', self.vol.get(),
            '-select', self.select.get()
        ]
        self.runJob(prog, args, cwd=self._getExtraPath())

        timeout = 800
        waited = 0
        interval = 2

        while not os.path.exists(outputFile):
            if waited >= timeout:
                raise RuntimeError(f"phase_find_common did not create expected file: {outputFile}")
            time.sleep(interval)
            waited += interval


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
        model.predictionsFile.set(os.path.join(outDir, "qsar_pred.csv"))
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
        if self.site.get() < 0:
            validations.append('Site score weight should be >=0.')
        if self.vect.get() < 0:
            validations.append('Vector score weight should be >=0.')
        if self.vol.get() < 0:
            validations.append('Volume score weight should be >=0.')
        if self.select.get() < 0:
            validations.append('Selectivity score weight should be >=0.')

        return validations

    def _warnings(self):
        warnings = []
        return warnings