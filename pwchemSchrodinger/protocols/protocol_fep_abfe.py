# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Absolute binding free energy (ABFE, "AB-FEP") for one ligand, via Schrodinger's FEP+
absolute-binding-FEP workflow: the ligand is alchemically decoupled from the receptor
(bound leg) and from bulk solvent (unbound/free leg), each with OPLS3/4 + FEP/REST2
enhanced sampling, and dG_bind is obtained from the two legs.

This is the FEP+ sibling of gromacs.protocols.protocol_pmx_abfe.GromacsPmxABFE - same aim
(absolute binding free energy for a single ligand, the alternative to RBFE when there is
no structurally similar partner to alchemically morph into/from), different underlying
engine. See claude/decisions/schrodinger/fep_plus.md for the full methodology research.

UPDATE (see claude/decisions/schrodinger/fep_plus.md §8): the binary-name guess this
docstring originally made ("abfep") was checked against a real, licensed Schrodinger
2024-4 install and was WRONG - the real binary is **`fep_absolute_binding`** (confirmed:
it exists, its `--help` is real and quoted below, and running it produced a real Python
traceback through `schrodinger.application.desmond.starter.generator.fep_absolute_binding`
- the actual, current module name, which differs from the older `abfep` module name found
in pre-2024 public API docs; Schrodinger evidently renamed it between releases). Running a
real `-prepare`-only job on this machine got as far as building the perturbation graph and
failed only on `ProductNotAvailableError: Cannot call
schrodinger.application.scisol.packages.fep.hot_atom.overwrite_hotatoms because the
associated product is not installed` - i.e. this specific machine's Schrodinger install has
core FEP+/Desmond but not the separately-licensed Absolute Binding FEP+ add-on. That is a
real product-licensing gap, not a bug in this protocol - see decisions doc §8 for the full
trace and what it does/doesn't prove.

Confirmed real facts (`$SCHRODINGER/fep_absolute_binding --help`, Schrodinger 2024-4):
- Positional input: "A Maestro pv structure file containing one or more ligand structures
  or a fmp generated for Absolute binding" - i.e. the same bare pose-viewer convention
  used by fep_plus (RBFE), and (like fep_plus) NOT restricted to needing a pre-built map.
- Confirmed real flags: "-ff {OPLS4|OPLS5}" (default OPLS4), "-md-sim-time <ps>" (initial
  MD relaxation before FEP, default 1000/min 500), "-fep-sim-time <ps>" (default 5000/min
  500 - analogous to fep_plus's "-time"), "-equilibration-time <ps>", "-time-complex"/
  "-time-solvent" (per-leg overrides), "-ensemble {muVT|NPT|NVT}" (default muVT),
  "-ligand-restraint"/"-adaptive-ligand-restraint" (dihedral restraints on the ligand -
  confirms the "has_ligand_restraint" parameter this docstring's public-API research had
  already found), "-skip-fep-cutoff-rmsd <real>" (default 4.0 A - skip FEP if the ligand
  drifts out of the pocket during the initial MD), "-extend <ligand-file>", "-OPLSDIR",
  "-ffbuilder"/"-ff-host", plus the same Job Control Options as fep_plus
  (-HOST/-SUBHOST/-WAIT/-RESTART/-checkpoint).
- **No "-water" flag exists for this binary** (unlike fep_plus) - an earlier draft of this
  protocol wrongly assumed it shared fep_plus's "-water" option; fixed by removing it.

Pipeline (single ligand, no edge/pair concept - ABFE evaluates one ligand's binding on its
own, same framing as GromacsPmxABFE):
  0. Prepare the receptor and the ligand (same quick PrepWizard/LigPrep helpers as the
     RBFE sibling protocol).
  1. Build a pose-viewer file (receptor + ligand) - the same real "jobname_pv.maegz"
     convention confirmed for fep_plus, and confirmed (see above) to be the documented
     input shape for fep_absolute_binding too.
  2. Run fep_absolute_binding.
  3. Export/parse the resulting ".fmp" the same way as the RBFE protocol (fmp2excel.py
     -cycle-closure), if it produces one in the same format - unconfirmed for ABFE
     specifically (never reached on this machine - see the product-licensing gap above).
     The raw output file is always kept regardless of parse success.
"""

import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb

from pwchem.utils import pdbqt2other, getBaseName

from .. import Plugin as schrodingerPlugin
from ..objects import SchrodingerFEPResult
from ..utils.fep_utils import (prepareReceptorMae, prepareLigandMae, buildPoseViewerFile,
                               runFmp2Excel, parseFmpCsv)

# Confirmed real binary name (Schrodinger 2024-4) - see module docstring §UPDATE.
progAbfep = schrodingerPlugin.getHome('fep_absolute_binding')

ENSEMBLES = ['muVT', 'NPT', 'NVT']
FORCE_FIELDS = ['OPLS4', 'OPLS5']


class ProtSchrodingerFepABFE(EMProtocol):
    """
    Absolute binding free energy for one ligand, via Schrodinger's FEP+ absolute-binding
    ("AB-FEP") workflow.
    """
    _label = 'absolute binding free energy (FEP+)'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False,
                      help='Set of docked/aligned molecules sharing a common receptor. The '
                           'ligand picked below must be genuinely co-located with the receptor '
                           '(a real docked/co-crystallized pose) - the restraint-atom selection '
                           'FEP+ builds the perturbation graph from (confirmed real via '
                           '"fep_absolute_binding --help": "-ligand-restraint"/'
                           '"-adaptive-ligand-restraint") needs real receptor atoms close enough '
                           'to the ligand to pick from, the same requirement the GROMACS/pmx '
                           'sibling protocol\'s own Boresch-restraint setup has, though the exact '
                           'restraint scheme FEP+ uses internally was not confirmed.')
        form.addParam('inputLigand', params.StringParam, allowsNull=False,
                      label='Ligand: ',
                      help='Ligand to compute the absolute binding free energy for, picked '
                           'from the set above.')

        group = form.addGroup('Preparation')
        group.addParam('prepareReceptor', params.BooleanParam, default=True, label='Prepare target: ',
                       help='Quick-prepare the receptor with PrepWizard (structure-fixing only) '
                            'if it is not already a Schrodinger-prepared structure. Run the full '
                            'ProtSchrodingerPrepWizard protocol upstream first if protonation/'
                            'refinement matter for this system.')
        group.addParam('prepareLigands', params.BooleanParam, default=True, label='Prepare ligand: ',
                       help='Quick-prepare the ligand with LigPrep if it is not already a '
                            'Schrodinger/LigPrep-produced pose.')

        form.addSection(label='FEP+ settings')
        group = form.addGroup('Simulation (real, confirmed fep_absolute_binding flags - see decisions doc §8)')
        group.addParam('forceField', params.EnumParam, choices=FORCE_FIELDS, default=0,
                       label='Force field: ',
                       help='"-ff" option. Default OPLS4.')
        group.addParam('ensemble', params.EnumParam, choices=ENSEMBLES, default=0,
                       label='Ensemble: ',
                       help='"-ensemble" option. muVT (default) is the grand-canonical ensemble '
                            'used for GCMC water sampling in the complex/bound leg.')
        group.addParam('mdSimTime', params.FloatParam, default=1000.0,
                       label='Initial MD relaxation time (ps): ',
                       help='"-md-sim-time" option: production MD before the FEP legs start. '
                            'Default 1000.0 ps; minimum 500.0.')
        group.addParam('fepSimTime', params.FloatParam, default=5000.0,
                       label='FEP production simulation time (ps): ',
                       help='"-fep-sim-time" option, analogous to fep_plus\'s "-time". Default '
                            '5000.0 ps (5 ns); minimum 500.0.')
        group.addParam('ligandRestraint', params.BooleanParam, default=False,
                       label='Apply ligand dihedral restraint: ',
                       help='"-ligand-restraint" option: apply dihedral restraints on the ligand '
                            'molecule during the decoupled/free-ligand leg (helps keep a flexible '
                            'ligand from sampling irrelevant conformations while decoupled).')
        group.addParam('oplsDir', params.FolderParam, allowsNull=True, expertLevel=params.LEVEL_ADVANCED,
                       label='Custom OPLS parameters directory (optional): ',
                       help='"-OPLSDIR" option - custom Force Field Builder parameters for the '
                            'ligand, if it has torsions not covered by the default OPLS '
                            'parameters. fep_absolute_binding can also run Force Field Builder '
                            'itself headlessly via "-ffbuilder -ff-host <host:njobs>" (confirmed '
                            'real, not wired up by this protocol).')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        recStep = self._insertFunctionStep(self.prepareReceptorStep)
        ligStep = self._insertFunctionStep(self.prepareLigandStep)
        pvStep = self._insertFunctionStep(self.buildPVStep, prerequisites=[recStep, ligStep])
        runStep = self._insertFunctionStep(self.runAbfeStep, prerequisites=[pvStep])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[runStep])

    def _getMolByName(self, name):
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == name:
                return mol.clone()
        raise ValueError(f'Ligand "{name}" not found in the input set of molecules')

    def getSpecifiedMol(self):
        return self._getMolByName(self.inputLigand.get())

    # -- 0. Receptor / ligand preparation --------------------------------------------
    def _getReceptorPdb(self):
        proteinFile = self.inputSetOfMols.get().getProteinFile()
        if proteinFile.endswith('.pdb'):
            return os.path.abspath(proteinFile)
        pdbFile = self._getTmpPath(getBaseName(proteinFile) + '.pdb')
        if proteinFile.endswith('.pdbqt'):
            pdbqt2other(self, proteinFile, pdbFile)
        else:
            toPdb(proteinFile, pdbFile)
        return os.path.abspath(pdbFile)

    def prepareReceptorStep(self):
        prepareReceptorMae(self, self._getReceptorPdb(), self._getExtraPath(),
                           prepare=self.prepareReceptor.get())

    def _receptorMaeFile(self):
        pdbFile = self._getReceptorPdb()
        targetName = os.path.splitext(os.path.basename(pdbFile))[0]
        return os.path.abspath(self._getExtraPath(targetName + '.maegz'))

    def prepareLigandStep(self):
        mol = self.getSpecifiedMol()
        prepareLigandMae(self, mol.getPoseFile(), self._getExtraPath(), mol.getMolName(),
                         prepare=self.prepareLigands.get())

    def _ligandMaeFile(self):
        mol = self.getSpecifiedMol()
        return os.path.abspath(self._getExtraPath(mol.getMolName() + '.maegz'))

    # -- 1. Pose-viewer assembly --------------------------------------------------------
    def _pvFile(self):
        molName = self.getSpecifiedMol().getMolName()
        return self._getExtraPath(f'{molName}_pv.maegz')

    def buildPVStep(self):
        buildPoseViewerFile(self, self._receptorMaeFile(), [self._ligandMaeFile()], self._pvFile())

    # -- 2. ABFE job ----------------------------------------------------------------------
    def _jobName(self):
        return self.getSpecifiedMol().getMolName()

    def runAbfeStep(self):
        jobName = self._jobName()
        pvFile = os.path.abspath(self._pvFile())
        args = f'-WAIT -HOST localhost -SUBHOST localhost -JOBNAME {jobName}'
        args += f' -ff {FORCE_FIELDS[self.forceField.get()]}'
        ensemble = ENSEMBLES[self.ensemble.get()]
        if ensemble != 'muVT':
            args += f' -ensemble {ensemble}'
        if self.mdSimTime.get() != 1000.0:
            args += f' -md-sim-time {self.mdSimTime.get()}'
        if self.fepSimTime.get() != 5000.0:
            args += f' -fep-sim-time {self.fepSimTime.get()}'
        if self.ligandRestraint.get():
            args += ' -ligand-restraint'
        if self.oplsDir.get():
            args += f' -OPLSDIR {os.path.abspath(self.oplsDir.get())}'
        args += f' {pvFile}'
        self.runJob(progAbfep, args, cwd=self._getExtraPath())

    def _fepResultFmp(self):
        jobName = self._jobName()
        for candidate in (f'{jobName}_out.fmp', f'{jobName}.fmp'):
            fPath = self._getExtraPath(candidate)
            if os.path.exists(fPath):
                return fPath
        return None

    # -- 3. Output -------------------------------------------------------------------------
    def createOutputStep(self):
        fmpFile = self._fepResultFmp()
        ligandName = self.inputLigand.get()

        dgByTitle = {}
        if fmpFile:
            csvFile = self._getExtraPath(f'{self._jobName()}_results.csv')
            if runFmp2Excel(self, fmpFile, csvFile, cwd=self._getExtraPath()):
                dgByTitle = parseFmpCsv(csvFile, ligandNames=[ligandName])

        outFEP = SchrodingerFEPResult(filename=fmpFile or self._getLogsPath('run.stdout'))
        outFEP.setLigands(ligandName)
        if dgByTitle:
            molName = self.getSpecifiedMol().getMolName()
            outFEP.setFreeEnergy(dgByTitle.get(molName, next(iter(dgByTitle.values()))))
        if fmpFile:
            outFEP.setFreeEnergyFile(fmpFile)

        self._defineOutputs(outputFEP=outFEP)
        self._defineSourceRelation(self.inputSetOfMols, outFEP)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        mols = self.inputSetOfMols.get()
        if mols is not None:
            names = [mol.__str__() for mol in mols]
            if self.inputLigand.get() not in names:
                errors.append(f'Ligand "{self.inputLigand.get()}" not found in the input set '
                              f'of molecules.')
        return errors

    def _warnings(self):
        return ['The real binary name ("fep_absolute_binding") and its flags are now '
               'confirmed (see this protocol\'s module docstring and '
               'claude/decisions/schrodinger/fep_plus.md §8), but a real run on the machine '
               'that verified this could not get past building the perturbation graph: '
               '"ProductNotAvailableError: ...schrodinger.application.scisol.packages.fep.'
               'hot_atom... because the associated product is not installed". If you see the '
               'same error, your Schrodinger license/install has core FEP+/Desmond but not '
               'the separately-licensed Absolute Binding FEP+ add-on - contact your '
               'Schrodinger administrator rather than assuming this protocol is broken. '
               'fmp2excel.py\'s exact CSV columns and this driver\'s exact output-file naming '
               'were never reached/confirmed for the same reason.']

    def _summary(self):
        summary = []
        if self.isFinished():
            out = getattr(self, 'outputFEP', None)
            if out is not None:
                dG = out.get().getFreeEnergy()
                summary.append(f'{self.inputLigand.get()}: '
                              f'{f"dG_bind = {dG:.2f} kcal/mol" if dG is not None else "see raw .fmp"}')
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append('Absolute binding free energy was calculated with Schrodinger '
                           'FEP+\'s absolute-binding ("AB-FEP") workflow: bound and free legs '
                           'with OPLS3/4 + FEP/REST2 enhanced sampling.')
        return methods
