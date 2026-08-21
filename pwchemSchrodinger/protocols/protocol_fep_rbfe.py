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
Relative binding free energy (RBFE) between two or more ligands docked/aligned in the
same pocket of a shared receptor, via Schrodinger's FEP+ (fep_plus): OPLS3/4 force field,
FEP/REST2 enhanced sampling (GCMC water sampling on by default), alchemical ligand A->B
transformation, Bennett Acceptance Ratio (BAR) free-energy estimation, and cycle-closure
error correction.

This is the FEP+ sibling of gromacs.protocols.protocol_pmx_rbfe.GromacsPmxRBFE - same aim
(relative binding free energy between congeneric ligands sharing a receptor pocket), same
"pick 2+ ligands from one SetOfSmallMolecules, chain 3+ by structural similarity" scope,
different underlying engine. See claude/decisions/schrodinger/fep_plus.md for the full
methodology research and, importantly, for exactly which parts of this implementation are
backed by real, confirmed FEP+ documentation versus best-effort inference.

UPDATE (see claude/decisions/schrodinger/fep_plus.md §8): this protocol was subsequently
checked against a real, licensed Schrodinger 2024-4 install (`$SCHRODINGER/fep_plus
--help`, plus a real `fep_plus -JOBNAME x -prepare <bare-2-ligand-pv>` run) - most of what
was originally "best-effort inference" below turned out correct or was fixed to match
reality. The single biggest open question at the time (whether fep_plus accepts a
from-scratch, map-less pose-viewer file) is now **empirically confirmed true**: it built a
real perturbation map and per-leg `.msj` files itself from a bare 2-ligand PV, no `.fmp`
required.

Real, confirmed facts this implementation relies on (see decisions doc for sources):
- File convention: "jobname_pv.maegz" (receptor + all ligands, pose-viewer format) is the
  input; "jobname.fmp" is the perturbation map/job file; "jobname_out.fmp" (+
  "jobname_out.fmpdb" for trajectories) is the results file FEP+ writes back.
- CLI shape: "$SCHRODINGER/fep_plus -HOST <host> -SUBHOST <subhost> -JOBNAME <name>
  <pv-or-fmp-file>" (confirmed real via `fep_plus --help`), accepting a bare pose-viewer
  file with no pre-built map (empirically confirmed, see decisions doc §8). Confirmed real
  optional flags: "-ff {OPLS4|OPLS5}" (default OPLS4), "-water <model>" (default SPC; real
  values SPC/SPCE/TIP3P/TIP3P_CHARMM/TIP4P/TIP4PEW/TIP4P2005/TIP5P/TIP4PD),
  "-ensemble <muVT|NPT|NVT>" (default muVT = GCMC on for the complex/bound leg),
  "-time <ps>" (production time, default 5000/min 500), "-lambda-windows <int>" (default
  12), "-OPLSDIR <dir>" (custom Force Field Builder parameters, which fep_plus can itself
  generate headlessly via "-ffbuilder -ff-host <host:njobs>" - confirmed real, not wired up
  by this protocol), "-RESTART -checkpoint <file>", "-WAIT" (Job Control Options, confirmed
  real - this protocol's blocking `runJob` calls do actually wait for the job this way).
- "fmp2excel.py" is a real, documented utility that exports a ".fmp"'s predicted
  affinities to CSV, with a "-cycle-closure" option for the cycle-closure-corrected ddG.

What is still NOT confirmed (flagged again in _warnings()):
- fmp2excel.py's exact output CSV column headers (the docs describe the content, not the
  header row) - see fep_utils.parseFmpCsv for the lenient, multi-candidate parsing this
  uses instead of hardcoding one guess.

Pipeline per edge (ligand A -> ligand B):
  0. Prepare the shared receptor once (PrepWizard quick structure-fix, same convention
     ProtSchrodingerGlideDocking/MMGBSA/DesmondSysPrep already use for a non-Maestro input).
  1. Prepare both ligands (LigPrep quick pass or plain structconvert) into Maestro format.
  2. Build "<edge>_pv.maegz" (receptor + ligand A + ligand B), the real documented FEP+
     input file convention.
  3. Run fep_plus.
  4. Export "<edge>_out.fmp" to CSV via fmp2excel.py -cycle-closure and parse it; the raw
     ".fmp"/".fmpdb" are always kept as this edge's output regardless of parse success.

With exactly 2 ligands selected, one edge runs (output 'outputFEP'). With 3+, ligands are
ordered into an A->B->C->... chain by RDKit Morgan/Tanimoto similarity (the same algorithm,
reused verbatim in spirit, as GromacsPmxRBFE._computeEdgeChain) and each edge runs, one
'outputFEP_edge<i>' output per edge - a simplified stand-in for a real perturbation network
built in the FEP+ Panel (which supports custom cores, intermediates and closed cycles for
exactly this reason - see decisions doc), same class of caveat as the GROMACS protocol's
own chain (see claude/decisions/amber/pmx_RBFE.md §15).
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
                               computeLigandEdgeChain, getEdgeName, runFmp2Excel, parseFmpCsv,
                               convertMaeToMol2)

progFepPlus = schrodingerPlugin.getHome('fep_plus')

LIGANDS_ALL, LIGANDS_SUBSET = 0, 1
# Real, confirmed values - read directly off "$SCHRODINGER/fep_plus --help" against a real
# Schrodinger 2024-4 install (see claude/decisions/schrodinger/fep_plus.md §8). Passed to
# "-water" verbatim (exact case, no hyphen) - an earlier draft of this protocol guessed a
# smaller, wrongly-cased/hyphenated list before this was verified.
WATER_MODELS = ['SPC', 'SPCE', 'TIP3P', 'TIP3P_CHARMM', 'TIP4P', 'TIP4PEW', 'TIP4P2005',
               'TIP5P', 'TIP4PD']
ENSEMBLES = ['muVT', 'NPT', 'NVT']
FORCE_FIELDS = ['OPLS4', 'OPLS5']


class ProtSchrodingerFepRBFE(EMProtocol):
    """
    Relative binding free energy via Schrodinger FEP+ (fep_plus), for two or more ligands
    picked from the same SetOfSmallMolecules, docked/aligned in the same pocket of a
    shared receptor.
    """
    _label = 'relative binding free energy (FEP+)'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False,
                      help='Set of docked/aligned molecules sharing a common receptor and a '
                           'common core (FEP+ needs the alchemically-transformed atoms to be '
                           'a small peripheral change, ideally <=10 heavy atoms per FEP+\'s '
                           'own validated range, not a full re-docking). Ligands are picked '
                           'from this same set below.')
        form.addParam('ligandSelection', params.EnumParam, choices=['All', 'Select subset'],
                      default=LIGANDS_SUBSET, display=params.EnumParam.DISPLAY_HLIST,
                      label='Ligands to use: ',
                      help='"All": every molecule in the set above is used. "Select subset": '
                           'pick specific ligands below.')
        form.addParam('selectedLigands', params.TextParam, width=70, allowsNull=False,
                      condition=f'ligandSelection=={LIGANDS_SUBSET}',
                      label='Ligands: ',
                      help='2 or more ligands, picked from the set above via the wizard '
                           '(Ctrl+Click or Shift+Click for multiple), one per line. With '
                           'exactly 2, a single A->B edge runs. With 3 or more, this protocol '
                           'orders them into an A->B->C->... chain by pairwise structural '
                           'similarity (RDKit Morgan/Tanimoto) and runs one edge per '
                           'consecutive pair - a simplified stand-in for a real perturbation '
                           'map built in the FEP+ Panel (which supports closed cycles, custom '
                           'cores and intermediates); see '
                           'claude/decisions/schrodinger/fep_plus.md.')
        form.addHidden('inputLigand', params.StringParam, label='Ligand A: ')
        form.addHidden('ligandB', params.StringParam, label='Ligand B: ')

        group = form.addGroup('Preparation')
        group.addParam('prepareReceptor', params.BooleanParam, default=True, label='Prepare target: ',
                       help='Quick-prepare the receptor with PrepWizard (structure-fixing only, '
                            'same "-noprotassign -noimpref -noepik" shortcut used elsewhere in '
                            'this plugin) if it is not already a Schrodinger-prepared structure. '
                            'Run the full ProtSchrodingerPrepWizard protocol upstream first if '
                            'protonation/refinement matter for this system - FEP+\'s own best '
                            'practices call this step critical to a reliable prediction.')
        group.addParam('prepareLigands', params.BooleanParam, default=True, label='Prepare ligands: ',
                       help='Quick-prepare each ligand with LigPrep if it is not already a '
                            'Schrodinger/LigPrep-produced pose. May subtly shift atom positions '
                            '(ionization/tautomerization) - disable to keep the exact input pose. '
                            'FEP+\'s own best practices recommend running the full '
                            'ProtSchrodingerLigPrep protocol (with Epik) upstream to exhaustively '
                            'enumerate stereoisomers/tautomers/protonation states first.')
        group.addParam('precomputedFmpFile', params.FileParam, allowsNull=True,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Precomputed perturbation map (.fmp, optional): ',
                       help='If you already built and exported a perturbation map in the FEP+ '
                            'Panel (Maestro GUI: Settings -> Write, or Export -> Perturbation '
                            'Map), point this at that ".fmp" file and it is passed to fep_plus '
                            'directly instead of the plain pose-viewer file this protocol would '
                            'otherwise build itself. Leave empty to let this protocol build a '
                            'from-scratch 2-ligand pose-viewer file per edge (see this '
                            'protocol\'s docstring for why that path is less certain to work '
                            'headless on every FEP+ version).')

        form.addSection(label='FEP+ settings')
        group = form.addGroup('Simulation (real, confirmed fep_plus flags - see decisions doc §8)')
        group.addParam('forceField', params.EnumParam, choices=FORCE_FIELDS, default=0,
                       label='Force field: ',
                       help='fep_plus "-ff" option. Default OPLS4 (this install\'s default); '
                            'OPLS5 is also available on Schrodinger 2024-4+.')
        group.addParam('waterModel', params.EnumParam, choices=WATER_MODELS, default=0,
                       label='Water model: ',
                       help='fep_plus "-water" option. Default SPC.')
        group.addParam('ensemble', params.EnumParam, choices=ENSEMBLES, default=0,
                       label='Ensemble: ',
                       help='fep_plus "-ensemble" option. muVT (default) is the grand-canonical '
                            'ensemble used for GCMC water sampling in the complex/bound leg - '
                            'recommended by FEP+\'s own documentation, which calls GCMC water '
                            'sampling important for buried/occluded binding pockets. NPT was the '
                            'ensemble used by older FEP+ versions before GCMC; NVT is also '
                            'supported. Switching away from muVT effectively turns GCMC off.')
        group.addParam('simTime', params.FloatParam, default=5000.0,
                       label='Production simulation time (ps): ',
                       help='fep_plus "-time" option. Default 5000.0 ps (5 ns, FEP+\'s own '
                            'documented default); minimum 500.0. Applies per lambda window/leg. '
                            'Raise this for edges that do not converge, as an alternative to '
                            'FEP+\'s own post-hoc "extend_jobname.sh" mechanism (which extends an '
                            'already-run edge rather than starting longer from scratch).')
        group.addParam('lambdaWindows', params.IntParam, default=12,
                       label='Lambda windows: ',
                       help='fep_plus "-lambda-windows" option. Default 12. FEP+\'s own troubleshooting '
                            'guidance suggests raising this (alongside longer sampling) for edges '
                            'with a high Bennett error, to improve phase-space overlap between '
                            'adjacent replicas.')
        group.addParam('oplsDir', params.FolderParam, allowsNull=True, expertLevel=params.LEVEL_ADVANCED,
                       label='Custom OPLS parameters directory (optional): ',
                       help='fep_plus "-OPLSDIR" option: directory with custom_release.opls '
                            'parameters from a Force Field Builder run, for ligand torsions not '
                            'covered by the default OPLS parameters. FEP+ itself can also run '
                            'Force Field Builder headlessly via its own "-ffbuilder"/"-ff-host" '
                            'flags (confirmed real, present in "fep_plus --help") - not wired up '
                            'by this protocol (out of scope for this version), so run it '
                            'separately (Maestro, or "fep_plus ... -ffbuilder -ff-host '
                            'HOST:MAX_JOBS") first if any ligand needs it, then point this param '
                            'at the resulting directory.')
        group.addParam('fepNotes', params.LabelParam,
                       label='Custom cores/linker-enumeration (fep_plus\'s real "-atom-mapping '
                             '<SMARTS>" flag) and per-leg time overrides ("-time-complex"/'
                             '"-time-solvent") are not exposed here - see '
                             'claude/decisions/schrodinger/fep_plus.md for the full, real flag '
                             'list ("fep_plus --help") if you need them.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        edges = self._computeEdgeChain()
        recStep = self._insertFunctionStep(self.prepareReceptorStep)

        if len(edges) <= 1:
            if edges:
                self.inputLigand.set(edges[0][0])
                self.ligandB.set(edges[0][1])
                self._store(self.inputLigand, self.ligandB)
            self._insertOneEdge(0, edges[0] if edges else (None, None), prerequisites=[recStep],
                               outputName='outputFEP')
            return

        for edgeIdx, (ligA, ligB) in enumerate(edges):
            self._insertOneEdge(edgeIdx, (ligA, ligB), prerequisites=[recStep],
                               outputName=f'outputFEP_edge{edgeIdx}')

    def _insertOneEdge(self, edgeIdx, edge, prerequisites, outputName):
        ligA, ligB = edge
        prepStep = self._insertFunctionStep(self.prepareEdgeLigandsStep, edgeIdx, ligA, ligB,
                                            prerequisites=prerequisites)
        pvStep = self._insertFunctionStep(self.buildPVStep, edgeIdx, ligA, ligB,
                                          prerequisites=[prepStep])
        runStep = self._insertFunctionStep(self.runFepStep, edgeIdx, ligA, ligB,
                                           prerequisites=[pvStep])
        self._insertFunctionStep(self.createOutputStep, edgeIdx, ligA, ligB, outputName,
                                 prerequisites=[runStep])

    # -- Multi-ligand chain construction (same algorithm as GromacsPmxRBFE) --------
    def _computeEdgeChain(self):
        names = self._getSelectedLigandNames()
        if len(names) < 2:
            return []
        return computeLigandEdgeChain(self, names, lambda n: self._molFileForFingerprint(n))

    def _getSelectedLigandNames(self):
        if self.ligandSelection.get() == LIGANDS_ALL:
            return [mol.__str__() for mol in self.inputSetOfMols.get()]
        return [line.strip() for line in self.selectedLigands.get().split('\n') if line.strip()]

    def _getMolByName(self, name):
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == name:
                return mol.clone()
        raise ValueError(f'Ligand "{name}" not found in the input set of molecules')

    def _molFileForFingerprint(self, name):
        """RDKit (via fep_utils.computeLigandEdgeChain) only reads .mol2/.sdf/.pdb - a
        real ligand pose here can be native Schrodinger .mae/.maegz (Glide docking output)
        or .cif (pwchem's own ProtExtractLigands output - confirmed real via the JNK1
        benchmark test: RDKit silently failed to fingerprint every ligand and the 2-ligand
        edge fell back to selection order, harmless here since 2 items only have one
        possible pairing anyway, but would matter for a real 3+-ligand chain). Neither
        format is RDKit-readable directly. Real gap found alongside the receptor .maegz
        bug above: without this conversion, computeLigandEdgeChain's RDKit parse silently
        fails for every ligand and the whole chain quietly degrades to plain selection
        order - never a crash (by design), but never using real structural similarity
        either. structconvert (used by convertMaeToMol2, despite its name) reads .cif
        directly and correctly, same as .mae/.maegz."""
        mol = self._getMolByName(name)
        molFile = mol.getPoseFile()
        if molFile.endswith('.pdbqt'):
            molFile = pdbqt2other(self, molFile, self._getTmpPath(getBaseName(molFile) + '.pdb'))
        elif molFile.endswith(('.mae', '.maegz', '.cif')):
            molFile = convertMaeToMol2(self, molFile, self._getTmpPath())
        return os.path.abspath(molFile)

    # -- 0. Receptor (shared across every edge) --------------------------------------
    def _getReceptorPdb(self):
        """Despite the name (kept for the internal _receptorMaeFile naming convention
        below), this returns whatever format the receptor already is in if Schrodinger's
        own tools can consume it directly - real bug found by actually running this
        protocol: a Glide-docked SetOfSmallMolecules' getProteinFile() is already a
        Maestro ".maegz" (prepwizard/ligprep/structconvert all accept that directly, per
        their own usage text), but pwem's toPdb() has no idea how to read Maestro format
        and silently produced no output file at all, so prepwizard then failed with
        "File does not exist" on the never-written .pdb. Only structures pwem's own
        converters actually understand (.pdb/.pdbqt) still go through a real conversion."""
        proteinFile = self.inputSetOfMols.get().getProteinFile()
        if proteinFile.endswith(('.mae', '.maegz', '.pdb')):
            return os.path.abspath(proteinFile)
        pdbFile = self._getTmpPath(getBaseName(proteinFile) + '.pdb')
        if proteinFile.endswith('.pdbqt'):
            pdbqt2other(self, proteinFile, pdbFile)
        else:
            toPdb(proteinFile, pdbFile)
        return os.path.abspath(pdbFile)

    def prepareReceptorStep(self):
        if self.precomputedFmpFile.get():
            return
        prepareReceptorMae(self, self._getReceptorPdb(), self._getExtraPath(),
                           prepare=self.prepareReceptor.get())

    def _receptorMaeFile(self):
        pdbFile = self._getReceptorPdb()
        targetName = os.path.splitext(os.path.basename(pdbFile))[0]
        return os.path.abspath(self._getExtraPath(targetName + '.maegz'))

    # -- 1. Per-edge ligand preparation ------------------------------------------------
    def _edgeDir(self, edgeIdx):
        edgeDir = self._getExtraPath(f'edge_{edgeIdx}')
        os.makedirs(edgeDir, exist_ok=True)
        return edgeDir

    def prepareEdgeLigandsStep(self, edgeIdx, ligA, ligB):
        if self.precomputedFmpFile.get():
            return
        edgeDir = self._edgeDir(edgeIdx)
        for name in (ligA, ligB):
            mol = self._getMolByName(name)
            prepareLigandMae(self, mol.getPoseFile(), edgeDir, mol.getMolName(),
                             prepare=self.prepareLigands.get())

    def _ligandMaeFile(self, edgeDir, name):
        mol = self._getMolByName(name)
        return os.path.abspath(os.path.join(edgeDir, mol.getMolName() + '.maegz'))

    # -- 2. Pose-viewer assembly ---------------------------------------------------------
    def buildPVStep(self, edgeIdx, ligA, ligB):
        if self.precomputedFmpFile.get():
            return
        edgeDir = self._edgeDir(edgeIdx)
        buildPoseViewerFile(self, self._receptorMaeFile(),
                           [self._ligandMaeFile(edgeDir, ligA), self._ligandMaeFile(edgeDir, ligB)],
                           self._pvFile(edgeIdx, ligA, ligB))

    def _pvFile(self, edgeIdx, ligA, ligB):
        """Real bug found by actually running this protocol: this must be absolute - see
        the identical fix (and full explanation) in ProtSchrodingerFepABFE._pvFile. This
        used to be built inline a second time (slightly differently) inside buildPVStep
        without the abspath wrapper either - both are now this one method."""
        return os.path.abspath(os.path.join(self._edgeDir(edgeIdx), f'{getEdgeName(ligA, ligB)}_pv.maegz'))

    # -- 3. fep_plus job ------------------------------------------------------------------
    def _fepPlusArgs(self):
        """Optional fep_plus flags common to this protocol's edges - every flag name/value
        here is real and confirmed against a real "fep_plus --help" (Schrodinger 2024-4,
        see claude/decisions/schrodinger/fep_plus.md §8), only default-valued ones are
        omitted to keep the command line minimal."""
        args = f' -ff {FORCE_FIELDS[self.forceField.get()]}'
        waterModel = WATER_MODELS[self.waterModel.get()]
        if waterModel != 'SPC':
            args += f' -water {waterModel}'
        ensemble = ENSEMBLES[self.ensemble.get()]
        if ensemble != 'muVT':
            args += f' -ensemble {ensemble}'
        if self.simTime.get() != 5000.0:
            args += f' -time {self.simTime.get()}'
        if self.lambdaWindows.get() != 12:
            args += f' -lambda-windows {self.lambdaWindows.get()}'
        if self.oplsDir.get():
            args += f' -OPLSDIR {os.path.abspath(self.oplsDir.get())}'
        return args

    def runFepStep(self, edgeIdx, ligA, ligB):
        edgeDir = self._edgeDir(edgeIdx)
        jobName = getEdgeName(ligA, ligB)
        inputFile = (os.path.abspath(self.precomputedFmpFile.get()) if self.precomputedFmpFile.get()
                    else os.path.abspath(self._pvFile(edgeIdx, ligA, ligB)))
        # "-HOST localhost -SUBHOST localhost -JOBNAME <name> <pv-or-fmp-file>" and "-WAIT"
        # are all real, confirmed flags (see decisions doc §8); a bare, from-scratch 2-ligand
        # pose-viewer file (no pre-built .fmp) was empirically confirmed to be accepted -
        # fep_plus builds the perturbation map itself.
        args = (f'-WAIT -HOST localhost -SUBHOST localhost -JOBNAME {jobName}'
               f'{self._fepPlusArgs()} {inputFile}')
        self.runJob(progFepPlus, args, cwd=edgeDir)

    def _fepResultFmp(self, edgeIdx, ligA, ligB):
        """Real, documented output naming: "jobname_out.fmp" (falls back to "jobname.fmp"
        if a version writes results in place instead, then None)."""
        edgeDir = self._edgeDir(edgeIdx)
        jobName = getEdgeName(ligA, ligB)
        for candidate in (f'{jobName}_out.fmp', f'{jobName}.fmp'):
            fPath = os.path.join(edgeDir, candidate)
            if os.path.exists(fPath):
                return fPath
        return None

    # -- 4. Output -------------------------------------------------------------------------
    def createOutputStep(self, edgeIdx, ligA, ligB, outputName):
        edgeDir = self._edgeDir(edgeIdx)
        fmpFile = self._fepResultFmp(edgeIdx, ligA, ligB)

        if fmpFile is None:
            # Real bug found by actually running this against a live install: fep_plus's
            # own multisim wrapper exits 0 (self.runJob raises nothing) even when an early
            # stage - fep_mapper, which builds the perturbation map before any real FEP
            # simulation runs - fails and every subsequent stage (including fep_launcher,
            # the actual production MD/FEP sampling) is skipped. On this specific machine
            # that happens with "ERROR: No compatible scisol installation found" (see
            # claude/decisions/schrodinger/fep_plus.md), the same underlying gap as ABFE's
            # ProductNotAvailableError, just surfaced as a soft multisim skip instead of a
            # hard crash. Silently substituting run.stdout as the output "file" here would
            # produce a SchrodingerFEPResult that looks real but carries no free-energy
            # data at all - raise instead of manufacturing a fake-looking success.
            multisimLog = os.path.join(edgeDir, f'{getEdgeName(ligA, ligB)}_multisim.log')
            raise RuntimeError(
                f'No .fmp result file was produced for edge {ligA} -> {ligB}. fep_plus '
                f'exited without error, but the real FEP simulation likely never ran - '
                f'check {multisimLog} for a stage that failed or was skipped '
                f'(e.g. "ERROR: No compatible scisol installation found" in fep_mapper).')

        csvFile = os.path.join(edgeDir, f'{getEdgeName(ligA, ligB)}_results.csv')
        dgByTitle = {}
        if runFmp2Excel(self, fmpFile, csvFile, cwd=edgeDir):
            dgByTitle = parseFmpCsv(csvFile, ligandNames=[ligA, ligB])

        outFEP = SchrodingerFEPResult(filename=fmpFile)
        outFEP.setLigands(ligA, ligB)
        if dgByTitle:
            # Relative dG of ligand B with respect to A is what an A->B edge computes;
            # prefer B's row if present (the "new" compound), else whatever was found.
            ligBMolName = self._getMolByName(ligB).getMolName()
            outFEP.setFreeEnergy(dgByTitle.get(ligBMolName, next(iter(dgByTitle.values()))))
        outFEP.setFreeEnergyFile(fmpFile)

        self._defineOutputs(**{outputName: outFEP})
        self._defineSourceRelation(self.inputSetOfMols, outFEP)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        names = self._getSelectedLigandNames()
        if len(names) < 2:
            errors.append('Select at least 2 ligands to compute a relative binding free energy.')
        allNames = [mol.__str__() for mol in self.inputSetOfMols.get()]
        for name in names:
            if name not in allNames:
                errors.append(f'Ligand "{name}" not found in the input set of molecules.')
        return errors

    def _warnings(self):
        return ['fmp2excel.py\'s exact CSV column names could not be verified (see '
               'claude/decisions/schrodinger/fep_plus.md §8) - most other details of this '
               'protocol\'s fep_plus invocation were confirmed against a real installation. '
               'Check the raw ".fmp"/".fmpdb" kept as this protocol\'s output (open it in '
               'the FEP+ Panel) rather than trusting the parsed free-energy value alone, and '
               'inspect the Bennett error and cycle-closure hysteresis for each edge before '
               'trusting a result - FEP+\'s own guidance flags a Bennett error above 0.3 '
               'kcal/mol as needing more sampling (custom cores, intermediates, or extended '
               'edges via '
               '"extend_jobname.sh").']

    def _summary(self):
        summary = []
        if self.isFinished():
            edges = self._computeEdgeChain() or [(self.inputLigand.get(), self.ligandB.get())]
            for edgeIdx, (ligA, ligB) in enumerate(edges):
                outputName = 'outputFEP' if len(edges) <= 1 else f'outputFEP_edge{edgeIdx}'
                out = getattr(self, outputName, None)
                if out is not None:
                    dG = out.getFreeEnergy()
                    summary.append(f'{ligA} -> {ligB}: '
                                  f'{f"ddG = {dG:.2f} kcal/mol" if dG is not None else "see raw .fmp"}')
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append('Relative binding free energy was calculated with Schrodinger FEP+ '
                           '(fep_plus): OPLS3/4 force field, FEP/REST2 enhanced sampling with '
                           'GCMC water sampling, BAR free-energy estimation and cycle-closure '
                           'error correction (Wang, Chambers & Abel 2019, in Bonomi & Camilloni '
                           'eds., Methods Mol Biol 2019:201-232).')
        return methods
