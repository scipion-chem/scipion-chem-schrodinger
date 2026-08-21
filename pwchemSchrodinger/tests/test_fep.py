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
Short smoke tests for ProtSchrodingerFepRBFE/ABFE, using the same real data and reuse
pattern as this codebase's GROMACS pmx RBFE/ABFE tests
(scipion-chem-gromacs/gromacs/tests/tests.py, branch ja_pmxFEP, class
TestGromacsPrepareSystem): a real published RBFE benchmark edge - ligands
18625-1/18626-1 against JNK1 (data/tests/smallMolecules/FEP/jnk1_18625-1_18626-1.pdb),
from pmx's own protLig_benchmark/ligand_tutorial.ipynb, grounded in Gapsys/Perez-Benito
et al. 2020, Chem. Sci. 11:1140 - extracted directly via pwchem's real
TestExtractLigand._runExtractLigand (ProtExtractLigands under the hood; real,
already-co-located poses baked into the PDB file itself, no docking involved).
`TestSchroFepPrepareLigands` below only overrides `_runImportPDB` to point at this file
instead of TestExtractLigand's own default (5ni1) - the same minimal-reuse pattern
scipion-chem-gromacs's TestGromacsPrepareSystem(TestExtractLigand) uses for 1uaz.

Real, licensed Schrodinger 2024-4 was available when this was last run (see
claude/decisions/schrodinger/fep_plus.md §8-§13). Both protocols' launches are correct -
they build receptor/ligand inputs and construct the real fep_plus/fep_absolute_binding
command lines with real, confirmed flags - but neither can produce a real .fmp result on
this specific machine: both hit the same underlying gap, a license/install tier that never
shipped the schrodinger.application.scisol.packages.fep submodules RBFE's fep_mapper stage
and ABFE's hot_atom call both depend on (§13). RBFE fails *softly* (multisim marks the
failed stage's job "completed" and exits 0 after skipping fep_mapper and everything after
it - no real free energy is ever computed) where ABFE fails *loudly* (a Python
ProductNotAvailableError traceback, §8.3). An earlier version of this test asserted a real
RBFE output and (wrongly, see §13) treated the fast, no-op "success" as a working end-to-end
run - it wasn't reading multisim.log, so it never saw fep_mapper's own
"ERROR: No compatible scisol installation found". Both tests here instead check that the
real command was actually launched with the expected arguments (see _launchedRealCommand /
_waitForLaunch), since asserting a real result would fail both tests for a reason outside
anyone's control in this environment, not a wiring problem.
"""

import os
import time

from pyworkflow.tests import DataSet
from pwem.protocols import ProtImportPdb
from pwchem import tests as pwchem_tests

from ..protocols import ProtSchrodingerFepRBFE, ProtSchrodingerFepABFE

# JNK1 chain A (data/tests/smallMolecules/FEP/jnk1_18625-1_18626-1.pdb), residues 1-358 -
# same real, already-verified value the GROMACS pmx RBFE test uses for this exact file.
JNK1_CHAIN_STR = '{"model": 0, "chain": "A", "residues": 358}'


def _launchedRealCommand(prot, *substrings):
    """True once `prot`'s own run.stdout log contains every one of `substrings` -
    confirms the real Schrodinger binary was actually invoked with the expected
    arguments, regardless of whether that job then succeeds or fails (both protocols fail
    here on the same real, external scisol licensing gap - see module docstring)."""
    logFile = prot._getLogsPath('run.stdout')
    if not logFile or not os.path.exists(logFile):
        return False
    with open(logFile) as f:
        text = f.read()
    return all(s in text for s in substrings)


def _waitForLaunch(testCase, prot, *substrings, attempts=30, sleepTime=10):
    """Polls until `prot`'s log shows the real command was launched (see
    _launchedRealCommand) or the protocol finishes/fails on its own, whichever comes
    first. Used instead of _waitOutput for both FEP+ protocols on this machine: neither
    can produce a real .fmp result here (missing scisol.packages.fep - see module
    docstring), so asserting a real output would fail for a reason outside anyone's
    control in this environment, not a wiring problem."""
    for _ in range(attempts):
        testCase.proj._updateProtocol(prot)
        if _launchedRealCommand(prot, *substrings) or prot.isFailed() or prot.isFinished():
            break
        time.sleep(sleepTime)


class TestSchroFepPrepareLigands(pwchem_tests.TestExtractLigand):
    """Reuses pwchem's own TestExtractLigand - its real `_runExtractLigand` classmethod
    and `setUpClass` (which calls `_runImportPDB`) - pointed at the real JNK1 RBFE
    benchmark PDB instead of the base class's own default (5ni1). Deliberately does NOT
    rely on TestExtractLigand's own `test()` method (it references a module-level
    `chainStr` that isn't defined in that file's own scope - a real, pre-existing bug in
    that base class, harmless here since both real test classes below define their own
    `test()` and never call the inherited one).

    Imported as `pwchem.tests` (module reference), not `from pwchem.tests import
    TestExtractLigand` directly - binding that name into this module's own namespace
    would make unittest's test discovery treat TestExtractLigand as a standalone test of
    *this* module too (the same real double-discovery issue already found and fixed once
    in this file's history for TestGlideDocking), and actually running it would hit the
    exact bug described above."""

    @classmethod
    def _runImportPDB(cls):
        ds = DataSet.getDataSet('smallMolecules')
        protImportPDB = cls.newProtocol(
            ProtImportPdb, inputPdbData=1,
            pdbFile=ds.getFile('FEP/jnk1_18625-1_18626-1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # Real bug found by actually running this concurrently with another test:
        # TestExtractLigand._runExtractLigand calls cls.proj.launchProtocol(protExtLig)
        # directly (not through BaseTest.launchProtocol's own wrapper, which defaults to
        # blocking) - not reliably synchronous. It "worked" in earlier sequential runs
        # purely because extraction is normally fast enough to finish before the next
        # line ran; under concurrent CPU contention (two scipion3 test processes at
        # once) it wasn't done yet. GromacsPmxRBFE's own test has the identical call
        # sequence and always follows it with exactly this _waitOutput - not a
        # defensive/optional step, a required one.
        cls.protExtract = cls._runExtractLigand(cls.protImportPDB, JNK1_CHAIN_STR)
        cls._waitOutput(cls.protExtract, 'outputSmallMolecules', sleepTime=5, timeOut=120)

    @classmethod
    def _molNamed(cls, tag):
        """ProtExtractLigands guesses the same molName for every ligand extracted from
        one source file (a real, pre-existing pwchem quirk - see
        claude/decisions/amber/pmx_RBFE.md §15) - the pose file path is what actually
        carries each ligand's real residue tag (e.g. "..._L25_900.cif" vs
        "..._L26_901.cif"), same real distinguishing trick the GROMACS pmx RBFE test uses."""
        return str(next(m for m in cls.protExtract.outputSmallMolecules if tag in m.getPoseFile()))

    # Suppresses the inherited TestExtractLigand.test (the real, pre-existing "chainStr
    # undefined" bug described above) from unittest's own discovery - same convention
    # this codebase's GROMACS pmx test suite already uses (TestGromacsPmxABFE.test2 =
    # None) to hide an unwanted inherited/placeholder test method.
    test = None


class TestSchroFepRBFE(TestSchroFepPrepareLigands):
    def test(self):
        molA = self._molNamed('L25')  # ligand 18625-1
        molB = self._molNamed('L26')  # ligand 18626-1

        protFep = self.newProtocol(ProtSchrodingerFepRBFE, ligandSelection=1)  # 1 = Select subset
        protFep.inputSetOfMols.set(self.protExtract)
        protFep.inputSetOfMols.setExtended('outputSmallMolecules')
        protFep.selectedLigands.set(f'{molA}\n{molB}')
        protFep.inputLigand.set(molA)
        protFep.ligandB.set(molB)
        self.proj.launchProtocol(protFep, wait=False)

        _waitForLaunch(self, protFep, 'fep_plus', '-JOBNAME')
        self.assertTrue(_launchedRealCommand(protFep, 'fep_plus', '-JOBNAME'),
                        'ProtSchrodingerFepRBFE never got to actually launching '
                        f'fep_plus - see {protFep._getLogsPath("run.stdout")}')


class TestSchroFepABFE(TestSchroFepPrepareLigands):
    def test(self):
        molA = self._molNamed('L25')  # ligand 18625-1

        protFep = self.newProtocol(ProtSchrodingerFepABFE, inputLigand=molA)
        protFep.inputSetOfMols.set(self.protExtract)
        protFep.inputSetOfMols.setExtended('outputSmallMolecules')
        self.proj.launchProtocol(protFep, wait=False)

        _waitForLaunch(self, protFep, 'fep_absolute_binding', '-JOBNAME')
        self.assertTrue(_launchedRealCommand(protFep, 'fep_absolute_binding', '-JOBNAME'),
                        'ProtSchrodingerFepABFE never got to actually launching '
                        f'fep_absolute_binding - see {protFep._getLogsPath("run.stdout")}')
