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
Short smoke tests for ProtSchrodingerFepRBFE/ABFE: reuse TestGlideDocking's own fixture
(4erf receptor prepared + 2 real docked ligands from the "smallMolecules" dataset) to wire
the two new FEP+ protocols up exactly the way a real workflow would, without duplicating
that setup.

Unlike the rest of this plugin's test suite, these do NOT wait for a real fep_plus/
fep_absolute_binding job to finish (a real FEP+ job runs for hours even with a license and
a GPU - see claude/decisions/schrodinger/fep_plus.md §8). Instead they wait only long
enough to see this protocol's own steps run for real (receptor/ligand prep, pose-viewer
assembly) and confirm the real Schrodinger binary actually got launched with the expected
arguments, by checking for the "Running command" line Scipion itself logs - the same
signal used to confirm every other real command in this session's manual verification
(see the decisions doc). This catches wiring mistakes (bad param names, a step that raises
before ever calling runJob, a wrong binary name) without requiring hours of runtime or,
for ABFE specifically, the separately-licensed Absolute Binding FEP+ product this session's
own test machine turned out not to have installed (a real, confirmed, environment-specific
gap - see protocol_fep_abfe.py's own module docstring - not a code bug, and not something
this test should be allowed to hide by asserting the job actually succeeded).
"""

import os
import time

from . import main_wf
from ..protocols import ProtSchrodingerFepRBFE, ProtSchrodingerFepABFE

# Imported via the module (not "from .main_wf import TestGlideDocking") so that unittest's
# module-level test discovery doesn't also pick up TestGlideDocking as a standalone test of
# *this* module (it already is one of main_wf's own) - confirmed for real: "scipion3 tests
# pwchemSchrodinger.tests.test_fep" re-ran TestGlideDocking.test() a second time otherwise,
# tripling this module's real runtime for no benefit.


def _waitForCommand(testCase, prot, *substrings, attempts=40, sleepTime=10):
    """Poll `prot`'s own run.stdout log until it contains every one of `substrings`, OR
    the protocol's own "<label> executed in a subproccess" line (Scipion prints this the
    instant runJob spawns the real command, before that command's own stdout has
    necessarily been flushed to the log - confirmed for real: a live run of
    ProtSchrodingerFepABFE showed this line appear while the *later* check for the exact
    binary name text in the log still hadn't matched after a 300s/30-attempt window, a
    false negative from being too strict/impatient, not a real wiring failure - see
    claude/decisions/schrodinger/fep_plus.md §8). Reloads the protocol's status each
    iteration via Project._updateProtocol (mutates `prot` in place - the same private
    helper pyworkflow's own BaseTest.launchProtocol/_waitOutput use internally; there is
    no public "updateProtocol" on BaseTest despite one being referenced - dead-code path -
    inside pyworkflow's own launchProtocol). Returns `prot` (updated in place)."""
    # Real, exact text pyworkflow itself logs the instant runJob spawns the real command
    # (pyworkflow/protocol/launch.py: "Protocol %s executed in a subproccess." %
    # protocol.getRunName()) - matched as a plain substring rather than reconstructing
    # getRunName()'s own formatting, since this file is already scoped to this one
    # protocol's own run.stdout.
    subprocessMarker = 'executed in a subproccess'
    for _ in range(attempts):
        testCase.proj._updateProtocol(prot)
        logFile = prot._getLogsPath('run.stdout')
        if os.path.exists(logFile):
            with open(logFile) as f:
                text = f.read()
            if all(s in text for s in substrings) or subprocessMarker in text:
                return prot
        if prot.isFailed() or prot.isFinished():
            return prot
        time.sleep(sleepTime)
    return prot


class TestSchroFepRBFE(main_wf.TestGlideDocking):
    def test(self):
        glideProt = self._runGlideDocking(self.ligProt, self.prepProt, mode=0)
        self._waitOutput(glideProt, 'outputSmallMolecules', sleepTime=10, timeOut=600)
        self.assertIsNotNone(getattr(glideProt, 'outputSmallMolecules', None))

        protFep = self.newProtocol(ProtSchrodingerFepRBFE, ligandSelection=0)  # 0 = All
        protFep.inputSetOfMols.set(glideProt)
        protFep.inputSetOfMols.setExtended('outputSmallMolecules')
        self.proj.launchProtocol(protFep, wait=False)

        protFep = _waitForCommand(self, protFep, 'fep_plus', '-JOBNAME')
        logFile = protFep._getLogsPath('run.stdout')
        logText = open(logFile).read() if os.path.exists(logFile) else ''
        self.assertTrue('fep_plus' in logText and '-JOBNAME' in logText,
                        'ProtSchrodingerFepRBFE never got to actually launching fep_plus - '
                        f'see {logFile}')


class TestSchroFepABFE(main_wf.TestGlideDocking):
    def test(self):
        glideProt = self._runGlideDocking(self.ligProt, self.prepProt, mode=0)
        self._waitOutput(glideProt, 'outputSmallMolecules', sleepTime=10, timeOut=600)
        self.assertIsNotNone(getattr(glideProt, 'outputSmallMolecules', None))
        ligandName = glideProt.outputSmallMolecules.getFirstItem().__str__()

        protFep = self.newProtocol(ProtSchrodingerFepABFE, inputLigand=ligandName)
        protFep.inputSetOfMols.set(glideProt)
        protFep.inputSetOfMols.setExtended('outputSmallMolecules')
        self.proj.launchProtocol(protFep, wait=False)

        protFep = _waitForCommand(self, protFep, 'fep_absolute_binding', '-JOBNAME')
        logFile = protFep._getLogsPath('run.stdout')
        logText = open(logFile).read() if os.path.exists(logFile) else ''
        self.assertTrue('fep_absolute_binding' in logText and '-JOBNAME' in logText,
                        'ProtSchrodingerFepABFE never got to actually launching '
                        f'fep_absolute_binding - see {logFile}')
