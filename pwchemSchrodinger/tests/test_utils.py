# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from pwchem.protocols import ProtChemImportSmallMolecules

from ..protocols import ProtSchrodingerConvert, ProtSchrodingerPrime, ProtSchrodingerSplitStructure
from ..protocols.protocol_convert import molChoices, targetChoices

from .main_wf import TestSchroProtPrep

chainStr = '{"model": 0, "chain": "C", "residues": 93}'
posStr = '{"index": "18-21", "residues": "QIPA"}'

sideChainList = '''{"model": 0, "chain": "C", "index": "18-23", "residues": "QIPASE"}
{"model": 0, "chain": "C", "index": "54-57", "residues": "LFYL"}
'''

class TestSchroConvert(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB, wait=False)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'), filesPattern='*9.mol2')
        cls.launchProtocol(cls.protImportSmallMols)
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)


    def _runConvert(self, inputProt, outFormat=0, mode=1, wait=False):
        protConv = self.newProtocol(
            ProtSchrodingerConvert,
            inputType=mode)

        if mode == 0:
            protConv.inputSmallMolecules.set(inputProt)
            protConv.inputSmallMolecules.setExtended('outputSmallMolecules')
            protConv.outputFormatSmall.set(outFormat)
            protConv.setObjLabel('Convert SmallMolecules to {}'.format(list(molChoices.keys())[outFormat]))

        elif mode == 1:
            protConv.inputStructure.set(inputProt)
            protConv.inputStructure.setExtended('outputPdb')
            protConv.outputFormatTarget.set(outFormat)
            protConv.setObjLabel('Convert AtomStruct to {}'.format(list(targetChoices.keys())[outFormat]))

        elif mode == 11:
            protConv.inputStructure.set(inputProt)
            protConv.inputType.set(1)
            protConv.inputStructure.setExtended('outputStructure')
            protConv.outputFormatTarget.set(outFormat)
            protConv.setObjLabel('Convert AtomStruct to {}'.format(list(targetChoices.keys())[outFormat]))

        self.proj.launchProtocol(protConv, wait=wait)
        return protConv

    def testSmallMols(self):
        protMols = []
        for i, outFormat in enumerate(molChoices.keys()):
            if 'mol2' not in outFormat.lower():
                protMols.append(self._runConvert(self.protImportSmallMols, outFormat=i, mode=0))
            else:
                mol2Mols = i

        for pm in protMols:
            self._waitOutput(pm, 'outputSmallMolecules', sleepTime=5)
        
        protMols.append(self._runConvert(protMols[-1], outFormat=mol2Mols, mode=0, wait=True))

        for pm in protMols:
            self.assertIsNotNone(getattr(pm, 'outputSmallMolecules', None))


    def testTargets(self):
        protTarget = []
        for i, outFormat in enumerate(targetChoices.keys()):
            protTarget.append(self._runConvert(self.protImportPDB, outFormat=i, mode=1))

        for pt in protTarget:
            self._waitOutput(pt, 'outputStructure', sleepTime=5)

        for pt in protTarget:
            self.assertIsNotNone(getattr(pt, 'outputStructure', None))

class TestPrimeSchro(TestSchroProtPrep):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.protPrepWizard = cls._runTargetPreparation()
        cls._waitOutput(cls.protPrepWizard, 'outputStructure', sleepTime=5)

    @classmethod
    def _runPrime(cls, targetProt, mode=0):
        protPrime = cls.newProtocol(
            ProtSchrodingerPrime, operation=mode,
            inputStructure=targetProt.outputStructure)

        if mode == 0:
            protPrime.residueList.set(sideChainList)
        elif mode == 2:
            protPrime.resChain.set(chainStr)
            protPrime.resPosition.set(posStr)

        cls.proj.launchProtocol(protPrime, wait=False)
        return protPrime

    def test(self):
        prots = []
        for i in range(3):
            prots.append(self._runPrime(self.protPrepWizard, mode=i))

        for p in prots:
            self._waitOutput(p, 'outputStructure', sleepTime=5)

        for p in prots:
            self.assertIsNotNone(getattr(p, 'outputStructure', None))


class TestSplitSchro(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='5ni1')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runSplitStructure(cls, mode=0):
        protSplit = cls.newProtocol(
            ProtSchrodingerSplitStructure,
            splitMode=mode)

        protSplit.inputStructure.set(cls.protImportPDB)
        protSplit.inputStructure.setExtended('outputPdb')

        cls.proj.launchProtocol(protSplit, wait=False)
        return protSplit

    def test(self):
        prots = []
        for i in range(3):
            prots.append(self._runSplitStructure(mode=i))

        outputs = {0: ['outputStructureA', 'outputStructureB', 'outputStructureC', 'outputStructureD'],
                   1: ['outputStructure1'],
                   2: ['outputStructure1', 'outputCofactors1']}

        for i, p in enumerate(prots):
            self._waitOutput(p, outputs[i][0], sleepTime=5)

        for i, p in enumerate(prots):
            for o in outputs[i]:
                self.assertIsNotNone(getattr(p, o, None))

