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
from pwchem.tests import TestExtractLigand
from pwchem.utils import assertHandle

from ..protocols import ProtSchrodingerSiteMap, ProtSchrodingerPrepWizard, \
    ProtSchrodingerLigPrep, ProtSchrodingerGrid, ProtSchrodingerGlideDocking, ProtSchrodingerMMGBSA

class TestSchroProtPrep(BaseTest):
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
            inputPdbData=0,
            pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runTargetPreparation(cls):
        kwargs = getPrepTargetWizardArgs()
        protPrepWizard = cls.newProtocol(
            ProtSchrodingerPrepWizard,
            cleanPDB=True, waters=False, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}',
            **kwargs)
        protPrepWizard.inputAtomStruct.set(cls.protImportPDB)
        protPrepWizard.inputAtomStruct.setExtended('outputPdb')
        cls.proj.launchProtocol(protPrepWizard, wait=False)
        return protPrepWizard

    def test(self):
        prepProt = self._runTargetPreparation()
        self._waitOutput(prepProt, 'outputStructure', sleepTime=5)
        self.assertIsNotNone(getattr(prepProt, 'outputStructure', None))


class TestSitemap(TestSchroProtPrep):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.protPrepWizard = cls._runTargetPreparation()
        cls._waitOutput(cls.protPrepWizard, 'outputStructure', sleepTime=5)

    @classmethod
    def _runSitemap(cls, targetProt):
        protSitemap = cls.newProtocol(
            ProtSchrodingerSiteMap,
            inputAtomStruct=targetProt.outputStructure)

        cls.proj.launchProtocol(protSitemap, wait=False)
        return protSitemap

    def test(self):
        protSitemap = self._runSitemap(self.protPrepWizard)
        self._waitOutput(protSitemap, 'outputStructROIs', sleepTime=5)
        self.assertIsNotNone(getattr(protSitemap, 'outputStructROIs', None))


class TestSchroLigPrep(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)
        cls._runImportSmallMols()

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'), filesPattern='*9.mol2')
        cls.launchProtocol(cls.protImportSmallMols, wait=True)

    @classmethod
    def _runLigandPreparation(cls):
        protPrepLigand = cls.newProtocol(
            ProtSchrodingerLigPrep,
            inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules,
            ionization=1)
        cls.launchProtocol(protPrepLigand, wait=True)
        return protPrepLigand

    def test(self):
        ligProt = self._runLigandPreparation()
        self.assertIsNotNone(getattr(ligProt, 'outputSmallMolecules', None))


class TestGridSchro(TestSitemap):
    @classmethod
    def _runGridDefinition(cls, filterProt):
        protGrid = cls.newProtocol(
            ProtSchrodingerGrid, manual=False,
            innerAction=1, diameterNin=0.8,
            outerAction=1, diameterNout=1.2)
        protGrid.inputStructROIs.set(filterProt)
        protGrid.inputStructROIs.setExtended("outputStructROIs")

        cls.launchProtocol(protGrid)
        return protGrid

    def test(self):
        protSitemap = self._runSitemap(self.protPrepWizard)
        self._waitOutput(protSitemap, 'outputStructROIs', sleepTime=5)
        gridProt = self._runGridDefinition(protSitemap)
        self.assertIsNotNone(getattr(gridProt, 'outputGrids', None))


class TestGlideDocking(TestGridSchro, TestSchroLigPrep):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

        cls.prepProt = cls._runTargetPreparation()
        cls.ligProt = cls._runLigandPreparation()
        cls._waitOutput(cls.prepProt, 'outputStructure', sleepTime=5)
        cls._waitOutput(cls.ligProt, 'outputSmallMolecules', sleepTime=5)

        cls.protSitemap = cls._runSitemap(cls.prepProt)
        cls._waitOutput(cls.protSitemap, 'outputStructROIs', sleepTime=5)
        
    def _runGlideDocking(self, ligProt, inputProt, mode=1):
        protGlide = self.newProtocol(
            ProtSchrodingerGlideDocking,
            posesPerLig=2, fromPockets=mode)

        if mode == 0:
            protGlide.inputAtomStruct.set(inputProt)
            protGlide.inputAtomStruct.setExtended('outputStructure')
            protGlide.radius.set(25)
        elif mode == 1:
            protGlide.inputStructROIs.set(inputProt)
            protGlide.inputStructROIs.setExtended('outputStructROIs')
        else:
            protGlide.inputGridSet.set(inputProt)
            protGlide.inputGridSet.setExtended('outputGrids')

        protGlide.inputLibrary.set(ligProt)
        protGlide.inputLibrary.setExtended('outputSmallMolecules')

        self.proj.launchProtocol(protGlide, wait=False)
        return protGlide

    def test(self):
        print('Docking on whole structure')
        glideProt0 = self._runGlideDocking(self.ligProt, self.prepProt, mode=0)

        print('Docking on StructROIs')
        glideProt1 = self._runGlideDocking(self.ligProt, self.protSitemap, mode=1)

        print('Docking on Schrodinger grids')
        gridProt = self._runGridDefinition(self.protSitemap)
        glideProt2 = self._runGlideDocking(self.ligProt, gridProt, mode=2)

        self._waitOutput(glideProt0, 'outputSmallMolecules', sleepTime=5)
        self._waitOutput(glideProt1, 'outputSmallMolecules', sleepTime=5)
        self._waitOutput(glideProt2, 'outputSmallMolecules', sleepTime=5)

        self.assertIsNotNone(getattr(glideProt0, 'outputSmallMolecules', None))
        self.assertIsNotNone(getattr(glideProt1, 'outputSmallMolecules', None))
        self.assertIsNotNone(getattr(glideProt2, 'outputSmallMolecules', None))


def getPrepTargetWizardArgs():
    args={'stage1':True, 'hydrogens': 1,
          'stage2':True, 'propKa':True,
          'stage3':True,
          'stage4':True, 'ms':True}
    return args

class TestMMGBSA(TestExtractLigand):
    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='4ERF')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runMMGBSA(cls, inputProt):
        protMMGBSA = cls.newProtocol(
            ProtSchrodingerMMGBSA,
            flexOption=1
        )

        protMMGBSA.inputSmallMolecules.set(inputProt)
        protMMGBSA.inputSmallMolecules.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(protMMGBSA)
        return protMMGBSA

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protMMGBSA = self._runMMGBSA(protExtract)
        self._waitOutput(protMMGBSA, 'outputSmallMolecules')
        assertHandle(self.assertIsNotNone, getattr(protMMGBSA, 'outputSmallMolecules', None),
                     cwd=protMMGBSA.getWorkingDir())


