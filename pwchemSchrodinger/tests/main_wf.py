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

import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from pwchem.protocols import ProtChemImportSmallMolecules
from pwchem.tests import TestExtractLigand
from pwchem.utils import assertHandle

from ..protocols import ProtSchrodingerSiteMap, ProtSchrodingerPrepWizard, \
    ProtSchrodingerLigPrep, ProtSchrodingerGrid, ProtSchrodingerGlideDocking, ProtSchrodingerMMGBSA, \
    ProtSchrodingerIFD, ProtSchrodingerConvert

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


# Smallest real IFD: initial Glide pose -> Prime refinement of nearby residues -> Glide redock -> IFD score
IFD_TEST_STEPS = [
    {"stageType": "Glide", "dockingPrecision": "Medium (SP)", "posesPerLig": 1, "selfDock": "False"},
    {"stageType": "Compile residues", "residuesCutoff": 5.0},
    {"stageType": "Residues refinement", "nMinPasses": 1},
    {"stageType": "Glide", "dockingPrecision": "Medium (SP)", "posesPerLig": 1, "selfDock": "True"},
    {"stageType": "Score poses", "scoreName": "r_psp_IFDScore",
     "scoreTerms": "  TERM 1.0,r_i_glide_gscore,0\n  TERM 0.05,r_psp_Prime_Energy,1"},
]


class TestSchroIFD(TestGlideDocking):
    @classmethod
    def _runSitemap(cls, targetProt):
        # One pocket -> one IFD job (IFD runs a job per pocket)
        protSitemap = cls.newProtocol(
            ProtSchrodingerSiteMap, inputAtomStruct=targetProt.outputStructure, maxsites=1)
        cls.proj.launchProtocol(protSitemap, wait=False)
        return protSitemap

    def test(self):
        protIFD = self.newProtocol(ProtSchrodingerIFD, fromPockets=1)
        protIFD.inputStructROIs.set(self.protSitemap)
        protIFD.inputStructROIs.setExtended('outputStructROIs')
        protIFD.inputLibrary.set(self.ligProt)
        protIFD.inputLibrary.setExtended('outputSmallMolecules')
        # Same stage format the "Default workflows" wizard writes
        protIFD.workFlowSteps.set(''.join(f'{protIFD.addDefaultForMissing(dict(s))}\n' for s in IFD_TEST_STEPS))
        self.launchProtocol(protIFD)
        self.assertIsNotNone(getattr(protIFD, 'outputSmallMolecules', None))
        self.assertGreater(protIFD.outputSmallMolecules.getSize(), 0, 'IFD output set is empty')

        # Docked (Maestro) poses and their receptors converted to PDB, e.g. for the PyMOL viewer
        protConv = self.newProtocol(ProtSchrodingerConvert, inputType=0, outputFormatSmall=1, usePose=True)
        protConv.inputSmallMolecules.set(protIFD)
        protConv.inputSmallMolecules.setExtended('outputSmallMolecules')
        self.proj.launchProtocol(protConv, wait=False)
        self._waitOutput(protConv, 'outputSmallMolecules', sleepTime=5)
        mol = protConv.outputSmallMolecules.getFirstItem()
        for fn in (mol.getPoseFile(), mol.getProteinFile(), protConv.outputSmallMolecules.getProteinFile()):
            self.assertTrue(fn.endswith('.pdb') and os.path.exists(fn), f'Not converted to PDB: {fn}')
            self.assertFalse(os.path.isabs(fn), f'Path should be relative to the project: {fn}')


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


