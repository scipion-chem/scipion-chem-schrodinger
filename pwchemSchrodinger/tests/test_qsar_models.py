# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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

from ..protocols import ProtSchrodingerQSARTest, ProtSchrodingerQSAR

class TestQSARModel(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)

    @classmethod
    def _runQSARmodel(cls):
        protQSAR = cls.newProtocol(
            ProtSchrodingerQSAR,
            input=0,
            chemblInput=True,
            ids='CHEMBL240',
            actFilter=7.80)

        cls.proj.launchProtocol(protQSAR, wait=True)
        return protQSAR

    def test(self):
        model = self._runQSARmodel()
        self.assertIsNotNone(getattr(model, 'SchrodingerQSARModel', None))

class TestQSARModelTesting(TestQSARModel):
    @classmethod
    def setUpClass(cls):
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)
        cls._runImportSmallMols()

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            defLibraries=True,
            choicesLibraries=0,
            choicesECBL=2
        )
        cls.launchProtocol(cls.protImportSmallMols, wait=True)

    @classmethod
    def _runQSARmodelTesting(cls, protQSAR):
        protQSARtest = cls.newProtocol(
            ProtSchrodingerQSARTest,
            model=protQSAR.SchrodingerQSARModel,
            inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules)

        cls.proj.launchProtocol(protQSARtest, wait=True)
        return protQSARtest

    def test(self):
        protQSAR = self._runQSARmodel()
        modelTest = self._runQSARmodelTesting(protQSAR)
        self.assertIsNotNone(getattr(modelTest, 'outputSmallMolecules', None))

