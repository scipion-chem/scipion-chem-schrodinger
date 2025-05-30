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
from ..protocols import ProtSchrodingerPrepWizard, ProtSchrodingerDesmondSysPrep, ProtSchrodingerDesmondMD

md_steps = '''{'simTime': 50.0, 'annealTemps': '[300, 0]', 'bondedT': 0.001, 'nearT': 0.001, 'farT': 0.003, 'velResamp': 1.0, 'glueSolute': True, 'trajInterval': 5.0, 'temperature': 10.0, 'deltaMax': 0.1, 'tempMDCons': 0.1, 'annealing': False, 'pressure': 1.01325, 'presMDCons': 2.0, 'surfTension': 0.0, 'restrainForce': 50.0, 'ensemType': 'Minimization (Brownian)', 'thermostat': 'Noose-Hover', 'barostat': 'Martyna-Tobias-Klein', 'coupleStyle': 'Isotropic', 'restrains': 'Solute_heavy_atom'}
{'simTime': 24.0, 'temperature': 300.0, 'ensemType': 'NPT', 'thermostat': 'Langevin', 'barostat': 'Langevin', 'presMDCons': 2.0}'''

class TestDesmondSysPrep(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

        cls.prepProt = cls._runTargetPreparation()
        cls._waitOutput(cls.prepProt, 'outputStructure', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB, wait=False)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runTargetPreparation(cls):
        kwargs = getPrepTargetWizardArgs()
        protPrepWizard = cls.newProtocol(
            ProtSchrodingerPrepWizard,
            cleanPDB=True, waters=False, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 141}',
            **kwargs)
        protPrepWizard.inputAtomStruct.set(cls.protImportPDB)
        protPrepWizard.inputAtomStruct.setExtended('outputPdb')
        cls.proj.launchProtocol(protPrepWizard, wait=False)
        return protPrepWizard

    @classmethod
    def _runSystemPrep(cls, inputProt, mode=0):
        protDesmondPrep = cls.newProtocol(
            ProtSchrodingerDesmondSysPrep,
            inputFrom=mode, prepareTarget=False, solCharge=4)

        if mode == 0:
            protDesmondPrep.inputStruct.set(inputProt)
            protDesmondPrep.inputStruct.setExtended('outputStructure')
        
        cls.proj.launchProtocol(protDesmondPrep, wait=False)
        return protDesmondPrep

    def test(self):
        sysProt = self._runSystemPrep(self.prepProt)
        self._waitOutput(sysProt, 'outputSystem', sleepTime=5)
        self.assertIsNotNone(getattr(sysProt, 'outputSystem', None))


class TestDesmondMD(TestDesmondSysPrep):
    @classmethod
    def _runMD(cls, inputProt):
        protDesmond= cls.newProtocol(
            ProtSchrodingerDesmondMD,
            workFlowSteps=md_steps)

        protDesmond.inputStruct.set(inputProt)
        protDesmond.inputStruct.setExtended('outputSystem')

        cls.proj.launchProtocol(protDesmond, wait=False)
        return protDesmond

    def test(self):
        sysProt = self._runSystemPrep(self.prepProt)
        self._waitOutput(sysProt, 'outputSystem', sleepTime=5)

        mdProt = self._runMD(sysProt)
        self._waitOutput(mdProt, 'outputSystem', sleepTime=5)
        self.assertIsNotNone(getattr(mdProt, 'outputSystem', None))



def getPrepTargetWizardArgs():
    args={'stage1':True, 'hydrogens': 1,
          'stage2':True, 'propKa':True,
          'stage3':True,
          'stage4':True, 'ms':True}
    return args
