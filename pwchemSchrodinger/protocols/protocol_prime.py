# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import os, json

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, EnumParam, FloatParam, LabelParam, TextParam
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from .. import Plugin
from ..objects import SchrodingerAtomStruct

progPrime = Plugin.getHome('prime')

class ProtSchrodingerPrime(EMProtocol):
    """Schrodinger's prime is a structure prediction program 

    User Documentation(AI GENERATED)
    The ProtSchrodingerPrime class is designed to interact with Schrodinger's Prime program, which is a tool for structure prediction. 
    The class provides a protocol for performing various operations such as side chain prediction, minimization of hydrogens, and loop prediction, which are essential for 
    refining and optimizing protein structures before further analysis or docking simulations.

    The protocol starts by allowing the user to define an input atomic structure, which is the protein structure to be processed. 
    The class includes an operation parameter that allows the user to choose between three main tasks: side chain prediction, minimization of all hydrogens, 
    and loop prediction. For side chain prediction, the user specifies which residues are of interest, and for loop prediction, the user defines the starting and 
    ending residues for the loop. Additionally, the user can adjust the radius of the residue sphere for loop prediction and set a minimum overlap parameter to refine the prediction.
    The class generates a job input file that contains the appropriate settings based on the user's choices, including the type of operation (e.g., side chain prediction or loop prediction) 
    and other related parameters. The job input file is then used to execute the Prime tool, and the results are saved to a .maegz file.
    After the process is complete, the class generates an output file that contains the refined structure, which is then linked back to the input structure for traceability. 
    This output file can be used for further simulations or analyses, such as molecular docking.
    
    The ProtSchrodingerPrime class is an important tool for protein structure refinement, as it automates the process of refining and predicting various aspects of the protein structure, 
    such as side chains, hydrogens, and loops, using Schrodinger's Prime software. This class is particularly useful in structural biology and drug discovery, where accurate protein models are essential 
    for understanding molecular interactions.
    """


    _label = 'target fixing (prime)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                       label='Atomic Structure:', allowsNull=False)
        form.addParam('operation', EnumParam,
                      choices=['Side chain prediction', 'Minimization of all hydrogens', 'Loop prediction'],
                      label='Operation', default=0)

        form.addParam('resChain', StringParam, label='Residues chain: ', condition='operation!=1',
                      help='Specify the protein chain where the residues of interest are')
        form.addParam('resPosition', StringParam, label='Residues positions: ', condition='operation!=1',
                      help='Specify the residues positions for side chain / loop prediction.'
                           'You can use Ctrl to select a range')
        form.addParam('addResidue', LabelParam, label='Add defined residues: ', condition='operation==0',
                      help='Click on the wizard to save the specified residues in the list below')
        form.addParam('residueList', TextParam, width=70, condition='operation==0',
                      label='List of residues:', default='',
                      help='List of residues where the side chain will be predicted')

        form.addParam('resSphere', FloatParam, condition='operation==2', expertLevel=LEVEL_ADVANCED,
                      label='Residue sphere (A):', default=7.5)
        form.addParam('minOverlap', FloatParam, condition='operation==2', expertLevel=LEVEL_ADVANCED,
                      label='Minimum overlap:', default=0.7)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('primeStep')
        self._insertFunctionStep('createOutput')

    def primeStep(self):

        fnInputStructure = self.inputStructure.get().getFileName()
        inExt = os.path.splitext(fnInputStructure)[1]
        fnIn = self._getPath("atomStructIn{}".format(inExt))

        createLink(fnInputStructure, fnIn)
        fnIn = os.path.split(fnIn)[1]

        fhJob = open(self._getPath('job.inp'), 'w')
        fhJob.write("STRUCT_FILE %s\n" % fnIn)
        if self.operation.get() == 0:
            # Side chain prediction
            fhJob.write("PRIME_TYPE  SIDE_PRED\n")
            fhJob.write("SELECT  pick\n")
            for i, residue in enumerate(self.getResidueList()):
                fhJob.write("RESIDUE_%d %s\n" % (i, residue.strip()))

        elif self.operation.get() == 1:
            # Minimization of all hydrogens
            fhJob.write("PRIME_TYPE  REAL_MIN\n")
            fhJob.write("SELECT  asl = (atom.ele H)\n")

        elif self.operation.get() == 2:
            # Loop prediction
            fRes, lRes = self.getFLResidues()
            fhJob.write("PRIME_TYPE  LOOP_BLD\n")
            fhJob.write("LOOP_0_RES_0 %s\n" % fRes)
            fhJob.write("LOOP_0_RES_1 %s\n" % lRes)
            fhJob.write("RES_SPHERE %f\n" % self.resSphere.get())
            fhJob.write("MIN_OVERLAP %f\n" % self.minOverlap.get())

        fhJob.write("USE_CRYSTAL_SYMMETRY no\n")
        fhJob.write("USE_RANDOM_SEED no\n")
        fhJob.write("SEED 0\n")
        fhJob.write("EXT_DIEL 80.00\n")
        fhJob.write("USE_MEMBRANE no\n")
        fhJob.close()
        args = '-WAIT -LOCAL job.inp'
        self.runJob(progPrime, args, cwd=self._getPath())

    def createOutput(self):
        fnMae = self._getPath('job-out.maegz')
        if os.path.exists(fnMae):
            maeFile = SchrodingerAtomStruct()
            maeFile.setFileName(fnMae)

            self._defineOutputs(outputStructure=maeFile)
            self._defineSourceRelation(self.inputStructure, maeFile)

    def _citations(self):
        return ['Jacobson2002', 'Jacobson2004']

    def getFLResidues(self):
        cDic, rDic = json.loads(self.resChain.get()), json.loads(self.resPosition.get())
        return '{}:{}'.format(cDic['chain'], rDic['index'].split('-')[0]), \
               '{}:{}'.format(cDic['chain'], rDic['index'].split('-')[1])

    def getResidueList(self):
        resList = []
        for line in self.residueList.get().split('\n'):
            if line.strip():
                rDic = json.loads(line.strip())
                chain, idxs = rDic['chain'], list(map(int, rDic['index'].split('-')))
                for i in range(idxs[0], idxs[1]+1):
                    resList.append('{}:{}'.format(chain, i))
        return resList

