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

import glob
import os

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, StringParam
from pwem.protocols import EMProtocol
from .. import Plugin
from ..objects import SchrodingerAtomStruct
from pwchem.objects import SmallMolecule

class ProtSchrodingerSplitStructure(EMProtocol):
    """Split a structure into different pieces
    User Documentation(AI GHENERATED)
    The ProtSchrodingerSplitStructure class is designed to split a structure into different components using Schrodinger's split_structure tool. 
    This protocol provides a way to decompose a given protein-ligand complex or other molecular structure into individual pieces, such as chains, ligands, cofactors, ions, and waters. 
    The resulting components can be further analyzed, visualized, or used for additional simulations or docking studies.
    The class includes several parameters that control how the splitting process is performed. The splitMode parameter allows users to select the mode of splitting: by chain, by ligand, 
    or by PDB. The splitMode=2 option splits the structure into receptor, individual ligands, non-metal ions, cofactors, and waters. Additionally, users can merge ligands or waters with 
    the closest chain by enabling the mergeLigands and mergeWaters options. The protocol also provides an option to group waters in the structure and to split cofactors and metals into 
    different structures using the splitAll parameter.

    The class also supports the use of ASL (Atom Specification Language) to define which residues or groups should be considered as ligands, cofactors, or ions. Users can specify custom 
    ASL expressions for ligands, cofactors, and ions through the ligandASL, cofactorASL, positiveASL, and negativeASL parameters. This allows for more flexibility and precision in how the structure is split.
    The execution of the protocol is structured into a single main step, splitStep, which calls the Schrodinger split_structure.py script with the appropriate arguments based on the user-defined parameters. 
    This step splits the input structure into the desired components, and the results are saved as separate files. After the splitting process, the protocol processes the output files and organizes them into
      different categories, such as receptor, chain, ligand, cofactor, and ion components. Each output component is defined as a separate structure and linked back to the input structure for traceability.
    The class also includes functionality for handling temporary files and ensuring that the split components are properly named and stored. The getNumber function is used to extract a unique identifier for 
    each component based on its file name, ensuring that each output structure is correctly labeled. The output is then defined using the _defineOutputs method, and the source relation between the input structure 
    and the output components is established using _defineSourceRelation.

    In terms of validation, the class ensures that the input structure is valid and that the parameters are correctly configured. The splitStep method takes care of invoking the Schrodinger tool with the 
    appropriate flags and arguments based on the user's choices. Once the splitting process is complete, the results are saved in the specified output directory, and the individual components are made available 
    for further analysis or simulations.
    The class also provides a summary method, _summary, which can be used to gather and report any relevant information or warnings generated during the execution of the protocol. However, in this case, the 
    method simply returns an empty list.
    
    Overall, the ProtSchrodingerSplitStructure class is a versatile tool for splitting molecular structures into their components, making it easier to analyze individual parts of a complex structure. 
    This can be especially useful in structural biology and drug discovery, where understanding the interactions between specific components, such as receptors, ligands, and cofactors, is crucial.
    """
    _label = 'split structure'
    _program = ""

    def _defineParams(self, form):
        # Defining param help texts
        aSLHelp = 'For help on ASL (Atom Specification Language), see Chap. 3 of https://shaker.umh.es/computing/Schrodinger_suites/maestro_command_reference.pdf'

        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                      label='Input structure:',)
        form.addParam('splitMode', EnumParam, default=0,
                       choices=["Chain", "Ligand", "PDB"],
                       label='Split mode',
                       help='Split structures by chain, ligand or '
                            'pdbMode=pdb will split the structure into 1) receptor, 2) each individual ligand, '
                            '3) all non-metal ions and cofactors, 4) all waters. By default')
        form.addParam('mergeLigands', BooleanParam, default=False, condition="splitMode==0",
                      label='Merge ligands with the closest chain')
        form.addParam('mergeWaters', BooleanParam, default=False, condition="splitMode==0",
                      label='Merge waters with the closest chain')
        form.addParam('keepProperties', BooleanParam, default=False,
                      label='Keep properties')
        form.addParam('groupWaters', BooleanParam, default=False, condition="splitMode!=2",
                      label='Group waters in structure')
        form.addParam('splitAll', BooleanParam, default=False,
                      label='Split cofactors and metals into different structures')
        form.addParam('ligandASL', StringParam, default="", expertLevel=LEVEL_ADVANCED,
                      label='ASL used to define ligand structures',
                      help=aSLHelp)
        form.addParam('cofactorASL', StringParam, default="", expertLevel=LEVEL_ADVANCED,
                      label='ASL used to define cofactor structures',
                      help=aSLHelp)
        form.addParam('positiveASL', StringParam, default="", expertLevel=LEVEL_ADVANCED, condition="splitMode==2",
                      label='ASL used to define non-metal positive ions',
                      help=aSLHelp)
        form.addParam('negativeASL', StringParam, default="", expertLevel=LEVEL_ADVANCED, condition="splitMode==2",
                      label='ASL used to define negative ions',
                      help=aSLHelp)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('splitStep')

    def splitStep(self):
        args=Plugin.getMMshareDir('python/common/split_structure.py')
        if self.splitMode.get()==0:
            args+=" -m chain"
            if self.mergeLigands.get():
                args+=" -merge_ligands_with_chain"
            if self.mergeWaters.get():
                args+=" -merge_waters_with_chain"
        elif self.splitMode.get()==1:
            args+=" -m ligand"
        elif self.splitMode.get()==2:
            args+=" -m pdb"
            if self.positiveASL.get()!="":
                args+=' -positive_ion_asl "%s"'%self.positiveASL.get()
            if self.negativeASL.get()!="":
                args+=' -negative_ion_asl "%s"'%self.negativeASL.get()
        if self.keepProperties.get():
            args+=" --keep_properties"
        if self.splitMode.get()!=3 and self.groupWaters.get():
            args+=" -groupwaters"
        if self.splitAll.get():
            args+=" -splitall"
        if self.ligandASL.get()!="":
            args+=' -ligand_asl "%s"'%self.ligandASL.get()
        if self.cofactorASL.get()!="":
            args+=' -cofactor_asl "%s"'%self.cofactorASL.get()
        args+=" -many_files %s %s"%(self.inputStructure.get().getFileName(),
                                    self._getExtraPath('output.maegz'))

        self.runJob(Plugin.getHome('run'), args)


        def getNumber(fn, suffix):
            fnBase = os.path.splitext(os.path.split(fn)[1])[0]
            tokens = fnBase.split(suffix)
            return tokens[1]

        for fn in glob.glob(self._getExtraPath("output*")):
            if "_receptor" in fn:
                target = SchrodingerAtomStruct(filename=fn)
                number = getNumber(fn, "_receptor")
                outputDict = {'outputStructure%s' % number: target}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, target)
            elif "_chain" in fn:
                target = SchrodingerAtomStruct(filename=fn)
                number = getNumber(fn, "_chain")
                outputDict = {'outputStructure%s' % number: target}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, target)
            elif "_ligand" in fn:
                ligand = SmallMolecule(smallMolFilename=fn, molName='guess')
                number = getNumber(fn, "_ligand")
                outputDict = {'outputLigand%s' % number: ligand}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, ligand)
            elif "_cof_ion" in fn:
                cof = SchrodingerAtomStruct(filename=fn)
                number = getNumber(fn, "_cof_ion")
                outputDict = {'outputCofactors%s' % number: cof}
                self._defineOutputs(**outputDict)
                self._defineSourceRelation(self.inputStructure, cof)
            else:
                print("Scipion: I don't know how to handle %s"%fn)

    def _summary(self):
        summary=[]
        return summary