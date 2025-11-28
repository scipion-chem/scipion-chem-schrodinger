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

from pyworkflow.protocol.params import PointerParam, StringParam, BooleanParam, FloatParam, \
  IntParam, EnumParam, LabelParam
from pyworkflow.utils.path import createLink
from pyworkflow.object import String

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import AtomicStructHandler

from pwchem.utils import cleanPDB, getBaseName

from .. import Plugin
from ..utils import mae2cif

class ProtSchrodingerPrepWizard(EMProtocol):
    """Calls the preparation wizard"""
    _label = 'target preparation (prepwizard)'
    _program = ""

    def _cleanStructureParams(self, form):
        clean = form.addGroup("Clean atomic structure", condition='not manual')
        clean.addParam('cleanPDB', BooleanParam, default=False, label='Clean PDB: ')
        clean.addParam("waters", BooleanParam,
                       label='Remove waters', condition='cleanPDB',
                       default=True, important=True,
                       help='Remove all waters molecules from a pdb file')

        clean.addParam("HETATM", BooleanParam,
                       label='Remove ligands HETATM',
                       default=True, important=True, condition='cleanPDB',
                       help='Remove all ligands and HETATM contained in the protein')

        clean.addParam("rchains", BooleanParam,
                       label='Remove redundant chains',
                       default=False, important=True, condition='cleanPDB',
                       help='Remove redundant chains in the proteins')

        clean.addParam("chain_name", StringParam,
                       label="Conserved chain",
                       important=True,
                       condition="cleanPDB and rchains==True",
                       help="Select the chain on which you want to carry out the "
                            "molecular docking. You must know the protein and structure "
                            "file that you loaded. \n\nFor example, the protein mdm2 "
                            "(4ERF) has a C1 symmetry, which indicates that its chains "
                            "are at least 95% equal, so you would write A, C or E "
                            "(names of the chains).")

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('manual', BooleanParam, default=False, label='Perform preparation manually: ',
                      help='Perform preparation manually using Maestro GUI. A tutorial can be found at '
                           'https://www.youtube.com/watch?v=YRFROyN88Fw&ab_channel=Schr%C3%B6dingerTV')

        form.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input Atomic Structure: ')

        self._cleanStructureParams(form)
        # Defining condition and label variables
        notManualCondition = 'not manual'
        stage1AndNotManualCondition = 'stage1 and not manual'
        stage2AndNotManualCondition = 'stage2 and not manual'
        stage3AndNotManualCondition = 'stage3 and not manual'
        stage4AndNotManualCondition = 'stage4 and not manual'
        optionsModifiableLabel = 'All options can be modified manually from Maestro GUI'

        form.addSection(label='Stage 1')
        form.addParam('manual1', LabelParam, label=optionsModifiableLabel,
                      condition='manual')
        form.addParam('stage1', BooleanParam, default=False, label='Stage 1 preparation:', condition=notManualCondition)

        form.addParam('fillSideChains', BooleanParam, default=True, condition=stage1AndNotManualCondition,
                      label='Fill side chains:',
                      help='Fill missing side chains Prime')
        form.addParam('fillLoops', BooleanParam, default=True, condition=stage1AndNotManualCondition,
                      label='Fill loops:',
                      help='Fill missing loops with Prime')
        form.addParam('disulfides', BooleanParam, default=True, condition=stage1AndNotManualCondition,
                      label='Create bonds to proximal Sulfurs:',
                      help='Delete hydrogens as needed')
        form.addParam('mse', BooleanParam, default=True, condition=stage1AndNotManualCondition,
                      label='Convert Selenomethionine residues to Methionines:')
        form.addParam('hydrogens', EnumParam, default=0, choices=["Don't add hydrogens", "Delete and re-add hydrogens"],
                      condition=stage1AndNotManualCondition, label='Hydrogen treatment:')
        form.addParam('glycosylation', BooleanParam, default=False, condition=stage1AndNotManualCondition,
                      label='Glycosylation:',
                      help='Create bonds to N-linked and O-linked sugars (delete hydrogens as needed)')
        form.addParam('palmitoylation', BooleanParam, default=False, condition=stage1AndNotManualCondition,
                      label='Palmitoylation:',
                      help='Create palmitoylation bonds even if not included in the CONNECT records ( delete hydrogens as needed)')
        form.addParam('captermini', BooleanParam, default=False, condition=stage1AndNotManualCondition,
                      label='Add cap termini:',
                      help='Add ACE and NME termini')
        form.addParam('keepFarWat', BooleanParam, default=False, condition=stage1AndNotManualCondition,
                      label='Keep far waters:',
                      help="Don't delete waters far from het groups")
        form.addParam('watdist', FloatParam, default=5, condition='stage1 and keepFarWat and not manual',
                      label='Water distance cutoff (A)',
                      help='Distance threshold for far waters')
        form.addParam('treatMetals', BooleanParam, default=True, condition=stage1AndNotManualCondition,
                      label='Treat metals:',
                      help="Adding zero-order bonds, etc")

        form.addSection(label='Stage 2. Protonation (Protassign)')
        form.addParam('manual2', LabelParam, label=optionsModifiableLabel,
                      condition='manual')
        form.addParam('stage2', BooleanParam, default=False, label='Stage 2 preparation:', condition=notManualCondition)
        form.addParam('sampleWaters', BooleanParam, default=True, condition=stage2AndNotManualCondition,
                      label='Sample waters:')
        form.addParam('xtal', BooleanParam, default=False, condition=stage2AndNotManualCondition,
                      label='Use crystal symmetry:')
        form.addParam('propKa', BooleanParam, default=False, condition=stage2AndNotManualCondition,
                      label='Adjust protonation to a specific pH:')
        form.addParam('propKapH', FloatParam, default=7, condition='stage2 and propKa and not manual',
                      label='pH:')
        form.addParam('minadjh', BooleanParam, default=True, condition=stage2AndNotManualCondition,
                      label='Minimize all adjustable hydrogens:',
                      help='Titratable hydrogens, water hydrogens, hydroxyls, thiols, and ASN/GLN carboxamide hydrogens')

        form.addSection(label='Stage 3. Restrained minimization (Impref)')
        form.addParam('manual3', LabelParam, label=optionsModifiableLabel,
                      condition='manual')
        form.addParam('stage3', BooleanParam, default=False, label='Stage 3 preparation:', condition=notManualCondition)
        form.addParam('rmsdD', FloatParam, default=0.3, condition=stage3AndNotManualCondition,
                      label='RMSD cutoff:')
        form.addParam('fix', BooleanParam, default=False, condition=stage3AndNotManualCondition,
                      label='Fix heavy atoms:', help='Minimize hydrogens only')
        form.addParam('force', EnumParam, default=1, choices=["2005", "3"],
                      condition=stage3AndNotManualCondition, label='OPLS force field:')

        form.addSection(label='Stage 4. Ionization and tautomerization (Epik)')
        form.addParam('manual4', LabelParam, label=optionsModifiableLabel,
                      condition='manual')
        form.addParam('stage4', BooleanParam, default=False, label='Stage 4 preparation: ', condition=notManualCondition)
        form.addParam('ms', BooleanParam, default=False, condition=stage4AndNotManualCondition,
                      label='Apply Max. States constrain:')
        form.addParam('msN', IntParam, default=1, condition='stage4 and ms and not manual',
                      label='Max. Number of states:')
        form.addParam('epikPh', FloatParam, default=7, condition=stage4AndNotManualCondition,
                      label='pH:')
        form.addParam('epikPht', FloatParam, default=2, condition=stage4AndNotManualCondition,
                      label='pH Range:')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.manual:
          self._insertFunctionStep('preparationManualStep')
          self._insertFunctionStep('createOutputManual')

        else:
          self._insertFunctionStep('preparationStep')
          self._insertFunctionStep('createOutput')

    def preparationManualStep(self):

        oriFile = self.inputAtomStruct.get().getFileName()
        _, inExt = os.path.splitext(oriFile)

        if inExt in ['.pdb', '.mae', '.maegz']:
            fnIn = os.path.abspath(self._getExtraPath("atomStructIn{}".format(inExt)))
            createLink(oriFile, fnIn)

        else:
            fnIn = os.path.abspath(self._getExtraPath("atomStructIn.pdb"))
            aStruct1 = AtomicStructHandler(oriFile)
            aStruct1.write(fnIn)

        self.runJob(Plugin.getHome('maestro'), " %s" % fnIn, cwd=self._getPath())


    def preparationStep(self):
        inputPDB = self.inputAtomStruct.get().getFileName()
        if self.cleanPDB and isinstance(self.inputAtomStruct.get(), AtomStruct):
            inputPDB = self.remove_WLHC()

        prog=Plugin.getHome('utilities/prepwizard')
        if isinstance(self.inputAtomStruct.get(), AtomStruct):
            fnIn = self._getExtraPath("atomStructIn.pdb")
            aStruct1 = AtomicStructHandler(inputPDB)
            aStruct1.write(fnIn)
            fnIn='extra/atomStructIn.pdb'
        else:
            fnIn = self._getExtraPath("atomStructIn.mae")
            createLink(inputPDB, fnIn)
            fnIn='extra/atomStructIn.mae'

        args='-WAIT'
        if self.stage1.get():
            if self.fillSideChains.get():
                args+=' -fillsidechains'
            if self.fillLoops.get():
                args+=' -fillloops'
            if self.disulfides.get():
                args+=' -disulfides'
            if self.mse.get():
                args+=' -mse'
            if self.hydrogens.get()==0:
                args+=" -nohtreat"
            elif self.hydrogens.get()==1:
                args+=" -rehtreat"
            if self.glycosylation.get():
                args+=" -glycosylation"
            if self.palmitoylation.get():
                args+=" -palmitoylation"
            if self.captermini.get():
                args+=" -captermini"
            if self.keepFarWat.get():
                args+=" -keepfarwat -watdist %f"%self.watdist.get()
            if not self.treatMetals.get():
                args+=" -nometaltreat"

        if self.stage2.get():
            if self.sampleWaters.get():
                args+=" -samplewater"
            if self.xtal.get():
                args+=" -xtal"
            if self.propKa.get():
                args+=" -propka_pH %f"%self.propKapH.get()
            else:
                args+=" -nopropka"
            if self.minadjh.get():
                args+=" -minimize_adj_h"
        else:
            args+=" -noprotassign"

        if self.stage3.get():
            args+=" -rmsd %f"%self.rmsdD.get()
            if self.fix.get():
                args+=" -fix"
            if self.force.get()==0:
                args+=" -f 2005"
            else:
                args+=" -f 3"
        else:
            args+=" -noimpref"

        if self.stage4.get():
            if self.ms.get():
                args+=" -ms %d"%self.msN.get()
            args+=" -epik_pH %f"%self.epikPh.get()
            args+=" -epik_pHt %f"%self.epikPht.get()
        else:
            args+=" -noepik"

        args+=' %s %s.maegz' % (fnIn, self.getJobName())
        self.runJob(prog, args, cwd=self._getPath())

    def createOutputManual(self):
        files = [self._getPath('prepwizard_workdir/%s'%file)
                 for file in os.listdir(self._getPath("prepwizard_workdir")) if (file.lower().endswith('out.mae'))]
        if len(files) > 0:
            files.sort(key=os.path.getmtime)
            filesSorted = sorted(files, key=os.path.getmtime)
            fnMae = filesSorted[-1]

            cifFile = self._getPath(getBaseName(fnMae) + '.cif')
            cifFile = mae2cif(fnMae, cifFile, self._getPath())

            cifAS = AtomStruct(filename=os.path.relpath(cifFile))
            cifAS._maeFile = String(fnMae)

            self._defineOutputs(outputStructure=cifAS)
            self._defineSourceRelation(self.inputAtomStruct, cifAS)

    def createOutput(self):
        fnMae = self._getPath(self.getJobName() + '.maegz')
        if os.path.exists(fnMae):
            cifFile = self._getPath(getBaseName(fnMae) + '.cif')
            cifFile = mae2cif(fnMae, cifFile, self._getPath())

            cifAS = AtomStruct(filename=os.path.relpath(cifFile))
            cifAS._maeFile = String(fnMae)

            self._defineOutputs(outputStructure=cifAS)
            self._defineSourceRelation(self.inputAtomStruct, cifAS)

    def getJobName(self):
      return self.inputAtomStruct.get().getFileName().split('/')[-1].split('.')[0]
    
    def remove_WLHC(self):
        """ Clean the pdb file from waters and ligands
        """
        # Get a PDB format file to the protein structure
        pdbIni = self.inputAtomStruct.get().getFileName()
        filename = os.path.splitext(os.path.basename(pdbIni))[0]
        fnPdb = self._getExtraPath('%s_clean.pdb' % filename)

        if self.rchains.get():
            chain = json.loads(self.chain_name.get())  # From wizard dictionary
            chainId = chain["chain"].upper().strip()
        else:
            chainId = None
        cleanedPDB = cleanPDB(self.inputAtomStruct.get().getFileName(), fnPdb,
                               self.waters.get(), self.HETATM.get(), chainId)
        return cleanedPDB


    def _citations(self):
        return ['Sastry2013']
