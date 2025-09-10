# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# Scipion em imports
from pyworkflow.protocol.params import EnumParam

# Scipion chem imports
from pwchem.viewers import SmallMoleculesViewer

# Plugin imports
from ..protocols import ProtSchrodingerGlideDocking, ProtSchrodingerMMGBSA
from ..viewers.viewers_data import MaestroView

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class ProtGlideDockingViewer(SmallMoleculesViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer glide schrodinger docking'
    _targets = [ProtSchrodingerGlideDocking, ProtSchrodingerMMGBSA]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        super()._defineParams(form)
        form.addSection(label='Maestro view')
        form.addParam('displayMaestroPocket', EnumParam,
                       choices=self.getChoices(vType=SET, pymol=False)[0], default=0,
                       label='Display ligands in Maestro: ',
                       help='Display output ligand set using Maestro')

    def _getVisualizeDict(self):
        visDic = super()._getVisualizeDict()
        visDic.update({'displayMaestroPocket': self._viewPocketMaestroDock})
        return visDic

    def getViewerOptions(self):
        return ['ChimeraX']

    def buildComplexDic(self, mols):
        '''Return a dic of the form: {receptorFile: [ligandFiles]}'''
        recDic = {}
        for mol in mols:
            molFile, recFile = mol.getPoseFile(), mol.getProteinFile()
            if recFile in recDic:
                recDic[recFile].append(molFile)
            else:
                recDic[recFile] = [molFile]
        return recDic

    def buildMaestroScript(self, mols, sFile):
        recDic = self.buildComplexDic(mols)
        with open(sFile, 'w') as f:
            ci, cr = 1, 1
            for recFile, molFiles in recDic.items():
                entrStr = '"\nentryimport "'.join([recFile, *molFiles])
                f.write(f'entryimport "{entrStr}"\n')

                f.write(f'entryselectonly entry {ci}\nentryandgroupcreatewithselectedattoplevel "Complex{cr}" {ci}\n')
                f.write(f'entrysetprop property=s_m_title value="Receptor{cr}" entry {ci}\n')
                for cii, molFile in enumerate(molFiles):
                    f.write(f'entryselectonly entry {ci+cii+1}\n'
                            f'entrygroupmoveselectionbyposition 1 parent_group="Complex{cr}"\n')
                ci += len(molFiles) + 1
                cr += 1
        return sFile


    def _viewPocketMaestroDock(self, e=None):
        ligandLabel = self.getEnumText('displayMaestroPocket')
        mols = self.getGroupMols(self.setLigandsDic, ligandLabel)
        maestroScript = self.buildMaestroScript(mols, self.protocol._getExtraPath(f'{ligandLabel}_mae.txt'))
        return [MaestroView(f'-c {maestroScript}')]


