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
import os, re, subprocess
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from .. import Plugin

from pwchem.constants import CIF_DEF_HEADER, CIF_DEF_COLS
from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import writePDBLine, splitPDBLine, cifFromASFile, getBaseName, natural_sort, writeCIFLine, \
  addCifCols, filterCifCols, writeCifBlocks

from pwchemSchrodinger import Plugin as schrodingerPlugin

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')

class ProtSchrodingerSiteMap(EMProtocol):
    """Calls sitemap to predict possible binding sites"""
    _label = 'binding site prediction (sitemap)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct", label='Atomic Structure:',
                      help='Input protein structure where the binding sites will be predicted')
        form.addParam('maxsites', IntParam, default=30, label='Number of predicted sites:', expertLevel=LEVEL_ADVANCED,
                      help='Maximum number of binding sites to be predicted on the structure')
        form.addParam('jobName', StringParam, label='Job Name:', default='', expertLevel=LEVEL_ADVANCED)


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertStep)
        self._insertFunctionStep(self.sitemapStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertStep(self):
      inAs = self.inputAtomStruct.get()
      if not hasattr(inAs, '_maeFile'):
          inFile = inAs.getFileName()
          cifFile = self._getCifFile()
          cifFromASFile(inFile, cifFile, atomStruct=inAs)

          maeFile = self.getInputMaeFile()
          prog = Plugin.getHome('utilities/prepwizard')
          args = ' -WAIT -noprotassign -noimpref -noepik {} {}'.\
            format(os.path.abspath(cifFile), os.path.abspath(maeFile))
          self.runJob(prog, args, cwd=self._getExtraPath())
      else:
        maeFile = inAs._maeFile.get()
        self.mae2cif(maeFile, self._getCifFile())

    def sitemapStep(self):
        prog=Plugin.getHome('sitemap')

        fnIn = os.path.abspath(self.getInputMaeFile())
        args = '-WAIT -prot %s -j %s -keepvolpts' % (fnIn, self.getJobName())
        args += f" -maxsites {self.maxsites.get()}"

        self.runJob(prog, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        fnBinding = self.getMaestroOutput()
        fnStructure = self.getInputMaeFile()
        fnLog = self.getOutputLogFile()

        setFile = self._getPath('structROIs.sqlite')
        if os.path.exists(setFile):
          os.remove(setFile)
        outPockets = SetOfStructROIs(filename=setFile)

        if os.path.exists(fnBinding):
            proteinFile, pocketFiles = self.createOutputStepCIFFile()

            for oFile in pocketFiles:
              pock = StructROI(oFile, proteinFile, fnLog, pClass='SiteMap')
              pock._maeFile = pwobj.String(fnStructure)
              outPockets.append(pock)

            outPockets.buildPDBhetatmFile()

        self._defineOutputs(outputStructROIs=outPockets)

    def _citations(self):
        return


########################## UTILS FUNCTIONS
    def getInputMaeFile(self):
        if not hasattr(self.inputAtomStruct.get(), '_maeFile'):
            maeFile = self._getExtraPath('inputReceptor.mae')
        else:
            maeFile = self.inputAtomStruct.get()._maeFile.get()
        return maeFile

    def getJobName(self):
      if self.jobName.get() != '':
        return self.jobName.get()
      else:
        return self._getInputName()

    def createOutputStepCIFFile(self):
      jobName = self.getJobName()
      cifFiles = self.buildPocketFiles(jobName, outDir=self._getExtraPath())
      proteinFile = self._getExtraPath(jobName + '.cif')
      return proteinFile, cifFiles

    def getMaestroOutput(self):
        return self._getExtraPath("{}_out.maegz".format(self.getJobName()))

    def getOutputLogFile(self):
        return self._getExtraPath('{}.log'.format(self.getJobName()))

    def buildPocketFiles(self, cifName, outDir):
      '''Build pocket files
      cifOut: name of the output (with or without .cif)'''
      cifFiles = []
      pdbFiles = self.searchOutPDBFiles(cifName, outDir)
      for pocketK, pFile in enumerate(pdbFiles):
        cifCols = '\n'.join(CIF_DEF_COLS)
        outStr = CIF_DEF_HEADER.format(cifCols)
        
        with open(pFile) as f:
          for i, line in enumerate(f):
            if line.startswith('HETATM'):
              sline = line.split()
              idx = (5, 8)
              if len(sline) < 11:
                sline = splitPDBLine(line)
                idx = (6, 9)
              coords = [float(c) for c in sline[idx[0]:idx[1]]]
              if coords:
                replacements = [str(i + 1), f'C{i + 1}', 'STP', 'C', 1, pocketK+1, *coords]
                outStr += writeCIFLine(*replacements)
        
        cifFile = pFile.replace('.pdb', '.cif')
        with open(cifFile, 'w') as fo:
          fo.write(outStr)
        cifFiles.append(cifFile)

      return cifFiles

    def searchOutPDBFiles(self, cifName, outDir):
      pdbFiles = []
      for outFile in os.listdir(outDir):
        if cifName in outFile and '.pdb' in outFile:
          pdbFiles.append(os.path.join(outDir, outFile))
      pdbFiles = natural_sort(self.getCompactPockets(pdbFiles))
      return pdbFiles

    def getCompactPockets(self, files):
      bestFiles = {}
      for f in files:
        match = re.search(r"site_(\d+)_volpts", f)
        if match:
          siteId = match.group(1)
          if siteId not in bestFiles or "compact" in f:
            bestFiles[siteId] = f

      return list(bestFiles.values())

    def getInputPath(self):
        return self.inputAtomStruct.get().getFileName()

    def getInputFileName(self):
        return self.getInputPath().split('/')[-1]

    def _getCifFile(self):
      return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

    def _getInputName(self):
        return getBaseName(self.getInputPath())

    def mae2cif(self, maeFile, cifFile):
      command = '{} {} {} -PDBx'.format(structConvertProg, os.path.abspath(maeFile), os.path.abspath(cifFile))
      subprocess.check_call(command, shell=True, cwd=self._getExtraPath())

      cifDic = MMCIF2Dict(cifFile)
      cifDic = filterCifCols(cifDic, CIF_DEF_COLS)
      cifDic = addCifCols(cifDic, '_atom_site.occupancy', 1)
      cifDic = addCifCols(cifDic, '_atom_site.B_iso_or_equiv', 1)
      cifDic = addCifCols(cifDic, '_atom_site.group_PDB', 'ATOM', 0)
      cifDic = addCifCols(cifDic, '_atom_site.pdbx_PDB_model_num', 1)
      cifDic = addCifCols(cifDic, '_atom_site.auth_seq_id', '_atom_site.label_seq_id', copyValues=True, position=13)

      with open(cifFile, 'w') as f:
        f.write(writeCifBlocks(cifDic) + '#\n')

      return cifFile








