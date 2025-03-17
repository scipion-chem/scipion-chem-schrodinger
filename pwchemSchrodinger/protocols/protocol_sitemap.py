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
import os, shutil
from subprocess import CalledProcessError

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam, StringParam
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from .. import Plugin

from pwchem import Plugin as pwchemPlugin
from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import writePDBLine, splitPDBLine
from pwchem.constants import OPENBABEL_DIC

class ProtSchrodingerSiteMap(EMProtocol):
    """Calls sitemap to predict possible binding sites"""
    _label = 'binding site prediction (sitemap)'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct", label='Atomic Structure:',
                      help='Input protein structure where the binding sites will be predicted')
        form.addParam('maxsites', IntParam, default=5, label='Number of predicted sites:',
                      help='Maximum numbe rof binding sites to be predicted on the structure')
        form.addParam('jobName', StringParam, label='Job Name:', default='', expertLevel=LEVEL_ADVANCED)


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('sitemapStep')
        self._insertFunctionStep('createOutput')

    def convertStep(self):
      if not hasattr(self.inputStructure.get(), '_maeFile'):
          inFile = self.inputStructure.get().getFileName()
          if inFile.endswith('.pdbqt'):
              outName, outDir = os.path.splitext(os.path.basename(inFile))[0], os.path.abspath(self._getTmpPath())
              args = ' -i "{}" -of pdb --outputDir "{}" --outputName {}'.format(os.path.abspath(inFile),
                                                                             os.path.abspath(outDir), outName)
              pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)
              pdbFile = os.path.abspath(os.path.join(outDir, '{}.pdb'.format(outName)))
          else:
              pdbFile = inFile

          maeFile = self.getInputMaeFile()
          prog = Plugin.getHome('utilities/prepwizard')
          print('Program: ', prog)
          args = ' -WAIT -noprotassign -noimpref -noepik {} {}'.\
            format(os.path.abspath(pdbFile), os.path.abspath(maeFile))
          self.runJob(prog, args, cwd=self._getExtraPath())

    def sitemapStep(self):
        prog=Plugin.getHome('sitemap')

        fnIn = os.path.abspath(self.getInputMaeFile())
        args='-WAIT -prot %s -j %s -keepvolpts' % (fnIn, self.getJobName())
        args+=" -maxsites %d"%self.maxsites.get()

        self.runJob(prog, args, cwd=self._getExtraPath())

    def createOutput(self):
        fnBinding = self.getMaestroOutput()
        fnStructure = self.getInputMaeFile()
        fnLog = self.getOutputLogFile()
        if os.path.exists(fnBinding):
            proteinFile, pocketFiles = self.createOutputPDBFile()
            outPockets = SetOfStructROIs(filename=self._getPath('structROIs.sqlite'))
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
        if not hasattr(self.inputStructure.get(), '_maeFile'):
            maeFile = self._getExtraPath('inputReceptor.mae')
        else:
            maeFile = self.inputStructure.get()._maeFile.get()
        return maeFile

    def getJobName(self):
      if self.jobName.get() != '':
        return self.jobName.get()
      else:
        return self.getInputFileName().split('/')[-1].split('.')[0]

    def createOutputPDBFile(self):
      maeFile = os.path.abspath(self.getMaestroOutput())
      pdbOutFile = self.getJobName()+'_out.pdb'
      pdbFiles = self.maestro2pdb(maeFile, pdbOutFile, outDir=self._getPath())

      # Creates a pdb with the HETATM corresponding to pocket points
      pdbFiles, proteinFile = self.mergePDBFiles(pdbFiles, pdbOutFile)
      return proteinFile, pdbFiles

    def getMaestroOutput(self):
        return self._getExtraPath("{}_out.maegz".format(self.getJobName()))

    def getInputFileName(self):
        return self.inputStructure.get().getFileName()

    def getOutputLogFile(self):
        return self._getExtraPath('{}.log'.format(self.getJobName()))

    def maestro2pdb(self, maeIn, pdbOut, outDir):
      '''Convert a maestro file (.mae) to a pdb file(s)
      maeIn: input maestro file (if contains several models, there will be several outputs
      pdbOut: name of the output (with or without .pdb)'''
      pdbOut = self.getPDBName(pdbOut)[1]

      prog = Plugin.getHome('utilities/structconvert')
      args = '{} {}'.format(maeIn, pdbOut)
      try:
          self.runJob(prog, args, cwd=self._getPath())
      except CalledProcessError as exception:
          # ask to Schrodinger why it returns a code 2 if it worked properly
          if exception.returncode != 2:
              raise exception

      pdbFiles = self.searchOutPDBFiles(pdbOut, outDir)
      if len(pdbFiles) > 1:
        return pdbFiles
      else:
        return pdbFiles[0]

    def formatPocketStrLine(self, line, numId):
      line = splitPDBLine(line)
      replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], '', 'Ve']
      pdbLine = writePDBLine(replacements)
      return pdbLine

    def mergePDBFiles(self, pdbFiles, pdbOutFile):
        atomLines, hetatmLines = '', ''
        idsDic = {}
        for pFile in pdbFiles:
            fileId = pFile.split('-')[-1].split('.')[0]
            with open(pFile) as fpdb:
                for line in fpdb:
                    if line.startswith('TITLE') and '_site_' in line:
                      numId = line.split('_site_')[1].strip()
                      idsDic[fileId] = numId
                    elif line.startswith('ATOM'):
                      atomLines += line
                    elif line.startswith('HETATM'):
                      newLine = self.formatPocketStrLine(line, numId)
                      hetatmLines += newLine

        with open(self._getPath(pdbOutFile), 'w') as f:
          f.write(atomLines)
          f.write(hetatmLines)
          f.write('\nTER')

        pdbFiles, proteinFile = self.renamePDBFiles(pdbFiles, idsDic)
        return pdbFiles, proteinFile

    def renamePDBFiles(self, pdbFiles, idsDic):
        tmpFiles = []
        for pFile in pdbFiles:
            fileId = pFile.split('-')[-1].split('.')[0]
            if fileId in idsDic:
                tmpFile = pFile.replace('-{}.pdb'.format(fileId), '-{}tmp.pdb'.format(idsDic[fileId]))
                shutil.move(pFile, tmpFile)
                tmpFiles.append(tmpFile)
            else:
                proteinFile = pFile.replace('out-{}.pdb'.format(fileId), 'protein.pdb')
                shutil.move(pFile, proteinFile)

        finalPDBFiles = []
        for tmpFile in tmpFiles:
            finalPDBFiles.append(tmpFile.replace('tmp.pdb', '.pdb'))
            shutil.move(tmpFile, finalPDBFiles[-1])
        return finalPDBFiles, proteinFile

    def getPDBName(self, pdbOut):
      if '.pdb' in pdbOut:
          pdbName = pdbOut.split('.pdb')[0]
      else:
          pdbName = pdbOut[:]
          pdbOut += '.pdb'
      return pdbName, pdbOut

    def searchOutPDBFiles(self, pdbOutFile, outDir):
      pdbFiles = []
      pdbName, _ = self.getPDBName(pdbOutFile)
      for outFile in os.listdir(outDir):
        if pdbName in outFile and '.pdb' in outFile:
          pdbFiles.append(os.path.join(outDir, outFile))
      return pdbFiles








