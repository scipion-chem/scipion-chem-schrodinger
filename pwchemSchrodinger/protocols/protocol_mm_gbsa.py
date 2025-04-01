# **************************************************************************
# *
# * Authors: Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# General imports
import os

# Scipion em imports
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb
from pyworkflow.protocol import params
import pyworkflow.object as pwobj

# Scipion chem imports
from pwchem.utils import pdbqt2other, convertToSdf, performBatchThreading, getBaseName, runInParallel
from pwchem.objects import SmallMolecule

# Plugin imports
from .. import Plugin as schrodingerPlugin
from ..utils import getNumberOfStructures

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
progLigPrep = schrodingerPlugin.getHome('ligprep')
progPrepWizard = schrodingerPlugin.getHome('utilities/prepwizard')
mmgbsaProg = schrodingerPlugin.getHome('prime_mmgbsa')

FNONE, FTHRES, FAUT = 0, 1, 2
MAEFILE_EXTENSION = '.maegz'

class ProtSchrodingerMMGBSA(EMProtocol):
    """Optimizes the docking position and calculates the binding energy"""
    _label = 'MM-GBSA'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ',
                       help="Input docked small molecules to be optimized and scored with MM-Prime GBSA")

        group = form.addGroup('Preparation')
        group.addParam('prepareLigand', params.BooleanParam, default=True, label='Prepare ligands: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Prepare target with Schrodinger LigPrep to ensure correct structure in case it was not '
                            'originally a Schrodinger molecule. It may subtly variate the ligand atom positions. '
                            'Might be necessary if docking was not performed in Schrodinger')
        group.addParam('prepareReceptor', params.BooleanParam, default=True, label='Prepare target: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Prepare receptor with Schrodinger PrepWizard to ensure correct structure in case it was '
                            'not originally a Schrodinger receptor. It may subtly variate the atom positions. '
                            'Might be necessary if docking was not performed in Schrodinger')

        group = form.addGroup('Target')
        group.addParam('flexOption', params.EnumParam, default=0, label='Flexible group: ',
                       choices=['None', 'Threshold distance', 'Flexibility estimation'],
                       help='Select the way to determine the target flexible residues.\n'
                            'None: target residues will remain rigid\n'
                            'Threshold distance: flexible residues will be determined by a distance to the ligand\n'
                            'Flexibility estimation: run a two-stage MMGBSA calculation where the second stage runs '
                            'with the subset of flexible protein residues identified by the first')
        group.addParam('flexDist', params.FloatParam, default=3.0, label='Flexible threshold (A): ',
                       condition='flexOption==1',
                       help='Treat all residues within this distance of the ligand as flexible. '
                            'By default the entire receptor is frozen')
        group.addParam('flexGroup', params.EnumParam, default=0, label='Flexible group: ',
                       choices=['Residue', 'Side', 'PolarH'], condition='flexOption==1',
                       help='Select a portion of the region defined with "flexible threshold" to be flexible.'
                            '\nResidue: Choose the entire residue.'
                            '\nSide: Choose the sidechain of each residue. '
                            '\nPolarH: Choose the polar hydrogens on each residue.')
        group.addParam('flexTargetCutOff', params.FloatParam, default=3.0, label='Flexible threshold (A): ',
                       condition='flexOption==2',
                       help='Cutoff for determining movement for target flexibility in Angstroms')

        group = form.addGroup('Ligand', expertLevel=params.LEVEL_ADVANCED,)
        group.addParam('useCharges', params.BooleanParam, default=False, label='Use ligand charges: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Use the partial charges in the input ligand file rather than from the force field')
        group.addParam('rigidLigand', params.BooleanParam, default=False, label='Rigid ligand: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Minimize the ligand as a rigid body')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        """ This function inserts all steps functions that will be run when running the protocol """
        convStep = self._insertFunctionStep('convertStep', prerequisites=[])
        mmStep = self._insertFunctionStep('runMMGBSAStep', prerequisites=[convStep])
        self._insertFunctionStep('createOutputStep', prerequisites=[mmStep])

    def convertStep(self):
        nt = self.numberOfThreads.get()
        outDir = self.getInputComplexDir()
        inMols = self.inputSmallMolecules.get()
        recMaeFile = self.prepareReceptorFile(inMols)

        fMol = inMols.getFirstItem()
        if '.mae' not in fMol.getPoseFile():
          print('Docking files are being converted to Maestro format')
          maeFiles = performBatchThreading(self.prepareLigandFile, inMols, nt)

          nSubMols, lastMol = self.getNMolsPerThread(len(maeFiles), nt), 0
          for i, nSubMol in enumerate(nSubMols):
              complexFile = os.path.join(outDir, os.path.basename(recMaeFile).
                                         replace('.mae', f'_{lastMol}-{lastMol+nSubMol-1}.mae'))
              self.runJob('zcat', '{} {} > {}'.format(recMaeFile,
                                                      ' '.join(maeFiles[lastMol:lastMol+nSubMol]), complexFile),
                          cwd=self._getTmpPath())
              lastMol = lastMol+nSubMol
        else:
          poseFiles = []
          for mol in inMols:
            molFn = mol.getPoseFile()
            molFn = molFn.split('@')[1] if '@' in molFn else molFn
            if not molFn in poseFiles:
              self.divideInThreads(molFn, nt, outDir)
              poseFiles.append(molFn)

    def runMMGBSAStep(self):
        """ This function runs the schrodinger binary file with the given params """
        nt = self.numberOfThreads.get()
        dockDir = self.getInputComplexDir()
        dockFiles = [os.path.join(dockDir, dockFile) for dockFile in os.listdir(dockDir)]
        performBatchThreading(self.performMMGBSA, dockFiles, nt, cloneItem=False)

    def createOutputStep(self):
      inMols = self.inputSmallMolecules.get()

      molDic, molIdxDic = {}, {}
      for mol in inMols:
        fnSmall = mol.getPoseFile()
        if mol.getMolClass() == 'Schrodinger':
          posId = fnSmall.split('@')[0]
          molIdxDic[int(posId)] = mol.clone()
        else:
          fnBase = getBaseName(fnSmall)
          molDic[fnBase] = mol.clone()

      outComplexFiles = {}
      for batchDir in self.getBatchDirs(complete=True):
        for file in os.listdir(batchDir):
          if file.endswith('-out.maegz'):
            maeFile = os.path.join(batchDir, file)
            csvFile = maeFile.replace(MAEFILE_EXTENSION, '.csv')
            outComplexFiles[maeFile] = csvFile

      oMols = []
      for maeFile, csvFile in outComplexFiles.items():
        csvBasesDic, csvIdxsDic = self.parseOutputCSVFile(csvFile)
        oMols += self.buildOutputMols(maeFile, molDic, csvBasesDic)
        oMols += self.buildOutputMols(maeFile, molIdxDic, csvIdxsDic)

      outputSet = inMols.createCopy(self._getPath(), copyInfo=True)
      for mol in oMols:
        outputSet.append(mol)
      self._defineOutputs(outputSmallMolecules=outputSet)


    def buildOutputMols(self, maeFile, molDic, csvDic):
      oMols = []
      for fnBase in csvDic:
        if fnBase in molDic:
          newMol = SmallMolecule()
          newMol.copy(molDic[fnBase], copyId=False)

          recFile, posFile = self.divideMaeComplex(maeFile)
          newMol.setProteinFile(recFile)
          newMol.setPoseFile(posFile)

          newMol._mmGBSA_BindEnergy = pwobj.Float(csvDic[fnBase])
          oMols.append(newMol.clone())
      return oMols


    def divideMaeComplex(self, maeFile, outDir=None):
      if not outDir:
        outDir = os.path.dirname(maeFile)

      molFile, recFile = os.path.join(outDir, getBaseName(maeFile) + '_ligand1.maegz'), \
                         os.path.join(outDir, getBaseName(maeFile) + '_receptor1.maegz')

      args = schrodingerPlugin.getMMshareDir('python/common/split_structure.py ')
      args += f'-m pdb -many_files {os.path.abspath(maeFile)} {os.path.abspath(maeFile)}'
      self.runJob(schrodingerPlugin.getHome('run'), args)
      return recFile, molFile

    # --------------------------- System functions --------------------------

    def _validate(self):
        """ Try to find warnings on define params. """
        validations = []
        pSet = self.inputSmallMolecules.get()
        if not pSet.isDocked():
            validations.append('Input molecules must be docked first.\nSet: {} has not been docked'.format(pSet))
        return validations

    # --------------------------- Parallelization functions --------------------
    def prepareLigandFile(self, mols, molLists, it):
      '''Return the corresponding MAE file for the molecule'''
      for mol in mols:
        molFile = mol.getPoseFile()
        baseName = os.path.splitext(os.path.basename(molFile))[0]
        if not molFile.endswith('.sdf'):
          # Use openbabel first to convert to sdf (usable by Schrodinger)
          molFile = convertToSdf(self, molFile)

        if self.prepareLigand.get():
          # Prepare sdf file using LigPrep
          tmpmaeFile = os.path.abspath(self._getTmpPath(baseName + '_tmp.maegz'))
          args = " -WAIT -R h -a -isd {} -omae {}".format(molFile, tmpmaeFile)
          self.runJob(progLigPrep, args, cwd=self._getTmpPath())
          maeFile = os.path.abspath(self._getTmpPath(baseName + MAEFILE_EXTENSION))
          os.rename(tmpmaeFile, maeFile)

        else:
          # Convert sdf to mae file
          maeFile = os.path.abspath(self._getTmpPath(mol.getUniqueName() + MAEFILE_EXTENSION))
          self.runJob(structConvertProg, '{} {}'.format(molFile, maeFile))

        molLists[it].append(maeFile)

    def getMMGBSAArgs(self):
      args = ''
      if self.flexOption.get() == FTHRES:
        args += f'-rflexdist {self.flexDist.get()} -rflexgroup {self.getEnumText("flexGroup").lower()} '
      elif self.flexOption.get() == FAUT:
        args += f'-target_flexibility -target_flexibility_cutoff {self.flexTargetCutOff.get()} '

      if self.useCharges.get():
        args += '-use_ligand_charges '
      if self.rigidLigand.get():
        args += '-rigid_body '
      return args

    def performMMGBSA(self, dockFiles, molLists, it):
      outDir = self._getExtraPath(f'batch_{it}')
      if not os.path.exists(outDir):
        os.mkdir(outDir)
      for dockFile in dockFiles:
        print('Performing prime mm-gbsa on : ', dockFile)
        args = f"-WAIT {os.path.abspath(dockFile)} -out_type COMPLEX "
        args += self.getMMGBSAArgs()
        self.runJob(mmgbsaProg, args, cwd=outDir)

    # --------------------------- Utils functions --------------------
    def getNMolsPerThread(self, nMols, nt):
      '''Return a list with the number of element to process per thread'''
      if nt > nMols:
        nt = nMols

      nSubMols = [0 for _ in range(nt)]
      for i in range(nMols):
          nSubMols[i % nt] += 1
      return nSubMols

    def divideInThreads(self, dockFile, nt, outDir):
      '''Divide the complex file into several containing the same receptor (first mol) so the following
      process of these complex files can be parallelized in nt threads'''
      nMols = getNumberOfStructures(dockFile) - 1
      if nMols < nt:
        nt = nMols

      if nt > 1:
        nSubMols = self.getNMolsPerThread(nMols, nt)
        firstMol, maeFiles = 2, []
        for nSubMol in nSubMols:
          lastMol = firstMol + nSubMol
          maeFile = os.path.basename(dockFile).replace('.mae', f'_{firstMol}-{lastMol-1}.mae')
          args = f' -n "1, {firstMol}:{lastMol-1}" {os.path.abspath(dockFile)} {maeFile}'
          self.runJob(structConvertProg, args, cwd=outDir)
          firstMol = lastMol
          maeFiles.append(os.path.join(outDir, maeFile))
      else:
        maeFiles = [os.path.join(outDir, os.path.basename(dockFile))]
        os.link(dockFile, maeFiles[0])
      return maeFiles

    def getMaestroDockGroups(self, molSet):
      '''Return the complex (receptor-molecules) files in a molSet'''
      gFiles = []
      for mol in molSet:
        g = mol.getPoseFile().split('@')[1]
        if g not in gFiles:
          gFiles.append(g)
      return gFiles

    def getInputComplexDir(self):
      cDir = self._getExtraPath('inputComplex')
      if not os.path.exists(cDir):
        os.mkdir(cDir)
      return os.path.abspath(cDir)

    def getInputComplexFiles(self):
      cFiles = []
      cDir = self.getInputComplexDir()
      for file in os.listdir(cDir):
        cFiles.append(os.path.join(cDir, file))
      return cFiles

    def getPdbFile(self, molSet):
      proteinFile = molSet.getProteinFile()
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))

      if inExt == '.pdb':
        return os.path.abspath(proteinFile)
      else:
        pdbFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))
        if inExt == '.pdbqt':
          pdbqt2other(self, proteinFile, pdbFile)
        else:
          toPdb(proteinFile, pdbFile)
      return os.path.abspath(pdbFile)

    def prepareReceptorFile(self, molSet):
      '''Return the corresponding MAE file for the receptor the molecule is docked to'''
      mol = molSet.getFirstItem()
      if hasattr(mol, 'structFile'):
        targetMaeFile = molSet.structFile
      else:
        pdbFile = self.getPdbFile(molSet)
        targetName = os.path.splitext(os.path.basename(pdbFile))[0]
        targetMaeFile = os.path.abspath(self._getTmpPath(targetName + MAEFILE_EXTENSION))
        if self.prepareReceptor.get():
          args = '-WAIT -noprotassign -noimpref -noepik '
          args += '%s %s' % (os.path.abspath(pdbFile), os.path.abspath(targetMaeFile))
          self.runJob(progPrepWizard, args, cwd=self._getTmpPath())
        else:
          self.runJob(structConvertProg, '{} {}'.format(pdbFile, targetMaeFile), cwd=self._getTmpPath())
      return targetMaeFile

    def getBatchDirs(self, complete=False):
      batchDirs = []
      for bdir in os.listdir(self._getExtraPath()):
        if bdir.startswith('batch_'):
          if not complete:
            batchDirs.append(self._getExtraPath(bdir))
          else:
            add = False
            for file in os.listdir(self._getExtraPath(bdir)):
              if file.endswith("-out.maegz"):
                add = True
                break

            if add:
              batchDirs.append(self._getExtraPath(bdir))
            else:
              print('No good poses found in batch ' + bdir.split('_')[1])
      return batchDirs

    def parseOutputCSVFile(self, csvFile):
      csvDic, csvIdxsDic = {}, {}
      posFIdx = csvFile.split('_')[-1].split('-')[0]
      with open(csvFile) as f:
        f.readline()
        for i, line in enumerate(f):
          fnBase, mmgbsaBindEnergy = line.split(',')[:2]
          fnBase = getBaseName(fnBase)
          csvDic[fnBase] = float(mmgbsaBindEnergy)
          csvIdxsDic[int(posFIdx) + i - 1] = float(mmgbsaBindEnergy)
      return csvDic, csvIdxsDic

    def getOriginalReceptorFile(self):
        return self.inputSmallMolecules.get().getProteinFile()
