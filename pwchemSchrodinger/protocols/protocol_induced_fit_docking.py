# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
import os, json, subprocess, shutil, time, gzip

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import Group, EnumParam, PointerParam, StringParam, FloatParam, \
  TextParam, IntParam, LEVEL_ADVANCED, BooleanParam, Line, STEPS_PARALLEL
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

from pwchem.utils import createMSJDic, calculate_centerMass, performBatchThreading, getBaseName, natural_sort
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC

# Plugin imports
from .. import Plugin as schrodingerPlugin
from ..protocols.protocol_glide_docking import ProtSchrodingerGlideDocking

prepWizardProg = schrodingerPlugin.getHome('utilities/prepwizard')
structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
propListerProg = schrodingerPlugin.getHome('utilities/proplister')

ifdProg = schrodingerPlugin.getHome('ifd')
runPath = schrodingerPlugin.getHome('run')

splitProg = schrodingerPlugin.getMMshareDir('python/common/split_structure.py')
closeScript = schrodingerPlugin.getPluginHome('scripts/getCloseResidues.py')


AS, POCKET, GRID = 0, 1, 2

#  0       1     2       3      4      5      6       7     8     9      10     11    12    13    14
COMPILE, GLIDE, IDOCK, PPREP, PFLEX, PENER, PHELIX, PLOOP, PMIN, PREF, PSIDE, SCORE, SORT, TRIM, VDW = list(range(15))

class ProtSchrodingerIFD(ProtSchrodingerGlideDocking):
  """Perform an Induced Fit Docking protocol using the user defined set of stages"""
  _label = 'Induced Fit Docking (IFD)'

  stageTypes = {'Compile residues': 'COMPILE_RESIDUE_LIST',
                'Glide': 'GLIDE_DOCKING2', 'Initial docking': 'INITIAL_DOCKING', 'Receptor preparation': 'PPREP',
                'Predict flexibility': 'PREDICT_FLEXIBILITY', 'Energy calculation': 'PRIME_ECALC',
                'Helix refinement': 'PRIME_HELIX', 'Loop refinement': 'PRIME_LOOP',
                'Residues minimization': 'PRIME_MINIMIZATION', 'Residues refinement': 'PRIME_REFINEMENT',
                'Side-chain prediction': 'PRIME_SIDECHAIN',
                'Score poses': 'SCORING', 'Sort and filter': 'SORT_AND_FILTER',
                'Trim side-chains': 'TRIM_SIDECHAINS', 'VDW Scaling': 'VDW_SCALING'}

  _omitParamNames = ['runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                     'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi',
                     'inputLibrary', 'fromPockets', 'inputAtomStruct', 'radius', 'inputStructROIs',
                     'inputGridSet']

  NONE, STAN, EXTE = 0, 1, 2

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = STEPS_PARALLEL

  def getStageParamsDic(self, type='All'):
    '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
    Type'''
    paramsDic = {}
    for paramName, param in self._definition.iterAllParams():
      if paramName not in self._omitParamNames and not isinstance(param, Group):
        if type == 'All' or (type == 'Enum' and isinstance(param, EnumParam)) or \
                (type == 'Normal' and not isinstance(param, EnumParam)):
          paramsDic[paramName] = param
    return paramsDic

  def getStageStr(self, sIdx):
    return list(self.stageTypes.keys())[sIdx]

  def _defineParams(self, form):
    # Defining costant variables
    chainListHelp = 'Use the wizard to select the chain you want'
    residueListHelp = 'Use the wizard to select the residues you want'
    cutoffDistanceLabel = 'Cutoff distance:'
    isGridStr = 'fromPockets!=2'

    form.addSection(label=Message.LABEL_INPUT)
    self._defineGlideReceptorParams(form)
    group = form.addGroup('Ligands')
    group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules:', help='Input small molecules to be docked with IFD')
    form.addParam('defSteps', EnumParam, label='Default MD workflows: ',
                  choices=['None', 'Standard', 'Extended sampling'], default=self.NONE,
                  help='Choose and add with the wizard some default IFD steps which will be reflected into '
                       'the summary')

    form = self._defineGridSection(form, condition=isGridStr)
    self._defineInnerGridParams(form, condition=isGridStr)
    self._defineOuterGridParams(form, condition=isGridStr)

    form.addSection('IFD stages')
    group = form.addGroup('Add Stage')
    group.addParam('stageType', EnumParam, default=0, label='Stage type: ', choices=list(self.stageTypes.keys()),
                   help='Type of the next stage to add to the protocol. For details about the different stages, '
                        'check the Schrodinger docs in {}'.
                   format(schrodingerPlugin.getHome('docs/inducedfit_command_reference/ifd_command_infile.htm')))

    # COMPILE_RESIDUE_LIST stage parameters
    group.addParam('residueCenter', StringParam, default='Z:999', label='Center to determine residues: ',
                   condition=f'stageType == {COMPILE}',
                   help='List of residues from which to measure the cutoff distance. Default: Z:999, which is the '
                        'default for the ligand.')
    group.addParam('residuesCutoff', FloatParam, label=cutoffDistanceLabel,
                   default=5.0, condition=f'stageType == {COMPILE}',
                   help='Cutoff distance (in angstroms) from the ligand pose, within which residues that have any '
                        'atoms are included in the refinement list')
    group.addParam('residuesAdd', TextParam, default='', width=120, label='Compile residues to add: ',
                   condition=f'stageType == {COMPILE}', expertLevel=LEVEL_ADVANCED,
                   help='Comma-separated list of residues to add outside of the cutoff to the refinement list')
    group.addParam('residuesOmit', TextParam, default='', width=120, label='Compile residues to omit: ',
                   condition=f'stageType == {COMPILE}', expertLevel=LEVEL_ADVANCED,
                   help='Comma-separated list of residues to omit inside of the cutoff to the refinement list')

    # GLIDE_DOCKING2 and INITIAL_DOCKING stage parameters
    group.addParam('selfDock', BooleanParam, label='Dock from previous dock: ',
                   default=False, condition=f'stageType in [{GLIDE}, {IDOCK}]',
                   help='Whether to dock from previously docked ligands or not (needed for second round docking)')
    self._defineGeneralGlideParams(form, condition=f'stageType in [{GLIDE}, {IDOCK}]')
    self._defineLigandParams(form, condition=f'stageType in [{GLIDE}, {IDOCK}]')
    self._defineOutputGlideParams(form, condition=f'stageType in [{GLIDE}, {IDOCK}]')

    # PPREP stage parameters
    group.addParam('convergenceRMSD', FloatParam, label='Convergence RMSD: ',
                   default=0.2, condition=f'stageType == {PPREP}',
                   help='RMSD value which specifies the convergence threshold for the constrained minimization of '
                        'the Glide protein preparation. The minimization ends when the RMSD is less than or equal to '
                        'value.')

    # PRIME_HELIX and PRIME_LOOP stage parameters
    group.addParam('residuesHelix', StringParam, default='',
                   label='Helix residues: ', condition=f'stageType in [{PHELIX}]',
                   help='Residue of the helix. It must be a subsequence of the residue selection.')
    group.addParam('helixLoopCutoff', FloatParam, label=cutoffDistanceLabel,
                   default=5.0, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='HELIX: Distance cutoff for prediction of side chains close to helix.\n'
                        'LOOP: Threshold for inclusion of residues for side-chain refinement. Any residues with atoms '
                        'within this distance are included')
    group.addParam('maxEnergyGap', FloatParam, label='Max. energy gap: ',
                   default=10000.0, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='Maximum energy gap for saved structures.'
                        'Energy threshold for predicted loop structures (in kcal/mol). Structures are discarded if '
                        'their energy is more than this amount above the lowest-energy structure.')
    group.addParam('maxStructures', IntParam, label='Max. number of structures: ',
                   default=1000, condition=f'stageType in [{PHELIX}, {PLOOP}]',
                   help='Maximum number of structures to store')
    group.addParam('includeResList', BooleanParam, label='Include compiled residues: ',
                   default=False, condition=f'stageType in [{PLOOP}]',
                   help='Include the residues from COMPILE_RESIDUE_LIST for side-chain refinement.')

    # PRIME_REFINEMENT and PRIME_MINIMIZATION stage parameters
    group.addParam('nMinPasses', IntParam, label='Number of refinement passes: ',
                   default=1, condition=f'stageType in [{PREF}]',
                   help='In Prime refinement, the side chains of the residue list compiled previously are optimized, '
                        'then the residues are minimized. This parameter controls the number of times to repeat those '
                        'two steps')

    # SCORING stage parameters
    group.addParam('scoreName', StringParam, default='',
                   label='Name of property to add to Maestro files: ', condition=f'stageType in [{SCORE}]',
                   help='Name of property to add to Maestro files')
    group.addParam('scoreTerms', TextParam, label='Scoring terms: ', default='', width=120,
                   condition=f'stageType in [{SCORE}]',
                   help='Add terms to the scoring function. You can include multiple TERM keywords to define the '
                        'scoring function (one per line). Format: coeff,property,stage, where coeff is the '
                        'coefficient, property is the property from the Maestro output file, and stage is the index '
                        'of the property-generating stage, counting backwards in the input file with 0 for the '
                        'previous stage. \nSee for details: '
                        f'{schrodingerPlugin.getHome("docs/inducedfit_command_reference/ifd_command_infile_scoring.htm")}')

    # SORT_AND_FILTER stage parameters
    group.addParam('sortType', EnumParam, label='Sort by: ', default=1, display=EnumParam.DISPLAY_HLIST,
                   choices=['Ligand', 'Pose'], condition=f'stageType == {SORT}',
                   help='Whether to sort and filter by Ligands or Poses')
    group.addParam('ligandFilter', StringParam, default='', label='Property to filter ligands: ',
                   condition=f'stageType in [{SORT}] and sortType == 0',
                   help='Name of Maestro property for filtering ligands, for example, r_psp_Prime_Energy')
    group.addParam('ligandKeep', StringParam, default='', label='Threshold on property for filtering ligands: ',
                   condition=f'stageType in [{SORT}] and sortType == 0',
                   help='Threshold on property for filtering ligands. The syntax is as follows: '
                        'n% Keep the n% of poses with the lowest property values'
                        'n# Keep the n poses with the lowest property values'
                        'n Keep poses with property values within n of the lowest value.')
    group.addParam('poseFilter', StringParam, default='', label='Property to filter poses: ',
                   condition=f'stageType in [{SORT}] and sortType == 1',
                   help='Name of Maestro property for filtering poses, for example, r_psp_Prime_Energy')
    group.addParam('poseKeep', StringParam, default='', label='Threshold on property for filtering poses: ',
                   condition=f'stageType in [{SORT}] and sortType == 1',
                   help='Threshold on property for filtering poses. The syntax is as follows: '
                        'n% Keep the n% of poses with the lowest property values'
                        'n# Keep the n poses with the lowest property values'
                        'n Keep poses with property values within n of the lowest value.')

    # TRIM_SIDECHAINS stage parameters
    group.addParam('trimResidues', EnumParam, label='Residues to mutate to Ala: ', default=1,
                   display=EnumParam.DISPLAY_HLIST,
                   choices=['Manual specification', 'Automatic'], condition=f'stageType == {TRIM}',
                   help='Residues for temporary mutation to alanine. Either a list of residues must be given, or AUTO.')
    group.addParam('trimMethod', EnumParam, label='Method for automatic specification of residues: ', default=0,
                   choices=['BFACTOR', 'FLEXIBILITY'], condition=f'stageType == {TRIM} and trimResidues == 1',
                   help='Method for automatic specification of residues. '
                        'BFACTOR: Mutate residues near the ligand with the highest B factors, subject to a cutoff on '
                        'the B-factor and the number of residues.\n'
                        'FLEXIBILITY: Mutate residues for which other rotamers exist that don’t have steric '
                        'clashes with the rest of the protein.')
    group.addParam('cutOffBFactor', FloatParam, label='B-factor cutoff: ',
                   default=40.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 0',
                   help='B-factor cutoff above which residues are selected for mutation.')
    group.addParam('maxBFactor', IntParam, label='Maximum number of residues to mutate: ',
                   default=3, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 0',
                   help='Maximum number of residues to mutate. If more than N residues have B-factors greater than the '
                        'value of BFACTOR_CUTOFF, the N residues with the highest B-factors are selected')
    group.addParam('maxFlexibility', FloatParam, label='Max flexibility cutoff (Å): ',
                   default=5.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 1',
                   help='If any heavy atom in the side chain is capable of moving more than this distance, the side '
                        'chain is considered flexible and is mutated')
    group.addParam('resolFlexibility', FloatParam, label='Granularity of rotamer search: ',
                   default=30.0, condition=f'stageType == {TRIM} and trimResidues == 1 and trimMethod == 1',
                   help='Granularity of the rotamer search for steric clashes, in degrees.')

    group.addParam('residuesTrim', TextParam, default='', width=120, label='Residues to trim: ',
                   condition=f'stageType == {TRIM} and trimResidues == 0',
                   help='List of residues to trim. Use the wizard to paste the selected residues.')

    group.addParam('extraParams', TextParam, label='Extra parameters: ', default='', width=120,
                   expertLevel=LEVEL_ADVANCED,
                   help='Any extra parameters to add in the current stage. Add each parameter in a different line as:\n'
                        'PARAMETER-NAME PARAMETER-VALUE')

    group = form.addGroup('Residue selection', condition=f'(stageType == {COMPILE}) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)')
    group.addParam('selChain', StringParam, default='',
                   label='Select chain: ',
                   condition=f'(stageType == {COMPILE}) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)',
                   help=chainListHelp)
    group.addParam('selResidue', StringParam, default='',
                   label='Select residues: ',
                   condition=f'(stageType == {COMPILE}) or (stageType in [{PHELIX}, {PLOOP}]) or (stageType == {TRIM} and trimResidues == 0)',
                   help=residueListHelp)

    group = form.addGroup('Summary')
    group.addParam('insertStep', StringParam, default='',
                   label='Insert stage number: ',
                   help='Insert the defined protocol stage into the workflow on the defined position.\n'
                        'The default (when empty) is the last position')
    group.addParam('summarySteps', TextParam, width=120, readOnly=True,
                   label='Summary of stages', help='Summary of the defined stages. \nManual modification will have no '
                                                   'effect, use the wizards to add / delete the stages')
    group.addParam('deleteStep', StringParam, default='',
                   label='Delete stage number: ',
                   help='Delete the stage of the specified index from the workflow.')
    group.addParam('watchStep', StringParam, default='',
                   label='Watch stage number: ',
                   help='''Watch the parameters stage of the specified index from the workflow..\n
                       This might be useful if you want to change some parameters of a predefined stage.\n
                       However, the parameters are not changed until you add the new stage (and probably\n
                       you may want to delete the previous unchanged stage)''')
    group.addParam('workFlowSteps', StringParam, label='User transparent', condition='False')

    form.addParallelSection(threads=4, mpi=1)

  ###################### SETTINGS #########################

  def _validate(self):
    return []

  def _summary(self):
    fnSummary = self._getExtraPath("summary.txt")
    if not os.path.exists(fnSummary):
        summary = ["No summary information yet."]
    else:
        fhSummary = open(fnSummary, "r")
        summary = []
        for line in fhSummary.readlines():
            summary.append(line.rstrip())
        fhSummary.close()
    return summary

  def getnThreads(self):
    '''Get the number of threads available for each pocket execution'''
    if self.fromPockets.get() == AS:
      nThreads = [self.numberOfThreads.get()]

    elif self.fromPockets.get() == POCKET:
      nP = len(self.inputStructROIs.get())
      nThreads = [0] * nP
      for i in range(self.numberOfThreads.get()):
        nThreads[i % nP] += 1

    elif self.fromPockets.get() == GRID:
      nP = len(self.inputGridSet.get())
      nThreads = [0] * nP
      for i in range(self.numberOfThreads.get()):
        nThreads[i % nP] += 1

    nThreads = list(map(lambda x: x if x != 0 else 1, nThreads))
    return nThreads

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    convStep = self._insertFunctionStep('convertStep', prerequisites=[])
    cSteps = [convStep]

    idfSteps, nt = [], self.getnThreads()
    if self.fromPockets.get() == AS:
      dStep = self._insertFunctionStep('ifdStep', self.inputAtomStruct.get(), nt[0], prerequisites=cSteps)
      idfSteps.append(dStep)

    elif self.fromPockets.get() == POCKET:
      for i, pocket in enumerate(self.inputStructROIs.get()):
        dStep = self._insertFunctionStep('ifdStep', pocket.clone(), nt[i], prerequisites=cSteps)
        idfSteps.append(dStep)

    elif self.fromPockets.get() == GRID:
      for i, pocket in enumerate(self.inputGridSet.get()):
        dStep = self._insertFunctionStep('ifdStep', pocket.clone(), nt[i], prerequisites=cSteps)
        idfSteps.append(dStep)

    self._insertFunctionStep('createOutputStep', prerequisites=idfSteps)

  def createLigandsFile(self, ligSet, molLists, it):
    curAllLigandsFile = self.getAllLigandsFile(suffix=it)
    with open(curAllLigandsFile, 'w') as fh:
      for small in ligSet:
        fnSmall = small.getFileName()
        if '.mae' not in fnSmall:
          fnSmall = self.convert2mae(fnSmall, it)
        with open(fnSmall) as fhLigand:
          fh.write(fhLigand.read())

  def convertStep(self):
    nt = self.numberOfThreads.get()
    # Convert receptor
    inFile = self.getOriginalReceptorFile()

    if '.mae' not in inFile:
      maeFile = self._getExtraPath('inputReceptor.maegz')
      args = ' -WAIT -noprotassign -noimpref -noepik {} {}'. \
        format(os.path.abspath(inFile), os.path.abspath(maeFile))
      self.runJob(prepWizardProg, args, cwd=self._getExtraPath())
    else:
      ext = os.path.splitext(inFile)[1]
      shutil.copy(inFile, self._getExtraPath('inputReceptor{}'.format(ext)))

    # Create ligand files
    ligSet = self.inputLibrary.get()
    performBatchThreading(self.createLigandsFile, ligSet, nt)

    allLigandsFile = self.getAllLigandsFile()
    with open(allLigandsFile, 'w') as f:
      for it in range(nt):
        lFile = self.getAllLigandsFile(suffix=it)
        if os.path.exists(lFile):
          with open(lFile) as fLig:
            f.write(fLig.read())

  def ifdStep(self, pocket, nt):
    pId = pocket.getObjId()
    fnGridDir = self.getGridDir(pId)
    os.mkdir(fnGridDir)

    msjStr = self.buildIFDStr(pocket)
    msjFile = self.getMSJFile(fnGridDir)
    with open(msjFile, 'w') as f:
      f.write(f'INPUT_FILE	{self.getInputReceptorFile()}\n\n{msjStr}')
    self.createGUISummary()

    args = f"-WAIT -SUBHOST localhost -NGLIDECPU {nt} -NPRIMECPU {nt} ifd.inp"
    self.runJob(ifdProg, args, cwd=fnGridDir)

  def performOutputParsing(self, gridDirs, molLists, it, smallDict):
    allSmallList = []
    scoreNames = self.getAllScoreNames()
    scoreNamesStr = '' if not scoreNames else f"-p {' -p '.join(scoreNames)} "
    for gridDir in gridDirs:
      gridId = gridDir.split('_')[-1]
      outFile, propFile = os.path.join(gridDir, 'ifd-out.maegz'), os.path.join(gridDir, 'report_properties.csv')

      subprocess.run(f'{propListerProg} -p "s_m_title" -p "r_i_docking_score" {scoreNamesStr}'
                     f'-c -o {propFile} {outFile}', cwd=gridDir, check=True, capture_output=True, text=True, shell=True)

      smallList, posDic = [], {}
      with open(propFile) as f:
        f.readline()
        for i, line in enumerate(f.readlines()):
          tokens = line.strip().split(',')
          title, glideEnergy = tokens[:2]
          posDic[title] = 1 if title not in posDic else posDic[title] + 1

          recFile, posFile = self.splitComplexFiles(outFile, gridDir, title, posDic[title])

          small = SmallMolecule()
          small.copy(smallDict[title], copyId=False)
          small.setMolClass('Schrodinger')
          small.setDockId(self.getObjId())
          small.setGridId(gridId)

          small._energy = pwobj.Float(glideEnergy)
          for j, sName in enumerate(scoreNames):
            setattr(small, sName, pwobj.Float(tokens[j+2]))
          small.setProteinFile(recFile)
          small.setPoseFile(posFile)
          small.setPoseId(posDic[title])

          smallList.append(small)

      allSmallList += smallList
    molLists[it] = allSmallList

  def createOutputStep(self):
      nt = self.numberOfThreads.get()
      smallDict = {}
      for mol in self.inputLibrary.get():
          fnBase = mol.getMolName()
          if fnBase not in smallDict:
              smallDict[fnBase] = mol.clone()

      gridDirs = self.getGridDirs(complete=True)
      allSmallList = performBatchThreading(self.performOutputParsing, gridDirs, nt, cloneItem=False,
                                           smallDict=smallDict)

      outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
      for small in allSmallList:
          outputSet.append(small)

      outputSet.setDocked(True)
      outputSet.proteinFile.set(self.getOriginalReceptorFile())
      outputSet.structFile = pwobj.String(self.getInputMaeFile())
      self._defineOutputs(outputSmallMolecules=outputSet)


  ##################### UTILS ###########################

  def getInputReceptorFile(self):
    for f in os.listdir(self._getExtraPath()):
      if 'inputReceptor.mae' in f:
        break
    return os.path.abspath(self._getExtraPath(f))

  def getAllLigandsFile(self, suffix=''):
      return os.path.abspath(self._getTmpPath('allMoleculesFile{}.mae'.format(suffix)))

  def convert2mae(self, fnSmall, it):
    baseName = os.path.splitext(os.path.basename(fnSmall))[0]
    outFile = os.path.abspath(self._getTmpPath('{}.mae'.format(baseName)))
    if fnSmall.endswith('.pdbqt'):
      # Manage files from autodock: 1) Convert to readable by schro (SDF). 2) correct preparation.
      # 3) Switch to mol2 to manage atom labels
      outDir = os.path.abspath(self._getTmpPath())
      args = ' -i "{}" -of sdf --outputDir "{}" --outputName {}_AD4'.format(os.path.abspath(fnSmall),
                                                                            os.path.abspath(outDir), baseName)
      pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir, popen=True)
      auxFile = os.path.abspath(os.path.join(outDir, '{}_AD4.sdf'.format(baseName)))
      fnSmall = auxFile.replace('_AD4.sdf', '_aux.sdf')

    args = "{} {}".format(fnSmall, outFile)
    subprocess.check_call([structConvertProg, *args.split()])

    while not os.path.exists(outFile):
      time.sleep(0.2)
    return outFile

  def getGridDirs(self, complete=False):
    gridDirs = []
    for dir in os.listdir(self._getExtraPath()):
      if dir.startswith('grid_'):
        if not complete:
          gridDirs.append(os.path.abspath(self._getExtraPath(dir)))
        else:
          gridId = dir.split('_')[1]
          if os.path.exists(self._getExtraPath(dir, "ifd-out.maegz")):
            gridDirs.append(os.path.abspath(self._getExtraPath(dir)))
          else:
            print('No output found in grid ' + gridId)
    return gridDirs

  def getMSJFile(self, fnGridDir):
    return os.path.join(fnGridDir, 'ifd.inp')

  def addDefaultForMissing(self, msjDic):
    '''Add default values for missing parameters in the msjDic'''
    paramDic = self.getStageParamsDic()
    for pName, pVal in paramDic.items():
      if pName not in msjDic and not isinstance(pVal, Line):
        if isinstance(pVal, EnumParam):
          msjDic[pName] = pVal.choices[int(pVal.default.get())]
        else:
          msjDic[pName] = pVal.default.get()
    return msjDic

  def buildIFDStr(self, pocket=None):
    '''Build the .msj (file used by IFD to specify the jobs performed by Schrodinger)
        defining the input parameters'''

    msjStr = ''
    if self.workFlowSteps.get() in ['', None]:
      msjDic = createMSJDic(self)
      msjStr += self.buildIFDStepStr(msjDic, pocket)
    else:
      workSteps = self.workFlowSteps.get().split('\n')
      if '' in workSteps:
        workSteps.remove('')
      for wStep in workSteps:
        msjDic = eval(wStep)
        msjStr += self.buildIFDStepStr(msjDic, pocket)

    return msjStr

  def countSteps(self):
    stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
    steps = stepsStr.split('\n')
    return len(steps) - 1

  def createGUISummary(self):
    with open(self._getExtraPath("summary.txt"), 'w') as f:
      if self.workFlowSteps.get():
        f.write(self.createSummary())
      else:
        f.write(self.createSummary(createMSJDic(self)))

  def createSummary(self, msjDic=None):
    '''Creates the displayed summary from the internal state of the steps'''
    sumStr = ''
    if not msjDic:
      for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
        if dicLine != '':
          msjDic = eval(dicLine)
          msjDic = self.addDefaultForMissing(msjDic)

          sumStr += f'{i+1}) {self.buildSummaryLine(msjDic)}'

    else:
      msjDic = self.addDefaultForMissing(msjDic)
      sumStr += self.buildSummaryLine(msjDic)
    return sumStr

  def buildSummaryLine(self, msjDic):
    sumStr = ''
    sType = msjDic["stageType"]
    sumStr += sType
    if sType == self.getStageStr(COMPILE):
      sumStr += f': {msjDic["residuesCutoff"]} Å from {msjDic["residueCenter"]}'

    elif sType in [self.getStageStr(IDOCK), self.getStageStr(GLIDE)]:
      sumStr += f': {msjDic["dockingMethod"]} with {msjDic["dockingPrecision"]} precision ' \
                f'up to {msjDic["posesPerLig"]} poses per ligand'

    elif sType == self.getStageStr(PPREP):
      sumStr += f': {msjDic["convergenceRMSD"]} convergence RMSD'

    elif sType in [self.getStageStr(PHELIX), self.getStageStr(PLOOP)]:
      sumStr += f': up to {msjDic["helixLoopCutoff"]} Å from {msjDic["selResidue"]}'

    elif sType in [self.getStageStr(PREF)]:
      sumStr += f': {msjDic["nMinPasses"]} min. passes'

    elif sType in [self.getStageStr(SORT)]:
      if msjDic["sortType"] == 0:
        sumStr += f': keep {msjDic["ligandKeep"]} for {msjDic["ligandFilter"]} by ligand'
      else:
        sumStr += f': keep {msjDic["poseKeep"]} for {msjDic["poseFilter"]} by pose'

    elif sType in [self.getStageStr(TRIM)]:
      if msjDic["trimResidues"] == 'Manual specification':
        sumStr += ': selected residues'
      else:
        sumStr += f': by {msjDic["trimMethod"]}'

    elif sType in [self.getStageStr(SCORE)]:
        sumStr += f': {msjDic["scoreName"]}'

    return sumStr + '\n'

  def getResidueIdxs(self, resSelectionJson):
    resDic = json.loads(resSelectionJson)
    chain, idxs = resDic['chain'], resDic['index'].split('-')
    return [f'{chain}:{idx}' for idx in range(int(idxs[0]), int(idxs[1])+1)]

  def buildIFDStepStr(self, msjDic, pocket=None):
    sType = msjDic["stageType"]
    idfStr = f'STAGE {self.stageTypes[sType]}\n'
    if sType == self.getStageStr(COMPILE):
      idfStr += f'  CENTER {msjDic["residueCenter"]}\n  DISTANCE_CUTOFF	{msjDic["residuesCutoff"]}\n'
      addStr, omitStr = msjDic['residuesAdd'], msjDic['residuesOmit']
      if addStr.strip():
        resIds = self.getResidueIdxs(addStr)
        idfStr += f'  RESIDUES_TO_ADD {", ".join(resIds)}\n'
      if omitStr.strip():
        resIds = self.getResidueIdxs(omitStr)
        idfStr += f'  RESIDUES_TO_OMIT {", ".join(resIds)}\n'

    elif sType in [self.getStageStr(GLIDE), self.getStageStr(IDOCK)]:
      usePrecision = sType == self.getStageStr(GLIDE)
      idfStr += self.getBindingSiteStr(msjDic, pocket)
      idfStr += self.getBoxDimensionsStr(msjDic, pocket)
      l2dock = 'self' if msjDic['selfDock'] in ['True', True] else 'all'
      idfStr += f'  LIGAND_FILE {self.getAllLigandsFile()}\n  LIGANDS_TO_DOCK {l2dock}\n'
      idfStr += self.getGridArgsStr(msjDic)
      idfStr += self.getDockArgsStr(msjDic, usePrecision)

    elif sType == self.getStageStr(PPREP):
      idfStr += f'  RMSD {msjDic["convergenceRMSD"]}\n'

    elif sType in [self.getStageStr(PFLEX)]:
      idfStr += self.getBindingSiteStr(msjDic, pocket)

    elif sType in [self.getStageStr(PHELIX), self.getStageStr(PLOOP)]:
      resIds = self.getResidueIdxs(msjDic["selResidue"])
      idfStr += f'  START_RESIDUE {resIds[0]}\n  END_RESIDUE {resIds[-1]}\n' \
                f'  DISTANCE_CUTOFF {msjDic["helixLoopCutoff"]}\n' \
                f'  MAX_ENERGY_GAP {msjDic["maxEnergyGap"]}\n' \
                f'  MAX_STRUCTURES {msjDic["maxStructures"]}\n'

      if sType in [self.getStageStr(PHELIX)]:
        resIds = self.getResidueIdxs(msjDic["residuesHelix"])
        idfStr += f'  START_HELIX_RESIDUE {resIds[0]}\n  END_HELIX_RESIDUE {resIds[-1]}\n'
      else:
        idfStr += f'  INCLUDE_RESIDUE_LIST {str(msjDic["includeResList"]).upper()}\n'

    elif sType in [self.getStageStr(PREF)]:
      idfStr += f'  NUMBER_OF_PASSES {msjDic["nMinPasses"]}\n'

    elif sType in [self.getStageStr(SCORE)]:
      termLines = msjDic["scoreTerms"].split("\n")
      termStr = '\n'.join(termLines)
      idfStr += f'  SCORE_NAME {msjDic["scoreName"]}\n{termStr}\n'

    elif sType in [self.getStageStr(SORT)]:
      if msjDic["sortType"] == 0:
        idfStr += f'  LIGAND_FILTER {msjDic["ligandFilter"]}\n  LIGAND_KEEP {msjDic["ligandKeep"]}\n'
      else:
        idfStr += f'  POSE_FILTER {msjDic["poseFilter"]}\n  POSE_KEEP {msjDic["poseKeep"]}\n'

    elif sType in [self.getStageStr(TRIM)]:
      if msjDic["trimResidues"] == 'Manual specification':
        resIds = self.getResidueIdxs(msjDic["residuesTrim"])
        idfStr += f'  RESIDUES {", ".join(resIds)}\n'
      else:
        idfStr += f'  RESIDUES AUTO\n  METHOD {msjDic["trimMethod"]}\n'
        if msjDic["trimMethod"] == 'BFACTOR':
          idfStr += f'  BFACTOR_CUTOFF {msjDic["cutOffBFactor"]}\n  MAX_RESIDUES {msjDic["maxBFactor"]}\n'
        else:
          idfStr += f'  MAX_FLEXIBILITY {msjDic["maxFlexibility"]}\n  RESOLUTION {msjDic["resolFlexibility"]}\n'

    elif sType in [self.getStageStr(VDW)]:
      idfStr += self.getBindingSiteStr(msjDic, pocket)


    extraParams = msjDic['extraParams']
    if extraParams.strip():
      idfStr += extraParams
    return idfStr + '\n'

  def getGridArgsStr(self, msjDic):
    gridDic = {}
    return self.dic2StrArgs(gridDic, toAdd='GRIDGEN_')

  def getDockArgsStr(self, msjDic, usePrecision):
    dockDic = {}
    dMethod = msjDic['dockingMethod'].split('(')[1].replace(')', '')
    dockDic['DOCKING_METHOD'] = dMethod if dMethod in ['confgen', 'rigid'] else 'confgen'

    dPrecision = msjDic['dockingPrecision'].split('(')[1].replace(')', '')
    if usePrecision:
      dockDic['PRECISION'] = dPrecision

    dockDic['POSES_PER_LIG'] = msjDic['posesPerLig']
    return self.dic2StrArgs(dockDic, toAdd='DOCKING_')

  def getBoxDimensionsStr(self, msjDic, pocket):
    gridDic = {}
    if self.fromPockets.get() == GRID:
      gridDic['INNERBOX'] = list(map(str, pocket.getInnerBox()))
      gridDic['OUTERBOX'] = list(map(str, pocket.getOuterBox()))
    else:
      if msjDic["innerAction"] == 0:
        gridDic['INNERBOX'] = f'{msjDic["innerX"], msjDic["innerY"], msjDic["innerZ"]}'
      else:
        diam = self.getDiameter(msjDic, pocket)
        gridDic['INNERBOX'] = f'{diam * float(msjDic["diameterNin"])}'

      if msjDic["outerAction"] == 0:
        gridDic['OUTERBOX'] = f'{msjDic["outerX"], msjDic["outerY"], msjDic["outerZ"]}'
      else:
        diam = self.getDiameter(msjDic, pocket)
        gridDic['OUTERBOX'] = f'{diam * float(msjDic["diameterNout"])}'
    return self.dic2StrArgs(gridDic)

  def getBindingSiteStr(self, msjDic, pocket=None):
    gridDic = {}
    if msjDic['selfDock'] in [True, 'True']:
      gridDic['BINDING_SITE'] = 'ligand Z:999'
    else:
      gridDic['BINDING_SITE'] = 'residues ' + self.getBindingSiteCenter(pocket)

    return self.dic2StrArgs(gridDic)

  def getCloseResidues(self, maeFile, center, pocketId, radius=5):
    closeFile = f"closeResidues_{pocketId}.txt"
    args = f' {maeFile} {radius} "{center}" {closeFile}'
    subprocess.run(f'{runPath} {closeScript} {args}', check=True, capture_output=True, text=True, shell=True,
                   cwd=self._getTmpPath())

    with open(self._getTmpPath(closeFile)) as f:
      f.readline()
      res = f.readline()
    return res

  def getBindingSiteCenter(self, pocket=None):
    if self.fromPockets.get() == AS:
      pocketId = None
      inAS = self.inputAtomStruct.get()
      pdbFile = inAS.getFileName()
      if not pdbFile.endswith('.pdb'):
        pdbFile = inAS.convert2PDB(self._getExtraPath('inputStructure.pdb'))

      _, x, y, z = calculate_centerMass(pdbFile)
    elif self.fromPockets.get() == POCKET:
      pocketId = pocket.getObjId()
      x, y, z = pocket.calculateMassCenter()
    elif self.fromPockets.get() == GRID:
      pocketId = pocket.getObjId()
      x, y, z = pocket.getCenter()

    closeResStr = self.getCloseResidues(self.getInputReceptorFile(), [x, y, z], pocketId)
    return closeResStr

  def getDiameter(self, msjDic, pocket=None):
    if self.fromPockets.get() == AS:
      d = 2 * msjDic['radius']
    elif self.fromPockets.get() == POCKET:
      d = pocket.getDiameter()
    return d

  def dic2StrArgs(self, argsDic, toAdd=''):
    iStr = ''
    for aKey, aVal in argsDic.items():
      iStr += f'  {toAdd}{aKey} {aVal}\n'
    return iStr

  def getScoreName(self, msjDic):
    if msjDic['stageType'] == 'Score poses':
      return msjDic['scoreName']

  def getAllScoreNames(self):
    scoreNames = []
    if self.workFlowSteps.get() in ['', None]:
      msjDic = createMSJDic(self)
      scoreName = self.getScoreName(msjDic)
      if scoreName:
        scoreNames.append(scoreName)

    else:
      workSteps = self.workFlowSteps.get().split('\n')
      if '' in workSteps:
        workSteps.remove('')
      for wStep in workSteps:
        msjDic = eval(wStep)
        scoreName = self.getScoreName(msjDic)
        if scoreName:
          scoreNames.append(scoreName)
    return scoreNames

  def splitComplexFiles(self, oFile, oDir, title, idx):
    complexFile = os.path.join(oDir, f'ifd_{title}_{idx}.maegz')
    args = f' -n {idx} {oFile} {complexFile}'
    subprocess.run(f'{structConvertProg} {args}', check=True, capture_output=True, text=True, shell=True, cwd=oDir)

    args = f' -m pdb -many_files {complexFile} {complexFile}'
    subprocess.run(f'{runPath} {splitProg} {args}', check=True, capture_output=True, text=True, shell=True,
                   cwd=oDir)

    os.remove(complexFile)
    recFile, molFile = complexFile.replace('.maegz', '_receptor1.maegz'), \
                       complexFile.replace('.maegz', '_ligand1.maegz')
    return recFile, molFile



