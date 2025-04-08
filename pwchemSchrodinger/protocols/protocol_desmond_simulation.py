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
import random as rd
import glob, os, re
from subprocess import check_call

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, FloatParam, BooleanParam, StringParam, TextParam, LEVEL_ADVANCED, Group, Line

# Scipion chem imports
from pwchem.utils import natural_sort, createMSJDic

# Plugin imports
from .. import Plugin as schrodinger_plugin
from ..constants import MSJ_SYSMD_INIT
from ..objects import SchrodingerSystem
from ..utils import buildSimulateStr, getJobName, setAborted, getSchJobId

multisimProg = schrodinger_plugin.getHome('utilities/multisim')
jobControlProg = schrodinger_plugin.getHome('jobcontrol')
mergeScript = schrodinger_plugin.getHome('internal/bin/trj_merge.py')

class ProtSchrodingerDesmondMD(EMProtocol):
    """Calls Desmond molecular dynamics for the preparation of the system via solvatation, the addition of ions
    and a force field"""
    _label = 'system molecular dynamics (desmond)'

    NONE, DESMOND_NPT = 0, 1

    _ensemTypes = ['NVE', 'NVT', 'NPT', 'NPAT', 'NPgT', 'Minimization (Brownian)']
    _thermoDic = {'Noose-Hover': 'NH', 'Langevin': 'Langevin', 'None': 'None'}
    _thermostats = ['Noose-Hover', 'Langevin', 'Brownie', 'None']

    _baroDic = {'Martyna-Tobias-Klein': 'MTK', 'Langevin': 'Langevin', 'None': 'None'}
    _barostats = ['Martyna-Tobias-Klein', 'Langevin', 'None']
    _coupleStyle = ['Isotropic', 'Semi-isotropic', 'Anisotropic', 'Constant area']
    _restrainTypes = ['None', 'Ligand', 'Protein', 'Solute_heavy_atom', 'Solute']

    _omitParamNames = ['runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                       'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi',
                       'inputStruct', 'defSteps']

    _defParams = {'simTime': 100, 'bondedT': 0.002, 'nearT': 0.002, 'farT': 0.006,
                  'annealTemps': '[300, 0]', 'annealing': False,
                  'velResamp': 1.0, 'glueSolute': True, 'trajInterval': 5.0,
                  'temperature': 300.0, 'deltaMax': 0.1, 'tempMDCons': 0.1, 'pressure': 1.01325,
                  'presMDCons': 2.0, 'surfTension': 0.0, 'restrainForce': 50.0,
                  'ensemType': 'NVE', 'thermostat': 'Noose-Hover', 'barostat': 'Martyna-Tobias-Klein',
                  'coupleStyle': 'Isotropic', 'restrains': 'None'}


    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStruct', PointerParam, pointerClass='SchrodingerSystem',
                      label='Input system to be prepared for MD:', allowsNull=False,
                      help='Atomic structure to be simulated')
        form.addParam('defSteps', EnumParam, label='Default MD workflows: ',
                      choices=['None', 'NPT desmond relaxation'], default=self.NONE,
                      help='Choose and add with the wizard some default MD steps which will be reflected into '
                           'the summary')

        form.addSection('Molecular dynamics simulation')
        group = form.addGroup('Simulation time')
        group.addParam('simTime', FloatParam, default=self._defParams['simTime'],
                       label='Simulation time (ps):',
                       help='Time of the simulation step (ps)')
        line = group.addLine('Time steps: ',
                             help='Time steps for the simulation (bonded / near / far) (ps)')
        line.addParam('bondedT', FloatParam, default=self._defParams['bondedT'],
                      label='Bonded: ')
        line.addParam('nearT', FloatParam, default=self._defParams['nearT'],
                      label='Near: ')
        line.addParam('farT', FloatParam, default=self._defParams['farT'],
                      label='Far: ')

        group.addParam('velResamp', FloatParam, default=self._defParams['velResamp'],
                       label='Velocity resampling (ps):', expertLevel=LEVEL_ADVANCED,
                       help='Velocities of the particles are resampled every x ps')

        group = form.addGroup('Trajectory')
        group.addParam('glueSolute', BooleanParam, default=self._defParams['glueSolute'],
                      label='Glue close solute molecules: ', expertLevel=LEVEL_ADVANCED,
                      help='Glue close solute molecules together. Does not affect energy, only for visualization'
                           ' of the trajectory')
        group.addParam('trajInterval', FloatParam, default=self._defParams['trajInterval'],
                       label='Interval time (ps):',
                       help='Time between each frame recorded in the simulation (ps)')

        group = form.addGroup('Ensemble')
        group.addParam('ensemType', EnumParam, default=0,
                       label='Simulation type: ', choices=self._ensemTypes,
                       help='Ensemble type of the simulation')

        line = group.addLine('Temperature settings: ', condition='ensemType!=0',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Simulation time constant for thermostat (ps)')
        line.addParam('annealing', BooleanParam, default=self._defParams['annealing'],
                      label='Annealing: ')
        line.addParam('temperature', FloatParam, default=self._defParams['temperature'], condition='not annealing',
                       label='Temperature: ')
        line.addParam('thermostat', EnumParam, default=0, condition='ensemType!=5',
                       label='Thermostat type: ', choices=self._thermostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('deltaMax', FloatParam, default=self._defParams['deltaMax'], condition='ensemType==5',
                      label='Max displacement: ', expertLevel=LEVEL_ADVANCED)
        line.addParam('tempMDCons', FloatParam, default=self._defParams['tempMDCons'],
                       label='Simulation temperature constant: ', expertLevel=LEVEL_ADVANCED)
        group.addParam('annealTemps', StringParam, default=self._defParams['annealTemps'], condition='annealing',
                      label='Annealing pairs [temperature, Starting Time]: ',
                      help='Temperatures and starting time (ps) of each of the annealing steps. '
                           'Temperature between steps is interpolated. '
                           'Format: pairs of temperature duration comma separated [300, 10], [280, 5]')

        line = group.addLine('Pressure settings: ', condition='ensemType not in [0, 1, 5]',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Simulation time constant for barostat (ps)')
        line.addParam('pressure', FloatParam, default=self._defParams['pressure'],
                       label='   Pressure:   ')
        line.addParam('barostat', EnumParam, default=0,
                       label='  Barostat type:   ', choices=self._barostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('presMDCons', FloatParam, default=self._defParams['presMDCons'],
                       label='Simulation pressure constant:   ', expertLevel=LEVEL_ADVANCED)
        line = group.addLine('Pressure coupling: ', condition='ensemType not in [0, 1, 5]',
                             help='Pressure coupling style', expertLevel=LEVEL_ADVANCED)
        line.addParam('coupleStyle', EnumParam, default=0,
                      label='  Coupling style:   ', choices=self._coupleStyle)

        line = group.addLine('Tension settings: ', condition='ensemType==4',
                             help='Surface tension during the simulation (bar·Å)')
        line.addParam('surfTension', FloatParam, default=self._defParams['surfTension'],
                      label='Surface tension: ')

        group = form.addGroup('Restrains')
        group.addParam('restrains', EnumParam, default=0,
                       label='Restrains: ', choices=self._restrainTypes,
                       help='Restrain movement of specific groups of atoms')
        group.addParam('restrainForce', FloatParam, default=self._defParams['restrainForce'],
                       label='Restrain force constant: ', condition='restrains!=0',
                       help='Restrain force applied to the selection (kcal/mol/Å2)')

        group = form.addGroup('Summary')
        group.addParam('insertStep', StringParam, default='',
                      label='Insert simulation step number: ',
                      help='Insert the defined simulation step into the workflow on the defined position.\n'
                           'The default (when empty) is the last position')
        group.addParam('summarySteps', TextParam, width=120, readOnly=True,
                      label='Summary of steps', help='Summary of the defined steps. \nManual modification will have no '
                                                     'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', StringParam, default='',
                       label='Delete simulation step number: ',
                       help='Delete the step of the specified index from the workflow.')
        group.addParam('watchStep', StringParam, default='',
                       label='Watch simulation step number: ',
                       help='''Watch the parameters step of the specified index from the workflow..\n
                       This might be useful if you want to change some parameters of a predefined step.\n
                       However, the parameters are not changed until you add the new step (and probably\n
                       you may want to delete the previous unchanged step)''')
        group.addParam('workFlowSteps', StringParam, label='User transparent', condition='False')

        #form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('simulationStep')
        self._insertFunctionStep('createOutputStep')

    def simulationStep(self):

        lastCheckFile = self.findLastCheckPoint()
        unmergedFile = 'unmerged-simulation.cms'
        self.jobName = self.getCurrentJobName()
        msjFile = self._getTmpPath('{}.msj'.format(self.jobName))

        if not lastCheckFile:
            maeFile = self.inputStruct.get().getFileName()

            msjStr = self.buildMSJStr()
            with open(msjFile, 'w') as f:
                f.write(msjStr)
            self.createGUISummary()

            args = ' -m {} {} -o {} -WAIT -JOBNAME {}'.format(os.path.abspath(msjFile), os.path.abspath(maeFile),
                                                              unmergedFile, self.jobName)

        else:
            lastTGZ = self.findLastDataTgz()
            args = ' -RESTART {} -d {} -o {} -WAIT -JOBNAME {}'.\
              format(os.path.abspath(lastCheckFile), os.path.abspath(lastTGZ),
                     unmergedFile, self.jobName)

        self.runJob(multisimProg, args, cwd=self._getTmpPath())


    def createOutputStep(self):
        #Managing trajectories and names
        maeFile = self.inputStruct.get().getFileName()
        sysName = maeFile.split('/')[-1].split('.')[0]
        unmergedFile = 'unmerged-simulation.cms'
        trjDirs = self.getTrjDirs()

        outBase = sysName + '-MD'
        outFile = os.path.abspath(self._getPath(outBase + '.cms'))
        outTrjDir = os.path.abspath(self._getPath(outBase + '_trj'))
        if len(trjDirs) > 1:
            mergedFile = self._getPath(sysName + '-out.cms')
            args = ' {} {} -o {} -concat 0 {}'.format(unmergedFile, ' '.join(trjDirs),
                                                      os.path.abspath(self._getPath(sysName)), self.trajInterval.get())
            schrodinger_plugin.runJobSchrodingerScript(self, mergeScript, args, cwd=os.path.abspath(self._getTmpPath()))

            cmsStruct = SchrodingerSystem()
            cmsStruct.setFileName(mergedFile)
            cmsStruct.changeTrajectoryDirName(outTrjDir, trjPath=self._getPath())
            cmsStruct.changeCMSFileName(os.path.relpath(outFile))

        else:
            cmsStruct = SchrodingerSystem()
            cmsStruct.setFileName(self._getTmpPath(unmergedFile))
            cmsStruct.changeTrajectoryDirName(outTrjDir, trjPath=self._getTmpPath())
            cmsStruct.changeCMSFileName(os.path.relpath(outFile))

        os.rename(self._getTmpPath('{}_multisim.log'.format(getJobName(self))),
                  self._getExtraPath(sysName+'_multisim.log'))
        os.rename(self._getTmpPath('{}.msj'.format(getJobName(self))), self._getPath('simulation.msj'))
        self._defineOutputs(outputSystem=cmsStruct)

#######################################

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

    def _validate(self):
        errors = []
        if not self.workFlowSteps.get():
            errors += ['You need to define some simulation. Try the default MD workflows in doubt.']
        else:
            workSteps = self.workFlowSteps.get().split('\n')
            if '' in workSteps:
                workSteps.remove('')
            for wStep in workSteps:
                msjDic = eval(wStep)
                msjDic = self.addDefaultForMissing(msjDic)
                errors += self.validateAnneal(msjDic)

        return errors

    def validateAnneal(self, msjDic):
        errors = []
        if msjDic['annealing']:
            annealTemps = msjDic['annealTemps'].replace('(', '[').replace(']', ']')
            annealStages = re.findall(r'\[\d+[.]?\d*, \d+\]', annealTemps)
            for aSt in annealStages:
                if eval(aSt)[1] > msjDic['simTime']:
                    errors.append('Annealing starting time {} should not be larger than total simulation time {}\n'.
                                  format(eval(aSt)[1], msjDic['simTime']))
        return errors

    def getStageParamsDic(self, type='All'):
      '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
      Type'''
      paramsDic = {}
      for paramName, param in self._definition.iterAllParams():
        if paramName not in self._omitParamNames and not isinstance(param, Group) and not isinstance(param, Line):
          if type == 'All' or (type == 'Enum' and isinstance(param, EnumParam)) or (type == 'Normal' and not isinstance(param, EnumParam)):
            paramsDic[paramName] = param
      return paramsDic

    ############# UTILS
    def addDefaultForMissing(self, msjDic):
      '''Add default values for missing parameters in the msjDic'''
      paramDic = self.getStageParamsDic()
      for pName in paramDic.keys():
        if pName not in msjDic:
          msjDic[pName] = self._defParams[pName]
      return msjDic

    def buildMSJStr(self):
        '''Build the .msj (file used by multisim to specify the jobs performed by Schrodinger)
        defining the input parameters'''
        msjStr = MSJ_SYSMD_INIT

        if self.workFlowSteps.get() in ['', None]:
            msjDic = createMSJDic(self)
            msjStr += buildSimulateStr(self, msjDic)
        else:
            workSteps = self.workFlowSteps.get().split('\n')
            if '' in workSteps:
                workSteps.remove('')
            for wStep in workSteps:
                msjDic = eval(wStep)
                msjStr += buildSimulateStr(self, msjDic)

        return msjStr

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
                    method, ensemType = self.getMethodEnsemType(msjDic)
                    sumStr += '{}) Sim. time: {} ps, {} ensemble with {} method'.\
                      format(i+1, msjDic['simTime'], ensemType, method)

                    if msjDic['restrains'] != 'None':
                        sumStr += ', restrain on {}'.format(msjDic['restrains'])
                    if msjDic['annealing']:
                        sumStr += ', annealing on'
                    else:
                        sumStr += ', {} K'.format(msjDic['temperature'])
                    sumStr += '\n'
        else:
            msjDic = self.addDefaultForMissing(msjDic)
            method, ensemType = self.getMethodEnsemType(msjDic)
            sumStr += 'Sim. time: {} ps, {} ensemble with {} method'. \
              format(msjDic['simTime'], ensemType, method)

            if msjDic['restrains'] != 'None':
              sumStr += ', restrain on {}'.format(msjDic['restrains'])
            if msjDic['annealing']:
              sumStr += ', annealing on'
            else:
              sumStr += ', {} K'.format(msjDic['temperature'])
            sumStr += '\n'
        return sumStr

    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def getMethodEnsemType(self, msjDic):
      method = self._thermoDic[msjDic['thermostat']]
      ensemType = msjDic['ensemType']
      if ensemType not in ['NVE', 'NVT', 'Minimization (Brownian)']:
        method = self._baroDic[msjDic['barostat']]

      if ensemType == 'Minimization (Brownian)':
        ensemType = 'NVT'
        method = 'Brownie'
      return method, ensemType


    def getTrjDirs(self):
        import re
        def atoi(text):
          return int(text) if text.isdigit() else text

        def naturalKeys(text):
          '''
          alist.sort(key=naturalKeys) sorts in human order
          http://nedbatchelder.com/blog/200712/human_sorting.html
          (See Toothy's implementation in the comments)
          '''
          return [atoi(c) for c in re.split(r'(\d+)', text)]

        files = os.listdir(self._getTmpPath())
        trjDirs = []
        for file in files:
            if file.endswith('_trj'):
                trjDirs.append(file)


        trjDirs.sort(key=naturalKeys)
        return trjDirs

    def parseAnnealing(self, annealTemps):
        tempArg = ' [\n'
        annealTemps = annealTemps.replace('(', '[').replace(']', ']')
        annealStages = re.findall(r'\[\d+[.]?\d*, \d+\]', annealTemps)
        for aSt in annealStages:
            tempArg += '    [{} {}]\n'.format(eval(aSt)[0], eval(aSt)[1])
        return tempArg + ']'
    
    def setAborted(self):
        super().setAborted()
        setAborted(getSchJobId(self), getJobName(self))

    def findLastCheckPoint(self):
        cpFiles = glob.glob(self._getTmpPath('*_checkpoint'))
        if cpFiles:
            return natural_sort(cpFiles)[-1]
        return cpFiles

    def findLastDataTgz(self):
        dataFiles = glob.glob(self._getTmpPath('*-out.tgz'))
        if dataFiles:
            return natural_sort(dataFiles)[-1]
        return dataFiles

    def getCurrentJobName(self):
        msjFile = glob.glob(self._getTmpPath('simulation*.msj'))
        if msjFile:
            return os.path.basename(msjFile[-1]).replace('.msj', '')
        else:
            return 'simulation_' + str(rd.randint(1000000, 9999999))
