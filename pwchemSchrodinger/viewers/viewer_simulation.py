# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params

from pwchem.viewers import VmdViewPopen

from ..protocols import ProtSchrodingerDesmondMD
from .. import Plugin as schPlugin
from ..objects import SchrodingerSystem
from ..constants import TCL_MD_STR

from .viewers_data import MaestroView

class DesmondSimulationViewer(pwviewer.ProtocolViewer):
    """ Visualize the output of Desmond simulation """
    _label = 'Viewer Desmond Simulation'
    _targets = [ProtSchrodingerDesmondMD]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        form.addSection(label='Visualization of Desmond Simulation')
        group = form.addGroup('Open Maestro GUI')
        group.addParam('displayMaestro', params.LabelParam,
                      label='Open simulation in Maestro: ',
                      help='Display Simulation in Maestro GUI.\n'
                           'in order to visualize the simulation, you need to load the trajectory by clicking on '
                           'the "T" next to the Maestro object, in the left panel.\n '
                           'https://www.schrodinger.com/kb/485'
                      )
        group.addParam('displayVMD', params.LabelParam,
                       label='Open simulation in VMD: ',
                       help='Display Simulation in VMD GUI.\n')

        group = form.addGroup('Launch Schrodinger Event Analysis')
        group.addParam('displayCompleteAnalysis', params.LabelParam,
                      label='Perform complete analysis in Maestro: ',
                      help='Analyze Simulation in Maestro and open its GUI. \n'
                           'This option will automatically run every analysis available and show the results in a '
                           'Schrodinger GUI.\n'
                           'Reports or images can be exported from this GUI'
                      )
        group.addParam('displayCustomAnalysis', params.LabelParam,
                      label='Perform custom analysis in Maestro: ',
                      help='Analyze Simulation in Maestro and open its GUI. \n'
                           'With this option, you can choose among the analysis available and it will show '
                           'their results in a Schrodinger GUI.\n'
                           'Reports or images can be exported from this GUI'
                      )


    def _getVisualizeDict(self):
        return {
          'displayMaestro': self._showMaestro,
          'displayCompleteAnalysis': self._showCompleteAnalysis,
          'displayCustomAnalysis': self._showCustomAnalysis,
          'displayVMD': self._showMdVMD
        }

    def _showMaestro(self, paramName=None):
        system = self.getMDSystem()
        return [MaestroView(os.path.abspath(system.getFileName()), cwd=self.protocol._getPath())]

    def _showMdVMD(self, paramName=None):
      system = self.getMDSystem()

      outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
      cmsFile = system.getFileName()
      clickMeFile = cmsFile.replace('.cms', '_trj/clickme.dtr')
      with open(outTcl, 'w') as f:
        f.write(TCL_MD_STR % (cmsFile, clickMeFile))
      args = '-e {}'.format(outTcl)

      return [VmdViewPopen(args)]


    def _showCompleteAnalysis(self, paramName=None):
        system = self.getMDSystem()
        #Generates the event analysis file (eaf) with instructions
        eventScript = schPlugin.getMMshareDir('python/scripts/event_analysis.py')
        baseName = system.getBaseName()
        inEAF = self.protocol._getExtraPath('{}-in.eaf'.format(baseName))
        outEAF = '{}_pl_complete.eaf'.format(baseName)

        inSys = os.path.abspath(system.getFileName())
        args = 'analyze {} -o extra/{}'.format(inSys, baseName)
        if not os.path.exists(inEAF):
            schPlugin.runSchrodingerScript(program=eventScript, args=args, cwd=system.getDirName())
            if not os.path.exists(inEAF):
                args = 'analyze {} -lig none -o extra/{}'.format(inSys, baseName)
                schPlugin.runSchrodingerScript(program=eventScript, args=args, cwd=system.getDirName())

        #Analyze the simulations with the instructions provided
        args = '{} {} extra/{} extra/{}-in.eaf'.format(inSys, system.getTrajectoryDirName(),
                                           outEAF, baseName)
        analysisScript = schPlugin.getHome('internal/bin/analyze_simulation.py')
        if not os.path.exists(self.protocol._getExtraPath(outEAF)):
            schPlugin.runSchrodingerScript(program=analysisScript, args=args, cwd=system.getDirName())

        # Run the analysis GUI
        args = 'gui ' + outEAF
        schPlugin.runSchrodingerScript(program=eventScript, args=args, cwd=self.protocol._getExtraPath(),
                                       popen=True)

    def _showCustomAnalysis(self, paramName=None):
        #Run the analysis GUI
        system = self.getMDSystem()
        inSys = os.path.abspath(system.getFileName())
        eventScript = schPlugin.getMMshareDir('python/scripts/event_analysis.py')
        args = 'gui ' + os.path.abspath(inSys)
        schPlugin.runSchrodingerScript(program=eventScript, args=args, cwd=self.protocol._getExtraPath(),
                                       popen=True)


    def getMDSystem(self, objType=SchrodingerSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        else:
            return self.protocol.outputSystem

