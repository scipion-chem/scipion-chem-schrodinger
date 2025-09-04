# **************************************************************************
# *
# * Authors: Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# General imports
import os, csv
from typing import Union, List

# Scipion em imports
from pwem.protocols import EMProtocol
from pwem.objects import String, Float, Integer
from pyworkflow.protocol.params import STEPS_PARALLEL, PointerParam, BooleanParam, FloatParam, IntParam, LEVEL_ADVANCED
from pyworkflow.utils import redStr, Message

# Scipion chem imports
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import removeElements

# Plugin imports
from .. import Plugin

# Output variable name
OUTPUTATTRIBUTE = "outputSmallMolecules"

class ProtSchrodingerQikprop(EMProtocol):
	""" Qikprop analyzes the properties of a given set of small molecules. 

	User Documentation(AI GENERATED)
	The ProtSchrodingerQikprop class is designed to facilitate the analysis of small molecules using Schrodinger's Qikprop tool. Qikprop is used to predict various molecular properties such as lipophilicity, 
	solubility, molecular weight, and other essential characteristics of molecules. The protocol allows users to input a set of small molecules and generate predictions about these properties. 
	Additionally, it can also generate a list of similar known drugs based on the properties of each molecule in the input set.
	This class offers several parameters that control how the Qikprop analysis is performed. Users can choose to run Qikprop in a faster mode by excluding certain calculations like dipole, 
	HOMO, and LUMO, or they can opt for a more comprehensive analysis. The fast parameter enables the fast mode, while the sim parameter allows the generation of a list of known drugs that 
	are similar to each processed molecule. Users can also specify how many similar drugs they want to report with the nsim parameter. The class includes advanced options, such as neutralizing 
	molecules before processing or using an alternative probe radius for SASA and PSA calculations.

	The execution of the protocol is parallelized, enabling efficient processing of multiple molecules at once. The process involves several steps: first, it runs Qikprop for each molecule in the 
	input set, applying the appropriate parameters. Then, if the cleanTmps parameter is selected, it cleans up any temporary files generated during execution. Finally, the results are processed, 
	and the properties calculated by Qikprop are added to each molecule in the output set. This output set contains the molecules with their newly calculated properties, ready for further analysis or use.
	Before running the protocol, the system performs several validation checks to ensure the correct parameters are set. For example, it checks that MPI is not selected (as only threads are supported), 
	ensures that if an alternative probe is used, its value is valid, and verifies that the input set contains at least one molecule. These checks prevent errors during execution and ensure that the protocol 
	runs smoothly.

	Once the Qikprop analysis is completed, the calculated properties are saved in a CSV file. This file is parsed, and the properties are added as attributes to each molecule in the output set. 
	The addCSVProperties method reads the CSV file, assigns the corresponding values to the molecules, and ensures that the properties are correctly stored. This detailed process allows users to analyze and 
	compare the properties of their small molecules with ease.
	
	The class also includes utility functions that assist in the execution of the protocol. The getQikpropBaseCmd method constructs the command string needed to run Qikprop, incorporating various 
	flags based on user input. The cleanTmpFiles function handles the deletion of temporary files, and the getCSVTextValue method ensures that the values read from the CSV file are correctly interpreted as 
	integers, floats, or strings. These utilities help streamline the execution of the protocol and improve its efficiency.
	In summary, the ProtSchrodingerQikprop class provides a comprehensive solution for analyzing the properties of small molecules. By running Qikprop and processing the results, users can gain valuable 
	insights into the molecular characteristics of their compounds. With customizable parameters, parallel execution, and built-in validation and cleanup, this protocol is a powerful tool for researchers 
	in drug discovery and molecular design.
	"""

	_label = 'qikprop'
	_possibleOutputs = {OUTPUTATTRIBUTE: SetOfSmallMolecules}

	# --------------------------- Class constructor --------------------------------------------
	def __init__(self, **args):
		# Calling parent class constructor
		super().__init__(**args)

		# Defining execution mode. Steps will take place in parallel now
		# Full tutorial on how to parallelize protocols can be read here:
		# https://scipion-em.github.io/docs/release-3.0.0/docs/developer/parallelization.html
		self.stepsExecutionMode = STEPS_PARALLEL

	def _defineParams(self, form):
		""" This function defines the params to be included in the protocol's form. """
		# Add parallel section
		form.addParallelSection(threads=4)

		# Create form
		form.addSection(label=Message.LABEL_INPUT)

		# Main params
		form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
			label='Input small molecules:', help='Input small molecules to be treated.')
		form.addParam('fast', BooleanParam, default=False, label='Run in fast mode:',
			help='Run job in fast mode (no dipole, homo or lumo).')
		form.addParam('sim', BooleanParam, default=True, label='Generate similar known drugs:',
			help='Generate a list of known drugs most similar to each processed molecule.')
		form.addParam('nsim', IntParam, default=5, label='Number of most similar drugs:', condition="sim==True",
			help="Number of most similar drug molecules to report.")
		
		# Additional optional params
		form.addParam('neut', BooleanParam, default=True, label='Neutralize molecules:', expertLevel=LEVEL_ADVANCED,
			help='Neutralize molecules in Maestro formatted files prior to processing.')
		form.addParam('altclass', BooleanParam, default=False, label='Alternative class:', expertLevel=LEVEL_ADVANCED,
			help='Run additional SASA, PSA calculations using an alternative class definition.')
		form.addParam('useAltprobe', BooleanParam, default=False, label='Use alternative probe:', expertLevel=LEVEL_ADVANCED,
			help='Run additional SASA, PSA calculations using an alternative probe radii (1.4A by default).')
		form.addParam('altprobe', FloatParam, default=1.4, label='Alternative probe:', expertLevel=LEVEL_ADVANCED, condition='useAltprobe==True')
		form.addParam('recap', BooleanParam, default=False, label='Replace CombiGlide with methyl:', expertLevel=LEVEL_ADVANCED,
			help='Replace the CombiGlide functional group with a methyl group prior to processing.\n'
				'When used in the CombiGlide reagent-preparation process, gives properties for the \'naked sidechain\'.')
		
		# Utils params
		form.addParam('cleanTmps', BooleanParam, default=True, label='Clean temporary files: ', expertLevel=LEVEL_ADVANCED,
			help='Clean temporary files after finishing the execution.\nThis is useful to reduce unnecessary disk usage.')

	# --------------------------- INSERT steps functions --------------------
	def _insertAllSteps(self):
		""" This function inserts all steps functions that will be run when running the protocol. """
		# Getting common command string for every call
		baseCommand = self.getQikpropBaseCmd()

		# For every SmallMolecule in the input set, build complete command
		deps = []
		for molecule in self.getInputFiles():
			deps.append(self._insertFunctionStep('runQikpropStep', baseCommand, molecule))
		
		# Clean tmp files if selected
		if self.cleanTmps.get():
			self._insertFunctionStep('cleanTmpFiles', prerequisites=deps)
		
		# Create output
		self._insertFunctionStep('createOutputStep', prerequisites=deps)

	def runQikpropStep(self, baseCommand, molecule):
		""" This function runs the schrodinger binary file with the given params. """
		self.runJob(baseCommand, f' {molecule}', cwd=self._getExtraPath())
	
	def cleanTmpFiles(self):
		""" This function removes the temporary files related to the execution of qikprop for all the molecules. """
		# Creating empty list to store all tmp files for all molecules
		tmpFileList = []

		# Generating tmp files for each molecule and adding them to the list
		for molecule in self.getInputFiles():
			moleculeBasePath = os.path.join(os.path.abspath(self._getExtraPath()), os.path.splitext(os.path.basename(molecule))[0])
			tmpFileList.append(moleculeBasePath + '.log')
			tmpFileList.append(moleculeBasePath + '.out')
			tmpFileList.append(moleculeBasePath + '-out.sdf')
			tmpFileList.append(moleculeBasePath + '.qpsa')
		
		# Deleting all tmp files
		removeElements(tmpFileList)
	
	def createOutputStep(self):
		""" This function generates the output of the protocol. """
		# Obtaining input small molecules
		inputMolecules = self.inputSmallMolecules.get()

		# Creating output small molecule set
		outputSmallMolecules = inputMolecules.createCopy(self._getPath(), copyInfo=True)

		# Add analyzed properties for each molecule
		for molecule in inputMolecules:
			outMol = self.addCSVProperties(molecule)
			outputSmallMolecules.append(outMol)
		
		# Generate output
		self._defineOutputs(**{OUTPUTATTRIBUTE: outputSmallMolecules})

	# --------------------------- INFO functions --------------------------------------------
	def _validate(self):
		"""
		The function of this hook is to add some validation before the protocol
		is launched to be executed. It should return a list of errors. If the list is
		empty the protocol can be executed.
		"""
		errors = []

		# Checking if MPI is selected (only threads are allowed)
		if self.numberOfMpi > 1:
			errors.append('MPI cannot be selected, because Scipion is going to drop support for it. Select threads instead.')
		
		# Checking if altprobe is used, and, if so, that the value is 0 or greater
		if self.useAltprobe.get() and self.altprobe.get() < 0:
			errors.append('Altprobe has to be a number equal or greater than 0.')
		
		# Checking if the number of elements of the input set is at least 1
		if len(self.inputSmallMolecules.get()) < 1:
			errors.append('The input set of molecules must at least contain 1 molecule.')

		return errors
	
	def _summary(self):
		"""
		This method dumps the contents of the .warning files produced by Qikprop into the summary text.
		"""
		summary = []

		# Getting a list of all .warning files
		warningFiles = []
		for extraFile in os.listdir(os.path.abspath(self._getExtraPath())):
			if extraFile.endswith(".warning"):
				warningFiles.append(os.path.abspath(self._getExtraPath(extraFile)))

		# Read each .warning file and dump the contents into summary
		for warningFile in warningFiles:
			with open(warningFile) as wf:
				summary.append(wf.read())

		return summary
	
	# --------------------------- Utils functions --------------------
	def getQikpropBaseCmd(self) -> str:
		""" This function returns the command string to run qikprop. """
		# Command starts with the executable file
		command = self.getQikpropBinaryFile()

		# Add permanent flags
		command += f' {self.getFastFlag()} {self.getSimFlag()} {self.getNeutFlag()} -LOCAL'

		# Add optional flags
		command += self.getNSim() + self.getAltClassFlag() + self.getAltProbeFlag() + self.getRecapFlag()

		# Return formatted command string
		return command

	def getQikpropBinaryFile(self) -> str:
		""" This function returns the location for the Schrodinger qikprop binary file. """
		# Getting path to the binary
		binaryPath = os.path.join(Plugin.getVar('SCHRODINGER_HOME'), 'qikprop')

		# If path exists, return it
		if os.path.exists(binaryPath):
			return binaryPath
		
		# If path was not found, raise exception
		raise FileNotFoundError(redStr(f"Path \"{binaryPath}\" not found. Is variable SCHRODINGER_HOME properly set within scipion.conf file?"))
	
	def getInputFiles(self) -> List[str]:
		""" This function returns a list with the full path to each one of the input files. """
		return [os.path.abspath(molecule.getFileName()) for molecule in self.inputSmallMolecules.get()]
	
	def getFastFlag(self) -> str:
		""" This function returns the flag string corresponding to the fast flag param. """
		return '-fast' if self.fast.get() else '-nofast'

	def getSimFlag(self) -> str:
		""" This function returns the flag string corresponding to the sim flag param. """
		return '-sim' if self.sim.get() else '-nosim'
	
	def getNSim(self) -> str:
		""" This function returns the flag string corresponding to the nsim flag param. """
		return f' -nsim {self.nsim.get()}' if self.sim.get() else ''
	
	def getNeutFlag(self) -> str:
		""" This function returns the flag string corresponding to the neut flag param. """
		return '-neut' if self.neut.get() else '-noneut'
	
	def getAltClassFlag(self) -> str:
		""" This function returns the flag string corresponding to the altclass flag param. """
		return ' -altclass' if self.altclass.get() else ''
	
	def getAltProbeFlag(self) -> str:
		""" This function returns the flag string corresponding to the altprobe flag param. """
		return f' -altprobe {self.altprobe.get()}' if self.useAltprobe.get() else ''
	
	def getRecapFlag(self) -> str:
		""" This function returns the flag string corresponding to the sim recap param. """
		return ' -recap' if self.recap.get() else ''

	def addCSVProperties(self, molecule : SmallMolecule) -> SmallMolecule:
		""" This function adds the properties dumped by Qikprop in a CSV file for a given molecule. """
		# Creating new molecule that will be a clone of the input one, with extra attributes
		outputMolecule = molecule.clone()

		# Getting CSV file path
		csvFile = os.path.splitext(os.path.abspath(self._getExtraPath(os.path.basename(molecule.getFileName()))))[0] + '.CSV'
		if os.path.isfile(csvFile):
			with open(csvFile, mode='r') as attrFile:
				rows = list(csv.reader(attrFile))

				# If number of headers does not match number of row columns, values will be empty
				# thus, not storing them. Also discard csv if it contains less than 2 columns (1st is id)
				if len(rows[0]) != len(rows[1]) or len(rows[1]) < 2:
					return outputMolecule
				
				# Setting info into output molecule
				for header, value in zip(rows[0], rows[1]):
					header = header.replace('.', '_') # Dots must be replaced with underscore to avoid errors
					value = self.getCSVTextValue(value)
					if header != 'molecule':
						setattr(outputMolecule, header, value)
			
		return outputMolecule

	def getCSVTextValue(self, text : str) -> Union[Integer, Float, String, None]:
		"""
		This function returns the value of the given text in the appropiate data type.
		Supported values are int, float, and str.
		"""
		# If string is empty, return None
		if not text:
			return None
		
		# If all characters are a number, it must be an integer
		if text.isnumeric():
			return Integer(text)
		
		# Try to convert to float.
		# If possible, it is float
		# If not, it is an str
		try:
			return Float(text)
		except ValueError:
			return String(text)
