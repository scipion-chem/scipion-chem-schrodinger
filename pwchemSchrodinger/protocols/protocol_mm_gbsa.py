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
import os, time, csv

# Scipion em imports
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb
from pyworkflow.protocol import params
import pyworkflow.object as pwobj

# Scipion chem imports
from pwchem.utils import pdbqt2other, convertToSdf, performBatchThreading, getBaseName
from pwchem.objects import SmallMolecule

# Plugin imports
from .. import Plugin as schrodingerPlugin
from ..utils import isMaeFile

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
structCatProg = schrodingerPlugin.getHome('utilities/structcat')
progLigPrep = schrodingerPlugin.getHome('ligprep')
progPrepWizard = schrodingerPlugin.getHome('utilities/prepwizard')
mmgbsaProg = schrodingerPlugin.getHome('prime_mmgbsa')

FNONE, FTHRES, FAUT = 0, 1, 2
JOB_TYPES = ['REAL_MIN', 'ENERGY', 'SIDE_PRED', 'SIDE_COMBI', 'SITE_OPT']
REAL_MIN, ENERGY = 0, 1
FORCE_FIELDS = ['S-OPLS', 'OPLS_2005']
MAEFILE_EXTENSION = '.maegz'
LICENSE_ERROR = 'Could not checkout licenses'
LICENSE_WAIT_MIN, LICENSE_WAIT_MAX = 60, 600  # seconds between license retries
BIND_COLUMN = 'r_psp_MMGBSA_dG_Bind'
ENERGIES_FILE = 'mmgbsa_energies.csv'

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
                       help='Prepare ligands with Schrodinger LigPrep to ensure a correct structure in case they were '
                            'not originally Schrodinger molecules. It may subtly vary the ligand atom positions. '
                            'Only the first LigPrep output structure of each pose is kept. '
                            'Maestro poses (e.g. from Glide) are always used as they are.')
        group.addParam('prepareReceptor', params.BooleanParam, default=True, label='Prepare target: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Prepare the receptor with Schrodinger PrepWizard in case it was not originally a '
                            'Schrodinger receptor: fix the structure, add hydrogens, assign protonation states at pH 7 '
                            '(PROPKA) and optimize the H-bond network (ProtAssign), then minimize the hydrogens only, '
                            'so heavy atoms stay where the ligands were docked.\n'
                            'Maestro receptors (e.g. from Glide) are always used as they are.')

        group = form.addGroup('Calculation')
        group.addParam('jobType', params.EnumParam, default=REAL_MIN, label='Sampling: ',
                       choices=['Local minimization', 'Single-point energy', 'Side-chain prediction',
                                'Combinatorial side-chain prediction', 'Binding-site optimization'],
                       help='How Prime treats the flexible parts (ligand and flexible receptor residues) '
                            'before computing the energies (prime_mmgbsa -job_type).\n\n'
                            '*Local minimization* (default)\n'
                            '  Minimizes the ligand and the flexible residues around the docked pose.\n'
                            '  Time: ~1 min per 2 poses.\n\n'
                            '*Single-point energy*\n'
                            '  No minimization: scores the docked pose exactly as it is.\n'
                            '  Fastest; good to rescore large libraries. Small clashes are penalized.\n'
                            '  Time: ~1 min per 2 poses (mostly setup).\n\n'
                            '*Side-chain prediction*\n'
                            '  Repacks the flexible side chains, then minimizes.\n'
                            '  Time: ~2 min per 2 poses.\n\n'
                            '*Combinatorial side-chain prediction*\n'
                            '  Exhaustive version of side-chain prediction.\n'
                            '  Only for fewer than 5 flexible side chains.\n\n'
                            '*Binding-site optimization*\n'
                            '  Side-chain predictions plus minimizations, to model induced-fit effects.\n'
                            '  Most expensive. Time: ~5 min per 2 poses.\n\n'
                            'Side-chain prediction, combinatorial side-chain prediction and binding-site '
                            'optimization need flexible receptor residues: set "Flexible residues" to an option '
                            'other than None.\n'
                            'Times measured on 1 core.')
        group.addParam('forceField', params.EnumParam, default=0, label='Force field: ',
                       choices=FORCE_FIELDS,
                       help='Force field used by Prime (-prime_opt OPLS_VERSION).\n'
                            'S-OPLS: Prime default, the newest OPLS force field (needs its license).\n'
                            'OPLS_2005: older public force field; use it to compare with older results or if '
                            'S-OPLS is not licensed.')

        group = form.addGroup('Target')
        group.addParam('flexOption', params.EnumParam, default=0, label='Flexible residues: ',
                       choices=['None', 'Threshold distance', 'Flexibility estimation'],
                       help='Select the way to determine the target flexible residues.\n'
                            'None: target residues will remain rigid\n'
                            'Threshold distance: flexible residues will be determined by a distance to the ligand\n'
                            'Flexibility estimation: run a two-stage MMGBSA calculation where the second stage runs '
                            'with the subset of flexible protein residues identified by the first')
        group.addParam('flexDist', params.FloatParam, default=3.0, label='Distance to ligand (A): ',
                       condition='flexOption==1',
                       help='Treat all residues within this distance of the ligand as flexible. '
                            'By default the entire receptor is frozen')
        group.addParam('flexGroup', params.EnumParam, default=0, label='Flexible group: ',
                       choices=['Residue', 'Side', 'PolarH'], condition='flexOption==1',
                       help='Select a portion of the region defined with "flexible threshold" to be flexible.'
                            '\nResidue: Choose the entire residue.'
                            '\nSide: Choose the sidechain of each residue. '
                            '\nPolarH: Choose the polar hydrogens on each residue.')
        group.addParam('flexTargetCutOff', params.FloatParam, default=3.0, label='Movement cutoff (A): ',
                       condition='flexOption==2',
                       help='Cutoff for determining movement for target flexibility in Angstroms')

        group = form.addGroup('Ligand', expertLevel=params.LEVEL_ADVANCED,)
        group.addParam('useCharges', params.BooleanParam, default=False, label='Use ligand charges: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Use the partial charges in the input ligand file rather than from the force field')
        group.addParam('rigidLigand', params.BooleanParam, default=False, label='Rigid ligand: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Minimize the ligand as a rigid body')

        form.addParallelSection(threads=4, mpi=1,
                                binThreadsHelp='Number of prime_mmgbsa jobs run at the same time (1 core each). '
                                            'Each job uses 8 Prime (PSP_PLOP) license tokens for its whole batch '
                                            'of poses. Jobs that cannot get licenses wait and retry until '
                                            'the tokens are freed.')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        """ This function inserts all steps functions that will be run when running the protocol """
        convStep = self._insertFunctionStep('convertStep', prerequisites=[])
        mmStep = self._insertFunctionStep('runMMGBSAStep', prerequisites=[convStep])
        self._insertFunctionStep('createOutputStep', prerequisites=[mmStep])

    def convertStep(self):
        '''Builds nt complex files (receptor + several ligand poses) and, next to each, the ids of its molecules'''
        nt = self.numberOfThreads.get()
        inMols = self.inputSmallMolecules.get()
        recMaeFile = self.prepareReceptorFile(inMols)
        molMaes = performBatchThreading(self.prepareLigandFile, inMols, nt)

        nSubMols, lastMol = self.getNMolsPerThread(len(molMaes), nt), 0
        for i, nSubMol in enumerate(nSubMols):
            batch = molMaes[lastMol:lastMol + nSubMol]
            lastMol += nSubMol
            complexFile = self.getComplexFile(i)
            args = f'-imae {recMaeFile} -imae {" ".join(maeFile for _, maeFile in batch)} -omae {complexFile}'
            self.runJob(structCatProg, args, cwd=self._getTmpPath())
            with open(complexFile.replace(MAEFILE_EXTENSION, '.ids'), 'w') as f:
                f.write('\n'.join(str(molId) for molId, _ in batch))

    def runMMGBSAStep(self):
        """ This function runs the schrodinger binary file with the given params """
        nt = self.numberOfThreads.get()
        performBatchThreading(self.performMMGBSA, self.getComplexFiles(), nt, cloneItem=False)

    def createOutputStep(self):
      inMols = self.inputSmallMolecules.get()
      outputSet = inMols.createCopy(self._getPath(), copyInfo=True)
      energyRows = []
      for complexFile in self.getComplexFiles():
        name = getBaseName(complexFile)
        outDir = self._getExtraPath(name)
        outMae, csvFile = os.path.join(outDir, f'{name}-out.maegz'), os.path.join(outDir, f'{name}-out.csv')
        with open(complexFile.replace(MAEFILE_EXTENSION, '.ids')) as f:
          molIds = [int(line) for line in f if line.strip()]
        results = self.parseOutputCSVFile(csvFile) if os.path.exists(csvFile) else []
        if len(results) != len(molIds):
          print(f'Skipping {name}: {len(results)} MM-GBSA results for {len(molIds)} input poses')
          continue

        self.splitMaeComplex(outMae)
        for k, (molId, row) in enumerate(zip(molIds, results), start=1):
          newMol = SmallMolecule()
          newMol.copy(inMols[molId], copyId=False)
          newMol.setProteinFile(os.path.join(outDir, f'{name}-out_receptor{k}.maegz'))
          newMol.setPoseFile(os.path.join(outDir, f'{name}-out_ligand{k}.maegz'))
          newMol._mmGBSA_BindEnergy = pwobj.Float(row[BIND_COLUMN])
          outputSet.append(newMol)
          energyRows.append({'molId': molId, 'poseFile': newMol.getPoseFile(), **row})

      self.writeEnergiesFile(energyRows)
      self._defineOutputs(outputSmallMolecules=outputSet)

    def writeEnergiesFile(self, rows):
      '''Writes every Prime energy term of all the scored poses in a single CSV'''
      if not rows:
        return
      with open(self.getEnergiesFile(), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    def getEnergiesFile(self):
      return self._getExtraPath(ENERGIES_FILE)

    def splitMaeComplex(self, maeFile):
      '''Splits each optimized complex of maeFile into <base>_receptor{k}.maegz and <base>_ligand{k}.maegz'''
      args = schrodingerPlugin.getMMshareDir('python/common/split_structure.py ')
      args += f'-m pdb -many_files {os.path.abspath(maeFile)} {os.path.abspath(maeFile)}'
      self.runJob(schrodingerPlugin.getHome('run'), args)

    # --------------------------- System functions --------------------------

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSmallMolecules'):
          summary.append(f'{len(self.outputSmallMolecules)} poses scored. Total binding energy (kcal/mol) stored '
                         f'as _mmGBSA_BindEnergy.')
          summary.append(f'Energy contributions of each pose (Coulomb, vdW, GB solvation, lipophilic, H-bond, '
                         f'covalent, packing, ligand strain, complex/receptor/ligand energies) are in: '
                         f'{self.getEnergiesFile()}')
        return summary

    def _validate(self):
        """ Try to find warnings on define params. """
        validations = []
        pSet = self.inputSmallMolecules.get()
        if not pSet.isDocked():
            validations.append('Input molecules must be docked first.\nSet: {} has not been docked'.format(pSet))
        if self.jobType.get() > ENERGY and self.flexOption.get() == FNONE:
            validations.append(f'"{self.getEnumText("jobType")}" samples flexible receptor residues: '
                               f'choose a "Flexible residues" option other than None')
        return validations

    # --------------------------- Parallelization functions --------------------
    def prepareLigandFile(self, mols, molLists, it):
      '''Appends (molId, MAE file) for each molecule. Maestro poses (e.g. from Glide) are used as they are'''
      for mol in mols:
        molFile = mol.getPoseFile()
        if isMaeFile(molFile):
          molLists[it].append((mol.getObjId(), os.path.abspath(molFile)))
          continue

        baseName = os.path.splitext(os.path.basename(molFile))[0]
        if not molFile.endswith('.sdf'):
          # Use openbabel first to convert to sdf (usable by Schrodinger)
          molFile = convertToSdf(self, molFile)

        if self.prepareLigand.get():
          # Prepare sdf file using LigPrep
          tmpmaeFile = os.path.abspath(self._getTmpPath(baseName + '_tmp.maegz'))
          args = " -WAIT -R h -a -isd {} -omae {}".format(molFile, tmpmaeFile)
          self.runJob(progLigPrep, args, cwd=self._getTmpPath())
          # LigPrep may output several states (tautomers, ionization) per pose: keep only the first one,
          # so each pose stays one entry in the complex file
          maeFile = os.path.abspath(self._getTmpPath(baseName + MAEFILE_EXTENSION))
          self.runJob(structConvertProg, f'-n 1 {tmpmaeFile} {maeFile}', cwd=self._getTmpPath())
          os.remove(tmpmaeFile)

        else:
          # Convert sdf to mae file
          maeFile = os.path.abspath(self._getTmpPath(mol.getUniqueName() + MAEFILE_EXTENSION))
          self.runJob(structConvertProg, '{} {}'.format(molFile, maeFile))

        molLists[it].append((mol.getObjId(), maeFile))

    def getMMGBSAArgs(self):
      args = ''
      if self.flexOption.get() == FTHRES:
        args += f'-rflexdist {self.flexDist.get()} -rflexgroup {self.getEnumText("flexGroup").lower()} '
      elif self.flexOption.get() == FAUT:
        args += f'-target_flexibility -target_flexibility_cutoff {self.flexTargetCutOff.get()} '

      if self.jobType.get() != REAL_MIN:
        args += f'-job_type {JOB_TYPES[self.jobType.get()]} '
      if self.forceField.get() != 0:
        args += f'-prime_opt OPLS_VERSION={FORCE_FIELDS[self.forceField.get()]} '

      if self.useCharges.get():
        args += '-use_ligand_charges '
      if self.rigidLigand.get():
        args += '-rigid_body '
      return args

    def performMMGBSA(self, complexFiles, molLists, it):
      for complexFile in complexFiles:
        name = getBaseName(complexFile)
        outDir = self._getExtraPath(name)
        if os.path.exists(os.path.join(outDir, f'{name}-out.csv')):
          print(f'{name}: already computed, skipping')
          continue
        os.makedirs(outDir, exist_ok=True)
        args = f"-WAIT {os.path.abspath(complexFile)} -out_type COMPLEX "
        args += self.getMMGBSAArgs()

        # prime_mmgbsa exits fine when it cannot get licenses, so check its log and wait for the tokens to be freed
        wait = LICENSE_WAIT_MIN
        while True:
          print('Performing prime mm-gbsa on : ', complexFile)
          self.runJob(mmgbsaProg, args, cwd=outDir)
          if os.path.exists(os.path.join(outDir, f'{name}-out.csv')) or not self.licenseFailed(outDir, name):
            break
          print(f'{name}: no Prime licenses available, retrying in {wait} s')
          time.sleep(wait)
          wait = min(2 * wait, LICENSE_WAIT_MAX)

    def licenseFailed(self, outDir, name):
      logFile = os.path.join(outDir, f'{name}.log')
      if not os.path.exists(logFile):
        return False
      with open(logFile) as f:
        return LICENSE_ERROR in f.read()

    # --------------------------- Utils functions --------------------
    def getNMolsPerThread(self, nMols, nt):
      '''Return a list with the number of element to process per thread'''
      if nt > nMols:
        nt = nMols

      nSubMols = [0 for _ in range(nt)]
      for i in range(nMols):
          nSubMols[i % nt] += 1
      return nSubMols

    def getInputComplexDir(self):
      cDir = self._getExtraPath('inputComplex')
      os.makedirs(cDir, exist_ok=True)
      return os.path.abspath(cDir)

    def getComplexFile(self, i):
      return os.path.join(self.getInputComplexDir(), f'complex_{i}{MAEFILE_EXTENSION}')

    def getComplexFiles(self):
      cDir = self.getInputComplexDir()
      return sorted(os.path.join(cDir, f) for f in os.listdir(cDir) if f.endswith(MAEFILE_EXTENSION))

    def getPdbFile(self, proteinFile):
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
      '''Return the corresponding MAE file for the receptor the molecule is docked to.
      Maestro receptors (e.g. from Glide/PrepWizard) are used as they are'''
      proteinFile = getattr(molSet, 'structFile', None)
      proteinFile = proteinFile.get() if proteinFile else molSet.getProteinFile()
      if isMaeFile(proteinFile):
        return os.path.abspath(proteinFile)

      pdbFile = self.getPdbFile(proteinFile)
      targetName = os.path.splitext(os.path.basename(pdbFile))[0]
      targetMaeFile = os.path.abspath(self._getTmpPath(targetName + MAEFILE_EXTENSION))
      if self.prepareReceptor.get():
        # Protonation + H-bond optimization, hydrogen-only minimization (-fix) keeps the docked heavy-atom frame
        args = '-WAIT -noepik -fix '
        args += '%s %s' % (os.path.abspath(pdbFile), os.path.abspath(targetMaeFile))
        self.runJob(progPrepWizard, args, cwd=self._getTmpPath())
      else:
        self.runJob(structConvertProg, '{} {}'.format(pdbFile, targetMaeFile), cwd=self._getTmpPath())
      return targetMaeFile

    def parseOutputCSVFile(self, csvFile):
      '''Returns the CSV rows (dicts of every Prime energy term), in the same order as the input poses'''
      with open(csvFile) as f:
        return list(csv.DictReader(f))

    def getOriginalReceptorFile(self):
        return self.inputSmallMolecules.get().getProteinFile()
