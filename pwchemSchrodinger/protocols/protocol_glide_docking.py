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
import os, shutil, threading, subprocess, time, glob

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam, IntParam, \
    STEPS_PARALLEL, LabelParam
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import relabelAtomsMol2, calculate_centerMass, getBaseName, performBatchThreading, \
    organizeThreads, insistentRun
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC

from .. import Plugin as schrodinger_plugin
from ..protocols.protocol_preparation_grid import ProtSchrodingerGrid
from ..utils.utils import putMolFileTitle, convertMAEMolSet

glideProg = schrodinger_plugin.getHome('glide')
progLigPrep = schrodinger_plugin.getHome('ligprep')
structConvertProg = schrodinger_plugin.getHome('utilities/structconvert')
structCatProg = schrodinger_plugin.getHome('utilities/structcat')
propListerProg = schrodinger_plugin.getHome('utilities/proplister')
maeSubsetProg = schrodinger_plugin.getHome('utilities/maesubset')

dockMethodDic = {0: 'confgen', 1: 'rigid', 2: 'mininplace', 3: 'inplace'}
dockPrecisionDic = {0: 'HTVS', 1: 'SP', 2: 'XP'}

class ProtSchrodingerGlideDocking(ProtSchrodingerGrid):
    """Calls glide to perform a docking of a set of compounds in a structural region defined by a grid.
       It is assumed that the input library of ligands is already prepared.

       The dockinsScore is the Glide Docking Score and it is measured in kcal/mol"""
    _label = 'docking (glide)'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

        self.tDic.update({'Flexible dock (confgen)': 'confgen', 'Rigid dock (rigid)': 'rigid',
                          'Refine (do not dock, mininplace)': 'mininplace',
                          'Score in place (do not dock, inplace)': 'inplace',
                          'Low (HTVS)': 'HTVS', 'Medium (SP)': 'SP', 'High (XP)': 'XP',
                          'Penalize': 'penal', 'Free': 'free', 'Fixed': 'fixed', 'Trans': 'trans',
                          'Generalized': 'gen'})

        self.paramsDic2 = {'posesPerLig': 'POSES_PER_LIG', 'canonicalize': 'CANONICALIZE',
                           'ligMaecharges': 'LIG_MAECHARGES', 'ligVScale': 'LIG_VSCALE', 'ligCCut': 'LIG_CCUT',
                           'sampleNinversions': 'SAMPLE_N_INVERSIONS', 'sampleRings': 'SAMPLE_RINGS',
                           'epikPenalties': 'EPIK_PENALTIES', 'skipMetalEpik': 'SKIP_EPIK_METAL_ONLY',
                           'expandedSampling': 'EXPANDED_SAMPLING', 'poseRMSD': 'POSE_RMSD',
                           'poseDisp': 'POSE_DISPLACEMENT', 'cvCutOffFilter': 'CV_CUTOFF',
                           'hbondCutOffFilter': 'HBOND_CUTOFF', 'metalCutOffFilter': 'METAL_CUTOFF',
                           'postDock': 'POSTDOCK', 'postDockN': 'POSTDOCK_NPOSE'}
        self.enumParamsDic2 = {'dockingMethod': 'DOCKING_METHOD', 'dockingPrecision': 'PRECISION',
                               'amideModel': 'AMIDE_MODE'}

    def _defineGlideReceptorParams(self, form):
        group = form.addGroup('Receptor')
        group.addParam('fromPockets', EnumParam, label='Dock on : ', default=1,
                       choices=['Whole protein', 'SetOfStructROIs', 'Schrodinger Grids'],
                       help='Whether to dock on a whole protein surface or on specific regions defines as StructROIs or'
                            ' Schrodinger Grids')

        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       condition='fromPockets==0', label='Input structure:',
                       help='Input atomic structure to perform the docking on. We suggest the use of Structural ROIs'
                            'to speed up the docking process instead of docking on the whole protein')
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='fromPockets == 0', allowsNull=False,
                       help='Radius of the Schrodinger grid for the whole protein')

        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                       condition='fromPockets==1', label='Input structural ROIs:',
                       help='Input structural ROIs defining the space where the docking will be performed')
        group.addParam('inputGridSet', PointerParam, pointerClass="SetOfSchrodingerGrids",
                       condition='fromPockets==2', label='Input grids:',
                       help='Input grids defining the space where the docking will be performed')
        return group

    def _defineGeneralGlideParams(self, form, condition='True'):
        group = form.addGroup('General params', condition=condition)
        group.addParam('dockingThreads', IntParam, default=2, label='No. of docking subjobs: ',
                       help='Maximum number of docking subjobs. Each subjobs will be run in a processor. \n'
                            'This number should generally be lower than the numbe rof threads, since each subjob will '
                            'require a license token.')
        group.addParam('dockingMethod', EnumParam, default=0, label='Docking method: ',
                       choices=['Flexible dock (confgen)', 'Rigid dock (rigid)', 'Refine (do not dock, mininplace)',
                                'Score in place (do not dock, inplace)'],
                       help='Glide method to use for docking')
        group.addParam('dockingPrecision', EnumParam, default=0, label='Docking precision: ',
                       choices=['Low (HTVS)', 'Medium (SP)', 'High (XP)'],
                       help='You may use a low to high strategy. HTVS takes about 2 s/ligand, SP about 10s, and XP '
                            'about 10 min.')
        group.addParam('posesPerLig', IntParam, default=5, label='No. Poses to report per ligand: ',
                       help='Maximum number of final poses to report per ligand')

        group.addParam('canonicalize', BooleanParam, default=True,
                       label='Canonicalize the input ligand structure :', expertLevel=LEVEL_ADVANCED,
                       help='Canonicalize the input structure by discarding the coordinates and regenerating the '
                            'structure from the connectivity and stereochemistry. Takes about 1 sec per ligand.')
        return group

    def _defineLigandParams(self, form, condition='True'):
        group = form.addGroup('Ligand params', condition=condition)
        group.addParam('ligMaecharges', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                       label='Use partial charges from the input ligand structure:',
                       help='Use partial charges from the input ligand structure.')
        group.addParam('ligVScale', FloatParam, default=0.8,
                       label='Scaling factor for van der Waals radii scaling: ',
                       help='Scaling factor for van der Waals radii scaling')
        group.addParam('ligCCut', FloatParam, default=0.15,
                       label='Partial charge cutoff for van der Waals radii scaling: ',
                       help='Partial charge cutoff for van der Waals radii scaling.')

        group.addParam('sampleNinversions', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample pyramid nitrogen inversions: ', expertLevel=LEVEL_ADVANCED,
                       help='Sample nitrogen inversions if set to True and DOCKING_METHOD is set to confgen.')
        group.addParam('sampleRings', BooleanParam, default=True, condition='dockingMethod==0',
                       label='Sample rings: ', expertLevel=LEVEL_ADVANCED,
                       help='Sample rings if set to True and DOCKING_METHOD is set to confgen.')
        group.addParam('amideModel', EnumParam, default=4, label='Amide model: ', expertLevel=LEVEL_ADVANCED,
                       choices=['Penalize', 'Free', 'Fixed', 'Trans', 'Generalized'],
                       help='Amide bond sampling mode.\n'
                            'Penalize: penalize nonplanar conformation\nFree: vary conformation\n'
                            'Fixed: retain original conformation\nTrans: allow trans conformation only\n'
                            'Generalized: use generalized torsion controls defined in torcontrol.txt.')

        group.addParam('epikPenalties', BooleanParam, default=False, label='Epik penalties: ',
                       expertLevel=LEVEL_ADVANCED,
                       help='Apply penalties for ionization or tautomeric states calculated by Epik')
        group.addParam('skipMetalEpik', BooleanParam, default=True, label='Skip Epik metal only: ',
                       expertLevel=LEVEL_ADVANCED,
                       help='Skip Epik-generated states of ligands that are designed for binding to metals. '
                            'This option is useful if the receptor has a metal but the ligand does not bind to it. '
                            'These states are skipped by default if the receptor does not have a metal.')
        group.addParam('expandedSampling', BooleanParam, default=False, label='Expanded sampling:',
                       expertLevel=LEVEL_ADVANCED,
                       help='Expand the sampling by bypassing the elimination of poses in the rough scoring stage. '
                            'Useful for fragment docking.')
        return group

    def _defineOutputGlideParams(self, form, condition='True'):
        group = form.addGroup('Output params', condition=condition)
        group.addParam('poseRMSD', FloatParam, default=0.5,
                       label='Maximum RMS deviation used in clustering (A): ',
                       help='RMS deviation used in clustering to discard poses, in angstroms')
        group.addParam('poseDisp', FloatParam, default=1.3, expertLevel=LEVEL_ADVANCED,
                       label='Maximum atomic displacement used in clustering (A): ',
                       help='Maximum atomic displacement used in clustering to discard poses, in angstroms.')
        group.addParam('cvCutOffFilter', FloatParam, default=0.0, label='Filter Coulomb-van der Waals cutoff: ',
                       expertLevel=LEVEL_ADVANCED,
                       help='Reject poses with Coulomb-van der Waals energy greater than cutoff kcal/mol.')
        group.addParam('hbondCutOffFilter', FloatParam, default=0.0, label='Filter H-bond score cutoff: ',
                       expertLevel=LEVEL_ADVANCED,
                       help='Reject poses with H-bond score greater than cutoff')
        group.addParam('metalCutOffFilter', FloatParam, default=10.0, label='Filter Metal score cutoff: ',
                       expertLevel=LEVEL_ADVANCED,
                       help='Reject poses with metal score greater than cutoff')

        group.addParam('postDock', BooleanParam, default=True, label='Use post docking minimization: ',
                       expertLevel=LEVEL_ADVANCED, help='Use post docking minimization')
        group.addParam('postDockN', IntParam, default=5, expertLevel=LEVEL_ADVANCED,
                       condition='postDock',
                       label='Number of poses to use in post-docking minimization: ',
                       help='Number of poses to use in post-docking minimization. Maestro sets this number to '
                            '10 for XP. Default: 5.')

        return group

    def _defineParams(self, form):
        notGridsCondition = 'fromPockets!=2'

        form.addSection(label='Input')
        self._defineGlideReceptorParams(form)
        self._defineInnerGridParams(form, condition=notGridsCondition)
        self._defineOuterGridParams(form, condition=notGridsCondition)

        group = form.addGroup('Ligands')
        group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules:', help='Input small molecules to be docked with Glide')

        group.addParam('convertOutput2Mol2', BooleanParam, label='Convert output to mol2: ', default=False,
                       help='Whether to convert to output molecules to mol2 files instead of keeping the mae files')

        form = self._defineGridSection(form, condition=notGridsCondition)

        form.addSection(label='Docking')
        self._defineGeneralGlideParams(form)
        self._defineLigandParams(form)

        form.addSection(label='Docking output')
        self._defineOutputGlideParams(form)

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        nt = self.dockingThreads.get()
        convStep = self._insertFunctionStep('convertStep', prerequisites=[])
        cSteps = [convStep]

        inputGrids = self.inputGridSet.get()
        if self.fromPockets.get() == 0:
            dStep = self._insertFunctionStep('gridStep', self.inputAtomStruct.get(), prerequisites=cSteps)
            cSteps, inputGrids = [dStep], [self.inputAtomStruct.get()]

        elif self.fromPockets.get() == 1:
            gridSteps = []
            for pocket in self.inputStructROIs.get():
                dStep = self._insertFunctionStep('gridStep', pocket.clone(), prerequisites=cSteps)
                gridSteps.append(dStep)
            cSteps, inputGrids = gridSteps, self.inputStructROIs.get()

        dockSteps = []
        for i, grid in enumerate(inputGrids):
            nThreads = organizeThreads(len(inputGrids), nt)
            dStep = self._insertFunctionStep('dockingStep', grid.clone(), nThreads[i], prerequisites=cSteps)
            dockSteps.append(dStep)
        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def getLigSetFormats(self, ligSet):
        formats = []
        for lig in ligSet:
            ext = os.path.splitext(lig.getFileName())[1]
            formats.append(ext)
        formats = list(set(formats))
        return formats

    def convertStep(self):
        nt = self.numberOfThreads.get()
        # Convert receptor
        inFile = self.getOriginalReceptorFile()

        if '.mae' not in inFile:
            maeFile = self._getExtraPath('inputReceptor.maegz')
            prog = schrodinger_plugin.getHome('utilities/prepwizard')
            args = ' -WAIT -noprotassign -noimpref -noepik {} {}'. \
                format(os.path.abspath(inFile), os.path.abspath(maeFile))
            self.runJob(prog, args, cwd=self._getExtraPath())
        else:
            ext = os.path.splitext(inFile)[1]
            shutil.copy(inFile, self._getExtraPath('inputReceptor{}'.format(ext)))

        # Create ligand files
        ligSet = self.inputLibrary.get()

        ligFormats = self.getLigSetFormats(ligSet)
        ligFiles = list(set([lig.getFileName() for lig in ligSet]))
        if len(ligFormats) > 1 or ligFormats[0] not in ['.mae', '.maegz', '.sdf', '.mol2']:
            allLigandsFile = self.getAllLigandsFile(format='.mol2')
            performBatchThreading(self.createLigandsFile, ligFiles, nt)
            with open(allLigandsFile, 'w') as f:
                for it in range(nt):
                    with open(self.getAllLigandsFile(suffix=it, format='.mol2')) as fLig:
                        f.write(fLig.read())
        else:
            allLigandsFile = self.getAllLigandsFile(format=ligFormats[0])
            with open(allLigandsFile, 'w') as f:
                for ligFile in ligFiles:
                    with open(ligFile) as fLig:
                        f.write(fLig.read())

    def gridStep(self, pocket):
        if self.fromPockets.get() == 0:
            inAS = self.inputAtomStruct.get()
            pdbFile = inAS.getFileName()
            if not pdbFile.endswith('.pdb'):
                pdbFile = inAS.convert2PDB(self._getExtraPath('inputStructure.pdb'))

            pocket = self.radius.get()
            _, x, y, z = calculate_centerMass(pdbFile)
            pId = 1
        else:
            x, y, z = pocket.calculateMassCenter()
            pId = pocket.getObjId()

        fnGridDir = self.getGridDir(pId)
        gridName = self.getGridName(pId)
        # print('Grid name: ', gridName)
        makePath(fnGridDir)

        fnJob = os.path.abspath(os.path.join(fnGridDir, gridName)) + '.inp'
        with open(fnJob, 'w') as fh:
            fh.write("GRIDFILE %s.zip\n" % gridName)
            fh.write("OUTPUTDIR %s\n" % fnGridDir)
            fh.write("RECEP_FILE %s\n" % os.path.abspath(self.getInputMaeFile()))
            fh.write("INNERBOX {},{},{}\n".format(*(self.getInnerBox(pocket))))
            fh.write("ACTXRANGE %d\n" % self.getOuterBox(pocket)[0])
            fh.write("ACTYRANGE %d\n" % self.getOuterBox(pocket)[1])
            fh.write("ACTZRANGE %d\n" % self.getOuterBox(pocket)[2])
            fh.write("OUTERBOX {},{},{}\n".format(*(self.getOuterBox(pocket))))
            fh.write("GRID_CENTER %s,%s,%s\n" % (x, y, z))

            argDic = self.getArgsDic(self.paramsDic, self.enumParamsDic)
            for glideArg, value in argDic.items():
                fh.write(f"{glideArg} {value}\n")

        args = f"-WAIT -LOCAL {gridName}.inp"
        insistentRun(self, glideProg, args, cwd=fnGridDir)

    def dockingStep(self, grid, nt):
        gridId = grid.getObjId()
        if self.fromPockets.get() == 0:
            gridId = 1
            gridDir = self.getGridDir(gridId)
        elif self.fromPockets.get() == 1:
            gridDir = self.getGridDir(grid.getObjId())

        elif self.fromPockets.get() == 2:
            gridDir = self._getExtraPath('grid_{}/'.format(gridId))

            makePath(gridDir)
            fnGrid = os.path.join(gridDir, "grid_{}.zip".format(gridId))
            if not os.path.exists(fnGrid):
                shutil.copy(grid.getFileName(), fnGrid)

        fnIn = os.path.join(gridDir, 'job_{}.inp'.format(gridId))
        if not os.path.exists(fnIn):
            with open(fnIn, 'w') as fhIn:
                fhIn.write("GRIDFILE %s\n" % ("grid_{}.zip".format(gridId)))
                fhIn.write("LIGANDFILE {}\n".format(os.path.abspath(self.getAllLigandsFile())))

                argDic = self.getArgsDic(self.paramsDic2, self.enumParamsDic2)
                for glideArg, value in argDic.items():
                    fhIn.write(f"{glideArg} {value}\n")

        args = "-WAIT -RESTART -LOCAL job_{}.inp".format(gridId)
        self.runJob(glideProg, args, cwd=gridDir)

        if os.path.exists(os.path.join(gridDir, "job_{}_pv.maegz".format(gridId))):
            self.runJob(propListerProg,
                        '-p "title" -p "docking score" -p "glide ligand efficiency" -p "glide ligand efficiency sa" '
                        '-p "glide ligand efficiency ln" -c -o %s %s'%\
                        ("job_{}_pv.csv".format(gridId), "job_{}_pv.maegz".format(gridId)),
                        cwd=gridDir)
        else:
            print('Failed to find ligands for grid {}'.format(gridId))

    def performOutputParsing(self, gridDirs, molLists, it, smallDict):
        allSmallList = []
        for gridDir in gridDirs:
            gridId = gridDir.split('_')[1]
            gridDir = 'grid_{}/'.format(gridId)

            smallList = []
            fnPv = self._getExtraPath(gridDir + 'job_{}_pv.maegz'.format(gridId))
            with open(self._getExtraPath(gridDir + 'job_{}_pv.csv'.format(gridId))) as fhCsv:
                fhCsv.readline()
                fhCsv.readline()
                for i, line in enumerate(fhCsv.readlines()):
                    tokens = line.split(',')
                    fnBase = getBaseName(tokens[0])

                    small = SmallMolecule()
                    small.copy(smallDict[fnBase], copyId=False)
                    small._energy = pwobj.Float(tokens[1])
                    small.ligandEfficiency = pwobj.Float(tokens[2])
                    small.ligandEfficiencySA = pwobj.Float(tokens[3])
                    small.ligandEfficiencyLn = pwobj.Float(tokens[4])
                    small.setMolClass('Schrodinger')
                    small.setDockId(self.getObjId())
                    small.setGridId(gridId)

                    _, posFile = self.divideMaeComplex(fnPv, posIdx=i+1)
                    small.setProteinFile(self.getInputMaeFile())
                    small.setPoseFile(posFile)
                    small.setPoseId(i + 1)

                    smallList.append(small)

            allSmallList += smallList

        molLists[it] = allSmallList

    def createOutputStep(self):
        nt = self.numberOfThreads.get()
        smallDict = {}
        for small in self.inputLibrary.get():
            fnSmall = small.getFileName()
            fnBase = getBaseName(fnSmall)
            if fnBase not in smallDict:
                smallDict[fnBase] = small.clone()
            molName = small.getMolName()
            if molName not in smallDict:
                smallDict[molName] = small.clone()

        fnStruct = self.getInputMaeFile()

        gridDirs = self.getGridDirs(complete=True)
        allSmallList = performBatchThreading(self.performOutputParsing, gridDirs, nt, cloneItem=False,
                                             smallDict=smallDict)

        if self.convertOutput2Mol2:
            outDir = os.path.abspath(self._getExtraPath('outputSmallMolecules'))
            os.mkdir(outDir)
            allSmallList = convertMAEMolSet(allSmallList, outDir, nt, updateSet=False, subset=False)

        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for small in allSmallList:
            outputSet.append(small)

        outputSet.setDocked(True)
        outputSet.proteinFile.set(self.getOriginalReceptorFile())
        outputSet.structFile = pwobj.String(fnStruct)
        if len(outputSet) > 0:
            self._defineOutputs(outputSmallMolecules=outputSet)
        else:
            print('No output docking files were generated or no poses were found')

    def divideMaeComplex(self, maeFile, posIdx=1, outDir=None):
      if not outDir:
        outDir = os.path.dirname(maeFile)
      molFile, recFile = os.path.join(outDir, getBaseName(maeFile) + f'_lig_{posIdx}.maegz'), \
                         os.path.join(outDir, getBaseName(maeFile) + '_rec.maegz')
      args = f' -n 1 {os.path.abspath(maeFile)} -o {os.path.abspath(recFile)}'
      subprocess.run(f'{maeSubsetProg} {args}', check=True, capture_output=True, text=True, shell=True, cwd=outDir)
      args = f' -n {posIdx+1} {os.path.abspath(maeFile)} -o {os.path.abspath(molFile)}'
      subprocess.run(f'{maeSubsetProg} {args}', check=True, capture_output=True, text=True, shell=True, cwd=outDir)

      return recFile, molFile

    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors

    ###########################  Utils functions #####################
    def convertOutputStep(self, curMolSet, molLists, it, outDir):
        for i, mol in enumerate(curMolSet):
            molName = mol.getUniqueName()
            fnOut = os.path.join(outDir, f'{molName}.mol2')

            if not os.path.exists(fnOut):
                posFile = os.path.abspath(mol.getPoseFile())
                args = f'{posFile} {os.path.abspath(fnOut)}'
                subprocess.call([structConvertProg, *args.split()])
                fnOut = relabelAtomsMol2(fnOut)

            mol.setPoseFile(fnOut)
            molLists[it].append(mol)


    def convert2mol2(self, fnSmall, it):
        baseName = getBaseName(fnSmall)
        outFile = os.path.abspath(self._getTmpPath('{}.mol2'.format(baseName)))
        args = f' -i "{os.path.abspath(fnSmall)}" -o {outFile} --outputDir {self._getTmpPath()}'
        pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=self._getTmpPath(), popen=True)

        while not os.path.exists(outFile):
            time.sleep(0.2)
        return outFile

    def getAllLigandsFile(self, suffix='', format=''):
        if format:
            ligFile = self._getTmpPath(f'allMoleculesFile{suffix}{format}')
        else:
            for file in os.listdir(self._getTmpPath()):
                if f'allMoleculesFile{suffix}' in file:
                    ligFile = self._getTmpPath(file)
                    break
        return ligFile

    def getGridDirs(self, complete=False):
        gridDirs = []
        for dir in os.listdir(self._getExtraPath()):
            if dir.startswith('grid_'):
                if not complete:
                    gridDirs.append(dir)
                else:
                    gridId = dir.split('_')[1]
                    if os.path.exists(self._getExtraPath(dir, "job_{}_pv.maegz".format(gridId))):
                        gridDirs.append(dir)
                    else:
                        print('No good poses found in grid ' + gridId)
        return gridDirs

    def mergeMAEfiles(self):
        maeFiles = []
        for gridDir in self.getGridDirs(complete=True):
            gridId = gridDir.split('_')[1]
            maeFiles.append(gridDir + '/job_{}_pv.maegz'.format(gridId))

        noRecMaeFiles = [maeFiles[0]]
        for mFile in maeFiles[1:]:
            mFile = os.path.abspath(self._getExtraPath(mFile))
            tFile = os.path.abspath(self._getTmpPath(getBaseName(mFile) + '.maegz'))
            args = f' -n 2: {mFile} -o {tFile}'
            self.runJob(maeSubsetProg, args, cwd=self._getTmpPath())
            noRecMaeFiles.append(tFile)

        outName = 'dockedMolecules.maegz'
        if len(noRecMaeFiles) > 1:
            command = 'zcat {} | gzip -c > {}'.format(' '.join(noRecMaeFiles), outName)
            self.runJob('', command, cwd=self._getExtraPath())
        else:
            os.symlink(noRecMaeFiles[0], self._getExtraPath(outName))
        return self._getExtraPath(outName)


    def createLigandsFile(self, ligFiles, molLists, it):
        curAllLigandsFile = self.getAllLigandsFile(suffix=it, format='.mol2')
        with open(curAllLigandsFile, 'w') as fh:
            for fnSmall in ligFiles:
                if not fnSmall.endswith('.mol2'):
                    fnSmall = self.convert2mol2(fnSmall, it)
                putMolFileTitle(fnSmall, ext='mol2')
                with open(fnSmall) as fhLigand:
                    fh.write(fhLigand.read())

    def getGridDir(self, pocketId):
        return self._getExtraPath(self.getGridName(pocketId))

    def getGridName(self, pocketId):
        return 'grid_{}'.format(pocketId)

    def getInputMaeFile(self):
        maeFile = glob.glob(self._getExtraPath('inputReceptor.mae*'))[0]
        return maeFile

    def getInnerBox(self, pocket):
        if self.fromPockets.get() == 0:
            return 3*[int(pocket * 2 * self.diameterNin.get())]
        else:
            if self.innerAction.get() == 0:
                return self.innerX.get(), self.innerY.get(), self.innerZ.get()
            else:
                diam = int(pocket.getDiameter() * self.diameterNin.get())
                return [diam, diam, diam]

    def getOuterBox(self, pocket):
        if self.fromPockets.get() == 0:
            return 3*[int(pocket * 2 * self.diameterNout.get())]
        else:
            if self.outerAction.get() == 0:
                return self.outerX.get(), self.outerY.get(), self.outerZ.get()
            else:
                diam = int(pocket.getDiameter() * self.diameterNout.get())
                return [diam, diam, diam]

    def getOriginalReceptorFile(self):
        if self.fromPockets.get() == 0:
            inAs = self.inputAtomStruct.get()
            if hasattr(inAs, '_maeFile'):
                return getattr(inAs, '_maeFile').get()
            else:
                return inAs.getFileName()

        elif self.fromPockets.get() == 1:
            inFile = self.inputStructROIs.get().getMAEFile()
            if not inFile:
                inFile = self.inputStructROIs.get().getProteinFile()

        elif self.fromPockets.get() == 2:
            inFile = self.inputGridSet.get().getProteinFile()
        return inFile



