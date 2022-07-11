# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
from schrodingerScipion.protocols.protocol_desmond_systemPrep import *
from schrodingerScipion.utils.utils import getChargeFromMAE

import pyworkflow.wizard as pwizard
from subprocess import check_call

from pwchem.wizards import SelectElementWizard

SelectElementWizard().addTarget(protocol=ProtSchrodingerDesmondSysPrep,
                               targets=['inputLigand'],
                               inputs=['inputSetOfMols'],
                               outputs=['inputLigand'])

class GetSoluteCharge(pwizard.Wizard):
    """Calculates the charge of the input atom structure"""
    _targets = [(ProtSchrodingerDesmondSysPrep, ['solCharge'])]

    def getSoluteFile(self, protocol):
        if protocol.inputFrom.get() == STRUCTURE:
            if type(protocol.inputStruct.get()) == SchrodingerAtomStruct:
                inSoluteFile = protocol.inputStruct.get().getFileName()
                if inSoluteFile.endswith('gz'):
                    structName = os.path.splitext(os.path.basename(inSoluteFile))[0]
                    soluteFile = os.path.join('/tmp', structName + '.mae')
                    check_call('zcat {} > {}'.format(os.path.abspath(inSoluteFile), soluteFile), shell=True)
            else:
                pdbFile = protocol.inputStruct.get().getFileName()
                structName = os.path.splitext(os.path.basename(pdbFile))[0]
                soluteFile = os.path.join('/tmp', structName + '.mae')
                if not os.path.exists(soluteFile):
                    check_call('{} {} {}'.format(structConvertProg, pdbFile, soluteFile), shell=True)

        elif protocol.inputFrom.get() == LIGAND:
            soluteFile = os.path.join('/tmp', 'complexSolute.mae')
            if not os.path.exists(soluteFile):
                mol = protocol.getSpecifiedMol()
                molMaeFile = os.path.join('/tmp', mol.getUniqueName() + '.maegz')
                check_call('{} {} {}'.format(structConvertProg, mol.getPoseFile(), molMaeFile), shell=True)

                if hasattr(mol, 'structFile'):
                    targetMaeFile = mol.structFile
                else:
                    targetFile = protocol.inputSetOfMols.get().getProteinFile()
                    targetName = os.path.splitext(os.path.basename(targetFile))[0]
                    targetMaeFile = os.path.join('/tmp', targetName + '.maegz')
                    check_call('{} {} {}'.format(structConvertProg, targetFile, targetMaeFile), shell=True)

                check_call('zcat {} {} > {}'.format(molMaeFile, targetMaeFile, soluteFile), shell=True)


        return soluteFile


    def show(self, form, *params):
        protocol = form.protocol
        try:
            soluteFile = self.getSoluteFile(protocol)
        except Exception as e:
            print("ERROR: ", e)
            return

        charge = getChargeFromMAE(soluteFile)
        form.setVar('solCharge', int(charge))

