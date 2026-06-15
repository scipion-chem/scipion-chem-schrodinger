# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
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
from pwchem.wizards import SelectChainWizardQT, SelectResidueWizardQT, AddResidueWizard, SelectElementWizard
from ..protocols import ProtSchrodingerPrepWizard, ProtSchrodingerPrime, ProtSchrodingerIFD
from ..viewers import ProtGlideDockingViewer

SelectChainWizardQT().addTarget(protocol=ProtSchrodingerPrepWizard,
                                targets=['chain_name'],
                                inputs=['inputAtomStruct'],
                                outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtSchrodingerPrime,
                                targets=['resChain'],
                                inputs=['inputStructure'],
                                outputs=['resChain'])

SelectResidueWizardQT().addTarget(protocol=ProtSchrodingerPrime,
                                  targets=['resPosition'],
                                  inputs=['inputStructure', 'resChain'],
                                  outputs=['resPosition'])

AddResidueWizard().addTarget(protocol=ProtSchrodingerPrime,
                             targets=['addResidue'],
                             inputs=['resChain', 'resPosition'],
                             outputs=['residueList'])

SelectChainWizardQT().addTarget(protocol=ProtSchrodingerIFD,
                                targets=['selChain'],
                                inputs=[{'fromPockets': ['inputAtomStruct', 'inputStructROIs']}],
                                outputs=['selChain'])

SelectResidueWizardQT().addTarget(protocol=ProtSchrodingerIFD,
                                  targets=['selResidue'],
                                  inputs=[{'fromPockets': ['inputAtomStruct', 'inputStructROIs']}, 'selChain'],
                                  outputs=['selResidue'])

AddResidueWizard().addTarget(protocol=ProtSchrodingerIFD,
                             targets=['residuesAdd'],
                             inputs=['selChain', 'selResidue'],
                             outputs=['residuesAdd'])

AddResidueWizard().addTarget(protocol=ProtSchrodingerIFD,
                             targets=['residuesOmit'],
                             inputs=['selChain', 'selResidue'],
                             outputs=['residuesOmit'])

AddResidueWizard().addTarget(protocol=ProtSchrodingerIFD,
                             targets=['residuesTrim'],
                             inputs=['selChain', 'selResidue'],
                             outputs=['residuesTrim'])

SelectResidueWizardQT().addTarget(protocol=ProtSchrodingerIFD,
                                  targets=['residuesHelix'],
                                  inputs=[{'fromPockets': ['inputAtomStruct', 'inputStructROIs']}, 'selChain'],
                                  outputs=['residuesHelix'])

SelectElementWizard().addTarget(protocol=ProtGlideDockingViewer,
                               targets=['displayMoleculeDock'],
                               inputs=['moleculeLabels'],
                               outputs=['displayMoleculeDock'])

SelectElementWizard().addTarget(protocol=ProtGlideDockingViewer,
                               targets=['displaySingleDock'],
                               inputs=['singleLabels'],
                               outputs=['displaySingleDock'])

SelectElementWizard().addTarget(protocol=ProtGlideDockingViewer,
                               targets=['displayPymolPLIP'],
                               inputs=['singleLabels'],
                               outputs=['displayPymolPLIP'])