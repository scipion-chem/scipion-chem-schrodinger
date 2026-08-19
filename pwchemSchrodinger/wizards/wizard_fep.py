# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""Wizards for the FEP+ RBFE/ABFE protocols. The multi-ligand picker is line-for-line the
same pattern gromacs/wizards/wizard_MDSimulation.py's SelectMultiLigandsWizard uses for
GromacsPmxRBFE - see that file for why it doesn't reuse pwchem's own
SelectMultiElementWizard.displayDialog() as-is."""

from pyworkflow.gui import ListTreeProviderString, dialog
import pyworkflow.object as pwobj
from pwem.objects import Pointer

from pwchem.wizards import SelectElementWizard, SelectMultiElementWizard

from ..protocols import ProtSchrodingerFepRBFE, ProtSchrodingerFepABFE

SelectElementWizard().addTarget(protocol=ProtSchrodingerFepABFE,
                               targets=['inputLigand'],
                               inputs=['inputSetOfMols'],
                               outputs=['inputLigand'])


class SelectMultiLigandsWizardFEP(SelectMultiElementWizard):
    """Same multi-select ListDialog as pwchem's own SelectMultiElementWizard (default
    selectmode='extended', Ctrl/Shift-click), additionally filling two more output vars
    with the first two selected names for ProtSchrodingerFepRBFE, so its existing
    single-ligand-A/B params (kept for the 2-ligand case) stay in sync with the
    multi-select list."""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        try:
            scipionSet = getattr(protocol, inputParam[0])
            if isinstance(scipionSet, Pointer):
                scipionSet = scipionSet.get()
            listOfElements = self.getListOfElements(protocol, scipionSet)
        except Exception as e:
            print("ERROR: ", e)
            return

        finalList = [pwobj.String(i) for i in listOfElements]
        provider = ListTreeProviderString(finalList)
        dlg = dialog.ListDialog(form.root, "Set items", provider,
                                "Select 2 or more items in the set\n"
                                "(Ctrl+Click or Shift+Click for multiple):")
        if not dlg.values:
            return
        selected = [v.get().strip() for v in dlg.values]
        form.setVar(outputParam[0], '\n'.join(selected))
        if len(outputParam) > 1 and len(selected) > 0:
            form.setVar(outputParam[1], selected[0])
        if len(outputParam) > 2 and len(selected) > 1:
            form.setVar(outputParam[2], selected[1])


SelectMultiLigandsWizardFEP().addTarget(protocol=ProtSchrodingerFepRBFE,
                                       targets=['selectedLigands'],
                                       inputs=['inputSetOfMols'],
                                       outputs=['selectedLigands', 'inputLigand', 'ligandB'])
