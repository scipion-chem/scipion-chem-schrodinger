# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
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
# Scipion em imports
import pyworkflow.protocol.params as params

# Scipion chem imports
from pwchem.viewers import SmallMoleculesViewer
import pyworkflow.viewer as pwviewer

# Plugin imports
from ..protocols import ProtSchrodingerQSAR
from ..viewers.viewers_data import MaestroView


class ProtSchrodingerQSARViewer(pwviewer.ProtocolViewer):
    _label = 'Viewer Schrodinger QSAR'
    _targets = [ProtSchrodingerQSAR, ProtSchrodingerQSARPharmacophore]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        form.addSection(label='Visualization of QSAR model')
        group = form.addGroup('Open Maestro GUI')
        group.addParam('displaySDF',
                       params.LabelParam,
                       label='Open molecules in Maestro',
                       help='Display the SDF file in Maestro GUI.')

        group = form.addGroup('Plots')
        group.addParam('displayPredictions',
                       params.LabelParam,
                       label='Plot predictions vs experimental',
                       help='Generate and display a scatter plot.')

    def _getVisualizeDict(self):
        visDic = super()._getVisualizeDict()
        visDic.update({
            'displaySDF': self._viewSDFMaestro,
            'displayPredictions': self._viewPredictionsPlot
        })
        return visDic

    def _viewSDFMaestro(self, paramName=None):
        model = self.protocol.SchrodingerQSARModel

        sdfFile = model.sdfFile.get()
        print(sdfFile)

        return [MaestroView(f'{sdfFile}')]

    def _viewPredictionsPlot(self, paramName=None):
        import pandas as pd
        import matplotlib.pyplot as plt
        import os

        model = self.protocol.SchrodingerQSARModel
        predFile = model.predictionsFile.get()

        df = pd.read_csv(os.path.abspath(predFile), sep=',')

        pred_cols = [col for col in df.columns if col.startswith('Pred(')]

        df['Predicted'] = df[pred_cols].mean(axis=1)

        plt.figure()

        plt.scatter(df['Activity'], df['Predicted'])

        plt.xlabel('Experimental Activity')
        plt.ylabel('Predicted Activity')
        plt.title('QSAR Model Performance')

        plt.grid()

        xmin = min(df['Activity'].min(), df['Predicted'].min())
        xmax = max(df['Activity'].max(), df['Predicted'].max())

        margin = (xmax - xmin) * 0.05 if xmax > xmin else 1.0

        plt.xlim(xmin - margin, xmax + margin)
        plt.ylim(xmin - margin, xmax + margin)
        plt.margins(0)

        plt.plot([xmin, xmax], [xmin, xmax], linestyle='--')

        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()
        plt.close()

        return []


