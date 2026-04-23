#!/usr/bin/env python3
"""
Convert Phase QSAR SDF output into a clean CSV.

Extracts:
- molecule name
- predicted activity (aggregated from r_phase_Pred_Activity(*))
"""

import sys
import pandas as pd
from rdkit import Chem


def sdfToCsv(inputSdf, outputCsv):

    suppl = Chem.SDMolSupplier(inputSdf)

    data = []

    for mol in suppl:
        if mol is None:
            continue

        name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if not name:
            continue

        preds = []
        for prop in mol.GetPropNames():
            if prop.startswith("r_phase_Pred_Activity"):
                try:
                    preds.append(float(mol.GetProp(prop)))
                except:
                    pass
        if preds:
            predictedActivity = sum(preds) / len(preds)   # mean
        else:
            predictedActivity = None

        data.append({
            "name": name,
            "predicted_activity": predictedActivity
        })

    df = pd.DataFrame(data)
    df.to_csv(outputCsv, index=False)



if __name__ == "__main__":

    sdfFile = sys.argv[1]
    csvFile = sys.argv[2]

    sdfToCsv(sdfFile, csvFile)