#!/usr/bin/env python3
"""
Convert CSV with SMILES and activity (pIC50) to 3D SDF for Phase QSAR.
Usage:
    python csvToSdf.py input.csv output.sdf
"""

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem

def csvToSdf(inputCsv, outputSdf, smilesCol="smiles", activityCol="pIC50", useActivity=True):
    df = pd.read_csv(inputCsv)

    # Add RDKit molecules
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol=smilesCol, molCol='ROMol', includeFingerprints=False)

    writer = Chem.SDWriter(outputSdf)
    for idx, row in df.iterrows():
        mol = row['ROMol']
        if mol is None:
            print(f"Skipping invalid SMILES at row {idx}: {row[smilesCol]}")
            continue

        # Add hydrogens
        molH = Chem.AddHs(mol)

        # Generate 3D coordinates
        try:
            AllChem.EmbedMolecule(molH, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(molH)
        except Exception as e:
            print(f"Failed to generate 3D for row {idx} ({row[smilesCol]}): {e}")
            continue

        molH.SetProp("_Name", str(row.get("name", f"mol_{idx}")))

        # Set activity property
        if useActivity and activityCol in row and pd.notna(row[activityCol]):
            molH.SetProp(activityCol, str(row[activityCol]))

        # Write molecule to SDF
        writer.write(molH)

    writer.close()
    print(f"Saved 3D SDF to: {outputSdf}")


if __name__ == "__main__":

    inputCsv = sys.argv[1]
    outputSdf = sys.argv[2]

    useActivity = True
    if len(sys.argv) > 3:
        useActivity = sys.argv[3].lower() == "true"

    csvToSdf(inputCsv, outputSdf, useActivity=useActivity)