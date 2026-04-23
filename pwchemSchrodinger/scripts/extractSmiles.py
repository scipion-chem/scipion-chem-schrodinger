import sys
import csv
import os
from rdkit import Chem

def loadMolecules(filePath):
    ext = os.path.splitext(filePath)[1].lower()

    if ext in [".sdf", ".sd"]:
        return Chem.SDMolSupplier(filePath)

    elif ext == ".mol2":
        mols = []
        with open(filePath, "r") as f:
            block = []
            for line in f:
                if line.startswith("@<TRIPOS>MOLECULE") and block:
                    mol = Chem.MolFromMol2Block("".join(block), sanitize=False)
                    if mol:
                        mols.append(mol)
                    block = []
                block.append(line)

            if block:
                mol = Chem.MolFromMol2Block("".join(block), sanitize=False)
                if mol:
                    mols.append(mol)

        return mols

    else:
        print(f"⚠️ Unsupported format: {filePath}")
        return []

def main():

    txtFile = sys.argv[1]
    outputCsv = sys.argv[2]

    allRows = []

    with open(txtFile, "r") as f:
        molFiles = [line.strip() for line in f if line.strip()]

    for molFile in molFiles:

        supplier = loadMolecules(molFile)

        for mol in supplier:
            if mol is None:
                continue

            try:
                Chem.SanitizeMol(mol)
                mol = Chem.RemoveHs(mol)
            except:
                continue

            name = mol.GetProp("_Name") if mol.HasProp("_Name") else os.path.basename(molFile)
            smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

            allRows.append({
                "name": name,
                "smiles": smiles
            })

    with open(outputCsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "smiles"])
        writer.writeheader()
        writer.writerows(allRows)

    print(f"Saved SMILES to: {outputCsv}")


if __name__ == "__main__":
    main()