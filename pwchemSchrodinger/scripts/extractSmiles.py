import sys
import csv
from rdkit import Chem


def main():

    txtFile = sys.argv[1]
    outputCsv = sys.argv[2]

    allRows = []

    with open(txtFile, "r") as f:
        sdfFiles = [line.strip() for line in f if line.strip()]

    for sdfFile in sdfFiles:

        supplier = Chem.SDMolSupplier(sdfFile)

        for mol in supplier:
            if mol is None:
                continue

            try:
                mol = Chem.RemoveHs(mol)
                Chem.SanitizeMol(mol)
            except:
                continue

            name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
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