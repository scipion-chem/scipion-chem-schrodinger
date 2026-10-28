'''Split every complex of an IFD output file into its receptor and its docked ligand.
IFD always places the docked ligand in chain Z, residue 999: selecting it explicitly avoids guessing ligands by
size, which takes small receptor fragments (e.g. a peptide segment after a chain break) as the ligand.
Usage: $SCHRODINGER/run splitIFDComplexes.py <ifd-out.maegz> <outDir>
Writes ifd_<title>_<n>_receptor.maegz and ifd_<title>_<n>_ligand.maegz, n being the pose number of each title
in file order.'''
import os
import sys

from schrodinger import structure
from schrodinger.structutils.analyze import evaluate_asl

LIGAND_ASL = 'chain.name Z and res.num 999'


def getSplitBase(outDir, title, poseIdx):
    return os.path.join(outDir, f'ifd_{title}_{poseIdx}')


def splitComplexes(complexFile, outDir):
    poseCounts = {}
    for st in structure.StructureReader(complexFile):
        poseCounts[st.title] = poseCounts.get(st.title, 0) + 1
        ligAtoms = evaluate_asl(st, LIGAND_ASL)
        if not ligAtoms:
            sys.exit(f'No docked ligand ({LIGAND_ASL}) found in the IFD complex of {st.title}')

        base = getSplitBase(outDir, st.title, poseCounts[st.title])
        st.extract(ligAtoms, copy_props=True).write(base + '_ligand.maegz')
        st.deleteAtoms(ligAtoms)
        st.write(base + '_receptor.maegz')


if __name__ == '__main__':
    splitComplexes(sys.argv[1], sys.argv[2])
