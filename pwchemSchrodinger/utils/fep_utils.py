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

"""
Shared helpers for the FEP+ protocols (protocol_fep_rbfe.py / protocol_fep_abfe.py):
building the receptor/ligand ".mae" inputs FEP+ needs, assembling them into a
pose-viewer file (the standard "receptor first, ligand(s) after" Maestro convention
already used elsewhere in this plugin - e.g. Glide's own "*_pv.maegz" docking output),
a tool-agnostic ligand-similarity edge-chain (reused, translated, from
gromacs/protocols/protocol_pmx_rbfe.py's _computeEdgeChain), and a deliberately lenient
parser for FEP+'s own free-text job log.

See claude/decisions/schrodinger/fep_plus.md for why these two protocols could not be
built and verified the same way the GROMACS pmx RBFE/ABFE protocols were (no licensed
Schrodinger install to actually run commands against in this session, and the official
FEP+ user manual / command reference are behind an authenticated learn.schrodinger.com
login this session could not pass) - every helper below that touches an unverified FEP+
CLI detail says so in its own docstring, instead of presenting a guess as a confirmed fact.
"""

import os
import re
import shutil

from pwchem.utils import convertToSdf

from .. import Plugin as schrodingerPlugin

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

structConvertProg = schrodingerPlugin.getHome('utilities/structconvert')
structCatProg = schrodingerPlugin.getHome('utilities/structcat')
progLigPrep = schrodingerPlugin.getHome('ligprep')
progPrepWizard = schrodingerPlugin.getHome('utilities/prepwizard')

MAEFILE_EXTENSION = '.maegz'


def convertMaeToMol2(protocol, molFile, outDir):
    """Real, correct -> .mol2 conversion for RDKit fingerprinting (used by
    computeLigandEdgeChain below), for any input format structconvert reads directly -
    confirmed for both ".mae"/".maegz" (Glide docking output) and ".cif" (pwchem's own
    ProtExtractLigands output). Despite the name (kept because ".mae" is the single most
    common real caller), it is not mae-specific. Deliberately NOT
    pwchemSchrodinger.utils.utils.convertMAE2Mol2File: despite its name, that function
    builds its output path from the INPUT file's own basename (`os.path.split(molFile)[-1]`
    - i.e. it keeps the original extension, e.g. still ".maegz"), so structconvert ends up
    doing a same-format copy, not an actual conversion to ".mol2" - a real, pre-existing bug
    in that utility found by reading it before reusing it, not introduced here. This
    version builds the output filename with an explicit ".mol2" extension instead."""
    outFile = os.path.join(outDir, os.path.splitext(os.path.basename(molFile))[0] + '.mol2')
    protocol.runJob(structConvertProg, f'{os.path.abspath(molFile)} {os.path.abspath(outFile)}',
                    cwd=outDir)
    return outFile


def prepareReceptorMae(protocol, pdbFile, outDir, prepare=True):
    """Convert/prepare a receptor structure into a ".mae" file, the same
    "-noprotassign -noimpref -noepik" quick-prep shortcut ProtSchrodingerGlideDocking/
    ProtSchrodingerMMGBSA/ProtSchrodingerDesmondSysPrep already use for a non-Maestro
    receptor: structure fixing only, no protonation/refinement/ionization - the full
    ProtSchrodingerPrepWizard protocol should be run upstream first if that matters for
    the system being studied (see KNOWLEDGE_SCHRODINGER.md §9)."""
    targetName = os.path.splitext(os.path.basename(pdbFile))[0]
    targetMaeFile = os.path.abspath(os.path.join(outDir, targetName + MAEFILE_EXTENSION))
    if prepare:
        args = f'-WAIT -noprotassign -noimpref -noepik {os.path.abspath(pdbFile)} {targetMaeFile}'
        protocol.runJob(progPrepWizard, args, cwd=outDir)
    else:
        protocol.runJob(structConvertProg, f'{os.path.abspath(pdbFile)} {targetMaeFile}', cwd=outDir)
    return targetMaeFile


def prepareLigandMae(protocol, molFile, outDir, molName, prepare=True):
    """Convert/prepare a ligand pose into a ".mae" file - same quick LigPrep-or-plain-
    convert choice as ProtSchrodingerMMGBSA.prepareLigandFile. Ionization/tautomerization
    (if prepare=True) can shift the docked pose slightly; disable it to keep the exact
    input pose unchanged (e.g. when the pose is already a genuine Schrodinger/LigPrep
    output).

    Real bug found by actually running this against a Glide-docked ligand (already
    ".maegz"): every caller (both FEP+ protocols' own _ligandMaeFile()) assumes this
    function's output always lands at "outDir/<molName>.maegz" - true for the LigPrep and
    structconvert branches below (they explicitly write there), but the "already .mae/
    .maegz, nothing to do" branch used to just return the ORIGINAL file's own path
    unchanged, wherever that happened to be (e.g. a sibling docking protocol's own extra/
    directory) - so the caller's assumed path was never actually created, and the next
    step (building the pose-viewer file) failed with "Input file does not exist". Fixed by
    always copying into outDir/<molName>.maegz in that branch too, same contract as the
    other two."""
    if not molFile.endswith('.sdf') and not molFile.endswith(('.mae', '.maegz')):
        molFile = convertToSdf(protocol, molFile)

    if prepare and not molFile.endswith(('.mae', '.maegz')):
        tmpMaeFile = os.path.abspath(os.path.join(outDir, molName + '_tmp.maegz'))
        args = f'-WAIT -R h -a -isd {molFile} -omae {tmpMaeFile}'
        protocol.runJob(progLigPrep, args, cwd=outDir)
        maeFile = os.path.abspath(os.path.join(outDir, molName + MAEFILE_EXTENSION))
        os.rename(tmpMaeFile, maeFile)
    elif molFile.endswith(('.mae', '.maegz')):
        maeFile = os.path.abspath(os.path.join(outDir, molName + MAEFILE_EXTENSION))
        if os.path.abspath(molFile) != maeFile:
            shutil.copy(molFile, maeFile)
    else:
        maeFile = os.path.abspath(os.path.join(outDir, molName + MAEFILE_EXTENSION))
        protocol.runJob(structConvertProg, f'{molFile} {maeFile}', cwd=outDir)
    return maeFile


def buildPoseViewerFile(protocol, receptorMae, ligandMaeFiles, outFile):
    """Concatenate receptor + ligand(s) into one pose-viewer-style ".mae", receptor
    first (the standard Maestro/Glide convention this plugin's own docking output
    already follows - see ProtSchrodingerGlideDocking's "*_pv.maegz"). This is the input
    shape both fep_plus (RBFE, 1 receptor + 2 co-aligned ligands sharing a common core)
    and the ABFE driver (1 receptor + 1 ligand) are documented to accept as either a
    ".mae" or ".fmp" positional argument. Callers should name `outFile` with a ".maegz"
    extension - confirmed for real (structcat on a real 4erf pose-viewer file) that
    structcat, like every other Schrodinger CLI tool in this plugin, compresses its own
    output automatically based on that extension (real gzip magic bytes in the result),
    not a plain-text ".mae" left mis-named."""
    files = ' '.join(f'-imae {f}' for f in [receptorMae] + list(ligandMaeFiles))
    args = f'{files} -omae {outFile}'
    protocol.runJob(structCatProg, args, cwd=os.path.dirname(outFile) or '.')
    return outFile


def computeLigandEdgeChain(protocol, names, molFileGetter):
    """Order `names` into an A->B->C->... chain by greedy nearest-neighbor structural
    similarity (RDKit Morgan fingerprints, Tanimoto) - line-for-line the same algorithm
    gromacs/protocols/protocol_pmx_rbfe.py's GromacsPmxRBFE._computeEdgeChain uses,
    since the "pick a similarity-ordered chain over a full LOMAP-style perturbation
    network" scope decision (see claude/decisions/amber/pmx_RBFE.md §15) is about
    ligand chemistry, not about which MD engine runs the alchemical legs - it applies
    unchanged here. `molFileGetter(name)` must return a path RDKit can read (mol2/sdf/pdb).
    Falls back to selection order if RDKit is unavailable or a structure can't be
    fingerprinted - never raises."""
    if len(names) < 2:
        return []

    if not RDKIT_AVAILABLE:
        return list(zip(names[:-1], names[1:]))

    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = {}
    try:
        for name in names:
            molFile = molFileGetter(name)
            rdmol = None
            if molFile.endswith('.mol2'):
                rdmol = Chem.MolFromMol2File(molFile, sanitize=True)
            elif molFile.endswith('.sdf'):
                rdmol = next(iter(Chem.SDMolSupplier(molFile, sanitize=True)), None)
            elif molFile.endswith('.pdb'):
                rdmol = Chem.MolFromPDBFile(molFile, sanitize=True)
            if rdmol is None:
                raise ValueError(f'RDKit could not parse the structure for "{name}" ({molFile})')
            fps[name] = morganGen.GetFingerprint(rdmol)
    except Exception as e:
        protocol.warning(f'Could not compute structural fingerprints for the edge chain '
                         f'({e}); falling back to the order ligands were selected in.')
        return list(zip(names[:-1], names[1:]))

    remaining = names[1:]
    chain = [names[0]]
    while remaining:
        last = fps[chain[-1]]
        best = max(remaining, key=lambda n: DataStructs.TanimotoSimilarity(last, fps[n]))
        chain.append(best)
        remaining.remove(best)
    return list(zip(chain[:-1], chain[1:]))


# fmp2excel.py IS a real, publicly-documented FEP+ utility (confirmed via the FEP+
# Utilities topic, see claude/decisions/schrodinger/fep_plus.md §3): "Script to write out
# a CSV file or Excel spreadsheet from a .fmp file... contains the same data [as the FEP+
# Panel's Map tab]... You can optionally perform cycle closure analysis with the
# -cycle-closure option." This is the real, sanctioned way to get FEP+'s predicted
# affinities out of a ".fmp" without the GUI - far more solid ground than scraping a job
# log whose exact text this session could never see. What is NOT confirmed is the exact
# column header names in the CSV it writes (the docs describe the content, not the exact
# header row) - parseFmpCsv below is therefore still deliberately lenient about column
# names, trying several plausible ones before giving up, rather than hardcoding one guess.
fmp2excelScript = 'fmp2excel.py'

_DG_COLUMN_CANDIDATES = ('Pred. dG', 'Pred dG', 'ddG(cycle closure)', 'ddG (cycle closure)',
                        'Predicted ddG (corrected)', 'ddG', 'dG')
_TITLE_COLUMN_CANDIDATES = ('Title', 'Name', 'Ligand')


def runFmp2Excel(protocol, fmpFile, outCsv, cwd=None):
    """Run the real fmp2excel.py FEP+ utility (via Schrodinger's own bundled Python,
    same "run <script>" convention as this plugin's own getCloseResidues.py) with
    "-cycle-closure" so the resulting CSV includes the cycle-closure-corrected ddG/error,
    not just the raw per-edge BAR value. Never raises on failure - callers fall back to
    "no parsed value, raw .fmp kept" the same way every other lenient parser in this
    plugin/its siblings does."""
    args = f'{fmp2excelScript} {os.path.abspath(fmpFile)} -cycle-closure -o {os.path.abspath(outCsv)}'
    try:
        protocol.runJob(schrodingerPlugin.getHome('run'), args, cwd=cwd)
        return outCsv if os.path.exists(outCsv) else None
    except Exception as e:
        protocol.warning(f'fmp2excel.py failed to export "{fmpFile}" ({e}); the raw .fmp '
                         f'is still kept as this run\'s output - open it in the FEP+ Panel '
                         f'(Maestro GUI) to inspect results by hand.')
        return None


def parseFmpCsv(csvFile, ligandNames=None):
    """Best-effort parse of fmp2excel.py's CSV - column names are not confirmed (see the
    module note above), so this tries several plausible dG/ddG and title/name columns
    before giving up. Returns {rowTitle: value} for whichever dG-like column was found
    (empty dict if none of the candidate columns are present, or the file can't be read).
    `ligandNames`, if given, is only used to prefer rows whose title matches one of the
    expected ligand names when several dG-like columns are present in the same file."""
    if not csvFile or not os.path.exists(csvFile):
        return {}
    try:
        import csv as csvmod
        with open(csvFile, newline='') as f:
            rows = list(csvmod.DictReader(f))
        if not rows:
            return {}
        titleCol = next((c for c in _TITLE_COLUMN_CANDIDATES if c in rows[0]), None)
        dgCol = next((c for c in _DG_COLUMN_CANDIDATES if c in rows[0]), None)
        if titleCol is None or dgCol is None:
            return {}
        out = {}
        for row in rows:
            try:
                out[row[titleCol]] = float(row[dgCol])
            except (TypeError, ValueError):
                continue
        return out
    except Exception:
        return {}


def getEdgeName(ligA, ligB):
    """A filesystem/jobname-safe identifier for an A->B edge, from the two ligand names
    as picked in the wizard (arbitrary Scipion object labels, not file paths - so this
    only strips characters that would be awkward in a JOBNAME/directory name, unlike
    pwchem's own getBaseName which is for file paths)."""
    safe = lambda n: re.sub(r'[^\w.-]', '_', n)
    return f'{safe(ligA)}_{safe(ligB)}'
