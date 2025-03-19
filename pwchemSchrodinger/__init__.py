# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

"""
This package contains protocols interfacing Schrodinger's Maestro
"""

# Useful reference: http://www.bi.cs.titech.ac.jp/mga_glide/xglide_mga.py

import glob
import os, fnmatch
from subprocess import Popen, run
import pwem
import pyworkflow.utils as pwutils
from .bibtex import _bibtexStr
from .constants import *

_logo = 'schrodinger_logo.png'
DEFAULT_VERSION = '2024-4'
SCHRODINGER_DIC = {'name': 'schrodinger', 'version': 'SCHRODINGER_VERSION', 'home': 'SCHRODINGER_HOME',
                   'SGL': 'MAESTRO_SGL'}

class Plugin(pwem.Plugin):
    _homeVar = SCHRODINGER_DIC['home']

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(SCHRODINGER_DIC['SGL'], False)
        cls._defineVar(SCHRODINGER_DIC['version'], DEFAULT_VERSION)

        schroVersion = cls.getVar(SCHRODINGER_DIC['version'])
        cls._defineEmVar(SCHRODINGER_DIC['home'], f'Schrodinger{schroVersion}')

    @classmethod
    def getEnviron(cls, schrodingerFirst=True):
        """ Create the needed environment for Schrodinger programs. """
        environ = pwutils.Environ(os.environ)
        pos = pwutils.Environ.BEGIN if schrodingerFirst else pwutils.Environ.END
        environ.update({
            'SCHRODINGER': cls.getHome(''),
            'PATH': cls.getHome(''),
            'LD_LIBRARY_PATH': '{}:{}:{}'.format(cls.getMMshareDir('lib/Linux-x86_64'), cls.getHome('internal/lib'),
                                                 cls.getHome('internal/lib/ssl')),
            'PYTHONPATH': cls.getSitePackages()
        }, position=pos)
        return environ

    @classmethod
    def getMMshareDir(cls, fn):
        fileList = glob.glob(cls.getHome('mmshare*'))
        if len(fileList)==0:
            return None
        else:
            return os.path.join(fileList[0],fn)

    @classmethod
    def getSchrodingerPython(cls):
        return cls.getHome('internal/bin/python3')

    @classmethod
    def getSitePackages(cls):
        fileList = glob.glob(cls.getHome('internal/lib/python3.*'))
        if len(fileList) == 0:
            return None
        else:
            return os.path.join(fileList[0], 'site-packages')

    @classmethod
    def getPluginHome(cls, path=""):
        import pwchemSchrodinger
        fnDir = os.path.split(pwchemSchrodinger.__file__)[0]
        return os.path.join(fnDir,path)

    @classmethod
    def runSchrodingerScript(cls, program, args, cwd=None, popen=False):
        """ Run Schrodinger python command from a given protocol. """
        runProg = cls.getHome('run')
        schArgs = ' {} {}'.format(program, args)
        if popen:
            Popen(runProg + schArgs, env=cls.getEnviron(), cwd=cwd, shell=True)
        else:
            run(runProg + schArgs, env=cls.getEnviron(), cwd=cwd, shell=True)

    @classmethod
    def runSchrodinger(cls, protocol, program, args, cwd=None):
        """ Run Schrodinger command from a given protocol. """
        protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runJobSchrodingerScript(cls, protocol, script, args, cwd=None):
        """ Run Schrodinger command from a given protocol. """
        runProg = cls.getHome('run')
        args = '{} {}'.format(script, args)
        protocol.runJob(runProg, args, env=cls.getEnviron(), cwd=cwd)

    # ---------------------------------- Utils functions  -----------------------
    @staticmethod
    def find(path, pattern):  # This function is analogous to linux find (case of folders)
        paths = []
        for root, dirs, files in os.walk(path):
            for name in dirs:
                if fnmatch.fnmatch(name, pattern):
                    paths.append(os.path.join(root, name))
        return paths