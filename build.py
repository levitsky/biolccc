from __future__ import print_function

import os, sys, shutil, argparse, subprocess

pjoin = os.path.join

parser = argparse.ArgumentParser(
    description='Build pyteomics.biolccc and the documentation.',
    epilog='',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-b', type=str, nargs=1, metavar='BUILD_PATH',
                    help='a path to the build directory',
                    required=True)
parser.add_argument(
    'tasks', type=str, metavar='TASKS', nargs='+',
    help=
        'build tasks. Possible values:\n'
        'pyteomics.biolccc - generate source files for pyteomics.biolccc\n'
        'doc - generate the documentation for libBioLCCC and\n'
        '      pyteomics.biolccc')

def _mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass

def _configure():
    global parser
    conf = {}
    args = parser.parse_args()
    if args.b == None:
        print('Please specify the build directory')
        sys.exit(1)
    conf['BUILD_PATH'] = os.path.abspath(args.b[0])
    conf['SRC_PATH'] = os.path.abspath(os.getcwd())
    conf['SWIG'] = os.environ.get('BIOLCCC_SWIG', 'swig')
    if conf['BUILD_PATH'] == conf['SRC_PATH']:
        print('Building in the source directory is forbidden.')
        print('Please specify another build directory.')
        sys.exit(1)
    _mkdir(conf['BUILD_PATH'])
    conf['VERSION'] = open('./VERSION').readline().strip()
    conf['TASKS'] = args.tasks

    return conf

def _swig(conf):
    _mkdir(pjoin(conf['BUILD_PATH'], 'pyteomics'))
    shutil.copy(pjoin('.', 'src', 'pyteomics', 'biolccc.i'),
                pjoin(conf['BUILD_PATH'], 'pyteomics'))
    print('Generate Python wrappers... ')
    subprocess.call('cd %s; %s -python -c++ -I%s biolccc.i' % (
            conf['SWIG'],
            pjoin(conf['BUILD_PATH'], 'pyteomics'),
            pjoin(conf['SRC_PATH'], 'include')),
        shell=True)
    print('Done!')

def _post_swig(conf):
    '''Inherit MutableMapping explicitly since SWIG does not allow to do that
    with Python < 3.0.
    '''
    modified_file = []
    py_module_path = pjoin(conf['BUILD_PATH'], 'pyteomics', 'biolccc.py')
    with open(py_module_path, 'r') as f:
        for line in f:
            modified_line = line
            modified_line = modified_line.replace(
                'class ChemicalGroup(_object)',
                'class ChemicalGroup(collections.MutableMapping, _object)')
            modified_line = modified_line.replace(
                'class ChemicalBasis(_object)',
                'class ChemicalBasis(collections.MutableMapping, _object)')
            modified_line = modified_line.replace(
                'class ChromoConditions(_object)',
                'class ChromoConditions(collections.MutableMapping, _object)')
            modified_line = modified_line.replace(
                'class GradientPoint(_object)',
                'class GradientPoint(collections.MutableMapping, _object)')
            modified_file.append(modified_line)

        newlines = f.newlines

    if not 'import collections' in modified_file[0]:
        modified_file.insert(0, 'VERSION = \"%s\"' % conf['VERSION']
                                + (newlines or '\n'))
        modified_file.insert(0, 'import collections' + (newlines or '\n'))

    with open(py_module_path, 'w') as f:
        f.writelines(modified_file)

def _configure_distutils(conf):
    print('Prepare sources for pyteomics.biolccc...')

    for filename in ['setup.py', 'MANIFEST.in', 'VERSION', 'README.rst']:
        shutil.copy(pjoin(conf['SRC_PATH'], filename),
                    pjoin(conf['BUILD_PATH'], filename))
    shutil.copy(pjoin(conf['SRC_PATH'], 'src', 'pyteomics', '__init__.py'),
                pjoin(conf['BUILD_PATH'], 'pyteomics', '__init__.py'))

    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'include')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'include'),
                        pjoin(conf['BUILD_PATH'], 'include'))
    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'src')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'src'),
                        pjoin(conf['BUILD_PATH'], 'src'))
    print('Done!')

def _sphinx_doc(conf):
    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'doc')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'doc'),
                        pjoin(conf['BUILD_PATH'], 'doc'),
                        symlinks=True)

    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'doc', 'sphinx',
                               'source', 'examples')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'src', 'examples'),
                        pjoin(conf['BUILD_PATH'], 'doc', 'sphinx',
                              'source', 'examples'))

    for filename in ['VERSION', 'README.rst', 'INSTALL', 'CHANGELOG']:
        shutil.copy(filename,
                    pjoin(conf['BUILD_PATH'], filename))

    subprocess.call('cd %s; make html' % (
            pjoin(conf['BUILD_PATH'], 'doc', 'sphinx'),),
        shell=True)

def _doxygen_doc(conf):
    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'doc')):
        _mkdir(pjoin(conf['BUILD_PATH'], 'doc'))

    if not os.path.isfile(pjoin(conf['BUILD_PATH'], 'doc', 'Doxyfile')):
        shutil.copy(pjoin(conf['SRC_PATH'], 'doc', 'Doxyfile'),
                    pjoin(conf['BUILD_PATH'], 'doc', 'Doxyfile'))

    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'include')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'include'),
                        pjoin(conf['BUILD_PATH'], 'include'))

    if not os.path.isdir(pjoin(conf['BUILD_PATH'], 'src')):
        shutil.copytree(pjoin(conf['SRC_PATH'], 'src'),
                        pjoin(conf['BUILD_PATH'], 'src'))

    subprocess.call('cd %s; doxygen %s' % (
        conf['BUILD_PATH'], pjoin('doc', 'Doxyfile'),), shell=True)

def _doc(conf):
    _sphinx_doc(conf)
    _doxygen_doc(conf)

    if os.path.isdir(pjoin(conf['BUILD_PATH'], 'doc', 'finished')):
        shutil.rmtree(pjoin(conf['BUILD_PATH'], 'doc', 'finished'))

    shutil.move(
        pjoin(conf['BUILD_PATH'], 'doc', 'sphinx', 'build', 'html'),
        pjoin(conf['BUILD_PATH'], 'doc', 'finished'))

    shutil.move(
        pjoin(conf['BUILD_PATH'], 'doc', 'doxygen', 'html'),
        pjoin(conf['BUILD_PATH'], 'doc', 'finished', 'API'))

    shutil.make_archive(
        pjoin(conf['BUILD_PATH'], 'doc', 'doc'),
        'zip',
        pjoin(conf['BUILD_PATH'], 'doc', 'finished'))

conf = _configure()
if 'pyteomics.biolccc' in conf['TASKS']:
    _swig(conf)
    _post_swig(conf)
    _configure_distutils(conf)

if 'doc' in conf['TASKS']:
    _doc(conf)

