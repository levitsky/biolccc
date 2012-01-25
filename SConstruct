import os
import inspect
import distutils.sysconfig

# Configuring the build.
#========================

Help("""
libBioLCCC and pyBioLCCC build script. 
Usage: scons -Y <path_to_repository> [OPTIONS] [TARGETS]
Note that building in a source directory is strictly prohibited.

Options:
    buildtype=<buildtype>    Build type. Possible values are release and debug.

Targets:
    all                      All targets.

    libBioLCCC_shared        Shared BioLCCC library. Compiled by default.
    pyBioLCCC                Preparations for pyBioLCCC package. The building
                             process is finished by invoking setup.py.
                             Compiled by default.
    libBioLCCC_static        Static BioLCCC library. 
    libgtest_static          Static Google Test library. 
    tests                    A test suite for libBioLCCC.
    test_pyBioLCCC           A test suite for pyBioLCCC.
    examples                 Examples for libBioLCCC.
    doc                      The documentation.
""")

# Get the mode flag from the command line.
# Default to 'release' if the user didn't specify.
buildtype = ARGUMENTS.get('buildtype', 'release')
if buildtype not in ['debug', 'release']:
   print "Error: expected 'debug' or 'release', found: %s" % (buildtype,)
   Exit(1)

# It is strictly prohibited to build the application in the source directory.
build_in_repository = not bool(GetOption('repository'))
if not build_in_repository:
    for dir in GetOption('repository'):
        if Dir(dir).abspath == Dir('.').abspath:
            build_in_repository = True
if build_in_repository:
    print 'Error:'
    print 'Avoid building the application in the source directory.'
    print 'Print \'scons -Y source_directory\' in the build directory.'
    Exit(1)

VariantDir('.', GetOption('repository'), duplicate=True)

# Setting the platform specific options.
platform = ARGUMENTS.get('OS', Platform())
version = open(str(File('./VERSION').srcnode())).readline().strip()

ccflags = ' -DVERSION=\"%s\"' % version
if platform.name in ['posix', 'linux', 'unix']:
    if buildtype=='release':
        ccflags += ' -O2'
    else:
        ccflags += ' -Wall -g'
    tools = 'default'

if platform.name in ['win32', 'windows']:
    if buildtype=='release':
        ccflags += ' -O2'
    else:
        ccflags += ' -Wall -g'
    tools = 'mingw'

env = Environment(
    PLATFORM=platform.name,
    CPPPATH=[Dir('include').abspath, distutils.sysconfig.get_python_inc()],
    #CPPPATH=[Dir('include').abspath, "/usr/include/python3.1"],
    CCFLAGS=ccflags,
    tools=[tools],
    BUILDTYPE=buildtype,
    LIBPATH=[os.path.join('#lib', 'static', platform.name, buildtype),
             os.path.join('#lib', 'shared', platform.name, buildtype)],
    ROOTBUILDDIR=Dir('.').abspath,
    )

# Building targets.
#===================

# Shared BioLCCC library.
#-------------------------
libBioLCCC_shared = SConscript(
    os.path.join('src', 'core', 'SConscript'),
    exports = {'env':env, 'libtype':'shared'},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'core', 'shared'), 
    duplicate=True
    )
Alias('libBioLCCC_shared', libBioLCCC_shared)

# Static BioLCCC library.
#-------------------------
libBioLCCC_static = SConscript(
    os.path.join('src', 'core', 'SConscript'),
    exports = {'env':env, 'libtype':'static'},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'core', 'static'), 
    duplicate=True
    )
Alias('libBioLCCC_static', libBioLCCC_static)

# Google test library.
#----------------------
libgtest_static = SConscript(
    os.path.join('src', 'gtest', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'gtest'), 
    duplicate=True,
    )
Alias('libgtest_static', libgtest_static)

# Test suite.
#--------------
tests = SConscript(
    os.path.join('src', 'apps', 'SConscript'),
    exports={'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'apps'), 
    duplicate=True,
    )

Requires(tests, libBioLCCC_static)
Requires(tests, libgtest_static)
Alias('tests', tests)

# pyBioLCCC package.
#-------------------
biolccc_py = SConscript(
    os.path.join('src', 'pyteomics', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'pyteomics'), 
    duplicate=True,
    )

# Copying source files required for the python source package.
build_dir_path = os.path.join(
    'build', platform.name, env['BUILDTYPE'], 'pyteomics')
src_dir_path = os.path.join(
        os.path.join(Dir('#.').abspath, 'src', 'pyteomics'))

env.AddPostAction(biolccc_py, 
    Copy(src_dir_path, os.path.join(build_dir_path, 'biolccc.py')))
env.AddPostAction(biolccc_py, 
    Copy(src_dir_path,  os.path.join(build_dir_path, 'biolccc_wrap.cc')))
env.AddPostAction(biolccc_py, 
    'python ' + os.path.join(src_dir_path, 'post_swig.py') + ' ' + version)

# Copying the remaining files for pyBioLCCC.
Depends(biolccc_py, Glob('src/pyteomics/*'))
Depends(biolccc_py, 'setup.py')
Depends(biolccc_py, 'MANIFEST.in')
# Copying the documentation to the build dir.
Depends(biolccc_py, 'VERSION')
Depends(biolccc_py, 'README')
Alias('biolccc_py', biolccc_py)

# A test suite for pyBioLCCC.
#----------------------------
test_pyBioLCCC = env.Install(
    os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings'),
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'test_pyBioLCCC.py'))
Depends(test_pyBioLCCC, 
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'test_pyBioLCCC.py'))
Alias('test_pyBioLCCC', test_pyBioLCCC)

# Examples.
#----------
examples = SConscript(
    os.path.join('src', 'examples', 'SConscript'),
    exports={'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'examples'), 
    duplicate=True,
    )

Requires(examples, libBioLCCC_static)
Alias('examples', examples)

# Doxygen documentation.
#-----------------------
doc_doxygen = env.Command('doc_doxygen', 'doc/Doxyfile', 
    [Mkdir('doc'), 'doxygen $SOURCE'])
Depends(doc_doxygen, 'doc/Doxyfile')
# Source code needs to be copied.
Depends(doc_doxygen, libBioLCCC_shared)

# Sphinx documentation.
#----------------------
doc_sphinx = env.Command('doc_sphinx', '', 
    'mkdir ./doc/sphinx/source/examples; '
    'cp src/examples/*.py ./doc/sphinx/source/examples; '
    'cp src/examples/*.cpp ./doc/sphinx/source/examples; '
    'cd ./doc/sphinx; make html; cd ../')
Depends(doc_sphinx, Glob('doc/sphinx/*'))
Depends(doc_sphinx, Glob('doc/sphinx/source/*'))
Depends(doc_sphinx, Glob('doc/sphinx/source/_static/*'))
Depends(doc_sphinx, Glob('doc/sphinx/source/_sphinxext/*.py'))
Depends(doc_sphinx, Glob('doc/sphinx/source/_templates/*'))
Depends(doc_sphinx, 'VERSION')
Depends(doc_sphinx, 'README')
Depends(doc_sphinx, 'INSTALL')
Depends(doc_sphinx, 'CHANGELOG')
Depends(doc_sphinx, Glob('src/examples/*.py'))
Depends(doc_sphinx, Glob('src/examples/*.cpp'))

# Complete documentation.
#------------------------
docs = env.Command('docs', '',
    #[''])
    ['mv doc/sphinx/build/html/* doc/',
     'mkdir doc/build',
     'mkdir doc/API',
     'mv doc/doxygen/html/* doc/API',
     'rm -r doc/doxygen',
     'rm -r doc/sphinx'])
Requires('docs', doc_doxygen)
Requires('docs', doc_sphinx)
# For some funny reason name 'doc' doesn't work, so we need to use an alias.
Alias('doc', docs)

# Final configuration of the build.
#===================================
env.Default([libBioLCCC_shared, biolccc_py])
Alias('all', 
    [libBioLCCC_shared, libBioLCCC_shared, libgtest_static, biolccc_py,
     tests, test_pyBioLCCC, examples, docs])

