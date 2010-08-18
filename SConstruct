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
    pyBioLCCC                pyBioLCCC package. Compiled by default.
    libBioLCCC_static        Static BioLCCC library. 
    libgtest_static          Static Google Test library. 
    test                     A test suite.
    doc                      A documentation for libBioLCCC and pyBioLCCC.
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

if platform.name in ['posix', 'linux', 'unix']:
    if buildtype=='release':
        ccflags = ' -O2'
    else:
        ccflags = ' -Wall -g'
    tools = 'default'

if platform.name in ['win32', 'windows']:
    if buildtype=='release':
        ccflags = ' -O2'
    else:
        ccflags = ' -Wall -g'
    tools = 'mingw'

env = Environment(
    PLATFORM=platform.name,
    CPPPATH=[Dir('include').abspath, distutils.sysconfig.get_python_inc()],
    CCFLAGS=ccflags,
    tools=[tools],
    BUILDTYPE=buildtype,
    LIBPATH=os.path.join('#lib', 'static', platform.name, buildtype),
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
#---------------------
pyBioLCCC = SConscript(
    os.path.join('src', 'bindings', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'bindings'), 
    duplicate=True,
    )

# Copying source files required for the python source package.
env.AddPostAction(pyBioLCCC, Copy(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'pyBioLCCC.py'),
    os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
        'pyBioLCCC.py')))
env.AddPostAction(pyBioLCCC, Copy(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'pyBioLCCC_wrap.cc'),
    os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
        'pyBioLCCC_wrap.cc')))
env.AddPostAction(pyBioLCCC, Touch(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', '__init__.py')))

# Copying the documentation to the build dir.
Depends(pyBioLCCC, 'setup.py')
Depends(pyBioLCCC, 'VERSION')
Depends(pyBioLCCC, 'MANIFEST.in')
Depends(pyBioLCCC, 'pyBioLCCC.README')
Alias('pyBioLCCC', pyBioLCCC)

# Doxygen documentation.
#------------------------
doc = env.Command('doc', 'Doxyfile', 'doxygen $SOURCE')
Depends(doc, 'Doxyfile')
# Source code needs to be copied.
Depends(doc, libBioLCCC_shared)

# Final configuration of the build.
#===================================
env.Default([libBioLCCC_shared, pyBioLCCC])
Alias('all', 
    [libBioLCCC_shared, libBioLCCC_shared, libgtest_static,
     tests, pyBioLCCC, doc])
