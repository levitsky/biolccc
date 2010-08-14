import os
import inspect
import distutils.sysconfig

# Configuring the build.

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
libBioLCCC_static = SConscript(
    os.path.join('src', 'core', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'core'), 
    duplicate=True
    )

pyBioLCCC_so = SConscript(
    os.path.join('src', 'bindings', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'bindings'), 
    duplicate=True,
    )

libgtest_static = SConscript(
    os.path.join('src', 'gtest', 'SConscript'),
    exports = {'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'gtest'), 
    duplicate=True,
    )

tests_app=SConscript(
    os.path.join('src', 'apps', 'SConscript'),
    exports={'env':env},
    variant_dir=os.path.join(
        'build', platform.name, env['BUILDTYPE'], 'apps'), 
    duplicate=True,
    )
Requires(tests_app, libBioLCCC_static)
Requires(tests_app, libgtest_static)

# Copying source files required for the python source package.
env.AddPostAction(pyBioLCCC_so, Copy(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'pyBioLCCC.py'),
    os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
        'pyBioLCCC.py')))
env.AddPostAction(pyBioLCCC_so, Copy(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', 'pyBioLCCC_wrap.cc'),
    os.path.join('build', platform.name, env['BUILDTYPE'], 'bindings',
        'pyBioLCCC_wrap.cc')))
env.AddPostAction(pyBioLCCC_so, Touch(
    os.path.join(Dir('#.').abspath, 'src', 'bindings', '__init__.py')))

# Copying the documentation to the build dir.
Depends(pyBioLCCC_so, 'setup.py')
Depends(pyBioLCCC_so, 'VERSION')
Depends(pyBioLCCC_so, 'MANIFEST.in')
Depends(pyBioLCCC_so, 'pyBioLCCC.README')

