import os
import inspect

# Configuring the build.

platform = ARGUMENTS.get('OS', Platform())

if platform.name in ['posix', 'linux', 'unix']:
    ccflags_debug = ' -Wall -g'
    ccflags_release = ' -O2'
    tools = 'default'

if platform.name in ['win32', 'windows']:
    ccflags_debug = ' -Wall -g'
    ccflags_release = ' -O2'
    tools = 'mingw'

debug = Environment(
    PLATFORM=platform.name,
    CPPPATH=[Dir('include').abspath],
    CCFLAGS=ccflags_debug,
    tools=[tools],
    BUILDTYPE='debug',
    LIBPATH=os.path.join('#lib', platform.name, 'debug')
    )

release = Environment(
    PLATFORM=platform.name,
    CPPPATH=[],
    CCFLAGS=ccflags_release,
    tools=[tools],
    BUILDTYPE='release'
    )

libBioLCCC=SConscript(
    os.path.join('src', 'core', 'SConscript'),
    exports = {'env':debug},
    variant_dir=os.path.join(
        'build', platform.name, debug['BUILDTYPE'], 'core'), 
    duplicate=True
    )

theorchromo_app=SConscript(
    os.path.join('src', 'apps', 'SConscript'),
    exports={'env':debug},
    variant_dir=os.path.join(
        'build', platform.name, debug['BUILDTYPE'], 'apps'), 
    duplicate=True,
    )

Requires(theorchromo_app, libBioLCCC)
