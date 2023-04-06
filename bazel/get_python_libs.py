from __future__ import print_function
import site
import os
try:
    input = raw_input
except NameError:
    pass
python_paths = []
if os.getenv('PYTHONPATH') is not None:
    python_paths = os.getenv('PYTHONPATH').split(':')
try:
    library_paths = site.getsitepackages()
except AttributeError as e:
    from distutils.sysconfig import get_python_lib
    library_paths = [get_python_lib()]
all_paths = set(python_paths + library_paths)
paths = []
for path in all_paths:
    if os.path.isdir(path):
        paths.append(path)
if len(paths) >= 1:
    print(paths[0])