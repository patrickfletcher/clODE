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
    import sysconfig
    library_paths = [sysconfig.get_path('purelib')]
all_paths = set(python_paths + library_paths)
paths = []
for path in all_paths:
    if os.path.isdir(path):
        paths.append(path)
if len(paths) >= 1:
    print(paths[0])