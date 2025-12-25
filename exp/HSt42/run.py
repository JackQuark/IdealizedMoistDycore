import sys, os
Pjoin = os.path.join

_dycore_pkg_path = Pjoin(os.environ.get('Dycore_BASE'), 'src', 'extra', 'python')
sys.path.append(_dycore_pkg_path)
from dycore import *


print(Experiment)