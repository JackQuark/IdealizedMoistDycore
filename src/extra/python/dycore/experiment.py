import os
Pjoin = os.path.join

from dycore.config import Config
from dycore.diagtable import DiagTable
from dycore import Dycore_BASE, Dycore_WORK, Dycore_DATA

class Experiment(object):
    def __init__(self, name, codebase, overwrite=False, workbase=Dycore_WORK, database=Dycore_DATA):
        super(Experiment, self).__init__()
        self.name = name
        self.codebase = codebase
        self.overwrite = overwrite
        
        self.workdir = Pjoin(workbase, 'experiment', self.name)
        self.datadir = Pjoin(database, self.name)
        self.restartdir = Pjoin(self.datadir, 'restart')
       
        self.workdir    = Pjoin(workbase, 'experiment', self.name)
        self.rundir     = Pjoin(self.workdir, 'run')      # temporary area an individual run will be performed
        self.datadir    = Pjoin(database, self.name)      # where run data will be moved to upon completion
        self.restartdir = Pjoin(self.datadir, 'restarts') # where restarts will be stored
        # self.template_dir = Pjoin(_module_directory, 'templates')
        
        self.config = Config()
        self.diagtable = DiagTable()

    def run(self):
        pass