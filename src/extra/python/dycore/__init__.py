"""
temporary docstring for src.extra.python.dycore
"""
import os, sys

from typing import Any

_module_directory = os.path.dirname(__file__)

try:
    Dycore_BASE        = os.environ['Dycore_BASE']
    Dycore_WORK        = os.environ['Dycore_WORK']
    Dycore_DATA        = os.environ['Dycore_DATA']
except Exception as e:
    # log.error('Environment variables Dycore_BASE, Dycore_WORK, Dycore_DATA must be set')
    raise ValueError('Environment variables Dycore_BASE, Dycore_WORK, Dycore_DATA must be set')

from dycore.experiment import Experiment
from dycore.config import Config
from dycore.diagtable import DiagTable