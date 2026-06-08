# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, A QIIME 2 Plugin Developer.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import metabat2, busco, filtering, utils, abundance, dereplication

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = ["metabat2", "busco", "filtering", "utils", "abundance", "dereplication"]
