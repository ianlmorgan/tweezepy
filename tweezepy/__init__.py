from .smmcalibration import AV
from .smmcalibration import PSD
from .simulations import simulate_trace
from .simulations import downsampled_trace


import os

from pkg_resources import DistributionNotFound,get_distribution

try:
    __version__ = get_distribution("tweezepy").version
except DistributionNotFound:
    __version__ = "unknown version"

#this_directory = os.path.abspath(os.path.dirname(__file__))
#pkginfo_path = os.path.join(this_directory,
#                            'tweezepy_info.json')
#pkginfo = json.load(open(pkginfo_path))

#__version__ = pkginfo["version"]
__all__ = [
           "AV",
           "PSD",
           "simulate_trace",
           "downsampled_trace",
           "__version__",
           ]