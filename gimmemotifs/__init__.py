import atexit as _atexit
import logging as _logging
import warnings as _warnings
from os import getpid as _getpid
from shutil import rmtree as _rmtree
from tempfile import mkdtemp as _mkdtemp
from warnings import warn as _warn


def _filtered_warn(*args, **kwargs):
    """warnings.filterwarnings does not work."""
    blacklist = [
        # maelstrom sklearn SVM
        "Liblinear failed to converge, increase the number of iterations.",
        # report.py (pandas internal DeprecationWarning)
        "The get_cmap function will be deprecated in a future version. "
        "Use ``matplotlib.colormaps[name]`` or "
        "``matplotlib.colormaps.get_cmap(obj)`` instead.",
        # report.py pandas FutureWarning, alternative works differently
        "this method is deprecated in favour of `Styler.to_html()`",
        # scanner/base.py set_threshold
        "DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, "
        "which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. "
        "To get a de-fragmented frame, use `newframe = frame.copy()`",
    ]
    message = str(args[0])
    if message in blacklist:
        pass
    else:
        _warn(*args, **kwargs)


_warnings.warn = _filtered_warn


def mytmpdir():
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        mytmpdir.dir = _mkdtemp(prefix=f"gimmemotifs.{_getpid()}.")
        _atexit.register(_rmtree, mytmpdir.dir)
    return mytmpdir.dir


# setup logger
logger = _logging.getLogger("gimme")
if logger.hasHandlers():
    logger.handlers.clear()

logger.setLevel(_logging.DEBUG)

# nice format
screen_formatter = _logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# log to screen
sh = _logging.StreamHandler()
sh.setLevel(_logging.INFO)
sh.setFormatter(screen_formatter)
logger.addHandler(sh)
del screen_formatter
del sh


# set __version__
from ._version import get_versions  # noqa: E402

__version__ = get_versions()["version"]
del get_versions


# import submodules to enable module autocomplete in editor

# fmt: off
# isort: off
from . import fasta            # noqa: F401
from . import plot             # noqa: F401
from . import rocmetrics       # noqa: F401
from . import config           # noqa: F401
from . import utils            # noqa: F401
# from . import commands       # not needed in API
# from . import preprocessing  # not needed in API
from . import motif            # noqa: F401
# from . import orthologs      # not needed in API
from . import comparison       # noqa: F401
from .motif import cluster     # noqa: F401
from . import background       # noqa: F401
from . import scanner          # noqa: F401
from . import stats            # noqa: F401
from . import report           # noqa: F401
# from . import tools          # not needed in API
from .motif import denovo      # noqa: F401
from . import maelstrom        # noqa: F401
# from . import commands       # not needed in API
# from . import cli            # not needed in API
# fmt: on
# isort: on
