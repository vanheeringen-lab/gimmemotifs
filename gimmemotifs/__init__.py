import atexit
from os import getpid
import shutil
from tempfile import mkdtemp
import logging


def mytmpdir():
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(prefix="gimmemotifs.{0}.".format(getpid()))
        atexit.register(shutil.rmtree, mytmpdir.dir)
    return mytmpdir.dir


class DuplicateFilter(logging.Filter):
    logs = {}
    # based on https://stackoverflow.com/questions/44691558

    def filter(self, record):
        if record.levelno != logging.INFO:
            return True

        if record.module not in self.logs:
            self.logs[record.module] = {}

        if record.msg in self.logs[record.module]:
            return False
        else:
            self.logs[record.module][record.msg] = 1
            return True


logger = logging.getLogger("gimme")
if logger.hasHandlers():
    logger.handlers.clear()

logger.setLevel(logging.DEBUG)
# logger.addFilter(DuplicateFilter('gimme'))

# nice format
screen_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# Log to screen
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(screen_formatter)
# sh.addFilter(DuplicateFilter('gimme'))
logger.addHandler(sh)


from ._version import get_versions  # noqa: E402

__version__ = get_versions()["version"]
del get_versions

# fmt: off
# easier import of gimme (config and cli left out)
from . import background  # noqa: F401
from . import cluster     # noqa: F401
from . import comparison  # noqa: F401
from . import denovo      # noqa: F401
from . import fasta       # noqa: F401
from . import maelstrom   # noqa: F401
from . import moap        # noqa: F401
from . import motif       # noqa: F401
from . import plot        # noqa: F401
from . import prediction  # noqa: F401
from . import rank        # noqa: F401
from . import report      # noqa: F401
from . import rocmetrics  # noqa: F401
from . import scanner     # noqa: F401
from . import shutils     # noqa: F401
from . import stats       # noqa: F401
from . import utils       # noqa: F401
from . import validation  # noqa: F401
# fmt: on
