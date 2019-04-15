import atexit
from os import getpid
import shutil
from tempfile import mkdtemp
import logging

def mytmpdir():
    if not hasattr(mytmpdir, 'dir') or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(prefix="gimmemotifs.{0}.".format(getpid()))
        atexit.register(shutil.rmtree, mytmpdir.dir)
    return mytmpdir.dir

class DuplicateFilter(logging.Filter):
    # https://stackoverflow.com/questions/44691558
    def filter(self, record):
        # add other fields if you need more granular comparison, depends on your app
        current_log = (record.module, record.levelno, record.msg)
        if current_log != getattr(self, "last_log", None):
            self.last_log = current_log
            return True
        return False

logger = logging.getLogger('gimme')
logger.setLevel(logging.DEBUG)
logger.propagate = 0

# nice format
screen_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# Log to screen
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(screen_formatter)
logger.addHandler(sh)
logger.addFilter(DuplicateFilter())

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
