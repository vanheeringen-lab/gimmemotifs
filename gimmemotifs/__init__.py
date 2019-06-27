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
