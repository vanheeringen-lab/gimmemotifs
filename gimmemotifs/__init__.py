import atexit
import logging
import shutil
from os import getpid
from tempfile import mkdtemp


def mytmpdir():
    if not hasattr(mytmpdir, "dir") or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(prefix=f"gimmemotifs.{getpid()}.")
        atexit.register(shutil.rmtree, mytmpdir.dir)
    return mytmpdir.dir


# setup logger
logger = logging.getLogger("gimme")
if logger.hasHandlers():
    logger.handlers.clear()

logger.setLevel(logging.DEBUG)

# nice format
screen_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# log to screen
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(screen_formatter)
logger.addHandler(sh)


from ._version import get_versions  # noqa: E402

__version__ = get_versions()["version"]
del get_versions
