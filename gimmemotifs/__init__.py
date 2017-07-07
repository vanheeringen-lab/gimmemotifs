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
