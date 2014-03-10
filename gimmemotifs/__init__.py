import atexit
from os import getpid
import shutil
from tempfile import mkdtemp

def mytmpdir():
    if not hasattr(mytmpdir, 'dir') or not mytmpdir.dir:
        mytmpdir.dir = mkdtemp(prefix="gimmemotifs.{0}.".format(getpid()))
        atexit.register(shutil.rmtree, mytmpdir.dir)
    return mytmpdir.dir

