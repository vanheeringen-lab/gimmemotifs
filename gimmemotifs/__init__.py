import atexit
import shutil
from tempfile import mkdtemp

mytmpdir = mkdtemp(prefix="gimmemotifs.")
atexit.register(shutil.rmtree, mytmpdir)
