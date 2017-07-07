import pkgutil
import six
import os

dirname = os.path.split(__file__)[0]

if six.PY3:
    level = 0
else:
    level = -1

# Dynamically load all commands
for importer, cmdname, _ in pkgutil.iter_modules([dirname]):
    m =  __import__("{0}.{1}".format(__name__, cmdname), globals(), locals(), [cmdname], level)
    globals()[cmdname] = getattr(m, cmdname)
