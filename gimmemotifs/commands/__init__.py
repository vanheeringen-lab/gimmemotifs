import pkgutil
import os

dirname = os.path.split(__file__)[0]

level = 0

# Dynamically load all commands
for _importer, cmdname, _ in pkgutil.iter_modules([dirname]):
    m = __import__(
        "{0}.{1}".format(__name__, cmdname), globals(), locals(), [cmdname], level
    )
    globals()[cmdname] = getattr(m, cmdname)
