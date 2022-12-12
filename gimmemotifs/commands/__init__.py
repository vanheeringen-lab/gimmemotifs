import os
import pkgutil

dirname = os.path.split(__file__)[0]

level = 0

# Dynamically load all commands
for _importer, cmdname, _ in pkgutil.iter_modules([dirname]):
    m = __import__(f"{__name__}.{cmdname}", globals(), locals(), [cmdname], level)
    globals()[cmdname] = getattr(m, cmdname)
