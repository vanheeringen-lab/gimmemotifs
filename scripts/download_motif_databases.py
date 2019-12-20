from gimmemotifs.db import MotifDb

for dbname in MotifDb.list_databases():
    print("Downloading:" < dbname)
    db = MotifDb.create(dbname)
    db.download()
