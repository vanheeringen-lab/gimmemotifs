# Release checklist

This is mainly for personal use at the moment.

## 0. (Re)install gimmemotifs

```shell
mamba activate base
mamba env create --force -n gimme -f requirements.yaml
mamba activate gimme
pip install --no-deps --no-cache-dir --use-pep517 -v -e .
```

## 1. Make sure all tests pass.

```shell
rm ~/.config/gimmemotifs/gimmemotifs.cfg
gimme -h
pytest -vvv --disable-pytest-warnings
```

## 2. Create release candidate with `git flow`:

```shell
new_version=0.0.0
echo ${new_version}

git flow release start ${new_version}
```

## 3. Make sure `__about__.py`, `pyproject.toml`, `CHANGELOG.md` are up-to-date.

    * set the new version in `__about__.py`
    * make sure all subpackages are listed in the `pyproject.toml`
    * add the new version & date to the header of `CHANGELOG.md`
    * link to the diff in the footer of `CHANGELOG.md`
    * add & commit the changes, but do not push

## 4. Test install using pip in fresh conda environment

```shell
pip wheel -w dist --no-deps --no-cache-dir --use-pep517 -v .
mamba env create -n test -f requirements.yaml
mamba activate test
pip install --no-deps --no-cache-dir --use-pep517 -v dist/gimmemotifs*.whl

rm ~/.config/gimmemotifs/gimmemotifs.cfg
gimme -h
pytest -vvv --disable-pytest-warnings

```

## 5. Upload to pypi testing server

```
# Check for warnings or errors
$ python setup.py check -r -s
# Create distribution
$ python setup.py sdist
# Upload to pypi testing server
$ twine upload -r testpypi dist/gimmemotifs-${version}.tar.gz
``` 

## 6. Finish release

```shell
git flow release finish ${new_version}
```


## 7. Push everything to github, including tags:

```shell
git push --follow-tags origin develop master
```

## 8. Upload to PyPi.

```shell
python setup.py sdist
twine upload dist/gimmemotifs-${new_version}.tar.gz
```

## 9. Finalize the release on Github.

Create a release. Download the tarball and then edit the release and attach the
tarball as binary (this binary's url is currently used by bioconda). 

## 10. Update the Bioconda package recipe

Wait for the automatic PR on the bioconda-recipes github, 
update the meta.yaml with the latest requirements, 
and approve/request approval.
