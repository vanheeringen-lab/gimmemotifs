# Release checklist

This is mainly for personal use at the moment.

## 0. (Re)install gimmemotifs

```shell
conda activate base
conda env remove -n gimme -yq
mamba env create -n gimme -f requirements.yaml
conda activate gimme
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
rm -rf dist
pip wheel -w dist --no-deps --no-cache-dir --use-pep517 -v .
conda env remove -n test -yq
mamba env create --force -n test -f requirements.yaml
conda activate test
pip install --no-deps --no-cache-dir --use-pep517 -v dist/gimmemotifs*.whl

rm ~/.config/gimmemotifs/gimmemotifs.cfg
gimme -h
pytest -vvv --disable-pytest-warnings

```

## 5. Upload to pypi testing server

For more info, see https://packaging.python.org/en/latest/tutorials/packaging-projects/#generating-distribution-archives

```shell
# Create the source distribution
python3 -m pip install --upgrade build
python3 -m build

# Test the source distribution
pip install --no-deps --no-cache-dir --use-pep517 -v dist/gimmemotifs*.tar.gz
rm ~/.config/gimmemotifs/gimmemotifs.cfg
gimme -h
pytest -vvv --disable-pytest-warnings

# Upload to pypi testing server
python3 -m pip install --upgrade twine
python3 -m twine upload --repository testpypi dist/*

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
twine upload dist/*
```

## 9. Finalize the release on Github.

Create a release. Download the tarball and then edit the release and attach the
tarball as binary (this binary's url is currently used by bioconda). 

## 10. Update the Bioconda package recipe

Wait for the automatic PR on the bioconda-recipes github, 
update the meta.yaml with the latest requirements, 
and approve/request approval.
