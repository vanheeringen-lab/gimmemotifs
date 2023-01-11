# Release checklist

This is mainly for personal use at the moment.

## 


1. Make sure all tests pass.

```shell
mamba env update -f environment.yml
pytest -vvv
```

2. Create release candidate with `git flow`:

```shell
new_version=0.0.0
echo ${new_version}

git flow release start ${new_version}
```

3. Make sure `CHANGELOG.md` is up-to-date.

    * add the new version & date to the header
    * link to the diff in the footer
    * add & commit the changes, but do not push

4. Test install using pip in fresh conda environment

```shell
python setup.py sdist
mamba create -n test python=3.9 pytest
mamba activate test
pip install dist/gimmemotifs*.tar.gz
pytest -vvv
```

5. Upload to pypi testing server

```
# Check for warnings or errors
$ python setup.py check -r -s
# Create distribution
$ python setup.py sdist
# Upload to pypi testing server
$ twine upload -r testpypi dist/gimmemotifs-${version}.tar.gz
``` 

6. Finish release

```shell
git flow release finish ${new_version}
```

7. Upload to PyPi.

```shell
python setup.py sdist
twine upload dist/gimmemotifs-${new_version}.tar.gz
```

8. Finalize the release on Github.

Create a release. Download the tarball and then edit the release and attach the
tarball as binary. 

9. Update the Bioconda package recipe
