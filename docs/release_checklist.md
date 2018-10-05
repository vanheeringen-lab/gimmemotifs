# Release checklist

This is mainly for personal use at the moment.

## 

1. Create release candidate with `git flow`.

```
$ git flow release start ${new_version} 
```

2. Update version in `gimmemotifs/config.py`

3. Make sure `CHANGELOG.md` is up-to-date.

4. Test install using pip in fresh conda environment

```
$ cd ${test_dir} 	# Not the gimmemotifs git directory
$ conda create -n testenv python=3 --file conda_env.txt 
$ conda activate testenv
$ pip install -e git+https://github.com/simonvh/gimmemotifs.git@release/${version}#egg=gimmemotifs
$ cd src/gimmemotifs
$ python run_tests.py
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

```
$ git flow release finish ${version}
```

7. Upload to PyPi.

```
$ python setup.py sdist
$ twine upload dist/gimmemotifs-${version}.tar.gz
```

8. Finalize the release on Github.

Create a release. Download the tarball and then edit the release and attach the
tarball as binary. 

8. Create bioconda package
