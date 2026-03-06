# hardcoded path
cd ~/git/gimmemotifs

# create python version specific yamls
mkdir test_py
cp requirements.yaml test_py/requirements_py3.9.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.9|' test_py/requirements_py3.9.yaml
cp requirements.yaml test_py/requirements_py3.10.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.10|' test_py/requirements_py3.10.yaml
cp requirements.yaml test_py/requirements_py3.11.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.11|' test_py/requirements_py3.11.yaml
cp requirements.yaml test_py/requirements_py3.12.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.12|' test_py/requirements_py3.12.yaml
cp requirements.yaml test_py/requirements_py3.13.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.13|' test_py/requirements_py3.13.yaml
cp requirements.yaml test_py/requirements_py3.14.yaml
sed -i -e 's|python >=3.9, <3.12|python ==3.14|' test_py/requirements_py3.14.yaml

mamba env list > ~/git/gimmemotifs/test_py/envs.txt
LOG=/home/$USER/git/gimmemotifs/test_py/log.txt
touch $LOG
for PY in 3.9 3.10 3.11 3.12 3.13 3.14; do
  MSG="create a blank python ${PY} conda environments without gimme installed"
  echo ""
  echo $MSG
  echo $MSG >> $LOG
  ENV="gimme_py${PY}_blank"
  FILE="test_py/requirements_py${PY}.yaml"
  if grep -wq "$ENV" test_py/envs.txt; then
    # delete env if already existing
    mamba env remove -n $ENV -yq > /dev/null 2>&1
  fi
  mamba env create -n $ENV -f $FILE -yq > /dev/null 2>&1

  # clone the conda environments for each install method
  MSG="default install"
  echo ""
  echo $MSG
  echo $MSG >> $LOG
  NAME="gimme_py${PY}_dflt"
  if grep -wq "$NAME" test_py/envs.txt; then
    # delete env if already existing
    mamba env remove -n $NAME -yq > /dev/null 2>&1
  fi
  echo "  - cloning environment"
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  mamba activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  echo "  - installing gimmemotifs"
  pip install -q --no-deps --no-cache-dir --use-pep517 .
  which gimme | tee >> $LOG 2>&1
  cd ..
  python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' | tee >> $LOG 2>&1
  cd gimmemotifs
  echo "  - uninstalling gimmemotifs"
  pip uninstall -qy gimmemotifs
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  conda deactivate
  echo "  - done"

  MSG="wheel install"
  echo ""
  echo $MSG
  echo $MSG >> $LOG
  NAME="gimme_py${PY}_whl"
  if grep -wq "$NAME" test_py/envs.txt; then
    # delete env if already existing
    mamba env remove -n $NAME -yq > /dev/null 2>&1
  fi
  echo "  - cloning environment"
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  mamba activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/dist
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  echo "  - installing gimmemotifs"
  pip wheel -q -w dist --no-deps --no-cache-dir --use-pep517 -v .
  pip install -q --no-deps --no-cache-dir --use-pep517 dist/gimmemotifs*.whl
  which gimme | tee >> $LOG 2>&1
  cd ..
  python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' | tee >> $LOG 2>&1
  cd gimmemotifs
  echo "  - uninstalling gimmemotifs"
  pip uninstall -qy gimmemotifs
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/dist
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  conda deactivate
  echo "  - done"

  MSG="pip install"
  echo ""
  echo $MSG
  echo $MSG >> $LOG
  NAME="gimme_py${PY}_pip"
  if grep -wq "$NAME" test_py/envs.txt; then
    # delete env if already existing
    mamba env remove -n $NAME -yq > /dev/null 2>&1
  fi
  echo "  - cloning environment"
  mamba create --name $NAME "python ==${PY}" "setuptools ==78.1.1" pip -yq > /dev/null 2>&1
  mamba activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  echo "  - installing gimmemotifs"
  pip install -q --no-cache-dir --use-pep517 .
  which gimme | tee >> $LOG 2>&1
  cd ..
  python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' | tee >> $LOG 2>&1
  cd gimmemotifs
  echo "  - uninstalling gimmemotifs"
  pip uninstall -qy gimmemotifs
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/build
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  conda deactivate
  echo "  - done"

  MSG="editable install"
  echo ""
  echo $MSG
  echo $MSG >> $LOG
  NAME="gimme_py${PY}_edit"
  if grep -wq "$NAME" test_py/envs.txt; then
    # delete env if already existing
    mamba env remove -n $NAME -yq > /dev/null 2>&1
  fi
  echo "  - cloning environment"
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  mamba activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  # editable specific:
  rm -f ~/git/gimmemotifs/gimmemotifs/c_metrics.cpython-3*-x86_64-linux-gnu.so
  rm -rf ~/git/gimmemotifs/gimmemotifs/included_tools/ChIPMunk
  rm -rf ~/git/gimmemotifs/gimmemotifs/included_tools/HMS
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/AMD.bin
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/ameme
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/BioProspector
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/clusterwd
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/CreateBackgroundModel
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/MDmodule
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/MotifSampler
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/posmo
  echo "  - installing gimmemotifs"
  pip install -q --no-deps --no-cache-dir --use-pep517 -e .
  which gimme 2>&1 | tee -a $LOG
  python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' | tee >> $LOG 2>&1

  echo "  - run tests"
  echo "" | tee >> $LOG 2>&1
  pytest -vvv --disable-pytest-warnings 2>&1 | tee -a $LOG
  echo "" | tee >> $LOG 2>&1
  echo "  - completed tests"

  echo "  - uninstalling gimmemotifs"
  pip uninstall -qy gimmemotifs
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf ~/git/gimmemotifs/gimmemotifs.egg-info
  # editable specific:
  rm -f ~/git/gimmemotifs/gimmemotifs/c_metrics.cpython-3*-x86_64-linux-gnu.so
  rm -rf ~/git/gimmemotifs/gimmemotifs/included_tools/ChIPMunk
  rm -rf ~/git/gimmemotifs/gimmemotifs/included_tools/HMS
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/AMD.bin
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/ameme
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/BioProspector
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/clusterwd
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/CreateBackgroundModel
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/MDmodule
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/MotifSampler
  rm -f ~/git/gimmemotifs/gimmemotifs/included_tools/posmo
  conda deactivate
  echo "  - done"
done
