#!/usr/bin/bash

# hardcoded path
cd ~/git/gimmemotifs
if ! [[ $(pwd) =~ /gimmemotifs ]]
then
  echo "working directory is not the gimmemotifs repo"
  exit 1
fi

# allow us to use conda
source $(dirname $(dirname $(which mamba)))/etc/profile.d/conda.sh
source $(dirname $(dirname $(which mamba)))/etc/profile.d/mamba.sh
if ! which mamba > /dev/null 2>&1
then
  echo "conda/mamba not found!"
  exit 1
fi

EXIT_CODE=0

# delete old runs
rm -rf test_py
rm -rf $(dirname $(dirname $(which mamba)))/envs/gimme_py*

# create python version specific yamls
mkdir test_py 2> /dev/null
cat requirements.yaml > test_py/template.yaml
echo "  - conda-forge::conda-ecosystem-user-package-isolation=1.0" >> test_py/template.yaml
cp test_py/template.yaml test_py/requirements_py3.9.yaml
sed -i -e 's|python >=3.9|python =3.9|' -e '/orthofinder/d' -e '/gffread/d' test_py/requirements_py3.9.yaml
cp test_py/template.yaml test_py/requirements_py3.10.yaml
sed -i -e 's|python >=3.9|python =3.10|' -e '/orthofinder/d' -e '/gffread/d' test_py/requirements_py3.10.yaml
cp test_py/template.yaml test_py/requirements_py3.11.yaml
sed -i -e 's|python >=3.9|python =3.11|' test_py/requirements_py3.11.yaml
cp test_py/template.yaml test_py/requirements_py3.12.yaml
sed -i -e 's|python >=3.9|python =3.12|' test_py/requirements_py3.12.yaml
cp test_py/template.yaml test_py/requirements_py3.13.yaml
sed -i -e 's|python >=3.9|python =3.13|' test_py/requirements_py3.13.yaml
cp test_py/template.yaml test_py/requirements_py3.14.yaml
sed -i -e 's|python >=3.9|python =3.14|' test_py/requirements_py3.14.yaml

LOG=test_py/log.txt
touch $LOG
echo "starting" 2>&1 | tee -a $LOG
for PY in 3.9 3.10 3.11 3.12; do  # scikit-learn is incompatible with >=3.13
  MSG="create a blank python ${PY} conda environments without gimme installed"
  echo -e "\n###########################################################################\n" 2>&1 | tee -a $LOG
  echo $MSG 2>&1 | tee -a $LOG
  ENV="gimme_py${PY}_blank"
  FILE="test_py/requirements_py${PY}.yaml"
  CONDA_LOG="test_py/log_conda_py${PY}.txt"
  if ! mamba env create -n $ENV -f $FILE -y > $CONDA_LOG 2>&1
  then
    echo "conda env ${ENV} creation failed! See ${CONDA_LOG} for details" 2>&1 | tee -a $LOG
    EXIT_CODE=1
    continue
  fi
  # log is only interesting upon a failure
  rm $CONDA_LOG

  # clone the conda environments for each install method
  MSG="default install - py${PY}"
  echo -e "\n###########################################################################\n" 2>&1 | tee -a $LOG
  echo $MSG 2>&1 | tee -a $LOG
  NAME="gimme_py${PY}_dflt"
  echo "  - cloning environment" 2>&1 | tee -a $LOG
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  conda activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf build
  rm -rf gimmemotifs.egg-info
  echo "  - installing gimmemotifs" 2>&1 | tee -a $LOG
  pip install -q --no-deps --no-cache-dir --use-pep517 .
  if which gimme > /dev/null 2>&1
  then
    echo "  - success" 2>&1 | tee -a $LOG
    echo "  - testing import:" 2>&1 | tee -a $LOG
    cd ..
    python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' 2>&1 | tee -a gimmemotifs/$LOG
    cd gimmemotifs
    echo "  - uninstalling gimmemotifs" 2>&1 | tee -a $LOG
    pip uninstall -qy gimmemotifs
  else
    echo "  - fail" 2>&1 | tee -a $LOG
    EXIT_CODE=1
  fi
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf build
  rm -rf gimmemotifs.egg-info
  conda deactivate
  echo "  - done" 2>&1 | tee -a $LOG

  MSG="wheel install - py${PY}"
  echo -e "\n###########################################################################\n" 2>&1 | tee -a $LOG
  echo $MSG 2>&1 | tee -a $LOG
  NAME="gimme_py${PY}_whl"
  echo "  - cloning environment" 2>&1 | tee -a $LOG
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  conda activate $NAME
  if ! command -v pip >/dev/null 2>&1
  then
    echo "  - installing pip manually" 2>&1 | tee -a $LOG
    mamba install --no-deps -y pip
  fi
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf dist
  rm -rf build
  rm -rf gimmemotifs.egg-info
  echo "  - installing gimmemotifs" 2>&1 | tee -a $LOG
  pip wheel -q -w dist --no-deps --no-cache-dir --use-pep517 -v .
  pip install -q --no-deps --no-cache-dir --use-pep517 dist/gimmemotifs*.whl
  if which gimme > /dev/null 2>&1
  then
    echo "  - success" 2>&1 | tee -a $LOG
    echo "  - testing import:" 2>&1 | tee -a $LOG
    cd ..
    python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' 2>&1 | tee -a gimmemotifs/$LOG
    cd gimmemotifs
    echo "  - uninstalling gimmemotifs" 2>&1 | tee -a $LOG
    pip uninstall -qy gimmemotifs
  else
    echo "  - fail" 2>&1 | tee -a $LOG
    EXIT_CODE=1
  fi
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf dist
  rm -rf build
  rm -rf gimmemotifs.egg-info
  conda deactivate
  echo "  - done" 2>&1 | tee -a $LOG

  MSG="pip install - py${PY}"
  echo -e "\n###########################################################################\n" 2>&1 | tee -a $LOG
  echo $MSG 2>&1 | tee -a $LOG
  NAME="gimme_py${PY}_pip"
  echo "  - creating blank environment" 2>&1 | tee -a $LOG
  mamba create --name $NAME "python =${PY}" "setuptools =78.1.1" pip pytest -yq > /dev/null 2>&1
  conda activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf build
  rm -rf gimmemotifs.egg-info
  echo "  - installing gimmemotifs" 2>&1 | tee -a $LOG
  pip install -q --no-cache-dir --use-pep517 .
  if which gimme > /dev/null 2>&1
  then
    echo "  - success" 2>&1 | tee -a $LOG
    echo "  - testing import:" 2>&1 | tee -a $LOG
    cd ..
    python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' 2>&1 | tee -a gimmemotifs/$LOG
    cd gimmemotifs
    echo "  - uninstalling gimmemotifs" 2>&1 | tee -a $LOG
    pip uninstall -qy gimmemotifs
  else
    echo "  - fail" 2>&1 | tee -a $LOG
    EXIT_CODE=1
  fi
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf build
  rm -rf gimmemotifs.egg-info
  conda deactivate
  echo "  - done" 2>&1 | tee -a $LOG

  MSG="editable install - py${PY}"
  echo -e "\n###########################################################################\n" 2>&1 | tee -a $LOG
  echo $MSG 2>&1 | tee -a $LOG
  NAME="gimme_py${PY}_edit"
  echo "  - cloning environment" 2>&1 | tee -a $LOG
  mamba create --name $NAME --clone $ENV -yq > /dev/null 2>&1
  conda activate $NAME
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf gimmemotifs.egg-info
  # editable specific:
  rm -f gimmemotifs/c_metrics.cpython-3*-x86_64-linux-gnu.so
  rm -rf gimmemotifs/included_tools/ChIPMunk
  rm -rf gimmemotifs/included_tools/HMS
  rm -f gimmemotifs/included_tools/AMD.bin
  rm -f gimmemotifs/included_tools/ameme
  rm -f gimmemotifs/included_tools/BioProspector
  rm -f gimmemotifs/included_tools/clusterwd
  rm -f gimmemotifs/included_tools/CreateBackgroundModel
  rm -f gimmemotifs/included_tools/MDmodule
  rm -f gimmemotifs/included_tools/MotifSampler
  rm -f gimmemotifs/included_tools/posmo
  echo "  - installing gimmemotifs" 2>&1 | tee -a $LOG
  pip install -q --no-deps --no-cache-dir --use-pep517 -e .
  if which gimme > /dev/null 2>&1
  then
    echo "  - success" 2>&1 | tee -a $LOG
    echo "  - testing import:" 2>&1 | tee -a $LOG
    python -c 'from gimmemotifs.config import MotifConfig; print(MotifConfig().bin("AMD"))' 2>&1 | tee -a $LOG
    echo "  - run tests" 2>&1 | tee -a $LOG
    echo "" 2>&1 | tee -a $LOG
    pytest -vvv --disable-pytest-warnings 2>&1 | tee -a $LOG
    echo "" 2>&1 | tee -a $LOG
    echo "  - completed tests" 2>&1 | tee -a $LOG
    echo "  - uninstalling gimmemotifs" 2>&1 | tee -a $LOG
    pip uninstall -qy gimmemotifs
  else
    echo "  - fail" 2>&1 | tee -a $LOG
    EXIT_CODE=1
  fi
  rm -f ~/.config/gimmemotifs/gimmemotifs.cfg
  rm -rf gimmemotifs.egg-info
  # editable specific:
  rm -f gimmemotifs/c_metrics.cpython-3*-x86_64-linux-gnu.so
  rm -rf gimmemotifs/included_tools/ChIPMunk
  rm -rf gimmemotifs/included_tools/HMS
  rm -f gimmemotifs/included_tools/AMD.bin
  rm -f gimmemotifs/included_tools/ameme
  rm -f gimmemotifs/included_tools/BioProspector
  rm -f gimmemotifs/included_tools/clusterwd
  rm -f gimmemotifs/included_tools/CreateBackgroundModel
  rm -f gimmemotifs/included_tools/MDmodule
  rm -f gimmemotifs/included_tools/MotifSampler
  rm -f gimmemotifs/included_tools/posmo
  conda deactivate
  echo "  - done" 2>&1 | tee -a $LOG

  echo "" 2>&1 | tee -a $LOG
done

exit $EXIT_CODE
