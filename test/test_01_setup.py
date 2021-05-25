import os
import pytest
import subprocess as sp

travis = os.environ.get("TRAVIS") == "true"


def test_black_formatting():
    if not travis:
        sp.check_call("black gimmemotifs/ setup.py scripts/", shell=True)
        pass

    try:
        sp.check_output(
            "black --check gimmemotifs/ setup.py scripts/", stderr=sp.STDOUT, shell=True
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = msg.split("\n")[:-3]
        msg = "\n".join(["Black output:"] + msg)
        pytest.fail(msg, False)


def test_flake8_formatting():
    try:
        sp.check_output(
            "flake8 setup.py gimmemotifs/ scripts/", stderr=sp.STDOUT, shell=True
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = "Flake8 output:\n" + msg
        pytest.fail(msg, False)
