import pytest
import subprocess as sp


def test_black_formatting():
    try:
        sp.check_output(
            "black --check gimmemotifs/ setup.py scripts/", stderr=sp.STDOUT, shell=True
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")
        msg = "\n".join(msg.split("\n")[:-3])  # multi line error
        print("Black output:")
        print(msg)

        msg = ", ".join(msg.split("\n"))[:-2]  # single line error
        pytest.fail(msg)


def test_flake8_formatting():
    try:
        sp.check_output(
            "flake8 setup.py gimmemotifs/ scripts/", stderr=sp.STDOUT, shell=True
        )
    except sp.CalledProcessError as e:
        msg = e.output.decode("utf-8")  # multi line error
        print("Flake8 output:")
        print(msg)

        msg = ", ".join(msg.split("\n"))[:-2]  # single line error
        pytest.fail(msg)
