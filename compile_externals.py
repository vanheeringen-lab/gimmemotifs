import os
import sys
from distutils import log
from subprocess import PIPE, Popen


def compile_simple(name, src_dir="src"):
    gcc = "gcc"
    if os.environ.get("GCC"):
        gcc = os.environ["GCC"]
    if os.environ.get("CC"):
        # implements conda patch
        # https://github.com/bioconda/bioconda-recipes/blob/
        # 7ced0315a5f4f7723d99285bfad16b609c660ae6/recipes/gimmemotifs/compile_fix.patch
        gcc = os.environ["CC"]

    path = os.path.join(src_dir, str(name))
    if not os.path.exists(path):
        return

    try:
        Popen(gcc, stdout=PIPE, stderr=PIPE).communicate()
    except Exception:
        return

    Popen([gcc, f"-o{name}", f"{name}.c", "-lm"], cwd=path, stdout=PIPE).communicate()
    if os.path.exists(os.path.join(path, name)):
        return True


# def compile_configmake(name, binary, configure=True, src_dir="src"):
#     path = os.path.join(src_dir, str(name))
#
#     if not os.path.exists(path):
#         return
#
#     if configure:
#         Popen(
#             ["chmod", "+x", "./configure"], cwd=path, stdout=PIPE, stderr=PIPE
#         ).communicate()
#         stdout, stderr = Popen(
#             ["./configure"], cwd=path, stdout=PIPE, stderr=PIPE
#         ).communicate()
#         print(stdout)
#         print(stderr)
#
#     stdout, stderr = Popen(
#         ["make -j 4"], cwd=path, stdout=PIPE, stderr=PIPE, shell=True
#     ).communicate()
#     print(stdout)
#     print(stderr)
#
#     if os.path.exists(os.path.join(path, binary)):
#         return True


def print_result(result):
    if not result:
        log.info("... failed")
    else:
        log.info("... ok")


def compile_all(src_dir="src"):
    sys.stderr.write("compiling BioProspector")
    sys.stderr.flush()
    result = compile_simple("BioProspector", src_dir=src_dir)
    print_result(result)

    sys.stderr.write("compiling MDmodule")
    sys.stderr.flush()
    result = compile_simple("MDmodule", src_dir=src_dir)
    print_result(result)

    return
