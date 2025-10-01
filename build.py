"""
This file is used during installation to compile all motif prediction tools,
and move them into place.
"""

import os
import subprocess as sp
import sys
from platform import processor

from setuptools.command.build_py import build_py

MOTIF_BINS = {
    "AMD": ["src/AMD/AMD.bin"],
    "BioProspector": ["src/BioProspector/BioProspector"],
    "MDmodule": ["src/MDmodule/MDmodule"],
    "Posmo": ["src/posmo/posmo", "src/posmo/clusterwd"],
}


class custom_build_py(build_py):  # noqa

    def initialize_options(self):
        super().initialize_options()
        self.build_base = None  # noqa
        self.build_lib = None  # noqa

    def finalize_options(self):
        super().finalize_options()
        self.set_undefined_options("build", ("build_base", "build_base"))
        self.set_undefined_options("build", ("build_lib", "build_lib"))

    def run(self):
        if not self.dry_run:
            # copy all tools to the build directory
            src_dir = os.path.join(self.build_base, "src")
            self.copy_tree("src/", src_dir)

            # compile tools where possible
            compile_tool("BioProspector", src_dir="src")
            compile_tool("MDmodule", src_dir="src")

            # copy tools binaries to the target_dir
            target_dir = os.path.join(self.build_lib, "gimmemotifs/included_tools")
            self.mkpath(target_dir)
            for exes in MOTIF_BINS.values():
                for exe in exes:
                    # if os.path.exists(exe):
                    #     self.copy_file(exe, target_dir)
                    exe = os.path.join(self.build_base, exe)
                    if os.path.exists(exe):
                        self.copy_file(exe, target_dir)

            # copy tool directories to the target_dir
            self.copy_tree(
                os.path.join(src_dir, "ChIPMunk"), os.path.join(target_dir, "ChIPMunk")
            )
            self.copy_tree(
                os.path.join(src_dir, "HMS"), os.path.join(target_dir, "HMS")
            )

            # copy pre-compiled tool binaries to the target_dir
            arch = processor()
            if arch in ["x86_64", "i386"]:
                self.copy_file(
                    os.path.join(src_dir, f"MotifSampler/MotifSampler_{arch}"),
                    os.path.join(target_dir, "MotifSampler"),
                )
                self.copy_file(
                    os.path.join(src_dir, f"MotifSampler/CreateBackgroundModel_{arch}"),
                    os.path.join(target_dir, "CreateBackgroundModel"),
                )
                self.copy_file(
                    os.path.join(src_dir, f"Improbizer/ameme_{arch}"),
                    os.path.join(target_dir, "ameme"),
                )

        super().run()


def compile_tool(name, src_dir="src"):
    sys.stderr.write(f"compiling {name}")
    sys.stderr.flush()

    # validate tool
    path = os.path.join(src_dir, name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find {path}")

    # select a compiler
    gcc = "gcc"
    if os.environ.get("GCC"):
        gcc = os.environ["GCC"]
    if os.environ.get("CC"):
        gcc = os.environ["CC"]
    try:
        sp.Popen(gcc, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    except Exception:  # noqa
        sys.stderr.write("... failed (compiler not found) \n")
        return

    # compile
    sp.Popen(
        [gcc, f"-o{name}", f"{name}.c", "-lm"], cwd=path, stdout=sp.PIPE
    ).communicate()
    if os.path.exists(os.path.join(path, name)):
        sys.stderr.write("... ok\n")
    else:
        sys.stderr.write("... failed (compilation failed)\n")
