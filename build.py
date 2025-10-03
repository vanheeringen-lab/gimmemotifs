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
    "AMD": [os.path.join("src", "AMD", "AMD.bin")],
    "Posmo": [
        os.path.join("src", "posmo", "posmo"),
        os.path.join("src", "posmo", "clusterwd"),
    ],
}


class custom_build_py(build_py):  # noqa
    def initialize_options(self):
        super().initialize_options()
        self.build_lib = None  # noqa

    def finalize_options(self):
        super().finalize_options()
        self.set_undefined_options("build", ("build_lib", "build_lib"))

    def run(self):
        if not self.dry_run:
            # # copy all tools to the build directory
            # src_dir = os.path.join(self.build_base, "src")
            # self.copy_tree("src", src_dir)
            src_dir = os.path.join(os.path.dirname(__file__), "src")
            target_dir = os.path.join(self.build_lib, "gimmemotifs", "included_tools")
            self.mkpath(target_dir)

            # compile tools where possible
            compile_tool("BioProspector", src_dir=src_dir, target_dir=target_dir)
            compile_tool("MDmodule", src_dir=src_dir, target_dir=target_dir)

            # copy tools binaries to the target_dir
            for tool, exes in MOTIF_BINS.items():
                for exe in exes:
                    # exe = os.path.join(self.build_base, exe)
                    self.copy_file(exe, target_dir, level=self.verbose)
                    # if os.path.exists(exe):
                    #     self.copy_file(exe, target_dir, level=self.verbose)
                    # else:
                    #     sys.stderr.write(f"Did not find {tool} binary {exe}\n")

            # copy tool directories to the target_dir
            self.copy_tree(
                os.path.join(src_dir, "ChIPMunk"),
                os.path.join(target_dir, "ChIPMunk"),
                level=self.verbose,
            )
            self.copy_tree(
                os.path.join(src_dir, "HMS"),
                os.path.join(target_dir, "HMS"),
                level=self.verbose,
            )

            # copy pre-compiled tool binaries to the target_dir
            arch = processor()
            if arch in ["x86_64", "i386"]:
                self.copy_file(
                    os.path.join(src_dir, "MotifSampler", f"MotifSampler_{arch}"),
                    os.path.join(target_dir, "MotifSampler"),
                    level=self.verbose,
                )
                self.copy_file(
                    os.path.join(
                        src_dir, "MotifSampler", f"CreateBackgroundModel_{arch}"
                    ),
                    os.path.join(target_dir, "CreateBackgroundModel"),
                    level=self.verbose,
                )
                self.copy_file(
                    os.path.join(src_dir, "Improbizer", f"ameme_{arch}"),
                    os.path.join(target_dir, "ameme"),
                    level=self.verbose,
                )

            # Copy the motif detection tools to the base directory.
            #
            # This is only needed when installing gimme in
            # development/editable mode (pip install -e)
            # TODO: automatically delete these files when uninstalling
            if self.editable_mode:
                self.copy_tree(
                    os.path.join(self.build_lib, "gimmemotifs", "included_tools"),
                    os.path.join(
                        os.path.dirname(__file__), "gimmemotifs", "included_tools"
                    ),
                    level=self.verbose,
                )

        super().run()


def compile_tool(name, src_dir, target_dir):
    # select a compiler
    gcc = "gcc"
    if os.environ.get("GCC"):
        gcc = os.environ["GCC"]
    if os.environ.get("CC"):
        gcc = os.environ["CC"]
    try:
        sp.Popen(gcc, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    except Exception:  # noqa
        raise FileNotFoundError(f"(G)CC compiler not found!")

    # compile
    infile = os.path.join(src_dir, name, f"{name}.c")
    outfile = os.path.join(target_dir, f"{name}")
    if not os.path.exists(infile):
        raise FileNotFoundError(f"Could not find {infile}")
    sp.Popen(
        [gcc, f"-o{outfile}", infile, "-lm"], stdout=sp.PIPE, stderr=sp.PIPE
    ).communicate()
    if os.path.exists(outfile):
        sys.stderr.write(f"compiled {infile} -> {outfile}\n")
    else:
        raise FileNotFoundError(f"{name} failed to compile\n")


# def uninstall():
#     """
#     Remove the motif detection tools installed in editable mode.
#     """
#     tool_dir = os.path.join(os.path.dirname(__file__), "gimmemotifs", "included_tools")
#     assert os.path.exists(tool_dir)
#     for tool in os.listdir(tool_dir):
#         if tool == "__init__.py":
#             continue
#         elif os.path.isdir(tool):
#             shutil.rmtree(tool)
#         else:
#             os.remove(tool)
