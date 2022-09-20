load(
    "//third_party/remote_config:common.bzl",
    "BAZEL_SH",
    "PYTHON_BIN_PATH",
    "PYTHON_LIB_PATH",
    "TF_PYTHON_CONFIG_REPO",
    "auto_config_fail",
    "config_repo_label",
    "execute",
    "get_bash_bin",
    "get_host_environ",
    "get_python_bin",
    "is_windows",
    "raw_exec",
    "read_dir",
)

def _genrule(src_dir, genrule_name, command, outs):
    """Returns a string with a genrule.
    Genrule executes the given command and produces the given outputs.
    """
    return (
        "genrule(\n" +
        '    name = "' +
        genrule_name + '",\n' +
        "    outs = [\n" +
        outs +
        "\n    ],\n" +
        '    cmd = """\n' +
        command +
        '\n   """,\n' +
        ")\n"
    )

def _norm_path(path):
    """Returns a path with '/' and remove the trailing slash."""
    path = path.replace("\\", "/")
    if path[-1] == "/":
        path = path[:-1]
    return path

def _symlink_genrule_for_dir(
        repository_ctx,
        src_dir,
        dest_dir,
        genrule_name,
        src_files = [],
        dest_files = []):
    """Returns a genrule to symlink(or copy if on Windows) a set of files.
    If src_dir is passed, files will be read from the given directory; otherwise
    we assume files are in src_files and dest_files
    """
    if src_dir != None:
        src_dir = _norm_path(src_dir)
        dest_dir = _norm_path(dest_dir)
        files = "\n".join(read_dir(repository_ctx, src_dir))

        # Create a list with the src_dir stripped to use for outputs.
        dest_files = files.replace(src_dir, "").splitlines()
        src_files = files.splitlines()
    command = []
    outs = []
    for i in range(len(dest_files)):
        if dest_files[i] != "":
            # If we have only one file to link we do not want to use the dest_dir, as
            # $(@D) will include the full path to the file.
            dest = "$(@D)/" + dest_dir + dest_files[i] if len(dest_files) != 1 else "$(@D)/" + dest_files[i]

            # Copy the headers to create a sandboxable setup.
            cmd = "cp -f"
            command.append(cmd + ' "%s" "%s"' % (src_files[i], dest))
            outs.append('        "' + dest_dir + dest_files[i] + '",')
    genrule = _genrule(
        src_dir,
        genrule_name,
        " && ".join(command),
        "\n".join(outs),
    )
    return genrule

def _get_python_lib(repository_ctx, python_bin):
    """Gets the python lib path."""
    python_lib = get_host_environ(repository_ctx, PYTHON_LIB_PATH)
    if python_lib != None:
        return python_lib

    # The interesting program to execute.
    print_lib = [
        "from __future__ import print_function",
        "import site",
        "import os",
        "python_paths = []",
        "if os.getenv('PYTHONPATH') is not None:",
        "  python_paths = os.getenv('PYTHONPATH').split(':')",
        "try:",
        "  library_paths = site.getsitepackages()",
        "except AttributeError:",
        "  from distutils.sysconfig import get_python_lib",
        "  library_paths = [get_python_lib()]",
        "all_paths = set(python_paths + library_paths)",
        "paths = []",
        "for path in all_paths:",
        "  if os.path.isdir(path):",
        "    paths.append(path)",
        "if len(paths) >=1:",
        "  print(paths[0])",
    ]

    # The below script writes the above program to a file
    # and executes it. This is to work around the limitation
    # of not being able to upload files as part of execute.
    cmd = "from os import linesep;"
    cmd += "f = open('script.py', 'w');"
    for line in print_lib:
        cmd += "f.write(\"%s\" + linesep);" % line
    cmd += "f.close();"
    cmd += "from subprocess import call;"
    cmd += "call([\"%s\", \"script.py\"]);" % python_bin

    result = execute(repository_ctx, [python_bin, "-c", cmd])
    return result.stdout.strip()

def _check_python_lib(repository_ctx, python_lib):
    """Checks the python lib path."""
    cmd = 'test -d "%s" -a -x "%s"' % (python_lib, python_lib)
    result = raw_exec(repository_ctx, [get_bash_bin(repository_ctx), "-c", cmd])
    if result.return_code == 1:
        auto_config_fail("Invalid python library path: %s" % python_lib)

def _check_python_bin(repository_ctx, python_bin):
    """Checks the python bin path."""
    cmd = '[[ -x "%s" ]] && [[ ! -d "%s" ]]' % (python_bin, python_bin)
    result = raw_exec(repository_ctx, [get_bash_bin(repository_ctx), "-c", cmd])
    if result.return_code == 1:
        auto_config_fail("--define %s='%s' is not executable. Is it the python binary?" % (
            PYTHON_BIN_PATH,
            python_bin,
        ))

def _get_python_include(repository_ctx, python_bin):
    """Gets the python include path."""
    result = execute(
        repository_ctx,
        [
            python_bin,
            "-Wignore",
            "-c",
            "import importlib; " +
            "import importlib.util; " +
            "print(importlib.import_module('distutils.sysconfig').get_python_inc() " +
            "if importlib.util.find_spec('distutils.sysconfig') " +
            "else importlib.import_module('sysconfig').get_path('include'))",
        ],
        error_msg = "Problem getting python include path.",
        error_details = ("Is the Python binary path set up right? " +
                         "(See ./configure or " + PYTHON_BIN_PATH + ".) " +
                         "Is distutils installed?"),
    )
    return result.stdout.splitlines()[0]

def _get_python_import_lib_name(repository_ctx, python_bin):
    """Get Python import library name (pythonXY.lib) on Windows."""
    result = execute(
        repository_ctx,
        [
            python_bin,
            "-c",
            "import sys;" +
            'print("python" + str(sys.version_info[0]) + ' +
            '      str(sys.version_info[1]) + ".lib")',
        ],
        error_msg = "Problem getting python import library.",
        error_details = ("Is the Python binary path set up right? " +
                         "(See ./configure or " + PYTHON_BIN_PATH + ".) "),
    )
    return result.stdout.splitlines()[0]
