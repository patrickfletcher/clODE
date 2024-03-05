from __future__ import annotations

import re


def _parse_ode_definition(xpp_line: str) -> list[tuple[str, str]]:
    line = re.sub("(par|p|init|aux) ", "", xpp_line)

    line = "=".join([a.strip() for a in line.split("=")])

    variable_pairs: list[tuple[str, str]] = []
    while line:
        first_variable = re.search(r"(\w+)(\s*=\s*)", line)
        if not first_variable:
            raise ValueError(f"Could not parse line '{xpp_line}'")
        name = first_variable.group(1)
        end_index = first_variable.end(0)
        second_variable = re.search(r"(\w+)(\s*=\s*)", line[end_index:])
        if second_variable:
            value = line[end_index : end_index + second_variable.start(0)]
            line = line[end_index + second_variable.start(0) :]
        else:
            value = line[end_index:]
            line = ""
        variable_pairs.append((name, value))

    return variable_pairs


def _parse_ode_definitions(store: dict[str, str], xpp_line: str) -> None:
    for name, value in _parse_ode_definition(xpp_line):
        store[name] = value.rstrip(",")


def read_ode_parameters(
    xpp_string: str,
) -> tuple[
    dict[str, str], dict[str, str], dict[str, str], dict[str, str], list[str], list[str]
]:
    parameters: dict[str, str] = dict()
    auxiliaries: dict[str, str] = dict()
    initial_values: dict[str, str] = dict()

    dx: dict[str, str] = dict()

    noise: list[str] = list()

    statements: list[str] = list()

    match_dx = "\w+'\s*=\s*"  # Match the name of the variable

    for line in xpp_string.splitlines():
        if line.startswith("par ") or line.startswith("p "):
            _parse_ode_definitions(parameters, line)
        elif line.startswith("aux "):
            _parse_ode_definitions(auxiliaries, line)
        elif line.startswith("init "):
            _parse_ode_definitions(initial_values, line)
        elif line.startswith("wiener "):
            noise.append(line[7:])
        elif re.match(match_dx, line):
            dx[line.split("'")[0]] = re.sub(match_dx, "", line)
        elif line.startswith("@"):
            pass  # TODO - handle configuration options
        elif line.strip() == "":
            pass  # Ignore empty lines
        else:
            statements.append(line)

    return parameters, auxiliaries, initial_values, dx, noise, statements


def format_opencl_rhs(
    parameters: dict[str, str],
    auxiliaries: dict[str, str],
    initial_values: dict[str, str],
    dx: dict[str, str],
    noise: list[str],
    statements: list[str],
) -> str:
    cl_file = """void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {

"""

    cl_file += "    /* State variables */\n"
    for index, name in enumerate(initial_values.keys()):
        cl_file += f"    realtype {name} = x_[{index}];\n"

    cl_file += "\n"

    cl_file += "    /* Parameters */\n"
    for index, name in enumerate(parameters):
        cl_file += f"    realtype {name} = p_[{index}];\n"

    cl_file += "\n"

    cl_file += "    /* Noise terms */\n"
    for index, name in enumerate(noise):
        cl_file += f"    realtype {name} = w_[{index}];\n"

    cl_file += "\n"

    cl_file += "    /* Core equations */\n"
    for line in statements:
        if line.find("=") != -1:
            non_ws_index = len(line) - len(line.lstrip())
            line = f"{line[:non_ws_index]}realtype {line[non_ws_index:]}"
        cl_file += f"    {line};\n"

    cl_file += "\n"

    cl_file += "    /* Auxiliary equations */\n"
    for name, value in auxiliaries.items():
        cl_file += f"    realtype {name} = {value};\n"

    cl_file += "\n"

    cl_file += "    /* Differential equations */\n"
    for name, value in dx.items():
        cl_file += f"    realtype d{name} = {value};\n"

    cl_file += "\n"

    cl_file += "    /* Auxiliary outputs */\n"
    for index, name in enumerate(auxiliaries.keys()):
        cl_file += f"    aux_[{index}] = {name};\n"

    cl_file += "\n"

    cl_file += "    /* Differential outputs */\n"

    for index, name in enumerate(initial_values.keys()):
        cl_file += f"    dx_[{index}] = d{name};\n"

    cl_file += "}"

    cl_file = re.sub("(\w+)\s*\^\s*2", r"\1*\1", cl_file, flags=re.MULTILINE)
    cl_file = re.sub("(\w+)\s*\^\s*3", r"\1*\1*\1", cl_file, flags=re.MULTILINE)
    cl_file = re.sub("(\w+)\s*\^\s*4", r"\1*\1*\1*\1", cl_file, flags=re.MULTILINE)
    cl_file = re.sub(
        "(\w+)\s*\^\s*([0-9]+)", r"pown(\1, \2)", cl_file, flags=re.MULTILINE
    )
    cl_file = re.sub(
        "(\w+)\s*\^\s*([-+]?(\d*\.*\d+))", r"pow(\1, \2)", cl_file, flags=re.MULTILINE
    )

    # Convert all floating numbers to single precision
    # using the suffix 'f'
    cl_file = re.sub("([-+]?(\d+\.\d*))", r"\1f", cl_file, flags=re.MULTILINE)

    cl_lines = cl_file.split("\n")

    cl_commented_lines: list[str] = []
    for line in cl_lines:
        comment_index = line.find("%")
        if comment_index != -1:
            commented_line = f"{line[:comment_index]}/* {line[comment_index+1:]} */"
            cl_commented_lines.append(commented_line)
        else:
            cl_commented_lines.append(line)

    return "\n".join(cl_commented_lines)


def convert_xpp_file(filename: str) -> str:
    with open(filename, "r") as f:
        xpp_string = f.read()

    (
        parameters,
        auxiliaries,
        initial_values,
        dx,
        noise,
        statements,
    ) = read_ode_parameters(xpp_string)

    cl_filename = filename[:-4] + ".cl"
    with open(cl_filename, "w") as f:
        f.write(
            format_opencl_rhs(
                parameters, auxiliaries, initial_values, dx, noise, statements
            )
        )
    return cl_filename
