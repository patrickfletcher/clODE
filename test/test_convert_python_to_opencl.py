from __future__ import annotations

from typing import List

import pytest

from clode import OpenCLConverter, convert_str_to_opencl


# Test error handling
# All the below cases should raise ValueError
@pytest.mark.parametrize(
    "error_description, python_code",
    [
        [
            "Unsupported type 'str' at line 1",
            """def unsupported_type(x: str) -> float:
                return x
            """,
        ],
        [
            "Argument 'x' has no type at line 1",
            """def no_annotation(x) -> float:
                return x
            """,
        ],
        [
            "Cannot change the type of variable 'a'",
            """def change_variable_type() -> int:
                a: int = 1
                a = 2.2
                return a
            """,
        ],
        [
            "Function 'no_return_type' must have a return type at line 1",
            """def no_return_type():
                a: int = 1
                return a
            """,
        ],
        [
            "Unsupported type 'str' at line 1",
            """def return_type_str() -> str:
                a: int = 1
                return a
            """,
        ],
        [
            "Unsupported type 'str' at line 2",
            """def variable_string() -> int:
                a: str = 2
                return 1
            """,
        ],
        [
            "Cannot assign multiple variables in one line at line 2",
            """def tuple_assign() -> int:
                a, b = 1, 2
                return a + b
            """,
        ],
        [
            "Cannot redeclare variable 'a' at line 3",
            """def redeclare_variable() -> int:
                a: int = 1
                a: int = 2
                return a
            """,
        ],
    ],
)
def test_type_errors(error_description: str, python_code: str):
    with pytest.raises(TypeError, match=error_description):
        convert_str_to_opencl(python_code)
        print("Error handling test passed")


class TestConversions:
    def test_convert_adder(self):
        def add_float(a: float, b: float) -> float:
            res: float = a + b
            return res

        converter = OpenCLConverter()

        opencl_code = converter.convert_to_opencl(add_float)
        expect_code = (
            "realtype add_float(const realtype a,\n"
            "                   const realtype b) {\n"
            "    realtype res = (a + b);\n"
            "    return res;\n"
            "}\n\n"
        )
        assert expect_code == opencl_code

    def test_convert_deinded_adder(self):
        def add_floats_in_list(lst_in: List[float], lst_out: List[float]) -> None:
            res: float = lst_in[0] + lst_in[1]
            lst_out[0] = res

        converter = OpenCLConverter()

        opencl_code = converter.convert_to_opencl(add_floats_in_list)
        expect_code = (
            "void add_floats_in_list(const realtype lst_in[],\n"
            "                        realtype lst_out[]) {\n"
            "    realtype res = (lst_in[0] + lst_in[1]);\n"
            "    lst_out[0] = res;\n"
            "}\n\n"
        )
        assert expect_code == opencl_code

    def test_all_operations(self):
        def all_operations(a: float, b: float, c: int, d: int) -> float:
            res1: float = a + b * c / d
            res2: float = ((a - b) % c) ** d
            res3: float = a**b
            res4: float = a**1
            res5: float = a**2
            res6: float = a**3
            res7: float = a**4
            res8: float = a**5
            res9: float = a**0
            res10: float = a**0.5
            sum_res: float = (
                res1 + res2 + -res3 + res4 + res5 + res6 + res7 + res8 + res9 + res10
            )
            return res1 + sum_res

        converter = OpenCLConverter()

        opencl_code = converter.convert_to_opencl(all_operations)

        expect_code = (
            "realtype all_operations(const realtype a,\n"
            "                        const realtype b,\n"
            "                        const int c,\n"
            "                        const int d) {\n"
            "    realtype res1 = (a + ((b * c) / d));\n"
            "    realtype res2 = pown(((a - b) % c), d);\n"
            "    realtype res3 = pow(a, b);\n"
            "    realtype res4 = a;\n"
            "    realtype res5 = (a * a);\n"
            "    realtype res6 = (a * a * a);\n"
            "    realtype res7 = (a * a * a * a);\n"
            "    realtype res8 = pown(a, 5);\n"
            "    realtype res9 = 1;\n"
            "    realtype res10 = pow(a, RCONST(0.5));\n"
            "    realtype sum_res = (((((((((res1 + res2) + (-res3)) + res4) + res5) + res6) + res7) + res8) + res9) + res10);\n"
            "    return (res1 + sum_res);\n"
            "}\n\n"
        )
        assert expect_code == opencl_code

    def test_convert_second_function_referencing_first(self):
        def add_float(a: float, b: float) -> float:
            res: float = a + b
            return res

        def get_rhs(var: List[float], derivatives: List[float]) -> None:
            res: float = add_float(var[0], var[1])
            derivatives[0] = res

        converter = OpenCLConverter()
        converter.convert_to_opencl(add_float)
        opencl_code = converter.convert_to_opencl(get_rhs)
        expect_code = (
            "realtype add_float(const realtype a,\n"
            "                   const realtype b) {\n"
            "    realtype res = (a + b);\n"
            "    return res;\n"
            "}\n"
            "\n"
            "void get_rhs(const realtype var[],\n"
            "             realtype derivatives[]) {\n"
            "    realtype res = add_float(var[0], var[1]);\n"
            "    derivatives[0] = res;\n"
            "}\n\n"
        )
        assert expect_code == opencl_code
