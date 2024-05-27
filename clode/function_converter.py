from __future__ import annotations

import ast
import inspect
import sys
import textwrap
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Union

OpenCLRhsEquation = Callable[
    [float, List[float], List[float], List[float], List[float], List[float]], None
]


class OpenCLType:
    name: str
    array: bool = False

    def __init__(self, name: str, array: bool = False) -> None:
        if name not in ["int", "realtype", "bool", "void"]:
            raise ValueError(f"Unsupported type '{name}'")
        self.name = name
        self.array = array

    def __str__(self) -> str:
        return self.name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, OpenCLType):
            return NotImplemented
        return self.name == other.name and self.array == other.array

    def __neq__(self, other: object) -> bool:
        return not self.__eq__(other)

    def casting_allowed(self, other: "OpenCLType") -> bool:
        if self.name == other.name:
            return self.array == other.array
        if self.name == "realtype" and other.name == "int":
            return self.array == other.array
        return False


class OpenCLBuiltin:
    arg_count: int
    name: str
    return_type: OpenCLType

    def __init__(self, arg_count: int, name: str, return_type: OpenCLType) -> None:
        self.arg_count = arg_count
        self.name = name
        self.return_type = return_type


_opencl_builtins = {
    "acos": OpenCLBuiltin(1, "acos", OpenCLType("realtype")),
    "acosh": OpenCLBuiltin(1, "acosh", OpenCLType("realtype")),
    "acospi": OpenCLBuiltin(1, "acospi", OpenCLType("realtype")),
    "asin": OpenCLBuiltin(1, "asin", OpenCLType("realtype")),
    "asinh": OpenCLBuiltin(1, "asinh", OpenCLType("realtype")),
    "asinpi": OpenCLBuiltin(1, "asinpi", OpenCLType("realtype")),
    "atan": OpenCLBuiltin(1, "atan", OpenCLType("realtype")),
    "atan2": OpenCLBuiltin(2, "atan2", OpenCLType("realtype")),
    "atanh": OpenCLBuiltin(1, "atanh", OpenCLType("realtype")),
    "atanpi": OpenCLBuiltin(1, "atanpi", OpenCLType("realtype")),
    "atan2pi": OpenCLBuiltin(2, "atan2pi", OpenCLType("realtype")),
    "cbrt": OpenCLBuiltin(1, "cbrt", OpenCLType("realtype")),
    "ceil": OpenCLBuiltin(1, "ceil", OpenCLType("realtype")),
    "copysign": OpenCLBuiltin(2, "copysign", OpenCLType("realtype")),
    "cos": OpenCLBuiltin(1, "cos", OpenCLType("realtype")),
    "cosh": OpenCLBuiltin(1, "cosh", OpenCLType("realtype")),
    "cospi": OpenCLBuiltin(1, "cospi", OpenCLType("realtype")),
    "erfc": OpenCLBuiltin(1, "erfc", OpenCLType("realtype")),
    "erf": OpenCLBuiltin(1, "erf", OpenCLType("realtype")),
    "exp": OpenCLBuiltin(1, "exp", OpenCLType("realtype")),
    "exp2": OpenCLBuiltin(1, "exp2", OpenCLType("realtype")),
    "exp10": OpenCLBuiltin(1, "exp10", OpenCLType("realtype")),
    "expm1": OpenCLBuiltin(1, "expm1", OpenCLType("realtype")),
    "abs": OpenCLBuiltin(1, "fabs", OpenCLType("realtype")),
    "fabs": OpenCLBuiltin(1, "fabs", OpenCLType("realtype")),
    "fdim": OpenCLBuiltin(2, "fdim", OpenCLType("realtype")),
    "fmod": OpenCLBuiltin(2, "fmod", OpenCLType("realtype")),
    "floor": OpenCLBuiltin(1, "floor", OpenCLType("realtype")),
    # "fma": OpenCLBuiltin(3, "fma", OpenCLType("realtype")), # OpenCL returns nan
    "max": OpenCLBuiltin(2, "fmax", OpenCLType("realtype")),
    "min": OpenCLBuiltin(2, "fmin", OpenCLType("realtype")),
    "mod": OpenCLBuiltin(2, "fmod", OpenCLType("realtype")),
    # "fract": OpenCLBuiltin(1, "fract", OpenCLType("realtype")), # OpenCL Compiler error
    # "frexp": OpenCLBuiltin(2, "frexp", OpenCLType("realtype")),
    "heaviside": OpenCLBuiltin(1, "heaviside", OpenCLType("realtype")),
    "hypot": OpenCLBuiltin(2, "hypot", OpenCLType("realtype")),
    "ilogb": OpenCLBuiltin(1, "ilogb", OpenCLType("int")),
    "ldexp": OpenCLBuiltin(2, "ldexp", OpenCLType("realtype")),
    "lgamma": OpenCLBuiltin(1, "lgamma", OpenCLType("realtype")),
    # "lgamma_r": OpenCLBuiltin(2, "lgamma_r", OpenCLType("realtype")), # Pointers not supported
    "log": OpenCLBuiltin(1, "log", OpenCLType("realtype")),
    "log2": OpenCLBuiltin(1, "log2", OpenCLType("realtype")),
    "log10": OpenCLBuiltin(1, "log10", OpenCLType("realtype")),
    "log1p": OpenCLBuiltin(1, "log1p", OpenCLType("realtype")),
    # "logb": OpenCLBuiltin(1, "logb", OpenCLType("realtype")), # OpenCL returns nan
    # "mad": OpenCLBuiltin(3, "mad", OpenCLType("realtype")), # OpenCL returns nan
    # "modf": OpenCLBuiltin(2, "modf", OpenCLType("realtype")), # Pointers not supported
    # "nan": OpenCLBuiltin(0, "nan", OpenCLType("realtype")),
    "nextafter": OpenCLBuiltin(2, "nextafter", OpenCLType("realtype")),
    "pow": OpenCLBuiltin(2, "pow", OpenCLType("realtype")),
    "pown": OpenCLBuiltin(2, "pown", OpenCLType("realtype")),
    "powr": OpenCLBuiltin(2, "powr", OpenCLType("realtype")),
    "remainder": OpenCLBuiltin(2, "remainder", OpenCLType("realtype")),
    # "remquo": OpenCLBuiltin(3, "remquo", OpenCLType("realtype")), # Pointers not supported
    "rint": OpenCLBuiltin(1, "rint", OpenCLType("realtype")),
    "rootn": OpenCLBuiltin(2, "rootn", OpenCLType("realtype")),
    "round": OpenCLBuiltin(1, "round", OpenCLType("realtype")),
    "rsqrt": OpenCLBuiltin(1, "rsqrt", OpenCLType("realtype")),
    "sin": OpenCLBuiltin(1, "sin", OpenCLType("realtype")),
    "sincos": OpenCLBuiltin(2, "sincos", OpenCLType("realtype")),
    "sinh": OpenCLBuiltin(1, "sinh", OpenCLType("realtype")),
    "sinpi": OpenCLBuiltin(1, "sinpi", OpenCLType("realtype")),
    "sqrt": OpenCLBuiltin(1, "sqrt", OpenCLType("realtype")),
    "tan": OpenCLBuiltin(1, "tan", OpenCLType("realtype")),
    "tanh": OpenCLBuiltin(1, "tanh", OpenCLType("realtype")),
    "tanpi": OpenCLBuiltin(1, "tanpi", OpenCLType("realtype")),
    "gamma": OpenCLBuiltin(1, "tgamma", OpenCLType("realtype")),
    "trunc": OpenCLBuiltin(1, "trunc", OpenCLType("realtype")),
}


def _convert_ast_op_to_cl_op(op: ast.operator) -> str:
    if isinstance(op, ast.Add):
        return "+"
    elif isinstance(op, ast.Sub):
        return "-"
    elif isinstance(op, ast.Mult):
        return "*"
    elif isinstance(op, ast.Div):
        return "/"
    elif isinstance(op, ast.Pow):
        return "pow"
    elif isinstance(op, ast.Mod):
        return "%"
    else:
        raise ValueError(f"Unsupported operator '{op}'  at line {op.lineno}")


class OpenCLExpression:
    def get_cl_type(self) -> OpenCLType:
        raise NotImplementedError()


class OpenCLInteger(OpenCLExpression):
    value: int

    def __init__(self, value: int) -> None:
        self.value = value

    def __str__(self) -> str:
        return f"{str(self.value)}"

    def get_cl_type(self) -> OpenCLType:
        return OpenCLType("int")


class OpenCLFloat(OpenCLExpression):
    value: float

    def __init__(self, value: float) -> None:
        self.value = value

    def __str__(self) -> str:
        return f"RCONST({str(self.value)})"

    def get_cl_type(self) -> OpenCLType:
        return OpenCLType("realtype")


class OpenCLVariable(OpenCLExpression):
    name: str
    cl_type: OpenCLType

    def __init__(self, name: str, cl_type: OpenCLType) -> None:
        self.name = name
        self.cl_type = cl_type

    def __str__(self) -> str:
        return self.name

    def get_cl_type(self) -> OpenCLType:
        return self.cl_type


class OpenCLFunctionCall(OpenCLExpression):
    name: str
    args: List[OpenCLExpression]
    cl_type: OpenCLType

    def __init__(
        self, name: str, args: List[OpenCLExpression], cl_type: OpenCLType
    ) -> None:
        if name in _opencl_builtins:
            # Verify that the number of arguments is correct
            if len(args) != _opencl_builtins[name].arg_count:
                raise ValueError(
                    (
                        f"Invalid number of arguments for function '{name}', expected "
                        f"'{_opencl_builtins[name].arg_count}', args = {args}"
                        "\n    Hint: OpenCL Arg count can differ from Python arg count."
                    )
                )
            self.name = _opencl_builtins[name].name
        else:
            self.name = name
        self.args = args
        self.cl_type = cl_type

    def __str__(self) -> str:
        if self.name == "int" or self.name == "float":
            name = "int" if self.name == "int" else "realtype"
            return f"({name})({', '.join(map(str, self.args))})"
        return f"{self.name}({', '.join(map(str, self.args))})"

    def get_cl_type(self) -> OpenCLType:
        return self.cl_type


class OpenCLUnaryOperation(OpenCLExpression):
    arg: OpenCLExpression
    op: str

    def __init__(self, op: ast.unaryop, arg: OpenCLExpression) -> None:
        self.arg = arg
        if isinstance(op, ast.USub):
            self.op = "-"
        else:
            raise ValueError(f"Unsupported unary operator {op} at line {op.lineno}")

    def __str__(self) -> str:
        return f"({self.op}{self.arg})"

    def get_cl_type(self) -> OpenCLType:
        return self.arg.get_cl_type()


class OpenCLArrayAccess(OpenCLExpression):
    name: str
    index: OpenCLExpression
    cl_type: OpenCLType

    def __init__(
        self, subscript: ast.Subscript, context: Dict[str, OpenCLType]
    ) -> None:
        if not isinstance(subscript.value, ast.Name):
            raise ValueError(
                f"Array access must be a variable at line {subscript.lineno}"
            )
        elif sys.version_info < (3, 9):
            if not isinstance(subscript.slice, ast.Index):
                raise ValueError(
                    f"Array access must be an index at line {subscript.lineno}, got '{type(subscript.slice)}'"
                )
            else:
                subscript_slice = subscript.slice.value
        else:
            subscript_slice = subscript.slice

        if not isinstance(subscript_slice, ast.Constant):
            raise ValueError(
                f"Array access must be an index at line {subscript.lineno}, got '{type(subscript.slice)}'"
            )

        self.name = subscript.value.id
        self.index = _convert_ast_expression_to_cl_expression(
            subscript_slice.value, context
        )
        # Return underlying type of self.name from context
        underlying_type = context[self.name]
        self.cl_type = OpenCLType(underlying_type.name, False)

    def __str__(self) -> str:
        return f"{self.name}[{self.index}]"

    def get_cl_type(self) -> OpenCLType:
        return self.cl_type


class OpenCLBinaryOperation(OpenCLExpression):
    left: OpenCLExpression
    right: OpenCLExpression
    op: str

    def __init__(
        self, left: OpenCLExpression, op: ast.operator, right: OpenCLExpression
    ) -> None:
        self.left = left
        self.right = right
        self.op = _convert_ast_op_to_cl_op(op)

    def __str__(self) -> str:
        # If the operation is a power, handle the pow special case
        if self.op == "pow":
            # If the exponent is an integer < 4, unroll the loop
            if isinstance(self.right, OpenCLInteger):
                if self.right.value == 0:
                    return "1"
                elif self.right.value == 1:
                    return str(self.left)
                elif self.right.value == 2:
                    return f"({self.left} * {self.left})"
                elif self.right.value == 3:
                    return f"({self.left} * {self.left} * {self.left})"
                elif self.right.value == 4:
                    return f"({self.left} * {self.left} * {self.left} * {self.left})"
                else:
                    return f"pown({self.left}, {self.right})"
            # Handle the case when the exponent is a float by using pow
            elif isinstance(self.right, OpenCLFloat):
                return f"pow({self.left}, {self.right})"
            else:
                # Get the type of the exponent and use pown or pow accordingly
                exponent_type = self.right.get_cl_type()
                if exponent_type.name == "int":
                    return f"pown({self.left}, {self.right})"
                elif exponent_type.name == "realtype":
                    return f"pow({self.left}, {self.right})"
                else:
                    raise TypeError(f"Invalid type '{exponent_type}' for exponent")
        else:
            return f"({self.left} {self.op} {self.right})"

    def get_cl_type(self) -> OpenCLType:
        # Returns the free-est type of the two operands
        left_type = self.left.get_cl_type()
        right_type = self.right.get_cl_type()
        if left_type == OpenCLType("realtype") or right_type == OpenCLType("realtype"):
            return OpenCLType("realtype")
        elif left_type == OpenCLType("int") or right_type == OpenCLType("int"):
            return OpenCLType("int")
        else:
            raise ValueError(
                f"Invalid types '{left_type}' and '{right_type}' for binary operation '{self.op}'"
            )


def _convert_ast_expression_to_cl_expression(
    expression: ast.expr, context: Dict[str, OpenCLType]
) -> OpenCLExpression:
    if isinstance(expression, int):
        return OpenCLInteger(expression)
    elif isinstance(expression, float):
        return OpenCLFloat(expression)
    elif isinstance(expression, ast.Constant):
        if isinstance(expression.value, int):
            return OpenCLInteger(expression.value)
        elif isinstance(expression.value, float):
            return OpenCLFloat(expression.value)
        else:
            raise ValueError(
                f"Unsupported constant '{expression.value}' at line {expression.lineno}"
            )
    elif isinstance(expression, ast.Name):
        # Check if variable exists in context and retrieve its type
        if expression.id not in context:
            raise ValueError(
                f"Variable '{expression.id}' not found at line {expression.lineno}"
            )
        cl_type = context[expression.id]
        return OpenCLVariable(expression.id, cl_type)
    elif isinstance(expression, ast.Call):
        args = [
            _convert_ast_expression_to_cl_expression(arg, context)
            for arg in expression.args
        ]
        # Check if function exists in context and retrieve its type
        function_name: str
        if isinstance(expression.func, ast.Name):
            if expression.func.id == "int":
                return_type = OpenCLType("int")
            elif expression.func.id == "float":
                return_type = OpenCLType("realtype")
            elif expression.func.id not in context:
                raise ValueError(
                    f"Function '{expression.func.id}' not found at line {expression.lineno}"
                )
            else:
                return_type = context[expression.func.id]
            function_name = expression.func.id
        # Check if func is a module (e.g. math.exp)
        elif isinstance(expression.func, ast.Attribute):
            if expression.func.attr not in context:
                module = (
                    expression.func.value.id
                    if isinstance(expression.func.value, ast.Name)
                    else ""
                )
                raise ValueError(
                    f"Function '{expression.func.attr}' from module '{module}' not found at line {expression.lineno}"
                )
            return_type = context[expression.func.attr]
            function_name = expression.func.attr
        else:
            raise ValueError(
                f"Unsupported function call '{expression.func}' at line {expression.lineno}"
            )
        return OpenCLFunctionCall(function_name, args, cl_type=return_type)
    elif isinstance(expression, ast.UnaryOp):
        return OpenCLUnaryOperation(
            expression.op,
            _convert_ast_expression_to_cl_expression(expression.operand, context),
        )
    elif isinstance(expression, ast.BinOp):
        return OpenCLBinaryOperation(
            _convert_ast_expression_to_cl_expression(expression.left, context),
            expression.op,
            _convert_ast_expression_to_cl_expression(expression.right, context),
        )
    elif isinstance(expression, ast.Subscript):
        return OpenCLArrayAccess(expression, context)
    else:
        raise ValueError(
            f"Unsupported expression '{expression}' at line {expression.lineno}"
        )


def _convert_ast_annotation_to_cl_type(
    name: Optional[str],
    annotation: Union[ast.Name, ast.Subscript, ast.BinOp, ast.Constant, None],
) -> OpenCLType:
    array = False
    if isinstance(annotation, ast.Subscript):
        if isinstance(annotation.value, ast.Name):
            if annotation.value.id == "list" or annotation.value.id == "List":
                array = True
            else:
                raise TypeError(
                    f"Variable '{name}' must be a list at line {annotation.lineno}"
                )
        else:
            raise TypeError(
                f"Variable '{name}' must be a list at line {annotation.lineno}"
            )

        annotation_slice: ast.expr
        # In Python 3.9+, annotation.slice is an ast.Name
        if sys.version_info < (3, 9):
            if isinstance(annotation.slice, ast.Index):
                annotation_slice = annotation.slice.value
            else:
                raise TypeError(
                    f"Variable '{name}' must be a list of strings or float at line {annotation.lineno}, got 'List[{type(annotation.slice)}]'"
                )
        else:
            annotation_slice = annotation.slice

        if isinstance(annotation_slice, ast.Name):
            if annotation_slice.id == "float":
                return OpenCLType("realtype", array)
            elif annotation_slice.id == "int":
                return OpenCLType("int", array)
            else:
                raise TypeError(
                    f"Variable '{name}' must be a list of floats or ints at line {annotation.lineno}, got '{annotation_slice.id}'"
                )
        else:
            raise TypeError(
                f"Variable '{name}' must be a list of strings or float at line {annotation.lineno}, got 'List[{type(annotation_slice)}]'"
            )
    elif isinstance(annotation, ast.Name):
        if annotation.id == "int":
            return OpenCLType("int", array)
        elif annotation.id == "float":
            return OpenCLType("realtype", array)
        else:
            raise TypeError(
                f"Unsupported type '{annotation.id}' at line {annotation.lineno}"
            )
    elif isinstance(annotation, ast.BinOp):
        # Handle the case when the right argument type is None
        if isinstance(annotation.right, ast.Constant):
            if annotation.right.value is None:
                if isinstance(annotation.left, ast.Name) or isinstance(
                    annotation.left, ast.Subscript
                ):
                    return _convert_ast_annotation_to_cl_type(name, annotation.left)
        raise TypeError(
            f"Unsupported type for variable '{name}' at line {annotation.lineno}"
        )
    # Handle the case when the annotation is None
    elif isinstance(annotation, ast.Constant) and annotation.value is None:
        return OpenCLType("void", False)
    else:
        raise TypeError(f"Unsupported type for function '{name}'")


class OpenCLArgument:
    function: str
    name: str
    cl_type: OpenCLType
    const: bool

    def __str__(self) -> str:
        arg = f"{self.cl_type} {self.name}"

        if self.cl_type.array:
            arg += "[]"
        if self.const:
            arg = f"const {arg}"
        return arg

    def __init__(self, fn_name: str, arg: ast.arg, const: bool = False) -> None:
        self.function = fn_name
        self.name = arg.arg
        if isinstance(arg.annotation, (ast.Name, ast.Subscript, ast.BinOp)):
            self.cl_type = _convert_ast_annotation_to_cl_type(self.name, arg.annotation)
        # Handle the case when the annotation is None
        # Other types are caught in the _convert_ast_annotation_to_cl_type function
        else:
            raise TypeError(f"Argument '{self.name}' has no type at line {arg.lineno}")
        self.const = const


class OpenCLExpressionType(Enum):
    DECLARE = "declare"
    DECLARE_AND_ASSIGN = "declare_and_assign"
    ASSIGN = "assign"
    RETURN = "return"


class OpenCLInstruction:
    target: Union[OpenCLArrayAccess, OpenCLVariable, None] = None
    cl_type: Optional[OpenCLType] = None
    expression_type: OpenCLExpressionType
    expression: Optional[OpenCLExpression] = None

    @property
    def target_name(self) -> Optional[str]:
        if self.target is not None:
            return self.target.name
        else:
            return None

    def __init__(
        self,
        fn_name: str,
        instruction: Union[ast.Assign, ast.AnnAssign, ast.Return],
        context: Dict[str, OpenCLType],
    ) -> None:
        if isinstance(instruction, ast.Assign) or isinstance(
            instruction, ast.AnnAssign
        ):
            if isinstance(instruction, ast.Assign):
                target = instruction.targets[0]
                self.expression_type = OpenCLExpressionType.ASSIGN
                if len(instruction.targets) > 1 or isinstance(
                    instruction.targets[0], ast.Tuple
                ):
                    raise TypeError(
                        f"Cannot assign multiple variables in one line at line {instruction.lineno}"
                    )
                if isinstance(target, ast.Name):
                    target_name = target.id
                elif isinstance(target, ast.Subscript):
                    if isinstance(target.value, ast.Name):
                        target_name = target.value.id
                    else:
                        raise ValueError(
                            f"Variable not declared at line {instruction.lineno}"
                        )
                else:
                    raise ValueError(
                        f"Variable not declared at line {instruction.lineno}"
                    )
                if target_name not in context:
                    if isinstance(target, ast.Name):
                        raise ValueError(
                            f"Variable '{target.id}' not declared at line {instruction.lineno}"
                        )
                    raise ValueError(
                        f"Variable not declared at line {instruction.lineno}"
                    )
            else:
                target = instruction.target
                if isinstance(instruction.annotation, (ast.Name, ast.Subscript)):
                    self.cl_type = _convert_ast_annotation_to_cl_type(
                        self.target_name, instruction.annotation
                    )
                else:
                    raise TypeError(
                        f"Unsupported type for variable '{self.target_name}' at line {instruction.lineno}"
                    )
                # Check if there is a target and set the type to declare
                # or declare_and_assign
                if instruction.value is None:
                    self.expression_type = OpenCLExpressionType.DECLARE
                else:
                    self.expression_type = OpenCLExpressionType.DECLARE_AND_ASSIGN

            # Check if the target is an Ast.Name
            if isinstance(target, ast.Name):
                self.target = OpenCLVariable(target.id, OpenCLType("void"))
            elif isinstance(target, ast.Subscript):
                self.target = OpenCLArrayAccess(target, context)
            else:
                raise TypeError(
                    f"Unsupported target '{target}' at line {instruction.lineno}"
                )

        elif isinstance(instruction, ast.Return):
            self.target = None
            self.expression_type = OpenCLExpressionType.RETURN

        # Check if there is an expression and convert it
        if instruction.value is not None:
            self.expression = _convert_ast_expression_to_cl_expression(
                instruction.value, context
            )
            if self.cl_type is not None:
                if not self.cl_type.casting_allowed(self.expression.get_cl_type()):
                    raise ValueError(
                        f"Type mismatch for variable '{self.target_name}' at line {instruction.lineno}"
                    )
            self.cl_type = self.expression.get_cl_type()

    def __str__(self) -> str:
        if self.expression_type == OpenCLExpressionType.DECLARE:
            return f"{self.cl_type} {self.target}"
        elif self.expression_type == OpenCLExpressionType.DECLARE_AND_ASSIGN:
            return f"{self.cl_type} {self.target} = {self.expression}"
        elif self.expression_type == OpenCLExpressionType.ASSIGN:
            return f"{self.target} = {self.expression}"
        elif self.expression_type == OpenCLExpressionType.RETURN:
            return f"return {self.expression}"
        else:
            raise ValueError(f"Unsupported expression type {self.expression_type}")


class OpenCLFunction:
    name: str
    args: List[OpenCLArgument]
    body: List[OpenCLInstruction]
    returns: OpenCLType
    declared_vars: Dict[str, OpenCLType]
    modified_vars: set[str]

    def __str__(self) -> str:
        fn: str = f"{self.returns.name} {self.name}("
        arg_indent = len(fn)
        for arg in self.args:
            fn += f"{str(arg)},\n{arg_indent * ' '}"
        fn = fn[: -2 - arg_indent] + ") "
        fn += "{\n"
        for instruction in self.body:
            fn += f"    {str(instruction)};\n"
        fn += "}\n"

        return fn

    def _convert_ast_args(self, fn_args: ast.arguments) -> None:
        for arg in fn_args.args:
            cl_arg = OpenCLArgument(self.name, arg)
            # One cannot redeclare an argument in Python
            self.declared_vars[cl_arg.name] = cl_arg.cl_type
            self.args.append(cl_arg)

    def __init__(
        self,
        fn_name: str,
        args: ast.arguments,
        body: List[ast.stmt],
        fn_returns: Union[ast.Name, ast.Constant],
        context: Dict[str, OpenCLType],
        mutable_args: Optional[Union[List[str], List[int]]] = None,
    ):
        self.args = []
        self.body = []
        self.name = fn_name
        self.declared_vars = {}
        self.returns = _convert_ast_annotation_to_cl_type(fn_name, fn_returns)
        self._convert_ast_args(args)
        arg_names = set([arg.name for arg in self.args])

        if isinstance(mutable_args, List):
            mutable_args_str: List[str] = []
            # Check if mutable_args is a list of ints and convert to a list of strings
            if isinstance(mutable_args[0], int):
                # Treat all arguments as ints and throw an error if they are not
                for mutable_arg in mutable_args:
                    if not isinstance(mutable_arg, int):
                        raise TypeError(
                            f"Mutable arguments must be a list of strings or integers, got '{mutable_args}'"
                        )
                    mutable_args_str.append(self.args[mutable_arg].name)
            elif isinstance(mutable_args[0], str):
                # Treat all arguments as strings and throw an error if they are not
                for mutable_arg in mutable_args:
                    if not isinstance(mutable_arg, str):
                        raise TypeError(
                            f"Mutable arguments must be a list of strings or integers, got '{mutable_args}'"
                        )
                    mutable_args_str.append(mutable_arg)
            else:
                # Throw an error if the list is not a list of strings or integers
                raise TypeError(
                    f"Mutable arguments must be a list of strings or integers, got '{mutable_args}'"
                )

            self.modified_vars = set(mutable_args_str)
        elif mutable_args is None:
            self.modified_vars = set()
        else:
            raise TypeError(
                f"Mutable arguments must be a list of strings or integers, got '{mutable_args}'"
            )
        if isinstance(body, list):
            for instruction in body:
                local_context = dict(context, **self.declared_vars)
                if not isinstance(instruction, (ast.Assign, ast.AnnAssign, ast.Return)):
                    raise TypeError(
                        f"Unsupported instruction type {type(instruction)} at line {instruction.lineno}"
                    )
                cl_arg = OpenCLInstruction(fn_name, instruction, context=local_context)

                # Record the variables that are modified in the function
                if cl_arg.target_name in arg_names:
                    self.modified_vars.add(cl_arg.target_name)

                # Check if the variable is already declared and if the type is changed
                if cl_arg.target_name in self.declared_vars:
                    if (
                        cl_arg.expression_type == OpenCLExpressionType.DECLARE
                        or cl_arg.expression_type
                        == OpenCLExpressionType.DECLARE_AND_ASSIGN
                    ):
                        raise TypeError(
                            f"Cannot redeclare variable '{cl_arg.target_name}' at line {instruction.lineno}"
                        )
                    if cl_arg.cl_type is not None:
                        if isinstance(cl_arg.target, OpenCLArrayAccess):
                            if self.declared_vars[cl_arg.target_name] != OpenCLType(
                                name=cl_arg.cl_type.name, array=True
                            ):
                                raise TypeError(
                                    f"Type mismatch for variable '{cl_arg.target_name}' at line {instruction.lineno}"
                                )
                        elif self.declared_vars[cl_arg.target_name] != cl_arg.cl_type:
                            raise TypeError(
                                f"Cannot change the type of variable '{cl_arg.target_name}' at line {instruction.lineno}"
                            )
                elif not isinstance(instruction, ast.Return):
                    if cl_arg.target_name is not None and cl_arg.cl_type is not None:
                        self.declared_vars[cl_arg.target_name] = cl_arg.cl_type
                    else:
                        # Something has gone wrong here
                        raise TypeError(
                            f"Cannot assign to variable '{cl_arg.target_name}' at line {instruction.lineno}"
                        )
                self.body.append(cl_arg)
        else:
            raise TypeError(
                f"Body of function '{fn_name}' must be a list of instructions"
            )

        # Set all the variables that are not modified to const
        for arg in self.args:
            if arg.name not in self.modified_vars:
                arg.const = True


class OpenCLSyntaxTree:
    functions: List[OpenCLFunction]

    def __init__(self) -> None:
        self.functions = []

    def add_function(
        self,
        fn: ast.FunctionDef,
        mutable_args: Optional[Union[List[str], List[int]]],
        function_name: Optional[str],
    ) -> None:
        if fn.returns is None:
            raise TypeError(
                f"Function '{fn.name}' must have a return type at line {fn.lineno}"
            )
        elif not isinstance(fn.returns, (ast.Name, ast.Constant)):
            raise TypeError(
                f"Unsupported return type {type(fn.returns)} for function '{fn.name}' at line {fn.lineno}"
            )
        context: Dict[str, OpenCLType] = {
            parsed_fn.name: parsed_fn.returns for parsed_fn in self.functions
        }
        context.update({key: val.return_type for key, val in _opencl_builtins.items()})

        if function_name is None:
            function_name = fn.name
        self.functions.append(
            OpenCLFunction(
                function_name,
                fn.args,
                fn.body,
                fn.returns,
                context,
                mutable_args=mutable_args,
            )
        )

    def __str__(self) -> str:
        syntax_tree: str = ""
        for fn in self.functions:
            syntax_tree += f"{str(fn)}\n"
        return syntax_tree


class OpenCLConverter(ast.NodeTransformer):
    syntax_tree: OpenCLSyntaxTree
    mutable_args: Optional[Union[List[str], List[int]]] = None
    function_name: Optional[str] = None

    def __init__(self, entry_function_name: str = "getRHS"):
        # Initialize any necessary variables
        self.entry_function_name = entry_function_name
        self.syntax_tree = OpenCLSyntaxTree()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> ast.FunctionDef:
        # Change the function signature for OpenCL kernel
        # For example, change 'def' to '__kernel void'
        # Add '__global' keyword for pointer arguments
        # Change argument types to 'realtype'
        # More modifications can be added here

        mutable_args = self.mutable_args
        self.syntax_tree.add_function(
            node, mutable_args=mutable_args, function_name=self.function_name
        )
        return node

    def visit_BinOp(self, node: ast.BinOp) -> ast.BinOp:
        # Convert binary operations to OpenCL syntax if needed
        # Example: Python's power operator '**' to OpenCL's 'pow' function
        self.generic_visit(node)
        return node

    def convert_to_opencl(
        self,
        python_fn: Union[Callable[[Any], Any], OpenCLRhsEquation] | str,
        dedent: bool = True,
        mutable_args: Optional[Union[List[str], List[int]]] = None,
        function_name: Optional[str] = None,
    ) -> str:
        # Convert a Python function to OpenCL
        # Example: 'def add_float(a: float, b: float) -> float:\n'
        #          '    res: float = a + b\n'
        #          '    return res'
        # to
        #          'realtype add_float(realtype a, realtype b) {\n'
        #          '    realtype res = a + b;\n'
        #          '    return res;\n'
        #          '}\n'
        # More modifications can be added here
        if isinstance(python_fn, str):
            python_source = python_fn
        else:
            python_source = inspect.getsource(python_fn)
            if dedent:
                python_source = textwrap.dedent(python_source)
        self.mutable_args = mutable_args
        self.function_name = function_name
        tree = ast.parse(python_source)
        self.visit(tree)
        self.function_name = None
        self.mutable_args = None
        return str(self.syntax_tree)


def convert_str_to_opencl(python_code: str) -> str:
    tree = ast.parse(python_code)
    converter = OpenCLConverter()
    converter.visit(tree)

    return str(converter.syntax_tree)
