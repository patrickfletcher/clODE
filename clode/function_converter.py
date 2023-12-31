import ast
import inspect
import textwrap
import typing
from enum import Enum


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
        return f"{str(self.value)}f"

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
    args: list[OpenCLExpression]
    cl_type: OpenCLType

    def __init__(
        self, name: str, args: list[OpenCLExpression], cl_type: OpenCLType
    ) -> None:
        self.name = name
        self.args = args
        self.cl_type = cl_type

    def __str__(self) -> str:
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
        self, subscript: ast.Subscript, context: dict[str, OpenCLType]
    ) -> None:
        if not isinstance(subscript.value, ast.Name):
            raise ValueError(
                f"Array access must be a variable at line {subscript.lineno}"
            )
        elif not isinstance(subscript.slice, ast.Constant):
            raise ValueError(
                f"Array access must be an index at line {subscript.lineno}"
            )

        self.name = subscript.value.id
        self.index = _convert_ast_expression_to_cl_expression(
            subscript.slice.value, context
        )
        # Return underlying type of self.name from context
        underlying_type = context[self.name]
        self.cl_type = OpenCLType(underlying_type.name, False)

    # def __init__(
    #     self, name: str, index: OpenCLExpression, context: dict[str, OpenCLType]
    # ) -> None:
    #     self.name = name
    #     self.index = index
    #     # Return underlying type of self.name from context
    #     underlying_type = context[self.name]
    #     self.cl_type = OpenCLType(underlying_type.name, False)

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
    expression: ast.expr, context: dict[str, OpenCLType]
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
        if isinstance(expression.func, ast.Name):
            if expression.func.id not in context:
                raise ValueError(
                    f"Function '{expression.func.id}' not found at line {expression.lineno}"
                )
            return_type = context[expression.func.id]
        else:
            raise ValueError(
                f"Unsupported function call '{expression.func}' at line {expression.lineno}"
            )
        return OpenCLFunctionCall(expression.func.id, args, cl_type=return_type)
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
    name, annotation: ast.Name | ast.Subscript | ast.BinOp | ast.Constant | None
) -> OpenCLType:
    array = False
    if isinstance(annotation, ast.Subscript):
        if isinstance(annotation.value, ast.Name):
            if annotation.value.id == "list":
                array = True
            else:
                raise TypeError(
                    f"Variable '{name}' must be a list at line {annotation.lineno}"
                )
        else:
            raise TypeError(
                f"Variable '{name}' must be a list at line {annotation.lineno}"
            )
        if isinstance(annotation.slice, ast.Name):
            if annotation.slice.id == "float":
                return OpenCLType("realtype", array)
            elif annotation.slice.id == "int":
                return OpenCLType("int", array)
            else:
                raise TypeError(
                    f"Variable '{name}' must be a list of floats or ints at line {annotation.lineno}"
                )
        else:
            raise TypeError(
                f"Variable '{name}' must be a list at line {annotation.lineno}"
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
    const: bool = False

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
        if isinstance(arg.annotation, (ast.Name | ast.Subscript | ast.BinOp)):
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
    target: OpenCLArrayAccess | OpenCLVariable | None = None
    cl_type: OpenCLType | None = None
    expression_type: OpenCLExpressionType
    expression: OpenCLExpression | None = None

    @property
    def target_name(self) -> str | None:
        if self.target is not None:
            return self.target.name
        else:
            return None

    def __init__(
        self,
        fn_name: str,
        instruction: ast.Assign | ast.AnnAssign | ast.Return,
        context: dict[str, OpenCLType],
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
                if isinstance(instruction.annotation, (ast.Name | ast.Subscript)):
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
    args: list[OpenCLArgument]
    body: list[OpenCLInstruction]
    returns: OpenCLType
    declared_vars: dict[str, OpenCLType]

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
        body: list[ast.stmt],
        fn_returns: ast.Name | ast.Constant,
        context: dict[str, OpenCLType],
    ):
        self.args = []
        self.body = []
        self.name = fn_name
        self.declared_vars = {}
        self.returns = _convert_ast_annotation_to_cl_type(fn_name, fn_returns)
        self._convert_ast_args(args)
        if isinstance(body, list):
            for instruction in body:
                local_context = dict(context, **self.declared_vars)
                if not isinstance(instruction, (ast.Assign, ast.AnnAssign, ast.Return)):
                    raise TypeError(
                        f"Unsupported instruction type {type(instruction)} at line {instruction.lineno}"
                    )
                cl_arg = OpenCLInstruction(fn_name, instruction, context=local_context)

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


class OpenCLSyntaxTree:
    functions: list[OpenCLFunction]

    def __init__(self) -> None:
        self.functions = []

    def add_function(self, fn: ast.FunctionDef):
        if fn.returns is None:
            raise TypeError(
                f"Function '{fn.name}' must have a return type at line {fn.lineno}"
            )
        elif not isinstance(fn.returns, (ast.Name, ast.Constant)):
            raise TypeError(
                f"Unsupported return type {type(fn.returns)} for function '{fn.name}' at line {fn.lineno}"
            )
        context: dict[str, OpenCLType] = {
            parsed_fn.name: parsed_fn.returns for parsed_fn in self.functions
        }
        self.functions.append(
            OpenCLFunction(fn.name, fn.args, fn.body, fn.returns, context)
        )

    def __str__(self) -> str:
        syntax_tree: str = ""
        for fn in self.functions:
            syntax_tree += f"{str(fn)}\n"
        return syntax_tree


class OpenCLConverter(ast.NodeTransformer):
    entry_function_seen: bool = False
    entry_function_name: str
    syntax_tree: OpenCLSyntaxTree

    def __init__(self, entry_function_name: str = "get_rhs"):
        # Initialize any necessary variables
        self.entry_function_name = entry_function_name
        self.syntax_tree = OpenCLSyntaxTree()

    def visit_FunctionDef(self, node):
        # Change the function signature for OpenCL kernel
        # For example, change 'def' to '__kernel void'
        # Add '__global' keyword for pointer arguments
        # Change argument types to 'realtype'
        # More modifications can be added here
        entrypoint = node.name == self.entry_function_name
        if entrypoint:
            entry_function_seen = True

        self.syntax_tree.add_function(node)
        return node

    def visit_BinOp(self, node):
        # Convert binary operations to OpenCL syntax if needed
        # Example: Python's power operator '**' to OpenCL's 'pow' function
        self.generic_visit(node)
        return node

    def convert_to_opencl(self, python_fn: typing.Callable, dedent: bool = True):
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
        python_source = inspect.getsource(python_fn)
        if dedent:
            python_source = textwrap.dedent(python_source)
        tree = ast.parse(python_source)
        self.visit(tree)
        return str(self.syntax_tree)


def convert_str_to_opencl(python_code):
    tree = ast.parse(python_code)
    converter = OpenCLConverter()
    converter.visit(tree)

    return str(converter.syntax_tree)
