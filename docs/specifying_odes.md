# Specifying ODEs

clODE supports many common use cases for specifying systems of ordinary
differential equations (ODEs). This document outlines the various ways to
specify ODEs in clODE.

## OpenCL

The ODE system can be specified directly as an OpenCL C-languange function, with the following signature:

``` c
void getRHS(const realtype t,
            const realtype variables[],
            const realtype parameters[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]);
```

## Python

To support conversion to OpenCL, an ODE system specified as a python function
must be **fully typed**, with the following signature:

```python
def getRHS(float t,
           list[float] variables,
           list[float] parameters,
           list[float] derivatives,
           list[float] aux,
           list[float] weiner) -> None
```

To support using the same function with Scipy's ```solve_ivp``, we provide a wrapper that
provides the signature expected there [TODO].

## Arithmetic operations

clODE supports the following arithmetic operations:

- Addition: `+`
- Subtraction: `-`
- Multiplication: `*`
- Division: `/`
- Exponentiation: `**`
- Modulo: `%`

```py run
from clode import OpenCLConverter

def get_rhs(t: float,
            variables: list[float],
            parameters: list[float],
            derivatives: list[float],
            aux: list[float],
            wiener: list[float]) -> None:
    x: float = variables[0]
    y: float = variables[1]
    
    derivatives[0] = x + y
    derivatives[1] = x - y
    aux[1] = x * y
    aux[2] = x / y
    aux[3] = x ** 3
    aux[4] = x % y
    
converter = OpenCLConverter()
print(converter.convert_to_opencl(get_rhs))
```

Note that the exponent operator `x ** y` is implicitly
converted to `x * ... * x` for integer exponents if x < 5.
If x >= 5, the exponent operator is converted to the
`pown` function.
If the exponent is a float, the `pow` function is used.

## Using builtins (sin, exp, etc.)

Certain builtin functions are available in clODE.
By convention, the floating point versions of
these functions are used.

The functions trigonometric currently available are:

| Sine   | Cosine | Tangent | Other    |
|--------|--------|---------|----------|
| sin    | cos    | tan     | atan2    |
| asin   | acos   | atan    | atan2pi  |
| sinh   | cosh   | tanh    |          |
| asinh  | acosh  | atanh   |          |
| sinpi  | cospi  | tanpi   |          |
| asinpi | acospi | atanpi  |          |

Additionally, the following functions are supported:

| Exponential & Logarithmic  | Power & Root  | Rounding & Remainder | Miscellaneous |
|----------------------------|---------------|----------------------|------------------------|
| exp                        | cbrt          | ceil                 | copysign               |
| exp2                       | pow           | floor                | fabs                   |
| exp10                      | pown          | fmod                 | fdim                   |
| expm1                      | powr          | remainder            | gamma                  |
| log                        | sqrt          | rint                 | hypot                  |
| log1p                      | rsqrt         | trunc                | ilogb                  |
| log2                       | rootn         |                      | ldexp                  |
| log10                      |               |                      | lgamma                 |
|                            |               |                      | nextafter              |
|                            |               |                      | erf                    |
|                            |               |                      | erfc                   |

You can change between int and float numbers by using
the Python builtins `int()` and `float()`.

**Note:** You can import arithmetic functions from
math, numpy or clode.

```py run
from clode import OpenCLConverter, exp

def get_rhs(
    t: float,
    variables: list[float],
    parameters: list[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    x: float = variables[0]
    y: float = variables[1]
    p1: int = int(parameters[0])

    derivatives[0] = (x - y) / exp(y)
    derivatives[1] = -x ** p1
    

converter = OpenCLConverter()
print(converter.convert_to_opencl(get_rhs))
```

## Specifying helper functions

It is possible to specify helper functions in Python and call them from the
main ODE function. This can be useful for breaking up complex ODEs into
smaller, more manageable pieces.

```py run
from clode import OpenCLConverter

def helper_function(x: float, y: float) -> float:
    return (x + y) / y

def get_rhs(t: float,
            variables: list[float],
            parameters: list[float],
            derivatives: list[float],
            aux: list[float],
            wiener: list[float]) -> None:
    x: float = variables[0]
    y: float = variables[1]
    
    derivatives[0] = helper_function(x, y)

converter = OpenCLConverter()
converter.convert_to_opencl(helper_function)
print(converter.convert_to_opencl(get_rhs))
```

To load a helper function into the OpenCL program, you must use
the list argument `supplementary_equations` in the Simulator initializers.

Note: The functions must be specified in the order that they are used.
clODE will not automatically resolve dependencies between functions.

```python
simulator = TrajectorySimulator(
    rhs_equation=get_rhs,
    variables={"x": 0.0, "y": 2.0},
    parameters={"a": 1.0},
    supplementary_equations=[helper_function],
    t_span=(0, 10),
    dt=0.01,
)
```
