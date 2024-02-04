import numpy as np
import clode
import time
from typing import List
import matplotlib.pyplot as plt

def lorenz(
    t: float,
    var: List[float],
    p: List[float],
    dvar: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    x: float = var[0]
    y: float = var[1]
    z: float = var[2]
    
    r: float = p[0],
    s: float = p[1],
    beta: float = p[2],

    dx: float = s*(y - x)
    dy: float = r*x - y - x*z
    dz: float = x*y - beta*z
    dvar[0] = dx
    dvar[1] = dy
    dvar[2] = dz

def test_lorenz_rk4(platform_id, device_id, num_pts, reps):
    
    parameters = {"r": 27.0, "s": 10.0, "b": 8.0/3.0}
    variables = {"x": 1.0, "y": 1.0, "z": 1.0}
    
    # convenience:
    variable_names = list(variables.keys())
    initial_state = list(variables.values())
    parameter_names = list(parameters.keys())
    default_parameters = list(parameters.values())

    t_span = (0.0, 10.0)

    src_file = "test/lorenz.cl"
    integrator = clode.Simulator(
        src_file=src_file,
        # rhs_equation=lorenz,
        variables=variables,
        parameters=parameters,
        single_precision=True,
        platform_id=platform_id,
        device_id=device_id,
        t_span=t_span,
        stepper=clode.Stepper.rk4,
        dt=0.01,
    )

    t_average = []
    print(f"Time for 1000 RK4 steps of the Lorenz system. {reps} repetitions.\n(N, min, median, max time)")
    for num in num_pts:
        times = []
        # Pars = np.tile(default_parameters, (num, 1))
        # X0 = np.tile(initial_state, (num, 1))
        # integrator.set_problem_data(X0, Pars)
        integrator.set_repeat_ensemble(num_repeats=num)

        #warm-up passes - seems to make output more reliable
        for _ in range(5):
            # integrator.set_problem_data(X0, Pars)
            integrator.transient(update_x0=False)
            integrator.get_final_state()

        for _ in range(reps):
            t0 = time.perf_counter()
            # integrator.initialize(X0, Pars) #include host->device transfers?
            integrator.transient(update_x0=False)
            # integrator.get_final_state() #include device->host transfers?
            times.append(time.perf_counter() - t0)

        t_mean = sum(times)/reps
        t_min = min(times)
        t_max = max(times)
        t_median = np.median(times)
        print(f"{num}\t {t_min:.3g}, {t_median:.3g}, {t_max:.3g} s")
        t_average.append(t_median)

    return t_average
    

# if using 'bazel test ...'
if __name__ == "__main__":
    # clode.set_log_level(clode.LogLevel.info)    
    # clode.runtime.print_opencl()
    # clode.set_log_level(clode.runtime.DEFAULT_LOG_LEVEL) 

    ocl_info = clode.runtime.query_opencl()
    ocl_info = ocl_info[:3]

    num_pts = [int(2**n) for n in np.arange(0,18)]
    reps = 50
    
    device_names = []
    for i, ocl in enumerate(ocl_info):
        if ocl.device_count==0:
            continue
        print("\n", ocl)
        for j, dev in enumerate(ocl.device_info):
            t_average = test_lorenz_rk4(platform_id=i, device_id=j, num_pts=num_pts, reps=reps)
            plt.plot(num_pts, t_average)
            device_names.append(dev.name)
    
    plt.legend(device_names)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()