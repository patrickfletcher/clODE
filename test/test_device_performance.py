import numpy as np
import clode
import time
import matplotlib.pyplot as plt


def test_lorenz_rk4(platform_id, device_id, num_pts, reps):
    
    parameters = {"r": 27.0, "s": 10.0, "b": 8.0/3.0}
    variables = {"x": 1.0, "y": 1.0, "z": 1.0}
    
    # convenience:
    variable_names = list(variables.keys())
    initial_state = list(variables.values())
    parameter_names = list(parameters.keys())
    default_parameters = list(parameters.values())

    tspan = (0.0, 100.0)

    src_file = "test/lorenz.cl"
    integrator = clode.Simulator(
        src_file=src_file,
        # rhs_equation=getRHS,
        variable_names=variable_names,
        parameter_names=parameter_names,
        aux=["aux1"],
        num_noise=1,
        single_precision=True,
        platform_id=platform_id,
        device_id=device_id,
        tspan=tspan,
        stepper=clode.Stepper.rk4,
        dt=0.1,
    )
    # first initialize with dummy parameters
    Pars = np.tile(default_parameters, (1, 1))
    X0 = np.tile(initial_state, (1, 1))
    integrator.initialize(X0, Pars)


    t_average = []
    print(f"Time for 1000 RK4 steps of the Lorenz system. Average of {reps} repetitions.\n(N, time)")
    for num in num_pts:
        tsum = 0.0
        Pars = np.tile(default_parameters, (num, 1))
        X0 = np.tile(initial_state, (num, 1))
        integrator.set_problem_data(X0, Pars)

        #warm-up pass
        for _ in range(3):
            integrator.transient(update_x0=False)
            # integrator.get_final_state() #include?

        for _ in range(reps):
            t0 = time.perf_counter()
            # integrator.initialize(X0, Pars) #include?
            integrator.transient(update_x0=False)
            # integrator.get_final_state() #include?
            tsum += time.perf_counter() - t0

        print(f"{num}\t{tsum/reps:0.4} s")
        t_average.append(tsum/reps)

    return t_average
    

# if using 'bazel test ...'
if __name__ == "__main__":
    # clode.set_log_level(clode.LogLevel.info)    
    # clode.runtime.print_opencl()
    # clode.set_log_level(clode.LogLevel.warn) 

    ocl_info = clode.runtime.query_opencl()
    ocl_info = ocl_info[:3]

    num_pts = [int(2**n) for n in np.arange(0,18)]
    reps = 30
    
    for i, ocl in enumerate(ocl_info):
        if ocl.device_count==0:
            continue
        print("\n", ocl)
        t_average = test_lorenz_rk4(platform_id=i, device_id=0, num_pts=num_pts, reps=reps)
        plt.plot(num_pts, t_average)
    
    plt.legend([ocl.name for ocl in ocl_info])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()