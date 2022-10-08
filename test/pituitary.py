import os.path

import warnings

import numpy as np
import pandas as pd
from collections import OrderedDict

from pyclode import Observer

peakutils = None

import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import pyclode as clode

#  OrderedDict remembers the key order. Required for python<3.9
default_ode_parameters = OrderedDict([
    ('g_CaL', (0, 4, 0)),
    ('g_CaT', (0, 4, 0)),
    ('g_K', (0, 10, 0)),
    ('g_SK', (0.05, 3.5, 0)),
    ('g_Kir', (0, 2, 0)),
    ('g_BK', (0, 4, 0)),
    ('g_NaV', (6, 16, 0)),
    ('g_A', (0, 100, 0)),
    ('g_leak', (0.05, 0.4, 0)),
    ('C_m', (4, 12, 10)),
    ('E_leak', (-75, -10, -50)),
    ('tau_m', (0.7, 1.3, 1)),
    ('tau_ht', (.7, 1.3, 1)),
    ('tau_n', (20, 40, 20)),
    ('tau_BK', (2, 10, 1)),
    ('tau_h', (10, 30, 1)),
    ('tau_hNa', (1.4, 2.6, 1)),
    ('k_c', (0.03, 0.21, 0.15))]
)


def pituitary_ori_ode_parameters():
    parameters = dict()

    # Maximal conductance
    parameters['g_CaL'] = np.random.uniform(0, 4)
    parameters['g_K'] = np.random.uniform(0., 10.)
    parameters['g_leak'] = np.random.uniform(0.05, 0.4)

    # Kinetic variables
    parameters['k_c'] = np.random.uniform(0.03, 0.21)
    parameters['tau_n'] = np.random.uniform(20, 40)

    return parameters


def pituitary_ori_ode_parameters_Isk():
    parameters = pituitary_ori_ode_parameters()

    # Maximal conductance
    parameters['g_SK'] = np.random.uniform(.5, 3.5)

    # Other structural parameters
    parameters['C_m'] = np.random.uniform(4, 12)

    return parameters


def pituitary_ori_ode_parameters_Isk_Ibk():
    parameters = pituitary_ori_ode_parameters_Isk()

    # Maximal conductance
    parameters['g_BK'] = np.random.uniform(0, 4)

    # Kinetic variables
    parameters['tau_m'] = np.random.uniform(.7, 1.3)
    parameters['tau_BK'] = np.random.uniform(2, 10)

    return parameters


def pituitary_ori_ode_parameters_Isk_Ibk_Ikir():
    parameters = pituitary_ori_ode_parameters_Isk_Ibk()

    # Maximal conductance
    parameters['g_Kir'] = np.random.uniform(0, 2)

    # Other structural parameters
    parameters['E_leak'] = np.random.uniform(-75, -10)

    return parameters


def pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat():
    parameters = pituitary_ori_ode_parameters_Isk_Ibk_Ikir()

    # Maximal conductance
    parameters['g_CaT'] = np.random.uniform(0, 4)

    # Kinetic variables
    parameters['tau_ht'] = np.random.uniform(.7, 1.3)

    return parameters


def pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat_Ia():
    parameters = pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat()

    # Maximal conductance
    parameters['g_A'] = np.random.uniform(0, 100)

    parameters['tau_h'] = np.random.uniform(10, 30)

    return parameters


def pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat_Ia_Inav():
    parameters = pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat_Ia()

    # Maximal conductance
    parameters['g_NaV'] = np.random.uniform(6, 16)

    # Kinetic variables
    parameters['tau_hNa'] = np.random.uniform(1.4, 2.6)

    return parameters


def generate_clode_pituitary(parameter_function, dt: float, num_simulations: int):

    tspan = (0.0, 1000.0)

    integrator = clode.CLODEFeatures(
        src_file="test/test.cl",
        variable_names=["V", "n", "m", "b", "h", "h_T", "h_Na", "c"],
        parameter_names=list(default_ode_parameters.keys()),
        num_noise=1,
        aux=["i_noise"],
        event_var="V",
        feature_var="V",
        observer=Observer.threshold_2,
        stepper=clode.Stepper.rk4,
        dt=0.1,
        dtmax=1,
        tspan=tspan,

    )
    #print(list(default_ode_parameters.keys()))

    #
    # pi = clode.problem_info("test/test.cl", 8, 18, 0, 1,
    #                         ["V", "n", "m", "b", "h", "h_T", "h_Na", "c"],
    #                         [], [])
    # stepper = "euler"
    # observer = "basic"
    # nReps = 1
    # nPts = 1000000
    # sp = clode.solver_params(0.1, 1.00, 1e-6, 1e-3, 10000000, 10000000, 50)
    # op = clode.observer_params(0, 7, 100, 1, 1, 0.01, 0.3, 0.2, 0, 0, 1e-7)
    # open_cl = clode.opencl_resource()
    # clode_features = clode.clode_features(pi, stepper, observer, True, open_cl, "src/")

    # Initial conditions
    V = -60.
    n = 0.1
    m = 0.1
    b = 0.1
    h = 0.01
    h_T = 0.01
    h_Na = 0.01
    c = 0.1

    x0 = np.array([(V, n, m, b, h, h_T, h_Na, c)])
    # x0_v = []
    # for index in range (len(x0)):
    #     for _ in range(nPts):
    #         x0_v.append(x0[index])
    num_simulations = 32

    x0_v = np.tile(x0, (num_simulations, 1))

    abserr = 1.0e-8
    relerr = 1.0e-6


    pars_v = []
    pars = ()
    for _ in range(num_simulations):
        parameters = parameter_function()

        # We need to add default parameters to the simulation so
        # our ODE function executes correctly
        simulation_parameters = parameters.copy()
        for key in default_ode_parameters.keys():
            if key not in simulation_parameters:
                simulation_parameters[key] = default_ode_parameters[key][2]

        # Note - the key order in default_ode_parameters is correct.
        # Disregard the key order in parameters
        pars = np.array([simulation_parameters[key] for key in default_ode_parameters.keys()])
        pars_v.append(pars)
        #print(pars)
        #break
    pars_v = np.array(pars_v)
    #pars = np.array((1.4, 0, 5, 0, 0, 0, 0, 0, 0.2, 10, -50, 1, 1, 39, 1, 1, 1, 0.03))

    #pars_v = np.tile(pars, (num_simulations, 1))
    # pars_v = []
    # for index in range (len(pars)):
    #     for _ in range(nPts):
    #         pars_v.append(pars[index])

    #pars_v = np.array(pars_v)
    #print(x0_v)
    integrator.initialize(x0_v, pars_v)

    integrator.transient()

    #integrator.features(initialize_observer=True)
    integrator.features()
    simulation_output = integrator.get_final_state()
    observer_output = integrator.get_observer_results()
    # for _ in range(nReps):
    #     clode_features.features()
    #
    # print(simulation_output)
    #
    # wsol = clode_features.get_f()
    #
    # F = clode_features.get_f()
    # print(pars)
    # print(len(F))
    # for i in range(clode_features.get_n_features()):
    #     print(F[i*nPts], end=" ")
    # print()
    # for i in range(8):
    #     print(clode_features.getXf()[i * nPts])

    return simulation_output, pars_v, observer_output


def find_pituitary_activation_event(wsol_trimmed, V_threshold, dV_max_threshold, dV_min_threshold, dVs):
    below_event_end_threshold = np.all(np.stack((wsol_trimmed[1:] < V_threshold, dVs > dV_min_threshold), axis=-1),
                                       axis=1)
    above_event_start_threshold = np.all(np.stack((wsol_trimmed[1:] > V_threshold, dVs > dV_max_threshold), axis=-1),
                                         axis=1)
    # Skip to the end of an event
    first_event_end_index = np.argmax(below_event_end_threshold)

    if first_event_end_index >= 20000:
        return 0, 20000

    # Skip to the start of the next event
    event_start_index = np.argmax(above_event_start_threshold[first_event_end_index:]) + first_event_end_index

    if event_start_index >= 20000:
        return 0, 20000

    # Skip to the end of the next event
    event_end_index = np.argmax(below_event_end_threshold[event_start_index:]) + event_start_index

    if event_end_index >= 20000:
        return 0, 20000

    return event_start_index, event_end_index


def classify_pituitary_ode(wsol, dt, recognise_one_burst_spiking=False):
    """
    Classifies the pituitary ODE as either spiking, bursting,
    depolarised, hyperpolarised or one-spike bursting.

    The classes are returned as numbers:
    0: hyperpolarised
    1: Depolarised
    2: Spiking
    3: Bursting
    4: One-spike bursting

    Arguments:
        wsol : The sequence to classify
        recognise_one_burst_spiking : Whether to recognise one-spike bursting
                                      as separate from bursting.
    """
    min_distance = 50.

    first_10_s = round(5000. / dt)
    wsol_trimmed = wsol[first_10_s:]

    highest_peak = np.max(wsol_trimmed)
    lowest_valley = np.min(wsol_trimmed)

    min_amplitude = 10.  # 20, 30?

    if highest_peak > lowest_valley + min_amplitude:
        dVs = np.diff(wsol_trimmed)

        highest_dV = np.max(dVs)
        lowest_dV = np.min(dVs)

        V_threshold = (highest_peak - lowest_valley) * .35 + lowest_valley
        dV_max_threshold = highest_dV * .25
        dV_min_threshold = lowest_dV * .25

        try:
            top_peaks = peakutils.indexes(wsol_trimmed, thres=0.8)
        except ValueError:
            # A small number of sequences generate ValueErrors.
            # Since there are only a few, we will suppress it.
            top_peaks = []
        if len(top_peaks) > 1:
            stride = max(round((top_peaks[1] - top_peaks[0]) / 50), 1)
        else:
            stride = 1

        event_start, event_end = find_pituitary_activation_event(wsol_trimmed, V_threshold, dV_max_threshold,
                                                                 dV_min_threshold, dVs)
        # Error handling
        if event_end >= 20000:
            return 0, (19000, 19500)
        elif event_start == event_end:
            return 0, (19000, 20000)

        min_dist = round((event_end - event_start) * 0.05)
        if min_dist < min_distance:
            min_dist = min_distance

        try:
            nearby_peaks = peakutils.indexes(wsol_trimmed[event_start:event_end], thres=0.3, min_dist=min_dist)
        except ValueError:
            # A small number of sequences generate ValueErrors.
            # Since there are only a few, we will suppress it.
            nearby_peaks = []

        if len(nearby_peaks) > 1:
            return 3, (event_start, event_end)

        else:
            event_amplitude = np.max(wsol_trimmed[event_start:event_end]) - np.min(wsol_trimmed[event_start:event_end])

            V_area = np.sum(wsol_trimmed[event_start:event_end]) - V_threshold * (event_end - event_start)
            if V_area > 3000 / dt and event_amplitude > 30:
                if recognise_one_burst_spiking:
                    return 4, (event_start, event_end)
                else:
                    return 3, (event_start, event_end)
            else:
                return 2, (event_start, event_end)

    else:
        if highest_peak > -30:
            return 1, (19000, 19500)
        else:
            return 0, (19000, 20000)


def generate_pitutary_dataframe(parameter_function, sample_id: int, trim_start: int, downsample_rate: int,
                                classify: bool, recognise_one_burst_spiking: bool, retain_trajectories: bool,
                                add_timesteps: bool,
                                compute_calcium_concentration: bool):
    #if not classify and not retain_trajectories:
    #    raise ValueError("Error! Generated samples will have no class and no trajectory!")

    dt = 0.5
    num_simulations=100
    pituitary_simulation, parameters, observer_output = generate_clode_pituitary(parameter_function, dt, num_simulations=num_simulations)
    # if retain_trajectories:
    #     df = pd.DataFrame(pituitary_simulation, columns=['V', 'n', 'm', 'b', 'h', 'h_T', 'h_Na', 'c'])
    # else:
    #     df = pd.DataFrame()
    # if add_timesteps:
    #     df['timesteps'] = np.arange(0, 10000, dt)
    # if classify:
    #     voltage_simulation = pituitary_simulation[:, 0]
    #     simulation_class, (_, _) = classify_pituitary_ode(voltage_simulation, dt,
    #                                                       recognise_one_burst_spiking=recognise_one_burst_spiking)
    #     df['class'] = simulation_class
    #     if retain_trajectories:
    #         df['class'] = simulation_class
    #     else:
    #         df['class'] = [simulation_class]
    # if retain_trajectories:
    #     df = df[trim_start:]
    #     df['ID'] = sample_id
    #     for key, value in parameters.items():
    #         df[key] = value
    #
    #     df = df.iloc[::downsample_rate, :]
    # else:
    #     df['ID'] = [sample_id]
    #     for key, value in parameters.items():
    #         df[key] = [value]
    #
    # if compute_calcium_concentration:
    #     df['calcium_concentration'] = sum(pituitary_simulation[16000:, 7] / 4000)

    avg_calcium = observer_output.get_var_mean('c')
    #print(pituitary_simulation.shape, parameters.shape, avg_calcium.shape)


    # The classes are returned as numbers:
    # 0: hyperpolarised
    # 1: Depolarised
    # 2: Spiking
    # 3: Bursting
    # 4: One-spike bursting

    max_v = observer_output.get_var_max('V')
    min_v = observer_output.get_var_min('V')

    max_period = observer_output.get_var_max('period')
    min_period = observer_output.get_var_min('period')

    active = max_v - min_v > 10

    inactive = np.logical_not(active)
    #print(inactive)
    depolarised = (max_v > -30) & inactive

    bursting = (min_period > max_period * 0.8) & active

    classes = np.zeros(active.shape)
    classes[depolarised] = 1
    classes[active] = 2

    #print("Period count", observer_output.get_var_count('period'))
    #print("Step count", observer_output.get_var_count('step'))
    #print(np.concatenate([active, bursting, max_period, min_period], axis=1))

    df_input = np.concatenate([pituitary_simulation, parameters, avg_calcium, classes], axis=1)
    columns = ['V', 'n', 'm', 'b', 'h', 'h_T', 'h_Na', 'c'] + list(default_ode_parameters.keys()) + ['calcium_concentration', 'class']
    df = pd.DataFrame(df_input, columns=columns)
    df.insert(0, 'ID', range(pituitary_simulation.shape[0]))

    return df


def generate_pituitary_dataset(parameter_function, num_samples, trim_start: int = 10000, downsample_rate: int = 20,
                               classify: bool = False, recognise_one_burst_spiking: bool = False,
                               retain_trajectories: bool = False, add_timesteps: bool = False,
                               compute_calcium_concentration: bool = False):
    """
    Computes a dataset of Traja dataframes representing
    pituitary gland simulations. The parameters are taken
    from Fletcher et al (2016) and slightly modified.

    To run this function, provide one of the parameter
    generating functions to create a dictionary of
    parameters. The full list is provided by the
    default_ode_parameters ordered dictionary.

    The parameter functions are of the format:
     * pituitary_ori_ode_parameters
    to
     * pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat_Ia_Inav
    with one ion channel (Ixx) being added in each step.

    Arguments:
        parameter_function :  Function generating a parameter
            dictionary
        num_samples        :  The number of samples to generate
        trim_start         :  How many samples to trim at the start
            of the sequence. Default is 20,000
        downsample_rate    :  The downsampling factor applied to each
            time series. Default is 20, meaning that there are
            100,000/20 = 5,000 steps in each sample.
        classify           :  Whether to classify the sequence as
            spiking, bursting, one-spike bursting (see next option),
            nonexcitable or depolarised. This operation is expensive
            and therefore disabled by default.
        recognise_one_burst_spiking : Whether to recognise one-spike
            bursting as a separate class or consider it a type
            of bursting. Disabling reduces the number of classes
            from 5 to 4. Usually one-spike bursting is less interesting
            when there are more ion channels.
        retain_trajectories : Whether to retain the trajectories in the
            dataframe. The dataframe will be significantly larger.
            Defaults to false.
        add_timesteps : Whether to add timesteps (in the column 'timesteps')
            in the dataframe. Defaults to false.
    """
    if recognise_one_burst_spiking and not classify:
        warnings.warn("Classification not requested but a classification option is set." +
                      "This is likely a mistake - please check the training options")

    df = generate_pitutary_dataframe(parameter_function, sample_id=0, trim_start=trim_start,
                                     downsample_rate=downsample_rate, classify=classify,
                                     recognise_one_burst_spiking=recognise_one_burst_spiking,
                                     retain_trajectories=retain_trajectories, add_timesteps=add_timesteps,
                                     compute_calcium_concentration=compute_calcium_concentration)

    return df


df_test = generate_pituitary_dataset(parameter_function=pituitary_ori_ode_parameters_Isk_Ibk_Ikir_Icat_Ia_Inav,
                                     num_samples=100,
                                     classify=True,
                                     retain_trajectories=False,
                                     compute_calcium_concentration=True)