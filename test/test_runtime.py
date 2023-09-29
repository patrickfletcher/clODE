import pytest

import clode


@pytest.mark.parametrize(
    "device_type, vendor, platform_id, device_id, device_ids",
    [
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            clode.cl_vendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            clode.cl_vendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            clode.cl_vendor.VENDOR_ANY,
            None,
            None,
            [0],
        ],
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            None,
            0,
            None,
            None,
        ],
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            None,
            None,
            0,
            None,
        ],
        [
            clode.cl_device_type.DEVICE_TYPE_CPU,
            None,
            None,
            None,
            [0],
        ],
        [
            None,
            clode.cl_vendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            None,
            clode.cl_vendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            None,
            clode.cl_vendor.VENDOR_ANY,
            None,
            None,
            [0],
        ],
        [
            None,
            None,
            0,
            0,
            [0],
        ],
    ],
)
def test_init_features_runtime_with_incorrect_config_fails(
    device_type, vendor, platform_id, device_id, device_ids
):
    input_file: str = "test/van_der_pol_oscillator.cl"

    tspan = (0.0, 1000.0)

    with pytest.raises(ValueError):
        _ = clode.CLODETrajectory(
            src_file=input_file,
            variable_names=["x", "y"],
            parameter_names=["mu"],
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            tspan=tspan,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )

    with pytest.raises(ValueError):
        _ = clode.CLODEFeatures(
            src_file=input_file,
            variable_names=["x", "y"],
            parameter_names=["mu"],
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            tspan=tspan,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )
