import pytest

import clode


@pytest.mark.parametrize(
    "device_type, vendor, platform_id, device_id, device_ids",
    [
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            clode.CLVendor.VENDOR_ANY,
            None,
            None,
            [0],
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            0,
            None,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            None,
            0,
            None,
        ],
        [
            clode.CLDeviceType.DEVICE_TYPE_CPU,
            None,
            None,
            None,
            [0],
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
            0,
            None,
            None,
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
            None,
            0,
            None,
        ],
        [
            None,
            clode.CLVendor.VENDOR_ANY,
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

    t_span = (0.0, 1000.0)

    with pytest.raises(ValueError):
        _ = clode.TrajectorySimulator(
            src_file=input_file,
            variables={"x": 1.0, "y": 1.0},
            parameters={"mu": 1.0},
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            t_span=t_span,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )

    with pytest.raises(ValueError):
        _ = clode.FeatureSimulator(
            src_file=input_file,
            variables={"x": 1.0, "y": 1.0},
            parameters={"mu": 1.0},
            num_noise=0,
            stepper=clode.Stepper.dormand_prince,
            t_span=t_span,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )
