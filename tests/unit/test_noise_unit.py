import numpy as np
from su2torus.lattice import HoneycombTorus
from su2torus.trotter import trotter_circuit
from su2torus import noise


def test_readout_rate_public_api():
    assert noise._readout_rate(noise.noise_model(1e-3, 3e-5, 2e-3)) == 2e-3
    assert noise._readout_rate(noise.noise_model(1e-3, 3e-5, 0.0)) == 0.0


def test_readout_corrected_echo():
    lat = HoneycombTorus(2, 2)
    circ = trotter_circuit(lat, 1.25, 0.7, 1, 1, False, None)
    init = [0, 4]
    only_readout = noise.noise_model(0.0, 0.0, 0.01)
    # readout alone: the echo returns exactly, so the corrected fidelity is exactly 1
    assert abs(noise.echo_fidelity_readout_corrected(circ, init, only_readout) - 1) < 1e-12
    assert noise.echo_fidelity(circ, init, only_readout) < 0.97
    nm = noise.noise_model(1.5e-3, 3e-5, 1e-3)
    fc = noise.echo_fidelity_readout_corrected(circ, init, nm)
    assert noise.echo_fidelity(circ, init, nm) < fc <= 1
    assert abs(fc / noise.predicted_fidelity(circ, 1.5e-3, 3e-5) - 1) < 0.25
