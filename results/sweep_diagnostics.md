# Sweep domain diagnosis

Independent Python equation audit, not MATLAB execution. Finite sign-change scans are not a proof of no root.

Source MAT SHA-256: `70cf2237a3b89975a7c6843977244b43733ff02408f386b2c91bf6d02a5ac174`.

CoolProp 7.2.0; 512 logarithmic samples per domain, plus phase boundaries. Flow searched to 5 kg/s. Configured cap: 250 K; relaxed cap: EOS limit 1000 K.

| Category | Points |
|---|---:|
| `temperature_cap_excluded` | 176 |
| `no_root_found_to_eos_limit` | 85 |
| `accepted` | 1339 |

## Representative unsolved points

| H [m] | Q [kW] | Diagnosis | Relaxed roots [kg/min] | Outlet T [K] |
|---:|---:|---|---|---|
| 0.010000 | 2.542041 | `no_root_found_to_eos_limit` | None found | Undefined |
| 0.010000 | 33.090612 | `no_root_found_to_eos_limit` | None found | Undefined |
| 0.202581 | 36.000000 | `no_root_found_to_eos_limit` | None found | Undefined |
| 0.010000 | 0.360000 | `temperature_cap_excluded` | 0.001902 | 776.783 |
| 0.330968 | 22.180408 | `temperature_cap_excluded` | 0.347706 | 261.858 |
| 0.844516 | 36.000000 | `temperature_cap_excluded` | 0.588354 | 251.629 |

## Interpretation

- `temperature_cap_excluded`: the denser scan finds no configured-domain root but finds a pressure-balanced root when the cap is relaxed. These out-of-domain states are diagnostic only; they are not accepted hardware operating points.
- `no_root_found_to_eos_limit`: neither scan found a sign-changing root. This does not establish physical stagnation; roots outside the EOS/flow bounds or missed by sampling are not ruled out.
- `root_found_on_dense_scan`: investigate sampling/bounds in the original solver before interpreting the saved failure.
- Source data, configured cap, and production solver remain unchanged.
