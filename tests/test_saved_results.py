"""Check the committed reference data without requiring a MATLAB license."""
from pathlib import Path
import unittest

import numpy as np
from scipy.io import loadmat
from CoolProp.CoolProp import PropsSI


ROOT = Path(__file__).resolve().parents[1]


class SavedResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = {
            name: loadmat(ROOT / "results" / f"{name}.mat", simplify_cells=True)
            for name in ("sweep_results", "power_decay_results",
                         "SOC_scenario_results", "constH_scenario_results")
        }

    def test_valid_states_and_failure_semantics(self):
        for name, s in self.data.items():
            with self.subTest(dataset=name):
                ok = np.asarray(s["ok"], dtype=bool)
                status = np.asarray(s["status"])
                self.assertEqual(int(s["failed_count"]), int(np.count_nonzero(~ok)))
                mdot = np.asarray(s.get("mdot_map", s.get("mdot")))
                temp = np.asarray(s.get("Tout_map", s.get("T_out")))
                quality = np.asarray(s.get("x_map", s.get("x_out")))
                self.assertEqual(mdot.shape, ok.shape)
                self.assertTrue(np.all(np.isnan(mdot[~ok])))
                solved = ok & (status == "solved")
                self.assertTrue(np.all(np.isfinite(mdot[solved])))
                self.assertTrue(np.all(mdot[solved] > 0))
                self.assertTrue(np.all(np.isfinite(temp[solved])))
                self.assertTrue(np.all(temp[solved] <= s["p"]["T_out_max"] + 1e-6))
                self.assertTrue(np.all((quality[solved] >= 0) & (quality[solved] <= 1)))
                self.assertTrue(np.all(np.abs(np.asarray(s["pressure_residual"])[solved]) < 1e-3))
                self.assertTrue(np.all(np.abs(np.asarray(s["energy_residual"])[solved]) < 1e-6))
                self.assertIn("coolprop_version", s["metadata"])

    def test_reference_point(self):
        r = self.data["sweep_results"]["sp"]
        self.assertTrue(r["ok"])
        self.assertAlmostEqual(60 * r["mdot"], 1.60928222448, places=7)
        self.assertAlmostEqual(r["T_out"], 61.6333127716, places=6)

    def test_default_sweep_coverage(self):
        s = self.data["sweep_results"]
        self.assertEqual(np.shape(s["ok"]), (32, 50))
        self.assertEqual(int(s["failed_count"]), 261)
        self.assertEqual(set(np.asarray(s["status"]).ravel()), {"solved", "no_root_in_domain"})
        np.testing.assert_allclose(np.asarray(s["H_list"])[[0, -1]], [.01, 2])
        np.testing.assert_allclose(np.asarray(s["Q_list"])[[0, -1]], [360, 36000])

    def test_shared_scenario_drivers_and_paired_means(self):
        soc = self.data["SOC_scenario_results"]
        const = self.data["constH_scenario_results"]
        for field in ("t", "Qdot"):
            np.testing.assert_array_equal(soc[field], const[field])
        for s in (soc, const):
            self.assertTrue(np.all(s["ok"]))
            for k in range(len(s["labels"])):
                paired = (np.asarray(s["ok"], bool)[:, k]
                          & np.isfinite(s["mdot"][:, k])
                          & np.isfinite(s["T_out"][:, k]))
                self.assertEqual(s["paired_count"][k], np.count_nonzero(paired))
                self.assertAlmostEqual(s["mean_mdot"][k], np.mean(s["mdot"][paired, k]))
                self.assertAlmostEqual(s["mean_T_out"][k], np.mean(s["T_out"][paired, k]))

    def test_inventory_conservation(self):
        s = self.data["power_decay_results"]
        self.assertTrue(np.all(s["ok"]))
        np.testing.assert_allclose(s["mass_L"] + s["mass_vaporized"], s["mass_L"][0], atol=1e-9, rtol=0)
        self.assertTrue(np.all(np.diff(s["mass_L"]) <= 0))
        np.testing.assert_allclose(-np.diff(s["mass_L"]), s["interval_vaporized"], atol=1e-9, rtol=0)
        np.testing.assert_allclose(s["vapor_transfer_rate"], s["mdot"] * s["x_out"], atol=1e-12, rtol=0)

    def test_scaled_calibration_target_is_inside_temperature_domain(self):
        p = self.data["sweep_results"]["p"]
        h_in = PropsSI("H", "P", p["P"], "T", p["T_L"], p["fluid"])
        h_out = h_in + p["Qdot_ref"] / (0.30 * 3.6 / 60)
        temp = PropsSI("T", "P", p["P"], "H", h_out, p["fluid"])
        self.assertLess(temp, p["T_out_max"])
        self.assertAlmostEqual(temp, 94.1517558349, places=6)


if __name__ == "__main__":
    unittest.main()
