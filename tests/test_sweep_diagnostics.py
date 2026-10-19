"""Check diagnostic provenance, domain classification, and relaxed roots."""
from pathlib import Path
import hashlib
import json
import sys
import unittest

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
from diagnose_sweep import EquationAudit


class SweepDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = json.loads((ROOT/'results'/'sweep_diagnostics.json').read_text())
        cls.s = loadmat(ROOT/'results'/'sweep_results.mat', simplify_cells=True)
        cls.audit = EquationAudit(cls.s['p'])

    def test_report_matches_source_and_failed_cells(self):
        source = ROOT/'results'/'sweep_results.mat'
        self.assertEqual(self.report['source_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(self.report['counts'], {'accepted': 1339, 'temperature_cap_excluded': 176,
                                                'no_root_found_to_eos_limit': 85})
        cells = {(r['row'], r['column']) for r in self.report['records']}
        self.assertEqual(cells, set(map(tuple, np.argwhere(~np.asarray(self.s['ok'], bool)))))

    def test_relaxed_roots_obey_balance_and_exceed_original_cap(self):
        for r in self.report['records']:
            with self.subTest(height=r['height_m'], power=r['power_W']):
                self.assertFalse(r['roots_at_configured_cap_kg_s'])
                if r['diagnosis'] == 'temperature_cap_excluded':
                    self.assertTrue(r['roots_at_eos_limit_kg_s'])
                    for m, recorded_temp in zip(r['roots_at_eos_limit_kg_s'], r['relaxed_root_temperatures_K']):
                        residual, temp = self.audit.state(m, r['power_W'])
                        residual += self.s['p']['g']*self.audit.rhoL*r['height_m']
                        self.assertLess(abs(float(residual[0])), 1e-3)
                        self.assertGreater(float(temp[0]), self.report['configured_temperature_cap_K'])
                        self.assertLessEqual(float(temp[0]), self.report['eos_temperature_limit_K'])
                        self.assertAlmostEqual(float(temp[0]), recorded_temp, places=6)
                else:
                    self.assertFalse(r['roots_at_eos_limit_kg_s'])

    def test_representative_points_with_denser_scan(self):
        for name in ('temperature_cap_excluded', 'no_root_found_to_eos_limit'):
            cases = [r for r in self.report['records'] if r['diagnosis'] == name]
            for idx in (0, len(cases)//2, len(cases)-1):
                r = cases[idx]
                with self.subTest(category=name, height=r['height_m'], power=r['power_W']):
                    current = self.audit.roots(r['height_m'], r['power_W'], 250, 1024)
                    relaxed = self.audit.roots(r['height_m'], r['power_W'], 1000, 1024)
                    self.assertFalse(current)
                    np.testing.assert_allclose(relaxed, r['roots_at_eos_limit_kg_s'], atol=1e-9, rtol=0)


if __name__ == '__main__':
    unittest.main()
