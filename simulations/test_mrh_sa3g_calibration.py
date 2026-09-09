"""Additional boundary and committed-artifact checks; no policy tuning."""

from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
import mrh_sa3g_calibration as experiment


class CalibrationChecks(unittest.TestCase):
    def test_nontrivial_near_threshold_exhaustive_paths(self):
        interior_count = 0
        for start, q, u in itertools.product((F(10), F(199, 10), F(20)),
                                            experiment.ASSUMED, experiment.RULES.values()):
            for p in (F(0), u, F(3, 4), F(1)):
                for r in experiment.curve(start, q, u, p, 8):
                    actual = experiment.brute_force(start, q, u, p, r['cap'])
                    self.assertEqual(actual, (r['crossing_probability'], r['expected_audits'],
                                              r['expected_stopped_evidence']))
                    interior_count += 0 < r['crossing_probability'] < 1
        self.assertGreater(interior_count, 0)

    def test_nontrivial_history_dependent_null_adversary(self):
        for start, q, u in itertools.product((F(10), F(199, 10)),
                                            experiment.ASSUMED, experiment.RULES.values()):
            for n in range(9):
                self.assertEqual(experiment.adversarial_null_crossing(start, q, u, n),
                                 experiment.curve(start, q, u, u, 8)[n]['crossing_probability'])

    def test_invalid_factor_models(self):
        for q, u in ((F(1, 2), F(1, 2)), (F(1, 2), F(3, 4)),
                     (F(1), F(1, 2)), (F(3, 4), F(0))):
            with self.assertRaises(ValueError):
                experiment.factors(q, u)

    def test_committed_artifact_schema_hashes_and_exact_replay(self):
        source = Path(experiment.__file__)
        result = json.loads(source.with_name('mrh_sa3g_calibration_results.json').read_text())
        self.assertEqual(result['source_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(result['registration_commit'], experiment.REGISTRATION)
        for filename, checksum in result['dependency_sha256'].items():
            self.assertEqual(checksum, hashlib.sha256(source.with_name(filename).read_bytes()).hexdigest())
        rows = result['results']['curves']
        self.assertEqual(len(rows), 96)
        keys = set()
        for row in rows:
            start, q, u, p = (F(row[k]) for k in ('start', 'assumed_alternative',
                                                  'null_ceiling', 'true_iid_agreement'))
            keys.add((start, q, u, p))
            self.assertEqual(row['null_bound_applicable'], p <= u)
            self.assertEqual(row['null_ceiling'], str(experiment.RULES[row['rule']]))
            self.assertEqual(len(row['caps']), 65)
            for saved, actual in zip(row['caps'], experiment.curve(start, q, u, p)):
                self.assertEqual(saved, experiment.serialize_record(actual))
        self.assertEqual(len(keys), 96)
        reports = result['results']['power_reports']
        self.assertEqual(len(reports), 12)
        for report in reports:
            start, q, u = (F(report[k]) for k in ('start', 'assumed_alternative', 'null_ceiling'))
            self.assertEqual(report['iid_alternative_interval'], ['3/4', '9/10'])
            self.assertEqual(len(report['caps']), 65)
            for n, row in enumerate(report['caps']):
                nominal = experiment.curve(start, q, u, q)[n]['crossing_probability']
                floor = experiment.curve(start, q, u, F(3, 4))[n]['crossing_probability']
                self.assertEqual(row['cap'], n)
                self.assertEqual(row['assumed_model_power'], experiment.number(nominal))
                self.assertEqual(row['interval_floor_power'], experiment.number(floor))
                self.assertEqual(row['assumed_model_supports_80_percent'], nominal >= F(4, 5))
                self.assertEqual(row['entire_interval_supports_80_percent'], floor >= F(4, 5))


if __name__ == '__main__':
    unittest.main()
