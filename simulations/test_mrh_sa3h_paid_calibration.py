"""Independent count/event factoring and artifact replay for SA-3H."""

from fractions import Fraction as F
import hashlib
import itertools
import json
from math import comb
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
import mrh_sa3h_paid_calibration as instrument


def direct_mass(n, k, p):
    return comb(n, k) * p ** k * (1 - p) ** (n - k)


class PaidCalibrationTests(unittest.TestCase):
    def test_exact_count_acceptance_rules(self):
        for n in (0, 8, 32):
            upper, lower = instrument.cutoffs(n)
            for k in range(n + 1):
                cdf = sum((direct_mass(n, j, instrument.U) for j in range(k + 1)), F(0))
                tail = sum((direct_mass(n, j, instrument.LOWER) for j in range(k, n + 1)), F(0))
                self.assertEqual(k <= upper, cdf <= instrument.DELTA)
                self.assertEqual(k >= lower, tail <= instrument.DELTA)

    def test_independent_calibration_count_pairs(self):
        for n, p0, p1 in itertools.product((0, 8, 32), instrument.NULL_RATES, instrument.SIGNAL_RATES):
            upper, lower = instrument.cutoffs(n)
            both = false_claim = F(0)
            for x, y in itertools.product(range(n + 1), repeat=2):
                mass = direct_mass(n, x, p0) * direct_mass(n, y, p1)
                a, b = x <= upper, y >= lower
                both += mass * (a and b)
                false_claim += mass * ((a and p0 > instrument.U) or (b and p1 < instrument.LOWER))
            row = instrument.configuration(n, p0, p1, 24)
            self.assertEqual(row['both_bounds_certified'], both)
            self.assertEqual(row['false_calibration_claim'], false_claim)

    def test_joint_bad_event_independent_indicator_enumeration(self):
        for n, p0, p1, cap in itertools.product((0, 8, 32, 128), instrument.NULL_RATES,
                                                instrument.SIGNAL_RATES, instrument.CAPS):
            row = instrument.configuration(n, p0, p1, cap)
            a, b, c = (row[k] for k in ('deployment_authorized', 'signal_floor_certified', 'conditional_null_crossing'))
            probability = F(0)
            for allow, certify, cross in itertools.product((False, True), repeat=3):
                mass = (a if allow else 1 - a) * (b if certify else 1 - b) * (c if cross else 1 - c)
                bad = (allow and p0 > instrument.U) or (certify and p1 < instrument.LOWER) or (allow and cross)
                probability += mass * bad
            self.assertEqual(row['false_claim_or_null_crossing'], probability)

    def test_artifact_replay_and_hashes(self):
        source = Path(instrument.__file__)
        saved = json.loads(source.with_name('mrh_sa3h_paid_calibration_results.json').read_text())
        self.assertEqual(saved['source_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(saved['registration_commit'], instrument.REGISTRATION)
        for name, checksum in saved['dependency_sha256'].items():
            self.assertEqual(checksum, hashlib.sha256(source.with_name(name).read_bytes()).hexdigest())
        actual = instrument.experiment()
        self.assertEqual(len(saved['results']['configurations']), 216)
        self.assertEqual(saved['results'], actual)
        self.assertEqual(sum(len(rows) for rows in actual['fresh_deployment_curves'].values()), 679)


if __name__ == '__main__':
    unittest.main()
