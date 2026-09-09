#!/usr/bin/env python3
"""SA-3H: exact paid calibration, abstention, and independent deployment."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
# Fixed n<=2048 produces exact joint probabilities with >4300-digit integers.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(20000)
import mrh_sa3g_calibration as previous

REGISTRATION = '1242022f'
PREVIOUS_HASH = '49b273ad2edbf1e0582ae02a9db12c097d0a686a023f6a905a1a3f1ab92cdf6a'
SIZES = (0, 8, 32, 128, 512, 2048)
NULL_RATES = (F(1, 2), F(11, 20), F(3, 5))
SIGNAL_RATES = (F(7, 10), F(3, 4), F(41, 50), F(9, 10))
CAPS = (24, 64, 96)
U, LOWER, Q = F(11, 20), F(3, 4), F(41, 50)
DELTA, ALPHA, TARGET = F(1, 200), F(1, 25), F(4, 5)
THRESHOLD, PRICE = F(25), F(3, 25)
UP, DOWN = Q / U, (1 - Q) / (1 - U)


@lru_cache(maxsize=32)
def binomial_weights(n, p):
    """Integer PMF numerators over one common denominator; no float tails."""
    if type(n) is not int or not 0 <= n <= 2048 or not 0 <= p <= 1:
        raise ValueError('invalid binomial sample size or rate')
    a, b = p.numerator, p.denominator
    denominator = b ** n
    if p == 1:
        return (0,) * n + (1,), 1
    value = (b - a) ** n
    weights = [value]
    for k in range(n):
        value = value * (n - k) * a // ((k + 1) * (b - a))
        weights.append(value)
    return tuple(weights), denominator


@lru_cache(maxsize=None)
def cutoffs(n):
    """X<=upper_cut certifies p0<=U; Y>=lower_cut certifies p1>=LOWER."""
    upper_cut, lower_cut = -1, n + 1
    weights, denominator = binomial_weights(n, U)
    cumulative = 0
    for k, mass in enumerate(weights):
        cumulative += mass
        if cumulative * DELTA.denominator <= DELTA.numerator * denominator:
            upper_cut = k
    weights, denominator = binomial_weights(n, LOWER)
    survival = denominator
    for k, mass in enumerate(weights):
        if survival * DELTA.denominator <= DELTA.numerator * denominator:
            lower_cut = k
            break
        survival -= mass
    return upper_cut, lower_cut


def certification_probabilities(n, p):
    upper_cut, lower_cut = cutoffs(n)
    weights, denominator = binomial_weights(n, p)
    return F(sum(weights[:upper_cut + 1]), denominator), F(sum(weights[lower_cut:]), denominator)


@lru_cache(maxsize=32)
def deployment_curve(p, cap=96):
    live = {0: F(1)}
    crossed = stopped_moment = expected_count = F(0)
    rows = []
    for n in range(cap + 1):
        survival = sum(live.values(), F(0))
        moment = stopped_moment + sum((mass * UP ** s * DOWN ** (n - s)
                                      for s, mass in live.items()), F(0))
        rows.append({'cap': n, 'crossing_probability': crossed,
                     'survival_probability': survival, 'expected_audits': expected_count,
                     'expected_spend': PRICE * expected_count, 'expected_stopped_evidence': moment})
        if n == cap:
            break
        expected_count += survival
        next_live = {}
        for s, mass in live.items():
            for bit, chance in ((1, p), (0, 1 - p)):
                amount, new_s = mass * chance, s + bit
                value = UP ** new_s * DOWN ** (n + 1 - new_s)
                if value >= THRESHOLD:
                    crossed += amount
                    stopped_moment += amount * value
                else:
                    next_live[new_s] = next_live.get(new_s, F(0)) + amount
        live = next_live
    return rows


def configuration(n, p0, p1, cap):
    authorize = certification_probabilities(n, p0)[0]
    certify_signal = certification_probabilities(n, p1)[1]
    null, signal = deployment_curve(p0)[cap], deployment_curve(p1)[cap]
    floor = deployment_curve(LOWER)[cap]['crossing_probability']
    both = authorize * certify_signal
    bad_null, bad_signal = p0 > U, p1 < LOWER
    false_claim = 1 - (1 - authorize * bad_null) * (1 - certify_signal * bad_signal)
    null_cross = authorize * null['crossing_probability']
    # Add null crossings not already included in either false-claim event.
    bad_union = false_claim + null_cross * (not bad_null) * (1 - certify_signal * bad_signal)
    fixed_cost = 2 * n * PRICE
    return {'n_per_calibration_source': n, 'deployment_cap': cap,
            'true_background_rate': str(p0), 'true_signal_rate': str(p1),
            'deployment_authorized': authorize, 'abstention': 1 - authorize,
            'signal_floor_certified': certify_signal, 'both_bounds_certified': both,
            'conditional_power_floor': floor,
            'reports_80_percent_power': both * (floor >= TARGET),
            'false_80_percent_report': both * (floor >= TARGET) * (signal['crossing_probability'] < TARGET),
            'false_calibration_claim': false_claim,
            'false_claim_or_null_crossing': bad_union,
            'end_to_end_null_crossing': null_cross,
            'end_to_end_signal_detection': authorize * signal['crossing_probability'],
            'conditional_null_crossing': null['crossing_probability'],
            'conditional_signal_detection': signal['crossing_probability'],
            'calibration_spend': fixed_cost,
            'expected_null_deployment_spend': authorize * null['expected_spend'],
            'expected_signal_deployment_spend': authorize * signal['expected_spend'],
            'expected_null_total_spend': fixed_cost + authorize * null['expected_spend'],
            'expected_signal_total_spend': fixed_cost + authorize * signal['expected_spend'],
            'maximum_total_audits': 2 * n + cap,
            'counterfactual_supplied_null_bound_valid': not bad_null,
            'counterfactual_supplied_signal_floor_valid': not bad_signal}


def controls():
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    check(hashlib.sha256(Path(previous.__file__).read_bytes()).hexdigest() == PREVIOUS_HASH,
          'inherited source frozen')
    check(U * UP + (1 - U) * DOWN == 1, 'boundary null martingale')
    check(2 * DELTA + ALPHA == F(1, 20), 'declared joint error budget')
    for n in SIZES:
        upper_cut, lower_cut = cutoffs(n)
        for bound, cutoff, upper in ((U, upper_cut, True), (LOWER, lower_cut, False)):
            weights, denominator = binomial_weights(n, bound)
            check(sum(weights) == denominator, 'boundary PMF mass')
            accepted = sum(weights[:cutoff + 1]) if upper else sum(weights[cutoff:])
            check(F(accepted, denominator) <= DELTA, 'exact boundary tail budget')
            if upper:
                check(F(sum(weights[:cutoff + 2]), denominator) > DELTA,
                      'upper cutoff maximal')
            else:
                check(F(sum(weights[max(0, cutoff - 1):]), denominator) > DELTA,
                      'lower cutoff minimal')
        old_upper, old_lower = F(1), F(0)
        for p in (F(k, 100) for k in range(101)):
            weights, denominator = binomial_weights(n, p)
            check(sum(weights) == denominator, 'PMF normalization across rate grid')
            a, b = certification_probabilities(n, p)
            check(a <= old_upper and b >= old_lower, 'monotone certificate probabilities')
            if p >= U:
                check(a <= DELTA, 'upper bound false-certificate control')
            if p <= LOWER:
                check(b <= DELTA, 'lower bound false-certificate control')
            old_upper, old_lower = a, b
        if n == 0:
            check(upper_cut == -1 and lower_cut == 1, 'empty data cannot certify')
    for n in range(9):
        for p in (F(0), F(1, 2), U, LOWER, Q, F(1)):
            empirical = [F(0)] * (n + 1)
            for path in itertools.product((0, 1), repeat=n):
                s = sum(path)
                empirical[s] += p ** s * (1 - p) ** (n - s)
            weights, denominator = binomial_weights(n, p)
            check(empirical == [F(v, denominator) for v in weights], 'independent binomial paths')
            for k in range(n + 1):
                check(sum(empirical[:k], F(0)) + sum(empirical[k:], F(0)) == 1,
                      'complementary tails')
    rates = sorted(set((*NULL_RATES, *SIGNAL_RATES)))
    for p in rates:
        last_probability = last_count = F(0)
        for row, old in zip(deployment_curve(p), previous.curve(F(4, 5), Q, U, p, 96)):
            probability = row['crossing_probability']
            check(probability + row['survival_probability'] == 1, 'deployment mass')
            check(last_probability <= probability <= 1, 'monotone crossing')
            check(last_count <= row['expected_audits'] <= row['cap'], 'audit cap')
            check(row['expected_spend'] == PRICE * row['expected_audits'], 'audit price')
            check(probability == old['crossing_probability'] and row['expected_audits'] == old['expected_audits'],
                  'threshold 25 equals SA-3G scaled start at threshold 20')
            check(row['expected_stopped_evidence'] * F(4, 5) == old['expected_stopped_evidence'],
                  'scaled stopped evidence regression')
            if p <= U:
                check(probability <= ALPHA and row['expected_stopped_evidence'] <= 1,
                      'null bound and stopped moment')
            if p == U:
                check(row['expected_stopped_evidence'] == 1, 'boundary optional stopping')
            if row['cap'] <= 8:
                values = previous.brute_force(F(4, 5), Q, U, p, row['cap'])
                check(values == (probability, row['expected_audits'], row['expected_stopped_evidence'] * F(4, 5)),
                      'independent deployment path enumeration')
            last_probability, last_count = probability, row['expected_audits']
    for low, high in zip(rates, rates[1:]):
        check(all(a['crossing_probability'] <= b['crossing_probability']
                  for a, b in zip(deployment_curve(low), deployment_curve(high))), 'monotone power floor')
    for n, p0, p1, cap in itertools.product(SIZES, NULL_RATES, SIGNAL_RATES, CAPS):
        row = configuration(n, p0, p1, cap)
        check(row['deployment_authorized'] + row['abstention'] == 1, 'authorization partition')
        check(row['false_calibration_claim'] <= 2 * DELTA, 'joint calibration guarantee')
        check(row['false_claim_or_null_crossing'] <= F(1, 20), 'joint bad-event guarantee')
        check(row['false_80_percent_report'] <= DELTA, 'false power report guarantee')
        check(row['reports_80_percent_power'] <= row['both_bounds_certified'], 'no unsupported power report')
        check(row['end_to_end_signal_detection'] <= row['conditional_signal_detection'], 'abstention paid in detection')
        for kind in ('null', 'signal'):
            check(row[f'expected_{kind}_total_spend'] == row['calibration_spend'] + row[f'expected_{kind}_deployment_spend'],
                  'cost accounting')
            check(2 * n * PRICE <= row[f'expected_{kind}_total_spend'] <= (2 * n + cap) * PRICE,
                  'all samples paid within allowance')
    return {'checks_passed': checks}


def number(value):
    return {'exact': str(value), 'decimal': float(value)}


def serialize(row):
    return {key: number(value) if isinstance(value, F) else value for key, value in row.items()}


def experiment():
    rows = [serialize(configuration(n, p0, p1, cap))
            for n, p0, p1, cap in itertools.product(SIZES, NULL_RATES, SIGNAL_RATES, CAPS)]
    calibration = []
    for n in SIZES:
        upper_cut, lower_cut = cutoffs(n)
        calibration.append({'n_per_source': n, 'null_count_at_most': upper_cut,
                            'signal_count_at_least': lower_cut, 'calibration_spend': number(2 * n * PRICE),
                            'null_certificate_probability': {str(p): number(certification_probabilities(n, p)[0]) for p in NULL_RATES},
                            'signal_certificate_probability': {str(p): number(certification_probabilities(n, p)[1]) for p in SIGNAL_RATES}})
    return {'calibration': calibration,
            'fresh_deployment_curves': {str(p): [serialize(row) for row in deployment_curve(p)]
                                        for p in sorted(set((*NULL_RATES, *SIGNAL_RATES)))},
            'configurations': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls-only', action='store_true')
    parser.add_argument('--output', type=Path, help='Create new artifact, refusing overwrite')
    args = parser.parse_args()
    result = {'protocol': 'SA-3H', 'registration_commit': REGISTRATION,
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'dependency_sha256': {Path(previous.__file__).name: PREVIOUS_HASH}, 'controls': controls()}
    if not args.controls_only:
        result['configuration'] = {'sizes_per_source': SIZES, 'deployment_caps': CAPS,
                                   'null_ceiling': str(U), 'signal_floor': str(LOWER),
                                   'assumed_alternative': str(Q), 'calibration_error_per_bound': str(DELTA),
                                   'test_error': str(ALPHA), 'threshold': str(THRESHOLD),
                                   'price_per_sample': str(PRICE), 'power_target': str(TARGET),
                                   'source_assumptions': 'labeled IID stationary sources; independent fresh deployment; no evidence reuse'}
        result['results'] = experiment()
    rendered = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n'
    if args.output:
        with args.output.open('x') as output:
            output.write(rendered)
        print(json.dumps({'output': str(args.output), 'controls': result['controls'],
                          'source_sha256': result['source_sha256']}, sort_keys=True))
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
