#!/usr/bin/env python3
"""Replay documented synthetic challenges; this is not assay validation."""
import argparse
import json
import tempfile
from pathlib import Path

from screen import VERSION, digest, run

ROOT = Path(__file__).resolve().parent
HITS = ROOT / 'data/demo_run/hits.tsv'
META = ROOT / 'data/demo_run/metadata.json'
RULES = ROOT / 'config/demo_rules.json'
REFERENCES = ROOT / 'config/synthetic_references.json'


def replay():
    cases = []

    def check(name, hits=HITS, metadata=META, expected='REVIEW_READY', count=None, qc=None):
        result = run(hits, metadata, RULES, REFERENCES)
        reasons = []
        if result['status'] != expected:
            reasons.append(f'status {result["status"]} instead of {expected}')
        if count is not None and len(result['worklist']) != count:
            reasons.append(f'worklist size {len(result["worklist"])} instead of {count}')
        if qc and result['qc'][qc] is not False:
            reasons.append(f'{qc} did not fail')
        if expected == 'QC_BLOCKED' and any(row['state'] != 'ON_HOLD_QC' for row in result['worklist']):
            reasons.append('candidate released despite QC block')
        cases.append({'case': name, 'passed': not reasons, 'observed_status': result['status'],
                      'worklist_count': len(result['worklist']), 'issues': reasons})
        return result

    nominal = check('nominal_ambiguous', count=2)
    if nominal['worklist'][1]['review_flags'] != ['NEGATIVE_BACKGROUND', 'LIMITED_REGION_SUPPORT', 'HOST_SIMILARITY']:
        cases[-1]['passed'] = False
        cases[-1]['issues'].append('ambiguous signal flags differ')
    check('negative_control_failure', ROOT / 'data/scenarios/negative_failure/hits.tsv',
          expected='QC_BLOCKED', qc='negative_control_clear')
    check('positive_control_failure', ROOT / 'data/scenarios/positive_failure/hits.tsv',
          expected='QC_BLOCKED', qc='expected_positive_detected')
    check('insufficient_depth', metadata=ROOT / 'data/scenarios/low_depth/metadata.json',
          expected='QC_BLOCKED', qc='minimum_depth_met')

    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp) / 'hits.tsv'
        baseline = HITS.read_text(encoding='utf-8')
        for name, coverage, count in [('threshold_equal', '500', 2), ('threshold_below', '499', 1)]:
            scratch.write_text(baseline.replace('24\t5\t1800\t10000\t98.2',
                                                f'3\t2\t{coverage}\t10000\t90'), encoding='utf-8')
            check(name, hits=scratch, count=count)
        scratch.write_text(baseline.replace('SYN-A1', 'SYN-X1'), encoding='utf-8')
        try:
            run(scratch, META, RULES, REFERENCES)
        except ValueError as exc:
            ok = 'unknown reference' in str(exc)
        else:
            ok = False
        cases.append({'case': 'unknown_reference_rejected', 'passed': ok,
                      'issues': [] if ok else ['unknown reference was not rejected']})

    repeat = run(HITS, META, RULES, REFERENCES)
    stable = nominal == repeat
    cases.append({'case': 'deterministic_replay', 'passed': stable,
                  'issues': [] if stable else ['identical inputs produced different results']})
    return {'title': 'Synthetic challenge replay (not assay validation)', 'software': VERSION,
            'all_passed': all(case['passed'] for case in cases), 'cases': cases,
            'input_sha256': {name: digest(path) for name, path in
                             [('hits', HITS), ('metadata', META), ('rules', RULES), ('references', REFERENCES)]},
            'limitations': 'Synthetic summarized evidence only; does not establish analytical performance, GMP compliance or lot release.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/validation_report.json')
    args = parser.parse_args()
    report = replay()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    for case in report['cases']:
        print(f'{"PASS" if case["passed"] else "FAIL"} {case["case"]}' +
              (f': {"; ".join(case["issues"])}' if case['issues'] else ''))
    print(f'Challenge report: {args.output}')
    return 0 if report['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
