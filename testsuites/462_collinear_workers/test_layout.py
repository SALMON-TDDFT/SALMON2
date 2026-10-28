"""Independent literal-loop oracle for the compiled integer layout selector."""
import argparse
import itertools
import random
import subprocess


def assignments(points, groups, block):
    counts = [0] * groups
    seen = []
    for base in range(0, points, groups * block):
        for group in range(groups):
            owned = list(range(base + group * block,
                               min(base + (group + 1) * block, points)))
            counts[group] += len(owned)
            seen.extend(owned)
    assert sorted(seen) == list(range(points)), "duplicate or missing grid index"
    return counts


def oracle(case):
    nfirst, ng, nk, nc, np, manual, divisor = case
    if min(nfirst, ng, nk, np) < 1 or nc not in (1, 2) or min(manual, divisor) < 0:
        return None
    pair_bytes = nk * nc * nc * 16
    if pair_bytes > 67108864 or manual and np % manual:
        return None
    choices = []
    for cols in range(1, np + 1):
        if np % cols or manual and cols != manual:
            continue
        rows = np // cols
        b = min(16, (nfirst + rows - 1) // rows, 67108864 // pair_bytes)
        divisions = divisor or cols
        chunk = min((ng + divisions - 1) // divisions, 67108864 // (b * pair_bytes))
        r = assignments(nfirst, rows, b)
        c = assignments(ng, cols, chunk)
        assert sum(r) * sum(c) == nfirst * ng
        peak = max(r) * max(c)
        active = sum(x > 0 for x in r) * sum(x > 0 for x in c)
        choices.append(((peak, -active, cols), (rows, cols, b, chunk, peak, active, 0)))
    return min(choices)[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    args = parser.parse_args()
    cases = [(164, 4096, 512, 1, 1280, 0, 0),
             (164, 4096, 4096, 1, 27648, 0, 0),
             (305, 8000, 512, 1, 2560, 0, 0),
             (164, 4096, 4096, 1, 27648, 48, 48),
             (8, 8, 3, 1, 4, 1, 16),
             (20, 25, 1048576, 2, 12, 0, 0),
             (2147483647, 2147483647, 1, 1, 1, 1, 2147483647)]
    # Large integer case has a closed-form independent expected result.
    large = cases.pop()
    rng = random.Random(20261028)
    for _ in range(160):
        np = rng.choice([1, 2, 3, 4, 7, 8, 12, 16, 24, 31, 64])
        groups = rng.choice([0, 1, np, 2, 3])
        cases.append((rng.randint(1, 49), rng.randint(1, 63),
                      rng.choice([1, 8, 4096, 262144, 1048576]),
                      rng.choice([1, 2]), np, groups, rng.choice([0, 1, 3, 16, 128])))
    cases += [(0, 8, 8, 1, 4, 0, 0), (8, 8, 8, 1, 4, -1, 0),
              (8, 8, 8, 1, 4, 0, -1), (8, 8, 8, 1, 4, 3, 0),
              (8, 8, 8, 3, 4, 0, 0), (8, 8, 4194305, 1, 4, 0, 0)]
    expected = [oracle(case) for case in cases]
    cases.append(large)
    # One rank processes every pair, with row batch 16 and one column per tile.
    expected.append((1, 1, 16, 1, 2147483647 ** 2, 1, 0))
    proc = subprocess.run([args.binary], input="".join(" ".join(map(str, c)) + "\n" for c in cases),
                          text=True, stdout=subprocess.PIPE, check=True)
    actual = [tuple(map(int, line.split())) for line in proc.stdout.splitlines()]
    assert len(actual) == len(cases)
    for case, want, got in zip(cases, expected, actual):
        if want is None:
            assert got[-1] != 0, (case, got)
        else:
            assert got == want, (case, want, got)
    assert actual[0][:4] == (16, 80, 11, 52)
    assert actual[1][:4] == (27, 1024, 7, 4)
    print("WORKER_LAYOUT_ORACLE_PASS", len(cases), "coverage, caps, overrides, overflow, rejection")


if __name__ == "__main__":
    main()
