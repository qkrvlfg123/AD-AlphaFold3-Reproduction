"""감별력 검정 엔진(src/detectability.py) 회귀 테스트.

README 의 결과 표에 실린 값이 저장소의 원자료(results/*.csv)에서 그대로 재현되는지 확인한다.
표준 라이브러리만 쓰므로 추가 설치 없이 실행할 수 있다.

    python tests/test_detectability.py      # 단독 실행
    python -m pytest tests                  # pytest 가 있으면
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from detectability import (matrix_from_rmsd_rows, mulberry32, permutation_p,  # noqa: E402
                           verdict_from_matrix, verdict_from_scalars)


def _read(name: str) -> list[dict]:
    with open(ROOT / "results" / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _structure_verdict(gene: str) -> dict:
    rows = [r for r in _read("rmsd_pairs.csv") if r["gene"] == gene]
    D, k = matrix_from_rmsd_rows(rows)
    return verdict_from_matrix(D, k, k)


def _binding_verdict(gene: str, a_label: str, b_label: str) -> dict:
    rows = [r for r in _read("boltz_iptm.csv") if r["gene"] == gene]
    a = [float(r["iptm"]) for r in rows if r["allele"] == a_label]
    b = [float(r["iptm"]) for r in rows if r["allele"] == b_label]
    return verdict_from_scalars(a, b)


# ---------------------------------------------------------------------------
# 구조 축 (README 결과 2): 조건당 모델 5개 -> 노이즈 20쌍, 신호 25쌍, 전수 열거 126 분할
# ---------------------------------------------------------------------------
STRUCTURE_EXPECTED = {
    #  gene     노이즈 중앙값  신호 중앙값   p        감별
    "CD33": (0.124, 0.120, 0.65, False),
    "TREM2": (0.168, 0.257, 0.0079, True),
    "PILRA": (0.100, 0.995, 0.0079, True),
}


def test_structure_axis_matches_readme():
    for gene, (mn, ms, p, det) in STRUCTURE_EXPECTED.items():
        v = _structure_verdict(gene)
        assert len(v["noise"]) == 20 and len(v["signal"]) == 25, gene
        assert v["method"] == "exact" and v["n_perm"] == 126, gene
        assert round(v["mn"], 3) == mn, (gene, v["mn"])
        assert round(v["ms"], 3) == ms, (gene, v["ms"])
        assert abs(v["p"] - p) < 0.005, (gene, v["p"])
        assert v["det"] is det, gene


def test_minimum_reachable_p_is_one_over_126():
    # 5 + 5 설계에서 전수 분할은 126 개이므로 p 의 바닥값은 1/126 이다.
    for gene in ("TREM2", "PILRA"):
        assert _structure_verdict(gene)["p"] == 1 / 126


# ---------------------------------------------------------------------------
# 결합 축 (README 결과 3): 분할 수가 많아 고정 시드 표본 20,000회
# ---------------------------------------------------------------------------
BINDING_EXPECTED = {
    #  gene     대립형질 A, B     노이즈쌍 신호쌍  노이즈  신호    비     p
    "TREM2": ("WT", "R62H", 210, 225, 0.046, 0.043, 0.95, 0.80),
    "PILRA": ("G78", "R78", 90, 100, 0.058, 0.046, 0.80, 0.97),
}


def test_binding_axis_matches_readme():
    for gene, (a, b, n_noise, n_signal, mn, ms, ratio, p) in BINDING_EXPECTED.items():
        v = _binding_verdict(gene, a, b)
        assert (len(v["noise"]), len(v["signal"])) == (n_noise, n_signal), gene
        assert v["method"] == "sample" and v["n_perm"] == 20000, gene
        assert round(v["mn"], 3) == mn, (gene, v["mn"])
        assert round(v["ms"], 3) == ms, (gene, v["ms"])
        assert round(v["ratio"], 2) == ratio, (gene, v["ratio"])
        assert abs(v["p"] - p) < 0.005, (gene, v["p"])
        assert v["det"] is False, gene


# ---------------------------------------------------------------------------
# 엔진 자체의 성질
# ---------------------------------------------------------------------------
def test_deterministic():
    # 고정 시드 PRNG 이므로 다시 실행해도 p 가 완전히 같아야 한다.
    assert _binding_verdict("PILRA", "G78", "R78")["p"] == _binding_verdict("PILRA", "G78", "R78")["p"]


def test_mulberry32_range_and_repeatability():
    r1, r2 = mulberry32(1), mulberry32(1)
    xs = [r1() for _ in range(1000)]
    assert xs == [r2() for _ in range(1000)]
    assert all(0.0 <= x < 1.0 for x in xs)
    assert len(set(xs)) > 990


def test_detects_clear_shift():
    # 두 그룹이 뚜렷이 떨어져 있으면 감별되어야 한다 (양성 대조).
    v = verdict_from_scalars([0.10, 0.11, 0.12, 0.13, 0.14], [0.80, 0.81, 0.82, 0.83, 0.84])
    assert v["det"] is True and v["p"] == 1 / 126


def test_no_signal_when_groups_identical():
    # 같은 분포에서 나온 두 그룹은 감별되면 안 된다 (음성 대조).
    v = verdict_from_scalars([0.50, 0.52, 0.54, 0.56, 0.58], [0.51, 0.53, 0.55, 0.57, 0.59])
    assert v["det"] is False


def test_label_swap_invariance():
    # 통계량은 A/B 라벨을 맞바꿔도 같아야 한다.
    a, b = [0.2, 0.25, 0.4, 0.41, 0.6], [0.3, 0.35, 0.5, 0.7, 0.9]
    assert verdict_from_scalars(a, b)["p"] == verdict_from_scalars(b, a)["p"]


def test_exact_and_sampled_agree():
    # 같은 자료를 전수 열거와 표본 순열로 계산하면 p 가 가까워야 한다.
    rows = [r for r in _read("rmsd_pairs.csv") if r["gene"] == "CD33"]
    D, k = matrix_from_rmsd_rows(rows)
    exact = permutation_p(D, k, k)
    sampled = permutation_p(D, k, k, exact_limit=0)
    assert exact["method"] == "exact" and sampled["method"] == "sample"
    assert abs(exact["p"] - sampled["p"]) < 0.02, (exact["p"], sampled["p"])


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {name}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
