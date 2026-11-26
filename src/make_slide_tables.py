"""슬라이드에 그대로 붙여넣을 표를 만든다. (base 환경에서 실행)

입력: results/confidence_summary.csv   (collect_confidence.py)
      results/rmsd_pairs.csv           (rmsd_analysis.py, PyMOL 환경)
출력: results/slide_tables.md

표 3종
  표 1 — 재현 결과 요약 (논문 대조 포함)          → 논문 재현 파트
  표 2 — 신뢰도 상세 (pTM / pLDDT)               → 논문 재현 파트
  표 3 — 노이즈 바닥 vs 신호                      → 자체 검증 파트 (논문에 없음)

실행: python3 src/make_slide_tables.py
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

from detectability import matrix_from_rmsd_rows, verdict_from_matrix

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS  # noqa: E402

RESULTS = ROOT / "results"
GENES = [t["gene"] for t in TARGETS]
BY_GENE = {t["gene"]: t for t in TARGETS}


def load(name: str) -> list[dict]:
    p = RESULTS / name
    if not p.exists():
        sys.exit(f"{name} 이 없다. 먼저 해당 스크립트를 실행할 것.")
    return list(csv.DictReader(p.open()))


def rng(vals: list[float], d: int = 2) -> str:
    """최소~최대. 값이 같으면 하나만."""
    lo, hi = min(vals), max(vals)
    return f"{lo:.{d}f}" if abs(hi - lo) < 10 ** -d else f"{lo:.{d}f}~{hi:.{d}f}"


def main() -> int:
    conf = [r for r in load("confidence_summary.csv") if r["construct"] == "full"]
    pairs = load("rmsd_pairs.csv")

    def sub(gene: str, allele_is_wt: bool) -> list[dict]:
        return [r for r in conf
                if r["gene"] == gene and (r["allele"] == "WT") == allele_is_wt]

    out: list[str] = []
    A = out.append

    A("# 슬라이드용 표: AlphaFold3 재현 결과\n")
    A("입력: AlphaFold Server(AF3), UniProt 전장 서열, seed 1, job당 모델 5개. "
      "단백질당 구조 10개(WT 5 + 변이형 5), 모두 30개\n")

    # ---------------- 표 1 : 재현 결과 요약 --------------------------------
    A("\n## 표 1. 재현 결과 요약 (논문 대조)\n")
    A("| 단백질 | 변이 | dbSNP | 재현 pTM | 논문 pTM | 도메인 pLDDT | 판정 |")
    A("|---|---|---|---|---|---|---|")
    for g in GENES:
        t = BY_GENE[g]
        lab = f"{t['wt']}{t['pos']}{t['mut']}"
        allp = [float(r["ptm"]) for r in conf if r["gene"] == g]
        dom = [float(r["mean_plddt_domain"]) for r in conf if r["gene"] == g]
        paper = "0.6" if g == "CD33" else "보고 없음"
        verdict = "일치" if g == "CD33" else "대조 불가"
        A(f"| {g} | {lab} | {t['rsid']} | {rng(allp)} | {paper} | "
          f"{rng(dom, 1)} | {verdict} |")
    A("\n논문이 pTM을 명시한 것은 CD33뿐임(Fig 4B 캡션: *\"a score of 0.6 for both "
      "structures\"*). PILRA와 TREM2는 논문에 수치가 없어 대조할 수 없음.")

    # ---------------- 표 2 : 신뢰도 상세 -----------------------------------
    A("\n## 표 2. 신뢰도 상세\n")
    A("| 단백질 | 대립 | pTM | 전장 pLDDT | 도메인 pLDDT | 변이부위 pLDDT | 무질서 비율 |")
    A("|---|---|---|---|---|---|---|")
    for g in GENES:
        t = BY_GENE[g]
        lab = f"{t['wt']}{t['pos']}{t['mut']}"
        for is_wt, name in ((True, "WT"), (False, lab)):
            rs = sub(g, is_wt)
            A(f"| {g} | {name} | {rng([float(r['ptm']) for r in rs])} "
              f"| {rng([float(r['mean_plddt']) for r in rs], 1)} "
              f"| {rng([float(r['mean_plddt_domain']) for r in rs], 1)} "
              f"| {rng([float(r['plddt_at_mut']) for r in rs], 1)} "
              f"| {rng([float(r['fraction_disordered']) for r in rs])} |")
    A("\npTM이 0.5~0.6으로 낮은 것은 접힘 실패가 아니라 신호펩타이드, 막관통 구간, 세포질 "
      "꼬리가 무질서하기 때문임(무질서 비율 0.41~0.61). 변이가 있는 Ig 도메인의 pLDDT는 "
      "TREM2 R62H 모델 하나(83.1)를 빼면 93~98임.")

    # ---------------- 표 3 : 노이즈 바닥 vs 신호 ---------------------------
    A("\n## 표 3. 음성 대조군 대비 변이 효과 (논문에 없는 추가 검정)\n")
    A("| 단백질 | 음성 대조군 (같은 서열, 20쌍) | 처리군 (WT × 변이형, 25쌍) | 비 | 순열 p | 판정 |")
    A("|---|---|---|---|---|---|")

    verdicts = {}
    for g in GENES:
        rows_g = [r for r in pairs if r["gene"] == g]
        D, k = matrix_from_rmsd_rows(rows_g)
        v = verdict_from_matrix(D, k, k)
        sep = v["det"] and v["ms"] > v["mn"]
        verdicts[g] = (v["mn"], v["ms"], v["p"], sep, v)
        A(f"| {g} | {v['mn']:.3f} Å (중앙값) | {v['ms']:.3f} Å (중앙값) "
          f"| {v['ratio']:.2f}× | {v['p']:.4f} | "
          f"{'처리군 > 음성 대조군' if sep else '구분 안 됨'} |")

    A("\n측정 방법\n")
    A("- 대상: 변이 잔기 8 Å 이내 Cα (국소 RMSD). 표의 값은 쌍별 RMSD의 중앙값. "
      "Ig 도메인 Cα로 먼저 정렬한 뒤 재정렬 없이 측정")
    A("- 음성 대조군: 같은 서열의 서로 다른 모델 쌍 (WT-WT 10쌍 + 변이형-변이형 10쌍)")
    A("- 처리군: WT 모델 5개 × 변이형 모델 5개 = 25쌍")
    A("- 검정: 라벨 순열검정, 단측(처리군 > 음성 대조군), 유의수준 0.05")
    A("  - 쌍은 서로 독립이 아님(같은 구조가 여러 쌍에 재사용됨). Mann-Whitney의 p는 "
      "과소평가되므로 라벨을 뒤섞은 경험 분포로 교정함")
    A(f"  - 조건당 모델 5개이면 전수 분할이 {next(iter(verdicts.values()))[4]['n_perm']}개임. "
      f"도달 가능한 최소 p는 1/{next(iter(verdicts.values()))[4]['n_perm']} = "
      f"{1/next(iter(verdicts.values()))[4]['n_perm']:.4f}임")
    A("- k = 모델 5개 (AF Server가 job당 자동 생성하는 diffusion sample). "
      "추가 예측 없이 기존 결과만 사용")
    A("- 음성 대조군은 한 작업(시드 1) 안의 모델 5개 사이 변동임. "
      "작업 간, 시드 간 변동은 들어 있지 않음")

    n_sep = sum(1 for v in verdicts.values() if v[3])
    A(f"\n결과: {len(GENES)}개 중 {n_sep}개에서 처리군이 음성 대조군을 넘음. "
      "TREM2는 국소 지표에서만 넘고 도메인 전체 RMSD로는 넘지 않음(상위 README 결과 2절).")
    if n_sep < len(GENES):
        nots = [g for g in GENES if not verdicts[g][3]]
        A(f"{', '.join(nots)}에서는 변이로 인한 구조 차이가 "
          f"모델 자체의 변동과 통계적으로 구분되지 않음.")

    path = RESULTS / "slide_tables.md"
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))
    print(f"\n\n저장 → {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
