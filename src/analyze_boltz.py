"""Boltz-2 결과 집계 — TREM2 ± Aβ42 결합 신뢰도 비교.

읽는 지표
  iptm            두 사슬(TREM2·Aβ42) 사이 인터페이스 신뢰도. **이게 결합 지표다**
  complex_iplddt  인터페이스 잔기들의 pLDDT
  ptm             복합체 전체
  complex_plddt   복합체 전체 pLDDT

검정 설계 (2026-09-04 전환 — 구조 축과 통일)
  ipTM은 구조당 스칼라이므로 쌍 차이 |ipTM_i − ipTM_j| 로 바꿔 쌍 설계에 맞춘다.
  같은 대립형질 안의 쌍 = 모델 진동(음성 대조군), 대립형질 간 쌍 = 변이 신호.
  통계량은 median(신호) − median(노이즈) 이고 **라벨 순열검정**으로 p 를 구한다.
  쌍은 서로 독립이 아니므로 Mann–Whitney 의 p 는 과소평가된다 — src/detectability.py 참조.

⚠️ affinity 수치는 없다. Boltz-2 affinity head는 저분자 전용이고 Aβ42는 펩타이드다.
⚠️ Bret et al. 2026(JCIM)은 Boltz-2가 결합부위 변이에 둔감하다고 보고했다.
   "구분 안 됨"이 나와도 그 자체가 보고할 결과다.

⚠️ CSV는 gene 칼럼을 포함한다 (다사례 확장 대비). 이 스크립트가 수집하지 않은 사례
   (예: PILRA ± NPDC1)의 행은 기존 CSV에서 읽어 그대로 보존한다 — 재실행해도 날아가지 않는다.

실행: python3 src/analyze_boltz.py
출력: results/boltz_iptm.csv, results/boltz_summary.md
"""

from __future__ import annotations

import csv
import json
import re
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detectability import verdict_from_scalars  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEARCH = ROOT / "notebooks" / "boltz_results" / "outputs"
RESULTS = ROOT / "results"

METRICS = ["iptm", "complex_iplddt", "ptm", "complex_plddt", "confidence_score"]


def collect() -> list[dict]:
    rows = []
    for f in sorted(SEARCH.rglob("confidence_*.json")):
        d = json.loads(f.read_text())
        # 경로 예: outputs/TREM2_WT_AB42_seed1/boltz_results_.../predictions/.../confidence_..._model_0.json
        m = re.search(r"TREM2_(WT|R62H)_AB42_seed(\d+)", str(f))
        if not m:
            continue
        allele, seed = m.group(1), int(m.group(2))
        model = int(f.stem.rsplit("_", 1)[-1])
        row = {"gene": "TREM2", "allele": allele, "seed": seed, "model": model}
        row.update({k: d.get(k) for k in METRICS})
        rows.append(row)
    return rows


FIELDNAMES = ["gene", "allele", "seed", "model"] + METRICS


def write_csv(out_csv: Path, rows: list[dict]) -> int:
    """수집한 행을 쓰되, 이번 실행이 다루지 않은 사례의 기존 행은 보존한다.

    이 스크립트는 TREM2만 수집한다. 단순 덮어쓰기를 하면 CSV에 들어 있는
    다른 사례(PILRA ± NPDC1 등)가 조용히 사라진다. 그래서 병합해서 쓴다.
    반환값은 보존한 행 수.
    """
    collected = {r["gene"] for r in rows}
    preserved: list[dict] = []
    if out_csv.exists():
        with out_csv.open(newline="", encoding="utf-8") as f:
            preserved = [r for r in csv.DictReader(f)
                         if r.get("gene") and r["gene"] not in collected]

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
        w.writerows(preserved)
    return len(preserved)


def describe(vals: list[float]) -> str:
    return (f"{st.median(vals):.3f}  "
            f"(min {min(vals):.3f} / max {max(vals):.3f}, n={len(vals)})")


ALLELES = {                      # gene -> (기준 대립형질, 비교 대립형질, 결합 파트너)
    "TREM2": ("WT", "R62H", "Aβ42"),
    "PILRA": ("G78", "R78", "NPDC1"),
}


def summarize(out_csv: Path) -> list[str]:
    """CSV 에 들어 있는 모든 사례로 요약 마크다운을 만든다."""
    rows = list(csv.DictReader(out_csv.open(encoding="utf-8")))
    lines: list[str] = []
    A = lines.append

    A("# Boltz-2 결과 · 결합 신뢰도 감별력 검정\n")
    A("**검정 설계** — 구조 축과 동일하다. 같은 대립형질 안의 쌍 = 모델 진동(음성 대조군),")
    A("대립형질 간 쌍 = 변이 신호. 통계량은 median(신호) − median(노이즈) 이고")
    A("**라벨 순열검정**으로 p 를 구한다 — 쌍이 서로 독립이 아니라 Mann–Whitney 는 과소평가한다.")
    A("`src/detectability.py` 와 `app/app.js` 가 같은 엔진·같은 PRNG 를 써 p 가 자릿수까지 일치한다.\n")

    for gene in sorted({r["gene"] for r in rows}):
        if gene not in ALLELES:
            continue
        a_lab, b_lab, lig = ALLELES[gene]
        g = [r for r in rows if r["gene"] == gene]
        na = sum(1 for r in g if r["allele"] == a_lab)
        nb = sum(1 for r in g if r["allele"] == b_lab)
        seeds = sorted({int(r["seed"]) for r in g})

        A(f"## {gene} ± {lig}\n")
        A(f"구조 {len(g)}개 ({a_lab} {na} + {b_lab} {nb}) · 시드 {seeds} × 모델 5개\n")
        A("| 지표 | 음성 대조군 (같은 대립형질 쌍) | 처리군 (대립형질 간 쌍) | 비 | 순열 p | 판정 |")
        A("|---|---|---|---|---|---|")

        iptm_v = None
        for k in METRICS:
            a = [float(r[k]) for r in g if r["allele"] == a_lab and r.get(k)]
            b = [float(r[k]) for r in g if r["allele"] == b_lab and r.get(k)]
            if len(a) < 2 or len(b) < 2:
                continue
            v = verdict_from_scalars(a, b)
            if k == "iptm":
                iptm_v = v
            A(f"| **{k}** | {v['mn']:.4f} (n={len(v['noise'])}쌍) "
              f"| {v['ms']:.4f} (n={len(v['signal'])}쌍) | {v['ratio']:.2f}× "
              f"| {v['p']:.4f} | {'**감별됨**' if v['det'] else '**구분 안 됨**'} |")

        ia = [float(r["iptm"]) for r in g if r["allele"] == a_lab]
        ib = [float(r["iptm"]) for r in g if r["allele"] == b_lab]
        A(f"\n원시 ipTM 중앙값 — {a_lab} {st.median(ia):.3f} "
          f"(범위 {min(ia):.3f}–{max(ia):.3f}) / "
          f"{b_lab} {st.median(ib):.3f} (범위 {min(ib):.3f}–{max(ib):.3f})\n")

        if iptm_v is None:
            continue
        A(f"검정 방식: {iptm_v['method']} · {iptm_v['n_perm']:,}"
          + (f" (도달 가능한 최소 p = {1/iptm_v['n_perm']:.4f})" if iptm_v["method"] == "exact" else "")
          + "\n")
        if iptm_v["det"]:
            A("→ 변이 쌍 차이가 모델 진동을 통계적으로 **초과했다.**\n")
        else:
            A(f"→ **Boltz-2 가 이 변이를 감별하지 못했다.** (비 {iptm_v['ratio']:.2f}×"
              + (", 1.0 미만 — 변이 쌍 차이가 같은 서열 쌍 차이보다도 작다)" if iptm_v["ratio"] < 1 else ")"))
            A("→ 결합 변화가 없다는 뜻이 아니라 **이 모델이 변이를 감지하지 못한다**는 뜻이다.\n")

    return lines


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    out_csv = RESULTS / "boltz_iptm.csv"

    # 원자료가 있으면 수집해 CSV 를 갱신하고, 없으면 기존 CSV 로 요약만 다시 만든다.
    if SEARCH.exists():
        rows = collect()
        if rows:
            n_kept = write_csv(out_csv, rows)
            print(f"수집 {len(rows)}행"
                  + (f" · 기존 CSV 에서 보존한 타 사례 행 {n_kept}개" if n_kept else ""))
        else:
            print("confidence json 을 못 찾았다 — 기존 CSV 로 요약만 갱신한다.")
    else:
        print(f"원자료 폴더 없음 ({SEARCH.relative_to(ROOT)}) — 기존 CSV 로 요약만 갱신한다.")

    if not out_csv.exists():
        print("results/boltz_iptm.csv 가 없다. 원자료 수집이 먼저다.")
        return 1

    lines = summarize(out_csv)
    A = lines.append

    A("\n## 반드시 함께 보고할 한계\n")
    A("1. **affinity 수치가 아니다.** Boltz-2 affinity head는 저분자 전용이고 "
      "Aβ42는 펩타이드다. ipTM은 결합력이 아니라 **인터페이스 신뢰도**다.")
    A("2. **Bret et al. 2026** (*J Chem Inf Model*) — Boltz-2 affinity가 결합부위 변이에 둔감. "
      "타겟을 바꿔도 분류가 잘 안 바뀌는 사례 보고.")
    A("3. **King et al. 2025** — Boltz-2를 단백질–단백질 친화도로 미세조정해도 "
      "서열 기반 모델보다 성능이 낮았다.")
    A("4. 실험적 검증 없음. 예측 대 예측 비교다.")
    A("5. **PILRA 는 글리칸 부재가 겹친다 (추정).** PILRA 의 리간드 인식은 시알산 매개인데 "
      "Boltz-2 단백질–단백질 예측에는 글리칸이 없다. 미감별이 모델 한계 때문인지 "
      "글리칸 부재 때문인지 이 결과만으로는 구분할 수 없다.")

    A("\n## 근거 문헌\n")
    A("- Zhao et al. 2018, *Neuron* — TREM2가 Aβ 올리고머에 나노몰 결합, AD 변이가 결합 감소")
    A("- Zhong et al. 2018, *Mol Neurodegener* — oAβ1-42 고친화도 결합, 결합 필수 잔기 31–91 (R62 포함)")
    A("- Yeh et al. 2016, *Neuron* — TREM2–APOE/CLU/LDL 결합, 질병 변이가 저해")
    A("- Passaro et al. 2025 — Boltz-2 (MIT, FEP 근접 성능)")
    A("- Bret et al. 2026, *JCIM* — Boltz-2 결합부위 변이 둔감성")
    A("- Rathore N, et al. 2018, *PLoS Genet* 14(11):e1007427 — "
      "PILRA G78R 이 시알산 결합 잔기를 바꿔 NPDC1 등 리간드 결합을 50% 이상 감소")

    (RESULTS / "boltz_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n저장 → {out_csv.relative_to(ROOT)}, results/boltz_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
