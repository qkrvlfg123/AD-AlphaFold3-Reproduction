"""슬라이드 11·13용 시각자료를 만든다.

  figures/slide11_targets.png  — 7개 인과 단백질 중 구조 예측 가능한 3개 (표 형식)
  figures/slide13_noise.png    — 노이즈 바닥 vs 변이 신호 분포 (작은 배수 3패널)

색은 검증된 카테고리 팔레트 슬롯 1·2를 쓴다 (validate_palette.js 전 항목 PASS).
  노이즈 = 파랑 #2a78d6 (배경 역할) · 신호 = 주황 #eb6834 (주목 대상)

실행: python3 src/make_slide_charts.py
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import matplotlib.font_manager as fm                # noqa: E402
from matplotlib.patches import Rectangle            # noqa: E402
import numpy as np                                  # noqa: E402
from detectability import (matrix_from_rmsd_rows,   # noqa: E402
                           verdict_from_matrix, verdict_from_scalars)

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS                       # noqa: E402

FIG_DIR = ROOT / "figures"
RESULTS = ROOT / "results"

# --- 색 역할 (검증된 팔레트) ------------------------------------------------
SURFACE = "#ffffff"        # 슬라이드 배경에 맞춤
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8a86"
GRID = "#e6e5e1"
C_NOISE = "#2a78d6"        # 슬롯 1 파랑
C_SIGNAL = "#eb6834"       # 슬롯 2 주황
C_NAVY = "#1e3a5f"         # 발표 타이틀색

# 한글 폰트 — 실행 환경에 있는 것을 고른다.
# macOS(AppleGothic) 하드코딩이었으나 Windows/Linux 에서 한글이 깨져 후보 목록으로 바꿨다.
_KO_FONTS = ["AppleGothic", "Malgun Gothic", "NanumGothic", "Noto Sans KR", "Gulim"]
_available = {f.name for f in fm.fontManager.ttflist}
_font = next((f for f in _KO_FONTS if f in _available), "DejaVu Sans")
if _font == "DejaVu Sans":
    print("  ⚠️ 한글 폰트를 찾지 못했다 — 라벨이 깨질 수 있다")

plt.rcParams.update({
    "font.family": _font,
    "axes.unicode_minus": False,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


# ===========================================================================
# 슬라이드 11 — 구조 예측 대상 표
# ===========================================================================
# 논문 Fig 4A(방향) + Fig 4B·S5(구조 유무) 기준
SEVEN = [
    # (단백질, 방향, 미스센스 변이 있음?, 변이 표기, 논문 그림)
    ("CD33",  "위험 ↑", True,  "R69G",  "Fig 4B"),
    ("PILRA", "위험 ↑", True,  "R78G",  "Fig S5(a)"),
    ("TREM2", "위험 ↓", True,  "R62H",  "Fig S5(b)"),
    ("PILRB", "위험 ↑", False, "—",     "없음"),
    ("RET",   "위험 ↑", False, "—",     "없음"),
    ("CD55",  "위험 ↓", False, "—",     "없음"),
    ("EPHA1", "위험 ↓", False, "—",     "없음"),
]


def slide11() -> Path:
    fig, ax = plt.subplots(figsize=(10, 5.4), dpi=200)
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 9.4)

    ax.text(0.15, 9.0, "MR-SPI가 지목한 인과 단백질 7개 중, 구조를 예측할 수 있는 것은 3개",
            fontsize=13.5, color=C_NAVY, fontweight="bold", va="center")
    ax.text(0.15, 8.45,
            "나머지 4개는 미스센스 변이 자체가 없어 변이형 서열을 만들 수 없다",
            fontsize=10, color=INK_2, va="center")

    cols = [0.35, 2.35, 4.35, 6.6, 8.4]
    heads = ["단백질", "AD 위험 방향", "미스센스 변이", "변이 표기", "논문 그림"]

    y0, dy = 7.5, 0.86
    for x, h in zip(cols, heads):
        ax.text(x, y0, h, fontsize=10, color=INK_2, fontweight="bold", va="center")
    ax.plot([0.15, 9.85], [y0 - 0.34] * 2, color=INK_2, lw=1.2)

    for i, (gene, direction, has_mis, mut, figref) in enumerate(SEVEN):
        y = y0 - dy * (i + 1)
        if has_mis:
            ax.add_patch(Rectangle((0.15, y - 0.35), 9.7, 0.7,
                                   facecolor="#fdf0e9", edgecolor="none", zorder=0))
        ink = INK if has_mis else INK_MUTED
        weight = "bold" if has_mis else "normal"

        ax.text(cols[0], y, gene, fontsize=11.5, color=ink, fontweight=weight, va="center")
        ax.text(cols[1], y, direction, fontsize=10.5, color=ink, va="center")
        ax.text(cols[2], y, "○" if has_mis else "×", fontsize=13,
                color=C_SIGNAL if has_mis else INK_MUTED,
                fontweight="bold", va="center")
        ax.text(cols[3], y, mut, fontsize=10.5, color=ink, fontweight=weight, va="center")
        ax.text(cols[4], y, figref, fontsize=10, color=ink, va="center")

        if i == 2:  # 3개와 4개 사이 구분선
            ax.plot([0.15, 9.85], [y - dy / 2] * 2, color=GRID, lw=1.4)

    ax.text(0.15, 0.55,
            "초록은 \"seven proteins with structural alterations\"로 기술되어 있어 "
            "7개 모두 구조가 보고된 것처럼 읽힌다.",
            fontsize=9, color=INK_2, va="center")
    ax.text(0.15, 0.15,
            "출처: Yao et al. 2024, Figure 4A · 4B · S5",
            fontsize=8, color=INK_MUTED, va="center")

    out = FIG_DIR / "slide11_targets.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return out


# ===========================================================================
# 슬라이드 13 — 노이즈 바닥 vs 신호
# ===========================================================================
def slide13() -> Path:
    rows = list((RESULTS / "rmsd_pairs.csv").open() and
                csv.DictReader((RESULTS / "rmsd_pairs.csv").open()))
    genes = [t["gene"] for t in TARGETS]

    data = {}
    for g in genes:
        rows_g = [r for r in rows if r["gene"] == g]
        D, k = matrix_from_rmsd_rows(rows_g)
        v = verdict_from_matrix(D, k, k)
        data[g] = (v["noise"], v["signal"], v["p"])

    xmax = max(max(v[0] + v[1]) for v in data.values()) * 1.18

    fig, axes = plt.subplots(len(genes), 1, figsize=(10, 6.4), dpi=200, sharex=True)
    fig.subplots_adjust(left=0.15, right=0.83, top=0.80, bottom=0.13, hspace=0.42)

    rng = np.random.default_rng(0)  # 지터 재현성 고정

    for ax, g in zip(axes, genes):
        noise, sig, p = data[g]
        sep = p < 0.05

        for vals, y, color, lab in ((noise, 1.0, C_NOISE, "노이즈"),
                                    (sig, 0.0, C_SIGNAL, "신호")):
            jit = rng.uniform(-0.17, 0.17, len(vals))
            ax.scatter(vals, np.full(len(vals), y) + jit, s=34, color=color,
                       alpha=0.55, linewidths=1.2, edgecolors=SURFACE, zorder=3)
            m = st.median(vals)
            ax.plot([m, m], [y - 0.32, y + 0.32], color=color, lw=2.6, zorder=4)
            ax.text(m, y + 0.42, f"{m:.3f}", fontsize=8.5, color=color,
                    ha="center", va="bottom", fontweight="bold")

        ax.set_ylim(-0.62, 1.72)
        ax.set_xlim(0, xmax)
        ax.set_yticks([1.0, 0.0])
        ax.set_yticklabels(["노이즈\n(같은 서열)", "신호\n(WT vs 변이형)"],
                           fontsize=9, color=INK_2)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", colors=INK_2, labelsize=9)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRID)

        mut = next(f"{t['wt']}{t['pos']}{t['mut']}" for t in TARGETS if t["gene"] == g)
        ax.text(-0.005, 1.98, f"{g}  {mut}", transform=ax.get_yaxis_transform(),
                fontsize=11.5, color=INK, fontweight="bold", va="top", ha="right")

        verdict = "분리됨" if sep else "구분 안 됨"
        vcolor = C_SIGNAL if sep else "#b02020"
        ax.text(1.02, 0.72, verdict, transform=ax.transAxes, fontsize=10,
                color=vcolor, fontweight="bold", va="center")
        ax.text(1.02, 0.30, f"순열 p = {p:.4f}", transform=ax.transAxes, fontsize=9,
                color=INK_2, va="center")

    axes[-1].set_xlabel("변이 잔기 8 Å 이내 Cα 국소 RMSD (Å)",
                        fontsize=10, color=INK_2, labelpad=8)

    fig.text(0.15, 0.955,
             "변이로 생긴 구조 차이가 모델 자체의 변동을 넘는가",
             fontsize=13.5, color=C_NAVY, fontweight="bold", ha="left")
    fig.text(0.15, 0.905,
             "노이즈 = 같은 서열의 모델 쌍 20개 · 신호 = WT × 변이형 모델 쌍 25개 · "
             "세로선 = 중앙값 · 라벨 순열검정 단측 (전수 126분할)",
             fontsize=9, color=INK_2, ha="left")
    fig.text(0.15, 0.032,
             "CD33은 논문이 본문 Figure 4B에 대표로 실은 사례다. "
             "AlphaFold Server(AF3), 전장 서열, seed 1, job당 모델 5개.",
             fontsize=8.5, color=INK_MUTED, ha="left")

    out = FIG_DIR / "slide13_noise.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ===========================================================================
# 슬라이드 14 — Boltz-2 결합 축 (ipTM)
# ===========================================================================
def slide14() -> Path:
    import csv as _csv
    rows = list(_csv.DictReader((RESULTS / "boltz_iptm.csv").open()))
    wt  = [float(r["iptm"]) for r in rows if r["allele"] == "WT"]
    mut = [float(r["iptm"]) for r in rows if r["allele"] == "R62H"]
    # 구조 축과 같은 음성 대조군 설계: 같은 대립형질 안의 쌍 = 노이즈, 대립형질 간 쌍 = 신호
    v = verdict_from_scalars(wt, mut)
    noise, signal, p = v["noise"], v["signal"], v["p"]

    fig, ax = plt.subplots(figsize=(10, 4.6), dpi=200)
    fig.subplots_adjust(left=0.17, right=0.80, top=0.68, bottom=0.22)
    rng = np.random.default_rng(0)

    for vals, y, color, lab in ((noise, 1.0, C_NOISE, "노이즈"),
                                (signal, 0.0, C_SIGNAL, "신호")):
        jit = rng.uniform(-0.17, 0.17, len(vals))
        ax.scatter(vals, np.full(len(vals), y) + jit, s=42, color=color,
                   alpha=0.6, linewidths=1.2, edgecolors=SURFACE, zorder=3)
        m = st.median(vals)
        ax.plot([m, m], [y - 0.30, y + 0.30], color=color, lw=2.8, zorder=4)
        ax.text(m, y + 0.40, f"{m:.3f}", fontsize=9.5, color=color,
                ha="center", va="bottom", fontweight="bold")

    ax.set_ylim(-0.62, 1.66)
    ax.set_yticks([1.0, 0.0])
    ax.set_yticklabels([f"기준 반복\n(같은 대립형질 쌍, {len(noise)})",
                        f"변이 비교\n(대립형질 간 쌍, {len(signal)})"],
                       fontsize=10, color=INK_2)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=INK_2, labelsize=9)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("ipTM 쌍 차이 |Δ| — 작을수록 두 구조의 인터페이스 신뢰도가 비슷하다",
                  fontsize=10, color=INK_2, labelpad=8)

    ax.text(1.03, 0.70, "구분 안 됨", transform=ax.transAxes, fontsize=11,
            color="#b02020", fontweight="bold", va="center")
    ax.text(1.03, 0.32, f"순열 p = {p:.2f}", transform=ax.transAxes, fontsize=9.5,
            color=INK_2, va="center")
    ax.text(1.03, 0.10, f"비 {v['ratio']:.2f}×", transform=ax.transAxes, fontsize=9.5,
            color=INK_2, va="center")

    fig.text(0.17, 0.93, "Boltz-2도 이 변이를 감별하지 못했다",
             fontsize=13.5, color=C_NAVY, fontweight="bold", ha="left")
    fig.text(0.17, 0.855,
             f"시드 3 × 모델 5 = 조건당 {len(wt)}개 구조 · 세로선 = 중앙값 · "
             f"라벨 순열검정 단측 ({v['method']}, {v['n_perm']:,})",
             fontsize=9, color=INK_2, ha="left")
    fig.text(0.17, 0.055,
             "복합체 예측 자체는 성공했다 (원시 ipTM 중앙값 0.85). 변이 감별력만 없다.  "
             "비 1.0 미만 = 변이 쌍 차이가 같은 서열 쌍 차이보다도 작다.",
             fontsize=8.5, color=INK_MUTED, ha="left")

    out = FIG_DIR / "slide14_boltz_iptm.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for f in (slide11(), slide13(), slide14()):
        print(f"  → {f.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
