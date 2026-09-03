"""Boltz-2 입력 YAML을 만든다 — PILRA ± NPDC1 복합체.

목적
  TREM2 확장과 같은 설계다. 구조 축에서 감별되는지와 별개로 **결합**이 변하는지를 본다.
  PILRA는 시알산 매개로 리간드와 결합하고, 78번 잔기가 바로 그 시알산 결합 부위에 있다.

근거
  · Rathore N, et al. (2018) PLoS Genetics 14(11):e1007427
    doi:10.1371/journal.pgen.1007427
    PILRA G78R이 시알산 결합 잔기를 바꿔 NPDC1 등 리간드 결합을 50% 이상 감소시킨다.

⚠️ 변이 방향 — 반드시 확인할 것
  문헌은 'G78R'이라 부르지만 UniProt 정본(Q9UKJ1)은 78번이 이미 **R**이다. 따라서
      정상형(G78) = 정본에 R78G 치환을 적용한 서열   ← 치환이 필요한 쪽
      변이형(R78) = 정본 그대로                       ← 치환이 없는 쪽
  기존 저장소 파일(inputs/seqs/PILRA_dom_WT.fasta = R78)과 라벨이 정반대이므로
  WT/MUT 표기를 쓰지 않고 잔기 이름(G78 / R78)으로만 명명한다.

구간 선택
  PILRA  32-150  Ig-like V-type 도메인. TARGETS의 construct 값을 그대로 쓴다.
                 결정구조 4NFB(78R) / 3WUZ(78G)가 이 구간과 정확히 일치한다.
  NPDC1  35-181  세포외 도메인. 신호펩타이드(1-34) 다음부터 막관통(182-202) 직전까지.
                 UniProt feature에서 유도하고 assert로 고정한다.

⚠️ 한계 (반드시 함께 보고할 것)
  · 글리칸 부재 — PILRA의 리간드 인식은 시알산 매개인데 Boltz-2 단백질-단백질
    예측에는 글리칸이 없다. 결과가 '구분 안 됨'이어도 그것이 모델의 감별 한계 때문인지
    글리칸이 빠졌기 때문인지 구분할 수 없다.
  · affinity 수치가 아니다 — Boltz-2 affinity head는 저분자 전용이다.
    얻는 값은 결합력이 아니라 ipTM(인터페이스 신뢰도)다.
  · NPDC1 세포외 도메인 안에 무질서 구간(138-175, 38 aa)이 있다. ipTM 분산이
    TREM2 대비 과도하면 35-137로 잘라 재생성하는 것을 검토한다.

실행: python src/make_boltz_npdc1_inputs.py
출력: inputs/boltz/PILRA_{G78,R78}_{NPDC1,alone}.yaml, inputs/boltz/README_NPDC1.md
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS, fetch_uniprot, apply_mutation  # noqa: E402

OUT = ROOT / "inputs" / "boltz"
CACHE_DIR = ROOT / "inputs" / "seqs" / "_uniprot"

NPDC1_ACC = "Q9NQX5"
NPDC1_LEN = 325
PILRA_LEN = 303
NPDC1_ECTO = (35, 181)     # UniProt feature에서 유도한 값. 아래 assert로 고정한다

SEEDS = [1, 2, 3]          # TREM2 실행과 동일 — 노이즈 바닥 검정용
DIFFUSION_SAMPLES = 5      # 조건당 15구조 (시드 3 x 샘플 5)


def fetch_uniprot_features(acc: str) -> dict:
    """UniProt JSON을 받는다. fetch_uniprot은 FASTA만 주므로 feature용으로 따로 받는다."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = CACHE_DIR / f"{acc}.json"
    if not cached.exists():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        print(f"  UniProt JSON 다운로드: {acc}")
        with urllib.request.urlopen(url, timeout=30) as r:
            cached.write_text(r.read().decode(), encoding="utf-8")
    else:
        print(f"  캐시 사용: {acc}.json")
    return json.loads(cached.read_text(encoding="utf-8"))


def ectodomain_range(data: dict, gene: str) -> tuple[int, int]:
    """신호펩타이드 끝 다음 ~ 막관통 시작 직전.

    app.js의 buildSequences가 화면에서 UniProt feature로 구간을 잡는 것과 같은 유도다.
    """
    feats = data.get("features", [])
    sig = next((f for f in feats if f["type"] == "Signal"), None)
    tm = next((f for f in feats if f["type"] == "Transmembrane"), None)
    if sig is None or tm is None:
        raise AssertionError(f"{gene}: Signal 또는 Transmembrane 주석이 없다 — 구간 유도 불가")
    return sig["location"]["end"]["value"] + 1, tm["location"]["start"]["value"] - 1


def yaml_for(chains: list[tuple[str, str]]) -> str:
    """Boltz-2 입력 YAML. chains = [(id, sequence), ...]

    affinity 블록은 넣지 않는다 — 단백질-단백질이라 affinity head 대상이 아니다.
    형식은 make_boltz_inputs.py와 동일하며, 그 형식으로 예측 30건이 성공했다.
    """
    lines = ["version: 1", "sequences:"]
    for cid, seq in chains:
        lines += ["  - protein:", f"      id: {cid}", f'      sequence: "{seq}"']
    return "\n".join(lines) + "\n"


def build_readme(cmds: str) -> str:
    return f"""# Boltz-2 실행 — PILRA ± NPDC1

## 무엇을 왜 하는가

TREM2 확장과 같은 설계다. 구조 축에서 감별되는지와 별개로 **결합**이 변하는지를 본다.
PILRA는 시알산 매개로 리간드와 결합하고, 78번 잔기가 그 시알산 결합 부위에 있다.

**근거** — Rathore N, et al. (2018) *PLoS Genetics* 14(11):e1007427
doi:[10.1371/journal.pgen.1007427](https://doi.org/10.1371/journal.pgen.1007427)
PILRA G78R이 시알산 결합 잔기를 바꿔 NPDC1 등 리간드 결합을 50% 이상 감소시킨다.

## ⚠️ 변이 방향 — 파일명을 읽기 전에 반드시 확인할 것

문헌은 'G78R'이라 부르지만 **UniProt 정본(Q9UKJ1)은 78번이 이미 R**이다.

| 이 디렉터리의 명칭 | 78번 잔기 | 생성 방법 | Rathore 2018 기준 |
|---|---|---|---|
| `PILRA_G78_*` | G | 정본에 R78G 치환 적용 | 정상형(common) |
| `PILRA_R78_*` | R | **정본 그대로** | 변이형(G78R의 산물) |

기존 저장소 파일과 라벨이 정반대다 — `inputs/seqs/PILRA_dom_WT.fasta`가 R78,
`PILRA_dom_R78G.fasta`가 G78이다. 혼동을 막기 위해 여기서는 WT/MUT 표기를 쓰지 않고
잔기 이름만으로 명명했다.

## 입력 구간

| 사슬 | 단백질 | UniProt | 구간 | 길이 | 근거 |
|---|---|---|---|---|---|
| A | PILRA | Q9UKJ1 | 32–150 | 119 aa | Ig-like V-type 도메인. 결정구조 4NFB(78R)/3WUZ(78G)와 정확히 일치 |
| B | NPDC1 | Q9NQX5 | 35–181 | 147 aa | 세포외 도메인. 신호펩타이드(1–34) 이후 ~ 막관통(182–202) 직전 |

복합체 = 266 aa (참고: TREM2 + Aβ42는 154 aa)

## 생성 파일

| 파일 | 사슬 | 용도 |
|---|---|---|
| `PILRA_G78_NPDC1.yaml` | A+B | 정상형 복합체 |
| `PILRA_R78_NPDC1.yaml` | A+B | 변이형 복합체 |
| `PILRA_G78_alone.yaml` | A | 단독 대조 |
| `PILRA_R78_alone.yaml` | A | 단독 대조 |

## 실행

조건당 시드 {len(SEEDS)}개 × diffusion sample {DIFFUSION_SAMPLES}개 = **15 구조**.
TREM2 실행과 동일한 설계이므로 같은 노이즈 바닥 검정을 그대로 적용할 수 있다.

```bash
{cmds}
```

## ⚠️ 반드시 함께 보고할 한계

1. **글리칸 부재 — 이번 예측의 핵심 한계.**
   PILRA의 리간드 결합은 **시알산 매개**인데 Boltz-2 단백질–단백질 예측에는
   **글리칸이 없다.** 따라서 결과가 '구분 안 됨'으로 나와도 그것이
   **모델의 감별 한계** 때문인지 **글리칸이 빠졌기 때문**인지 구분할 수 없다.
   → 이는 이 도구의 핵심 논지를 오히려 강화하는 사례다. 계산 예측이 '차이 없음'을
   냈을 때 그것은 생물학적 결론이 아니라 **모델이 답할 수 없는 영역**이라는 뜻이다.
   실험 문헌(Rathore 2018)은 50% 이상의 결합 감소를 보고했다.

2. **affinity 수치가 아니다.** Boltz-2 affinity head는 저분자 전용이다.
   얻는 값은 결합력이 아니라 **ipTM(인터페이스 신뢰도)** 다.

3. **Bret et al. 2026** (*J Chem Inf Model*) — Boltz-2가 결합부위 변이에 둔감하다는 보고.
   ipTM 차이도 노이즈 바닥 검정을 거쳐야 한다.

4. **NPDC1 무질서 구간.** 세포외 도메인 35–181 안에 무질서 구간 138–175(38 aa)가 있다.
   시드 간 흔들림이 ipTM 분산을 키울 수 있다. TREM2 대비 분산이 과도하면
   35–137(103 aa)로 잘라 재생성하는 것을 검토한다.
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # --- PILRA 도메인 G78 / R78 -------------------------------------------
    t = next(x for x in TARGETS if x["gene"] == "PILRA")
    s, e = t["construct"]                       # 32-150
    canon = fetch_uniprot(t["acc"])             # 정본 = 78R
    assert len(canon) == PILRA_LEN, f"PILRA 길이 {len(canon)} — {PILRA_LEN}을 기대"
    assert canon[t["pos"] - 1] == "R", "PILRA 78번이 R이 아니다 — isoform 확인 필요"

    # 정본이 R78이므로 '정상형 G78'을 만들려면 치환이 필요하다 (변이형은 정본 그대로)
    g78_full = apply_mutation(canon, t["pos"], t["wt"], t["mut"], t["gene"])
    dom = {"R78": canon[s - 1 : e], "G78": g78_full[s - 1 : e]}

    assert len(dom["R78"]) == len(dom["G78"]) == e - s + 1
    assert sum(a != b for a, b in zip(dom["R78"], dom["G78"])) == 1, "상이 잔기가 1개가 아님"
    idx = t["pos"] - s
    assert dom["R78"][idx] == "R" and dom["G78"][idx] == "G", "구간 내 78번 위치가 어긋남"
    print(f"PILRA {s}-{e}  {len(dom['R78'])} aa  G78(정상형) / R78(변이형·정본)  OK")

    # --- NPDC1 세포외 도메인 ----------------------------------------------
    npdc1_full = fetch_uniprot(NPDC1_ACC)
    assert len(npdc1_full) == NPDC1_LEN, f"NPDC1 길이 {len(npdc1_full)} — {NPDC1_LEN}을 기대"
    ns, ne = ectodomain_range(fetch_uniprot_features(NPDC1_ACC), "NPDC1")
    assert (ns, ne) == NPDC1_ECTO, f"NPDC1 세포외 구간 {ns}-{ne} — {NPDC1_ECTO}를 기대 (UniProt 갱신 확인)"
    npdc1 = npdc1_full[ns - 1 : ne]
    assert len(npdc1) == ne - ns + 1 == 147
    print(f"NPDC1 {ns}-{ne}  {len(npdc1)} aa  세포외 도메인 (막관통 182-202 제외)  OK")

    # --- YAML 생성 ---------------------------------------------------------
    jobs = []
    for allele in ("G78", "R78"):
        for with_lig, tag in ((True, "NPDC1"), (False, "alone")):
            chains = [("A", dom[allele])] + ([("B", npdc1)] if with_lig else [])
            name = f"PILRA_{allele}_{tag}"
            (OUT / f"{name}.yaml").write_text(yaml_for(chains), encoding="utf-8")
            jobs.append((name, sum(len(c[1]) for c in chains)))

    for n, ln in jobs:
        print(f"  {n:24s} {ln:>4d} aa")

    # --- 실행 안내 ---------------------------------------------------------
    cmds = "\n".join(
        f"boltz predict inputs/boltz/{n}.yaml \\\n"
        f"    --out_dir outputs/boltz --use_msa_server \\\n"
        f"    --diffusion_samples {DIFFUSION_SAMPLES} --seed {seed} \\\n"
        f"    --output_format mmcif"
        for n, _ in jobs if n.endswith("NPDC1")
        for seed in SEEDS
    )
    (OUT / "README_NPDC1.md").write_text(build_readme(cmds), encoding="utf-8")

    print(f"\n출력: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
