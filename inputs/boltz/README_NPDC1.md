# Boltz-2 실행 — PILRA ± NPDC1

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

조건당 시드 3개 × diffusion sample 5개 = **15 구조**.
TREM2 실행과 동일한 설계이므로 같은 노이즈 바닥 검정을 그대로 적용할 수 있다.

```bash
boltz predict inputs/boltz/PILRA_G78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/PILRA_G78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/PILRA_G78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
boltz predict inputs/boltz/PILRA_R78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/PILRA_R78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/PILRA_R78_NPDC1.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
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
