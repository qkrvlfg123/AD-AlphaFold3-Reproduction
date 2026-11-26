# Boltz-2 실행: TREM2 ± Aβ42

## 목적과 배경

논문은 변이의 구조 변화를 봄. 재현 결과 TREM2 R62H의 구조 차이는 작았고(국소 RMSD 중앙값 차 0.09 Å) 판정이 지표에 따라 갈렸음.
실험 문헌은 TREM2 변이가 Aβ 결합을 떨어뜨린다고 보고함.
그래서 결합 축을 따로 봄. Boltz-2로 TREM2 + Aβ42 복합체를 예측하고 인터페이스 신뢰도(ipTM)를 비교함.

## 입력

| 파일 | 내용 | 길이 |
|---|---|---|
| `TREM2_WT_AB42.yaml` | TREM2 야생형 + Aβ42 | 154 aa |
| `TREM2_R62H_AB42.yaml` | TREM2 R62H + Aβ42 | 154 aa |
| `TREM2_WT_alone.yaml` | TREM2 야생형 단독 (대조) | 112 aa |
| `TREM2_R62H_alone.yaml` | TREM2 R62H 단독 (대조) | 112 aa |

- TREM2는 Ig 도메인 19~130 구간을 씀(결정구조 5UD8, 5ELI 구간이며 결합 필수 잔기 31~91을 포함).
- Aβ42는 APP P05067 672~713임.
- 단독 입력 2개는 대조용으로 만든 것이며, 아래 실행 명령과 저장소의 결과에는 복합체 입력만 포함됨.

## 읽을 지표

- ipTM (0~1): 두 사슬 사이 인터페이스 신뢰도. 결합 축에서 비교하는 지표임
- interface pLDDT: 접촉면 잔기들의 확신도
- affinity 수치는 나오지 않음. Boltz-2 affinity head는 저분자 전용이고 Aβ42는 펩타이드임. `properties: affinity` 블록을 넣지 않았음

## 실행 (Colab T4 권장)

```bash
!pip install -q "boltz[cuda]==2.2.1"
# 런타임 재시작 후

boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
```

## 해석할 때 함께 보고할 것

- Bret et al. 2026 (J Chem Inf Model)은 Boltz-2 affinity가 결합부위 변이에 둔감하다고 보고함. 따라서 WT와 R62H의 ipTM 차이도 음성 대조군과 비교하는 검정을 거쳐야 함.
- 시드 [1, 2, 3] × diffusion sample 5개로 분포를 만들고, 구조 축과 같은 라벨 순열검정(`src/detectability.py`)을 적용함.
- 처음에는 Mann-Whitney U 검정을 썼고 2026-09-03에 순열검정으로 바꿈. 쌍이 서로 독립이 아니어서 Mann-Whitney의 p가 과소평가되기 때문임.
- 검정 결과는 `results/boltz_summary.md`에 있음.

## 근거 문헌

- Zhao et al. 2018, *Neuron*: TREM2가 Aβ 올리고머에 나노몰 결합, AD 변이가 결합 감소
- Zhong et al. 2018, *Mol Neurodegener*: oAβ1-42 고친화도 결합, 필수 잔기 31~91
- Yeh et al. 2016, *Neuron*: TREM2와 APOE, CLU, LDL의 결합, 질병 변이가 저해
- Passaro et al. 2025, *bioRxiv*: Boltz-2 (MIT 라이선스, FEP 근접 성능)
- Bret et al. 2026, *JCIM*: Boltz-2의 결합부위 변이 둔감성 (한계)
