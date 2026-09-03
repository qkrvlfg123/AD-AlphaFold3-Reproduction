# PROGRESS — 교차검증 종합 판정 추가

**상태: 1단계 완료. 2단계 진행 중 — PILRA 결합 축 확보(2026-09-03), CD33만 미측정.**
작성일 2026-09-02 · 대상 저장소 `AD-AlphaFold3-Reproduction`

변경 파일: `app/app.js` · `app/index.html` (데이터·CSV·파이썬 파이프라인 무변경)

---

## 1. 목표

구조 축(AlphaFold3)과 결합 축(Boltz-2)의 판정 결과를 **하나의 종합 판정으로 교차검증**한다.

현재는 두 축이 각각 독립된 화면에서 따로 판정되고, 리포트 화면(s5)에서도 축별 행으로
평면 나열될 뿐이다. 사례 단위로 "구조는 어떻고 결합은 어떤가"를 묶어 읽을 수 있는 자리가 없다.

- 구조 축 — 변이 잔기 8 A 이내 Ca 국소 RMSD, 재현 노이즈 바닥 대비 초과 여부
- 결합 축 — TREM2 + Ab42 복합체 ipTM, 재현 변동 대비 차이 유무

두 축이 서로를 검증한다: 구조가 안 변해도 결합이 변할 수 있고, 그 반대도 가능하다.
어느 축에서도 감별되지 않는 경우는 **"효과 없음"이 아니라 "계산 모델의 감별 한계"** 로 보고한다.

---

## 2. 방향 — s5 리포트 화면에 종합 판정 (2x2 규칙)

### 붙이는 위치

**새 화면을 만들지 않고 기존 s5(리포트)를 확장한다.**

- `app/index.html:398` 의 `<div class="panel"><h3>판정 요약</h3>` **위에** 종합 판정 패널을 신설
- 기존 `#rTable`(축별 상세 4행)은 근거 표로 그대로 유지
- 위계: 상단 = 사례별 종합 결론 / 하단 = 축별 근거 상세

선택 근거

- s5 는 `RMSD` 와 `IPTM` 두 원자료를 모두 로드하는 유일한 화면이다 (`app/app.js:427-428`)
- 화면 설명이 이미 "전 사례 판정 요약" 이다
- 네비 흐름 `홈 - 분석 - 구조 검증 - 결합 인터페이스 검증 - 리포트 - 방법론` 이
  "입력 -> 축별 검증 -> 종합 -> 방법" 서사다. s4/s5 사이 삽입이나 s6 뒤 추가는 이 흐름을 깬다
- 라우팅 훅 `if(id==='s5') renderReport();` (`app/app.js:452, 466`) 를 그대로 재사용 -> 라우팅 변경 없음
- 새 화면 s7 신설은 비채택: 네비 7개로 늘고 `app/slides/` 스크린샷 10장 및
  `DOC_CALLOUTS` 구성과 어긋나며, s5 와 역할이 대부분 겹친다

### 2x2 판정 규칙

구조 축 판정 x 결합 축 판정의 네 조합.

| | 결합 축 — 변동성 초과 | 결합 축 — 변동성 내 |
|---|---|---|
| **구조 축 — 변동성 초과** | 양축 감별 | 구조만 감별 |
| **구조 축 — 변동성 내** | 결합만 감별 | 양축 미감별 |

5번째 상태로 **결합 축 미측정** 을 둔다 (CD33 / PILRA). 빈 값을 "차이 없음" 으로
렌더하지 않는다 — 미측정과 미감별은 다른 결론이다.

### 규칙 적용 시 반드시 다뤄야 할 축 간 비대칭

|  | 구조 축 | 결합 축 |
|---|---|---|
| 대상 | 3개 사례 (인자로 전환) | TREM2 하드코딩, 인자 없음 |
| 검정 | `mannWhitneyU(sig, noise, 'greater')` **단측** | `mannWhitneyU(mut, wt, 'two-sided')` **양측** |
| 귀무가설 | 변이 편차가 재현 변동을 초과하지 않는다 | 두 군의 ipTM 분포가 같다 |
| 단위 | A (클수록 신호) | ipTM (무차원, 방향 없음) |
| 효과크기 | `ms/mn` 비율 | `spread` 대비 `d` |

두 축의 `det` 불리언을 단순 AND/OR 하면 **서로 다른 귀무가설을 뭉개게 된다.**
2x2 는 두 판정을 합성하지 않고 **나란히 배치해 읽는 매트릭스** 로 설계한다.

---

## 3. 단계 계획

### 1단계 — TREM2 (데이터 있음, 먼저) — 완료

TREM2 R62H 는 **세 사례 중 유일하게 양축 데이터가 갖춰져 있다.**

| 사례 | 구조 축 | 결합 축 |
|---|---|---|
| CD33 R69G | 있음 — noise 20 / signal 25 쌍 | 없음 |
| PILRA R78G | 있음 — noise 20 / signal 25 쌍 | 없음 |
| **TREM2 R62H** | **있음 — noise 20 / signal 25 쌍** | **있음 — WT 15 / R62H 15** |

추가 예측 없이 기존 CSV 만으로 종합 판정을 구현할 수 있다. 비용 0.

작업 항목

1. 계산부를 순수 함수로 분리 — 현재 동일 계산이 **3곳에 중복** 되어 있다
   (`renderStructureAxis` `app/app.js:240-295` / `renderBindingAxis` `app/app.js:296-336`
   / `renderReport` `app/app.js:426-453`). 종합 판정을 그냥 추가하면 4번째 복사본이 생긴다
   - `computeStructureVerdict(gene)` -> `{p, det, mn, ms, ratio, n}`
   - `computeBindingVerdict()` -> `{p, det, mw, mm, d, spread, n}`
2. 2x2 매트릭스 패널을 `#rTable` 위에 렌더
3. 결합 축 미측정 상태를 CD33 / PILRA 에 명시적으로 표시
4. `DOC_CALLOUTS.s5` (`app/app.js:660`) 에 새 패널 셀렉터 추가

**실측 판정 (앱 구현으로 재계산, `results/slide_tables.md` scipy 값과 일치)**

| 사례 | 구조 축 p | 판정 | 결합 축 p | 판정 | 결합 파트너 |
|---|---|---|---|---|---|
| CD33 R69G | 0.53 | 변동성 내 | — | 미측정 | 미정 |
| PILRA R78 | 6.0e-09 | 변동성 초과 | 0.85 | 변동성 내 | NPDC1 |
| TREM2 R62H | 0.0054 | 변동성 초과 | 0.30 | 변동성 내 | Aβ42 |

PILRA 결합 축 p 는 앱의 정규근사값이다. 시드 1–2 (n = 10 + 10) 기준으로
정규근사 0.8501 / 정확검정 0.8534 (U = 47.0, G78 기준). 시드 1만 있던 단계에서는
정규근사 0.8345 / 정확검정 0.8413 (U = 11.0) 이었다. 네 경우 모두 판정은 변동성 내다.

TREM2 는 **구조 축 초과 x 결합 축 변동성 내** 칸에 배치된다. 단 효과크기는 작다 —
Δ +0.089 Å (1.53배), noise 범위 0.069-1.250 Å 와 signal 범위 0.163-1.187 Å 가 크게 겹친다.

참고: TREM2 는 결합 축에서 ipTM WT 0.851 vs R62H 0.824, p = 0.3 으로
**변동성 내** 다 (`results/boltz_summary.md`). Bret et al. 2026 의 결합부위 변이 둔감성 보고와 일치.

### 2단계 — CD33 / PILRA (나중)

두 사례는 결합 축 데이터가 **아예 없다.** 만들려면 예측 실행이 선행된다.

선결 과제

1. **결합 파트너 조사** — TREM2 는 Ab42 라는 문헌 근거가 확립돼 있다
   (Zhao 2018 Neuron / Zhong 2018 Mol Neurodegener). CD33 / PILRA 는 상응하는
   파트너를 문헌에서 확정해야 한다. 임의 선택은 안 된다
2. **Boltz-2 예측 실행** — 사례당 WT/변이 x seed 3 x model 5 = 30 구조. GPU 필요
3. **입력 구간 결정** — TREM2 는 Ig 도메인 19-130 을 썼다 (`src/make_boltz_inputs.py`).
   각 사례별로 인터페이스 지표를 오염시키지 않는 구간을 따로 정해야 한다

스키마 변경 필요

`results/boltz_iptm.csv` 와 `app/data/boltz_iptm.csv` 는 현재 `allele` 값이 `WT` / `R62H` 뿐이고
**`gene` 칼럼이 없다.** 단일 사례 전제 스키마이므로 다사례 확장 시 컬럼 추가가 불가피하다.
`CASES` 상수(`app/app.js:15-25`)에도 결합 축 관련 필드가 없다
(`wtCif` / `mutCif` / `img` 는 모두 구조 축용).

---

## 4. 현황 요약 — 이미 있는 것

조사 결과 결합 축은 데이터 / 계산 / 화면 / 문구까지 전부 구현된 상태다. 신규 구축이 아니라 기존 자산 연결 작업이다.

| 자산 | 위치 |
|---|---|
| 결합 축 화면 (s4) | `app/index.html:364-392` — `#bScore` `#bPlot` `#bViewer` |
| 결합 축 계산 | `app/app.js:296-336` `renderBindingAxis()` |
| 구조 축 화면 (s3) | `app/index.html:334-362` — `#sScore` `#sTitle` `#sCond` `#sPlot` `#sViewer` |
| 구조 축 계산 | `app/app.js:240-295` `renderStructureAxis(gene)` |
| 통계 (JS 자체 구현) | `app/app.js:32` `mannWhitneyU` + `app/app.js:60` `normCdf` |
| 원자료 (구조) | `app/data/rmsd_pairs.csv` — 135행, `results/` 사본과 동일 |
| 원자료 (결합) | `app/data/boltz_iptm.csv` — 30행, `results/` 사본과 동일 |
| Boltz 원시 출력 | `data/boltz2/TREM2_{WT,R62H}_AB42_seed{1,2,3}/` — 60 파일 |
| AF3 원시 출력 | `data/alphafold3/` 6개 잡 — 66 파일 |

### 알려진 사항

- 판정 결과가 **어디에도 저장되지 않는다.** `det` `p` `mn` `ms` 전부 렌더 함수 지역 변수이고
  HTML 문자열로 즉시 소비된 뒤 폐기된다. 함수 밖으로 나가는 것은 `window._curGene` 하나뿐이다.
  종합 판정을 붙이려면 이 부분을 먼저 정리해야 한다
- 업로드 경로(`window._userResult`, `app/app.js:591`)는 유일하게 결과를 객체로 보존한다.
  `{baseline, signal, p, detected, nFit, usedCore, pos, label, acc, ...}` — 상태 설계 시 선례로 참고
- 업로드 경로에는 구조 축만 있고 결합 축이 없다. 종합 판정 적용 범위에서 일단 제외
- `app/data/*.csv` 3개는 `results/` 의 수동 사본이다. 파이프라인 재실행 시 동기화 필요
- 파이썬(scipy)과 웹앱(JS 자체 구현)이 같은 검정을 각각 구현하고 있다
- `app/app.js:52` 의 `const U = alternative==='greater' ? U1 : U1;` 은 양쪽 분기가 같은 무의미한 줄
  (기능에는 영향 없음, 미수정)

---

## 5. 적용 결과 (2026-09-02)

**리팩터링**
- `computeStructureVerdict(gene)` / `computeBindingVerdict()` 순수 함수 신설 (`app/app.js:244, 256`)
- `renderStructureAxis` / `renderBindingAxis` / `renderReport` 세 곳이 모두 이 함수만 호출
- `mannWhitneyU` 호출부 **4곳 -> 2곳** (두 compute 함수 내부에만 존재). 검정 로직 복사본 0개

**회귀 확인** — node vm 하니스로 각 화면이 대입하는 HTML 문자열 19개를 캡처해
`git show HEAD:app/app.js` (리팩터링 전) 결과와 바이트 단위 비교. **diff 비었음 = 회귀 없음.**
대상: s3 3개 사례 x {#sTitle #sScore #sCond #sPlot #sViewer} + s4 {#bScore #bPlot #bViewer} + s5 {#rTable}

**매트릭스** — `#rMatrix` 패널을 s5 판정 요약 위에 추가.
칸 배치는 실측 판정에서 유도하며 위치를 하드코딩하지 않는다 (`renderVerdictMatrix`, `app/app.js:484`).
결합 예측이 추가되거나 판정이 바뀌면 칸도 따라 움직인다.

## 6. 2단계 진행 (2026-09-03) — PILRA + NPDC1

**입력** — `src/make_boltz_npdc1_inputs.py` 로 YAML 4개 생성
(`PILRA_{G78,R78}_{NPDC1,alone}.yaml`). PILRA 32-150 + NPDC1 35-181 = 266 aa.
⚠️ UniProt 정본이 R78 이므로 정상형이 G78 이다 — 기존 `PILRA_dom_WT.fasta`(R78)와 라벨이 반대다.

**결과** — ipTM 조건당 10개(seed 1-2 x model 0-4). G78 중앙값 0.9088 / R78 0.8935 (차 -0.015).
**변동성 내 = 미감별.** 시드 1만 있던 n = 5 단계에서는 차가 +0.017 로 **부호가 반대**였다 —
표본을 늘리자 방향이 뒤집힌다는 것은 이 차이가 노이즈라는 증거다.

**스키마 확장** — `boltz_iptm.csv` 에 `gene` 칼럼 추가, 40행(TREM2 30 + PILRA 10).
`app/data/` 와 `results/` 두 사본 모두 갱신. PILRA 는 `iptm` 만 있고 나머지 4개 지표는 빈칸이다.

**`src/analyze_boltz.py` 병합화** — `write_csv()` 헬퍼 신설. 이 스크립트는 TREM2 만 수집하므로
단순 덮어쓰기를 하면 PILRA 행이 사라진다. 수집하지 않은 사례의 기존 행을 읽어 보존하도록 고쳤다.
(검증: 가짜 TREM2 30행으로 재실행 시 PILRA 10행이 값 그대로 유지됨)

**미감할 원인 추정** — 매트릭스 칸에 '추정' 으로 명시했다. PILRA 결합은 시알산 매개인데
Boltz-2 단백질-단백질 예측에는 글리칸이 없다. 모델의 감별 한계와 글리칸 부재를
이 결과만으로는 구분할 수 없다. 근거: Rathore N, et al. (2018) *PLoS Genet* 14(11):e1007427.

**범위 밖으로 둔 것** — s4(결합 인터페이스 검증) 화면은 TREM2 전용 유지.
PILRA 는 ipTM 수치만 있고 복합체 구조(.cif)가 없어 Mol* 뷰어에 띄울 대상이 없다.

## 7. 미결 사항

- CD33 의 결합 파트너 — 문헌 조사 선행 필요 (남은 미측정 사례 하나)
- PILRA 결합 축 표본은 n = 10 (시드 1-2). TREM2(n = 15)보다 적어 시드 3 추가 여지가 있다
- `#rTable` 구조 행은 `PILRA R78G`(CASES 라벨), 인터페이스 행은 `PILRA R78`(BIND_ALLELES 라벨)로
  같은 표 안에서 표기가 갈린다. 구조 행 라벨을 바꾸면 s3 출력도 바뀌므로 별도 판단 필요
- 업로드 경로(`window._userResult`)는 구조 축만 있어 매트릭스 적용 범위 밖
