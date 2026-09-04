"""감별력 검정 엔진 — 음성 대조군 쌍 설계 + 라벨 순열검정.

두 축이 같은 논리를 쓴다.

  노이즈 = 같은 조건 안의 쌍   (같은 서열을 반복 예측한 모델끼리 = 모델 진동)
  신호   = 조건 간의 쌍        (야생형 모델 x 변이형 모델)
  통계량 = median(신호) - median(노이즈)

  구조 축: 쌍 값이 RMSD 로 이미 주어진다 (거리 행렬을 그대로 쓴다)
  결합 축: ipTM 은 구조당 스칼라이므로 |ipTM_i - ipTM_j| 로 쌍을 만든다

왜 순열검정인가
  쌍은 서로 독립이 아니다. 같은 구조가 여러 쌍에 재사용되므로
  Mann-Whitney 의 p 는 과소평가된다 (예: PILRA 구조 축 6e-9 -> 0.0079).
  라벨을 뒤섞어 얻은 경험 분포로 교정한다.

  분할 수 <= PERM_EXACT_LIMIT 이면 전수 열거(정확), 아니면 고정 시드 표본.
  둘 다 결정적이므로 재실행해도 같은 값이 나온다.

⚠️ app/app.js 의 순열검정 엔진과 **같은 알고리즘·같은 PRNG(mulberry32)** 를 쓴다.
   두 구현의 p 는 자릿수까지 일치해야 한다 — 어긋나면 한쪽이 틀린 것이다.

⚠️ n = 5 + 5 (조건당 모델 5개) 에서 전수 분할은 126 개다.
   즉 도달 가능한 최소 p 는 1/126 = 0.0079 다. 그보다 작은 p 를 주장하려면
   모델 수를 늘려야 한다.
"""

from __future__ import annotations

import itertools
import math
from statistics import median

PERM_EXACT_LIMIT = 20000   # 분할 수가 이하이면 전수 열거
PERM_N = 20000             # 초과하면 표본 순열 횟수
PERM_SEED = 1

M32 = 0xFFFFFFFF


def mulberry32(seed: int):
    """JS 구현과 난수열이 동일한 결정적 PRNG.

    원본(JS):
        a = a + 0x6D2B79F5 | 0;
        let t = Math.imul(a ^ a>>>15, 1|a);
        t = t + Math.imul(t ^ t>>>7, 61|t) ^ t;
        return ((t ^ t>>>14) >>> 0) / 4294967296;
    """
    state = seed & M32

    def rnd() -> float:
        nonlocal state
        state = (state + 0x6D2B79F5) & M32
        a = state
        t = ((a ^ (a >> 15)) * (1 | a)) & M32
        t = ((t + (((t ^ (t >> 7)) * (61 | t)) & M32)) & M32) ^ t
        return ((t ^ (t >> 14)) & M32) / 4294967296

    return rnd


def matrix_from_scalars(values: list[float]) -> list[list[float]]:
    """스칼라 지표(ipTM 등) -> 쌍 차이 |x_i - x_j| 행렬."""
    n = len(values)
    return [[abs(values[i] - values[j]) for j in range(n)] for i in range(n)]


def matrix_from_rmsd_rows(rows, value_key: str = "rmsd_local8A"):
    """rmsd_pairs.csv 의 쌍 목록 -> 거리 행렬.

    0..k-1 = 야생형 모델, k..2k-1 = 변이형 모델.
    반환: (행렬, k). 쌍이 하나라도 비면 AssertionError 를 낸다.
    app/app.js 의 computeStructureVerdict 와 같은 인덱싱이다.
    """
    k = max(int(r["model_b"]) for r in rows) + 1
    N = 2 * k
    D = [[0.0] * N for _ in range(N)]
    filled = 0
    for r in rows:
        i, j, v = int(r["model_a"]), int(r["model_b"]), float(r[value_key])
        if r["comparison"] == "noise_WT":
            a, b = i, j
        elif r["comparison"] == "noise_MUT":
            a, b = k + i, k + j
        else:
            a, b = i, k + j
        if D[a][b] == 0.0:
            filled += 1
        D[a][b] = D[b][a] = v
    assert filled == N * (N - 1) // 2, f"쌍 행렬이 불완전하다 ({filled}/{N*(N-1)//2})"
    return D, k


def split_pairs(D, N: int, in_a) -> tuple[list[float], list[float]]:
    """라벨에 따라 쌍을 노이즈(같은 조건)/신호(조건 간)로 가른다."""
    noise, signal = [], []
    for i in range(N):
        for j in range(i + 1, N):
            (noise if in_a[i] == in_a[j] else signal).append(D[i][j])
    return noise, signal


def pair_stat(D, N: int, in_a) -> float:
    noise, signal = split_pairs(D, N, in_a)
    return median(signal) - median(noise)


def permutation_p(D, n1: int, n2: int,
                  exact_limit: int = PERM_EXACT_LIMIT,
                  nperm: int = PERM_N,
                  seed: int = PERM_SEED) -> dict:
    """라벨 순열검정. -> {'p', 'method', 'n_perm'}"""
    N = n1 + n2
    in_a = [1] * n1 + [0] * n2
    obs = pair_stat(D, N, in_a)

    if math.comb(N, n1) <= exact_limit:
        # 통계량이 A/B 라벨 교환에 불변이므로 0번을 A 에 고정해 여집합 대칭을 없앤다
        cnt = tot = 0
        for pick in itertools.combinations(range(N), n1):
            if 0 not in pick:
                continue
            lab = [0] * N
            for i in pick:
                lab[i] = 1
            tot += 1
            if pair_stat(D, N, lab) >= obs:
                cnt += 1
        return {"p": cnt / tot, "method": "exact", "n_perm": tot}

    rnd = mulberry32(seed)
    ord_ = list(range(N))
    cnt = 0
    for _ in range(nperm):
        for i in range(N - 1, 0, -1):                 # JS 와 동일한 Fisher-Yates
            j = int(rnd() * (i + 1))
            ord_[i], ord_[j] = ord_[j], ord_[i]
        lab = [0] * N
        for i in range(n1):
            lab[ord_[i]] = 1
        if pair_stat(D, N, lab) >= obs:
            cnt += 1
    return {"p": (cnt + 1) / (nperm + 1), "method": "sample", "n_perm": nperm}


def verdict_from_scalars(a: list[float], b: list[float], alpha: float = 0.05) -> dict:
    """스칼라 두 그룹(예: ipTM) -> 감별력 판정."""
    D = matrix_from_scalars(list(a) + list(b))
    N = len(a) + len(b)
    in_a = [1] * len(a) + [0] * len(b)
    noise, signal = split_pairs(D, N, in_a)
    perm = permutation_p(D, len(a), len(b))
    mn, ms = median(noise), median(signal)
    return {"noise": noise, "signal": signal, "mn": mn, "ms": ms,
            "ratio": ms / mn if mn else float("inf"),
            "p": perm["p"], "method": perm["method"], "n_perm": perm["n_perm"],
            "det": perm["p"] < alpha}


def verdict_from_matrix(D, n1: int, n2: int, alpha: float = 0.05) -> dict:
    """이미 쌍 값(RMSD 등)이 있는 경우 -> 감별력 판정."""
    N = n1 + n2
    in_a = [1] * n1 + [0] * n2
    noise, signal = split_pairs(D, N, in_a)
    perm = permutation_p(D, n1, n2)
    mn, ms = median(noise), median(signal)
    return {"noise": noise, "signal": signal, "mn": mn, "ms": ms,
            "ratio": ms / mn if mn else float("inf"),
            "p": perm["p"], "method": perm["method"], "n_perm": perm["n_perm"],
            "det": perm["p"] < alpha}
