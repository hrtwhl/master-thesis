"""
mvo.py
======
Transaction-cost-aware long-only mean-variance optimisation.

    max_w  mu' w  -  lam * w' Sigma w  -  tc * ||w - w_prev||_1
    s.t.   1' w = 1,   w >= 0,   ||w||_inf <= w_max

Exactly the formulation in Boukardagha (2026) §5.6 / §7.1.23 and the
notebook's `solve_mvo` function.  We try OSQP first, then SCS as a
fallback (also matching the notebook).
"""
import numpy as np
import cvxpy as cp

from config import LAM, TC, W_MAX


def solve_mvo(mu: np.ndarray,
              Sigma: np.ndarray,
              w_prev: np.ndarray,
              lam: float = LAM,
              tc: float = TC,
              w_max: float = W_MAX,
              solver: str = "OSQP") -> np.ndarray:
    mu = np.asarray(mu).reshape(-1)
    n  = mu.shape[0]

    Sigma = np.asarray(Sigma)
    Sigma = 0.5 * (Sigma + Sigma.T) + 1e-8 * np.eye(n)

    w = cp.Variable(n)
    objective = cp.Maximize(
        mu @ w
        - lam * cp.quad_form(w, cp.psd_wrap(Sigma))
        - tc  * cp.norm1(w - w_prev)
    )
    constraints = [
        cp.sum(w) == 1,
        w >= 0,
        cp.norm_inf(w) <= w_max,
    ]
    prob = cp.Problem(objective, constraints)

    try:
        prob.solve(solver=getattr(cp, solver), verbose=False)
        if w.value is None:
            prob.solve(solver=cp.SCS, verbose=False)
    except Exception:
        try:
            prob.solve(solver=cp.SCS, verbose=False)
        except Exception:
            return w_prev.copy()

    if w.value is None:
        return w_prev.copy()

    w_sol = np.asarray(w.value).reshape(-1)
    w_sol = np.maximum(w_sol, 0.0)
    s = w_sol.sum()
    if s <= 1e-12:
        return w_prev.copy()
    return w_sol / s
