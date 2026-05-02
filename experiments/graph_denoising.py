"""
Graph spectral denoising on a pixel grid with similarity weights.

Mirrors src/denoise_experiment.m so results can be reproduced without MATLAB.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def build_pixel_graph(img: np.ndarray, sigma_r: float) -> np.ndarray:
    """4-neighbor adjacency with weights exp(-(I_i-I_j)^2 / (2 sigma_r^2))."""
    h, w = img.shape
    n = h * w
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    inv_den = 1.0 / (2.0 * sigma_r**2)

    for r in range(h):
        for c in range(w):
            i = r * w + c
            vi = img[r, c]
            if c + 1 < w:
                j = i + 1
                vj = img[r, c + 1]
                ww = float(np.exp(-((vi - vj) ** 2) * inv_den))
                rows.extend([i, j])
                cols.extend([j, i])
                vals.extend([ww, ww])
            if r + 1 < h:
                j = i + w
                vj = img[r + 1, c]
                ww = float(np.exp(-((vi - vj) ** 2) * inv_den))
                rows.extend([i, j])
                cols.extend([j, i])
                vals.extend([ww, ww])

    a = np.zeros((n, n), dtype=float)
    for ri, ci, v in zip(rows, cols, vals, strict=True):
        a[ri, ci] += v
    return a


def combinatorial_laplacian(a: np.ndarray) -> np.ndarray:
    d = np.sum(a, axis=1)
    return np.diag(d) - a


def ideal_low_pass(lam: np.ndarray, cutoff: float) -> np.ndarray:
    return (lam <= cutoff).astype(float)


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def checkerboard(h: int, w: int, block: int = 12) -> np.ndarray:
    rr, cc = np.meshgrid(np.arange(w), np.arange(h))
    return ((np.floor(rr / block) + np.floor(cc / block)) % 2 == 0).astype(float)


def run_experiment(
    h: int = 48,
    w: int = 48,
    block: int = 12,
    sigma_n: float = 0.25,
    sigma_r: float = 0.2,
    cutoff_frac: float = 0.15,
    seed: int = 42,
) -> dict:
    rng = np.random.default_rng(seed)
    clean = checkerboard(h, w, block)
    noisy = clean + sigma_n * rng.standard_normal((h, w))

    a = build_pixel_graph(noisy, sigma_r)
    lmat = combinatorial_laplacian(a)
    lam, v = np.linalg.eigh(lmat)
    lam = np.real(lam)
    lam_max = float(np.max(lam))
    cutoff = cutoff_frac * lam_max
    h_resp = ideal_low_pass(lam, cutoff)

    x = noisy.reshape(-1)
    x_clean = clean.reshape(-1)
    spec = v.T @ x
    y = v @ (h_resp * spec)

    return {
        "clean": clean,
        "noisy": noisy,
        "filtered": y.reshape(h, w),
        "rmse_noisy": rmse(x, x_clean),
        "rmse_filtered": rmse(y, x_clean),
        "lam_max": lam_max,
        "cutoff": cutoff,
        "n": h * w,
    }


def save_figure_triptych(out_path: Path, result: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    for ax, img, title in zip(
        axes,
        [result["clean"], result["noisy"], result["filtered"]],
        ["Clean", "Noisy", "Filtered"],
        strict=True,
    ):
        im = ax.imshow(img, cmap="gray", vmin=0, vmax=1)
        ax.set_title(title)
        ax.axis("off")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(
        f"RMSE noisy={result['rmse_noisy']:.4f}, filtered={result['rmse_filtered']:.4f}"
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("figures/denoising_triptych.png"))
    args = p.parse_args()

    r = run_experiment()
    print(f"RMSE noisy: {r['rmse_noisy']:.4f}")
    print(f"RMSE filtered (ideal LPF): {r['rmse_filtered']:.4f}")
    save_figure_triptych(args.out, r)
    print(f"Wrote {args.out.resolve()}")


if __name__ == "__main__":
    main()
