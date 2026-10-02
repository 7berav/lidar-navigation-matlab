"""타깃 메시 레이캐스팅 → 커팅 입력 점군 (.mat)

모드
  complete : 피보나치 구면 시점 합집합 + voxel 균일화  (가림 없는 정합 데이터, L1)
  s1       : CW 접근 hold point 누적 점군               (시나리오 S1, L2)
             타깃 중심 LVLH (x: R-bar 반경, y: V-bar 진행, z: H-bar 궤도법선), 자세 고정 (--att)

거리/노이즈/voxel 은 타깃 경계구 반경 R 대비 상대값으로 지정 → 타깃 크기와 무관하게 같은 규칙

출력: out_cutting_sim/<yyyyMMdd_HHmmss>_<target>/
  complete.mat, s1_<dir>_h<k>.mat   (k번째 hold point 까지 누적)
  preview_complete.png, preview_s1.png, summary.json
  .mat 필드: P P_true N label range incidence view_id sensor_pos part_names meta_json
  (label 1부터 = part_names 순서, P 는 타깃 중심 LVLH 좌표 [m], body 변환은 meta_json)

사용 예
  .venv/Scripts/python sim/raycast_target.py sim/meshes/cubesat6u.glb
"""
import argparse
import json
import os
from datetime import datetime

import numpy as np
import scipy.io as sio
from scipy.spatial import cKDTree

from lidar import LidarParams, TargetScene, fibonacci_sphere, merge, voxel_downsample
from targets import bounding_sphere, load_target

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ATTITUDES = {                        # R_b2l: 열 = body 축을 LVLH 로 표현
    "identity": np.eye(3),
    # 날개 윗면(body +y) → 천정(+R-bar), 버스 장축(body z) → V-bar, 날개 스팬(body x) → H-bar
    "zenith": np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float),
}

S1_DIRS = {                          # 추적선이 있는 쪽 (LVLH), 타깃을 향해 접근
    "+vbar":   (0, 1, 0),
    "-vbar":   (0, -1, 0),
    "-rbar":   (-1, 0, 0),           # 지구 쪽(아래)에서 접근
    "oblique": (-0.4, 0.8, 0.45),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("glb")
    ap.add_argument("--modes", nargs="+", default=["complete", "s1"], choices=["complete", "s1"])
    ap.add_argument("--res", type=int, default=256)
    ap.add_argument("--fov", type=float, default=40.0)
    ap.add_argument("--noise-rel", type=float, default=0.002, help="거리 노이즈 sigma / R")
    ap.add_argument("--voxel-rel", type=float, default=0.01, help="voxel 크기 / R (target-n=0 일 때)")
    ap.add_argument("--target-n", type=int, default=20000,
                    help="완전 데이터 목표 점 수. voxel 크기를 여기에 맞춰 정하고 S1 에도 같은 값을 쓴다 (0=끔)")
    ap.add_argument("--complete-views", type=int, default=150)
    ap.add_argument("--complete-dist", type=float, default=4.0, help="완전 데이터 촬영 거리 / R")
    ap.add_argument("--s1-dirs", nargs="+", default=["+vbar", "-rbar", "oblique"], choices=list(S1_DIRS))
    ap.add_argument("--s1-holds", nargs="+", type=float, default=[10, 5, 3], help="hold point 거리 / R")
    ap.add_argument("--att", default="zenith", choices=list(ATTITUDES),
                    help="타깃 자세 (body→LVLH). 점군은 타깃 중심 LVLH 좌표로 저장")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-root", default=os.path.join(REPO, "out_cutting_sim"))
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    target = os.path.splitext(os.path.basename(a.glb))[0]
    parts = load_target(a.glb)
    names = [n for n, _ in parts]
    center_body, R = bounding_sphere(np.vstack([m.vertices for _, m in parts]))
    # body → 타깃 중심 LVLH:  p_l = R_b2l (p_b - center_body)
    Rb2l = ATTITUDES[a.att]
    T = np.eye(4); T[:3, :3] = Rb2l; T[:3, 3] = -Rb2l @ center_body
    parts = [(n, m.copy().apply_transform(T)) for n, m in parts]
    center = np.zeros(3)
    h = a.voxel_rel * R
    lp = LidarParams(res=a.res, fov_deg=a.fov, noise_sigma=a.noise_rel * R)
    ts = TargetScene(parts)

    out_dir = os.path.join(a.out_root, datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + target)
    os.makedirs(out_dir, exist_ok=True)
    base_meta = dict(target=target, glb=os.path.abspath(a.glb), R=float(R),
                     attitude=a.att, R_b2l=Rb2l.tolist(), center_body=center_body.tolist(),
                     voxel=float(h), lidar=lp.to_dict(), seed=a.seed,
                     frame="target-centered LVLH (x R-bar, y V-bar, z H-bar); "
                           "p_body = R_b2l' * p + center_body")
    print(f"타깃 {target}: 부품 {len(names)}개 {names}, R={R:.4f} m, voxel h={h * 1000:.2f} mm, "
          f"noise sigma={lp.noise_sigma * 1000:.2f} mm")
    summary = dict(meta=base_meta, clouds={})

    # ---- complete: 구면 다시점 합집합
    ref_keys = None
    if "complete" in a.modes:
        views = center + a.complete_dist * R * fibonacci_sphere(a.complete_views)
        frames = []
        for v, c in enumerate(views):
            fr = ts.cast(c, center, lp, rng)
            fr["view_id"] = np.full(fr["label"].size, v + 1, np.uint32)
            frames.append(fr)
        raw = merge(frames)
        if a.target_n > 0:                       # 표면이므로 N ∝ h^-2 → 두 번 보정하면 충분
            for _ in range(2):
                n_now = voxel_downsample(raw, h, np.random.default_rng(0))["label"].size
                h *= (n_now / a.target_n) ** 0.5
            base_meta["voxel"] = float(h)
            print(f"  voxel 자동 조정: h={h * 1000:.2f} mm (h/R={h / R:.4f})")
        cloud = voxel_downsample(raw, h, rng)
        small = [f"{nm}({c})" for nm, c in zip(names, np.bincount(cloud["label"], minlength=len(names) + 1)[1:]) if c < 200]
        if small:
            print("  [주의] 200점 미만 부품:", ", ".join(small))
        ref_keys = {"all": voxel_keys(cloud["P_true"], h)}
        ref_keys.update({k: voxel_keys(cloud["P_true"][cloud["label"] == k], h)
                         for k in range(1, len(names) + 1)})
        meta = dict(base_meta, mode="complete", n_views=a.complete_views,
                    dist_R=a.complete_dist, n_raw=int(raw["label"].size))
        s = save_cloud(out_dir, "complete", cloud, views, names, meta, ref_keys, h)
        summary["clouds"]["complete"] = s
        preview_complete(cloud, views, names, os.path.join(out_dir, "preview_complete.png"))

    # ---- s1: hold point 누적
    if "s1" in a.modes:
        grid = {}
        for d in a.s1_dirs:
            u = np.asarray(S1_DIRS[d], float); u /= np.linalg.norm(u)
            frames, sensors = [], []
            for k, dist in enumerate(a.s1_holds):
                c = center + dist * R * u
                fr = ts.cast(c, center, lp, rng)
                fr["view_id"] = np.full(fr["label"].size, k + 1, np.uint32)
                frames.append(fr); sensors.append(c)
                cloud = voxel_downsample(merge(frames), h, rng)
                name = f"s1_{d.replace('+', 'p').replace('-', 'm')}_h{k + 1}"
                meta = dict(base_meta, mode="s1", direction=d, dir_vec=u.tolist(),
                            holds_R=a.s1_holds[:k + 1])
                s = save_cloud(out_dir, name, cloud, np.array(sensors), names, meta, ref_keys, h)
                summary["clouds"][name] = s
                grid[(d, k)] = cloud
        preview_s1(grid, a.s1_dirs, a.s1_holds, names, os.path.join(out_dir, "preview_s1.png"))

    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"\n저장: {out_dir}")


# ---------------------------------------------------------------- 저장/진단
def voxel_keys(P, h):
    return set(map(tuple, np.floor(P / h).astype(np.int64)))


def save_cloud(out_dir, name, cloud, sensors, names, meta, ref_keys, h):
    n = cloud["label"].size
    cnt = np.bincount(cloud["label"], minlength=len(names) + 1)[1:]
    # 밀도 불균일: 10-NN 거리 분위수 비 → 밀도비 ~ (d95/d5)^2
    d10 = cKDTree(cloud["P"]).query(cloud["P"], k=11)[0][:, 10]
    q5, q50, q95 = np.percentile(d10, [5, 50, 95])
    s = dict(n=int(n), part_counts=dict(zip(names, cnt.tolist())),
             density_ratio_5_95=float((q95 / q5) ** 2), knn10_median=float(q50))
    if ref_keys is not None and meta["mode"] != "complete":
        # 가시율: 완전 데이터가 차지한 voxel 중 이 점군이 덮은 비율 (전체 / 부품별)
        s["coverage"] = len(voxel_keys(cloud["P_true"], h) & ref_keys["all"]) / len(ref_keys["all"])
        s["part_coverage"] = {
            nm: len(voxel_keys(cloud["P_true"][cloud["label"] == k], h) & ref_keys[k]) / max(len(ref_keys[k]), 1)
            for k, nm in enumerate(names, start=1)}
    meta = dict(meta, stats=s)
    sio.savemat(os.path.join(out_dir, name + ".mat"), dict(
        P=cloud["P"], P_true=cloud["P_true"], N=cloud["N"], label=cloud["label"],
        range=cloud["range"], incidence=cloud["incidence"], view_id=cloud["view_id"],
        sensor_pos=np.asarray(sensors, float), part_names=np.array(names, dtype=object),
        meta_json=json.dumps(meta, ensure_ascii=False)), do_compression=True)
    cov = (f"  coverage={s['coverage']:.2f} parts="
           f"{[round(v, 2) for v in s['part_coverage'].values()]}") if "coverage" in s else ""
    print(f"  {name:<18s} N={n:<6d} parts={cnt.tolist()}  density_ratio(5-95%)={s['density_ratio_5_95']:.1f}x{cov}")
    return s


def _scatter(ax, P, lab, names, nmax=15000, s=0.6):
    import matplotlib.pyplot as plt
    idx = np.random.default_rng(0).permutation(len(P))[:nmax]
    cm = plt.get_cmap("tab10")
    for k in range(1, len(names) + 1):
        m = idx[lab[idx] == k]
        ax.scatter(P[m, 0], P[m, 1], P[m, 2], s=s, color=cm((k - 1) % 10), label=names[k - 1])
    lo, hi = P.min(0), P.max(0)
    ax.set_box_aspect(np.maximum(hi - lo, 1e-3 * (hi - lo).max()))
    ax.set_xlabel("x (R-bar)"); ax.set_ylabel("y (V-bar)"); ax.set_zlabel("z (H-bar)")


def preview_complete(cloud, views, names, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 5))
    for i, (el, az) in enumerate([(25, -60), (90, -90), (0, 0)]):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        _scatter(ax, cloud["P"], cloud["label"], names)
        ax.view_init(el, az)
        ax.set_title(f"complete N={cloud['label'].size} (el={el}, az={az})")
    ax.legend(markerscale=10, loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def preview_s1(grid, dirs, holds, names, path):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(4.2 * len(holds), 3.8 * len(dirs)))
    for i, d in enumerate(dirs):
        for k, dist in enumerate(holds):
            c = grid[(d, k)]
            ax = fig.add_subplot(len(dirs), len(holds), i * len(holds) + k + 1, projection="3d")
            _scatter(ax, c["P"], c["label"], names, s=1.0)
            ax.view_init(25, -60)
            ax.set_title(f"{d}  ~{dist:g}R  N={c['label'].size}", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
