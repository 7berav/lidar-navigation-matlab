"""절차적 타깃 메시 생성 → sim/meshes/<이름>.glb + <이름>_parts.json

사용 예
  .venv/Scripts/python sim/make_target.py all                 # 전체 생성 + 한 장짜리 미리보기
  .venv/Scripts/python sim/make_target.py cubesat6u
  .venv/Scripts/python sim/make_target.py cubesat6u --out sim/meshes/cubesat6u_a30.glb --param wing_angle_deg=30
"""
import argparse
import json
import os

import numpy as np

from targets import BUILDERS, save_target

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name", choices=sorted(BUILDERS) + ["all"])
    ap.add_argument("--out", default=None, help="출력 GLB 경로 (기본 sim/meshes/<name>.glb)")
    ap.add_argument("--param", action="append", default=[],
                    help="빌더 파라미터 덮어쓰기 key=value (value 는 JSON, 예: bus=[0.2,0.1,0.3])")
    a = ap.parse_args()

    kw = {}
    for kv in a.param:
        k, v = kv.split("=", 1)
        kw[k] = json.loads(v)

    names = list(BUILDERS) if a.name == "all" else [a.name]
    built = []
    for name in names:
        parts, params = BUILDERS[name](**kw)
        out = (a.out if a.out and a.name != "all" else os.path.join(HERE, "meshes", name + ".glb"))
        info = save_target(parts, out, name, params)
        built.append((name, parts, info))
        print(f"{name:<16s} K={len(parts):<2d} R={info['R']:7.3f} m  "
              + ", ".join(f"{p['name']}({p['area']:.2f})" for p in info["parts"]))

    if a.name == "all":
        path = os.path.join(HERE, "meshes", "targets_overview.png")
        overview(built, path)
        print(f"미리보기: {path}")


def overview(built, path):
    """전체 타깃을 부품별 색으로 한 장에 그림"""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    cm = plt.get_cmap("tab10")
    light = np.array([0.5, -0.6, 0.62]); light /= np.linalg.norm(light)
    nc = 4; nr = int(np.ceil(len(built) / nc))
    fig = plt.figure(figsize=(4.6 * nc, 4.2 * nr))
    for i, (name, parts, info) in enumerate(built):
        ax = fig.add_subplot(nr, nc, i + 1, projection="3d")
        for k, (_, m) in enumerate(parts):
            lit = 0.45 + 0.55 * np.abs(m.face_normals @ light)        # 면 법선으로 직접 음영
            fc = np.column_stack([np.outer(lit, cm(k % 10)[:3]), np.ones(len(lit))])
            ax.add_collection3d(Poly3DCollection(m.triangles, facecolors=fc, linewidths=0))
        V = np.vstack([m.vertices for _, m in parts])
        lo, hi = V.min(0), V.max(0)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(np.maximum(hi - lo, 0.02 * (hi - lo).max()))
        ax.view_init(24, -58); ax.set_axis_off()
        ax.set_title(f"{name}  (K={len(parts)}, R={info['R']:.2f} m)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
