"""플래시 라이다 레이캐스팅 (Open3D RaycastingScene, Embree)

센서 모델 v0
- 해상도 res x res, FOV fov_deg x fov_deg 핀홀 격자
- 거리 노이즈: 레이 방향 가우시안 sigma [m]
- 입사각 > grazing_max_deg 인 점 제거, 나머지에서 무작위 dropout 비율만큼 제거
- 제외(v1 이후): 가장자리 mixed pixel, 재질 반사율
"""
from dataclasses import dataclass, asdict

import numpy as np
import open3d as o3d


@dataclass
class LidarParams:
    res: int = 256
    fov_deg: float = 40.0
    noise_sigma: float = 0.0        # [m] 레이 방향 거리 노이즈
    grazing_max_deg: float = 80.0
    dropout: float = 0.01
    max_range: float = np.inf

    def to_dict(self):
        d = asdict(self)
        d["max_range"] = float(d["max_range"]) if np.isfinite(d["max_range"]) else -1.0
        return d


class TargetScene:
    """부품별 메시를 geometry 하나씩 등록 → 레이가 맞은 geometry id = 부품 인덱스"""

    def __init__(self, parts):
        self.scene = o3d.t.geometry.RaycastingScene()
        self.names = [n for n, _ in parts]
        self.gid2label = {}
        for k, (_, m) in enumerate(parts):
            gid = self.scene.add_triangles(
                o3d.core.Tensor(np.asarray(m.vertices, np.float32)),
                o3d.core.Tensor(np.asarray(m.faces, np.uint32)))
            self.gid2label[int(gid)] = k + 1               # 1-based 부품 라벨

    def cast(self, origin, target, lp: LidarParams, rng, up=(0, 0, 1)):
        """센서 위치 origin 에서 target 점을 바라보고 한 프레임 촬영

        반환 dict: P(노이즈), P_true, N(센서 쪽 법선), label, range, incidence[deg]
        """
        o = np.asarray(origin, float)
        D = look_at_rays(o, np.asarray(target, float), lp.res, lp.fov_deg, up)
        rays = np.hstack([np.broadcast_to(o, D.shape), D]).astype(np.float32)
        ans = self.scene.cast_rays(o3d.core.Tensor(rays))

        t = ans["t_hit"].numpy().astype(float)
        gid = ans["geometry_ids"].numpy().astype(np.int64)
        n = ans["primitive_normals"].numpy().astype(float)
        hit = np.isfinite(t) & (t <= lp.max_range)
        t, gid, n, D = t[hit], gid[hit], n[hit], D[hit]

        n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
        flip = np.einsum("ij,ij->i", n, D) > 0              # 뒷면 → 센서 쪽으로
        n[flip] *= -1
        inc = np.degrees(np.arccos(np.clip(-np.einsum("ij,ij->i", n, D), -1, 1)))

        keep = (inc <= lp.grazing_max_deg) & (rng.random(t.size) >= lp.dropout)
        t, gid, n, D, inc = t[keep], gid[keep], n[keep], D[keep], inc[keep]

        t_noisy = t + lp.noise_sigma * rng.standard_normal(t.size)
        label = np.array([self.gid2label[g] for g in gid], np.uint32)
        return dict(P=o + t_noisy[:, None] * D, P_true=o + t[:, None] * D, N=n,
                    label=label, range=t, incidence=inc)


def look_at_rays(origin, target, res, fov_deg, up=(0, 0, 1)):
    """핀홀 격자 레이 방향 (res^2 x 3, 단위벡터)"""
    f = target - origin
    f /= np.linalg.norm(f)
    up = np.asarray(up, float)
    if abs(f @ up) > 0.99:                                 # 시선과 up 이 거의 평행
        up = np.array([1.0, 0, 0]) if abs(f[0]) < 0.9 else np.array([0, 1.0, 0])
    r = np.cross(f, up); r /= np.linalg.norm(r)
    u = np.cross(r, f)
    h = np.tan(np.deg2rad(fov_deg) / 2)
    s = (np.arange(res) + 0.5) / res * 2 - 1               # 픽셀 중심 [-1, 1]
    X, Y = np.meshgrid(s * h, -s * h)
    D = f + X.reshape(-1, 1) * r + Y.reshape(-1, 1) * u
    return D / np.linalg.norm(D, axis=1, keepdims=True)


def fibonacci_sphere(n):
    """구면 위 거의 균일한 n개 방향"""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])


def voxel_downsample(cloud, h, rng):
    """voxel 당 대표점 1개 (무작위 선택). 평균을 내지 않으므로 라벨/법선이 섞이지 않는다"""
    n = cloud["P"].shape[0]
    perm = rng.permutation(n)
    key = np.floor(cloud["P"][perm] / h).astype(np.int64)
    _, first = np.unique(key, axis=0, return_index=True)
    idx = np.sort(perm[first])
    return {k: v[idx] for k, v in cloud.items()}


def merge(clouds):
    return {k: np.concatenate([c[k] for c in clouds]) for k in clouds[0]}
