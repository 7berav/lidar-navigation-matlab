"""타깃 메시 생성/로드

- 부품 하나 = GLB 노드 하나. 노드 이름이 부품 이름, 부품 ID는 1부터 (0 = 배경/미할당, MATLAB 규약)
- 단위 m, 타깃 body 좌표계. 공통 규약: 주축 z, 날개 스팬 ±x, 태양전지판 법선 ±y
  (raycast_target.py 의 'zenith' 자세에서 y→천정, z→V-bar, x→H-bar)
- build_* 함수는 부품 메시 목록 [(name, trimesh.Trimesh)] 과 파라미터 dict 를 돌려준다
- 가는 붐/요크/마스트는 점이 거의 안 찍히므로 독립 부품으로 두지 않고 매달린 부품에 합친다
"""
import json
import os

import numpy as np
import trimesh


# ---------------------------------------------------------------- 기본 도형
def _place(m, center=(0, 0, 0), axis=(0, 0, 1)):
    """로컬 z축을 axis 로 돌린 뒤 center 로 이동"""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    if not np.allclose(axis, (0, 0, 1)):
        m.apply_transform(trimesh.geometry.align_vectors((0, 0, 1), axis))
    m.apply_translation(center)
    return m


def _box(ext, center=(0, 0, 0)):
    return _place(trimesh.creation.box(extents=ext), center)


def _cyl(r, h, center=(0, 0, 0), axis=(0, 0, 1), sections=48, r_top=None):
    """원기둥 (sections=6 이면 육각기둥). r_top 을 주면 원뿔대 (로컬 +z 쪽 반지름)"""
    m = trimesh.creation.cylinder(radius=r, height=h, sections=sections)
    if r_top is not None:
        v = m.vertices.copy()
        v[v[:, 2] > 0, :2] *= max(r_top, 1e-4 * r) / r
        m.vertices = v
    return _place(m, center, axis)


def _sphere(r, center=(0, 0, 0), scale=(1, 1, 1)):
    m = trimesh.creation.icosphere(subdivisions=3, radius=r)
    m.apply_scale(scale)
    return _place(m, center)


def _dish(R, depth, thick, center=(0, 0, 0), axis=(0, 0, 1), n=16, sections=48):
    """포물면 접시 (두께 있는 껍질). 꼭짓점이 center, 열린 쪽이 axis 방향"""
    r = np.linspace(0, R, n)
    prof = np.vstack([np.column_stack([r, depth * (r / R) ** 2]),
                      np.column_stack([r[::-1], depth * (r[::-1] / R) ** 2 + thick])])
    m = trimesh.creation.revolve(prof, sections=sections)
    if m.volume < 0:
        m.invert()
    return _place(m, center, axis)


def _panel(root, span_dir, span, width, thick=0.02, normal=(0, 1, 0)):
    """판: root 에서 span_dir 로 span 만큼 뻗고, 법선 normal, 폭은 나머지 축"""
    x = np.asarray(span_dir, float); x /= np.linalg.norm(x)
    y = np.asarray(normal, float); y -= (y @ x) * x; y /= np.linalg.norm(y)
    T = np.eye(4)
    T[:3, :3] = np.column_stack([x, y, np.cross(x, y)])
    T[:3, 3] = np.asarray(root, float) + x * span / 2
    m = trimesh.creation.box(extents=(span, thick, width))
    m.apply_transform(T)
    return m


def _boom(p0, p1, r=0.03):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    return _cyl(r, np.linalg.norm(p1 - p0), (p0 + p1) / 2, p1 - p0, sections=12)


def _merge(*ms):
    return trimesh.util.concatenate(list(ms))


# ---------------------------------------------------------------- A. 단순 (K=1~5)
def build_cubesat6u(bus=(0.2263, 0.100, 0.366),
                    wing_len=0.44, wing_width=0.34, wing_thick=0.0025,
                    hinge_gap=0.002, wing_angle_deg=0.0, wing_y=1.0):
    """6U 큐브위성 + ±x 방향 이중 전개 태양전지판 2장 (K=3)

    body 좌표: x = 2U 폭, y = 1U 두께, z = 3U 길이, 원점 = 버스 중심
    날개: ±x 모서리에 힌지, 힌지 축은 z. wing_y = 힌지 높이 (1 = +y 면과 같은 높이, 0 = 옆면 중간)
          wing_angle_deg = 0 이면 xz 평면과 평행, 양수면 +y 쪽으로 접혀 올라감
    wing_y=1, angle=0 (기본)은 날개가 버스 윗면과 동일 평면·같은 폭 → 병목이 없는 어려운 형상
    """
    bx, by, bz = bus
    parts = [("bus", trimesh.creation.box(extents=bus))]

    for sign, name in ((+1, "wing_px"), (-1, "wing_nx")):
        # 힌지 원점 로컬 좌표: x 가 바깥쪽, 판 윗면이 y=0
        w = trimesh.creation.box(extents=(wing_len, wing_thick, wing_width))
        w.apply_translation((hinge_gap + wing_len / 2, -wing_thick / 2, 0.0))
        th = np.deg2rad(wing_angle_deg)
        w.apply_transform(trimesh.transformations.rotation_matrix(th, (0, 0, 1)))
        if sign < 0:                                     # x 거울상 (법선 방향 유지 위해 면 순서 뒤집기)
            v = w.vertices.copy(); v[:, 0] *= -1; w.vertices = v
            w.invert()
        w.apply_translation((sign * bx / 2, wing_y * by / 2, 0.0))
        parts.append((name, w))

    params = dict(bus=list(bus), wing_len=wing_len, wing_width=wing_width, wing_thick=wing_thick,
                  hinge_gap=hinge_gap, wing_angle_deg=wing_angle_deg, wing_y=wing_y)
    return parts, params


def build_cubesat6u_mid():
    """6U, 날개를 옆면 중간 높이에 장착 (T자 접합 → 법선이 달라지는 쉬운 형상)"""
    return build_cubesat6u(wing_y=0.0, hinge_gap=0.004)


def build_cubesat3u_petal(bus=(0.10, 0.10, 0.34), petal_len=0.34, petal_width=0.083, thick=0.002):
    """3U 버스 + 한쪽 끝에서 꽃잎처럼 펼친 패널 4장 (K=5)"""
    bx, by, bz = bus
    z0 = -bz / 2 + thick
    parts = [("bus", _box(bus))]
    for name, d in (("petal_px", (1, 0, 0)), ("petal_nx", (-1, 0, 0)),
                    ("petal_py", (0, 1, 0)), ("petal_ny", (0, -1, 0))):
        root = np.array(d, float) * (bx / 2 + 0.002) + (0, 0, z0)
        parts.append((name, _panel(root, d, petal_len, petal_width, thick, normal=(0, 0, 1))))
    return parts, dict(bus=list(bus), petal_len=petal_len, petal_width=petal_width)


def build_hexsat(r=0.6, h=1.2, wing_span=2.4, wing_width=1.0, dish_r=0.45):
    """육각기둥 버스 + 요크 달린 날개 2장 + 상단 접시 안테나 (K=4)"""
    parts = [("bus", _cyl(r, h, sections=6))]
    for s, name in ((1, "wing_px"), (-1, "wing_nx")):
        parts.append((name, _merge(
            _boom((s * r * 0.9, 0, 0), (s * (r + 0.35), 0, 0), 0.025),
            _panel((s * (r + 0.35), 0, 0), (s, 0, 0), wing_span, wing_width, 0.02))))
    parts.append(("dish", _merge(_boom((0, 0, h / 2 - 0.02), (0, 0, h / 2 + 0.27), 0.04),
                                 _dish(dish_r, 0.12, 0.015, (0, 0, h / 2 + 0.25)))))
    return parts, dict(r=r, h=h, wing_span=wing_span, wing_width=wing_width, dish_r=dish_r)


def build_spin_drum(r=1.08, h=2.2, dish_r=0.6):
    """스핀 안정 드럼형 위성 + 상단 디스펀 안테나 (K=2)"""
    parts = [("drum", _cyl(r, h)),
             ("antenna", _merge(_boom((0, 0, h / 2 - 0.02), (0, 0, h / 2 + 0.62), 0.06),
                                _dish(dish_r, 0.18, 0.02, (0, 0.05, h / 2 + 0.75), axis=(0, 1, 0.25))))]
    return parts, dict(r=r, h=h, dish_r=dish_r)


# ---------------------------------------------------------------- B. 중형 (K=2~7)
def build_geo_comsat(bus=(2.0, 2.2, 3.0), wing_span=8.0, wing_width=2.2, yoke=1.2, dish_r=1.1):
    """GEO 통신위성: 박스 버스 + 긴 날개 2장 + 측면 반사판 2개 (K=5)"""
    bx, by, bz = bus
    parts = [("bus", _box(bus))]
    for s, name in ((1, "wing_px"), (-1, "wing_nx")):
        x0 = s * (bx / 2 + yoke)
        parts.append((name, _merge(_boom((s * bx / 2 * 0.95, 0, 0), (x0, 0, 0), 0.04),
                                   _panel((x0, 0, 0), (s, 0, 0), wing_span, wing_width, 0.03))))
    for s, name in ((1, "reflector_py"), (-1, "reflector_ny")):
        c = np.array((0, s * (by / 2 + 1.0), 0.9))
        parts.append((name, _merge(_boom((0, s * by / 2 * 0.95, 0.2), c, 0.04),
                                   _dish(dish_r, 0.25, 0.03, c, axis=(0, s * 0.35, 0.94)))))
    return parts, dict(bus=list(bus), wing_span=wing_span, wing_width=wing_width, dish_r=dish_r)


def build_eo_sat(bus=(1.8, 1.8, 4.0), wing_span=7.0, wing_width=2.4):
    """지구관측 위성 (Landsat 유형): 박스 버스 + 한쪽 대형 패널 + 망원경 경통 + HGA (K=4, 비대칭)"""
    bx, by, bz = bus
    x0 = -(bx / 2 + 1.1)
    parts = [
        ("bus", _box(bus)),
        ("wing", _merge(_boom((-bx / 2 * 0.95, 0, 0.5), (x0, 0, 0.5), 0.05),
                        _panel((x0, 0, 0.5), (-1, 0, 0), wing_span, wing_width, 0.03))),
        ("instrument", _cyl(0.55, 1.6, (0, -(by / 2 + 0.75), -0.8), axis=(0, -1, 0))),
        ("hga", _merge(_boom((0.6, by / 2 * 0.95, 1.6), (0.6, by / 2 + 1.0, 1.6), 0.03),
                       _dish(0.4, 0.1, 0.015, (0.6, by / 2 + 1.0, 1.6), axis=(0, 1, 0)))),
    ]
    return parts, dict(bus=list(bus), wing_span=wing_span, wing_width=wing_width)


def build_telescope():
    """우주망원경 (Hubble 유형): 후방 장비부 + 전방 경통 + 덮개 + 날개 2장 + HGA 2개 (K=7)"""
    r_aft, r_fwd = 2.14, 1.53
    parts = [
        ("aft_shroud", _cyl(r_aft, 3.5, (0, 0, 1.75))),
        ("fwd_tube", _cyl(r_fwd, 9.0, (0, 0, 8.0))),
        ("door", _cyl(r_fwd, 0.08, (0, r_fwd + 0.45, 12.5 + 1.35), axis=(0, 0.966, 0.259))),
    ]
    for s, name in ((1, "wing_px"), (-1, "wing_nx")):
        parts.append((name, _merge(_boom((s * r_fwd * 0.95, 0, 8.0), (s * 2.6, 0, 8.0), 0.06),
                                   _box((2.6, 0.03, 7.1), (s * 3.9, 0, 8.0)))))
    for s, name in ((1, "hga_py"), (-1, "hga_ny")):
        c = (0, s * 4.2, 6.0)
        parts.append((name, _merge(_boom((0, s * r_fwd * 0.95, 6.0), c, 0.04),
                                   _dish(0.65, 0.15, 0.02, c, axis=(0, s, 0)))))
    return parts, {}


def build_soyuz():
    """Soyuz 유형: 기계선(원뿔대) + 귀환선(종형) + 궤도선(구) + 날개 2장 (K=5)"""
    parts = [
        ("service_module", _cyl(1.36, 2.3, (0, 0, 1.15), r_top=1.1)),
        ("descent_module", _cyl(1.1, 2.1, (0, 0, 3.35), r_top=0.62)),
        ("orbital_module", _merge(_sphere(1.1, (0, 0, 5.55), scale=(1, 1, 1.15)),
                                  _cyl(0.25, 0.6, (0, 0, 6.95)))),
    ]
    for s, name in ((1, "wing_px"), (-1, "wing_nx")):
        parts.append((name, _panel((s * 1.2, 0, 1.2), (s, 0, 0), 4.1, 1.0, 0.03)))
    return parts, {}


def build_capsule():
    """캡슐 + 트렁크 (Crew Dragon 유형): 매끈하게 이어진 접합 (K=2, 병목 없음)"""
    r = 1.85
    parts = [
        ("trunk", _merge(_cyl(r, 3.7, (0, 0, 1.85)),
                         _box((0.6, 0.06, 2.0), (r + 0.25, 0, 1.2)),
                         _box((0.6, 0.06, 2.0), (-r - 0.25, 0, 1.2)))),
        ("capsule", _merge(_cyl(r, 3.2, (0, 0, 5.3), r_top=0.75),
                           _sphere(0.75, (0, 0, 6.9), scale=(1, 1, 0.6)))),
    ]
    return parts, {}


# ---------------------------------------------------------------- C. 대형 / D. 잔해
def build_station_t():
    """T자형 우주정거장 (톈궁 유형): 코어 + 노드 + 실험모듈 2 + 날개 6장 (K=10)

    코어 축 z, 실험모듈 축 ±x, 패널 법선 y
    """
    parts = [
        ("core", _merge(_cyl(2.1, 8.0, (0, 0, -5.4)), _cyl(1.4, 5.0, (0, 0, -11.9)),
                        _cyl(2.1, 0.8, (0, 0, -9.4), r_top=1.4, axis=(0, 0, -1)))),
        ("node", _sphere(1.5)),
    ]
    for s, lab in ((1, "lab_px"), (-1, "lab_nx")):
        parts.append((lab, _cyl(2.1, 13.0, (s * 7.8, 0, 0), axis=(s, 0, 0))))
        xe = s * 15.2
        truss = _boom((s * 14.3, 0, 0), (xe, 0, 0), 0.25)
        for t, side in ((1, "a"), (-1, "b")):
            w = _panel((xe, 0, t * 0.9), (0, 0, t), 12.0, 5.0, 0.04)
            parts.append((f"wing_{lab}_{side}", _merge(w, truss) if t == 1 else w))
    for s, name in ((1, "wing_core_px"), (-1, "wing_core_nx")):
        parts.append((name, _merge(_boom((s * 1.3, 0, -12.5), (s * 2.2, 0, -12.5), 0.08),
                                   _panel((s * 2.2, 0, -12.5), (s, 0, 0), 9.0, 3.0, 0.04))))
    return parts, {}


def build_upper_stage():
    """로켓 상단 (H-IIA 2단 유형 잔해): 탱크 + 엔진 노즐 + 탑재체 어댑터 (K=3)"""
    parts = [
        ("tank", _merge(_cyl(2.0, 9.2), _sphere(2.0, (0, 0, -4.6), scale=(1, 1, 0.45)))),
        ("nozzle", _cyl(1.1, 2.6, (0, 0, -6.6), r_top=0.35)),
        ("adapter", _cyl(2.0, 1.4, (0, 0, 5.3), r_top=0.8)),
    ]
    return parts, {}


BUILDERS = {
    "cubesat6u": build_cubesat6u,
    "cubesat6u_mid": build_cubesat6u_mid,
    "cubesat3u_petal": build_cubesat3u_petal,
    "hexsat": build_hexsat,
    "spin_drum": build_spin_drum,
    "geo_comsat": build_geo_comsat,
    "eo_sat": build_eo_sat,
    "telescope": build_telescope,
    "soyuz": build_soyuz,
    "capsule": build_capsule,
    "station_t": build_station_t,
    "upper_stage": build_upper_stage,
}


# ---------------------------------------------------------------- 저장/로드
def save_target(parts, out_glb, target_name, params=None):
    """부품 목록을 GLB(부품별 노드) + <이름>_parts.json 으로 저장"""
    scene = trimesh.Scene()
    for name, m in parts:
        scene.add_geometry(m, node_name=name, geom_name=name)
    os.makedirs(os.path.dirname(out_glb) or ".", exist_ok=True)
    scene.export(out_glb)

    V = np.vstack([m.vertices for _, m in parts])
    center, R = bounding_sphere(V)
    info = dict(target=target_name, unit="m",
                parts=[dict(id=i + 1, name=n, n_faces=int(len(m.faces)),
                            area=float(m.area)) for i, (n, m) in enumerate(parts)],
                center=center.tolist(), R=float(R), params=params or {})
    with open(os.path.splitext(out_glb)[0] + "_parts.json", "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    return info


def load_target(glb_path, parts_map=None):
    """GLB → 부품별 메시 목록 [(name, Trimesh)] (노드 변환 적용된 body 좌표)

    parts_map : {부품이름: [노드이름, ...]} 를 주면 여러 노드를 한 부품으로 묶는다 (ISS 등).
                None 이면 노드 하나 = 부품 하나. <이름>_parts.json 이 있으면 그 순서를 따른다.
    """
    scene = trimesh.load(glb_path, force="scene")
    node_meshes = {}
    for node in scene.graph.nodes_geometry:
        T, gname = scene.graph[node]
        g = scene.geometry[gname]
        if not isinstance(g, trimesh.Trimesh):
            continue
        node_meshes[node] = g.copy().apply_transform(T)

    if parts_map is None:
        order = list(node_meshes)
        js = os.path.splitext(glb_path)[0] + "_parts.json"
        if os.path.exists(js):
            with open(js, encoding="utf-8") as f:
                names = [p["name"] for p in json.load(f)["parts"]]
            order = [n for n in names if n in node_meshes] + [n for n in order if n not in names]
        parts_map = {n: [n] for n in order}

    parts = []
    for pname, nodes in parts_map.items():
        ms = [node_meshes[n] for n in nodes if n in node_meshes]
        if ms:
            parts.append((pname, trimesh.util.concatenate(ms)))
    return parts


def bounding_sphere(V):
    """AABB 중심 기준 경계구 (타깃 스케일 R 정의용, 최소 경계구는 아님)"""
    center = (V.min(0) + V.max(0)) / 2
    R = np.linalg.norm(V - center, axis=1).max()
    return center, R
