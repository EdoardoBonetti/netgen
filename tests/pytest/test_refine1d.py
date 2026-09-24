"""Adaptive bisection of 1D meshes (Bisect1D in libsrc/meshing/bisect.cpp)."""
import pytest
from netgen.meshing import Mesh, MeshPoint, Element1D, Element0D, Pnt


def make_1d_mesh(n, a=0.0, b=1.0, regions=True):
    m = Mesh(dim=1)
    pids = [m.Add(MeshPoint(Pnt(a + (b - a) * i / n, 0, 0))) for i in range(n + 1)]
    if regions:   # as ngsolve.meshes.Make1DMesh does
        idx = m.AddRegion("dom", dim=1)
        il, ir = m.AddRegion("left", dim=0), m.AddRegion("right", dim=0)
    else:         # bare indices, no region descriptors
        idx, il, ir = 1, 1, 2
    for i in range(n):
        m.Add(Element1D([pids[i], pids[i + 1]], index=idx))
    m.Add(Element0D(pids[0], index=il))
    m.Add(Element0D(pids[-1], index=ir))
    return m


def test_refine_without_edge_regions():
    # a hand-built mesh whose segments carry an index with no region
    # descriptor must still refine (linear midpoint), not crash
    m = make_1d_mesh(4, regions=False)
    for i, seg in enumerate(m.Elements1D()):
        seg.refine = i % 2 == 0
    m.Refine(adaptive=True)
    assert len(list(m.Elements1D())) == 6


def segments(m):
    return [(seg.vertices[0].nr, seg.vertices[1].nr) for seg in m.Elements1D()]


def test_refine_marked_segments():
    m = make_1d_mesh(8)
    for i, seg in enumerate(m.Elements1D()):
        seg.refine = i < 3
    m.Refine(adaptive=True)
    assert len(list(m.Elements1D())) == 11
    xs = sorted(p.p[0] for p in m.Points())
    for x in (1 / 16, 3 / 16, 5 / 16):
        assert any(abs(x - y) < 1e-14 for y in xs)
    # boundary points untouched
    assert len(list(m.Elements0D())) == 2


def test_refine_none_is_noop():
    m = make_1d_mesh(5)
    before = segments(m)
    m.Refine(adaptive=True)
    assert segments(m) == before


def test_refine_all_equals_uniform():
    m1 = make_1d_mesh(6)
    for seg in m1.Elements1D():
        seg.refine = True
    m1.Refine(adaptive=True)
    m2 = make_1d_mesh(6)
    m2.Refine(adaptive=False)
    xs1 = sorted(p.p[0] for p in m1.Points())
    xs2 = sorted(p.p[0] for p in m2.Points())
    assert xs1 == pytest.approx(xs2)


def test_repeated_refinement_and_children_are_conforming():
    m = make_1d_mesh(4)
    for _ in range(3):
        for i, seg in enumerate(m.Elements1D()):
            seg.refine = (i % 2 == 0)
        m.Refine(adaptive=True)
    segs = segments(m)
    coords = [p.p[0] for p in m.Points()]
    base = min(min(a, b) for a, b in segs)          # PointId.nr is base-dependent
    pts = {base + i: c for i, c in enumerate(coords)}
    # every interior point is used by exactly two segments, endpoints by one
    use = {}
    for a, b in segs:
        use[a] = use.get(a, 0) + 1
        use[b] = use.get(b, 0) + 1
    ends = [k for k, v in use.items() if v == 1]
    assert len(ends) == 2
    assert all(v == 2 for k, v in use.items() if k not in ends)
    # segments tile [0,1] without overlap
    lengths = sorted(abs(pts[b] - pts[a]) for a, b in segs)
    assert sum(lengths) == pytest.approx(1.0)
