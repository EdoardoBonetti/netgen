# ngmesh.SplitAlfeld puts the new vertex at the barycentre: the split mesh has
# the same volume and no inverted elements; SplitPowellSabin (2D) is conforming
import numpy as np
from netgen.csg import unit_cube
from netgen.geom2d import unit_square
import pytest


def tet_volumes(mesh):
    p = np.array([q.p for q in mesh.Points()])
    vols = []
    for el in mesh.Elements3D():
        v = p[[pi.nr - 1 for pi in el.vertices]]
        vols.append(np.linalg.det(np.array([v[1] - v[0], v[2] - v[0], v[3] - v[0]])) / 6)
    return np.array(vols)


def trig_areas(mesh):
    p = np.array([q.p for q in mesh.Points()])[:, :2]
    areas = []
    for el in mesh.Elements2D():
        v = p[[pi.nr - 1 for pi in el.vertices]]
        areas.append(np.linalg.det(np.array([v[1] - v[0], v[2] - v[0]])) / 2)
    return np.array(areas)


def edge_counts(elements):
    count = {}
    for el in elements:
        vs = [pi.nr for pi in el.vertices]
        for i in range(len(vs)):
            e = tuple(sorted((vs[i], vs[(i + 1) % len(vs)])))
            count[e] = count.get(e, 0) + 1
    return count


def test_split_alfeld_3d():
    mesh = unit_cube.GenerateMesh(maxh=0.5)
    vol0 = tet_volumes(mesh)
    ne = len(mesh.Elements3D())
    mesh.SplitAlfeld()
    vol = tet_volumes(mesh)
    assert len(vol) == 4 * ne
    assert abs(vol.sum() - vol0.sum()) < 1e-12
    assert (np.sign(vol) == np.sign(vol0[0])).all()   # no inverted elements


@pytest.mark.parametrize("split", ["SplitAlfeld", "SplitPowellSabin"])
def test_split_2d(split):
    mesh = unit_square.GenerateMesh(maxh=0.3)
    area0 = trig_areas(mesh)
    ne = len(mesh.Elements2D())
    getattr(mesh, split)()
    area = trig_areas(mesh)
    assert len(area) == (3 if split == "SplitAlfeld" else 6) * ne
    assert abs(area.sum() - area0.sum()) < 1e-12
    assert (np.sign(area) == np.sign(area0[0])).all()
    # conforming: every edge is shared by two triangles, or lies on the boundary
    count = edge_counts(mesh.Elements2D())
    bnd = edge_counts(mesh.Elements1D())
    for e, c in count.items():
        assert c == (1 if e in bnd else 2)
