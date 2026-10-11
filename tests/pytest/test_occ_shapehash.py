# TopoDS_Shape.__hash__ and the global shape maps use NgShapeMapHasher, which
# also works with OpenCascade < 7.8
import pytest
occ = pytest.importorskip("netgen.occ")


def test_shape_hash():
    box = occ.Box((0, 0, 0), (1, 1, 1))
    faces = box.faces
    assert len({hash(f) for f in faces}) == 6
    assert hash(box.faces[0]) == hash(faces[0])
    assert faces[0] == box.faces[0]
    s = {f for f in faces} | {f for f in box.faces}
    assert len(s) == 6


def test_shape_properties():
    # names and maxh are stored in the global shape-property map
    box = occ.Box((0, 0, 0), (1, 1, 1))
    box.faces.Min(occ.X).name = "left"
    box.faces.Min(occ.X).maxh = 0.1
    assert box.faces.Min(occ.X).name == "left"
    assert [f.name for f in box.faces].count("left") == 1
    mesh = occ.OCCGeometry(box).GenerateMesh(maxh=0.5)
    assert "left" in mesh.GetRegionNames(dim=2)
