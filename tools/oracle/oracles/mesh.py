"""trimesh as the geometry oracle for OBJ and PLY (#340).

Like `gemmi`/`meeko`/`mdanalysis`, this does not implement the `Oracle`
dataclass -- it answers a geometric question, not a SMILES-shaped one. The
bar here is the issue's own: "geometry either matches or does not." A
closed, watertight fixture (a box) makes that concrete -- vertex count,
face count and volume, compared within a small float tolerance rather than
exactly, since neither OBJ's nor PLY's text encoding is required to round a
coordinate the same way twice.
"""

from pathlib import Path
from typing import NamedTuple, Optional

import trimesh


class Summary(NamedTuple):
    vertex_count: int
    face_count: int
    volume: float


def summarize(path: Path) -> Optional[Summary]:
    """Reads a mesh file and summarises what trimesh saw, or `None` if
    trimesh could not read it, or read something with no faces at all.
    """
    try:
        mesh = trimesh.load(str(path), process=False, force="mesh")
    except Exception:
        return None
    if len(mesh.faces) == 0:
        return None
    return Summary(len(mesh.vertices), len(mesh.faces), float(mesh.volume))


def write_fixture(path: Path) -> None:
    """Authors a small, closed (watertight) fixture -- a unit cube -- via
    trimesh's own writer. Format is inferred from `path`'s extension,
    matching how trimesh's own API expects to be called.
    """
    trimesh.creation.box(extents=[2.0, 2.0, 2.0]).export(str(path))


def load_trimesh() -> None:
    """Confirms trimesh can write and read back its own OBJ before its
    answers are trusted -- the same discipline as `load_gemmi`/
    `load_meeko`/`load_mdanalysis`.
    """
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "sanity.obj"
        write_fixture(path)
        summary = summarize(path)
    if summary is None or summary.vertex_count == 0 or summary.face_count == 0:
        raise RuntimeError(
            "trimesh cannot write and read back its own OBJ cube, so it is "
            "broken rather than strict — refusing to report its answers as "
            "findings"
        )
