"""Writes `nmr-ensemble.bcif` from `../mmcif/nmr-ensemble.cif` (#404).

`chem convert` can't produce this file: its BinaryCIF writer encodes one
model per block, and the point of the fixture is three models in one
`_atom_site` category. So the bytes are assembled here, stdlib only, with a
deliberately different encoding choice from the crate's writer (plain
`Int32` integers, `Int32` string indices and offsets) so the decoder is
also checked against bytes it didn't write.

Run from anywhere: `python3 generate_nmr_ensemble.py`.
"""

import struct
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "mmcif" / "nmr-ensemble.cif"
TARGET = HERE / "nmr-ensemble.bcif"

FLOAT_COLUMNS = {"Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"}
INT_COLUMNS = {"id", "auth_seq_id", "pdbx_PDB_model_num"}


def pack(value) -> bytes:
    """The MessagePack subset BinaryCIF needs: map, array, str, bin, int."""
    if isinstance(value, dict):
        n = len(value)
        head = bytes([0x80 | n]) if n < 16 else b"\xde" + struct.pack(">H", n)
        return head + b"".join(pack(k) + pack(v) for k, v in value.items())
    if isinstance(value, list):
        n = len(value)
        head = bytes([0x90 | n]) if n < 16 else b"\xdc" + struct.pack(">H", n)
        return head + b"".join(pack(v) for v in value)
    if isinstance(value, str):
        raw = value.encode()
        n = len(raw)
        if n < 32:
            head = bytes([0xA0 | n])
        elif n < 256:
            head = b"\xd9" + struct.pack(">B", n)
        else:
            head = b"\xda" + struct.pack(">H", n)
        return head + raw
    if isinstance(value, bytes):
        n = len(value)
        if n < 256:
            head = b"\xc4" + struct.pack(">B", n)
        else:
            head = b"\xc5" + struct.pack(">H", n)
        return head + value
    if isinstance(value, int):
        return b"\xd2" + struct.pack(">i", value)
    raise TypeError(type(value))


def byte_array(type_code: int) -> dict:
    return {"kind": "ByteArray", "type": type_code}


def encode(name: str, values: list[str]) -> dict:
    if name in FLOAT_COLUMNS:
        data = b"".join(struct.pack("<d", float(v)) for v in values)
        encoding = [byte_array(33)]
    elif name in INT_COLUMNS:
        data = b"".join(struct.pack("<i", int(v)) for v in values)
        encoding = [byte_array(3)]
    else:
        table: list[str] = []
        for v in values:
            if v not in table:
                table.append(v)
        offsets = [0]
        for s in table:
            offsets.append(offsets[-1] + len(s.encode()))
        data = b"".join(struct.pack("<i", table.index(v)) for v in values)
        encoding = [
            {
                "kind": "StringArray",
                "dataEncoding": [byte_array(3)],
                "stringData": "".join(table),
                "offsetEncoding": [byte_array(3)],
                "offsets": b"".join(struct.pack("<i", o) for o in offsets),
            }
        ]
    return {"name": name, "data": {"encoding": encoding, "data": data}}


def main() -> None:
    lines = [l for l in SOURCE.read_text().splitlines() if l.strip()]
    header = lines[0].removeprefix("data_")
    tags = [l.removeprefix("_atom_site.") for l in lines if l.startswith("_atom_site.")]
    rows = [l.split() for l in lines if l.startswith(("ATOM", "HETATM"))]
    columns = [encode(tag, [row[i] for row in rows]) for i, tag in enumerate(tags)]
    document = {
        "version": "0.3.0",
        "encoder": "generate_nmr_ensemble.py",
        "dataBlocks": [
            {
                "header": header,
                "categories": [
                    {"name": "_atom_site", "rowCount": len(rows), "columns": columns}
                ],
            }
        ],
    }
    TARGET.write_bytes(pack(document))


if __name__ == "__main__":
    main()
