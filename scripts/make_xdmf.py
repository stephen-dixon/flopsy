#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
import h5py


def write_xdmf_for_flopsy_h5(h5_path: Path, xdmf_path: Path | None = None) -> Path:
    h5_path = Path(h5_path)
    if xdmf_path is None:
        xdmf_path = h5_path.with_suffix(".xdmf")
    else:
        xdmf_path = Path(xdmf_path)

    with h5py.File(h5_path, "r") as h5:
        times = h5["time"][:]
        x = h5["mesh/x"][:]
        field_names = sorted(list(h5["fields"].keys()))

    nt = len(times)
    nx = len(x)
    h5_ref = h5_path.name  # assume same directory as .xdmf

    with open(xdmf_path, "w", encoding="utf-8") as f:
        f.write("""<?xml version="1.0" ?>\n""")
        f.write("""<Xdmf Version="3.0">\n""")
        f.write("""  <Domain>\n""")
        f.write("""    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n""")

        for it, tval in enumerate(times):
            f.write(f"""      <Grid Name="step_{it}" GridType="Uniform">\n""")
            f.write(f"""        <Time Value="{float(tval)}"/>\n""")

            f.write(f"""        <Topology TopologyType="1DSMesh" NumberOfElements="{nx}"/>\n""")
            f.write("""        <Geometry GeometryType="X">\n""")
            f.write(
                f"""          <DataItem Format="HDF" Dimensions="{nx}" NumberType="Float" Precision="8">{h5_ref}:/mesh/x</DataItem>\n"""
            )
            f.write("""        </Geometry>\n""")

            for name in field_names:
                f.write(
                    f"""        <Attribute Name="{name}" AttributeType="Scalar" Center="Node">\n"""
                )
                f.write(
                    f"""          <DataItem ItemType="HyperSlab" Dimensions="{nx}" Type="HyperSlab">\n"""
                )
                f.write("""            <DataItem Dimensions="3 2" Format="XML">\n""")
                f.write(f"""              {it} 0\n""")
                f.write("""              1 1\n""")
                f.write(f"""              1 {nx}\n""")
                f.write("""            </DataItem>\n""")
                f.write(
                    f"""            <DataItem Format="HDF" Dimensions="{nt} {nx}" NumberType="Float" Precision="8">{h5_ref}:/fields/{name}</DataItem>\n"""
                )
                f.write("""          </DataItem>\n""")
                f.write("""        </Attribute>\n""")

            f.write("""      </Grid>\n""")

        f.write("""    </Grid>\n""")
        f.write("""  </Domain>\n""")
        f.write("""</Xdmf>\n""")

    return xdmf_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate an XDMF companion file for a Flopsy HDF5 field output file."
    )
    parser.add_argument("h5", type=Path, help="Input Flopsy .h5 file")
    parser.add_argument(
        "--xdmf",
        type=Path,
        default=None,
        help="Output .xdmf path (defaults to input basename with .xdmf)",
    )
    args = parser.parse_args()

    out = write_xdmf_for_flopsy_h5(args.h5, args.xdmf)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
