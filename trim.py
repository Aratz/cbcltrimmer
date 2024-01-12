""" A primitive tool to remove all tiles but the first one from a cbcl file.

Usage:
    python3 trim.py <input_file> <output_file>

Documentation for .cbcl format (p7):
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
"""

import sys

with open(sys.argv[1], "rb") as f:
    version = int.from_bytes(f.read(2), "little")
    header_size = int.from_bytes(f.read(4), "little")
    bits_per_bc = int.from_bytes(f.read(1), "little")
    bits_per_qscore = int.from_bytes(f.read(1), "little")

    n_bins = int.from_bytes(f.read(4), "little")
    q_val_mapping = {
        int.from_bytes(f.read(4), "little"):
            int.from_bytes(f.read(4), "little")
        for _ in range(n_bins)
    }

    n_tile_records = int.from_bytes(f.read(4), "little")

    gzip_virtual_file_offset = {
        int.from_bytes(f.read(4), "little"):
            {
                "n_clusters": int.from_bytes(f.read(4), "little"),
                "uncompressed_size": int.from_bytes(f.read(4), "little"),
                "compressed_size": int.from_bytes(f.read(4), "little"),
            }
        for _ in range(n_tile_records)
    }

    non_PF_flag = f.read(1)

    (tile, metadata) = next(iter(gzip_virtual_file_offset.items()))
    gzip_data = f.read(metadata["compressed_size"])

with open(sys.argv[2], "wb") as f:
    new_header_size = 2 + 4 + 1 + 1 + 4 + 2*4*n_bins + 4 + 4*4 + 1

    f.write(version.to_bytes(2, byteorder="little"))
    f.write(new_header_size.to_bytes(4, byteorder="little"))
    f.write(bits_per_bc.to_bytes(1, byteorder="little"))
    f.write(bits_per_qscore.to_bytes(1, byteorder="little"))

    f.write(n_bins.to_bytes(4, byteorder="little"))
    for fr, to in q_val_mapping.items():
        f.write(fr.to_bytes(4, byteorder="little"))
        f.write(to.to_bytes(4, byteorder="little"))

    f.write((1).to_bytes(4, byteorder="little"))  # n_tile_records

    f.write(tile.to_bytes(4, byteorder="little"))
    f.write(metadata["n_clusters"].to_bytes(4, byteorder="little"))
    f.write(metadata["uncompressed_size"].to_bytes(4, byteorder="little"))
    f.write(metadata["compressed_size"].to_bytes(4, byteorder="little"))

    f.write(non_PF_flag)

    f.write(gzip_data)
