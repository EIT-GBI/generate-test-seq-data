"""Small shared helpers."""

import gzip
import io


def gzip_text(path, compresslevel=9):
    """Open a .gz for text writing, reproducibly.

    gzip.open() stamps the current time into the header, which makes the
    checksum of otherwise identical output change on every run. mtime=0
    removes the only non-deterministic byte.
    """
    raw = gzip.GzipFile(filename="", mode="wb", compresslevel=compresslevel,
                        fileobj=open(path, "wb"), mtime=0)
    return io.TextIOWrapper(raw, encoding="utf-8")
