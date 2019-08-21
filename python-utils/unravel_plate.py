import sys
from pathlib import Path
import lxml

from lxml.etree import parse, tostring, fromstring

with Path(sys.argv[1]).open("rb") as original:
    file_tree = parse(original)
for plate in file_tree.findall(".//plate"):
    print(tostring(plate).decode("utf-8"))
    response = input()
    if not response or response.strip().lower().startswith("n"):
        continue
    range = plate.get("range").split(",")
    var = plate.get("var")
    inner = "".join(
        [plate.text or ""] + [tostring(x).decode("utf-8") for x in plate[:]])
    print(inner)
    parent = plate.getparent()
    index = parent.index(plate)
    replacement = fromstring("<wrapper>" + ("".join(
        inner.replace("$({:})".format(var), val)
        for val in range)) + "</wrapper>")
    parent[index:index + 1] = replacement[:]

with Path(sys.argv[-1]).open("wb") as new:
    file_tree.write(
        new,
        pretty_print=True,
        encoding='utf-8')
