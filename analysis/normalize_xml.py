import sys
from pathlib import Path
import lxml

from lxml.etree import parse, tostring
with Path(sys.argv[1]).open("rb") as original:
    file_tree = parse(original)
for element in file_tree.iter():
    if element.tail and not element.tail.strip():
        element.tail = None
with Path(sys.argv[-1]).open("wb") as new:
    file_tree.write(
        new,
        pretty_print=True,
        encoding='utf-8')
