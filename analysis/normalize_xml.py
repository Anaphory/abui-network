import sys
from pathlib import Path
from bs4 import BeautifulSoup

with Path(sys.argv[1]).open("r") as original:
    file_tree = BeautifulSoup(original, "xml")
print(file_tree)
with Path(sys.argv[-1]).open("w") as new:
    new.write(file_tree.prettify())
