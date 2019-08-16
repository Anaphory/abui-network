import sys
from pathlib import Path
import lxml

from lxml.etree import parse, tostring, fromstring, SubElement


species = {
    "A-Atimelang": ["atim1239", "abui1241-fuime"],
    "A-Takalelang": ["abui1241-takal"],
    "A-Petleng": ["abui1241-petle"],
    "A-Ulaga": ["abui1241-ulaga"],
    "Kafoa": ["kafo1240"],
    "Kamang": ["kama1365"],
    "Kui": ["kuii1253"],
    "Kiraman": ["kira1248"],
    "Klon": ["kelo1247-bring", "kelo1247-hopte"],
}


with Path("../beastling/abui-neighbours-spdc.xml").open("rb") as original:
    file_tree = parse(original)

# STEP 1: Update TaxonSet
taxonset = file_tree.find("./taxonset")
print(tostring(taxonset).decode("utf-8"))
print("====")
taxa = taxonset[0].get("range").split(",")
taxonset[:] = []
for species_name, individuals in species.items():
    species_tag = SubElement(taxonset, "taxon", id=sp, spec="TaxonSet")
    for identifier in individuals:
        SubElement(species_tag, "taxon", id=identifier, spec="TaxonSet")
        taxa.remove(identifier)
print(tostring(taxonset).decode("utf-8"))
assert not taxa

# STEP 2: Update BeastLing comment
comment = file_tree.xpath("//comment()")[0]
comment.text = comment.text[:comment.text.index("Please")] + """
This configuration had network capabilities added
using the script 'make_network_config.py'."""

if False:
    response = input()
    if not response or response.strip().lower().startswith("n"):
        continue
    inner = "".join(
        [plate.text or ""] + [tostring(x).decode("utf-8") for x in plate[:]])
    print(inner)
    parent = plate.getparent()
    index = parent.index(plate)
    replacement = fromstring("<wrapper>" + ("".join(
        inner.replace("$({:})".format(var), val)
        for val in range)) + "</wrapper>")
    parent[index:index + 1] = replacement[:]

for element in file_tree.iter():
    if element.tail and not element.tail.strip():
        element.tail = None
with Path("../beastling/NW-abui-neighbours-spdc.xml").open("wb") as new:
    file_tree.write(
        new,
        pretty_print=True,
        encoding='utf-8')
