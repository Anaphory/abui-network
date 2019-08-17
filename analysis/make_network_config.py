import sys
from pathlib import Path
from bs4 import BeautifulSoup as BS


species = {
    "A-Atimelang": ["atim1239", "abui1241-fuime"],
    "A-Takalelang": ["abui1241-takal"],
    "A-Petleng": ["abui1241-petle"],
    "A-Ulaga": ["abui1241-ulaga"],
    "Papuna": ["abui1241-papun"],
    "Tiyai": ["tiay1238"],
    "Kamang": ["kama1365"],
    "Suboo": ["sibo1242", "sibo1242-atii"],
    "Kui": ["kuii1253"],
    "Kabola": ["kabo1247"],
    "Hamap": ["hama1240"],
    "Kafoa": ["kafo1240"],
    "Klon": ["kelo1247-bring", "kelo1247-hopte"],
}

n_trees = 5


with Path("../beastling/abui-neighbours-spdc.xml").open("r") as original:
    file_tree = BS(original, "xml")

# STEP 1: Update TaxonSet
taxonset = file_tree.find("taxonset")
taxa = taxonset.plate["range"].split(",")
taxonset.string = ""
for species_name, individuals in species.items():
    species_tag = file_tree.new_tag("taxon", id=species_name, spec="TaxonSet")
    taxonset.append(species_tag)
    for identifier in individuals:
        species_tag.append(file_tree.new_tag("taxon", id=identifier, spec="Taxon"))
        taxa.remove(identifier)
assert not taxa

# STEP 2: Change Comment
...


tree_index = []

# STEP 7: Use gene trees in the likelihoods
for likelihood in file_tree.find_all(spec="TreeLikelihood"):
    del likelihood["tree"]
    basic_id = likelihood["id"]
    likelihood.insert(0, BS("""
     <tree id="{:}.t" spec="speciesnetwork.TreeSelector">
      <plate range="{:}" var="tree">
       <tree idref="tree:gene$(tree)"/>
      </plate>
      <index idref="index:{:}"/>
     </tree>
    """.format(basic_id,
               ",".join(map(str, range(n_trees))),
               basic_id), "xml"))
    tree_index.append("index:" + basic_id)

# STEP 4: Define the tree selector choice parameter
state = file_tree.find("state")
state_plate = state.find_all("plate")[0]
assert state_plate["var"] == "rate"
state_plate.append(BS("""
    <stateNode id="index:featureLikelihood:vocabulary:$(rate)" lower="0" spec="beast.core.parameter.IntegerParameter" upper="{:d}">
     0
    </stateNode>""".format(n_trees - 1), "xml"))

# STEP 3: Change state concerning tree
tree = state.tree
print(tree)
tree["id"] = "network:species"
tree["spec"] = "speciesnetwork.NetworkParser"
tree["tree"] = "@newick:species"
del tree["name"]
tree.name = "stateNode"
print(tree)

# STEP 4: Change the initialization
init = file_tree.find("init")
init_parent, init_spot = init.parent, init.parent.index(init)
init.extract()
init_parent.insert(init_spot, BS("""
  <init estimate="false" id="SNI" method="random" origin="@originTime:species" spec="speciesnetwork.SpeciesNetworkInitializer" speciesNetwork="@network:species">
   <plate range="{0:}" var="tree">
    <geneTree idref="tree:gene$(tree)"/>
   </plate>
   <rebuildEmbedding id="initEmbed" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="0.0">
    <plate range="{0:}" var="tree">
     <geneTree idref="tree:gene$(tree)"/>
    </plate>
   </rebuildEmbedding>
  </init>
""".format(",".join(map(str, range(n_trees)))), "xml"))

def turn_into_tree(l):
    if len(l) == 1:
        return "{:}".format(l[0]), 0
    else:
        i = len(l) // 2
        t1, h1 = turn_into_tree(l[:i])
        t2, h2 = turn_into_tree(l[i:])
        h = max(h1, h2) + 1
        return "({:s}:{:d},{:s}:{:d})".format(
            t1, h - h1,
            t2, h - h2), h

starting_tree, _ = turn_into_tree(list(species))

init_parent.insert(init_spot, BS("""
  <init IsLabelledNewick="true" adjustTipHeights="false" id="newick:species" newick="{:}" spec="beast.util.TreeParser"/>
""".format(starting_tree), "xml"))

state.insert(3, BS("""
   <plate range="{:}" var="tree">
    <stateNode id="tree:gene$(tree)" spec="speciesnetwork.EmbeddedTree">
     <taxonset alignment="@data_vocabulary" id="taxonset:gene$(tree)" spec="TaxonSet"/>
    </stateNode>
   </plate>
""".format(",".join(map(str, range(n_trees)))), "xml"))
state.insert(4, BS("""
   <parameter id="turnOverRate:species" lower="0.0" name="stateNode" upper="1.0">
    0.5
   </parameter>
""", "xml"))
state.insert(4, BS("""
   <parameter id="originTime:species" lower="0.0" name="stateNode">
    0.1
   </parameter>
""", "xml"))

# STEP 5: Add the multispecies coalescent
posterior = file_tree.find(id="posterior")
posterior.insert(0, BS("""
   <distribution id="coalescent" spec="speciesnetwork.MultispeciesCoalescent" speciesNetwork="@network:species">
    <plate range="{:}" var="tree">
     <geneTreeWithin geneTree="@tree:gene$(tree)" id="geneTree:gene$(tree)" spec="speciesnetwork.GeneTreeInSpeciesNetwork" speciesNetwork="@network:species"/>
    </plate>
    <populationModel alpha="5.0" id="popModel" mean="2.0" spec="speciesnetwork.ConstantPopIntegrated"/>
   </distribution>
""".format(",".join(map(str, range(n_trees)))), "xml"))

# STEP 6: Add the relevant bits to the prior
tree_prior = file_tree.find(id="YuleModel.t:beastlingTree")
tree_prior_parent, tree_prior_spot = tree_prior.parent, tree_prior.parent.index(tree_prior)
tree_prior.extract()
tree_prior_parent.insert(init_spot, BS("""
    <prior id="turnOverPrior" name="distribution" x="@turnOverRate:species">
     <Beta alpha="1.0" beta="10.0" id="betadistr.01" name="distr"/>
    </prior>
""", "xml"))

tree_prior_parent.insert(init_spot, BS("""
    <prior id="networkOrigin" name="distribution" x="@originTime:species">
     <LogNormal M="6.907755278982137" S="2.302585092994046" id="exponential.0" name="distr" offset="0.0"/>
    </prior>
""", "xml"))

tree_prior_parent.insert(init_spot, BS("""
    <distribution id="networkPrior" netDiversification="@birthRate.t:beastlingTree" network="@network:species" spec="speciesnetwork.BirthHybridizationModel" turnOver="@turnOverRate:species"/>
""", "xml"))

# STEP 7: Wrap the tree operators
operators = []
for tree_operator in file_tree.find_all(name="operator", tree="@Tree.t:beastlingTree"):
    wt = tree_operator["weight"]
    tree_operator["weight"] = "0.0"
    tree_operator["tree"] = "@tree:gene$(tree)"
    tree_operator["id"] = tree_operator["id"].replace("beastlingTree", "gene$(tree)")
    operators.append((tree_operator.extract(), wt))
operator_plate = file_tree.new_tag("plate", range=",".join(map(str, range(n_trees))), var="tree")
run = file_tree.find("run")
run.insert(5, operator_plate)
for operator, weight in operators:
    rebuild = file_tree.new_tag("op", spec="speciesnetwork.operators.RebuildEmbedding", speciesNetwork="@network:species", taxonset="@taxa", weight=weight)
    rebuild["name"] = "operator"
    rebuild.append(file_tree.new_tag("geneTree", idref="tree:gene$(tree)"))
    rebuild.append(operator)
    operator_plate.append(rebuild)

# STEP 8: Add more operators
for operator in [
        """<operator id="turnOverScale:species" parameter="@turnOverRate:species" scaleFactor="0.5" spec="ScaleOperator" weight="10.0"/>""",
        """<operator id="gammaProbUniform:species" spec="speciesnetwork.operators.GammaProbUniform" speciesNetwork="@network:species" weight="20.0"/>""",
        """<operator id="gammaProbRndWalk:species" spec="speciesnetwork.operators.GammaProbRndWalk" speciesNetwork="@network:species" weight="20.0"/>""",
        """<operator id="originMultiplier:species" origin="@originTime:species" spec="speciesnetwork.operators.OriginMultiplier" speciesNetwork="@network:species" weight="5.0"/>""",
        """<operator id="networkMultiplier:species" origin="@originTime:species" spec="speciesnetwork.operators.NetworkMultiplier" speciesNetwork="@network:species" weight="10.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="coordNodeUniformAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="80.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="coordNodeUniform:species" spec="speciesnetwork.operators.CoordinatedNodeUniform" speciesNetwork="@network:species" weight="0.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        </operator>
        </operator>""".format(",".join(map(str, range(n_trees))), ",".join(map(str, range(n_trees)))),
        """<operator id="nodeSliderAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="80.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="nodeSlider:species" isNormal="true" origin="@originTime:species" sigma="0.005" spec="speciesnetwork.operators.NodeSlider" speciesNetwork="@network:species" weight="0.0"/>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="relocateBranchAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="150.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="relocateBranch:species" spec="speciesnetwork.operators.RelocateBranch" speciesNetwork="@network:species" weight="0.0"/>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="flipReticulationAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="15.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="flipReticulation:species" spec="speciesnetwork.operators.FlipReticulation" speciesNetwork="@network:species" weight="0.0"/>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="addReticulationAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="80.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="addReticulation:species" spec="speciesnetwork.operators.AddReticulation" speciesNetwork="@network:species" weight="0.0"/>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="deleteReticulationAndEmbed:species" spec="speciesnetwork.operators.RebuildEmbedding" speciesNetwork="@network:species" taxonset="@taxa" weight="80.0">
        <plate range="{:}" var="tree">
        <geneTree idref="tree:gene$(tree)"/>
        </plate>
        <operator id="deleteReticulation:species" spec="speciesnetwork.operators.DeleteReticulation" speciesNetwork="@network:species" weight="0.0"/>
        </operator>""".format(",".join(map(str, range(n_trees)))),
        """<operator id="clockScaler.c:default" parameter="@clockRate.c:default" scaleFactor="0.5" spec="ScaleOperator" weight="3.0"/>"""]:
    run.append(BS(operator, "xml"))
file_tree.find(id="UpDown").extract()

# STEP 9: Add switch operators
run.append(BS("""
  <plate range="{:}" var="gene">
   <operator id="changeTree:$(gene)" parameter="@$(gene)" spec="beast.evolution.operators.UniformOperator" weight="100"/>
  </plate>""".format(",".join(tree_index)), "xml"))


# STEP 10: Fix loggers
file_tree.find(id="treeStats").extract()
file_tree.find(id="treeLogger").extract()

run.append(BS("""
  <logger fileName="networks.newick" id="specieslog" logEvery="100000" mode="tree">
   <log id="networkLogger:species" spec="speciesnetwork.NetworkWithMetaDataLogger" speciesNetwork="@network:species"/>
  </logger>""", "xml"))
run.append(BS("""
  <plate range="{:}" var="tree">
   <logger fileName="genetrees$(tree).trees" id="treelog:gene$(tree)" logEvery="100000" mode="tree">
    <log dp="4" id="treeLogger:gene$(tree)" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@tree:gene$(tree)"/>
   </logger>
""".format(",".join(map(str, range(n_trees)))), "xml"))

# Output the results
with Path("../beastling/NW-abui-neighbours-spdc.xml").open("w") as new:
    new.write(file_tree.prettify())
