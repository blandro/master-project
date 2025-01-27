import argparse
import pandas as pd
import rdflib

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",type=str)
args = parser.parse_args()
taxa = args.input

u=rdflib.Graph()
u.parse(f"{taxa}.rdf",format="xml")
print(f"{taxa} rdf loaded")

query_results=dict()
query_dict=dict()

query_dict["uniReactions"]="""
    PREFIX up: <http://purl.uniprot.org/core/>
    SELECT ?reaction ?protein ?direction
    WHERE {
       ?protein rdf:type up:Protein .
       ?protein up:annotation ?b .
       ?b up:catalyticActivity ?a .
       ?a up:catalyzedReaction ?reaction .
       OPTIONAL { ?b up:catalyzedPhysiologicalReaction ?direction }      
        }
    """

for key in query_dict.keys():
        print(key)
        result = u.query(query_dict[key])
        query_results[key] = pd.DataFrame(result.bindings).applymap(str).rename(columns=str)
        query_results[key] = query_results[key].drop_duplicates()
        query_results[key].to_csv(key + ".tsv", sep="\t", index=False)