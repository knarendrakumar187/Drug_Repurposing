// Removes ONLY what write_new_links_to_graph.py added
MATCH ()-[r:EXTERNAL_ACTS_ON]->() DELETE r;
MATCH (x:ExternalDrug) WHERE NOT x:Drug AND NOT (x)--() DELETE x;
