import requests

test_url = "https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecules/5j7l"

response = requests.get(test_url)

print(response.text)


