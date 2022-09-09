from models import ChainInfo
from database import db_session
import csv
import itertools
import operator
from collections import defaultdict


#query = db_session.query(User.id, User.name).filter(User.id.in_([123,456]))

def get_nr_proteins():

    test_list = ['5JTE', '6W7N', '6VWN', '6W7M', '6VYS', '6OSK', '6OT3', '6WDE', '6WDD', '6WDG', '6WDF', '6WDA', '6WDB', '5J5B', '7JSZ', '6VYQ', '6WDI', '6WDK', '6WDJ', '7D6Z', '4V54', '6OSQ', '6XDQ', '4V57', '5U4J', '5MDZ', '6ENF', '7JSW', '6NQB', '5UYN', '7JT3', '7JT2', '7JT1', '5UYQ', '5UYP', '5MDY', '5WDT', '7OE0', '6GXN', '7JSS', '6WD6', '6H4N', '5UYK', '4V89', '5UYM', '5UYL', '5WE6', '6WD5', '6WD4', '7N30', '7N31', '6WD1', '6WD0', '6WD3', '5WE4', '5KCS', '4V4Q', '6WD9', '6WD8', '4V4H', '5KCR', '5IT8', '7OE1', '7ACJ', '6X7K', '7AFI', '4WOI', '6C4I', '5JC9', '5L3P', '5KPW', '5JU8', '6WDM', '6WNV', '6LKQ', '6WDL', '7K00', '7ACR', '3J9Y', '3J9Z', '6WD7', '6VYR', '6BU8', '7OJ0', '6W7W', '5WF0', '4V55', '4WWW', '6TBV', '4V56', '4V50', '4V53', '4V52', '4WF1', '5H5U', '7K52', '4V5B', '7N2U', '7N2V', '6SZS', '5WFS', '5O2R', '7BOI', '7ABZ', '5WFK', '7NAS', '6X6T', '7N2C', '7BOF', '6Y69', '5LZD', '7AFH', '5LZA', '7O5H', '6YSR', '6YSS', '6O9J', '6O9K', '6YST', '6YSU', '6WD2', '7AC7', '6OGF', '4V85', '5MDV', '5MDW', '4U27', '4U26', '4U25', '4U24', '6OGI', '4U20', '5KPS', '6TC3', '6GWT', '6W77', '6GXM', '5KPX', '5IQR', '3JBU', '4V9P', '6ORL', '3JBV', '6Q9A', '6DNC', '4U1V', '7BOG', '6GXO', '4U1U', '5NWY', '4V9C', '5U4I', '7AFN', '4V9D', '7AFA', '6X7F', '7OIZ', '4V9O', '7AFO', '6WNW', '7AF8', '5MGP', '7AF3', '7AF5', '7OI0', '6Q97', '6XZA', '6XZB', '7BOD', '3JCJ', '7O19', '7NSP', '7NSQ', '5J91', '6VU3', '3JCD', '3JCE', '6I7V', '5LZE', '6OG7', '6W6K', '7AFD', '7NSO', '7N1P', '4V64', '5J7L', '5AFI', '6BY1', '6ENU', '6VWM', '4V7V', '4V7U', '4V7T', '4V7S', '3JA1', '7LV0', '6OUO', '6ENJ', '6ZU1', '7K51', '7K50', '7K53', '5U9F', '7K55', '7K54', '7B5K', '5J8A', '6VWL', '5NP6', '6OFX', '6ZTP', '4YBB', '7O1A', '7O1C', '6OM6', '6ORE', '5U9G', '5J88', '4V6D', '4V6E', '6ZTJ', '4V6C', '6ZTL', '6ZTM', '6ZTN', '6ZTO']
    nr_proteins = set()

    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(ChainInfo.pdb_id.in_(test_list)) \
                       .filter(ChainInfo.entity_macromolecule_type == 'Polypeptide(L)') \

    for row in query:
        nr_proteins.add((row.pdb_id, row.compound))

    return list(nr_proteins)

def group_elements(l):
    it = itertools.groupby(l, operator.itemgetter(1))
    for k, v in it:
        yield k, [item[0] for item in v]


proteins_list = get_nr_proteins()

substr = ['ribosomal protein', 'RIBOSOMAL PROTEIN']

filter_data = [x for x in proteins_list if all(y not in x[1] for y in substr)]

#filter_data = [x for x in proteins_list if substr not in x]

d = defaultdict(list)

for k, v in filter_data:
    d[v].append(k)

with open('dict.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in d.items():
       writer.writerow([key, value])



