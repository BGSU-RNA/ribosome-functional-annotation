from database import db_session
from models import UnitPairsDistances, UnitInfo, ChainInfo
from sqlalchemy import tuple_
import csv
import time

test_units = ['5J7L|1|AA|A|964', '5J7L|1|AA|A|968', '5J7L|1|AA|A|969', '5J7L|1|AA|C|970', '5J7L|1|AA|C|972', '5J7L|1|AA|G|963', '5J7L|1|AA|G|971', '5J7L|1|AA|U|965']

def get_protein_contacts(units):
    """A function that returns the protein residues that are located
       within 8.5 Angstroms of the loop nts

    Args:
        units (list of strings): A list containing the unitid of loop nts

    Returns:
        list of tuples: A list of 3-elements tuple in the form (rna_unit_id, protein_unit_id, distance)
    """

    contacts_list = []
    with db_session() as session:
        query = session.query(UnitPairsDistances) \
                       .join(UnitInfo, UnitPairsDistances.unit_id_2 == UnitInfo.unit_id) \
                       .filter(UnitInfo.unit_type_id == 'aa') \
                       .filter(UnitPairsDistances.distance <= 8.5) \
					   .filter(UnitPairsDistances.unit_id_1.in_(units))

        return [(row.unit_id_1, row.unit_id_2, row.distance) for row in query]

def get_bound_protein_chains(rna_chain):
    """A function that returns all the protein chains bound to a particular RNA chain

    Args:
        rna_chain (string): RNA chain

    Returns:
        set of tuples: A set of 2-elements tuple in the form of (pdbid, chain)
    """

    bound_protein_chains = set()
    search_keyword = rna_chain + "%"
    with db_session() as session:
        query = session.query(UnitPairsDistances) \
                       .join(UnitInfo, UnitPairsDistances.unit_id_2 == UnitInfo.unit_id) \
                       .filter(UnitPairsDistances.unit_id_1.like(search_keyword)) \
                       .filter(UnitInfo.unit_type_id == 'aa') \
                       .filter(UnitPairsDistances.distance <= 8.5) 
    
        for row in query:
            protein_info_list = row.unit_id_2.split("|")
            pdb, _, chain, _ = protein_info_list[0], protein_info_list[1], protein_info_list[2], protein_info_list[3:]
            protein_info = (pdb, chain)
            bound_protein_chains.add(protein_info)

        return list(bound_protein_chains)

def get_protein_names(protein_info):
    """A function that returns the names of protein chains

    Args:
        protein_info (set of tuples): A set of 2-elements tuple in the form of (pdbid, chain)

    Returns:
        dictionary: The key is the chain_id while the value is the compound name
    """

    chain_names = {}
    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(tuple_(ChainInfo.pdb_id, ChainInfo.chain_name) \
					   .in_(protein_info))

        for row in query:
            chain_names[row.chain_name] = row.compound
            #chain_names['pdb'] = row.pdb_id
            #chain_names['chain'] = row.chain_name
            #chain_names['compound'] = row.compound

        return chain_names

def lookup(dct, *args):
    for keyword in args:
        dct = {key: value for key, value in dct.items() if keyword not in value}
    return dct

filter_string = ['ribosomal protein']

def run(test_chain):

    pdb_id = test_chain.split('|')[0]

    protein_chains = get_bound_protein_chains(test_chain)
    protein_names = get_protein_names(protein_chains)
    filtered_dict = lookup(protein_names, "ribosomal protein", "RIBOSOMAL PROTEIN")

    if filtered_dict:
        with open('protein_chain_annotation.csv', 'a') as f:  
            #w = csv.DictWriter(f, filtered_dict.keys())
            #w.writerow(filtered_dict)

            writer = csv.writer(f)
            for key, value in filtered_dict.items():
                writer.writerow([pdb_id, key, value])


if __name__ == "__main__":
    test_chains = ['4V4Q|1|BB', '4V4Q|1|DB', '4V50|1|BB', '4V50|1|DB', '4V5B|1|AB', '4V5B|1|CB', '4V9D|1|CA', '4V9D|1|DA', '4V9O|1|AA', '4V9O|1|CA', '4V9O|1|EA', '4V9O|1|GA', '4V9P|1|AA', '4V9P|1|CA', '4V9P|1|EA', '4V9P|1|GA', '4YBB|1|DA', '4YBB|1|CA', '5IT8|1|DA', '5IT8|1|CA', '5J5B|1|DA', '5J5B|1|CA', '5J7L|1|DA', '5J7L|1|CA', '5J88|1|DA', '5J88|1|CA', '5J8A|1|DA', '5J8A|1|CA', '5J91|1|DA', '5J91|1|CA', '5JC9|1|DA', '5JC9|1|CA', '5MDZ|1|1', '6BU8|1|01', '6GWT|1|A', '6GXM|1|A', '6GXN|1|A', '6GXO|1|A', '6I7V|1|CA', '6I7V|1|DA', '4U1U|1|BA', '4U1U|1|DA', '4U1V|1|BA', '4U1V|1|DA', '4U20|1|BA', '4U20|1|DA', '4U24|1|BA', '4U24|1|DA', '4U25|1|BA', '4U25|1|DA', '4U26|1|BA', '4U26|1|DA', '4U27|1|BA', '4U27|1|DA', '4V4H|1|BB', '4V4H|1|DB', '4V52|1|BB', '4V52|1|DB', '4V53|1|BB', '4V53|1|DB', '4V54|1|BB', '4V54|1|DB', '4V55|1|BB', '4V55|1|DB', '4V56|1|BB', '4V56|1|DB', '4V57|1|BB', '4V57|1|DB', '4V64|1|BB', '4V64|1|DB', '4V7S|1|BA', '4V7S|1|DA', '4V7T|1|BA', '4V7T|1|DA', '4V7U|1|BA', '4V7U|1|DA', '4V7V|1|BA', '4V7V|1|DA', '4V9C|1|BA', '4V9C|1|DA', '4WF1|1|BA', '4WF1|1|DA', '4WOI|1|BA', '4WOI|1|CA', '4WWW|1|RA', '4WWW|1|YA', '5KCR|1|1A', '5KCS|1|1A', '3J9Y|1|A', '3JCE|1|A', '5AFI|1|A', '5UYK|1|01', '5UYL|1|01', '5UYM|1|01', '5UYN|1|01', '5UYP|1|01', '5UYQ|1|01', '5WDT|1|A', '5WE4|1|A', '5WE6|1|A', '5WF0|1|A', '5WFK|1|A', '5WFS|1|A', '3J9Z|1|LA', '3JA1|1|LA', '3JCJ|1|A', '6DNC|1|B', '5NP6|1|Y', '6H4N|1|A', '5H5U|1|A', '5MDV|1|1', '5MDW|1|1', '5MDY|1|1', '5MGP|1|A', '5U4I|1|A', '5U4J|1|A', '5U9F|1|01', '5U9G|1|01', '6ENF|1|A', '6ENJ|1|A', '6ENU|1|A', '6C4I|1|A', '3JBU|1|b', '3JBV|1|b', '5JTE|1|BA', '5JU8|1|BA', '5NWY|1|N', '5O2R|1|A', '5LZA|1|A', '5LZD|1|A', '5IQR|1|1', '5KPS|1|28', '5KPW|1|27', '5KPX|1|27', '5L3P|1|A', '4V85|1|BA', '4V89|1|BA', '4V6C|1|BA', '4V6C|1|DA', '4V6D|1|BA', '4V6D|1|DA', '4V6E|1|BA', '4V6E|1|DA', '3JCD|1|A', '6O9J|1|B', '6O9K|1|A', '6OFX|1|1', '6OGI|1|1', '6OGF|1|1', '6OG7|1|1', '6BY1|1|DA', '6BY1|1|CA', '6ORE|1|1', '6OSQ|1|1', '6ORL|1|1', '6OUO|1|1', '6OT3|1|1', '6OSK|1|1', '6Q97|1|1', '6Q9A|1|1', '6SZS|1|A', '6TC3|1|23S1', '5LZE|1|A']
    #test_chains = ['4V54|1|BB']
    failed_chains = []

    start_time = time.time()
    for chain in test_chains:
        try:
            run(chain)
        except:
            "Cannot produce data for the 23S chain " + chain
            failed_chains.append(chain)
            continue


    print "The run is complete"
    if failed_chains:
        print "Cannot produce data for the following chains:"
        for chain in failed_chains:
            print chain
    print "--- {} seconds ---".format(time.time() - start_time)





