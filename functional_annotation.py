from collections import OrderedDict
import process_annotation as rcs
import csv
import time

# This doesn't need to be ordered. But I might want to print this as a csv, so it's good to keep the columns ordered
ribosome_components = OrderedDict()
ribosome_components['SSU_16S'] = None
ribosome_components['LSU_23S'] = None
ribosome_components['LSU_5S'] = None
ribosome_components['mRNA'] = None
ribosome_components['aminoacyl-trna'] = None
ribosome_components['peptidyl-trna'] = None
ribosome_components['exit-trna'] = None

# This is the reference dictionary that stores nts that make contacts between specific pair of ribosomal RNA chains
reference_units = {
    'ssu_nts_lsu': ['5J7L|1|AA|A|1408', '5J7L|1|AA|A|1418', '5J7L|1|AA|A|1483'],
    'ssu_nts_mrna': ['5J7L|1|AA|G|926', '5J7L|1|AA|4OC|1402', '5J7L|1|AA|C|1403'],
    'ssu_nts_ptrna': ['5J7L|1|AA|G|1338', '5J7L|1|AA|A|1339', '5J7L|1|AA|C|1400'],
    'ssu_nts_atrna': ['5J7L|1|AA|G|530', '5J7L|1|AA|A|1492', '5J7L|1|AA|A|1493'],
    'ssu_nts_etrna': ['5J7L|1|AA|G|693', '5J7L|1|AA|A|694'],
    'lsu_nts_atrna': ['5J7L|1|DA|G|2553', '5J7L|1|DA|G|2583', '5J7L|1|DA|U|2585'],
    'lsu_nts_etrna': ['5J7L|1|DA|G|2112', '5J7L|1|DA|G|2421', '5J7L|1|DA|C|2422'],
    'lsu_nts_ptrna': ['5J7L|1|DA|OMG|2251', '5J7L|1|DA|G|2252', '5J7L|1|DA|G|2253'],
    'lsu_nts_5S': ['5J7L|1|DA|A|861', '5J7L|1|DA|G|862', '5J7L|1|DA|A|917']
}


def run(test_ife):

    ribosome_components['SSU_16S'] = test_ife

    # test_ssu is a tuple of (pdb, chain) of the ssu
    # components_ife is a list of rna_chains in the pdbid that is associated with the ssu
    test_ssu, components_ife = rcs.get_components_ife(test_ife)

    # Returns a list of unit_ids that should interact with the 23S in the test_ife based on the reference_units
    ssu_nts_lsu_corr = rcs.get_units_correspondence(reference_units['ssu_nts_lsu'], test_ssu)

    # Infer 23S ife
    ribosome_components['LSU_23S'], components_ife = rcs.infer_interacting_ife(ssu_nts_lsu_corr,
                                                                                components_ife)

    # Get the pdb & chain information from the LSU ife
    lsu_pdb, _, lsu_chain = ribosome_components['LSU_23S'].split('|')

    # Build a tuple of (pdb, chain)
    test_lsu = (lsu_pdb, lsu_chain)

    # Get correspondence for all the reference nts involved in chain-chain pairing
    test_units_corr = rcs.units_correspondence(reference_units, test_ssu, test_lsu)

    # Infer mRNA ife
    ribosome_components['mRNA'], components_ife = rcs.infer_interacting_ife(test_units_corr['ssu_mrna'],
                                                                            components_ife)

    # Infer 5S ife
    ribosome_components['LSU_5S'], components_ife = rcs.infer_interacting_ife(test_units_corr['lsu_5S'],
                                                                            components_ife)

    # Infer P-tRNA ife
    ribosome_components['peptidyl-trna'], components_ife = rcs.infer_interacting_ife(test_units_corr['ssu_ptrna'],
                                                                                    components_ife)
    # Infer A-tRNA ife
    ribosome_components['aminoacyl-trna'], components_ife = rcs.infer_interacting_ife(test_units_corr['ssu_atrna'],
                                                                                    components_ife,
                                                                                    ribosome_components,
                                                                                    'atrna')
    # Infer E-tRNA ife
    ribosome_components['exit-trna'], components_ife = rcs.infer_interacting_ife(test_units_corr['lsu_etrna'],
                                                                                components_ife,
                                                                                ribosome_components, 'etrna')
    
    # Infer A-tRNA state
    atrna_state = rcs.infer_tRNA_state(ribosome_components, 'aminoacyl-trna', test_units_corr)

    # Infer P-tRNA state
    ptrna_state = rcs.infer_tRNA_state(ribosome_components, 'peptidyl-trna', test_units_corr)

    # Infer E-tRNA state
    etrna_state = rcs.infer_tRNA_state(ribosome_components, 'exit-trna', test_units_corr)

    #trna_states = (atrna_state, ptrna_state, etrna_state)

    ribosome_components['atrna_state'] = atrna_state
    ribosome_components['ptrna_state'] = ptrna_state
    ribosome_components['etrna_state'] = etrna_state
    

    with open('ribosome_chain_annotation_new.csv', 'a') as f:  
        w = csv.DictWriter(f, ribosome_components.keys())
        w.writerow(ribosome_components)
    

if __name__ == "__main__":
    test_chains = ['5UYM|1|A', '5WDT|1|a', '3JBV|1|A', '3JCJ|1|g', '4V4Q|1|AA', '4V50|1|AA']
    test_chains = ['4V4Q|1|AA', '4V4Q|1|CA', '4V50|1|AA', '4V50|1|CA', '4V5B|1|BA', '4V5B|1|DA', '4V9D|1|AA', '4V9D|1|BA', '4V9O|1|BA', '4V9O|1|DA', '4V9O|1|FA', '4V9O|1|HA', '4V9P|1|BA', '4V9P|1|DA', '4V9P|1|FA', '4V9P|1|HA', '4YBB|1|AA', '4YBB|1|BA', '5IT8|1|AA', '5IT8|1|BA', '5J5B|1|AA', '5J5B|1|BA', '5J7L|1|AA', '5J7L|1|BA', '5J88|1|AA', '5J88|1|BA', '5J8A|1|AA', '5J8A|1|BA', '5J91|1|AA', '5J91|1|BA', '5JC9|1|AA', '5JC9|1|BA', '5MDZ|1|2', '6BU8|1|A', '6GWT|1|a', '6GXM|1|a', '6GXN|1|a', '6GXO|1|a', '6I7V|1|AA', '6I7V|1|BA', '4U1U|1|AA', '4U1U|1|CA', '4U1V|1|AA', '4U1V|1|CA', '4U20|1|AA', '4U20|1|CA', '4U24|1|AA', '4U24|1|CA', '4U25|1|AA', '4U25|1|CA', '4U26|1|AA', '4U26|1|CA', '4U27|1|AA', '4U27|1|CA', '4V4H|1|AA', '4V4H|1|CA', '4V52|1|AA', '4V52|1|CA', '4V53|1|AA', '4V53|1|CA', '4V54|1|AA', '4V54|1|CA', '4V55|1|AA', '4V55|1|CA', '4V56|1|AA', '4V56|1|CA', '4V57|1|AA', '4V57|1|CA', '4V64|1|AA', '4V64|1|CA', '4V7S|1|AA', '4V7S|1|CA', '4V7T|1|AA', '4V7T|1|CA', '4V7U|1|AA', '4V7U|1|CA', '4V7V|1|AA', '4V7V|1|CA', '4V9C|1|AA', '4V9C|1|CA', '4WF1|1|AA', '4WF1|1|CA', '4WOI|1|AA', '4WOI|1|DA', '4WWW|1|QA', '4WWW|1|XA', '5KCR|1|1a', '5KCS|1|1a', '3J9Y|1|a', '3JCD|1|a', '3JCE|1|a', '5AFI|1|a', '5UYK|1|A', '5UYL|1|A', '5UYM|1|A', '5UYN|1|A', '5UYP|1|A', '5UYQ|1|A', '5WDT|1|a', '5WE4|1|a', '5WE6|1|a', '5WF0|1|a', '5WFK|1|a', '5WFS|1|a', '3J9Z|1|SA', '3JA1|1|SA', '3JCJ|1|g', '6DNC|1|A', '5NP6|1|D', '6H4N|1|a', '5H5U|1|h', '5MDV|1|2', '5MDW|1|2', '5MDY|1|2', '5MGP|1|a', '5U4I|1|a', '5U4J|1|a', '5U9F|1|A', '5U9G|1|A', '6ENF|1|a', '6ENJ|1|a', '6ENU|1|a', '6C4I|1|a', '3JBU|1|A', '3JBV|1|A', '5JTE|1|AA', '5JU8|1|AA', '5NWY|1|0', '5O2R|1|a', '5LZA|1|a', '5LZD|1|a', '5IQR|1|2', '5KPS|1|27', '5KPW|1|26', '5KPX|1|26', '5L3P|1|a', '4V85|1|AA', '4V89|1|AA', '4V6C|1|AA', '4V6C|1|CA', '4V6D|1|AA', '4V6D|1|CA', '4V6E|1|AA', '4V6E|1|CA']
    #test_chains = ['5LZE|1|a']
    failed_chains = []

    start_time = time.time()
    for chain in test_chains:
        try:
            run(chain)
        except:
            "Cannot produce data for the 16S chain " + chain
            failed_chains.append(chain)
            continue


    print "The run is complete"
    if failed_chains:
        print "Cannot produce data for the following chains:"
        for chain in failed_chains:
            print chain
    print "--- {} seconds ---".format(time.time() - start_time)


"""
# Infer A-tRNA state
atrna_state = rcs.infer_tRNA_state(ribosome_components, 'aminoacyl-trna', test_units_corr)

# Infer P-tRNA state
ptrna_state = rcs.infer_tRNA_state(ribosome_components, 'peptidyl-trna', test_units_corr)

# Infer E-tRNA state
etrna_state = rcs.infer_tRNA_state(ribosome_components, 'exit-trna', test_units_corr)

trna_states = (atrna_state, ptrna_state, etrna_state)
"""




