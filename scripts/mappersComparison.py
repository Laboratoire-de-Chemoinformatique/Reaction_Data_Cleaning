#############################################################################
# Chemical reactions data curation best practices
#############################################################################
# GNU LGPL https://www.gnu.org/licenses/lgpl-3.0.en.html
#############################################################################
# Corresponding Author: Timur Madzhidov and Alexandre Varnek
# Corresponding Author email: tmadzhidov@gmail.com and varnek@unistra.fr
# Main contributors: Arkadii Lin, Natalia Duybankova, Ramil Nugmanov, Rail Suleymanov and Timur Madzhidov
# Copyright: Copyright 2020,
#            MaDeSmart, Machine Design of Small Molecules by AI
#            VLAIO project HBC.2018.2287
# Credits: Kazan Federal University, Russia
#          University of Strasbourg, France
#          University of Linz, Austria
#          University of Leuven, Belgium
#          Janssen Pharmaceutica N.V., Beerse, Belgium
#          Rail Suleymanov, Arcadia, St. Petersburg, Russia
# License: GNU LGPL https://www.gnu.org/licenses/lgpl-3.0.en.html
# Version: 00.01
#############################################################################


from typing import Dict, Tuple
from CGRtools import SMILESRead
from CGRtools.files import RDFRead
import logging
import pickle


def parse_reactions(input_file: str, id_tag: str, do_basic_standardization: bool = True) -> Dict:
    """
    Reads an RDF/SMILES file. Returns a dictionary where the key is the reaction ID, and the value is the
    ReactionContainer.
    :param do_basic_standardization: bool
    :param input_file: str
    :param id_tag: str
    :return: Dict
    """
    data = {}
    if input_file.endswith('.rdf'):
        with RDFRead(input_file) as f:
            while True:
                try:
                    r = next(f)
                    key = r.meta[id_tag]
                    if do_basic_standardization:
                        r.thiele()
                        r.standardize()
                    data[key] = r
                    logging.info(f'Reaction {key} passed..')
                except StopIteration:
                    break
    elif input_file.endswith('.smi') or input_file.endswith('.smiles'):
        with SMILESRead(input_file, ignore=True, store_log=True, remap=False, header=True) as ifile, \
                open(input_file) as meta_searcher:
            id_tag_position = meta_searcher.readline().strip().split().index(id_tag)
            if id_tag_position is None or id_tag_position == 0:
                logging.critical(f'No reaction ID tag was found in the header!')
                raise ValueError(f'No reaction ID tag was found in the header!')
            for reaction in ifile._data:
                if isinstance(reaction, tuple):
                    meta_searcher.seek(reaction.position)
                    line = meta_searcher.readline().strip().split()
                    if len(line) <= id_tag_position:
                        logging.critical(f'No reaction ID tag was found in line {reaction.number}!')
                        raise ValueError(f'No reaction ID tag was found in line {reaction.number}!')
                    r_id = line[id_tag_position]
                    logging.critical(f'Reaction {r_id}: Parser has returned an error message\n{reaction.log}')
                    continue
                if do_basic_standardization:
                    reaction.thiele()
                    reaction.standardize()
                key = reaction.meta[id_tag]
                data[key] = reaction
                logging.info(f'Reaction {key} passed..')
    return data


def filter_data(referenced_aam: dict, generated_aam: dict) -> Tuple[Dict, Dict]:
    """
    To prevent any errors, only common reactions are selected, and others are removed.
    :param referenced_aam: dict
    :param generated_aam: dict
    :return: Tuple[Dict, Dict]
    """
    common_IDs = set(referenced_aam.keys()) & set(generated_aam.keys())
    new_ref, new_gen = {}, {}
    for key in common_IDs:
        new_ref[key] = referenced_aam[key]
        new_gen[key] = generated_aam[key]
    return new_ref, new_gen


def make_cgr(data: dict) -> Dict:
    """
    Transforms the reactions to CGRs.
    :param data: dict
    :return: Dict
    """
    result = {}
    for key, r in data.items():
        cgr = ~ r
        result[key] = cgr

    return result


def __config_log(log_file):
    logging.basicConfig(filename=log_file, level=logging.DEBUG, filemode='w', format='%(asctime)s: %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')


def main(reference_file: str, generated_file: str, log_file: str, id_tag: str, archive_file: str,
         ignore_basic_stand: bool):
    """
    Computes statistics on the comparison of the referenced and generated mappings.
    :param ignore_basic_stand: bool
    :param reference_file: str
    :param generated_file: str
    :param log_file: str
    :param archive_file: str
    :param id_tag: str
    """
    __config_log(log_file=log_file)

    print('Loading the reference dataset..')
    logging.info('Loading the reference dataset..')
    ref_mapping = parse_reactions(input_file=reference_file, id_tag=id_tag,
                                  do_basic_standardization=(not ignore_basic_stand))
    print(f'{len(ref_mapping)} reactions were found..')

    print("Loading the generated dataset..")
    logging.info("Loading the generated dataset..")
    gen_mapping = parse_reactions(input_file=generated_file, id_tag=id_tag,
                                  do_basic_standardization=(not ignore_basic_stand))
    print(f'{len(gen_mapping)} reactions were found..')

    print("Filtering data..")
    ref_mapping, gen_mapping = filter_data(referenced_aam=ref_mapping, generated_aam=gen_mapping)
    print(f'{len(ref_mapping)} reactions stayed..')

    print('Convert reference data to CGR..')
    ref_cgr = make_cgr(data=ref_mapping)

    print('Convert generated data to CGR..')
    gen_cgr = make_cgr(data=gen_mapping)

    statistics = dict(same_mapping=0, differentCD=0, differentCGRs=0, not_equal_reactions=0)
    reactions = dict(good_mapping=[], bad_mapping_diff_CD=[], bad_mapping_same_CD_but_diff_CGRs=[],
                     not_equal_reactions=[])

    for key in sorted(ref_cgr.keys()):
        if ref_mapping[key] != gen_mapping[key]:
            statistics['not_equal_reactions'] += 1
            reactions['not_equal_reactions'].append(key)
        elif ref_cgr[key] != gen_cgr[key]:
            chemical_distance = abs(len(ref_cgr[key].center_atoms) + len(ref_cgr[key].center_bonds)
                                    - len(gen_cgr[key].center_atoms) - len(gen_cgr[key].center_bonds))
            if chemical_distance > 0:
                statistics['differentCD'] += 1
                reactions['bad_mapping_diff_CD'].append(key)
            else:
                statistics['differentCGRs'] += 1
                reactions['bad_mapping_same_CD_but_diff_CGRs'].append(key)
        else:
            statistics['same_mapping'] += 1
            reactions['good_mapping'].append(key)

    pickle.dump(reactions, open(archive_file, 'wb'))

    print('----------------------------------------------------------------------')
    print('**********************************************************************')
    print(f'Reactions mapped identically: {statistics["same_mapping"]}')
    print(f'Reactions mapped differently (CD>0): {statistics["differentCD"]}')
    print(f'Reactions mapped differently (CD=0 but CGRs are different): {statistics["differentCGRs"]}')
    print(f'Not equal reactions: {statistics["not_equal_reactions"]}')
    print('**********************************************************************')
    print('----------------------------------------------------------------------')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="This tool is done to compare AAMs.",
                                     epilog="Arkadii Lin, Strasbourg/Kazan 2020", prog="Standardizer")
    parser.add_argument("-rd", "--reference_data", type=str, required=True,
                        help="Input RDF/SMILES file containing the referenced AAM.")
    parser.add_argument("-md", "--mapped_data", type=str, required=True,
                        help="Input RDF/SMILES file containing the generated AAM.")
    parser.add_argument("-id", "--id_tag", type=str, required=True,
                        help="The property field that stores the reaction ID.")
    parser.add_argument("--ignore_basic_stand", action="store_true", help="Will ignore basic standardization.")
    parser.add_argument("--ids_archive", type=str, help="A pickle archive containing a dictionary of four lists: "
                                                        "good_mapping - reactions' IDs that were mapped in the same "
                                                        "way; bad_mapping_diff_CD - reactions that were mapped "
                                                        "differently with the chemical distance above zero; "
                                                        "bad_mapping_same_CD_but_diff_CGRs - reactions that were "
                                                        "mapped differently with the same chemical distance; "
                                                        "not_equal_reactions - reactions that possess the same ID "
                                                        "but they are chemically different.")
    parser.add_argument("--logFile", type=str, default='comparison.log', help="Log file name.")
    args = parser.parse_args()

    main(reference_file=args.reference_data, generated_file=args.mapped_data, log_file=args.logFile,
         id_tag=args.id_tag, archive_file=args.ids_archive, ignore_basic_stand=args.ignore_basic_stand)
