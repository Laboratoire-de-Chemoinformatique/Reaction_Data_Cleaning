#############################################################################
# Chemical reactions data curation best practices
# including optimized RDTool
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

from CGRtools.files import RDFRead, RDFWrite, SDFWrite, SDFRead, SMILESRead
from CGRtools.containers import MoleculeContainer, ReactionContainer
import logging
from ordered_set import OrderedSet
from typing import Tuple
import os
import io
import pathlib
from pathlib import PurePosixPath


class Standardizer:
    def __init__(self, skip_errors=False, log_file=None, keep_unbalanced_ions=False, id_tag='Reaction_ID',
                 action_on_isotopes=0, keep_reagents=False, logger=None, ignore_mapping=False, jvm_path=None,
                 rdkit_dearomatization=False, remove_unchanged_parts=True, skip_tautomerize=False,
                 jchem_path=None) -> None:
        if logger is None:
            self.logger = self._config_log(log_file, logger_name='logger')
        else:
            self.logger = logger
        self._skip_errors = skip_errors
        self._keep_unbalanced_ions = keep_unbalanced_ions
        self._id_tag = id_tag
        self._action_on_isotopes = action_on_isotopes
        self._keep_reagents = keep_reagents
        self._ignore_mapping = ignore_mapping
        self._remove_unchanged_parts_flag = remove_unchanged_parts
        self._skip_tautomerize = skip_tautomerize
        self._dearomatize_by_rdkit = rdkit_dearomatization
        if not skip_tautomerize:
            if jvm_path:
                os.environ['JDK_HOME'] = jvm_path
                os.environ['JAVA_HOME'] = jvm_path
                os.environ['PATH'] += f';{PurePosixPath(jvm_path).joinpath("bin").joinpath("server")};' \
                                      f'{PurePosixPath(jvm_path).joinpath("bin").joinpath("server")};'
            if jchem_path:
                import jnius_config
                jnius_config.add_classpath(jchem_path)
            from jnius import autoclass
            Standardizer = autoclass('chemaxon.standardizer.Standardizer')
            self._Molecule = autoclass('chemaxon.struc.Molecule')
            self._MolHandler = autoclass('chemaxon.util.MolHandler')
            self._standardizer = Standardizer('tautomerize')

    def standardize_file(self, input_file=None) -> OrderedSet:
        """
        Standardize a set of reactions in a file. Returns an ordered set of ReactionContainer objects passed the
        standardization protocol.
        :param input_file: str
        :return: OrderedSet
        """
        if pathlib.Path(input_file).suffix == '.rdf':
            data = self._read_RDF(input_file)
        elif pathlib.Path(input_file).suffix == '.smi' or pathlib.Path(input_file).suffix == '.smiles':
            data = self._read_SMILES(input_file)
        else:
            raise ValueError('Data format is not recognized!')

        print("{0} reactions passed..".format(len(data)))
        return data

    def _read_RDF(self, input_file) -> OrderedSet:
        """
        Reads an RDF file. Returns an ordered set of ReactionContainer objects passed the standardization protocol.
        :param input_file: str
        :return: OrderedSet
        """
        data = OrderedSet()
        self.logger.info('Start..')
        with RDFRead(input_file, ignore=self._ignore_mapping, store_log=True, remap=self._ignore_mapping) as ifile, \
                open(input_file) as meta_searcher:
            for reaction in ifile._data:
                if isinstance(reaction, tuple):
                    meta_searcher.seek(reaction.position)
                    flag = False
                    for line in meta_searcher:
                        if flag and '$RFMT' in line:
                            self.logger.critical(f'Reaction id extraction problem rised for the reaction '
                                                 f'#{reaction.number + 1}: a reaction id was expected but $RFMT line '
                                                 f'was found!')
                        if flag:
                            self.logger.critical(f'Reaction {line.strip().split()[1]}: Parser has returned an error '
                                                 f'message\n{reaction.log}')
                            break
                        elif '$RFMT' in line:
                            self.logger.critical(f'Reaction #{reaction.number + 1} has no reaction id!')
                        elif f'$DTYPE {self._id_tag}' in line:
                            flag = True
                    continue
                standardized_reaction = self.standardize(reaction)
                if standardized_reaction:
                    if standardized_reaction not in data:
                        data.add(standardized_reaction)
                    else:
                        i = data.index(standardized_reaction)
                        if 'Extraction_IDs' not in data[i].meta:
                            data[i].meta['Extraction_IDs'] = ''
                        data[i].meta['Extraction_IDs'] = ','.join(data[i].meta['Extraction_IDs'].split(',') +
                                                                  [reaction.meta[self._id_tag]])
                        self.logger.info('Reaction {0} is a duplicate of the reaction {1}..'
                                         .format(reaction.meta[self._id_tag], data[i].meta[self._id_tag]))
        return data

    def _read_SMILES(self, input_file) -> OrderedSet:
        """
        Reads a SMILES file. Returns an ordered set of ReactionContainer objects passed the standardization protocol.
        :param input_file: str
        :return: OrderedSet
        """
        data = OrderedSet()
        self.logger.info('Start..')
        with SMILESRead(input_file, ignore=True, store_log=True, remap=self._ignore_mapping, header=True) as ifile, \
                open(input_file) as meta_searcher:
            id_tag_position = meta_searcher.readline().strip().split().index(self._id_tag)
            if id_tag_position is None or id_tag_position == 0:
                self.logger.critical(f'No reaction ID tag was found in the header!')
                raise ValueError(f'No reaction ID tag was found in the header!')
            for reaction in ifile._data:
                if isinstance(reaction, tuple):
                    meta_searcher.seek(reaction.position)
                    line = meta_searcher.readline().strip().split()
                    if len(line) <= id_tag_position:
                        self.logger.critical(f'No reaction ID tag was found in line {reaction.number}!')
                        raise ValueError(f'No reaction ID tag was found in line {reaction.number}!')
                    r_id = line[id_tag_position]
                    self.logger.critical(f'Reaction {r_id}: Parser has returned an error message\n{reaction.log}')
                    continue

                standardized_reaction = self.standardize(reaction)
                if standardized_reaction:
                    if standardized_reaction not in data:
                        data.add(standardized_reaction)
                    else:
                        i = data.index(standardized_reaction)
                        if 'Extraction_IDs' not in data[i].meta:
                            data[i].meta['Extraction_IDs'] = ''
                        data[i].meta['Extraction_IDs'] = ','.join(data[i].meta['Extraction_IDs'].split(',') +
                                                                  [reaction.meta[self._id_tag]])
                        self.logger.info('Reaction {0} is a duplicate of the reaction {1}..'
                                         .format(reaction.meta[self._id_tag], data[i].meta[self._id_tag]))
        return data

    def standardize(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Standardization protocol: transform functional groups, kekulize, remove explicit hydrogens,
        check for radicals (remove if something was found), check for isotopes, regroup ions (if the total charge
        of reactants and/or products is not zero, and the 'keep_unbalanced_ions' option is False which is by default,
        such reactions are removed; if the 'keep_unbalanced_ions' option is set True, they are kept), check valences
        (remove if something is wrong), aromatize (thiele method), fix mapping (for symmetric functional groups) if
        such is in, remove unchanged parts.
        :param reaction: ReactionContainer
        :return: ReactionContainer
        """
        self.logger.info('Reaction {0}..'.format(reaction.meta[self._id_tag]))
        try:
            reaction.standardize()
        except:
            self.logger.exception(
                'Reaction {0}: Cannot standardize functional groups..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception(
                    'Reaction {0}: Cannot standardize functional groups..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            reaction.kekule()
        except:
            self.logger.exception('Reaction {0}: Cannot kekulize..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot kekulize..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            if self._check_valence(reaction):
                self.logger.info(
                    'Reaction {0}: Bad valence: {1}'.format(reaction.meta[self._id_tag], reaction.meta['mistake']))
                return
        except:
            self.logger.exception('Reaction {0}: Cannot check valence..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                self.logger.critical('Stop the algorithm!')
                raise Exception('Reaction {0}: Cannot check valence..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            if not self._skip_tautomerize:
                reaction = self._tautomerize(reaction)
        except:
            self.logger.exception('Reaction {0}: Cannot tautomerize..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot tautomerize..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            reaction.implicify_hydrogens()
        except:
            self.logger.exception(
                'Reaction {0}: Cannot remove explicit hydrogens..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot remove explicit hydrogens..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            if self._check_radicals(reaction):
                self.logger.info('Reaction {0}: Radicals were found..'.format(reaction.meta[self._id_tag]))
                return
        except:
            self.logger.exception('Reaction {0}: Cannot check radicals..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot check radicals..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            if self._action_on_isotopes == 1 and self._check_isotopes(reaction):
                self.logger.info('Reaction {0}: Isotopes were found..'.format(reaction.meta[self._id_tag]))
                return
            elif self._action_on_isotopes == 2 and self._check_isotopes(reaction):
                reaction.clean_isotopes()
                self.logger.info('Reaction {0}: Isotopes were removed but the reaction was kept..'.format(
                    reaction.meta[self._id_tag]))
        except:
            self.logger.exception('Reaction {0}: Cannot check for isotopes..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot check for isotopes..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            reaction, return_code = self._split_ions(reaction)
            if return_code == 1:
                self.logger.info('Reaction {0}: Ions were split..'.format(reaction.meta[self._id_tag]))
            elif return_code == 2:
                self.logger.info('Reaction {0}: Ions were split but the reaction is imbalanced..'.format(
                    reaction.meta[self._id_tag]))
                if not self._keep_unbalanced_ions:
                    return
        except:
            self.logger.exception('Reaction {0}: Cannot group ions..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot group ions..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            reaction.thiele()
        except:
            self.logger.exception('Reaction {0}: Cannot aromatize..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot aromatize..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            reaction.fix_mapping()
        except:
            self.logger.exception('Reaction {0}: Cannot fix mapping..'.format(reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot fix mapping..'.format(reaction.meta[self._id_tag]))
            else:
                return
        try:
            if self._remove_unchanged_parts_flag:
                reaction = self._remove_unchanged_parts(reaction)
                if not reaction.reactants:
                    self.logger.info('Reaction {0}: Reactants are empty..'.format(reaction.meta[self._id_tag]))
                    return
                if not reaction.products:
                    self.logger.info('Reaction {0}: Products are empty..'.format(reaction.meta[self._id_tag]))
                    return
        except:
            self.logger.exception('Reaction {0}: Cannot remove unchanged parts or the reaction is empty..'.format(
                reaction.meta[self._id_tag]))
            if not self._skip_errors:
                raise Exception('Reaction {0}: Cannot remove unchanged parts or the reaction is empty..'.format(
                    reaction.meta[self._id_tag]))
            else:
                return
        self.logger.debug('Reaction {0} is done..'.format(reaction.meta[self._id_tag]))
        return reaction

    def write(self, output_file: str, data: OrderedSet) -> None:
        """
        Dump a set of reactions.
        :param data: OrderedSet
        :param output_file: str
        :return: None
        """
        with RDFWrite(output_file) as out:
            for r in data:
                out.write(r)

    def _check_valence(self, reaction: ReactionContainer) -> bool:
        """
        Checks valences.
        :param reaction: ReactionContainer
        :return: bool
        """
        mistakes = []
        for molecule in (reaction.reactants + reaction.products + reaction.reagents):
            valence_mistakes = molecule.check_valence()
            if valence_mistakes:
                mistakes.append(("|".join([str(num) for num in valence_mistakes]),
                                          "|".join([str(molecule.atom(n)) for n in valence_mistakes]), str(molecule)))
        if mistakes:
            message = ",".join([f'{atom_nums} at {atoms} in {smiles}' for atom_nums, atoms, smiles in mistakes])
            reaction.meta['mistake'] = f'Valence mistake: {message}'
            return True
        return False

    def _config_log(self, log_file: str, logger_name: str):
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(fmt='%(asctime)s: %(message)s', datefmt='%d/%m/%Y %H:%M:%S')
        logger.handlers.clear()
        fileHandler = logging.FileHandler(filename=log_file, mode='w')
        fileHandler.setFormatter(formatter)
        fileHandler.setLevel(logging.DEBUG)
        logger.addHandler(fileHandler)
        # logging.basicConfig(filename=log_file, level=logging.info, filemode='w', format='%(asctime)s: %(message)s',
        #                     datefmt='%d/%m/%Y %H:%M:%S')
        return logger

    def _check_radicals(self, reaction: ReactionContainer) -> bool:
        """
        Checks radicals.
        :param reaction: ReactionContainer
        :return: bool
        """
        for molecule in (reaction.reactants + reaction.products + reaction.reagents):
            for n, atom in molecule.atoms():
                if atom.is_radical:
                    return True
        return False

    def _calc_charge(self, molecule: MoleculeContainer) -> int:
        """Computing charge of molecule
        :param: molecule: MoleculeContainer
        :return: int
        """
        return sum(molecule._charges.values())

    def _group_ions(self, reaction: ReactionContainer):
        """
        Ungroup molecules recorded as ions, regroup ions. Returns a tuple with the corresponding ReactionContainer and
        return code as int (0 - nothing was changed, 1 - ions were regrouped, 2 - ions are unbalanced).
        :param reaction: current reaction
        :return: tuple[ReactionContainer, int]
        """
        meta = reaction.meta
        reaction_parts = []
        return_codes = []
        for molecules in (reaction.reactants, reaction.reagents, reaction.products):
            divided_molecules = [x for m in molecules for x in m.split('.')]

            if len(divided_molecules) == 0:
                reaction_parts.append(())
                continue
            elif len(divided_molecules) == 1 and self._calc_charge(divided_molecules[0]) == 0:
                return_codes.append(0)
                reaction_parts.append(molecules)
                continue
            elif len(divided_molecules) == 1:
                return_codes.append(2)
                reaction_parts.append(molecules)
                continue

            new_molecules = []
            cations, anions, ions = [], [], []
            total_charge = 0
            for molecule in divided_molecules:
                mol_charge = self._calc_charge(molecule)
                total_charge += mol_charge
                if mol_charge == 0:
                    new_molecules.append(molecule)
                elif mol_charge > 0:
                    cations.append((mol_charge, molecule))
                    ions.append((mol_charge, molecule))
                else:
                    anions.append((mol_charge, molecule))
                    ions.append((mol_charge, molecule))

            if len(cations) == 0 and len(anions) == 0:
                return_codes.append(0)
                reaction_parts.append(tuple(new_molecules))
                continue
            elif total_charge != 0:
                return_codes.append(2)
                reaction_parts.append(tuple(divided_molecules))
                continue
            else:
                salt = MoleculeContainer()
                for ion_charge, ion in ions:
                    salt = salt.union(ion)
                    total_charge += ion_charge
                    if total_charge == 0:
                        new_molecules.append(salt)
                        salt = MoleculeContainer()
                if total_charge != 0:
                    new_molecules.append(salt)
                    return_codes.append(2)
                    reaction_parts.append(tuple(new_molecules))
                else:
                    return_codes.append(1)
                    reaction_parts.append(tuple(new_molecules))
        return ReactionContainer(reactants=reaction_parts[0], reagents=reaction_parts[1], products=reaction_parts[2],
                                 meta=meta), max(return_codes)

    def _split_ions(self, reaction: ReactionContainer):
        """
        Split ions in a reaction. Returns a tuple with the corresponding ReactionContainer and
        a return code as int (0 - nothing was changed, 1 - ions were split, 2 - ions were split but the reaction
        is imbalanced).
        :param reaction: current reaction
        :return: tuple[ReactionContainer, int]
        """
        meta = reaction.meta
        reaction_parts = []
        return_codes = []
        for molecules in (reaction.reactants, reaction.reagents, reaction.products):
            divided_molecules = [x for m in molecules for x in m.split('.')]

            total_charge = 0
            ions_present = False
            for molecule in divided_molecules:
                mol_charge = self._calc_charge(molecule)
                total_charge += mol_charge
                if mol_charge != 0:
                    ions_present = True

            if ions_present and total_charge:
                return_codes.append(2)
            elif ions_present:
                return_codes.append(1)
            else:
                return_codes.append(0)

            reaction_parts.append(tuple(divided_molecules))

        return ReactionContainer(reactants=reaction_parts[0], reagents=reaction_parts[1], products=reaction_parts[2],
                                 meta=meta), max(return_codes)

    def _remove_unchanged_parts(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Ungroup molecules, remove unchanged parts from reactants and products.
        :param reaction: current reaction
        :return: ReactionContainer
        """
        meta = reaction.meta
        reactants = [m for m in reaction.reactants]
        new_reactants, new_reagents, new_products = [m for m in reaction.reactants], [m for m in reaction.reagents], \
                                                    [m for m in reaction.products]
        for reactant in reactants:
            if reactant in new_products:
                # if self._ignore_mapping:
                new_reagents.append(reactant)
                new_reactants.remove(reactant)
                new_products.remove(reactant)
                # elif self._confirm_equivalence_by_mapping(reactant, new_products):
                #     new_reagents.append(reactant)
                #     new_reactants = [m for m in new_reactants if
                #                      not self._confirm_equivalence_by_mapping(reactant, tuple([m]))]
                #     new_products = [m for m in new_products if
                #                     not self._confirm_equivalence_by_mapping(reactant, tuple([m]))]
        if not self._keep_reagents:
            new_reagents = []
        return ReactionContainer(reactants=tuple(new_reactants), reagents=tuple(new_reagents),
                                 products=tuple(new_products), meta=meta)

    def _check_isotopes(self, reaction: ReactionContainer) -> bool:
        for molecules in (reaction.reactants, reaction.products):
            for molecule in molecules:
                for _, atom in molecule.atoms():
                    if atom.isotope:
                        return True
        return False

    def _confirm_equivalence_by_mapping(self, molecule: MoleculeContainer,
                                        molecules_list: Tuple[ReactionContainer, ...]) -> bool:
        """
        Checks if the molecule really corresponds to one of the molecules in the molecules_list
        taking into account their mapping if such exists.
        :param molecule: molecule for search, molecules_list: where to search
        :return: bool
        """
        analogs = [retreived_molecule for retreived_molecule in molecules_list if retreived_molecule == molecule]
        reference_mapping = molecule.atoms_numbers
        for analog in analogs:
            analog_mapping = analog.atoms_numbers
            if len(set(reference_mapping) & set(analog_mapping)) == len(reference_mapping):
                return True
        return False

    def _tautomerize(self, reaction: ReactionContainer) -> ReactionContainer:
        """
        Perform ChemAxon tautomerization.
        :param reaction: reaction that needs to be tautomerized
        :return: ReactionContainer
        """
        new_molecules = []
        for part in [reaction.reactants, reaction.reagents, reaction.products]:
            tmp = []
            for mol in part:
                with io.StringIO() as f, SDFWrite(f) as i:
                    i.write(mol)
                    sdf = f.getvalue()
                mol_handler = self._MolHandler(sdf)
                mol_handler.clean(True, '2')
                molecule = mol_handler.getMolecule()
                self._standardizer.standardize(molecule)
                new_mol_handler = self._MolHandler(molecule)
                new_sdf = new_mol_handler.toFormat('SDF')
                with io.StringIO('\n  ' + new_sdf.strip()) as f, SDFRead(f, remap=False) as i:
                    new_mol = next(i)
                tmp.append(new_mol)
            new_molecules.append(tmp)
        return ReactionContainer(reactants=tuple(new_molecules[0]), reagents=tuple(new_molecules[1]),
                                 products=tuple(new_molecules[2]), meta=reaction.meta)

    # def _dearomatize_by_RDKit(self, reaction: ReactionContainer) -> ReactionContainer:
    #     """
    #     Dearomatizes by RDKit (needs in case of some mappers, such as RXNMapper).
    #     :param reaction: ReactionContainer
    #     :return: ReactionContainer
    #     """
    #     with io.StringIO() as f, RDFWrite(f) as i:
    #         i.write(reaction)
    #         s = '\n'.join(f.getvalue().split('\n')[3:])
    #         rxn = rdChemReactions.ReactionFromRxnBlock(s)
    #         reactants, reagents, products = [], [], []
    #         for mol in rxn.GetReactants():
    #             try:
    #                 Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_KEKULIZE, catchErrors=True)
    #             except Chem.rdchem.KekulizeException:
    #                 return reaction
    #             with io.StringIO(Chem.MolToMolBlock(mol)) as f2, SDFRead(f2, remap=False) as sdf_i:
    #                 reactants.append(next(sdf_i))
    #         for mol in rxn.GetAgents():
    #             try:
    #                 Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_KEKULIZE, catchErrors=True)
    #             except Chem.rdchem.KekulizeException:
    #                 return reaction
    #             with io.StringIO(Chem.MolToMolBlock(mol)) as f2, SDFRead(f2, remap=False) as sdf_i:
    #                 reagents.append(next(sdf_i))
    #         for mol in rxn.GetProducts():
    #             try:
    #                 Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_KEKULIZE, catchErrors=True)
    #             except Chem.rdchem.KekulizeException:
    #                 return reaction
    #             with io.StringIO(Chem.MolToMolBlock(mol)) as f2, SDFRead(f2, remap=False) as sdf_i:
    #                 products.append(next(sdf_i))
    #
    #     new_reaction = ReactionContainer(reactants=tuple(reactants), reagents=tuple(reagents), products=tuple(products),
    #                                      meta=reaction.meta)
    #
    #     return new_reaction


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="This is a tool for reaction standardization.",
                                     epilog="Arkadii Lin, Strasbourg/Kazan 2020", prog="Standardizer")
    parser.add_argument("-i", "--input", type=str, help="Input RDF file.")
    parser.add_argument("-o", "--output", type=str, help="Output RDF file.")
    parser.add_argument("-id", "--idTag", default='Reaction_ID', type=str, help="ID tag in the RDF file.")
    parser.add_argument("--skipErrors", action="store_true", help="Skip errors.")
    parser.add_argument("--keep_unbalanced_ions", action="store_true", help="Will keep reactions with unbalanced ions.")
    parser.add_argument("--action_on_isotopes", type=int, default=0, help="Action performed if an isotope is "
                                                                          "found: 0 - to ignore isotopes; "
                                                                          "1 - to remove reactions with isotopes; "
                                                                          "2 - to clear isotopes' labels.")
    parser.add_argument("--keep_reagents", action="store_true", help="Will keep reagents from the reaction.")
    parser.add_argument("--ignore_mapping", action="store_true", help="Will ignore the initial mapping in the file.")
    parser.add_argument("--keep_unchanged_parts", action="store_true", help="Will keep unchanged parts in a reaction.")
    parser.add_argument("--logFile", type=str, default='logFile.txt', help="Log file name.")
    parser.add_argument("--skip_tautomerize", action="store_true", help="Will skip generation of the major tautomer.")
    parser.add_argument("--rdkit_dearomatization", action="store_true", help="Will kekulize the reaction using RDKit "
                                                                             "facilities.")
    parser.add_argument("--jvm_path", type=str,
                        help="JVM path (e.g. C:\\Program Files\\Java\\jdk-13.0.2).")
    parser.add_argument("--jchem_path", type=str, help="JChem path (e.g. C:\\Users\\user\\JChemSuite\\lib\\jchem.jar).")
    args = parser.parse_args()

    standardizer = Standardizer(skip_errors=args.skipErrors, log_file=args.logFile,
                                keep_unbalanced_ions=args.keep_unbalanced_ions, id_tag=args.idTag,
                                action_on_isotopes=args.action_on_isotopes, keep_reagents=args.keep_reagents,
                                ignore_mapping=args.ignore_mapping, skip_tautomerize=args.skip_tautomerize,
                                remove_unchanged_parts=(not args.keep_unchanged_parts), jvm_path=args.jvm_path,
                                jchem_path=args.jchem_path, rdkit_dearomatization=args.rdkit_dearomatization)
    data = standardizer.standardize_file(input_file=args.input)
    standardizer.write(output_file=args.output, data=data)
