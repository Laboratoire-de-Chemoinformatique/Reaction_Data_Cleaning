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
# Sources: [1] ReactionDecoder (or RDTool)
#              author: Syed Asad Rahman https://github.com/asad
#              repository: https://github.com/asad/ReactionDecoder
#############################################################################

from multiprocessing import Pool
from scripts.standardizer import Standardizer


def __run(i):
    standardizer = Standardizer(skip_errors=True, keep_unbalanced_ions=False, id_tag='Reaction_ID',
                                log_file='./Log_' + str(i) + '.txt', keep_reagents=False, ignore_mapping=True,
                                action_on_isotopes=2)
    data = standardizer.standardize_file(input_file='./chunk_' + str(i) + '.smi', )
    standardizer.write(output_file='./clean_' + str(i) + '.rdf',
                       data=data)
    standardizer.logger.handlers.pop().close()


if __name__ == '__main__':
    import jnius_config

    jnius_config.add_classpath("/home/opt/chemaxon/jchemsuite/lib/jchem.jar")
    n_threads = 45
    p = Pool(processes=n_threads)  # Number of threads that can be used in parallel
    n_chunks = 14843
    p.map(__run, [i for i in range(0, n_chunks)])  # Enumerate the chunks
