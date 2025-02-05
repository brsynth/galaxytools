import argparse
import logging
import sys
import time

import pandas as pd
from molsig.enumerate_signature import enumerate_molecule_from_morgan
from molsig.SignatureAlphabet import load_alphabet
from rdkit import Chem
from rdkit.Chem import AllChem


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-smiles-str", required=True, help="SMILES string")
    parser.add_argument("--input-alphabet-npz", required=True, help="Alphabet file")
    parser.add_argument("--output-data-tsv", required=True, help="Output file")
    args = parser.parse_args()

    # Init
    logging.basicConfig(
        level=logging.DEBUG,  # Set the minimum log level (DEBUG, INFO, WARNING, etc.)
        format="%(asctime)s - %(levelname)s - %(message)s",  # Log format
        handlers=[logging.StreamHandler()],  # Log to standard output
    )
    logger = logging.getLogger(__name__)

    # Load Alphabet
    logging.info("Load alphabet")
    Alphabet = load_alphabet(args.input_alphabet_npz, verbose=True)

    # Create ECFP
    logging.info("Create ECFP")
    smiles = args.input_smiles_str
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error(f"SMILES is not validated by rdkit: {smiles}")
        sys.exit(1)
    fpgen = AllChem.GetMorganGenerator(
        radius=Alphabet.radius,
        fpSize=Alphabet.nBits,
        includeChirality=Alphabet.use_stereo,
    )
    morgan = fpgen.GetCountFingerprint(mol).ToList()

    logging.info("Enumerate molecules - start")
    start = time.time()
    Ssig, Smol, Nsig, thresholds_reached, computational_times = (
        enumerate_molecule_from_morgan(morgan, Alphabet)
    )
    end = round(time.time() - start, 2)
    logging.info("Enumerate molecules - end")
    logging.info(f"Time computed: {end}")

    sthresholds_reached = ", ".join([str(x) for x in thresholds_reached])
    scomputational_times = ",".join([str(x) for x in computational_times])

    logging.info(f"Thresholds_reached: {sthresholds_reached}")
    logging.info(f"Computational times: {scomputational_times}")

    df = pd.DataFrame(list(Smol), columns=["SMILES"])
    df.to_csv(args.output_data_tsv, sep="\t", index=False)
