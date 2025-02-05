import argparse
import logging
import sys

import pandas as pd
from molsig.Signature import MoleculeSignature
from rdkit import Chem


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-smiles-str", required=True, help="SMILES string")
    parser.add_argument("--output-data-tsv", required=True, help="Output file")
    args = parser.parse_args()

    # Init
    logging.basicConfig(
        level=logging.DEBUG,  # Set the minimum log level (DEBUG, INFO, WARNING, etc.)
        format="%(asctime)s - %(levelname)s - %(message)s",  # Log format
        handlers=[logging.StreamHandler()],  # Log to standard output
    )
    logger = logging.getLogger(__name__)

    # Build signature
    logging.info("Load SMILES into rdkit")
    smiles = args.input_smiles_str
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.error(f"SMILES is not validated by rdkit: {smiles}")
        sys.exit(1)

    logging.info("Build signature")
    mol_sig = MoleculeSignature(mol)

    data = dict(SMILES=args.input_smiles_str, signature=mol_sig.to_list())
    df = pd.DataFrame([data])
    logging.info(f"Save signature to: {args.output_data_tsv}")
    df.to_csv(args.output_data_tsv, sep="\t", index=False)
