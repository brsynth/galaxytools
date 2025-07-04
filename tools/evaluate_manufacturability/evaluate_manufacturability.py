import argparse
import os
from Bio import SeqIO
import dnacauldron
import dnachisel
import dnachisel.reports.constraints_reports as cr


def evaluate_manufacturability(files_to_evaluate, file_name_mapping, output_tsv, output_pdf, outdir_gb, use_file_names_as_id,
                               avoid_patterns, hairpin_constraints, gc_constraints, kmer_size):
    """Evaluate manufacturability of DNA sequences based on user-defined constraints."""

    files_to_evaluate = files_to_evaluate.split(',')

 #   try:
 #       if pdf_evaluation:
 #           print(f'PDF evaluation file path: {pdf_evaluation}')
 #       else:
 #          print('PDF evaluation file path is empty.')

 #       if excel_evaluation:
 #           print(f'Excel evaluation file path: {excel_evaluation}')
 #       else:
 #           print('Excel evaluation file path is empty.')
 #   except Exception as e:
 #       print(f'An error occurred while checking file paths: {e}')

    records_to_evaluate = dnacauldron.biotools.load_records_from_files(
        files=files_to_evaluate,
        folder=None,
        use_file_names_as_ids=use_file_names_as_id
    )
 #   #try:
 #       #if not records_to_evaluate:
 #           #print('records to evaluate: is empty')
 #       #else:
 #           #for record in records_to_evaluate:
 #               #print(f'records to evaluate: {record}')
 #   #except Exception as e:
 #       #print(f'An error occurred: {e}')

    length_cutoff = 100
    for part in records_to_evaluate:
        if len(part) < length_cutoff:
            raise ValueError(f"Error: Sequence '{part.id}' is too short (length {len(part)}). Minimum required: {length_cutoff}")

    constraint_list = []

    for pattern in avoid_patterns:
        constraint_list.append(dnachisel.AvoidPattern(pattern))
        print(f"AvoidPattern constraint: {pattern}")

    for constraint in hairpin_constraints:
        try:
            stem_size, hairpin_window = map(int, constraint.split(';'))
            constraint_list.append(dnachisel.AvoidHairpins(stem_size=stem_size, hairpin_window=hairpin_window))
            print(f"Hairpin constraint: stem_size={stem_size}, hairpin_window={hairpin_window}")
        except ValueError:
            print(f"Skipping invalid hairpin constraint: {constraint}")

    for constraint in gc_constraints:
        try:
            mini, maxi, window = map(float, constraint.split(';'))
            constraint_list.append(dnachisel.EnforceGCContent(mini=mini, maxi=maxi, window=int(window)))
            print(f"GC constraint: mini={mini}, maxi={maxi}, window={window}")
        except ValueError:
            print(f"Skipping invalid GC constraint: {constraint}")

    for k_size in kmer_size:
        try:
            constraint_list.append(dnachisel.UniquifyAllKmers(k=int(k_size)))
            print(f"k-mer size is: {k_size}")
        except ValueError:
            print(f"Skipping invalid k-mer size: {k_size}")

    dataframe = cr.constraints_breaches_dataframe(constraint_list, records_to_evaluate)
    if isinstance(file_name_mapping, str):
        file_name_mapping = dict(
            item.split(":") for item in file_name_mapping.split(",")
        )
    dataset_to_real_name = {
        key.split("/")[-1].replace(".dat", ""): value.replace(".gb", "")
        for key, value in file_name_mapping.items()
    }
    # rename sequences with file name mapped
    if "sequence" in dataframe.columns:
        sequences = dataframe["sequence"]
        dataframe[sequences] = dataframe[sequences].map(dataset_to_real_name)
    else:
        dataframe.index = dataframe.index.map(dataset_to_real_name)
 #   try:
 #       if dataframe.empty:
 #           print('dataframe: is empty')
 #       else:
 #           print(f'dataframe is:\n{dataframe}')
 #           print(dataframe.columns)
 #   except Exception as e:9x2mer
 #       print(f'An error occurred: {e}')
    dataframe.to_csv(output_tsv, sep='\t')
 #   try:
 #       if os.path.exists(output_tsv):
 #           if os.path.getsize(output_tsv) == 0:
 #               print('excel_evaluation file is EMPTY after writing.')
 #           else:
 #               print(f'excel_evaluation file has been created successfully: {output_tsv}')

 #               read_back = pd.read_csv(output_tsv, sep='\t')
 #               print("Contents of the excel_evaluation file:")
 #               print(read_back)
 #       else:
 #           print('excel_evaluation file does NOT exist after writing.')
 #   except Exception as e:
 #       print(f'An error occurred while saving or reading the file: {e}')
    records_annotated = cr.records_from_breaches_dataframe(dataframe, records_to_evaluate)
 #  rename records id with file name mapped
    for record in records_annotated:
        if record.id in dataset_to_real_name:
            record.id = dataset_to_real_name[record.id]
 #   # try:
 #       # if not records_annotated:
 #           # print('records_annotated is empty')
 #       # else:
 #          # print(f'records_annotated is:\n{records_annotated}')
 #    except Exception as e:
 #        print(f'An error occurred: {e}')

    # generate annotated genbank files
    output_dir = outdir_gb
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for record in records_annotated:
        file_path = os.path.join(output_dir, f"{record.id}.gb")
        with open(file_path, "w") as output_handle:
            SeqIO.write(record, output_handle, "genbank")

    # generate PDF report
    cr.breaches_records_to_pdf(records_annotated, output_pdf)

    return output_tsv, output_pdf, output_dir


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Evaluate manufacturability of DNA sequences.")

    parser.add_argument("--files_to_evaluate", required=True,
                        help="List of GenBank files (Comma-separated)")
    parser.add_argument('--file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original filenames')
    parser.add_argument("--output_tsv", required=True, help="Excel file name")
    parser.add_argument("--output_pdf", required=True, help="PDF file name")
    parser.add_argument("--outdir_gb", required=True, help="DIR for annotated GenBank files")
    parser.add_argument("--use_file_names_as_id", type=lambda x: x.lower() == 'true', default=True,
                        help="Use file names as IDs (True/False)")
    parser.add_argument("--avoid_patterns", required=True,
                        help="List of patterns to avoid (comma-separated, e.g., 'BsaI_site,BsmBI_site')")
    parser.add_argument("--hairpin_constraints", required=True,
                        help="Hairpin constraints as 'stem_size;window_size' (space-separated, e.g., '20;200 30;250')")
    parser.add_argument("--gc_constraints", required=True,
                        help="GC content constraints as 'min;max;window' (space-separated, e.g., '0.3;0.7;100 0.1;0.3;100')")
    parser.add_argument("--kmer_size", required=True,
                        help="K-mer uniqueness size (e.g., '15')")

    return parser.parse_args()


def extract_constraints_from_args(args):
    """Extract constraints directly from the command-line arguments."""

    avoid_patterns = args.avoid_patterns.split(',')

    hairpin_constraints = [hp.strip() for hp in args.hairpin_constraints.split(' ')]

    gc_constraints = [gc.strip() for gc in args.gc_constraints.split(' ')]

    kmer_size = [int(k) for k in args.kmer_size.split(',')]

    return avoid_patterns, hairpin_constraints, gc_constraints, kmer_size


if __name__ == "__main__":
    args = parse_command_line_args()

    avoid_patterns, hairpin_constraints, gc_constraints, kmer_size = extract_constraints_from_args(args)

    evaluate_manufacturability(
        args.files_to_evaluate, args.file_name_mapping, args.output_tsv, args.output_pdf,
        args. outdir_gb, args.use_file_names_as_id, avoid_patterns,
        hairpin_constraints, gc_constraints, kmer_size
    )
