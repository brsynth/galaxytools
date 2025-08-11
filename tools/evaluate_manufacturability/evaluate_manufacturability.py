import argparse
import sys
import json
import os
import re
from Bio import SeqIO
import dnacauldron
import dnachisel
import dnachisel.reports.constraints_reports as cr


def smart_number(val):
        try:
            float_val = float(val)
            if float_val.is_integer():
                return int(float_val)
            else:
                return float_val
        except ValueError:
            raise ValueError(f"Invalid number: {val}")
        

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
 #   try:
 #       if not records_to_evaluate:
 #           print('records to evaluate: is empty')
 #       else:
 #           for record in records_to_evaluate:
 #               print(f'records to evaluate: {record}')
 #   except Exception as e:
 #       print(f'An error occurred: {e}')

    length_cutoff = 100
    for part in records_to_evaluate:
        if len(part) < length_cutoff:
            raise ValueError(f"Error: Sequence '{part.id}' is too short (length {len(part)}). Minimum required: {length_cutoff}")

    constraint_list = []

 # avoid_patterns pearsing
    for pattern in avoid_patterns:
        constraint_list.append(dnachisel.AvoidPattern(pattern))
        print(f"AvoidPattern constraint: {pattern}")

 # hairpin_constraints pearsing
    print (f'hairpin_constraints is {hairpin_constraints}')
    for constraint in hairpin_constraints:
        constraints = [c.strip() for c in constraint.split('  ') if c.strip()] 
        for line in constraints:
            if not line:
                continue
            print(f"Hairpin constraint: {line}")

            try:
                pairs = [kv.strip() for kv in line.split(',')]
                params = {}
                for pair in pairs:
                    key, val = pair.split('=')
                    key = key.strip()
                    val = val.strip()
                    params[key] = smart_number(val)
                    
                constraint_list.append(dnachisel.AvoidHairpins(**params))

            except Exception as e:
                print(f"Skipping invalid hairpin_constraints: {line} ({e})")
 # gc_constraints pearsing
    for constraint in gc_constraints:
        constraints = [c.strip() for c in constraint.split('  ') if c.strip()] 
        for line in constraints:
            if not line:
                continue
            print(f"GC constraint: {line}")

            try:
                pairs = [kv.strip() for kv in line.split(',')]
                params = {}
                for pair in pairs:
                    key, val = pair.split('=')
                    key = key.strip()
                    val = val.strip()
                    params[key] = smart_number(val)

                constraint_list.append(dnachisel.EnforceGCContent(**params))

            except Exception as e:
                print(f"Skipping invalid gc_constraints: {line} ({e})")

 # k_size pearsing
    for k_size in kmer_size:
        try:
            constraint_list.append(dnachisel.UniquifyAllKmers(k=int(k_size)))
            print(f"k-mer size is: {k_size}")
        except ValueError:
            print(f"Skipping invalid k-mer size: {k_size}")

 #   print(f'constraint_list is{constraint_list}')

 # constraint_list apply
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
    try:
        if dataframe.empty:
            print('dataframe: is empty')
        else:
            print(f'dataframe is:\n{dataframe}')
            print(dataframe.columns)
    except Exception as e:
        print(f'An error occurred: {e}')
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

    parser.add_argument("--files_to_evaluate", required=False,
                        help="List of GenBank files (Comma-separated)")
    parser.add_argument('--file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original filenames')
    parser.add_argument('--DB_file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original DB filenames')
    parser.add_argument("--output_tsv", required=True, help="Excel file name")
    parser.add_argument("--output_pdf", required=True, help="PDF file name")
    parser.add_argument("--outdir_gb", required=True, help="DIR for annotated GenBank files")
    parser.add_argument("--use_file_names_as_id", type=lambda x: x.lower() == 'true', default=True,
                        help="Use file names as IDs (True/False)")
    parser.add_argument("--avoid_patterns", required=False,
                        help="List of patterns to avoid (comma-separated, e.g., 'BsaI_site,BsmBI_site')")
    parser.add_argument("--hairpin_constraints", required=False,
                        help="Hairpin constraints as 'stem_size;window_size' (space-separated, e.g., '20;200 30;250')")
    parser.add_argument("--gc_constraints", required=False,
                        help="GC content constraints as 'min;max;window' (space-separated, e.g., '0.3;0.7;100 0.1;0.3;100')")
    parser.add_argument("--kmer_size", required=False,
                        help="K-mer uniqueness size (e.g., '15')")
    parser.add_argument("--json_params", required=False,
                        help="JSON params for the tool")
    parser.add_argument("--use_json_param", required=True,
                            help="If use JSON as param source")
    parser.add_argument("--mode", required=True,
                            help="mode d'utilisation: standard ou workflow")
    parser.add_argument("--DB_report", required=False,
                            help="In wkf mode")
    parser.add_argument("--DB_genbank_files", required=False,
                            help="IN wkf mode")


    return parser.parse_args()



def extract_constraints_from_args(args):
    """Extract constraints directly from the command-line arguments."""

    split_pattern = r'(?:__cn__|\s{2,})'

    # 1. Avoid patterns (split by any whitespace)
    avoid_patterns = re.split(split_pattern, args.avoid_patterns.strip())

    # 2. Hairpin constraint: one dictionary (print as string later)
    hairpin_constraints = re.split(split_pattern,args.hairpin_constraints.strip())

    # 3. GC constraints: split by 2+ spaces or newlines
    gc_constraints = re.split(split_pattern, args.gc_constraints.strip())

    # 4. k-mer size: single value or list
    kmer_size = [int(k.strip()) for k in args.kmer_size.strip().split(',') if k.strip()]
    
    return avoid_patterns, hairpin_constraints, gc_constraints, kmer_size


def load_params_from_json(json_path):
    with open(json_path, 'r') as f:
        params = json.load(f)
    return params


def load_constraints_from_json(json_path):
    with open(json_path, 'r') as f:
        params = json.load(f)

    def split_lines(val):
        if isinstance(val, str):
            return [line.strip() for line in val.strip().split('\n') if line.strip()]
        return val

    avoid_patterns = split_lines(params.get("avoid_patterns", ""))
    hairpin_constraints = split_lines(params.get("hairpin_constraints", ""))
    gc_constraints = split_lines(params.get("gc_constraints", ""))
    kmer_size = [int(k.strip()) for k in str(params.get("kmer_size", "")).split(',') if k.strip()]

    return {
        "avoid_patterns": avoid_patterns,
        "hairpin_constraints": hairpin_constraints,
        "gc_constraints": gc_constraints,
        "kmer_size": kmer_size
    }


if __name__ == "__main__":

    args = parse_command_line_args()

    ###
    if "--mode" in sys.argv:
        mode_index = sys.argv.index("--mode") + 1
        mode = sys.argv[mode_index].strip()

    skip_evaluation = False
    use_DB_files = False
    DB_genbank_files = []

    if mode == "wkf":

        if "--DB_report" not in sys.argv:
            print("ERROR: --DB_report is required in wkf mode.")
            sys.exit(1)
        db_index = sys.argv.index("--DB_report") + 1
        db_report_path = sys.argv[db_index]

        if "--DB_genbank_files" in sys.argv:
            db_gb_index = sys.argv.index("--DB_genbank_files") + 1
            DB_genbank_files = sys.argv[db_gb_index].split(",")
        else:
            DB_genbank_files = []

        if not os.path.isfile(db_report_path):
            print(f"ERROR: DB report file not found at {db_report_path}")
            sys.exit(1)

        with open(db_report_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        if not lines:
            print("ERROR: DB_report file is empty.")
            sys.exit(1)

        header = lines[0]
        if header.startswith("Missing fragment in DB:"):

            missing_fragments = lines[1:]

           # Parse file_name_mapping
            if isinstance(args.file_name_mapping, str):
                mapping_dict = dict(item.split(":") for item in args.file_name_mapping.split(","))
            else:
                mapping_dict = {}

            # Logical names
            provided_filenames = [os.path.splitext(v)[0] for v in mapping_dict.values()]

            # print(f'provided_filenames is : {provided_filenames}')

            unmatched = [
                frag for frag in missing_fragments
                if os.path.splitext(frag)[0] not in provided_filenames
            ]

            if unmatched:
                print(f"ERROR: The following missing fragment(s) must be provided as .gb files: {', '.join(unmatched)}")
                sys.exit(1)
            else:
                use_DB_files = True  # Append after evaluation

        elif header.startswith("NO missing fragments in DB"):
            skip_evaluation = True
        else:
            print(f"ERROR: Invalid header in DB_report: '{header}'")
            sys.exit(1)
    ###

    # Default values from command-line
    avoid_patterns, hairpin_constraints, gc_constraints, kmer_size = extract_constraints_from_args(args)

    
    # Check if the flag --use_json_param is present and set to true
    if "--use_json_param" in sys.argv:
        use_json_index = sys.argv.index("--use_json_param") + 1
        use_json = sys.argv[use_json_index].lower() == "true"
    else:
        use_json = False

    # Now only check --json_params if use_json is True
    if use_json:
        if "--json_params" in sys.argv:
            json_index = sys.argv.index("--json_params") + 1
            json_file = sys.argv[json_index]
            if json_file.lower() != "none":
                json_constraints = load_constraints_from_json(json_file)
                avoid_patterns = json_constraints["avoid_patterns"]
                hairpin_constraints = json_constraints["hairpin_constraints"]
                gc_constraints = json_constraints["gc_constraints"]
                kmer_size = json_constraints["kmer_size"]

    params = {
        "files_to_evaluate": args.files_to_evaluate,
        "file_name_mapping": args.file_name_mapping,
        "output_tsv": args.output_tsv,
        "output_pdf": args.output_pdf,
        "outdir_gb": args.outdir_gb,
        "use_file_names_as_id": args.use_file_names_as_id,
        "avoid_patterns": avoid_patterns,
        "hairpin_constraints": hairpin_constraints,
        "gc_constraints": gc_constraints,
        "kmer_size": kmer_size
    }

    if not skip_evaluation:
        evaluate_manufacturability(
            params["files_to_evaluate"], params["file_name_mapping"],
            params["output_tsv"], params["output_pdf"], params["outdir_gb"],
            params["use_file_names_as_id"], params["avoid_patterns"],
            params["hairpin_constraints"], params["gc_constraints"],
            params["kmer_size"]
        )

    if mode == "wkf" and (skip_evaluation or use_DB_files):
        if DB_genbank_files:
            print(f"DB_genbank_files is: {DB_genbank_files}")
            print("Adding DB GenBank files to output collection using DB_file_name_mapping...")

            os.makedirs(params["outdir_gb"], exist_ok=True)

            # mapping real DB gb file name
            if isinstance(args.DB_file_name_mapping, str):
                print (f'DB_file_name_mapping is: {args.DB_file_name_mapping}')
                DB_mapping_dict = dict(item.split(":") for item in args.DB_file_name_mapping.split(","))
            else:
                DB_mapping_dict={}

            for path in DB_genbank_files:
                basename = os.path.basename(path)
                logical_name = DB_mapping_dict.get(path) or DB_mapping_dict.get(basename)

                if not logical_name:
                    print(f"WARNING: No mapping found for DB GenBank file: {path}. Skipping.")
                    continue

                output_filename = os.path.splitext(logical_name)[0] + ".gb"
                dest_path = os.path.join(params["outdir_gb"], output_filename)

                try:
                    with open(path, 'r') as src, open(dest_path, 'w') as dst:
                        dst.write(src.read())
                    print(f"Copied and renamed: {path} → {dest_path}")
                except Exception as e:
                    print(f"ERROR: Failed to copy {path} → {dest_path}: {e}")
        else:
            print("No DB GenBank files to append, continuing without error.")


