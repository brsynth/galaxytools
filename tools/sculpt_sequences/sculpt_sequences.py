import argparse
import os
import re
import dnacauldron
import dnachisel
from Bio import SeqIO
import proglog


def smart_number(val):
        try:
            float_val = float(val)
            if float_val.is_integer():
                return int(float_val)
            else:
                return float_val
        except ValueError:
            raise ValueError(f"Invalid number: {val}")


def sculpt_sequances(files_to_sculpt, file_name_mapping, outdir_scul, outdir_unscul, use_file_names_as_id,
                     avoid_patterns, enforce_gc_content, DnaOptimizationProblemClass,
                     kmer_size,hairpin_constraints):

    files_to_sculpt = files_to_sculpt.split(',')

    records_to_sculpt = dnacauldron.biotools.load_records_from_files(
        files=files_to_sculpt,
        folder=None,
        use_file_names_as_ids=use_file_names_as_id
    )

    constraint_list = []

 # avoid_patterns pearsing
    for pattern in avoid_patterns:
        constraint_list.append(dnachisel.AvoidPattern(pattern))
        print(f"AvoidPattern constraint: {pattern}")

 # gc_constraints pearsing
    for constraint in enforce_gc_content:
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

 # hairpin_constraints pearsing
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

 # k_size pearsing
    for k_size in kmer_size:
        try:
            constraint_list.append(dnachisel.UniquifyAllKmers(k=int(k_size)))
            print(f"k-mer size is: {k_size}")
        except ValueError:
            print(f"Skipping invalid k-mer size: {k_size}")
    
 #   refine the real record name dict
    if isinstance(file_name_mapping, str):
        file_name_mapping = dict(
            item.split(":") for item in file_name_mapping.split(",")
        )
    real_names = {
        os.path.splitext(os.path.basename(k))[0]: v.replace(".gb", "")
        for k, v in file_name_mapping.items()
    } 

    counter = 0
 #   not_divisible_seqs = []
    for record in records_to_sculpt:

 #       Look up the "real name" from the simplified dict (file_name_mapping)
        real_name = real_names.get(record.id, record.id)
        if real_name is None:
            raise ValueError(f"Error: {record.id} not found in the fragment names dictionary!")

        if DnaOptimizationProblemClass == 'DnaOptimizationProblem':
            problem = dnachisel.DnaOptimizationProblem.from_record(record, extra_constraints=constraint_list)
        elif DnaOptimizationProblemClass == 'CircularDnaOptimizationProblem':
            problem = dnachisel.CircularDnaOptimizationProblem.from_record(record, extra_constraints=constraint_list)

        if not problem.all_constraints_pass():
            counter += 1
            print(real_name)
            len_seq = len(problem.sequence)
            print("Length:", len_seq)
            target_path = os.path.join(outdir_scul, real_name + ".zip")
            problem.optimize_with_report(target=target_path)
            print()
            print(problem.constraints_text_summary())
            print()
        else:  # just save the genbank
            genbank_filename = real_name + ".gb"
            with open(os.path.join(outdir_unscul, genbank_filename), "w") as output_handle:
                SeqIO.write(record, output_handle, "genbank")
    print()
    print("Counter:", counter)
 #   print("Contents of outdir_scul:")
 #   for f in os.listdir(outdir_scul):
 #       print(" -", f)

 #   Unzip all files in scul
    import zipfile
    for file in os.listdir(outdir_scul):
        if file.endswith(".zip"):
            zip_path = os.path.join(outdir_scul, file)
            extract_dir = os.path.join(outdir_scul, os.path.splitext(file)[0])
            os.makedirs(extract_dir, exist_ok=True)
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
 #   print("Contents of outdir_scul:")
 #   for f in os.listdir(outdir_scul):
 #       print(" -", f)
    reports = os.listdir(outdir_scul)
 #   print (f'report is {reports}')
    from shutil import copyfile
    for report in reports:
        if report.endswith(".zip"):
            continue
        gb_file = os.path.join(outdir_scul, report, "final_sequence.gb")

        target_filename = report + ".gb"
        target_file = os.path.join(outdir_unscul, target_filename)
        print(gb_file, target_file)
        copyfile(gb_file, target_file)
 #   print("Contents of outdir_unscul:")
 #   for f in os.listdir(outdir_unscul):
 #         print(" -", f)

    return outdir_scul, outdir_unscul


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Evaluate manufacturability of DNA sequences.")

    parser.add_argument("--files_to_sculpt", required=True,
                        help="List of GenBank files (Comma-separated)")
    parser.add_argument('--file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original filenames')
    parser.add_argument("--outdir_scul", required=True, help="scul file dir")
    parser.add_argument("--outdir_unscul", required=True, help="unscul file dir")
    parser.add_argument("--use_file_names_as_id", type=lambda x: x.lower() == 'true', default=True,
                        help="Use file names as IDs (True/False)")
    parser.add_argument("--avoid_patterns", required=True,
                        help="List of patterns to avoid (comma-separated, e.g., 'BsaI_site,BsmBI_site')")
    parser.add_argument("--gc_constraints", required=True,
                        help="GC content constraints as 'min;max;window' (space-separated, e.g., '0.3;0.7;100 0.1;0.3;100')")
    parser.add_argument("--DnaOptimizationProblemClass", required=True,
                        help="the class to use for DnaOptimizationProblem")
    parser.add_argument("--hairpin_constraints", required=True,
                        help="Hairpin constraints as 'stem_size;window_size' (space-separated, e.g., '20;200 30;250')")
    parser.add_argument("--kmer_size", required = True,
                        help="K-mer uniqueness size (e.g., '15')")

    return parser.parse_args()


def extract_constraints_from_args(args):
    """Extract constraints directly from the command-line arguments."""

    split_pattern = r'(?:__cn__|\s{2,})'

    # 1. Avoid patterns (split by any whitespace)
    avoid_patterns = re.split(split_pattern, args.avoid_patterns.strip()) if args.avoid_patterns.strip() else []

    # 2. Hairpin constraint: one dictionary (print as string later)
    hairpin_constraints = re.split(split_pattern,args.hairpin_constraints.strip()) if args.hairpin_constraints.strip() else []

    # 3. GC constraints: split by 2+ spaces or newlines
    gc_constraints = re.split(split_pattern, args.gc_constraints.strip()) if args.gc_constraints.strip() else []

    # 4. k-mer size: single value or list
    kmer_size = [int(k.strip()) for k in args.kmer_size.strip().split(',') if k.strip()] if args.kmer_size.strip() else []
    
    return avoid_patterns, hairpin_constraints, gc_constraints, kmer_size



if __name__ == "__main__":
    args = parse_command_line_args()

    avoid_patterns, hairpin_constraints, gc_constraints, kmer_size, = extract_constraints_from_args(args)

    sculpt_sequances(
        args.files_to_sculpt, args.file_name_mapping,
        args.outdir_scul, args.outdir_unscul,
        args.use_file_names_as_id, avoid_patterns,
        gc_constraints, args.DnaOptimizationProblemClass,
        kmer_size, hairpin_constraints
    )
