import argparse
import os
import zipfile
import pandas
import dnacauldron
import genedom
import proglog
import shutil
# proglog.notebook()
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def domestication(files_to_domestication, csv_file,file_name_mapping, use_file_names_as_id,
                  allow_edits, output_dom, gb_dom):

    file_to_domestication = files_to_domestication.split(',')

    records_to_domesticate = dnacauldron.biotools.load_records_from_files(
        files=file_to_domestication,
        folder=None,
        use_file_names_as_ids=use_file_names_as_id
    )

    # refine the real record name dict
    if isinstance(file_name_mapping, str):
        file_name_mapping = dict(
            item.split(":") for item in file_name_mapping.split(",")
        )
    real_names = {
        os.path.splitext(os.path.basename(k))[0]: v.replace(".gb", "")
        for k, v in file_name_mapping.items()
    }
    updated_records = []
    for record in records_to_domesticate:
        original_id = record.id
        if original_id in real_names:
            new_id = real_names[original_id]
            record.id = new_id
        updated_records.append(record)

    df = pandas.read_csv(csv_file)
    EMMA_PLUS = genedom.GoldenGateDomesticator.standard_from_spreadsheet(dataframe=df)
    genedom.batch_domestication(
        records=updated_records,
        standard=EMMA_PLUS,
        allow_edits=allow_edits,
        target=output_dom)

    # Check if any names were truncated:
    if isinstance(output_dom, str) and output_dom.endswith(".zip"):
        with zipfile.ZipFile(output_dom, "r") as zipf:
            with zipf.open("order_ids.csv") as f:
                order_ids = pandas.read_csv(f)
    else:
        order_ids = pandas.read_csv(os.path.join(output_dom, "order_ids.csv"))

    any_truncated = False
    for index, row in order_ids.iterrows():
        if row["sequence"] != row["order_id"]:
            any_truncated = True
            print("Changed names:", end=" ")
            print(" --> ".join(row))
    if not any_truncated:
        print("Part names were not truncated")

    # zip compressing
    if os.path.isdir(output_dom):
        zip_path = output_dom.rstrip("/\\") + ".zip"
        shutil.make_archive(output_dom, 'zip', output_dom)
        shutil.move(zip_path, output_dom)

    temp_extract_dir = 'temp_extracted'
    os.makedirs(temp_extract_dir, exist_ok=True)

    with zipfile.ZipFile(output_dom, 'r') as zip_ref:
        zip_ref.extractall(temp_extract_dir)
    # Navigate into domesticated_genbanks and copy .gb files
    dom_gb_dir = os.path.join(temp_extract_dir, 'domesticated_genbanks')
    if os.path.isdir(dom_gb_dir):
        for file in os.listdir(dom_gb_dir):
            if file.endswith('.gb'):
                shutil.copy(os.path.join(dom_gb_dir, file), gb_dom)
                # print(f"Copied: {file}")
    else:
        print("Folder 'domesticated_genbanks' not found inside zip.")

    # Optional: clean up
    shutil.rmtree(temp_extract_dir)

    return output_dom


def methylation_protection(domestication_target, output_methprot):

    if domestication_target.endswith(".zip"):
        extracted_dir = "extracted_genbanks"
        with zipfile.ZipFile(domestication_target, "r") as zipf:
            os.makedirs(extracted_dir, exist_ok=True)
            for member in zipf.namelist():
                if member.startswith("domesticated_genbanks/") and member.endswith(".gb"):
                    zipf.extract(member, path=extracted_dir)
        gb_folder = os.path.join(extracted_dir, "domesticated_genbanks")
    else:
        gb_folder = os.path.join(domestication_target, "domesticated_genbanks")

    records_to_protect = dnacauldron.biotools.load_records_from_files(folder=gb_folder, use_file_names_as_ids=True)

    for record in records_to_protect:
        new_seqrecord = SeqRecord("TTC") + record + SeqRecord("GAA")  # these sequences are designed to prevent Dcm recognition site formation
        new_seqrecord.id = record.id
        new_seqrecord.name = record.name
        new_seqrecord.annotations = {"molecule_type": "DNA"}

        target_file = os.path.join(output_methprot, new_seqrecord.id + ".gb")
        with open(target_file, "w") as output_handle:
            SeqIO.write(new_seqrecord, output_handle, "genbank")


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Domestication")

    parser.add_argument("--files_to_domestication", required=True,
                        help="List of GenBank files (Comma-separated)")
    parser.add_argument("--csv_file", required=True,
                        help="csv file")
    parser.add_argument('--file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original filenames')
    parser.add_argument("--use_file_names_as_id", type=lambda x: x.lower() == 'true', default=True,
                        help="Use file names as IDs (True/False)")
    parser.add_argument("--allow_edits", type=lambda x: x.lower() == 'true', default=True,
                        help="Allow sequence edits")
    parser.add_argument("--output_dom", required=True,
                        help="zip output for domestication results")
    parser.add_argument("--output_methprot", required=True,
                        help="gb output for methylation protection")
    parser.add_argument("--output_gb_dom", required=True,
                        help="gb output for mdomesticated gb files")
    parser.add_argument("--methylation_protection", type=lambda x: x.lower() == 'true', default=False,
                        help="Enable methyl protection (true/false)")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line_args()

    domestication(
        args.files_to_domestication, args.csv_file,
        args.file_name_mapping, args.use_file_names_as_id,
        args.allow_edits, args.output_dom, args.output_gb_dom
    )

    if args.methylation_protection:
        methylation_protection(
            args.output_dom, args.output_methprot
        )
