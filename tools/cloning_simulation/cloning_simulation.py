import argparse
import os
import json
import zipfile
import pandas
import dnacauldron


def cloning_simulation(files_to_assembly, domesticated_list,
                       csv_file, assembly_type, topology,
                       file_name_mapping, file_name_mapping_dom,
                       use_file_names_as_id,
                       outdir_simulation, output_simulation, enzyme, outdir_gb):

    files_to_assembly = files_to_assembly.split(',')

    repository = dnacauldron.SequenceRepository()
    repository.import_records(files=files_to_assembly,
                              use_file_names_as_ids=use_file_names_as_id,
                              topology=topology)
    if domesticated_list:
       domesticated_files = domesticated_list.split(',')
       repository.import_records(files=domesticated_files,
                                 use_file_names_as_ids=use_file_names_as_id,
                                 topology=topology)

    # refine the real record name dict
    if isinstance(file_name_mapping, str):
        file_name_mapping = dict(
            item.split(":") for item in file_name_mapping.split(",")
        )
    real_names = {
        os.path.splitext(os.path.basename(k))[0]: v.replace(".gb", "")
        for k, v in file_name_mapping.items()
    }

    # refine the real record name dict_dom
    if file_name_mapping_dom == "":
        file_name_mapping_dom = {}
    else:
        if isinstance(file_name_mapping_dom, str):
            file_name_mapping_dom = dict(
                item.split(":") for item in file_name_mapping_dom.split(",")
            )
        dom_real_names = {
            os.path.splitext(os.path.basename(k))[0]: v.replace(".gb", "")
            for k, v in file_name_mapping_dom.items()
        }
        real_names.update(dom_real_names)

    # update the records

    for key, record in list(repository.collections["parts"].items()):
        current_id = record.id
        if current_id in real_names:
            new_id = real_names[current_id]
            record.id = new_id
            record.name = new_id
            record.description = new_id
            repository.collections["parts"][new_id] = repository.collections["parts"].pop(key)
    ########################################################
    # print (f"repo: {vars(repository)}")
    # any(pandas.read_csv(csv_file, index_col=0, header=None).duplicated())
    df = pandas.read_csv(csv_file, index_col=0, header=None)
    if df.duplicated().any():
        raise ValueError("Duplicate rows found in the data!")

    if assembly_type == "Type2sRestrictionAssembly":
        assembly_class = dnacauldron.Type2sRestrictionAssembly
    elif assembly_type == "GibsonAssembly":
        assembly_class = dnacauldron.GibsonAssembly
    elif assembly_type == "BASICAssembly":
        assembly_class = dnacauldron.BASICAssembly
    elif assembly_type == "BioBrickStandardAssembly":
        assembly_class = dnacauldron.BioBrickStandardAssembly
    elif assembly_type == "OligoPairAnnealin":
        assembly_class = dnacauldron.OligoPairAnnealin
    elif assembly_type == "LigaseCyclingReactionAssembly":
        assembly_class = dnacauldron.LigaseCyclingReactionAssembly
    else:
        raise ValueError(f"Unsupported assembly type: {assembly_type}")

    new_csvname = "assambly.csv"
    os.rename(csv_file, new_csvname)

    assembly_plan = dnacauldron.AssemblyPlan.from_spreadsheet(
        name="auto_from_filename",
        path=new_csvname,
        dataframe=None,
        header=None,
        assembly_class=assembly_class
    )
    if enzyme != 'auto':
        for assembly in assembly_plan.assemblies:
            assembly.enzyme = enzyme

    simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = simulation.compute_stats()
    print(stats)

    report_writer = dnacauldron.AssemblyReportWriter(
        include_mix_graphs=True,
        include_assembly_plots=True,
        show_overhangs_in_graph=True,
        annotate_parts_homologies=True,
        include_pdf_report=True,
    )
    simulation.write_report(outdir_simulation, assembly_report_writer=report_writer)

    # Append report files to .dat (ZIP)
    with zipfile.ZipFile(output_simulation, mode='a', compression=zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(outdir_simulation):
            for file in files:
                full_path = os.path.join(root, file)
                arcname = os.path.relpath(full_path, outdir_simulation)
                zipf.write(full_path, arcname)
 #       print("Files in the zip archive:")
 #       for info in zipf.infolist():
 #           print(info.filename)
        for member in zipf.namelist():
            # Only extract actual files inside 'all_construct_records/' (not subfolders)
            if member.startswith("assambly_simulation/all_construct_records/") and not member.endswith("/"):
                # Get the file name only (strip folder path)
                filename = os.path.basename(member)
                if not filename:
                    continue  # skip any edge cases

                # Destination path directly in outdir_dir
                target_path = os.path.join(outdir_gb, filename)

                # Write the file content
                with zipf.open(member) as source, open(target_path, "wb") as target:
                    target.write(source.read())

    return output_simulation, outdir_gb


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Domestication")

    parser.add_argument("--parts_files", required=True,
                        help="List of GenBank files (Comma-separated)")
    parser.add_argument("--domesticated_seq", required=True,
                        help="output of domestication (ganbank list)")
    parser.add_argument("--assembly_csv", required=True,
                        help="csv assembly")
    parser.add_argument('--assembly_plan_name', type=str, required=False,
                        help='type of assembly')
    parser.add_argument('--topology', type=str, required=False,
                        help='"circular" or "linear"')
    parser.add_argument('--file_name_mapping', type=str,
                        help='Mapping of Galaxy filenames to original filenames')
    parser.add_argument('--file_name_mapping_dom', type=str,
                        help='Mapping of Galaxy filenames to original domestication filenames')
    parser.add_argument("--use_file_names_as_id", type=lambda x: x.lower() == 'true', default=True,
                        help="Use file names as IDs (True/False)")
    parser.add_argument("--outdir_simulation", required=True,
                        help="dir output for cloning simulation results")
    parser.add_argument("--output_simulation", required=True,
                        help="zip output for cloning simulation results")
    parser.add_argument('--enzyme', type=str,required=False,
                        help='enzyme to use')
    parser.add_argument("--outdir_gb", required=True,
                        help="dir output constructs gb files")
    parser.add_argument("--use_json_paramers", required=True,
                         help="Use parameters from JSON: true/false")
    parser.add_argument("--json_conf", required=False,
                         help="JSON config file with DB parameters")
 
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line_args()

    #json param checking 
    config_params = {}
    use_json = args.use_json_paramers == 'true'
    if use_json:
        if not args.json_conf:
            raise ValueError("You must provide --json_conf when --use_json_paramers is 'true'")
        with open(args.json_conf, "r") as f:
            config_params = json.load(f)
    else:
        config_params = {
            "assembly_plan_name": args.assembly_plan_name,
            "topology": args.topology,
            "enzyme": args.enzyme
        }
    assembly_plan_name = config_params["assembly_plan_name"]
    topology = config_params["topology"]
    enzyme = config_params["enzyme"]

    cloning_simulation(
        args.parts_files, args.domesticated_seq,
        args.assembly_csv, assembly_plan_name, topology,
        args.file_name_mapping, args.file_name_mapping_dom,
        args.use_file_names_as_id, args.outdir_simulation,
        args.output_simulation, enzyme, args.outdir_gb
    )
