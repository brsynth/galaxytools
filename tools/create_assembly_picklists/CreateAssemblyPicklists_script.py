#!/usr/bin/env python
# coding: utf-8
# Code copied from CUBA backend tools.py and create_assembly_picklists/CreateAssemblyPicklistsView.py
# Code modified for running in a script in Galaxy.
##############################################################################
##############################################################################
# App code
## EGF Galaxy Create assembly picklists -- script

##############################################################################
# IMPORTS
import argparse
import os
from io import StringIO, BytesIO
import re
from base64 import b64encode, b64decode
from copy import deepcopy
import sys

from collections import OrderedDict
from fuzzywuzzy import process
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import bandwagon as bw
import crazydoc
from dnachisel.biotools import sequence_to_biopython_record
import dnacauldron
import flametree
from plateo import AssemblyPlan
from plateo.parsers import plate_from_content_spreadsheet
from plateo.containers import Plate4ti0960
from plateo.exporters import AssemblyPicklistGenerator, picklist_to_assembly_mix_report
from plateo.exporters import (
    picklist_to_labcyte_echo_picklist_file,
    picklist_to_tecan_evo_picklist_file,
    plate_to_platemap_spreadsheet,
    PlateTextPlotter,
)
from plateo.tools import human_volume
from snapgene_reader import snapgene_file_to_seqrecord


##############################################################################
# FUNCTIONS

def fix_and_rename_paths(paths):
    fixed_paths = []
    for path in paths:
        new_path = path.replace("__sq__", "'")
        if new_path != path:
            os.rename(path, new_path)
        fixed_paths.append(new_path)
    return fixed_paths


def parse_optional_float(x):
    if x == '':
        return None
    return float(x)


def did_you_mean(name, other_names, limit=5, min_score=50):  # test
    results = process.extract(name, list(other_names), limit=limit)
    return [e for (e, score) in results if score >= min_score]


def fix_ice_genbank(genbank_txt):
    lines = genbank_txt.splitlines()
    lines[0] += max(0, 80 - len(lines[0])) * " "
    return "\n".join(lines)


def write_record(record, target, fmt="genbank"):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    if fmt == "genbank":
        if isinstance(record, (list, tuple)):
            for r in record:
                r.name = r.name[:20]
        else:
            record.name = record.name[:20]
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, fmt)


def autoname_genbank_file(record):
    return record.id.replace(".", "_") + ".gb"


def string_to_records(string):
    """Convert a string of a fasta, genbank... into a simple ATGC string.

    Can also be used to detect a format.
    """
    matches = re.match("([ATGC][ATGC]*)", string)
    # print("============", len(matches.groups()[0]), len(string))
    # print (matches.groups()[0] == string)
    if (matches is not None) and (matches.groups()[0] == string):
        return [SeqRecord(Seq(string))], "ATGC"

    for fmt in ("fasta", "genbank"):
        if fmt == "genbank":
            string = fix_ice_genbank(string)
        try:
            stringio = StringIO(string)
            records = list(SeqIO.parse(stringio, fmt))
            if len(records) > 0:
                return (records, fmt)
        except:
            pass
    try:
        record = snapgene_file_to_seqrecord(filecontent=StringIO(string))
        return [record]
    except:
        pass
    raise ValueError("Invalid sequence format")


def file_to_filelike_object(file_, type="byte"):
    content = file_.content.split("base64,")[1]
    filelike = BytesIO if (type == "byte") else StringIO
    return filelike(b64decode(content))


def spreadsheet_file_to_dataframe(filedict, header="infer"):
    filelike = file_to_filelike_object(filedict)
    if filedict.name.endswith(".csv"):
        return pandas.read_csv(filelike, header=header)
    else:
        return pandas.read_excel(filelike, header=header)


def records_from_zip_file(zip_file, use_file_names_as_ids=False):
    zip_name = zip_file.name
    zip_file = flametree.file_tree(file_to_filelike_object(zip_file))
    records = []
    for f in zip_file._all_files:
        ext = f._extension.lower()
        if ext in ["gb", "gbk", "fa", "dna"]:
            try:
                new_records, fmt = string_to_records(f.read())
                if not isinstance(new_records, list):
                    new_records = [new_records]
            except:
                content_stream = BytesIO(f.read("rb"))
                try:
                    record = snapgene_file_to_seqrecord(fileobject=content_stream)
                    new_records, fmt = [record], "snapgene"
                except:
                    try:
                        parser = crazydoc.CrazydocParser(
                            ["highlight_color", "bold", "underline"]
                        )
                        new_records = parser.parse_doc_file(content_stream)
                        fmt = "doc"
                    except:
                        raise ValueError("Format not recognized for file " + f._path)

            single_record = len(new_records) == 1
            for i, record in enumerate(new_records):
                name = record.id
                if name in [
                    None,
                    "",
                    "<unknown id>",
                    ".",
                    " ",
                    "<unknown name>",
                ]:
                    number = "" if single_record else ("%04d" % i)
                    name = f._name_no_extension.replace(" ", "_") + number
                record.id = name
                record.name = name
                record.file_name = f._name_no_extension
                record.zip_file_name = zip_name
                if use_file_names_as_ids and single_record:
                    basename = os.path.basename(record.file_name)
                    basename_no_extension = os.path.splitext(basename)[0]
                    record.id = basename_no_extension
            records += new_records
    return records


def records_from_data_file(data_file):
    content = b64decode(data_file.content.split("base64,")[1])
    try:
        records, fmt = string_to_records(content.decode("utf-8"))
    except:
        try:
            record = snapgene_file_to_seqrecord(fileobject=BytesIO(content))
            records, fmt = [record], "snapgene"
        except:
            try:
                parser = crazydoc.CrazydocParser(
                    ["highlight_color", "bold", "underline"]
                )
                records = parser.parse_doc_file(BytesIO(content))
                fmt = "doc"
            except:
                try:
                    df = spreadsheet_file_to_dataframe(data_file, header=None)
                    records = [
                        sequence_to_biopython_record(sequence=seq, id=name, name=name)
                        for name, seq in df.values
                    ]
                    fmt = "spreadsheet"
                except:
                    raise ValueError("Format not recognized for file " + data_file.name)
    if not isinstance(records, list):
        records = [records]
    return records, fmt


def record_to_formated_string(record, fmt="genbank", remove_descr=False):
    if remove_descr:
        record = deepcopy(record)
        if isinstance(record, (list, tuple)):
            for r in record:
                r.description = ""
        else:
            record.description = ""
    fileobject = StringIO()
    write_record(record, fileobject, fmt)
    return fileobject.getvalue().encode("utf-8")


def records_from_data_files(data_files, use_file_names_as_ids=False):
    records = []
    for file_ in data_files:
        circular = ("circular" not in file_) or file_.circular
        if file_.name.lower().endswith("zip"):
            records += records_from_zip_file(
                file_, use_file_names_as_ids=use_file_names_as_ids
            )
            continue
        recs, fmt = records_from_data_file(file_)
        single_record = len(recs) == 1
        for i, record in enumerate(recs):
            record.circular = circular
            record.linear = not circular
            name_no_extension = "".join(file_.name.split(".")[:-1])
            name = name_no_extension + ("" if single_record else ("%04d" % i))
            name = name.replace(" ", "_")
            UNKNOWN_IDS = [
                "None",
                "",
                "<unknown id>",
                ".",
                "EXPORTED",
                "<unknown name>",
                "Exported",
            ]
            # Sorry for this parts, it took a lot of "whatever works".
            # keep your part names under 20c and pointless, and everything
            # will be good
            if str(record.id).strip() in UNKNOWN_IDS:
                record.id = name
            if str(record.name).strip() in UNKNOWN_IDS:
                record.name = name
            record.file_name = name_no_extension
            if use_file_names_as_ids and single_record:
                basename = os.path.basename(record.source_file)
                basename_no_extension = os.path.splitext(basename)[0]
                record.id = basename_no_extension
        records += recs
    return records


def data_to_html_data(data, datatype, filename=None):
    """Data types: zip, genbank, fasta, pdf"""
    datatype = {
        "zip": "application/zip",
        "genbank": "application/genbank",
        "fasta": "application/fasta",
        "pdf": "application/pdf",
        "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    }.get(datatype, datatype)
    datatype = "data:%s;" % datatype
    data64 = "base64,%s" % b64encode(data).decode("utf-8")
    headers = ""
    if filename is not None:
        headers += "headers=filename%3D" + filename + ";"
    return datatype + headers + data64


def zip_data_to_html_data(data):
    return data_to_html_data(data, "application/zip")


LADDERS = {"100_to_4k": bw.ladders.LADDER_100_to_4k}


def matplotlib_figure_to_svg_base64_data(fig, **kwargs):
    """Return a string of the form 'data:image/svg+xml;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    output = BytesIO()
    fig.savefig(output, format="svg", **kwargs)
    svg_txt = output.getvalue().decode("utf-8")
    svg_txt = "\n".join(svg_txt.split("\n")[4:])
    svg_txt = "".join(svg_txt.split("\n"))

    content = b64encode(svg_txt.encode("utf-8"))
    result = (b"data:image/svg+xml;base64," + content).decode("utf-8")

    return result


def matplotlib_figure_to_bitmap_base64_data(fig, fmt="png", **kwargs):
    """Return a string of the form 'data:image/png;base64,XXX' where XXX
    is the base64-encoded svg version of the figure."""
    output = BytesIO()
    fig.savefig(output, format=fmt, **kwargs)
    bitmap = output.getvalue()
    content = b64encode(bitmap)
    result = (b"data:image/%s;base64,%s" % (fmt.encode("utf-8"), content)).decode(
        "utf-8"
    )
    return result


def figures_to_pdf_report_data(figures, filename="report.pdf"):
    pdf_io = BytesIO()
    with PdfPages(pdf_io) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
    return {
        "data": (
            "data:application/pdf;base64,"
            + b64encode(pdf_io.getvalue()).decode("utf-8")
        ),
        "name": filename,
        "mimetype": "application/pdf",
    }


def csv_to_list(csv_string, sep=","):
    return [
        element.strip()
        for line in csv_string.split("\n")
        for element in line.split(sep)
        if len(element.strip())
    ]


def set_record_topology(record, topology):
    """Set the Biopython record's topology, possibly passing if already set.

    This actually sets the ``record.annotations['topology']``.The ``topology``
    parameter can be "circular", "linear", "default_to_circular" (will default
    to circular if ``annotations['topology']`` is not already set) or
    "default_to_linear".
    """
    valid_topologies = [
        "circular",
        "linear",
        "default_to_circular",
        "default_to_linear",
    ]
    if topology not in valid_topologies:
        raise ValueError(
            "topology (%s) should be one of %s."
            % (topology, ", ".join(valid_topologies))
        )
    annotations = record.annotations
    default_prefix = "default_to_"
    if topology.startswith(default_prefix):
        if "topology" not in annotations:
            annotations["topology"] = topology[len(default_prefix) :]
    else:
        annotations["topology"] = topology


##############################################################################
def main():

    parser = argparse.ArgumentParser(description="Generate picklist for DNA assembly.")
    parser.add_argument("--parts_files", help="Directory with parts data or file with part sizes")
    parser.add_argument("--picklist", type=str, help="Path to the assembly plan CSV or Excel file")
    parser.add_argument("--source_plate", help="Source plate file (CSV or Excel)")
    parser.add_argument("--backbone_name", required=False, help="Name of the backbone")
    parser.add_argument("--result_zip", help="Name of the output zip file")
    parser.add_argument("--part_backbone_ratio", type=parse_optional_float, required=False, help="Part to backbone molar ratio")
    parser.add_argument("--quantity_unit", choices=["fmol", "nM", "ng"], help="Quantity unit")
    parser.add_argument("--part_quantity", type=float, help="Quantity of each part")
    parser.add_argument("--buffer_volume", type=float, help="Buffer volume in µL")
    parser.add_argument("--total_volume", type=float, help="Total reaction volume in µL")
    parser.add_argument("--dispenser", choices=["labcyte_echo", "tecan_evo"], help="Dispenser machine")

    args = parser.parse_args()

    # Parameters:
    picklist = args.picklist  # assembly plan
    # directory or can be a csv/Excel with part sizes
    if isinstance(args.parts_files, str):
        args.parts_files = args.parts_files.split(",")
    parts_dir = fix_and_rename_paths(args.parts_files)
    source_plate_path = args.source_plate
    backbone_name = args.backbone_name
    part_backbone_ratio = args.part_backbone_ratio
    result_zip_file = args.result_zip  # output file name "picklist.zip"
    ##############################################################################
    # Defaults:
    destination_plate = None
    destination_type = "new"  # this parameter is not actually used
    destination_size = 96  # this parameter is not actually used
    fill_by = "column"  # this parameter is not actually used
    quantity_unit = args.quantity_unit
    part_quantity = args.part_quantity # 1.3
    buffer_volume = args.buffer_volume # 0.3  # (µL)
    total_volume = args.total_volume # 1  # (µL)
    dispenser_machine = args.dispenser
    dispenser_min_volume = 0.5  # (nL), this parameter is not actually used
    dispenser_max_volume = 5  # (µL), this parameter is not actually used
    dispenser_resolution = 2.5  # (nL), this parameter is not actually used
    dispenser_dead_volume = 8  # (µL), this parameter is not actually used
    use_file_names_as_ids = True

    # CODE
    if picklist.endswith(".csv"):
        csv = picklist.read().decode()
        rows = [line.split(",") for line in csv.split("\n") if len(line)]
    else:
        dataframe = pandas.read_excel(picklist)
        rows = [row for i, row in dataframe.iterrows()]

    assembly_plan = AssemblyPlan(
        OrderedDict(
            [
                (
                    row[0],
                    [
                        str(e).strip()
                        for e in row[1:]
                        if str(e).strip() not in ["-", "nan", ""]
                    ],
                )
                for row in rows
                if row[0] not in ["nan", "Construct name", "constructs", "construct"]
            ]
        )
    )
    for assembly, parts in assembly_plan.assemblies.items():
        assembly_plan.assemblies[assembly] = [part.replace(" ", "_") for part in parts]

    # Reading part infos
    if not isinstance(parts_dir, list):
        if parts_dir.endswith((".csv", ".xls", ".xlsx")):  # part sizes specified in table
            if parts_dir.endswith(".csv"):
                dataframe = pandas.read_csv(parts_dir)
            else:
                dataframe = pandas.read_excel(parts_dir)
            parts_data = {row.part: {"size": row["size"]} for i, row in dataframe.iterrows()}
    else:  # input records
        records = dnacauldron.biotools.load_records_from_files(
            files=parts_dir, use_file_names_as_ids=use_file_names_as_ids
        )
        parts_data = {rec.id.replace(" ", "_").lower(): {"record": rec} for rec in records}
        #parts_data = process_parts_with_mapping(records, args.file_name_mapping)
    assembly_plan.parts_data = parts_data
    parts_without_data = assembly_plan.parts_without_data()
    if len(parts_without_data):
        print("success: False")
        print("message: Some parts have no provided record or data.")
        print("missing_parts: ", parts_without_data)
        sys.exit()
    # Reading protocol
    if quantity_unit == "fmol":
        part_mol = part_quantity * 1e-15
        part_g = None
    if quantity_unit == "nM":
        part_mol = part_quantity * total_volume * 1e-15
        part_g = None
    if quantity_unit == "ng":
        part_mol = None
        part_g = part_quantity * 1e-9
        # Backbone:part molar ratio calculation is not performed in this case.
        # This ensures no change regardless of form input:
        part_backbone_ratio = 1
    print("Generating picklist")
    picklist_generator = AssemblyPicklistGenerator(
        part_mol=part_mol,
        part_g=part_g,
        complement_to=total_volume * 1e-6,  # convert uL to L
        buffer_volume=buffer_volume * 1e-6,
        volume_rounding=2.5e-9,  # not using parameter from form
        minimal_dispense_volume=5e-9,  # Echo machine's minimum dispense -
    )
    if backbone_name != '' and backbone_name != 'Non':
        backbone_name_list = backbone_name.split(",")
    source_plate = plate_from_content_spreadsheet(source_plate_path)

    for well in source_plate.iter_wells():
        if well.is_empty:
            continue
        quantities = well.content.quantities
        part, quantity = list(quantities.items())[0]
        quantities.pop(part)
        quantities[part.replace(" ", "_")] = quantity

        if backbone_name != '' and backbone_name != 'Non':
            if part in backbone_name_list:
                # This section multiplies the backbone concentration with the
                # part:backbone molar ratio. This tricks the calculator into making
                # a picklist with the desired ratio.
                # For example, a part:backbone = 2:1 will multiply the
                # backbone concentration by 2, therefore half as much of it will be
                # added to the well.
                quantities[part.replace(" ", "_")] = quantity * part_backbone_ratio
            else:
                quantities[part.replace(" ", "_")] = quantity

    source_plate.name = "Source"
    if destination_plate:
        dest_filelike = file_to_filelike_object(destination_plate)
        destination_plate = plate_from_content_spreadsheet(destination_plate)
    else:
        destination_plate = Plate4ti0960("Mixplate")
    destination_wells = (
        well for well in destination_plate.iter_wells(direction="column") if well.is_empty
    )
    picklist, picklist_data = picklist_generator.make_picklist(
        assembly_plan,
        source_wells=source_plate.iter_wells(),
        destination_wells=destination_wells,
    )
    if picklist is None:
        print("success: False")
        print("message: Some parts in the assembly plan have no corresponding well.")
        print("picklist_data: ", picklist_data)
        print("missing_parts:", picklist_data.get("missing_parts", None))
        sys.exit()

    future_plates = picklist.simulate(inplace=False)


    def text(w):
        txt = human_volume(w.content.volume)
        if "construct" in w.data:
            txt = "\n".join([w.data["construct"], txt])
        return txt


    plotter = PlateTextPlotter(text)
    ax, _ = plotter.plot_plate(future_plates[destination_plate], figsize=(20, 8))

    ziproot = flametree.file_tree(result_zip_file, replace=True)

    # MIXPLATE MAP PLOT
    ax.figure.savefig(
        ziproot._file("final_mixplate.pdf").open("wb"),
        format="pdf",
        bbox_inches="tight",
    )
    plt.close(ax.figure)
    plate_to_platemap_spreadsheet(
        future_plates[destination_plate],
        lambda w: w.data.get("construct", ""),
        filepath=ziproot._file("final_mixplate.xls").open("wb"),
    )

    # ASSEMBLY REPORT
    print("Writing report...")
    picklist_to_assembly_mix_report(
        picklist,
        ziproot._file("assembly_mix_picklist_report.pdf").open("wb"),
        data=picklist_data,
    )
    assembly_plan.write_report(ziproot._file("assembly_plan_summary.pdf").open("wb"))

    # MACHINE PICKLIST

    if dispenser_machine == "labcyte_echo":
        picklist_to_labcyte_echo_picklist_file(
            picklist, ziproot._file("ECHO_picklist.csv").open("w")
        )
    else:
        picklist_to_tecan_evo_picklist_file(
            picklist, ziproot._file("EVO_picklist.gwl").open("w")
        )
    # We'll not write the input source plate.
    # raw = file_to_filelike_object(source_plate_path).read()
    # f = ziproot.copy(source_plate_path)
    # f.write(raw, mode="wb")
    ziproot._close()
    print("success: True")
    

if __name__ == "__main__":
    main()
