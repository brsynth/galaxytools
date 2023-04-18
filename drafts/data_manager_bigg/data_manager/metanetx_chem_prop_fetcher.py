import argparse
import json
import os
import shutil
import sys
import tempfile
# import pandas as pd
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen


METANETX_URL = "https://www.metanetx.org/ftp/"


def url_download(url, path):
    try:
        with urlopen(Request(url)) as fod:
            with open(path, "wb") as dst:
                while True:
                    chunk = fod.read(2**10)
                    if chunk:
                        dst.write(chunk)
                    else:
                        break
    except Exception as e:
        sys.exit(str(e))


def clean_metanetx_file(path):
    ftmp = tempfile.NamedTemporaryFile()
    isHeaderFound = False
    with open(path) as fid, open(ftmp.name, 'w') as fod:
        for line in fid:
            if line.startswith("#"):
                last_line = line
            else:
                if not isHeaderFound:
                    last_line = last_line.replace("#", "")
                    fod.write(last_line)
                    isHeaderFound = True
                fod.write(line)
    shutil.copyfile(ftmp.name, path)


# def records_chem_prop_pandas(path):
#    df = pd.read_csv(path, sep="\t")
#    df["name"] = df.apply(lambda x: "%s: %s (%s)" % (x["ID"], x["name"], x["formula"]), axis=1)
#    df.drop(columns=["reference", "formula", "charge", "mass", "InChIKey", "SMILES"], inplace=True)
#    df.rename(columns={"ID": "value", "InChI": "inchi"})
#    return df.to_dict('records')


def records_chem_prop(path):
    records = []
    with open(path) as fid:
        for ix, line in enumerate(fid):
            if ix == 0:
                continue
            line = line.split("\t")
            if line[7] != '':
                records.append({
                    "value": line[0],
                    "name": "%s: %s (%s)" % (line[0], line[1], line[3]),
                    "inchi": line[7],
                })
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    pinput = parser.add_mutually_exclusive_group(required=True)
    pinput.add_argument("--version", help="Version to download")
    parser.add_argument("--out-file", help="JSON output file")
    args = parser.parse_args()

    # Init.
    data_manager_json = {"data_tables": {}}
    with open(args.out_file) as fh:
        params = json.load(fh)

    workdir = params["output_data"][0]["extra_files_path"]
    os.makedirs(workdir)

    # Load models and models metadata.
    ftmp = tempfile.NamedTemporaryFile()
    url = '/'.join([METANETX_URL, args.version, 'chem_prop.tsv'])
    url_download(url, ftmp.name)

    # Clean header
    clean_metanetx_file(ftmp.name)

    # Select records.
    records = records_chem_prop(ftmp.name)

    # Write data.
    data_manager_json["data_tables"]["metanetx_chem_prop"] = records
    with open(args.out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)
