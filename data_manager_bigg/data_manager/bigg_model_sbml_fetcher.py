import argparse
import json
import os
import sys
import time
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen


MODEL_URL = "http://bigg.ucsd.edu/static/models/"
MODEL_DETAIL_URL = "http://bigg.ucsd.edu/api/v2/models/"


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


def url_json(url):
    data = {}
    try:
        with urlopen(Request(url)) as fod:
            data = fod.read().decode("utf-8")
        data = json.loads(data)
    except Exception as e:
        sys.exit(str(e))
    return data


def download_entries(model_ids, id2org, workdir):
    for ix, model_id in enumerate(model_ids):
        model_filename = model_id + ".xml"
        path = os.path.abspath(os.path.join(workdir, model_filename))

        url_download(MODEL_URL + model_filename, path)
        data_manager_entry = {}
        data_manager_entry["value"] = model_id
        data_manager_entry["name"] = model_id + " - " + id2org[model_id]
        data_manager_entry["path"] = path

        # Make sure that less than 10 requests per second, as required by host (http://bigg.ucsd.edu/data_access)
        if ix % 5 == 0:
            time.sleep(1)
        yield data_manager_entry


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    pinput = parser.add_mutually_exclusive_group(required=True)
    pinput.add_argument("--model-id", help="Model BIGG id")
    pinput.add_argument("--model-all", action="store_true", help="Download all models")
    parser.add_argument("--out-file", help="JSON output file")
    args = parser.parse_args()

    # Init.
    data_manager_json = {"data_tables": {}}
    with open(args.out_file) as fh:
        params = json.load(fh)

    workdir = params["output_data"][0]["extra_files_path"]
    os.makedirs(workdir)

    # Load models and models metadata.
    models = url_json(MODEL_DETAIL_URL)
    id2org = {}
    for result in models.get("results", []):
        ident = result["bigg_id"]
        ident = ident.replace(" ", "")
        id2org[ident] = result["organism"]

    # Select model_ids.
    model_ids = []
    if args.model_id:
        if args.model_id not in id2org.keys():
            sys.exit("Model id: %s is not available" % (args.model_id,))
        model_ids.append(args.model_id)
    else:
        model_ids.extend(list(id2org.keys()))

    # Download.
    entries = list(download_entries(model_ids, id2org, workdir))

    # Write data.
    data_manager_json["data_tables"]["bigg_model_sbml"] = entries
    with open(args.out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)
