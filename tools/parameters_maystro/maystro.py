import argparse
import tempfile
import os
import json
import shutil


def parse_command_line_args():
    parser = argparse.ArgumentParser(description="Maystro JSON Handler")

    parser.add_argument("--distribute_json", required=True, help="true or false")
    parser.add_argument("--json_from_workflow", required=False, nargs='+', help="JSON files from tools", default=[])
    parser.add_argument("--json_from_user", required=False, help="User-provided JSON")
    parser.add_argument("--json_name_mapping", required=True, help="map the real json name")
    parser.add_argument("--output_workflow", required=True, help="JSON output for next workflow steps")
    parser.add_argument("--output_user", required=True, help="Final JSON output to user")

    return parser.parse_args()


def parse_file_name_mapping(mapping_str):
    mapping = {}
    if mapping_str:
        for pair in mapping_str.split(','):
            stored, original = pair.strip().split(':', 1)
            # Strip .json from original
            real_name = os.path.splitext(original)[0]
            mapping[os.path.basename(stored)] = real_name
    return mapping


def handle_distribute_json_false(args):
    temp_dir = tempfile.mkdtemp(prefix="maystro_merge_")
    print(f"[INFO] Watching temp dir for new JSONs: {temp_dir}")

    try:
        # Collect JSONs from json_from_workflow
        initial_jsons = list(filter(os.path.isfile, args.json_from_workflow))
        print(f"[INFO] Initial JSONs from workflow: {initial_jsons}")

        # Parse filename mapping if provided
        filename_mapping = parse_file_name_mapping(getattr(args, 'json_name_mapping', ''))

        # Merge all together
        merged = {}
        for file_path in initial_jsons:
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                basename = os.path.basename(file_path)
                real_name = filename_mapping.get(basename, basename)  # fallback if not in mapping
                merged[real_name] = data
                print(f"[INFO] Added data under key: {real_name}")
            except json.JSONDecodeError as e:
                print(f"[WARN] Skipping invalid JSON file {file_path}: {e}")

        with open(args.output_user, "w") as f:
            json.dump(merged, f, indent=2)
        print(f"[INFO] Merged JSON written to: {args.output_user}")

    finally:
        print(f"[INFO] Cleaning up: {temp_dir}")
        shutil.rmtree(temp_dir)


def merge_json_files(paths):
    merged = {}
    for path in paths:
        try:
            with open(path, "r") as f:
                data = json.load(f)
                merged.update(data)
        except Exception as e:
            print(f"[WARN] Skipping {path}: {e}")
    return merged


def handle_distribute_json_true(args):
    if not args.json_from_user:
        raise ValueError("json_from_user is required when distribute_json is true")

    with open(args.json_from_user, 'r') as in_f:
        user_data = json.load(in_f)

    with open(args.output_workflow, 'w') as out_f:
        json.dump(user_data, out_f, indent=2)


def main():
    args = parse_command_line_args()

    if args.distribute_json.lower() == 'false':
        handle_distribute_json_false(args)
    else:
        handle_distribute_json_true(args)

if __name__ == "__main__":
    main()
