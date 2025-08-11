import subprocess
import argparse
import time
import os
import socket
import re
import json
from sqlalchemy import create_engine, inspect
from sqlalchemy.engine.url import make_url
from sqlalchemy.sql import text
from sqlalchemy.exc import OperationalError


def fix_db_uri(uri):
    """Replace __at__ with @ in the URI if needed."""
    return uri.replace("__at__", "@")


def is_port_in_use(uri):
    """Check if a TCP port is already in use on host."""
    url = make_url(uri)
    host = url.host
    port = url.port
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(2)
        return s.connect_ex((host, port)) == 0

def extract_db_name(uri):
    """Extract the database name from the SQLAlchemy URI."""
    url = make_url(uri)
    return url.database


# this fuction is to activate the Docker id the DB is in container. BUT IT IS NOT USED IN MAIN()
def start_postgres_container(db_name):
    """Start a PostgreSQL container with the given database name as the container name."""
    container_name = db_name

    # Check if container is already running
    container_running = subprocess.run(
        f"docker ps -q -f name={container_name}", shell=True, capture_output=True, text=True
    )

    if container_running.stdout.strip():
        print(f"Container '{container_name}' is already running.")
        return

    # Check if container exists (stopped)
    container_exists = subprocess.run(
        f"docker ps -a -q -f name={container_name}", shell=True, capture_output=True, text=True
    )

    if container_exists.stdout.strip():
        print(f"Starting existing container '{container_name}'...")
        subprocess.run(f"docker start {container_name}", shell=True)
        print(f"PostgreSQL Docker container '{container_name}' activated.")
        return

    # If container does not exist, create and start a new one
    port = 5432 if not is_port_in_use(5432) else 5433
    postgres_password = os.getenv("POSTGRES_PASSWORD", "RK17")

    start_command = [
        "docker", "run", "--name", container_name,
        "-e", f"POSTGRES_PASSWORD={postgres_password}",
        "-p", f"{port}:5432",
        "-d", "postgres"
    ]

    try:
        subprocess.run(start_command, check=True)
        print(f"PostgreSQL Docker container '{container_name}' started on port {port}.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to start Docker container: {e}")


def wait_for_db(uri, timeout=60):
    """Try connecting to the DB until it works or timeout."""
    engine = create_engine(uri)
    start_time = time.time()
    while time.time() - start_time < timeout:
        try:
            with engine.connect():
                print("Connected to database.")
                return
        except OperationalError:
            print("Database not ready, retrying...")
            time.sleep(2)
    raise Exception("Database connection failed after timeout.")


def push_gb_annotations(gb_files, sequence_column, annotation_column, db_uri, table_name, fragment_column_name, output, file_name_mapping):
    """Push GenBank file content into the database if the fragment is not already present."""
    db_uri = fix_db_uri(db_uri)
    engine = create_engine(db_uri)
    inserted_fragments = []

    try:
        # Parse the file_name_mapping string into a dictionary {base_file_name: fragment_name}
        file_name_mapping_dict = {
            os.path.basename(path): os.path.splitext(fragment_name)[0]
            for mapping in file_name_mapping.split(",")
            for path, fragment_name in [mapping.split(":")]
        }

        #print("File name mapping dictionary:")
        #print(file_name_mapping_dict)  # Debugging: Print the mapping dictionary

        with engine.begin() as connection:
            inspector = inspect(engine)
            columns = [col['name'] for col in inspector.get_columns(table_name)]

            if fragment_column_name not in columns:
                raise ValueError(f"Fragment column '{fragment_column_name}' not found in table '{table_name}'.")

            # Get existing fragments
            all_rows = connection.execute(text(f"SELECT {fragment_column_name} FROM {table_name}")).fetchall()
            existing_fragments = {row[0] for row in all_rows}

            insert_rows = []

            for gb_file in gb_files:
                # Extract base file name (just the file name, not the full path)
                real_file_name = os.path.basename(gb_file)
                fragment_name = file_name_mapping_dict.get(real_file_name)

                print(f"Processing file: {real_file_name}({fragment_name})")  # Debugging: Log the current file

                # Get the corresponding fragment name from the mapping
                fragment_name = file_name_mapping_dict.get(real_file_name)

                if not fragment_name:
                    raise ValueError(f"Fragment name not found for file '{real_file_name}' in file_name_mapping.")

                # If the fragment is already in the DB, raise an error and stop the process
                if fragment_name in existing_fragments:
                    raise RuntimeError(f"Fatal Error: Fragment '{fragment_name}' already exists in DB. Stopping the process.")

                with open(gb_file, "r") as f:
                    content = f.read()

                origin_match = re.search(r"^ORIGIN.*$", content, flags=re.MULTILINE)
                if not origin_match:
                    raise ValueError(f"ORIGIN section not found in file: {gb_file}")

                origin_start = origin_match.start()
                annotation_text = content[:origin_start].strip()
                sequence_text = content[origin_start:].strip()

                values = {}
                values[fragment_column_name] = fragment_name
                values[annotation_column] = annotation_text
                values[sequence_column] = sequence_text

                insert_rows.append(values)
                inserted_fragments.append(fragment_name)

            # Insert the rows into the database
            for values in insert_rows:
                col_names = ", ".join(values.keys())
                placeholders = ", ".join([f":{key}" for key in values.keys()])
                insert_stmt = text(f"INSERT INTO {table_name} ({col_names}) VALUES ({placeholders})")

                # print(f"Inserting into DB: {values}")  # Debugging print statement
                connection.execute(insert_stmt, values)

                # print(f"Insert result: {result.rowcount if hasattr(result, 'rowcount') else 'N/A'}")  # Debugging the row count

            print(f"Inserted {len(insert_rows)} fragments.")

            # Write inserted fragment names to a text file
            with open(output, "w") as log_file:
                for frag in inserted_fragments:
                    log_file.write(f"{frag}\n")
            print(f"Fragment names written to '{output}'.")

    except Exception as e:
        print(f"Error during GB file insertion: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(description="Fetch annotations from PostgreSQL database and save as JSON.")
    parser.add_argument("--input", required=True, help="Input gb files")
    parser.add_argument("--use_json_paramers", required=False, help="Use parameters from JSON: true/false")
    parser.add_argument("--sequence_column", required=False, help="DB column contains sequence for GenBank file")
    parser.add_argument("--annotation_column", required=False, help="DB column contains head for GenBank file")
    parser.add_argument("--db_uri", required=False, help="Database URI connection string")
    parser.add_argument("--table", required=False, help="Table name in the database")
    parser.add_argument("--fragment_column", required=False, help="Fragment column name in the database")
    parser.add_argument("--output", required=False, help="Text report")
    parser.add_argument("--file_name_mapping", required=False, help="Real fragment names")
    parser.add_argument("--json_conf", required=False, help="JSON config file with DB parameters")
    parser.add_argument("--execution_enable", required=False, help="Enable or disable execution")
    parser.add_argument("--json_generating", required=False, help="Generate JSON: true/false")
    parser.add_argument("--json_output", required=False, help="Output path for generated JSON")

    args = parser.parse_args()

    if args.execution_enable == 'false':
        print("Execution disabled. 'Send Request to DB' is set to 'false'")
        return

    config_params = {}
    use_json = args.use_json_paramers == 'true'
    generate_json = args.json_generating == 'true'

    if use_json:
        if not args.json_conf:
            raise ValueError("You must provide --json_conf when --use_json_paramers is 'true'")
        with open(args.json_conf, "r") as f:
            config_params = json.load(f)
        if config_params.get("execution", "") == "false":
            print("Execution blocked by config (execution = false)")
            return
    else:
        config_params = {
            "table": args.table,
            "sequence_column": args.sequence_column,
            "annotation_column": args.annotation_column,
            "fragment_column": args.fragment_column,
            "db_uri": args.db_uri,
            "execution": args.execution_enable
        }

    if generate_json:
        if not args.json_output:
            raise ValueError("You must provide --json_output when --json_generating is 'true'")
        with open(args.json_output, "w") as f:
            json.dump(config_params, f, indent=2)
        print(f"JSON configuration written to: {args.json_output}")

    # Extract final resolved parameters
    table = config_params["table"]
    sequence_column = config_params["sequence_column"]
    annotation_column = config_params["annotation_column"]
    fragment_column = config_params["fragment_column"]
    db_uri = fix_db_uri(config_params["db_uri"])

    gb_file_list = [f.strip() for f in args.input.split(",") if f.strip()]

    # Connect to DB
    MAX_RETRIES = 3
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            wait_for_db(db_uri)
            break
        except Exception as e:
            if attempt == MAX_RETRIES:
                print(f"Attempt {attempt} failed: Could not connect to database at {db_uri}.")
                raise e
            time.sleep(2)

    # Push annotations
    push_gb_annotations(
        gb_file_list,
        sequence_column,
        annotation_column,
        db_uri,
        table,
        fragment_column,
        args.output,
        args.file_name_mapping
    )


if __name__ == "__main__":
    main()
