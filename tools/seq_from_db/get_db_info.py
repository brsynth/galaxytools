import subprocess
import argparse
import time
import os
import socket
import re
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
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


def fetch_annotations(csv_file, sequence_column, annotation_columns, db_uri, table_name, fragment_column_name, output):
    """Fetch annotations from the database and save the result as GenBank files."""
    db_uri = fix_db_uri(db_uri)
    df = pd.read_csv(csv_file, sep=',', header=None)

    engine = create_engine(db_uri)
    connection = engine.connect()

    annotated_data = []

    try:
        with connection:
            inspector = inspect(engine)
            columns = [column['name'] for column in inspector.get_columns(table_name)]

            # Fetch all fragments from the table once
            if fragment_column_name not in columns:
                raise ValueError(f"Fragment column '{fragment_column_name}' not found in table '{table_name}'.")

            fragment_column_index = columns.index(fragment_column_name)
            all_rows = connection.execute(text(f"SELECT * FROM {table_name}")).fetchall()
            fragment_map = {row[fragment_column_index]: row for row in all_rows}

            # Compare fragments between CSV and DB
            csv_fragments = set()
            all_ids = set(df[0].dropna().astype(str))
            for _, row in df.iterrows():
                for col in df.columns:
                    if col != 0:
                        fragment = row[col]
                        if pd.notna(fragment):
                            fragment_str = str(fragment)
                            if fragment_str not in all_ids:
                                csv_fragments.add(fragment_str)

            db_fragments = set(fragment_map.keys())
            missing_fragments = sorted(list(csv_fragments - db_fragments))
            if missing_fragments:
                raise ValueError(
                    f" Missing fragments in DB: {', '.join(missing_fragments)}"
                )

            # === CONTINUE WITH GB FILE CREATION ===
            for _, row in df.iterrows():
                annotated_row = {"Backbone": row[0], "Fragments": []}
                for col in df.columns:
                    if col != 0:
                        fragment = row[col]
                        if fragment not in csv_fragments:
                            continue
                        db_row = fragment_map.get(fragment)

                        if db_row:
                            fragment_data = {"id": fragment}
                            for i, column_name in enumerate(columns[1:]):  # skip ID column
                                fragment_data[column_name] = db_row[i + 1]
                        else:
                            fragment_data = {"id": fragment, "metadata": "No data found"}

                        annotated_row["Fragments"].append(fragment_data)

                annotated_data.append(annotated_row)

    except Exception as e:
        print(f"Error occurred during annotation: {e}")
        raise  # Ensures the error exits the script

    # GenBank file generation per fragment
    try:
        for annotated_row in annotated_data:
            backbone_id = annotated_row["Backbone"]
            for fragment in annotated_row["Fragments"]:
                fragment_id = fragment["id"]
                sequence = fragment.get(sequence_column, "")
                annotation = fragment.get(annotation_columns, "")

                # Create the SeqRecord
                record = SeqRecord(
                    Seq(sequence),
                    id=fragment_id,
                    name=fragment_id,
                    description=f"Fragment {fragment_id} from Backbone {backbone_id}"
                )

                # Add annotations to GenBank header
                record.annotations = {
                    k: str(fragment[k]) for k in annotation_columns if k in fragment
                }

                # LOCUS line extraction from annotation (copy-paste the LOCUS from annotation)
                locus_line_match = re.search(r"LOCUS\s+.+", annotation)
                if locus_line_match:
                    locus_line = locus_line_match.group()
                else:
                    print(f"LOCUS info missing for fragment {fragment_id}")
                    locus_line = f"LOCUS       {fragment_id: <20} {len(sequence)} bp    DNA     linear   UNK 01-JAN-2025"

                # Format sequence as per GenBank standards (with ORIGIN and line breaks)
                if "ORIGIN" in sequence:
                    origin_block = sequence.strip()
                else:
                    # Format sequence as per GenBank standards (with ORIGIN and line breaks)
                    formatted_sequence = "ORIGIN\n"
                    seq_str = str(record.seq)
                    for i in range(0, len(seq_str), 60):  # 60 bases per line
                        line_seq = seq_str[i:i + 60]
                        formatted_sequence += f"{str(i + 1).rjust(9)} { ' '.join([line_seq[j:j+10] for j in range(0, len(line_seq), 10)]) }\n"
                    origin_block = formatted_sequence.strip()

                # Find and copy the FEATURES section directly from annotation
                features_section = ""
                features_start = annotation.find("FEATURES")
                if features_start != -1:
                    features_section = annotation[features_start:]

                # Writing the GenBank file
                if not os.path.exists(output):
                    os.makedirs(output)

                gb_filename = os.path.join(output, f"{fragment_id}.gb")
                with open(gb_filename, "w") as f:
                    # Write the LOCUS line
                    f.write(locus_line + "\n")
                    # Write DEFINITION, ACCESSION, and other annotations
                    f.write(f"DEFINITION  {record.description}\n")
                    f.write(f"ACCESSION   {record.id}\n")
                    f.write(f"VERSION     DB\n")
                    f.write(f"KEYWORDS    .\n")
                    f.write(f"SOURCE      .\n")
                    # Write the FEATURES section directly from annotation
                    f.write(features_section)
                    # Write the ORIGIN section
                    f.write(origin_block + "\n")
                    f.write("//\n")

    except Exception as e:
        print(f"Error saving GenBank files: {e}")
        return


def main():
    parser = argparse.ArgumentParser(description="Fetch annotations from PostgreSQL database and save as JSON.")
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--sequence_column", required=True, help="DB column contains sequence for ganbank file")
    parser.add_argument("--annotation_columns", required=True, help="DB column contains head for ganbank file")
    parser.add_argument("--db_uri", required=True, help="Database URI connection string")
    parser.add_argument("--table", required=True, help="Table name in the database")
    parser.add_argument("--fragment_column", required=True, help="Fragment column name in the database")
    parser.add_argument("--output", required=True, help="Output dir for gb files")
    args = parser.parse_args()

    # Wait until the database is ready
    db_uri = fix_db_uri(args.db_uri)
    # db_name = extract_db_name(db_uri)
    # start_postgres_container(db_name)
    MAX_RETRIES = 3
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            wait_for_db(db_uri)
            break  # Success
        except Exception as e:
            if attempt == MAX_RETRIES:
                print(f"Attempt {attempt} failed: Could not connect to database at {db_uri}.")
                raise e
            else:
                time.sleep(2)

    # Fetch annotations from the database and save as gb
    fetch_annotations(args.input, args.sequence_column, args.annotation_columns, db_uri, args.table, args.fragment_column, args.output)


if __name__ == "__main__":
    main()
