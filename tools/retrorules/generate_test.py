import hashlib
import subprocess
import tempfile
from typing import List


def compute_md5(path: str):
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def run_cmd(subcommand: str, sub_cmd: List):
    with tempfile.NamedTemporaryFile() as fd:
        cmd = ["python", "query.py"]
        cmd.append(subcommand)
        cmd.append("--output-data-json")
        cmd.append(fd.name)
        cmd.extend(sub_cmd)

        print("Cmd:", " ".join(cmd))
        ret = subprocess.run(cmd)
        assert ret.returncode == 0

        return compute_md5(path=fd.name)


if __name__ == "__main__":

    print("=" * 50)
    print("EC number")
    ec_number = "1.1.1.1"
    cmd = ["--input-ec-number-str", ec_number]
    value = run_cmd(subcommand="ec-number", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd = ["--input-ec-number-str", ec_number, "--input-min-diameter-int", "16"]
    value = run_cmd(subcommand="ec-number", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    print("Substrate")
    substrate = "pyruvate"
    cmd = ["--input-substrate-str", substrate]
    value = run_cmd(subcommand="substrate", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd = ["--input-substrate-str", substrate, "--input-min-diameter-int", "16"]
    value = run_cmd(subcommand="substrate", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    print("Reaction id")
    reaction_id = "MNXR104443"
    repository = "mnx"
    cmd = ["--input-reaction-id-str", reaction_id, "--input-repository-str", repository]
    value = run_cmd(subcommand="reaction-id", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd.extend(["--input-min-diameter-int", "1"])
    value = run_cmd(subcommand="reaction-id", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    print("Inchi")
    inchi = "InChI=1S/C4H4O5/c5-2(4(8)9)1-3(6)7/h1H2,(H,6,7)(H,8,9)/p-2"
    cmd = ["--input-inchi-str", inchi]
    value = run_cmd(subcommand="inchi", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd.extend(["--input-min-diameter-int", "1"])
    value = run_cmd(subcommand="inchi", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    print("Repository")
    repository = "kegg"
    cmd = ["--input-repository-str", repository]
    value = run_cmd(subcommand="repository", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd.extend(["--input-min-diameter-int", "16"])
    value = run_cmd(subcommand="repository", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    print("Smarts id")
    smarts_id = "RR00239878"
    cmd = ["--input-smarts-id-str", smarts_id]
    value = run_cmd(subcommand="smarts-id", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    cmd.extend(["--input-min-diameter-int", "16"])
    value = run_cmd(subcommand="smarts-id", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    smarts_ids = ["RR00239877", "RR00239878"]
    cmd = ["--input-smarts-id-str"] + smarts_ids
    value = run_cmd(subcommand="smarts-id", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
