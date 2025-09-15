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
        if ret.returncode != 0:
            msg = f"Command failed: {' '.join(cmd)}\n"
            if ret.stdout:
                msg += "Stdout:", ret.stdout + "\n"
            if ret.stderr:
                msg += "Stderr:", ret.stderr + "\n"
            raise Exception()

        return compute_md5(path=fd.name)


if __name__ == "__main__":

    print("=" * 50)
    print("templates")
    smarts_str = "[O]-[C](=[O])"
    cmd = ["--input-smarts-str", smarts_str, "--input-limit-int", "5"]
    value = run_cmd(subcommand="templates", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)

    print("templates-summary")
    template_id = "RR:03-27BC85-19184A-A71018"
    cmd = ["--input-template-id-str", template_id]
    value = run_cmd(subcommand="templates-summary", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    
    print("=" * 50)

    print("templates-sources")
    template_id = "RR:03-27BC85-19184A-A71018"
    cmd = ["--input-template-id-str", template_id]
    value = run_cmd(subcommand="templates-sources", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    
    print("=" * 50)
    
    print("templates-count")
    smarts_str = "[O]-[C](=[O])"
    cmd = ["--input-smarts-str", smarts_str]
    value = run_cmd(subcommand="templates", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)

    print("=" * 50)
    """
    print("templates-export")
    smarts_str = "[O]-[C](=[O])"
    cmd = ["--input-smarts-str", smarts_str]
    value = run_cmd(subcommand="templates-export", sub_cmd=cmd)
    print("Test:", " ".join(cmd), "md5:", value)
    """