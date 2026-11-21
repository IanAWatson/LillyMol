#!/usr/bin/env python3
"""Exercise CLI file safety and error propagation without Python bindings."""
from pathlib import Path
import os
import subprocess
import sys

exe = Path(sys.argv[1]).resolve()


def run(options, inputs=("input.smi",), success=True):
    result = subprocess.run([str(exe), *options, *map(str, inputs)],
                            text=True, capture_output=True, timeout=30)
    assert (result.returncode == 0) == success, (options, result.returncode, result.stderr)
    return result


source = Path("input.smi")
original = "CC example\n"
source.write_text(original)
# No output is opened until all destinations are checked.
for options in (["-S", "input"], ["-Y", "B4F=input.smi", "-S", "-"],
                ["-L", "input.smi", "-S", "-"],
                ["-Y", "B4F=auxiliary", "-S", "input"]):
    run(options, success=False)
    assert source.read_text() == original
assert not Path("auxiliary.smi").exists()

# Protect raw rejection-log names as well as molecule outputs.
raw = Path("input.txt")
raw.write_text(original)
run(["-i", "smi", "-L", "input.txt", "-S", "-"], [raw], success=False)
assert raw.read_text() == original

# F: file lists must not hide an input from overwrite checks.
Path("inputs.list").write_text("# comment\n\ninput.smi\n")
run(["-S", "input"], ["F:inputs.list"], success=False)
assert source.read_text() == original

# Equivalent paths, symlinks, and hard links must be treated as the same file.
run(["-S", "./input"], success=False)
Path("alias.smi").symlink_to(source.resolve())
run(["-S", "alias"], success=False)
os.link(source, "hardlink.smi")
run(["-S", "hardlink"], success=False)
assert source.read_text() == original

# A numbered output stem must not overwrite the current input either.
Path("numbered0.smi").write_text(original)
run(["-G", "-S", "numbered"], ["numbered0.smi"], success=False)
assert Path("numbered0.smi").read_text() == original

# Converting one input must not overwrite a later input of another type.
run(["-o", "sdf", "-S", "pair"])
Path("pair.smi").write_text(original)
run([], ["pair.sdf", "pair.smi"], success=False)
assert Path("pair.smi").read_text() == original

# Named auxiliary and main outputs must also be distinct from each other.
run(["-Y", "B4F=collision", "-S", "collision"], success=False)
assert not Path("collision.smi").exists()

for options in (["-B", "bad"], ["-Y", "rpt=bad"], ["-Y", "rpt=0"],
                ["-Y", "maxplen=bad"]):
    run([*options, "-S", "-"], success=False)

run(["-Y", "rpt=1", "-S", "-"])
run(["-Y", "B4F=missing/destination", "-S", "-"], success=False)

# Rejected structures and explanations must use different files.
run(["-c", "3", "-L", "rejected.smi", "-S", "-"])
assert Path("rejected.smi").read_text().startswith("CC example REJ ")
assert "REASON too few atoms" in Path("rejected.smi.log").read_text()
run(["-c", "3", "-L", "rejects", "-S", "-"])
assert Path("rejects.smi").read_text().startswith("CC example REJ ")
assert "REASON too few atoms" in Path("rejects").read_text()
assert run(["-c", "3", "-L", "-", "-S", "accepted"]).stdout.startswith("CC example REJ ")

# A failed flush must affect the exit status, including small buffered outputs.
if Path("/dev/full").exists():
    Path("full.smi").symlink_to("/dev/full")
    run(["-Y", "B4F=full.smi", "-S", "-"], success=False)
    run(["-S", "full"], success=False)
    with open("/dev/full", "w") as output:
        result = subprocess.run([str(exe), "-S", "-", str(source)],
                                stdout=output, stderr=subprocess.PIPE, text=True, timeout=30)
    assert result.returncode != 0, result.stderr
    Path("reject_full.smi").symlink_to("/dev/full")
    run(["-c", "3", "-L", "reject_full.smi", "-S", "-"], success=False)

print("Fileconv correctness regressions passed")
