# Python reaction examples

Each script is standalone and checks its expected products. Its reaction and
reagent files are located relative to the script, so it can be run from any
working directory. Use a virtual environment matching the installed LillyMol
Python ABI, and run via the repository wrapper:

```shell
"$LILLYMOL_HOME/run_python.sh" "$LILLYMOL_HOME/src/nanobind/examples/reactions/in_place.py"
```

| Script | Workflow |
| --- | --- |
| `in_place.py` | Change all chlorine sites in the input molecule. |
| `per_embedding.py` | Generate a separate product at each scaffold match. |
| `explicit_reagents.py` | Supply scaffold and sidechain molecules directly. |
| `enumerate_reagents.py` | Enumerate multiple scaffolds against stored reagents. |

See the [reaction guide](../../../../docs/python/reactions.md) for explanations
and the API reference. Its full example listings are checked against these files.
Each script is also a Bazel test, named
`//nanobind:reaction_<script-name>_example_test`.
