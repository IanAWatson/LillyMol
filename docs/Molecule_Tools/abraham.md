# abraham

Computes additive Abraham solvation descriptors, McGowan volume (`abr_Vx`), and
a derived `abr_Logp` value.

```shell
abraham.sh molecules.smi > molecules.abraham
```

For a direct invocation, supply the model files:

```shell
abraham -l -g all -F "$LILLYMOL_HOME/data/queries/abraham/Abraham" \
  -P "$LILLYMOL_HOME/data/queries/abraham/Alpha2H" molecules.smi
```

Descriptor values use platform-independent formatting with four decimal places
and trailing zeros omitted. Negative values and values above one use the same
formatting. Calculations retain their original precision; rounding affects only
the text output. Fingerprint generation (`-J`) retains its existing encoding.

`-Y wpipe` includes SMILES in descriptor output. `-Y rpipe` reads a headered
SMILES/identifier/descriptor table, appending Abraham descriptors while preserving
existing descriptor columns. Combine both options to retain the SMILES column.
