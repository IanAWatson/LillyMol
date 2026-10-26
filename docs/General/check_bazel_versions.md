# check_bazel_versions

`contrib/bin/check_bazel_versions.py` reports whether the versions declared in a
`MODULE.bazel` file are the latest eligible versions in the Bazel Central
Registry. It uses only the Python standard library; LillyMol Python bindings
and Bazel itself are not required.

From `src`, check against a local registry clone:

```shell
python3 ../contrib/bin/check_bazel_versions.py \
  --registry /d1/ian/bazel-central-registry
```

Pass a module file explicitly to run from another directory:

```shell
python3 /path/to/LillyMol/contrib/bin/check_bazel_versions.py \
  /path/to/project/MODULE.bazel --registry /path/to/bazel-central-registry
```

A local clone is optional. Omit `--registry` to read metadata directly from
`https://bcr.bazel.build` over HTTPS. `--timeout SECONDS` controls each request's
timeout (default 15 seconds).

The script prints the declared version, latest eligible version and comparison
status for each dependency. Newer versions are informational: the script never
edits files, runs Bazel, updates the clone or resolves transitive dependencies.
A local report reflects only the metadata currently present in the clone.

## Version selection

Versions are compared using Bazel's relaxed SemVer ordering, including date
versions, `.bcr.N` registry patches, numeric segments and prerelease identifiers.
Build metadata after `+` does not affect ordering. Numeric segments are not
sorted as strings, so `1.10` follows `1.9`.

Yanked versions are excluded from the latest-version candidates. A yanked
currently declared version is reported with its reason. Missing declared
versions and module deprecation notices are also shown.

By default, prereleases are excluded. Use `--include-prereleases` to consider
them. Prerelease detection follows Bazel's hyphen syntax with one deliberate
exception: valid ISO date releases such as `2025-11-05` and
`2025-11-05.bcr.1` are treated as releases. This accommodates RE2's naming
scheme. Other unusual upstream preview naming conventions are not inferred.

## Declaration handling and status

Only direct, top-level `bazel_dep` calls with literal string names and versions
are checked. Comments, embedded BUILD-file strings and repository aliases do
not create dependencies. The script parses the declaration syntax without
executing Starlark. Computed names/versions and `include()` produce an explicit
error; included files can be checked separately.

Dependencies with local, Git, archive, single-version or multiple-version
overrides are shown as skipped, along with the override type. Versionless
requirements are also skipped. Override-selected versions are not compared.

Statuses include `CURRENT`, `UPDATE_AVAILABLE`, `DECLARED_NEWER`, `YANKED`,
`NO_ELIGIBLE_VERSION`, `SKIPPED` and `ERROR`. Finding newer or yanked versions
does not make the command fail. A parsing, metadata or retrieval error produces
exit status 1; a completed check returns 0. Errors for individual modules are
reported while checking the remaining dependencies.

A newer version may require build or source changes. This report compares the
declared versions only; Bazel's resolved dependency graph can select different
versions through transitive requirements.

## Tests

```shell
python3 contrib/bin/test/check_bazel_versions_test.py
```

Run from the repository root. These tests use fixture registries and mocked
HTTPS responses; they do not require network access.
