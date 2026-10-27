"""Check generated script protos without loading LillyMol's native bindings."""
import importlib
import shlex
import sys

MODULES = ('xgboost_model_pb2', 'random_forest_model_pb2',
           'class_label_translation_pb2', 'linear_scaling_pb2')


def runtime_needs_update(exc):
    # Avoid importing VersionError: old runtimes may not provide it at all.
    if (type(exc).__module__ == 'google.protobuf.runtime_version'
            and type(exc).__name__ == 'VersionError'):
        return True
    if isinstance(exc, ModuleNotFoundError):
        name = exc.name or ''
        return name == 'google' or name == 'google.protobuf' or name.startswith('google.protobuf.')
    return (isinstance(exc, ImportError) and 'google.protobuf' in str(exc)
            and 'runtime_version' in str(exc))


def main():
    try:
        for name in MODULES:
            importlib.import_module('lillymol_protos.' + name)
    # VersionError is not available in older protobuf runtimes. Catch import
    # failures here, including that version-specific exception.
    except Exception as exc:
        print(f'Cannot load LillyMol script protos with {sys.executable}: {exc}', file=sys.stderr)
        if runtime_needs_update(exc):
            command = shlex.join([sys.executable, '-m', 'pip', 'install', '--upgrade', 'protobuf'])
            uv_command = shlex.join(['uv', 'pip', 'install', '--python', sys.executable,
                                     '--upgrade', 'protobuf'])
            print('Update protobuf in this Python environment, then rerun:', file=sys.stderr)
            print(f'  {command}', file=sys.stderr)
            print('Or, using uv:', file=sys.stderr)
            print(f'  {uv_command}', file=sys.stderr)
        else:
            print('Build //lillymol_protos:generated_protos and run src/copy_script_protos.sh.',
                  file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
