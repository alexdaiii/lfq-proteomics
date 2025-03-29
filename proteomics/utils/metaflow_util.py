from pathlib import Path

from metaflow import FlowSpec, current


def get_run_output(
    flowspec: FlowSpec,
) -> Path:
    return (
        Path(flowspec._datastore._storage_impl.datastore_root)
        / current.pathspec
    ).parent.parent


def get_task_output(
    flowspec: FlowSpec,
) -> Path:
    return (
        Path(flowspec._datastore._storage_impl.datastore_root)
        / current.pathspec
    )
