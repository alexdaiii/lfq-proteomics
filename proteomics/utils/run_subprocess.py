import selectors
import subprocess
from pathlib import Path
from subprocess import SubprocessError


def _run_subprocess(command: list[str], *, cwd: Path = None) -> list[str]:
    """
    Run a subprocess and capture the output

    Args:
        command:  list[str]: The command to run
        cwd:  Path: The current working directory to run the command in

    Raises:
        SubprocessError: If the subprocess returns a non-zero exit code

    Returns: stdout if the command was successful

    """
    stdout = []
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd
    ) as process:
        sel = selectors.DefaultSelector()
        sel.register(process.stdout, selectors.EVENT_READ)
        sel.register(process.stderr, selectors.EVENT_READ)

        while True:
            for key, _ in sel.select():
                data = key.fileobj.read1().decode()
                if not data:
                    return stdout
                if key.fileobj is process.stdout:
                    print(data, end="")
                    stdout.append(data)
                else:
                    print(data, end="")

    return stdout


def add_metaflow_items(command: list[str]):
    try:
        from metaflow import current
        from metaflow.cards import Markdown

        print("Metaflow successfully imported")
    except ModuleNotFoundError:
        print("Metaflow not found. Skipping.")
        return

    # check if current has card attribute
    if hasattr(current, "card"):
        current.card.append(
            Markdown(f"""
```bash
{" \\ \n".join(command)}
```
    """)
        )


def run_command(
    command: list[str], *, cwd: Path = None, raise_if_fail: bool = True
) -> list[str] | None:
    """
    Runss the command and captures the output

    Args:
        command:  list[str]: The command to run
        cwd:  Path: The current working directory to run the command in

    Returns: stdout if the command was successful, None otherwise

    """
    command = [str(c) for c in command]

    add_metaflow_items(command)

    try:
        stdout = _run_subprocess(command, cwd=cwd)
    except SubprocessError as e:
        # Handle errors in the subprocess
        print(f"Error running command: {command}")
        if raise_if_fail:
            raise e

        return None

    return stdout
