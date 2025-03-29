import json
from datetime import datetime
from pathlib import Path
from typing import TypeVar
from pydantic import BaseModel, FilePath, DirectoryPath, NewPath
import hashlib


__all__ = [
    "BaseParams",
    "SampleInfo",
    "fix_dir_params",
    "is_xlsx_file",
    "get_project_root",
]


class SampleInfo(BaseModel):
    input_file: Path


# Define a generic type variable
T = TypeVar("T", bound="SampleInfo")


def get_project_root() -> Path:
    return Path(__file__).parent.parent.parent.resolve()


class BaseParams(BaseModel):
    hash_str_length: int = 12
    output_dir: NewPath | DirectoryPath | None = None

    def get_run_output_dir(self) -> Path:
        if self.output_dir is None:
            raise ValueError("output_dir is not set")

        return self.output_dir / self.hash_params()

    def hash_params(self) -> str:
        # Get the model dump as a JSON string
        model_dump_str = self.model_dump_json(exclude={"get_samples"})

        # Hash the JSON string
        hash_str = hashlib.sha256(model_dump_str.encode()).hexdigest()

        return hash_str[: self.hash_str_length]

    def save_config(self) -> None:
        """
        Save the config to a file in the output directory.
        """
        # Get the JSON string with indentation
        json_str = json.loads(self.model_dump_json())
        json_str["run_name"] = self.hash_params()
        json_str["run_date"] = datetime.now().isoformat()
        json_str = json.dumps(json_str, indent=2)

        # Save the JSON string to a file
        output_file = self.get_run_output_dir() / "config.json"

        if not self.get_run_output_dir().exists():
            self.get_run_output_dir().mkdir(parents=True)
            print(f"Created output directory {self.get_run_output_dir()}")
        else:
            print(
                f"Output directory {self.get_run_output_dir()} already exists."
            )

        with open(output_file, "w") as f:
            f.write(json_str)



def fix_dir_params(
    config: dict, *, group_name: str | None, config_keys: list[str]
) -> dict:
    """Converts a list of config keys to Path objects that are relative
    to the project root.

    Args:
        config (dict): The config loaded from a yaml or json file.
        group_name (str): The key in the config dict that contains the
            config keys to be converted to Path objects.
        config_keys (list[str]): The keys in the config dict that should
            be converted to Path objects.


    Returns:
        dict: The config dict with the specified keys converted to Path
            objects.

    """
    for key in config_keys:
        config[group_name][key] = get_project_root() / config.get(
            group_name, {}
        ).get(key, "")

    return config


def is_xlsx_file(value: FilePath) -> FilePath:
    if not value.__str__().endswith(".xlsx"):
        raise ValueError("File must be an excel file")
    return value
