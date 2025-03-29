import abc
import re
from pathlib import Path
from typing import TypeVar, Generic

from pydantic import FilePath, BaseModel


__all__ = ["RConfig", "RunRMixin", "RResultType"]

from proteomics.utils.run_subprocess import run_command


class RConfig(BaseModel):
    rscript_bin: FilePath = "/opt/conda/bin/Rscript"


RResultType = TypeVar("RResultType")


class RunRMixin(RConfig, abc.ABC, Generic[RResultType]):
    @abc.abstractmethod
    def get_r_script(self) -> Path:
        pass

    @abc.abstractmethod
    def process_stdout(self, stdout: list[str]) -> RResultType:
        pass

    def ignored_linters(self) -> set[str]:
        return {"object_usage_linter"}

    def raise_if_warning(self, stdout: list[str]):
        regex_str = r"R:\d+:\d+: (?P<lint_type>\w+): \[(?P<lint_name>\w+)\]"

        for line in stdout:
            matches = re.findall(regex_str, line, re.DOTALL)

            for match in matches:
                if (match[0] == "warning") and match[
                    1
                ] not in self.ignored_linters():
                    raise ValueError(line)
                else:
                    print(f"Ignoring linter: {match}")

    def lint_r_script(self):
        """
        Runs R -e "lint(filename = 'self.get_r_script()')" on the R script
        """
        command = [
            self.rscript_bin,
            "-e",
            f"lintr::lint(filename = '{self.get_r_script()}')",
        ]
        result = run_command(command)
        self.raise_if_warning(result)

    def create_command(self) -> list[str]:
        r_script = self.get_r_script()
        # run r_script
        print(f"Running R script: {r_script}")

        command = [self.rscript_bin, r_script]

        for arg_flag, argument in self.model_dump(
            exclude={"rscript_bin", "hash_str_length"}
        ).items():
            command.extend([f"--{arg_flag}", argument])

        command = [str(arg) for arg in command]

        return command

    def run_analysis(self) -> RResultType | None:
        self.lint_r_script()

        command = self.create_command()
        print(" ".join(command))

        result = run_command(command, cwd=self.get_r_script().parent)

        if not result:
            print(f"Error running command: {result}")
            return None

        try:
            return self.process_stdout(result)
        except Exception as e:
            print(f"Error processing stdout: {e}")
            print(result)
