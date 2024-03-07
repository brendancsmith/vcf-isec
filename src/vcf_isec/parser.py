import types
from pathlib import Path
from typing import Iterable, cast

import click
import pysam
from pysam import VariantRecord, bcftools


class VCFPreparer:
    """Builder pattern for preparing VCF files for intersection operation."""

    def __init__(self, path: Path | str, prefix: Path | str = "output"):
        self.provided_path = Path(path)
        self.prefix = Path(prefix)

        self.gz_path = None
        self.gz_is_new = False

        self.tbi_path = None

    def prepare(self) -> Path:
        if self.provided_path.suffix == ".vcf":
            self.compress()
        elif self.provided_path.suffix == ".gz":
            self.gz_path = self.provided_path

        self.index()

        if self.gz_path is None:
            raise RuntimeError("BGZipped VCF file not set.")

        return self.gz_path

    def compress(self):
        """
        Compresses the provided VCF file to a BGZipped VCF file.

        -
        -
        - If the user chooses to overwrite the existing file or if no
            existing file is found, a new BGZipped VCF file is created.

        Returns:
            None
        """

        # If a BGZipped VCF file already exists in the same directory...
        gz_path = self.provided_path.with_suffix(".vcf.gz")
        if gz_path.exists():
            # ... the user is prompted to use it or overwrite it.
            if click.confirm(f"BGZipped VCF found at {gz_path}. Use it?"):
                self.gz_path = gz_path
                return

            if click.confirm(f"Overwrite {gz_path}?"):
                self.gz_path = gz_path

        # Otherwise, the user is prompted to create a new BGZipped VCF file.
        elif click.confirm(f"Create BGZipped VCF at {gz_path}?"):
            self.gz_path = gz_path

        # If the user declines, we use the prefix directory
        if self.gz_path is None:
            self.gz_path = self.prefix / self.provided_path.with_suffix(".vcf.gz").name

        # Sort input VCF and output to BGZipped VCF
        bcftools.sort(
            "-Oz", "-o", str(self.gz_path), str(self.provided_path), catch_stdout=False
        )
        self.gz_is_new = True

    def index(self):
        # TODO: I think bcftools.isec expects index to be in the same directory,
        #      so we may need to copy the VCF to the prefix directory and write
        #      the index there too

        # Ensure that we have a compressed VCF file to index
        if self.gz_path is None:
            import inspect

            method_name = cast(types.FrameType, inspect.currentframe()).f_code.co_name
            raise RuntimeError(
                f"{method_name} should not have been called before {self.gz_path.__qualname__} is set."
            )

        # Construct the path for the tabix index file
        tbi_path = Path(f"{self.gz_path}.tbi")

        # If the prefix directory is the expected parent of the index, set tbi_path
        if self.prefix and self.prefix == tbi_path.parent:
            self.tbi_path = tbi_path

        elif not self.gz_is_new:
            # If the compressed VCF is not new, we can check to use existing index
            if tbi_path.exists():
                if click.confirm("Tabix index found at {tbi_path}. Use it?"):
                    self.tbi_path = tbi_path
                    return

            # If the tbi file does not exist, ask user if they want to create it
            elif click.confirm(f"Create tabix index at {tbi_path}?"):
                self.tbi_path = tbi_path

            # If user does not want to create tbi file, set the tbi_path to the prefix
            else:
                self.tbi_path = Path(self.prefix.name) / f"{self.gz_path.name}.tbi"

        # If compressed VCF is new, we try to overwrite old index if it exists
        elif tbi_path.exists():
            if click.confirm(f"Tabix index found at {tbi_path}. Overwrite it?"):
                self.tbi_path = tbi_path
                return

            self.tbi_path = Path(self.prefix.name) / f"{self.gz_path.name}.tbi"

        # If the tbi file does not exist, ask user if they want to create it
        elif click.confirm(f"Create tabix index at {tbi_path}?"):
            self.tbi_path = tbi_path
        else:
            # Otherwise, set the tbi_path to the prefix directory
            self.tbi_path = Path(self.prefix.name) / f"{self.gz_path.name}.tbi"

        # Create index
        pysam.tabix_index(
            str(self.gz_path),
            force=True,
            preset="vcf",
            index=str(self.tbi_path),
            keep_original=True,  # not strictly necessary since VCF is already compressed
        )


def parse_header(vcf_path: Path | str) -> str:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return str(vcf.header)


def parse_variants(vcf_path: Path | str) -> Iterable[VariantRecord]:
    """
    Extracts variants from a VCF file using pysam.VariantFile.

    Args:
        vcf_path (str): Path to the VCF file.

    Returns:
        Iterable[str]: Iterable of variants.
    """

    with pysam.VariantFile(str(vcf_path)) as vcf:
        yield from vcf.fetch()
