## vcf_comparator.py
import shutil
import tempfile
from pathlib import Path
from typing import Iterator, List, NamedTuple

import click
import pysam
from pysam import bcftools
from pysam.utils import SamtoolsError

from vcf_isec.errors import FileFormatError, IntersectError
from vcf_isec.parser import parse_variants


class ISecOutput(NamedTuple):
    """
    Represents the output of the intersection operation between two sets of variant records.
    """

    intersection: List[pysam.VariantRecord]
    setdiff1: List[pysam.VariantRecord]
    setdiff2: List[pysam.VariantRecord]


class VCFIntersection:
    """
    Class for comparing two VCF files and extracting shared and unique variants.

    The provided vcf files need to be compressed and indexed (i.e. csi or tbi files exist in the same directory).

    Args:
        vcf1 (Path | str): Path to the first VCF file.
        vcf2 (Path | str): Path to the second VCF file.
        prefix (Path | str | None, optional): Path to the output directory. Defaults to None.

    Raises:
        FileFormatError: If either isec_dir or out_dir are files (expect directories or nonexistant).
    """

    def __init__(
        self,
        vcf1: Path | str,
        vcf2: Path | str,
        out_dir: Path | str,
        isec_dir: Path | str | None = None,
    ) -> None:
        self.vcf1 = Path(vcf1)
        self.vcf2 = Path(vcf2)

        self.out_dir = Path(out_dir)
        if self.out_dir.is_file():
            raise FileFormatError(f"{self.out_dir} is not a directory.")

        self.isec_dir = Path(isec_dir) if isec_dir else None
        if self.isec_dir and self.isec_dir.is_file():
            raise FileFormatError(f"{self.isec_dir} is not a directory.")

    def intersect(self) -> ISecOutput:
        """
        Compares two VCF files and returns a tuple containing lists of shared variants,
        unique variants in file1, and unique variants in file2.

        Returns:
            Tuple[List[str], List[str], List[str]]: shared variants, unique to file1, unique to file2

        Raises:
            IntersectError: If there is an error in the pysam library bcftools.isec function.
            FileFormatError: If there is an error reading the file.
        """

        out_exists = self.out_dir.is_dir()
        isec_exists = self.isec_dir and self.isec_dir.is_dir()

        if out_exists and isec_exists:
            click.confirm(
                f"Directories {self.out_dir} and {self.isec_dir} will be overwritten. Continue?",
                default=True,
                abort=True,
            )
        elif out_exists:
            click.confirm(
                f"Directory {self.out_dir} will be overwritten. Continue?",
                default=True,
                abort=True,
            )
        elif isec_exists:
            click.confirm(
                f"Directory {self.isec_dir} will be overwritten. Continue?",
                default=True,
                abort=True,
            )

        if isec_exists:
            shutil.rmtree(self.isec_dir)  # type: ignore (can't be None)

        if self.isec_dir:
            self.isec_dir.mkdir(parents=True, exist_ok=True)
            isec_files = self._isec(self.isec_dir)
        else:
            with tempfile.TemporaryDirectory() as tmpdirname:
                isec_files = self._isec(Path(tmpdirname))

        out_files = self._merge_output(isec_files)

        results = self._gather_variant_results(out_files)

        return results

    def _isec(self, isec_dir: Path) -> Iterator[Path]:
        """
        Runs the isec operation and produces the intermediary files.
        """

        try:
            bcftools.isec(
                str(self.vcf1),
                str(self.vcf2),
                f"--prefix={str(isec_dir)}",
                # "--nfiles=2",
                "-w 1,2",
                capture_stdout=False,
            )
        except SamtoolsError as e:
            raise IntersectError(
                f"Error in pysam library bcftools.isec function: {e}"
            ) from e
        except OSError as e:
            raise FileFormatError(f"Error reading file: {e.filename}") from e

        for i in range(4):
            yield isec_dir / f"000{i}.vcf"

    def _merge_output(self, isec_files: Iterator[Path]) -> Iterator[Path]:
        """
        Copies the output of the intersection operation to the specified directory.

        Returns:
            Iterator[Path]: The paths to the output files.
        """

        if self.out_dir.exists():
            shutil.rmtree(self.out_dir)

        self.out_dir.mkdir(parents=True, exist_ok=True)

        # copy the two setdiff files to the output directory
        setdiff1 = self.out_dir / self.vcf1.stem.replace(".vcf", "_setdiff.vcf")
        setdiff2 = self.out_dir / self.vcf2.stem.replace(".vcf", "_setdiff.vcf")

        shutil.copy(next(isec_files), setdiff1)
        shutil.copy(next(isec_files), setdiff2)

        intersection = self.out_dir / "intersection.vcf"
        intermediate1 = intersection.with_stem(intersection.stem + "-1").with_suffix(
            ".vcf.gz"
        )
        intermediate2 = intersection.with_stem(intersection.stem + "-2").with_suffix(
            ".vcf.gz"
        )
        bcftools.view(
            str(next(isec_files)),
            "-Oz",
            "-o",
            str(intermediate1),
            "--write-index",
            catch_stdout=False,
        )

        bcftools.view(
            str(next(isec_files)),
            "-Oz",
            "-o",
            str(intermediate2),
            "--write-index",
            catch_stdout=False,
        )

        # merge the two intersection files into one
        bcftools.merge(
            str(intermediate1),
            str(intermediate2),
            "-Ov",
            "-o",
            str(intersection),
            "--force-samples",
            catch_stdout=False,
        )

        for f in self.out_dir.glob(f"{intersection.stem}-*"):
            f.unlink()

        yield from [setdiff1, setdiff2, intersection]

    def _gather_variant_results(self, out_files) -> ISecOutput:
        """
        Gathers the results of the intersection operation.

        Returns:
            ISecOutput: A named tuple containing the shared variants, unique to file1, and unique to file2.
        """

        unique1 = list(parse_variants(next(out_files)))
        unique2 = list(parse_variants(next(out_files)))
        shared = list(parse_variants(next(out_files)))

        return ISecOutput(shared, unique1, unique2)
