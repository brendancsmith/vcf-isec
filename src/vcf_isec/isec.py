## vcf_comparator.py
import tempfile
from pathlib import Path
from typing import List, NamedTuple

from pysam import bcftools
from pysam.utils import SamtoolsError

from vcf_isec.errors import FileFormatError, IntersectError
from vcf_isec.parser import Variant, parse_variants


class ISecOutput(NamedTuple):
    intersection: List[Variant]
    complement1: List[Variant]
    complement2: List[Variant]


class VCFIntersection:
    def __init__(
        self,
        file1_path: Path | str,
        file2_path: Path | str,
        prefix: Path | str | None = None,
    ) -> None:
        self.file1_path = file1_path
        self.file2_path = file2_path

        if prefix:
            self.prefix = Path(prefix)
            if self.prefix.is_file():
                raise FileFormatError(f"{self.prefix} is not a directory.")
        else:
            self.prefix = None

    def compare_vcf(self) -> ISecOutput:
        """
        Compares two VCF files and returns a tuple containing lists of shared variants,
        unique variants in file1, and unique variants in file2.

        Returns:
            Tuple[List[str], List[str], List[str]]: shared variants, unique to file1, unique to file2
        """

        try:
            if self.prefix:
                self.prefix.mkdir(parents=True, exist_ok=True)
                results = self._compare_vcf(self.prefix)
            else:
                with tempfile.TemporaryDirectory() as tmpdirname:
                    results = self._compare_vcf(tmpdirname)

        except SamtoolsError as e:
            raise IntersectError(
                f"Error in pysam library bcftools.isec function: {e}"
            ) from e
        except OSError as e:
            raise FileFormatError(f"Error reading file: {e.filename}") from e

        return results

    def _compare_vcf(self, dirname: Path | str) -> ISecOutput:
        dirname = Path(dirname)

        bcftools.isec(
            str(self.file1_path),
            str(self.file2_path),
            f"--prefix={str(dirname)}",
            # "--nfiles=2",
            "-w 1,2",
            capture_stdout=False,
        )

        complements1_path = dirname / "0000.vcf"
        complements2_path = dirname / "0001.vcf"
        intersection1_path = dirname / "0002.vcf"
        intersection2_path = dirname / "0003.vcf"

        # TODO: merge or use sitemap
        intersection1 = list(parse_variants(intersection1_path))
        intersection2 = list(parse_variants(intersection2_path))
        # trunk-ignore(bandit/B101)
        assert len(intersection1) == len(intersection2)

        complement1 = list(parse_variants(complements1_path))
        complement2 = list(parse_variants(complements2_path))

        return ISecOutput(intersection1, complement1, complement2)
