## test_vcf_comparator.py

import random
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import IO

import pytest
from vcf_isec.parser import parse_variants

SAMPLE_VCF = Path(__file__).parent.parent / "resources/sample.vcf"


@dataclass
class VCFFixture:
    file1: IO[str]
    file2: IO[str]
    shared_count: int
    complement1_count: int
    complement2_count: int


class TestVCFIntersection:
    @pytest.fixture(scope="session")
    def vcf_files(self):
        """
        Fixture to create temporary VCF files for testing.
        """

        variants = parse_variants(SAMPLE_VCF)

        def tmp():
            return tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)

        with tmp() as f1, tmp() as f2:
            for variant in variants:
                # trunk-ignore(bandit/B311)
                in_f1, in_f2 = random.choices([True, False], k=2)

                if in_f1:
                    f1.write(str(variant))
                if in_f2:
                    f2.write(str(variant))

                yield f1.file, f2.file

    def test_compare_vcf(self, vcf_files):
        """
        Generalized test for comparing VCF files.
        """

        # comparator, should_pass, shared_count, unique_count = setup_vcf_comparator
        # if not should_pass:
        #     with pytest.raises(VCFError):  # Corrected exception handling
        #         shared, unique1, unique2 = comparator.compare_vcf()
        # else:
        #     shared, unique1, unique2 = comparator.compare_vcf()
        #     if shared_count is not None:
        #         assert len(shared) == shared_count, f"Expected {shared_count} shared variants."
        #     if unique_count is not None:
        #         # Separately checking the counts of unique variants in each file
        #         assert len(unique1) == unique_count, f"Expected {unique_count} unique variants in file1."
        #         assert len(unique2) == unique_count, f"Expected {unique_count} unique variants in file2."
