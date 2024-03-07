## test_vcf_comparator.py

import random
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest
from vcf_isec.isec import VCFIntersection
from vcf_isec.parser import VCFPreparer, parse_header, parse_variants

RESOURCES = Path(__file__).parent.parent / "resources"

SAMPLE_VCF = RESOURCES / "sample.vcf"
SIMPLE_CMP_VCFS = [RESOURCES / "cmp-test-a.vcf", RESOURCES / "cmp-test-b.vcf"]


@dataclass
class VCFFixture:
    file1: Path
    file2: Path
    intersected_count: int
    complement1_count: int
    complement2_count: int


class TestVCFIntersection:
    @pytest.fixture(scope="session")
    def vcf_fixture(self):
        """
        Fixture to create temporary VCF files for testing.
        """

        header = parse_header(SAMPLE_VCF)
        variants = parse_variants(SAMPLE_VCF)

        intersected_count = 0
        complement1_count = 0
        complement2_count = 0

        with tempfile.TemporaryDirectory() as prefix:
            prefix = Path(prefix)
            fpath1 = prefix / "file1.vcf"
            fpath2 = prefix / "file2.vcf"
            with open(fpath1, "w") as f1, open(fpath2, "w") as f2:
                f1.write(header)
                f2.write(header)

                for variant in variants:
                    # trunk-ignore(bandit/B311)
                    in_f1, in_f2 = random.choices([True, False], k=2)

                    if in_f1 and in_f2:
                        intersected_count += 1
                    elif in_f1:
                        complement1_count += 1
                    elif in_f2:
                        complement2_count += 1

                    if in_f1:
                        f1.write(str(variant))
                    if in_f2:
                        f2.write(str(variant))

            yield VCFFixture(
                fpath1,
                fpath2,
                intersected_count,
                complement1_count,
                complement2_count,
            )

    def test_isec_simple(self, monkeypatch: pytest.MonkeyPatch):
        """
        Test the VCFIntersection class with simple data.
        """

        monkeypatch.setattr("click.confirm", lambda *args, **kwargs: False)

        file1, file2 = SIMPLE_CMP_VCFS
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf1 = VCFPreparer(file1, tmp_dir).prepare()
            vcf2 = VCFPreparer(file2, tmp_dir).prepare()

            prefix = Path(tmp_dir) / "isec"
            isec = VCFIntersection(vcf1, vcf2, prefix)
            intersected, unique1, unique2 = isec.compare_vcf()
            assert len(intersected) == 5, "Expected 5 intersected variants."
            assert len(unique1) == 1, "Expected no unique variants in file1."
            assert len(unique2) == 1, "Expected no unique variants in file2."

    def test_randomized_isec(
        self, monkeypatch: pytest.MonkeyPatch, vcf_fixture: VCFFixture
    ):
        """
        Generalized test for comparing VCF files.
        """

        monkeypatch.setattr("click.confirm", lambda *args, **kwargs: False)

        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf1 = VCFPreparer(vcf_fixture.file1, tmp_dir).prepare()
            vcf2 = VCFPreparer(vcf_fixture.file2, tmp_dir).prepare()

            prefix = Path(tmp_dir) / "isec"
            isec = VCFIntersection(vcf1, vcf2, prefix)
            intersected, unique1, unique2 = isec.compare_vcf()
            assert (
                len(intersected) == vcf_fixture.intersected_count
            ), f"Expected {vcf_fixture.intersected_count} intersected variants, found {len(intersected)}."
            assert (
                len(unique1) == vcf_fixture.complement1_count
            ), f"Expected {vcf_fixture.complement1_count} unique variants in file1, found {len(unique1)}."
            assert (
                len(unique2) == vcf_fixture.complement2_count
            ), f"Expected {vcf_fixture.complement2_count} unique variants in file2, found {len(unique2)}."
