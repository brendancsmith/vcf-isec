import tempfile
from pathlib import Path

import click

from vcf_isec import parser
from vcf_isec.errors import ErrorHandler
from vcf_isec.isec import VCFIntersection


@click.command()
@click.option(
    "-p",
    "--prefix",
    type=click.Path(
        exists=False,
        file_okay=False,
        dir_okay=True,
        allow_dash=False,  # TODO: allow writing to stdout
        writable=True,
        path_type=Path,
    ),
    default=None,
    help="The directory to write isec outputs. Will create a tempory directory if not set.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(
        exists=False,
        file_okay=False,
        dir_okay=True,
        allow_dash=False,  # TODO: allow writing to stdout
        writable=True,
        path_type=Path,
    ),
    default="output",
    help="The directory to write the output files.",
)
@click.argument(
    "file1",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        allow_dash=False,
        readable=True,
        path_type=Path,
    ),
)
@click.argument(
    "file2",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        allow_dash=False,
        readable=True,
        path_type=Path,
    ),
)
def main(prefix, output, file1, file2):
    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf1 = parser.VCFPreparer(file1, tmp_dir).prepare()
            vcf2 = parser.VCFPreparer(file2, tmp_dir).prepare()

            isec = VCFIntersection(vcf1, vcf2, output, prefix or Path(tmp_dir) / "isec")
            shared, unique1, unique2 = isec.intersect()

        print(f"Unique to {file1.name}: {len(unique1)}")
        print(f"Unique to {file1.name}: {len(unique2)}")
        print(f"Shared Variants: {len(shared)}")
    except KeyboardInterrupt as e:
        click.echo("Received termination signal. Exiting now!")
        raise click.Abort(1) from e
    except Exception as e:
        ErrorHandler().handle(e)
        raise click.Abort(1) from e


if __name__ == "__main__":
    main()
