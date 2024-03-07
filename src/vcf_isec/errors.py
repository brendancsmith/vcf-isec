import logging
import sys
from datetime import datetime
from pathlib import Path

import click

log = logging.getLogger(__name__)


class FileFormatError(Exception):
    """
    An exception to be raised when a file is not in the expected format.
    """

    def __init__(self, filename: str):
        self.filename = filename
        super().__init__(f"File is unknown or unsupported format: {filename}")


class IntersectError(Exception):
    """
    An exception that wraps/decouples the isec depenency.
    """

    def __init__(self, *args):
        if not len(args):
            args = [
                f"An error occurred while running the intersection operation: {self.__cause__}."
            ]
            super().__init__(*args)


class ClickStreamHandler(logging.StreamHandler):
    def emit(self, record):
        try:
            msg = self.format(record)

            if record.levelno >= logging.ERROR:
                msg = click.style(msg, fg="red", bold=True)
            elif record.levelno >= logging.WARNING:
                msg = click.style(msg, fg="yellow")

            click.echo(msg)

        except Exception:
            self.handleError(record)


class ErrorHandler:
    """
    A class to handle errors during VCF file parsing and comparison.
    """

    LOG_DIR = Path("logs")

    def __init__(self):
        self.configure_logging()

    def configure_logging(self) -> None:
        """
        Configures logging to ensure that it is set up before handling any errors.
        """
        pkg_logger = logging.getLogger(__package__)
        if pkg_logger.hasHandlers():
            return

        self.LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_file = self.LOG_DIR / datetime.now().strftime("%Y-%m-%d %H-%M-%S.log")
        log_fmt = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        # Log all messages to a file
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(log_fmt)
        file_handler.setLevel(logging.DEBUG)
        pkg_logger.addHandler(file_handler)

        # Log only errors and warnings to stderr
        stderr_handler = ClickStreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        pkg_logger.addHandler(stderr_handler)

    @staticmethod
    def handle(exc: Exception) -> bool:
        """
        Handles exceptions by logging an error message to a file.
        Ensures that logging is configured before attempting to log the error.

        Parameters:
            error (Exception): The exception that occurred.

        Returns:
            bool: Unrecoverable if True, otherwise False.
        """

        log.debug("Caught exception:", exc_info=True)
        if isinstance(exc, IntersectError):
            log.error(exc)
            return True
        elif isinstance(exc, FileNotFoundError):
            log.error(f"File not found: {exc.filename}")
        elif isinstance(exc, FileFormatError):
            log.error(f"File is unknown/unsupported format: {exc.filename}")
        elif isinstance(exc, OSError):
            log.error(f"Error reading file: {exc.filename}")
        else:
            log.critical("Unhandled exception: %s" % exc)
            return True

        return False
        return False
