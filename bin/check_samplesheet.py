#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()

class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_ALIGN_FORMATS = (
        ".bam",
        ".cram"
    )

    VALID_INDEX_FORMATS = (
            ".bai",
            ".crai"
        )

    VALID_PED_FORMAT = (
            ".ped",
        )

    def __init__(
        self,
        sample_col="sample_id",
        bam_col="bam_or_cram",
        bai_col="bai_or_crai",
        batch_col="sample_set_id",
        cohort_col="sample_set_set_id",
        ped_col="cohort_ped",
        **kwargs,
    ):
        """
            Initialize the row checker with the expected column names.

            Args:
                sample_col (str): The name of the column that contains the sample name
                    (default "sample_id").
                bam_col (str): The name of the column that contains the BAM or CRAM file
                    path (default "bam_or_cram").
                bai_col (str): The name of the column that contains the BAI or CRAI file
                    path (default "bai_or_crai").
                batch_col (str): The name of the new column that contains the batch 
                    identifier (default "sample_set_id").
                cohort_col (str): The name of the new column that contains the cohort 
                    identifier (default "sample_set_set_id").
                ped_col (str): The name of the new column that contains the PED file 
                    path with per-sample family structure information. For more details,
                    please check https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format 
                    (default "cohort_ped").
        """
        super().__init__()
        self._sample_col = sample_col
        self._bam_col = bam_col
        self._bai_col = bai_col
        self._batch_col = batch_col
        self._cohort_col = cohort_col
        self._ped_col = ped_col
        self._seen = set()
        self.modified = []
        
    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_bam(row)
        self._validate_bai(row)
        self._validate_alignment_index_pair(row)
        self._validate_batch(row)
        self._validate_cohort(row)
        self._validate_ped(row)
        self._seen.add((row[self._sample_col], row[self._bam_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        if row[self._sample_col].find(" ") != -1:
            logger.warning(f"Spaces have been replaced by underscores for sample: {row[self._sample_col]}")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_bam(self, row):
        """Assert that the BAM or CRAM entry is non-empty and has the right format."""
        if len(row[self._bam_col]) <= 0:
            raise AssertionError("BAM OR CRAM file is required.")
        self._validate_bam_format(row[self._bam_col])
    
    def _validate_bai(self, row):
        """Assert that the BAI or CRAI entry is non-empty and has the right format."""
        if len(row[self._bai_col]) <= 0:
            raise AssertionError("BAI or CRAI file is required.")
        self._validate_bai_format(row[self._bai_col])

    def _validate_alignment_index_pair(self, row):
        """Assert that alignment and index extensions are a valid pair."""
        alignment = row[self._bam_col].lower()
        index = row[self._bai_col].lower()
        if alignment.endswith(".bam") and not index.endswith(".bai"):
            raise AssertionError(
                f"Invalid pair: BAM alignment must use BAI index. "
                f"Got alignment={row[self._bam_col]} index={row[self._bai_col]}."
            )
        if alignment.endswith(".cram") and not index.endswith(".crai"):
            raise AssertionError(
                f"Invalid pair: CRAM alignment must use CRAI index. "
                f"Got alignment={row[self._bam_col]} index={row[self._bai_col]}."
            )

    def _validate_bam_format(self, filename):
        """Assert that a given filename has one of the expected alignment extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_ALIGN_FORMATS):
            raise AssertionError(
                f"The alignment file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_ALIGN_FORMATS)}"
            )
    
    def _validate_bai_format(self, filename):
        """Assert that a given filename has one of the expected index extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_INDEX_FORMATS):
            raise AssertionError(
                f"The index file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_INDEX_FORMATS)}"
            )
    
    def _validate_batch(self, row):
        """Assert that the batch identifier exists and convert spaces to underscores."""
        if len(row[self._batch_col]) <= 0:
            raise AssertionError("Batch ID input is required.")
        if row[self._batch_col].find(" ") != -1:
            logger.warning(f"Spaces have been replaced by underscores for batch: {row[self._batch_col]}")
            # Sanitize batch identifier(s) slightly.
            row[self._batch_col] = row[self._batch_col].replace(" ", "_")
    
    def _validate_cohort(self, row):
        """Assert that the cohort identifier exists and convert spaces to underscores."""
        if len(row[self._cohort_col]) <= 0:
            raise AssertionError("Cohort ID input is required.")
        if row[self._cohort_col].find(" ") != -1:
            logger.warning(f"Spaces have been replaced by underscores for cohort: {row[self._cohort_col]}")
            # Sanitize cohort identifier(s) slightly.
            row[self._cohort_col] = row[self._cohort_col].replace(" ", "_")
    
    def _validate_ped(self, row):
        """Assert that the PED entry is non-empty and has the right format."""
        if len(row[self._ped_col]) <= 0:
            raise AssertionError("PED file is required.")
        self._validate_ped_format(row[self._ped_col])
    
    def _validate_ped_format(self, filename):
        """Assert that a given filename has the expected extension."""
        if not any(filename.endswith(extension) for extension in self.VALID_PED_FORMAT):
            raise AssertionError(
                f"The cohort ped file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_PED_FORMAT)}"
            )

    def validate_unique_samples(self):
        """Assert that the combination of sample name and BAM/CRAM filename is unique."""
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and BAM/CRAM must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1

def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect

def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by GermGatkSVonHPC pipeline.

    Validate the general shape of the table, expected columns, and each row.

    Args:
    file_in (pathlib.Path): The given tabular samplesheet. The format can be either
        CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
    file_out (pathlib.Path): Where the validated and transformed samplesheet should
        be created; always in CSV format.
    
    Example:
    This function checks that the samplesheet follows the following structure,

        sample_id,bam_or_cram,bai_or_crai,sample_set_id,sample_set_set_id,cohort_ped
        SAMPLE_ID1,SAMPLE_ID1.bam,SAMPLE_ID1.bai,BATCH1,COHORT1,SAMPLE_ID1.ped
        SAMPLE_ID2,SAMPLE_ID2.bam,SAMPLE_ID2.bai,BATCH1,COHORT1,SAMPLE_ID2.ped

        For an example see:
        https://github.com/balayev1/GermGatkSVonHPC/tree/main/data/samplesheet.csv
    """
    required_columns = {"sample_id", "bam_or_cram", "bai_or_crai", "sample_set_id", "sample_set_set_id", "cohort_ped"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
