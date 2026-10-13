import os
import sys
import csv
import json
import tempfile
from pathlib import Path
import pytest

# Ensure bin/ is on Python path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

from generate_methylation_dashboard import parse_methbat_tsv, calculate_stats


class TestSamplesheetIngestion:
    """Validate samplesheet format, column requirements, and URI schemes."""

    REQUIRED_HEADERS = ["run_id", "input_path", "input_type", "kit", "model", "modified_bases_model"]
    VALID_INPUT_TYPES = ["pod5", "fast5", "bam", "fastq", "pacbio_bam", "pacbio_hifi", "hifi", "ubam", "aligned_bam"]

    def test_pacbio_samplesheet_structure(self):
        csv_path = REPO_ROOT / "tests" / "csv" / "samplesheet_pacbio_test.csv"
        assert csv_path.exists(), "PacBio test samplesheet must exist"

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == self.REQUIRED_HEADERS
            rows = list(reader)
            assert len(rows) >= 1
            for row in rows:
                assert row["input_type"] in self.VALID_INPUT_TYPES
                assert row["input_path"].startswith("http://") or row["input_path"].startswith("https://") or Path(row["input_path"]).name.endswith(".bam")

    def test_nanopore_samplesheet_structure(self):
        csv_path = REPO_ROOT / "tests" / "samplesheet_test.csv"
        assert csv_path.exists(), "Nanopore test samplesheet must exist"

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == self.REQUIRED_HEADERS
            rows = list(reader)
            assert len(rows) >= 1
            for row in rows:
                assert row["input_type"] in self.VALID_INPUT_TYPES


class TestMethylationDashboardGenerator:
    """Validate parser, calculations, and HTML payload injection for the clinical dashboard."""

    def test_parse_methbat_tsv_with_mock_data(self):
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
            tf.write("chrom\tstart\tend\tcpg_label\texpected_category\tsummary_label\tmean_combined_methyl\tmean_hap1_methyl\tmean_hap2_methyl\tmean_meth_delta\tqc_warnings\n")
            tf.write("chr1\t1000\t2000\tH19-DMR\tImprinting\tPASS\t0.520\t0.980\t0.040\t0.940\tnone\n")
            tf.write("chr7\t5000\t6000\tPEG10\tImprinting\tPASS\t0.490\t0.950\t0.020\t0.930\tlow_coverage\n")
            tf.flush()
            temp_tsv = Path(tf.name)

        try:
            rows = parse_methbat_tsv(temp_tsv)
            assert len(rows) == 2
            assert rows[0]["cpg_label"] == "H19-DMR"
            assert rows[0]["qc_status"] == 0  # PASS
            assert rows[1]["qc_status"] == 1  # Warning
            assert float(rows[0]["mean_combined_methyl"]) == 0.520
            assert float(rows[0]["mean_hap1_methyl"]) == 0.980

            stats = calculate_stats(rows, "test_sample")
            assert stats["passCount"] == 1
            assert stats["cov"] == "Phased-5mC"
            assert "96.5%" in stats["h1Avg"]  # Average of 98.0 and 95.0
        finally:
            if temp_tsv.exists():
                os.remove(temp_tsv)

    def test_parse_empty_tsv_fallback(self):
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
            tf.write("chrom\tstart\tend\n")
            tf.flush()
            temp_tsv = Path(tf.name)

        try:
            rows = parse_methbat_tsv(temp_tsv)
            assert rows == []
            stats = calculate_stats(rows, "empty_sample")
            assert stats["passCount"] == 0
            assert stats["cov"] == "N/A"
        finally:
            if temp_tsv.exists():
                os.remove(temp_tsv)


class TestConfigAndModuleIntegrity:
    """Validate nextflow configs and subworkflow bindings."""

    def test_nextflow_profiles_exist(self):
        nf_config = REPO_ROOT / "nextflow.config"
        content = nf_config.read_text()
        for expected_profile in ["test_ont", "test_pacbio", "human", "eukaryote", "prokaryote", "virus", "fraka", "version_lock"]:
            assert expected_profile in content, f"Profile '{expected_profile}' missing in nextflow.config"

    def test_assets_dashboard_template_present(self):
        template = REPO_ROOT / "assets" / "methylation_dashboard" / "dashboard-template.html"
        assert template.exists(), "dashboard-template.html must be present in assets"
        content = template.read_text(errors="ignore")
        assert "{{SAMPLE_DATA_JSON}}" in content
        assert "{{SAMPLE_STATS_JSON}}" in content
