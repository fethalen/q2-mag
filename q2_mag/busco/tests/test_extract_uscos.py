# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from pytest import approx
from q2_types.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_mag.busco.extract_orthologs import (
    CaseMode,
    DuplicateMode,
    FragmentMode,
    SequenceType,
    _append_uscos,
    _extract_uscos,
    _filter_by_mode,
    _filter_full_table_df,
    _get_corresponding_files,
    _normalize_fasta_header,
    _parse_full_table_file,
)


class TestBUSCOExtractOrthologs(TestPluginBase):
    package = "q2_mag.busco.tests"

    def setUp(self):
        super().setUp()
        self.full_table_fp_1 = self.get_data_path("full_tables/full_table_1.tsv")
        self.full_table_fp_2 = self.get_data_path("full_tables/full_table_2.tsv")

        self.full_table_df_1 = _parse_full_table_file(self.full_table_fp_1)
        self.full_table_df_2 = _parse_full_table_file(self.full_table_fp_2)

        self.busco_dir = Path(self.get_data_path("busco_results"))

    def test_is_nucleotide_returns_true_for_nucleotide(self):
        assert SequenceType.NUCLEOTIDE.is_nucleotide() is True
        assert SequenceType.NUCLEOTIDE.is_protein() is False

    def test_is_protein_returns_true_for_protein(self):
        assert SequenceType.PROTEIN.is_protein() is True
        assert SequenceType.PROTEIN.is_nucleotide() is False

    def test_enum_values_are_correct(self):
        assert SequenceType.NUCLEOTIDE.value == "nucleotide"
        assert SequenceType.PROTEIN.value == "protein"

    def test_parse_full_table_file(self):
        df = _parse_full_table_file(self.full_table_fp_1)

        assert isinstance(df, pd.DataFrame)

        expected_cols = [
            "busco_id",
            "status",
            "sequence",
            "gene_start",
            "gene_end",
            "strand",
            "score",
            "length",
        ]
        assert list(df.columns) == expected_cols

        assert len(df) == 8
        assert df.loc[0, "busco_id"] == "44125at2759"
        assert df.loc[0, "score"] == approx(1653.2)
        assert df.loc[1, "strand"] == "-"
        assert df.loc[1, "length"] == 210
        assert df.loc[7, "status"] == "Missing"

    def test_filter_by_mode_skip_duplicates(self):
        df = _filter_by_mode(self.full_table_df_1, DuplicateMode.SKIP)
        assert df.empty

    def test_filter_by_mode_skip_fragments(self):
        df = _filter_by_mode(self.full_table_df_1, FragmentMode.SKIP)
        assert df.empty

    def test_filter_by_mode_best_scoring_duplicates(self):
        df = _filter_by_mode(self.full_table_df_1, DuplicateMode.BEST_SCORE)
        assert len(df) == 1
        assert df.loc[0, "score"] == approx(714.9)

    def test_filter_by_mode_best_scoring_fragments(self):
        df = _filter_by_mode(self.full_table_df_1, FragmentMode.BEST_SCORE)
        assert len(df) == 1
        assert df.loc[0, "score"] == approx(195.3)

    def test_filter_by_mode_longest_duplicates(self):
        df = _filter_by_mode(self.full_table_df_1, DuplicateMode.LONGEST)
        assert len(df) == 1
        assert df.loc[0, "length"] == 483

    def test_filter_by_mode_longest_fragments(self):
        df = _filter_by_mode(self.full_table_df_1, FragmentMode.LONGEST)
        assert len(df) == 1
        assert df.loc[0, "length"] == 354

    def test_filter_by_mode_unexpected_mode(self):
        df = pd.DataFrame()
        fake_mode = SequenceType.NUCLEOTIDE

        with self.assertRaisesRegex(ValueError, "Unexpected mode type"):
            _filter_by_mode(df, fake_mode)

    def test_filter_by_mode_empty_subset(self):
        result = _filter_by_mode(self.full_table_df_2, FragmentMode.BEST_SCORE)
        assert result.empty

    def test_filter_full_table_df_min_len_and_score(self):
        result = _filter_full_table_df(
            self.full_table_df_1, min_len=1200, min_score=1600
        )
        assert len(result) == 1
        assert result["busco_id"].iloc[0] == "44125at2759"

    def test_filter_full_table_df_drop_missing(self):
        df = _filter_full_table_df(
            full_table_df=self.full_table_df_2, drop_missing=True
        )
        assert all(df["status"] == "Complete")

    def test_filter_full_table_df_keep_missing(self):
        df = _filter_full_table_df(
            full_table_df=self.full_table_df_2, drop_missing=False
        )
        assert set(df["status"]) == {"Complete", "Missing"}

    def test_get_corresponding_files_nucl(self):
        busco_seq_dir = Path("/fake/busco_dir")
        df = pd.DataFrame(
            {
                "busco_id": ["A1", "B2", "C3"],
                "status": ["Complete", "Duplicated", "Fragmented"],
            }
        )
        paths = _get_corresponding_files(
            df, busco_seq_dir, seq_type=SequenceType.NUCLEOTIDE
        )
        expected = [
            busco_seq_dir / "single_copy_busco_sequences" / "A1.fna",
            busco_seq_dir / "multi_copy_busco_sequences" / "B2.fna",
            busco_seq_dir / "fragmented_busco_sequences" / "C3.fna",
        ]
        assert paths == expected

    def test_get_corresponding_files_prot(self):
        busco_seq_dir = Path("/fake/busco_dir")
        df = pd.DataFrame(
            {
                "busco_id": ["A1", "B2", "C3"],
                "status": ["Complete", "Duplicated", "Fragmented"],
            }
        )
        paths = _get_corresponding_files(
            df, busco_seq_dir, seq_type=SequenceType.PROTEIN
        )
        expected = [
            busco_seq_dir / "single_copy_busco_sequences" / "A1.faa",
            busco_seq_dir / "multi_copy_busco_sequences" / "B2.faa",
            busco_seq_dir / "fragmented_busco_sequences" / "C3.faa",
        ]
        assert paths == expected

    def test_get_corresponding_files_empty_df(self):
        busco_seq_dir = Path("/fake/busco_dir")
        df = pd.DataFrame()
        paths = _get_corresponding_files(
            df, busco_seq_dir, seq_type=SequenceType.NUCLEOTIDE
        )
        expected = []
        assert paths == expected

    def test_normalize_fasta_header_replaces_spaces(self):
        result = _normalize_fasta_header(
            "gene 1 with spaces",
        )
        assert result == "gene_1_with_spaces\n"

    def test_normalize_fasta_header_replaces_special(self):
        result = _normalize_fasta_header("gene@#!^&|1")
        assert result == "gene______1\n"

    def test_normalize_fasta_header_with_species_tag(self):
        result = _normalize_fasta_header("gene_1", "Escherichia_coli")
        assert result == "Escherichia_coli|gene_1\n"

    def test_normalize_fasta_header_custom_separator(self):
        result = _normalize_fasta_header("gene_1", "Escherichia_coli", ":")
        assert result == "Escherichia_coli:gene_1\n"

    def test_normalize_fasta_header_disallows_separator(self):
        result = _normalize_fasta_header("gene|1")
        assert result == "gene_1\n"

    def test_normalize_fasta_header_custom_replacement_char(self):
        result = _normalize_fasta_header(header="gene 1", replacement_char=";")
        assert result == "gene;1\n"

    def test_normalize_fasta_header_custom_allowed_char(self):
        result = _normalize_fasta_header(
            header="gene@1#test$", allowed_chars=r"A-Za-z0-9_.\-\@\$"
        )
        assert result == "gene@1_test$\n"

    def test_append_uscos_returns_correct_type(self):
        result_nucl = _append_uscos(GenesDirectoryFormat(), [], SequenceType.NUCLEOTIDE)
        result_prot = _append_uscos(ProteinsDirectoryFormat(), [], SequenceType.PROTEIN)
        self.assertIsInstance(result_nucl, GenesDirectoryFormat)
        self.assertIsInstance(result_prot, ProteinsDirectoryFormat)

    def test_append_uscos_interleaved_seqs(self):
        fp = Path(self.get_data_path("usco_sequences/9802219at2.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(), usco_paths, SequenceType.NUCLEOTIDE, wrap_column=80
        )
        contents = (Path(str(out_dir)) / "9802219at2.fasta").read_text()
        assert len(contents.split("\n", 2)[1]) == 80

    def test_append_uscos_non_interleaved_seqs(self):
        fp = Path(self.get_data_path("usco_sequences/9802219at2.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(), usco_paths, SequenceType.NUCLEOTIDE, wrap_column=0
        )
        contents = (Path(str(out_dir)) / "9802219at2.fasta").read_text()
        assert contents.split("\n", 3)[2].startswith(">")

    def test_append_uscos_lower(self):
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(),
            usco_paths,
            SequenceType.NUCLEOTIDE,
            case_mode=CaseMode.LOWER,
        )
        contents = (Path(str(out_dir)) / "tiny_gene_1.fasta").read_text()
        assert contents.split("\n", 2)[1].islower()

    def test_append_uscos_upper(self):
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(),
            usco_paths,
            SequenceType.NUCLEOTIDE,
            case_mode=CaseMode.UPPER,
        )
        contents = (Path(str(out_dir)) / "tiny_gene_1.fasta").read_text()
        assert contents.split("\n", 2)[1].isupper()

    def test_append_uscos_preserve(self):
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(),
            usco_paths,
            SequenceType.NUCLEOTIDE,
            case_mode=CaseMode.PRESERVE,
        )
        contents = (Path(str(out_dir)) / "tiny_gene_1.fasta").read_text()
        assert contents.split("\n", 2)[1].isupper()

    def test_append_uscos_invalid_case_mode_raises(self):
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        with self.assertRaisesRegex(ValueError, "Invalid case_mode:"):
            _append_uscos(
                GenesDirectoryFormat(),
                usco_paths,
                seq_type=SequenceType.NUCLEOTIDE,
                case_mode="invalid",
            )

    def test_append_uscos_creates_new_file(self):
        out_dir = GenesDirectoryFormat()
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        result = _append_uscos(out_dir, usco_paths, seq_type=SequenceType.NUCLEOTIDE)
        output_file = Path(str(result)) / "tiny_gene_1.fasta"
        assert output_file.exists()

    def test_append_uscos_appends_to_existing_file(self):
        out_dir = GenesDirectoryFormat()
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        _append_uscos(out_dir, usco_paths, seq_type=SequenceType.NUCLEOTIDE)
        _append_uscos(out_dir, usco_paths, seq_type=SequenceType.NUCLEOTIDE)
        output_file = Path(str(out_dir)) / "tiny_gene_1.fasta"
        contents = output_file.read_text()
        assert contents.count(">") == 4

    def test_append_uscos_species_tag_and_separator(self):
        fp = Path(self.get_data_path("usco_sequences/tiny_gene_1.fna"))
        usco_paths = [fp]
        out_dir = _append_uscos(
            GenesDirectoryFormat(),
            usco_paths,
            seq_type=SequenceType.NUCLEOTIDE,
            species_tag="Escherichia_coli",
            species_separator="|",
        )
        contents = (Path(str(out_dir)) / "tiny_gene_1.fasta").read_text()
        assert contents.startswith(">Escherichia_coli|")

    def test_append_uscos_with_description(self):
        fp = Path(
            self.get_data_path("uscos_with_description/usco_with_description_1.fna")
        )

        usco_dir = GenesDirectoryFormat()

        out_dir = _append_uscos(
            usco_dir,
            [fp],
            seq_type=SequenceType.NUCLEOTIDE,
            species_tag="Escherichia_coli",
            species_separator="|",
        )

        contents = (
            (Path(str(out_dir)) / "usco_with_description_1.fasta")
            .read_text()
            .splitlines()
        )
        header_lines = [line for line in contents if line.startswith(">")]

        assert header_lines[0].startswith(
            ">Escherichia_coli|gene1_some_description_text"
        )

    def test_extract_uscos_nonexisting_dir(self):
        with self.assertRaisesRegex(
            NotADirectoryError, "BUSCO results directory does not exist or"
        ):
            _extract_uscos(
                busco_results_dir=Path("not_a_dir"), lineage_dataset="lineage_1"
            )

    @patch("q2_mag.busco.extract_orthologs._parse_full_table_file")
    @patch("q2_mag.busco.extract_orthologs._filter_full_table_df")
    @patch("q2_mag.busco.extract_orthologs._get_corresponding_files")
    def test_extract_uscos_basic(
        self,
        mock_get_files,
        mock_filter_df,
        mock_parse_file,
    ):
        df_mock = pd.DataFrame(
            [
                {
                    "busco_id": "gene1",
                    "status": "Complete",
                    "Score": 250,
                    "Length": 500,
                },
                {
                    "busco_id": "gene2",
                    "status": "Fragmented",
                    "Score": 100,
                    "Length": 200,
                },
            ]
        )
        mock_parse_file.return_value = df_mock

        mock_filter_df.return_value = df_mock[df_mock["status"] == "Complete"]

        mock_get_files.return_value = [Path("/fake/path/gene1.fna")]

        result = _extract_uscos(
            busco_results_dir=self.busco_dir,
            lineage_dataset="bacteria_odb10",
            seq_type=SequenceType.NUCLEOTIDE,
            min_len=0,
            min_score=0,
            fragment_mode=FragmentMode.SKIP,
            duplicate_mode=DuplicateMode.SKIP,
        )

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], Path("/fake/path/gene1.fna"))

        analysis_run_dir = self.busco_dir / "run_bacteria_odb10"
        full_table_fp = analysis_run_dir / "full_table.tsv"
        busco_sequences_dir = analysis_run_dir / "busco_sequences"

        mock_parse_file.assert_called_once_with(full_table_fp)
        mock_filter_df.assert_called_once_with(
            df_mock,
            min_len=0,
            min_score=0,
            drop_missing=True,
            fragment_mode=FragmentMode.SKIP,
            duplicate_mode=DuplicateMode.SKIP,
        )
        mock_get_files.assert_called_once_with(
            mock_filter_df.return_value,
            busco_sequences_dir=busco_sequences_dir,
            seq_type=SequenceType.NUCLEOTIDE,
        )

    def test_filter_full_table_accepts_action_parameter_strings(self):
        observed = _filter_full_table_df(
            self.full_table_df_1, fragment_mode="best_score", duplicate_mode="longest"
        )
        expected = _filter_full_table_df(
            self.full_table_df_1,
            fragment_mode=FragmentMode.BEST_SCORE,
            duplicate_mode=DuplicateMode.LONGEST,
        )
        pd.testing.assert_frame_equal(observed, expected)

    def test_selected_hits_are_the_only_records_written(self):
        root = Path(self.temp_dir.name)
        run = root / "run_test"
        run.mkdir()
        table = run / "full_table.tsv"
        for status, folder, option in (
            ("Duplicated", "multi_copy_busco_sequences", "duplicate_mode"),
            ("Fragmented", "fragmented_busco_sequences", "fragment_mode"),
        ):
            sequence_dir = run / "busco_sequences" / folder
            sequence_dir.mkdir(parents=True)
            table.write_text(
                "# BUSCO version\n# lineage\n"
                "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand"
                "\tScore\tLength\n"
                f"marker\t{status}\tcontig\t1\t12\t+\t100\t4\n"
                f"marker\t{status}\tcontig\t30\t13\t-\t50\t6\n"
            )
            for seq_type, extension, short, long in (
                (SequenceType.NUCLEOTIDE, "fna", "ATGC", "ATGCAT"),
                (SequenceType.PROTEIN, "faa", "MACK", "MACKLG"),
            ):
                for direct_ids in (False, True):
                    # BUSCO can retain the full hit ID in Sequence, or split
                    # the contig name and coordinates into separate columns.
                    contents = table.read_text()
                    if direct_ids:
                        contents = contents.replace(
                            "\tcontig\t1\t", "\tcontig:1-12|+\t1\t"
                        ).replace("\tcontig\t30\t", "\tcontig:13-30|-\t30\t")
                        table.write_text(contents)
                    (sequence_dir / f"marker.{extension}").write_text(
                        f">contig:1-12|+ first copy\n{short}\n"
                        f">contig:13-30|- second copy\n{long}\n"
                    )
                    for mode, expected, header in (
                        ("longest", long, "contig_13-30_-"),
                        ("best_score", short, "contig_1-12__"),
                        ("skip", None, None),
                    ):
                        with self.subTest(
                            status=status, seq_type=seq_type, mode=mode,
                            direct_ids=direct_ids,
                        ):
                            selected = _extract_uscos(
                                root, "test", seq_type=seq_type, **{option: mode}
                            )
                            output = (ProteinsDirectoryFormat() if seq_type.is_protein()
                                      else GenesDirectoryFormat())
                            _append_uscos(output, selected, seq_type, species_tag="MAG")
                            if expected is None:
                                self.assertEqual(list(Path(str(output)).iterdir()), [])
                            else:
                                lines = (Path(str(output)) / "marker.fasta").read_text(
                                ).splitlines()
                                self.assertEqual(sum(line.startswith(">")
                                                     for line in lines), 1)
                                self.assertTrue(lines[0].startswith(f">MAG|{header}"))
                                self.assertEqual("".join(lines[1:]), expected)
                    if direct_ids:
                        table.write_text(contents.replace(
                            "\tcontig:1-12|+\t", "\tcontig\t"
                        ).replace("\tcontig:13-30|-\t", "\tcontig\t"))

    def test_selected_hit_requires_a_unique_fasta_match(self):
        from q2_mag.busco.extract_orthologs import _SelectedUSCO

        path = Path(self.temp_dir.name) / "marker.fna"
        for contents, count in (
            (">other:1-12\nATGC\n", 0),
            (">contig:1-12\nATGC\n>contig:1-12\nATGC\n", 2),
        ):
            path.write_text(contents)
            selected = _SelectedUSCO(path, "contig", 1, 12, "+")
            with self.subTest(count=count), self.assertRaisesRegex(
                ValueError, f"Expected exactly one FASTA record.*found {count}"
            ):
                _append_uscos(GenesDirectoryFormat(), [selected])
