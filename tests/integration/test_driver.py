import tempfile
import shutil
from pathlib import Path
from ms_matching import driver
import sqlite3

def test_run_index_peptides_from_fasta():
    # Paths
    base = Path(__file__).parent.parent.parent
    config_path = base / "config" / "test.yaml"
    fasta_path = base / "data" / "testdata" / "ref" / "test.fa"  # substitute with a test FASTA if needed

    # Create a temp directory for DB file
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        #outname = Path(tmpdir) / "out.tsv"

        # Run driver
        fragments_info = driver.run(
            config_path=str(config_path),
            db_path=str(db_path),
            fasta_path=str(fasta_path),
            #outname=str(outname)
        )
        info1 = next(fragments_info)
        # import pdb; pdb.set_trace()
        1+1
        for key in ('fragments', 'peptide', 'aa_mass', 'mod_aliases', 'bin_width'):
            assert key in info1

        # # Assertions
        # assert db_path.exists()

        # # Check DB contents
        # with sqlite3.connect(db_path) as conn:
        #     cur = conn.cursor()
        #     cur.execute("SELECT COUNT(*) FROM peptide_index")
        #     count = cur.fetchone()[0]
        #     import pdb; pdb.set_trace()
        #     assert count > 0, "Expected peptides to be indexed in peptide_index"

        #     cur.execute("SELECT COUNT(*) FROM mod_definitions")
        #     mod_count = cur.fetchone()[0]
        #     assert mod_count > 0, "Expected mod definitions to be present"

        # print(f"Test ran successfully with {count} peptides and {mod_count} mods.")

def _xxtest_driver_run_end_to_end():
    # Paths
    base = Path(__file__).parent.parent.parent
    config_path = base / "config" / "test.yaml"
    fasta_path = base / "data" / "testdata" / "ref" / "test.fa"  # substitute with a test FASTA if needed

    # Create a temp directory for DB file
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        outname = Path(tmpdir) / "out.tsv"

        # Run driver
        driver.run(
            config_path=str(config_path),
            db_path=str(db_path),
            fasta_path=str(fasta_path),
            outname=str(outname)
        )

        # # Assertions
        # assert db_path.exists()

        # # Check DB contents
        # with sqlite3.connect(db_path) as conn:
        #     cur = conn.cursor()
        #     cur.execute("SELECT COUNT(*) FROM peptide_index")
        #     count = cur.fetchone()[0]
        #     import pdb; pdb.set_trace()
        #     assert count > 0, "Expected peptides to be indexed in peptide_index"

        #     cur.execute("SELECT COUNT(*) FROM mod_definitions")
        #     mod_count = cur.fetchone()[0]
        #     assert mod_count > 0, "Expected mod definitions to be present"

        # print(f"Test ran successfully with {count} peptides and {mod_count} mods.")

