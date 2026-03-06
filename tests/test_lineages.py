from speciator import lineages
from speciator.lineages import Lineage, TaxonKit


def test_run_taxonkit_klebsiella_pneumoniae():
    """Test run_taxonkit with Klebsiella pneumoniae returns expected lineage. Assumes taxonkit is available."""
    taxonkit = TaxonKit(exe="taxonkit")
    result = taxonkit.run_taxonkit({"573"})

    assert len(result) == 1
    assert "573" in result

    lineage = result["573"]
    assert lineage.LinkCode == "573"
    assert lineage.SuperkingdomName == "Bacteria"
    assert lineage.PhylumName == "Pseudomonadota"
    assert lineage.ClassName == "Gammaproteobacteria"
    assert lineage.OrderName == "Enterobacterales"
    assert lineage.FamilyName == "Enterobacteriaceae"
    assert lineage.GenusName == "Klebsiella"
    assert lineage.SpeciesName == "Klebsiella pneumoniae"
    assert lineage.OrganismName == "Klebsiella pneumoniae"
    assert lineage.SuperkingdomId == "2"
    assert lineage.PhylumId == "1224"
    assert lineage.ClassId == "1236"
    assert lineage.OrderId == "91347"
    assert lineage.FamilyId == "543"
    assert lineage.GenusId == "570"
    assert lineage.SpeciesId == "573"
    assert lineage.OrganismId == "573"


def test_convert_names_to_taxid():
    """Test conversion of taxon names to taxids. Assumes taxonkit is available."""
    taxonkit = TaxonKit(exe="taxonkit")
    result = taxonkit.convert_names_to_taxid(
        {"Klebsiella pneumoniae", "Klebsiella", "Klebsiella new", "Odd ball"}
    )
    assert len(result) == 4
    assert result["Klebsiella pneumoniae"] == "573"
    assert result["Klebsiella"] == "570"
    assert result["Klebsiella new"] == Lineage.digest("Klebsiella new")
    assert result["Odd ball"] == lineages._UNKNOWN_SPECIES_CODE


def test_extract_lineages():
    taxonkit: TaxonKit = TaxonKit(exe="taxonkit")
    result: dict[str, Lineage] = taxonkit.extract_lineages({"573", "54291", "123456789012"}, add_unknown=True, clean_raoultella=True)
    assert len(result) == 3
    assert result["573"].SpeciesName == "Klebsiella pneumoniae"
    assert result["54291"].GenusName == "Klebsiella"
    assert result["54291"].GenusId == "570"
    assert result["54291"].SpeciesName == "Klebsiella ornithinolytica"
    assert result[lineages._UNKNOWN_SPECIES_CODE].LinkCode == Lineage.build_unknown_species().LinkCode
    assert "123456789012" not in result

