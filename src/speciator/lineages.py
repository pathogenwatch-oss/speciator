from __future__ import annotations

import dataclasses
import logging
import os
import subprocess

# Initialize logger
logger = logging.getLogger(__name__)

_UNKNOWN_CODE = "99999999"
_UNKNOWN_SPECIES_CODE = "32644"

@dataclasses.dataclass
class Lineage:
    """A lineage."""

    LinkCode: str
    SuperkingdomName: str
    PhylumName: str
    ClassName: str
    OrderName: str
    FamilyName: str
    GenusName: str
    SpeciesName: str
    OrganismName: str
    SuperkingdomId: str
    PhylumId: str
    ClassId: str
    OrderId: str
    FamilyId: str
    GenusId: str
    SpeciesId: str
    OrganismId: str

    @staticmethod
    def build(link_code: str, lineage_names: str, lineage_codes: str) -> "Lineage":
        """Builds a standard lineage object with the organism code set to the species code."""
        names = lineage_names.split(";")
        codes = lineage_codes.split(";")
        return Lineage.build_full(link_code, names, codes, names[6], codes[6])

    @staticmethod
    def species_only(name: str, taxon_code: str):
        return Lineage(
            taxon_code,
            None,
            None,
            None,
            None,
            None,
            name.split()[0],
            name,
            name,
            None,
            None,
            None,
            None,
            None,
            None,
            taxon_code,
            taxon_code,
        )

    @staticmethod
    def build_full(
        link_code: str,
        lineage_names: list[str],
        lineage_codes: list[str],
        organism_name: str,
        organism_code,
    ) -> "Lineage":
        lineage_names.append(organism_name)
        lineage_codes.append(organism_code)
        return Lineage(link_code, *lineage_names, *lineage_codes)

    @staticmethod
    def build_unknown_species() -> "Lineage":
        """Represent a novel species that isn't in the taxonomy database."""
        return Lineage(
            _UNKNOWN_SPECIES_CODE,
            "Unidentified",
            "Unknown",
            "Unknown",
            "Unknown",
            "Unknown",
            "Unknown",
            "Unknown",
            "Unknown",
            _UNKNOWN_SPECIES_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
            _UNKNOWN_CODE,
        )

    @staticmethod
    def replace_species(
        base_lineage: Lineage, new_species: Lineage, link_code: str
    ) -> Lineage:
        """This method is used to replace the Raoultella lineage with the Klebsiella lineage"""
        return Lineage(
            link_code,
            base_lineage.SuperkingdomName,
            base_lineage.PhylumName,
            base_lineage.ClassName,
            base_lineage.OrderName,
            base_lineage.FamilyName,
            base_lineage.GenusName,
            new_species.SpeciesName.replace(
                new_species.GenusName, base_lineage.GenusName
            ),
            new_species.OrganismName.replace(
                new_species.GenusName, base_lineage.GenusName
            ),
            base_lineage.SuperkingdomId,
            base_lineage.PhylumId,
            base_lineage.ClassId,
            base_lineage.OrderId,
            base_lineage.FamilyId,
            base_lineage.GenusId,
            new_species.SpeciesId,
            new_species.OrganismId,
        )

    @staticmethod
    def digest(taxon_name: str) -> str:
        """Returns a unique integer representation of the taxon name for use as a taxon code."""
        import hashlib

        hash_bytes = hashlib.sha256(taxon_name.encode("utf-8")).digest()
        hash_int = int.from_bytes(
            hash_bytes[:6], "big"
        )  # Use first 6 bytes for 12 digits
        return str(hash_int).zfill(12)


@dataclasses.dataclass
class TaxonKit:
    exe: str

    def extract_lineages(
        self,
        taxon_ids: set[str],
        add_unknown: bool = True,
        clean_raoultella: bool = True,
    ) -> dict[str, Lineage]:
        """Extracts lineages for the given taxon IDs. If add_unknown is True, a record for an unknown taxon ID will
        be added. Also, the Raoultella lineages is replaced by Klebsiella if clean_raoultella is True."""
        lineages: dict[str, Lineage] = self.run_taxonkit(taxon_ids)
        if add_unknown:
            lineages[str(_UNKNOWN_SPECIES_CODE)] = Lineage.build_unknown_species()
        if clean_raoultella:
            klebsiella_lineage = next(iter(self.run_taxonkit({"573"}).values()))
            replacement_lineages: dict[str, Lineage] = {}
            for taxon_id, lineage in lineages.items():
                if lineage.GenusName == "Raoultella":
                    logger.debug(f"Replacing Raoultella lineage for {taxon_id}")
                    replacement_lineages[taxon_id] = Lineage.replace_species(
                        klebsiella_lineage, lineage, taxon_id
                    )
            lineages = lineages | replacement_lineages
        return lineages

    def run_taxonkit(self, taxon_ids: set[str]) -> dict[str, Lineage]:
        codes_file = "taxon_codes.txt"
        with open(codes_file, "w") as temp:
            temp.write("\n".join(taxon_ids) + "\n")

        tax_lines = subprocess.run(
            f"{self.exe} lineage {codes_file} | awk '$2!=\"\"' | taxonkit reformat2 -t ",
            shell=True,
            capture_output=True,
            text=True,
        )
        lineages: dict[str, Lineage] = {}
        for line in tax_lines.stdout.splitlines():
            (taxon_id, full_lineage, lineage_names, lineage_codes) = line.split("\t")
            lineages[taxon_id] = Lineage.build(taxon_id, lineage_names, lineage_codes)
        os.remove(codes_file)
        return lineages

    def convert_names_to_taxid(self, names: set[str]) -> dict[str, str]:
        """Uses taxonkit to convert the names to taxids."""
        with subprocess.Popen(
            [self.exe, "name2taxid", "--line-buffered"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ) as p:
            for name in names:
                p.stdin.write(f"{name}\n".encode("utf-8"))
            p.stdin.close()
            p.wait()
            mapping: dict[str, str] = {}
            for line in p.stdout.readlines():
                parts = line.decode("utf-8").strip().split("\t")
                if len(parts) == 1:
                    logger.warning(f"{parts[0]} not a recognised species.")
                    # Deal with novel Klebsiella species from Kleborate
                    if parts[0].startswith("Klebsiella"):
                        mapping[parts[0]] = str(Lineage.digest(parts[0]))
                    else:
                        mapping[parts[0]] = _UNKNOWN_SPECIES_CODE
                else:
                    mapping[parts[0]] = parts[1]
        return mapping
