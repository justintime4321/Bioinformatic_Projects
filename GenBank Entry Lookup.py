from Bio import Entrez

Entrez.email = "justintime4321@gmail.com"

def count_genbank_entries(genus_name, start_date, end_date):
    """
    Counts the number of Nucleotide GenBank entries for a given genus
    published between the specified dates.

    Args:
        genus_name (str): The genus name to search for.
        start_date (str): The start date in YYYY/M/D format.
        end_date (str): The end date in YYYY/M/D format.

    Returns:
        int: The number of GenBank entries found.
    """
    term = f"{genus_name}[Organism] AND {start_date}[PDAT] : {end_date}[PDAT]"
    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    return int(record["Count"])

# Example usage with your provided dates:
genus = "Gyalecta"
start = "2004/08/10"
end = "2011/10/31"

count = count_genbank_entries(genus, start, end)
print(f"Number of GenBank entries for {genus} published between {start} and {end}: {count}")