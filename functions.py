import matplotlib.pyplot as plt
from Bio import SeqIO

def PatternCount(a: str, b: str) -> int:
    """Counts occurences of a pattern in a string
    Parameters:
        a : str -> DNA sequence

        b : str -> pattern to count

        Returns: -> int, count of occurences in the sequence
    """
    count = 0
    for i in range(len(a) - len(b) + 1):
        if a[i:i+len(b)] == b:
            count += 1
    return count       

def FrequentWords(text: str, k: str) -> list[str]:
    """Find the most frequent k-mers in a given DNA sequence
    text : str -> DNA sequence
    k : int -> length of k-mer to analyze

    Returns: -> list of most frequent k-mers

    """
    FrequentPatterns = set() #stores unique freq. k-mers
    count = []  #array
    for i in range(len(text) - k+1):
        Pattern = text[i:i+k]   #extract k-mers
        count.append(PatternCount(text, Pattern)) #stores k-mer in arrray
    maxCount = max(count) #finds maximum count
    for i in range(len(text) - k + 1):
        if count[i] == maxCount:     #identifying highest count
            FrequentPatterns.add(text[i:i+k]) #add k-mer to set

    return list(FrequentPatterns)

def FrequenceTable(text: str, k: int) -> dict[str, int]:
    """Generates a table of all k-mers in a DNA string

    text : str -> DNA sequence
    k : int -> length of k-mers to analyze

    Returns -> dictionary with kmers as keys and frequency as  value
    """
    freqMap = {}
    for i in range(len(text) - k + 1):
        Pattern = text[i:i+k]
        if Pattern in freqMap:
            freqMap[Pattern] += 1
        else:
            freqMap[Pattern] = 1

    return freqMap


def MaxMap(freqMap: dict[str, int]) -> int: 
    """ Finds maximum frequency value in a given k-mer frequency table
    freqMap : dict -> k-mer as key, count as value

    Returns -> the highest frequency found in the table (int)

    """
    return max(freqMap.values()) if freqMap else 0 

def ImprovedFrequentWords(text: str, k: int) -> list[str]:
    freqMap = FrequenceTable(text,k)  #Build table
    maxCount = MaxMap(freqMap)  #find max counts
    FrequentPattern = [pattern for pattern, count in freqMap.items() if count == maxCount] #collect patterns

    return FrequentPattern

def Complementreverse(Pattern: str) -> str:
   
    complement = {'a':'t', 't':'a', 'c':'g', 'g':'c'} #dictionary to map complement
    Patternrc =''.join(complement[base] for base in reversed(Pattern))  #complement each n and reverse result
    return Patternrc

def patternmatching(Pattern: str, Genome: str) -> list[int]:
    """Finds all starting positions of a given pattern within a sequence
        pattern : str -> pattern to search for
        genome : str -> DNA sequence

        Returns -> a list of starting indices where the pattern appears
    """
    position = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == Pattern:
            position.append(i)   #finds position of chosen pattern

    reverse_complement = Complementreverse(Pattern) #finds reverse complement
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == reverse_complement:
            position.append(i)
            
    return position

# day 2

def FindClumps(genome: str, k: int, L: int, t: int) -> list[str]:
    """ Find all k-mers that form (L, t)-clumps in a given genome
    genome : str -> DNA sequece
    k : str -> length of k-mer to search for
    L : int -> window length in which to search for clumps
    t : int -> min. number of times a k-mer has to appear in the window

    Returns -> List of distinct k-mers that form (L, t)-clumps.
    """
    patterns = set()
    n = len(genome)

    for i in range(n - L + 1):          #window slides over genome
        window = genome[i:i + L]
        freqMap = FrequenceTable(window, k)
        for kmer, count in freqMap.items(): #if k-mer inside window is there t-times
            if count >= t:
                patterns.add(kmer)          #it gets added to a list
    return list(patterns)

def read_genome_file(filename: str) -> str:
    """Reads a genome fromn .txt file and returns it as a single string
    filename : str -> path to genome file

    Returns -> sequence as a single string
    """
    with open(filename, 'r') as file:
        genome = file.read().replace('\n', '').replace('\r', '')  #opens text file and deletes empty lines or spaces left -> important for reading text documents created on dos to unix for example
    return genome

def print_clumps(clumps: list[str], per_line: int = 10) -> None:
    """
    Prints k-mers found in clumps in a nicely formatted way.

    Args:
        clumps (list[str]): List of k-mers found as clumps.
        per_line (int): Number of k-mers to print per line (default: 10).
    """
    print("\nFound k-mers forming clumps:\n")
    for i in range(0, len(clumps), per_line):  #lists the found k-mers in lines of 10
        print(" ".join(clumps[i:i+per_line]))

#day 3

def skew_compute(genome: str) -> list[int]:
    """Computes skew array for given genome

    genome : str -> DNA sequence

    Returns -> List of skew values at each position
    """
    skew = [0]
    for base in genome:
        if base == "G":
            skew.append(skew[-1] + 1)
        elif base == "C":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

def skew_minimum(genome: str) -> list[int]:
    """ Finds all positions where the skew is at minimum value
    genome : str -> DNA sequence
    Returns -> List of positions with minimal skew
    """
    skew_array = skew_compute(genome)
    min_skew = min(skew_array)
    return [i for i, value in enumerate(skew_array) if value == min_skew]
    
def plot_skew(genome: str, output_file: str):
    """Creates a skew diagram and saves it as a PNG file."""
    # Calculate the skew array
    skew_array = skew_compute(genome)
    # Plot the skew values
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(skew_array)), skew_array, label="Skew", color='blue')
    # Mark minimum skew positions
    min_positions = skew_minimum(genome)
    plt.scatter(min_positions, [skew_array[i] for i in min_positions], color='red', zorder=5, label="Min Skew Positions")
    # Labels and title
    plt.title("Skew Diagram")
    plt.xlabel("Position in Genome")
    plt.ylabel("Skew Value")
    plt.legend(loc="upper right")
    # Save the plot to a PNG file
    plt.savefig(output_file)
    plt.close()
    print(f"Skew diagram saved to {output_file}")

def plot_skew_pdf(fasta_file: str, output_pdf: str):
    """Creates skew diagrams for all sequences in the FASTA file and saves them in a single PDF."""
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Create a PdfPages object to save the plots to a PDF
    with PdfPages(output_pdf) as pdf:
        # Parse the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            genome = str(record.seq)  # Get the genome as a string
            skew_array = skew_compute(genome)
            plt.figure(figsize=(10, 6))
            plt.plot(range(len(skew_array)), skew_array, label="Skew", color='blue')
            min_positions = skew_minimum(genome)
            plt.scatter(min_positions, [skew_array[i] for i in min_positions], color='red', zorder=5, label="Min Skew Positions")
            plt.title(f"Skew Diagram for {record.id}")
            plt.xlabel("Position in Genome")
            plt.ylabel("Skew Value")
            plt.legend(loc="upper right")
            pdf.savefig()
            plt.close()
    print(f"Skew diagrams saved to {output_pdf}")

# day 4

def hamming_distance(p: str, q: str) -> int:
    """Computes the hamming distance between two equal length strings
    p : str -> first dna sequence
    q: str -> second dna sequence
    Returns -> number of differing positions between two sequences
    """
    if len(p) !=len(q):
        raise ValueError("Strings must be equal length")
    
    return sum(1 for i in range(len(p)) if p[i] != q[i])

