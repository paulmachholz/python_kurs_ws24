def PatternCount(a, b):
    """
    Counts the number of times a pattern appears in a given DNA sequence.

    Args:
        a (str): The DNA sequence.
        b (str): The pattern to count.

    Returns:
        int: The count of occurrences of the pattern in the sequence.
    """
    count = 0
    for i in range(len(a) - len(b) + 1):
        if a[i:i+len(b)] == b:
            count += 1
    return count
        

def FrequentWords(text, k):
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

def FrequenceTable(text, k):
    """
    Generates a frequency table of all k-mers in a given DNA string.

    Args:
        text (str): The DNA sequence.
        k (int): The length of k-mers to analyze.

    Returns:
        dict: A dictionary where keys are k-mers and values are their frequency in the text.
    """
    freqMap = {}
    for i in range(len(text) - k + 1):
        Pattern = text[i:i+k]
        if Pattern in freqMap:
            freqMap[Pattern] += 1
        else:
            freqMap[Pattern] = 1

    return freqMap


def MaxMap(freqMap):
    
    return max(freqMap.values()) if freqMap else 0 



def ImprovedFrequentWords(text, k):
    freqMap = FrequenceTable(text,k)  #Build table
    maxCount = MaxMap(freqMap)  #find max counts
    FrequentPattern = [pattern for pattern, count in freqMap.items() if count == maxCount] #collect patterns

    return FrequentPattern


def Complementreverse(Pattern):
    """
    Computes the reverse complement of a DNA sequence.

    Args:
        pattern (str): The DNA sequence.

    Returns:
        str: The reverse complement of the input DNA sequence.
    """
    complement = {'a':'t', 't':'a', 'c':'g', 'g':'c'} #dictionary to map complement
    Patternrc =''.join(complement[base] for base in reversed(Pattern))  #complement each n and reverse result
    return Patternrc


def patternmatching(Pattern, Genome):
     """
    Finds all starting positions of a given pattern within a DNA sequence.

    Args:
        pattern (str): The pattern to search for.
        genome (str): The DNA sequence.

    Returns:
        list: A list of starting indices where the pattern appears in the genome.
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

