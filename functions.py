def PatternCount(a, b):
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
    complement = {'a':'t', 't':'a', 'c':'g', 'g':'c'} #dictionary to map complement
    Patternrc =''.join(complement[base] for base in reversed(Pattern))  #complement each n and reverse result
    return Patternrc