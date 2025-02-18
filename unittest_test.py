import unittest
import functions

class TestFunctions(unittest.TestCase):
    def test_complementreverse(self):
        """Test reverse complement function."""
        self.assertEqual(functions.Complementreverse("atgc"), "gcat")
        self.assertEqual(functions.Complementreverse("ggcc"), "ccgg")  # Palindromic test

    def test_PatternCount(self):
        """Test PatternCount function."""
        self.assertEqual(functions.PatternCount("atgatgatg", "atg"), 3)
        self.assertEqual(functions.PatternCount("", "a"), 0)  # Edge case: Empty genome

    def test_patternmatching(self):
        """Test pattern matching function."""
        self.assertEqual(functions.patternmatching("atg", "atgatgatg"), [0, 3, 6])

    def test_FrequenceTable(self):
        """Test FrequenceTable function."""
        self.assertEqual(functions.FrequenceTable("atgatgatg", 3), {"atg": 3, "tga": 2, "gat": 2})
        self.assertEqual(functions.FrequenceTable("", 3), {})  # Edge case: Empty input

    def test_FindClumps(self):
        """Test FindClumps function."""
        genome = "atgatcaaggatgatcaag"
        k, L, t = 3, 10, 2
        self.assertEqual(set(functions.FindClumps(genome, k, L, t)), {"atg", "gat", "tga"})

    def test_MaxMap(self):
        """Test MaxMap function."""
        self.assertEqual(functions.MaxMap({"atg": 3, "tga": 2, "gat": 2}), 3)
        self.assertEqual(functions.MaxMap({}), 0)  # Edge case: Empty dictionary

    def test_ImprovedFrequentWords(self):
        """Test ImprovedFrequentWords function."""
        self.assertEqual(set(functions.ImprovedFrequentWords("atgatgatg", 3)), {"atg"})

    
if __name__ == "__main__":
    unittest.main()

