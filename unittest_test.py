import unittest
import functions

class TestFunctions(unittest.TestCase):
    def test_reverse_complement(self):
        """Test reverse complement function."""
        self.assertEqual(functions.reverse_complement("ATGC"), "GCAT")
        self.assertEqual(functions.reverse_complement("CCGG"), "CCGG")  # Palindromic test

    def test_pattern_count(self):
        """Test PatternCount function."""
        self.assertEqual(functions.PatternCount("ATGATGATG", "ATG"), 3)
        self.assertEqual(functions.PatternCount("", "A"), 0)  # Edge case: Empty genome

    def test_pattern_matching(self):
        """Test pattern matching function."""
        self.assertEqual(functions.pattern_matching("ATG", "ATGATGATG"), [0, 3, 6])
        self.assertEqual(functions.pattern_matching("GGT", "AAAAAA"), [])  # No match case

    def test_frequence_table(self):
        """Test FrequenceTable function."""
        self.assertEqual(functions.FrequenceTable("ATGATGATG", 3), {"ATG": 3, "TGA": 2, "GAT": 2})
        self.assertEqual(functions.FrequenceTable("", 3), {})  # Edge case: Empty input

    def test_find_clumps(self):
        """Test FindClumps function."""
        genome = "ATGATCAAGGATGATCAAG"
        k, L, t = 3, 10, 2
        self.assertEqual(set(functions.FindClumps(genome, k, L, t)), {"ATG", "GAT", "TGA"})

    def test_max_map(self):
        """Test MaxMap function."""
        self.assertEqual(functions.MaxMap({"ATG": 3, "TGA": 2, "GAT": 2}), 3)
        self.assertEqual(functions.MaxMap({}), 0)  # Edge case: Empty dictionary

    def test_improved_frequent_words(self):
        """Test ImprovedFrequentWords function."""
        self.assertEqual(set(functions.ImprovedFrequentWords("ATGATGATG", 3)), {"ATG"})

    
if __name__ == "__main__":
    unittest.main()

