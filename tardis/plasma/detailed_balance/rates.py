class RadiativeRatesSolver:
    def __init__(self, einstein_coefficients):
        assert einstein_coefficients.index.names == [
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ]
        assert {"A_ul", "B_ul", "B_lu", "nu"} - set(
            einstein_coefficients.columns
        ) == set()
        self.einstein_coefficients = einstein_coefficients.sort_index()

    def solve(self, mean_intensities):
        pass
