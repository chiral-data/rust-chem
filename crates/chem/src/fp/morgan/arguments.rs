use crate::fp::fingerprint::FingerprintArguments;

/// Morgan fingerprint specific arguments
pub struct MorganArguments {
    pub radius: u32,
    pub count_simulation: bool,
    pub include_chirality: bool,
    pub only_nonzero_invariants: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub include_redundant_environments: bool,
    pub use_bond_types: bool,
}

impl MorganArguments {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius: u32,
        count_simulation: bool,
        include_chirality: bool,
        only_nonzero_invariants: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        include_redundant_environments: bool,
        use_bond_types: bool,
    ) -> Self {
        Self {
            radius,
            count_simulation,
            include_chirality,
            only_nonzero_invariants,
            count_bounds,
            fp_size,
            include_redundant_environments,
            use_bond_types,
        }
    }
}

impl FingerprintArguments for MorganArguments {
    fn info_string(&self) -> String {
        format!(
            "MorganArguments onlyNonzeroInvariants={} radius={}",
            self.only_nonzero_invariants, self.radius
        )
    }

    fn fp_size(&self) -> u32 {
        self.fp_size
    }

    fn count_simulation(&self) -> bool {
        self.count_simulation
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}
