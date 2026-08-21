fn main() {
    // `option_env!("CHEM_BUILD_ID")` is baked in when the crate compiles, and
    // cargo has no way to know the value changed — so without this, an e2e
    // build would happily reuse a crate compiled under a different id and show
    // the wrong one. Which is precisely the stale-build confusion the id exists
    // to prevent.
    println!("cargo:rerun-if-env-changed=CHEM_BUILD_ID");
}
