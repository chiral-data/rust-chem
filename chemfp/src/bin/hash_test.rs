fn main() {
    // Simulate GPU hash for layer 1, atom 0
    let mut invar = 0u32; // layer 0

    // Hash with current invariant
    chemfp::hash::hash_combine(&mut invar, &37324u32);
    println!("After hash_combine(layer, 37324): {}", invar);

    // Two neighbors, both (12, 37324)
    chemfp::hash::hash_pair(&mut invar, 12i32, 37324u32);
    println!("After first neighbor (12, 37324): {}", invar);

    chemfp::hash::hash_pair(&mut invar, 12i32, 37324u32);
    println!("After second neighbor (12, 37324): {}", invar);

    println!("\nFinal: {}", invar);
    println!("Bit position: {}", invar % 2048);
}
