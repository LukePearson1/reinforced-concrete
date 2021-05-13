// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This file contains the circuit implementation of the
//! zelbet hash function

use crate::constants::{MONTGOMERY_FOUR, MONTGOMERY_THREE, MONTGOMERY_TWO};
use dusk_plonk::constraint_system::Variable;
use dusk_plonk::prelude::*;

/// This function computes the in-circuit brick function,
/// as part of the hashing gadget
pub fn brick_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
) -> [Variable; 3] {
    let two = composer.add_witness_to_circuit_description(MONTGOMERY_TWO);

    let x_squared = composer.mul(
        BlsScalar::one(),
        state[0],
        state[0],
        BlsScalar::one(),
        None,
    );

    let var_one = composer.big_mul(
        BlsScalar::one(),
        x_squared,
        x_squared,
        Some((BlsScalar::one(), state[0])),
        BlsScalar::one(),
        None,
    );

    let tuple = composer.big_add(
        (BlsScalar::one(), x_squared),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::one(),
        None,
    );

    let var_two =
        composer.mul(BlsScalar::one(), state[1], tuple, BlsScalar::one(), None);

    let three = composer.add_witness_to_circuit_description(MONTGOMERY_THREE);

    let four = composer.add_witness_to_circuit_description(MONTGOMERY_FOUR);

    // x3 ·(x2 +α2 ·x2 +β2))

    let tuple_one =
        composer.mul(BlsScalar::one(), three, state[1], BlsScalar::one(), None);
    let tuple_two = composer.big_add(
        (BlsScalar::one(), tuple_one),
        (BlsScalar::one(), four),
        Some((BlsScalar::one(), state[1])),
        BlsScalar::one(),
        None,
    );

    let var_three = composer.mul(
        BlsScalar::one(),
        tuple_two,
        state[2],
        BlsScalar::one(),
        None,
    );

    [var_one, var_two, var_three]
}

#[cfg(test)]
mod tests {
    use super::*;
    use dusk_plonk::constraint_system::StandardComposer;

    #[test]
    fn test_bricks_gadget() {
        let mut composer = StandardComposer::new();
        let one = composer.add_witness_to_circuit_description(BlsScalar::one());
        let output = brick_gadget(&mut composer, &[one, one, one]);
        composer.constrain_to_constant(
            output[2],
            // This Bls is taken from a print of the 
            // same value in modular_hashing. This 
            // is performed as it makes it easier
            // to compare to the obfuscated `variable`
            BlsScalar::from(16384),
            Some(BlsScalar::zero()),
        );
    }
}
