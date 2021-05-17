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

/// In-circuit concrete function as part of the Zelbet hashing
/// gadget with t = 3 and MDS matrix M = circ(2, 1, 1).
pub fn concrete_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    constants: &[Variable; 3],
) -> [Variable; 3] {
    // out0 = 2*u[0] + u[1] + u[2] + c[0];
    let a0 = composer.big_add(
        (BlsScalar::from(2), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        None,
    );
    let out0 = composer.add(
        (BlsScalar::one(), a0),
        (BlsScalar::one(), constants[0]),
        BlsScalar::zero(),
        None,
    );

    // out1 = u[0] + 2*u[1] + u[2] + c[1];
    let a1 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::from(2), state[1]),
        Some((BlsScalar::one(), state[2])),
        BlsScalar::zero(),
        None,
    );
    let out1 = composer.add(
        (BlsScalar::one(), a1),
        (BlsScalar::one(), constants[1]),
        BlsScalar::zero(),
        None,
    );

    // out2 = u[0] + u[1] + 2*u[2] + c[2];
    let a2 = composer.big_add(
        (BlsScalar::one(), state[0]),
        (BlsScalar::one(), state[1]),
        Some((BlsScalar::from(2), state[2])),
        BlsScalar::zero(),
        None,
    );
    let out2 = composer.add(
        (BlsScalar::one(), a2),
        (BlsScalar::one(), constants[2]),
        BlsScalar::zero(),
        None,
    );

    return [out0, out1, out2];
}

// TODO: verify all functions against python outputs
// for hashing the same values.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::gadget_tester;
    use crate::hashing::zelbet::brick;

    #[ignore]
    #[test]
    fn test_bricks_gadget() {
        let res = gadget_tester(
            |composer| {
                let one = composer
                    .add_witness_to_circuit_description(BlsScalar::one());
                let output = brick_gadget(composer, &[one, one, one]);
                let output_1 = brick([
                    BlsScalar::one(),
                    BlsScalar::one(),
                    BlsScalar::one(),
                ]);
                println!("1 {:?}", output_1[0]);
                println!("2 {:?}", output_1[1]);
                composer.constrain_to_constant(
                    output[0],
                    // This Bls is taken from a print of the
                    // same value in modular_hashing. This
                    // is performed as it makes it easier
                    // to compare to the obfuscated `variable`
                    output_1[0],
                    Some(BlsScalar::zero()),
                );
            },
            32,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_concrete_gadget() {
        let res = gadget_tester(
            |composer| {
                let two = composer
                    .add_witness_to_circuit_description(BlsScalar::from(2));
                let a0 = composer.big_add(
                    (BlsScalar::from(2), two),
                    (BlsScalar::one(), two),
                    Some((BlsScalar::one(), two)),
                    BlsScalar::zero(),
                    None,
                );
                composer.constrain_to_constant(a0, BlsScalar::from(8), None);
            },
            32,
        );
        assert!(res.is_ok());
    }
}
