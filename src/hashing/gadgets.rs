// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This file contains the circuit implementation of the
//! zelbet hash function

use crate::constants::MONTGOMERY_TWO;
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
        BlsScalar::zero(),
        None,
    );
    let x_fourth = composer.mul(
        BlsScalar::one(),
        x_squared,
        x_squared,
        BlsScalar::zero(),
        None,
    );
    let var_one = composer.mul(
        BlsScalar::one(),
        x_fourth,
        state[0],
        BlsScalar::zero(),
        None,
    );

    let tuple = composer.big_add(
        (BlsScalar::one(), x_squared),
        (BlsScalar::one(), state[0]),
        Some((BlsScalar::one(), two)),
        BlsScalar::zero(),
        None,
    );
    let var_two = composer.mul(
        BlsScalar::one(),
        state[1],
        tuple,
        BlsScalar::zero(),
        None,
    );

    // x3 ·(x2^2 +α2 ·x2 +β2))
    let y_squared_plus_4 = composer.mul(
        BlsScalar::one(),
        state[1],
        state[1],
        BlsScalar::from(4),
        None,
    );
    let tuple_one = composer.big_add(
        (BlsScalar::one(), y_squared_plus_4),
        (BlsScalar::from(3), state[1]),
        None,
        BlsScalar::zero(),
        None,
    );
    let var_three = composer.mul(
        BlsScalar::one(),
        tuple_one,
        state[2],
        BlsScalar::zero(),
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

    [out0, out1, out2]
}

// TODO: verify all functions against python outputs
// for hashing the same values.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{constants::MATRIX_BLS, hashing::zelbet::brick};
    use crate::{gadget_tester, hashing::zelbet::concrete};

    #[test]
    fn test_bricks_gadget() {
        let res = gadget_tester(
            |composer| {
                let two = composer
                    .add_witness_to_circuit_description(BlsScalar::from(2));
                let three = composer
                    .add_witness_to_circuit_description(BlsScalar::from(3));
                let four = composer
                    .add_witness_to_circuit_description(BlsScalar::from(4));
                let output = brick_gadget(composer, &[two, three, four]);
                let output_1 = brick([
                    BlsScalar::from(2),
                    BlsScalar::from(3),
                    BlsScalar::from(4),
                ]);
                // Check in circuit result against out of circuit result
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        output_1[i],
                        None,
                    );
                }

                // Check in circuit result against expected result
                let expected_result = [
                    BlsScalar::from(32),
                    BlsScalar::from(24),
                    BlsScalar::from(88),
                ];
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        expected_result[i],
                        None,
                    );
                }
            },
            32,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_concrete_gadget() {
        let res = gadget_tester(
            |composer| {
                let one = composer
                    .add_witness_to_circuit_description(BlsScalar::one());
                let two = composer
                    .add_witness_to_circuit_description(BlsScalar::from(2));
                let three = composer
                    .add_witness_to_circuit_description(BlsScalar::from(3));
                let output =
                    concrete_gadget(composer, &[one, two, three], &[two; 3]);
                let output_1 = concrete(
                    [BlsScalar::one(), BlsScalar::from(2), BlsScalar::from(3)],
                    MATRIX_BLS,
                    [BlsScalar::from(2); 3],
                );

                // Check in circuit result against out of circuit result
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        output_1[i],
                        None,
                    );
                }

                // Check in circuit result against expected result
                let expected_result = [
                    BlsScalar::from(9),
                    BlsScalar::from(10),
                    BlsScalar::from(11),
                ];
                for i in 0..3 {
                    composer.constrain_to_constant(
                        output[i],
                        expected_result[i],
                        None,
                    );
                }
            },
            32,
        );
        assert!(res.is_ok());
    }
}
