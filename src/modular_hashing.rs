// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module contains the fundamental functions required for hashing,
//! using lookups. These are the three main functions of Zelbet:
//! Bricks, Concrete and Bars.

extern crate dusk_bls12_381 as BLS;

use crate::constants::{self, decomposition_s_i};
use bigint::U256 as u256;
use BLS::BlsScalar as Scalar;

use constants::{decomposition_inverses_mont, vu256, SboxBLS};

const DECOMPOSITION_LEN: usize = 27;

// Convert representation from tuple in (Z_{s1} x ... x Z_{s_n}) to single
// element
fn compute_whole_representation(
    decomposition: [u256; DECOMPOSITION_LEN],
) -> Scalar {
    // Note that decomposition[53] is s_1, not s_n, so decomposition[1] is s_n
    Scalar::from_raw(
        (0..DECOMPOSITION_LEN)
            .rev()
            .fold(u256::zero(), |single, k| match k > 0 {
                true => (single + decomposition[k]) * decomposition_s_i[k - 1],
                false => single + decomposition[k],
            })
            .0,
    )
}

// S-box used in bar function
fn small_s_box(x: u256) -> u256 {
    match x < vu256 {
        true => SboxBLS[x.as_u32() as usize],
        false => x,
    }
}

// Lookup-table-based Sbox
fn bar(mut state: [Scalar; 3]) {
    let mut nibbles = [u256::zero(); 27];

    for scalar in state.iter_mut() {
        // 1. Decomposition
        // Get state value that we are decomposing in non-Montgomery form (comes
        // in Montgomery form by default due to BLS library; but the
        // modular operations won't work as intended if left like this)
        let mut intermediate = u256(scalar.reduce().0);
        let mut value = u256::zero();

        (0..27).for_each(|k| {
            value = intermediate % decomposition_s_i[k];

            // Reduce intermediate representation
            match k < 26 {
                true => {
                    // Convert to BLS scalar form to make use of fast modular
                    // multiplication (rather than dividing)
                    let intermediate_scalar: Scalar =
                        Scalar::from_raw((intermediate - value).0)
                            * decomposition_inverses_mont[k];
                    intermediate = u256(intermediate_scalar.reduce().0);
                }
                false => value = intermediate,
            };

            // 2. S-box
            nibbles[k] = small_s_box(value);
        });

        // 3. Composition
        *scalar = compute_whole_representation(nibbles);
    }
}

// Element-wise power function
fn brick(state: [Scalar; 3]) -> [Scalar; 3] {
    let mut new_state = [Scalar::zero(); 3];
    let two = Scalar([
        17179869180,
        12756850513266774020,
        3681868479150465002,
        3479420709561305823,
    ]);
    let x_squared = state[0] * state[0];
    new_state[0] = x_squared * x_squared * state[0];
    new_state[1] = state[1] * (&x_squared + state[0] + &two);
    new_state[2] = state[2]
        * ((state[1] * state[1])
            + Scalar([
                25769803770,
                688531696190609414,
                14746174755580473312,
                5219131064341958734,
            ]) * state[1]
            + Scalar([
                34359738360,
                7066956952823996424,
                7363736958300930005,
                6958841419122611646,
            ]));
    new_state
}

// Apply affine transformation to state via MDS matrix multiplication
fn concrete(
    state: [Scalar; 3],
    matrix: [[Scalar; 3]; 3],
    constants: [Scalar; 3],
) -> [Scalar; 3] {
    let mut new_state = constants;

    // matrix multiplication
    for i in 0..3 {
        for j in 0..3 {
            new_state[i] += matrix[i][j] * state[j];
        }
    }

    new_state
}

// Reinforced concrete hash function, taking in the hash parameters and
// three-element item to be hashed, and outputting the hash value (three BLS
// scalar elements)
pub fn zelbet_hash(
    scalar_inputs: [Scalar; 3],
    matrix: [[Scalar; 3]; 3],
    constants: [[Scalar; 3]; 6],
) -> [Scalar; 3] {
    let mut new_state = concrete(scalar_inputs, matrix, constants[0].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[1].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[2].clone());
    bar(new_state);
    new_state = concrete(new_state, matrix, constants[3].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[4].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[5].clone());
    new_state
}

#[cfg(test)]
mod tests {
    use crate::constants::BLS_scalar_decomposition;

    use super::*;
    use constants::{constantsBLS, decomposition_inverses_mont, matrixBLS};

    #[test]
    fn decomposition_inverses_correct() {
        for k in 0..27 {
            let product = Scalar::from_raw(decomposition_s_i[k].0)
                * decomposition_inverses_mont[k];
            assert_eq!(product, Scalar::one());
        }
    }

    #[test]
    fn test_bricks() {
        let input = [Scalar::from(4), Scalar::from(3), Scalar::from(2)];
        let output = brick(input);

        let two = Scalar([
            17179869180,
            12756850513266774020,
            3681868479150465002,
            3479420709561305823,
        ]);
        let three = Scalar([
            25769803770,
            688531696190609414,
            14746174755580473312,
            5219131064341958734,
        ]);
        let four = Scalar([
            34359738360,
            7066956952823996424,
            7363736958300930005,
            6958841419122611646,
        ]);
        let quadratic_x = Scalar::from(4) * Scalar::from(4);
        let quadratic_y = Scalar::from(3) * Scalar::from(3);
        let element_1 = quadratic_x * quadratic_x * Scalar::from(4);
        let element_2 = Scalar::from(3) * (quadratic_x + Scalar::from(4) + two);
        let element_3 =
            Scalar::from(2) * (quadratic_y + (three * Scalar::from(3)) + four);
        let calculated_output = [element_1, element_2, element_3];
        println!("{:?}", output);
        println!("{:?}", calculated_output);
        assert_eq!(output[0], calculated_output[0]);
        assert_eq!(output[1], calculated_output[1]);
        assert_eq!(output[2], calculated_output[2]);
    }

    #[test]
    fn test_concrete() {
        let state = [Scalar::from(4), Scalar::from(3), Scalar::from(2)];
        let matrix = matrixBLS;
        let constants = constantsBLS;
        let output = concrete(state, matrix, constants[0]);

        let copy_matrix = [
            [Scalar::from(2), Scalar::from(1), Scalar::from(1)],
            [Scalar::from(1), Scalar::from(2), Scalar::from(1)],
            [Scalar::from(1), Scalar::from(1), Scalar::from(2)],
        ];

        let mut new_state = constants[0];
        for i in 0..3 {
            for j in 0..3 {
                new_state[i] += copy_matrix[i][j] * state[j];
            }
        }

        assert_eq!(new_state, output);
    }
}
