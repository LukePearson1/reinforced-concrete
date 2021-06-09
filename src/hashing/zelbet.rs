// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module contains the fundamental functions required for hashing,
//! using lookups. These are the three main functions of Zelbet:
//! Bricks, Concrete and Bars.

mod bar;
mod brick;
mod concrete;

pub use bar::bar;
pub use brick::brick;
pub use concrete::concrete;
// use dusk_bls12_381::BlsScalar;
use super::gadgets::*;
use crate::constants::{
    CONSTANTS_BLS, MONTGOMERY_TWO, S_I_DECOMPOSITION_MONTGOMERY,
};
use dusk_plonk::{
    constraint_system::{StandardComposer, Variable},
    prelude::BlsScalar,
};

/// In circuit Zelbet hash
pub fn zelbet(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
) -> [Variable; 3] {
    let one = composer.add_witness_to_circuit_description(BlsScalar::one());
    let two = composer.add_witness_to_circuit_description(MONTGOMERY_TWO);
    let mut constants_vector = [one; 3];

    // Constants vector used for each concrete round must be added as variables
    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[0][k]);
    });
    let mut item = concrete_gadget(composer, state, &constants_vector);
    item = brick_gadget(composer, &item, two);

    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[1][k]);
    });
    item = concrete_gadget(composer, &item, &constants_vector);
    item = brick_gadget(composer, &item, two);

    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[2][k]);
    });
    item = concrete_gadget(composer, &item, &constants_vector);
    // Initialise the s_i values as variables
    let mut s_i_decomposition = [state[0]; 27];
    (0..27).for_each(|k| {
        s_i_decomposition[k] = composer.add_witness_to_circuit_description(
            S_I_DECOMPOSITION_MONTGOMERY[k],
        );
    });
    // Apply bar function to each entry
    (0..3).for_each(|k| {
        item[k] = bar_gadget(composer, item[k], s_i_decomposition, one, two);
    });

    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[3][k]);
    });
    item = concrete_gadget(composer, &item, &constants_vector);
    item = brick_gadget(composer, &item, two);

    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[4][k]);
    });
    item = concrete_gadget(composer, &item, &constants_vector);
    item = brick_gadget(composer, &item, two);

    (0..3).for_each(|k| {
        constants_vector[k] = composer.add_input(CONSTANTS_BLS[5][k]);
    });
    item = concrete_gadget(composer, &item, &constants_vector);

    item
}

/// Reinforced concrete hash function, taking in the hash parameters and
/// three-element item to be hashed, and outputting the hash value (three BLS
/// scalar elements)
pub fn zelbet_out_of_circuit(
    scalar_inputs: [BlsScalar; 3],
    matrix: [[BlsScalar; 3]; 3],
    constants: [[BlsScalar; 3]; 6],
) -> [BlsScalar; 3] {
    let mut new_state = concrete(scalar_inputs, matrix, constants[0].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[1].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[2].clone());
    bar(&mut new_state);
    new_state = concrete(new_state, matrix, constants[3].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[4].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, matrix, constants[5].clone());
    new_state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gadget_tester;
    use dusk_plonk::plookup::PlookupTable4Arity;

    #[test]
    fn test_zelbet_circuit() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let one = composer.add_input(BlsScalar::one());
                let _result = zelbet(composer, &[one; 3]);

                let one_eight_seven = composer.add_input(BlsScalar::from(187));
                let zero = composer.add_input(BlsScalar::from(0));

                // Plookup is designed to not pass if the number of plookup
                // checks is much smaller than the size of the lookup table, so
                // these extra checks are to increase the number of plookup
                // gates
                (0..500).for_each(|_| {
                    composer.plookup_gate(
                        one,
                        zero,
                        one_eight_seven,
                        Some(one),
                        BlsScalar::zero(),
                    );
                    composer.plookup_gate(
                        one,
                        one,
                        one,
                        Some(one),
                        BlsScalar::zero(),
                    );
                });
            },
            3000,
        );
        assert!(res.is_ok());
    }
}
