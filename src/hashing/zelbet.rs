// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module contains the fundamental functions required for hashing,
//! using lookups. These are the three main functions of Zelbet:

mod bar;
mod brick;
mod concrete;

use super::gadgets::*;
use crate::constants::CONSTANTS_BLS;
pub use bar::bar;
pub use brick::brick;
pub use concrete::concrete;
use dusk_plonk::{
    constraint_system::{StandardComposer, Variable},
    prelude::BlsScalar,
};

/// In circuit Zelbet hash
pub fn zelbet_gadget(
    composer: &mut StandardComposer,
    state: &[Variable; 3],
    s_i_decomposition: [Variable; 27],
    constants_for_rounds: [Variable; 18],
    zero: Variable,
    one: Variable,
    two: Variable,
) -> [Variable; 3] {
    let mut constants_vector = [one; 3];

    // Round 1
    let mut round1_constants = [one; 3];
    round1_constants.copy_from_slice(&constants_for_rounds[0..3]);
    let mut item = concrete_gadget(composer, state, &round1_constants);
    item = brick_gadget(composer, &item, two);

    // Round 2
    let mut round2_constants = [one; 3];
    round2_constants.copy_from_slice(&constants_for_rounds[3..6]);
    item = concrete_gadget(composer, &item, &round2_constants);
    item = brick_gadget(composer, &item, two);

    // Round 3
    let mut round3_constants = [one; 3];
    round3_constants.copy_from_slice(&constants_for_rounds[6..9]);
    item = concrete_gadget(composer, &item, &round3_constants);
    // Apply bar function to each entry
    (0..3).for_each(|k| {
        item[k] = bar_gadget(composer, item[k], s_i_decomposition, zero, one, two);
    });

    // Round 4
    let mut round4_constants = [one; 3];
    round4_constants.copy_from_slice(&constants_for_rounds[9..12]);
    item = concrete_gadget(composer, &item, &round4_constants);
    item = brick_gadget(composer, &item, two);

    // Round 5
    let mut round5_constants = [one; 3];
    round5_constants.copy_from_slice(&constants_for_rounds[12..15]);
    item = concrete_gadget(composer, &item, &round5_constants);
    item = brick_gadget(composer, &item, two);

    // Final concrete
    let mut round6_constants = [one; 3];
    round6_constants.copy_from_slice(&constants_for_rounds[15..18]);
    item = concrete_gadget(composer, &item, &round6_constants);

    item
}

/// Reinforced concrete hash function, taking in the hash parameters and
/// three-element item to be hashed, and outputting the hash value (three BLS
/// scalar elements)
pub fn zelbet_out_of_circuit(scalar_inputs: [BlsScalar; 3]) -> [BlsScalar; 3] {
    let mut new_state = concrete(scalar_inputs, CONSTANTS_BLS[0].clone());
    new_state = brick(new_state);
    new_state = concrete(new_state, CONSTANTS_BLS[1]);
    new_state = brick(new_state);
    new_state = concrete(new_state, CONSTANTS_BLS[2]);
    bar(&mut new_state);
    new_state = concrete(new_state, CONSTANTS_BLS[3]);
    new_state = brick(new_state);
    new_state = concrete(new_state, CONSTANTS_BLS[4]);
    new_state = brick(new_state);
    new_state = concrete(new_state, CONSTANTS_BLS[5]);
    new_state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{constants::S_I_DECOMPOSITION_MONTGOMERY, gadget_tester};
    use dusk_plonk::plookup::PlookupTable4Arity;

    #[test]
    fn test_zelbet_gadget_circuit() {
        let res = gadget_tester(
            |composer| {
                let hash_table = PlookupTable4Arity::create_hash_table();
                composer.append_lookup_table(&hash_table);
                let zero = composer.add_input(BlsScalar::zero());
                let one = composer.add_input(BlsScalar::one());
                let two = composer.add_input(BlsScalar::from(2));
                let mut s_i_decomposition = [one; 27];
                (0..27).for_each(|k| {
                    s_i_decomposition[k] = composer
                        .add_witness_to_circuit_description(
                            S_I_DECOMPOSITION_MONTGOMERY[k],
                        );
                });
                let mut constants_for_rounds = [one; 18];
                (0..6).for_each(|k| {
                    (0..3).for_each(|j| {
                        constants_for_rounds[3 * k + j] =
                            composer.add_input(CONSTANTS_BLS[k][j]);
                    })
                });
                let _result = zelbet_gadget(
                    composer,
                    &[one; 3],
                    s_i_decomposition,
                    constants_for_rounds,
                    zero,
                    one,
                    two,
                );

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
