// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constants::{DECOMPOSITION_S_I, INVERSES_S_I, VU_256, SBOX_BLS};
use dusk_bls12_381::BlsScalar as Scalar;
use bigint::U256 as u256;

const DECOMPOSITION_LEN: usize = 27;

// Convert representation from tuple in (Z_{s_n} x ... x Z_{s_1}) to single
// element
fn compute_whole_representation(
    decomposition: [u256; DECOMPOSITION_LEN],
) -> Scalar {
    // Note that decomposition_s_i[26] is s_1, so decomposition_s_i[0] is s_27
    Scalar::from_raw(
        (0..DECOMPOSITION_LEN)
            .rev()
            .fold(u256::zero(), |single, k| match k > 0 {
                true => (single + decomposition[k]) * DECOMPOSITION_S_I[k - 1],
                false => single + decomposition[k],
            })
            .0,
    )
}

// S-box used in bar function
fn small_s_box(x: u256) -> u256 {
    match x < VU_256 {
        true => SBOX_BLS[x.as_u32() as usize],
        false => x,
    }
}

// Lookup-table-based Sbox
pub fn bar(state: &mut [Scalar; 3]) {
    let mut nibbles = [u256::zero(); 27];

    for scalar in state.iter_mut() {
        // 1. Decomposition
        // Get state value that we are decomposing in non-Montgomery form (comes
        // in Montgomery form by default due to BLS library; but the
        // modular operations won't work as intended if left like this)
        let mut intermediate = u256(scalar.reduce().0);
        let mut value = u256::zero();

        (0..27).for_each(|k| {
            value = intermediate % DECOMPOSITION_S_I[k];

            // Reduce intermediate representation
            match k < 26 {
                true => {
                    // Convert to BLS scalar form to make use of fast modular
                    // multiplication (rather than dividing)
                    let intermediate_scalar: Scalar =
                        Scalar((intermediate - value).0) * INVERSES_S_I[k];
                    intermediate = u256(intermediate_scalar.0);
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

mod tests {
    use super::*;
    
    #[test]
    fn test_bar() {
        let mut input = [Scalar::one(); 3];
        bar(&mut input);
        let mut breakdown = [u256([248, 0, 0, 0]); 27];
        breakdown[0] = u256([131, 0, 0, 0]);
        let composed = compute_whole_representation(breakdown);
        assert_eq!(input[0], composed);
    }
}