// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constants::{DECOMPOSITION_S_I, INVERSES_S_I, SBOX_U256, VU_256};
use bigint::U256 as u256;
use dusk_plonk::prelude::BlsScalar as Scalar;

const DECOMPOSITION_LEN: usize = 27;

/// Convert representation from tuple in (Z_{s_n} x ... x Z_{s_1}) to single
/// scalar element in Montgomery form (out of circuit)
fn compute_whole_representation(
    decomposition: [u256; DECOMPOSITION_LEN],
) -> Scalar {
    // Note that decomposition_s_i[26] is s_1, so decomposition_s_i[0] is s_27
    Scalar::from_raw(
        (0..DECOMPOSITION_LEN)
            .rev()
            .fold(u256::zero(), |single, k| match k > 0 {
                true => {
                    (single + decomposition[k])
                        * u256(DECOMPOSITION_S_I[k - 1].0)
                }
                false => single + decomposition[k],
            })
            .0,
    )
}

/// S-box used in bar function (out of circuit)
fn small_s_box(x: u256) -> u256 {
    match x < VU_256 {
        true => SBOX_U256[x.as_u32() as usize],
        false => x,
    }
}

/// Bar function (out of circuit)
pub fn bar(state: &mut [Scalar; 3]) {
    let mut nibbles = [u256::zero(); 27];

    for scalar in state.iter_mut() {
        // println!("input is {:?}", scalar.reduce().0);
        // 1. Decomposition
        // Get state value that we are decomposing in non-Montgomery form (comes
        // in Montgomery form by default due to BLS library; but the
        // modular operations can't be done if left like this)
        let mut intermediate = u256(scalar.reduce().0);
        let mut remainder = u256::zero();

        (0..27).for_each(|k| {
            // Reduce intermediate representation
            match k < 26 {
                true => {
                    remainder = intermediate % u256(DECOMPOSITION_S_I[k].0);
                    // Ensure s_i inverses are in Montgomery form, because BLS
                    // scalar multiplication removes a
                    // factor of 2^512
                    let intermediate_scalar: Scalar =
                        Scalar((intermediate - remainder).0) * INVERSES_S_I[k];
                    intermediate = u256(intermediate_scalar.0);
                }
                false => remainder = intermediate,
            };

            // 2. S-box
            nibbles[k] = small_s_box(remainder);
            print!("{}, ", nibbles[k]);
        });
        println!(" ");

        // 3. Composition
        *scalar = compute_whole_representation(nibbles);
    }
}

mod tests {
    use crate::constants::BLS_SCALAR_REAL;

    use super::*;

    #[test]
    fn test_bar() {
        let mut input = [Scalar::one(); 3];
        bar(&mut input);
        let mut breakdown = [u256([15, 0, 0, 0]); 27];
        breakdown[0] = u256([187, 0, 0, 0]);
        let composed = compute_whole_representation(breakdown);
        assert_eq!(input[0], composed);

        // Check whether -5 is dealt with correctly
        let mut input2 = [-Scalar::from(5), -Scalar::from(3), -Scalar::from(1)];
        bar(&mut input2);

        assert_eq!(
            input2[0],
            Scalar([
                18446742991377793276,
                7975156907413507843,
                7849958771744875838,
                2157424182152352530
            ])
        );
        assert_eq!(
            input2[1],
            Scalar([
                18446741084412314296,
                12364043610437066055,
                4990927207653396091,
                3323350968747990091
            ])
        );
        assert_eq!(
            input2[2],
            Scalar([
                18446744060824649731,
                18102478225614246908,
                11073656695919314959,
                6613806504683796440
            ])
        );
    }

    #[test]
    fn test_compute_whole() {
        // Check if -5 is composed correctly
        let expected_breakdown = [
            656, 660, 673, 663, 674, 682, 687, 683, 669, 684, 672, 666, 680,
            662, 686, 668, 661, 678, 692, 686, 689, 660, 690, 687, 683, 674,
            678, 658, 660, 673, 663, 674, 682, 687, 683, 669, 684, 672, 666,
            680, 662, 686, 668, 661, 678, 692, 686, 689, 660, 690, 687, 683,
            674, 678, 658, 660, 673, 663, 674, 682, 687, 683, 669, 684, 672,
            666, 680, 662, 686, 668, 661, 678, 692, 686, 689, 660, 690, 687,
            683, 674, 678,
        ];
        let mut expected = [u256::zero(); 27];
        (0..27).for_each(|k| {
            expected[k] = u256::from(expected_breakdown[k]);
        });
        let composition = compute_whole_representation(expected);
        assert_eq!(composition, -Scalar::from(5));
    }

    #[test]
    fn test_s_box() {
        let six_five_eight = u256::from(658);
        let six_five_nine = u256::from(659);
        let thirty = u256::from(30);
        let six_seventy = u256::from(670);
        assert_eq!(small_s_box(six_five_eight), u256::from(346));
        assert_eq!(small_s_box(six_five_nine), u256::from(659));
        assert_eq!(small_s_box(thirty), u256::from(179));
        assert_eq!(small_s_box(six_seventy), u256::from(670));
    }

    #[test]
    fn test_field_size_decomposition() {
        let size = compute_whole_representation(BLS_SCALAR_REAL);
        assert_eq!(size + Scalar::one(), Scalar::zero());
    }
    #[test]
    fn test_bar_minus_one() {
        let mut input = [-Scalar::one(); 3];
        bar(&mut input);
        let breakdown = [
            u256([660, 0, 0, 0]),
            u256([660, 0, 0, 0]),
            u256([673, 0, 0, 0]),
            u256([663, 0, 0, 0]),
            u256([674, 0, 0, 0]),
            u256([682, 0, 0, 0]),
            u256([687, 0, 0, 0]),
            u256([683, 0, 0, 0]),
            u256([669, 0, 0, 0]),
            u256([684, 0, 0, 0]),
            u256([672, 0, 0, 0]),
            u256([666, 0, 0, 0]),
            u256([680, 0, 0, 0]),
            u256([662, 0, 0, 0]),
            u256([686, 0, 0, 0]),
            u256([668, 0, 0, 0]),
            u256([661, 0, 0, 0]),
            u256([678, 0, 0, 0]),
            u256([692, 0, 0, 0]),
            u256([686, 0, 0, 0]),
            u256([689, 0, 0, 0]),
            u256([660, 0, 0, 0]),
            u256([690, 0, 0, 0]),
            u256([687, 0, 0, 0]),
            u256([683, 0, 0, 0]),
            u256([674, 0, 0, 0]),
            u256([678, 0, 0, 0]),
        ];
        let composed = compute_whole_representation(breakdown);
        assert_eq!(input[0], composed);
    }

    #[test]
    fn test_inverses() {
        for k in 0..27 {
            let product = Scalar(DECOMPOSITION_S_I[k].0) * (INVERSES_S_I[k]);
            assert_eq!(Scalar::from_raw(product.0), Scalar::one());
        }
    }
}
